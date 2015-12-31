#! /usr/bin/env python
""" Basic parser for the Rhim Thunderbolts manifest file. """

from __future__ import print_function, absolute_import, division

import ampliconsoftclipper
import ampliconsoftclipper.cigar as cigar
from collections import defaultdict
import csv
from datetime import datetime
import natsort
import sys
import traceback
import os
import pysam

__VERSION__ = ampliconsoftclipper.__version__

_SAM_FLAG_NEGATIVE_STRAND = 0x10

# I would rather just say pysam.index(...), but since that global is
# added dynamically, Eclipse flags this as a compilation problem. So
# instead we connect directly to the pysam.SamtoolsDispatcher.
PYSAM_INDEX = pysam.SamtoolsDispatcher("index", None).__call__
PYSAM_SORT = pysam.SamtoolsDispatcher("sort", None).__call__

def _is_positive_strand(read):
    return read.flag & _SAM_FLAG_NEGATIVE_STRAND == 0

def _build_read_transform_key(read):
    return (read.query_name,
            _is_positive_strand(read),
            read.reference_name,
            read.reference_start)

def log(msg_format, *args):
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    print("{}|{}".format(timestamp, msg_format).format(*args),
          file=sys.stderr)
    sys.stderr.flush()

class NullPrimerPair(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.target_id = "PRIMER_NOT_RECOGNIZED"
        self.chrom = "."
        self.sense_start = 0

    @staticmethod
    def softclip_primers(old_cigar):
        return old_cigar


class PrimerPair(object):
    '''Sense start and antisense regions should be start and end (exclusive),
    zero-based, positive genomic index.'''

    _all_primers = {}
    NULL_PRIMER_PAIR = NullPrimerPair()

    def __init__(self,
                 target_id,
                 chrom,
                 sense_primer_region,
                 antisense_primer_region):
        self.target_id = target_id
        self.chrom = chrom #TODO: test
        self.sense_start = sense_primer_region[0] #TODO: test
        self._query_region_start = sense_primer_region[1]
        self._query_region_end =  antisense_primer_region[0]
        self._add_primer(chrom,
                         sense_primer_region[0],
                         antisense_primer_region[1])

    def _add_primer(self, chrom, sense_start, antisense_start):
        sense_key = (chrom, sense_start, '+')
        antisense_key = (chrom, antisense_start, '-')
        self._all_primers[sense_key] = self
        self._all_primers[antisense_key] = self

    @staticmethod
    def _key_for_read(read):
        if _is_positive_strand(read):
            return (read.reference_name, read.reference_start, '+')
        else:
            return (read.reference_name, read.reference_end, '-')

    @staticmethod
    def get_primer_pair(read):
        try:
            read_key = PrimerPair._key_for_read(read)
            primer_pair = PrimerPair._all_primers[read_key]
            return primer_pair
        except KeyError:
            return PrimerPair.NULL_PRIMER_PAIR

    def softclip_primers(self, old_cigar):
        return old_cigar.softclip_target(self._query_region_start,
                                         self._query_region_end)

class ReadHandler(object):
    def begin(self):
        pass
    def handle(self, read):
        pass
    def end(self):
        pass


#TODO: Add PG and CO header lines for tags
class WriteReadHandler(ReadHandler):
    def __init__(self, input_bam_filename, output_bam_filename):
        self._input_bam_filename = input_bam_filename
        output_dir = os.path.dirname(output_bam_filename)
        output_basename = "tmp_unsorted_" \
                + os.path.basename(output_bam_filename)
        self._tmp_bam_filename = os.path.join(output_dir, output_basename)
        self._output_bam_filename = output_bam_filename
        self._bamfile = None

    def begin(self):
        #pylint: disable=no-member
        input_bam = None
        try:
            input_bam = pysam.AlignmentFile(self._input_bam_filename, "rb")
            self._bamfile = pysam.AlignmentFile(self._tmp_bam_filename,
                                                "wb",
                                                template=input_bam)
        finally:
            if input_bam:
                input_bam.close()

    def handle(self, read):
        self._bamfile.write(read)

    #TODO: test sort/index
    def end(self):
        self._bamfile.close()
        self._bamfile = None
        output_root = os.path.splitext(self._output_bam_filename)[0]
        log("Sorting BAM")
        PYSAM_SORT(self._tmp_bam_filename, output_root)
        log("Indexing BAM")
        PYSAM_INDEX(self._output_bam_filename)
        os.remove(self._tmp_bam_filename)


class TransformReadHandler(ReadHandler):
    def __init__(self,
                 read_transformations,
                 key_builder=_build_read_transform_key):
        self._read_transformations = read_transformations
        self._key_builder = key_builder

    #TODO: test
    def handle(self, read):
        key = self._key_builder(read)
        (dummy,
         new_reference_start,
         new_cigar_string) = self._read_transformations[key]
        read.reference_start = new_reference_start
        read.cigarstring = new_cigar_string


class TagReadHandler(ReadHandler):
    def __init__(self,
                 read_transformations,
                 key_builder=_build_read_transform_key):
        self._read_transformations = read_transformations
        self._key_builder = key_builder

    #TODO: test
    def handle(self, read):
        key = self._key_builder(read)
        primer_pair = self._read_transformations[key][0]
        read.set_tag("X0", primer_pair.target_id, "Z")
        read.set_tag("X1", read.cigarstring, "Z")
        read.set_tag("X2", read.reference_start, "i")
        read.set_tag("X3", read.reference_end, "i")

class CaptureStatsHandler(ReadHandler):
    def __init__(self,
                 read_transformations,
                 key_builder=_build_read_transform_key):
        self._read_transformations = read_transformations
        self._total_read_count = len(read_transformations)
        self._key_builder = key_builder
        self._primer_stats = defaultdict(int)

    #TODO: test
    def handle(self, read):
        read_key = self._key_builder(read)
        primer_pair = self._read_transformations[read_key][0]
        stat_key = (primer_pair, _is_positive_strand(read))
        self._primer_stats[stat_key] += 1

    #TODO: test
    def end(self):
        stat_format = "{}|{}|{}|{}|{}|{}%|{}%"
        primers = set([primer for (primer, _) in self._primer_stats])
        log(stat_format, "chrom", "target_id", "sense_start", "sense", "antisense",
            "sense", "antisense")
        sort_key = lambda x:(x.chrom, x.target_id, x.sense_start)
        for primer in natsort.natsorted(primers, key=sort_key):  #TODO: sort
            log(stat_format,
                primer.chrom,
                primer.target_id,
                primer.sense_start,
                self._primer_stats[primer, True],
                self._primer_stats[primer, False],
                int(100 * self._primer_stats[primer, True] / self._total_read_count),
                int(100 * self._primer_stats[primer, False] / self._total_read_count))


#TODO: reset mate pos

def _build_read_transformations(read_iter):
    read_transformations = {}
    read_count = 0
    for read in read_iter:
        primer_pair = PrimerPair.get_primer_pair(read)
        old_cigar = cigar.cigar_factory(read)
        new_cigar = primer_pair.softclip_primers(old_cigar)
        key = _build_read_transform_key(read)
        read_transformations[key] = (primer_pair,
                                     new_cigar.reference_start,
                                     new_cigar.cigar)
        read_count += 1
    log("Read [{}] alignments", read_count)
    return read_transformations

def _handle_reads(read_iter, read_handlers):
    for handler in read_handlers:
        handler.begin()
    for read in read_iter:
        for handler in read_handlers:
            handler.handle(read)
    for handler in read_handlers:
        handler.end()

#TODO: test
def _read_primer_pairs(base_reader):
    dict_reader = csv.DictReader(base_reader, delimiter='\t')
    for row in dict_reader:
        sense_start = int(row["Sense Start"]) - 1
        sense_end = sense_start + len(row["Sense Sequence"])
        antisense_start = int(row["Antisense Start"])
        antisense_end = antisense_start - len(row["Antisense Sequence"])
        PrimerPair(row["Customer TargetID"],
                   "chr" + row["Chr"],  #TODO: this prefix seems hackish?
                   (sense_start, sense_end),
                   (antisense_end, antisense_start))

#TODO: test
def main2(input_primer_manifest_filename,
          input_bam_filename,
          output_bam_filename):
    #pylint: disable=no-member
    log("Read primer pairs from [{}]", input_primer_manifest_filename)
    with open(input_primer_manifest_filename, "r") as input_primer_manifest:
        _read_primer_pairs(input_primer_manifest)
    log("Read [{}] primer pairs", len(PrimerPair._all_primers))
    input_bamfile = None
    try:
        log("Reading alignments from BAM [{}]", input_bam_filename)
        input_bamfile = pysam.AlignmentFile(input_bam_filename,"rb")
        read_iter = input_bamfile.fetch()
        read_transformations = _build_read_transformations(read_iter)

        log("Writing transformed alignments to [{}]", output_bam_filename)
        handlers = [CaptureStatsHandler(read_transformations),
                    TagReadHandler(read_transformations),
                    TransformReadHandler(read_transformations),
                    WriteReadHandler(input_bam_filename, output_bam_filename),]
        read_iter = input_bamfile.fetch()
        _handle_reads(read_iter, handlers)
    except Exception as exception:
        log("ERROR: An unexpected error occurred")
        log(traceback.format_exc())
        exit(1)
    finally:
        if input_bamfile:
            input_bamfile.close()
    log("Done")

#TODO: Maybe pass a hash?
class PrimerPairRecord(object):
    """Represents a primer pair derived from a MiSeq manifest. """
    def __init__(self, id, target_id, primer_set, chrom, original_start, original_end, converted_end, converted_start, genome_build, sense_start,
                 antisense_start, sense_sequence, antisense_sequence, sense_sequence_tailed_illumina, antisense_sequence_tailed_illumina):
        self.id = id
        self.target_id = target_id
        self.primer_set = primer_set
        #TODO: careful with this pattern
        self.chrom = 'chr'+chrom
        self.orginal_start = int(original_start)
        self.original_end = int(original_end)
        self.converted_start = int(converted_start)
        self.converted_end = int(converted_end)
        self.genome_build = genome_build
        self.sense_start = int(sense_start)
        self.antisense_start = int(antisense_start)
        self.sense_sequence = sense_sequence
        self.antisense_sequence = antisense_sequence
        self.sense_sequence_tailed_illumina = sense_sequence_tailed_illumina
        self.antisense_sequence_tailed_illumina = antisense_sequence_tailed_illumina

#    @property
#     def to_fasta(self):
#         fa_header1 = ">{0}_{1}_{2}_F".format(self.id, self.target_id, self.sense_sequence)
#         outline1 = "{0}\n{1}\n".format(fa_header1, self.sense_sequence.upper())
#         fa_header2 = ">{0}_{1}_{2}_R".format(self.id, self.target_id, self.sense_sequence)
#         outline1 = "{0}\n{1}\n".format(fa_header1, self.sense_sequence.upper())
# 


#TODO: CSVReader?
def parse_thunderbolts_manifest(filename):
    """ Parses a Rhim manifest file to pull out primer data. Returns a dict. """
    #primer_l = []
    primer_pos_d = {}
    datafile = file(filename)
    header = datafile.readline()
    for line in datafile.readlines():
        bits = line.strip().split("\t")
        ppr = PrimerPairRecord(*bits[0:15])
        primer_pos_d[(ppr.chrom, ppr.sense_start-1, '+')] = ppr
        primer_pos_d[(ppr.chrom, ppr.antisense_start, '-')] = ppr
        #primer_l.append(ppr)
    return primer_pos_d


# def get_region(primer_pair_record):
#     return (primer_pair_record.chrom,
#             primer_pair_record.sense_start,
#             primer_pair_record.antisense_start)

# def get_reads_with_coords(alignment_file, chrom, start, end):
#     read_l = []
#     read_total = 0
#     forw_total = 0
#     rev_total = 0
#     fmatch_total = 0
#     rmatch_total = 0
#     print("coords: {0}\t{1}\n".format(start, end))
#     for read in alignment_file.fetch(chrom, start, end):
#         read_total += 1
#         if (read.flag & 16 == 0):
#             forw_total += 1
#             if (read.reference_start == start):
#                 fmatch_total += 1
#                 read_l.append(read)
#         else:
#             rev_total += 1
#             if (read.reference_end == end):
#                 rmatch_total += 1
#                 read_l.append(read)
#     print("{0}\t{1}\t{2}\t{3}\t{4}".format(read_total, forw_total, rev_total, fmatch_total, rmatch_total))
#     return read_l


def create_bam_from_reads_and_template(read_list, template_bam, new_fn):
    out_bam = pysam.AlignmentFile(new_fn, "wh", template=template_bam)
    for read in read_list:
        out_bam.write(read)
    out_bam.close()


def get_key_from_read(read):
    if (read.flag & 16 == 0):
        return (read.reference_name, read.reference_start, '+')
    else:
        return (read.reference_name, read.reference_end, '-')


#TODO: argparse
def main():
    if len(sys.argv) != 4:
        print("usage: {0} [thunderbolt_manifest] [input_bam] [output_bam]".format(os.path.basename(sys.argv[0])))
        sys.exit()

    main2(sys.argv[1], sys.argv[2], sys.argv[3])
# 
#     bam_fn = sys.argv[2]
#     outfolder = sys.argv[3]
# 
#     thunderbolt_manifest = sys.argv[1]
#     primer_pos_d = parse_thunderbolts_manifest(thunderbolt_manifest)
# 
#     bamfile = pysam.AlignmentFile(bam_fn, 'rb')
#     outfile = outfolder+"/"+"subset.sam"
#     out_bam = pysam.AlignmentFile(outfile, "wh", template=bamfile)
#     count = 0
#     hit_count = 0
#     sense_count = 0
#     indels_in_primer_region = 0
#     for read in bamfile.fetch():
#         cigar_util = cigar.cigar_factory(read.reference_start, read.cigarstring)
#         if (cigar_util.is_null()):
#             continue
#         count += 1
#         read_key = get_key_from_read(read)
#         try:
#             primer_pair_record = primer_pos_d[read_key]
#             read.set_tag("X0", primer_pair_record.target_id, "Z")
#             hit_count += 1
#             if read_key[2] == '+':
#                 sense_clip_len = len(primer_pair_record.sense_sequence)
# ##                 if cigar.edge_indels("front", sense_clip_len):
# ##                     indels_in_primer_region += 1
# ##                     continue
#                 new_cigar = cigar_util.softclip(primer_pair_record.sense_start + sense_clip_len,
#                                                 primer_pair_record.antisense_start - antisense_clip_len)
#                 read.cigarstring = new_cigar.cigar
#                 read.reference_start = new_cigar.pos
#             elif read_key[2] == '-':
#                 antisense_clip_len = len(primer_pair_record.antisense_sequence)
#                 if cigar.edge_indels("back", antisense_clip_len):
#                     indels_in_primer_region += 1
#                     continue
#                 new_cigar = cigar.softclip_back(antisense_clip_len)
#                 read.cigarstring = new_cigar.cigar
#                 read.reference_start = new_cigar.pos
# 
#             out_bam.write(read)
#         except KeyError:
#             #This read not an amplicon
#             pass
#     out_bam.close()
#     print("total: {}\nhit: {}\nsense: {}\nindel in primer region: {}".format(count, hit_count, sense_count, indels_in_primer_region))
#     sys.exit()
# 
# 
# #     selected_reads = get_reads_with_coords(bamfile, chrom, start-1, end)
# #     print("Total reads: {0}".format(len(selected_reads)))
#     
#     create_bam_from_reads_and_template(selected_reads, bamfile, outfile)
# 
#     print("done.")

if __name__ == '__main__':
    main()

