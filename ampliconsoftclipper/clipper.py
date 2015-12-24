#! /usr/bin/env python
""" Basic parser for the Rhim Thunderbolts manifest file. """

from __future__ import print_function, absolute_import, division

import ampliconsoftclipper.cigar as cigar
#from ampliconsoftclipper import __version__
import sys
import os
import pysam

class PrimerPair(object):
    '''Sense start and antisense start should be zero-based, positive genomic
    index.'''

    _all_primers = {}
    SAM_FLAG_NEGATIVE_STRAND = 0x10
    TAG_NAME = "X0"
    TAG_TYPE = "Z"
    TAG_VALUE_NOT_FOUND = "PRIMER_NOT_FOUND"

    def __init__(self,
                 target_id,
                 chrom,
                 sense_start,
                 antisense_start,
                 sense_sequence,
                 antisense_sequence):
        self._target_id = target_id
        self._query_region_start = sense_start + len(sense_sequence)
        self._query_region_end =  antisense_start - len(antisense_sequence)
        self._add_primer(chrom, sense_start, antisense_start)

    def _add_primer(self, chrom, sense_start, antisense_start):
        sense_key = (chrom, sense_start, '+')
        antisense_key = (chrom, antisense_start, '-')
        self._all_primers[sense_key] = self
        self._all_primers[antisense_key] = self

    @staticmethod
    def _is_positive_strand(read):
        return read.flag & PrimerPair.SAM_FLAG_NEGATIVE_STRAND == 0

    @staticmethod
    def _key_for_read(read):
        if PrimerPair._is_positive_strand(read):
            return (read.reference_name, read.reference_start, '+')
        else:
            return (read.reference_name, read.reference_end, '-')

    #TODO: extend to handle KeyError
    @staticmethod
    def softclip_read(read, cigar_util):
        try:
            read_key = PrimerPair._key_for_read(read)
            primer_pair = PrimerPair._all_primers[read_key]
            new_cigar = cigar_util.softclip_target(primer_pair._query_region_start,
                                                   primer_pair._query_region_end)
            read.reference_start = new_cigar.reference_start
            read.cigarstring = new_cigar.cigar
            tag_value = "{}|{}|{}|{}".format(read_key[0],
                                             str(read_key[1]),
                                             read_key[2],
                                             primer_pair._target_id)
            read.set_tag(PrimerPair.TAG_NAME, tag_value, PrimerPair.TAG_TYPE)
        except KeyError:
            read.set_tag(PrimerPair.TAG_NAME,
                         PrimerPair.TAG_VALUE_NOT_FOUND,
                         PrimerPair.TAG_TYPE)
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



def main():
    if (len(sys.argv) != 4):
        print("usage: {0} [thunderbolts manifest] [bam] [outfolder]".format(os.path.basename(sys.argv[0])))
        sys.exit()

    bam_fn = sys.argv[2]
    outfolder = sys.argv[3]

    thunderbolt_manifest = sys.argv[1]
    primer_pos_d = parse_thunderbolts_manifest(thunderbolt_manifest)
    #print(primer_pos_d)

    bamfile = pysam.AlignmentFile(bam_fn, 'rb')
    outfile = outfolder+"/"+"subset.sam"
    out_bam = pysam.AlignmentFile(outfile, "wh", template=bamfile)
    count = 0
    hit_count = 0
    sense_count = 0
    indels_in_primer_region = 0
    for read in bamfile.fetch():
        cigar_util = cigar.cigar_factory(read.reference_start, read.cigarstring)
        if (cigar_util.is_null()):
            continue
        count += 1
        read_key = get_key_from_read(read)
        try:
            primer_pair_record = primer_pos_d[read_key]
            read.set_tag("X0", primer_pair_record.target_id, "Z")
            hit_count += 1
            if read_key[2] == '+':
                sense_clip_len = len(primer_pair_record.sense_sequence)
#                 if cigar.edge_indels("front", sense_clip_len):
#                     indels_in_primer_region += 1
#                     continue
                new_cigar = cigar_util.softclip(primer_pair_record.sense_start + sense_clip_len,
                                                primer_pair_record.antisense_start - antisense_clip_len)
                read.cigarstring = new_cigar.cigar
                read.reference_start = new_cigar.pos
            elif read_key[2] == '-':
                antisense_clip_len = len(primer_pair_record.antisense_sequence)
                if cigar.edge_indels("back", antisense_clip_len):
                    indels_in_primer_region += 1
                    continue
                new_cigar = cigar.softclip_back(antisense_clip_len)
                read.cigarstring = new_cigar.cigar
                read.reference_start = new_cigar.pos

            out_bam.write(read)
        except KeyError:
            #This read not an amplicon
            pass
    out_bam.close()
    print("total: {}\nhit: {}\nsense: {}\nindel in primer region: {}".format(count, hit_count, sense_count, indels_in_primer_region))
    sys.exit()


#     selected_reads = get_reads_with_coords(bamfile, chrom, start-1, end)
#     print("Total reads: {0}".format(len(selected_reads)))
    
    create_bam_from_reads_and_template(selected_reads, bamfile, outfile)

    print("done.")

if __name__ == '__main__':
    main()

