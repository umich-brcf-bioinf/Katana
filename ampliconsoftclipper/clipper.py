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

def _log(msg_format, *args):
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    try:
        print("{}|{}".format(timestamp, msg_format).format(*args),
              file=sys.stderr)
    except IndexError:
        print(args)
    sys.stderr.flush()


class Read(object):
    def __init__(self, aligned_segment):
        self.aligned_segment = aligned_segment

    #TODO: Simplify using is_reverse?
    def _is_positive_strand(self):
        return self.aligned_segment.flag & _SAM_FLAG_NEGATIVE_STRAND == 0

    @property
    def key(self):
        return (self.aligned_segment.query_name,
                _is_positive_strand(self.aligned_segment),
                self.aligned_segment.reference_name,
                self.aligned_segment.reference_start)

    def mate_key(self):
        if not self.aligned_segment.is_paired \
                or self.aligned_segment.mate_is_unmapped:
            return None
        return (self.aligned_segment.query_name,
                not self.aligned_segment.mate_is_reverse,
                self.aligned_segment.next_reference_name,
                self.aligned_segment.next_reference_start)

    @property
    def reference_name(self):
        return self.aligned_segment.reference_name

    @property
    def reference_start(self):
        return self.aligned_segment.reference_start

    @reference_start.setter
    def reference_start(self, value):
        self.aligned_segment.reference_start = value

    @property
    def cigarstring(self):
        return self.aligned_segment.cigarstring

    @cigarstring.setter
    def cigarstring(self, value):
        self.aligned_segment.cigarstring = value

    def set_tag(self, tag_name, tag_value, tag_type):
        self.aligned_segment.set_tag(tag_name, tag_value, tag_type)

class _NullPrimerPair(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.target_id = "PRIMER_NOT_RECOGNIZED"
        self.chrom = "."
        self.sense_start = 0

    @staticmethod
    def softclip_primers(old_cigar):
        return old_cigar


class _PrimerPair(object):
    '''Sense start and antisense regions should be start and end (exclusive),
    zero-based, positive genomic index.'''

    _all_primers = {}
    NULL_PRIMER_PAIR = _NullPrimerPair()

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
        sense_key = (chrom, sense_start, True)
        antisense_key = (chrom, antisense_start, False)
        self._all_primers[sense_key] = self
        self._all_primers[antisense_key] = self

    @staticmethod
    def _key_for_read(read):
        if _is_positive_strand(read):
            return (read.reference_name, read.reference_start, True)
        else:
            return (read.reference_name, read.reference_end, False)

    @staticmethod
    def get_primer_pair(read):
        try:
            read_key = _PrimerPair._key_for_read(read)
            primer_pair = _PrimerPair._all_primers[read_key]
            return primer_pair
        except KeyError:
            return _PrimerPair.NULL_PRIMER_PAIR

    def softclip_primers(self, old_cigar):
        return old_cigar.softclip_target(self._query_region_start,
                                         self._query_region_end)


class _ReadHandler(object):
    def begin(self):
        pass
    def handle(self, read):
        pass
    def end(self):
        pass


#TODO: Add PG and CO header lines for tags
class _WriteReadHandler(_ReadHandler):
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
        _log("Sorting BAM")
        PYSAM_SORT(self._tmp_bam_filename, output_root)
        _log("Indexing BAM")
        PYSAM_INDEX(self._output_bam_filename)
        os.remove(self._tmp_bam_filename)


class _TransformReadHandler(_ReadHandler):
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


class _TagReadHandler(_ReadHandler):
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

#TODO: Capture paired
#TODO: Also emit stats file
class _CaptureStatsHandler(_ReadHandler):
    def __init__(self,
                 read_transformations,
                 key_builder=_build_read_transform_key):
        self._read_transformations = read_transformations
        self._total_read_count = len(read_transformations)
        self._key_builder = key_builder
        self._primer_stats = defaultdict(int)

    def _percent(self, primer, is_sense):
        read_count = self._primer_stats[primer, is_sense]
        return int(100 * read_count / self._total_read_count)

    #TODO: test
    def handle(self, read):
        read_key = self._key_builder(read)
        primer_pair = self._read_transformations[read_key][0]
        stat_key = (primer_pair, _is_positive_strand(read))
        self._primer_stats[stat_key] += 1

    #TODO: test
    def end(self):
        header = "chrom|target_id|sense_start|sense|antisense|sense|antisense"
        stat_format = "SUMMARY|{}|{}|{}|{}|{}|{}%|{}%"
        primers = set([primer for (primer, _) in self._primer_stats])
        _log(stat_format, *header.split("|"))
        sort_key = lambda x:(x.chrom, x.target_id, x.sense_start)
        for primer in natsort.natsorted(primers, key=sort_key):
            _log(stat_format,
                primer.chrom,
                primer.target_id,
                primer.sense_start,
                self._primer_stats[primer, True],
                self._primer_stats[primer, False],
                self._percent(primer, True),
                self._percent(primer, False))


def _build_read_transformations(read_iter):
    read_transformations = {}
    read_count = 0
    for read in read_iter:
        primer_pair = _PrimerPair.get_primer_pair(read)
        old_cigar = cigar.cigar_factory(read)
        new_cigar = primer_pair.softclip_primers(old_cigar)
        key = _build_read_transform_key(read)
        read_transformations[key] = (primer_pair,
                                     new_cigar.reference_start,
                                     new_cigar.cigar)
        read_count += 1
    _log("Read [{}] alignments", read_count)
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
        _PrimerPair(row["Customer TargetID"],
                   "chr" + row["Chr"],  #TODO: this prefix seems hackish?
                   (sense_start, sense_end),
                   (antisense_end, antisense_start))

#TODO: test
#TODO: guard if input bam missing index
#TODO: guard if input bam regions disjoint with primer regions
#TODO: guard if less than 5% reads transformed
#TODO: argparse
def main(command_line_args=None):

    if not command_line_args:
        command_line_args = sys.argv
    if len(command_line_args) != 4:
        usage = "usage: {0} [thunderbolt_manifest] [input_bam] [output_bam]"
        print(usage.format(os.path.basename(command_line_args[0])),
              file=sys.stderr)
        sys.exit(1)
    (input_primer_manifest_filename,
     input_bam_filename,
     output_bam_filename) = command_line_args[1:]

    #pylint: disable=no-member
    _log("Reading primer pairs from [{}]", input_primer_manifest_filename)
    with open(input_primer_manifest_filename, "r") as input_primer_manifest:
        _read_primer_pairs(input_primer_manifest)
    _log("Read [{}] primer pairs", len(_PrimerPair._all_primers))
    input_bamfile = None
    try:
        _log("Reading alignments from BAM [{}]", input_bam_filename)
        input_bamfile = pysam.AlignmentFile(input_bam_filename,"rb")
        read_iter = input_bamfile.fetch()
        read_transformations = _build_read_transformations(read_iter)

        _log("Writing transformed alignments to [{}]", output_bam_filename)
        handlers = [_CaptureStatsHandler(read_transformations),
                    _TagReadHandler(read_transformations),
                    _TransformReadHandler(read_transformations),
                    _WriteReadHandler(input_bam_filename, output_bam_filename),]
        read_iter = input_bamfile.fetch()
        _handle_reads(read_iter, handlers)
    except Exception as exception:
        _log("ERROR: An unexpected error occurred")
        _log(traceback.format_exc())
        exit(1)
    finally:
        if input_bamfile:
            input_bamfile.close()
    _log("Done")


if __name__ == '__main__':
    main(sys.argv)
