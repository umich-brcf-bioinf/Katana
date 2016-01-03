"""Softclips primers at edge of aligned reads"""
#TODO: elaborate module doc
from __future__ import print_function, absolute_import, division

import os
import pysam


# I would rather just say pysam.index(...), but since that global is
# added dynamically, Eclipse flags this as a compilation problem. So
# instead we connect directly to the pysam.SamtoolsDispatcher.
PYSAM_INDEX = pysam.SamtoolsDispatcher("index", None).__call__
PYSAM_SORT = pysam.SamtoolsDispatcher("sort", None).__call__


class _BaseReadHandler(object):
    def begin(self):
        pass
    def handle(self, read):
        pass
    def end(self):
        pass


#TODO: Add PG and CO header lines for tags
class _WriteReadHandler(_BaseReadHandler):
    '''Writes reads to a BAM file (ultimately sorting and indexing the BAM).'''
    def __init__(self,
                 input_bam_filename,
                 output_bam_filename,
                 log_method):
        self._input_bam_filename = input_bam_filename
        output_dir = os.path.dirname(output_bam_filename)
        output_basename = "tmp_unsorted_" \
                + os.path.basename(output_bam_filename)
        self._tmp_bam_filename = os.path.join(output_dir, output_basename)
        self._output_bam_filename = output_bam_filename
        self._bamfile = None
        self._log = log_method

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
        self._bamfile.write(read.aligned_segment)

    def end(self):
        self._bamfile.close()
        self._bamfile = None
        output_root = os.path.splitext(self._output_bam_filename)[0]
        self._log("Sorting BAM")
        PYSAM_SORT(self._tmp_bam_filename, output_root)
        self._log("Indexing BAM")
        PYSAM_INDEX(self._output_bam_filename)
        os.remove(self._tmp_bam_filename)


class _TransformReadHandler(_BaseReadHandler):
    '''Updates reference_start and cigar string.'''
    def __init__(self, read_transformations):
        self._read_transformations = read_transformations

    def handle(self, read):
        (dummy,
         new_reference_start,
         new_cigar_string) = self._read_transformations[read.key]
        read.reference_start = new_reference_start
        read.cigarstring = new_cigar_string


class _AddTagsReadHandler(_BaseReadHandler):
    '''Adds original read values and other explanatory tags.'''
    def __init__(self, read_transformations):
        self._read_transformations = read_transformations

    def handle(self, read):
        primer_pair = self._read_transformations[read.key][0]
        read.set_tag("X0", primer_pair.target_id, "Z")
        read.set_tag("X1", read.cigarstring, "Z")
        read.set_tag("X2", read.reference_start, "i")
        read.set_tag("X3", read.reference_end, "i")


#TODO: test
class _UnmatchedPrimerHandler(_BaseReadHandler):
    '''Optionally excludes unmatched primers.'''
    def __init__(self, read_transformations,
                 suppress_unmatched_primer=False,
                 log_method=None):
        self._read_transformations = read_transformations
        self._suppress_unmatched_primer = suppress_unmatched_primer
        self._unmatched_count = 0
        self._log = log_method

    def handle(self, read):
        primer_pair = self._read_transformations[read.key][0]
        if primer_pair.is_unmatched():
            self._unmatched_count += 1
            if self._suppress_unmatched_primer:
                raise StopIteration()

    def end(self):
        if self._suppress_unmatched_primer:
            msg = "excluded from"
        else:
            msg = "preserved in"
        self._log("[{}] reads did not match primers and were {} output",
                  self._unmatched_count,
                  msg)

class _StatsHandler(_BaseReadHandler):
    '''Processes reads and primers connecting PrimerStats and
    PrimerStatsDumper'''
    def __init__(self, read_transformations, primer_stats, primer_stats_dumper):
        self._read_transformations = read_transformations
        self._primer_stats = primer_stats
        self._primer_stats_dumper = primer_stats_dumper

    def handle(self, read):
        primer_pair = self._read_transformations[read.key][0]
        self._primer_stats.add_read_primer(read, primer_pair)

    def end(self):
        self._primer_stats_dumper.dump(self._primer_stats)
