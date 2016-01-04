"""Various classes that handle reads"""
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
    def handle(self, read, read_transformation, mate_transformation):
        pass
    def end(self):
        pass


class _AddTagsReadHandler(_BaseReadHandler):
    '''Adds original read values and other explanatory tags.'''
    def handle(self, read, read_transformation, mate_transformation):
        primer_pair = read_transformation[0]
        read.set_tag("X0", primer_pair.target_id, "Z")
        read.set_tag("X1", read.cigarstring, "Z")
        read.set_tag("X2", read.reference_start, "i")
        read.set_tag("X3", read.reference_end, "i")

class _ExcludeNonMatchedReadHandler(_BaseReadHandler):
    '''Excludes reads from further processing'''
    STOP_ITERATION_EXCEPTION = StopIteration()
    def __init__(self, log_method):
        self._unmatched_count = 0
        self._broken_pair_count = 0
        self._neither_match_count = 0
        self._log_method = log_method

    def handle(self, read, read_transformation, mate_transformation):
        primer_pair = read_transformation[0]
        mate_primer_pair = mate_transformation[0]
        if primer_pair.is_unmatched:
            self._unmatched_count += 1
            raise self.STOP_ITERATION_EXCEPTION
        if read.mate_is_mapped and mate_primer_pair.is_unmatched:
            self._broken_pair_count += 1
            read.mate_is_mapped = False

    def end(self):
        msg = ("EXCLUDE|[{}] alignments did not match a primer and will be "
               "excluded from the output")
        self._log_method(msg, self._unmatched_count)
        msg = ("EXCLUDE|[{}] alignment pairs were broken because their mates "
               "did not match a primer")
        self._log_method(msg, self._broken_pair_count)

class _StatsHandler(_BaseReadHandler):
    '''Processes reads and primers connecting PrimerStats and
    PrimerStatsDumper'''
    def __init__(self, primer_stats, primer_stats_dumper):
        self._primer_stats = primer_stats
        self._primer_stats_dumper = primer_stats_dumper

    def handle(self, read, read_transformation, mate_transformation):
        primer_pair = read_transformation[0]
        self._primer_stats.add_read_primer(read, primer_pair)

    def end(self):
        self._primer_stats_dumper.dump(self._primer_stats)


class _TransformReadHandler(_BaseReadHandler):
    '''Updates reference_start, cigar string, and mate_reference_start.'''
    def handle(self, read, read_transformation, mate_transformation):
        (new_reference_start, new_cigar_string) = read_transformation[1:]
        read.reference_start = new_reference_start
        read.cigarstring = new_cigar_string
        (mate_primer_pair, mate_reference_start) = mate_transformation[0:2]
        if not mate_primer_pair.is_unmatched:
            read.mate_reference_start = mate_reference_start


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
        self._read_count = 0

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

    def handle(self, read, read_transformation, mate_transformation):
        self._read_count += 1
        self._bamfile.write(read.aligned_segment)

    def end(self):
        self._bamfile.close()
        self._bamfile = None
        output_root = os.path.splitext(self._output_bam_filename)[0]
        self._log("WRITE_BAM|sorting BAM")
        PYSAM_SORT(self._tmp_bam_filename, output_root)
        self._log("WRITE_BAM|indexing BAM")
        PYSAM_INDEX(self._output_bam_filename)
        os.remove(self._tmp_bam_filename)
        self._log("WRITE_BAM|wrote [{}] alignments to [{}]",
                  self._read_count,
                  self._output_bam_filename)
