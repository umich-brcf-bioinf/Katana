"""Various classes that handle reads"""
#TODO: elaborate module doc
from __future__ import print_function, absolute_import, division

from collections import defaultdict
import os
import re
import sys

import pysam
from katana.util import KatanaException
from katana import readhandler

PYSAM_ADAPTER = None
class _Pysam8(object):
    _SUPPORTED = re.match(r"^0\.8\.*", pysam.__version__)

    def pysam_index(input_filename):
        pysam.index(input_filename, catch_stdout=False)

    def pysam_sort(input_filename, output_prefix):
        pysam.sort(input_filename, output_prefix, catch_stdout=False)
    
    def pysam_view(input_filename):
        stdout_orig = sys.stdout
        try:
            sys.stdout = sys.__stdout__
            return [x for x in pysam.view(input_filename)]
        finally:
            sys.stdout = stdout_orig

class _Pysam9_10_11_12(object):
    _SUPPORTED = re.match(r"^0\.(9|10|11|12)\.*", pysam.__version__)

    def pysam_index(input_filename):
        pysam.index(input_filename, catch_stdout=False)

    def pysam_sort(input_filename, output_prefix):
        pysam.sort(input_filename, '-o', output_prefix + '.bam', catch_stdout=False)
    
    def pysam_view(input_filename):
        stdout_orig = sys.stdout
        try:
            sys.stdout = sys.__stdout__
            view = pysam.view(input_filename)
            try:
                view = view.decode("utf-8")
            except AttributeError:
                pass # already a string
            return view.strip().split('\n')
        finally:
            sys.stdout = stdout_orig

for pysam_adpater in [_Pysam8, _Pysam9_10_11_12]:
    if pysam_adpater._SUPPORTED:
        PYSAM_ADAPTER=pysam_adpater
        break
else:
    raise KatanaException('Unsupported version of pysam; please review config and install instructions.')

class _BaseReadHandler(object):
    def begin(self):
        pass
    def handle(self, read, read_transformation, mate_transformation):
        pass
    def end(self):
        pass

class AddTagsReadHandler(_BaseReadHandler):
    '''Adds original read values and other explanatory tags.'''
    
    _SANITIZED_REGEX = re.compile('[^A-Za-z0-9-_.]+')

    #pylint: disable=unused-parameter
    def handle(self, read, read_transformation, mate_transformation):
        primer_pair = read_transformation.primer_pair
        sanitized_target_id = self._sanitize(primer_pair.target_id)
        read.set_tag("X0", sanitized_target_id, "Z")
        read.set_tag("X1", read.cigarstring, "Z")
        read.set_tag("X2", read.reference_start, "i")
        read.set_tag("X3", read.reference_end, "i")
        if read_transformation.filters:
            filter_string = ",".join(read_transformation.filters)
            read.set_tag("X4", filter_string, "Z")

    def _sanitize(self, value):
         return re.sub(self._SANITIZED_REGEX, '_', value)

class ExcludeNonMatchedReadHandler(_BaseReadHandler):
    '''Excludes reads from further processing'''
    _STOP_ITERATION_EXCEPTION = StopIteration()
    def __init__(self, log_method):
        self._log_method = log_method
        self._all_exclusions = defaultdict(int)

    def handle(self, read, read_transformation, mate_transformation):
        if mate_transformation.filters:
            read.is_paired = False
        if read_transformation.filters:
            self._all_exclusions[read_transformation.filters] += 1
            raise self._STOP_ITERATION_EXCEPTION

    def end(self):
        for (filters, count) in self._all_exclusions.items():
            msg = "EXCLUDE|{} alignments were excluded because: {}"
            self._log_method(msg, count, ",".join(filters))


class StatsHandler(_BaseReadHandler):
    '''Processes reads and primers connecting PrimerStats and
    PrimerStatsDumper'''
    def __init__(self, primer_stats, primer_stats_dumper):
        self._primer_stats = primer_stats
        self._primer_stats_dumper = primer_stats_dumper

    def handle(self, read, read_transformation, mate_transformation):
        self._primer_stats.add_read_primer(read,
                                           read_transformation.primer_pair)

    def end(self):
        self._primer_stats_dumper.dump(self._primer_stats)


class TransformReadHandler(_BaseReadHandler):
    '''Updates reference_start, cigar string, and mate_reference_start.'''
    def handle(self, read, read_transformation, mate_transformation):
        read.reference_start = read_transformation.reference_start
        read.cigarstring = read_transformation.cigar
        if read.is_paired:
            if mate_transformation.is_unmapped:
                read.mate_cigar = None
            else:
                read.next_reference_start = mate_transformation.reference_start
                read.mate_cigar = mate_transformation.cigar

#TODO: Add PG and CO header lines for tags
class WriteReadHandler(_BaseReadHandler):
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
        readhandler.PYSAM_ADAPTER.pysam_sort(self._tmp_bam_filename, output_root)
        self._log("WRITE_BAM|indexing BAM")
        readhandler.PYSAM_ADAPTER.pysam_index(self._output_bam_filename)
        os.remove(self._tmp_bam_filename)
        self._log("WRITE_BAM|wrote [{}] alignments to [{}]",
                  self._read_count,
                  self._output_bam_filename)
