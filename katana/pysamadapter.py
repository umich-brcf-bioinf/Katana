from __future__ import print_function, absolute_import, division

import re
import sys

import pysam

from katana.util import KatanaException

PYSAM_ADAPTER = None

class _Pysam(object):
    @staticmethod
    def alignment_file(bam_filename, mode='rb', template=None, header=None):
        return pysam.AlignmentFile(bam_filename, mode, template=template, header=header)

    @staticmethod
    def index(input_filename):
        pysam.index(input_filename, catch_stdout=False)

    @staticmethod
    def aligned_segment():
        return pysam.AlignedSegment()

class _Pysam8(_Pysam):
    _SUPPORTED = re.match(r"^0\.8\.*", pysam.__version__)
    
    @staticmethod
    def sort(input_filename, output_prefix):
        pysam.sort(input_filename, output_prefix, catch_stdout=False)
    
    @staticmethod
    def view(input_filename):
        stdout_orig = sys.stdout
        try:
            sys.stdout = sys.__stdout__
            return [x for x in pysam.view(input_filename)]
        finally:
            sys.stdout = stdout_orig

class _Pysam9_10_11_12(_Pysam):
    _SUPPORTED = re.match(r"^0\.(9|10|11|12)\.*", pysam.__version__)

    @staticmethod
    def sort(input_filename, output_prefix):
        pysam.sort(input_filename, '-o', output_prefix + '.bam', catch_stdout=False)
    
    @staticmethod
    def view(input_filename):
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
        PYSAM_ADAPTER=pysam_adpater()
        break
else:
    raise KatanaException('Unsupported version of pysam; please review config and install instructions.')
