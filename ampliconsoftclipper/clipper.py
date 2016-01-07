#! /usr/bin/env python
"""Softclips primers at edge of aligned reads based on primer locations; emits
new BAM file with clipped reads optionally excluding alignments that did not
match primer locations."""
#TODO: Add TODO file
#TODO: elaborate module doc
#TODO: Add README.rst
#TODO: Add ability to parse Illumna primer files
#TODO: Emit primer bed file
#TODO: Add to travis
#TODO: Test in Py3
#TODO: Add setup.py
#TODO: Add to PyPI

##   Copyright 2014 Bioinformatics Core, University of Michigan
##
##   Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##
##       http://www.apache.org/licenses/LICENSE-2.0
##
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

from __future__ import print_function, absolute_import, division

import argparse
import csv
from datetime import datetime
import pysam
import sys
import traceback

import ampliconsoftclipper
import ampliconsoftclipper.cigar as cigar
import ampliconsoftclipper.readhandler as readhandler
from ampliconsoftclipper.util import PrimerStats, PrimerStatsDumper,\
        PrimerPair, Read

__version__ = ampliconsoftclipper.__version__

DESCRIPTION=\
'''Match each alignment in input BAM to primer, softclipping the primer region.

This helps identify variants that land in the primer region because sequences
which start or end with the primer region are only measuring the efficacy of
the sample prep and the presence of those alignments tend to overwhelm true
variants which occur in that region measured by overlapping primers. The output
is conceptually similar to clipping the primers from the FASTQ reads but
preserves the primers during alignment to improve alignment quality.
'''

class _ClipperUsageError(Exception):
    """Raised for malformed command or invalid arguments."""
    def __init__(self, msg, *args):
        super(_ClipperUsageError, self).__init__(msg, *args)


class _ClipperArgumentParser(argparse.ArgumentParser):
    """Argument parser that raises UsageError instead of exiting."""
    #pylint: disable=too-few-public-methods
    def error(self, message):
        '''Suppress default exit behavior'''
        raise _ClipperUsageError(message)


#TODO: make this a logger object that writes to file and console and supports debug and info calls
def _log(msg_format, *args):
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    try:
        print("{}|{}".format(timestamp, msg_format).format(*args),
              file=sys.stderr)
    except IndexError:
        print(args)
    sys.stderr.flush()

def _build_read_transformations(read_iter):
    read_transformations = {}
    read_count = 0
    for read in read_iter:
        primer_pair = PrimerPair.get_primer_pair(read)
        old_cigar = cigar.cigar_factory(read)
        new_cigar = primer_pair.softclip_primers(old_cigar)
        read_transformations[read.key] = (primer_pair,
                                          new_cigar.reference_start,
                                          new_cigar.cigar)
        read_count += 1
    _log("Built transforms for [{}] alignments", read_count)
    return read_transformations

def _handle_reads(read_handlers, read_iter, read_transformations):
    null_transformation = (PrimerPair.NULL_PRIMER_PAIR, 0, "")
    for handler in read_handlers:
        handler.begin()
    for read in read_iter:
        read_transformation = read_transformations[read.key]
        mate_transformation = read_transformations.get(read.mate_key,
                                                       null_transformation)
        try:
            for handler in read_handlers:
                handler.handle(read, read_transformation, mate_transformation)
        except StopIteration:
            pass
    for handler in read_handlers:
        handler.end()

def _initialize_primer_pairs(base_reader):
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

def _build_handlers(input_bam_filename,
                    output_bam_filename,
                    include_unmatched_reads):
    stats = readhandler.StatsHandler(PrimerStats(),
                                      PrimerStatsDumper(log_method=_log))
    exclude = readhandler.ExcludeNonMatchedReadHandler(log_method=_log)
    tag = readhandler.AddTagsReadHandler()
    transform = readhandler.TransformReadHandler()
    write = readhandler.WriteReadHandler(input_bam_filename,
                                          output_bam_filename,
                                          log_method=_log)
    handlers = [stats, exclude, tag, transform, write]
    if include_unmatched_reads:
        handlers.remove(exclude)
    return handlers

#TODO: test
def _parse_command_line_args(arguments):
    parser = _ClipperArgumentParser( \
        formatter_class=argparse.RawTextHelpFormatter,
        usage="clipper primer_manifest input_bam output_bam",
        description=(DESCRIPTION))

    parser.add_argument("-V",
                        "--version",
                        action='version',
                        version=__version__)
    parser.add_argument('primer_manifest',
                        help="path to primer manifest (tab-separated text)")
    parser.add_argument('input_bam',
                        help="path to input BAM")
    parser.add_argument('output_bam',
                        help="path to output BAM")
    parser.add_argument("--include_unmatched",
                        action="store_true",
                        help=("Preserve all incoming alignments (even if they"
                              "cannot be matched with primers"))
    args = parser.parse_args(arguments)
    return args

#TODO: test
#TODO: deal if input bam missing index
#TODO: deal if input bam regions disjoint with primer regions
#TODO: warn/stop if less than 5% reads transformed
def main(command_line_args=None):
    '''Clipper entry point.'''
    try:
        if not command_line_args:
            command_line_args = sys.argv
        args = _parse_command_line_args(command_line_args[1:])

        _log("Reading primer pairs from [{}]", args.primer_manifest)
        with open(args.primer_manifest, "r") as input_primer_manifest:
            _initialize_primer_pairs(input_primer_manifest)
        _log("Read [{}] primer pairs", len(PrimerPair._all_primers))
        input_bamfile = None

        _log("Building transformations from BAM [{}]", args.input_bam)
        #pylint: disable=no-member
        input_bamfile = pysam.AlignmentFile(args.input_bam,"rb")
        aligned_segment_iter = input_bamfile.fetch()
        read_iter = Read.iter(aligned_segment_iter)
        read_transformations = _build_read_transformations(read_iter)

        _log("Writing transformed alignments to [{}]", args.output_bam)
        handlers = _build_handlers(args.input_bam,
                                   args.output_bam,
                                   args.include_unmatched)
        aligned_segment_iter = input_bamfile.fetch()
        read_iter = Read.iter(aligned_segment_iter)
        _handle_reads(handlers, read_iter, read_transformations)

        _log("Done")
    except _ClipperUsageError as usage_error:
        message = "clipper usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'clipper --help'.", file=sys.stderr)
        sys.exit(1)
    except Exception: #pylint: disable=broad-except
        _log("ERROR: An unexpected error occurred")
        _log(traceback.format_exc())
        exit(1)
    finally:
        try:
            if input_bamfile:
                input_bamfile.close()
        except NameError:
            pass


if __name__ == '__main__':
    main(sys.argv)
