#! /usr/bin/env python
"""Softclips primers at edge of aligned reads based on primer locations; emits
new BAM file with clipped reads optionally excluding alignments that did not
match primer locations."""
#TODO: elaborate module doc
#TODO: see TODO.rst
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
import resource
import sys
import time
import traceback

import katana
import katana.cigar as cigar
import katana.pysamadapter as pysamadapter
import katana.readhandler as readhandler
from katana.util import KatanaException, PrimerStats, PrimerStatsDumper, \
    PrimerPair, Read, ReadTransformation

__version__ = katana.__version__

DESCRIPTION=\
'''Match each alignment in input BAM to primer, softclipping the primer region.

Katana matches each read to its corresponding primer pair based on start
position of the read. Katana then soft-clips the primer region from the edge of
the read sequence, rescuing the signal of true variants measured by overlapping
amplicons. The output is conceptually similar to hard-clipping the primers from
the original FASTQ reads based on sequence identity but with the advantage that
retaining the primers during alignment improves alignment quality.
'''

class _KatanaUsageError(Exception):
    """Raised for malformed command or invalid arguments."""
    def __init__(self, msg, *args):
        super(_KatanaUsageError, self).__init__(msg, *args)


class _KatanaArgumentParser(argparse.ArgumentParser):
    """Argument parser that raises UsageError instead of exiting."""
    #pylint: disable=too-few-public-methods
    def error(self, message):
        '''Suppress default exit behavior'''
        raise _KatanaUsageError(message)


#TODO: make this a logger object that writes to file and console and
# supports debug and info calls
def _log(msg_format, *args):
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    try:
        print("{}|{}".format(timestamp, msg_format).format(*args),
              file=sys.stderr)
    except IndexError:
        print(args)
    sys.stderr.flush()

def _filter_builder(read_transformation):
    filters = []
    if read_transformation.is_unmapped:
        filters.append("UNMAPPED_ALIGNMENT")
    else:
        if read_transformation.primer_pair.is_unmatched:
            filters.append("UNMATCHED_PRIMER_PAIR")
        if not read_transformation.is_cigar_valid:
            filters.append("INVALID_CIGAR")
    return filters

#TODO: refactor to expedite testing (e.g. clipped_cigar_provider,
#  cached_clipped_cigar_provider)
def _build_read_transformations(read_iter, filter_builder):
    read_transformations = {}
    read_count = 0
    cigar_cache={}
    for read in read_iter:
        try:
            read_count += 1
            primer_pair = PrimerPair.get_primer_pair(read)
            key = (primer_pair, read.reference_start, read.cigarstring)
            if not cigar_cache.get(key):
                old_cigar = cigar.cigar_factory(read)
                new_cigar = primer_pair.softclip_primers(old_cigar)
                cigar_cache[key] = new_cigar
            new_cigar = cigar_cache[key]
            transform = ReadTransformation(read,
                                           primer_pair,
                                           new_cigar,
                                           filter_builder)
            read_transformations[read.key] = transform
        except Exception as exception:
            msg = "Problem with read {} [line {}] and primer pair {}: {}"
            raise KatanaException(msg.format(read.query_name,
                                             read_count,
                                             primer_pair.target_id,
                                             exception))
    _log("Built transforms for [{}] alignments", read_count)
    return read_transformations


def _handle_reads(read_handlers, read_iter, read_transformations):
    for handler in read_handlers:
        handler.begin()
    for read in read_iter:
        read_transformation = read_transformations[read.key]
        mate_transformation = read_transformations.get(read.mate_key,
                                                       ReadTransformation.NULL)
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
    handlers = [stats, tag, transform, exclude, write]
    if include_unmatched_reads:
        handlers.remove(exclude)
    return handlers

def _parse_command_line_args(arguments):
    parser = _KatanaArgumentParser( \
        formatter_class=argparse.RawTextHelpFormatter,
        usage="katana primer_manifest input_bam output_bam",
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
    parser.add_argument("--preserve_all_alignments",
                        action="store_true",
                        help=("Preserve all incoming alignments (even if they "
                              "are unmapped, cannot be matched with primers, "
                              "result in invalid CIGARs, etc.)"))
    args = parser.parse_args(arguments)
    return args

def _peak_memory():
    peak_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_memory_mb = peak_memory/1024
    if sys.platform == 'darwin':
        peak_memory_mb /= 1024
    return int(peak_memory_mb)


#TODO: test
#TODO: check input files exist
def main(command_line_args=None):
    '''Katana entry point.'''
    try:
        start_time = time.time()
        if not command_line_args:
            command_line_args = sys.argv
        args = _parse_command_line_args(command_line_args[1:])

        _log("Reading primer pairs from [{}]", args.primer_manifest)
        with open(args.primer_manifest, "r") as input_primer_manifest:
            _initialize_primer_pairs(input_primer_manifest)
        _log("Read [{}] primer pairs", len(PrimerPair._all_primers))
        input_bamfile = None

        _log("Building transformations from BAM [{}]", args.input_bam)
        input_bamfile = pysamadapter.PYSAM_ADAPTER.alignment_file(args.input_bam)
        aligned_segment_iter = input_bamfile.fetch()
        read_iter = Read.iter(aligned_segment_iter, input_bamfile)
        read_transformations = _build_read_transformations(read_iter,
                                                            _filter_builder)

        _log("Writing transformed alignments to [{}]", args.output_bam)
        handlers = _build_handlers(args.input_bam,
                                   args.output_bam,
                                   args.preserve_all_alignments)
        aligned_segment_iter = input_bamfile.fetch()
        read_iter = Read.iter(aligned_segment_iter, input_bamfile)
        _handle_reads(handlers, read_iter, read_transformations)

        elapsed_time = int(time.time() - start_time)
        _log("Done ({} seconds, {}mb peak memory)",
             elapsed_time,
             _peak_memory())
    except _KatanaUsageError as usage_error:
        message = "katana usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'katana --help'.", file=sys.stderr)
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
    #import cProfile
    #cProfile.run('main()')
    main(sys.argv)
