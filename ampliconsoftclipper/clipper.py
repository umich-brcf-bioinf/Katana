#! /usr/bin/env python
"""Softclips primers at edge of aligned reads"""
#TODO: elaborate module doc
from __future__ import print_function, absolute_import, division

import csv
from datetime import datetime
import os
import pysam
import sys
import traceback
import ampliconsoftclipper
import ampliconsoftclipper.cigar as cigar
import ampliconsoftclipper.readhandler as readhandler
from ampliconsoftclipper.util import PrimerStats, PrimerStatsDumper,\
        PrimerPair, Read

__VERSION__ = ampliconsoftclipper.__version__

#TODO: rename private classes to omit underscore


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
                    exclude_unmatched_reads):
    stats = readhandler.StatsHandler(PrimerStats(),
                                      PrimerStatsDumper(log_method=_log))
    exclude = readhandler.ExcludeNonMatchedReadHandler(log_method=_log)
    tag = readhandler.AddTagsReadHandler()
    transform = readhandler.TransformReadHandler()
    write = readhandler.WriteReadHandler(input_bam_filename,
                                          output_bam_filename,
                                          log_method=_log)
    handlers = [stats, exclude, tag, transform, write]
    if not exclude_unmatched_reads:
        handlers.remove(exclude)
    return handlers

#TODO: test
#TODO: allow suppress/divert unmatched reads
#TODO: deal if input bam missing index
#TODO: deal if input bam regions disjoint with primer regions
#TODO: warn/stop if less than 5% reads transformed
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
        _initialize_primer_pairs(input_primer_manifest)
    _log("Read [{}] primer pairs", len(PrimerPair._all_primers))
    input_bamfile = None
    try:
        _log("Reading alignments from BAM [{}]", input_bam_filename)
        input_bamfile = pysam.AlignmentFile(input_bam_filename,"rb")
        aligned_segment_iter = input_bamfile.fetch()
        read_iter = Read.iter(aligned_segment_iter)
        _log("Building transformations")
        read_transformations = _build_read_transformations(read_iter)

        _log("Writing transformed alignments to [{}]", output_bam_filename)
        handlers = _build_handlers(input_bam_filename,
                                   output_bam_filename,
                                   True)
        aligned_segment_iter = input_bamfile.fetch()
        read_iter = Read.iter(aligned_segment_iter)
        _handle_reads(handlers, read_iter, read_transformations)
    except Exception: #pylint: disable=broad-except
        _log("ERROR: An unexpected error occurred")
        _log(traceback.format_exc())
        exit(1)
    finally:
        if input_bamfile:
            input_bamfile.close()
    _log("Done")


if __name__ == '__main__':
    main(sys.argv)
