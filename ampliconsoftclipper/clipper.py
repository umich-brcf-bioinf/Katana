#! /usr/bin/env python
"""Softclips primers at edge of aligned reads"""
#TODO: elaborate module doc
from __future__ import print_function, absolute_import, division

from collections import defaultdict
import csv
from datetime import datetime
import itertools
import os
import sys
import traceback

import natsort
import pysam

import ampliconsoftclipper
import ampliconsoftclipper.cigar as cigar
import ampliconsoftclipper.readhandler as readhandler


__VERSION__ = ampliconsoftclipper.__version__


#TODO: make this a logger object that writes to file and console and supports debug and info calls
def _log(msg_format, *args):
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    try:
        print("{}|{}".format(timestamp, msg_format).format(*args),
              file=sys.stderr)
    except IndexError:
        print(args)
    sys.stderr.flush()

#TODO: To avoid circular dependencies, move this and other classes to util module (and reorg Mocks)
#TODO: Capture mapped pairs for each primer
#TODO: Capture overall mapped pairs
#TODO: Capture unmatched primers
class _PrimerStats(object):
    '''Collects simple counts for overall reads and reads per primer.'''
    STAT_KEYS = ["chrom", "target_id", "sense_start", "sense_count",
                 "antisense_count", "sense_percent", "antisense_percent"]
    def __init__(self):
        self.total_read_count = 0
        self._primer_stats = defaultdict(int)


    def _percent(self, primer, is_sense):
        read_count = self._primer_stats[primer, is_sense]
        return int(100 * read_count / self.total_read_count)

    def add_read_primer(self, read, primer_pair):
        self.total_read_count +=1
        stat_key = (primer_pair, read.is_positive_strand)
        self._primer_stats[stat_key] += 1

    @property
    def primer_pairs(self):
        sort_key = lambda x:(x.chrom, x.sense_start, x.target_id)
        primers = set([primer for (primer, _) in self._primer_stats])
        return natsort.natsorted(primers, key=sort_key)

    def stats(self, primer_pair):
        values = [primer_pair.chrom,
                  primer_pair.target_id,
                  primer_pair.sense_start,
                  self._primer_stats[primer_pair, True],
                  self._primer_stats[primer_pair, False],
                  self._percent(primer_pair, True),
                  self._percent(primer_pair, False)]
        return dict(itertools.izip(_PrimerStats.STAT_KEYS, values))


#TODO: Extend log method to emit subset to screen and full set to file
class _PrimerStatsDumper():
    '''Prints PrimerStats object to log.'''
    def __init__(self, log_method=_log):
        self._log_method = log_method

    def dump(self, primer_stats):
        # pylint: disable=star-args
        stat_format = "PRIMER_STATS" + "|{}" * len(primer_stats.STAT_KEYS)
        self._log_method(stat_format.format(*primer_stats.STAT_KEYS))
        for primer_pair in primer_stats.primer_pairs:
            stat_dict = primer_stats.stats(primer_pair)
            stats = [stat_dict[key] for key in primer_stats.STAT_KEYS]
            self._log_method(stat_format.format(*stats))


class _Read(object):
    '''Lightweight wrapper around AlignedSegment'''
    def __init__(self, aligned_segment):
        self.aligned_segment = aligned_segment

    @property
    def is_positive_strand(self):
        return not self.aligned_segment.is_reverse

    @property
    def key(self):
        return (self.aligned_segment.query_name,
                self.is_positive_strand,
                self.aligned_segment.reference_name,
                self.aligned_segment.reference_start)

    @property
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
    def reference_end(self):
        return self.aligned_segment.reference_end

    @property
    def cigarstring(self):
        return self.aligned_segment.cigarstring

    @cigarstring.setter
    def cigarstring(self, value):
        self.aligned_segment.cigarstring = value

    def set_tag(self, tag_name, tag_value, tag_type):
        self.aligned_segment.set_tag(tag_name, tag_value, tag_type)

    @staticmethod
    def iter(aligned_segment_iter):
        for aligned_segment in aligned_segment_iter:
            yield _Read(aligned_segment)

class _NullPrimerPair(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.target_id = "PRIMER_NOT_RECOGNIZED"
        self.chrom = "."
        self.sense_start = 0

    @staticmethod
    def softclip_primers(old_cigar):
        return old_cigar

    @staticmethod
    def is_unmatched():
        return True

class _PrimerPair(object):
    '''Sense start and antisense regions should be start and end (inclusive and
    exclusive respectively), zero-based, positive genomic index.'''
    _all_primers = {}
    NULL_PRIMER_PAIR = _NullPrimerPair()

    def __init__(self,
                 target_id,
                 chrom,
                 sense_primer_region,
                 antisense_primer_region):
        self.target_id = target_id
        self.chrom = chrom
        self.sense_start = sense_primer_region[0]
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
        if read.is_positive_strand:
            return (read.reference_name, read.reference_start, True)
        else:
            return (read.reference_name, read.reference_end, False)

    @staticmethod
    def get_primer_pair(read):
        '''Matches a read to a primer pair based on exact start position of the
        read; returns NullPrimerPair if no match.

        Note this exact approach will disqualify some otherwise valid reads,
        e.g. if a read adapter is not trimmed precisely, those reads will have
        an "early start" and will not be correctly matched with their primer
        pair.
        '''
        try:
            read_key = _PrimerPair._key_for_read(read)
            primer_pair = _PrimerPair._all_primers[read_key]
            return primer_pair
        except KeyError:
            return _PrimerPair.NULL_PRIMER_PAIR

    @staticmethod
    def is_unmatched():
        return False

    def softclip_primers(self, old_cigar):
        return old_cigar.softclip_target(self._query_region_start,
                                         self._query_region_end)

def _build_read_transformations(read_iter):
    read_transformations = {}
    read_count = 0
    for read in read_iter:
        primer_pair = _PrimerPair.get_primer_pair(read)
        old_cigar = cigar.cigar_factory(read)
        new_cigar = primer_pair.softclip_primers(old_cigar)
        read_transformations[read.key] = (primer_pair,
                                          new_cigar.reference_start,
                                          new_cigar.cigar)
        read_count += 1
    _log("Built transforms for [{}] alignments", read_count)
    return read_transformations

def _handle_reads(read_handlers, read_iter, read_transformations):
    for handler in read_handlers:
        handler.begin()
    for read in read_iter:
        read_transformation = read_transformations[read.key]
        try:
            for handler in read_handlers:
                handler.handle(read, read_transformation)
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
        _PrimerPair(row["Customer TargetID"],
                   "chr" + row["Chr"],  #TODO: this prefix seems hackish?
                   (sense_start, sense_end),
                   (antisense_end, antisense_start))

def _build_handlers(input_bam_filename,
                    output_bam_filename,
                    exclude_unmatched_reads):
    stats = readhandler._StatsHandler(_PrimerStats(),
                                      _PrimerStatsDumper(log_method=_log))
    exclude = readhandler._ExcludeNonMatchedReadHandler(log_method=_log)
    tag = readhandler._AddTagsReadHandler()
    transform = readhandler._TransformReadHandler()
    write = readhandler._WriteReadHandler(input_bam_filename,
                                          output_bam_filename,
                                          log_method=_log)
    handlers = [stats, exclude, tag, transform, write]
    if not exclude_unmatched_reads:
        handlers.remove(exclude)
    return handlers

#TODO: test
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
    _log("Read [{}] primer pairs", len(_PrimerPair._all_primers))
    input_bamfile = None
    try:
        _log("Reading alignments from BAM [{}]", input_bam_filename)
        input_bamfile = pysam.AlignmentFile(input_bam_filename,"rb")
        aligned_segment_iter = input_bamfile.fetch()
        read_iter = _Read.iter(aligned_segment_iter)
        _log("Building transformations")
        read_transformations = _build_read_transformations(read_iter)

        _log("Writing transformed alignments to [{}]", output_bam_filename)
        handlers = _build_handlers(input_bam_filename,
                                   output_bam_filename,
                                   True)
        aligned_segment_iter = input_bamfile.fetch()
        read_iter = _Read.iter(aligned_segment_iter)
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
