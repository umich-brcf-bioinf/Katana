"""Classes common to several modules"""
#TODO: elaborate module doc
from __future__ import print_function, absolute_import, division

from collections import defaultdict

import natsort


try:
    from itertools import izip
except ImportError:  #python3.x
    izip = zip

class KatanaException(Exception):
    """Flagging cases that we can not process at this time."""
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(KatanaException, self).__init__(error_msg)


class PrimerStats(object):
    '''Collects simple counts for overall reads and reads per primer.'''
    STAT_KEYS = ["chrom", "target_id", "sense_start", "sense_count",
                 "antisense_count", "sense_percent"]
    def __init__(self):
        self.total_read_count = 0
        self._primer_stats = defaultdict(int)


    def _sense_percent(self, primer):
        sense_read_count = self._primer_stats[primer, True]
        antisense_read_count = self._primer_stats[primer, False]
        return int(100 * sense_read_count / \
                   (antisense_read_count + sense_read_count))

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
                  self._sense_percent(primer_pair)]
        return dict(izip(PrimerStats.STAT_KEYS, values))


#TODO: Extend log method to emit subset to screen and full set to file
class PrimerStatsDumper(object):
    '''Prints PrimerStats object to log.'''
    def __init__(self, log_method):
        self._log_method = log_method

    def dump(self, primer_stats):
        # pylint: disable=star-args
        stat_format = "PRIMER_STATS" + "|{}" * len(primer_stats.STAT_KEYS)
        self._log_method(stat_format.format(*primer_stats.STAT_KEYS))
        for primer_pair in primer_stats.primer_pairs:
            stat_dict = primer_stats.stats(primer_pair)
            stats = [stat_dict[key] for key in primer_stats.STAT_KEYS]
            self._log_method(stat_format.format(*stats))

#TODO: sort methods
class Read(object):
    '''Lightweight wrapper around AlignedSegment'''
    def __init__(self, aligned_segment, input_bamfile):
        self.aligned_segment = aligned_segment
        self.input_bamfile = input_bamfile

    @property
    def is_positive_strand(self):
        return not self.aligned_segment.is_reverse

    def _key(self, mate=False):
        key = None
        if not mate:
            key= (self.aligned_segment.query_name,
                  self.is_positive_strand,
                  self.input_bamfile.getrname(self.aligned_segment.reference_id),
                  self.aligned_segment.reference_start,
                  self.aligned_segment.is_read1)
        elif self.aligned_segment.is_paired and \
                (not self.aligned_segment.mate_is_unmapped):
            key = (self.aligned_segment.query_name,
                   not self.aligned_segment.mate_is_reverse,
                   self.input_bamfile.getrname(self.aligned_segment.next_reference_id),
                   self.aligned_segment.next_reference_start,
                   not self.aligned_segment.is_read1)
        return key

    @property
    def key(self):
        return self._key(False)

    @property
    def next_reference_start(self):
        return self.aligned_segment.next_reference_start

    @next_reference_start.setter
    def next_reference_start(self, value):
        self.aligned_segment.next_reference_start = value


    @property
    def mate_key(self):
        return self._key(True)

    @property
    def is_unmapped(self):
        return self.aligned_segment.is_unmapped

    @property
    def mate_is_paired(self):
        return self.aligned_segment.is_paired

    @mate_is_paired.setter
    def is_paired(self, value):
        self.aligned_segment.is_paired = value
        if not value:
            self.aligned_segment.mate_is_unmapped = False
            self.aligned_segment.mate_is_reverse = False
            self.aligned_segment.is_proper_pair = False
            self.aligned_segment.is_read1 = False
            self.aligned_segment.is_read2 = False
            self.aligned_segment.next_reference_id = -1
            self.aligned_segment.next_reference_start = -1

    @property
    def query_name(self):
        return self.aligned_segment.query_name

    @property
    def reference_name(self):
        return self.input_bamfile.getrname(self.aligned_segment.reference_id)
#        return self.aligned_segment.rname

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
    def iter(aligned_segment_iter, input_bamfile):
        for aligned_segment in aligned_segment_iter:
            yield Read(aligned_segment, input_bamfile)

    @property
    def mate_cigar(self):
        return self.aligned_segment.get_tag('MC')

    @mate_cigar.setter
    def mate_cigar(self, value):
        self.aligned_segment.set_tag('MC',value)


class _NullPrimerPair(object):
    #pylint: disable=too-few-public-methods
    def __init__(self):
        self.target_id = "PRIMER_NOT_RECOGNIZED"
        self.chrom = "."
        self.sense_start = 0

    @staticmethod
    def softclip_primers(old_cigar):
        return old_cigar

    @property
    def is_unmatched(self):
        #pylint: disable=no-self-use
        return True

class PrimerPair(object):
    '''Sense start and antisense regions should be start and end (inclusive and
    exclusive respectively), zero-based, positive genomic index.'''
    _all_primers = {}
    NULL = _NullPrimerPair()

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
            read_key = PrimerPair._key_for_read(read)
            primer_pair = PrimerPair._all_primers[read_key]
            return primer_pair
        except KeyError:
            return PrimerPair.NULL

    @property
    def is_unmatched(self):
        #pylint: disable=no-self-use
        return False

    def softclip_primers(self, old_cigar):
        return old_cigar.softclip_target(self._query_region_start,
                                         self._query_region_end)

class _NullReadTransformation(object):
    def __init__(self):
        self.is_paired = False
        self.is_unmapped = False
        self.primer_pair = PrimerPair.NULL
        self.reference_start = 0
        self.cigar = ""
        self.is_cigar_valid = False
        self.filters=("NULL",)


class ReadTransformation(object):
    '''Lightweight container for changes to be made to a read.'''

    NULL = _NullReadTransformation()

    def __init__(self,
                 read,
                 primer_pair,
                 new_cigar,
                 filter_builder=None):
        self.is_paired = read.is_paired
        self.is_unmapped = read.is_unmapped
        self.primer_pair = primer_pair
        self.reference_start = new_cigar.reference_start
        self.cigar = new_cigar.cigar
        self.is_cigar_valid = new_cigar.is_valid
        if filter_builder == None:
            self.filters = ()
        else:
            self.filters = tuple(filter_builder(self))

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return ("{}("
                "primer_pair={}, "
                "reference_start={}, "
                "cigar='{}', "
                "is_paired={}, "
                "filters={})").format(self.__class__,
                                      self.primer_pair.target_id,
                                      self.reference_start,
                                      self.cigar,
                                      self.is_paired,
                                      self.filters)

