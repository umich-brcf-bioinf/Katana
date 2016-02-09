"""Basic CIGAR manipulation and querying. """

from __future__ import print_function, absolute_import, division

import itertools
import re

import katana.util as util


class CigarUtil(object):
    _QUERY_CONSUMING = set(list("MIS=X"))
    _REF_CONSUMING = set(list("MDNS=X"))
    _REGEX_CIGAR = re.compile("([0-9]+)([MIDNSHP=X])")
    _REGEX_MATCHING_OP = re.compile("[MX=]")
    _REGEX_NON_HARDCLIP = re.compile("[^H]")
    _REGEX_REF_CONSUMING = re.compile("[MDNS=X]")
    _REGEX_REQUIRED_OPS = re.compile("[MIDN=X]")
    _REGEX_QUERY_CONSUMING = re.compile("[MIS=X]")
    _REGEX_QUERY_NON_CONSUMING = re.compile("[DNP]")

    def __init__(self, reference_start, cigar=None, cigar_profile=None):
        self.reference_start = reference_start
        self.cigar = ""
        self.cigar_profile = ""
        if cigar:
            self.cigar = cigar
        elif cigar_profile:
            self.cigar = self._collapse_cigar_profile(cigar_profile)
        if cigar_profile:
            self.cigar_profile = cigar_profile
        else:
            self.cigar_profile = self._expand_cigar(self.cigar)
        self.query_length = \
                len(self._REGEX_QUERY_CONSUMING.findall(self.cigar_profile))
        self.is_valid = self._REGEX_REQUIRED_OPS.search(self.cigar) is not None

    def __repr__(self):
        return ("{}(reference_start={}, "
                "cigar='{}')").format(self.__class__,
                                      self.reference_start,
                                      self.cigar)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def _is_ref_consuming(self, cigar_op):
        return cigar_op in self._REF_CONSUMING

    def _is_match(self, cigar_op):
        return self._REGEX_MATCHING_OP.search(cigar_op) != None

    def _softclip(self, cigar_profile):
        profile = self._REGEX_QUERY_NON_CONSUMING.sub("", cigar_profile)
        return self._REGEX_NON_HARDCLIP.sub("S", profile)

    def _cut_at_first_match(self, profile):
        cut_index = 0
        first_match = self._REGEX_MATCHING_OP.search(profile)
        if first_match:
            cut_index = first_match.start(0)
        return (profile[0:cut_index], profile[cut_index:])

    def _softclip_to_first_match(self, old_pos, old_profile):
        new_pos = old_pos
        (pre_match, post_match) = self._cut_at_first_match(old_profile)
        new_pos += len(re.findall(self._REGEX_REF_CONSUMING, pre_match))
        pre_match = self._softclip(pre_match)
        return (new_pos, pre_match + post_match)

    def _expand_cigar(self, cigar_string):
        expanded_cigar = []
        for cigar_tuple in self._REGEX_CIGAR.findall(cigar_string):
            expanded_cigar.append(int(cigar_tuple[0]) * cigar_tuple[1])
        return "".join(expanded_cigar)

    @staticmethod
    def _collapse_cigar_profile(profile):
        op_strings = ["".join(g) for _, g in itertools.groupby(profile)]
        new_cigar = [str(len(op)) + op[0] for op in op_strings]
        return "".join(new_cigar)

    def _pos_profiles(self, profile):
        '''Returns a list of tuple of first match index and a list of profiles.
        Each list element is the cigar profile for that reference position; i.e.
        the length of this list is the length of consumed reference bases.
        '''
        pos_profiles = list()
        pos = 0
        pos_profiles.append([])
        first_match_index = -1
        for cigar_op in profile:
            if pos >= len(pos_profiles):
                pos_profiles.append([])
            pos_profiles[pos].append(cigar_op)
            if self._is_match(cigar_op) and first_match_index == -1:
                first_match_index = pos
            if self._is_ref_consuming(cigar_op):
                pos += 1

        profiles = ["".join(x) for x in pos_profiles]

        return (first_match_index, profiles)

    def _partition_cigar(self, region_start, region_end):
        '''Split the profile into 3-tuple: before region, in region,
        and after region; each tuple defines the start coordinate and
        profile fragment.

        ref_start, region start, and region_end are genomic coordinates.

        Regions can extend outside the range of the read but the resulting
        start coordinates are constrained to the range of the read.  For
        example, if region_end overhangs the end of the read, the "after region"
        start coordinate will be "pushed back" to the trailing edge of the
        read.'''
        ref_start = self.reference_start
        profile = self.cigar_profile
        (first_match_index, pos_profile) = self._pos_profiles(profile)

        read_start = ref_start - first_match_index
        read_end = read_start + len(pos_profile)
        constrain = lambda x: max(read_start, min(read_end, x))

        region_start_index = max(region_start - read_start, 0)
        region_end_index = region_end - read_start

        before_profile = "".join(pos_profile[0:region_start_index])
        region_profile="".join(pos_profile[region_start_index:region_end_index])
        after_profile="".join(pos_profile[region_end_index:])

        return (CigarUtil(read_start,
                          cigar_profile=before_profile),
                CigarUtil(constrain(region_start),
                          cigar_profile=region_profile),
                CigarUtil(constrain(region_end),
                          cigar_profile=after_profile))

    def _assert_query_lengths_match(self, new_cigar):
        if self.query_length != new_cigar.query_length:
            msg = ("Old CIGAR query length [{}] ({}) != new CIGAR length"
                    "[{}] ({})").format(self.cigar,
                                        self.query_length,
                                        new_cigar.cigar,
                                        new_cigar.query_length)
            raise util.KatanaException(msg)

    def softclip_target(self, target_start, target_end):
        (pre_target, target, post_target) = self._partition_cigar(target_start,
                                                                  target_end)
        pre_profile = self._softclip(pre_target.cigar_profile)
        (new_pos,
         target_profile) = self._softclip_to_first_match(target.reference_start,
                                                         target.cigar_profile)
        post_profile = self._softclip(post_target.cigar_profile)
        new_profile = pre_profile + target_profile + post_profile
        new_cigar = CigarUtil(new_pos, cigar_profile = new_profile)
        self._assert_query_lengths_match(new_cigar)
        return  new_cigar


class NullCigarUtil(object):
#pylint: disable=unused-argument
    def __init__(self, reference_start):
        self.reference_start = reference_start
        self.cigar = "*"
        self.is_valid = True
        self.query_length = 0

    def softclip_target(self, target_start, target_end):
        return self


def cigar_factory(read):
    if not read.cigarstring or read.cigarstring == "*":
        return NullCigarUtil(read.reference_start)
    else:
        return CigarUtil(read.reference_start, read.cigarstring)

