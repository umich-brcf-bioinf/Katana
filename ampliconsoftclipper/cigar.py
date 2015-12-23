""" Super-basic CIGAR manipulation and querying. """

from __future__ import print_function, absolute_import, division
import itertools
import re


class IndelException(Exception):
    """Flagging cases that we can not process at this time."""
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(IndelException, self).__init__(error_msg)


class CigarUtil(object):

    _REGEX_QUERY_NON_CONSUMING = re.compile("[DNP]") #Preserve H?
    _REGEX_NON_HARDCLIP = re.compile("[^H]")
    _REGEX_MATCHING_OP = re.compile("[MX=]")
    _REGEX_CIGAR = re.compile("([0-9]+)([MIDNSHP=X])")
    _QUERY_CONSUMING = set(list("MIS=X"))
    _REF_CONSUMING = set(list("MDNS=X")) #S is debatable, but works better
    _REGEX_REF_CONSUMING = re.compile("[MDNS=X]")

    def __init__(self):
        pass

#     def in_query_consuming(self, cigar_op):
#         return cigar_op in self._QUERY_CONSUMING

    def _is_ref_consuming(self, cigar_op):
        return cigar_op in self._REF_CONSUMING

    def _is_match(self, cigar_op):
        return self._REGEX_MATCHING_OP.search(cigar_op) != None

    #TODO move softclip_target inside and make this private
    def softclip(self, cigar_profile):
        profile = self._REGEX_QUERY_NON_CONSUMING.sub("", cigar_profile)
        return self._REGEX_NON_HARDCLIP.sub("S", profile)

    def _cut_at_first_match(self, profile):
        cut_index = 0
        first_match = self._REGEX_MATCHING_OP.search(profile)
        if first_match:
            cut_index = first_match.start(0)
        return (profile[0:cut_index], profile[cut_index:])

    #TODO move softclip_target inside and make this private
    def softclip_to_first_match(self, old_pos, old_profile):
        new_pos = old_pos
        (pre_match, post_match) = self._cut_at_first_match(old_profile)
        new_pos += len(re.findall(self._REGEX_REF_CONSUMING, pre_match))
        pre_match = CIGAR_UTIL.softclip(pre_match)
        return (new_pos, pre_match + post_match)

    #TODO: transition CigarEditor into this and make this private
    def expand_cigar(self, cigar_string):
        expanded_cigar = []
        for cigar_tuple in self._REGEX_CIGAR.findall(cigar_string):
            expanded_cigar.append(int(cigar_tuple[0]) * cigar_tuple[1])
        return "".join(expanded_cigar)

    @staticmethod
    def collapse_cigar_profile(profile):
        op_strings = ["".join(g) for _, g in itertools.groupby(profile)]
        new_cigar = [str(len(op)) + op[0] for op in op_strings]
        return "".join(new_cigar)

    def _pos_profiles(self, profile):
        '''Returns a list of tuple of first match index and a list of profiles.
        The list elements match the genomic coordinates of the profile; i.e.
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

    #TODO: transition CigarEditor to CigarUtil and make this private
    def partition_cigar(self, ref_start, profile, region_start, region_end):
        '''Split the profile into 3-tuple: before region, in region,
        and after region; each tuple defines the start coordinate and
        profile fragment.

        ref_start, region start, and region_end are genomic coordinates.

        Regions can extend outside the range of the read but the resulting
        start coordinates are constrained to the range of the read.  For
        example, if region_end overhangs the end of the read, the "after region"
        start coordinate will be "pushed back" to the trailing edge of the
        read.'''
        (first_match_index, pos_profile) = self._pos_profiles(profile)

        read_start = ref_start - first_match_index
        read_end = read_start + len(pos_profile)
        constrain = lambda x: max(read_start, min(read_end, x))

        region_start_index = region_start - read_start
        region_end_index = region_end - read_start

        before = (read_start,
                  "".join(pos_profile[0:region_start_index]))
        region = (constrain(region_start),
                  "".join(pos_profile[region_start_index:region_end_index]))
        after = (constrain(region_end),
                 "".join(pos_profile[region_end_index:]))
        return (before, region, after)


CIGAR_UTIL = CigarUtil()
INDEL_PATTERN = re.compile("[ID]")

class CigarEditor(object):
    def __init__(self, reference_start, cigar_string):
        self.reference_start = reference_start
        self.cigar = cigar_string
        self.cigar_profile = CIGAR_UTIL.expand_cigar(cigar_string)


    def indels_in_region(self, start, end):
        (_, region_profile) = CIGAR_UTIL.partition_cigar(self.reference_start,
                                                         self.cigar_profile,
                                                         start,
                                                         end)[1]
        return self._contains_indel(region_profile)



    @staticmethod
    def _contains_indel(profile):
        return INDEL_PATTERN.search(profile) != None


    def edge_indels(self, front_or_back, length):
        if length > len(self.cigar_profile):
            msg = "Requested edge indel length [{}] > sequence length [{}]"
            raise ValueError(msg.format(length, len(self.cigar_profile)))
        if front_or_back == 'front':
            cigar_sequence = self.cigar_profile[0:length]
        elif front_or_back == 'back':
            cigar_sequence = self.cigar_profile[-length:]
        else:
            raise ValueError("Parameter 'front_or_back' must "
                             "be 'front' or 'back'")
        return INDEL_PATTERN.search(cigar_sequence) != None

    #TODO normalize names e.g. pre_target <> before_region, pos<>ref_start
    def softclip_target(self, target_start, target_end):
#         (pre_target,
#          target,
#          post_target) = self._partition_cigar(target_start,
#                                                      target_end)
        #TODO: Yuk
        ((_dummy, pre_target),
         (target_pos, target),
         (_dummy, post_target)) = CIGAR_UTIL.partition_cigar(self.reference_start,
                                                   self.cigar_profile,
                                                   target_start,
                                                    target_end)


        pre_profile = CIGAR_UTIL.softclip(pre_target)

        (new_pos,
         new_target_profile) = CIGAR_UTIL.softclip_to_first_match(target_pos,
                                                                 target)

        post_profile = CIGAR_UTIL.softclip(post_target)
        new_cigar = CIGAR_UTIL.collapse_cigar_profile(pre_profile \
                                                 + new_target_profile \
                                                 + post_profile)
        return CigarEditor(new_pos, new_cigar)

    def softclip_front(self, length):
        if length > len(self.cigar_profile):
            msg = "Requested softclip length [{}] > sequence length [{}]"
            raise ValueError(msg.format(length, len(self.cigar_profile)))
        if self.edge_indels('front', length):
            msg = "Indel found in first [{}] bases of cigar [{}]"
            raise IndelException(msg.format(length, self.cigar))
        front = self.cigar_profile[:length]
        back = self.cigar_profile[length:]
        new_front = re.sub("[^H]", "S", front)
        new_profile = new_front + back
        new_first_match = new_profile.find("M")
        old_first_match = self.cigar_profile.find("M")
        new_pos = self.reference_start + (new_first_match - old_first_match)
        new_cigar = CigarEditor(new_pos,
                                CIGAR_UTIL.collapse_cigar_profile(new_profile))
        return new_cigar

    def softclip_back(self, length):
        if length > len(self.cigar_profile):
            msg = "Requested softclip length [{}] > sequence length [{}]"
            raise ValueError(msg.format(length, len(self.cigar_profile)))
        if self.edge_indels('back', length):
            msg = "Indel found in last [{}] bases of cigar [{}]"
            raise IndelException(msg.format(length, self.cigar))
        front = self.cigar_profile[:-length]
        back = self.cigar_profile[-length:]
        new_back = re.sub("[^H]", "S", back)
        new_profile = front + new_back
        new_cigar = CigarEditor(self.reference_start,
                                CIGAR_UTIL.collapse_cigar_profile(new_profile))
        return new_cigar

    @staticmethod
    def is_null():
        return False


class NullCigarEditor(object):
#pylint: disable=unused-argument
    def __init__(self, pos):
        self.cigar = "*"
        self.reference_start = pos

    @staticmethod
    def edge_indels(front_or_back, length):
        return False

    @staticmethod
    def indels_in_back(length):
        return False

    def softclip_front(self, length):
        return self

    def softclip_back(self, length):
        return self

    @staticmethod
    def is_null():
        return True

def cigar_editor_factory(pos, cigar_string):
    if not cigar_string or cigar_string == "*":
        return NullCigarEditor(pos)
    else:
        return CigarEditor(pos, cigar_string)

