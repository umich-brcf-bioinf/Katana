""" Super-basic CIGAR manipulation and querying. """

from __future__ import print_function, absolute_import, division
import collections
import itertools

class IndelException(Exception):
    """Flagging cases that we can not process at this time."""
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(IndelException, self).__init__(error_msg)

import re

INDEL_PATTERN = re.compile("[ID]")
REF_CONSUMING_OPS = re.compile("[MDN=XS]")
REF_CONSUMING_OP_LIST = list("MDN=XS")

class CigarEditor(object):
# // The Consume values for each of the CigarOpTypes is as follows:
# //
# //                    Query  Reference
# //  CigarMatch          1        1
# //  CigarInsertion      1        0
# //  CigarDeletion       0        1
# //  CigarSkipped        0        1
# //  CigarSoftClipped    1        0
# //  CigarHardClipped    0        0
# //  CigarPadded         0        0
# //  CigarEqual          1        1
# //  CigarMismatch       1        1
# //  CigarBack           0       -1


    def __init__(self, reference_start, cigar_string):
        self.reference_start = reference_start
        self.cigar = cigar_string
        self.cigar_profile = self._expand_cigar(cigar_string)
        self.reference_end = self._reference_end(self.reference_start,
                                                 self.cigar_profile)

    #TODO: This wrongly assumes the first cigar is a match; let's make this go away?
    def _reference_end(self, reference_start, cigar_profile):
        return reference_start +  len(REF_CONSUMING_OPS.findall(cigar_profile)) - 1

    def indels_in_region(self, start, end):
        return self._contains_indel(self._get_profile_overlap_tuple(start, end)[1])

    @staticmethod
    def _contains_indel(profile):
        return INDEL_PATTERN.search(profile) != None

    def _get_cigar_editor_tuple(self, reference_start, reference_end):
        pos_profiles = collections.defaultdict(str)

        preamble = self.cigar_profile[:self.cigar_profile.find("M")]
        pos = self.reference_start - 1
        for cigar_op in preamble[::-1]:
            pos_profiles[pos] += cigar_op
            if cigar_op in REF_CONSUMING_OP_LIST:
                pos -= 1
        min_pos = pos

        pos = self.reference_start
        for cigar_op in self.cigar_profile[self.cigar_profile.find("M"):]:
            pos_profiles[pos] += cigar_op
            max_pos= pos
            if cigar_op in REF_CONSUMING_OP_LIST:
                pos += 1

        pre_target = []
        target = []
        post_target = []
        target_start = -1
        post_start = -1
        for pos in xrange(min_pos, max_pos + 1):
            ops = pos_profiles[pos]
            if pos < reference_start:
                pre_target.append(ops)
            elif reference_start <= pos < reference_end:
                #TODO: this smells bad
                target_start = target_start if target_start > 0 else pos
                target.append(ops)
            elif pos >= reference_end:
                post_start = post_start if post_start > 0 else pos
                post_target.append(ops)
            else:
                raise IndexError("incorrect position for profile")
        pre_cigar_editor = CigarEditor(min_pos,
                                       self._collapse_cigar_profile("".join(pre_target)))
        target_cigar_editor = CigarEditor(target_start,
                                          self._collapse_cigar_profile("".join(target)))
        post_cigar_editor = CigarEditor(post_start,
                                        self._collapse_cigar_profile("".join(post_target)))
        return (pre_cigar_editor, target_cigar_editor, post_cigar_editor)

    def _get_profile_overlap_tuple(self, reference_start, reference_end):
        (a, b, c) = self._get_cigar_editor_tuple(reference_start, reference_end)
        return (a.cigar_profile, b.cigar_profile, c.cigar_profile)


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

    @staticmethod
    def _expand_cigar(cigar_string):
        pattern = re.compile("([0-9]+)([MIDNSHP=X])")
        expanded_cigar = []
        for cigar_tuple in pattern.findall(cigar_string):
            expanded_cigar.append(int(cigar_tuple[0]) * cigar_tuple[1])
        return "".join(expanded_cigar)

    def softclip_target(self, target_start, target_end):
        (pre_target,
         target,
         post_target) = self._get_cigar_editor_tuple(target_start,
                                                     target_end)
        pre_profile = re.sub("[DNP]", "", pre_target.cigar_profile)
        pre_profile = re.sub("[^H]", "S", pre_profile)
        first_match = target.cigar_profile.find("M")
        new_pos = target.reference_start
        if first_match > -1:
            new_pos += first_match
        target_prematch = re.sub("[DNP]",
                                 "",
                                 target.cigar_profile[0:first_match])
        target_prematch = re.sub("[^H]",
                                 "S",
                                 target_prematch)
        target_profile = target_prematch + target.cigar_profile[first_match:]
        post_profile = re.sub("[DNP]", "", post_target.cigar_profile)
        post_profile = re.sub("[^H]", "S", post_profile)
        new_cigar = self._collapse_cigar_profile(pre_profile \
                                                 + target_profile \
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
                                self._collapse_cigar_profile(new_profile))
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
                                self._collapse_cigar_profile(new_profile))
        return new_cigar

    @staticmethod
    def _collapse_cigar_profile(profile):
        op_strings = ["".join(g) for _, g in itertools.groupby(profile)]
        new_cigar = [str(len(op)) + op[0] for op in op_strings]
        return "".join(new_cigar)

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

