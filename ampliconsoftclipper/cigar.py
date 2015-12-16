""" Super-basic CIGAR manipulation and querying. """

from __future__ import print_function, absolute_import, division

class IndelException(Exception):
    """Flagging cases that we can not process at this time."""
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(IndelException, self).__init__(error_msg)

import re
from itertools import groupby

#CigarEditor?
class Cigar(object):
    INDEL_PATTERN = re.compile("[ID]")

    def __init__(self, pos, cigar_string):
        self.pos = pos
        self.cigar = cigar_string
        self.cigar_profile = self._expand_cigar(cigar_string)

    def edge_indels(self, front_or_back, length):
        if length > len(self.cigar_profile):
            msg = "Requested edge indel length [{}] > sequence length [{}]"
            raise ValueError(msg.format(length, len(self.cigar_profile)))
        if front_or_back == 'front':
            cigar_sequence = self.cigar_profile[0:length]
        elif front_or_back == 'back':
           cigar_sequence = self.cigar_profile[-length:]
        else:
            raise ValueError("Parameter 'front_or_back' must be 'front' or 'back'")
        return self.INDEL_PATTERN.search(cigar_sequence) != None

    def _expand_cigar(self, cigar_string):
        try:
            pattern = re.compile("([0-9]+)([MIDNSHP=X])")
            expanded_cigar = []
            for cigar_tuple in pattern.findall(cigar_string):
                expanded_cigar.append(int(cigar_tuple[0]) * cigar_tuple[1])
            return "".join(expanded_cigar)
        except TypeError as e:
            print (e)
            print(cigar_string)
            raise

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
        new_first_M = new_profile.find("M")
        old_first_M = self.cigar_profile.find("M")
        new_pos = self.pos + (new_first_M - old_first_M)
        new_cigar = Cigar(new_pos, self._collapse_cigar_profile(new_profile))
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
        new_cigar = Cigar(self.pos, self._collapse_cigar_profile(new_profile))
        return new_cigar

    def _collapse_cigar_profile(self, profile):
        op_strings = []
        current_op = profile[0]
        op_string = current_op
        for op in profile[1:]:
            if op != current_op:
                op_strings.append(op_string)
                op_string = op
                current_op = op
            else:
                op_string += op
        op_strings.append(op_string)
        new_cigar = ""
        for op_string in op_strings:
            new_cigar += str(len(op_string)) + op_string[0]
        return new_cigar

    def is_null(self):
        return False


class NullCigar(object):
    def __init__(self, pos):
        self.cigar = "*"
        self.pos = pos

    def edge_indels(self, front_or_back, length):
        return False

    def indels_in_back(self, length):
        return False

    def softclip_front(self, length):
        return self

    def softclip_back(self, length):
        return self
        
    def is_null(self):
        return True

def cigar_factory(pos, cigar_string):
    if not cigar_string or cigar_string == "*":
        return NullCigar(pos)
    else:
        return Cigar(pos, cigar_string)

