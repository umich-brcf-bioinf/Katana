""" Super basic CIGAR manipulation and querying. """

from __future__ import print_function, absolute_import, division

import re
from itertools import groupby

class Cigar(object):
    def __init__(self, cigar_string):
        self.cigar = cigar_string
        self.cigar_profile = self._expand_cigar(cigar_string)

    def indels_in_front(self, length):
        indel_pattern = re.compile("[ID]")
        return indel_pattern.search(self.cigar_profile[0:length]) != None

    def indels_in_back(self, length):
        indel_pattern = re.compile("[ID]")
        return indel_pattern.search(self.cigar_profile[-length:]) != None

    def _expand_cigar(self, cigar_string):
        pattern = re.compile("([0-9]+)([MIDNSHP=X])")
        expanded_cigar = []
        for cigar_tuple in pattern.findall(cigar_string):
            expanded_cigar.append(int(cigar_tuple[0]) * cigar_tuple[1])
        return "".join(expanded_cigar)

    
    def softclip_front(self, length=5):
        new_profile = "S" * length + self.cigar_profile[length:]
        new_cigar = Cigar(self._collapse_cigar_profile(new_profile))
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

#5M1I4S
#MMMMMISSSS
#

#     read_consuming_ops = ("M", "I", "S", "=", "X")
#     ref_consuming_ops = ("M", "D", "N", "=", "X")
# 
#     def __init__(self, cigar_string):
#         self.cigar = cigar_string
# 
#     def items(self):
#         if self.cigar == "*":
#             yield (0, None)
#             raise StopIteration
#         cig_iter = groupby(self.cigar, lambda c: c.isdigit())
#         for g, n in cig_iter:
#             yield int("".join(n)), "".join(next(cig_iter)[1])
# 
#     def __str__(self):
#         return self.cigar
# 
#     def __repr__(self):
#         return "Cigar('%s')" % self
# 
#     def __len__(self):
#         """
#         sum of MIS=X ops shall equal the sequence length.
#         """
#         return sum(l for l, op,in self.items() \
#                                if op in Cigar.read_consuming_ops)
# 
#     def reference_length(self):
#         return sum(l for l, op in self.items() \
#                                if op in Cigar.read_consuming_ops)
# 
#     def mask_left(self, n_seq_bases, mask="S"):
#         """
#         Return a new cigar with cigar string where the first `n_seq_bases` are
#         soft-masked unless they are already hard-masked.
#         """
#         cigs = list(self.items())
#         new_cigs = []
# 
#         c, cum_len  = self.cigar, 0
#         for i, (l, op) in enumerate(cigs):
#             if op in Cigar.read_consuming_ops:
#                 cum_len += l
#             if op == "H":
#                 cum_len += l
#                 new_cigs.append(cigs[i])
#             elif cum_len < n_seq_bases:
#                 new_cigs.append(cigs[i])
#             else:
#                 # the current cigar element is split by the masking.
#                 right_extra = cum_len - n_seq_bases
#                 new_cigs.append((l - right_extra, 'S'))
#                 if right_extra != 0:
#                     new_cigs.append((right_extra, cigs[i][1]))
#             if cum_len >= n_seq_bases: break
#         else:
#             pass
# 
#         new_cigs[:i] = [(l, op if op in "HS" else "S") for l, op in
#                 new_cigs[:i]]
#         new_cigs.extend(cigs[i + 1:])
#         return Cigar(Cigar.string_from_elements(new_cigs))
# 
# 
#     @classmethod
#     def string_from_elements(self, elements):
#         return "".join("%i%s" % (l, op) for l, op in elements if l !=0)
# 
# 
#     def mask_right(self, n_seq_bases, mask="S"):
#         """
#         Return a new cigar with cigar string where the last `n_seq_bases` are
#         soft-masked unless they are already hard-masked.
#         """
#         return Cigar(Cigar(self._reverse_cigar()).mask_left(n_seq_bases, mask)._reverse_cigar())
# 
# 
#     def _reverse_cigar(self):
#         return Cigar.string_from_elements(list(self.items())[::-1])
