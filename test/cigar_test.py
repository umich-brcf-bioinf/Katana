#pylint: disable=invalid-name, too-few-public-methods
from __future__ import print_function, absolute_import
import unittest
from ampliconsoftclipper.cigar import Cigar

class CigarTestCase(unittest.TestCase):
    def test_init(self):
        cigar = Cigar("10M")
        self.assertEquals("10M", cigar.cigar)

    def test_expand_cigar(self):
        cigar = Cigar("5M1I4S")
        self.assertEquals("MMMMMISSSS", cigar.cigar_profile)

    def test_expand_cigar_longer(self):
        cigar = Cigar("10M1I1D10H")
        self.assertEquals("MMMMMMMMMMIDHHHHHHHHHH", cigar.cigar_profile)

    def test_indels_in_front_falseIfNoIndels(self):
        cigar = Cigar("10M")
        self.assertEquals(False, cigar.indels_in_front(length=10))

    def test_indels_in_front_falseIfNoIndels(self):
        cigar = Cigar("5M5S")
        self.assertEquals(False, cigar.indels_in_front(length=10))

    def test_indels_in_front_trueIfInsert(self):
        cigar = Cigar("5M1I4M")
        self.assertEquals(True, cigar.indels_in_front(length=10))

    def test_indels_in_front_falseIfInsertOutsideRange(self):
        cigar = Cigar("5M1I4M")
        self.assertEquals(False, cigar.indels_in_front(length=5))

    def test_indels_in_front_trueIfDelete(self):
        cigar = Cigar("5M1D4M")
        self.assertEquals(True, cigar.indels_in_front(length=10))

    def test_indels_in_front_falseIfDeleteOutsideRange(self):
        cigar = Cigar("5M1D4M")
        self.assertEquals(False, cigar.indels_in_front(length=5))

    def test_indels_in_back_falseIfNoIndels(self):
        cigar = Cigar("10M")
        self.assertEquals(False, cigar.indels_in_back(length=10))

    def test_indels_in_back_falseIfNoIndels(self):
        cigar = Cigar("5M5S")
        self.assertEquals(False, cigar.indels_in_back(length=10))

    def test_indels_in_back_trueIfInsert(self):
        cigar = Cigar("4M1I5M")
        self.assertEquals(True, cigar.indels_in_back(length=10))

    def test_indels_in_back_falseIfInsertOutsideRange(self):
        cigar = Cigar("4M1I5M")
        self.assertEquals(False, cigar.indels_in_back(length=5))

    def test_indels_in_back_trueIfDelete(self):
        cigar = Cigar("4M1D5M")
        self.assertEquals(True, cigar.indels_in_back(length=10))

    def test_indels_in_back_falseIfDeleteOutsideRange(self):
        cigar = Cigar("4M1D5M")
        self.assertEquals(False, cigar.indels_in_back(length=5))

    def test_softclip_front(self):
        cigar = Cigar("10M")
        self.assertEquals("5S5M", cigar.softclip_front(length=5).cigar)

    #test softclip_front with indel
    
    #test cigar = *
    #test malformed cigar
