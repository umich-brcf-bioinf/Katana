#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods

from __future__ import print_function, absolute_import
import unittest
from ampliconsoftclipper import clipper


class MicroMock(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class MockRead(MicroMock):
    def __init__(self, **kwargs):
        self._tags={}
        super(MockRead, self).__init__(**kwargs)

    def set_tag(self, tag_name, tag_value, tag_type):
        self._tags[tag_name] = "{}:{}:{}".format(tag_name, tag_type, tag_value)

class MockCigarUtil(MicroMock):
    def __init__(self, **kwargs):
        self._softclip_targets_calls=[]
        self._softclip_target_return = None
        super(MockCigarUtil, self).__init__(**kwargs)

    def softclip_target(self, target_start, target_end):
        self._softclip_targets_calls.append((target_start, target_end))
        return self._softclip_target_return


class PrimerPairTestCase(unittest.TestCase):
    def test_init(self):
        primer_pair = clipper.PrimerPair(target_id="target_1",
                                         chrom="chr42",
                                         sense_start=100,
                                         antisense_start=150,
                                         sense_sequence="AAACCCGGGT",
                                         antisense_sequence="TTTCCCGGGA")
        self.assertEquals("target_1", primer_pair._target_id)
        self.assertEquals(110, primer_pair._query_region_start)
        self.assertEquals(140, primer_pair._query_region_end)

    def test_all_primers(self):
        primer_pair1 = clipper.PrimerPair(target_id="target_1",
                                          chrom="chr42",
                                          sense_start=100,
                                          antisense_start=150,
                                          sense_sequence="AACCGGTT",
                                          antisense_sequence="CTTTA")
        primer_pair2 = clipper.PrimerPair(target_id="target_1",
                                          chrom="chr42",
                                          sense_start=200,
                                          antisense_start=250,
                                          sense_sequence="AACCGGTT",
                                          antisense_sequence="CTTTA")
        actual_primers = clipper.PrimerPair._all_primers
        self.assertEquals(primer_pair1, actual_primers[('chr42', 100, '+')])
        self.assertEquals(primer_pair1, actual_primers[('chr42', 150, '-')])
        self.assertEquals(primer_pair2, actual_primers[('chr42', 200, '+')])
        self.assertEquals(primer_pair2, actual_primers[('chr42', 250, '-')])
        self.assertEquals(4, len(actual_primers))

    def test_softclip_read_positivePrimer(self):
        clipper.PrimerPair._all_primers = {}
        clipper.PrimerPair(target_id="target_1",
                           chrom="chr42",
                           sense_start=100,
                           antisense_start=150,
                           sense_sequence="AACCGGTT",
                           antisense_sequence="CTTTA")
        read = MockRead(reference_name="chr42",
                        flag=0,
                        reference_start=100,
                        reference_end=200,
                        cigarstring="75M")
        clipped_cigar_util = MockCigarUtil(reference_start=42, cigar="75X")
        mock_cigar_util = MockCigarUtil(_softclip_target_return=clipped_cigar_util)

        clipper.PrimerPair.softclip_read(read, mock_cigar_util)

        self.assertEquals(42, read.__dict__["reference_start"])
        self.assertEquals("75X", read.__dict__["cigarstring"])
        self.assertEquals({"X0" : "X0:Z:chr42|100|+|target_1"}, read._tags)

    def test_softclip_read_negativePrimer(self):
        clipper.PrimerPair._all_primers = {}
        clipper.PrimerPair(target_id="target_1",
                           chrom="chr42",
                           sense_start=100,
                           antisense_start=150,
                           sense_sequence="AACCGGTT",
                           antisense_sequence="CTTTA")
        read = MockRead(reference_name="chr42",
                        flag=16,
                        reference_start=0,
                        reference_end=150,
                        cigarstring="75M")
        clipped_cigar_util = MockCigarUtil(reference_start=42, cigar="75X")
        mock_cigar_util = MockCigarUtil(_softclip_target_return=clipped_cigar_util)

        clipper.PrimerPair.softclip_read(read, mock_cigar_util)

        self.assertEquals(42, read.__dict__["reference_start"])
        self.assertEquals("75X", read.__dict__["cigarstring"])
        self.assertEquals({"X0" : "X0:Z:chr42|150|-|target_1"}, read._tags)


    def test_softclip_read_readNotRecognized(self):
        clipper.PrimerPair._all_primers = {}
        clipper.PrimerPair(target_id="target_1",
                           chrom="chr42",
                           sense_start=100,
                           antisense_start=150,
                           sense_sequence="AACCGGTT",
                           antisense_sequence="CTTTA")
        read = MockRead(reference_name="chr42",
                        flag=0,
                        reference_start=42,
                        reference_end=242,
                        cigarstring="75M")
#        clipped_cigar_util = MockCigarUtil(reference_start=42, cigar="75X")
        mock_cigar_util = MockCigarUtil()

        clipper.PrimerPair.softclip_read(read, mock_cigar_util)

        self.assertEquals(42, read.__dict__["reference_start"])
        self.assertEquals("75M", read.__dict__["cigarstring"])
        self.assertEquals({"X0" : "X0:Z:PRIMER_NOT_FOUND"}, read._tags)
