#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods, unused-parameter
from __future__ import print_function, absolute_import

import sys
import unittest

from nose.exc import SkipTest
import pysam

from katana import util
from katana.util import ReadTransformation
from katana import readhandler

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class MicroMock(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def make_bam_file(filename, reads, header=None):
    #pylint: disable=no-member,too-many-arguments
    if header is None:
        header = { 'HD': {'VN': '1.0'},
                  'SQ': [{'LN': 1575, 'SN': 'chr1'},
                         {'LN': 1584, 'SN': 'chr2'}] }
    outfile = pysam.AlignmentFile(filename, "wb", header=header)
    for read in reads:
        outfile.write(read.aligned_segment)
    outfile.close()
    readhandler.PYSAM_ADAPTER.pysam_index(filename)

def build_aligned_segment(query_name = "read_28833_29006_6945",
               query_sequence="AGCTTAGCTA",
               flag = 99,
               reference_id = 0,
               reference_start = 32,
               mapping_quality = 20,
               cigar = None,
               next_reference_id = 0,
               next_reference_start=199,
               template_length=167,
               query_qualities = None):
    #pylint: disable=no-member,too-many-arguments
    a = pysam.AlignedSegment()
    a.query_name = query_name
    a.query_sequence = query_sequence
    a.flag = flag
    a.reference_id = reference_id
    a.reference_start = reference_start
    a.mapping_quality = mapping_quality
    if cigar is None:
        a.cigar = ((0, len(query_sequence)),)
    else:
        a.cigar = cigar
    a.next_reference_id = next_reference_id
    a.next_reference_start = next_reference_start
    a.template_length = template_length
    if query_qualities is None:
        a.query_qualities = [27] * len(query_sequence)
    return a

def build_read(query_name = "read_28833_29006_6945",
               query_sequence="AGCTTAGCTA",
               flag = 99,
               reference_id = 0,
               reference_start = 32,
               mapping_quality = 20,
               cigar = None,
               next_reference_id = 0,
               next_reference_start=199,
               template_length=167,
               query_qualities = None):
    #pylint: disable=no-member,too-many-arguments
    a = build_aligned_segment(query_name, query_sequence, flag, reference_id, reference_start, mapping_quality, cigar, next_reference_id, next_reference_start, template_length, query_qualities)
    return MicroMock(aligned_segment=a)


class MockLog(object):
    def __init__(self):
        self._log_calls = []

    def log(self, msg_format, *args):
        self._log_calls.append((msg_format, args))


class MockRead(MicroMock):
    def __init__(self, **kwargs):
        self._tags={}
        super(MockRead, self).__init__(**kwargs)

    def set_tag(self, tag_name, tag_value, tag_type):
        self._tags[tag_name] = "{}:{}:{}".format(tag_name, tag_type, tag_value)


class MockAlignedSegment(MicroMock):
    def __init__(self, **kwargs):
        self._tags={}
        super(MockAlignedSegment, self).__init__(**kwargs)

    def set_tag(self, tag_name, tag_value, tag_type):
        self._tags[tag_name] = "{}:{}:{}".format(tag_name, tag_type, tag_value)


class MockCigarUtil(MicroMock):
    def __init__(self, **kwargs):
        self._softclip_targets_calls=[]
        self._softclip_target_return = None
        self.reference_start=100
        self.cigar="10M"
        self.is_valid=True
        super(MockCigarUtil, self).__init__(**kwargs)

    def softclip_target(self, target_start, target_end):
        self._softclip_targets_calls.append((target_start, target_end))
        return self._softclip_target_return


class MockReadHandler(MicroMock):
    def __init__(self, **kwargs):
        self.begin_calls=0
        self.handle_calls = []
        self.end_calls = 0
        self.handle_raise = lambda x,y: False
        super(MockReadHandler, self).__init__(**kwargs)

    def begin(self):
        self.begin_calls += 1

    def handle(self, read, read_transformation, mate_transformation):
        self.handle_calls.append((read,
                                  read_transformation,
                                  mate_transformation))
        if self.handle_raise and self.handle_raise(read, read_transformation):
            raise StopIteration()

    def end(self):
        self.end_calls += 1


class MockPrimerPair(MicroMock):
    def __init__(self, **kwargs):
        self._softclip_primers_calls=[]
        self._softclip_primers_return = None
        super(MockPrimerPair, self).__init__(**kwargs)

    def add_tags(self, read):
        pass

    def softclip_primers(self, old_cigar):
        self._softclip_primers_calls.append(old_cigar)
        return self._softclip_primers_return


class KatanaBaseTestCase(unittest.TestCase):
    def setUp(self):
        util.PrimerPair._all_primers = {}
        unittest.TestCase.setUp(self)
        self.stderr = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.stderr

    def tearDown(self):
        self.stderr.close()
        sys.stderr = self.saved_stderr
        unittest.TestCase.tearDown(self)

    @staticmethod
    def check_sysout_safe():
        try:
            # nosetest and pysam fight over stdout; run nosetest -s to be safe
            sys.stdout.fileno() #pylint disable=pointless-statement
        except Exception: #pylint disable=broad=except
            raise SkipTest("sysout unsafe to test: run nosetest with -s option")

class PrimerStatsTestCase(KatanaBaseTestCase):
    def test_stat_keys(self):
        self.assertEquals(6, len(util.PrimerStats.STAT_KEYS))

    def test_total_read_count(self):
        read_sense = MockRead(is_positive_strand=True)
        primer_pair1 = MockPrimerPair(target_id="target_2",
                                      chrom="chr2",
                                      sense_start=222)
        stats = util.PrimerStats()
        self.assertEquals(0, stats.total_read_count)
        stats.add_read_primer(read_sense, primer_pair1)
        stats.add_read_primer(read_sense, primer_pair1)
        stats.add_read_primer(read_sense, primer_pair1)
        self.assertEquals(3, stats.total_read_count)

    def test_primer_pairs(self):
        read = MockRead(is_positive_strand=True)
        primer_pair1 = MockPrimerPair(target_id="target_1",
                                     chrom="chr1",
                                     sense_start=222)
        stats = util.PrimerStats()
        stats.add_read_primer(read, primer_pair1)
        self.assertEquals([primer_pair1], stats.primer_pairs)
        stats.add_read_primer(read, primer_pair1)
        self.assertEquals([primer_pair1], stats.primer_pairs)
        primer_pair2 = MockPrimerPair(target_id="target_1",
                                     chrom="chr1",
                                     sense_start=555)
        stats.add_read_primer(read, primer_pair2)
        self.assertEquals([primer_pair1, primer_pair2], stats.primer_pairs)

    def test_primer_pairs_sortedChromStartTarget(self):
        read = MockRead(is_positive_strand=True)
        expected_primer_pairs=[]
        expected_primer_pairs.append(MockPrimerPair(target_id="target_2",
                                                   chrom="chr2",
                                                   sense_start=222))
        expected_primer_pairs.append(MockPrimerPair(target_id="target_10",
                                                   chrom="chr2",
                                                   sense_start=222))
        expected_primer_pairs.append(MockPrimerPair(target_id="target_1",
                                                   chrom="chr2",
                                                   sense_start=333))
        expected_primer_pairs.append(MockPrimerPair(target_id="target_1",
                                                   chrom="chr10",
                                                   sense_start=111))
        expected_primer_pairs.append(MockPrimerPair(target_id="target_1",
                                                   chrom="chr10",
                                                   sense_start=222))
        stats = util.PrimerStats()
        for primer_pair in expected_primer_pairs[::-1]:
            stats.add_read_primer(read, primer_pair)
        self.assertEquals(expected_primer_pairs,
                          stats.primer_pairs)

    def test_stats(self):
        read_sense = MockRead(is_positive_strand=True)
        read_antisense = MockRead(is_positive_strand=False)
        primer_pair1 = MockPrimerPair(target_id="target_2",
                                      chrom="chr2",
                                      sense_start=222)
        primer_pair2 = MockPrimerPair(target_id="target_10",
                                      chrom="chr4",
                                      sense_start=444)
        stats = util.PrimerStats()
        stats.add_read_primer(read_sense, primer_pair1)
        stats.add_read_primer(read_sense, primer_pair2)
        stats.add_read_primer(read_sense, primer_pair2)
        stats.add_read_primer(read_antisense, primer_pair2)

        self.assertEquals(2, len(stats.primer_pairs))

        stats1 = stats.stats(primer_pair1)
        self.assertEquals(len(util.PrimerStats.STAT_KEYS), len(stats1))
        self.assertEquals("target_2", stats1["target_id"])
        self.assertEquals("chr2", stats1["chrom"])
        self.assertEquals(222, stats1["sense_start"])
        self.assertEquals(1, stats1["sense_count"])
        self.assertEquals(0, stats1["antisense_count"])
        self.assertEquals(100, stats1["sense_percent"])

        stats2 = stats.stats(primer_pair2)
        self.assertEquals(len(util.PrimerStats.STAT_KEYS), len(stats2))
        self.assertEquals("target_10", stats2["target_id"])
        self.assertEquals("chr4", stats2["chrom"])
        self.assertEquals(444, stats2["sense_start"])
        self.assertEquals(2, stats2["sense_count"])
        self.assertEquals(1, stats2["antisense_count"])
        self.assertEquals(66, stats2["sense_percent"])


class PrimerStatsDumperTestCase(KatanaBaseTestCase):
    def test_dump(self):
        mock_log = MockLog()
        dumper = util.PrimerStatsDumper(log_method=mock_log.log)
        stat_dict = {"primer1":{"statA":"A1", "statB":"B1"},
                     "primer2":{"statA":"A2", "statB":"B2"}}
        primer_stats = MicroMock(STAT_KEYS=["statA", "statB"],
                                       primer_pairs=["primer1", "primer2"],
                                       stats=lambda x: stat_dict[x])
        dumper.dump(primer_stats)
        self.assertEquals(3, len(mock_log._log_calls))
        self.assertEquals('PRIMER_STATS|statA|statB', mock_log._log_calls[0][0])
        self.assertEquals('PRIMER_STATS|A1|B1', mock_log._log_calls[1][0])
        self.assertEquals('PRIMER_STATS|A2|B2', mock_log._log_calls[2][0])


class MockAlignmentFile(object):
    def __init__(self, reference_id_to_name):
        self.reference_id_to_name = reference_id_to_name

    def getrname(self, reference_id):
        return self.reference_id_to_name[reference_id]

class ReadTestCase(KatanaBaseTestCase):
    def test_init(self):
        mock_alignment_file = MockAlignmentFile({1:"chr1"})
        mock_aligned_segment = MockAlignedSegment(query_name="readA",
                                                  reference_id=1,
                                                  reference_start=100,
                                                  reference_end=110,
                                                  cigarstring="10M",
                                                  is_paired=True)
        read = util.Read(mock_aligned_segment, mock_alignment_file)
        self.assertEquals("readA", read.query_name)
        self.assertEquals("chr1", read.reference_name)
        self.assertEquals(100, read.reference_start)
        self.assertEquals(110, read.reference_end)
        self.assertEquals("10M", read.cigarstring)
        self.assertEquals(True, read.is_paired)

    def test_mutatorsPassThroughToAlignedSegment(self):
        mock_alignment_file = MockAlignmentFile({1:"chr1"})
        aligned_segment = build_aligned_segment(query_name="read1",
                                                  reference_id=1,
                                                  reference_start=100,
                                                  cigar=((0, 10),),
                                                  next_reference_start=150)
        read = util.Read(aligned_segment, mock_alignment_file)
        read.reference_start = 142
        read.cigarstring = "10S"
        read.next_reference_start = 200
        self.assertEquals(142,
                          aligned_segment.reference_start)
        self.assertEquals("10S",
                          aligned_segment.cigarstring)
        self.assertEquals(200,
                          aligned_segment.next_reference_start)

    def test_isPaired_mutatorPassedThroughToAlignedSegment(self):
        mock_alignment_file = MockAlignmentFile({1:"chr1"})
        aligned_segment = build_aligned_segment(query_name="read1",
                                                reference_id=1,
                                                reference_start=100,
                                                cigar=((0, 10),),
                                                next_reference_id = 1,
                                                next_reference_start=200)
        aligned_segment.is_paired=True
        aligned_segment.mate_is_unmapped = True
        aligned_segment.mate_is_reverse = True
        aligned_segment.is_proper_pair = True
        aligned_segment.is_read1 = True
        original_flag = aligned_segment.flag

        read = util.Read(aligned_segment, mock_alignment_file)

        read.is_paired = False

        self.assertEquals(107, original_flag)
        self.assertEquals(0, aligned_segment.flag)
        self.assertEquals(False,
                          aligned_segment.is_paired)
        self.assertEquals(False,
                          aligned_segment.mate_is_unmapped)
        self.assertEquals(False,
                          aligned_segment.is_proper_pair)
        self.assertEquals(False,
                          aligned_segment.is_read1)
        self.assertEquals(-1,
                          aligned_segment.next_reference_id)
        self.assertEquals(-1,
                          aligned_segment.next_reference_start)

    def test_is_positive(self):
        mock_alignment_file = MockAlignmentFile({1:"chr1"})
        read = util.Read(MockAlignedSegment(is_reverse=False),
                         mock_alignment_file)
        self.assertEquals(True, read.is_positive_strand)
        read = util.Read(MockAlignedSegment(is_reverse=True),
                         mock_alignment_file)
        self.assertEquals(False, read.is_positive_strand)

    def test_key(self):
        mock_alignment_file = MockAlignmentFile({1:"chr1"})
        mock_aligned_segment = MockAlignedSegment(query_name="readA",
                                                  reference_id=1,
                                                  reference_start=100,
                                                  is_reverse=False,
                                                  is_read1=True)
        read = util.Read(mock_aligned_segment, mock_alignment_file)
        expected_key = ("readA", True, "chr1", 100, True)
        self.assertEquals(expected_key, read.key)

    def test_mate_key(self):
        mock_alignment_file = MockAlignmentFile({2:"chr2"})
        mock_aligned_segment = MockAlignedSegment(query_name="readA",
                                                  is_paired=True,
                                                  mate_is_unmapped=False,
                                                  mate_is_reverse=True,
                                                  next_reference_id=2,
                                                  next_reference_start=200,
                                                  is_read1=True)
        read = util.Read(mock_aligned_segment, mock_alignment_file)
        expected_key = ("readA", False, "chr2", 200, False)
        self.assertEquals(expected_key, read.mate_key)

    def test_mate_key_noneWhenNoMate(self):
        mock_alignment_file = MockAlignmentFile({2:"chr2"})
        #pylint: disable=attribute-defined-outside-init
        mock_aligned_segment = MockAlignedSegment(query_name="read1",
                                                  is_paired=True,
                                                  mate_is_unmapped=False,
                                                  mate_is_reverse=True,
                                                  next_reference_id=2,
                                                  next_reference_start=200)
        mock_aligned_segment.is_paired=False
        read = util.Read(mock_aligned_segment, mock_alignment_file)
        self.assertEquals(None, read.mate_key)
        mock_aligned_segment.is_paired=True

        mock_aligned_segment.mate_is_unmapped=True
        read = util.Read(mock_aligned_segment, mock_alignment_file)
        self.assertEquals(None, read.mate_key)

    def test_set_tag(self):
        mock_alignment_file = MockAlignmentFile({2:"chr2"})
        mock_aligned_segment = MockAlignedSegment()
        read = util.Read(mock_aligned_segment, mock_alignment_file)
        read.set_tag("name", "value", "type")
        self.assertEquals("name:type:value", mock_aligned_segment._tags["name"])

    def test_iter(self):
        mock_alignment_file = MockAlignmentFile({1:"chr1"})
        aligned_segment1 = MockAlignedSegment(query_name="read1")
        aligned_segment2 = MockAlignedSegment(query_name="read2")
        aligned_segment_iter = iter([aligned_segment1, aligned_segment2])
        actual_iter = util.Read.iter(aligned_segment_iter, mock_alignment_file)
        actual_reads = [read for read in actual_iter]
        self.assertEquals(2, len(actual_reads))
        self.assertIsInstance(actual_reads[0], util.Read)
        self.assertEquals(aligned_segment1, actual_reads[0].aligned_segment)
        self.assertIsInstance(actual_reads[1], util.Read)
        self.assertEquals(aligned_segment2, actual_reads[1].aligned_segment)

    def test_mate_cigar(self):
        mock_alignment_file = MockAlignmentFile({2:"chr2"})
        aligned_segment = build_aligned_segment()
        aligned_segment.set_tag('MC', '10M')
        read = util.Read(aligned_segment, mock_alignment_file)

        self.assertEqual('10M', read.mate_cigar)
        read.mate_cigar = '2S8M'
        self.assertEqual('2S8M', aligned_segment.get_tag('MC'))

class ReadTransformationTestCase(KatanaBaseTestCase):
    '''Lightweight container of what we need to update the read.'''
    def test_init(self):
        read = MockRead(is_unmapped=False, is_paired=True)
        new_cigar_util = MockCigarUtil(reference_start=42,
                                       cigar="10M",
                                       is_valid=True)
        mock_primer_pair = MockPrimerPair(is_unmatched=False)
        filter_builder = lambda x: ["filter1"]

        transform = util.ReadTransformation(read,
                                            mock_primer_pair,
                                            new_cigar_util,
                                            filter_builder)
        self.assertEquals(mock_primer_pair, transform.primer_pair)
        self.assertEquals(42, transform.reference_start)
        self.assertEquals("10M", transform.cigar)
        self.assertEquals(False, transform.is_unmapped)
        self.assertEquals(True, transform.is_cigar_valid)
        self.assertEquals(True, transform.is_paired)
        self.assertEquals(("filter1",), transform.filters)

    def test_NULL(self):
        transform = ReadTransformation.NULL
        self.assertEquals(util.PrimerPair.NULL, transform.primer_pair)
        self.assertEquals(0, transform.reference_start)
        self.assertEquals("", transform.cigar)
        self.assertEquals(False, transform.is_unmapped)
        self.assertEquals(False, transform.is_cigar_valid)
        self.assertEquals(False, transform.is_paired)
        self.assertEquals(("NULL",), transform.filters)


    def test_eq(self):
        read = MockRead(is_unmapped=False, is_paired=True)
        primer_pair = MockPrimerPair(is_unmatched=False)
        new_cigar_util = MockCigarUtil(cigar="10M")
        filter_builder = lambda x: ["filter1"]
        base = util.ReadTransformation(read,
                                       primer_pair,
                                       new_cigar_util,
                                       filter_builder)
        same = util.ReadTransformation(read,
                                       primer_pair,
                                       new_cigar_util,
                                       filter_builder)
        self.assertEquals(base, same)
        read2 = MockRead(is_unmapped=True, is_paired=False)
        diff_read = util.ReadTransformation(read2,
                                            primer_pair,
                                            new_cigar_util,
                                            filter_builder)
        self.assertNotEquals(base, diff_read)
        primer_pair2 = MockPrimerPair(is_unmatched=False)
        diff_primer = util.ReadTransformation(read,
                                              primer_pair2,
                                              new_cigar_util,
                                              filter_builder)
        self.assertNotEquals(base, diff_primer)
        new_cigar_util2 = MockCigarUtil(cigar="20M")
        diff_cigar = util.ReadTransformation(read,
                                             primer_pair,
                                             new_cigar_util2,
                                             filter_builder)
        self.assertNotEquals(base, diff_cigar)
        filter_builder2 = lambda x: ["filter2"]
        diff_filter = util.ReadTransformation(read,
                                              primer_pair,
                                              new_cigar_util,
                                              filter_builder2)
        self.assertNotEquals(base, diff_filter)


class NullPrimerPairTestCase(KatanaBaseTestCase):
    def test_init(self):
        null_primer_pair = util._NullPrimerPair()
        old_cigar = MockCigarUtil()
        new_cigar = null_primer_pair.softclip_primers(old_cigar)
        self.assertEquals(old_cigar, new_cigar)
        self.assertEquals(True, null_primer_pair.is_unmatched)

class PrimerPairTestCase(KatanaBaseTestCase):
    def test_init(self):
        primer_pair = util.PrimerPair(target_id="target_1",
                                         chrom="chr42",
                                         sense_primer_region=(100,110),
                                         antisense_primer_region=(140,150))
        self.assertEquals("target_1", primer_pair.target_id)
        self.assertEquals("chr42", primer_pair.chrom)
        self.assertEquals(100, primer_pair.sense_start)
        self.assertEquals(110, primer_pair._query_region_start)
        self.assertEquals(140, primer_pair._query_region_end)
        self.assertEquals(False, primer_pair.is_unmatched)

    def test_all_primers(self):
        primer_pair1 = util.PrimerPair(target_id="target_1",
                                          chrom="chr42",
                                          sense_primer_region=(100,108),
                                          antisense_primer_region=(145,150))
        primer_pair2 = util.PrimerPair(target_id="target_1",
                                          chrom="chr42",
                                          sense_primer_region=(200,208),
                                          antisense_primer_region=(245,250))
        actual_primers = util.PrimerPair._all_primers
        self.assertEquals(primer_pair1, actual_primers[('chr42', 100, True)])
        self.assertEquals(primer_pair1, actual_primers[('chr42', 150, False)])
        self.assertEquals(primer_pair2, actual_primers[('chr42', 200, True)])
        self.assertEquals(primer_pair2, actual_primers[('chr42', 250, False)])
        self.assertEquals(4, len(actual_primers))

    def test_get_primer_pair_matchPositiveStrand(self):
        util.PrimerPair._all_primers = {}
        util.PrimerPair(target_id="target_1",
                           chrom="chr42",
                           sense_primer_region=(100,110),
                           antisense_primer_region=(140,150))
        read = MockRead(reference_name="chr42",
                        is_positive_strand=True,
                        reference_start=100,
                        reference_end=242,
                        cigarstring="75M")
        actual_primer_pair = util.PrimerPair.get_primer_pair(read)
        self.assertEquals("target_1", actual_primer_pair.target_id)

    def test_get_primer_pair_matchNegativeStrand(self):
        util.PrimerPair._all_primers = {}
        util.PrimerPair(target_id="target_2",
                           chrom="chr42",
                           sense_primer_region=(100,110),
                           antisense_primer_region=(140,150))
        read = MockRead(reference_name="chr42",
                        is_positive_strand=False,
                        reference_start=42,
                        reference_end=150,
                        cigarstring="75M")
        actual_primer_pair = util.PrimerPair.get_primer_pair(read)
        self.assertEquals("target_2", actual_primer_pair.target_id)

    def test_get_primer_pair_noMatch(self):
        util.PrimerPair._all_primers = {}
        util.PrimerPair(target_id="target_1",
                           chrom="chr42",
                           sense_primer_region=(100,110),
                           antisense_primer_region=(140,150))
        read = MockRead(reference_name="chrX",
                        is_positive_strand=True,
                        reference_start=42,
                        reference_end=142,
                        cigarstring="75M")
        actual_primer_pair = util.PrimerPair.get_primer_pair(read)
        self.assertIsInstance(actual_primer_pair, util._NullPrimerPair)

    def test_softclip_primers(self):
        primer_pair = util.PrimerPair(target_id="target_1",
                                         chrom="chr42",
                                         sense_primer_region=(100,110),
                                         antisense_primer_region=(140,150))

        expect_clipped_cigar = MockCigarUtil(reference_start=42, cigar="75X")
        mock_cigar = MockCigarUtil(_softclip_target_return=expect_clipped_cigar)
        actual_clipped_cigar = primer_pair.softclip_primers(mock_cigar)
        self.assertEquals([(110, 140)], mock_cigar._softclip_targets_calls)
        self.assertEquals(expect_clipped_cigar, actual_clipped_cigar)
