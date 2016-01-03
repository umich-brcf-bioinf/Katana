#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods

from __future__ import print_function, absolute_import
from ampliconsoftclipper import clipper
import os.path
import pysam
from testfixtures.tempdirectory import TempDirectory
import unittest
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class ClipperBaseTestCase(unittest.TestCase):
    def setUp(self):
        clipper._PrimerPair._all_primers = {}
        unittest.TestCase.setUp(self)
        self.stderr = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.stderr

    def tearDown(self):
        self.stderr.close()
        sys.stderr = self.saved_stderr
        unittest.TestCase.tearDown(self)

class MicroMock(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

class MockPrimerStatsDumper(MicroMock):
    def __init__(self, **kwargs):
        self._dump_calls=[]
        super(MockPrimerStatsDumper, self).__init__(**kwargs)

    def dump(self, primer_stats):
        self._dump_calls.append(primer_stats)

class MockLog(object):
    def __init__(self):
        self._log_calls = []

    def log(self, msg_format, *args):
        self._log_calls.append((msg_format, args))

class MockPrimerStats(MicroMock):
    def __init__(self, **kwargs):
        self._add_read_primer_calls=[]
        super(MockPrimerStats, self).__init__(**kwargs)

    def add_read_primer(self, read, primer):
        self._add_read_primer_calls.append((read, primer))

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
        super(MockCigarUtil, self).__init__(**kwargs)

    def softclip_target(self, target_start, target_end):
        self._softclip_targets_calls.append((target_start, target_end))
        return self._softclip_target_return


class MockReadHandler(MicroMock):
    def __init__(self, **kwargs):
        self.begin_calls=0
        self.handle_calls = []
        self.end_calls = 0
        super(MockReadHandler, self).__init__(**kwargs)

    def begin(self):
        self.begin_calls += 1

    def handle(self, read):
        self.handle_calls.append(read)

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


class MockWriter(object):
    def __init__(self):
        self._write_calls = []

    def write(self, obj):
        self._write_calls.append(obj)

class ClipperTestCase(ClipperBaseTestCase):
    def test_build_read_transformations(self):
        primer_pair = clipper._PrimerPair(target_id="target_1",
                                         chrom="chr42",
                                         sense_primer_region=(100,102),
                                         antisense_primer_region=(148,150))
        read1 = MockRead(query_name="read1-sense",
                        is_positive_strand=True,
                        reference_name="chr42",
                        reference_start=100,
                        reference_end=150,
                        cigarstring="10M",
                        key=111)
        read2 = MockRead(query_name="read2-antisense",
                        is_positive_strand=False,
                        reference_name="chr42",
                        reference_start=140,
                        reference_end=150,
                        cigarstring="10M",
                        key=222)
        read3 = MockRead(query_name="read3",
                        is_positive_strand=True,
                        reference_name="chr42",
                        reference_start=333,
                        reference_end=343,
                        cigarstring="10M",
                        key=333)

        read_iter = iter([read1, read2, read3])
        actual_read_transforms = clipper._build_read_transformations(read_iter)
        self.assertEquals(3, len(actual_read_transforms))
        self.assertEquals((primer_pair, 102, "2S8M"),
                          actual_read_transforms[111])
        self.assertEquals((primer_pair, 140, "8M2S"),
                          actual_read_transforms[222])
        self.assertEquals((clipper._PrimerPair.NULL_PRIMER_PAIR, 333, "10M"),
                          actual_read_transforms[333])

    def test_handle_reads(self):
        handler1 = MockReadHandler()
        handler2 = MockReadHandler()
        handler3 = MockReadHandler()
        handlers = [handler1, handler2, handler3]
        read1 = MockRead(query_name="read1-sense")
        read2 = MockRead(query_name="read2-antisense")
        read3 = MockRead(query_name="read3")
        read_iter = iter([read1, read2, read3])
        clipper._handle_reads(read_iter, handlers)
        self.assertEquals(1, handler1.begin_calls)
        self.assertEquals(1, handler2.begin_calls)
        self.assertEquals(1, handler3.begin_calls)
        self.assertEquals([read1, read2, read3], handler1.handle_calls)
        self.assertEquals([read1, read2, read3], handler2.handle_calls)
        self.assertEquals([read1, read2, read3], handler3.handle_calls)
        self.assertEquals(1, handler1.end_calls)
        self.assertEquals(1, handler2.end_calls)
        self.assertEquals(1, handler3.end_calls)

    def test_initialize_primer_pairs(self):
        clipper._PrimerPair._all_primers = {}
        #pylint: disable=line-too-long
        input_manifest = (\
'''id|Customer TargetID|Chr|Genome Build|Sense Start|Antisense Start|Sense Sequence|Antisense Sequence
1|NRAS_3|1|hg19 dbSNP137|111|222|ATCCGACAA|AGTACAAACT
2|IDH1_1|2|hg19 dbSNP137|333|444|TGCCAACAT|TGGAAATCAC
''').replace("|", "\t")
        clipper._initialize_primer_pairs(StringIO(input_manifest))
        actual_primers = clipper._PrimerPair._all_primers
        self.assertEquals(4, len(actual_primers))
        targets = set([primer.target_id for primer in actual_primers.values()])
        self.assertEquals(set(["IDH1_1", "NRAS_3"]), targets)

class AddTagsReadHandlerTestCase(ClipperBaseTestCase):
    def test_handle(self):
        #pylint: disable=no-member,too-many-arguments
        original_reference_start = 100
        original_reference_end = 110
        original_cigar_string = "10M"
        read = MockRead(key=42,
                        reference_start=original_reference_start,
                        reference_end=original_reference_end,
                        cigarstring=original_cigar_string)
        primer_pair_target = "target_1"
        primer_pair = MockPrimerPair(target_id=primer_pair_target)

        new_reference_start = 102
        new_cigar_string = "2S8M"
        transformations = {42: (primer_pair,
                                new_reference_start,
                                new_cigar_string)}
        handler = clipper._AddTagsReadHandler(transformations)

        handler.handle(read)

        self.assertEquals("X0:Z:" + primer_pair_target,
                          read._tags["X0"])
        self.assertEquals("X1:Z:" + original_cigar_string,
                          read._tags["X1"])
        self.assertEquals("X2:i:" + str(original_reference_start),
                          read._tags["X2"])
        self.assertEquals("X3:i:" + str(original_reference_end),
                          read._tags["X3"])

class PrimerStatsTestCase(ClipperBaseTestCase):
    def test_stat_keys(self):
        self.assertEquals(7, len(clipper._PrimerStats.STAT_KEYS))

    def test_total_read_count(self):
        read_sense = MockRead(is_positive_strand=True)
        primer_pair1 = MockPrimerPair(target_id="target_2",
                                      chrom="chr2",
                                      sense_start=222)
        stats = clipper._PrimerStats()
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
        stats = clipper._PrimerStats()
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
        stats = clipper._PrimerStats()
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
        stats = clipper._PrimerStats()
        stats.add_read_primer(read_sense, primer_pair1)
        stats.add_read_primer(read_sense, primer_pair2)
        stats.add_read_primer(read_sense, primer_pair2)
        stats.add_read_primer(read_antisense, primer_pair2)

        self.assertEquals(2, len(stats.primer_pairs))

        stats1 = stats.stats(primer_pair1)
        self.assertEquals(len(clipper._PrimerStats.STAT_KEYS), len(stats1))
        self.assertEquals("target_2", stats1["target_id"])
        self.assertEquals("chr2", stats1["chrom"])
        self.assertEquals(222, stats1["sense_start"])
        self.assertEquals(1, stats1["sense_count"])
        self.assertEquals(0, stats1["antisense_count"])
        self.assertEquals(25, stats1["sense_percent"])
        self.assertEquals(0, stats1["antisense_percent"])

        stats2 = stats.stats(primer_pair2)
        self.assertEquals(len(clipper._PrimerStats.STAT_KEYS), len(stats2))
        self.assertEquals("target_10", stats2["target_id"])
        self.assertEquals("chr4", stats2["chrom"])
        self.assertEquals(444, stats2["sense_start"])
        self.assertEquals(2, stats2["sense_count"])
        self.assertEquals(1, stats2["antisense_count"])
        self.assertEquals(50, stats2["sense_percent"])
        self.assertEquals(25, stats2["antisense_percent"])


class PrimerStatsDumperTestCase(ClipperBaseTestCase):
    def test_dump(self):
        mock_log = MockLog()
        dumper = clipper._PrimerStatsDumper(log_method=mock_log.log)
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


class StatsReadHandlerTestCase(ClipperBaseTestCase):
    def test_handle(self):
        read1 = MockRead(key=42, is_positive_strand=True)
        read2 = MockRead(key=42, is_positive_strand=True)
        primer_pair = MockPrimerPair(target_id="target_1",
                                     chrom="chr1",
                                     sense_start=242)
        new_reference_start = 102
        new_cigar_string = "2S8M"
        transformations = {42: (primer_pair,
                                new_reference_start,
                                new_cigar_string)}
        mock_primer_stats = MockPrimerStats()
        mock_primer_stats_dumper = MockPrimerStatsDumper()
        handler = clipper._StatsHandler(transformations,
                                               mock_primer_stats,
                                               mock_primer_stats_dumper)

        handler.handle(read1)
        handler.handle(read2)
        expected_calls = [(read1, primer_pair), (read2, primer_pair)]
        self.assertEquals(expected_calls,
                          mock_primer_stats._add_read_primer_calls)

        handler.end()
        self.assertEquals([mock_primer_stats],
                          mock_primer_stats_dumper._dump_calls)


class TransformReadHandlerTestCase(ClipperBaseTestCase):
    def test_handle(self):
        #pylint: disable=no-member,too-many-arguments
        read = MockRead(key=42, reference_start=100, cigarstring="10M")
        primer_pair = None
        new_reference_start = 102
        new_cigar_string = "2S8M"
        transformations = {42: (primer_pair,
                                new_reference_start,
                                new_cigar_string)}
        handler = clipper._TransformReadHandler(transformations)
        handler.handle(read)
        self.assertEquals(102, read.reference_start)
        self.assertEquals("2S8M", read.cigarstring)

class WriteReadHandlerTestCase(ClipperBaseTestCase):
    #pylint: disable=no-member,too-many-arguments
    @staticmethod
    def make_bam_file(filename, reads, header=None):
        if header is None:
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
        outfile = pysam.AlignmentFile(filename, "wb", header=header)
        for read in reads:
            outfile.write(read.aligned_segment)
        outfile.close()
        clipper.PYSAM_INDEX(filename)

    @staticmethod
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
                   query_qualities = None,
                   tags = None):
        a = pysam.AlignedSegment()
        a.query_name = query_name
        a.query_sequence = query_sequence
        a.flag = flag
        a.reference_id = reference_id
        a.reference_start = reference_start
        a.mapping_quality = mapping_quality
        if cigar is None:
            a.cigar = ((0,10), (2,1), (0,25))
        a.next_reference_id = next_reference_id
        a.next_reference_start = next_reference_start
        a.template_length = template_length
        if query_qualities is None:
            a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<")
        if tags is None:
            a.tags = (("NM", 1),
                      ("RG", "L1"))
        return clipper._Read(a)

    def test_handle_sortsAndIndexes(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_bam_filename = os.path.join(input_dir.path, "input.bam")
            self.make_bam_file(input_bam_filename, [self.build_read()])
            output_bam_filename = os.path.join(output_dir.path, "output.bam")
            handler = clipper._WriteReadHandler(input_bam_filename,
                                               output_bam_filename)
            read1 = self.build_read(query_name="read1",
                                    reference_id=0,
                                    reference_start=20)
            read2 = self.build_read(query_name="read2",
                                    reference_id=0,
                                    reference_start=10)

            handler.begin()
            handler.handle(read1)
            handler.handle(read2)
            handler.end()

            actual_files = sorted(os.listdir(output_dir.path))
            self.assertEquals(["output.bam", "output.bam.bai"], actual_files)
            actual_bam = pysam.AlignmentFile(output_bam_filename, "rb")
            actual_reads = [read for read in actual_bam.fetch()]
            actual_bam.close()

        self.assertEquals(2, len(actual_reads))
        self.assertEquals("read2", actual_reads[0].query_name)
        self.assertEquals("read1", actual_reads[1].query_name)

class ReadTestCase(ClipperBaseTestCase):
    def test_init(self):
        mock_aligned_segment = MockAlignedSegment(query_name="read1",
                                                  reference_name="chr1",
                                                  reference_start=100,
                                                  reference_end=110,
                                                  cigarstring="10M")
        read = clipper._Read(mock_aligned_segment)
        self.assertEquals("chr1", read.reference_name)
        self.assertEquals(100, read.reference_start)
        self.assertEquals(110, read.reference_end)
        self.assertEquals("10M", read.cigarstring)

    def test_mutatorsPassThroughToAlignedSegment(self):
        mock_aligned_segment = MockAlignedSegment(query_name="read1",
                                                  reference_name="chr1",
                                                  reference_start=100,
                                                  cigarstring="10M")
        read = clipper._Read(mock_aligned_segment)
        read.reference_start = 142
        read.cigarstring = "10S"
        self.assertEquals(142, mock_aligned_segment.__dict__['reference_start'])
        self.assertEquals("10S", mock_aligned_segment.__dict__['cigarstring'])

    def test_is_positive(self):
        read = clipper._Read(MockAlignedSegment(is_reverse=False))
        self.assertEquals(True, read.is_positive_strand)
        read = clipper._Read(MockAlignedSegment(is_reverse=True))
        self.assertEquals(False, read.is_positive_strand)

    def test_key(self):
        mock_aligned_segment = MockAlignedSegment(query_name="read1",
                                                  reference_name="chr1",
                                                  reference_start=100,
                                                  is_reverse=False)
        read = clipper._Read(mock_aligned_segment)
        expected_key = ("read1", True, "chr1", 100)
        self.assertEquals(expected_key, read.key)

    def test_mate_key(self):
        mock_aligned_segment = MockAlignedSegment(query_name="read1",
                                                  is_paired=True,
                                                  mate_is_unmapped=False,
                                                  mate_is_reverse=True,
                                                  next_reference_name="chr2",
                                                  next_reference_start=200)
        read = clipper._Read(mock_aligned_segment)
        expected_key = ("read1", False, "chr2", 200)
        self.assertEquals(expected_key, read.mate_key)

    def test_mate_key_noneWhenNoMate(self):
        #pylint: disable=attribute-defined-outside-init
        mock_aligned_segment = MockAlignedSegment(query_name="read1",
                                                  is_paired=True,
                                                  mate_is_unmapped=False,
                                                  mate_is_reverse=True,
                                                  next_reference_name="chr2",
                                                  next_reference_start=200)
        mock_aligned_segment.is_paired=False
        self.assertEquals(None, clipper._Read(mock_aligned_segment).mate_key)
        mock_aligned_segment.is_paired=True

        mock_aligned_segment.mate_is_unmapped=True
        self.assertEquals(None, clipper._Read(mock_aligned_segment).mate_key)
        mock_aligned_segment.mate_is_unmapped=True

    def test_set_tag(self):
        mock_aligned_segment = MockAlignedSegment()
        read = clipper._Read(mock_aligned_segment)
        read.set_tag("name", "value", "type")
        self.assertEquals("name:type:value", mock_aligned_segment._tags["name"])

    def test_iter(self):
        aligned_segment1 = MockAlignedSegment(query_name="read1")
        aligned_segment2 = MockAlignedSegment(query_name="read2")
        aligned_segment_iter = iter([aligned_segment1, aligned_segment2])
        actual_iter = clipper._Read.iter(aligned_segment_iter)
        actual_reads = [read for read in actual_iter]
        self.assertEquals(2, len(actual_reads))
        self.assertIsInstance(actual_reads[0], clipper._Read)
        self.assertEquals(aligned_segment1, actual_reads[0].aligned_segment)
        self.assertIsInstance(actual_reads[1], clipper._Read)
        self.assertEquals(aligned_segment2, actual_reads[1].aligned_segment)

class PrimerPairTestCase(ClipperBaseTestCase):
    def test_init(self):
        primer_pair = clipper._PrimerPair(target_id="target_1",
                                         chrom="chr42",
                                         sense_primer_region=(100,110),
                                         antisense_primer_region=(140,150))
        self.assertEquals("target_1", primer_pair.target_id)
        self.assertEquals("chr42", primer_pair.chrom)
        self.assertEquals(100, primer_pair.sense_start)
        self.assertEquals(110, primer_pair._query_region_start)
        self.assertEquals(140, primer_pair._query_region_end)

    def test_all_primers(self):
        primer_pair1 = clipper._PrimerPair(target_id="target_1",
                                          chrom="chr42",
                                          sense_primer_region=(100,108),
                                          antisense_primer_region=(145,150))
        primer_pair2 = clipper._PrimerPair(target_id="target_1",
                                          chrom="chr42",
                                          sense_primer_region=(200,208),
                                          antisense_primer_region=(245,250))
        actual_primers = clipper._PrimerPair._all_primers
        self.assertEquals(primer_pair1, actual_primers[('chr42', 100, True)])
        self.assertEquals(primer_pair1, actual_primers[('chr42', 150, False)])
        self.assertEquals(primer_pair2, actual_primers[('chr42', 200, True)])
        self.assertEquals(primer_pair2, actual_primers[('chr42', 250, False)])
        self.assertEquals(4, len(actual_primers))

    def test_get_primer_pair_matchPositiveStrand(self):
        clipper._PrimerPair._all_primers = {}
        clipper._PrimerPair(target_id="target_1",
                           chrom="chr42",
                           sense_primer_region=(100,110),
                           antisense_primer_region=(140,150))
        read = MockRead(reference_name="chr42",
                        is_positive_strand=True,
                        reference_start=100,
                        reference_end=242,
                        cigarstring="75M")
        actual_primer_pair = clipper._PrimerPair.get_primer_pair(read)
        self.assertEquals("target_1", actual_primer_pair.target_id)

    def test_get_primer_pair_matchNegativeStrand(self):
        clipper._PrimerPair._all_primers = {}
        clipper._PrimerPair(target_id="target_2",
                           chrom="chr42",
                           sense_primer_region=(100,110),
                           antisense_primer_region=(140,150))
        read = MockRead(reference_name="chr42",
                        is_positive_strand=False,
                        reference_start=42,
                        reference_end=150,
                        cigarstring="75M")
        actual_primer_pair = clipper._PrimerPair.get_primer_pair(read)
        self.assertEquals("target_2", actual_primer_pair.target_id)

    def test_get_primer_pair_noMatch(self):
        clipper._PrimerPair._all_primers = {}
        clipper._PrimerPair(target_id="target_1",
                           chrom="chr42",
                           sense_primer_region=(100,110),
                           antisense_primer_region=(140,150))
        read = MockRead(reference_name="chrX",
                        is_positive_strand=True,
                        reference_start=42,
                        reference_end=142,
                        cigarstring="75M")
        actual_primer_pair = clipper._PrimerPair.get_primer_pair(read)
        self.assertIsInstance(actual_primer_pair, clipper._NullPrimerPair)

    def test_softclip_primers(self):
        primer_pair = clipper._PrimerPair(target_id="target_1",
                                         chrom="chr42",
                                         sense_primer_region=(100,110),
                                         antisense_primer_region=(140,150))

        expected_clipped_cigar = MockCigarUtil(reference_start=42, cigar="75X")
        mock_cigar = MockCigarUtil(_softclip_target_return=expected_clipped_cigar)
        actual_clipped_cigar = primer_pair.softclip_primers(mock_cigar)
        self.assertEquals([(110, 140)], mock_cigar._softclip_targets_calls)
        self.assertEquals(expected_clipped_cigar, actual_clipped_cigar)
