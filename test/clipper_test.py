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
        return clipper.Read(a)

    def test_handle(self):
        with TempDirectory() as tmp_dir:
            input_bam_filename = os.path.join(tmp_dir.path, "input.bam")
            output_bam_filename = os.path.join(tmp_dir.path, "output.bam")
            self.make_bam_file(input_bam_filename, [self.build_read()])
            handler = clipper._WriteReadHandler(input_bam_filename,
                                               output_bam_filename)
            read1 = self.build_read(query_name="read1")
            read2 = self.build_read(query_name="read2")

            handler.begin()
            handler.handle(read1)
            handler.handle(read2)
            handler.end()

            clipper.PYSAM_INDEX(output_bam_filename)
            actual_bam = pysam.AlignmentFile(output_bam_filename, "rb")
            actual_reads = [read for read in actual_bam.fetch()]
            actual_bam.close()

        self.assertEquals(2, len(actual_reads))
        self.assertEquals("read1", actual_reads[0].query_name)
        self.assertEquals("read2", actual_reads[1].query_name)

class ReadTestCase(ClipperBaseTestCase):
    def test_init(self):
        aligned_segment = MockRead(query_name="read1",
                                   reference_name="chr1",
                                   reference_start=100,
                                   reference_end=110,
                                   cigarstring="10M")
        read = clipper.Read(aligned_segment)
        self.assertEquals("chr1", read.reference_name)
        self.assertEquals(100, read.reference_start)
        self.assertEquals(110, read.reference_end)
        self.assertEquals("10M", read.cigarstring)

    def test_mutatorsPassThroughToAlignedSegment(self):
        aligned_segment = MockRead(query_name="read1",
                                   reference_name="chr1",
                                   reference_start=100,
                                   cigarstring="10M")
        read = clipper.Read(aligned_segment)
        read.reference_start = 142
        read.cigarstring = "10S"
        self.assertEquals(142, aligned_segment.__dict__['reference_start'])
        self.assertEquals("10S", aligned_segment.__dict__['cigarstring'])

    def test_is_positive(self):
        read = clipper.Read(MockRead(is_reverse=False))
        self.assertEquals(True, read.is_positive_strand)
        read = clipper.Read(MockRead(is_reverse=True))
        self.assertEquals(False, read.is_positive_strand)

    def test_key(self):
        aligned_segment = MockRead(query_name="read1",
                                   reference_name="chr1",
                                   reference_start=100,
                                   is_reverse=False)
        read = clipper.Read(aligned_segment)
        expected_key = ("read1", True, "chr1", 100)
        self.assertEquals(expected_key, read.key)

    def test_mate_key(self):
        aligned_segment = MockRead(query_name="read1",
                                   is_paired=True,
                                   mate_is_unmapped=False,
                                   mate_is_reverse=True,
                                   next_reference_name="chr2",
                                   next_reference_start=200)
        read = clipper.Read(aligned_segment)
        expected_key = ("read1", False, "chr2", 200)
        self.assertEquals(expected_key, read.mate_key)

    def test_mate_key_noneWhenNoMate(self):
        #pylint: disable=attribute-defined-outside-init
        aligned_segment = MockRead(query_name="read1",
                                   is_paired=True,
                                   mate_is_unmapped=False,
                                   mate_is_reverse=True,
                                   next_reference_name="chr2",
                                   next_reference_start=200)
        aligned_segment.is_paired=False
        self.assertEquals(None, clipper.Read(aligned_segment).mate_key)
        aligned_segment.is_paired=True

        aligned_segment.mate_is_unmapped=True
        self.assertEquals(None, clipper.Read(aligned_segment).mate_key)
        aligned_segment.mate_is_unmapped=True

    def test_set_tag(self):
        aligned_segment = MockRead()
        read = clipper.Read(aligned_segment)
        read.set_tag("name", "value", "type")
        self.assertEquals("name:type:value", aligned_segment._tags["name"])

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
