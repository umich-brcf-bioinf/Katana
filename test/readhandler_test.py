#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
from __future__ import print_function, absolute_import
from ampliconsoftclipper import readhandler
import os.path
import pysam
from testfixtures.tempdirectory import TempDirectory
from test.util_test import ClipperBaseTestCase, MicroMock, MockRead,\
        MockPrimerPair, MockLog


class MockPrimerStatsDumper(MicroMock):
    def __init__(self, **kwargs):
        self._dump_calls=[]
        super(MockPrimerStatsDumper, self).__init__(**kwargs)

    def dump(self, primer_stats):
        self._dump_calls.append(primer_stats)


class MockPrimerStats(MicroMock):
    def __init__(self, **kwargs):
        self._add_read_primer_calls=[]
        super(MockPrimerStats, self).__init__(**kwargs)

    def add_read_primer(self, read, primer):
        self._add_read_primer_calls.append((read, primer))


class AddTagsReadHandlerTestCase(ClipperBaseTestCase):
    def test_handle(self):
        #pylint: disable=no-member
        original_reference_start = 100
        original_reference_end = 110
        original_cigar_string = "10M"
        read = MockRead(reference_start=original_reference_start,
                        reference_end=original_reference_end,
                        cigarstring=original_cigar_string)
        primer_pair_target = "target_1"
        primer_pair = MockPrimerPair(target_id=primer_pair_target)

        new_reference_start = 102
        new_cigar_string = "2S8M"
        transformation = (primer_pair,
                          new_reference_start,
                          new_cigar_string)
        mate_transformation = None
        handler = readhandler.AddTagsReadHandler()

        handler.handle(read, transformation, mate_transformation)

        self.assertEquals("X0:Z:" + primer_pair_target,
                          read._tags["X0"])
        self.assertEquals("X1:Z:" + original_cigar_string,
                          read._tags["X1"])
        self.assertEquals("X2:i:" + str(original_reference_start),
                          read._tags["X2"])
        self.assertEquals("X3:i:" + str(original_reference_end),
                          read._tags["X3"])


class ExcludeReadHandlerTestCase(ClipperBaseTestCase):
    def test_handle_readAndMateMatchPrimers(self):
        #pylint: disable=no-member
        read = MockRead(is_paired=True, mate_is_mapped=True, is_unmapped=False)
        primer_pair = MockPrimerPair(is_unmatched=False)
        transformation = (primer_pair, None, "1M")
        mate_transformation = (primer_pair, None, "1M")
        mock_log = MockLog()
        handler = readhandler.ExcludeNonMatchedReadHandler(log_method=mock_log)
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(True, read.mate_is_mapped)

    def test_handle_raisesIfUnmatchedPrimer(self):
        #pylint: disable=no-member
        read = MockRead(is_paired=True, mate_is_mapped=True, is_unmapped=False)
        primer_pair = MockPrimerPair(is_unmatched=True)
        transformation = (primer_pair, None, "1M")
        mate_transformation = (primer_pair, None, "1M")
        mock_log = MockLog()
        handler = readhandler.ExcludeNonMatchedReadHandler(log_method=mock_log)
        self.assertRaises(StopIteration,
                          handler.handle,
                          read,
                          transformation,
                          mate_transformation)

    def test_handle_unpairsReadIfMateUnmatchedPrimer(self):
        #pylint: disable=no-member
        read = MockRead(is_paired=True, is_unmapped=False, mate_is_mapped=True)
        primer_pair = MockPrimerPair(is_unmatched=False)
        mate_primer_pair = MockPrimerPair(is_unmatched=True)
        transformation = (primer_pair, None, "1M")
        mate_transformation = (mate_primer_pair, None, None)
        mock_log = MockLog()
        handler = readhandler.ExcludeNonMatchedReadHandler(log_method=mock_log)
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(False, read.is_paired)


class StatsReadHandlerTestCase(ClipperBaseTestCase):
    def test_handle(self):
        read1 = MockRead(is_positive_strand=True)
        read2 = MockRead(is_positive_strand=True)
        primer_pair = MockPrimerPair(target_id="target_1",
                                     chrom="chr1",
                                     sense_start=242)
        new_reference_start = 102
        new_cigar_string = "2S8M"
        transformation = (primer_pair,
                          new_reference_start,
                          new_cigar_string)
        mate_transformation = None
        mock_primer_stats = MockPrimerStats()
        mock_primer_stats_dumper = MockPrimerStatsDumper()
        handler = readhandler.StatsHandler(mock_primer_stats,
                                            mock_primer_stats_dumper)

        handler.handle(read1, transformation, mate_transformation)
        handler.handle(read2, transformation, mate_transformation)
        expected_calls = [(read1, primer_pair), (read2, primer_pair)]
        self.assertEquals(expected_calls,
                          mock_primer_stats._add_read_primer_calls)

        handler.end()
        self.assertEquals([mock_primer_stats],
                          mock_primer_stats_dumper._dump_calls)


class TransformReadHandlerTestCase(ClipperBaseTestCase):
    def test_handle(self):
        #pylint: disable=no-member
        read = MockRead(reference_start=100, cigarstring="10M")
        primer_pair = None
        new_reference_start = 102
        new_cigar_string = "2S8M"
        transformation = (primer_pair,
                          new_reference_start,
                          new_cigar_string)
        mate_transformation = (MockPrimerPair(is_unmatched=True),
                               0,
                               "")
        handler = readhandler.TransformReadHandler()
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(102, read.reference_start)
        self.assertEquals("2S8M", read.cigarstring)

    def test_handle_adjustsMateIfMatchedPrimerPair(self):
        #pylint: disable=no-member
        read = MockRead(mate_reference_start=100)
        transformation = (None, 0, "")
        mate_transformation = (MockPrimerPair(is_unmatched=False),
                               222,
                               "")
        handler = readhandler.TransformReadHandler()
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(222, read.next_reference_start)

    def test_handle_ignoresMateIfUnmatchedPrimerPair(self):
        #pylint: disable=no-member
        read = MockRead(mate_reference_start=100)
        transformation = (None, 0, "")
        mate_transformation = (MockPrimerPair(is_unmatched=True),
                               222,
                               "")
        handler = readhandler.TransformReadHandler()
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(100, read.mate_reference_start)


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
        readhandler.PYSAM_INDEX(filename)

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
        return MicroMock(aligned_segment=a)

    def test_end_sortsAndIndexes(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_bam_filename = os.path.join(input_dir.path, "input.bam")
            self.make_bam_file(input_bam_filename, [self.build_read()])
            output_bam_filename = os.path.join(output_dir.path, "output.bam")
            mock_log = MockLog()
            handler = readhandler.WriteReadHandler(input_bam_filename,
                                                    output_bam_filename,
                                                    log_method=mock_log.log)
            read1 = self.build_read(query_name="read1",
                                    reference_id=0,
                                    reference_start=20)
            read2 = self.build_read(query_name="read2",
                                    reference_id=0,
                                    reference_start=10)

            handler.begin()
            handler.handle(read1, None, None)
            handler.handle(read2, None, None)
            handler.end()

            actual_files = sorted(os.listdir(output_dir.path))
            self.assertEquals(["output.bam", "output.bam.bai"], actual_files)
            actual_bam = pysam.AlignmentFile(output_bam_filename, "rb")
            actual_reads = [read for read in actual_bam.fetch()]
            actual_bam.close()

        self.assertEquals(2, len(actual_reads))
        self.assertEquals("read2", actual_reads[0].query_name)
        self.assertEquals("read1", actual_reads[1].query_name)

