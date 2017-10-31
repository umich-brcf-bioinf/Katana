#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
from __future__ import print_function, absolute_import

import os.path

import pysam
from testfixtures.tempdirectory import TempDirectory

from katana import readhandler
from katana.readhandler import ExcludeNonMatchedReadHandler
from test.util_test import KatanaBaseTestCase, MicroMock, MockRead, \
    MockPrimerPair, MockLog, build_read, make_bam_file


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


class AddTagsReadHandlerTestCase(KatanaBaseTestCase):
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

        transformation = MicroMock(primer_pair=primer_pair,
                                   filters=())
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

    def test_handle_addsFiltersWhenPresent(self):
        #pylint: disable=no-member
        read = MockRead(reference_start=100,
                        reference_end=110,
                        cigarstring="10M")
        primer_pair_target = "target_1"
        primer_pair = MockPrimerPair(target_id=primer_pair_target)

        transformation = MicroMock(primer_pair=primer_pair,
                                   filters=("filter1", "filter2"))
        mate_transformation = None
        handler = readhandler.AddTagsReadHandler()

        handler.handle(read, transformation, mate_transformation)

        self.assertEquals("X4:Z:" + "filter1,filter2", read._tags["X4"])

    def test_handle_sanitizes_tags(self):
        #pylint: disable=no-member
        original_reference_start = 100
        original_reference_end = 110
        original_cigar_string = "10M"
        read = MockRead(reference_start=100,
                        reference_end=110,
                        cigarstring="10M")
        primer_pair_target = "123!@#$%target \t\t A \n thing-_."
        primer_pair = MockPrimerPair(target_id=primer_pair_target)

        transformation = MicroMock(primer_pair=primer_pair,
                                   filters=())
        mate_transformation = None
        handler = readhandler.AddTagsReadHandler()

        handler.handle(read, transformation, mate_transformation)

        self.assertEquals("X0:Z:123_target_A_thing-_.",
                          read._tags["X0"])


class ExcludeReadHandlerTestCase(KatanaBaseTestCase):
    def test_handle_noExceptionIfNoFilters(self):
        read = MockRead()
        transformation = MicroMock(filters=())
        mate_transformation = MicroMock(filters=())
        mock_log = MockLog()
        handler = ExcludeNonMatchedReadHandler(log_method=mock_log.log)
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(True, True)

    def test_handle_raisesIfTransformHasFilter(self):
        #pylint: disable=no-member
        read = MockRead()
        transformation = MicroMock(filters=("filterA"))
        mate_transformation = MicroMock(filters=())
        mock_log = MockLog()
        handler = ExcludeNonMatchedReadHandler(log_method=mock_log.log)
        self.assertRaises(StopIteration,
                          handler.handle,
                          read,
                          transformation,
                          mate_transformation)

    def test_handle_logsExclusions(self):
        #pylint: disable=no-member
        read = MockRead()
        mate_transform = MicroMock(filters=())
        mock_log = MockLog()
        handler = ExcludeNonMatchedReadHandler(log_method=mock_log.log)
        try:
            handler.handle(read, MicroMock(filters=("filterA")), mate_transform)
        except StopIteration:
            pass
        try:
            handler.handle(read, MicroMock(filters=("filterA")), mate_transform)
        except StopIteration:
            pass
        try:
            handler.handle(read, MicroMock(filters=("filterB")), mate_transform)
        except StopIteration:
            pass
        self.assertEquals(2, len(handler._all_exclusions))
        self.assertEquals(2, handler._all_exclusions["filterA"])
        self.assertEquals(1, handler._all_exclusions["filterB"])
        handler.end()
        self.assertEquals(2, len(mock_log._log_calls))

    def test_handle_unpairsReadIfMateHasFilter(self):
        #pylint: disable=no-member
        read = MockRead()
        transformation = MicroMock(filters=())
        mate_transformation = MicroMock(filters=("filterA"))
        mock_log = MockLog()
        handler = readhandler.ExcludeNonMatchedReadHandler(log_method=mock_log)
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(False, read.is_paired)


class StatsReadHandlerTestCase(KatanaBaseTestCase):
    def test_handle(self):
        read1 = MockRead(is_positive_strand=True)
        read2 = MockRead(is_positive_strand=True)
        primer_pair = MockPrimerPair(target_id="target_1",
                                     chrom="chr1",
                                     sense_start=242)
        new_reference_start = 102
        new_cigar_string = "2S8M"
        transformation = MicroMock(primer_pair=primer_pair,
                                   reference_start=new_reference_start,
                                   cigar=new_cigar_string)
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


class TransformReadHandlerTestCase(KatanaBaseTestCase):
    def test_handle(self):
        #pylint: disable=no-member
        read = MockRead(reference_start=100,
                        cigarstring="10M",
                        next_reference_start=150,
                        is_paired = True,
                        mate_cigar='10M')
        new_reference_start = 102
        new_cigar_string = "2S8M"
        new_mate_cigar_string = '8M2S'
        new_next_ref_start = 200
        transformation = MicroMock(primer_pair=None,
                                   reference_start=new_reference_start,
                                   cigar=new_cigar_string)
        mate_transformation = MicroMock(primer_pair=None,
                                        reference_start=new_next_ref_start,
                                        cigar=new_mate_cigar_string,
                                        is_unmapped=False)
        handler = readhandler.TransformReadHandler()
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(new_reference_start, read.reference_start)
        self.assertEquals(new_cigar_string, read.cigarstring)
        self.assertEquals(new_next_ref_start, read.next_reference_start)
        self.assertEquals(new_mate_cigar_string, read.mate_cigar)

    def test_handle_ignoresMateIfUnpaired(self):
        #pylint: disable=no-member
        reference_start = 150
        read = MockRead(reference_start=100,
                        cigarstring="10M",
                        next_reference_start=reference_start,
                        is_paired = False)
        new_reference_start = 102
        new_cigar_string = "2S8M"
        transformation = MicroMock(primer_pair=None,
                                   reference_start=new_reference_start,
                                   cigar=new_cigar_string)
        mate_transformation = None
        handler = readhandler.TransformReadHandler()
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(reference_start, read.next_reference_start)

    def test_handle_doesNotResetNextStartIfMateUnmapped(self):
        #pylint: disable=no-member
        reference_start = 111
        original_next_reference_start = 333
        read = MockRead(reference_start=reference_start,
                        cigarstring="10M",
                        next_reference_start=original_next_reference_start,
                        is_paired = True)
        new_reference_start = 222
        transformation = MicroMock(reference_start=new_reference_start,
                                   cigar="10M")
        mate_transformation = MicroMock(is_unmapped=True,
                                        reference_start=333)
        handler = readhandler.TransformReadHandler()
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(original_next_reference_start,
                          read.next_reference_start)

    def test_handle_RemovesMateCigarTagIfMateUnmapped(self):
        #pylint: disable=no-member
        read = MockRead(reference_start=111,
                        cigarstring="10M",
                        next_reference_start=222,
                        is_paired = True,
                        mate_cigar='10M')
        transformation = MicroMock(reference_start=111,
                                   cigar="10M")
        mate_transformation = MicroMock(is_unmapped=True,
                                        reference_start=333)
        handler = readhandler.TransformReadHandler()
        handler.handle(read, transformation, mate_transformation)
        self.assertEquals(None, read.mate_cigar)


class WriteReadHandlerTestCase(KatanaBaseTestCase):

    def test_end_sortsAndIndexes(self):
        #pylint: disable=no-member
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_bam_filename = os.path.join(input_dir.path, "input.bam")
            make_bam_file(input_bam_filename, [build_read()])
            output_bam_filename = os.path.join(output_dir.path, "output.bam")
            mock_log = MockLog()
            handler = readhandler.WriteReadHandler(input_bam_filename,
                                                    output_bam_filename,
                                                    log_method=mock_log.log)
            read1 = build_read(query_name="read1",
                               reference_id=0,
                               reference_start=20)
            read2 = build_read(query_name="read2",
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

