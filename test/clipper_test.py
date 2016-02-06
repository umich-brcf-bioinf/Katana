#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
from __future__ import print_function, absolute_import
from katana import clipper
from katana import readhandler
from test.util_test import ClipperBaseTestCase, MockPrimerPair, MockRead, \
        MockReadHandler, MockCigarUtil, MicroMock
import katana.util as util
from katana.util import ReadTransformation
from katana.clipper import _build_read_transformations
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class ClipperTestCase(ClipperBaseTestCase):
    def test_build_handlers_excludeUnmatchedReads(self):
        actual_handlers = clipper._build_handlers("input_bam_filename",
                                                  "output_bam_filename",
                                                  False)
        actual_handler_classes = [x.__class__ for x in actual_handlers]
        self.assertEquals([readhandler.StatsHandler,
                           readhandler.AddTagsReadHandler,
                           readhandler.TransformReadHandler,
                           readhandler.ExcludeNonMatchedReadHandler,
                           readhandler.WriteReadHandler],
                          actual_handler_classes)

    def test_build_handlers_includeUnmatchedReads(self):
        actual_handlers = clipper._build_handlers("input_bam_filename",
                                                  "output_bam_filename",
                                                  True)
        actual_handler_classes = [x.__class__ for x in actual_handlers]
        self.assertEquals([readhandler.StatsHandler,
                           readhandler.AddTagsReadHandler,
                           readhandler.TransformReadHandler,
                           readhandler.WriteReadHandler],
                          actual_handler_classes)

    def test_filter_builder_noFilters(self):
        transform = MicroMock(primer_pair=MockPrimerPair(is_unmatched=False),
                              is_unmapped=False,
                              is_cigar_valid=True)
        self.assertEquals([], clipper._filter_builder(transform))

    def test_filter_builder_unmatchedPrimerPair(self):
        transform = MicroMock(primer_pair=MockPrimerPair(is_unmatched=True),
                              is_unmapped=False,
                              is_cigar_valid=True)
        self.assertEquals(["UNMATCHED_PRIMER_PAIR"],
                          clipper._filter_builder(transform))

    def test_filter_builder_unmappedAlignment(self):
        transform = MicroMock(primer_pair=MockPrimerPair(is_unmatched=False),
                              is_unmapped=True,
                              is_cigar_valid=True)
        self.assertEquals(["UNMAPPED_ALIGNMENT"],
                          clipper._filter_builder(transform))

    def test_filter_builder_invalidCigar(self):
        transform = MicroMock(primer_pair=MockPrimerPair(is_unmatched=False),
                              is_unmapped=False,
                              is_cigar_valid=False)
        self.assertEquals(["INVALID_CIGAR"],
                          clipper._filter_builder(transform))

    def test_filter_builder_multiple(self):
        transform = MicroMock(primer_pair=MockPrimerPair(is_unmatched=True),
                              is_unmapped=False,
                              is_cigar_valid=False)
        self.assertEquals(["UNMATCHED_PRIMER_PAIR","INVALID_CIGAR"],
                          clipper._filter_builder(transform))
        transform = MicroMock(primer_pair=MockPrimerPair(is_unmatched=True),
                              is_unmapped=True,
                              is_cigar_valid=False)
        self.assertEquals(["UNMAPPED_ALIGNMENT"],
                          clipper._filter_builder(transform))

    def test_build_read_transformations(self):
        primer_pair = clipper.PrimerPair(target_id="target_1",
                                         chrom="chr42",
                                         sense_primer_region=(100,102),
                                         antisense_primer_region=(148,150))
        read1 = MockRead(query_name="read1-sense",
                        is_positive_strand=True,
                        reference_name="chr42",
                        reference_start=100,
                        reference_end=150,
                        cigarstring="10M",
                        is_paired=True,
                        is_unmapped=False,
                        key=111)
        read2 = MockRead(query_name="read2-antisense",
                        is_positive_strand=False,
                        reference_name="chr42",
                        reference_start=140,
                        reference_end=150,
                        cigarstring="10M",
                        is_paired=False,
                        is_unmapped=False,
                        key=222)
        read3 = MockRead(query_name="read3",
                        is_positive_strand=True,
                        reference_name="chr42",
                        reference_start=333,
                        reference_end=343,
                        cigarstring="10M",
                        is_paired=True,
                        is_unmapped=False,
                        key=333)
        filter_builder = lambda transform: []

        read_iter = iter([read1, read2, read3])
        actual_read_transforms = _build_read_transformations(read_iter,
                                                             filter_builder)
        self.assertEquals(3, len(actual_read_transforms))
        filter_builder=lambda x: []
        transform1 = ReadTransformation(read1,
                                        primer_pair,
                                        MockCigarUtil(reference_start=102,
                                                      cigar="2S8M"),
                                        filter_builder)
        transform2 = ReadTransformation(read2,
                                        primer_pair,
                                        MockCigarUtil(reference_start=140,
                                                      cigar="8M2S"),
                                        filter_builder)
        transform3 = ReadTransformation(read3,
                                        clipper.PrimerPair.NULL,
                                        MockCigarUtil(reference_start=333,
                                                      cigar="10M"),
                                        filter_builder)
        self.assertEquals(transform1,
                          actual_read_transforms[111])
        self.assertEquals(transform2,
                          actual_read_transforms[222])
        self.assertEquals(transform3,
                          actual_read_transforms[333])

    def test_build_read_transformations2_addsFilters(self):
        clipper.PrimerPair(target_id="target_1",
                           chrom="chr42",
                           sense_primer_region=(100,102),
                           antisense_primer_region=(148,150))
        clipper.PrimerPair(target_id="target_2",
                           chrom="chr42",
                           sense_primer_region=(200,202),
                           antisense_primer_region=(248,250))
        read1 = MockRead(query_name="readA",
                        is_positive_strand=True,
                        reference_name="chr42",
                        reference_start=100,
                        reference_end=150,
                        cigarstring="50M",
                        is_paired=True,
                        is_unmapped=False,
                        key=111)
        read2 = MockRead(query_name="readB",
                        is_positive_strand=False,
                        reference_name="chr42",
                        reference_start=200,
                        reference_end=250,
                        cigarstring="50M",
                        is_paired=False,
                        is_unmapped=False,
                        key=222)
        read_iter = iter([read1, read2])
        filter_builder = lambda transform: [str(transform.reference_start),
                                            transform.primer_pair.target_id]

        actual_read_transforms = _build_read_transformations(read_iter,
                                                             filter_builder)
        self.assertEquals(2, len(actual_read_transforms))
        self.assertEquals(("102", "target_1"),
                          actual_read_transforms[111].filters)
        self.assertEquals(("202", "target_2"),
                          actual_read_transforms[222].filters)



    def test_build_read_transformations2_reraisesExceptions(self):
        clipper.PrimerPair(target_id="target_1",
                           chrom="chr42",
                           sense_primer_region=(100,102),
                           antisense_primer_region=(148,150))
        read1 = MockRead(query_name="readA-sense",
                        is_positive_strand=True,
                        reference_name="chr42",
                        reference_start=100,
                        reference_end=150,
                        cigarstring="10M",
                        is_paired=True,
                        is_unmapped=False,
                        key=111)
        read2 = MockRead(query_name="readB-antisense",
                        key=222)
        read_iter = iter([read1, read2])
        filter_builder = lambda transform: []

        self.assertRaisesRegexp(util.ClipperException,
                                (r"Problem with read readB-antisense "
                                 r"\[line 2\] and primer pair target_1: .*"),
                                clipper._build_read_transformations,
                                read_iter,
                                filter_builder)

    def test_handle_reads(self):
        handler1 = MockReadHandler()
        handler2 = MockReadHandler()
        handler3 = MockReadHandler()
        handlers = [handler1, handler2, handler3]
        read1 = MockRead(key=1, mate_key=11)
        read2 = MockRead(key=2, mate_key=12)
        read3 = MockRead(key=3, mate_key=13)
        transform1 = ()
        mate_transform1=()
        transform2 = ()
        mate_transform2=()
        transform3 = ()
        mate_transform3=()

        read_transformations = {1: transform1,
                                2: transform2,
                                3:transform3,
                                11: mate_transform1,
                                12: mate_transform2,
                                13: mate_transform3}
        read_iter = iter([read1, read2, read3])
        clipper._handle_reads(handlers, read_iter, read_transformations)
        self.assertEquals(1, handler1.begin_calls)
        self.assertEquals(1, handler2.begin_calls)
        self.assertEquals(1, handler3.begin_calls)

        expected_calls = [(read1, transform1, mate_transform1),
                          (read2, transform2, mate_transform2),
                          (read3, transform3, mate_transform3)]
        self.assertEquals(expected_calls, handler1.handle_calls)
        self.assertEquals(expected_calls, handler2.handle_calls)
        self.assertEquals(expected_calls, handler3.handle_calls)
        self.assertEquals(1, handler1.end_calls)
        self.assertEquals(1, handler2.end_calls)
        self.assertEquals(1, handler3.end_calls)

    def test_handle_reads_stopIterationSkipsHandlers(self):
        read1 = MockRead(key=1, mate_key=11)
        read2 = MockRead(key=2, mate_key=12)
        read3 = MockRead(key=3, mate_key=13)
        handler1 = MockReadHandler()
        handler2 = MockReadHandler(handle_raise=lambda read,_: read.key==2)
        handler3 = MockReadHandler()
        handlers = [handler1, handler2, handler3]
        transform1 = ()
        mate_transform1=()
        transform2 = ()
        mate_transform2=()
        transform3 = ()
        mate_transform3=()
        read_transformations = {1: transform1,
                                2: transform2,
                                3:transform3,
                                11: mate_transform1,
                                12: mate_transform2,
                                13: mate_transform3}
        read_iter = iter([read1, read2, read3])

        clipper._handle_reads(handlers, read_iter, read_transformations)
        self.assertEquals(1, handler1.begin_calls)
        self.assertEquals(1, handler2.begin_calls)
        self.assertEquals(1, handler3.begin_calls)

        expected_calls12 = [(read1, transform1, mate_transform1),
                            (read2, transform2, mate_transform2),
                            (read3, transform3, mate_transform3)]
        expected_calls3 = [(read1, transform1, mate_transform1),
                           (read3, transform3, mate_transform3)]
        self.assertEquals(expected_calls12, handler1.handle_calls)
        self.assertEquals(expected_calls12, handler2.handle_calls)
        self.assertEquals(expected_calls3, handler3.handle_calls)
        self.assertEquals(1, handler1.end_calls)
        self.assertEquals(1, handler2.end_calls)
        self.assertEquals(1, handler3.end_calls)


    def test_handle_reads_unmatedReadOk(self):
        read1 = MockRead(key=1, mate_key=None)
        handler1 = MockReadHandler()
        handlers = [handler1]
        transform1 = ()
        read_transformations = {1: transform1}
        read_iter = iter([read1])

        clipper._handle_reads(handlers, read_iter, read_transformations)
        self.assertEquals(1, handler1.begin_calls)

        self.assertEquals(1, len(handler1.handle_calls))
        actual_calls = handler1.handle_calls[0]
        self.assertEquals(read1, actual_calls[0])
        self.assertEquals(transform1, actual_calls[1])
        self.assertEquals(ReadTransformation.NULL, actual_calls[2])

        self.assertEquals(1, handler1.end_calls)


    def test_initialize_primer_pairs(self):
        clipper.PrimerPair._all_primers = {}
        #pylint: disable=line-too-long
        input_manifest = (\
'''id|Customer TargetID|Chr|Genome Build|Sense Start|Antisense Start|Sense Sequence|Antisense Sequence
1|NRAS_3|1|hg19 dbSNP137|111|222|ATCCGACAA|AGTACAAACT
2|IDH1_1|2|hg19 dbSNP137|333|444|TGCCAACAT|TGGAAATCAC
''').replace("|", "\t")
        clipper._initialize_primer_pairs(StringIO(input_manifest))
        actual_primers = clipper.PrimerPair._all_primers
        self.assertEquals(4, len(actual_primers))
        targets = set([primer.target_id for primer in actual_primers.values()])
        self.assertEquals(set(["IDH1_1", "NRAS_3"]), targets)

