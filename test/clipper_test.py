#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
from __future__ import print_function, absolute_import
from ampliconsoftclipper import clipper
from ampliconsoftclipper import readhandler
from test.util_test import ClipperBaseTestCase, MockRead, MockReadHandler
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
                           readhandler.ExcludeNonMatchedReadHandler,
                           readhandler.AddTagsReadHandler,
                           readhandler.TransformReadHandler,
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
        self.assertEquals((clipper.PrimerPair.NULL_PRIMER_PAIR, 333, "10M"),
                          actual_read_transforms[333])

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
        null_transform = (clipper.PrimerPair.NULL_PRIMER_PAIR, 0, "")
        self.assertEquals(null_transform, actual_calls[2])

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

