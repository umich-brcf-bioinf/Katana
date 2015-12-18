#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
from __future__ import print_function, absolute_import
import unittest
from ampliconsoftclipper.cigar import IndelException
from ampliconsoftclipper import cigar

class CigarEditorTestCase(unittest.TestCase):
    def test_init(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertEquals("10M", c.cigar)
        self.assertEquals(42, c.reference_start)

    def test_init_reference_end(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertEquals(51, c.reference_end)
        c = cigar.cigar_editor_factory(42, "5M" "5I" "5S" "5H" "5M")
        self.assertEquals(51, c.reference_end)
        c = cigar.cigar_editor_factory(42, "5S" "5M" "5I" "3D" "5M" "5S")
        self.assertEquals(54, c.reference_end)

    def test_cigar_profile(self):
        c = cigar.cigar_editor_factory(42, "5M1I4S")
        self.assertEquals("MMMMMISSSS", c.cigar_profile)

    def test_cigar_profile_longer(self):
        c = cigar.cigar_editor_factory(42, "10M1I1D10H")
        self.assertEquals("MMMMMMMMMMIDHHHHHHHHHH", c.cigar_profile)

    def test_get_profile_overlap(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertEquals("MMMMMMMMMM", c._get_profile_overlap(42,52))

    def test_get_profile_overlap_refConsumingOps(self):
        c = cigar.cigar_editor_factory(42, "5M5X")
        self.assertEquals("MMMMM", c._get_profile_overlap(42,47))
        self.assertEquals("XXXXX", c._get_profile_overlap(47, 52))

    def test_get_profile_overlap_refNonconsumingOps(self):
        c = cigar.cigar_editor_factory(42, "5M4I1M")
        self.assertEquals("MMMMM", c._get_profile_overlap(42,47))
        self.assertEquals("IIIIM", c._get_profile_overlap(47, 48))

    def test_get_profile_overlap_refNonconsumingOpsWithFlankingSoftclips(self):
        c = cigar.cigar_editor_factory(42, "5S5M4I1M5S")
        self.assertEquals("SSSSS", c._get_profile_overlap(37,42))
        self.assertEquals("MMMMM", c._get_profile_overlap(42,47))
        self.assertEquals("IIIIM", c._get_profile_overlap(47, 48))
        self.assertEquals("SSSSS", c._get_profile_overlap(48,53))

    def test_get_profile_overlap_outOfBounds(self):
        c = cigar.cigar_editor_factory(42, "2S6M2X")
        self.assertEquals("SS", c._get_profile_overlap(30,42))
        self.assertEquals("XX", c._get_profile_overlap(48,60))
        self.assertEquals("", c._get_profile_overlap(20,30))

#     def test_partition(self):
#         c = cigar.cigar_editor_factory(42, "10M")
#         self.assertEquals((42,"2M"), (44,"6M"), (50, "2M"), c.partition(44,50))

    def test_indels_in_region_falseIfNoIndels(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertEquals(False, c.indels_in_region(42,52))

    def test_indels_in_region_true(self):
        c = cigar.cigar_editor_factory(42, "5M5I")
        self.assertEquals(True, c.indels_in_region(42,48))

    def test_indels_in_region_falseOutsideRegion(self):
        c = cigar.cigar_editor_factory(42, "5M5I")
        self.assertEquals(False, c.indels_in_region(42,47))

#     def test_indels_in_region_falseRegionBeforeSequence(self):
#         c = cigar.cigar_editor_factory(42, "5M5I")
#         self.assertEquals(False, c.indels_in_region(30,40))

    def test_indels_in_region_falseRegionAfterSequence(self):
        c = cigar.cigar_editor_factory(42, "5M5I")
        self.assertEquals(False, c.indels_in_region(53,63))

    def test_indels_in_region_falseRegionAfterComplexSequence(self):
        c = cigar.cigar_editor_factory(42, "5M5I5M")
        self.assertEquals(False, c.indels_in_region(53,63))

    def test_indels_in_region_falseRegionWithLeadingSoftclips(self):
        c = cigar.cigar_editor_factory(42, "5S5M5I")
        self.assertEquals(False, c.indels_in_region(42,47))
        self.assertEquals(True, c.indels_in_region(47,48))

    def test_indels_in_front_falseIfNoIndels(self):
        c = cigar.cigar_editor_factory(42, "5M5S")
        self.assertEquals(False, c.edge_indels('front', length=10))

    def test_indels_in_front_trueIfInsert(self):
        c = cigar.cigar_editor_factory(42, "5M1I4M")
        self.assertEquals(True, c.edge_indels('front',length=10))

    def test_indels_in_front_falseIfInsertOutsideRange(self):
        c = cigar.cigar_editor_factory(42, "5M1I4M")
        self.assertEquals(False, c.edge_indels('front',length=5))

    def test_indels_in_front_trueIfDelete(self):
        c = cigar.cigar_editor_factory(42, "5M1D4M")
        self.assertEquals(True, c.edge_indels('front',length=10))

    def test_indels_in_front_falseIfDeleteOutsideRange(self):
        c = cigar.cigar_editor_factory(42, "5M1D4M")
        self.assertEquals(False, c.edge_indels('front',length=5))

    def test_indels_in_front_outOfBounds(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertRaises(ValueError, c.edge_indels, 'front', 15)

    def test_indels_in_back_falseIfNoIndels(self):
        c = cigar.cigar_editor_factory(42, "5M5S")
        self.assertEquals(False, c.edge_indels('back',length=10))

    def test_indels_in_back_trueIfInsert(self):
        c = cigar.cigar_editor_factory(42, "4M1I5M")
        self.assertEquals(True, c.edge_indels('back',length=10))

    def test_indels_in_back_falseIfInsertOutsideRange(self):
        c = cigar.cigar_editor_factory(42, "4M1I5M")
        self.assertEquals(False, c.edge_indels('back',length=5))

    def test_indels_in_back_trueIfDelete(self):
        c = cigar.cigar_editor_factory(42, "4M1D5M")
        self.assertEquals(True, c.edge_indels('back',length=10))

    def test_indels_in_back_falseIfDeleteOutsideRange(self):
        c = cigar.cigar_editor_factory(42, "4M1D5M")
        self.assertEquals(False, c.edge_indels('back',length=5))

    def test_indels_in_back_outOfBounds(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertRaises(ValueError, c.edge_indels, 'back', 15)

    def test_softclip_front(self):
        c = cigar.cigar_editor_factory(42, "10M")
        new_cigar = c.softclip_front(length=5)
        self.assertEquals("5S5M", new_cigar.cigar)
        self.assertEquals(47, new_cigar.reference_start)

    def test_softclip_front_existingSoftClips(self):
        c = cigar.cigar_editor_factory(42, "3S7M")
        new_cigar = c.softclip_front(length=5)
        self.assertEquals("5S5M", new_cigar.cigar)
        self.assertEquals(44, new_cigar.reference_start)

    def test_softclip_front_existingHardClips(self):
        c = cigar.cigar_editor_factory(42, "3H7M")
        new_cigar = c.softclip_front(length=5)
        self.assertEquals("3H2S5M", new_cigar.cigar)
        self.assertEquals(44, new_cigar.reference_start)

    def test_softclip_front_withIndel(self):
        c = cigar.cigar_editor_factory(42, "3M2I5M")
        self.assertRaisesRegexp(IndelException,
                                (r"Indel found in first \[5\] bases of "
                                 r"cigar \[3M2I5M\]"),
                                c.softclip_front,
                                5)
        c = cigar.cigar_editor_factory(42, "3M2D5M")
        self.assertRaisesRegexp(IndelException,
                                (r"Indel found in first \[5\] bases of "
                                 r"cigar \[3M2D5M\]"),
                                c.softclip_front,
                                5)

    def test_softclip_front_outOfBounds(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertRaises(ValueError, c.softclip_front, 15)

    def test_softclip_back(self):
        c = cigar.cigar_editor_factory(42, "10M")
        new_cigar = c.softclip_back(length=5)
        self.assertEquals("5M5S", new_cigar.cigar)
        self.assertEquals(42, new_cigar.reference_start)

    def test_softclip_back_existingSoftClips(self):
        c = cigar.cigar_editor_factory(42, "7M3S")
        new_cigar = c.softclip_back(length=5)
        self.assertEquals("5M5S", new_cigar.cigar)
        self.assertEquals(42, new_cigar.reference_start)

    def test_softclip_back_existingHardClips(self):
        c = cigar.cigar_editor_factory(42, "7M3H")
        new_cigar = c.softclip_back(length=5)
        self.assertEquals("5M2S3H", new_cigar.cigar)
        self.assertEquals(42, new_cigar.reference_start)

    def test_softclip_back_withIndel(self):
        c = cigar.cigar_editor_factory(42, "5M2I3M")
        self.assertRaisesRegexp(IndelException,
                                (r"Indel found in last \[5\] bases of "
                                 r"cigar \[5M2I3M\]"),
                                c.softclip_back,
                                5)
        c = cigar.cigar_editor_factory(42, "5M2D3M")
        self.assertRaisesRegexp(IndelException,
                                (r"Indel found in last \[5\] bases of "
                                 r"cigar \[5M2D3M\]"),
                                c.softclip_back,
                                5)

    def test_softclip_back_outOfBounds(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertRaises(ValueError, c.softclip_back, 15)

    def test_null_cigar_false(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertEquals(False, c.is_null())

    def test_null_cigar(self):
        c = cigar.cigar_editor_factory(42, "*")
        self.assertEquals("*", c.cigar)
        self.assertEquals(False, c.edge_indels('front',42))
        self.assertEquals(False, c.edge_indels('back',42))
        self.assertEquals("*", c.softclip_front(42).cigar)
        self.assertEquals("*", c.softclip_back(42).cigar)
        self.assertEquals(True, c.is_null())
    #test malformed cigar
    #test softclips eliminate all matches
