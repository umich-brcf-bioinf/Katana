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
        c = cigar.cigar_editor_factory(42, "5M" "5I" "5M" "5S" "5H")
        self.assertEquals(56, c.reference_end)
        c = cigar.cigar_editor_factory(42, "5S" "5M" "5I" "3D" "5M" "5S")
        self.assertEquals(64, c.reference_end)

    def test_cigar_profile(self):
        c = cigar.cigar_editor_factory(42, "5M1I4S")
        self.assertEquals("MMMMMISSSS", c.cigar_profile)

    def test_cigar_profile_longer(self):
        c = cigar.cigar_editor_factory(42, "10M1I1D10H")
        self.assertEquals("MMMMMMMMMMIDHHHHHHHHHH", c.cigar_profile)

    def test_get_profile_overlap_tuple(self):
        c = cigar.cigar_editor_factory(42, "10M")
        self.assertEquals(("","MMMMMMMMMM",""),
                          c._get_profile_overlap_tuple(42,52))

    def test_get_profile_overlap_refConsumingOps(self):
        c = cigar.cigar_editor_factory(42, "5M5X")
        self.assertEquals(("", "MMMMM", "XXXXX"),
                          c._get_profile_overlap_tuple(42,47))
        self.assertEquals(("MMMMM", "XXXXX", ""),
                          c._get_profile_overlap_tuple(47, 52))

    def test_get_profile_overlap_refNonconsumingOps(self):
        c = cigar.cigar_editor_factory(42, "5M4I1M")
        self.assertEquals(("", "MMMMM", "IIIIM"),
                          c._get_profile_overlap_tuple(42,47))
        self.assertEquals(("MMMMM", "IIIIM", ""),
                          c._get_profile_overlap_tuple(47, 48))

    def test_get_profile_overlap_refNonconsumingOpsWithFlankingSoftclips(self):
        c = cigar.cigar_editor_factory(42, "5S" "5M" "4I" "1M" "5S")
        self.assertEquals(("", "SSSSS", "MMMMMIIIIMSSSSS"),
                          c._get_profile_overlap_tuple(37,42))
        self.assertEquals(("SSSSS", "MMMMM", "IIIIMSSSSS"),
                          c._get_profile_overlap_tuple(42,47))
        self.assertEquals(("SSSSSMMMMM","IIIIM", "SSSSS"),
                          c._get_profile_overlap_tuple(47, 48))
        self.assertEquals(("SSSSSMMMMMIIIIM", "SSSSS", ""),
                          c._get_profile_overlap_tuple(48,53))

    def test_get_profile_overlap_outOfBounds(self):
        c = cigar.cigar_editor_factory(42, "2S" "6M" "2X")
        self.assertEquals(("", "SS","MMMMMMXX"),
                          c._get_profile_overlap_tuple(30,42))
        self.assertEquals(("SSMMMMMM", "XX",""),
                          c._get_profile_overlap_tuple(48,60))
        self.assertEquals(("","","SSMMMMMMXX"),
                          c._get_profile_overlap_tuple(20,30))
        self.assertEquals(("SSMMMMMMXX", "", ""),
                          c._get_profile_overlap_tuple(100,110))

    def test_softclip_target(self):
        c = cigar.cigar_editor_factory(42, "10M")
        new_cigar_editor = c.softclip_target(44,50)
        self.assertEquals("2S6M2S", new_cigar_editor.cigar)
        self.assertEquals(44, new_cigar_editor.reference_start)

    def test_softclip_target_flankingSoftclips(self):
        c = cigar.cigar_editor_factory(42, "2S" "6M" "2S")
        #4444444444
        #0123456789
        #SSMMMMMMSS
        #SSSMMMMSSS
        new_cigar_editor = c.softclip_target(43,47)
        self.assertEquals("3S4M3S", new_cigar_editor.cigar)
        self.assertEquals(43, new_cigar_editor.reference_start)

    def test_softclip_target_flankingHardclips(self):
        c = cigar.cigar_editor_factory(42, "2H1S" "4M" "1S2H")
        #3444444444
        #9012345678
        #HHSMMMMSHH
        #HHSSMMSSHH
        new_cigar_editor = c.softclip_target(43,45)
        self.assertEquals("2H2S" "2M" "2S2H", new_cigar_editor.cigar)
        self.assertEquals(43, new_cigar_editor.reference_start)

    def test_softclip_target_deletes(self):
        c = cigar.cigar_editor_factory(42, "2M3D" "4M" "1S")
        #4444444455
        #2345678901
        #MMDDD
        #     MMMM
        #         S
        #SS   MMMMS
        new_cigar_editor = c.softclip_target(47,51)
        self.assertEquals("2S" "4M" "1S", new_cigar_editor.cigar)
        self.assertEquals(47, new_cigar_editor.reference_start)

#    TODO the edge adjustment in softclipping doesn't know that inserts don't take a genomic space (fix this)
    def test_softclip_target_edgeInsert(self):
        c = cigar.cigar_editor_factory(42, "3M" "1I4M" "2X")
        #444 444445
        #234 567890
        #ATAAACGTAC
        #MMMI
        #    MMMM
        #        XX
        #SSSSMMMMSS
        new_cigar_editor = c.softclip_target(45,49)
        self.assertEquals("4S" "4M" "2S", new_cigar_editor.cigar)
        self.assertEquals(45, new_cigar_editor.reference_start)

    def test_softclip_target_edgeDelete(self):
        c = cigar.cigar_editor_factory(42, "3M" "1D5M" "2S")
        #44444444555
        #23456789012
        #ATA ACGTACG
        #MMM
        #   DMMMM
        #         MSS
        #SSS MMMMSSS
        new_cigar_editor = c.softclip_target(45,50)
        self.assertEquals("3S" "4M" "3S", new_cigar_editor.cigar)
        self.assertEquals(46, new_cigar_editor.reference_start)


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
