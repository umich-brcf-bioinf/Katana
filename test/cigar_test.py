#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
from __future__ import print_function, absolute_import

import re
import unittest

from katana import cigar
import katana.util
from test.util_test import MicroMock


class CigarUtilTestCase(unittest.TestCase):
    def test_assert_query_lengths_match(self):
        base = cigar.CigarUtil(42, "10M")
        trivially_matching = cigar.CigarUtil(42, "10M")
        base._assert_query_lengths_match(trivially_matching)

        subtly_matching = cigar.CigarUtil(42, "2M15N8S")
        base._assert_query_lengths_match(subtly_matching)

        not_matching = cigar.CigarUtil(42, "2S8M2S")
        self.assertRaisesRegexp(katana.util.KatanaException,
                               r"Old CIGAR query length.*10M.*10.*2S8M2S.*12.*",
                               base._assert_query_lengths_match,
                               not_matching)

    def test_ref_consuming(self):
        util = cigar.CigarUtil(0, "1M")
        self.assertEquals(True, util._is_ref_consuming("M"))
        self.assertEquals(False, util._is_ref_consuming("I"))
        self.assertEquals(True, util._is_ref_consuming("S"))
        self.assertEquals(True, util._is_ref_consuming("="))
        self.assertEquals(True, util._is_ref_consuming("X"))
        self.assertEquals(True, util._is_ref_consuming("D"))
        self.assertEquals(True, util._is_ref_consuming("N"))
        self.assertEquals(False, util._is_ref_consuming("H"))
        self.assertEquals(False, util._is_ref_consuming("P"))

    def test_init(self):
        util = cigar.CigarUtil(42, "10M")
        self.assertEquals("10M", util.cigar)
        self.assertEquals(42, util.reference_start)

    def test_init_valid(self):
        util = cigar.CigarUtil(42, "1M1I1D1N1=1X" "1H1S")
        self.assertEquals(True, util.is_valid)

    def test_init_invalid(self):
        util = cigar.CigarUtil(42, "1H1S")
        self.assertEquals(False, util.is_valid)

    def test_init_withProfile(self):
        util = cigar.CigarUtil(42, cigar_profile="XXMMM")
        self.assertEquals("2X3M", util.cigar)
        self.assertEquals("XXMMM", util.cigar_profile)
        self.assertEquals(42, util.reference_start)

    def test_init_inconsistentCigarProfileDoesNotFail(self):
        util = cigar.CigarUtil(42, cigar="10M", cigar_profile="XXMMM")
        self.assertEquals("10M", util.cigar)
        self.assertEquals("XXMMM", util.cigar_profile)
        self.assertEquals(42, util.reference_start)

    def test_init_with_neitherCigarNorProfile(self):
        util = cigar.CigarUtil(42, cigar="", cigar_profile="")
        self.assertEquals("", util.cigar)
        self.assertEquals("", util.cigar_profile)
        self.assertEquals(42, util.reference_start)
        util = cigar.CigarUtil(42, cigar=None, cigar_profile=None)
        self.assertEquals("", util.cigar)
        self.assertEquals("", util.cigar_profile)
        self.assertEquals(42, util.reference_start)
        util = cigar.CigarUtil(42)
        self.assertEquals("", util.cigar)
        self.assertEquals("", util.cigar_profile)
        self.assertEquals(42, util.reference_start)

    def test_eq(self):
        base = cigar.CigarUtil(42, "10M")
        self.assertEquals(base, cigar.CigarUtil(42, "10M"))
        self.assertNotEquals(base, cigar.CigarUtil(10, "10M"))
        self.assertNotEquals(base, cigar.CigarUtil(42, "1M"))

    def test_repr(self):
        base = cigar.CigarUtil(42, "10M")
        match = re.match("<.*(cigar.*?)'>(.*)", repr(base))
        eval_repr = match.group(1) + match.group(2)
        #pylint: disable=eval-used
        self.assertEquals(base, eval(eval_repr))

    def test_cigar_profile(self):
        util = cigar.CigarUtil(42, "5M1I4S")
        self.assertEquals("MMMMMISSSS", util.cigar_profile)

    def test_cigar_profile_longer(self):
        util = cigar.CigarUtil(42, "10M1I1D10H")
        self.assertEquals("MMMMMMMMMMIDHHHHHHHHHH", util.cigar_profile)

    def test_softclip_replacesConsumingOpsWithSoftClips(self):
        util = cigar.CigarUtil(0, "1M")
        self.assertEquals("HS" "SSSS" "SH",
                          util._softclip(cigar_profile="HS" "MI=X" "SH"))

    def test_softclip_removesNonConsumingOps(self):
        util = cigar.CigarUtil(0, "1M")
        self.assertEquals("HS" "SH",
                          util._softclip(cigar_profile="HS" "DP" "SH"))

    def test_softclip_to_first_match(self):
        util = cigar.CigarUtil(0, "1M")
        self.assertEquals((45, "MM"),
                          util._softclip_to_first_match(45, "MM"))

    def test_softclip_to_first_noMatchesPassesThrough(self):
        util = cigar.CigarUtil(0, "1M")
        self.assertEquals((45, "SS"),
                          util._softclip_to_first_match(45, "SS"))

    def test_softclip_to_first_match_posIsAdjusted(self):
        util = cigar.CigarUtil(0, "1M")
        self.assertEquals((45, "SS" "MPNDISH"),
                          util._softclip_to_first_match(42, "SDNIP" "MPNDISH"))

    def test_softclip_to_first_match_uncommonMatchOpOk(self):
        util = cigar.CigarUtil(0, "1M")
        self.assertEquals((45, "SS" "XPNDISH"),
                          util._softclip_to_first_match(42, "SDNIP" "XPNDISH"))
        util = cigar.CigarUtil(0, "1M")
        self.assertEquals((45, "SS" "=PNDISH"),
                          util._softclip_to_first_match(42, "SDNIP" "=PNDISH"))

    def test_pos_profiles(self):
        util = cigar.CigarUtil(0, "1M")
        actual = util._pos_profiles("MMM")
        self.assertEquals((0, ["M", "M", "M"]), actual)

    def test_pos_profiles_uncommonMatchOpOk(self):
        util = cigar.CigarUtil(0, "1M")
        actual = util._pos_profiles("=MM")
        self.assertEquals((0, ["=", "M", "M"]), actual)

    def test_pos_profiles_leadingSoftclips(self):
        util = cigar.CigarUtil(0, "1M")
        actual = util._pos_profiles("SMM")
        self.assertEquals((1, ["S", "M", "M"]), actual)

    def test_pos_profiles_insertOk(self):
        util = cigar.CigarUtil(0, "1M")
        actual = util._pos_profiles("SIM")
        self.assertEquals((1, ["S", "IM"]), actual)

    def test_pos_profiles_deleteOk(self):
        util = cigar.CigarUtil(0, "1M")
        actual = util._pos_profiles("SDM")
        self.assertEquals((2, ["S", "D", "M"]), actual)

    def test_pos_profiles_trailingHardclips(self):
        util = cigar.CigarUtil(0, "1M")
        actual = util._pos_profiles("SDMHH")
        self.assertEquals((2, ["S", "D", "M", "HH"]), actual)


    def test_partition_cigar(self):
        util = cigar.CigarUtil(42, "2S" "6M" "2S")
        (c1, c2, c3) = util._partition_cigar(42, 48)
        self.assertEquals(cigar.CigarUtil(40,"2S"), c1)
        self.assertEquals(cigar.CigarUtil(42,"6M"), c2)
        self.assertEquals(cigar.CigarUtil(48,"2S"), c3)

    def test_partition_cigar_refConsumingOps(self):
        util = cigar.CigarUtil(42, "5M5X")
        # | 4444444455
        # | 2345678901
        # | MMMMM
        # |      XXXXX
        (c1, c2, c3) = util._partition_cigar(42, 47)
        self.assertEquals(cigar.CigarUtil(42,""), c1)
        self.assertEquals(cigar.CigarUtil(42,"5M"), c2)
        self.assertEquals(cigar.CigarUtil(47,"5X"), c3)
        (c1, c2, c3) = util._partition_cigar(47, 52)
        self.assertEquals(cigar.CigarUtil(42,"5M"), c1)
        self.assertEquals(cigar.CigarUtil(47,"5X"), c2)
        self.assertEquals(cigar.CigarUtil(52,""), c3)

    def test_partition_cigar_refNonconsumingOps(self):
        util = cigar.CigarUtil(42, "5M4I1M")
        (c1, c2, c3) = util._partition_cigar(42, 47)
        self.assertEquals(cigar.CigarUtil(42,""), c1)
        self.assertEquals(cigar.CigarUtil(42,"5M"), c2)
        self.assertEquals(cigar.CigarUtil(47,"4I1M"), c3)
        (c1, c2, c3) = util._partition_cigar(47, 48)
        self.assertEquals(cigar.CigarUtil(42,"5M"), c1)
        self.assertEquals(cigar.CigarUtil(47,"4I1M"), c2)
        self.assertEquals(cigar.CigarUtil(48,""), c3)

    def test_partition_cigar_refNonconsumingOpsWithFlankingSoftclips(self):
        util = cigar.CigarUtil(42, "5S" "5M" "4I" "1M" "5S")
        (c1, c2, c3) = util._partition_cigar(42, 47)
        self.assertEquals(cigar.CigarUtil(37,"5S"), c1)
        self.assertEquals(cigar.CigarUtil(42,"5M"), c2)
        self.assertEquals(cigar.CigarUtil(47,"4I1M5S"), c3)
        (c1, c2, c3) = util._partition_cigar(47, 48)
        self.assertEquals(cigar.CigarUtil(37,"5S5M"), c1)
        self.assertEquals(cigar.CigarUtil(47,"4I1M"), c2)
        self.assertEquals(cigar.CigarUtil(48,"5S"), c3)
        (c1, c2, c3) = util._partition_cigar(48, 53)
        self.assertEquals(cigar.CigarUtil(37,"5S5M4I1M"), c1)
        self.assertEquals(cigar.CigarUtil(48,"5S"), c2)
        self.assertEquals(cigar.CigarUtil(53,""), c3)

    def test_partition_outOfBounds(self):
        # ref   | 4444444444
        #       | 0123456789
        # read  | SSMMMMMXX
        # index | 210123456
        util = cigar.CigarUtil(42, "2S" "6M" "2X")
        (c1, c2, c3) = util._partition_cigar(30, 42)
        self.assertEquals(cigar.CigarUtil(40,""), c1)
        self.assertEquals(cigar.CigarUtil(40,"2S"), c2)
        self.assertEquals(cigar.CigarUtil(42,"6M2X"), c3)
        (c1, c2, c3) = util._partition_cigar(48, 60)
        self.assertEquals(cigar.CigarUtil(40,"2S6M"), c1)
        self.assertEquals(cigar.CigarUtil(48,"2X"), c2)
        self.assertEquals(cigar.CigarUtil(50,""), c3)
        (c1, c2, c3) = util._partition_cigar(20, 30)
        self.assertEquals(cigar.CigarUtil(40,""), c1)
        self.assertEquals(cigar.CigarUtil(40,""), c2)
        self.assertEquals(cigar.CigarUtil(40,"2S6M2X"), c3)
        (c1, c2, c3) = util._partition_cigar(100, 110)
        self.assertEquals(cigar.CigarUtil(40,"2S6M2X"), c1)
        self.assertEquals(cigar.CigarUtil(50,""), c2)
        self.assertEquals(cigar.CigarUtil(50,""), c3)

    def test_partition_cigar_targetRegionBoundedByShortReverseRead(self):
        util = cigar.CigarUtil(100, "20M")
        (c1, c2, c3) = util._partition_cigar(90, 120)
        self.assertEquals(100, c1.reference_start)
        self.assertEquals("", c1.cigar)
        self.assertEquals(100, c2.reference_start)
        self.assertEquals("20M", c2.cigar)
        self.assertEquals(120, c3.reference_start)
        self.assertEquals("", c3.cigar)



    def test_softclip_target(self):
        util = cigar.CigarUtil(42, "10M")
        new_util = util.softclip_target(44,50)
        self.assertEquals("2S6M2S", new_util.cigar)
        self.assertEquals(44, new_util.reference_start)

    def test_softclip_target_validatesLength(self):
        base = cigar.CigarUtil(42, "10M")
        base.query_length = 100
        self.assertRaises(katana.util.KatanaException,
                          base.softclip_target,
                          44,
                          50)


    def test_softclip_target_flankingSoftclips(self):
        util = cigar.CigarUtil(42, "2S" "6M" "2S")
        # 4444444444
        # 0123456789
        # SSMMMMMMSS
        # SSSMMMMSSS
        new_util = util.softclip_target(43,47)
        self.assertEquals("3S4M3S", new_util.cigar)
        self.assertEquals(43, new_util.reference_start)

    def test_softclip_target_flankingHardclips(self):
        util = cigar.CigarUtil(42, "2H1S" "4M" "1S2H")
        #3444444444
        #9012345678
        #HHSMMMMSHH
        #HHSSMMSSHH
        new_util = util.softclip_target(43,45)
        self.assertEquals("2H2S" "2M" "2S2H", new_util.cigar)
        self.assertEquals(43, new_util.reference_start)

    def test_softclip_target_deletes(self):
        util = cigar.CigarUtil(42, "2M3D" "4M" "1S")
        #4444444455
        #2345678901
        #MMDDD
        #     MMMM
        #         S
        #SS   MMMMS
        new_util = util.softclip_target(47,51)
        self.assertEquals("2S" "4M" "1S", new_util.cigar)
        self.assertEquals(47, new_util.reference_start)

    def test_softclip_target_edgeInsert(self):
        util = cigar.CigarUtil(42, "3M" "1I4M" "2X")
        #444 444445
        #234 567890
        #ATAAACGTAC
        #MMMI
        #    MMMM
        #        XX
        #SSSSMMMMSS
        new_util = util.softclip_target(45,49)
        self.assertEquals("4S" "4M" "2S", new_util.cigar)
        self.assertEquals(45, new_util.reference_start)

    def test_softclip_target_edgeDelete(self):
        util = cigar.CigarUtil(42, "3M" "1D5M" "2S")
        #44444444555
        #23456789012
        #ATA ACGTACG
        #MMM
        #   DMMMM
        #         MSS
        #SSS MMMMSSS
        new_util = util.softclip_target(45,50)
        self.assertEquals("3S" "4M" "3S", new_util.cigar)
        self.assertEquals(46, new_util.reference_start)


class NullCigarTestCase(unittest.TestCase):
    def test_null_cigar(self):
        c = cigar.NullCigarUtil(reference_start=42)
        self.assertEquals(42, c.reference_start)
        self.assertEquals("*", c.cigar)
        self.assertEquals(0, c.query_length)
        self.assertEquals(True, c.is_valid)
        self.assertEquals("*", c.softclip_target(0,100).cigar)
        self.assertEquals(42, c.softclip_target(0,100).reference_start)


class CigartestCase(unittest.TestCase):
    def test_cigar_factory(self):
        mock_read = MicroMock(reference_start=42, cigarstring="10M")
        c = cigar.cigar_factory(mock_read)
        self.assertIsInstance(c, cigar.CigarUtil)
        self.assertEquals(42, c.reference_start)
        self.assertEquals("10M", c.cigar)

    def test_cigar_factory_nullCigar(self):
        mock_read = MicroMock(reference_start=42, cigarstring="*")
        c = cigar.cigar_factory(mock_read)
        self.assertIsInstance(c, cigar.NullCigarUtil)
        self.assertEquals(42, c.reference_start)
        self.assertEquals("*", c.cigar)

        mock_read = MicroMock(reference_start=42, cigarstring="")
        c = cigar.cigar_factory(mock_read)
        self.assertIsInstance(c, cigar.NullCigarUtil)
        self.assertEquals(42, c.reference_start)
        self.assertEquals("*", c.cigar)
