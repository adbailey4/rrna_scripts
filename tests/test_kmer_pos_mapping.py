#!/usr/bin/env python
"""Tests for kmer_pos_mapping.py"""
########################################################################
# File: kmer_pos_mapping.py
#  executable: kmer_pos_mapping.py
#
# Author: Andrew Bailey
# History: Created 12/03/20
########################################################################
import unittest

from rrna_analysis.kmer_pos_mapping import *


class TestKmerPosMapping(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestKmerPosMapping, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])
        cls.test_reference = os.path.join(cls.HOME, "tests/test_files/yeast_25S_18S.fa")
        cls.test_positions = os.path.join(cls.HOME, "tests/test_files/test.positions")
        cls.test_mod_file = os.path.join(cls.HOME, "tests/test_files/test_mod_file.csv")
        cls.mod_file = os.path.join(cls.HOME, "tests/test_files/mod_file.csv")
        cls.pos_file = os.path.join(cls.HOME, "tests/test_files/yeast_18S_25S_variants.positions")
        cls.kpm = KmerPosMapping(cls.test_reference, cls.pos_file, cls.mod_file)

        cls.alt_c_pos = os.path.join(cls.HOME, "tests/test_files/alt_c_yeast_18S_25S_variants.positions")
        cls.alt_c_kpm = KmerPosMapping(cls.test_reference, cls.alt_c_pos, cls.mod_file)

    def test_create_mod_handler(self):
        self.assertSequenceEqual(list(self.kpm.create_mod_handler().columns),
                                 ['contig', 'position', 'strand', 'change_from', 'change_to', 'mod',
                                  'pos', 'percent', 'reference_index', 'delta1_below', 'delta1_above', 'delta2_below',
                                  'delta2_above', 'delta3_below', 'delta3_above', 'delta4_below', 'delta4_above',
                                  'delta', 'in_2prime', 'in_pseudo', 'pseudo_in_other',
                                  '2prime_in_other', 'in_unknown'])

    def test_get_kmers_from_seq(self):
        kmers = get_kmers_from_seq("ATGCA", kmer_length=5)
        self.assertSequenceEqual(["ATGCA"], kmers)
        self.assertRaises(AssertionError, get_kmers_from_seq, "AT")
        kmers = get_kmers_from_seq("AT*0A", kmer_length=2)
        self.assertSequenceEqual(["AT", "T*", "*0", "0A"], kmers)

    def test_get_kmer_permutations(self):
        kmers = get_kmer_permutations(["AT", "GC"])
        self.assertSequenceEqual(["AG", "AC", "TG", "TC"], kmers)
        self.assertRaises(AssertionError, get_kmer_permutations, ["AA"])

    def test_get_covered_kmers(self):
        rh = ReferenceHandler(self.test_reference)
        kmers = get_covered_kmers(rh, "RDN18-1", "+", [4], ["TE"], kmer_length=2)
        self.assertSequenceEqual([{"CT", "CE"}, {"TG", "EG"}], kmers)
        kmers = get_covered_kmers(rh, "RDN18-1", "+", [27], ["Aa"], kmer_length=5, rna=True)
        self.assertSequenceEqual([{'ATACT', 'aTACT'},
                                  {'TATAC', 'TaTAC'},
                                  {'GTaTA', 'GTATA'},
                                  {'CGTAT', 'CGTaT'},
                                  {'TCGTa', 'TCGTA'}], kmers)
        self.assertRaises(AssertionError, get_covered_kmers, rh, "RDN18-1", "+", [4], ["AE"])
        kmers = get_covered_kmers(rh, "RDN18-1", "+", [4, 5], ["TE", "GF"], kmer_length=2)
        self.assertSequenceEqual([{'CT', 'CE'}, {'EG', 'TG', 'EF', 'TF'}, {'GG', 'FG'}], kmers)

    def test_get_covered_bases(self):
        csp = self.kpm.contig_strand_position(contig="RDN18-1", strand="+", position=1186)
        self.assertSetEqual(self.kpm.pos_2_kmers[csp], {'TCAGl', 'gCAGT', 'TCAGT', 'gCAGl'})
        self.assertSequenceEqual(self.kpm.kmer_2_pos["gCAGl"], [csp])
        self.assertSequenceEqual(self.kpm.pos_2_covered_kmers[csp],
                                 [{'lTTAA', 'TTTAA'},
                                  {'GlTTA', 'GTTTA'},
                                  {'AGTTT', 'AGlTT'},
                                  {'CAGlT', 'CAGTT'},
                                  {'TCAGl', 'gCAGT', 'TCAGT', 'gCAGl'}])
        self.assertSequenceEqual(self.kpm.pos_2_overlap_pos[csp],
                                 [self.kpm.contig_strand_position(contig="RDN18-1", strand="+", position=1190)])
        # csp = self.kpm.contig_strand_position(contig="RDN25-1", strand="+", position=648)
        # print(self.kpm.pos_2_covered_kmers[csp])

    def test_alt_c_get_covered_bases(self):
        csp = self.alt_c_kpm.contig_strand_position(contig="RDN18-1", strand="+", position=1186)
        self.assertSetEqual(self.alt_c_kpm.pos_2_kmers[csp], {"TCAGT", 'qCAGl', 'gCAGq', 'qCAGq', 'gCAGl'})
        self.assertSequenceEqual(self.alt_c_kpm.kmer_2_pos["gCAGl"], [csp])
        self.assertSequenceEqual(self.alt_c_kpm.pos_2_covered_kmers[csp],
                                 [{'lTTAA', 'qTTAA'},
                                  {'GlTTA', 'GqTTA'},
                                  {'AGlTT', 'AGqTT'},
                                  {'CAGlT', 'CAGqT'},
                                  {'qCAGl', 'gCAGq', 'qCAGq', 'gCAGl'}])
        self.assertSequenceEqual(self.alt_c_kpm.pos_2_overlap_pos[csp],
                                 [self.alt_c_kpm.contig_strand_position(contig="RDN18-1", strand="+", position=1190)])
        # csp = self.kpm.contig_strand_position(contig="RDN25-1", strand="+", position=648)
        # print(self.kpm.pos_2_covered_kmers[csp])


    def test_read_in_mod_data(self):
        mod_data = KmerPosMapping.read_in_mod_data(self.test_mod_file)
        self.assertSequenceEqual(list(mod_data.columns),
                                 ["contig", "mod", "pos", "percent", "strand", "reference_index"])

    def test_get_positions_covering_kmer(self):
        positions = self.kpm.get_positions_covering_kmer("gCAGl")
        self.assertSequenceEqual([self.kpm.contig_strand_position(contig="RDN18-1", strand="+", position=1186)],
                                 positions)

    def test_get_kmers_covering_mod(self):
        kmers = self.kpm.get_kmers_covering_mod(
            contig="RDN18-1", strand="+", position=1186)
        self.assertSequenceEqual([{'lTTAA', 'TTTAA'},
                                  {'GlTTA', 'GTTTA'},
                                  {'AGTTT', 'AGlTT'},
                                  {'CAGlT', 'CAGTT'},
                                  {'TCAGl', 'gCAGT', 'TCAGT', 'gCAGl'}],
                                 kmers)

    def test_get_reference_position_kmers(self):
        kmers = self.kpm.get_reference_position_kmers(contig="RDN18-1", strand="+", position=1186)
        self.assertSetEqual({'TCAGl', 'gCAGT', 'TCAGT', 'gCAGl'}, kmers)
        kmers = self.kpm.get_reference_position_kmers(contig="RDN18-1", strand="+", position=1800)
        self.assertFalse(kmers)
        kmers = self.kpm.get_reference_position_kmers(contig="RDN18-1", strand="+", position=4)
        self.assertSetEqual({"TTGGT"}, kmers)

    def test_get_neighboring_mods(self):
        mods = self.kpm.get_neighboring_mods(contig="RDN18-1", strand="+", position=1186)
        self.assertSequenceEqual(mods,
                                 [self.kpm.contig_strand_position(contig="RDN18-1", strand="+", position=1190)])
        mods = self.kpm.get_neighboring_mods(contig="RDN18-1", strand="+", position=1)
        self.assertFalse(mods)


if __name__ == "__main__":
    unittest.main()
