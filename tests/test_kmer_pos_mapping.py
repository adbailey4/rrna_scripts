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

        cls.cpm = KmerPosMapping(cls.test_reference, cls.pos_file, cls.mod_file)

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
        self.assertTrue(self.cpm.get_covered_bases())
        csp = self.cpm.contig_strand_position(contig="RDN18-1", strand="+", position=1186)
        self.assertSetEqual(self.cpm.pos_2_kmers[csp], {'TCAGl', 'gCAGT', 'TCAGT', 'gCAGl'})
        self.assertSequenceEqual(self.cpm.kmer_2_pos["gCAGl"], [csp])
        self.assertSequenceEqual(self.cpm.pos_2_covered_kmers[csp],
                                 [{'lTTAA', 'TTTAA'},
                                  {'GlTTA', 'GTTTA'},
                                  {'AGTTT', 'AGlTT'},
                                  {'CAGlT', 'CAGTT'},
                                  {'TCAGl', 'gCAGT', 'TCAGT', 'gCAGl'}])
        self.assertSequenceEqual(self.cpm.pos_2_overlap_pos[csp],
                                 [self.cpm.contig_strand_position(contig="RDN18-1", strand="+", position=1190)])

    def test_read_in_mod_data(self):
        mod_data = KmerPosMapping.read_in_mod_data(self.test_mod_file)
        self.assertSequenceEqual(list(mod_data.columns),
                                 ["contig", "mod", "pos", "percent", "strand", "reference_index"])

    # TODO
    def test_index_kmers_in_reference(self):
        pass

    def test_get_positions_covering_kmer(self):
        positions = self.cpm.get_positions_covering_kmer("gCAGl")
        self.assertSequenceEqual([self.cpm.contig_strand_position(contig="RDN18-1", strand="+", position=1186)],
                                 positions)

    def test_get_kmers_covering_mod(self):
        kmers = self.cpm.get_kmers_covering_mod(
            contig="RDN18-1", strand="+", position=1186)
        self.assertSequenceEqual([{'lTTAA', 'TTTAA'},
                                  {'GlTTA', 'GTTTA'},
                                  {'AGTTT', 'AGlTT'},
                                  {'CAGlT', 'CAGTT'},
                                  {'TCAGl', 'gCAGT', 'TCAGT', 'gCAGl'}],
                                 kmers)

    def test_get_reference_position_kmers(self):
        kmers = self.cpm.get_reference_position_kmers(contig="RDN18-1", strand="+", position=1186)
        self.assertSetEqual({'TCAGl', 'gCAGT', 'TCAGT', 'gCAGl'}, kmers)
        kmers = self.cpm.get_reference_position_kmers(contig="RDN18-1", strand="+", position=1800)
        self.assertFalse(kmers)
        kmers = self.cpm.get_reference_position_kmers(contig="RDN18-1", strand="+", position=4)
        self.assertSetEqual({"TTGGT"}, kmers)

    def test_get_neighboring_mods(self):
        mods = self.cpm.get_neighboring_mods(contig="RDN18-1", strand="+", position=1186)
        self.assertSequenceEqual(mods,
                                 [self.cpm.contig_strand_position(contig="RDN18-1", strand="+", position=1190)])
        mods = self.cpm.get_neighboring_mods(contig="RDN18-1", strand="+", position=1)
        self.assertFalse(mods)


if __name__ == "__main__":
    unittest.main()
