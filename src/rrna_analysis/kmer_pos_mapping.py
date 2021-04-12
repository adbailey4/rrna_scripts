#!/usr/bin/env python
"""class for easy kmer to position to modification mapping access"""
########################################################################
# File: kmer_pos_mapping.py
#  executable: kmer_pos_mapping.py
#
# Author: Andrew Bailey
# History: Created 12/03/20
########################################################################

import os
from collections import namedtuple, defaultdict

import numpy as np
import pandas as pd
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement
from signalalign.utils.sequenceTools import CustomAmbiguityPositions, reverse_complement


def get_kmers_from_seq(sequence, kmer_length=5):
    """Get list of kmers given a sequence and kmer length
    :param sequence: string
    :param kmer_length: kmer length (default 5)
    """
    assert kmer_length <= len(sequence), \
        "kmer length cannot be larger than sequence length: {} > {}".format(kmer_length, len(sequence))
    kmers = []
    x = 0
    while x <= len(sequence) - kmer_length:
        kmers.append(sequence[x:kmer_length + x])
        x += 1
    return kmers


def get_kmer_permutations(bases: list):
    """Given a list of strings, generate a list of permutations where each string are all the nucleotide options
    for each position

    ex: input: ["AT", "GC"] output-> ["AG", "AC", "TG", "TC"]
    :param bases: list of bases for each position
    """
    switch = [1]
    multiply = 1
    for x in bases[::-1]:
        assert len(set(x)) == len(x), f"Each string cannot have any duplicate characters. len(set(x)) != len(x): x={x}"
        switch.insert(0, multiply)
        multiply *= len(x)
    num_kmers = multiply
    i = 0
    kmers = []
    kmer = None
    while i < num_kmers:
        kmer = ""
        num = i
        for x, chars in enumerate(bases):
            # print(num, switch[x], num // switch[x], end=" ")
            # print(chars[num // switch[x]])
            kmer += chars[num // switch[x]]
            num -= (num // switch[x]) * switch[x]
        kmers.append(kmer)
        i += 1
    assert len(set(kmers)) == num_kmers, "Made wrong number of kmers: {} != {}".format(len(set(kmer)), num_kmers)
    return kmers


def get_covered_kmers(reference_handler: ReferenceHandler, chromosome_name: str, strand: str, pos: list,
                      variant_chars: list,
                      rna: bool = False,
                      kmer_length: int = 5,
                      canonical_replacement={"A": "n", "T": 'q', "G": "p", "C": "o"}) -> list:
    """Get a list of lists of kmers covering a certain position. The first list will always be canonical,
    any subsequent kmer lists are the same kmers with the specified position replaced with non
    reference characters
    :param canonical_replacement: if there are replacement bases for canonical, observe those instead
    :param reference_handler: ReferenceHandler object
    :param chromosome_name: reference chromosome name
    :param strand: strand
    :param pos: 0 index position of reference
    :param variant_chars: all modified or canonical nucleotides
    :param rna: boolean option (default False)
    :param kmer_length: kmer length (default 5)
    """
    min_pos = min(pos)
    max_pos = max(pos)

    sequence = reference_handler.get_sequence(chromosome_name, max(0, min_pos - (kmer_length - 1)),
                                              max_pos + kmer_length)
    if strand == "-":
        sequence = reverse_complement(sequence)
    sequence = list(sequence)
    for pos, chars in zip(pos, variant_chars):
        c_base = sequence[kmer_length - 1 + (pos - min_pos)]
        c_rep = canonical_replacement[c_base]
        if len(chars) > 1:
            assert c_rep in chars or c_base in chars, \
                f"Reference base is not in variant characters: pos{pos}:{c_base} not in {chars} " \
                f"and pos{pos}:{c_rep} not in {chars}"
        sequence[kmer_length - 1 + (pos - min_pos)] = chars
    if rna:
        sequence = sequence[::-1]
    kmers = get_kmer_permutations(sequence)
    all_kmers = []
    for x in kmers:
        all_kmers.append(get_kmers_from_seq(x, kmer_length))

    all_kmers = [set(x) for x in zip(*all_kmers)]
    if rna:
        all_kmers = all_kmers[::-1]
    return all_kmers


class KmerPosMapping(object):
    contig_strand_position = namedtuple('contig_strand_position', ['contig', 'strand', 'position'])

    def __init__(self, reference, positions, mods_csv, kmer_length=5,
                 canonical_replacement={"A": "n", "T": 'q', "G": "p", "C": "o"}):
        assert os.path.exists(reference), f"reference path does not exist:{reference}"
        assert os.path.exists(positions), f"positions path does not exist:{positions}"
        assert os.path.exists(mods_csv), f"mods_csv path does not exist:{mods_csv}"
        self.reference = reference
        self.positions = positions
        self.mods_csv = mods_csv
        self.canonical_replacement = canonical_replacement
        self.mod_data = self.read_in_mod_data(self.mods_csv)
        self.subset_mod_data = self.mod_data[["contig", "reference_index", "percent"]]
        self.kmer_length = kmer_length
        self.rna = True
        self.strands = ["+"]
        self.pos_2_covered_kmers = defaultdict()
        self.pos_2_overlap_pos = defaultdict()
        self.pos_2_kmers = defaultdict(set)
        self.kmer_2_pos = defaultdict(list)
        self.ref_handler = ReferenceHandler(self.reference)
        self.positions_data = CustomAmbiguityPositions.parseAmbiguityFile(self.positions)
        self._index_kmers_in_reference()
        self._get_covered_bases()
        self.mod_handler = self.create_mod_handler()

    def create_mod_handler(self):
        all_pos_data = pd.merge(self.positions_data, self.mod_data, left_on=["contig", "strand", "position"],
                                right_on=["contig", "strand", "reference_index"]).sort_values(by=['contig', 'position'])

        all_pos_data['delta1_below'] = all_pos_data.reference_index.diff().shift(-1)
        all_pos_data['delta1_below'] = np.abs(all_pos_data.reference_index.diff(periods=1))
        all_pos_data['delta1_above'] = np.abs(all_pos_data.reference_index.diff(periods=-1))

        all_pos_data['delta2_below'] = np.abs(all_pos_data.reference_index.diff(periods=2))
        all_pos_data['delta2_above'] = np.abs(all_pos_data.reference_index.diff(periods=-2))
        all_pos_data['delta3_below'] = np.abs(all_pos_data.reference_index.diff(periods=3))
        all_pos_data['delta3_above'] = np.abs(all_pos_data.reference_index.diff(periods=-3))
        all_pos_data['delta4_below'] = np.abs(all_pos_data.reference_index.diff(periods=4))
        all_pos_data['delta4_above'] = np.abs(all_pos_data.reference_index.diff(periods=-4))

        all_pos_data['delta'] = all_pos_data[["delta1_below",
                                              "delta1_above",
                                              "delta2_below",
                                              "delta2_above",
                                              "delta3_below",
                                              "delta3_above",
                                              "delta4_below",
                                              "delta4_above"]].min(axis=1)

        all_pos_data["in_2prime"] = (((all_pos_data.change_to.shift(1).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                       (all_pos_data.delta1_below <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(-1).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                       (all_pos_data.delta1_above <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(2).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                       (all_pos_data.delta2_below <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(-2).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                       (all_pos_data.delta2_above <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(3).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                       (all_pos_data.delta3_below <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(-3).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                       (all_pos_data.delta3_above <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(4).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                       (all_pos_data.delta4_below <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(-4).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                       (all_pos_data.delta4_above <= self.kmer_length))
                                      ) & (~all_pos_data.change_to.isin(["Aa", "Cb", "Gc", "Td", "Tdm"])))
        all_pos_data["in_pseudo"] = (((all_pos_data.change_to.shift(1).isin(["Tl"]) &
                                       (all_pos_data.delta1_below <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(-1).isin(["Tl"]) &
                                       (all_pos_data.delta1_above <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(2).isin(["Tl"]) &
                                       (all_pos_data.delta2_below <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(-2).isin(["Tl"]) &
                                       (all_pos_data.delta2_above <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(3).isin(["Tl"]) &
                                       (all_pos_data.delta3_below <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(-3).isin(["Tl"]) &
                                       (all_pos_data.delta3_above <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(4).isin(["Tl"]) &
                                       (all_pos_data.delta4_below <= self.kmer_length)) |
                                      (all_pos_data.change_to.shift(-4).isin(["Tl"]) &
                                       (all_pos_data.delta4_above <= self.kmer_length))
                                      ) & (~all_pos_data.change_to.isin(["Tl"])))

        all_pos_data["pseudo_in_other"] = (((~all_pos_data.change_to.shift(1).isin(["Tl"]) &
                                             (all_pos_data.delta1_below <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(-1).isin(["Tl"]) &
                                             (all_pos_data.delta1_above <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(2).isin(["Tl"]) &
                                             (all_pos_data.delta2_below <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(-2).isin(["Tl"]) &
                                             (all_pos_data.delta2_above <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(3).isin(["Tl"]) &
                                             (all_pos_data.delta3_below <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(-3).isin(["Tl"]) &
                                             (all_pos_data.delta3_above <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(4).isin(["Tl"]) &
                                             (all_pos_data.delta4_below <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(-4).isin(["Tl"]) &
                                             (all_pos_data.delta4_above <= self.kmer_length))
                                            ) & (all_pos_data.change_to.isin(["Tl"])))

        all_pos_data["2prime_in_other"] = (((~all_pos_data.change_to.shift(1).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                             (all_pos_data.delta1_below <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(-1).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                             (all_pos_data.delta1_above <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(2).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                             (all_pos_data.delta2_below <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(-2).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                             (all_pos_data.delta2_above <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(3).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                             (all_pos_data.delta3_below <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(-3).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                             (all_pos_data.delta3_above <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(4).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                             (all_pos_data.delta4_below <= self.kmer_length)) |
                                            (~all_pos_data.change_to.shift(-4).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                             (all_pos_data.delta4_above <= self.kmer_length))
                                            ) & (all_pos_data.change_to.isin(["Aa", "Cb", "Gc", "Td", "Tdm"])))

        all_pos_data["in_unknown"] = (all_pos_data["in_pseudo"] | all_pos_data["in_2prime"])
        return all_pos_data

    @staticmethod
    def read_in_mod_data(mods_csv):
        """
        Read in a mod_csv file with the header [contig,mod,pos,percent,strand]
        :param mods_csv: path to mod_csv file
        :return: pandas df handle
        """
        mods_df = pd.read_csv(mods_csv)
        mods_df["reference_index"] = mods_df["pos"] - 1
        return mods_df

    def _get_covered_bases(self):
        """Get all covered kmers by position.
        """
        contig = None
        pos = None
        strand = None
        variant_bases = None
        all_pos = []
        all_variant_bases = []
        for i, row in self.positions_data.iterrows():
            if i > 0:
                all_pos.append(pos)
                all_variant_bases.append(variant_bases)
                if pos + self.kmer_length <= row[1] or strand != row[2] or contig != row[0]:
                    self._update_kmer_pos_maps(contig, strand, all_pos, all_variant_bases)
                    all_pos = []
                    all_variant_bases = []
            contig = row[0]
            pos = row[1]
            strand = row[2]
            # base = row[3]
            variant_bases = row[4]
        all_pos.append(pos)
        all_variant_bases.append(variant_bases)
        self._update_kmer_pos_maps(contig, strand, all_pos, all_variant_bases)

        return True

    def _update_kmer_pos_maps(self, contig, strand, all_pos, all_variant_bases):
        kmers = get_covered_kmers(self.ref_handler, contig, strand, all_pos, all_variant_bases,
                                  self.rna, self.kmer_length, canonical_replacement=self.canonical_replacement)
        for i, p in enumerate(range(min(all_pos) - (self.kmer_length - 1), max(all_pos) + 1)):
            csp = self.contig_strand_position(contig=contig, strand=strand, position=p)
            self.pos_2_kmers[csp] = self.pos_2_kmers[csp] | kmers[i]
            for k in self.pos_2_kmers[csp]:
                self.kmer_2_pos[k].append(csp)

        for i, p in enumerate(all_pos):
            csp = self.contig_strand_position(contig=contig, strand=strand, position=p)
            offset = p - all_pos[0]
            self.pos_2_covered_kmers[csp] = kmers[offset:self.kmer_length + offset]
            missing_pos = all_pos[:i] + all_pos[i + 1:]
            self.pos_2_overlap_pos[csp] = [self.contig_strand_position(contig=contig,
                                                                       strand=strand,
                                                                       position=p) for p in missing_pos]

    def _index_kmers_in_reference(self):
        """Find all kmers in reference sequence including possibly modified kmers and keep track of index"""
        contigs = self.ref_handler.fasta.references
        rc = ReverseComplement()
        for contig in contigs:
            sequence = self.ref_handler.get_sequence(contig, None, None)
            for strand in self.strands:
                if strand == "-":
                    sequence = rc.complement(sequence)
                ref_kmers = get_kmers_from_seq(sequence, kmer_length=self.kmer_length)
                for p, kmer in enumerate(ref_kmers):
                    if self.rna:
                        kmer = kmer[::-1]
                    csp = self.contig_strand_position(contig=contig, strand=strand, position=p)
                    self.kmer_2_pos[kmer].append(csp)
                    self.pos_2_kmers[csp] = self.pos_2_kmers[csp] | {kmer}

    def get_positions_covering_kmer(self, kmer):
        """Given a kmer, get all reference positions covering that kmer"""
        assert len(kmer) == self.kmer_length, \
            f"Kmer does not match class kmer length. len(kmer) == self.kmer_length. len({kmer}) != {self.kmer_length}"
        if kmer in self.kmer_2_pos.keys():
            return self.kmer_2_pos[kmer]
        else:
            return []

    def get_kmers_covering_mod(self, contig, strand, position):
        """Get all kmers covering a reference position even if there are modifications present
        Note: the kmer length is set by the class instantiation
        :param contig: name of contig
        :param strand: strand ("+", "-")
        :param position: reference 0 based index position

        kmer_length = 2
        There is a mod at position 3 showing C or E could be present
        >chr1
        AAACGGAAAA
        ex: input: "chr1", "+", 4
            output-> [{"CG", "EG"}, {"GG"}]
        """
        csp = self.contig_strand_position(contig, strand, position)
        if csp in self.pos_2_covered_kmers:
            return self.pos_2_covered_kmers[csp]
        else:
            return False

    def get_reference_position_kmers(self, contig, strand, position):
        """Get all kmers covering a reference position even if there are modifications present
        Note: the kmer length is set by the class instantiation
        :param contig: name of contig
        :param strand: strand ("+", "-")
        :param position: reference 0 based index position

        kmer_length = 2
        There is a mod at position 3 showing C or E could be present
        >chr1
        AAACGGAAAA
        ex: input: "chr1", "+", 4
            output-> {"GG"}
        """
        csp = self.contig_strand_position(contig, strand, position)
        if csp in self.pos_2_kmers:
            return self.pos_2_kmers[csp]
        else:
            return False

    def get_neighboring_mods(self, contig, strand, position):
        """Return positions that are within a kmer distance of the given modification
        :param contig: name of contig
        :param strand: strand ("+", "-")
        :param position: reference 0 based index position
        """
        csp = self.contig_strand_position(contig, strand, position)
        if csp in self.pos_2_overlap_pos:
            return self.pos_2_overlap_pos[csp]
        else:
            return False
