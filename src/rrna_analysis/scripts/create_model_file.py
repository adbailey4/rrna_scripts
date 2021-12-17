#!/usr/bin/env python
"""Create model file from reference sequence and positions"""
########################################################################
# File: create_model_file.py
#  executable: create_model_file.py
#
# Author: Andrew Bailey
# History: 12/16/21 Created
########################################################################

import os
from argparse import ArgumentParser
from pathlib import Path
from collections import defaultdict
from itertools import zip_longest

from py3helpers.seq_tools import ReverseComplement, ReferenceHandler
from py3helpers.utils import time_it, merge_lists

from signalalign.hiddenMarkovModel import HmmModel
from signalalign.utils.sequenceTools import CustomAmbiguityPositions, reverse_complement

from rrna_analysis.kmer_pos_mapping import get_kmers_from_seq, get_kmer_permutations

rc = ReverseComplement()


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--reference', action='store',
                        dest='reference', required=True, type=str, default=None,
                        help="Path to reference")
    parser.add_argument('--original_model', action='store',
                        dest='original_model', required=True, type=str, default=None,
                        help="Path to original_model")
    parser.add_argument('--positions', action='store',
                        dest='positions', required=True, type=str, default=None,
                        help="Path to positions of modifications")
    parser.add_argument('--output_dir', action='store',
                        dest='output_dir', required=True, type=str, default=None,
                        help="Path to output directory")
    parser.add_argument('--kmer_length', action='store',
                        dest='kmer_length', required=False, type=int, default=5,
                        help="kmer length")
    parser.add_argument('--debug', action='store_true',
                        dest='debug', required=False,
                        help="Print debug messages")

    args = parser.parse_args()
    return args


def get_mod_positions(alt_c_positions, k=5):
    mods_data = CustomAmbiguityPositions.parseAmbiguityFile(alt_c_positions)
    mod_positions = {}
    grouped_positions = {}
    for i, data in mods_data.groupby(["contig", "strand"]):
        contig = i[0]
        strand = i[1]
        #     print(data)
        mod_positions[contig + strand] = data["position"].values
        values = sorted(data["position"].values)
        recent = values[0]
        groups = []
        group = [recent]
        for x in values[1:]:
            if x - recent <= k:
                group.append(x)
            else:
                groups.append(group)
                group = [x]
            recent = x
        grouped_positions[contig + strand] = groups
    return mod_positions, grouped_positions


def create_rrna_model(original_rna_model, output_dir, noise=0, new_variants="abc",
                      replacement_bases="AAA"):
    alphabet = "ACGT"
    if noise > 0:
        tmp_rna_model = os.path.join(output_dir, "rna_r94_5mer_{}_noise_" + str(noise) + ".model")
    else:
        tmp_rna_model = os.path.join(output_dir, "rna_r94_5mer_{}.model")

    current_model_path = original_rna_model
    new_variants = sorted(new_variants)
    for i, (variant, rep_base) in enumerate(zip(new_variants, replacement_bases)):
        alphabet += variant
        rna_model = HmmModel(current_model_path, rna=True)
        # print(alphabet)
        rna_model.write_new_model(tmp_rna_model.format(alphabet), alphabet, rep_base, noise=noise)
        if i > 0:
            os.remove(current_model_path)
        current_model_path = tmp_rna_model.format(alphabet)
    return current_model_path


def get_list_full_sequence(rh, positions_data, edit=True):
    full_sequence = {}
    for i, data in positions_data.groupby(["contig", "strand"]):
        contig = i[0]
        strand = i[1]
        sequence = rh.get_sequence(contig, None, None)
        #     sequence = self.ref_handler.get_sequence(contig, None, None)
        if strand == "-":
            sequence = rc.complement(sequence)
        list_seq = list(sequence)
        if edit:
            for i, row in data.iterrows():
                assert list_seq[row.position] == row.change_from, f"ref base:{list_seq[row.position]} != position " \
                                                                  f"change_from:{row.change_from} "
                list_seq[row.position] = row.change_to
        full_sequence[contig + strand] = list_seq
    return full_sequence


def get_kmers_from_list(sequence, index, k, rna):
    sub_seq = sequence[index - k + 1:index + k]
    if rna:
        sub_seq = sub_seq[::-1]
    kmers = get_kmer_permutations(sub_seq)
    all_kmers = []
    for x in kmers:
        all_kmers.append(get_kmers_from_seq(x, k))
    all_kmers = [set(x) for x in zip(*all_kmers)]
    if rna:
        all_kmers = all_kmers[::-1]
    return all_kmers


def get_kmers_from_full_sequence(full_sequence, mod_positions, k=5, rna=True):
    contig_strand_position_kmers = {}
    for contig_strand, ref_indices in mod_positions.items():
        position_kmers = {}
        contig_sequence = full_sequence[contig_strand]
        for index in ref_indices:
            all_kmers = get_kmers_from_list(contig_sequence, index, k, rna)
            position_kmers[str(index)] = all_kmers
        contig_strand_position_kmers[contig_strand] = position_kmers
    return contig_strand_position_kmers


def get_expected_total(group, k=5):
    total = 0
    curr_window = []
    prev_min = group[0]
    for index in range(group[0] - k + 1, group[-1] + 1):
        if index + (k - 1) in group:
            curr_window.append(index + (k - 1))
        if len(curr_window) > 0 and index > curr_window[0]:
            curr_window.pop(0)
        kmers = 2 ** len(curr_window)
        total += kmers
    #         print(index, kmers, curr_window)
    return total


def get_contig_strand_group_kmers(contig_strand_position_kmers, grouped_positions, k=5):
    contig_strand_group_kmers = {}
    seen_kmers = {}
    for contig_strand, position_kmers in contig_strand_position_kmers.items():
        gp = grouped_positions[contig_strand]
        group_kmers_list = []
        for group in gp:
            #         print("NEW GROUP", group)
            group_kmers = []
            #         total = 0
            for index in group:
                #             list_of_pos = list(itertools.chain.from_iterable([list(x) for x in position_kmers[str(index)]]))
                #             total += len(list_of_pos)
                group_kmers.extend(set.union(*position_kmers[str(index)]))
            group_kmers = set(group_kmers)
            expected_total = get_expected_total(group, k=k)
            #         print(len(group), group)
            #         print(expected_total, len(group_kmers), len(set(group_kmers)))
            #             print(group_kmers)
            assert expected_total == len(
                group_kmers), f"expected_total != len(group_kmers), {expected_total} != {len(group_kmers)}"
            group_kmers_list.append(group_kmers)
        contig_strand_group_kmers[contig_strand] = group_kmers_list
    return contig_strand_group_kmers


def get_fix_positions(contig_strand_group_kmers):
    all_kmers_set = set()
    for cstrand in contig_strand_group_kmers.keys():
        all_kmers_set |= set.union(*contig_strand_group_kmers[cstrand])
    all_kmers_dict = {key: [] for key in all_kmers_set}
    fix_positions = defaultdict(list)
    for cstrand in contig_strand_group_kmers.keys():
        for i, group in enumerate(contig_strand_group_kmers[cstrand]):
            for kmer in group:
                all_kmers_dict[kmer].append(cstrand + str(i))
                if len(all_kmers_dict[kmer]) > 1:
                    fix_positions[cstrand + str(i)] += [kmer]
    return fix_positions


def get_max_pos_max_len(fix_positions, skip=0):
    max_pos = -1
    max_len = -1
    prev_bests = []
    for pos, k_list in fix_positions.items():
        k_list_len = len(k_list)
        if k_list_len > max_len:
            max_len = k_list_len
            max_pos = pos
            prev_bests.insert(0, [max_pos, max_len])
    if skip >= len(prev_bests):
        max_pos = prev_bests[-1][0]
        max_len = prev_bests[-1][1]
        over = True
    else:
        max_pos = prev_bests[skip][0]
        max_len = prev_bests[skip][1]
        over = False

    return max_pos, max_len, over

def get_kmer_mapping(og_sequence, contig_strand_position_kmers, k=5, rna=True):
    kmer_mapping = []
    for contig_strand, pos_kmer_map in contig_strand_position_kmers.items():
        for pos, kmer_list in pos_kmer_map.items():
            #         print(contig_strand, pos, kmer_list)
            index = int(pos)
            curr_kmers = get_kmers_from_list(og_sequence[contig_strand], index, k, rna)
            #         print(curr_kmers)
            for kmer_pair in zip(curr_kmers, kmer_list):
                k1 = list(kmer_pair[0])[0]
                k2 = list(kmer_pair[1])
                #             print(k1, k2)
                kmer_mapping.append([k1, k2])
    return kmer_mapping


def get_change_base_index(max_pos, fix_positions, grouped_positions, contig_strand_position_kmers, rna=True, k=5,
                          skip=0):
    contig_strand = max_pos[:8]
    pos = int(max_pos[8:])
    kmers = set(fix_positions[max_pos])
    #     print(grouped_positions)
    #     print(contig_strand)
    #     print(grouped_positions[contig_strand])
    #     print(max_pos, pos, len(grouped_positions[contig_strand]))
    target_positions = grouped_positions[contig_strand][pos]
    best_pos = None
    best_pos_index = None
    best_kmers = None
    curr_max = 0
    for x in target_positions:
        kmer_sets = contig_strand_position_kmers[contig_strand][str(x)]
        for i, kmer_set in enumerate(kmer_sets):
            overlap_kmers = kmers & kmer_set
            n_kmers = len(overlap_kmers)
            if n_kmers > curr_max:
                curr_max = n_kmers
                best_pos = x
                best_pos_index = i
                best_kmers = overlap_kmers
    #     print(best_pos, best_pos_index, best_kmers)
    sub_pos = None
    canonical = {"A", "T", "G", "C"}
    sub_pos_list = []
    for i, x in enumerate(zip_longest(*best_kmers)):
        bases = set(x)
        canonical & bases
        n_bases = len(set(x))
        if n_bases == 1 and len(canonical & bases) == 1:
            sub_pos = i
            sub_pos_list.insert(0, sub_pos)

    if skip >= len(sub_pos_list):
        sub_pos = sub_pos_list[-1]
        over = True
    else:
        sub_pos = sub_pos_list[skip]
        over = False
    if rna:
        sub_index = k - sub_pos - 1
    change_base_index = (best_pos - k + 1) + best_pos_index + sub_index
    return change_base_index, over


def create_efficient_positions_file(reference, ambig_positions, output_path, k=5, rna=True, debug=False):
    rh = ReferenceHandler(reference)
    mod_positions, grouped_positions = get_mod_positions(ambig_positions, k=k)
    positions_data = CustomAmbiguityPositions.parseAmbiguityFile(ambig_positions)
    positions_data.loc[:, "change_to"] = "ab"
    # positions_data = CustomAmbiguityPositions.parseAmbiguityFile(path)
    full_sequence = get_list_full_sequence(rh, positions_data)

    fix_positions = [1, 1]
    iterration = 0
    kmer_mapping = []

    pos_skip = 0
    base_skip = 0
    while len(fix_positions) > 0 and iterration < 1000:
        if debug:
            print("ANOTHER ITER", len(fix_positions))
        contig_strand_position_kmers = get_kmers_from_full_sequence(full_sequence, mod_positions, k=k, rna=rna)
        contig_strand_group_kmers = get_contig_strand_group_kmers(contig_strand_position_kmers,
                                                                  grouped_positions=grouped_positions, k=k)
        fix_positions = get_fix_positions(contig_strand_group_kmers)
        if len(fix_positions) == 0:
            break
        max_pos, max_len, pos_over = get_max_pos_max_len(fix_positions, skip=pos_skip)
        change_base_index, base_over = get_change_base_index(max_pos, fix_positions, grouped_positions,
                                                             contig_strand_position_kmers, k=k, rna=rna, skip=base_skip)

        contig_strand = max_pos[:8]
        contig = max_pos[:7]
        position = change_base_index
        strand = max_pos[7]
        change_from = full_sequence[contig_strand][change_base_index]
        change_to = "c"

        try:
            full_sequence[contig_strand][change_base_index] = change_to
            contig_strand_position_kmers = get_kmers_from_full_sequence(full_sequence, mod_positions, k=k, rna=rna)
            contig_strand_group_kmers = get_contig_strand_group_kmers(contig_strand_position_kmers,
                                                                      grouped_positions=grouped_positions, k=k)
            #  will not run if causes error

            new_row = {"contig": max_pos[:7],
                       "position": change_base_index,
                       "strand": max_pos[7],
                       "change_from": change_from,
                       "change_to": change_to}
            positions_data = positions_data.append(new_row, ignore_index=True)
            if debug:
                print("PASS", max_pos, change_base_index, pos_skip, base_skip)
            fix_positions = get_fix_positions(contig_strand_group_kmers)
            if len(fix_positions) == 0:
                print("DONE")

        except AssertionError:
            full_sequence[contig_strand][change_base_index] = change_from
            if pos_over:
                base_skip += 1
            if not pos_over:
                pos_skip += 1
            if pos_over and base_over:
                raise AssertionError("ERROR", contig_strand, change_base_index, pos_skip, base_skip)
            if debug:
                print("SKIPPING", contig_strand, change_base_index, pos_skip, base_skip)

        iterration += 1
    positions_data.to_csv(output_path, sep="\t", header=False, index=False)
    #     get kmer mapping
    og_sequence = get_list_full_sequence(rh, positions_data, edit=False)
    kmer_mapping = get_kmer_mapping(og_sequence, contig_strand_position_kmers, k=k, rna=rna)

    return kmer_mapping


def edit_model_with_kmer_mapping(model, kmer_mapping, rna=True):
    rna_model = HmmModel(model, rna=rna)
    for k, change in kmer_mapping:
        normal_mean, normal_sd = rna_model.get_event_mean_gaussian_parameters(k)
        for k2 in change:
            rna_model.set_kmer_event_mean_params(k2, normal_mean, normal_sd)
    rna_model.normalized = True
    rna_model.write(model)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = os.path.join(output_dir, "small_model_" + os.path.basename(args.positions))
    print(f"Removing matching kmers... ")

    kmer_mapping = create_efficient_positions_file(args.reference, args.positions, output_path,
                                                   k=args.kmer_length, rna=True, debug=args.debug)
    new_chars = "".join(set(''.join(map(str, (merge_lists([x[1] for x in kmer_mapping]))))) - {'A', 'C', 'G', 'T'})
    final_model_path = create_rrna_model(args.original_model, output_dir,
                                         noise=0, new_variants=new_chars, replacement_bases="A"*len(new_chars))
    edit_model_with_kmer_mapping(final_model_path, kmer_mapping, rna=True)
    print(f"Final Model: {final_model_path}")
    print(f"Final Positions: {output_path}")


if __name__ == '__main__':
    ret, time = (time_it(main))
    print(f"Running Time: {time} s")
