#!/usr/bin/env python
"""Plot per position kmer distribution wrapper on training runs"""
########################################################################
# File: run_embed_plot_wrapper.py
#  executable: run_embed_plot_wrapper.py
#
# Author: Andrew Bailey
# History: 06/15/20 Created
########################################################################

import os
import sys
from argparse import ArgumentParser
import numpy as np
from plot_per_position_kmer_distributions import plot_per_position_kmer_distributions
from subprocess import check_call

from py3helpers.utils import create_dot_dict, load_json, time_it, list_dir, merge_lists


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to config file")
    args = parser.parse_args()
    return args


def check_config(config):
    assert os.path.exists(config.save_fig_dir), "save_fig_dir {} does not exist".format(config.save_fig_dir)
    assert os.path.exists(config.positions), "positions {} does not exist".format(config.positions)
    assert os.path.exists(config.reference), "reference {} does not exist".format(config.reference)
    assert os.path.exists(config.ambig_model), "ambig_model {} does not exist".format(config.ambig_model)
    assert isinstance(config.min_prob, float), "min_prob parameter needs to be an float"
    assert isinstance(config.num_threads, int), "Num threads parameter needs to be an int"
    assert isinstance(config.debug, bool), "debug parameter needs to be an boolean"
    assert isinstance(config.rna, bool), "rna parameter needs to be an boolean"
    assert isinstance(config.mod_only, bool), "mod_only parameter needs to be an boolean"
    assert isinstance(config.alphabet, str), "alphabet parameter needs to be an string"
    assert isinstance(config.force_pos_hist, bool), "force_pos_hist parameter needs to be an boolean"
    assert isinstance(config.force_event, bool), "force_event parameter needs to be an boolean"

    for data in config.iterations:
        hmm_model = data["hmm_model"]
        hdp_model = data["hdp_model"]
        if hmm_model:
            assert os.path.exists(hmm_model), "{} does not exist".format(hmm_model)
        if hdp_model:
            assert os.path.exists(hdp_model), "{} does not exist".format(hdp_model)
        for directory in data["sa_output_dirs"]:
            assert os.path.exists(directory), "{} does not exist".format(directory)


def run_embed_plot_wrapper(config):
    config = create_dot_dict(config)
    check_config(config)

    duplicate_config = config
    for i, data in enumerate(config.iterations):
        split_by_pos = "embed_main split_by_position --reference {} --threads {} " \
                       "--alphabet {}".format(config.reference, config.num_threads, config.alphabet)
        if config.rna:
            split_by_pos += " --rna"
        kmer_distributions = "embed_main kmer_distributions --reference {} --ambig_model {} " \
                             "--positions_file {} --min_prob {}".format(config.reference, config.ambig_model,
                                                                     config.positions, config.min_prob)
        for alignment_dir in data["sa_output_dirs"]:
            split_by_pos += " --alignment_files {}".format(alignment_dir)
        output_file = os.path.join(os.path.split(data["sa_output_dirs"][0])[0],
                                   "_".join([os.path.split(x)[-1] for x in data["sa_output_dirs"]])+".event")
        split_by_pos += " --output {}".format(output_file)
        # generate .event file
        if not os.path.exists(output_file):
            print(split_by_pos)
            check_call(split_by_pos.split())
        elif config.force_event:
            os.remove(output_file)
            print(split_by_pos)
            check_call(split_by_pos.split())

        # get kmer distributions
        kmer_distributions += " --event_file {}".format(output_file)
        pos_histograms = os.path.join(os.path.split(data["sa_output_dirs"][0])[0], "position_histograms")
        kmer_distributions += " --output {}".format(pos_histograms)
        if not os.path.isdir(pos_histograms):
            os.mkdir(pos_histograms)
            print(pos_histograms)
            check_call(kmer_distributions.split())
        elif config.force_pos_hist:
            print(kmer_distributions)
            check_call(kmer_distributions.split())
        duplicate_config.iterations[i]["kmer_hist_dir"] = pos_histograms
    plot_per_position_kmer_distributions(duplicate_config)


def main():
    args = parse_args()
    assert os.path.exists(args.config), "{} does not exist".format(args.config)
    run_embed_plot_wrapper(load_json(args.config))


if __name__ == '__main__':
    print(time_it(main)[1])
