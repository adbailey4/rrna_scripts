#!/usr/bin/env python
"""Train model, test model, gather accuracy metrics from tests"""
########################################################################
# File: train_test_accuracy_wrapper.py
#  executable: train_test_accuracy_wrapper.py
#
# Author: Andrew Bailey
# History: 07/23/20 Created
########################################################################

import os
import re
import shutil
from argparse import ArgumentParser
from itertools import zip_longest
from subprocess import check_call

from py3helpers.utils import list_dir, load_json, create_dot_dict
from rrna_analysis.scripts.re_run_signalalign import re_run_signalalign
from rrna_analysis.scripts.run_embed_plot_wrapper import run_embed_plot_wrapper
from signalalign.visualization.plot_multiple_variant_accuracy import plot_multiple_variant_accuracy


# * Train
# * trainModels run --config /home/ubuntu/rRNA2/sa_workspace/supervised/less_data_testing_07_03_20/mod_only_nopgal_80_prob/rrna_trainModels_config.json
#
# * Test
# * [Copy model files bash helper](#Copy-model-files-bash-helper)
# * python /home/ubuntu/rRNA2/re_run_signalalign.py --dir /home/ubuntu/rRNA2/sa_workspace/supervised/less_data_testing_07_03_20/mod_only_nopgal_80_prob/training_models --output_dir /home/ubuntu/rRNA2/sa_workspace/supervised/less_data_testing_07_03_20/mod_only_nopgal_80_prob/testing --base_model /home/ubuntu/rRNA2/sa_workspace/supervised/less_data_testing_07_03_20/mod_only_nopgal_80_prob/re_runSignalAlign_testing_config.json --variants BDEFHIJKLMOPQR --rna
#
#                                                                                                                                                                                                                                                                                                                                                                                                                                                          * Plot Distributions
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                 * [Generate iterations helper](#Generate-iterations-helper)
#     * python /home/ubuntu/rRNA2/run_embed_plot_wrapper.py -c /home/ubuntu/rRNA2/sa_workspace/supervised/less_data_testing_07_03_20/mod_only_nopgal_80_prob/run_embed_plot_wrapper2.config.json
#       * Time:
#
# * Plot accuracy over time
#                      * [Run plot_mulitple_variant_accuracy.py many times](#Re run plot_multiple_variant_accuracy)
#     * [Get all variant calls to single directory](Get all variant calls to single directory)
# * Plot via jupyter notebook


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to config file")
    args = parser.parse_args()
    return args


def re_run_plot_variant_accuracy(output_dir, testing_dir, positions_files, names, suffixes, threshold=0.5):
    assert os.path.isdir(testing_dir), "testing_dir does not exist: {}".format(testing_dir)
    assert os.path.isdir(output_dir), "output_dir does not exist: {}".format(output_dir)
    dirs = [os.path.join(testing_dir, x) for x in os.listdir(testing_dir) if not x.endswith("created_models")]
    assert len(names) == len(suffixes) == len(positions_files), \
        "The length of suffixes, names and positions files must be equal"
    for testing_dir in dirs:
        test_output_dir = os.path.join(output_dir, os.path.basename(testing_dir))
        if os.path.exists(test_output_dir):
            print("{} exists, continuing".format(test_output_dir))
            continue
        samples = []
        print(testing_dir)
        for name, suffix, positions_file in zip_longest(names, suffixes, positions_files):
            print(os.path.join(testing_dir, suffix))
            samples.append({"name": name,
                            "variant_csv": os.path.join(testing_dir, suffix),
                            "positions_file": positions_file})
        plot_multiple_variant_accuracy(test_output_dir, samples, threshold)


def copy_variant_call_files(testing_dir, suffixes, output_dir):
    """Copy variant call files from testing into their own directories"""
    dirs = [os.path.join(testing_dir, x) for x in os.listdir(testing_dir) if not x.endswith("created_models")]
    for vc_dir in dirs:
        test_output_dir = os.path.join(output_dir, os.path.basename(vc_dir))
        os.mkdir(test_output_dir)
        for suffix in suffixes:
            shutil.copy2(os.path.join(vc_dir, suffix), test_output_dir)


def get_trailing_number(s):
    m = re.search(r'\d+$', s)
    return int(m.group()) if m else None


def main():
    args = parse_args()
    assert os.path.exists(args.config), "{} does not exist".format(args.config)
    # plot_per_position_kmer_distributions(load_json(args.config))
    config_dict = create_dot_dict(load_json(args.config))
    print(config_dict)
    training_dir = os.path.join(config_dict.base_directory, "training")
    training_model_dir = os.path.join(config_dict.base_directory, "training_models")
    testing_dir = os.path.join(config_dict.base_directory, "testing")
    all_variant_calls = os.path.join(config_dict.base_directory, "all_variant_calls")

    testing_accuracy_dir = os.path.join(config_dict.base_directory, "testing_accuracy")
    training_accuracy_dir = os.path.join(config_dict.base_directory, "training_accuracy")
    testing_accuracy_csvs_dir = os.path.join(config_dict.base_directory, "testing_accuracy_csvs")
    training_accuracy_csvs_dir = os.path.join(config_dict.base_directory, "training_accuracy_csvs")
    training_distributions_dir = os.path.join(config_dict.base_directory, "training_distributions")

    if config_dict.options["train"]:
        assert not os.path.exists(training_dir), "Directory should not exist: {}. Check base_directory".format(training_dir)
        os.mkdir(training_dir)
        assert not os.path.exists(training_model_dir), \
            "Directory should not exist: {}. Check base_directory".format(training_model_dir)
        os.mkdir(training_model_dir)
    assert os.path.exists(training_dir), "Directory should exist: {}. Check base_directory".format(training_dir)
    assert os.path.exists(training_model_dir), \
        "Directory should exist: {}. Check base_directory".format(training_model_dir)

    if config_dict.options["test"]:
        assert not os.path.exists(testing_dir), "Directory should not exist: {}. Check base_directory".format(testing_dir)
        os.mkdir(testing_dir)
        assert not os.path.exists(all_variant_calls), "Directory should not exist: {}. Check base_directory".format(all_variant_calls)
        os.mkdir(all_variant_calls)
    assert os.path.exists(testing_dir), "Directory should exist: {}. Check base_directory".format(testing_dir)
    assert os.path.exists(all_variant_calls), "Directory should exist: {}. Check base_directory".format(all_variant_calls)

    if config_dict.options["plot_accuracies"]:
        assert not os.path.exists(testing_accuracy_dir), \
            "Directory should not exist: {}. Check base_directory".format(testing_accuracy_dir)
        assert not os.path.exists(testing_accuracy_csvs_dir), \
            "Directory should not exist: {}. Check base_directory".format(testing_accuracy_csvs_dir)
        os.mkdir(testing_accuracy_dir)
        os.mkdir(testing_accuracy_csvs_dir)
        assert not os.path.exists(training_accuracy_dir), \
            "Directory should not exist: {}. Check base_directory".format(training_accuracy_dir)
        assert not os.path.exists(training_accuracy_csvs_dir), \
            "Directory should not exist: {}. Check base_directory".format(training_accuracy_csvs_dir)
        os.mkdir(training_accuracy_dir)
        os.mkdir(training_accuracy_csvs_dir)

    if config_dict.options["plot_distributions"]:
        assert not os.path.exists(training_distributions_dir), \
            "Directory should not exist: {}. Check base_directory".format(training_distributions_dir)
        os.mkdir(training_distributions_dir)

    # Train models
    if config_dict.options["train"]:
        print("TRAINING", flush=True)
        execute = "trainModels.py run --config {}"
        assert os.path.exists(config_dict.train_model), "train_model path does not exist {}".format(config_dict.train_model)
        train_model_config = load_json(config_dict.train_model)
        assert training_dir == train_model_config["output_dir"], \
            "train_model path is not same as output_dir: {} != {}".format(training_dir, train_model_config["output_dir"])
        check_call(execute.format(config_dict.train_model).split())
        # copy training models
        print("COPYING MODELS", flush=True)
        shutil.copy2(os.path.join(training_dir, "tempFiles_trainModels/template_hmm1.model"), training_model_dir)
        for i in range(2, train_model_config["training"]["em_iterations"]+1):
            shutil.copy2(os.path.join(training_dir,
                                      "tempFiles_trainModels_{}/template_hmm{}.model".format(i, i)),
                         training_model_dir)

    if config_dict.options["test"]:
        print("TESTING", flush=True)
        # test model
        re_run_signalalign(training_model_dir,
                           testing_dir,
                           config_dict.test_model["base_model"],
                           config_dict.test_model["variants"],
                           config_dict.test_model["rna"])
        test_config = load_json(config_dict.test_model["base_model"])
        names = [sample["name"] for sample in test_config["samples"]]
        suffixes = ["variant_calls/{}.csv".format(name) for name in names]
        copy_variant_call_files(testing_dir, suffixes, all_variant_calls)

    if config_dict.options["plot_accuracies"]:
        print("PLOTTING ACCURACIES:TEST", flush=True)
        # plot accuracies
        positions_files = config_dict.plot_accuracies.test["positions_files"]
        names = config_dict.plot_accuracies.test["names"]
        suffixes = ["variant_calls/{}.csv".format(name) for name in names]
        re_run_plot_variant_accuracy(testing_accuracy_dir, testing_dir, positions_files, names, suffixes, threshold=0.5)
        # copy into one directory
        print("COPYING ACCURACY DATA:TEST", flush=True)
        for i in range(1, len(list_dir(training_model_dir))+1):
            shutil.copy2(os.path.join(testing_accuracy_dir,
                                      "template_hmm{}/per_position/per_position_data_0.5.csv".format(i)),
                         os.path.join(testing_accuracy_csvs_dir, "{}_per_position_data_0.5.csv".format(i)))

        print("PLOTTING ACCURACIES:TRAIN", flush=True)
        # plot accuracies
        positions_files = config_dict.plot_accuracies.train["positions_files"]
        names = config_dict.plot_accuracies.train["names"]
        suffixes = ["variant_calls/{}.csv".format(name) for name in names]
        re_run_plot_variant_accuracy(training_accuracy_dir, testing_dir, positions_files, names, suffixes, threshold=0.5)
        # copy into one directory
        print("COPYING ACCURACY DATA:TRAIN", flush=True)
        for i in range(1, len(list_dir(training_model_dir))+1):
            shutil.copy2(os.path.join(training_accuracy_dir,
                                      "template_hmm{}/per_position/per_position_data_0.5.csv".format(i)),
                         os.path.join(training_accuracy_csvs_dir, "{}_per_position_data_0.5.csv".format(i)))

    if config_dict.options["plot_distributions"]:
        print("PLOTTING DISTRIBUTIONS", flush=True)
        model_path = config_dict.plot_distributions
        assert(os.path.exists(training_dir))
        assert(os.path.exists(model_path))
        model_json = load_json(model_path)
        model_json["save_fig_dir"] = training_distributions_dir
        dirs_in_training = [os.path.join(training_dir, x) for x in os.listdir(training_dir) if
                            os.path.exists(os.path.join(training_dir, x)) and
                            os.path.isdir(os.path.join(training_dir, x))]
        assert len(dirs_in_training) > 0
        iterations = []
        for directory in dirs_in_training:
            number = get_trailing_number(directory)
            if number is None:
                hmm_model0 = list_dir(directory, ext="hmm")[0]
                hmm_model1 = list_dir(directory, ext="model")[0]
                sa_output_dirs = [os.path.join(directory, x) for x in os.listdir(directory) if os.path.isdir(os.path.join(directory, x))]
                iterations.append({"name": "iteration0", "hmm_model": hmm_model0, "hdp_model": None, "sa_output_dirs": sa_output_dirs})
                iterations.append({"name": "iteration1", "hmm_model": hmm_model1, "hdp_model": None, "sa_output_dirs": sa_output_dirs})
            else:
                name = "iteration{}".format(number)
                hmm_model = list_dir(directory, ext="model")
                if len(hmm_model) == 1:
                    hmm_model = hmm_model[0]
                    sa_output_dirs = [os.path.join(directory, x) for x in os.listdir(directory) if os.path.isdir(os.path.join(directory, x))]
                    iterations.append({"name": name, "hmm_model": hmm_model, "hdp_model": None, "sa_output_dirs": sa_output_dirs})

        model_json["iterations"] = iterations
        run_embed_plot_wrapper(model_json)


if __name__ == '__main__':
    main()
