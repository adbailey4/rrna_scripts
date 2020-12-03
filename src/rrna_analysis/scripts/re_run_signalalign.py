#!/usr/bin/env python
"""Re-run signalalign multiple times for each model"""
########################################################################
# File: re_run_signalalign.py
#  executable: re_run_signalalign.py
#
# Author: Andrew Bailey
# History: 05/04/20 Created
########################################################################

import os
from argparse import ArgumentParser
from py3helpers.utils import list_dir, load_json, time_it, save_json
from subprocess import check_call


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--dir', '-d', action='store',
                        dest='dir', required=True, type=str, default=None,
                        help="Path to directory of model files")
    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True, type=str, default=None,
                        help="Path to directory to output signalalign runs")
    parser.add_argument('--base_model', '-m', action='store',
                        dest='base_model', required=True, type=str, default=None,
                        help="Base model file for signalalign")
    parser.add_argument('--variants', '-v', action='store',
                        dest='variants', required=True, type=str, default="Y",
                        help="variants to analyze")
    parser.add_argument('--rna', action='store_true', dest='rna',
                        help="set if rna reads")

    args = parser.parse_args()
    return args


def re_run_signalalign(directory, output_dir, base_model, variants, rna):
    assert os.path.isdir(directory), "{} is not a directory".format(directory)
    assert os.path.isdir(output_dir), "{} is not a directory".format(output_dir)
    assert os.path.exists(base_model), "{} does not exist".format(base_model)
    models = list_dir(directory, ext="model")
    sa_base_model = load_json(base_model)
    created_models_dir = os.path.join(output_dir, "created_models")
    if not os.path.exists(created_models_dir):
        os.mkdir(created_models_dir)
    execute = "runSignalAlign.py run --config {}"
    embed_main_execute = "embed_main sa2bed -d {}/tempFiles_alignment/{}/ -a {} -o {}/{}.bed -t {} -c {} --overwrite"
    if rna:
        embed_main_execute += " --rna"
    for model in models:
        # run sa
        try:
            base_name = os.path.splitext(os.path.basename(model))[0]
            round_output_dir = os.path.join(output_dir, base_name)
            if not os.path.exists(round_output_dir):
                sa_base_model["template_hmm_model"] = model
                sa_base_model["output_dir"] = round_output_dir
                new_model_file = os.path.join(created_models_dir, base_name+".json")
                save_json(sa_base_model, new_model_file)
                check_call(execute.format(new_model_file).split())
            else:
                print(round_output_dir, " exists already: continuing")
                continue
            #       run sa2bed
            variants_dir = os.path.join(round_output_dir, "variant_calls")
            if not os.path.exists(variants_dir):
                os.mkdir(variants_dir)
            for sample in sa_base_model["samples"]:
                check_call(embed_main_execute.format(round_output_dir, sample["name"], sa_base_model["ambig_model"], variants_dir, sample["name"], sa_base_model["job_count"], variants).split())
        except Exception as e:
            print(e)
            continue


def main():
    args = parse_args()
    print(args)
    re_run_signalalign(args.dir, args.output_dir, args.base_model, args.variants, args.rna)


if __name__ == '__main__':
    print(time_it(main)[1])
