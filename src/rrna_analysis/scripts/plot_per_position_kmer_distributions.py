#!/usr/bin/env python
"""Plot per position kmer distributions"""
########################################################################
# File: plot_per_position_kmer_distributions.py
#  executable: plot_per_position_kmer_distributions.py
#
# Author: Andrew Bailey
# History: 06/04/20 Created
########################################################################

import os
import sys
from collections import defaultdict
from itertools import zip_longest
from argparse import ArgumentParser
import platform
from multiprocessing import Manager, Process, current_process

import pandas as pd
import numpy as np
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
import matplotlib as mpl

if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
if platform.system() == "Darwin":
    mpl.use("macosx")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.animation import FuncAnimation

from py3helpers.utils import create_dot_dict, load_json, time_it, list_dir, merge_lists
from py3helpers.multiprocess import run_service
from py3helpers.seq_tools import ReverseComplement

from signalalign.utils.sequenceTools import CustomAmbiguityPositions
from signalalign.hiddenMarkovModel import HmmModel


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to config file")
    args = parser.parse_args()
    return args


def get_kde_from_data(data, bandwidth=1):
    x = np.asarray(data).reshape(len(data), 1)
    x_plot = np.linspace(30, 200, 1000)[:, np.newaxis]
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(x)
    log_dens = kde.score_samples(x_plot)
    x = x_plot[:, 0]
    y = np.exp(log_dens)
    return [x, y]


def parse_hist_output(hist_path):
    with open(hist_path, 'r') as fh:
        header = [float(x) for x in fh.readline().rstrip().split(",")]
        x = list(np.linspace(header[0], header[1], int(header[2])))
        y = [[x[index]] * int(count) for index, count in enumerate(fh.readline().rstrip().split(",")) if count != "0"]
        data = merge_lists(y)
    return data


def parse_hist_output2(hist_path):
    data = defaultdict()
    with open(hist_path, 'r') as fh:
        header = [float(x) for x in fh.readline().rstrip().split(",")]
        x = list(np.linspace(header[0], header[1], int(header[2])))
        n_lines = int(header[4])
        min_prob = header[3]
        for _ in range(n_lines):
            name = fh.readline().rstrip()
            l = merge_lists([[x[index]] * int(count) for index, count in enumerate(fh.readline().rstrip().split(",")) if count != "0"])
            data[name] = l
    return data, min_prob


def get_number_of_lines(hist_path):
    with open(hist_path, 'r') as fh:
        header = [float(x) for x in fh.readline().rstrip().split(",")]
    return int(header[4])


def get_paired_colors(num_needed):
    cmap1 = mpl.cm.get_cmap("tab20")
    cmap2 = mpl.cm.get_cmap("tab20b")
    colors = [cmap1(x) for x in np.linspace(0, 1, num=20)]
    colors.extend([cmap2(x) for x in np.linspace(0, 1, num=10)])
    if num_needed > 15:
        cmap3 = mpl.cm.get_cmap("hsv")
        new_colors = [cmap3(x) for x in np.linspace(0, 1, num=num_needed - 15)]
        faded_new_colors = [x[:3] + tuple([0.6]) for x in new_colors]
        # colors.extend([j for i in zip(new_colors, faded_new_colors) for j in i])
        colors.extend([j for i in zip(new_colors, faded_new_colors) for j in i])
    return colors


def plot_kmer_gif2(per_pos_hist_data, hmm_model_params, output_file, strand="t", mod_only=False):
    hmm_models = []
    for hmm_model_params in hmm_model_params:
        hmm_model = None
        if hmm_model_params is not None:
            hmm_model = HmmModel(**hmm_model_params)
        hmm_models.append(hmm_model)
    plot_kmer_gif(per_pos_hist_data, hmm_models, output_file, strand=strand, mod_only=mod_only)


def plot_kmer_gif(pos_hist_data, hmm_models, output_file, strand="t", mod_only=False):
    n_images = len(pos_hist_data)
    n_models = len(hmm_models)
    canonical_bases = set("ATGC")
    assert n_images == n_models, "n_images != n_models"
    n_lines = get_number_of_lines(pos_hist_data[0][2])
    # print("n_lines: ", n_lines)
    colors = get_paired_colors(2*n_lines)

    min_x = 1000
    max_x = 0
    # min_x = 0
    # max_x = 200
    rc = ReverseComplement()

    position_main = []
    kmer_main = []
    model_main = []
    for hist_data, model in zip_longest(pos_hist_data, hmm_models):
        model_list = []
        kmer_list = []
        position_list = []
        n_positions = 0
        n_kmers = 0
        n_plot_models = 0
        data, min_prob = parse_hist_output2(hist_data[2])
        for i, (pos_kmer, kmer_data) in enumerate(data.items()):
            # print(pos_kmer.split("_"))
            # reset min_x max_x
            tmp_min_x = min([kmer_data[0]])
            tmp_max_x = max([kmer_data[-1]])
            if min_x > tmp_min_x:
                min_x = tmp_min_x
            if max_x < tmp_max_x:
                max_x = tmp_max_x
            # reset params

            kmer_x, kmer_y, kmer_name = None, None, None
            model_x, model_y, model_name = None, None, None
            pos_x, pos_y, pos_name = None, None, None
            # if kmer dist
            if len(pos_kmer.split("_")) == 1:
                kmer_x, kmer_y = get_kde_from_data(kmer_data)
                kmer_name = pos_kmer + ": N=" + str(len(kmer_data))
                n_kmers += 1
                if model is not None:
                    n_plot_models += 1
                    nuc_type = "RNA" if model.rna else "DNA"
                    strand = "t" if strand is None else strand

                    normal_mean, normal_sd = model.get_event_mean_gaussian_parameters(pos_kmer)
                    model_name = model.name + ": " + "_".join([nuc_type, strand, pos_kmer]) + " N(" + str(
                        round(normal_mean, ndigits=2)) + ", " + str(round(normal_sd, ndigits=2)) + ")"

                    # plot ont normal distribution
                    model_x = np.linspace(normal_mean - 4 * normal_sd, normal_mean + 4 * normal_sd, 200)
                    model_y = norm.pdf(model_x, normal_mean, normal_sd)
                    # plot HDP distribution
                    if model.has_hdp_model:
                        model_x = model.linspace
                        kmer_id = model.get_kmer_index(pos_kmer)
                        model_y = model.all_posterior_pred[kmer_id]
                        model_name = f"{model.name} HDP: {'_'.join([nuc_type, strand, pos_kmer])} "
                    tmp_min_x = normal_mean - (5 * normal_sd)
                    tmp_max_x = normal_mean + (5 * normal_sd)
                    if min_x > tmp_min_x:
                        min_x = tmp_min_x
                    if max_x < tmp_max_x:
                        max_x = tmp_max_x
            else:
                if not mod_only or not set(pos_kmer.split("_")[-1]) <= canonical_bases:
                    pos_x, pos_y = get_kde_from_data(kmer_data)
                    pos_name = pos_kmer + ": N=" + str(len(kmer_data))
                    n_positions += 1

            kmer_color = colors[n_kmers * 2]
            model_color = colors[n_kmers * 2]
            pos_color = colors[(n_kmers * 2) + 1]

            model_list.append([model_x, model_y, model_color, model_name])
            kmer_list.append([kmer_x, kmer_y, kmer_color, kmer_name])
            position_list.append([pos_x, pos_y, pos_color, pos_name])

        position_main.append(position_list)
        kmer_main.append(kmer_list)
        model_main.append(model_list)

    fig = plt.figure(figsize=(10, 10))
    # panel1 = plt.axes([0.1, 0.1, .8, .8])
    panel1 = plt.axes([0.13, 0.5, .8, .45])
    # panel1 = plt.axes()
    panel1.set_xlabel('pA')
    panel1.set_ylabel('Density')
    panel1.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)
    panel1.xaxis.set_major_locator(ticker.AutoLocator())
    panel1.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    panel1.set_xlim(min_x, max_x)
    panel1.set_ylim(0, 0.5)
    panel1.set_title("Kmer distribution comparisons")
    linestyle_tuple = [
        (0, (1, 10)),
        (0, (1, 1)),
        (0, (5, 10)),
        (0, (5, 5)),
        (0, (5, 1)),
        (0, (3, 10, 1, 10)),
        (0, (3, 5, 1, 5)),
        (0, (3, 1, 1, 1)),
        (0, (3, 5, 1, 5, 1, 5)),
        (0, (3, 10, 1, 10, 1, 10)),
        (0, (3, 1, 1, 1, 1, 1))]

    pos_lines = [panel1.plot([], [], linestyle=linestyle_tuple[i % 11], lw=1, label="label")[0] for i in range(n_positions)]
    model_lines = [panel1.plot([], [], "--", lw=2, label="label")[0] for _ in range(n_plot_models)]
    kmer_lines = [panel1.plot([], [], lw=2, label="label")[0] for _ in range(n_kmers)]

    # legend = [panel1.legend(loc='upper left')]
    legend = panel1.legend(loc='upper left', bbox_to_anchor=(0.0, -0.1),
                            handles=merge_lists([pos_lines, model_lines, kmer_lines]))
    # legend = [panel1.legend(loc='upper left', bbox_to_anchor=(0.0, -0.1), handles=pos_lines),
    #           panel1.legend(loc='upper left', handles=kmer_lines),
    #           panel1.legend(loc='upper right', handles=model_lines)]

    ax = plt.gca().add_artist(legend)
    # model_x, model_y, model_color, model_name

    def animate(i):
        pos_index = 0
        model_index = 0
        kmer_index = 0
        for lnum in range(n_lines):
            model_x, model_y, model_color, model_name = model_main[i][lnum]
            kmer_x, kmer_y, kmer_color, kmer_name = kmer_main[i][lnum]
            pos_x, pos_y, pos_color, pos_name = position_main[i][lnum]
            if pos_x is not None:
                pos_lines[pos_index].set_data(pos_x, pos_y)  # set data for each line separately.
                pos_lines[pos_index].set_color(pos_color)
                pos_lines[pos_index].set_label(pos_name)
                legend.texts[pos_index].set_text(pos_name)
                legend.legendHandles[pos_index].set_color(pos_color)
                pos_index += 1
            if model_x is not None:
                model_lines[model_index].set_data(model_x, model_y)  # set data for each line separately.
                model_lines[model_index].set_color(model_color)
                model_lines[model_index].set_label(model_name)
                legend.texts[n_positions+model_index].set_text(model_name)
                legend.legendHandles[n_positions+model_index].set_color(model_color)
                model_index += 1
            if kmer_x is not None:
                kmer_lines[kmer_index].set_data(kmer_x, kmer_y)  # set data for each line separately.
                kmer_lines[kmer_index].set_color(kmer_color)
                kmer_lines[kmer_index].set_label(kmer_name)
                legend.texts[n_kmers+n_positions+kmer_index].set_text(kmer_name)
                legend.legendHandles[n_kmers+n_positions+kmer_index].set_color(kmer_color)
                kmer_index += 1
        return pos_lines + model_lines + kmer_lines + [legend]
    # for i in range(n_models):
    #     animate(i)

    anim = FuncAnimation(fig, animate, frames=list(range(n_models)), interval=400, blit=False)

    if output_file:
        anim.save(output_file, writer='imagemagick', fps=2, dpi=200)
        # anim.save(output_file, writer="imagemagick")
    else:
        plt.show()
    plt.close(fig)


def get_hist_data(hist_dir, contig, position, strand):
    positions_files = list_dir(os.path.join(hist_dir, "positions", contig + strand + str(position)))
    positions = {int(positions_files[i].split("_")[-1].split(".")[0]): [] for i in range(len(positions_files))}
    for pos_file in positions_files:
        split_file_path = pos_file.split("_")
        tmp_pos = int(split_file_path[-1].split(".")[0])
        kmer = split_file_path[-1][:-4]
        assert (os.path.exists(pos_file))
        positions[tmp_pos].append((contig, kmer, pos_file))

    return positions


def check_config(config):
    assert os.path.exists(config.save_fig_dir), "{} does not exist".format(config.save_fig_dir)
    assert os.path.exists(config.positions), "{} does not exist".format(config.positions)
    assert isinstance(config.num_threads, int), "Num threads parameter needs to be an int"
    assert isinstance(config.debug, bool), "debug parameter needs to be an boolean"
    assert isinstance(config.rna, bool), "rna parameter needs to be an boolean"
    assert isinstance(config.mod_only, bool), "mod_only parameter needs to be an boolean"

    for data in config.iterations:
        hmm_model = data["hmm_model"]
        hdp_model = data["hdp_model"]
        kmer_hist_dir = data["kmer_hist_dir"]
        if hmm_model:
            assert os.path.exists(hmm_model), "{} does not exist".format(hmm_model)
        if hdp_model:
            assert os.path.exists(hdp_model), "{} does not exist".format(hdp_model)
        assert os.path.exists(kmer_hist_dir), "{} does not exist".format(kmer_hist_dir)


def run_plotting_worker(output_file, per_pos_hist_data, hmm_models, strand, mod_only):
    # plot_kmer_gif2(per_pos_hist_data, hmm_model_params, output_file, strand=strand, mod_only=mod_only)
    print(output_file)
    plot_kmer_gif(per_pos_hist_data, hmm_models, output_file, strand=strand, mod_only=mod_only)


class BasicService2(object):

    def __init__(self, function, models, service_name="example_service"):
        self.function = function
        self.service_name = service_name
        self.models = models

    def run(self, work_queue, done_queue):

        """
        Service used by the multiprocess module for an example of what should be created for multiprocessing
        :param work_queue: arguments to be done
        :param done_queue: errors and returns to be put
        :param service_name: name of the service
        """
        # prep
        total_handled = 0
        failure_count = 0

        # catch overall exceptions
        try:
            for f in iter(work_queue.get, 'STOP'):
                # catch exceptions on each element
                try:
                    return_value = self.function(hmm_models=self.models, **f)
                    done_queue.put(return_value)
                except Exception as e:
                    # get error and log it
                    message = "{}:{}".format(type(e), str(e))
                    error = "{} '{}' failed with: {}".format(self.service_name, current_process().name, message)
                    print("[{}] ".format(self.service_name) + error)
                    done_queue.put(error)
                    failure_count += 1

                # increment total handling
                total_handled += 1

        except Exception as e:
            # get error and log it
            message = "{}:{}".format(type(e), str(e))
            error = "{} '{}' critically failed with: {}".format(self.service_name, current_process().name, message)
            print("[{}] ".format(self.service_name) + error)
            done_queue.put(error)

        finally:
            # logging and final reporting
            print("[%s] '%s' completed %d calls with %d failures"
                  % (self.service_name, current_process().name, total_handled, failure_count))
            done_queue.put("{}:{}".format("total", total_handled))
            done_queue.put("{}:{}".format("failure", failure_count))


def plot_per_position_kmer_distributions(config_dict):
    config = create_dot_dict(config_dict)
    check_config(config)
    print("Output dir: ", config.save_fig_dir)
    print("Positions: ", config.positions)
    print("RNA: ", config.rna)
    print("Model/Data: ")

    pos_data = CustomAmbiguityPositions.parseAmbiguityFile(config.positions)
    # contig         pUC19
    # position         709
    # strand             +
    # change_from        G
    # change_to          G
    strand = "t"

    hmm_models = []
    for iteration_data in config.iterations:
        hmm_model = None
        if iteration_data["hmm_model"]:
            # if not multiprocess:
            hmm_model = HmmModel(iteration_data["hmm_model"], hdp_model_file=iteration_data["hdp_model"],
                                 rna=config.rna, name=iteration_data["name"])
        # if not multiprocess:
        hmm_models.append(hmm_model)

    pos_hist_data = defaultdict(list)
    for i in range(len(pos_data)):
        pos_row = pos_data.iloc[i]
        output_dir = os.path.join(config.save_fig_dir, "_".join([pos_row.contig,
                                                                 str(pos_row.position),
                                                                 pos_row.strand]))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for iteration_data in config.iterations:
            hist_data = get_hist_data(iteration_data["kmer_hist_dir"], pos_row.contig, pos_row.position, pos_row.strand)
            for position, tuple_data in hist_data.items():
                assert len(tuple_data) == 1, "Something is wrong with Histogram directory"
                pos_hist_data[os.path.join(output_dir, str(position) + ".gif")].append(tuple_data[0])

    if not config.debug:
        extra_args = {"strand": strand,
                      "mod_only": config.mod_only}

        output_files = list(pos_hist_data.keys())
        per_pos_hist_data = list(pos_hist_data.values())
        service = BasicService2(run_plotting_worker, hmm_models, service_name="multiprocess_run_plotting_worker")
        run_service(service.run, zip_longest(output_files, per_pos_hist_data), extra_args,
                    ["output_file", "per_pos_hist_data"], config.num_threads)
    else:
        for output_file, each_iteration_data in pos_hist_data.items():
            print(output_file)
            plot_kmer_gif(each_iteration_data, hmm_models, output_file, strand=strand, mod_only=config.mod_only)
            break


def main():
    args = parse_args()
    assert os.path.exists(args.config), "{} does not exist".format(args.config)
    plot_per_position_kmer_distributions(load_json(args.config))


if __name__ == '__main__':
    print(time_it(main)[1])
