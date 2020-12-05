#!/usr/bin/env python
"""class for plotting accuracy of multiple model training runs"""
########################################################################
# File: multiple_model_accuracy.py
#  executable: multiple_model_accuracy.py
#
# Author: Andrew Bailey
# History: Created 12/03/20
########################################################################

import re
import os
from scipy.stats import norm, entropy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from collections import defaultdict

from py3helpers.utils import list_dir
from signalalign.hiddenMarkovModel import HmmModel
from rrna_analysis.kmer_pos_mapping import KmerPosMapping


def get_first_int(string):
    return int(re.search(r'\d+', string).group())


def get_kmer_kl_divergence(model, kmer1, kmer2):
    """Get Kullback–Leibler divergence between the HDP and ONT models for a specific kmer"""
    normal_mean1, normal_sd1 = model.get_event_mean_gaussian_parameters(kmer1)
    normal_mean2, normal_sd2 = model.get_event_mean_gaussian_parameters(kmer2)

    ont_normal_dist1 = norm.pdf(model.linspace, normal_mean1, normal_sd1)
    ont_normal_dist2 = norm.pdf(model.linspace, normal_mean2, normal_sd2)

    kl_divergence = entropy(pk=ont_normal_dist2, qk=ont_normal_dist1, base=2)

    return kl_divergence


def get_kmer_mean_delta(model, kmer1, kmer2):
    """Get Kullback–Leibler divergence between the HDP and ONT models for a specific kmer"""
    normal_mean1, normal_sd1 = model.get_event_mean_gaussian_parameters(kmer1)
    normal_mean2, normal_sd2 = model.get_event_mean_gaussian_parameters(kmer2)

    return np.abs(normal_mean1 - normal_mean2)


# Index(['contig', 'reference_index', 'strand', 'variants', 'accuracy',
#        'precision', 'negative_predictive_value', 'recall', 'specificity',
#        'positive_likelihood_ratio', 'negative_likelihood_ratio',
#        'diagnostic_odds_ratio', 'f1_score', 'prevalence', 'aucroc',
#        'avg_precision', 'brier_score'],
#       dtype='object')


def preprocess_accuracy_csv(path, mod_data):
    assert (os.path.exists(path))
    accuracy_csv = pd.read_csv(path).sort_values(by=['contig', 'reference_index'])

    accuracy_csv['delta1'] = accuracy_csv.reference_index.diff().shift(-1)
    accuracy_csv['delta2'] = np.abs(accuracy_csv.reference_index.diff().shift(0))
    accuracy_csv['delta'] = accuracy_csv[["delta1", "delta2"]].min(axis=1)
    accuracy_csv["in_2prime"] = (((accuracy_csv.variants.shift().isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                   (accuracy_csv.delta2 <= 5)) |
                                  (accuracy_csv.variants.shift(-1).isin(["Aa", "Cb", "Gc", "Td", "Tdm"]) &
                                   (accuracy_csv.delta1 <= 5))) & (
                                     ~accuracy_csv.variants.isin(["Aa", "Cb", "Gc", "Td"])))
    accuracy_csv["in_pseudo"] = (((accuracy_csv.variants.shift().isin(["Tl"]) &
                                   (accuracy_csv.delta2 <= 5)) |
                                  (accuracy_csv.variants.shift(-1).isin(["Tl", "Tdm"]) &
                                   (accuracy_csv.delta1 <= 5))) &
                                 (~accuracy_csv.variants.isin(["Tl"])))
    accuracy_csv["in_unknown"] = (accuracy_csv["in_pseudo"] | accuracy_csv["in_2prime"])
    accuracy_csv = pd.merge(accuracy_csv, mod_data, on=["contig", "reference_index"])
    return accuracy_csv


def csv_model_mapping(csv_dir, model_dir):
    """Map accuracy csv file with correct model
    :param csv_dir: path to accuracy csv dir
    :param model_dir: path to model dir
    """
    csvs = list_dir(csv_dir, "csv")
    models = list_dir(model_dir, "model")
    assert len(csvs) == len(models), f"len(csvs) != len(models): {len(csvs)} != {len(models)}"
    csv_model_map = defaultdict(list)
    for c in csvs:
        csv_model_map[get_first_int(os.path.basename(c))].append(c)
    for m in models:
        csv_model_map[get_first_int(os.path.basename(m))].append(m)
    return csv_model_map


def sort_dir(csv_dir, ext="csv"):
    """Sort files by the first number in the basename"""
    files = list_dir(csv_dir, ext)
    csv_model_map = []
    for f in files:
        csv_model_map.append((get_first_int(os.path.basename(f)), f))
    return [x for (i, x) in sorted(csv_model_map, key=lambda x: x[0])]


def plot_accuracy_vs_delta_and_accuracy_over_time(kmer_pos_mapping, directory, model_dir, model_n, high_percent=100,
                                                  low_percent=0, low_delta=0, high_delta=np.inf, key="accuracy",
                                                  max_delta=False, aot=True, avd=True):
    """Plot accuracy vs model delta from canonical"""
    csv_model_map = csv_model_mapping(directory, model_dir)

    plot_me2 = defaultdict(list)
    plot_me = defaultdict(defaultdict)
    model = None
    for i, (x, model_path) in csv_model_map.items():
        assert (os.path.exists(x))
        round_data = preprocess_accuracy_csv(x, kmer_pos_mapping.subset_mod_data)
        round_number = get_first_int(os.path.basename(x))
        #     print(round_number)
        if i + 1 == model_n:
            model = HmmModel(model_path, rna=True)
        for group, data in round_data.groupby(["contig", "reference_index", "strand"]):
            mod = "_".join([str(x) for x in group])
            if (low_percent <= data["percent"].iloc[0] <= high_percent and
                    low_delta <= data["delta"].iloc[0] <= high_delta):
                plot_me[mod][round_number] = data[key].iloc[0]
                # other plot
                if i + 1 == model_n:
                    kmer_pairs = kmer_pos_mapping.get_kmers_covering_mod(contig=str(group[0]),
                                                                         strand=str(group[2]),
                                                                         position=int(group[1]))
                    deltas = []
                    assert kmer_pairs, f"missing contig, strand, position: contig, strand, position: " \
                                       f"{group[0]}, {group[2]}, {group[1]}"
                    for kmer_pair in kmer_pairs:
                        kmer_pair = list(kmer_pair)
                        if len(kmer_pair) > 2:
                            subset_kmer_deltas = []
                            subset_kmer_pairs = list(combinations(kmer_pair, 2))
                            for k1, k2 in subset_kmer_pairs:
                                subset_kmer_deltas.append(get_kmer_mean_delta(model, k1, k2))
                            deltas.append(max(subset_kmer_deltas))
                        else:
                            deltas.append(get_kmer_mean_delta(model, kmer_pair[0], kmer_pair[1]))
                    if max_delta:
                        delta = max(deltas)
                    else:
                        delta = sum(deltas)
                    plot_me2[mod] = [data[key].iloc[0], delta]
    if aot:
        plt, keys = plot_fig_1(plot_me, key)
    if avd:
        plt, keys = plot_me_fig2(plot_me2, key, max_delta)
    return plt, keys


def plot_fig_1(plot_me, key):
    #     fig 1
    fig = plt.figure(figsize=[10, 15])
    panel1 = plt.axes([0.1, 0.5, .6, .4])
    panel1.set_title('{} over training runs'.format(key))
    panel1.set_xlabel('Training Round')
    panel1.set_ylabel(key)

    lines = []
    for label, data_dict in plot_me.items():
        sorted_items = sorted(data_dict.items())
        l, = panel1.plot([x[0] for x in sorted_items], [x[1] for x in sorted_items], label=label)
        lines.append(l)

    panel1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
                  fancybox=True, shadow=True, ncol=4)
    annot = panel1.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(line, idx):
        posx, posy = [line.get_xdata()[idx], line.get_ydata()[idx]]
        annot.xy = (posx, posy)
        text = f'{line.get_label()}: {posx:.2f}-{posy:.2f}'
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == panel1:
            for line in lines:
                cont, ind = line.contains(event)
                if cont:
                    update_annot(line, ind['ind'][0])
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)
    return plt, list(plot_me.keys())


def plot_me_fig2(plot_me, key, max_delta):
    fig = plt.figure(figsize=[10, 15])
    panel1 = plt.axes([0.1, 0.5, .6, .4])
    if max_delta:
        panel1.set_title(f'{key} vs Max delta between modified and canonical')
        panel1.set_xlabel("Max delta between modified and canonical")
    else:
        panel1.set_title(f'{key} vs Sum deltas between modified and canonical')
        panel1.set_xlabel("Sum deltas between modified and canonical")

    panel1.set_ylabel(key)

    lines = []
    accuracy_values = []
    delta_values = []
    for label, (acc, delta) in plot_me.items():
        accuracy_values.append(acc)
        delta_values.append(delta)
        l, = panel1.plot(delta, acc, 'o', label=label)
        lines.append(l)

    m, b = np.polyfit(delta_values, accuracy_values, 1)
    panel1.plot(delta_values, m * np.array(delta_values) + b,
                label=f"slope={round(m, 4)} \n intercept={round(b, 4)}")

    panel1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
                  fancybox=True, shadow=True, ncol=4)
    annot = panel1.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(line, idx):
        posx, posy = [line.get_xdata()[idx], line.get_ydata()[idx]]
        annot.xy = (posx, posy)
        text = f'{line.get_label()}: {posx:.2f}-{posy:.2f}'
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == panel1:
            for line in lines:
                cont, ind = line.contains(event)
                if cont:
                    update_annot(line, ind['ind'][0])
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)

    return plt, list(plot_me.keys())
