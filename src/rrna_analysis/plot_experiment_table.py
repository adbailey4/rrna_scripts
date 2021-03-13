#!/usr/bin/env python
"""class for plotting accuracy of multiple model training runs"""
########################################################################
# File: plot_experiment_table.py
#  executable: plot_experiment_table.py
#
# Author: Andrew Bailey
# History: Created 12/10/20
########################################################################
import matplotlib
matplotlib.use('SVG')

import os
import re
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from itertools import zip_longest
from py3helpers.utils import list_dir, get_all_sub_directories
from py3helpers.aws import AwsS3
from rrna_analysis.multiple_model_accuracy import sort_dir, \
    preprocess_accuracy_csv
from rrna_analysis.kmer_pos_mapping import KmerPosMapping

# from IPython.display import set_matplotlib_formats
#
# set_matplotlib_formats('svg')


def find_experiments(top_dirs):
    """Find experiments in any subdirectory. An experiment must only have a 'testing_accuracy_csvs' directory

    :param top_dirs: list of paths or path to directory to search for experiments
    """
    dirs = []
    if not isinstance(top_dirs, list):
        top_dirs = [top_dirs]
    for top_dir in top_dirs:
        assert os.path.exists(top_dir), f"top_dir path does not exist:{top_dir}"
        for x in get_all_sub_directories(top_dir):
            if x.endswith("testing_accuracy_csvs") and x not in dirs:
                dirs.append(x)
    assert len(dirs) > 0, "No experiment directory found. len(dirs) == 0"
    return dirs


def get_acc_data_from_experiments(experiments, key, round_n, kpm):
    """Given a list of directories and a key, create a dataframe with the modified position information
    and key accuracy metric
    """
    data = []
    for rn, experiment in zip(round_n, experiments):
        accuracy_csv = sort_dir(experiment)[rn - 1]
        acc = preprocess_accuracy_csv(accuracy_csv, kpm.mod_data)[key].T
        # if len(acc) != 110:
        #     continue
        data.append(acc)

    split_experiments = [x.split("/") for x in experiments]
    i = 0
    for i, n in enumerate(zip_longest(*split_experiments)):
        if len(set(n)) != 1:
            break
    experiment_names = ["\n".join(["\n".join(i.split("_")) for i in x[i:-1]]) for x in split_experiments]

    final_data_frame = pd.merge(kpm.mod_handler, pd.DataFrame(np.array(data).T, columns=experiment_names),
                                left_index=True, right_index=True)
    return final_data_frame


def plot_heatmap_of_experiment(plot_df, key, urls=None, show_numbers=True, savefig=None, cmap="viridis", norm=None):
    cbarlabel = f'{key}'
    # cmap = "RdYlGn"
    # cmap = "viridis"

    not_columns = ['contig', 'position', 'strand', 'change_from', 'change_to', 'mod',
                   'pos', 'percent', 'reference_index', 'delta1_below', 'delta1_above',
                   'delta2_below', 'delta2_above', 'delta3_below', 'delta3_above',
                   'delta4_below', 'delta4_above', 'delta', 'in_2prime', 'in_pseudo',
                   'in_unknown']
    x_labels = [x for x in plot_df.columns if x not in not_columns]
    y_labels = ["_".join([str(x) for x in plot_df[["contig", "reference_index"]].loc[i]]) for i in plot_df.index]
    data = plot_df[x_labels]

    fig, ax = plt.subplots(figsize=(8, 20))
    plt.subplots_adjust(top=0.95, bottom=0.2, right=0.98, left=0.15)

    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm)
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_xticklabels(x_labels)
    ax.set_yticklabels(y_labels)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    if show_numbers:
        for x, x_label in enumerate(data.columns):
            for y, y_label in enumerate(data.index):
                #         print(data.loc[y_label, x_label])
                mod_name = "_".join([str(x) for x in plot_df[["contig", "position", "strand"]].loc[y_label]])
                url = None
                if urls:
                    url = urls[mod_name]
                text = ax.text(x, y, round(data.loc[y_label, x_label], 3),
                               ha="center", va="center", color="w",
                               url=url)

    ax.set_title("Data")
    # option to save figure or just show it
    if savefig is not None:
        plt.savefig(savefig, dpi=300)
    else:
        plt.show()
    plt.close(fig)


def check_s3_training_distributions(client, experiments):
    s3_bucket_dists = "s3://bailey-k8s/rrna_experiments/"
    check_buckets = [
        os.path.join(s3_bucket_dists, "/".join(experiment.split("/")[experiment.split("/").index("rrna_kube_testing")+1:-1]),
                     "training_distributions") for experiment in experiments]
    for bucket in check_buckets:
        #     print(bucket)
        assert (client.folder_exists(bucket[5:]))

    s3_bucket_dists = "https://bailey-k8s.s3-us-west-2.amazonaws.com/rrna_experiments/"
    experiment_buckets = [
        os.path.join(s3_bucket_dists, "/".join(experiment.split("/")[experiment.split("/").index("rrna_kube_testing")+1:-1]),
                     "training_distributions") for experiment in experiments]
    names = ["-".join(x.split("/")[-4:-1]) for x in experiment_buckets]
    return experiment_buckets, names


def create_html_file_for_mod(client, mod_title, mod_key, pos_list, experiment_buckets, names,
                             html_output_bucket="bailey-k8s/html/experiment_plotting_htmls/", overwrite=False):
    soup = BeautifulSoup("", 'html.parser')
    table = soup.new_tag("table")
    #     table["style"] = "width:100%"
    soup.insert(1, table)

    title = soup.new_tag("strong")
    title["style"] = "font-size:30px"
    title.string = mod_title
    table.insert(1, title)

    title_row = soup.new_tag("tr")
    table.insert(2, title_row)

    image_rows = []
    for i in range(2, len(pos_list) + 2):
        image_row = soup.new_tag("tr")
        table.insert(i, image_row)
        image_rows.append(image_row)

    #     print(experiment_buckets)
    for i, bucket in enumerate(experiment_buckets):
        title_tag = soup.new_tag("th")
        extra_rows = 0
        for z, row in enumerate(names[i].split("-")):
            if len(row) > 40:
                split_row = row.split("_")
                new_line = []
                row_len = 0
                counter = 0
                for split_name in split_row:
                    new_line += [split_name]
                    row_len += len(split_name)
                    if row_len > 30 or split_name == split_row[-1]:
                        par_tag = soup.new_tag("p")
                        par_tag.string = "_".join(new_line)
                        title_tag.insert(z+1+extra_rows+counter, par_tag)
                        new_line = []
                        row_len = 0
                        counter += 1
                extra_rows += (counter - 1)
            else:
                par_tag = soup.new_tag("p")
                par_tag.string = row
                title_tag.insert(z+1+extra_rows, par_tag)

        #         title_tag.string = names[i]
        title_row.insert(i + 1, title_tag)
        for j, pos in enumerate(pos_list):
            row_tag = soup.new_tag("td")
            row_tag["style"] = "border: 1px solid black"
            image_rows[j].insert(i + 1, row_tag)

            par_tag = soup.new_tag("p")
            par_tag.string = str(pos)
            row_tag.insert(1, par_tag)

            image_tag = soup.new_tag("img")
            image_tag["width"] = 500
            image_tag['src'] = os.path.join(bucket, mod_key.replace("+", "%2B"), str(pos) + ".gif")
            row_tag.insert(2, image_tag)

    with open("output.html", "w") as fh:
        print(soup.prettify(), file=fh)
    #         print(soup.prettify())
    out_file_name = "_".join([str(hash("_".join(names))), mod_key]) + ".html"
    s3_file = os.path.join(html_output_bucket, out_file_name)
    if not client.object_exists(s3_file) or overwrite:
        _ = client.upload_object(file_path="output.html", destination=s3_file, use_original_name=False,
                                 ExtraArgs={'ContentType': "text/html", 'ACL': "public-read"})
    #     need to replace +
    s3_gif_url = s3_file.replace("bailey-k8s/", "https://bailey-k8s.s3-us-west-2.amazonaws.com/")
    s3_gif_url = s3_gif_url.replace("+", "%2B")
    os.remove("output.html")
    return mod_key, s3_gif_url


def create_html_distributions(experiments, plot_df, client=None, overwrite=False):
    if client is None:
        client = AwsS3()
    try:
        experiment_buckets, names = check_s3_training_distributions(client, experiments)
    except AssertionError as e:
        print(e)
        return None
    url_dict = {"_".join([str(x) for x in list(plot_df[["contig", "position", "strand"]].loc[i])]): "" for i in
                plot_df.index}
    mod_title_dict = {"_".join([str(x) for x in list(plot_df[["contig", "position", "strand"]].loc[i])]):
                          "_".join([str(x) for x in list(plot_df[["contig", "position", "strand", "mod",
                                                                  "percent", "delta"]].loc[i])]) for i in
                      plot_df.index}

    for mod_key in url_dict.keys():
        mod_pos = int(mod_key.split("_")[1])
        pos_list = list(range(mod_pos - 4, mod_pos + 1))
        mod_key, s3_gif_url = create_html_file_for_mod(client, mod_title_dict[mod_key], mod_key, pos_list,
                                                       experiment_buckets, names, overwrite=overwrite)
        url_dict[mod_key] = s3_gif_url
    return url_dict


def plot_acc_heatmap_for_experiment(dirs, key, kpm, min_percent=90, max_percent=100, max_delta=np.inf, min_delta=6,
                                    round_n=30, show_numbers=True, client=None, overwrite=False, savefig=None,
                                    cmap="viridis", norm=None):
    experiments = find_experiments(dirs)
    if not isinstance(round_n, list):
        round_n = [round_n for i in range(len(experiments))]
    assert len(round_n) == len(experiments), f"len(round_n) != len(experiments): {len(round_n)} != {len(experiments)}"
    final_data_frame = get_acc_data_from_experiments(experiments, key, round_n, kpm)
    print(f"Got final_data_frame: len={len(final_data_frame)}")
    plot_df = final_data_frame[(final_data_frame["percent"] >= min_percent) &
                               (final_data_frame["percent"] <= max_percent) &
                               (final_data_frame["delta"] >= min_delta) &
                               (final_data_frame["delta"] <= max_delta)]
    print(f"Got plot_df: len={len(plot_df)}")
    mod_s3_urls = create_html_distributions(experiments, plot_df, client=client, overwrite=overwrite)
    if mod_s3_urls is not None:
        print(f"Got mod_s3_urls: len={len(mod_s3_urls)}")
    else:
        print("No mod_s3_urls")
    plot_heatmap_of_experiment(plot_df, key, mod_s3_urls, show_numbers=show_numbers, savefig=savefig, cmap=cmap,
                               norm=norm)
