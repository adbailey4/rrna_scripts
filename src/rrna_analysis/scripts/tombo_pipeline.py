#!/usr/bin/env python
"""Run end to end tombo pipeline for yeast rRNA analysis"""
########################################################################
# File: tombo_pipeline.py
#  executable: tombo_pipeline.py
#
# Author: Andrew Bailey
# History: 06/28/21 Created
########################################################################

import os
from argparse import ArgumentParser
from subprocess import check_call, Popen, check_output, PIPE
import urllib.request
from pathlib import Path
import shutil

from py3helpers.utils import list_dir, load_json, time_it, save_json, list_dir_recursive


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--fastq', action='store',
                        dest='fastq', required=True, type=str, default=None,
                        help="Path to fastq")
    parser.add_argument('--fast5', action='store',
                        dest='fast5', required=False, type=str, default=None,
                        help="Path to fast5 directory")
    parser.add_argument('--reference', '-r', action='store',
                        dest='reference', required=True, type=str, default=None,
                        help="Path to reference fasta")
    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=False, type=str, default=".",
                        help="Path to output_dir")
    parser.add_argument('--threads', '-t', action='store',
                        dest='threads', required=False, type=int, default=1,
                        help="Number of threads")
    parser.add_argument('--seq_summary', action='store',
                        dest='seq_summary', required=False, type=str, default=None,
                        help="Path to sequence summary file")
    parser.add_argument('--name', action='store',
                        dest='name', required=True, type=str,
                        help="Name of experiment")
    parser.add_argument('--wt_fast5s', action='store',
                        dest='wt_split_fast5s', required=True, type=str,
                        help="Wild type fast5 files to compare with experiment")

    args = parser.parse_args()
    return args


def align_and_filter(fastq, reference, threads=1):
    """Align reads using minimap and filter bam files"""
    out_bam = os.path.splitext(fastq)[0] + ".bam"
    filtered_sorted_bam = os.path.splitext(fastq)[0] + ".2308.sorted.bam"
    p1 = Popen(f"minimap2 --MD -t {threads} -ax map-ont {reference} {fastq}".split(), stdout=PIPE)
    p2 = Popen(f"samtools sort -@ {threads} ".split(), stdin=p1.stdout, stdout=PIPE)
    p1.stdout.close()
    p3 = Popen(f"samtools view -h -@ {threads} -o {out_bam}".split(), stdin=p2.stdout, stdout=PIPE)
    p2.wait()
    p2.stdout.close()
    p3.wait()
    rcode = p3.returncode
    assert rcode == 0, "Return code is not 0, check input paths and if both minimap2 and samtools are installed"
    check_call(f"samtools index -@ {threads} {out_bam}".split())
    p4 = Popen(f"samtools view -@ {threads} -bSF 2308 {out_bam} -o {filtered_sorted_bam}".split(), stdout=PIPE)
    p4.wait()
    rcode = p4.returncode
    assert rcode == 0, "Return code is not 0, check input paths and if both minimap2 and samtools are installed"
    check_call(f"samtools index -@ {threads} {filtered_sorted_bam}".split())
    return out_bam, filtered_sorted_bam


def run_qc(seq_summary, bam, html):
    p1 = Popen(f"pycoQC -f {seq_summary} -a {bam} -o {html}".split())
    p1.wait()
    rcode = p1.returncode
    assert rcode == 0, "Return code is not 0, check input paths and if pycoQC is installed"
    return html


def split_fast5s(fast5_dir, output_dir, threads=1):
    check_call(f"multi_to_single_fast5 --t {threads} "
               f"--input_path {fast5_dir} --save_path {output_dir}".split())


def concat_fastq_files(list_of_paths, output_file_path):
    """Concat all fastq files together"""
    with open(output_file_path, 'wb') as wfd:
        for f in list_of_paths:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    return output_file_path


def main():
    args = parse_args()
    print("Align and filter BAM")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    name = args.name
    if os.path.isdir(args.fastq):
        fastqs = [x for x in list_dir_recursive(args.fastq, ext="fastq")]
        output_file_path = os.path.join(output_dir, args.name + ".fastq")
        assert len(fastqs) >= 1, f"No fastqs found in {args.fastq}"
        fastq = concat_fastq_files(fastqs, output_file_path)
    else:
        fastq = args.fastq
    out_bam, filtered_sorted_bam = align_and_filter(fastq, args.reference, threads=args.threads)
    outpath = os.path.join(output_dir, "tombo_output")
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    if args.seq_summary is not None:
        print("pycoQC")
        html = os.path.join(output_dir, name + ".html")
        run_qc(args.seq_summary, out_bam, html)

    print("Split Fast5s")
    split_fast5s_path = os.path.join(outpath, "split_fast5s")
    if not os.path.exists(split_fast5s_path):
        os.mkdir(split_fast5s_path)
    split_fast5s(args.fast5, split_fast5s_path, threads=args.threads)

    print("Run Tombo")
    check_call(f"tombo preprocess annotate_raw_with_fastqs --fast5-basedir {split_fast5s_path} "
               f"--fastq-filenames {fastq} --processes {args.threads} --overwrite".split())

    check_call(f"tombo resquiggle {split_fast5s_path} {args.reference} --processes {args.threads} "
               f"--num-most-common-errors 5 --rna --overwrite".split())

    check_call(f"tombo detect_modifications level_sample_compare --fast5-basedirs {args.wt_split_fast5s}"
               f"--alternate-fast5-basedirs {split_fast5s_path}"
               f"--statistics-file-basename WT_vs_{name}.level_compare_sample --processes {args.threads}".split())


if __name__ == '__main__':
    ret, time = (time_it(main))
    print(f"Running Time: {time} s")
