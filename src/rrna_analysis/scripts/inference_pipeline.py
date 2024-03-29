#!/usr/bin/env python
"""Run end to end signalAlign pipeline for yeast rRNA analysis"""
########################################################################
# File: inference_pipeline.py
#  executable: inference_pipeline.py
#
# Author: Andrew Bailey
# History: 04/07/21 Created
########################################################################

import os
from argparse import ArgumentParser
from subprocess import check_call, Popen, check_output, PIPE
import importlib.resources as resources
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
    parser.add_argument('--path_to_bin', action='store',
                        dest='path_to_bin', required=True, type=str, default=None,
                        help="Path to signalalign bin folder")
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
    parser.add_argument('--embed', action='store_true',
                        dest='embed', required=False,
                        help="Run MEA during signalAlign and embed into fast5s")
    parser.add_argument('--debug', action='store_true',
                        dest='debug', required=False,
                        help="Run signalalign in debug mode")

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


run_config = {"signal_alignment_args": {
    "target_regions": None,
    "track_memory_usage": False,
    "threshold": 0.1,
    "event_table": None,
    "embed": False,
    "delete_tmp": True,
    "output_format": "full"},
    "samples": [
        {
            "positions_file": None,
            "fast5_dirs": ["."],
            "bwa_reference": None,
            "fofns": [],
            "readdb": None,
            "fw_reference": None,
            "bw_reference": None,
            "kmers_from_reference": False,
            "motifs": None,
            "name": None,
            "probability_threshold": 0.7,
            "number_of_kmer_assignments": 10000,
            "alignment_file": False,
            "recursive": False,
            "assignments_dir": None
        }
    ],
    "path_to_bin": None,
    "complement_hdp_model": None,
    "template_hdp_model": None,
    "complement_hmm_model": None,
    "template_hmm_model": None,
    "job_count": None,
    "debug": False,
    "two_d": False,
    "output_dir": None,
    "constraint_trim": None,
    "diagonal_expansion": None,
    "traceBackDiagonals": 150,
    "filter_reads": 0,
    "perform_kmer_event_alignment": True,
    "overwrite": True,
    "rna": True,
    "ambig_model": None,
    "built_alignments": None,
    "delete_alignments": False
}


def create_config(outpath, bam, name, path_to_bin, readdb, fast5_dir, embed=False, threads=1, debug=False):
    template_hmm_model = "yeast_rrna_ivt_wt_trained_071521.model"
    ambig_model = "small_variants.model"
    variants = "yeast_18S_25S_variants.positions"
    reference = "yeast_25S_18S.fa"
    ref_index = "yeast_25S_18S.fa.fai"
    for x in [template_hmm_model, ambig_model, variants, reference, ref_index]:
        with resources.path("rrna_analysis.data", x) as p:
            shutil.copy(p, os.path.join(outpath, x))
            assert os.path.exists(os.path.join(outpath, x)), f"{os.path.join(outpath, x)} doesnt exist"

    run_config["samples"][0]["positions_file"] = os.path.join(outpath, variants)
    run_config["samples"][0]["bwa_reference"] = os.path.join(outpath, reference)
    run_config["samples"][0]["name"] = name
    run_config["samples"][0]["alignment_file"] = bam
    run_config["samples"][0]["readdb"] = readdb
    run_config["samples"][0]["fast5_dirs"] = [fast5_dir]

    run_config["signal_alignment_args"]["embed"] = embed

    run_config["ambig_model"] = os.path.join(outpath, ambig_model)
    run_config["template_hmm_model"] = os.path.join(outpath, template_hmm_model)
    run_config["output_dir"] = outpath
    run_config["job_count"] = threads
    run_config["path_to_bin"] = path_to_bin
    run_config["debug"] = debug

    return run_config


def split_fast5s(fast5_dir, output_dir, threads=1):
    temp_dir = os.path.join(output_dir, "tmp")
    check_call(f"compress_fast5 -t {threads} "
               f"-i {fast5_dir} -s {temp_dir} -c gzip".split())
    check_call(f"multi_to_single_fast5 --t {threads} "
               f"--input_path {temp_dir} --save_path {output_dir}".split())
    shutil.rmtree(temp_dir)


def index_reads(directory, fastq):
    check_call(f"embed_main index --directory {directory} {fastq}".split())
    readdb_file = fastq + ".index.readdb"
    assert os.path.exists(readdb_file), f"{readdb_file} doesnt exist"
    return readdb_file


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
        output_file_path = os.path.join(output_dir, args.name+".fastq")
        assert len(fastqs) >= 1, f"No fastqs found in {args.fastq}"
        fastq = concat_fastq_files(fastqs, output_file_path)
    else:
        fastq = args.fastq
    out_bam, filtered_sorted_bam = align_and_filter(fastq, args.reference, threads=args.threads)
    outpath = os.path.join(output_dir, "signalalign_output")
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    if args.seq_summary is not None:
        print("pycoQC")
        html = os.path.join(output_dir, name + ".html")
        run_qc(args.seq_summary, out_bam, html)

    print("Split and Index Fast5s")
    split_fast5s_path = os.path.join(outpath, "split_fast5s")
    if not os.path.exists(split_fast5s_path):
        os.mkdir(split_fast5s_path)
    split_fast5s(args.fast5, split_fast5s_path, threads=args.threads)
    readdb_path = index_reads(split_fast5s_path, fastq)

    print("Running SignalAlign")
    run_config_dict = create_config(outpath, filtered_sorted_bam, name, args.path_to_bin, readdb_path,
                                    split_fast5s_path, args.embed,
                                    threads=args.threads, debug=args.debug)
    config_path = os.path.join(outpath, "sa_run_config.json")
    save_json(run_config_dict, config_path)

    check_call(f"runSignalAlign.py run --config {config_path}".split())
    variant_call_path = os.path.join(outpath, "variant_calls")
    if not os.path.exists(variant_call_path):
        os.mkdir(variant_call_path)
    print("Running sa2bed")
    check_call(f"embed_main sa2bed -d {outpath}/tempFiles_alignment/{name}/ -a {run_config_dict['ambig_model']} "
               f"-o {variant_call_path}/{name}.bed -t {args.threads} -c B --overwrite --rna".split())


if __name__ == '__main__':
    ret, time = (time_it(main))
    print(f"Running Time: {time} s")
