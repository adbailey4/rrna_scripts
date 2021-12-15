# rrna_scripts
Scripts to run analysis of yeast rRNA data 


#### Inference pipeline

* If the installation was performed correctly you should have `inference_pipeline.py` in your path.
```
usage: inference_pipeline.py [-h] --fastq FASTQ [--fast5 FAST5] --reference
                             REFERENCE --path_to_bin PATH_TO_BIN
                             [--output_dir OUTPUT_DIR] [--threads THREADS]
                             [--seq_summary SEQ_SUMMARY] --name NAME

Run end to end signalAlign pipeline for yeast rRNA analysis

optional arguments:
  -h, --help            show this help message and exit
  --fastq FASTQ         Path to fastq
  --fast5 FAST5         Path to fast5 directory
  --reference REFERENCE, -r REFERENCE
                        Path to reference fasta
  --path_to_bin PATH_TO_BIN
                        Path to signalalign bin folder
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Path to output_dir
  --threads THREADS, -t THREADS
                        Number of threads
  --seq_summary SEQ_SUMMARY
                        Path to sequence summary file
  --name NAME           Name of experiment
  --embed               Run MEA during signalAlign and embed into fast5s

```