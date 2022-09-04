#!/usr/bin/env python
from setuptools import setup, find_packages  # , Extension


setup(name="rrna_scripts",
      version="0.0.1",
      description="Scripts and helper functions for data analysis of yeast rrna signalAlign sequencing data",
      author="Andrew Bailey",
      author_email="andbaile@ucsc.edu",
      url="https://github.com/adbailey4/rrna_scripts",
      package_dir={"": "src"},
      # ext_modules=extensions,
      packages=find_packages("src"),
      scripts=["src/rrna_analysis/scripts/re_run_signalalign.py",
               "src/rrna_analysis/scripts/plot_per_position_kmer_distributions.py",
               "src/rrna_analysis/scripts/run_embed_plot_wrapper.py",
               "src/rrna_analysis/scripts/train_test_accuracy_wrapper.py",
               "src/rrna_analysis/scripts/inference_pipeline.py",
               "src/rrna_analysis/scripts/tombo_pipeline.py"],
      install_requires=["numpy>=1.9.2",
                        "pysam>=0.8.2.1,<0.16.0",
                        "pandas>=0.23.1",
                        "matplotlib>=2.0.2",
                        "signalAlign>=0.3.0",
                        "py3helpers[seq_tools]>=0.5.0",
                        'embed>=0.0.5'],
      package_data={'': ['data/*']},
      )
