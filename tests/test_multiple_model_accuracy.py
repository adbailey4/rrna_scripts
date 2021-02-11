#!/usr/bin/env python
"""Tests for kmer_pos_mapping.py"""
########################################################################
# File: kmer_pos_mapping.py
#  executable: kmer_pos_mapping.py
#
# Author: Andrew Bailey
# History: Created 12/03/20
########################################################################
import unittest
import numpy as np
from rrna_analysis.kmer_pos_mapping import *
from rrna_analysis.multiple_model_accuracy import *
from rrna_analysis.scripts.train_test_accuracy_wrapper import get_hdp_models


class TestKmerPosMapping(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestKmerPosMapping, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])
        cls.test_reference = os.path.join(cls.HOME, "tests/test_files/yeast_25S_18S.fa")
        cls.test_positions = os.path.join(cls.HOME, "tests/test_files/test.positions")
        cls.test_mod_file = os.path.join(cls.HOME, "tests/test_files/test_mod_file.csv")
        cls.mod_file = os.path.join(cls.HOME, "tests/test_files/mod_file.csv")
        cls.pos_file = os.path.join(cls.HOME, "tests/test_files/yeast_18S_25S_variants.positions")

        cls.kpm = KmerPosMapping(cls.test_reference, cls.pos_file, cls.mod_file)

        cls.hdp_dir = os.path.join(cls.HOME, "tests/test_files/hdp_models")

    # def test_sort_dir(self):
    #     dir_path = "/Volumes/gdrive/rrna_kube_testing/supervised/probability_sweep/" \
    #                "train_500_test_500_prob_0.5_em_iterations_30/testing_accuracy_csvs"
    #     print(sort_dir(dir_path))

    # def test_plot(self):
    #     key = "accuracy"
    #     dir_path = "/Volumes/gdrive/rrna_kube_testing/supervised/probability_sweep/" \
    #                "train_500_test_500_prob_0.5_em_iterations_30/testing_accuracy_csvs"
    #     model_dir = "/Volumes/gdrive/rrna_kube_testing/supervised/probability_sweep/" \
    #                 "train_500_test_500_prob_0.5_em_iterations_30/training_models"
    #     model_n = 30
    #     p, k = plot_accuracy_vs_delta_and_accuracy_over_time(self.kpm, dir_path, model_dir, model_n, high_percent=100,
    #                                                          low_percent=90, low_delta=6, high_delta=np.inf, key=key,
    #                                                          max_delta=False)

    def test_get_hdp_models(self):
        models = get_hdp_models(self.hdp_dir)
        for x in models:
            self.assertTrue(os.path.exists(x))
        self.assertEqual(len(models), 2)


if __name__ == "__main__":
    unittest.main()
