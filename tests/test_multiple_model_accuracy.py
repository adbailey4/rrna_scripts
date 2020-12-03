#!/usr/bin/env python
"""Tests for multiple_model_accuracy.py"""
########################################################################
# File: multiple_model_accuracy.py
#  executable: multiple_model_accuracy.py
#
# Author: Andrew Bailey
# History: Created 12/02/20
########################################################################
import unittest
import os

# from py3helpers.utils import merge_dicts, captured_output
# from signalalign.alignedsignal import *
# from signalalign.signalAlignment import SignalAlignment, create_signalAlignment_args


class MultipleModelAccuracy(unittest.TestCase):

    tmp_directory = None

    @classmethod
    def setUpClass(cls):
        super(MultipleModelAccuracy, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])
        cls.test_file = os.path.join(cls.HOME, "src/rrna_analysis/__init__.py")

    def test_test(self):
        self.assertTrue(os.path.exists(self.test_file))


if __name__ == "__main__":
    unittest.main()
