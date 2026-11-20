#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
# Original Author: Angeline G. Burrell, NRL (September 2026)
# -----------------------------------------------------------------------------
"""Tests for functions in `utils.calc`."""

import numpy as np
import unittest

from pySpyEPI.utils import calc


class TestCalc(unittest.TestCase):
    """Tests for fitlers functions."""

    def setUp(self):
        """Set up the test runs."""
        self.in_arr = np.array([0, 6, 6, 7, 8, 50, np.nan])
        self.out = None
        return

    def tearDown(self):
        """Tear down the test environment."""
        del self.in_arr, self.out
        return

    def test_detect_outliers(self):
        """Test success of detect outliers."""
        # Find outliers
        self.out = calc.detect_outliers(self.in_arr)

        # Test output, outliers must be numbers
        self.assertListEqual(list(self.in_arr[self.out]), [0, 50])
        return

    def test_detect_outliers_single_val(self):
        """Test success of detect outliers."""
        # Find outliers
        self.out = calc.detect_outliers(self.in_arr[0])

        # Test output, can't have outliers with a single value
        self.assertEqual(len(self.out), 0)
        return

    def test_detect_outliers_nan(self):
        """Test detect outliers routine with only NaN."""
        # Find outliers
        self.in_arr = np.full(shape=self.in_arr.shape, fill_value=np.nan)
        self.out = calc.detect_outliers(self.in_arr)

        # Test no outliers, as no real values
        self.assertEqual(len(self.out), 0)
        return

    def test_moving_median(self):
        """Test moving median routine."""
        self.out = calc.moving_median(self.in_arr, 3)

        out = np.array([ 3., 6., 6., 7., 8., 29., 50.])

        # Test
        self.assertTrue(np.isfinite(self.out).all())
        self.assertTrue(np.all(out == self.out))
        self.assertEqual(len(self.out), len(self.in_arr))

    def test_detect_moving_median_nan(self):
        """Test detect outliers routine with only NaN."""
        # Find outliers
        self.in_arr = np.full(shape=self.in_arr.shape, fill_value=np.nan)
        self.out = calc.moving_median(self.in_arr, 3)

        # Test no outliers, as no real values
        self.assertTrue(np.isnan(self.out).all())
        self.assertEqual(len(self.out), len(self.in_arr))
        return




