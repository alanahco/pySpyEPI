#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
# Original Author: Angeline G. Burrell, NRL (September 2026)
# -----------------------------------------------------------------------------
"""Tests for functions in `pySpyEPI.barrel`."""

import numpy as np
import unittest

from pySpyEPI.barrel import barrel_roll as br


class TestBarrelRollFuncs(unittest.TestCase):
    """Tests for fitlers functions."""

    def setUp(self):
        """Set up the test runs."""
        self.in_arr = np.array([0, 6, 6, 7, 8, 50, np.nan])
        self.outx = None
        self.outy = None
        return

    def tearDown(self):
        """Tear down the test environment."""
        del self.in_arr, self.outx, self.outy
        return

    def test_simple_barrel_forward(self):
        # Set a second array of the same shape and range, but different values
        xvar = np.linspace(np.nanmin(self.in_arr), np.nanmax(self.in_arr),
                           self.in_arr.shape[0])


        self.outx, self.outy = br.simple_barrel(
            xvar, self.in_arr, 30, direction='forward')

        # Test the output
        self.assertTrue(np.isfinite(self.outy).all())
        self.assertLessEqual(self.outy.max(), np.nanmax(self.in_arr))
        self.assertGreaterEqual(self.outy.min(),
                                np.nanmin(self.in_arr))
        return

    def test_simple_barrel_no_roi_forward(self):
        """Test success of `simple_barrel` with a small barrel."""
        # Set a second array of the same shape and range, but different values
        xvar = np.linspace(np.nanmin(self.in_arr), np.nanmax(self.in_arr),
                           self.in_arr.shape[0])

        self.outx, self.outy = br.simple_barrel(
            xvar, self.in_arr, 1, direction='forward')

        # Test the output
        self.assertTrue(np.isfinite(self.outy).all())
        self.assertLessEqual(self.outy.max(), np.nanmax(self.in_arr))
        self.assertGreaterEqual(self.outy.min(),
                                np.nanmin(self.in_arr))
        return

    def test_simple_barrel_backward(self):
        """Test success of `simple_barrel` with a backward barrel."""
        # Set a second array of the same shape and range, but different values
        xvar = np.linspace(np.nanmin(self.in_arr), np.nanmax(self.in_arr),
                           self.in_arr.shape[0])


        self.outx, self.outy = br.simple_barrel(
            xvar, self.in_arr, 30, direction='backward')

        # Test the output
        self.assertLessEqual(np.nanmax(self.outy), np.nanmax(self.in_arr))
        self.assertGreaterEqual(np.nanmin(self.outy),
                                np.nanmin(self.in_arr))
        self.assertGreater(self.outx[0], self.in_arr[0])
        return
