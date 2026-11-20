"""Test for functions in 'barrel.barrel_utils.py'
"""

import numpy as np
import unittest

from pySpyEPI.barrel import barrel_utils as bu


class TestBarrelUtils(unittest.TestCase):
    """Tests for barrel_utils."""

    def setUp(self):
        """ Set up the test runs."""
        pass

    def tearDown(self):
        pass

    def test_find_all_gaps_present(self):
        """Test find_all_gaps with gaps present."""
        arr = np.array([1, 2, 3, 4, 5, 8, 9, 14, 16, 20, 21, 22, 24])

        out = np.array([0, 5, 7, 8, 9, 12])

        result = np.all(out == np.array(bu.find_all_gaps(arr)))

        self.assertTrue(result)

    def test_find_all_gaps_absent(self):
        """Test find_all_gaps without gaps present."""
        arr = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

        out = np.array([0, 9])

        result = np.all(out == np.array(bu.find_all_gaps(arr)))

        self.assertTrue(result)


    def test_find_nan_ranges_some(self):
        """Test find_nan_ranges when some are nan."""
        arr = np.array([3, 4, np.nan, np.nan, 9, 14, 16, np.nan, np.nan, 22])
        out = np.array([(2, 4), (7, 9)])

        result = np.all(out == np.array(bu.find_nan_ranges(arr)))

        self.assertTrue(result)


    def test_find_nan_ranges_none(self):
        """Test find_nan_ranges when none are nan."""
        arr = np.array([3, 4, 11, 12, 20, 22])
        out = np.array([])

        result = np.all(out == np.array(bu.find_nan_ranges(arr)))

        self.assertTrue(result)

    def test_find_nan_ranges_all(self):
        """Test find_nan_ranges when all are nan."""
        arr = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
        out = np.array([(0, 5)])

        result = np.all(out == np.array(bu.find_nan_ranges(arr)))

        self.assertTrue(result)
