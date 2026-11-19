"""Test for functions in 'detection.midlat_trough.py'
"""

import numpy as np
import unittest

from pySpyEPI.detection import midlat_trough as mdt

class TestZeroCrossings(unittest.TestCase):
    """Tests test_candidates."""

    def setUp(self):
        """ Set up the test runs."""
        pass

    def tearDown(self):
        pass

    def test_zero_cross(self):
        """Test zero_crossings."""

        lat = np.array([35, 40, 45])
        arr = np.array([10, 0, -10])
        zlats = np.array([40])
        slope = np.array([-1])
        cross = np.array([1])
        zero_lat, m_sign, icross = mdt.zero_crossings(arr, lat)

        result = np.all(zlats == zero_lat)
        self.assertTrue(result)

        result = np.all(slope == m_sign)
        self.assertTrue(result)

        result = np.all(cross == icross)
        self.assertTrue(result)