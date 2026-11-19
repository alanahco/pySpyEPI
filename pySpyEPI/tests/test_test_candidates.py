"""Test for functions in 'detection.test_candidates.py'
"""

import numpy as np
import unittest

from pySpyEPI.detection import test_candidates as tcan


class TestTestCandidates(unittest.TestCase):
    """Tests test_candidates."""

    def setUp(self):
        """ Set up the test runs."""
        pass

    def tearDown(self):
        pass

    def test_no_doubbubs_present(self):
        """Test no_doubbubs with doubles present."""
        peaks = np.array([7, 12, 13, 21, 23, 30])
        properties = {}
        properties["left_ips"] = np.array([5.2, 10.2, 11.1, 17.7, 20.8, 29.3])
        properties["right_ips"] = np.array([8.8, 14.3, 14.3, 24.9, 23.7, 32.1])

        out = np.array([0, 0, 1, 0, 1, 0])

        result = np.all(out == tcan.no_doubbubs(peaks, properties))
        self.assertTrue(result)

    def test_no_doubbubs_absent(self):
        """Test no_doubbubs without doubles present."""
        peaks = np.array([7, 12, 16, 22, 27, 30])
        properties = {}
        properties["left_ips"] = np.array([5.2, 10.2, 15.1, 20.5, 26.4, 29.3])
        properties["right_ips"] = np.array([8.8, 14.3, 16.9, 24.9, 28.1, 32.1])

        out = np.array([0, 0, 0, 0, 0, 0])

        result = np.all(out == tcan.no_doubbubs(peaks, properties))
        self.assertTrue(result)
