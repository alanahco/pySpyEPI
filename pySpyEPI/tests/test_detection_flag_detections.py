"""Test for functions in 'detection.flag_detections.py'
"""

import numpy as np
import unittest

from pySpyEPI.detection import flag_detections as fdet

class TestDetectionFlags(unittest.TestCase):
    """Tests for detection flags."""

    def setUp(self):
        """ Set up the test runs."""
        self.mlat_50 = np.array([-41, -10, 10, 41])

        self.mlat_51 = np.array([-41, -36, -30, -10, 10, 30, 36, 41])

        self.mlat_52 = np.array([-41, -30, -10, 10, 30, 41])

        self.mlat_53 = np.array([-41, -36, -10, 10, 36, 41])

        self.mlat_54 = np.array([-41, -36, -30, -10, 10, 30, 36, 41])

    def tearDown(self):
        pass

    def test_mlat_flag_low_one(self):
        """Test mlat flags when low lats have 1 flag."""
        # 5.0
        flag = np.array([1, 1, 1, 1])
        o1 = np.array([5.0, 0, 0, 5.0])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_50, flag))

        self.assertEqual(result, True)

        # 5.1
        flag = np.array([1] * 8)
        o1 = np.array([5.1, 4.1, 3, 0, 0, 3, 4.1, 5.1])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_51, flag))

        self.assertEqual(result, True)

        # 5.2
        flag = np.array([1] * 6)
        o1 = np.array([5.2, 3, 0, 0, 3, 5.2])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_52, flag))

        self.assertEqual(result, True)

        # 5.3
        flag = np.array([1] * 6)
        o1 = np.array([5.3, 4.0, 0, 0, 4.0, 5.3])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_53, flag))

        self.assertEqual(result, True)

        # 5.4
        flag = np.array([1, 9, 9, 1, 1, 9, 9, 1])
        o1 = np.array([5.4, 4.2, 3, 0, 0, 3, 4.2, 5.4])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_54, flag))

        self.assertEqual(result, True)


    def test_mlat_flag_low_none(self):
        """Test mlat flags when low lat has none."""
        # 10.0
        mlat = self.mlat_50[abs(self.mlat_50) > 25]
        flag = np.array([1, 1])
        o1 = np.array([10.0, 10.0])

        result = np.all(o1 == fdet.mlat_flag(mlat, flag))

        self.assertEqual(result, True)

        # 10.1
        mlat = self.mlat_51[abs(self.mlat_51) > 25]
        flag = np.array([1] * 6)
        o1 = np.array([10.1, 8.1, 6, 6, 8.1, 10.1])

        result = np.all(o1 == fdet.mlat_flag(mlat, flag))

        self.assertEqual(result, True)

         # 10.2
        mlat = self.mlat_52[abs(self.mlat_52) > 25]
        flag = np.array([1] * 4)
        o1 = np.array([10.2, 6, 6, 10.2])

        result = np.all(o1 == fdet.mlat_flag(mlat, flag))

        self.assertEqual(result, True)

          # 10.3
        mlat = self.mlat_53[abs(self.mlat_53) > 25]
        flag = np.array([1] * 4)
        o1 = np.array([10.3, 8.0, 8.0, 10.3])

        result = np.all(o1 == fdet.mlat_flag(mlat, flag))

        self.assertEqual(result, True)

        # 10.4
        mlat = self.mlat_54[abs(self.mlat_54) > 25]
        flag = np.array([1, 9, 9, 9, 9, 1])
        o1 = np.array([10.4, 8.2, 6, 6, 8.2, 10.4])

        result = np.all(o1 == fdet.mlat_flag(mlat, flag))

        self.assertEqual(result, True)

    def test_mlat_flag_low_nonone(self):
        """Test mlat flags when low lat is non-1 flag."""
        # 15.0
        flag = np.array([1, 9, 9, 1])
        o1 = np.array([15.0, 0, 0, 15.0])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_50, flag))

        self.assertEqual(result, True)

        # 15.1
        flag = np.array([1, 1, 1, 9, 9, 1, 1, 1])
        o1 = np.array([15.1, 12.1, 9, 0, 0, 9, 12.1,15.1])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_51, flag))

        self.assertEqual(result, True)

        # 15.2
        flag = np.array([1, 1, 9, 9, 1, 1])
        o1 = np.array([15.2, 9, 0, 0, 9, 15.2])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_52, flag))

        self.assertEqual(result, True)

        # 15.3
        flag = np.array([1, 1, 9, 9, 1, 1])
        o1 = np.array([15.3, 12.0, 0, 0, 12.0, 15.3])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_53, flag))

        self.assertEqual(result, True)

        # 15.4
        flag = np.array([1, 9, 9, 9, 9, 9, 9, 1])
        o1 = np.array([15.4, 12.2, 9, 0, 0, 9, 12.2, 15.4])

        result = np.all(o1 == fdet.mlat_flag(self.mlat_54, flag))

        self.assertEqual(result, True)

    def test_mlat_flag_mix(self):
        """Test mlat flags for all low latitude 1 flags."""

        i1 = np.array([-42, -30, -25, -20, -10, 5, 15, 35])
        i2 = np.array([1] * 8)
        o1 = np.array([5.2, 3, 0, 0, 0, 0, 0, 4.0])

        result = np.all(o1 == fdet.mlat_flag(i1, i2))

        self.assertEqual(result, True)

        i1 = np.array([-42, -30, -26, 5, 15, 35])
        i2 = np.array([1] * 6)
        o1 = np.array([10.2, 6, 6, 0, 0, 4.0])

        result = np.all(o1 == fdet.mlat_flag(i1, i2))

        self.assertEqual(result, True)

        i1 = np.array([-42, -30, -26, -5])
        i2 = np.array([1] * 4)
        o1 = np.array([5.2, 3, 3, 0])

        result = np.all(o1 == fdet.mlat_flag(i1, i2))

        self.assertEqual(result, True)

        i1 = np.array([5, 26, 32, 40])
        i2 = np.array([1] * 4)
        o1 = np.array([0, 3, 3, 5.2])

        result = np.all(o1 == fdet.mlat_flag(i1, i2))

        self.assertEqual(result, True)










