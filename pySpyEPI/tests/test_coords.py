#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
# Original Author: Angeline G. Burrell, NRL (September 2026)
# -----------------------------------------------------------------------------
"""Tests for functions in `utils.coords`."""

import datetime as dt
import numpy as np
import unittest

from pySpyEPI.utils import coords


class TestTimeFuncs(unittest.TestCase):
    """Tests for time-handling functions."""

    def setUp(self):
        """Set up the test runs."""
        self.dtime = dt.datetime(1999, 2, 11)
        self.lon = 10.0
        self.out = None
        return

    def tearDown(self):
        """Tear down the test environment."""
        del self.dtime, self.lon, self.out
        return

    def evaluate_offset(self):
        """Evaluate the offset between UT and local time."""
        # Get the time difference in seconds regardless of the timezone
        sec = (dt.datetime.strptime(
            self.out[0].strftime("%Y-%m-%d %H:%M:%S:%u"),
            "%Y-%m-%d %H:%M:%S:%u") - self.dtime).total_seconds()
        self.assertEqual(self.lon, sec / 240.0)
        return

    def test_longitude_to_local_time_list(self):
        """Test success for datetime casting with list inputs."""
        # Cycle through potential time formats
        for in_time in [dt.date(self.dtime.year, self.dtime.month,
                                self.dtime.day), self.dtime,
                        np.datetime64(self.dtime.strftime('%Y-%m-%d')),
                        dt.datetime(self.dtime.year, self.dtime.month,
                                    self.dtime.day, tzinfo=dt.timezone.utc)]:
            with self.subTest(in_time=in_time):
                # Convert the time
                self.out = coords.longitude_to_local_time([self.lon], [in_time])
                self.assertTupleEqual(self.out.shape, (1,))
                self.evaluate_offset()
        return

    def test_longitude_to_local_time_value(self):
        """Test success for datetime casting with single value inputs."""
        # Cycle through potential time formats
        for in_time in [dt.date(self.dtime.year, self.dtime.month,
                                self.dtime.day), self.dtime,
                        np.datetime64(self.dtime.strftime('%Y-%m-%d')),
                        dt.datetime(self.dtime.year, self.dtime.month,
                                    self.dtime.day, tzinfo=dt.timezone.utc)]:
            with self.subTest(in_time=in_time):
                # Convert the time
                self.out = [coords.longitude_to_local_time(self.lon, in_time)]
                self.evaluate_offset()
        return

    def test_longitude_to_local_time_array(self):
        """Test success for datetime casting with array inputs."""
        # Cycle through potential time formats
        for in_time in [dt.date(self.dtime.year, self.dtime.month,
                                self.dtime.day), self.dtime,
                        np.datetime64(self.dtime.strftime('%Y-%m-%d')),
                        dt.datetime(self.dtime.year, self.dtime.month,
                                    self.dtime.day, tzinfo=dt.timezone.utc)]:
            with self.subTest(in_time=in_time):
                # Convert the time
                self.out = coords.longitude_to_local_time(
                    np.array([self.lon]), np.array([in_time]))
                self.assertTupleEqual(self.out.shape, (1,))
                self.evaluate_offset()
        return

    def test_longitude_to_local_time_mult_time(self):
        """Test success for datetime casting with multiple times."""
        # Cycle through different unequal length combinations
        for lon_len, time_len in [[1, 5], [3, 1], [3, 5]]:
            lon_in = np.full(shape=(lon_len,), fill_value=self.lon)
            ut_in = np.full(shape=(time_len,), fill_value=self.dtime)
            with self.subTest(lon_len=lon_len, time_len=time_len):
                with self.assertRaisesRegex(
                        ValueError, 'cannot add indices of unequal length'):
                    coords.longitude_to_local_time(lon_in, ut_in)
        return

class TestLocFuncs(unittest.TestCase):
    """Tests for location-handling functions."""

    def setUp(self):
        """Set up the test runs."""
        self.dtime = dt.datetime(1999, 2, 11)
        self.lon = 10.0
        self.lat = 0.0
        self.ref_alt = 6378137.0
        self.out = None
        return

    def tearDown(self):
        """Tear down the test environment."""
        del self.dtime, self.lon, self.lat, self.out
        return

    def evaluate_mag_loc(self, mag_type):
        """Evaluate the magnetic location using the specified mag type.

        Parameters
        ----------
        mag_type : str
            ApexPy coordinate type string

        """
        # Define the expected outputs
        out_lat = {6378137.0: {'geo': self.lat, 'qd': -4.0133490562438965,
                               'apex': -45.15626907348633},
                   0.0: {'geo': self.lat, 'qd': -12.547586441040039,
                         'apex': -12.547586441040039}}
        out_lon = {6378137.0: {'geo': self.lon, 'qd': 81.84486389160156,
                               'apex': 81.84486389160156},
                   0.0: {'geo': self.lon, 'qd': 83.14351654052734,
                         'apex': 83.14351654052734}}

        # Ensure the values to check exist
        self.assertTrue(self.ref_alt in out_lat.keys())
        self.assertTrue(mag_type in out_lat[self.ref_alt].keys())

        # Ensure the values are the same
        self.assertAlmostEqual(self.out[0], out_lat[self.ref_alt][mag_type])
        self.assertAlmostEqual(self.out[1], out_lon[self.ref_alt][mag_type])
        return


    def test_earth_radius_list(self):
        """Test the default earth radius calculation with list inputs."""
        self.out = coords.earth_radius([self.lat, self.lat])
        self.assertTrue((self.out == self.ref_alt).all())
        self.assertTupleEqual(self.out.shape, (2,))
        return

    def test_earth_radius_array(self):
        """Test the default earth radius calculation with list inputs."""
        self.out = coords.earth_radius(np.full(shape=(3, 4),
                                               fill_value=self.lat))
        self.assertTrue((self.out == self.ref_alt).all())
        self.assertTupleEqual(self.out.shape, (3, 4))
        return

    def test_compute_magnetic_coords_mag_type(self):
        """Test the mag coord calculation for different coord systems."""
        # Define the heights
        heights = [self.ref_alt, 0.0]

        # Cycle through the valid combinations
        for mag_type in ['apex']:
            for height in heights:
                self.ref_alt = height
                with self.subTest(mag_type=mag_type, height=self.ref_alt):
                    self.out = coords.compute_magnetic_coords(
                        self.lat, self.lon, self.ref_alt / 1000.0, [self.dtime])

                    self.evaluate_mag_loc(mag_type)
        return
