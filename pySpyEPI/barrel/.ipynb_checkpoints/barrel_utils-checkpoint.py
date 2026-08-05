"""Barrel Rolling utility functions
Created 19 February 2026
Updated 29 June 2026
alanahco
List Functions
--------------
find_nan_ranges
find_all_gaps
"""
import numpy as np


def find_nan_ranges(arr):
    """Identify continuous ranges of NaN values in an array.

    Parameters
    ----------
    arr : array-like
        array with nans
    Returns
    -------
    nan_list : array-like
        List of (start_idx, end_idx) for each continuous NaN section.
    """
    # Get continuous ranges of nan values
    isnan = np.isnan(arr)
    edges = np.diff(isnan.astype(int))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0] + 1

    if isnan[0]:
        starts = np.insert(starts, 0, 0)
    if isnan[-1]:
        ends = np.append(ends, len(arr))
    nan_list = list(zip(starts, ends))

    return nan_list


def find_all_gaps(arr):
    """Find gap indices.

    Parameters
    ----------
    arr : array-like
        array of indices
    Returns
    -------
    gap_indices : array-like
        indices of gap start and end and first and last index
    Notes
    -----
    e.g. in an array 2,3,5,6,7,8
    find_all_gaps will return index 2 indicating no gap from 0:2 and no gap
    includes 0 and last index, so return would be [0, 2, len(arr) - 1]
    """
    gap_indices = []
    # Iterate through the array and find where the gaps start
    for i in range(len(arr) - 1):
        if arr[i + 1] != arr[i] + 1:
            gap_indices.append(i + 1)  # Append the index where the gap starts

    if len(arr) != 0:
        gap_indices.insert(0, 0)
        gap_indices.insert(len(gap_indices), len(arr) - 1)
        gap_indices = np.unique(gap_indices)

    return gap_indices
