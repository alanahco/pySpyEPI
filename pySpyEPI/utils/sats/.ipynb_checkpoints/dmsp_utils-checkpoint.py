"""DMSP data utility functions.
Updated 21 February 2026
created by grkwon
updated by alanahco

List Functions
--------------
rem_flags
find_intervals
rot_sat2field
join_magvel
"""

# created by Grace Kwon
import numpy as np
import pandas as pd
from pySpyEPI.utils import coords


def rem_flags(target_data, flag_data):
    """Processes data to replace flagged values with NaN.

    Parameters
    ----------
    target_data: array-like
        Given physical parameter array
    flag_data: array-like
        Flags for given parameter array

    Returns
    -------
    flagged: array-like
        The modified data with flagged indices replaced by NaN.

    Notes
    -----
    Created by Grace Kwon (grkwon@umich.edu)
    """

    # Do some reassigning of values (for F17)
    flag_data[flag_data == 6] = 1
    flag_data[flag_data == 7] = 2
    flag_data[flag_data == 8] = 3

    # Identify the indices where flag values are 3 or greater
    flagged_indices = flag_data >= 3

    # Copy the target data
    # Ensure float type to support NaN
    flagged = np.array(target_data, dtype=float)

    # Replace the flagged locations with NaN
    flagged[flagged_indices] = np.nan

    return flagged


def find_intervals(mlat, val=50):
    """Locates data within a certain latitude range.

    Parameters
    ----------
    mlat : array-like
        Magnetic latitude array

    Returns
    -------
    intervals : list-like
        of [start, end] of intervals included within the latitude range

    Notes
    -----
    Created by Grace Kwon (grkwon@umich.edu)
    """

    # Pick mlat ranges to include
    mask = (mlat >= -val) & (mlat <= val)  # For equatorial region

    # Find where satellite moves in/out of range
    # +1 to shift to the right indices
    change_indices = np.where(np.diff(mask))[0] + 1

    intervals = []
    start_index = None

    # Handle the first interval if it starts with mlat > 40
    if mask[0]:
        start_index = 0

    for idx in change_indices:
        if mask[idx]:  # Start of a new interval
            if start_index is None:
                start_index = idx
        else:  # End of an interval
            if start_index is not None:
                intervals.append([start_index, idx - 1])
                start_index = None

    # If the last interval ends with True
    # add the final interval (ensure it's valid)
    if start_index is not None and start_index != len(mlat) - 1:
        intervals.append([start_index, len(mlat) - 1])

    return intervals


def rot_sat2field(time, gdlat, glon, gdalt, alongV, crossV, vertV, alongB,
                  crossB, vertB):
    """Locates data within a certain latitude range.

    Parameters
    -----------
    time : array-like
        satellite time
    gdlat, glon, gdalt : array-like
        satellite location in geodetic coords
    alongV, crossV, vertV : array-like
        ion velocity in satellite coords
    alongB, crossB, vertB : array-like
        magnetic field in satellite coords

    Returns
    -------
    ionVelMFA[:,0], ionVelMFA[:,1], ionVelMFA[:,2]:
        ion velocity in field-aligned coordinates (parallel, zonal, meridional)

    Notes
    -----
    Created by Grace Kwon (grkwon@umich.edu)
    """

    satN, satE, satC = coords.ecef2ned(time, gdlat, glon, gdalt)

    # Rotation matrix--sat to NEC
    rotSat2NEC = coords.track_to_nec_rotation(satN, satE, satC, C_is_up=False)

    # Ion velocity from sat to NEC
    ionVel_NEC = coords.coordinate_transform(alongV, crossV, vertV, rotSat2NEC)

    # Magnetic field from sat to NEC
    mag_NEC = coords.coordinate_transform(alongB, crossB, vertB, rotSat2NEC)

    # Rotation matrix-- NEC to field-aligned
    mag_MFA, rotNEC2MFA = coords.NEC_2_MFA_Rotation(mag_NEC)

    # Ion velocity from NEC to field-aligned
    ionVelMFA = coords.coordinate_transform(ionVel_NEC[:, 0], ionVel_NEC[:, 1],
                                            ionVel_NEC[:, 2], rotNEC2MFA)

    return ionVelMFA[:, 0], ionVelMFA[:, 1], ionVelMFA[:, 2]


def join_magvel(df):
    """Joins velocity data and magnetic field data (slightly different times).

    Parameters
    ----------
    df: DataFrame
        velocity data
    Returns
    -------
    dfFull : DataFrame
        dataframe containing time, mlt, mlat, lat, lon, alt, ni, te, ti,
        alongV, crossV, and vertV

    Notes
    -----
    Created by Grace Kwon (grkwon@umich.edu)
    """
    dfFull = df
    time = dfFull.time.to_numpy()
    mlt = dfFull.mlt.to_numpy()
    mlat = dfFull.mlat.to_numpy()
    lat = dfFull.lat.to_numpy()
    lon = dfFull.lon.to_numpy()
    alt = dfFull.alt.to_numpy()
    ni = dfFull.ni.to_numpy()
    ti = dfFull.ti.to_numpy()
    te = dfFull.te.to_numpy()
    alongV = dfFull.alongV.to_numpy()
    crossV = dfFull.crossV.to_numpy()
    vertV = dfFull.vertV.to_numpy()

    dfFull = pd.DataFrame({
        "time": time,
        "mlt": mlt,
        "mlat": mlat,
        "lat": lat,
        "lon": lon,
        "alt": alt,
        "ni": ni,
        "ti": ti,
        "te": te,
        "alongV": alongV,
        "crossV": crossV,
        "vertV": vertV,
    })

    return dfFull
