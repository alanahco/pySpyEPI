"""Detect Midlatitude Trough
Adapted from Liu and Xiong et al 2020
Updated 21 February 2026

List Functions
--------------
zero_crossings
find_trough_lat
"""
import pandas as pd
import numpy as np
from pySpyEPI.utils.calc import moving_average_same_length


def zero_crossings(arr, lat):
    """Find indices where the array crosses zero (sign change).
    Parameters
    ----------
    arr : array-like
        array to be analyzed for zero crossings

    Returns
    -------
    zero_lat : array-like
        latitudes at zero crossing found using linear fit
    m_sign : array-like
        signs of slope of linear fit used for zero crossing lats
    icross : array-like
        Indices just before the crossing
    """
    signs = np.sign(arr)
    icross = np.where(np.diff(signs))[0]

    # make sure that it is not the first or last value in the array
    icross = icross[icross != 0]
    icross = icross[icross != (len(arr) - 1)]

    # Get the locations of sign changes
    # ichange = np.where(arr[1:] != arr[:-1])[0]

    # Use a linear fit to estimate the latitude of the sign change
    zero_lat = list()

    # record if the slope is negative or positive
    m_sign = list()
    for cind in icross:
        ind_bf = cind - 1
        ind_af = cind + 1
        slope = (arr[ind_bf] - arr[ind_af]) / (lat[ind_bf] - lat[ind_af])
        intercept = arr[cind] - slope * lat[cind]
        zero_lat.append(-intercept / slope)
        m_sign.append(np.sign(slope))

    zero_lat = np.array(zero_lat)
    m_sign = np.array(m_sign)

    return zero_lat, m_sign, icross


def find_trough_lat(data, equator_bound=35, auroral_bound=75, set_lat=55,
                    ne_str='ne', mlat_str='mlat', freq=1):
    """Find latitude of trough for the northern and southern hemisphere.

    Parameters
    ----------
    data : DataFrame
        IRR data inclduing magnetic latitude and Ne
    equator_bound : float kwarg (35)
        equatorward magnetic latitude to search down to
        set as a positive number, will be used as +/-equator_bound
    auroral_bound : float kwarg (75)
        autroalward magnetic latitude to search up to
        set as a positive number, will be used as +/-auroral_bound
    set_lat : float kwarg (55)
        if a trough cannot be found, set_lat will be used in its place
        set as a positive number, will be used as +/-set_lat
    ne_str, mlat_str : string
        default 'ne' and 'mlat'
        density and mlat column names
    freq : int
        frequency of data
        defualt 1 Hz for IRR
    Returns
    -------
    lat_retrun : array of length 2
        northern and southern hemispehre midlatitude trough location
    """

    # northern Hemisphere
    data_north = data[((data[mlat_str] >= equator_bound)
                       & (data[mlat_str] <= auroral_bound))]

    # southern hemsiphere
    data_south = data[((data[mlat_str] >= -1 * auroral_bound)
                       & (data[mlat_str] <= -1 * equator_bound))]

    # Electron density to be used to determine trough location
    # iterate for norhtern and southern hemisphere
    lat_return = []
    for i in range(2):
        if i == 0:

            # Northern
            data_lat_filt = data_north.copy()

        else:
            # Southern
            data_lat_filt = data_south.copy()

        # iterate through the length of the data
        if len(data_lat_filt) > 0:

            # use a large moving average (320s) to get the overall trend
            window_size = 320 * freq
            filt320 = moving_average_same_length(data_lat_filt[ne_str].values,
                                                 window_size)

            # Small moving average (10 s) to remove minor variations)
            filt10 = moving_average_same_length(data_lat_filt[ne_str].values,
                                                10 * freq)

            # Detrend the data 10s filter - 320 second filter
            detrended = filt10 - filt320

            # Get the zero crossings mlats using zero_crossings function
            zs, ms, zi = zero_crossings(detrended,
                                        data_lat_filt[mlat_str].values)

            # Create a zero crossing DataFrame
            zero_cross_df = pd.DataFrame({"mlats": zs})

            zero_cross_df["m_signs"] = ms
            zero_cross_df["idx"] = zi

            # set the trough to nans
            zero_cross_df["trough"] = np.nan

            # find crossings where slope is negative
            z_first = zero_cross_df[zero_cross_df["m_signs"] < 0]

            # find crossings where slope is positive
            z_second = zero_cross_df[zero_cross_df["m_signs"] > 0]

            # get the indices of z_first and z_second
            z2_index = z_second.index

            # iterate through the crossings assuming there is more than 1
            if len(zero_cross_df["mlats"] > 1):
                for r in range(len(z_first)):
                    # get index of first crossing
                    z_f = z_first.index[r]

                    # get second index using first index plus 1
                    z_s = z2_index[z_f + 1 == z2_index]

                    # check that second index is present
                    if len(z_s) > 0:

                        # Get Ne min between zero crossings
                        z1_ml = zero_cross_df["mlats"].loc[z_f]
                        z2_ml = zero_cross_df["mlats"].loc[z_s[0]]

                        # mask data using magnetic latitudes of each crossing
                        if z1_ml < z2_ml:
                            mask = ((data_lat_filt[mlat_str] >= z1_ml)
                                    & (data_lat_filt[mlat_str] <= z2_ml))
                        else:
                            mask = ((data_lat_filt[mlat_str] >= z2_ml)
                                    & (data_lat_filt[mlat_str] <= z1_ml))

                        # get minimum Ne in the Zero crossing width
                        # first check if masked array is longer than 0
                        # if it is then use that to get zmin
                        # if not try using the indices (2 values)
                        if len(data_lat_filt[ne_str][mask]) != 0:
                            zmin = min(data_lat_filt[ne_str][mask])
                        else:
                            z1_id = zero_cross_df["idx"].loc[z_f]
                            z2_id = zero_cross_df["idx"].loc[z_s[0]]
                            z1 = data_lat_filt[ne_str].iloc[z1_id]
                            z2 = data_lat_filt[ne_str].iloc[z2_id]
                            z_arr = np.array([z1, z2])
                            zmin = min(z_arr)

                        # get closest ne values
                        z1_id = zero_cross_df["idx"].loc[z_f]
                        z2_id = zero_cross_df["idx"].loc[z_s[0]]
                        z1 = data_lat_filt[ne_str].iloc[z1_id]
                        z2 = data_lat_filt[ne_str].iloc[z2_id]

                        # calculate the depth of both sides returns a percent
                        dr1 = (z1 - zmin) / z1 * 100
                        dr2 = (z2 - zmin) / z2 * 100

                        # if either side is > 40%, then set that as the trough
                        if (dr1 >= 40) | (dr2 >= 40):
                            zfr = z_first.index[r]
                            zero_cross_df.loc[zfr, "trough"] = 1
                            zero_cross_df.loc[zfr + 1, "trough"] = 1

                        else:
                            zfr = z_first.index[r]
                            zero_cross_df.loc[zfr, "trough"] = 0
                            zero_cross_df.loc[zfr + 1, "trough"] = 0

                # is there a candidate for trough location?
                if len(zero_cross_df["mlats"][zero_cross_df["trough"]
                       == 1]) == 0:

                    # If no candidates, set trogh as set_lat
                    if i == 0:
                        tr_lat = set_lat
                    else:
                        tr_lat = -1 * set_lat

                else:
                    # if there are candidates, set the trough as the
                    # most equatorward location
                    tr_lat = min(
                        zero_cross_df["mlats"][zero_cross_df["trough"] == 1],
                        key=abs)

            else:

                # if there are no zero crossings, set trough lat to set-lat
                if i == 0:
                    tr_lat = set_lat
                else:
                    tr_lat = -1 * set_lat

            # get sign of trough lat
            s0 = np.sign(tr_lat)

            # if tr_lat is greater than set_lat, use set_lat * sign of tr_lat
            if abs(tr_lat) > set_lat:
                tr_lat = set_lat * s0

            lat_return.append(tr_lat)

        else:

            # if there is not enough data, set tr_lat to set_lat
            if i == 0:
                tr_lat = set_lat
            else:
                tr_lat = -1 * set_lat
            lat_return.append(tr_lat)

    if np.sign(lat_return[0]) == np.sign(lat_return[1]):
        lat_return = [set_lat, -1 * set_lat]

    return lat_return
