"""Barrel Rolling codes
Created 19 February 2026
Updated 29 June 2026
alanahco

List Functions
--------------
simple_barrel
full_barrel
get_barrel_info
triple_barrel
just_barrel
"""
import numpy as np
import pandas as pd
import math
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d
from pySpyEPI.barrel import barrel_utils


def simple_barrel(x_sc, ne_sc, r_sc, direction='forward'):
    """Roll barrel over array to get detrended tec/ne (1 direciton, 1 radius).

    Parameters
    ----------
    x_sc: array-like
        scaled x array
        could be latitude or time
    ne_sc: array-like
        scaled electron density array
    r_sc : double
        radius of barrel (same unit as x_sc)
    direction : str kwarg
        forward (default) and backward depending on rolling direction

    Returns
    -------
    contact_xs : array-like
        array of x data points that were touched by barrel
    contact_ys : array-like
        array of y data points that were touched by barrel
    """

    if direction == 'forward':
        # Keep the original direction
        ne = ne_sc
        lat = x_sc
    elif direction == 'backward':
        # Reverse the x and y values
        ne = ne_sc[::-1]
        lat = x_sc[::-1]

    # starting x and y values for barrel rolling
    strt_con_y = ne[0]
    strt_con_x = lat[0]

    # empty array of contact points
    f_con_xs = []
    f_con_ys = []

    # set j to 0 to keep track of index of array of next contact point
    j = 0

    # while index is less than last point in the array
    # keep lookin for contact points
    while (j < len(ne) - 1):
        if j < len(ne) - 1:

            # Forward Rolling only
            f_con_xs .append(strt_con_x)
            f_con_ys.append(strt_con_y)

            # get the regions of interest (within barrel view) FORWARD ROLLING
            # > than start time and less than + diameter (same units as x)

            if direction == 'forward':
                x_roi = lat[((lat > strt_con_x)
                             & (lat <= strt_con_x + 2 * r_sc))]
                y_roi = ne[(lat > strt_con_x) & (lat <= strt_con_x + 2 * r_sc)]
            elif direction == 'backward':
                x_roi = lat[((lat < strt_con_x)
                             & (lat >= strt_con_x - 2 * r_sc))]
                y_roi = ne[(lat < strt_con_x) & (lat >= strt_con_x - 2 * r_sc)]

            # Calcualte angular distance delta for each delta = beta - theta
            deltas = []

            # iterate through x region of interest
            for i in range(len(x_roi)):

                # Calculate change in x and y
                if direction == 'forward':
                    del_x = x_roi[i] - strt_con_x
                    del_y = y_roi[i] - strt_con_y
                elif direction == 'backward':
                    # For backward, change direction of del_x to keep del_x
                    # positive
                    # if del_y is flipped, then we have the same vector flipped
                    # by 180 degrees
                    del_x = strt_con_x - x_roi[i]
                    del_y = y_roi[i] - strt_con_y

                # Calculate Theta
                theta = math.atan(del_y / del_x)
                if (2 * r_sc) >= (((del_x) ** 2 + (del_y) ** 2) ** 0.5):
                    beta = math.asin((((del_x) ** 2 + (del_y) ** 2) ** 0.5)
                                     / (2 * r_sc))
                else:
                    beta = math.pi / 2
                delta = beta - theta
                deltas.append(delta * 180 / math.pi)

            if len(x_roi) != 0:

                # FORWARD CONTACTS
                strt_con_y = y_roi[deltas.index(min(deltas))]  # minimum delta
                strt_con_x = x_roi[deltas.index(min(deltas))]
                j = np.where(strt_con_x == lat)[0]  # update j for while loop
            else:
                # Append last value if there is no region of interest
                strt_con_y = ne[len(ne) - 1]
                strt_con_x = lat[len(ne) - 1]
                j = len(ne)
    contact_x = np.array(f_con_xs)
    contact_y = np.array(f_con_ys)

    return contact_x, contact_y


def full_barrel(x_sc, ne_sc, barrel_start, envelope=True, envelope_lower=0.06,
                envelope_upper=0.02, small_large=False, num_barrels=10,
                filter_barrel=True, svg_window=5, svg_poly=2):
    """ Roll barrel over an array to get the forward and backward rollings
    over 10 barrel radius sizes

    Parameters
    ----------
    x_sc: array-like
        scaled x array
        could be latitude or time
    ne_sc: array-like
        scaled density array
    barrel_start : double
        starting radius of the barrel
    envelope : bool kwarg
        default true, envelope will be used
        if false, no envelope will be used
    envelope_lower : double kwarg
        lower limit of envelope
        default 0.06 (.6%) of min value from contact points
    envelope_upper : double kwarg
        upper limit of envelope
        default 0.02 (.2%) of max value from contact points
    small_large : bool
        if True, smallest barrel is used first leading to largest
        False opposite (default)
    num_barrels : int
        number of barrel radii for rolling
        default is 10 (10 sizes/ 10 rollings)
    filter_barrel : boolean
        if True (default), barrel will be filtered after all rollings
    svg_window : int
        savitzky golay filter window length
    svg_poly : int
        savitzky golay filter poly order
    Returns
    -------
    filt_y : array-like
        trend, filtered scaled ne with same length as ne_sc
    gap_ind_st_ed : array-like
        indices of the gap locations (where barrel did not roll)
    """

    # Create array of zeros with length of x
    track = [0] * len(x_sc)

    # df2 keeps track of the gaps
    df2 = pd.DataFrame()

    # Initialize parameters
    # xs is the original x array
    df2['xs'] = x_sc

    # counter tracks of the number of times the contact points were contacted
    df2['counter'] = track

    # zs keeps track of how many times the values were 0 (not contacted),
    # resets every new barrel size, changes to 1 if contacted
    df2['zs'] = track

    # zs_previous keeps track of the previous z value
    # so that we know if it is consecutive or not
    df2['zs_previous'] = track

    # keeps track of how many consecutive times the zs is 0 and not 1,
    # z_counter = 3 means 2 surrounding points were set in stone
    # reset if it is not consecutive
    df2['z_counter'] = track

    # so if they are maintained for 3 consecutive rollings,
    # they get a 1 flag
    df2['x_flags'] = track

    # keeps track of previous iteration of x_flag
    df2['x_flags_old'] = track

    # iterate through 10 barrel sizes
    for k in range(num_barrels):

        # zs keeps track of how many times the values were 0,
        # resets every new barrel size
        df2['zs'] = track

        # Call simple_barrel
        if small_large:  # start small go large
            r_sc = barrel_start * (k + 1)
        else:  # start large go small
            r_sc = barrel_start * (num_barrels - k)

        f_con_xs, f_con_ys = simple_barrel(x_sc, ne_sc, r_sc,
                                           direction='forward')
        b_con_xs, b_con_ys = simple_barrel(x_sc, ne_sc, r_sc,
                                           direction='backward')

        # unite forward and backward rolling contact points
        x_b_f = np.concatenate((f_con_xs, b_con_xs))

        # Sort combined data by x values and take out repeats
        bf_sorted_indices = np.argsort(x_b_f)  # in order
        x_sort = x_b_f[bf_sorted_indices]

        # the repeat values should be the same for both xs and ys
        unique_indices = np.unique(x_sort, return_index=True)[1]
        contact_xs = x_sort[unique_indices]

        for xx in contact_xs:
            df2.loc[df2['xs'] == xx, 'counter'] += 1
            df2.loc[df2['xs'] == xx, 'zs'] = 1

        # if zs_previous == 0. then z_counter +=1 else, z_counter = 0 again
        # for zz in range(len(df2['zs'])):
        # this counts how many times the gaps are seen. if not successive,
        # set to 0. if successive +1 for z_counter (zero counter)
        for zz in range(len(df2['xs'])):

            # round to nearest whole number
            zz = round(zz)

            # if previous and the current rolling are 0s, then add to counter
            if (df2['zs_previous'].loc[zz] == 0) & (df2['zs'].loc[zz] == 0):
                df2.loc[zz, 'z_counter'] += 1

            # if the previous is not a gap (1) and the current is a gap (0),
            # reset counter at 1
            elif (df2['zs_previous'].loc[zz] == 1) & (df2['zs'].loc[zz] == 0):
                df2.loc[zz, 'z_counter'] = 1

            # if they both rollings are not gaps, set to 0, so not counted
            elif (df2['zs_previous'].loc[zz] == 1) & (df2['zs'].loc[zz] == 1):
                df2.loc[zz, 'z_counter'] = 0

            # if the previous is a gap, but the current is not, set to 0.
            elif (df2['zs_previous'].loc[zz] == 0) & (df2['zs'].loc[zz] == 1):
                df2.loc[zz, 'z_counter'] = 0

        # this loop is flagging gaps that have a z_counter of at least 3.
        for zz in range(len(df2['xs'])):
            zz = round(zz)
            if df2['z_counter'].loc[zz] >= 3:

                # not previously a gap
                # No need to check 1 since it is already done
                if df2['x_flags'].loc[zz] == 0:

                    # this tells us we are not merging two gaps to be 1 big gap
                    # because the previous iteration shows the value before
                    # OR after are not gaps.
                    if ((df2['x_flags_old'].loc[zz - 1] == 0)
                            or (df2['x_flags_old'].loc[zz + 1] == 0)):
                        df2.loc[zz, 'x_flags'] = 1
                        # if the previous AND next values =1.
                        # then we leave it as a 0,
                        # so we are not merging 2 flagged regions

        # now equals previous z values going into next one
        df2['zs_previous'] = df2['zs']
        df2['x_flags_old'] = df2['x_flags']

    # ONLY INCLUDE TIMES of x_flags == 0
    # and continue on from here outside of the loop
    # x_flags == 1 are gaps, and x_flags == 0 is contact points
    fin_contact_xs = x_sc[df2['x_flags'] == 0]
    fin_contact_ys = ne_sc[df2['x_flags'] == 0]

    int_y = np.interp(np.setdiff1d(x_sc, fin_contact_xs),
                      fin_contact_xs, fin_contact_ys)  # linear interpolation

    x_combined = np.concatenate((fin_contact_xs,
                                 np.setdiff1d(x_sc, fin_contact_xs)))
    y_combined = np.concatenate((fin_contact_ys, int_y))

    # Sort combined data by x values
    sorted_indices = np.argsort(x_combined)
    x_combined = x_combined[sorted_indices]
    y_combined = y_combined[sorted_indices]

    # collect all data points between envelopes,
    # to do this code wise i need to close all gaps in code
    # and have same length for original curve and other ones.
    ne_cont = np.array(y_combined)

    if envelope:  # use an envelope
        BRC_upper = ne_cont + envelope_upper * max(ne_cont)
        BRC_lower = ne_cont - envelope_lower * min(ne_cont)
        a = np.empty((len(y_combined)))
        a[:] = np.nan

        # this is a combination of nan and values inside envelope
        a[(ne_sc < BRC_upper) & (ne_sc > BRC_lower)] = ne_sc[(
            (ne_sc < BRC_upper) & (ne_sc > BRC_lower))]

        # Check to see if there are nan and non-nan values
        if len(x_sc[~np.isnan(a)]) != 0:
            # Check to see if there are NaN values left from a
            int_y2 = np.interp(x_sc[np.isnan(a)],
                               x_sc[~np.isnan(a)], a[~np.isnan(a)])

            # replace nan values with interpolated values
            a[np.isnan(a)] = int_y2
        else:
            a = np.array(y_combined)
    else:
        a = np.array(y_combined)

    # if the data is long enough apply a savitzky-golay filter

    if filter_barrel:
        if len(a) > 4:
            filt_y = savgol_filter(a, svg_window, svg_poly)
        else:
            filt_y = a

    # Find the gap ranges
    gap_ind_st_ed = barrel_utils.find_nan_ranges(a)

    return filt_y, gap_ind_st_ed


def get_barrel_info(x_sec, ne, barrel_start=8, det_filt=5, width=3,
                    envelope=True, envelope_lower=0.06, envelope_upper=0.02,
                    small_large=False, num_barrels=10, up=True,
                    prom_height=None, get_peaks=True, scale_n=True,
                    exp_inc=True, ex_scale=2, filter_barrel=True, svg_window=5,
                    svg_poly=2, x_type='time'):
    """Get Barrel dataframe, peaks, and properties.
    Parameters
    ----------
    x_sec: array-like
        non-scaled time in seconds (not dataframe)
    ne : array-like
        non-scaled density array (not dataframe)
        already clean from "bad data points"
    barrel_start : int kwarg
        starting radius of the barrel, default is
        8 (for 1 Hz data, 8 seconds). The barrel will
        range from 8-80 seconds
    det_filt : int kwarg
        filter for unifrom_filter1d for detrended ne
        default is 5 data points (5 seconds for 1 Hz)
    width : int kwarg
        minimum width of a peak for find_peaks
        a peak has to last at least width xs
        to be considered a peak
        default is 3 (1 Hz data, 3 seconds or 22.5 km)
    envelope_lower : double kwarg
        lower limit of envelope
        default 0.06 (.6%) of min value from contact points
    envelope_upper : double kwarg
        upper limit of envelope
        default 0.02 (.2%) of max value from contact points
    small_large : bool
        if True, smallest barrel is used first leading to largest
        False opposite (default)
    num_barrels : int
        number of barrel radii for rolling
        default is 10 (10 sizes/ 10 rollings)
    up : boolean
        True if this is the upper barrel (default)
        False if this is the lower barrel
    prom_height : NoneType or float
        height of prominence for find_peaks
        if None, no minimum will be used
    get_peaks : boolean
        if True (default), peaks and proeprties will be returned
        if False, only barrel_df will be returned
    scale_n : boolean
        if True (default), density will be scaled
            by x_span / (np.nanmax(ne) * 10**exp)
    exp_inc : boolean
        if True (default), ne scaled by x_span / (np.nanmax(ne) * 10**exp)
            where exp = int(math.log10(x_span)) - ex_scale
        if False, ne scaled by x_span / (np.nanmax(ne) * 10**0)
    ex_scale : float
        scaling param for the exponent.
        exp = int(math.log10(x_span)) - ex_scale
    filter_barrel : boolean
        if True (default), barrel will be filtered after all rollings
    svg_window : int
        savitzky golay filter window length
    svg_poly : int
        savitzky golay filter poly order
    x_type : string
        default 'time', will start at 0
        if not time but rather latitude or longitude, no conversion

    Returns
    ------
    barrel_df : dataframe
        dataframe that includes: ne, filt_ne, detrended_ne, and
        filt_det_ne
        all nan values have been put back into all of those
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    """
    # if lower barrel, turn the data upside down
    ne_copy = ne.copy()
    if not up:
        ne = np.nanmax(ne_copy) - ne_copy + np.nanmin(ne_copy)

    # Interpolate nans
    nans = np.isnan(ne)
    ne[nans] = np.interp(np.flatnonzero(nans), np.flatnonzero(~nans),
                         ne[~nans])

    # x_sec and ne should be an array not a dataframe
    # Scale x to make it start at 0
    if x_type.lower() == 'time':
        x_sc = x_sec - x_sec[0]

        # adjust axes if one day leads into another
        if x_sec[0] > x_sec[-1]:
            x_sc = np.linspace(0, len(x_sec) - 1, len(x_sec))
    else:
        x_sc = x_sec

    # Scale the density ---
    if scale_n:
        x_span = np.nanmax(x_sc) - np.nanmin(x_sc)
        if exp_inc:
            exp = int(math.log10(x_span)) - ex_scale
        else:
            exp = 0
        scal_param = x_span / (np.nanmax(ne) * 10**exp)
        ne_sc = ne * scal_param
    else:
        ne_sc = ne
        scal_param = 1

    # Roll the barrel ---
    filt_y, gap_ind_st_ed = full_barrel(
        x_sc, ne_sc, barrel_start, envelope=envelope,
        envelope_lower=envelope_lower, envelope_upper=envelope_upper,
        small_large=small_large, num_barrels=num_barrels,
        filter_barrel=filter_barrel, svg_window=svg_window, svg_poly=svg_poly)

    # Rescale y
    filt_ne = filt_y / scal_param

    # upper or lower barrel
    if up:  # upper
        # Detrend Ne where valleys are now peaks and apply a filter
        detrended_ne = filt_ne - ne
        ne_save = ne
        # Filter
        if det_filt != 0:
            filt_det_ne = uniform_filter1d(detrended_ne, size=det_filt)
        else:
            filt_det_ne = detrended_ne
    else:  # lower
        # flip data right side up
        filt_ne = np.nanmax(ne_copy) - filt_ne + np.nanmin(ne_copy)
        ne_save = np.nanmax(ne_copy) - ne + np.nanmin(ne_copy)
        detrended_ne = filt_ne - ne_save

        # Filter
        if det_filt != 0:
            filt_det_ne = uniform_filter1d(detrended_ne, size=det_filt)
        else:
            filt_det_ne = detrended_ne

    # Replace nans
    filt_ne[nans] = np.nan
    detrended_ne[nans] = np.nan
    filt_det_ne[nans] = np.nan

    # create dataframe to return
    barrel_df = pd.DataFrame()

    # replace ne nans
    ne[nans] = np.nan
    barrel_df['ne'] = ne_save
    barrel_df['barrel_ne'] = filt_ne
    barrel_df['dne_dd'] = detrended_ne
    barrel_df['dne_dd_filt'] = filt_det_ne
    barrel_df['ne_scale'] = ne_sc
    barrel_df['x_scale'] = x_sc

    if get_peaks:
        # use find_peaks
        peaks, properties = find_peaks(filt_det_ne, width=width,
                                       prominence=prom_height)
        return barrel_df, peaks, properties
    else:
        return barrel_df


def triple_barrel(data_df, barrel_start=8, det_filt=5, peak_width=3,
                  num_barrels=5, prom_height=5000, x_type='time', x_str='time',
                  ne_str='ne', e_up_large=0.1, e_lo_large=0.05, e_up_small=0.5,
                  e_lo_small=0.1, upper_weight=2, lower_weight=1, scale_n=True,
                  exp_inc=True, ex_scale=2, filter_barrel=True, svg_window=5,
                  svg_poly=2, nflag_mask=None):
    """Roll multiple barrels over/under the data.

    Parameters
    ----------
    data_df : DataFrame
        Containing density and time or latitude
    barrel_start : double
        starting radius of the barrel
        Largest rolling : 2*barrel_start * num_barrels
            decreasing to 2* barrel_start
        Smallest rolling : barrel_start
            increasing to barrel_start * num_barrels
    det_filt : int kwarg
        filter for unifrom_filter1d for detrended ne
        default is 5 data points (5 seconds for 1 Hz)
    peak_width : int kwarg
        minimum width of a peak for find_peaks
        a peak has to last at least width xs
        to be considered a peak
        default is 3 (1 Hz data, 3 seconds or 22.5 km for Swarm)
    num_barrels : int
        number of barrel radii for rolling
        default is 5 (5 sizes/ 5 rollings)
    prom_height : NoneType or float
        height of prominence for find_peaks
        if None, no minimum will be used
        default 5000
    x_type : string
        default 'time', will convert to seconds
        if not time but rather latitude or longitude, no conversion
    x_str : string
        x axis column name string for rolling over
        default : 'time'
    ne_str : string
        density string
        default : 'ne'
    e_up_large, e_lo_large : double kwarg
        upper and lower limit of envelope for large upper barrel
        default 0.1 (1%) and 0.05 (0.5%) of min value from contact points
    e_lo_small, e_up_small : double kwarg
        upper and lower limit of envelope for small upper and lower barrel
        default 0.5 (5%) and 0.05 (0.5%) of min value from contact points
    upper_weight : float
        weight to put on upper barrel
        default 2
    lower_weight : float
        weight to put on lower barrel
        default 1
        (upper_weight * up_barrel + lower_weight * lo_barrel) / (sum(weights))
    scale_n : boolean
        if True (default), density will be scaled
            by x_span / (np.nanmax(ne) * 10**exp)
    exp_inc : boolean
        if True (default), ne scaled by x_span / (np.nanmax(ne) * 10**exp)
            where exp = int(math.log10(x_span)) - ex_scale
        if False, ne scaled by x_span / (np.nanmax(ne) * 10**0)
    ex_scale : float
        scaling param for the exponent.
        exp = int(math.log10(x_span)) - ex_scale
    filter_barrel : boolean
        if True (default), barrel will be filtered after all rollings
    svg_window : int
        savitzky golay filter window length
    svg_poly : int
        savitzky golay filter poly order
    nflag_mask : array-like or NoneType
        density flags mask
        where True values are good data
        False, bad data
        default None, no density flag assumed
    Returns
    -------
    barrel_df : DataFrame
        Barrel info for the largest upper barrel with added columns
        for narlo_barrel, narup_barrel, dne_wdr, weight_barrel
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    """

    # Set up the x data
    if x_type.lower() == 'time':
        x_secdf = (data_df[x_str].dt.hour * 3600
                   + data_df[x_str].dt.minute * 60
                   + data_df[x_str].dt.second
                   + data_df[x_str].dt.microsecond / 10 ** 6)
        x_sec = x_secdf.values
    else:
        x_sec = data_df[x_str].values

    # set up the y data
    ne = data_df[ne_str].values

    # MAIN BARREL--------------------------------------------------------------
    # Upper barrel moving from largest radius to smallest
    barrel_large = barrel_start * 2
    e_up = 0.1
    e_lo = 0.05

    # Upper barrel
    barrel_up, peaks, properties = get_barrel_info(
        x_sec, ne, barrel_start=barrel_large, det_filt=det_filt,
        width=peak_width, envelope_lower=e_lo_large, envelope_upper=e_up_large,
        small_large=False, num_barrels=num_barrels, up=True,
        prom_height=prom_height, get_peaks=True, scale_n=scale_n,
        exp_inc=exp_inc, ex_scale=ex_scale, filter_barrel=filter_barrel,
        svg_window=svg_window, svg_poly=svg_poly, x_type=x_type)

    barrel_up['x'] = data_df[x_str].values

    if nflag_mask is not None:
        barrel_up['ne_flag'] = nflag_mask

    barrel_up = barrel_up.copy()
    barrel_up['div_det_ne'] = (barrel_up['ne'].values
                               / barrel_up['barrel_ne'].values)

    # SMALL UPPER BARREL-------------------------------------------------------
    # moving from smallest radius to largest
    # Barrel Roll
    barrel_up2 = get_barrel_info(
        x_sec, ne, barrel_start=barrel_start, det_filt=det_filt,
        width=peak_width, envelope_lower=e_lo_small, envelope_upper=e_up_small,
        small_large=True, num_barrels=num_barrels, up=True,
        prom_height=prom_height, get_peaks=False, scale_n=scale_n,
        exp_inc=exp_inc, ex_scale=ex_scale, filter_barrel=filter_barrel,
        svg_window=svg_window, svg_poly=svg_poly)

    # SMALL LOWER BARREL-------------------------------------------------------
    # moving from smallest radius to largest
    # Barrel Roll
    barrel_lo = get_barrel_info(
        x_sec, ne, barrel_start=barrel_start, det_filt=det_filt,
        width=peak_width, envelope_lower=e_lo, envelope_upper=e_up,
        small_large=True, num_barrels=num_barrels, up=False, get_peaks=False,
        scale_n=scale_n, exp_inc=exp_inc, ex_scale=ex_scale,
        filter_barrel=filter_barrel, svg_window=svg_window, svg_poly=svg_poly)

    # Barrel dataframe with weighted barrel and lower density
    # weighted Barrel
    barrel_df = barrel_up.copy()
    total_weight = upper_weight + lower_weight
    barrel_df['weight_barrel'] = ((upper_weight
                                   * barrel_up2['barrel_ne'].values
                                   + lower_weight
                                   * barrel_lo['barrel_ne'].values)
                                  / total_weight)
    barrel_df['dne_wdr'] = (data_df[ne_str].values
                            / barrel_df['weight_barrel'].values)
    barrel_df['narlo_barrel'] = barrel_lo['barrel_ne'].values
    barrel_df['narup_barrel'] = barrel_up2['barrel_ne'].values

    return barrel_df, peaks, properties


def just_barrel(IRR_data, barrel_start=8, det_filt=5, peak_width=3,
                x_type='time', ne_str='ne', x_str='time', scale_n=True,
                exp_inc=True, ex_scale=2, filter_barrel=True, svg_window=5,
                svg_poly=2, nflag_mask=None):
    """ Get the peaks, properties, and flags in one step
    Parameters
    ----------
    IRR_data : DataFrame
        Swarm IRR data (1 second cadence)
    barrel_start : int kwarg
        starting radius of the barrel, default is
        8 (for 1 Hz data, 8 seconds). The barrel will
        range from 8-80 seconds
    det_filt : int kwarg
        filter for unifrom_filter1d for detrended ne
        default is 5 data points (5 seconds for 1 Hz)
    width : int kwarg
        minimum width of a peak for find_peaks
        a peak has to last at least width xs
        to be considered a peak
        default is 3 (1 Hz data, 3 seconds or 22.5 km)
    x_type : string
        default 'time', will convert to seconds
        if not time but rather latitude or longitude, no conversion
    ne_str, x_str : strings
        column names of density, flag params, and time param
        defaults: 'ne' and 'time' respectively
    scale_n : boolean
        if True (default), density will be scaled
            by x_span / (np.nanmax(ne) * 10**exp)
    exp_inc : boolean
        if True (default), ne scaled by x_span / (np.nanmax(ne) * 10**exp)
            where exp = int(math.log10(x_span)) - ex_scale
        if False, ne scaled by x_span / (np.nanmax(ne) * 10**0)
    ex_scale : float
        scaling param for the exponent.
        exp = int(math.log10(x_span)) - ex_scale
    filter_barrel : boolean
        if True (default), barrel will be filtered after all rollings
    svg_window : int
        savitzky golay filter window length
    svg_poly : int
        savitzky golay filter poly order
    nflag_mask : array-like or NoneType
        density flags mask
        where True values are good data
        False, bad data
        default None, no density flag assumed
    Retruns
    -------
    barrel_df : dataframe
        dataframe that includes: ne, filt_ne, detrended_ne, filt_det_ne,
        ne_scale, x_scale, time, ne_flag
        all nan values have been put back into all of those
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    """
    # Set up the x data
    x_secdf = (IRR_data[x_str].dt.hour * 3600
               + IRR_data[x_str].dt.minute * 60
               + IRR_data[x_str].dt.second
               + IRR_data[x_str].dt.microsecond / 10 ** 6)
    x_sec = x_secdf.values

    # set up the y data
    ne_copy = IRR_data[ne_str].copy()
    ne = ne_copy.values

    # Barrel Roll
    barrel_df, peaks, properties = get_barrel_info(
        x_sec, ne, barrel_start=barrel_start, det_filt=det_filt,
        width=peak_width, scale_n=scale_n, exp_inc=exp_inc, ex_scale=ex_scale,
        filter_barrel=filter_barrel, svg_window=svg_window, svg_poly=svg_poly,
        x_type=x_type)

    barrel_df['x'] = IRR_data[x_str].values
    if nflag_mask is not None:
        barrel_df['ne_flag'] = nflag_mask

    return barrel_df, peaks, properties
