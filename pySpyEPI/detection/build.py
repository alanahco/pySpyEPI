"""Detection building functions
created 24 February 2026
alanahco

List Functions
--------------
find_cand
"""
import pandas as pd
from pySpyEPI.detection import midlat_trough as mid_tr
from pySpyEPI.barrel import barrel_roll
from pySpyEPI.detection.stats import cand_info


def find_cand(
        data_df, satellite, sat_num, trough=True, equator_bound=35,
        auroral_bound=70, set_lat=50, num_barrels=5, barrel_start=8,
        det_filt=5, peak_width=6, x_type='time', x_str='time',
        e_up_large=0.1, e_lo_large=0.05, e_up_small=0.1, e_lo_small=0.05,
        upper_weight=2, lower_weight=1, scale_n=True, exp_inc=True, ex_scale=2,
        filter_barrel=True, svg_window=5, svg_poly=2, freq=1, fft_test=False,
        nest_test=True, prom_test=True, percent_test=True, wdr_test=True,
        utest_mask=None, min_height=10**4, min_perc=10, wdr50_set=0.995,
        wdr20_set=0.95, eia_state='unknown', pass_id=None, prom_height=5000,
        time_str='time', ne_str='ne', glat_str='lat', glon_str='lon',
        mlat_str='mlat', mlon_str=None, alt_str=None, lshell_str=None,
        lt_str='lt', nflag_mask=None):
    """Find EPI candidates.

    Parameters
    ----------
    time, lat, lon, mlat, mlon, lt, lshell, alt : array-like
        arrays of time, geographic latitude, geographic longitude,
        magnetic latitude, magnetic longitude, local time, L shell,
        and altitude
    satellite : string
        observation satellite ('swarm', 'dmsp')
    sat_num : string
        satellite number or letter
    equator_bound : float kwarg (35)
        equatorward magnetic latitude to search down to
        set as a positive number, will be used as +/-equator_bound
    auroral_bound : float kwarg (75)
        autroalward magnetic latitude to search up to
        set as a positive number, will be used as +/-auroral_bound
    set_lat : float kwarg (55)
        if a trough cannot be found, set_lat will be used in its place
        set as a positive number, will be used as +/-set_lat
    num_barrels : int
        number of barrel radii for rolling
        default is 5 (5 sizes/ 5 rollings)
    barrel_start : int kwarg
        starting radius of the barrel, default is
        8 (for 1 Hz data, 8 seconds). The barrel will
        range from 8-80 seconds
    det_filt : int kwarg
        filter for unifrom_filter1d for detrended ne
        default is 5 data points (5 seconds for 1 Hz)
    peak_width : int kwarg
        minimum width of a peak for find_peaks
        a peak has to last at least width xs
        to be considered a peak
        default is 6 (1 Hz data, 6 seconds or 45 km)
    x_type : string
        default 'time', will convert to seconds
        if not time but rather latitude or longitude, no conversion
    x_str : string
        x axis column name string for rolling over
        default : 'time'
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
    freq : int
        data frequency
        default 1 Hz
    fft_test : bool
        default False
        fft of delta ne 20s for data will be calculated and erroneous
        data will be flagged (SWARM specific)
    nest_test : bool
        default True
        if True, test for nested bubbles
    prom_test : bool
        default True
        if True, test for a minimum prominence
    percent_test : bool
        if True, test for a minimum percent depth
    wdr_test : bool
        if True, test the weighted density ratio
    utest_mask : array-like or NoneType
        Additional test from user as a mask
        True values are assumed to be good and False values assumed to be for
        removal
        array must be same length as ne
    min_height: float
        minimum height for a prominence to be considered a prominence
        default 10**4
    min_perc : float
        minimum percent depth
        default 10%
    wdr50_set : float
        default 0.995
        minimum wdr for median of EPI WDR
    wdr20_set : float
        default 0.95
        if 25th percentile WDR of curve is below wdr20_set, use wdr20_set to
        test 20th percentile of EPI WDR
    eia_state : string
        eia category/orientation string from pyValEIA
    pass_id : string
        satelltie pass identifier
    prom_height : NoneType or float
        height of prominence for find_peaks
        if None, no minimum will be used
        default 5000
    time_str, ne_str, glat_str, glon_str, mlat_str, mlon_str, alt_str,
        lshell_st, lt_str : strings
        column names from data_df for time, density, geographic lat and lon,
        magnetic lat and lon, altitude, lshell, and lt_str
        alt_str, mlon_str, and lshell_str can all be NoneTypes as they are
        not requried for analysis.
    nflag_mask : array-like or NoneType
        density flags mask
        where True values are good data
        False, bad data
        default None, no density flag assumed
    Returns
    -------
    barrel_df : dataframe
        dataframe that includes: ne, filt_ne, detrended_ne, filt_det_ne,
        ne_scale, x_scale, time
        all nan values have been put back into all of those
    peaks, properties: array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    trough_lats : array of length 2
        northern and southern hemispehre midlatitude trough location
    pass_id : string
        pass identifier of format '%Y%b%d_%H%M_sat{sat}' for future reference
        time is based on 0 degrees mangetic latitude
    trough_lats : array-like
        array of length 2 with the northern and southern trough bounds
    """
    # if pass_id is not given
    if pass_id is None:
        # name pass ybd_HM_sat{sat} based on 0 degrees mlat
        z_id = abs(data_df[mlat_str].values).argmin()
        ztime = data_df[time_str].iloc[z_id]
        t_id = pd.Timestamp(ztime).to_pydatetime().strftime('%Y%m%d_%H%M')
        pass_id = t_id + f'_{satellite}{sat_num}'

    if trough:
        # Segment the data using the found trough latitudes -------------------
        trough_lats = mid_tr.find_trough_lat(
            data_df, equator_bound=equator_bound, auroral_bound=auroral_bound,
            set_lat=set_lat, ne_str=ne_str, mlat_str=mlat_str)
    else:
        trough_lats = [abs(set_lat), -1 * abs(set_lat)]

    # magnetic latitudinally filtered data ------------------------------------
    lat_mask = ((data_df[mlat_str] >= min(trough_lats))
                & (data_df[mlat_str] <= max(trough_lats)))
    data_df_lat = data_df[lat_mask]

    # Roll barrel and get peak info -------------------------------------------
    barrel_df, peaks, properties = barrel_roll.triple_barrel(
        data_df_lat, barrel_start=barrel_start, det_filt=det_filt,
        peak_width=peak_width, num_barrels=num_barrels,
        prom_height=prom_height, x_type=x_type, x_str=x_str, ne_str=ne_str,
        e_up_large=e_up_large, e_lo_large=e_lo_large, e_up_small=e_up_small,
        e_lo_small=e_lo_small, upper_weight=upper_weight,
        lower_weight=lower_weight, scale_n=scale_n, exp_inc=exp_inc,
        ex_scale=ex_scale, filter_barrel=filter_barrel, svg_window=svg_window,
        svg_poly=svg_poly, nflag_mask=nflag_mask)

    # separate the dataframe for epi_stats ------------------------------------
    time2 = data_df_lat[time_str].values
    time2 = time2.astype('datetime64[us]')
    lat2 = data_df_lat[glat_str].values
    lon2 = data_df_lat[glon_str].values
    mlat2 = data_df_lat[mlat_str].values
    lt2 = data_df_lat[lt_str].values
    lt2 = lt2.astype('datetime64[us]')
    ne2 = data_df_lat[ne_str].values

    # for parameters that are nice to have but not required, None is an option
    if lshell_str is not None:
        lshell2 = data_df_lat[lshell_str].values
    else:
        lshell2 = None
    if mlon_str is not None:
        mlon2 = data_df_lat[mlon_str].values
    else:
        mlon2 = None
    if alt_str is not None:
        alt2 = data_df_lat[alt_str].values
    else:
        alt2 = None

    if nflag_mask is not None:
        nflag_lat = nflag_mask[lat_mask]
    else:
        nflag_lat = None

    # get candidate info ------------------------------------------------------
    stats_df = cand_info.epi_stats(
        time2, lat2, lon2, mlat2, mlon2, lt2, lshell2, alt2, peaks, properties,
        barrel_df, sat_num, nflag_mask=nflag_lat, freq=freq, fft_test=fft_test,
        nest_test=nest_test, prom_test=prom_test, percent_test=percent_test,
        wdr_test=wdr_test, utest_mask=utest_mask, min_height=min_height,
        min_perc=min_perc, wdr50_set=wdr50_set, wdr20_set=wdr20_set)

    # save pass info ----------------------------------------------------------
    plen = len(peaks)
    pass_df = cand_info.pass_info(time2, mlat2, ne2, plen, pass_id,
                                  eia_state=eia_state)

    # combine stats_df and pass_df, ignoring index and align row by row
    cand_df = pd.concat(
        [stats_df.reset_index(drop=True),
         pass_df.reset_index(drop=True)], axis=1)

    return barrel_df, cand_df, peaks, properties, pass_id, trough_lats
