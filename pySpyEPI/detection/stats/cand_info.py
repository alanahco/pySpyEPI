"""Save peak info.
Created 19 February 2026
alanahco
List Functions
--------------
center_stats
"""
import pandas as pd
import numpy as np
from pySpyEPI.detection.test_candidates import test_peaks
from pySpyEPI.detection.flag_detections import mlat_flag


def epi_stats(time, lat, lon, mlat, mlon, lt, lshell, alt, peaks, properties,
              barrel_df, satellite, nflag_mask=None, freq=1, fft_test=False,
              nest_test=True, prom_test=True, percent_test=True, wdr_test=True,
              utest_mask=None, min_height=10**4, min_perc=10, wdr50_set=0.995,
              wdr20_set=0.95):
    """Collect information about candidate EPIs.
    Parameters
    ----------
    time, lat, lon, mlat, mlon, lt, lshell, alt : array-like
        arrays of time, geographic latitude, geographic longitude,
        magnetic latitude, magnetic longitude, local time, L shell,
        and altitude
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    barrel_df : dataframe
        dataframe created by barrel.barrel_roll.get_barrel_info
        and updated by barrel.barrel_roll.just_barrel
    satellite: str
        swarm 'A', 'B' or 'C'
    nflag_mask : array-like or NoneType
        density flags mask
        where True values are good data
        False, bad data
        default None, no density flag assumed
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
    Returns
    ------
    stats_df : dataframe
        includes info about candidates and candidate tests
        time of center and edges
        plasma density (np) of center and edges
        center prominence
        background density at center (from barrel trend)
        magnetic latitude of center and edges
        magnetic longitude of center and edges
        geographic latitude of center and edges
        geographic longitude of center and edges
        lshell of center and edges
        local time of center and edges
        barrel slope
        altitude
        satellite
        mlat_flag
    """

    # ititiate columns and create dataframe
    columns = ['center_time', 'edge1_time', 'edge2_time', 'center_np',
               'edge1_np', 'edge2_np', 'center_prom',
               'center_bnp', 'center_mlat', 'edge1_mlat', 'edge2_mlat',
               'center_mlon', 'edge1_mlon', 'edge2_mlon', 'center_lat',
               'edge1_lat', 'edge2_lat', 'center_lon', 'edge1_lon',
               'edge2_lon', 'center_lshell', 'edge1_lshell', 'edge2_lshell',
               'center_lt', 'edge1_lt', 'edge2_lt', 'barrel_slope',
               'barrel_slope_left', 'barrel_slope_right', 'altitude',
               'satellite', 'fpeak_prom']

    df = pd.DataFrame(columns=columns)

    # get the edges of the peaks
    xleft = np.trunc(properties["left_ips"]).astype(int)
    xright = np.ceil(properties["right_ips"]).astype(int)
    fpeak_prom = properties['prominences']
    # iterate through peaks
    for i, xp in enumerate(peaks):
        xr = xright[i]
        xl = xleft[i]

        # get times of peaks
        l_time = time[xl]
        c_time = time[xp]
        r_time = time[xr]

        df.at[i, 'center_time'] = pd.Timestamp(c_time).to_pydatetime()
        df.at[i, 'edge1_time'] = pd.Timestamp(l_time).to_pydatetime()
        df.at[i, 'edge2_time'] = pd.Timestamp(r_time).to_pydatetime()

        # density
        df.at[i, 'center_np'] = barrel_df['ne'].iloc[xp]
        df.at[i, 'edge1_np'] = barrel_df['ne'].iloc[xl]
        df.at[i, 'edge2_np'] = barrel_df['ne'].iloc[xr]

        # save prominence
        df.at[i, 'center_prom'] = barrel_df["dne_dd_filt"].iloc[xp]
        df.at[i, 'center_bnp'] = barrel_df["barrel_ne"].iloc[xp]

        # Coordinates at peak
        df.at[i, 'center_mlat'] = mlat[xp]
        df.at[i, 'center_lat'] = lat[xp]
        df.at[i, 'center_lon'] = lon[xp]
        df.at[i, 'center_lt'] = pd.Timestamp(lt[xp]).to_pydatetime()

        # coordinates at edges
        df.at[i, 'edge1_mlat'] = mlat[xl]
        df.at[i, 'edge1_lat'] = lat[xl]
        df.at[i, 'edge1_lon'] = lon[xl]
        df.at[i, 'edge1_lt'] = pd.Timestamp(lt[xl]).to_pydatetime()

        df.at[i, 'edge2_mlat'] = mlat[xr]
        df.at[i, 'edge2_lat'] = lat[xr]
        df.at[i, 'edge2_lon'] = lon[xr]
        df.at[i, 'edge2_lt'] = pd.Timestamp(lt[xr]).to_pydatetime()

        # values that are nice to have but not required
        if mlon is not None:
            df.at[i, 'center_mlon'] = mlon[xp]
            df.at[i, 'edge1_mlon'] = mlon[xl]
            df.at[i, 'edge2_mlon'] = mlon[xr]
        if lshell is not None:
            df.at[i, 'center_lshell'] = lshell[xp]
            df.at[i, 'edge1_lshell'] = lshell[xl]
            df.at[i, 'edge2_lshell'] = lshell[xr]
        if alt is not None:
            df.at[i, 'altitude'] = np.nanmedian(alt[xl:xr + 1])

        # barrel slope
        rise = barrel_df['ne'].iloc[xr] - barrel_df['ne'].iloc[xl]
        run = xr - xl
        df.at[i, 'barrel_slope'] = rise / run

        rise = barrel_df['ne'].iloc[xp] - barrel_df['ne'].iloc[xl]
        run = xp - xl
        df.at[i, 'barrel_slope_left'] = rise / run

        rise = barrel_df['ne'].iloc[xr] - barrel_df['ne'].iloc[xp]
        run = xr - xp
        df.at[i, 'barrel_slope_right'] = rise / run

        df.at[i, 'satellite'] = satellite
        df.at[i, 'fpeak_prom'] = fpeak_prom[i]

    # get peak flag info from process_peaks
    flag_df = test_peaks(barrel_df['ne'].values, peaks, properties,
                         barrel_df, nflag_mask=nflag_mask, freq=freq,
                         fft_test=fft_test, nest_test=nest_test,
                         prom_test=prom_test, percent_test=percent_test,
                         wdr_test=wdr_test, utest_mask=utest_mask,
                         min_height=min_height, min_perc=min_perc,
                         wdr50_set=wdr50_set, wdr20_set=wdr20_set)

    # flag by magnetic latitude
    ml_flag = mlat_flag(df['center_mlat'].values,
                        flag_df['epi_flag'].values)

    # combine dataframes and add mlat flag
    stats_df = pd.concat([df, flag_df], axis=1)
    stats_df['mlat_flag'] = ml_flag

    # save max and min Ne for different latitudinal regions of barrel Ne
    mlat_lims = [-50, -40, -25, 0, 25, 40, 50]

    # set up lengths
    plen = len(peaks)

    # set up strings
    p_strings = ['barrel_high_south', 'barrel_mid_south', 'barrel_eq_south',
                 'barrel_eq_north', 'barrel_mid_north', 'barrel_high_north']

    for mi in range(len(mlat_lims) - 1):
        m1 = mlat_lims[mi]
        m2 = mlat_lims[mi + 1]

        bne = barrel_df['barrel_ne'].values
        ml_lim = mlat[(mlat >= m1) & (mlat <= m2)]
        ne_lim = bne[(mlat >= m1) & (mlat <= m2)]

        ne_s_max = p_strings[mi] + 'ne_max'
        ne_s_min = p_strings[mi] + 'ne_min'
        ml_s_max = p_strings[mi] + 'mlat_max'
        ml_s_min = p_strings[mi] + 'mlat_min'

        nan_check = len(ne_lim[np.isnan(ne_lim)])

        if (len(ne_lim) > 0) & (nan_check != len(ne_lim)):
            stats_df[ne_s_max] = [np.nanmax(ne_lim)] * plen
            stats_df[ne_s_min] = [np.nanmin(ne_lim)] * plen
            stats_df[ml_s_max] = [ml_lim[np.nanargmax(ne_lim)]] * plen
            stats_df[ml_s_min] = [ml_lim[np.nanargmin(ne_lim)]] * plen
        else:
            stats_df[ne_s_max] = np.nan
            stats_df[ne_s_min] = np.nan
            stats_df[ml_s_max] = np.nan
            stats_df[ml_s_min] = np.nan

    return stats_df


def pass_info(time, mlat, ne, plen, pass_id, eia_state='unknown'):
    """Save pass info in a dataframe.

    Parameters
    ----------
    time : array-like
        time array
    mlat : array-like
        magnetic latitudes
    ne : array-like
        densities
    plen : int
        length of final pass_df
    pass_id : str
        pass_id created in build.find_cand or swarm_detections.swarm_detect
    eia_state : string
        eia category/orientation string from pyValEIA
        default : 'unknown'

    Returns
    -------
    pass_df : pandas DataFrame
        contains pass_id, starting and ending times of pass,
        starting and ending magnetic latitudes of pass,
        max and minimum densities and what magnetic lats they occurred at,
        time when pass crossed 0 degrees magnetic latitude,
        and eia_state
    """
    # save pass info including num, mlat vals, and EIA_type
    pass_df = pd.DataFrame()
    pass_df['pass_id'] = [pass_id] * plen
    pass_df['pass_time1'] = [pd.Timestamp(time[0]).to_pydatetime()] * plen
    pass_df['pass_time2'] = [pd.Timestamp(time[-1]).to_pydatetime()] * plen
    pass_df['pass_mlat1'] = [mlat[0]] * plen
    pass_df['pass_mlat2'] = [mlat[-1]] * plen

    # save maximum and minimum plasma densities
    pass_df['pass_max_np'] = [np.nanmax(ne)] * plen
    pass_df['pass_min_np'] = [np.nanmin(ne)] * plen
    pass_df['pass_mlat_max_np'] = [mlat[np.nanargmax(ne)]] * plen
    pass_df['pass_mlat_min_np'] = [mlat[np.nanargmin(ne)]] * plen

    # save when equator was 0
    z_id = np.nanargmin(abs(mlat))
    pass_df['pass_eqtime'] = [pd.Timestamp(time[z_id]).to_pydatetime()] * plen

    # save eia state
    pass_df['eia_state'] = [eia_state] * plen

    return pass_df
