"""Build Monthly files for EPI detection from Swarm.
Created: 24 February 2026
alanahco

List Functions
--------------
swarm_detect_daily
open_dayb4_check
rewrite_files
efi_file_only
"""
import datetime as dt
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from pySpyEPI.barrel import barrel_utils as bar_uts
from pySpyEPI.detection.stats import swarm_stats
from pySpyEPI.detection import build
from pySpyEPI.io.load import swarm_load as load
from pySpyEPI.detection.stats import eia_stats
from pySpyEPI.plotting import irr_plot
from pySpyEPI.io.write import write
from pySpyEPI.io.load import det_load


def swarm_detect_daily(date_array, satellite, fdir, file_dir, fig_dir=None,
                       remake=False, lt_skip=None,
                       equator_bound=35, auroral_bound=70, set_lat=45,
                       barrel_start=8, det_filt=5, peak_width=6,
                       fig_on=True, fft_test=False, nest_test=True,
                       prom_test=True, percent_test=True, wdr_test=True,
                       min_height=10**4, min_perc=10, wdr50_set=0.995,
                       wdr20_set=0.95,
                       num_barrels=5, prom_height=5000, e_up_large=0.1,
                       e_lo_large=0.05, e_up_small=0.1, e_lo_small=0.05,
                       upper_weight=2, lower_weight=1, scale_n=True,
                       exp_inc=True, ex_scale=2, filter_barrel=True,
                       svg_window=5, svg_poly=2, eia_mlat=30):
    """Create monthly detection files for bubble detection.

    Parameters
    ----------
    date_array : array-like
        date range for files and plots, one month at a time
        pandas date range
    satellite : str
        Swarm Satellite for files and plots
    fdir : string
        file directory for Swarm files
    fildir : string
        directory for saving detection files
    figdir : string or NoneType
        directory for saving detection files
    equator_bound : float
        equatorward magnetic latitude for midlat trough detection
        default 35 degrees maglat (+/-)
    auroral_bound : float
        auroralward magnetic latitude for midlat trough detection
        default 70 degrees maglat (+/-)
    set_lat : float
        highest magnetic latitude acceptable for trough detection
        default 50 degrees maglat (+/-)
    fig_on : boolean
        If True (default), figures will be made
    remake : boolean
        If True, file will be remade even if they already exist
        if False, files will not be remade unless they do not exist
    lt_skip : NoneType or array-like
        if specified, local times between lt_skip[0] and lt_skip[1] will not
        be recorded
    equator_bound : float kwarg (35)
        equatorward magnetic latitude to search down to
        set as a positive number, will be used as +/-equator_bound
    auroral_bound : float kwarg (75)
        autroalward magnetic latitude to search up to
        set as a positive number, will be used as +/-auroral_bound
    set_lat : float kwarg (55)
        if a trough cannot be found, set_lat will be used in its place
        set as a positive number, will be used as +/-set_lat
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
    obs : string
        observation satellite ('swarm', 'dmsp')
        default 'sat'
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
    Np : string
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
    time_str, ne_str, glat_str, glon_str, mlat_str, mlon_str, alt_str,
        lshell_st, lt_str : strings
        column names from data_df for time, density, geographic lat and lon,
        magnetic lat and lon, altitude, lshell, and lt_str
    eia_mlat : float
        maximum search magnetic latitude range for EIA type.
    Returns
    -------
    epi_all : list of DataFrame
        containing epi information for whole date_array
    eia_all : list of DataFrame
        containing eia information for whole date_array

    Notes
    -----
    Also saves a file and figure (if specified)
    Only makes EPI files if it does not already exists
    If EPI file exists, IRR and EFI files will not be made
    If EIA file does not exist, will be remade
    """

    if fig_dir is None:
        fig_on = False
        print('Figure will not be created or saved')

    # create arrays
    epi_all = []
    eia_all = []

    obs = 'swarm'
    # go through each day in month
    for day_st in date_array:

        # create empty daily arrays
        epi_day = []
        eia_day = []

        # use start and end days
        day_ed = day_st + dt.timedelta(days=1)

        # check if files need to be created
        if remake:
            re_epi = True
            re_eia = True
        else:
            re_epi, re_eia = rewrite_files(day_st, obs, file_dir, satellite)

        # if the files already exists, skip to next day
        if not re_epi and not re_eia:
            continue

        # open only IRR first
        data, IRR_bool, EFI_bool = load.load_all(
            day_st, day_ed, satellite, fdir, EFI=False)

        # if IRR_bool is false, skip the run
        if IRR_bool:
            IRR = data['IRR']
        else:
            continue

        # see if day before also needs to be opened based on mlat
        IRR_update, EFI_update, EFI_bool = open_dayb4_check(
            IRR, fdir, satellite, mlat_thresh=auroral_bound)

        # Limit by auroral upper bound
        IRR_lat = IRR_update[abs(IRR_update['mlat']) < auroral_bound]

        # use find gaps to get mlat segments
        sw_idx = IRR_lat.index.values
        lat_gaps = bar_uts.find_all_gaps(sw_idx)

        # iterate through the gaps
        for mgap in range(len(lat_gaps) - 1):

            # get each segment
            IRR_auroral = IRR_lat.iloc[lat_gaps[mgap]:lat_gaps[mgap + 1]]

            # make sure there is enough data to analyze
            # also ensures that we are not double analyzing
            min_mlat = IRR_auroral['mlat'].min()
            max_mlat = IRR_auroral['mlat'].max()
            if ((min_mlat > -30) | (max_mlat < 30)):
                continue

            zero_id = abs(IRR_auroral['mlat']).argmin()
            time_id = IRR_auroral['time'].iloc[zero_id].strftime(
                '%Y%m%d_%H%M')
            pass_id = f'{time_id}_{obs}{satellite}'

            # skip irrelevant local time hours
            if lt_skip is not None:
                lt_id = IRR_auroral['lt_hour'].iloc[zero_id]
                if (lt_id >= lt_skip[0]) & (lt_id <= lt_skip[1]):
                    continue

            time = IRR_auroral['time'].values
            time = time.astype('datetime64[us]')
            ne = IRR_auroral['ne'].values
            lat = IRR_auroral['lat'].values
            lon = IRR_auroral['lon'].values
            mlat = IRR_auroral['mlat'].values
            mlon = IRR_auroral['mlon'].values
            lt = IRR_auroral['lt'].values
            lt = lt.astype('datetime64[us]')
            alt = IRR_auroral['altitude'].values
            ne_flag = IRR_auroral['ne_flag'].values

            # Make sure that the majority of points are not Nans
            ne_check = ne[abs(mlat) <= set_lat]
            nans = np.isnan(ne_check)
            nperc = len(ne_check[nans]) / len(ne_check) * 100
            if nperc > 50:
                continue

            # eia detect
            print(pass_id)
            eia_df, eia_state, plats = eia_stats.eia_info(
                pass_id, time, lat, lon, mlat, mlon, lt, alt, ne,
                "Ne", mlat_val=eia_mlat, filt='barrel_average',
                interpolate=1, barrel_envelope=True, envelope_lower=0.6,
                envelope_upper=0.2, barrel_radius=3, window_lat=2)

            # if EPI file does not need to be remade, don't
            # if EPI file does need to be remade, EIA still has to be detected
            if re_epi:
                # get barrel, trough, and peak values
                (barrel_df, cand_df, peaks, properties,
                 pass_id, trough_lats) = build.find_cand(
                     IRR_auroral, obs, satellite, ne_flag=ne_flag, trough=True,
                     equator_bound=equator_bound, auroral_bound=auroral_bound,
                     set_lat=set_lat, barrel_start=barrel_start,
                     det_filt=det_filt, peak_width=peak_width, freq=1,
                     fft_test=fft_test, nest_test=nest_test,
                     prom_test=prom_test, percent_test=percent_test,
                     wdr_test=wdr_test, min_height=min_height,
                     min_perc=min_perc, wdr50_set=wdr50_set,
                     wdr20_set=wdr20_set, eia_state=eia_state, pass_id=pass_id,
                     num_barrels=num_barrels, prom_height=prom_height,
                     e_up_large=e_up_large, e_lo_large=e_lo_large,
                     e_up_small=e_up_small, e_lo_small=e_lo_small,
                     upper_weight=upper_weight, lower_weight=lower_weight,
                     scale_n=scale_n, exp_inc=exp_inc, ex_scale=ex_scale,
                     filter_barrel=filter_barrel, svg_window=svg_window,
                     svg_poly=svg_poly)

                # create figure
                if fig_on:

                    # limit by trough lats
                    trough_mask = ((IRR_auroral['mlat'] >= trough_lats[1])
                                   & (IRR_auroral['mlat'] <= trough_lats[0]))

                    irr_use = IRR_auroral[trough_mask]

                    if EFI_bool:
                        t1A = irr_use['time'].iloc[0]
                        t2A = irr_use['time'].iloc[-1]

                        # limit EFI by time (EFI update is fine)
                        efi_use = EFI_update[((EFI_update['time'] > t1A)
                                              & (EFI_update['time'] < t2A))]
                    else:
                        efi_use = None

                    # create filename for figure
                    date_dir, figname = write.build_figname(
                        day_st, pass_id, fig_dir, obs, satellite)
                    fig_file = os.path.join(date_dir, figname)

                    # check if figure already exists
                    if not os.path.exists(fig_file):
                        fig = irr_plot.ne_2barrel(
                            irr_use, satellite, barrel_df, peaks,
                            properties, cand_df, EFI_data=efi_use,
                            blo_bool=True, hspace=0.4, fig_size=None,
                            leg_loc='center', fs=14, flag_mark=True,
                            eia_state=eia_state, eia_lats=plats)

                        write.write_fig(day_st, fig, pass_id, obs,
                                        satellite, fig_dir, date_dir=date_dir,
                                        figname=figname)
                        plt.close(fig)

                # save info in daily files
                if not cand_df.empty:
                    epi_day.append(cand_df)
            if not eia_df.empty:
                eia_day.append(eia_df)

        if re_epi:
            if len(epi_day) != 0:
                epi_day_df = pd.concat(epi_day, ignore_index=True)

                # get swarm specific statistics DataFrames
                irr_day = swarm_stats.IRR_stats(IRR_update, epi_day_df)

                # if EFI bool exists
                if EFI_bool:
                    efi_day = swarm_stats.EFI_stats(EFI_update, epi_day_df)
                    write.write_stats(efi_day, day_st, obs, file_dir, 'efi',
                                      satellite, daily=True)

                # write daily files
                write.write_stats(epi_day_df, day_st, obs, file_dir,
                                  'detection', satellite, daily=True)
                write.write_stats(irr_day, day_st, obs, file_dir, 'irr',
                                  satellite, daily=True)

                # save DataFrames
                if not epi_day_df.empty:
                    epi_all.append(epi_day_df)

        if re_eia:
            # convert epi day to dataframe
            if len(eia_day) != 0:
                eia_day_df = pd.concat(eia_day, ignore_index=True)
                write.write_stats(eia_day_df, day_st, obs, file_dir, 'eia',
                                  satellite, daily=True)
                # save DataFrame
                if not eia_day_df.empty:
                    eia_all.append(eia_day_df)

    return epi_all, eia_all


def open_dayb4_check(IRR, fdir, satellite, mlat_thresh=70, IRR_only=False):
    """Open file of day before and incldue it in original dataframes

    Parameters
    ----------
    IRR : pd.DataFrame
        dataframe of IRR data
    fdir : str
        file directory
    satelltie: string
        indicates which satellite to use
    mlat_thresh : float
        magnetic latitude threshold indicating what values to concatenate
        default 70 degrees mlat
    IRR_only : bool
        False default will open all files
        True will open only IRR file
    Returns
    -------
    IRR_update : pd.DataFrame
        IRR with previous day data included up to mlat_thresh
    EFI_update : pd.DataFrame
        EFI with previous day data included up to mlat_thresh based on IRR time
    """
    # open previous day as well if necessary
    if abs(IRR['mlat'].iloc[0]) < mlat_thresh:
        # get day before
        st = IRR['time'].iloc[0].replace(hour=0, minute=0)
        st_b4 = st - dt.timedelta(days=1)
        ed_b4 = st_b4 + dt.timedelta(days=1)

        # open day before of IRR only
        datab4, IRR_boolb4, EFI_boolb4 = load.load_all(
            st_b4, ed_b4, satellite, fdir, EFI=False)

        # if IRR_bool is false, skip the run
        if IRR_boolb4:
            IRRb4 = datab4['IRR']
        else:
            IRR_update = IRR

        # if there is data present on the day before
        if IRR_boolb4:
            IRR_lat = IRRb4[abs(IRRb4['mlat']) < mlat_thresh]
            b4_gaps = bar_uts.find_all_gaps(IRR_lat.index.values)

            # limit IRR b4 to just b
            IRR_lim = IRR_lat.iloc[b4_gaps[-1]:-1]

            # update IRR
            IRR_update = pd.concat([IRR_lim, IRR], ignore_index=True)

    else:
        IRR_update = IRR

    # Open MAG and EFI files using time from IRR_update

    if IRR_only:
        return IRR_update
    else:
        st_open = IRR_update["time"].iloc[0]
        ed_open = IRR_update["time"].iloc[-1]

        data_update, IRR_bool, EFI_bool = load.load_all(
            st_open, ed_open, satellite, fdir, IRR=False, EFI=True)

        # Get EFI data, make empty array if false
        if EFI_bool:
            EFI_update = data_update['EFI']
        else:
            EFI_update = pd.DataFrame()

        return IRR_update, EFI_update, EFI_bool


def rewrite_files(stime, obs, file_dir, sat, daily=True):
    """See if file exists or needs to be created.

    Parameters
    ----------
    stime : datetime
        day of desired file
    obs : str
        Name of data set e.g., 'SWARM', 'MADRIGAL')
    file_dir : str
        File directory
    sat : string
        sat number or letter
    daily : boolean
        creates name for daily files if True
        creates name for monthly files if false
        default True
    Returns
    -------
    re_epi : bool
        True if file does not exist
        False if file already exists
    re_eia : bool
        True if file does not exist
        False if file already exists
    """
    # EPI FILE
    date_dir, fname = write.build_stats_filename(stime, obs, file_dir,
                                                 'detection', sat, daily=daily)
    save_file = os.path.join(date_dir, fname)

    if not os.path.exists(save_file):
        re_epi = True
    else:
        re_epi = False

    # EIA FILE
    date_dir, fname = write.build_stats_filename(stime, obs, file_dir,
                                                 'eia', sat, daily=daily)
    save_file = os.path.join(date_dir, fname)

    if not os.path.exists(save_file):
        re_eia = True
    else:
        re_eia = False

    return re_epi, re_eia


def efi_file_only(date_array, sat, fdir, epi_dir):
    """Create EFI detection files only.
    Parameters
    ----------
    date_array : array-like
        date range for files and plots, one month at a time
        pandas date range
    sat : str
        Swarm Satellite for files and plots
    fdir : string
        file directory for Swarm files
    epi_dir : string
        directory for saving detection files and opening detection files

    Returns
    -------
    creates daily efi files
    """

    # iterate through dates
    for dayt in date_array:
        d_str = dayt.strftime("%Y-%m-%d")
        date_range = pd.date_range(start=d_str, end=d_str)
        epi_df = det_load.load_detections(date_range, 'swarm', epi_dir,
                                          'detection', sat)

        if len(epi_df) == 0:
            continue

        # load efi file
        st = epi_df['edge1_time'].iloc[0] - dt.timedelta(minutes=5)
        ed = epi_df['edge2_time'].iloc[-1] + dt.timedelta(minutes=5)

        data, IRR_bool, EFI_bool = load.load_all(
            st, ed, sat, fdir, IRR=False, EFI=True)

        # if EFI file exists, use it, else continue
        if EFI_bool:
            EFI_data = data['EFI']
        else:
            continue

        # pass_id = epi_df['pass_id'].iloc[0]
        efi_df = swarm_stats.EFI_stats(EFI_data, epi_df)
        write.write_stats(efi_df, dayt, 'swarm', epi_dir, 'efi', sat,
                          daily=True)
