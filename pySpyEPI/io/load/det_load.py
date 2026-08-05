"""Load detection files.
Updated 14 March 2026
alanahco

List Functions
--------------
load_detections
epi_only
lt_mask
epi_te
omni_bubble
pass_maxmin
update_eianorm
update_eianormAC
"""
import os
import pandas as pd
from pySpyEPI.io.write import write
import numpy as np


def load_detections(date_array, obs, file_dir, filetype, sat, daily=True):
    """Load detection files from daily save files.

    Parameters
    ----------
    date_array : array-like
        date range for files and plots, one month at a time
        pandas date range
    obs : str
        Name of data set e.g., 'SWARM', 'MADRIGAL')
    file_dir : str
        File directory
    filetype : str
        detection, irr, mag, efi, tii, eia
    sat : string
        sat number or letter
    daily : boolean
        creates name for daily files if True
        creates name for yearly files if false
        default True

    Returns
    -------
    df_full : DataFrame
        contains all daily file info for specified filetype for specified
        date range
    """

    if filetype == 'detection':
        time_strs = ['center_time', 'center_lt', 'edge1_time', 'edge1_lt',
                     'edge2_time', 'edge2_lt', 'pass_eqtime']
        lt_str = 'center_lt'
    elif filetype == 'eia':
        time_strs = ['time_start', 'time_end', 'lt']
        lt_str = 'lt'
    elif filetype == 'irr':
        time_strs = ['center_time']
        lt_str = None
    elif filetype == 'efi':
        time_strs = ['center_time']
        lt_str = None
    elif filetype == 'solarwind':
        time_strs = ['time_start', 'time_end', 'time_min']
        lt_str = None
    elif filetype == 'omni':
        time_strs = ['center_time', 'time_sym_dip_b4', 'time_symh_min_b4',
                     'time_sym_rec_b4', 'time_sym_dip_af', 'time_symh_min_af',
                     'time_sym_rec_af', 'time_al_dip_b4', 'time_al_min_b4',
                     'time_al_rec_b4', 'time_al_dip_af', 'time_al_min_af',
                     'time_al_rec_af']
        lt_str = None

    df_list = []

    # iterate through date array
    for day_st in date_array:

        # get filename
        date_dir, fname = write.build_stats_filename(day_st, obs, file_dir,
                                                     filetype, sat,
                                                     daily=daily)
        save_file = os.path.join(date_dir, fname)

        # Check if file exists
        if os.path.exists(save_file):

            # get dataframe
            df = pd.read_csv(save_file, sep='\t')
            df.columns = df.columns.str.lstrip('#')

            df = df.copy()

            for t_str in time_strs:
                df[t_str] = pd.to_datetime(df[t_str], format='ISO8601')

            if lt_str is not None:
                # add lt hour
                df["lt_hr"] = df[lt_str].apply(
                    lambda tim: tim.hour + tim.minute / 60 + tim.second / 3600)

            # add it to dataframe list
            df_list.append(df)

    # concat df_list to make one large dataframe
    if len(df_list) != 0:
        df_full = pd.concat(df_list, ignore_index=True)
    else:
        df_full = pd.DataFrame()

    return df_full


def epi_only(epi_df, lt_limit=None, mlat_limit=None):
    """Filter EPB events from detection files.

    Parameters
    ----------
    info_df : dataframe
        detection file as a dataframe made by open_epb_file function
    lt_limit : NoneType or array-like
        default None, local time limit will be set to 18 to 10 LT
    mlat_limit : None Type or int
        limit the data by +/- mlat limit
    wdr_thresh : float
        maximum value for to reintroduce some 9 flags
        default 0.95
    Returns
    -------
    bubs : dataframe
        updated dataframe with non-epbs removed, limited by local time,
        month added, and unique id added, if te is True,
        limited by Te availability and quality as well
    """
    bubs = epi_df.copy()

    # Update WDR
    wdr_mask = ((bubs['epi_flag'] == 9) & (bubs['med_epi_l'] <= 0.995)
                & (bubs['med_epi_r'] <= 0.995)
                & (bubs['p20_epi_r'] < 0.95)
                & (bubs['p20_epi_l'] < 0.95)
                & (bubs['bub_in_bub'] == 0))
    bubs.loc[wdr_mask, 'epi_flag'] = 1

    # bubble flags
    bubs = bubs[bubs['epi_flag'] == 1]

    mlat_mask = ((bubs['mlat_flag'] < 4) | (bubs['mlat_flag'] == 4.1)
                 | (bubs['mlat_flag'] == 5.1))
    bubs = bubs[mlat_mask]

    # create lt_hr column
    # limit by local times
    mask = lt_mask(bubs, lt_str="lt_hr", lt_limit=lt_limit)
    bubs = bubs[mask]

    # limit by magnetic latitude
    if mlat_limit is not None:
        bubs = bubs[abs(bubs['center_mlat']) <= abs(mlat_limit)]

    # add a month val for plotting purposes
    bubs['month'] = bubs['center_time'].dt.month

    ct = bubs['center_time']

    dm = ct.dt.days_in_month

    bubs['decimal_month'] = (ct.dt.month + (ct.dt.day - 1 + ct.dt.hour / 24
                                            + ct.dt.minute / 1440) / dm)

    return bubs


def lt_mask(dat, lt_str=None, lt_limit=None, max_inclusive=True,
            min_inclusive=True):
    """Create local time mask.

    Parameters
    ----------
    dat : dataframe or array
        dataframe with local time as a param or an array of local times
    lt_str : NoneType or string
        if None, assumed to be array
        if string, assumed to be dataframe with lt_str as local time column
    lt_limit : NoneType or array-like
        default None, local time limit will be set to 18 to 10 LT
    max_inclusive, min_inclusive : boolean
        if True, mask will be applied as >= and <= instead of > or <
        default True
    Returns
    -------
    mask : boolean mask
        representing data within lt limit
    """

    # set limit if none is provided
    if lt_limit is None:
        lt_limit = [18, 10]

    # get min and max lt
    min_lt = min(lt_limit)
    max_lt = max(lt_limit)

    # get only lt data
    if lt_str is None:
        lt_dat = dat
    else:
        lt_dat = dat[lt_str].values

    # create mask
    if (min_lt <= 12) & (max_lt >= 12):
        # e.g. 23 and 6
        if min_inclusive:
            min_mask = (lt_dat >= max_lt)
        else:
            min_mask = (lt_dat > max_lt)
        if max_inclusive:
            max_mask = (lt_dat <= min_lt)
        else:
            max_mask = (lt_dat < min_lt)
        mask = max_mask | min_mask
    else:
        # e.g. 0 and 6 or 21 and 23
        if min_inclusive:
            min_mask = (lt_dat >= min_lt)
        else:
            min_mask = (lt_dat > min_lt)
        if max_inclusive:
            max_mask = (lt_dat <= max_lt)
        else:
            max_mask = (lt_dat < max_lt)
        mask = max_mask & min_mask

    return mask


def epi_te(bub_df, efi_df, perc_nan=1, nan_test=False, jump_test=True):
    """Create a DataFrame for only good temperature bubbles from efi_df.

    Parameters
    ----------
    bub_df : DataFrame
        DataFrame of confirmed bubbles from det_load.epi_only()
    efi_df : DataFrame
        loaded DataFrame of all efi data from det_load.load_detections()
    perc_nan : float
        percent of nan values acceptable for temperature data
        default below 1%
    nan_test : boolean
        For higher perc_nan, it is good to do a nan_test to make sure that
        the the edges and center te are not nan
        default is False
    jump_test : boolean
        Remove Jump Te = 1 if True (default)
    Returns
    -------
    efi_bub : DataFrame
        DataFrame of bubble info with only good temperature data
        also includes lt_hr and te_dif_err, which is abs(edge1 te - edge2 te)/2
    """
    # get indices of bubbles by center time
    bub_idx = []
    lt_hr = []
    lons = []
    mlats = []
    perc_depths = []
    lats = []
    eq_time = []
    proms = []
    edge1 = []
    edge2 = []
    mlons = []

    for ct in bub_df['center_time'].values:
        if len(np.where(efi_df['center_time'] == ct)[0]) > 0:
            bub_idx.append(np.where(efi_df['center_time'] == ct)[0][0])
            lt_hr.append(
                bub_df[bub_df['center_time'] == ct]['lt_hr'].values[0])
            lons.append(
                bub_df[bub_df['center_time'] == ct]['center_lon'].values[0])
            mlats.append(
                bub_df[bub_df['center_time'] == ct]['center_mlat'].values[0])
            perc_depths.append(
                bub_df[bub_df['center_time'] == ct]['percent_depth'].values[0])
            lats.append(
                bub_df[bub_df['center_time'] == ct]['center_lat'].values[0])
            eq_time.append(
                bub_df[bub_df['center_time'] == ct]['pass_eqtime'].values[0])
            proms.append(
                bub_df[bub_df['center_time'] == ct]['center_prom'].values[0])
            edge1.append(
                bub_df[bub_df['center_time'] == ct]['edge1_time'].values[0])
            edge2.append(
                bub_df[bub_df['center_time'] == ct]['edge2_time'].values[0])
            mlons.append(
                bub_df[bub_df['center_time'] == ct]['center_mlon'].values[0])

    # only pick bubbles
    efi_bub = efi_df.copy()
    efi_bub = efi_bub.iloc[bub_idx]
    efi_bub['lt_hr'] = lt_hr
    efi_bub['center_mlat'] = mlats
    efi_bub['center_mlon'] = mlons
    efi_bub['center_lon'] = lons
    efi_bub['center_lat'] = lats
    efi_bub['percent_depth'] = perc_depths
    efi_bub['pass_eqtime'] = eq_time
    efi_bub['edge1_time'] = edge1
    efi_bub['edge2_time'] = edge2
    efi_bub['center_prom'] = proms

    # filter temperatures
    efi_bub = efi_bub[efi_bub['lp_flag'] == 1]
    if jump_test:
        efi_bub = efi_bub[efi_bub['te_jump'] == 0]
    efi_bub = efi_bub[efi_bub['te_percnan'] < perc_nan]
    efi_bub = efi_bub[efi_bub['time_lpsweep'] != 0]

    # nan test if requested
    if nan_test:
        mask_nan = ((efi_bub['edge1_te'].notna())
                    & (efi_bub['edge2_te'].notna())
                    & (efi_bub['center_te'].notna()))
        efi_bub = efi_bub[mask_nan]

    # create te error using difference between
    te_dif_err = abs(efi_bub['edge1_te'] - efi_bub['edge2_te']) / 2

    efi_bub['te_dif_err'] = te_dif_err

    return efi_bub


def omni_bubble(bub_df, omni_df):
    """Create a DataFrame for only good temperature bubbles from efi_df.

    Parameters
    ----------
    bub_df : DataFrame
        DataFrame of confirmed bubbles from det_load.epi_only()
    omni_df : DataFrame
        loaded DataFrame of all omni data from det_load.load_detections()

    Returns
    -------
    omni_bub : DataFrame
        DataFrame of bubble info with additional params:
        'hr_sym_dip_b4', 'hr_sym_min_b4', 'hr_sym_rec_b4', 'hr_sym_dip_af',
        'hr_sym_min_af', 'hr_sym_rec_af', 'hr_al_dip_b4', 'hr_al_min_b4',
        'hr_al_rec_b4', 'hr_al_dip_af', 'hr_al_min_af', 'hr_al_rec_af'
    """
    # get indices of bubbles by center time
    bub_idx = []
    lt_hr = []
    lons = []
    mlats = []
    perc_depths = []
    id = []

    for ct in bub_df['center_time'].values:

        if len(np.where(omni_df['center_time'] == ct)[0]) > 0:
            bub_idx.append(np.where(omni_df['center_time'] == ct)[0][0])
            lt_hr.append(
                bub_df[bub_df['center_time'] == ct]['lt_hr'].values[0])
            lons.append(
                bub_df[bub_df['center_time'] == ct]['center_lon'].values[0])
            mlats.append(
                bub_df[bub_df['center_time'] == ct]['center_mlat'].values[0])
            perc_depths.append(
                bub_df[bub_df['center_time'] == ct]['percent_depth'].values[0])
            id.append(
                bub_df[bub_df['center_time'] == ct]['pass_id'].values[0])

    # bubble params
    omni_bub = omni_df.copy()
    omni_bub = omni_bub.iloc[bub_idx]
    omni_bub['lt_hr'] = lt_hr
    omni_bub['mlat'] = mlats
    omni_bub['lon'] = lons
    omni_bub['percent_depth'] = perc_depths
    omni_bub['pass_id'] = id

    # Add params
    # SYM H
    omni_bub['hr_sym_dip_b4'] = (
        (omni_bub['center_time'] - omni_bub['time_sym_dip_b4'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_sym_min_b4'] = (
        (omni_bub['center_time'] - omni_bub['time_symh_min_b4'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_sym_rec_b4'] = (
        (omni_bub['center_time'] - omni_bub['time_sym_rec_b4'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_sym_dip_af'] = (
        (omni_bub['center_time'] - omni_bub['time_sym_dip_af'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_sym_min_af'] = (
        (omni_bub['center_time'] - omni_bub['time_symh_min_af'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_sym_rec_af'] = (
        (omni_bub['center_time'] - omni_bub['time_sym_rec_af'])
        / pd.Timedelta('1 hour'))

    # if recover time before is < 0, make it equal to 0
    omni_bub.loc[omni_bub['hr_sym_rec_b4'] < 0, 'hr_sym_rec_b4'] = 0

    # AL INDEX
    omni_bub['hr_al_dip_b4'] = (
        (omni_bub['center_time'] - omni_bub['time_al_dip_b4'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_al_min_b4'] = (
        (omni_bub['center_time'] - omni_bub['time_al_min_b4'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_al_rec_b4'] = (
        (omni_bub['center_time'] - omni_bub['time_al_rec_b4'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_al_dip_af'] = (
        (omni_bub['center_time'] - omni_bub['time_al_dip_af'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_al_min_af'] = (
        (omni_bub['center_time'] - omni_bub['time_al_min_af'])
        / pd.Timedelta('1 hour'))
    omni_bub['hr_al_rec_af'] = (
        (omni_bub['center_time'] - omni_bub['time_al_rec_af'])
        / pd.Timedelta('1 hour'))

    # if recover time before is < 0, make it equal to 0
    omni_bub.loc[omni_bub['hr_al_rec_b4'] < 0, 'hr_al_rec_b4'] = 0

    return omni_bub


def pass_maxmin(bub, param='mlat', min_max='max', absolute=True):
    """Keep only pass max or min based on param.
    bub : DataFrame
        bubble info with mlat as a column. Any det_load with mlat as a param
    param : str
        parameter for limiting bubble DF
        default 'mlat'
    min_max : str
        determines if min or max will be kept
        default 'max'
    absolute : bool
        if True, absolute value of min/max will be used
        if False, true min/max will be used
        default True
    Returns
    -------
    bub_pass : DataFrame
        limited bubble dataframe with only param max or min

    Notes
    -----
    DataFrame must contain pass_id
    """
    # save pass id and unique id
    uq_id = bub['pass_id'].unique()

    bub_list = []
    # iterate through Pass IDs
    for id in uq_id:

        # limit by ID
        bub_lim = bub.copy()
        bub_lim = bub_lim[bub_lim['pass_id'] == id]

        # get max/min index and value
        i_min = bub_lim[param].argmin()
        v_min = bub_lim[param].min()
        i_max = bub_lim[param].argmin()
        v_max = bub_lim[param].max()

        # get index to save
        # use absolute value if absolute is True
        if absolute:
            if min_max == 'min':
                if abs(v_min) < abs(v_max):
                    im = i_min
                else:
                    im = i_max
            elif min_max == 'max':
                if abs(v_min) > abs(v_max):
                    im = i_min
                else:
                    im = i_max
        else:
            if min_max == 'min':
                im = i_min
            else:
                im = i_max

        # append to bub_list
        bub_list.append(bub_lim.iloc[[im]])

    # concatenate bub_list
    bub_pass = pd.concat(bub_list, ignore_index=True)

    return bub_pass


def efi_combine_AC(efi_bubA, efi_bubC, t_dif=1):
    """Combine A and C satellites for EFI values.

    Parameters
    ----------
    efi_bubA, efi_bubC : pd.DataFrames
        created by det_load.epi_te
        A and C satellite EFI info
    t_dif : float
        comparison time difference in minutes between pass eqautorial times
        of A and C

    Returns
    -------
    efi_bubAC : pd.DataFrame
        combined A and C data with only best data kept, ordered by center_time
    AC_check : pd.DataFrame
        info about how satellites were chosen including A and C center times,
        mean nan percents, lengths, and mean te_dif_err. Also included is
        'sat_keep' for which satelltie was kept.
    Notes
    -----
    First, find matching equatorial passes for Swarm A and C
    Second, calculate the average te_nan_perc, average te_dif_err, and
        array length corresponding to each Swarm A and C pass
    Third, compare values
        Start by comparing nan_perc, keep pass with smallest averge nan_perc
        if  they are equal, compare length
            keep pass with most EPIs
        if they are the same length,
            keep the pass with the smallest average te_dif_err
        if they are the same (unlikely),
            keep A
    """
    utime_A = efi_bubA['pass_eqtime'].unique()
    utime_C = efi_bubC['pass_eqtime'].unique()

    new_tC = []
    new_tA = []
    dif_A = []
    dif_C = []
    nan_A = []
    nan_C = []
    len_C = []
    len_A = []

    i_remA = []
    i_remC = []
    sat_use = []

    for uA in utime_A:
        td = abs(utime_C - uA)

        if td.min() < pd.Timedelta(minutes=t_dif):
            uC = utime_C[td.argmin()]
            new_tC.append(uC)
            new_tA.append(uA)

            # nan percent and array length Swarm C ----------------------------
            nan_listC = efi_bubC['te_percnan'][efi_bubC['pass_eqtime'] == uC]
            nC = np.nanmean(nan_listC)
            nan_C.append(np.nanmean(nC))
            lC = len(nan_listC)
            len_C.append(lC)

            # nan percent and array length Swarm A ----------------------------
            nan_listA = efi_bubA['te_percnan'][efi_bubA['pass_eqtime'] == uA]
            nA = np.nanmean(nan_listA)
            nan_A.append(nA)
            lA = len(nan_listC)
            len_A.append(lA)

            # Te Dif Error Swarm C --------------------------------------------
            dif_listC = efi_bubC['te_dif_err'][efi_bubC['pass_eqtime'] == uC]
            dC = np.nanmean(dif_listC)
            dif_C.append(dC)

            # Te Dif Error Swarm A --------------------------------------------
            dif_listA = efi_bubA['te_dif_err'][efi_bubA['pass_eqtime'] == uA]
            dA = np.nanmean(dif_listA)
            dif_A.append(dA)

            idA = efi_bubA[efi_bubA['pass_eqtime'] == uA].index.tolist()
            idC = efi_bubC[efi_bubC['pass_eqtime'] == uC].index.tolist()

            # check percent nan first
            if nC < nA:
                # keep
                sat_use.append('C')

                i_remA.append(idA)
            elif nA < nC:
                # keep
                sat_use.append('A')

                # remove
                i_remC.append(idC)
            else:
                # if they are equal, do more testing!
                # keep one with more EPIs, otherwise keep lowest te_dif_err
                if lC > lA:
                    # keep
                    sat_use.append('C')

                    # remove
                    i_remA.append(idA)

                elif lA > lC:
                    # keep
                    sat_use.append('A')

                    # remove
                    i_remC.append(idC)

                else:
                    # if both arrays are the same length, check te_dif_err
                    if dC < dA:
                        # keep
                        sat_use.append('C')

                        # remove
                        i_remA.append(idA)
                    elif dA < dC:
                        # keep
                        sat_use.append('A')

                        # remove
                        i_remC.append(idC)
                    else:
                        # if te_dif_err are equal, (unlikely) use A
                        sat_use.append('A')

                        # remove
                        i_remC.append(idC)

    AC_check = pd.DataFrame()
    AC_check['timeA'] = new_tA
    AC_check['timeC'] = new_tC

    AC_check['nanA'] = nan_A
    AC_check['nanC'] = nan_C

    AC_check['lenA'] = len_A
    AC_check['lenC'] = len_C

    AC_check['errA'] = dif_A
    AC_check['errC'] = dif_C

    AC_check['sat_keep'] = sat_use

    # create new clean dataframes
    efi_cleanC = efi_bubC.copy()
    efi_cleanC['sat'] = ['C'] * len(efi_cleanC)
    efi_cleanA = efi_bubA.copy()
    efi_cleanA['sat'] = ['A'] * len(efi_cleanA)

    if len(i_remA) > 0:
        i_remA = np.concatenate(i_remA)
        efi_cleanA = efi_cleanA.drop(index=i_remA)

    if len(i_remC) > 0:
        i_remC = np.concatenate(i_remC)
        efi_cleanC = efi_cleanC.drop(index=i_remC)

    efi_bubAC = [efi_cleanA, efi_cleanC]
    efi_bubAC = pd.concat(efi_bubAC, ignore_index=True)

    efi_bubAC = efi_bubAC.sort_values(by='center_time')

    return efi_bubAC, AC_check


def update_eianorm(eia_df, efi_df):
    """Update EIA DataFrame for normalizing local time counts.

    Parameters
    ----------
    eia_df : DataFrame
        eia dataframe from det_load.load_detections filetype='eia'
    efi_df : DataFrame
        Bubble dataframe from det_load.load_detections filetype='efi'
        Limiting dataframe updated by det_load.epi_te
    Returns
    -------
    eia_return : DataFrame
        eia_df to include only the same local times as efi_df and only the same
        days
    """
    # Limit by local time based on bub_df
    # using largest value greater than 12 and largest value less than 12
    l1 = int(efi_df['lt_hr'][efi_df['lt_hr'] > 12].min())
    l2 = np.round(efi_df['lt_hr'][efi_df['lt_hr'] < 12].max())

    # mask EIA by same local times
    mask = lt_mask(eia_df, lt_str='lt_hr', lt_limit=[l1, l2])
    eia_LT = eia_df[mask]

    # Only want days that we have Te data available
    efi_day = efi_df['center_time'].dt.strftime('%d%m%Y')
    eia_day = eia_LT['time_start'].dt.strftime('%d%m%Y')

    # only look at unique days
    day_uq = efi_day.unique()

    # create new array for dataframes
    eia_new = []

    for d in day_uq:
        eia_update = eia_LT.copy()
        eia_update = eia_update[eia_day == d]

        eia_new.append(eia_update)

    eia_return = pd.concat(eia_new, ignore_index=True)

    return eia_return


def update_eianormAC(eia_dfA, eia_dfC, efi_dfAC):
    """Update EIA DataFrame for normalizing local time counts.

    Parameters
    ----------
    eia_dfA : DataFrame
        eia dataframe from det_load.load_detections filetype='eia' for sat A
    eia_dfC : DataFrame
        eia dataframe from det_load.load_detections filetype='eia' for sat C
    efi_dfAC : DataFrame
        Bubble dataframe from det_load.load_detections filetype='efi'
        updated by det_load.epi_te and det_load.efi_combineAC
    Returns
    -------
    eia_AC : boolean
        if True, efi_df is a combined A and C dataframe instead of only A
    """
    # Limit by local time based on bub_df
    # using largest value greater than 12 and largest value less than 12
    l1 = int(efi_dfAC['lt_hr'][efi_dfAC['lt_hr'] > 12].min())
    l2 = np.round(efi_dfAC['lt_hr'][efi_dfAC['lt_hr'] < 12].max())

    # mask EIA by same local times
    mask = lt_mask(eia_dfA, lt_str='lt_hr', lt_limit=[l1, l2])
    eia_LTA = eia_dfA[mask]
    mask = lt_mask(eia_dfC, lt_str='lt_hr', lt_limit=[l1, l2])
    eia_LTC = eia_dfC[mask]

    # only want days where Te data is available
    efi_day = efi_dfAC['center_time'].dt.strftime('%d%m%Y')
    eia_dayA = eia_LTA['time_start'].dt.strftime('%d%m%Y')
    eia_dayC = eia_LTC['time_start'].dt.strftime('%d%m%Y')

    # only look at unique days
    day_uq = efi_day.unique()
    # create new array for dataframes of mixed A and C eia_df
    eia_new = []

    for d in day_uq:
        # if eia info available for A, use it
        eia_updateA = eia_LTA.copy()
        eia_updateA = eia_updateA[eia_dayA == d]

        # if no EIA available for A, use C
        if len(eia_updateA) == 0:
            eia_updateC = eia_LTC.copy()
            eia_updateC = eia_updateC[eia_dayC == d]
            eia_new.append(eia_updateC)
        else:
            eia_new.append(eia_updateA)

    eia_AC = pd.concat(eia_new, ignore_index=True)

    return eia_AC
