"""Validate pySpyEPI by finding conjunctions with other data sets.
Created 02 July 2026
By Alanah Cardenas-O'Toole

Functions
---------
goldcon
reduce_conj
gold_by_lon
"""
from pySpyEPI.utils.sats import gold_utils as gu
import pandas as pd
import numpy as np
from pySpyEPI.barrel import barrel_utils
import datetime as dt


def goldcon(time, mlat, lon, max_mindif=60, max_mlat=40):
    """Find conjunction between any satellite and GOLD satellite.

    Parameters
    ----------
    time : array-like
        Universal Times as DateTimes associated with the satellite
        can also be numpy.datetime64 array
    mlat : array-like
        magnetic latitude of the satellite pass
    lon : array-like
        geographic longitude of the satellite pass
    max_mindif : float
        maximum minute difference between satellite time and GOLD time
        default 60
    max_mlat : float
        maximum magnetic latitude for search
        default 40 degrees
    Returns
    -------
    conj_df : pd.DataFrame
        conjunction starting time (sat_t1) and ending time (sat_t2), satellite
        starting (sat_lon1) and ending (sat_lon2) longitudes,
        gold time (gold_time), and gold starting (gold_lon1) and ending
        (gold_lon2) longitudes
    """

    # make an index array for sat info
    id = np.arange(len(time))

    # get gold time and lon ranges
    gtl = gu.time2lonRange()

    # convert max_min_dif to fractional hour
    mdif = max_mindif / 60

    # convert to datetime then fractional hour for GOLD
    time_filters = {dt.datetime.strptime(t, '%H:%M').time():
                    lon_range for t, lon_range in gtl.items()}
    gold_hr = []
    gl1 = []
    gl2 = []
    for time_key, lon_range in time_filters.items():
        gold_hr.append(time_key.hour + time_key.minute / 60)
        gl1.append(lon_range[0])
        gl2.append(lon_range[-1])

    # get fractional hour of satellite
    # try it as if it is pandas dataframe
    # otherwise treat it as a numpy.datetime64
    try:
        sat_hr = time.dt.hour + time.dt.minute / 60
        sat_hr = np.asarray(sat_hr)
    except AttributeError:
        time_arr = pd.to_datetime(time)
        sat_hr = time_arr.hour + time_arr.minute / 60

    # conver to array
    lon = np.asarray(lon)
    mlat = np.asarray(mlat)
    time = np.asarray(time)

    # set lists
    t1 = []
    t2 = []
    gtime = []
    glon1 = []
    glon2 = []
    sat_lon1 = []
    sat_lon2 = []

    # iterate through GOLD info and mask by longitude and time
    for gt, l1, l2 in zip(gold_hr, gl1, gl2):
        # mask by longitude and time
        hr_mask = (sat_hr >= (gt - mdif)) & (sat_hr <= (gt + mdif))
        if l1 <= l2:
            lon_mask = (lon >= l1) & (lon <= l2)
        else:
            lon_mask = (lon >= l1) | (lon <= l2)

        # Mask sat by mlat
        mlat_mask = (abs(mlat) <= max_mlat)

        mask = hr_mask & lon_mask & mlat_mask

        time_new = time[mask]
        lon_new = lon[mask]
        id_new = id[mask]

        if len(time_new) == 0:
            continue
        gaps = barrel_utils.find_all_gaps(id_new)

        # go through gap list and save data
        for i in range(len(gaps) - 1):
            gi = gaps[i]
            gj = gaps[i + 1]

            # satellite
            t1.append(time_new[gi])
            t2.append(time_new[gj])
            sat_lon1.append(lon_new[gi])
            sat_lon2.append(lon_new[gj])
            # gold
            glon1.append(l1)
            glon2.append(l2)
            gtime.append(gt)

    # save as dataframe
    conj_df = pd.DataFrame()
    conj_df['sat_t1'] = t1
    conj_df['sat_t2'] = t2
    conj_df['sat_lon1'] = sat_lon1
    conj_df['sat_lon2'] = sat_lon2
    conj_df['gold_t'] = gtime
    conj_df['gold_lon1'] = glon1
    conj_df['gold_lon2'] = glon2

    # sort by sat_t1
    conj_df = conj_df.sort_values(by='sat_t1')

    return conj_df


def reduce_conj(conj_df, lon_edge=1):
    """Reduce the conjunction DataFrame from goldcon to exclude repeats.

    Parameters
    ----------
    conj_df : pd.DataFrame
        conjunction starting time (sat_t1) and ending time (sat_t2), satellite
        starting (sat_lon1) and ending (sat_lon2) longitudes,
        gold time (gold_time), and gold starting (gold_lon1) and ending
        (gold_lon2) longitudes
    lon_edge : float
        how far the satellite lon 1 should be from GOLD lon edge
        default 1
        If GOLD extends from -60 to -20 and the satellite starts at -20.5,
        lon_edge means it needs to be between -59 and -21
        or else it would be excluded if other options are available
    Returns
    -------
    con_use : pd.DataFrame
        dataframe containing only the non-repeated and latest in time rows
        based on sat_t1 and sat_t2
    """
    # check starting time
    id = np.linspace(0, len(conj_df) - 1, len(conj_df))
    keep_ilocs = []

    # look for conjunctions with same starting time and pick longest
    # and latest
    for tu in conj_df['sat_t1'].unique():
        umask = (conj_df['sat_t1'] == tu)
        con_mini = conj_df[umask].copy()
        m_id = np.array(id[umask])
        lon_mask = ((con_mini['sat_lon1'] > (con_mini['gold_lon1'] + lon_edge))
                    & (con_mini['sat_lon1'] < (con_mini['gold_lon2']
                                               - lon_edge)))
        con_lon = con_mini[lon_mask]
        l_id = m_id[lon_mask]

        if len(con_lon) != 0:
            cm_tdif = con_lon['sat_t2'] - con_lon['sat_t1']
            ci = np.where(cm_tdif == np.max(cm_tdif))[0]
            keep_ilocs.append(int(l_id[ci[-1]]))
        else:
            cm_tdif = con_mini['sat_t2'] - con_mini['sat_t1']
            ci = np.where(cm_tdif == np.max(cm_tdif))[0]
            keep_ilocs.append(int(m_id[ci[-1]]))

    con2 = conj_df.iloc[keep_ilocs].copy()

    # check ending time
    id = np.linspace(0, len(con2) - 1, len(con2))
    keep_ilocs = []

    # look for conjunctions with same ending time and pick longest
    # and latest
    for tu in con2['sat_t2'].unique():
        umask = (con2['sat_t2'] == tu)
        con_mini = con2[umask].copy()
        m_id = np.array(id[umask])

        lon_mask = ((con_mini['sat_lon1'] > (con_mini['gold_lon1'] + lon_edge))
                    & (con_mini['sat_lon1'] < (con_mini['gold_lon2']
                                               - lon_edge)))
        con_lon = con_mini[lon_mask]
        l_id = m_id[lon_mask]

        if len(con_lon) != 0:
            cm_tdif = con_lon['sat_t2'] - con_lon['sat_t1']
            ci = np.where(cm_tdif == np.max(cm_tdif))[0]
            keep_ilocs.append(int(l_id[ci[-1]]))
        else:
            cm_tdif = con_mini['sat_t2'] - con_mini['sat_t1']
            ci = np.where(cm_tdif == np.max(cm_tdif))[0]
            keep_ilocs.append(int(m_id[ci[-1]]))

    # reduce dataframe
    con_use = con2.iloc[keep_ilocs].copy()

    # sort by sat_t1
    con_use = con_use.sort_values(by='sat_t1')
    return con_use


def gold_by_lon(gdc, sat_lat, sat_lon):
    """Get GOLD OI conjunctions for each latitude based on satellite longitude.

    Parameters
    ----------
    gdc : dictionary
        GOLD data including columns 'oi', 'lat', and 'lon'
    sat_lat : array-like
        geograhic latitudes of satellite
    sat_lon : array-like
        geogrpahic longitudes of satellite

    Returns
    -------
    oi_list : list-like
        GOLD OI at nearest lat and lon
    mlat_list : list-like
        GOLD magnetic latitude at nearest lat and lon
    """
    # make sure they are arrays and not dataframes
    sat_lat = np.asarray(sat_lat)
    sat_lon = np.asarray(sat_lon)

    # go through each geogrpahic latitude (200) and get oi at closest lon
    oi_list = []
    mlat_list = []
    for ul in np.unique(gdc['lat']):

        # get closest longitude to each lat
        slon = sat_lon[np.argmin(abs(sat_lat - ul))]
        loni = np.argmin(abs(gdc['lon'] - slon))
        lati = np.where(ul == gdc['lat'])[0]
        oi_list.append(np.unique(gdc['oi'][lati, loni])[0])
        mlat_list.append(np.unique(gdc['mlat'][lati, loni])[0])

    return oi_list, mlat_list
