"""GOLD loading and file creation
Adapted from code created by Maria Smirnova
Created 19 February 2026
alanahco

List Functions
--------------
GOLD_barrel
GOLD_no_barrel
open_GOLD_nc
open_GOLD_barrel
open_GOLD_NObarrel
"""
import xarray as xr
import numpy as np
import datetime as dt
import pandas as pd
import scipy
import os
from pySpyEPI.utils.sats import gold_utils as go_ut


def GOLD_barrel(st, output_data_folder, lon_range=None, closest=False):
    """ create a barrel rolled file for desired swarm/ GOLD conjunction

    Parameters
    ----------
    st : datetime
        time of conjunction
    output_data_folder : str
        folder name for data output
    lon_range: None type or list-like
        range of longitudes to search
        default None, all longitudes will be used
    closest : bool kwarg False
        if True, only the closest in time and longitude file will be created
        otherwise, it will be all files that fit the longitude criteria
        if lon_range is None, this should be False

    Returns
    -------
    mat file of format name barrel_{save_time}.mat
    including interpolated oi, interpolated baseline, and depletion matrices

    Notes
    -----
    Based on code by Maria Smirnova
    """
    # Need to open GOLD file from day before if st.hour is between 0 and 3 (6
    # is used to be safe) The GOLD files range from 20 to 00:22 the next day
    # So if we have an early morning date, it needs the file before
    if st.hour < 12:
        # Get GOLD files for day before
        st_use = st - dt.timedelta(days=1)
        dcA, dcB, timeA, timeB = open_GOLD_nc(st_use)
    else:
        # Get GOLD files for same day
        dcA, dcB, timeA, timeB = open_GOLD_nc(st)

    # get index of conjunction
    min_ind, list_inds = go_ut.get_conj_ind(st, lon_range=lon_range)

    if closest:
        ind_use = [min_ind]
        if lon_range is None:
            print('Warning: need to supply lon range!')
    else:
        ind_use = list_inds

    # iterate through ind_use
    for ii in ind_use:

        # make sure it is an integer
        i = int(ii)

        # Check time to make sure that we are dealing with the same indices
        if len(dcA['time']) != 18:
            print('Check dcA time length')

        # get A and B times
        tA = dcA['time'][i]
        tB = dcB['time'][i]
        diff_sec = (tB - tA).total_seconds()

        # Save the final data, including the new interpolated oi,
        # lat, lon, baseline, and depletion
        save_time = tA.strftime('%Y%m%d_%H%M')
        mnth_fold = tA.strftime('%Y%m%d')

        # create extract folder if necessary
        extract_folder = os.path.join(
            output_data_folder, tA.strftime('%Y'), mnth_fold)
        os.makedirs(extract_folder, exist_ok=True)

        # Out folder name complete
        output_data_filename = os.path.join(
            extract_folder, f'barrel_{save_time}.mat')

        # check if file already exists
        if os.path.exists(output_data_filename):
            continue

        # if file does not exist, continue on
        # Check that they are close in time
        if diff_sec > 240:
            print('diff_sec > 240')

        # double check that they are in different hemispheres
        if dcA['hemisphere'][i] == dcB['hemisphere'][i]:
            print('same hemispheres')

        # combine A and B together after those checks
        # Get necessary Params A and B
        latA = dcA["lat"][i, :, :]
        lonA = dcA["lon"][i, :, :]
        oiA = dcA["oi"][i, :, :]
        szaA = dcA["sza"][i, :, :]
        mlatA = dcA["mlat"][i, :, :]
        mlonA = dcA["mlon"][i, :, :]

        # Filter data
        validA = ((mlatA > -40) & (mlatA < 40)
                  & go_ut.get_longitude_filter(lonA, tA))
        latA, lonA, oiA, szaA, mlatA, mlonA = (latA[validA], lonA[validA],
                                               oiA[validA], szaA[validA],
                                               mlatA[validA], mlonA[validA])

        # remove data from A below the equator
        lat_mask = latA >= 0
        latA = latA[lat_mask]
        lonA = lonA[lat_mask]
        oiA = oiA[lat_mask]
        szaA = szaA[lat_mask]
        mlonA = mlonA[lat_mask]
        mlatA = mlatA[lat_mask]

        latB = dcB["lat"][i, :, :]
        lonB = dcB["lon"][i, :, :]
        oiB = dcB["oi"][i, :, :]
        szaB = dcA["sza"][i, :, :]
        mlatB = dcB["mlat"][i, :, :]
        mlonB = dcB["mlon"][i, :, :]

        # Limit magnetic latitudinally
        validB = ((mlatB > -40) & (mlatB < 40)
                  & go_ut.get_longitude_filter(lonB, tB))
        latB, lonB, oiB, szaB, mlatB, mlonB = (latB[validB], lonB[validB],
                                               oiB[validB], szaB[validB],
                                               mlatB[validB], mlonB[validB])

        lat = np.concatenate((latA, latB), axis=0)
        lon = np.concatenate((lonA, lonB), axis=0)
        oi = np.concatenate((oiA, oiB), axis=0)

        # Remove and nan values form oiA and other params
        valid_idx = ~np.isnan(oi)
        lat_valid = lat[valid_idx]
        lon_valid = lon[valid_idx]
        oi_valid = oi[valid_idx]

        # Define regular grid for lat and lon
        lat_reg = np.linspace(min(lat_valid), max(lat_valid), 200)
        lon_reg = np.linspace(min(lon_valid), max(lon_valid), 100)
        lon_grid, lat_grid = np.meshgrid(lon_reg, lat_reg)

        # Interpolate 'oi' values on the regular grid
        oi_grid = scipy.interpolate.griddata((lon_valid, lat_valid), oi_valid,
                                             (lon_grid, lat_grid),
                                             method='linear')

        # Adjusting scales
        x = lon_reg / 12
        y = lat_reg / 5

        # Finite values for scaling
        g = 24 + np.min(oi_grid[np.isfinite(oi_grid)])
        G = 0.012
        z = np.log10((oi_grid + g) / G)

        # Start at the highest point on the terrain
        max_index = np.unravel_index(np.nanargmax(z), z.shape)
        x0, y0, z0 = x[max_index[1]], y[max_index[0]], z[max_index]

        # Initialize the bearing angle in radians
        phi = np.radians(np.random.uniform(0, 360))
        n = len(z.flatten())  # Number of points
        visited_points_matrix = np.full(z.shape, np.nan)
        visited_points_matrix[max_index] = z0

        # Initialize delta phi
        delta_phi = 20  # degrees

        for num_rolls in range(n):
            # Calculate the number of rolls as a % of the total grid points
            roll_percentage = (num_rolls / n) * 100

            # Update Δφ based on the percentage of executed rolls
            if roll_percentage >= 80:
                delta_phi = 60
            elif roll_percentage >= 60:
                delta_phi = 50
            elif roll_percentage >= 40:
                delta_phi = 40
            elif roll_percentage >= 20:
                delta_phi = 30

            # Initialize variables to find the minimum delta
            min_delta = float('inf')
            min_index = None

            # Iterate over all points in the z matrix to find the candidate
            for i in range(z.shape[0]):  # Iterate over the rows (y)
                for j in range(z.shape[1]):  # Iterate over the columns (x)

                    # Skip if the current point is the starting point
                    if (x0 == x[j] and y0 == y[i] and z0 == z[i, j]
                            or z[i, j] == np.nan):
                        continue

                    # Is point (i, j) inside hit zone?
                    x_center = x0 + np.sin(phi)  # radians
                    y_center = y0 + np.cos(phi)
                    distance = np.sqrt((x[j] - x_center) ** 2
                                       + (y[i] - y_center) ** 2)

                    if distance > 1:
                        continue

                    # Calculate delta for this candidate point
                    delta = go_ut.calc_delta(x0, y0, z0, x[j], y[i], z[i, j],
                                             phi, distance)

                    # Update the minimum delta and corresponding candidate
                    if np.abs(delta) < np.abs(min_delta):
                        min_delta = delta
                        min_index = (i, j)

            # If no valid candidate is found, break the loop
            if min_index is None:
                phi += np.radians(np.random.uniform(-delta_phi, delta_phi))
                continue

        # Update the current point to the one with the minimum delta
            i, j = min_index
            x0, y0, z0 = x[j], y[i], z[i, j]
            # Store the visited point's z value in the matrix
            visited_points_matrix[i, j] = z0

            # Vary the bearing angle φ by a random variable within ±Δφ
            phi += np.radians(np.random.uniform(-delta_phi, delta_phi))

            # Update max_index for the new current position
            max_index = min_index

        # Rolling Ball 2D Rastering --------------------
        # Convert visited points' z values back to oi
        visited_oi_matrix = G * (10 ** visited_points_matrix) - g

        # Step 2: Interpolate NaN values in the oi matrix
        # Get the indices of the non-NaN values
        non_nan_indices = np.argwhere(~np.isnan(visited_oi_matrix))
        non_nan_x = non_nan_indices[:, 1]  # X coordinates
        non_nan_y = non_nan_indices[:, 0]  # Y coordinates

        # OI values at these points
        non_nan_values = visited_oi_matrix[non_nan_y, non_nan_x]

        # Create a grid for interpolation
        grid_x, grid_y = np.meshgrid(np.arange(visited_oi_matrix.shape[1]),
                                     np.arange(visited_oi_matrix.shape[0]))

        # Interpolate NaN values using griddata
        interpolated_oi_matrix = scipy.interpolate.griddata((non_nan_x,
                                                             non_nan_y),
                                                            non_nan_values,
                                                            (grid_x, grid_y),
                                                            method='linear')

        # Net depletion
        baseline_subtracted_matrix = oi_grid - interpolated_oi_matrix

        # Save time list as a string
        time_matA = np.array(tA.strftime('%Y-%m-%d/%H:%M:%S')).astype(object)
        time_matB = np.array(tB.strftime('%Y-%m-%d/%H:%M:%S')).astype(object)

        # save as .mat file
        scipy.io.savemat(output_data_filename, {
            'data': {'timeA': time_matA, 'timeB': time_matB, 'lat': lat_grid,
                     'lon': lon_grid, 'oi': oi_grid,
                     'baseline': interpolated_oi_matrix,
                     'depletion': baseline_subtracted_matrix}})


def GOLD_no_barrel(st, output_data_folder, lon_range=None, closest=False):
    """create a non-barrel file for desired swarm/ GOLD conjunction.

    Parameters
    ----------
    st : datetime
        time of conjunction
    output_data_folder : str
        folder name for data output
    lon_range: None type or list-like
        range of longitudes to search
        default None, all longitudes will be used
    closest : bool kwarg False
        if True, only the closest in time and longitude file will be created
        otherwise, it will be all files that fit the longitude criteria
        if lon_range is None, this should be False

    Returns
    -------
    mat file of format name barrel_{save_time}.mat
    including interpolated oi, interpolated baseline, and depletion matrices

    Notes
    -----
    Based on code by Maria Smirnova
    """
    # Need to open GOLD file from day before if st.hour is between 0 and 3 (6
    # is used to be safe) The GOLD files range from 20 to 00:22 the next day
    # So if we have an early morning date, it needs the file before
    if st.hour < 12:
        # Get GOLD files for day before
        st_use = st - dt.timedelta(days=1)
        dcA, dcB, timeA, timeB = open_GOLD_nc(st_use)
    else:
        # Get GOLD files for same day
        dcA, dcB, timeA, timeB = open_GOLD_nc(st)

    # get index of conjunction
    min_ind, list_inds = go_ut.get_conj_ind(st, lon_range=lon_range)

    if closest:
        ind_use = [min_ind]
        if lon_range is None:
            print('Warning: need to supply lon range!')
    else:
        ind_use = list_inds

    # iterate through ind_use
    for ii in ind_use:

        # make sure it is an integer
        i = int(ii)

        # Check time to make sure that we are dealing with the same indices
        if len(dcA['time']) != 18:
            print('Check dcA time length')

        # get A and B times
        tA = dcA['time'][i]
        tB = dcB['time'][i]
        diff_sec = (tB - tA).total_seconds()

        # Save the final data, including the new interpolated oi,
        # lat, lon, baseline, and depletion
        save_time = tA.strftime('%Y%m%d_%H%M')
        mnth_fold = tA.strftime('%Y%m%d')

        # create extract folder if necessary
        extract_folder = os.path.join(
            output_data_folder, tA.strftime('%Y'), mnth_fold)
        os.makedirs(extract_folder, exist_ok=True)

        # Out folder name complete
        output_data_filename = os.path.join(
            extract_folder, f'NObarrel_{save_time}.mat')

        # check if file already exists
        if os.path.exists(output_data_filename):
            continue

        # if file does not exist, continue on
        # Check that they are close in time
        if diff_sec > 240:
            print('diff_sec > 240')

        # double check that they are in different hemispheres
        if dcA['hemisphere'][i] == dcB['hemisphere'][i]:
            print('same hemispheres')

        # combine A and B together after those checks
        # Get necessary Params A and B
        latA = dcA["lat"][i, :, :]
        lonA = dcA["lon"][i, :, :]
        oiA = dcA["oi"][i, :, :]
        szaA = dcA["sza"][i, :, :]
        mlatA = dcA["mlat"][i, :, :]
        mlonA = dcA["mlon"][i, :, :]

        # Filter data
        validA = ((mlatA > -40) & (mlatA < 40)
                  & go_ut.get_longitude_filter(lonA, tA))
        latA, lonA, oiA, szaA, mlatA, mlonA = (latA[validA], lonA[validA],
                                               oiA[validA], szaA[validA],
                                               mlatA[validA], mlonA[validA])

        # remove data from A below the equator
        lat_mask = latA >= 0
        latA = latA[lat_mask]
        lonA = lonA[lat_mask]
        oiA = oiA[lat_mask]
        szaA = szaA[lat_mask]
        mlonA = mlonA[lat_mask]
        mlatA = mlatA[lat_mask]

        latB = dcB["lat"][i, :, :]
        lonB = dcB["lon"][i, :, :]
        oiB = dcB["oi"][i, :, :]
        szaB = dcA["sza"][i, :, :]
        mlatB = dcB["mlat"][i, :, :]
        mlonB = dcB["mlon"][i, :, :]

        # Limit magnetic latitudinally
        validB = ((mlatB > -40) & (mlatB < 40)
                  & go_ut.get_longitude_filter(lonB, tB))
        latB, lonB, oiB, szaB, mlatB, mlonB = (latB[validB], lonB[validB],
                                               oiB[validB], szaB[validB],
                                               mlatB[validB], mlonB[validB])

        lat = np.concatenate((latA, latB), axis=0)
        lon = np.concatenate((lonA, lonB), axis=0)
        oi = np.concatenate((oiA, oiB), axis=0)
        mlat = np.concatenate((mlatA, mlatB), axis=0)
        mlon = np.concatenate((mlonA, mlonB), axis=0)

        # Remove and nan values form oiA and other params
        valid_idx = ~np.isnan(oi)
        lat_valid = lat[valid_idx]
        lon_valid = lon[valid_idx]
        mlat_valid = mlat[valid_idx]
        mlon_valid = mlon[valid_idx]
        oi_valid = oi[valid_idx]

        # Define regular grid for lat and lon
        lat_reg = np.linspace(min(lat_valid), max(lat_valid), 200)
        lon_reg = np.linspace(min(lon_valid), max(lon_valid), 100)
        lon_grid, lat_grid = np.meshgrid(lon_reg, lat_reg)

        mlat_reg = np.linspace(min(mlat_valid), max(mlat_valid), 200)
        mlon_reg = np.linspace(min(mlon_valid), max(mlon_valid), 100)
        mlon_grid, mlat_grid = np.meshgrid(mlon_reg, mlat_reg)

        # Interpolate 'oi' values on the regular grid
        oi_grid = scipy.interpolate.griddata((lon_valid, lat_valid), oi_valid,
                                             (lon_grid, lat_grid),
                                             method='linear')

        # Save time list as a string
        time_matA = np.array(tA.strftime('%Y-%m-%d/%H:%M:%S')).astype(object)
        time_matB = np.array(tB.strftime('%Y-%m-%d/%H:%M:%S')).astype(object)

        # save as .mat file
        scipy.io.savemat(output_data_filename, {
            'data': {'timeA': time_matA, 'timeB': time_matB, 'lat': lat_grid,
                     'lon': lon_grid, 'oi': oi_grid, 'mlat': mlat_grid,
                     'mlon': mlon_grid}})


def open_GOLD_nc(st, fdir):
    """ open and process GOLD nc files

    Parameters
    ----------
    st : datetime
        day to open file
    fdir : string
        location of file before /Y

    Returns
    -------
    dcA : dictionary
        A channel data dictionary
    dcB : dictionary
        B channel data dictionary
    time : array-like
        time array

    Notes
    -----
    Based on code by Maria Smirnova
    """

    # Establish filename
    fil = f"GOLD_L2_NMAX_{st.strftime('%Y_%j')}_v05_r01_c01.nc"

    fname = os.path.join(fdir, st.strftime('%Y'), fil)
    ds = xr.open_dataset(fname)

    # turn into dictionary
    data = {var: ds[var].values for var in ds.data_vars}

    # Extract Parameters
    # get channel, time, and hemsiphere [Shape: (nscans)]
    channel = data['channel'].astype(str)
    time = pd.to_datetime(data['scan_start_time'].astype(str), utc=True)
    hemisphere = data['hemisphere'].astype(str)

    # Get lat, lon, oi, and sza [Shape: (nscans, nlats, nlons)]
    lat = data['latitude']
    lon = data['longitude']
    oi = data['radiance_oi_1356']
    sza = data['solar_zenith_angle']

    # Calculate geomagnetic parameters
    mlat, mlon, lshells = go_ut.convert_to_mag(lat, lon, time)

    # Separate Measurements into Channels A and B (North and South)
    ch_str = 'CHA'
    channelA = channel[channel == ch_str]
    timeA = time[channel == ch_str]
    hemisphereA = hemisphere[channel == ch_str]
    latA = lat[channel == ch_str, :, :]
    lonA = lon[channel == ch_str, :, :]
    oiA = oi[channel == ch_str, :, :]
    szaA = sza[channel == ch_str, :, :]
    mlatA = mlat[channel == ch_str, :, :]
    mlonA = mlon[channel == ch_str, :, :]

    # create a dictionary of A data
    dcA = {
        "time": timeA, "oi": oiA, "lat": latA, "lon": lonA, "mlat": mlatA,
        "mlon": mlonA, "sza": szaA, "hemisphere": hemisphereA,
        "channel": channelA}

    # B southern hemisphere
    ch_str = 'CHB'
    channelB = channel[channel == ch_str]
    timeB = time[channel == ch_str]
    hemisphereB = hemisphere[channel == ch_str]
    latB = lat[channel == ch_str, :, :]
    lonB = lon[channel == ch_str, :, :]
    oiB = oi[channel == ch_str, :, :]
    szaB = sza[channel == ch_str, :, :]
    mlatB = mlat[channel == ch_str, :, :]
    mlonB = mlon[channel == ch_str, :, :]

    # create a dictionary of B data
    dcB = {
        "time": timeB, "oi": oiB, "lat": latB, "lon": lonB, "mlat": mlatB,
        "mlon": mlonB, "sza": szaB, "hemisphere": hemisphereB,
        "channel": channelB}

    return dcA, dcB, timeA, timeB


def open_GOLD_barrel(st, fdir):

    """ Open GOLD files created by GOLD.GOLD_fxns.GOLD_barrel
    Parameters
    ----------
    st : datetime
        datetime of mat file to open
        will open closest file to time provided
    fdir : string
        file directory
    Returns
    -------
    Gdc : dictionary
        dictionary of GOLD data containing timeA, timeB, oi, lat, lon, base_oi,
        deplete_oi, mlat, mlon, lshell
    """
    # Separate time range array
    time_to_lon_range = go_ut.time2lonRange()

    # looking for file of closest time based on time_to_lon_range
    time_filters = [dt.datetime.strptime(t, '%H:%M').time()
                    for t, lon_range in time_to_lon_range.items()]

    # From the list of times, pick the closest time to st
    # initiate a starting difference
    min_ind = 100
    t_dif = 100
    t_use = 100

    # convert st to fractional hour
    st_hr = st.hour + st.minute / 60
    # update the hour so that 0-2 am are 24-26
    if st_hr < 12:
        st_hr = st_hr + 24

    # set a second fractional hour to check against t_use
    st_use = st.hour + st.minute / 60

    # find the closest index
    for s, ts in enumerate(time_filters):
        ts_hr = ts.hour + ts.minute / 60
        t_ch = ts.hour + ts.minute / 60
        # update the hour so that 0-2 am are 24-26
        if ts_hr < 12:
            ts_hr = ts_hr + 24

        # update t_dif and min_ind
        if abs(st_hr - ts_hr) < t_dif:
            min_ind = int(s)
            t_use = t_ch
            t_dif = abs(st_hr - ts_hr)

    # use min_ind to open the closest file
    # if they are on different days
    if (st_use < 12) & (t_use > 12):

        # make it the day before
        st_update = st - dt.timedelta(days=1)

        st_open = dt.datetime(st_update.year, st_update.month, st_update.day,
                              time_filters[min_ind].hour,
                              time_filters[min_ind].minute)
    elif (st_use > 12) & (t_use < 12):

        # make it the next day
        st_update = st + dt.timedelta(days=1)

        st_open = dt.datetime(st_update.year, st_update.month, st_update.day,
                              time_filters[min_ind].hour,
                              time_filters[min_ind].minute)
    else:
        # same day
        st_open = dt.datetime(st.year, st.month, st.day,
                              time_filters[min_ind].hour,
                              time_filters[min_ind].minute)

    fname = os.path.join(fdir, st_open.strftime('%Y'),
                         st_open.strftime('%Y%m%d'),
                         f"barrel_{st_open.strftime('%Y%m%d_%H%M')}.mat")

    Gdc = {}

    # Try original time, then +8 min, then -8 min
    offsets = [0]
    for i in range(1, 9):
        offsets.extend([-i, i])

    for offset in offsets:
        if offset == 0:
            st_new = st_open
        else:
            st_new = st_open + dt.timedelta(minutes=offset)
            # print(f"Trying {st_new.strftime('%Y%m%d_%H%M')}")

        fname = os.path.join(
            fdir,
            st_new.strftime('%Y'),
            st_new.strftime('%Y%m%d'),
            f"barrel_{st_new.strftime('%Y%m%d_%H%M')}.mat"
        )

        if not os.path.exists(fname):
            continue

        # ---------------- Open GOLD file ----------------
        GOLD = scipy.io.loadmat(fname)
        GOLD_data = GOLD['data']

        oi_grid = GOLD_data['oi'][0][0]
        lat_grid = GOLD_data['lat'][0][0]
        lon_grid = GOLD_data['lon'][0][0]
        base_oi = GOLD_data['baseline'][0][0]
        deplete_oi = GOLD_data['depletion'][0][0]
        timeA_GOLD = GOLD_data['timeA'][0][0]
        timeB_GOLD = GOLD_data['timeB'][0][0]

        format_string = "%Y-%m-%d/%H:%M:%S"
        tA = dt.datetime.strptime(timeA_GOLD[0][0][0], format_string)
        tB = dt.datetime.strptime(timeB_GOLD[0][0][0], format_string)

        mlat_grid, mlon_grid, gold_ls = go_ut.convert_to_mag(
            lat_grid, lon_grid, tA
        )

        Gdc = {
            "timeA": tA,
            "timeB": tB,
            "oi": oi_grid,
            "lat": lat_grid,
            "lon": lon_grid,
            "base_oi": base_oi,
            "deplete_oi": deplete_oi,
            "mlat": mlat_grid,
            "mlon": mlon_grid,
            "lshell": gold_ls,
        }

        break

    if not Gdc:
        print("Attempt unsuccessful")

    return Gdc


def open_GOLD_NObarrel(st, fdir):

    """ Open GOLD files created by GOLD.GOLD_fxns.GOLD_barrel
    Parameters
    ----------
    st : datetime
        datetime of mat file to open
        will open closest file to time provided
    fdir : string
        file directory
    Returns
    -------
    Gdc : dictionary
        dictionary of GOLD data containing timeA, timeB, oi, lat, lon, base_oi,
        deplete_oi, mlat, mlon, lshell

    """
    # Separate time range array
    time_to_lon_range = go_ut.time2lonRange()

    # looking for file of closest time based on time_to_lon_range
    time_filters = [dt.datetime.strptime(t, '%H:%M').time()
                    for t, lon_range in time_to_lon_range.items()]

    # From the list of times, pick the closest time to st
    # initiate a starting difference
    min_ind = 100
    t_dif = 100
    t_use = 100

    # convert st to fractional hour
    st_hr = st.hour + st.minute / 60

    # update the hour so that 0-2 am are 24-26
    if st_hr < 12:
        st_hr = st_hr + 24

    # set a second fractional hour to check against t_use
    st_use = st.hour + st.minute / 60

    # find the closest index
    for s, ts in enumerate(time_filters):
        ts_hr = ts.hour + ts.minute / 60
        t_ch = ts.hour + ts.minute / 60
        # update the hour so that 0-2 am are 24-26
        if ts_hr < 12:
            ts_hr = ts_hr + 24

        # update t_dif and min_ind
        if abs(st_hr - ts_hr) < t_dif:
            min_ind = int(s)
            t_use = t_ch
            t_dif = abs(st_hr - ts_hr)

    # use min_ind to open the closest file
    # if they are on different days
    if (st_use < 12) & (t_use > 12):

        # make it the day before
        st_update = st - dt.timedelta(days=1)

        st_open = dt.datetime(st_update.year, st_update.month, st_update.day,
                              time_filters[min_ind].hour,
                              time_filters[min_ind].minute)
    elif (st_use > 12) & (t_use < 12):

        # make it the next day
        st_update = st + dt.timedelta(days=1)

        st_open = dt.datetime(st_update.year, st_update.month, st_update.day,
                              time_filters[min_ind].hour,
                              time_filters[min_ind].minute)
    else:
        # same day
        st_open = dt.datetime(st.year, st.month, st.day,
                              time_filters[min_ind].hour,
                              time_filters[min_ind].minute)

    fname = os.path.join(fdir, st_open.strftime('%Y'),
                         st_open.strftime('%Y%m%d'),
                         f"NObarrel_{st_open.strftime('%Y%m%d_%H%M')}.mat")

    Gdc = {}

    # Try original time, then +8 min, then -8 min
    offsets = [0]
    for i in range(1, 9):
        offsets.extend([-i, i])

    for offset in offsets:
        if offset == 0:
            st_new = st_open
        else:
            st_new = st_open + dt.timedelta(minutes=offset)
            # print(f"Trying {st_new.strftime('%Y%m%d_%H%M')}")

        fname = os.path.join(
            fdir,
            st_new.strftime('%Y'),
            st_new.strftime('%Y%m%d'),
            f"NObarrel_{st_new.strftime('%Y%m%d_%H%M')}.mat"
        )

        if not os.path.exists(fname):
            continue

        # ---------------- Open GOLD file ----------------
        GOLD = scipy.io.loadmat(fname)
        GOLD_data = GOLD['data']
        oi_grid = GOLD_data['oi'][0][0]
        lat_grid = GOLD_data['lat'][0][0]
        lon_grid = GOLD_data['lon'][0][0]
        timeA_GOLD = GOLD_data['timeA'][0][0]
        timeB_GOLD = GOLD_data['timeB'][0][0]

        format_string = "%Y-%m-%d/%H:%M:%S"
        tA = dt.datetime.strptime(timeA_GOLD[0][0][0], format_string)
        tB = dt.datetime.strptime(timeB_GOLD[0][0][0], format_string)

        mlat_grid, mlon_grid, gold_ls = go_ut.convert_to_mag(
            lat_grid, lon_grid, tA
        )

        Gdc = {
            "timeA": tA,
            "timeB": tB,
            "oi": oi_grid,
            "lat": lat_grid,
            "lon": lon_grid,
            "mlat": mlat_grid,
            "mlon": mlon_grid,
            "lshell": gold_ls,
        }

        break

    if not Gdc:
        print("Attempt unsuccessful")

    return Gdc
