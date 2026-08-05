"""GOLD utilility functions.
Created 19 February 2026
alanahco

List Functions
--------------
time2lonRange
convert_to_mag
get_longitude_filter
calc_delta
get_conj_ind
time_update
myround_ud
"""
import numpy as np
import apexpy
import datetime as dt


def time2lonRange():
    """ return array of time to longitude range

    Parameters
    ----------
    NONE

    Returns
    -------
    time_to_lon_range : array-like
        array of times and their corresponding longitude ranges
        for GOLD

    Notes
    -----
    time_to_lon_range established by Mary Smirnova
    marysgk@umich.edu
    """
    time_to_lon_range = {
        '20:10': (-15, 9),
        '20:22': (-20, 8),
        '20:40': (-25, 8),
        '20:52': (-27, 6),
        '21:10': (-30, 5),
        '21:22': (-30, 5),
        '21:40': (-35, 3),
        '21:52': (-35, 0),
        '22:10': (-40, -5),
        '22:22': (-45, -10),
        '22:40': (-45, -12),
        '22:52': (-50, -15),
        '23:10': (-55, -20),
        '23:22': (-60, -25),
        '23:40': (-60, -25),
        '23:52': (-65, -30),
        '00:10': (-70, -35),
        '00:22': (-75, -40)
    }
    return time_to_lon_range


def convert_to_mag(lat, lon, time):
    """ Convert GOLD geographic to magnetic coordinates

    Parameters
    ----------
    lat : array-like
        2D array of geographic latitudes
    lon : array-like
        2D array of geographic longitudes
    time : datetime
        time for apex conversion

    Returns
    -------
    mlat : array-like
        2D array of magnetic latitudes
    mlon : array-like
        2D array of magnetic longitudes
    L_shells : array-like
        array of apex L_shells

    Notes
    -----
    function developed by Maria Smirnova

    Change made by alanahco@umich.edu on 09/29/2025
        instead of apex at 0km, using quasi-dipole
    """
    # Establish Apex
    apex = apexpy.Apex(date=time.year)

    # find valid indices (not NaN)
    valid_indices = ~np.isnan(lat) & ~np.isnan(lon)

    # Create nan array similar to mlat and mlon
    mlat = np.full_like(lat, np.nan)
    mlon = np.full_like(lon, np.nan)

    # iterate through valid indices and fill up mlat and mlon arrays
    if np.any(valid_indices):
        mlat[valid_indices], mlon[valid_indices] = apex.convert(
            lat[valid_indices], lon[valid_indices], 'geo', 'qd')

    # Get L shells
    aalt = apex.get_apex(mlat, 160)
    L_shells = 1.0 + aalt / apex.RE

    return mlat, mlon, L_shells


def get_longitude_filter(lon, file_time):
    """Determine longitude filter based on time difference within 5 minutes.

    Parameters
    ----------
    lon : array-like
        2D array of geographic longitudes
    filt_time : datetime
        time of data file

    Returns
    -------
    valid_lon : array-like
        mask of valid longitudes

    Notes
    -----
    If there are no valid longitudes
    an array of all True values is returned

    function developed by Maria Smirnova
    """

    # initialize time_to_lon_range
    time_to_lon_range = time2lonRange()

    # Convert dictionary keys to datetime.time objects for comparison
    time_filters = {dt.datetime.strptime(t, '%H:%M').time():
                    lon_range for t, lon_range in time_to_lon_range.items()}

    # Extract the time part from file_time (ignore the date)
    file_time_of_day = file_time.time()

    # Look for a time in the dictionary within 5 minutes of file's time of day
    for time_key, lon_range in time_filters.items():

        # Calculate the time difference
        # by converting both times to datetime objects
        file_time_dt = dt.datetime.combine(dt.datetime.min, file_time_of_day)
        time_key_dt = dt.datetime.combine(dt.datetime.min, time_key)

        # Difference in minutes
        time_diff = abs((file_time_dt - time_key_dt).total_seconds() / 60)

        # If within 5 minutes, apply the longitude filter
        if time_diff <= 10:
            valid_lon = (lon >= lon_range[0]) & (lon <= lon_range[1])
            return valid_lon

    # If no time matches, return a mask that keeps all data (no filtering)
    return np.ones_like(lon, dtype=bool)


def calc_delta(x0, y0, z0, x, y, z, phi, distance):
    """ calculate delta (angle between barrel and next point)

    Parameters
    ----------
    x0 : float
        starting x point
    y0 : float
        starting y point
    z0 : float
        starting z point
    x : float
        next x point
    y : float
        next y point
    z : float
        next z point
    phi : float
        bearing angle in radians
    distance : float
        distance between (x,y) and (x_center, y_center)
        where
        x_center = x0 + np.sin(phi) # radians
        y_center = y0 + np.cos(phi)
    Returns
    -------
    delta : float
        angle between barrel and next point in radians
    Notes
    -----
    If the data point is actually located outside the barrel, return pi

    function developed by Mary Smirnova
    marysgk@umich.edu

    Based on:
        https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023JA031963
        https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2015JA021723
    """
    # Establish dx, dy, and dz
    dx = x - x0
    dy = y - y0
    dz = z - z0

    # Calculate numerator
    num = dx**2 + dy**2 + dz**2

    # Calculate denominator
    denom = 2 * np.sqrt((dx * np.sin(phi) + dy * np.cos(phi))**2 + dz**2)

    # if num/denom is greater then 1, return pi (outside barrel)
    if abs(num / denom) > 1:
        return np.pi

    # calculate delta
    delta = np.arcsin(num / denom) - np.arctan(dz / (dx * np.sin(phi)
                                                     + dy * np.cos(phi)))
    return delta


def get_conj_ind(st, lon_range=None):
    """Get the index of the closest conjunction between swarm and GOLD

    Parameters
    ----------
    st : datetime
        time of conjunction
    lon_range: None type or list-like
        range of longitudes to search
        default None, -75 to 9 will be used (whole range)

    Returns
    -------
    min_ind : int
        index of conjunction closest in time and longitude
    list_inds : list
        list of indices that fit the longitude criteria without time criteria

    """
    # Get time to lon array:
    time_to_lon_range = time2lonRange()

    # use list from lon_range
    if lon_range is None:

        # set min_ind to none and list_inds to all
        min_ind = None
        list_inds = np.linspace(
            0, len(time_to_lon_range) - 1, len(time_to_lon_range))

    else:
        lon_min = min(lon_range)
        lon_max = max(lon_range)

        # Separate times, lons, and indices
        time_filters = [dt.datetime.strptime(t, '%H:%M').time()
                        for t, lon_range in time_to_lon_range.items()]

        lon_filters = [lon_range for t, lon_range in time_to_lon_range.items()]
        inds = np.linspace(0, 17, 18)

        if len(inds) != len(lon_filters):
            print('Warning: Wrong region may be chosen!!')

        swarm_time = st.time()

        # pick out the potential longitudes
        lons_new = []
        times_new = []
        list_inds = []
        for li in range(len(lon_filters)):
            l1 = lon_filters[li][0]
            l2 = lon_filters[li][1]

            # pick lon ranges that cover lon_min and lon_max
            if (lon_min > l1) & (lon_max < l2):
                lons_new.append(lon_filters[li])
                times_new.append(time_filters[li])
                list_inds.append(int(inds[li]))

        # From the new list of lons and times,
        # pick the closest time to Swarm time
        # initiate a starting difference
        t_dif = 100
        min_ind = 100
        for s, ts in enumerate(times_new):
            ts_hr = ts.hour + ts.minute / 60

            # update the hour so that 0-2 am are 24-26
            if ts_hr < 12:
                ts_hr = ts_hr + 24

            # get swarm hour
            sw_hr = swarm_time.hour + swarm_time.minute / 60

            # update the hour so that 0-2 am are 24-26
            if sw_hr < 12:
                sw_hr = sw_hr + 24

            # update t_dif and min_ind
            if abs(sw_hr - ts_hr) < t_dif:
                t_dif = abs(sw_hr - ts_hr)
                min_ind = int(list_inds[s])

    return min_ind, list_inds


def time_update(swarm_time, gold_time):
    """ Get file opening time from swarm time and gold time
    Parameters
    ----------
    swarm_time : datetime
        datetime of swarm pass
    gold_time : datetime HH:MM
        datetime of GOLD
    Returns
    -------
    st_open : datetime
        updated time for GOLD file opening purposes to be used in
        open_GOLD_barrel
    """
    # get swarm time and gold time fractional minutes

    st_hr = swarm_time.hour + swarm_time.minute / 60

    gt_hr = gold_time.hour + gold_time.minute / 60

    # if they are on different days and Swarm is > 0 and less than 6
    if (st_hr < 12) & (gt_hr > 12):
        # make it the day before
        st_update = swarm_time - dt.timedelta(days=1)

    # if they are different days and swarm is > 6
    elif (st_hr > 12) & (gt_hr < 12):
        # make it the next day
        st_update = swarm_time + dt.timedelta(days=1)

    # if they are on the same day (both > 6 or both < 6)
    else:
        # same day
        st_update = swarm_time

    st_open = dt.datetime(st_update.year, st_update.month, st_update.day,
                          gold_time.hour, gold_time.minute)

    return st_open


def myround_ud(x, base=5, updown='nearest'):
    """ Round array to the nearest base
    Parameters
    ----------
    x : array-like
        array of values to be rounded
    base : int kwarg
        base to be rounded to, 5 default
    updown : str
        string inicating whether rounding 'up', 'down', or 'nearest' (Default)
    Returns
    -------
    rounded array to nearest base
    """
    rounded_array = []

    for xx in x:
        if updown == 'nearest':
            rounded_array.append(int(base * np.round(float(xx) / base)))
        elif updown == 'up':
            rounded_array.append(int(base * np.ceil(float(xx) / base)))
        elif updown == 'down':
            rounded_array.append(int(base * np.floor(float(xx) / base)))

    return np.array(rounded_array)
