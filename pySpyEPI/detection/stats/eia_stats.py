"""EIA detection info from pyValEIA.
updated 25 February 2026
alanahco

List Functions
--------------
eia_info
"""
import pandas as pd
from pyValEIA.eia.detection import eia_complete
import numpy as np


def eia_info(pass_id, time, lat, lon, mlat, mlon, lt, alt, density,
             den_type, mlat_val=30, filt='barrel_average', interpolate=1,
             barrel_envelope=True, envelope_lower=0.6,
             envelope_upper=0.2, barrel_radius=3, window_lat=2):
    """Detect and classify EIAs in plasma density data.

    Parameters
    ----------
    pass_id : str
        pass_id created in build.find_cand
    time, lat, lon, mlat, mlon, lt, alt : array-like
        arrays of time, geographic latitude and longitude, magnetic latitude
        and longitude, local time, and altitude
    density : array-like
        Plasma density data, e.g., TEC, electron density, or ion density
    den_type : str
        String specifying 'tec' if density is TEC or 'ne' for ion or electron
        density
    mlat_val : float
        how far to limit magnetic latitude
        default 30 degrees mlat
    filt : str
        Filter method(s) for density. An empty string means no filtering, and
        and underscore combines two methods in the order they are specified.
        Valid methods include 'barrel', 'median', 'mean', and 'average'
        default 'barrel_average'
    interpolate : int
        Interpolate data to a higher resolution; the integer determines the
        number of data points in the interpolated output (e.g.,
        len(`density`) * `interpolate), so a value of one or less means there
        will not be any interpolation
        default 1
    barrel_envelope : bool kwarg
        if True, barrel roll will include points inside an
        envelope, if False, no envelope will be used
        default True
    envelope_lower : double kwarg
        lower limit of envelope
        default 0.6 (6%) of min value from contact points
    envelope_upper : double kwarg
        upper limit of envelope
        default 0.2 (2%) of max value from contact points
    barrel_radius : double kwarg
        latitudinal radius of barrel
        default 3
    window_lat : double kwarg
        latitudinal width of moving window
        default 2

    Returns
    -------
    eia_df : DataFrame
        length of 1 containing columns:
            'pass_id', 'time_start', 'time_end', 'mlat_start', 'mlat_end',
               'glon_start', 'glon_end', 'glat_start', 'glat_end', 'lt',
               'eia_type', 'crest1_mlat', 'crest1_np', 'crest2_mlat',
               'crest2_np', 'crest3_mlat', 'crest3_np', 'alt'
    eia_state : string
        eia category/orientation string
    plats : array
        crest latitudes
    """
    # limit to +/- 30 degrees
    mlat_mask = (abs(mlat) < 30)
    mlat_lim = mlat[mlat_mask]
    mlon_lim = mlon[mlat_mask]
    time_lim = time[mlat_mask]
    lat_lim = lat[mlat_mask]
    lon_lim = lon[mlat_mask]
    lt_lim = lt[mlat_mask]
    alt_lim = alt[mlat_mask]
    density_lim = density[mlat_mask]

    lt_center = lt_lim[abs(mlat_lim).argmin()]

    # try to run eia_complete
    try:
        (sort_lat, den_filt2, eia_state,
         z_lat, plats, p3lats) = eia_complete(
             mlat_lim, density_lim, "Ne", filt=filt, interpolate=interpolate,
             barrel_envelope=barrel_envelope, envelope_lower=envelope_lower,
             envelope_upper=envelope_upper, barrel_radius=barrel_radius,
             window_lat=window_lat)
    except Exception:
        eia_state = 'unknown'
        plats = [None]

    if eia_state != 'unknown':
        # Find trough lat between peaks by using filtered density and sort lats
        if len(plats) == 2:
            tr_n = den_filt2[(sort_lat > min(plats)) & (sort_lat < max(plats))]
            tr_l = sort_lat[(sort_lat > min(plats)) & (sort_lat < max(plats))]

            tr_i = np.argmin(tr_n)
            tr_lat1 = [tr_l[tr_i]]
            tr_lat2 = np.nan
        elif len(plats) == 3:
            pl_min = min(plats)
            pl_max = max(plats)
            pl_med = plats[(plats != pl_min) & (plats != pl_max)]

            tr_n1 = den_filt2[(sort_lat > pl_min) & (sort_lat < pl_med)]
            tr_l1 = sort_lat[(sort_lat > pl_min) & (sort_lat < pl_med)]
            tr_i1 = np.argmin(tr_n1)
            tr_lat1 = [tr_l1[tr_i1]]

            tr_n2 = den_filt2[(sort_lat > pl_med) & (sort_lat < pl_max)]
            tr_l2 = sort_lat[(sort_lat > pl_med) & (sort_lat < pl_max)]
            tr_i2 = np.argmin(tr_n2)
            tr_lat2 = [tr_l2[tr_i2]]
        else:
            tr_lat1 = np.nan
            tr_lat2 = np.nan
    else:
        tr_lat1 = np.nan
        tr_lat2 = np.nan

    # EIA DataFrame
    eia_df = pd.DataFrame()
    eia_df['pass_id'] = [pass_id]
    eia_df['time_start'] = [pd.Timestamp(time_lim[0]).to_pydatetime()]
    eia_df['time_end'] = [pd.Timestamp(time_lim[-1]).to_pydatetime()]
    eia_df['mlat_start'] = [mlat_lim[0]]
    eia_df['mlat_end'] = [mlat_lim[-1]]
    eia_df['mlon_start'] = [mlon_lim[-1]]
    eia_df['mlon_end'] = [mlon_lim[-1]]
    eia_df['glon_start'] = [lon_lim[0]]
    eia_df['glon_end'] = [lon_lim[-1]]
    eia_df['glat_start'] = [lat_lim[0]]
    eia_df['glat_end'] = [lat_lim[-1]]

    eia_df['lt'] = [pd.Timestamp(lt_center).to_pydatetime()]
    eia_df['alt'] = [alt_lim[0]]
    eia_df['eia_type'] = [eia_state]

    # nan place holder
    for pi in range(3):
        ml_str = f'crest{pi + 1}_mlat'
        np_str = f'crest{pi + 1}_np'
        eia_df[ml_str] = [np.nan]
        eia_df[np_str] = [np.nan]

    # if plats is not None, replace nans by plats
    for pi, p in enumerate(plats):
        if p:
            mlat_loc = (abs(p - mlat_lim).argmin())
            ml_str = f'crest{pi + 1}_mlat'
            np_str = f'crest{pi + 1}_np'
            eia_df.at[0, ml_str] = [mlat_lim[mlat_loc]]
            eia_df.at[0, np_str] = [density_lim[mlat_loc]]

    eia_df['min_lat1'] = tr_lat1
    eia_df['min_lat2'] = tr_lat2

    return eia_df, eia_state, plats
