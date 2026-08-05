"""Functions to aid swarm  loading.
Created 19 February 2026
alanahco

List Functions
--------------
"""
import numpy as np
import cdflib
from pySpyEPI.utils import calc
from pySpyEPI.barrel import barrel_utils as bar_uts


def load_cdf_data(file_path, variable_names):
    """Load CDF File data from Swarm.

    Parameters
    ----------
    file_path : string
        path of file including filename
    varible_names: list
        variables in cdf file

    Returns
    -------
    dictionary containing variable_names from cdf file
    """
    cdf = cdflib.CDF(file_path)
    return {var: cdf.varget(var) for var in variable_names}, cdf


# CHange below function to have a radius input
def extract_common_coords(cdf_file, satellite):
    """ Extract data from cdf_file
    Parameters
    ----------
    cdf_file : dictionary
        cdf file dictionary returned from load_cdf_data
    satellite : str
        'A', 'B', or 'C' for Swarm
    rad : radius
    Returns
    -------
        lat : array-like
            geo latitudes
        lon : array-like
            geo longitudes
        time : array-like
            time
        rad : array-like
            satellite altitude in km
    """
    lat = cdf_file.varget('Latitude')
    lon = cdf_file.varget('Longitude')
    epoch = cdf_file.varget('Timestamp')
    time = cdflib.cdfepoch.to_datetime(epoch)

    return lat, lon, time


def create_ne_flag(data, Np='ne'):
    """ Create a new ne flag
    Parameters
    ----------
    data : dataframe
        IRR dataset
    Np : str
        default 'Ne'
        plasma density column name in DataFrame
    Returns
    -------
    flg : array of integers
        0 : okay to use
        5 : don't use, Ne Values that are too low
            Ne < 10^4
            OR
            bNe < 10^4
    Notes
    -----
    based on

    Jin Y, Xiong C, Clausen L, Spicher A, Kotova D, et al. 2020.
    Ionospheric plasma irregularities based on in situ measurements
    from the swarm satellites. J Geophys Res 125(7): e2020JA028103.
    https://doi.org/10.1029/2020JA028103

    Kotova D, Jin Y & Miloch W 2022. Interhemispheric variability
    of the electron density and derived parameters by the Swarm
    satellites during different solar activity. J. Space Weather
    Space Clim. 12, 12. https://doi.org/10.1051/swsc/2022007.
    """

    # Initiate an array of Zeros
    flg = np.zeros(len(data[Np]))

    Ne_mask = (data[Np] < 10**4)

    flg[Ne_mask] = 5

    return flg


def Te_calibrate_2021(EFI_data, Satellite):
    """Calibrate the electron temperature based on calibration
        from Lomizde et al 2021.

    Parameters
    ----------
    EFI_data : pandas DataFrame
        DataFrame of Swarm EFI data created by Load_Swarm.py
    Satellite : str
        'A', 'B', or 'C'

    Returns
    -------
    Te_adj : array-like
        calibrated electron temperature
    dA : array-like
        error propagated from te_error (random error) and calibration error

    Notes
    -----
    error based on df = |df/dx|*dx + |df/dy|* dy
    """
    Te_og = EFI_data["te"]
    Te_adj = Te_og.copy()
    LP_flag = EFI_data['lp_flag']
    dT = abs(EFI_data['te_error'])
    dA = dT.copy()

    # Adjust by Satellite
    if Satellite == "A":

        # high gain
        C1_high = 1.2844
        dC1_high = 0.0464
        C2_high = 1083
        dC2_high = 127

        # low gain
        C_low = 723
        dC_low = 35

    elif Satellite == "B":
        # high gain
        C1_high = 1.1626
        dC1_high = 0.0622
        C2_high = 827
        dC2_high = 190

        # low gain
        C_low = 698
        dC_low = 41

    elif Satellite == "C":

        # high gain
        C1_high = 1.2153
        dC1_high = 0.0410
        C2_high = 916
        dC2_high = 111

        # low gain
        C_low = 682
        dC_low = 32

    # high gain
    Te_adj[LP_flag != 3] = C1_high * Te_og[LP_flag != 3] - C2_high

    # error
    dA[LP_flag != 3] = (Te_og[LP_flag != 3] * dC1_high
                        + C1_high * dT[LP_flag != 3] + dC2_high)

    # if all dTe is Nan, use error from formula
    if np.all(np.isnan(dT[LP_flag != 3])):
        dA[LP_flag != 3] = (Te_og[LP_flag != 3] * dC1_high + dC2_high)

    # low gain
    Te_adj[LP_flag == 3] = Te_og[LP_flag == 3] - C_low

    # error
    dA[LP_flag == 3] = (dT[LP_flag == 3] + dC_low)

    # if all dTe is Nan, use error from formula
    if np.all(np.isnan(dT[LP_flag == 3])):
        dA[LP_flag == 3] = dC_low

    return Te_adj, dA


def process_Te(EFI_data, satellite, lp_s=30, freq=2, mlat_cut=55, jump=500):
    """ process EFI electron temperature by calibrating the Te, removing the
        bad flags, removing the outliers, and then median filtering
    Parameters
    ----------
    EFI_data : pandas DataFrame
        DataFrame of Swarm EFI data created by Load_Swarm.py
    satellite : str
        Swarm Satellite for calibration purposes
    lp_s : int kwarg
        number of seconds after lp sweep to remove
        default 30s
        30 seconds decided observationally
        and confirmed by Shasha on 10/17/2025
    mlat_cut : float
        need to detect outliers for limited lattiudinal ranges, mlat_cut is the
        magnetic latitude cutoff. beyond that range, te_meds_all will be np.nan
        default is 55 degrees mlat (meaning within +/- 55)
    jump : float
        maximum jump allowed between one point and another in temperature data
        defualt is 500 K
    Returns
    -------
    te_meds_all : array-like
        median filtered Te that is calibrated, flag checked,
        and outlier removed that is same length as EFI_data['Te']
        outside of mlat_cut is np.nan
    dte_adj_all : array-like
        error propagated from te_error (random error)
        and calibration error from Te_calibrate_2021 (not te_meds_all)
    te_sweep_all : array-like
        te values at lp sweep for plotting purposes
    jump_all : array-like
        0 if there was no jump of > jump (default 500K)
        1 if there was a jump from 1 point to the next (not gradual)
    te_jump_b4 : array-like
        te difference between each point for jump from cleaned te
    te_jump_af : array-like
        te difference between each point from median filtered te
    NOTES
    -----
    Consider propagating the error through the rolling median
    Need to use mlat_cut because outliers need to be handled with limited
        segments
    Updated so that flag 21 is also included
        Note: no flag 21 for version 0602, but nominal flag no error computed
        for version 0701
    """

    # Limit by mlat_cut
    efi_lim = EFI_data[abs(EFI_data['mlat']) < mlat_cut]

    # use find gaps to get mlat segments
    sw_idx = efi_lim.index.values
    lat_gaps = bar_uts.find_all_gaps(sw_idx)

    # Make nan lists of same length as EFI_data
    te_meds_all = np.array([np.nan] * len(EFI_data))
    dte_adj_all = np.array([np.nan] * len(EFI_data))
    te_sweep_all = np.array([np.nan] * len(EFI_data))
    jump_all = np.array([np.nan] * len(EFI_data))
    te_jump_b4 = np.array([np.nan] * len(EFI_data))
    te_jump_af = np.array([np.nan] * len(EFI_data))

    # iterate through the gaps
    for mgap in range(len(lat_gaps) - 1):

        # get each segment
        efi_lat = efi_lim.iloc[lat_gaps[mgap]:lat_gaps[mgap + 1]]
        idx_all = sw_idx[lat_gaps[mgap]:lat_gaps[mgap + 1]]

        # calibrate -----------------------------------------------------------
        te_cal = efi_lat["te"].values.copy()
        te_cal, dte_adj = Te_calibrate_2021(efi_lat, satellite)

        # update full list
        dte_adj_all[idx_all] = dte_adj

        # sweep ---------------------------------------------------------------
        not_sweep = (efi_lat["lp_flag"] != 9)
        te_sweep = te_cal.copy()
        te_sweep[not_sweep] = np.nan

        # update full list
        te_sweep_all[idx_all] = te_sweep

        # remove bad flags ----------------------------------------------------
        te_flagged = te_cal.copy()
        te_rem = ((efi_lat["te_flag"] != 10) & (efi_lat["te_flag"] != 19)
                  & (efi_lat["te_flag"] != 20) & (efi_lat["te_flag"] != 21))
        te_flagged[te_rem] = np.nan

        # remove outliers -----------------------------------------------------
        out_ind = calc.detect_outliers(te_flagged.values)
        clean_te = te_flagged.copy()
        if len(out_ind) > 0:
            clean_te.iloc[out_ind] = np.nan

        # Remove 30s after LP sweep (20 data points) --------------------------
        sweep_idx = np.where(efi_lat['lp_flag'] == 9)[0]
        lp_pt = lp_s * freq + 1
        for sp in sweep_idx:
            if sp + lp_pt < len(clean_te):
                clean_te.iloc[sp: sp + lp_pt] = np.nan
            else:
                clean_te.iloc[sp:] = np.nan

        # flag te jumps -------------------------------------------------------
        te_jump = [0] * len(clean_te)
        te_jv = np.array([np.nan]* len(clean_te))
        for t in range(len(clean_te) - 1):
            t_og = clean_te.iloc[t]
            t_next = clean_te.iloc[t + 1]
            if abs(t_og - t_next) > 500:
                te_jump[t] = 1
            te_jv[t] = abs(t_og - t_next)

        jump_all[idx_all] = te_jump
        te_jump_b4[idx_all] = te_jv

        # save nan info -------------------------------------------------------

        nan_mask = np.isnan(clean_te)

        # Median Filter -------------------------------------------------------
        window = 10 * freq  # 10 seconds, 2 Hz
        te_meds = calc.moving_median(clean_te.values, window)

        # replace nans --------------------------------------------------------
        te_meds[nan_mask.values] = np.nan

        # update full list
        te_meds_all[idx_all] = te_meds

        te_jv_af = np.array([np.nan]* len(te_meds))
        for t in range(len(clean_te) - 1):
            t_og = te_meds[t]
            t_next = te_meds[t + 1]
            te_jv_af[t] = abs(t_og - t_next)

        te_jump_af[idx_all] = te_jv_af



    return (te_meds_all, dte_adj_all, te_sweep_all, jump_all,
            te_jump_b4, te_jump_af)


def process_Te_OLD(EFI_data, satellite, lp_s=30, freq=2):
    """ process EFI electron temperature by calibrating the Te, removing the
        bad flags, removing the outliers, and then median filtering
    Parameters
    ----------
    EFI_data : pandas DataFrame
        DataFrame of Swarm EFI data created by Load_Swarm.py
    satellite : str
        Swarm Satellite for calibration purposes
    lp_s : int kwarg
        number of seconds after lp sweep to remove
        default 30s
        30 seconds decided observationally
        and confirmed by Shasha on 10/17/2025
    te_jump
    Returns
    -------
    te_meds : array-like
        median filtered Te that is calibrated, flag checked,
        and outlier removed that is same length as EFI_data['Te']
    NOTES
    -----
    Something to be considered is removing data that is +10s after LP sweep
    """

    # Calibrate Te
    Te_adj, dte_adj = Te_calibrate_2021(EFI_data, satellite)

    # Remove Bad Flags
    te_flagged = Te_adj
    te_rem = ((EFI_data["te_flag"] != 10) & (EFI_data["te_flag"] != 19)
              & (EFI_data["te_flag"] != 20))
    te_flagged[te_rem] = np.nan

    # Remove Outliers
    out_ind = calc.detect_outliers(te_flagged.values)

    clean_te = te_flagged.copy()
    if len(out_ind) > 0:
        clean_te.iloc[out_ind] = np.nan

    # Remove 10s after LP sweep (20 data points)
    sweep_idx = np.where(EFI_data['lp_flag'] == 9)[0]
    lp_pt = lp_s * freq + 1
    for sp in sweep_idx:
        if sp + lp_pt < len(clean_te):
            clean_te.iloc[sp: sp + lp_pt] = np.nan
        else:
            clean_te.iloc[sp:] = np.nan

    # save nan info
    nan_mask = np.isnan(clean_te)

    # Median Filter
    window = 10 * freq  # 10 seconds, 2 Hz
    te_meds = calc.moving_median(clean_te.values, window)

    not_sweep = (EFI_data["lp_flag"] != 9)
    te_sweep = te_meds.copy()
    te_sweep[not_sweep] = np.nan

    # replace nans
    te_meds[nan_mask] = np.nan

    return te_meds, dte_adj, te_sweep
