"""Swarm Loading Functions.
Upated 25 March 2026
alanahco

List Functions
--------------
load_IRR
load_EFI
load_all
quick_load
"""
import pandas as pd
import numpy as np
import glob
import os
from pySpyEPI.utils import calc
from pySpyEPI.utils import coords
from pySpyEPI.utils.sats import swarm_utils as swu


def load_IRR(start_date, end_date, satellite, fdir):
    """Load Swarm IRR dataset.

    Parameters
    ----------
    start_date : datetime
        starting time
    end_date : dateteim
        ending time
    satellite : str
        'A', 'B', or 'C'
    fdir : str
        File directory string

    Returns
    -------
    df2 : Dataframe
        Dataframe of swarm data
        columns: Time, Ne, bNe, fNe, Grad_Ne_at_100km, Grad_Ne_at_50km,
        Grad_Ne_at_20km, delta_Ne10s, delta_Ne20s, delta_Ne40s, ROD, RODI10,
        RODI20, mVTEC, mROT, Te, Ne_flag, Ionosphere_region_flag, PCP_flag,
        mROTI10, mROTI20, IBI_flag, IPIR_Index, TEC_STD,
        Latitude, Longitude, Mag_Lat, Mag_Lon, and LT
    Notes
    -----
    The returned delta_NeXs are scaled by Ne
    """
    # Set Directory
    base_path = f"{fdir}/IPD/Sat_{satellite}"
    date_yr = start_date.strftime('%Y')
    date_str = start_date.strftime('%Y%m%d')

    dir_path = f"{base_path}/{date_yr}/{date_str}/"
    search_pattern = os.path.join(dir_path, "*.cdf")
    filename = 'swarm_file'

    if len(glob.glob(search_pattern)) > 0:
        filename = glob.glob(search_pattern)[0]

    if os.path.exists(filename):
        variables = [
            "Timestamp", "Latitude", "Longitude", "Ne", "Background_Ne",
            "Foreground_Ne", "Grad_Ne_at_100km", "Grad_Ne_at_50km",
            "Grad_Ne_at_20km", "delta_Ne10s", "delta_Ne20s", "delta_Ne40s",
            "ROD", "RODI10s", "RODI20s", "mVTEC", "mROT", "Te",
            "Ne_quality_flag", "Ionosphere_region_flag", "PCP_flag",
            "mROTI10s", "mROTI20s", "IBI_flag", "IPIR_index", "TEC_STD",
            "Radius"
        ]

        data, cdf_file = swu.load_cdf_data(filename, variables)
        lat, lon, time = swu.extract_common_coords(cdf_file, satellite)
        in_time = pd.DataFrame(time)

        local_time = coords.longitude_to_local_time(lon, time)

        # calculate altitude and convert to km
        alt = data['Radius'] / 1000 - coords.earth_radius(lat) / 1000
        rad_km = np.nanmedian(alt)

        # Also compute qd coords
        glat, glon, l_shells = coords.compute_qd_coords(
            lat, lon, rad_km, in_time[0])

        # set up dataframe
        df = pd.DataFrame({
            "time": time, "ne": data["Ne"], "bne": data["Background_Ne"],
            "fne": data["Foreground_Ne"],
            "grad_ne_at_100km": data["Grad_Ne_at_100km"],
            "grad_ne_at_50km": data["Grad_Ne_at_50km"],
            "grad_ne_at_20km": data["Grad_Ne_at_20km"], "rod": data["ROD"],
            "rodi10": data["RODI10s"], "rodi20": data["RODI20s"],
            "mvtec": data["mVTEC"], "mrot": data["mROT"], "te": data["Te"],
            "ionosphere_region_flag": data["Ionosphere_region_flag"],
            "pcp_flag": data["PCP_flag"], "mroti10": data["mROTI10s"],
            "mroti20": data["mROTI20s"], "ibi_flag": data["IBI_flag"],
            "ipir_index": data["IPIR_index"], "tec_std": data["TEC_STD"],
            "lat": lat, "lon": lon, "mlat": glat, "mlon": glon,
            "lt": local_time, "l_shell": l_shells, "altitude": alt
        })

        # Create Ne Flag
        # First add LT_hour
        df["lt_hour"] = df["lt"].apply(
            lambda tim: tim.hour + tim.minute / 60 + tim.second / 3600)
        flg = swu.create_ne_flag(df)
        df['ne_flag'] = flg

        # calculate ne delta nes with flag removed
        # calculate your own delta_NeXs instead of using what is given
        df_calc = df.copy()
        df_calc.loc[df_calc['ne_flag'] > 2, 'ne'] = np.nan
        delta_Ne10s = calc.calc_delta_ne(df_calc["ne"], t=10, freq=1)
        delta_Ne20s = calc.calc_delta_ne(df_calc["ne"], t=20, freq=1)
        delta_Ne40s = calc.calc_delta_ne(df_calc["ne"], t=40, freq=1)

        df['delta_ne10s'] = delta_Ne10s
        df['delta_ne20s'] = delta_Ne20s
        df['delta_ne40s'] = delta_Ne40s

        # set delta ne bad ne flags to nan
        df.loc[df['ne_flag'] > 2, 'delta_ne10s'] = np.nan
        df.loc[df['ne_flag'] > 2, 'delta_ne20s'] = np.nan
        df.loc[df['ne_flag'] > 2, 'delta_ne40s'] = np.nan
        df.loc[df['ne_flag'] > 2, 'ne'] = np.nan

        # cut off by time
        df2 = df[(df['time'] > start_date) & (df['time'] < end_date)]
    else:
        df2 = pd.DataFrame()

    return df2


def load_EFI(start_date, end_date, satellite, fdir):
    """Load Swarm EFI dataset.

    Parameters
    ----------
    start_date : datetime
        starting time
    end_date : dateteim
        ending time
    satellite : str
        'A', 'B', or 'C'
    fdir : str
        File directory string

    Returns
    -------
    df2 : Dataframe
        Dataframe of swarm data
        columns: Time, Ne, Ne_error, Ne_flag, Te, Te_error, Te_flag, LP_flag,
        Latitude, Longitude, Mag_Lat, Mag_Lon, and LT
    """
    # Set Directory
    base_path = f"{fdir}/EFI/Sat_{satellite}"

    date_yr = start_date.strftime('%Y')
    date_str = start_date.strftime('%Y%m%d')

    dir_path = f"{base_path}/{date_yr}/{date_str}/"
    search_pattern = os.path.join(dir_path, "*.cdf")
    filename = 'swarm_file'
    if len(glob.glob(search_pattern)) > 0:
        filename = glob.glob(search_pattern)[0]

    if os.path.exists(filename):
        if ('_0701' in filename) | ('_0702' in filename):  # updated EFI files
            variables = ['Timestamp', 'Latitude', 'Longitude', 'N_ion',
                         'N_elec', 'T_elec', 'dN_ion', 'N_ion_error',
                         'Flags_T_elec', 'N_elec_error', 'T_elec_error',
                         'Flags_N_elec', 'Flags_N_ion', 'dT_elec', 'Flagbits1',
                         'Flagbits2', 'Gamma1', 'Gamma2', 'Flags_LP', 'Radius',
                         'Flags_Vs']
            n_str = 'N_ion'
            n_err_str = 'N_ion_error'
            f_n_str = 'Flags_N_ion'
            t_str = 'T_elec'
            t_err_str = 'T_elec_error'
            f_t_str = 'Flags_T_elec'
            f_lp_str = 'Flags_LP'
            if ('_0701' in filename):
                ftype = '0701'
            elif ('_0702' in filename):
                ftype = '0702'
        else:  # old file format
            variables = ["Timestamp", "Latitude", "Longitude", "Ne",
                         "Ne_error", "Te", "Te_error", "Flags_Ne", "Flags_Te",
                         "Flags_LP", "Radius", 'Flags_Vs']
            n_str = 'Ne'
            n_err_str = 'Ne_error'
            f_n_str = 'Flags_Ne'
            t_str = 'Te'
            t_err_str = 'Te_error'
            f_t_str = 'Flags_Te'
            f_lp_str = 'Flags_LP'
            ftype = '0602'

        # load file
        data, cdf_file = swu.load_cdf_data(filename, variables)

        lat, lon, time = swu.extract_common_coords(cdf_file, satellite)
        in_time = pd.DataFrame(time)

        alt = data['Radius'] / 1000 - coords.earth_radius(lat) / 1000
        rad_km = np.nanmedian(alt)

        # Also compute qd coords
        glat, glon, L_shells = coords.compute_qd_coords(lat, lon, rad_km,
                                                        in_time[0])
        local_time = coords.longitude_to_local_time(lon, time)

        file_form = [ftype] * len(glat)

        df = pd.DataFrame({
            "time": time, "ne": data[n_str], "ne_error": data[n_err_str],
            "flag_ne": data[f_n_str], "te": data[t_str],
            "te_error": data[t_err_str], "te_flag": data[f_t_str],
            "lp_flag": data[f_lp_str], "lat": lat, "lon": lon,
            "mlat": glat, "mlon": glon, "lt": local_time, "fnum": file_form
        })

        # process and adjust Te
        te_pro, dt_adj, te_sweep, te_jump, tb4, taf = swu.process_Te(df,
                                                                     satellite)
        df['te_adj'] = te_pro

        # save calibration error
        df['dte_adj'] = dt_adj

        # save sweep tes
        df['te_sweep'] = te_sweep

        # save jump flag
        df['te_jump'] = te_jump

        # save jump flag
        df['te_dif_b4'] = tb4

        # save jump flag
        df['te_dif_af'] = taf

        # limit by time
        df2 = df[(df['time'] > start_date) & (df['time'] < end_date)]

    else:
        df2 = pd.DataFrame()

    return df2


def load_all(start_date, end_date, satellite, fdir, IRR=True, EFI=True):
    """Load all four datasets.

    Parameters
    ----------
    start_date : datetime
        starting time
    end_date : dateteim
        ending time
    satellite : str
        'A', 'B', or 'C'
    fdir : str
        File directory string
    IRR : bool kwarg
        True default- IRR data will be retreived
    EFI : bool kwarg
        True default- EFI data will be retreived

    Returns
    -------
    data : Dictionary
         dictionary containing all True datasets
         as long as data is present
    IRR : bool
        if IRR data is present in data
    EFI : bool
        if EFI data is present in data

    Notes
    -----
    Printed Warnings if a dataset does not have any parameters present
    for desired time period
    """

    data = {}

    # Establish midnight of end day
    midnight = end_date.replace(hour=0, minute=0, second=0, microsecond=0)
    # open IRR if True
    if IRR:
        try:
            IRR_data = load_IRR(start_date, end_date, satellite, fdir=fdir)

            # if start_date.day != end_date.day
            # Get next day too-- midnight to end_date
            if start_date.day != end_date.day:
                if (midnight != end_date):
                    IRR_data2 = load_IRR(midnight, end_date, satellite,
                                         fdir=fdir)
                    IRR_data = pd.concat([IRR_data, IRR_data2],
                                         ignore_index=True)

            # Check length of IRR data
            len_ch = len(IRR_data['time'])

            if len_ch > 0:
                # Add to dictionary
                data['IRR'] = IRR_data
            else:
                IRR = False
                print(f'Warning: No IRR data for {start_date}')
        except Exception:
            IRR = False
            print(f'Warning: No IRR File for {start_date}')

    # open EFI if True
    if EFI:
        try:
            EFI_data = load_EFI(start_date, end_date, satellite, fdir=fdir)

            # if start_date.day != end_date.day
            # Get next day too-- midnight to end_date
            if start_date.day != end_date.day:
                if (midnight != end_date):
                    EFI_data2 = load_EFI(midnight, end_date, satellite,
                                         fdir=fdir)
                    EFI_data = pd.concat([EFI_data, EFI_data2],
                                         ignore_index=True)

            # Check length of EFI data
            len_ch = len(EFI_data['time'])

            if len_ch > 0:
                # Add to dictionary
                data['EFI'] = EFI_data
            else:
                EFI = False
                print(f'Warning: No EFI data for {start_date}')
        except Exception:
            EFI = False
            print(f'Warning: No EFI File for {start_date}')

    return data, IRR, EFI


def quick_load(start_date, end_date, satellite, fdir):
    """Load Swarm IRR dataset.

    Parameters
    ----------
    start_date : datetime
        starting time
    end_date : dateteim
        ending time
    satellite : str
        'A', 'B', or 'C'
    fdir : str
        File directory string

    Returns
    -------
    df2 : Dataframe
        Dataframe of swarm data
        columns: Time, Ne, bNe, fNe, Grad_Ne_at_100km, Grad_Ne_at_50km,
        Grad_Ne_at_20km, delta_Ne10s, delta_Ne20s, delta_Ne40s, ROD, RODI10,
        RODI20, mVTEC, mROT, Te, Ne_flag, Ionosphere_region_flag, PCP_flag,
        mROTI10, mROTI20, IBI_flag, IPIR_Index, TEC_STD,
        Latitude, Longitude, Mag_Lat, Mag_Lon, and LT
    Notes
    -----
    The returned delta_NeXs are scaled by Ne
    """
    # Set Directory
    base_path = f"{fdir}/IPD/Sat_{satellite}"
    date_yr = start_date.strftime('%Y')
    date_str = start_date.strftime('%Y%m%d')

    dir_path = f"{base_path}/{date_yr}/{date_str}/"
    search_pattern = os.path.join(dir_path, "*.cdf")
    filename = 'swarm_file'

    if len(glob.glob(search_pattern)) > 0:
        filename = glob.glob(search_pattern)[0]

    if os.path.exists(filename):
        variables = [
            "Timestamp", "Latitude", "Longitude", "Ne", "Background_Ne",
            "Foreground_Ne", "Grad_Ne_at_100km", "Grad_Ne_at_50km",
            "Grad_Ne_at_20km", "delta_Ne10s", "delta_Ne20s", "delta_Ne40s",
            "ROD", "RODI10s", "RODI20s", "mVTEC", "mROT", "Te",
            "Ne_quality_flag", "Ionosphere_region_flag", "PCP_flag",
            "mROTI10s", "mROTI20s", "IBI_flag", "IPIR_index", "TEC_STD",
            "Radius"
        ]

        data, cdf_file = swu.load_cdf_data(filename, variables)
        lat, lon, time = swu.extract_common_coords(cdf_file, satellite)

        local_time = coords.longitude_to_local_time(lon, time)

        # calculate altitude and convert to km
        alt = data['Radius'] / 1000 - coords.earth_radius(lat) / 1000

        # set up dataframe
        df = pd.DataFrame({
            "time": time, "ne": data["Ne"], "bne": data["Background_Ne"],
            "fne": data["Foreground_Ne"],
            "grad_ne_at_100km": data["Grad_Ne_at_100km"],
            "grad_ne_at_50km": data["Grad_Ne_at_50km"],
            "grad_ne_at_20km": data["Grad_Ne_at_20km"], "rod": data["ROD"],
            "rodi10": data["RODI10s"], "rodi20": data["RODI20s"],
            "mvtec": data["mVTEC"], "mrot": data["mROT"], "te": data["Te"],
            "ionosphere_region_flag": data["Ionosphere_region_flag"],
            "pcp_flag": data["PCP_flag"], "mroti10": data["mROTI10s"],
            "mroti20": data["mROTI20s"], "ibi_flag": data["IBI_flag"],
            "ipir_index": data["IPIR_index"], "tec_std": data["TEC_STD"],
            "lat": lat, "lon": lon,
            "lt": local_time, "altitude": alt
        })

        # Create Ne Flag
        # First add LT_hour
        df["lt_hour"] = df["lt"].apply(
            lambda tim: tim.hour + tim.minute / 60 + tim.second / 3600)
        flg = swu.create_ne_flag(df)
        df['ne_flag'] = flg

        # calculate ne delta nes with flag removed
        # calculate your own delta_NeXs instead of using what is given
        df_calc = df.copy()
        df_calc.loc[df_calc['ne_flag'] > 2, 'ne'] = np.nan
        delta_Ne10s = calc.calc_delta_ne(df_calc["ne"], t=10, freq=1)
        delta_Ne20s = calc.calc_delta_ne(df_calc["ne"], t=20, freq=1)
        delta_Ne40s = calc.calc_delta_ne(df_calc["ne"], t=40, freq=1)

        df['delta_ne10s'] = delta_Ne10s
        df['delta_ne20s'] = delta_Ne20s
        df['delta_ne40s'] = delta_Ne40s

        # set delta ne bad ne flags to nan
        df_calc.loc[df_calc['ne_flag'] > 2, 'delta_ne10s'] = np.nan
        df_calc.loc[df_calc['ne_flag'] > 2, 'delta_ne20s'] = np.nan
        df_calc.loc[df_calc['ne_flag'] > 2, 'delta_ne40s'] = np.nan
        df_calc.loc[df_calc['ne_flag'] > 2, 'ne'] = np.nan

        # cut off by time
        df2 = df[(df['time'] > start_date) & (df['time'] < end_date)]
    else:
        df2 = pd.DataFrame()

    return df2
