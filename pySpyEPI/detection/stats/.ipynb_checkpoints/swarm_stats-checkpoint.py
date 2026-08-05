"""Save EPI info Swarm data.
Created 19 February 2026
alanahco
List Functions
--------------
IRR_stats
EFI_stats
"""
import pandas as pd
import numpy as np
from pySpyEPI.utils import calc


def IRR_stats(IRR_data, stats_df):
    """Stats from IRR files for detections.
    Parameters
    ----------
    IRR_data : DataFrame
        Swarm data from IRR
    stats_df : DataFrame
        dataframe containing EPI info created by epi_stats
    pass_id : str
        pass_id created in build.find_cand or swarm_detections.swarm_detect
    Returns
    ------
    df : dataframe
        dataframe containing IRR dataset info for detections
        Includes maximum, minimum, median, and standard deviation of the
        following:
            Grad Ne (100 km, 50 km, 20 km)
            Delta Ne (10s, 20s, 40s)
            ROD, RODI (10s, 20s)
            mVTEC
            mROT, mROTI (10s, 20s)
            TEC_STD
        Additionally:
            IBI (max value)
            IPIR (max value)
            Ionosphere Region Flag (IRF) at center and edges
            pass_id
    """

    # iterate through columns and create dataframe
    columns = ['center_time', 'gne100_max', 'gne100_min', 'gne100_med',
               'gne100_std', 'gne50_max', 'gne50_min', 'gne50_med',
               'gne50_std', 'gne20_max', 'gne20_min', 'gne20_med', 'gne20_std',
               'dne10_max', 'dne10_min', 'dne10_med', 'dne10_std', 'dne20_max',
               'dne20_min', 'dne20_med', 'dne20_std', 'dne40_max', 'dne40_min',
               'dne40_med', 'dne40_std', 'rod_max', 'rod_min', 'rod_med',
               'rod_std', 'rodi10_max', 'rodi10_min', 'rodi10_med',
               'rodi10_std', 'rodi20_max', 'rodi20_min', 'rodi20_med',
               'rodi20_std', 'mvtec_max', 'mvtec_min', 'mvtec_med',
               'mvtec_std', 'mrot_max', 'mrot_min', 'mrot_med', 'mrot_std',
               'mroti10_max', 'mroti10_min', 'mroti10_med', 'mroti10_std',
               'mroti20_max', 'mroti20_min', 'mroti20_med', 'mroti20_std',
               'tec_std_max', 'tec_std_min', 'tec_std_med', 'tec_std_std',
               'ibi_max', 'ipir_max', 'irf_center', 'irf_edge1', 'irf_edge2',
               'pass_id']

    df = pd.DataFrame(columns=columns)

    c_time = stats_df['center_time'].values
    l_time = stats_df['edge1_time'].values
    r_time = stats_df['edge2_time'].values

    # set up a lis tof variables from IRR dataset to iterate through
    var_strs = ['grad_ne_at_100km', 'grad_ne_at_50km',
                'grad_ne_at_20km', 'delta_ne10s', 'delta_ne20s', 'delta_ne40s',
                'rod', 'rodi10', 'rodi20', 'mvtec', 'mrot',
                'mroti10', 'mroti20', 'tec_std']

    save_strs = ['gne100', 'gne50', 'gne20', 'dne10', 'dne20', 'dne40', 'rod',
                 'rodi10', 'rodi20', 'mvtec', 'mrot', 'mroti10', 'mroti20',
                 'tec_std']

    # iterate through peaks
    for i in range(len(c_time)):

        # get positional arguments of edges and center
        xl = (IRR_data['time'] - l_time[i]).abs().argmin()
        xp = (IRR_data['time'] - c_time[i]).abs().argmin()
        xr = (IRR_data['time'] - r_time[i]).abs().argmin()

        df.at[i, 'center_time'] = IRR_data['time'].iloc[xp]

        # Variables for max, min, median, and standard deviation --------------
        for v_str, s_str in zip(var_strs, save_strs):
            var = IRR_data[v_str].iloc[xl:xr + 1]

            # set save strings
            max_str = f'{s_str}_max'
            min_str = f'{s_str}_min'
            med_str = f'{s_str}_med'
            std_str = f'{s_str}_std'

            # save stats
            df.at[i, max_str] = var.max()
            df.at[i, min_str] = var.min()
            df.at[i, med_str] = var.median()
            df.at[i, std_str] = var.std()

        # IBI -----------------------------------------------------------------
        ibi = IRR_data['ibi_flag'].iloc[xl:xr + 1]
        df.at[i, 'ibi_max'] = ibi.max()

        # IPIR ----------------------------------------------------------------
        ipir = IRR_data['ipir_index'].iloc[xl:xr + 1]
        df.at[i, 'ipir_max'] = ipir.max()

        # Ionosphere Region Flag ----------------------------------------------
        df.at[i, 'irf_center'] = IRR_data['ionosphere_region_flag'].iloc[xp]
        df.at[i, 'irf_edge1'] = IRR_data['ionosphere_region_flag'].iloc[xl]
        df.at[i, 'irf_edge2'] = IRR_data['ionosphere_region_flag'].iloc[xr]

    return df


def EFI_stats(EFI_data, stats_df):
    """Stats from EFI files for detections.

    Parameters
    ----------
    EFI_data : DataFrame
        Swarm data from EFI
    stats_df : DataFrame
        dataframe containing EPI info created by epi_stats
    pass_id : str
        pass_id created in build.find_cand or swarm_detections.swarm_detect

    Returns
    ------
    df : dataframe
        dataframe containing EFI dataset info for detections
        Includes:
            center te (3 point averge) and standard deviation
            edge 1 te (3 point averge) and standard deviation
            edge 2 te (3 point averge) and standard deviation
            te median for whole span
            te standard deviation for whole span
            percent nan for whole span
            lp flag (indicating high and low gains)
            time_lpsweep (time from last lp sweep and edge 1)
            te difference
            pass_id
            file number ('0701' or '0602')
    """

    # iterate through columns and create dataframe
    columns = ['center_time', 'center_te', 'center_te_std', 'edge1_te',
               'edge1_te_std', 'edge2_te', 'edge2_te_std', 'te_med', 'te_std',
               'te_percnan', 'lp_flag', 'te_jump', 'te_dif_b4', 'te_dif_af',
               'time_lpsweep', 'center_dte_max', 'edge1_dte_max',
               'edge2_dte_max', 'center_dte_min', 'edge1_dte_min',
               'edge2_dte_min', 'center_dte_med', 'edge1_dte_med',
               'edge2_dte_med', 'fnum', 'pass_id']

    df = pd.DataFrame(columns=columns)

    c_time = stats_df['center_time'].values
    l_time = stats_df['edge1_time'].values
    r_time = stats_df['edge2_time'].values
    passid = stats_df['pass_id'].values

    # get LP sweep indices
    lp_id = np.where(EFI_data["lp_flag"] == 9)[0]

    # iterate through peaks
    for i in range(len(c_time)):
        t_l = l_time[i]
        t_r = r_time[i]
        t_c = c_time[i]
        pass_id_use = passid[i]

        efi_peak_mask = ((EFI_data["time"] >= t_l)
                         & (EFI_data["time"] <= t_r))

        idx = np.where(efi_peak_mask)[0]

        if len(idx) != 0:

            # Include one before and one after if possible
            start = max(idx[0] - 1, 0)
            end = min(idx[-1] + 1, len(EFI_data) - 1)

            efi_use = EFI_data.iloc[start:end + 1]
            te_use = efi_use['te_adj'].values
            te_err = efi_use['dte_adj'].values

            # closest to center time
            p_id = (efi_use["time"] - t_c).abs().argmin()

            # time ------------------------------------------------------------
            df.at[i, 'center_time'] = efi_use['time'].iloc[p_id]
            df.at[i, 'pass_id'] = pass_id_use

            # Temp at peak, edges, and avg ------------------------------------
            df.at[i, 'center_te'] = calc.safe_nanmean(
                te_use[p_id - 1: p_id + 2])
            df.at[i, 'center_te_std'] = calc.safe_nanstd(
                te_use[p_id - 1: p_id + 2])
            df.at[i, 'center_dte_max'] = te_err[p_id - 1: p_id + 2].max()
            df.at[i, 'center_dte_min'] = te_err[p_id - 1: p_id + 2].min()
            df.at[i, 'center_dte_med'] = calc.safe_nanmedian(
                te_err[p_id - 1: p_id + 2])

            df.at[i, 'edge1_te'] = calc.safe_nanmean(te_use[0:3])
            df.at[i, 'edge1_te_std'] = calc.safe_nanstd(te_use[0:3])
            df.at[i, 'edge1_dte_max'] = te_err[0:3].max()
            df.at[i, 'edge1_dte_min'] = te_err[0:3].min()
            df.at[i, 'edge1_dte_med'] = calc.safe_nanmedian(te_err[0:3])

            df.at[i, 'edge2_te'] = calc.safe_nanmean(te_use[[-3, -2, -1]])
            df.at[i, 'edge2_te_std'] = calc.safe_nanstd(te_use[[-3, -2, -1]])
            df.at[i, 'edge2_dte_max'] = te_err[[-3, -2, -1]].max()
            df.at[i, 'edge2_dte_min'] = te_err[[-3, -2, -1]].min()
            df.at[i, 'edge2_dte_med'] = calc.safe_nanmedian(
                te_err[[-3, -2, -1]])

            df.at[i, 'te_med'] = calc.safe_nanmedian(te_use)
            df.at[i, 'te_std'] = calc.safe_nanstd(te_use)

            # LP flag for high gain and low gain probes------------------------
            lp_flag = efi_use['lp_flag']
            if np.any(lp_flag == 9):  # LP Sweep
                df.at[i, 'lp_flag'] = 9
            elif np.all(lp_flag == 1):  # high gain only
                df.at[i, 'lp_flag'] = 1
            elif np.all(lp_flag == 3):  # Low gain, may not be in 0702 files
                df.at[i, 'lp_flag'] = 3
            elif np.any(lp_flag == 7):  # only original, with sweep or error
                df.at[i, 'lp_flag'] = 7
            else:
                # LP flag = 5 in 0702 file indicates mixed probe
                # also if some are 1 and some are 3, then 2
                df.at[i, 'lp_flag'] = 2

            # check te_jump ---------------------------------------------------
            te_jump = efi_use['te_jump']

            if np.any(te_jump == 1):
                df.at[i, 'te_jump'] = 1
            else:
                df.at[i, 'te_jump'] = 0

            # get max jump before and corresponding after filtering
            # using maximum at b4 to find after

            if np.all(np.isnan(efi_use['te_dif_b4'])):
                df.at[i, 'te_dif_b4'] = np.nan
            else:
                df.at[i, 'te_dif_b4'] = np.nanmax(efi_use['te_dif_b4'])

            # Same index from b4
            if np.all(np.isnan(efi_use['te_dif_b4'])):
                df.at[i, 'te_dif_af'] = np.nan
            else:
                af_iloc = np.nanargmax(efi_use['te_dif_b4'])
                df.at[i, 'te_dif_af'] = efi_use['te_dif_af'].iloc[af_iloc]

            # percent nan -----------------------------------------------------
            df.at[i, 'te_percnan'] = np.isnan(te_use).mean() * 100

            # file  number
            df.at[i, 'fnum'] = efi_use['fnum'].iloc[p_id]

            # time between last lp sweep and depletion ------------------------
            lp_check = lp_id[(lp_id < (idx[0] - 1))]

            if len(lp_check) != 0:
                taflp = ((idx[0] - 1) - lp_check[-1]) / 2  # 2 hz cadence
            else:
                taflp = np.nan

            if np.any(efi_use['lp_flag'] == 9):
                taflp = 0

            df.at[i, 'time_lpsweep'] = taflp
        else:
            for cname in columns:
                df.at[i, cname] = np.nan
    # calculate and save te difference
    df['te_dif'] = df['center_te'] - (df['edge1_te'] + df['edge2_te']) / 2

    return df


def MAG_stats(MAG_dic, stats_df, pass_id):
    """Stats from MAG files for detections.

    Parameters
    ----------
    MAG_dic : DataFrame
        Swarm data from MAG files
    stats_df : DataFrame
        dataframe containing EPI info created by epi_stats
    pass_id : str
        pass_id created in build.find_cand or swarm_detections.swarm_detect

    Returns
    ------
    df : dataframe
        dataframe containing MAG dataset info for detections
        Includes maximum, minimum, median, and standard deviation of the
        following:
            B parallel (original, 2s and 60s HP)
            B zonal (original, 2s and 60s HP)
            B meridional (origininal, 2s and 60s HP)
    """

    # ititiate columns and create dataframe
    columns = ['center_time',
               'bpar_max', 'bpar_min', 'bpar_med', 'bpar_std',
               'bzon_max', 'bzon_min', 'bzon_med', 'bzon_std',
               'bmer_max', 'bmer_min', 'bmer_med', 'bmer_std',
               'bpar2s_max', 'bpar2s_min', 'bpar2s_med', 'bpar2s_std',
               'bzon2s_max', 'bzon2s_min', 'bzon2s_med', 'bzon2s_std',
               'bmer2s_max', 'bmer2s_min', 'bmer2s_med', 'bmer2s_std',
               'bpar60s_max', 'bpar60s_min', 'bpar60s_med', 'bpar60s_std',
               'bzon60s_max', 'bzon60s_min', 'bzon60s_med', 'bzon60s_std',
               'bmer60s_max', 'bmer60s_min', 'bmer60s_med', 'bmer60s_std',
               'pass_id']

    df = pd.DataFrame(columns=columns)

    c_time = stats_df['center_time'].values
    l_time = stats_df['edge1_time'].values
    r_time = stats_df['edge2_time'].values

    var_strs = ['B_MFA_hp2_par', 'B_MFA_hp2_zon', 'B_MFA_hp2_mer',
                'B_MFA_hp60_par', 'B_MFA_hp60_zon', 'B_MFA_hp60_mer']
    save_strs = ['bpar2s', 'bzon2s', 'bmer2s', 'bpar60s', 'bzon60s', 'bmer60s']

    # iterate through peaks
    for i in len(c_time):
        t_l = l_time[i]
        t_r = r_time[i]
        t_c = c_time[i]

        # closest to center time
        p_id = (MAG_dic["Time"] - t_c).abs().argmin()

        # time ----------------------------------------------------------------
        df.at[i, 'center_time'] = MAG_dic['Time'][p_id]

        # Mask MAG data
        MAG_mask = ((MAG_dic["Time"] >= t_l)
                    & (MAG_dic["Time"] <= t_r))
        filt_mag = (
            {key: value[MAG_mask] for key, value in MAG_dic.items()
             if isinstance(value, np.ndarray)})

        B_MFA = filt_mag['B_MFA']
        mfa_str = ['bpar', 'bzon', 'bmer']

        # Non filtered B MFA --------------------------------------------------
        for mf in range(3):
            b_var = B_MFA[:, mf]

            # set save strings
            max_str = f'{mfa_str[mf]}_max'
            min_str = f'{mfa_str[mf]}_min'
            med_str = f'{mfa_str[mf]}_med'
            std_str = f'{mfa_str[mf]}_std'

            if len(b_var) > 0:
                # save stats
                df.at[i, max_str] = b_var.max()
                df.at[i, min_str] = b_var.min()
                df.at[i, med_str] = np.nanmedian(b_var)
                df.at[i, std_str] = b_var.std()
            else:
                df.at[i, max_str] = np.nan
                df.at[i, min_str] = np.nan
                df.at[i, med_str] = np.nan
                df.at[i, std_str] = np.nan

        # Variables for max, min, median, and standard deviation---------------
        for v_str, s_str in zip(var_strs, save_strs):
            var = filt_mag[v_str]

            # set save strings
            max_str = f'{s_str}_max'
            min_str = f'{s_str}_min'
            med_str = f'{s_str}_med'
            std_str = f'{s_str}_std'

            if len(var) > 0:
                # save stats
                df.at[i, max_str] = var.max()
                df.at[i, min_str] = var.min()
                df.at[i, med_str] = np.nanmedian(var)
                df.at[i, std_str] = var.std()
            else:
                df.at[i, max_str] = np.nan
                df.at[i, min_str] = np.nan
                df.at[i, med_str] = np.nan
                df.at[i, std_str] = np.nan

    plen = len(stats_df)
    df['pass_id'] = [pass_id] * plen
    return df
