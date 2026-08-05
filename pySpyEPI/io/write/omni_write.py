"""Write solar wind files.
Created : 14 March 2026
alanahco
List Functions
--------------
sw_year_save
sw_epi
nearest

Future updates:
test sw_year_save
make a if made function for sw_epi
"""

from dateutil.relativedelta import relativedelta
from pySpyEPI.io.load import omni_load
from pySpyEPI.io.write import write
import pandas as pd
import datetime as dt
import numpy as np
from pySpyEPI.detection.stats.solarwind_stats import swind_year
from pySpyEPI.io.load import det_load


def sw_year_save(day_st, day_ed, file_dir, param='SYM_H', mval=-30, below=True,
                 param_lab='symh'):
    """Get yearly solar wind data and save based on solarwind_stats.swind_year.
    day_st : datetime
        start of period
    day_ed : datetime
        end of period
    param : string
        parameter of focus
        default SYM_H
    mval : float
        value above or below which to save data on
        default -20
    below : boolean
        True: param <= mval (default)
        False : param >= mval
    param_lab : string
        label for dataframe
    """

    # get solar wind data
    stat_data = swind_year(day_st, day_ed, param=param, mval=mval, below=below,
                           param_lab=param_lab)

    # save solar wind data
    sat = '_omni'
    obs = param_lab
    filetype = 'solarwind'
    write.write_stats(stat_data, day_st, obs, file_dir, filetype, sat,
                      daily=False)
    return stat_data


def sw_epi(date_array, sat, epi_dir, year_dir):
    """Iterate through EPI files and save info such as how long after dip
    in AL or SYM H was the EPI. What was the min? How close was the time to the
    minimum
    Also use swind_year data for this and include B field values where available
    And HP30 if possible
    """

    # ALL columns
    columns = ['center_time', 'epi_by_gsm', 'epi_bz_gsm', 'epi_flow_speed',
               'epi_vx', 'epi_vy', 'epi_vz', 'epi_proton_density',
               'epi_pressure', 'epi_efield', 'epi_plasma_beta', 'epi_bsn_x',
               'epi_bsn_y', 'epi_bsn_z', 'epi_ae_index', 'epi_al_index',
               'epi_au_index', 'epi_symh', 'time_sym_dip_b4',
               'time_symh_min_b4', 'symh_min_b4', 'time_sym_rec_b4',
               'time_sym_dip_af', 'time_symh_min_af', 'symh_min_af',
               'time_sym_rec_af', 'time_al_dip_b4', 'time_al_min_b4',
               'al_min_b4', 'time_al_rec_b4', 'time_al_dip_af',
               'time_al_min_af', 'al_min_af', 'time_al_rec_af']

    # Download OMNI data for full date range
    day_st = date_array[0] - dt.timedelta(days=1)
    day_ed = date_array[-1] + dt.timedelta(days=1)
    OMNI_data = omni_load.get_OMNI(day_st, day_ed)

    # Update dataframe
    om_df = OMNI_data.reset_index()          # moves index into a column
    om_df = om_df.rename(columns={'index': 'Epoch'})  # rename if needed
    om_df['Epoch'] = pd.to_datetime(om_df['Epoch'])

    omni_str = ['BY_GSM', 'BZ_GSM', 'flow_speed', 'Vx', 'Vy', 'Vz',
                'proton_density', 'Pressure', 'E', 'Beta', 'BSN_x',
                'BSN_y', 'BSN_z', 'AE_INDEX', 'AL_INDEX', 'AU_INDEX', 'SYM_H']
    s_str = ['epi_by_gsm', 'epi_bz_gsm', 'epi_flow_speed', 'epi_vx', 'epi_vy',
             'epi_vz', 'epi_proton_density', 'epi_pressure', 'epi_efield',
             'epi_plasma_beta', 'epi_bsn_x', 'epi_bsn_y', 'epi_bsn_z',
             'epi_ae_index', 'epi_al_index', 'epi_au_index', 'epi_symh']

    # load yearly omni files
    date_1 = date_array[0] - relativedelta(years=1)
    date_2 = date_array[-1] + relativedelta(years=1)

    d1_str = date_1.strftime("%Y-%m-%d")
    d2_str = date_2.strftime("%Y-%m-%d")

    date_omni = pd.date_range(start=d1_str, end=d2_str, freq='YS')

    sym = det_load.load_detections(date_omni, 'symh', year_dir, 'solarwind',
                                   '_omni', daily=False)
    al_ind = det_load.load_detections(date_omni, 'al_index', year_dir,
                                      'solarwind', '_omni', daily=False)

    # iterate through dates
    for dayt in date_array:

        sw_save = pd.DataFrame(columns=columns)

        # open files one day at a time
        d_str = dayt.strftime("%Y-%m-%d")
        date_range = pd.date_range(start=d_str, end=d_str)

        # load the epi file
        epi_df = det_load.load_detections(date_range, 'swarm', epi_dir,
                                          'detection', sat)

        if len(epi_df) == 0:
            continue

        for ei in range(len(epi_df)):
            ct = epi_df['center_time'].iloc[ei]
            sw_save.at[ei, 'center_time'] = ct

            # use yearly omni files to get info on symh and alindex
            # SYM H
            sb4 = nearest(sym['time_start'], ct, direction='before')
            saf = nearest(sym['time_start'], ct, direction='after')
            b4_iloc = np.where(sb4 == sym['time_start'])[0][0]
            af_iloc = np.where(saf == sym['time_start'])[0][0]

            # before
            sw_save.at[ei, 'time_sym_dip_b4'] = sb4
            sw_save.at[ei, 'time_symh_min_b4'] = sym['time_min'].iloc[b4_iloc]
            sw_save.at[ei, 'time_sym_rec_b4'] = sym['time_end'].iloc[b4_iloc]
            sw_save.at[ei, 'symh_min_b4'] = sym['symh_min'].iloc[b4_iloc]

            # after
            sw_save.at[ei, 'time_sym_dip_af'] = saf
            sw_save.at[ei, 'time_symh_min_af'] = sym['time_min'].iloc[af_iloc]
            sw_save.at[ei, 'time_sym_rec_af'] = sym['time_end'].iloc[af_iloc]
            sw_save.at[ei, 'symh_min_af'] = sym['symh_min'].iloc[af_iloc]

            # AL
            alb4 = nearest(al_ind['time_start'], ct, direction='before')
            alaf = nearest(al_ind['time_start'], ct, direction='after')
            b4_iloc = np.where(alb4 == al_ind['time_start'])[0][0]
            af_iloc = np.where(alaf == al_ind['time_start'])[0][0]

            # before
            sw_save.at[ei, 'time_al_dip_b4'] = alb4
            sw_save.at[ei, 'time_al_min_b4'] = al_ind['time_min'].iloc[b4_iloc]
            sw_save.at[ei, 'time_al_rec_b4'] = al_ind['time_end'].iloc[b4_iloc]
            sw_save.at[ei, 'al_min_b4'] = al_ind['al_index_min'].iloc[b4_iloc]

            # after
            sw_save.at[ei, 'time_al_dip_af'] = alaf
            sw_save.at[ei, 'time_al_min_af'] = al_ind['time_min'].iloc[af_iloc]
            sw_save.at[ei, 'time_al_rec_af'] = al_ind['time_end'].iloc[af_iloc]
            sw_save.at[ei, 'al_min_af'] = al_ind['al_index_min'].iloc[af_iloc]

            # limit omni data by halfwidth
            om_epi = om_df.copy()
            om_near = nearest(om_epi['Epoch'], ct, direction='both')
            i_om = np.where(om_near == om_epi['Epoch'])[0][0]
            om_epi = om_epi.iloc[i_om]

            # iterate through all omni save strings
            for oi, ostring in enumerate(omni_str):
                # save max and min of each value
                sw_save.at[ei, s_str[oi]] = om_epi[ostring]

        # save sw_save
        write.write_stats(sw_save, dayt, 'swarm', epi_dir, 'omni', sat,
                          daily=True)


def nearest(items, pivot, direction="both"):
    """Find nearest time to pivot.

    Parameters
    ----------
    items : datetime list
    pivot : datetime
    direction: string
        'before', 'after', or 'both' if you want the closest time before
        pivot, after pivot, or nearest time in whole list
    """
    if direction == "before":
        items = [x for x in items if x <= pivot]
        return min(items, key=lambda x: pivot - x)

    elif direction == "after":
        items = [x for x in items if x >= pivot]
        return min(items, key=lambda x: x - pivot)

    else:  # both (original behavior)
        return min(items, key=lambda x: abs(x - pivot))
