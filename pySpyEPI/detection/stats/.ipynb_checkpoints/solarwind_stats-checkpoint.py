"""Save Solar Wind stats for each day.
Updated on 27 February 2026
alanahco

List Functions
--------------
swind_day
swind_epi
"""
import pandas as pd
from pySpyEPI.barrel import barrel_utils as bar_uts
from pySpyEPI.depi_io.load import omni_load


def swind_year(day_st, day_ed, param='SYM_H', mval=-30, below=True,
               param_lab='sym'):
    """Create solar wind df for specified param when dips below certain value.

    Parameters
    ----------

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
    OMNI_data = omni_load.get_OMNI(day_st, day_ed)

    # Update dataframe
    df = OMNI_data.reset_index()          # moves index into a column
    df = df.rename(columns={'index': 'Epoch'})  # rename if needed
    df['Epoch'] = pd.to_datetime(df['Epoch'])

    if below:
        df_lim = df[df[param] <= mval]
        meas = 'min'
    else:
        df_lim = df[df[param] >= mval]
        meas = 'max'

    # get gap indices
    sw_idx = df_lim.index.values
    gaps = bar_uts.find_all_gaps(sw_idx)

    columns = ['time_start', 'time_end', f'time_{meas}', f'{param_lab}_{meas}']
    sw_save = pd.DataFrame(columns=columns)

    for i in range(len(gaps) - 1):

        # segment dataframe
        sw_df = df_lim.iloc[gaps[i]: gaps[i + 1]]

        # save variables
        sw_save.at[i, 'time_start'] = sw_df['Epoch'].iloc[0]
        sw_save.at[i, 'time_end'] = sw_df['Epoch'].iloc[-1]

        # save min or max depending on below
        if meas == 'min':
            mi = sw_df[param].argmin()
            mv = sw_df[param].min()
        else:
            mi = sw_df[param].argmax()
            mv = sw_df[param].max()

        sw_save.at[i, f'time_{meas}'] = sw_df['Epoch'].iloc[mi]
        sw_save.at[i, f'{param_lab}_{meas}'] = mv

    return sw_save
