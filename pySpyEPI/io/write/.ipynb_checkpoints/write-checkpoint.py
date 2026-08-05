"""Write monthly files.
Updated 25 February 2026
alanahco

List Functions
--------------
build_stats_filename
write_stats
build_figname
write_fig
"""
import os
from pathlib import Path


def build_stats_filename(stime, obs, file_dir, filetype, sat, daily=True):
    """Build the filename and directory for daily EIA stat files.

    Parameters
    ----------
    stime : datetime
        day of desired file
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
        creates name for yearly files if False
        default True
    Returns
    -------
    date_dir : str
        Directory path in which file should exist
        for daily
            (e.g., 'file_dir/filetype/Y/m')
        for monthly
            (e.g., 'file_dir/{filetype}_monthly/Y/m')
    fname : str
        Filename without directory

    """
    # Build the directory path
    if daily:
        date_dir = os.path.join(file_dir, filetype, stime.strftime('%Y'),
                                stime.strftime('%m'))
    else:
        date_dir = os.path.join(file_dir, filetype)

    # Build the filename
    if daily:
        t_str = stime.strftime('%Y%m%d')
    else:
        t_str = stime.strftime('%Y')

    fname = f'{filetype}_{t_str}_{obs}{sat}_ascii.txt'

    return date_dir, fname


def write_stats(stat_data, stime, obs, file_dir, filetype, sat,
                daily=True):
    """Write the daily statistics file for model-data comparisons.

    Parameters
    ----------
    stat_data : pd.DataFrame
        DataFrame that includes all statistics to save
    stime : datetime
        day of desired file
    model : str
        Case-sensitive name of model requested (e.g., 'NIMO', 'PyIRI').
    obs : str
        Name of data set requested (e.g., 'SWARM', 'MADRIGAL')
    file_dir : str
        File directory, if it does not exist it will use the current directory
    kwargs : dict
        Optional kwargs by data type.  Includes 'mad_lon', which expects
        longitudes of either -90 deg E or 60 deg E for Madrigal data.

    """
    # Test the output directory
    Path(file_dir).mkdir(parents=True, exist_ok=True)

    # Build the output directory path and filename
    date_dir, fname = build_stats_filename(stime, obs, file_dir, filetype, sat,
                                           daily=daily)

    # Ensure the directory exists
    Path(date_dir).mkdir(parents=True, exist_ok=True)
    save_file = os.path.join(date_dir, fname)

    # Create the custom header row with a hashtag
    header_line = '#{:s}\n'.format('\t'.join(stat_data.columns))

    # Write the header to the file
    with open(save_file, 'w') as fout:
        fout.write(header_line)

    # Append the DataFrame data without the header and index
    stat_data.to_csv(save_file, sep='\t', index=False,
                     na_rep='NaN', header=False, mode='a', encoding='ascii')

    return


def build_figname(stime, pass_id, fig_dir, obs, sat):
    """Build figure name.
    Parameters
    ----------
    stime : datetime
        day of desired file
    pass_id : str
        pass_id created in build.find_cand
    fig_dir : string
        figure directory for saving
    obs : str
        Name of data set requested (e.g., 'SWARM', 'MADRIGAL')
    sat : string
        satellite letter or number
    """

    date_dir = os.path.join(fig_dir, stime.strftime('%Y'),
                            stime.strftime('%m'), stime.strftime('%Y%m%d'))

    # {obs}{sat}detection_{pass_id}.png
    figname = f'{sat}_{pass_id}.png'

    return date_dir, figname


def write_fig(stime, fig, pass_id, obs, sat, fig_dir, date_dir=None,
              figname=None):
    """Write figure png files.
    Parameters
    ----------
    stime : datetime
        day of desired file
    fig : figure handle
        figure to save
    pass_id : str
        pass_id created in build.find_cand
    fig_dir : string
        figure directory for saving
    obs : str
        Name of data set requested (e.g., 'SWARM', 'MADRIGAL')
    sat : string
        satellite letter or number
    """

    if date_dir is None:
        date_dir, figname = build_figname(stime, pass_id, fig_dir, obs, sat)

    # make the directory
    Path(date_dir).mkdir(parents=True, exist_ok=True)

    save_file = os.path.join(date_dir, figname)

    # if file path does not exist, create a new figure
    if not os.path.exists(save_file):
        # save figure in save figure directory as png
        fig.savefig(save_file, bbox_inches='tight', dpi=300)
    return
