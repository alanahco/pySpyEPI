"""Plotting functions for barrel rolling.
Updated : 23 February 2026
alanahco

List Functions
--------------
barrel_ax
barrel2_ax
"""
import numpy as np


def barrel_ax(barrel_df, peaks, properties, ax, detrended=False, marker='o',
              mark_color="green", np_flag_str="ne_flag"):
    """ plot barrel info
    Parameters
    ----------
    barrel_df : dataframe
        dataframe created by barrel_roll.get_barrel_info and updated by
        barre_roll.just_barrel
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    ax : matplotlib axis
        ax for plotting
    detrended : kwarg bool
        default False, will plot ne with filtered ne
        True will plot detrended ne with peaks/prominences
    marker : kwarg str
        marker for peaks
        'o' is good for keep and 'x' is good for remove
    mark_color : kwarg str
        color string
        'green' is good for kept peaks
    np_flag_str : string
        column name in barrel_df for electron density flag
    Returns
    -------
    ax : axis with data plotted
    """
    # if we are not plotting the detrended and just ne with trend
    if not detrended:
        ax.plot(barrel_df['time'], barrel_df['ne'], color='#1f77b4')
        ax.plot(barrel_df['time'], barrel_df['barrel_ne'], color='orange')
        ax.scatter(barrel_df["time"][barrel_df[np_flag_str] > 2],
                   barrel_df["ne"][barrel_df[np_flag_str] > 2],
                   color="red", lw=1)
    else:
        # if we are plotting the detrended ne
        ax.plot(barrel_df['time'], barrel_df['dne_dd_filt'], color='#1f77b4',
                zorder=1)

        # plotting peaks will need to change when
        # we have confrimed and unconfirmed peaks
        ax.plot(barrel_df['time'].iloc[peaks],
                barrel_df['dne_dd_filt'].iloc[peaks], marker, color=mark_color,
                zorder=2)

        # plot prominences
        ymin = (barrel_df['dne_dd_filt'].iloc[peaks]
                - properties['prominences'])
        ymax = barrel_df['dne_dd_filt'].iloc[peaks]
        ax.vlines(x=barrel_df['time'].iloc[peaks], ymin=ymin, ymax=ymax,
                  color="C1", zorder=0)

        # plot widths
        xl = np.ceil(properties["left_ips"]).astype(int)
        xr = np.trunc(properties["right_ips"]).astype(int)
        ax.hlines(y=properties["width_heights"],
                  xmin=barrel_df["time"].iloc[xl],
                  xmax=barrel_df["time"].iloc[xr], color="C1", zorder=0)

    return ax


def barrel2_ax(barrel_up, barrel_lo, peaks, properties, ax, detrended=False,
               marker='o', mark_color="green", np_flag_str="ne_flag"):
    """Plot barrel upper and lower envelope info
    Parameters
    ----------
    barrel_up : dataframe
        dataframe created by barrel_roll.get_barrel_info and updated by
        barre_roll.just_barrel
        upper barrel
    barrel_lo : dataframe
        dataframe created by barrel_roll.get_barrel_info and updated by
        barre_roll.just_barrel
        lower barrel
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    ax : matplotlib axis
        ax for plotting
    detrended : kwarg bool
        default False, will plot ne with filtered ne
        True will plot detrended ne with peaks/prominences
    marker : kwarg str
        marker for peaks
        'o' is good for keep and 'x' is good for remove
    mark_color : kwarg str
        color string
        'green' is good for kept peaks
    np_flag_str : string
        column name in barrel_df for electron density flag
    Returns
    -------
    ax : axis with data plotted
    """
    # if we are not plotting the detrended and just ne with trend
    if not detrended:
        ax.plot(barrel_up['time'], barrel_up['ne'], color='gray')
        ax.scatter(barrel_up["time"][barrel_up[np_flag_str] > 2],
                   barrel_up["ne"][barrel_up[np_flag_str] > 2],
                   color="red", lw=1)
        ax.plot(barrel_up['time'], barrel_up['barrel_ne'], color='C0',
                linestyle='-')
        ax.plot(barrel_lo['time'], barrel_lo['barrel_ne'], color='C1',
                linestyle='--')
        ax.set_xlim([barrel_up['time'].iloc[0], barrel_up['time'].iloc[-1]])
    else:
        # if we are plotting the detrended ne
        ax.plot(barrel_up['time'], barrel_up['dne_dd_filt'], color='C0',
                zorder=1)

        # plotting peaks will need to change when
        # we have confrimed and unconfirmed peaks
        ax.plot(barrel_up['time'].iloc[peaks],
                barrel_up['dne_dd_filt'].iloc[peaks], marker, color=mark_color,
                zorder=2)

        # plot prominences
        ymin = (barrel_up['dne_dd_filt'].iloc[peaks]
                - properties['prominences'])
        ymax = barrel_up['dne_dd_filt'].iloc[peaks]
        ax.vlines(x=barrel_up['time'].iloc[peaks], ymin=ymin, ymax=ymax,
                  color="purple", zorder=0)

        # plot widths
        xl = np.ceil(properties["left_ips"]).astype(int)
        xr = np.trunc(properties["right_ips"]).astype(int)
        ax.hlines(y=properties["width_heights"],
                  xmin=barrel_up["time"].iloc[xl],
                  xmax=barrel_up["time"].iloc[xr], color="purple", zorder=0)

        ax.plot(barrel_lo['time'], barrel_lo['dne_dd_filt'], linestyle='--',
                color='C1')
        ax.set_xlim([barrel_up['time'].iloc[0], barrel_up['time'].iloc[-1]])

        ax.hlines(y=0, xmin=barrel_up["time"].iloc[0],
                  xmax=barrel_up["time"].iloc[-1], color="gray", zorder=0,
                  alpha=0.3)
    return ax
