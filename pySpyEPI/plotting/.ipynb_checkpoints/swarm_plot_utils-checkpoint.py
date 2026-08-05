"""Utility functions for plotting swarm
Updated 20 February 2026
alanahco

List of Functions
-----------------
swarm_panel
swarm_barrel
"""
import numpy as np
import math
from pySpyEPI.plotting import plotL1
from pySpyEPI.plotting. barrel_plot import barrel_ax
from pySpyEPI.plotting.detection_plotting import flag_markers


def swarm_panel(ax, data_df, y_strs, t_str='time', cols=None, ylab='y data',
                modes=None, styles=None, sizes=None, leg_ax=None,
                leg_labs=None, leg_loc='center left', mlat0=False, mlat25=True,
                x_on=True, time_bool=True, sci_y=True, x_ax_str='time',
                xlab='Universal Time', fs=14):
    """Plot a single panel of swarm data.
    Parameters
    ----------
    ax : axis handle
        axis of figure to plot
    data_df : pandas DataFrame
        whole dataset as a dataframe
    y_strs : list of stirngs
        data_df column labels to plot on y axis
    t_str : string
        data_df time column label
        default "time"
    cols : NoneType or list-like
        line or scatter colors
        list of same length as y_strs
        if None (default), '#1f77b4' is used for all
    ylab : string
        label for y axis of plot
        default 'y data'
    modes : NoneType or list-like
        strings of either "line" or "scatter" or "shading"
        if shading, y_strs must be a list within the list of 2 values
        list of same length as y_strs
        if None (default), "line" will be assumed for all y_strs
    styles : NoneType or list-like
        linestyles
        list of same length as y_strs
        if None (default), '-' for line and '.' for scatter
    sizes : NoneType or list-like
        list of floats for the linewidth (mode=line)
        or markersize (mode=scatter)
        default NoneType, 1 is used.
    leg_ax : NoneType or axis handle
        axis handle for legend
        if None (default), no legend will be included
    leg_labs : NoneType or list-like
        list of same length as y_strs
        if None (default) and leg_ax is not None, ylabs will be used
    leg_loc : string
        legend location
        default "center left"
    mlat0, mlat25 : bool
        if True, 0 and +/- 25 degree mlat lines will be added
        default mlat0, False
        default mlat25, True
    x_on : bool
        if True (default), xlabels will be plotted
    time_bool : boolean
        True for plotting the x axis as a time array
        default True
    sci_y : bool
        y axis labels in scientific notation
        default True
    x_ax_str : string
        parameter for x axis labelling
        default 'time'
    xlab : string
        x axis label
        default "Universal Time"
    fs : int
        fontsize
        default 14
    Return
    ------
    ax : axis handle
        updated figure axis
    """

    # initialize parameters ---------------------------------------------------

    # colors
    # make all blue
    if cols is None:
        cols = ['#1f77b4'] * len(y_strs)

    # modes
    # assume line  plots
    if modes is None:
        modes = ['line'] * len(y_strs)

    # styles
    # '-' for line
    # '.' for scatter
    if styles is None:
        styles = []
        for m in modes:
            if m == 'line':
                styles.append('-')
            elif m == 'scatter':
                styles.append('.')

    # sizes
    if sizes is None:
        sizes = [2] * len(y_strs)

    # set time array
    t_data = data_df[t_str]
    # -------------------------------------------------------------------------

    # Set lanbel locations ----------------------------------------------------
    x_lick = len(t_data) - 1
    xtick_loc = [math.ceil(x_lick * i / 10) for i in range(11)]

    leg_alphas = []
    # Iterate through y strings and plot --------------------------------------
    for y, col, style, mode, size in zip(y_strs, cols, styles, modes, sizes):

        # initialize y data
        y_data = data_df[y]

        if mode == "line":
            ax.plot(t_data, y_data, color=col, linestyle=style, linewidth=size)
            leg_alphas.append(1)
        elif mode == "scatter":
            ax.scatter(t_data, y_data, color=col, marker=style, s=size)
            leg_alphas.append(1)
        elif mode == "shading":
            y1 = data_df[y[0]]
            y2 = data_df[y[1]]
            if size > 1:
                alpha = 0.8
            else:
                alpha = size
            ax.fill_between(t_data, y1, y2, color=col, alpha=alpha)
            leg_alphas.append(alpha)

        else:
            raise ValueError(f"Unknown mode: {mode}")
    # -------------------------------------------------------------------------

    # adjust y axes and add lat lines -----------------------------------------
    x_data = data_df[x_ax_str].values
    plotL1.adjust_axes_time(x_data, ax, t_data.values,
                            xlab=xlab, ylab=ylab, time_bool=time_bool,
                            leg=False, fs=fs, x_on=x_on, sci_y=sci_y)

    # Add lattiude lines
    if mlat0:
        plotL1.add_lat_lines(data_df, ax, 0, lat_type='mlat')
    if mlat25:
        plotL1.add_lat_lines(data_df, ax, 25, lat_type='mlat')
    # -------------------------------------------------------------------------

    # add legend if leg_ax is given -------------------------------------------
    if leg_ax is not None:

        # specify leg_labs if none provided
        if leg_labs is None:
            leg_labs = y_strs

        plotL1.make_legend(leg_ax, leg_labs=leg_labs, leg_cols=cols,
                           leg_styles=styles, modes=modes,
                           leg_alphas=leg_alphas, loc=leg_loc, fontsize=fs)

    return ax


def swarm_barrel(ax, leg_ax, data_df, barrel_df, peaks, properties, stats_df,
                 xtick_loc, detrended=False, marker='o', mark_color="green",
                 Np_flag="ne_flag", x_on=True, flag_mark=False,
                 leg_loc='center left', xlab='Universal Time', eia_state=None,
                 eia_lats=None, barrel_weight=False, fs=14):
    """Plot swarm barrel panel.

    Parameters
    ----------
    ax : axis handle
        ax for plotting
    leg_ax : axis handle
        axis handle for legend
    data_df : dataframe
        dataframe of swarm data
    barrel_df : dataframe
        dataframe created by barrel_roll.get_barrel_info and updated by
        barre_roll.just_barrel
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    stats_df :
        dataframe containing peak info such as flag
        from cand_info
    xtick_loc : array-like
        array of tick mark locs
    detrended : kwarg bool
        default False, will plot ne with filtered ne
        True will plot detrended ne with peaks/prominences
    marker : kwarg str
        marker for peaks
        'o' is good for keep and 'x' is good for remove
    mark_color : kwarg str
        color string
        'green' is good for kept peaks
    Np_flag : string
        np flag parameter string
    x_on : bool
        if True (default), xlabels will be plotted
    flag_mark : bool
        if True, flag markers will be plotted based on stats_df
        if false (default) only red X's and green O's will be plotted
    leg_loc : string
        legend location
        default "center left"
    xlab : string
        x axis label
        default 'Univeral Time'
    eia_state : NoneType or string
        if None (default), EIA type will not be added to legend of density
    eia_lats : NoneType or list-like
        if provided, EIA creat latitudes will be plotted as vertical lines in
        density plot
        default None
    fs : int
        fontsize
        default 14
    Returns
    -------
    ax : axis with data plotted
    """

    if not detrended:
        # Ne with trend -------------------------------------------------------
        barrel_ax(barrel_df, peaks, properties, ax, detrended=False)

        if barrel_weight:
            ax.plot(barrel_df['time'], barrel_df['weight_barrel'],
                    color='#C64B8C', linestyle='dashed')

        plotL1.adjust_axes_time(data_df['time'].values, ax,
                                data_df['time'].values, xlab=xlab,
                                ylab=r'$N_{p}$ ($cm^{-3}$)', time_bool=True,
                                leg=False, fs=fs, sci_y=True, x_on=x_on)

        # add lines for EIA peaks
        if eia_lats is not None:
            for pp in eia_lats:
                if pp:
                    plotL1.add_lat_lines(data_df, ax, pp, lat_type='mlat',
                                         col='green', dir='N', style=':')
        # add latitude lines
        plotL1.add_lat_lines(data_df, ax, 0, lat_type='mlat')
        plotL1.add_lat_lines(data_df, ax, 25, lat_type='mlat')

        # Add legend on side axis

        if np.any(data_df["ne_flag"] > 2):
            leg_labs = [r'$N_{p}$', r'Barrel $N_{p}$', 'flagged']
            leg_cols = ['#1f77b4', 'orange', 'red']
            modes = ['line', 'line', 'scatter']
            leg_styles = ['-', '-', 'o']
        else:
            leg_labs = [r'$N_{p}$', r'Barrel $N_{p}$']
            leg_cols = ['#1f77b4', 'orange']
            modes = ['line', 'line']
            leg_styles = ['-', '-']

        if barrel_weight:
            leg_labs.append('Weighted Barrel')
            leg_cols.append('#C64B8C')
            modes.append('line')
            leg_styles.append('--')

        if eia_state is not None:
            leg_labs.append(eia_state)
            leg_cols.append('green')
            modes.append('line')
            leg_styles.append(':')

        plotL1.make_legend(leg_ax, leg_labs=leg_labs, leg_cols=leg_cols,
                           leg_styles=leg_styles, modes=modes, loc=leg_loc,
                           fontsize=fs)
    if detrended:
        # Panel 4 Detrended Ne ------------------------------------------------

        # Plot kept peaks
        mask_keep = (stats_df['epi_flag'].values == 1)
        peaks_keep = peaks[mask_keep]
        prop_keep = (
            {key: value[mask_keep] for key, value in properties.items()
             if isinstance(value, np.ndarray)})

        if not flag_mark:

            barrel_ax(barrel_df, peaks_keep, prop_keep, ax, detrended=True,
                      marker='o', mark_color="green")

            # Plot Discarded peaks
            mask_rem = (stats_df['epi_flag'].values != 1)
            peaks_rem = peaks[mask_rem]
            prop_rem = (
                {key: value[mask_rem] for key, value in properties.items()
                 if isinstance(value, np.ndarray)})
            barrel_ax(barrel_df, peaks_rem, prop_rem, ax, detrended=True,
                      marker='x', mark_color="red")

            plotL1.adjust_axes_time(data_df['time'].values, ax,
                                    data_df['time'].values, xlab=xlab,
                                    ylab=r'Detrended $N_{p}$ ($cm^{-3}$)',
                                    time_bool=True, leg=False,
                                    fs=fs, sci_y=True, x_on=x_on)

            # add latitude lines
            plotL1.add_lat_lines(data_df, ax, 0, lat_type='mlat')
            plotL1.add_lat_lines(data_df, ax, 25, lat_type='mlat')

            # Add legend on side axis
            leg_labs = [r'Barrel - $N_{p}$', 'EPI', 'Discard']
            leg_cols = ['#1f77b4', 'green', 'red']
            modes = ['line', 'scatter', 'scatter']
            leg_styles = ['-', 'o', 'x']

            plotL1.make_legend(leg_ax, leg_labs=leg_labs, leg_cols=leg_cols,
                               leg_styles=leg_styles, modes=modes, loc=leg_loc,
                               fontsize=fs)
        else:
            f_marks, f_labs, f_cols, f_mode, f_flag = flag_markers(
                ax, peaks, properties, stats_df,
                barrel_df)

            plotL1.adjust_axes_time(data_df['time'].values, ax,
                                    data_df['time'].values, xlab=xlab,
                                    ylab=r'Detrended $N_p$ ($cm^{-3}$)',
                                    time_bool=True, leg=False, fs=fs,
                                    sci_y=True, x_on=x_on)

            # add latitude lines
            plotL1.add_lat_lines(data_df, ax, 0, lat_type='mlat')
            plotL1.add_lat_lines(data_df, ax, 25, lat_type='mlat')

            # Add legend on side axis
            leg_labs = [r'Barrel - $N_{p}$']
            leg_cols = ['#1f77b4']
            modes = ['line']
            leg_styles = ['-']

            # add flag labels
            if len(f_marks) > 0:
                leg_labs += f_labs.tolist()
                leg_cols += f_cols.tolist()
                leg_styles += f_marks.tolist()
                modes += f_mode.tolist()

            plotL1.make_legend(leg_ax, leg_labs=leg_labs, leg_cols=leg_cols,
                               leg_styles=leg_styles, modes=modes, loc=leg_loc,
                               fontsize=fs)
    return ax
