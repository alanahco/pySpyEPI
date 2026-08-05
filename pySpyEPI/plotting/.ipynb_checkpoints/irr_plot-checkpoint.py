"""Plotting functions for IRR data
Updated : 23 February 2026
alanahco

List Functions
--------------
ne_ibarrel
ne_2barrel
fft_ax
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
from pySpyEPI.plotting import plotL1
from pySpyEPI.plotting import swarm_plot_utils as spu
import matplotlib.ticker as mticker
from matplotlib.colors import Normalize


def ne_ibarrel(IRR_data, satellite, barrel_df, peaks,
               properties, stats_df, EFI_data=None, dne_bool=True,
               hspace=0.4, fig_size=None, leg_loc='center left', fs=14,
               flag_mark=False, eia_state=None, eia_lats=None, y_sup=0.96,
               x_sup=0.5):
    """Plot coordinates, Ne, deltaNe 40s, and barrel with barrel info as input.

    Parameters
    ----------
    IRR_data : dataframe
        dataframe of IRR data
    satellite : str
        Swarm A, B, or C for Plotting Purposes
    barrel_df : DataFrame
        containing barrel info from get_barrel_info plus added columns:
        'time', 'ne_og', and 'Ne_flag'
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    stats_df : DataFrame
        df containing info about the peaks and if they are bubbles
        from cand_info
    EFI_data : dataframe kwarg empty
        EFI data fro Te if te_bool is True
    te_bool : bool kwarg False
        will plot EFI Te if True
    dne_bool : bool kwarg True
        will plot delta ne if True
    hspace : float
        height space between each panel
        default is 0.4
    fig_size : NoneType or list-like
        if None, figsize wil be 18,15 for te_bool and 18,13 for not te_bool
        list of 2 numbers width, height
    leg_loc : string
        location of legend default (center)
    fs : kwarg int
        font size default 14
    flag_mark : bool
        if True, flag markers will be plotted based on stats_df
        if false (default) only red X's and green O's will be plotted
    eia_state : NoneType or string
        if None (default), EIA type will not be added to legend of density
    eia_lats : NoneType or list-like
        if provided, EIA creat latitudes will be plotted as vertical lines in
        density plot
        default None
    y_sup, x_sup : float
        y and x placement of super title
        default 0.96
    Returns
    -------
    Fig : figure handle
        2-4 panel figure
        Coordinate axis,
        Ne, and detrended Ne
    Notes
    -----
    """
    # Initialize Figure size --------------------------------------------------

    if EFI_data is not None:
        t_add = 3
    else:
        t_add = 0

    if dne_bool:
        d_add = 3
    else:
        d_add = 0

    h_fig = 9 + d_add + t_add
    h_pans = int(h_fig / 3)

    if fig_size is None:
        fig = plt.figure(figsize=(18, h_fig))
    else:
        fig = plt.figure(figsize=(fig_size[0], fig_size[1]))

    gs = gridspec.GridSpec(nrows=h_pans, ncols=4, width_ratios=[1, 1, 1, 0.05],
                           hspace=hspace)

    plt.rcParams.update({'font.size': fs})

    # set up ax_list for shading purposes
    ax_list = []

    # Panel 1: Top plot (with colorbar)----------------------------------------
    ax0 = fig.add_subplot(gs[0, 0:3])
    cax = fig.add_subplot(gs[0, 3])
    plotL1.plot_maglat_panel(IRR_data, ax=ax0, cax=cax, fs=fs)

    # Initialize Tick locations based off of consistent time intervals
    xlab = ''
    next_pan = 1
    last_pan = h_pans - 1
    if dne_bool:
        # Panel 2: Delta Ne----------------------------------------------------
        ax1 = fig.add_subplot(gs[next_pan, 0:3])
        leg_ax = fig.add_subplot(gs[next_pan, 3])

        if next_pan % 2 == 0:
            x_on = True
        else:
            x_on = False

        if next_pan == last_pan:
            xlab = 'Universal Time'
            x_on = True
        else:
            xlab = ''

        y_strs = ['delta_ne40s']
        cols = ['#1f77b4']
        ylab = r'$\Delta$$N_{p}$/$N_{p}$'
        modes = ['line']
        styles = ['-']
        leg_labs = ['40s']

        spu.swarm_panel(ax1, IRR_data, y_strs, t_str='time', cols=cols,
                        ylab=ylab, modes=modes, styles=styles, leg_ax=leg_ax,
                        leg_labs=leg_labs, leg_loc='center left', mlat0=False,
                        mlat25=True, x_on=x_on, time_bool=True, sci_y=True,
                        x_ax_str='time', xlab=xlab, fs=fs)

        next_pan += 1
        ax_list.append(ax1)

    # Panel 3: IRI Ne----------------------------------------------------------
    ax2 = fig.add_subplot(gs[next_pan, 0:3])
    leg_ax = fig.add_subplot(gs[next_pan, 3])

    x_lick = len(IRR_data['time']) - 1
    xtick_loc = [math.ceil(x_lick * i / 10) for i in range(11)]

    if next_pan % 2 == 0:
        x_on = True
    else:
        x_on = False

    if next_pan == last_pan:
        xlab = 'Universal Time'
        x_on = True
    else:
        xlab = ''

    spu.swarm_barrel(ax2, leg_ax, IRR_data, barrel_df, peaks, properties,
                     stats_df, xtick_loc, detrended=False, marker='o',
                     mark_color="green", Np_flag="ne_flag", x_on=x_on,
                     flag_mark=False, xlab=xlab, eia_state=eia_state,
                     eia_lats=eia_lats, fs=fs)
    next_pan += 1
    ax_list.append(ax2)

    # Panel 4 Detrended Ne ----------------------------------------------------
    ax3 = fig.add_subplot(gs[next_pan, 0:3])
    leg_ax = fig.add_subplot(gs[next_pan, 3])

    if next_pan % 2 == 0:
        x_on = True
    else:
        x_on = False

    if next_pan == last_pan:
        xlab = 'Universal Time'
        x_on = True
    else:
        if EFI_data is not None:
            if len(EFI_data) == 0:
                xlab = 'Universal Time'
                x_on = True
            else:
                xlab = ''
        else:
            xlab = ''

    spu.swarm_barrel(ax3, leg_ax, IRR_data, barrel_df, peaks, properties,
                     stats_df, xtick_loc, detrended=True, marker='o',
                     mark_color="green", Np_flag="ne_flag", x_on=x_on,
                     flag_mark=flag_mark, xlab=xlab, fs=fs)

    next_pan += 1
    ax_list.append(ax3)

    # Panel 5 Te --------------------------------------------------------------
    if EFI_data is not None:
        if len(EFI_data) > 0:
            ax4 = fig.add_subplot(gs[next_pan, 0:3])
            leg_ax = fig.add_subplot(gs[next_pan, 3])
            if next_pan % 2 == 0:
                x_on = True
            else:
                x_on = False

            if next_pan == last_pan:
                xlab = 'Universal Time'
                x_on = True
            else:
                xlab = ''

            y_strs = ['te_adj']
            cols = ['#1f77b4']
            ylab = r'$T_{e}$ (K)'
            modes = ['scatter']
            styles = ['.']
            sizes = [8]
            leg_labs = [r'10s Filt $T_e$']

            spu.swarm_panel(ax4, EFI_data, y_strs, t_str='time', cols=cols,
                            ylab=ylab, modes=modes, styles=styles, sizes=sizes,
                            leg_ax=leg_ax, leg_labs=leg_labs, mlat0=False,
                            mlat25=True, x_on=x_on, time_bool=True, sci_y=True,
                            x_ax_str='time', xlab=xlab, fs=fs)
            ax_list.append(ax4)

    # add shading to figure
    mask_keep = (stats_df['epi_flag'].values == 1)
    prop_keep = (
        {key: value[mask_keep] for key, value in properties.items()
         if isinstance(value, np.ndarray)})
    plotL1.add_peak_shading(ax_list, IRR_data['time'].values, prop_keep,
                            col='gray', alpha=0.1)

    # Super Title -------------------------------------------------------------
    title_date = IRR_data['time'].iloc[0].strftime('%d %b %Y')
    titl = title_date + ' Swarm ' + satellite
    plt.suptitle(titl, y=y_sup, x=x_sup, fontsize=fs)

    return fig


def ne_2barrel(IRR_data, satellite, barrel_df, peaks,
               properties, stats_df, EFI_data=None, blo_bool=True,
               hspace=0.4, fig_size=None, leg_loc='center left', fs=14,
               flag_mark=False, eia_state=None, eia_lats=None, y_sup=0.96,
               x_sup=0.5):
    """Plot coordinates, Ne, deltaNe 40s, and barrel with barrel info as input.

    Parameters
    ----------
    IRR_data : dataframe
        dataframe of IRR data
    satellite : str
        Swarm A, B, or C for Plotting Purposes
    barrel_df : DataFrame
        containing barrel info from get_barrel_info plus added columns:
        'time', 'ne_og', and 'Ne_flag'
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    stats_df : DataFrame
        df containing info about the peaks and if they are bubbles
        from cand_info
    EFI_data : dataframe kwarg empty
        EFI data fro Te if te_bool is True
    te_bool : bool kwarg False
        will plot EFI Te if True
    dne_bool : bool kwarg True
        will plot delta ne if True
    hspace : float
        height space between each panel
        default is 0.4
    fig_size : NoneType or list-like
        if None, figsize wil be 18,15 for te_bool and 18,13 for not te_bool
        list of 2 numbers width, height
    leg_loc : string
        location of legend default (center)
    fs : kwarg int
        font size default 14
    flag_mark : bool
        if True, flag markers will be plotted based on stats_df
        if false (default) only red X's and green O's will be plotted
    eia_state : NoneType or string
        if None (default), EIA type will not be added to legend of density
    eia_lats : NoneType or list-like
        if provided, EIA creat latitudes will be plotted as vertical lines in
        density plot
        default None
    y_sup, x_sup : float
        y and x placement of super title
        default 0.96
    Returns
    -------
    Fig : figure handle
        2-4 panel figure
        Coordinate axis,
        Ne, and detrended Ne
    Notes
    -----
    """
    # Initialize Figure size --------------------------------------------------

    if EFI_data is not None:
        t_add = 3
    else:
        t_add = 0

    if blo_bool:
        d_add = 3
    else:
        d_add = 0

    h_fig = 9 + d_add + t_add
    h_pans = int(h_fig / 3)

    if fig_size is None:
        fig = plt.figure(figsize=(18, h_fig))
    else:
        fig = plt.figure(figsize=(fig_size[0], fig_size[1]))

    gs = gridspec.GridSpec(nrows=h_pans, ncols=4, width_ratios=[1, 1, 1, 0.05],
                           hspace=hspace)

    plt.rcParams.update({'font.size': fs})

    # set up ax_list for shading purposes
    ax_list = []

    # Panel 1: Top plot (with colorbar)----------------------------------------
    ax0 = fig.add_subplot(gs[0, 0:3])
    cax = fig.add_subplot(gs[0, 3])
    plotL1.plot_maglat_panel(IRR_data, ax=ax0, cax=cax, fs=fs)

    # Initialize Tick locations based off of consistent time intervals
    xlab = ''
    next_pan = 1
    last_pan = h_pans - 1
    if blo_bool:
        # Panel 2: Delta Ne----------------------------------------------------
        ax1 = fig.add_subplot(gs[next_pan, 0:3])
        leg_ax = fig.add_subplot(gs[next_pan, 3])

        if next_pan % 2 == 0:
            x_on = True
        else:
            x_on = False

        if next_pan == last_pan:
            xlab = 'Universal Time'
            x_on = True
        else:
            xlab = ''

        IRR_copy = IRR_data.copy()
        IRR_copy['dne_wdr'] = barrel_df['dne_wdr'].values
        y_strs = ['dne_wdr']
        cols = ['#C64B8C']
        ylab = r'$N_{p}$/B$N_p$'
        modes = ['line']
        styles = ['-']
        leg_labs = ['Weighted Density Ratio']

        spu.swarm_panel(ax1, IRR_copy, y_strs, t_str='time', cols=cols,
                        ylab=ylab, modes=modes, styles=styles, leg_ax=leg_ax,
                        leg_labs=leg_labs, leg_loc='center left', mlat0=False,
                        mlat25=True, x_on=x_on, time_bool=True, sci_y=True,
                        x_ax_str='time', xlab=xlab, fs=fs)

        next_pan += 1
        ax_list.append(ax1)

    # Panel 3: IRI Ne----------------------------------------------------------
    ax2 = fig.add_subplot(gs[next_pan, 0:3])
    leg_ax = fig.add_subplot(gs[next_pan, 3])

    x_lick = len(IRR_data['time']) - 1
    xtick_loc = [math.ceil(x_lick * i / 10) for i in range(11)]

    if next_pan % 2 == 0:
        x_on = True
    else:
        x_on = False

    if next_pan == last_pan:
        xlab = 'Universal Time'
        x_on = True
    else:
        xlab = ''

    spu.swarm_barrel(ax2, leg_ax, IRR_data, barrel_df, peaks, properties,
                     stats_df, xtick_loc, detrended=False, marker='o',
                     mark_color="green", Np_flag="ne_flag", x_on=x_on,
                     flag_mark=False, xlab=xlab, eia_state=eia_state,
                     eia_lats=eia_lats, barrel_weight=True, fs=fs)

    next_pan += 1
    ax_list.append(ax2)

    # Panel 4 Detrended Ne ----------------------------------------------------
    ax3 = fig.add_subplot(gs[next_pan, 0:3])
    leg_ax = fig.add_subplot(gs[next_pan, 3])

    if next_pan % 2 == 0:
        x_on = True
    else:
        x_on = False

    if next_pan == last_pan:
        xlab = 'Universal Time'
        x_on = True
    else:
        if EFI_data is not None:
            if len(EFI_data) == 0:
                xlab = 'Universal Time'
                x_on = True
            else:
                xlab = ''
        else:
            xlab = ''

    spu.swarm_barrel(ax3, leg_ax, IRR_data, barrel_df, peaks, properties,
                     stats_df, xtick_loc, detrended=True, marker='o',
                     mark_color="green", Np_flag="ne_flag", x_on=x_on,
                     flag_mark=flag_mark, xlab=xlab, fs=fs)

    next_pan += 1
    ax_list.append(ax3)

    # Panel 5 Te --------------------------------------------------------------
    if EFI_data is not None:
        if len(EFI_data) > 0:
            ax4 = fig.add_subplot(gs[next_pan, 0:3])
            leg_ax = fig.add_subplot(gs[next_pan, 3])
            if next_pan % 2 == 0:
                x_on = True
            else:
                x_on = False

            if next_pan == last_pan:
                xlab = 'Universal Time'
                x_on = True
            else:
                xlab = ''

            y_strs = ['te_adj']
            cols = ['#1f77b4']
            ylab = r'$T_{e}$ (K)'
            modes = ['scatter']
            styles = ['.']
            sizes = [8]
            leg_labs = [r'10s Filt $T_e$']

            spu.swarm_panel(ax4, EFI_data, y_strs, t_str='time', cols=cols,
                            ylab=ylab, modes=modes, styles=styles, sizes=sizes,
                            leg_ax=leg_ax, leg_labs=leg_labs, mlat0=False,
                            mlat25=True, x_on=x_on, time_bool=True, sci_y=True,
                            x_ax_str='time', xlab=xlab, fs=fs)
            ax_list.append(ax4)

    # add shading to figure
    mask_keep = (stats_df['epi_flag'].values == 1)
    prop_keep = (
        {key: value[mask_keep] for key, value in properties.items()
         if isinstance(value, np.ndarray)})
    plotL1.add_peak_shading(ax_list, IRR_data['time'].values, prop_keep,
                            col='gray', alpha=0.1)

    # Super Title -------------------------------------------------------------
    title_date = IRR_data['time'].iloc[0].strftime('%d %b %Y')
    titl = title_date + ' Swarm ' + satellite
    plt.suptitle(titl, y=y_sup, x=x_sup, fontsize=fs)

    return fig


def fft_ax(spec, time_axis, freq_axis, ax, cax, t1, fs=14,
           amp_lab=r'($cm^{-3}$)'):
    """
    Plot the dynamic spectrum as a heatmap in a single panel
    Parameters
    ----------
    spec : 2D array
        FFT magnitudes (time x frequency) from compute_dynamic_spectra
    time_axis: array-like
        center times (can be in seconds)
    freq_axis : array-like
        frequency bins (Hz)
    ax : axis
        for heatmap
    cax : axis
        for colorbar
    t1 : datetime
        starting time, for now, not used, but may be if I want to change from
        distance into datetime for the x axis
    fs : kwarg int
        fontsize
    Returns
    -------
    ax : axis handle
        axis with data heatmap
    """

    # Plot spectra
    # Transpose spec so that freq is y-axis and time is x-axis
    # Plot pcolormesh
    ax.pcolormesh(time_axis, freq_axis, spec.T, shading='auto', cmap='cividis')

    # Add labels
    ax.set_xlabel("Distance Traveled (km)", fontsize=fs)
    ax.set_ylabel("Frequency (Hz)", fontsize=fs)
    ax.grid(True)
    ax.set_xlim([0, time_axis[-1]])

    ax.xaxis.set_major_locator(mticker.LinearLocator(11))
    ax.xaxis.set_minor_locator(mticker.LinearLocator(31))
    # set tick lengths
    ax.tick_params(axis='both', which='major', labelsize=fs, width=1.5,
                   length=10)
    ax.tick_params(axis='both', which='minor', labelsize=fs, width=1.5,
                   length=5)

    # Change tick labels to distance
    # Get the ticks and labels
    xticks = ax.get_xticks()
    ax.set_xticks(xticks)
    x_labels = [label.get_text() for label in ax.get_xticklabels()]

    # Convert to distance
    x_ints = [int(s) for s in x_labels]
    x_ints = np.array(x_ints)
    x_dist = x_ints * 7.5  # Swarm travels ~7.5 km/s
    ax.set_xticklabels(x_dist)

    # Add color bar
    sm = plt.cm.ScalarMappable(
        cmap='cividis', norm=Normalize(vmin=spec.T.min(), vmax=spec.T.max()))
    cb = ax.figure.colorbar(sm, cax=cax)
    cb.set_label("Amplitude " + amp_lab, fontsize=fs, rotation=270,
                 labelpad=20)

    return ax
