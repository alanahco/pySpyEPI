
"""GOLD plotting Functions
alanahco
Created 14 April 2026

List Functions
--------------
gold_ax
zdrift_plot
gold_swarm_single
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import Normalize
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import numpy as np
from pySpyEPI.plotting import plotL1
from pySpyEPI.utils.calc import moving_median


def gold_ax(axg, Gdc1, heat_param, vmin, vmax, cax=None,
            mlat_lim=None, lon_lim=None, cax_bool=True, cmap='gray',
            pan='both', side='left', t_color='k', t_loc='inside',
            lon_type='lon', lat_type='mlat', facecolor=None,
            tformat='%d %b %H:%M UT', fs=12, cax_loc='right'):
    """ Plot GOLD axis panel

    Parameters
    ----------
    axg : axis handle
        axis for first GOLD map
    Gdc1 : dictionary
        GOLD dictionary for first pannel made by open_GOLD_barrel
    heat_param : str
        which parameter to plot
        either 'oi' or 'deplete oi'
    vmin : float
        minimum value for heatmap
    vmax : float
        maximum value for heatmap
    cax : axis handle or NoneType
        axis for color map
        default (None) nocolor bar axis plotted
    mlat_lim : NoneType or list of 2
        magnetic latitude limit for y axis
        default None, [-40, 40] will be used
    lon_lim : list of 2 values kwarg
        longitude limit for x axis
        default None, assumes GOLD span from Gdc1
    cax_bool: bool kwarg True
        plot color bar if true
    cmap : str kwarg
        colormap colors
        default is gray
    pan : str kwarg default 'both'
        pannel type to determine where to place axes information
        options: 'top', 'bottom', 'middle', 'both'
    side : str kwarg default left
        determines what side info will be plotted on left or right
    t_color : str kwarg default black
        color for text
    t_loc : string
        location for time string, 'inside' default
        if 'top', it will be the axis title
    lon_type : str
        determines if geographic or mangeitc lon will be plotted
        default 'lon' for geographic lon
        'mlon' for magnetic
    lat_type : str
        determines if geographic or mangeitc lat will be plotted
        'lat' for geographic lon
        default 'mlat' for magnetic
    facecolor : None Type or string
        if not none, plot facecolor will be set as color provided
    tformat : string
        format for strftime
        default '%d %b %H:%M UT'
    fs : int
        fontsize default 12
    Returns
    -------
    plotted panels with GOLD and swarm trajectories if desired

    """

    # if lon_lim and mlat_lim are not provided, use
    if lon_lim is None:
        lon_grid = Gdc1[lon_type]
        lon_min = lon_grid.min()
        lon_max = lon_grid.max()
        lon_lim = [lon_min, lon_max]
    else:
        lon_min = min(lon_lim)
        lon_max = max(lon_lim)

    if mlat_lim is None:
        mlat_lim = [-40, 40]

    # set labels based on lon and lat types
    if lon_type == 'lon':
        lon_lab = 'Geographic Longitude'
    else:
        lon_lab = 'Magnetic Longitude'

    if lat_type == 'lat':
        lat_lab = 'Geographic Latitude'
    else:
        lat_lab = 'Magnetic Latitude'

    # Panel Left GOLD ---------------------------------------------------------
    axg.set_global()
    axg.pcolormesh(Gdc1[lon_type], Gdc1[lat_type], Gdc1[heat_param], vmin=vmin,
                   vmax=vmax, cmap=cmap, transform=ccrs.PlateCarree())
    axg.set_xlim(lon_lim)
    axg.set_ylim(mlat_lim)

    # add grid lines
    gl = axg.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                       alpha=0.5)

    gl.xlocator = mticker.MultipleLocator(5)
    gl.ylocator = mticker.MultipleLocator(10)

    gl.right_labels = False  # Optional: Turn off right labels
    axg.set_aspect('auto', adjustable='box')

    # Add labels based on limits provided
    if side == 'left':
        axg.text(-0.2, 0.5, lat_lab, transform=axg.transAxes, ha='right',
                 va='center', rotation=90, fontsize=fs)
    else:
        gl.left_labels = False

    # turn off bottom labels if this is a top panel or middle panel
    if (pan == 'bottom') | (pan == 'both'):
        axg.text(0.5, -0.06, lon_lab, transform=axg.transAxes, ha='center',
                 va='top', fontsize=fs)
        gl.top_labels = False

    # if it is a middle panel turn offf top and bottom labels
    elif (pan == 'middle'):
        gl.top_labels = False
        gl.bottom_labels = False
    else:
        gl.top_labels = True
        gl.bottom_labels = False
        axg.text(0.5, -0.1, lon_lab, transform=axg.transAxes, ha='center',
                 va='top', fontsize=fs)

    # Add time to plot
    t_str = Gdc1['timeA'].strftime(tformat)

    if t_loc != 'top':
        axg.text(0.5, 0.95, t_str, transform=axg.transAxes, ha='center',
                 va='top', color=t_color, fontsize=fs, fontweight='bold')
    else:
        axg.set_title(t_str, fontsize=fs)

    if facecolor is not None:
        axg.set_facecolor(facecolor)

    # Color Bar panel ---------------------------------------------------------
    # create a color bar
    if cax_bool:
        sm = plt.cm.ScalarMappable(cmap=cmap,
                                   norm=Normalize(vmin=vmin, vmax=vmax))

        if cax_loc == 'right':
            cb = axg.figure.colorbar(sm, cax=cax)

            # set heatmap label
            if heat_param == 'oi':
                cb.set_label("OI 135.6 nm Brightness (R)", rotation=270,
                             fontsize=fs, labelpad=30)
            else:
                cb.set_label("Radiance: Net Depletion (R)", rotation=270,
                             fontsize=fs, labelpad=30)
        elif cax_loc == 'left':
            cb = axg.figure.colorbar(sm, cax=cax)

            # set heatmap label
            if heat_param == 'oi':
                cb.set_label("OI 135.6 nm Brightness (R)",
                             fontsize=fs, labelpad=30)
            else:
                cb.set_label("Radiance: Net Depletion (R)",
                             fontsize=fs, labelpad=30)
        elif cax_loc == 'top':
            cb = axg.figure.colorbar(sm, cax=cax, orientation='horizontal')

            # set heatmap label
            if heat_param == 'oi':
                cb.set_label("OI 135.6 nm Brightness [R]",
                             fontsize=fs, labelpad=10)
            else:
                cb.set_label("Radiance: Net Depletion [R]",
                             fontsize=fs, labelpad=10)

            cb.ax.xaxis.set_label_position('top')
            cb.ax.xaxis.set_ticks_position('top')
        elif cax_loc == 'bottom':
            cb = axg.figure.colorbar(sm, cax=cax, orientation='horizontal')

            # set heatmap label
            if heat_param == 'oi':
                cb.set_label("OI 135.6 nm Brightness (R)",
                             fontsize=fs, labelpad=10)
            else:
                cb.set_label("Radiance: Net Depletion (R)",
                             fontsize=fs, labelpad=10)


def zdrift_plot(all_gold, lon_type, oi_type, drift_df, gcols, xlon, fs=12,
                o_smooth=2, plot2=False):
    """Plot Zonal Drifts for different times like in Karan et al 2020 paper.

    Parameters:
    all_gold : list-like
        list of GOLD dictionaries
    lon_type : str
        determines if geographic or mangeitc lon will be plotted
        default 'mlon' for magnetic lon
        'lon' for geographic
     oi_type : str
        which parameter to plot
        either 'oi' or 'deplete oi'
    drift_df : DataFrame
        dataframe of drift and bubble info made by drift_calc_GOLD.zdrift_df
    gcols : list-like
        list of colors for super groups plotting
    xlon : array-like
        longitude limits
    fs : float
        fontsize default 12
    Returns
    -------
    fig : figure handle
        3 columns, first column: oi vs lon, second: lon vs time, third: legend
    Notes
    -----
    Warning: mlon or oi_deplete may break this function!!!
    """

    if plot2 is False:
        # remove doubles before plotting
        sing_gr = []

        for gr in drift_df['group'].unique():
            gr_len = len(drift_df['group'][drift_df['group'] == gr])
            if gr_len == 2:
                sing_gr.append(gr)

        for si in sing_gr:
            drift_df = drift_df.copy()
            drift_df = drift_df[drift_df['group'] != si]

    # set mlat_bins based on drift_df
    mlat_bins = np.linspace(drift_df['mlat1'].min(), drift_df['mlat2'].max(),
                            len(drift_df['mlat2'].unique()) + 1)

    # create figure, set panels based on lat bins
    pan_nums = len(mlat_bins) - 1
    fig = plt.figure(figsize=(16, 20))
    gs = gridspec.GridSpec(nrows=pan_nums, ncols=3, width_ratios=[1, 1, 0.1],
                           wspace=0.3, hspace=0.3)
    plt.rcParams.update({'font.size': fs})

    # create lists of figure handles
    ax_list = []
    ax2_list = []
    for pn in range(pan_nums):

        # oi lon plots
        ax = fig.add_subplot(gs[pn, 0])
        ax_list.append(ax)

        # lon time plots
        ax2 = fig.add_subplot(gs[pn, 1])
        ax2_list.append(ax2)

    lsty = ['-', '-.', '--', ':']
    scat_mark = ['^', 'o', 'P', 's']
    t_labs = []
    Gdc1 = all_gold[0]
    t_start = Gdc1['timeA'].hour * 60 + Gdc1['timeA'].minute

    for g in range(len(all_gold)):
        # get GOLD data and time label
        Gdc = all_gold[g]
        t_lab = Gdc['timeA'].strftime('%H:%M UT')
        t_labs.append(t_lab)
        t_g = Gdc['timeA'].hour * 60 + Gdc['timeA'].minute - t_start

        # get basic grids
        oi_grid = Gdc[oi_type]
        mlat_grid = Gdc['mlat']
        lon_grid = Gdc['lon']
        mlon_grid = Gdc['mlon']

        for i in range(len(mlat_bins) - 1):
            adj_ax = len(mlat_bins) - 2
            ax = ax_list[adj_ax - i]
            ax2 = ax2_list[adj_ax - i]
            m1 = mlat_bins[i]
            m2 = mlat_bins[i + 1]

            # get bubble info for scattering
            drift_copy = drift_df.copy()
            drift_copy = drift_copy[drift_copy['time'] == t_g]
            drift_copy = drift_copy[drift_copy['mlat1'] == m1]

            # mask magnetic latitude by mlat bins
            mlat_mask = (mlat_grid < m2) & (mlat_grid > m1)

            # use latitude mask to sum the oxygen emissions
            o_bin = []
            mlon_bin = []

            for j in range(100):
                o_bin.append(np.nansum(oi_grid[mlat_mask[:, j], j]))

                # get median magnetic longitudes
                mlon_bin.append(np.median(mlon_grid[mlat_mask[:, j], j]))

            # convert lists to arrays
            o_bin = np.array(o_bin)
            mlon_bin = np.array(mlon_bin)

            # smooth summed emissions
            o_med = moving_median(o_bin, o_smooth)

            # plot emissions and remove each end, which = 0
            o_plot = o_med[o_med != 0]
            lon_plot = lon_grid[0, o_med != 0]
            mlon_plot = mlon_bin[o_med != 0]

            if lon_type == 'mlon':
                ax.plot(mlon_plot, o_plot, linestyle=lsty[g], color='gray',
                        zorder=0)
                lon_lab = 'Magnetic Longitude'

            else:
                ax.plot(lon_plot, o_plot, linestyle=lsty[g], color='gray',
                        zorder=0)
                for si, sg in enumerate(drift_copy['super_group']):
                    col_use = gcols[sg]
                    dc_lon = drift_copy['lon'][drift_copy['super_group'] == sg]
                    dc_oi = drift_copy[oi_type][drift_copy['super_group'] == sg]
                    ax.scatter(dc_lon, dc_oi,
                               color=col_use, marker=scat_mark[g])

                    for lp in dc_lon:
                        ax2.scatter(t_g, lp, color=col_use,
                                    marker=scat_mark[g])
                lon_lab = 'Geographic Longitude'

            ax.set_ylabel("Summed Brightness")

            # set title of each MLat Sector
            if g == 0:
                deg = r'$^\circ$'

                # set direction of Mlat
                if m1 < 0:
                    dir1 = 'S'
                elif m1 > 0:
                    dir1 = 'N'
                else:
                    dir1 = ''

                if m2 < 0:
                    dir2 = 'S'
                elif m2 > 0:
                    dir2 = 'N'
                else:
                    dir2 = ''

                la1 = f"{int(abs(m1))}{deg}{dir1}"
                la2 = f"{int(abs(m2))}{deg}{dir2}"

                if (dir1 == 'S') | (dir2 == 'S'):
                    ax.set_title(f"{la2} to {la1} MLat", fontsize=fs)
                else:
                    ax.set_title(f"{la1} to {la2} MLat", fontsize=fs)

            # add x axis grid
            # remove xlabels for all but bottom row
            if i != 0:
                ax.tick_params(axis='x', labelbottom=False)
            else:
                plotL1.format_longitude_labels(ax, xy='x')
                ax.set_xlabel(lon_lab)

            # add x limits and grid
            ax.set_xlim(xlon)
            if g == 3:
                ax.xaxis.grid()

            # right plots params
            ax2.set_xlabel('Time (minutes)')

            # add x axis grid
            # remove xlabels for all but bottom row
            plotL1.format_longitude_labels(ax2, xy='y')
            ax2.set_ylabel(lon_lab)

            # add x limits and grid
            ax2.set_ylim(xlon)
            if g == 3:
                ax2.yaxis.grid()

    # Add legend on side axis
    leg_ax = fig.add_subplot(gs[0, 2])

    leg_cols = ['k', 'k', 'k', 'k']
    modes = ['line+scatter', 'line+scatter', 'line+scatter', 'line+scatter']

    plotL1.make_legend(leg_ax, leg_labs=t_labs, leg_cols=leg_cols,
                       leg_styles=lsty, modes=modes, ls_marker=scat_mark,
                       loc="center left", fontsize=fs)

    # get x axis loc for text
    max_t = drift_df['time'].max() + 3
    min_t = drift_df['time'].min() - 3

    # second axis
    # for each latitude sector
    for i in range(len(mlat_bins) - 1):
        m1 = mlat_bins[i]
        m2 = mlat_bins[i + 1]
        ax2 = ax2_list[adj_ax - i]

        bub_lat = drift_df.copy()
        bub_lat = bub_lat[bub_lat["mlat1"] == m1]

        rg = bub_lat["group"].values
        for ir in np.unique(rg):
            gcol = bub_lat["super_group"].iloc[rg == ir]

            if min(gcol) != max(gcol):
                print('WARNING')

            gcol = np.mean(gcol)
            gcol = gcols[int(gcol)]
            ax2.plot(bub_lat["time"].iloc[rg == ir],
                     bub_lat["lon"].iloc[rg == ir], color=gcol, linewidth=1)
            ax2.set_xlim([min_t, max_t - 1])
            # Returns [slope, intercept]
            if len(bub_lat["time"].iloc[rg == ir]) > 2:
                coefficients = np.polyfit(bub_lat["time"].iloc[rg == ir],
                                          bub_lat["lon"].iloc[rg == ir], 1)
                slope, intercept = coefficients

                # convert slope to degrees/hour from degrees/minute
                slope_dhr = slope * 60
                y_text = bub_lat["lon"].iloc[rg == ir]
                y_text = y_text.iloc[-1]
                f_text = f"{np.round(slope_dhr, 2)} {deg}/h"
                ax2.text(max_t, y_text, f_text)

    return fig


def gold_swarm_single(gdc, irr, prop_keep, peaks_keep, satellite='SAT',
                      heat_param='oi', vmin=0, vmax=200, mlat_lim=None,
                      lon_lim=None, cmap='civids', facecolor=None,
                      tformat='%d %b %H:%M UT', swarm_color="#F48FB1",
                      add_trajectory=False, fs=12):
    """Create a figure with the Swarm Np and detections and GOLD Heat Map.

    Parameters
    ----------
    gdc : dictionary
        GOLD data in dictionary format
    irr : DataFrame
        Swarm Data
    prop_keep, peaks_keep : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation.
        Only EPIs as determined by test_peaks function
    satellite : string
        Satellite letter for Swarm, 'A', 'B', or 'C'
        defualt 'SAT'
    heat_param : str
        which parameter to plot
        either 'oi' or 'deplete oi' (for barrel only)
        default 'oi'
    vmin : float
        minimum value for heatmap
        default 0
    vmax : float
        maximum value for heatmap
        default 200
    mlat_lim : NoneType or list of 2
        magnetic latitude limit for y axis
        default None, [-40, 40] will be used
    lon_lim : list of 2 values kwarg
        longitude limit for x axis
        default None, assumes GOLD span from gdc
    cmap : str kwarg
        colormap colors
        default is cividis
    facecolor : None Type or string
        if not none, plot facecolor will be set as color provided
    tformat : string
        format for strftime
        default '%d %b %H:%M UT'
    swarm_color : string
        color for Swarm trajectory and EPIs
        default "#F48FB1" (light rose)
    add_trajectory : boolean
        if True swarm trajectory will be added to heat map
        if Flase (default), the trajectory will not be added
    fs : int
        fontsize
        default 12
    """
    # create figure and gridspecs based on the number of columns
    fig = plt.figure(figsize=(9, 7))
    gs = gridspec.GridSpec(nrows=1, ncols=4, width_ratios=[1.5, 1, 0.1, 0.1],
                           wspace=0)
    plt.rcParams.update({'font.size': fs})

    if mlat_lim is None:
        mlat_lim = [-40, 40]

    # Swarm Panel -------------------------------------------------------------
    axs = fig.add_subplot(gs[0, 0])
    axs.plot(irr['ne'], irr['mlat'], color=swarm_color)

    # set axes
    axs.set_ylim(mlat_lim)
    axs.grid(axis='y')
    plotL1.format_latitude_labels(axs, xy='y', dec=0)
    axs.set_ylabel('Magnetic Latitude')
    axs.set_xlabel(r'$N_p$ (cm$^{-3}$)')
    axs.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))

    # set subtitle
    t1 = irr['time'].iloc[0].strftime('%d %b %H:%M')
    t2 = irr['time'].iloc[-1].strftime('%H:%M UT')

    dash = r'$-$'
    sw_title = f'Swarm {satellite} ({t1} {dash} {t2})'
    axs.set_title(sw_title)

    # add EPI shading
    ax_list = [axs]
    plotL1.add_peak_shading(ax_list, irr['mlat'].values, prop_keep, col='gray',
                            alpha=0.15, axis='y')

    # GOLD panel --------------------------------------------------------------
    axg = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
    cax = fig.add_subplot(gs[0, 3])
    gold_ax(axg, gdc, heat_param, vmin, vmax, cax=cax, mlat_lim=mlat_lim,
            lon_lim=lon_lim, cax_bool=True, cmap=cmap, pan='both', side='right',
            t_color='w', lon_type='lon', lat_type='mlat', facecolor=facecolor,
            tformat=tformat, fs=fs)

    # add EPIs
    axg.scatter(irr['lon'].iloc[peaks_keep], irr['mlat'].iloc[peaks_keep],
                color=swarm_color)

    if add_trajectory:
        axg.plot(irr['lon'], irr['mlat'], color=swarm_color)

    return fig
