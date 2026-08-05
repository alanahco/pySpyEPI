"""Plotting Conjunction figures for GOLD and arbitrary satellite.
Created : 07 July 2026
By Alanah Cardenas-O'Toole
List Functions
--------------
gold_conj_pans
gold_sat_plot
gold_sat_dailies
"""
import datetime as dt
import numpy as np
from pySpyEPI.validation import conjunctions
from pySpyEPI.io.load import gold_load
import matplotlib.pyplot as plt
from pySpyEPI.plotting import gold_plot as gplot
import cartopy.crs as ccrs
import cmocean
import os
from pathlib import Path


def gold_conj_pans(map_ax, cbar_ax, line_ax, gdc, data_lim, center_time=None,
                   edge1_mlat=None, edge2_mlat=None, cmap=None, vmin=None,
                   vmax=None, mlat_lim=None, lon_lim=None, sat_time_str='time',
                   sat_lon_str='lon', sat_lat_str='lat', sat_mlat_str='mlat',
                   sat_n_str='ne', sat_name='Satellite', gold_line=True,
                   fs=12):
    """GOLD satellite conjunction map and line plot panels.

    Parameters
    ----------
    map_ax, cbar_ax, line_ax : axes handles
        map axis handle fr heatmap
        color bar axis
        line axis for satellite density and GOLD OI
    gdc : dictionary
        GOLD data including columns 'oi', 'lat', and 'lon'
    data_lim : dataframe
        satellite data for designated period
    center_time : array-like or NoneType
        EPI center Time
        Default : None, no EPIs plotted
    edge1_mlat : array-like  or NoneType
        EPI edge 1 magnetic latitude
        Default : None, no EPI width plotted
    edge2_mlat : array-like or NoneType
        EPI edge 2 magnetic latitude
        Default : None, no EPI width plotted
    cmap : color bar type or string or NoneType
        color bar for heatmap
        if None (default), cmap = cmocean.cm.algae_r
    vmin, vmax : float or NoneType
        color range
        if None (default), vmin = 5 and vmax = maximum of GOLD OI
    mlat_lim, lon_lim : array-like or NoneType
        magnetic latitude limit for plot default is None [-40, 40]
        longiutde limit, default is derived from gdc['lon']
    sat_{}_str : strings
        satellite parameter column string names for time, lon, lat, mlat, and n
    sat_name : string
        Satellite name for labeling
    gold_line : boolean
        if True (default), GOLD line plot will be added to line_ax
    fs : float
        fontsize
        default 12
    """

    # Establish None params
    if cmap is None:
        cmap = cmocean.cm.algae_r
    if vmin is None:
        vmin = 5
    if vmax is None:
        vmax = np.nanmax(gdc['oi'])
    if mlat_lim is None:
        mlat_lim = [-40, 40]

    # GOLD axis ---------------------------------------------------------------
    gplot.gold_ax(
        map_ax, gdc, 'oi', vmin, vmax, cax=cbar_ax, mlat_lim=mlat_lim,
        lon_lim=lon_lim, cax_bool=True, cmap=cmap, pan='both', side='left',
        t_color='w', lon_type='lon', lat_type='mlat', cax_loc='top', fs=fs)
    map_ax.set_facecolor('gray')

    # add sat trajectory ------------------------------------------------------
    map_ax.plot(data_lim[sat_lon_str], data_lim[sat_mlat_str], color="#F48FB1",
                linewidth=2.5)

    # add EPI detections
    if center_time is not None:
        epi_clim = center_time[(
            (center_time >= data_lim[sat_time_str].iloc[0])
            & (center_time <= data_lim[sat_time_str].iloc[-1]))]
        elons = []
        emlats = []
        for ctime in epi_clim:
            ei = np.argmin(np.abs(ctime - data_lim[sat_time_str]))
            elons.append(data_lim[sat_lon_str].iloc[ei])
            emlats.append(data_lim[sat_mlat_str].iloc[ei])

        map_ax.scatter(elons, emlats, color="C1", zorder=2)

    # Line_plot ---------------------------------------------------------------

    # plot satellite density
    if gold_line:  # add gold
        sx = line_ax.twiny()
        sx.plot(data_lim[sat_n_str], data_lim[sat_mlat_str], color="#F48FB1")

        gold_ois, gold_mlats = conjunctions.gold_by_lon(gdc,
                                                        data_lim[sat_lat_str],
                                                        data_lim[sat_lon_str])
        line_ax.plot(gold_ois, gold_mlats, color='#09834c', linestyle='--')

        line_ax.set_ylim(mlat_lim)
        line_ax.tick_params(left=True, labelleft=False)
        line_ax.yaxis.grid()
        sx.set_xlabel(fr'{sat_name} $N_p$ [$cm^{-3}$]', color="#F06292")
        line_ax.set_xlabel('GOLD OI (R)', color="#09834c")
        sx.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))

        # add sat time text
        mini = np.argmin(data_lim[sat_mlat_str])
        south_t = data_lim[sat_time_str].iloc[mini].strftime('%H:%M UT')
        maxi = np.argmax(data_lim[sat_mlat_str])
        north_t = data_lim[sat_time_str].iloc[maxi].strftime('%H:%M UT')

        sx.text(0.5, 0.98, north_t, transform=sx.transAxes, ha='center',
                va='top', color="#F06292", fontsize=fs, fontweight='bold')
        sx.text(0.5, 0.05, south_t, transform=sx.transAxes, ha='center',
                va='top', color="#F06292", fontsize=fs, fontweight='bold')

        sx.set_title(data_lim[sat_time_str].iloc[0].strftime('%d %B %Y'),
                     fontsize=fs)
        # add EPI width shading
        if (edge1_mlat is not None) and (edge2_mlat is not None):
            epi_llim = edge1_mlat[(
                (center_time >= data_lim[sat_time_str].iloc[0])
                & (center_time <= data_lim[sat_time_str].iloc[-1]))]
            epi_rlim = edge2_mlat[(
                (center_time >= data_lim[sat_time_str].iloc[0])
                & (center_time <= data_lim[sat_time_str].iloc[-1]))]
            for xl, xr in zip(epi_llim, epi_rlim):
                sx.axhspan(xl, xr, color='gray', alpha=0.15, zorder=0)
    else:
        line_ax.plot(data_lim[sat_n_str], data_lim[sat_mlat_str],
                     color="#F48FB1")

        line_ax.set_ylim(mlat_lim)
        line_ax.tick_params(left=True, labelleft=False)
        line_ax.yaxis.grid()
        unit = '$cm^{-3}$'
        line_ax.set_xlabel(f'{sat_name} $N_p$ [{unit}]')
        line_ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))

        # add sat time text
        mini = np.argmin(data_lim[sat_mlat_str])
        south_t = data_lim[sat_time_str].iloc[mini].strftime('%H:%M UT')
        maxi = np.argmax(data_lim[sat_mlat_str])
        north_t = data_lim[sat_time_str].iloc[maxi].strftime('%H:%M UT')

        line_ax.text(0.5, 0.98, north_t, transform=line_ax.transAxes,
                     ha='center', va='top', color="#F06292", fontsize=fs,
                     fontweight='bold')
        line_ax.text(0.5, 0.05, south_t, transform=line_ax.transAxes,
                     ha='center', va='top', color="#F06292", fontsize=fs,
                     fontweight='bold')

        line_ax.set_title(data_lim[sat_time_str].iloc[0].strftime('%d %b %Y'),
                          fontsize=fs)

        # add EPI width shading
        if (edge1_mlat is not None) and (edge2_mlat is not None):
            epi_llim = edge1_mlat[(
                (center_time >= data_lim[sat_time_str].iloc[0])
                & (center_time <= data_lim[sat_time_str].iloc[-1]))]
            epi_rlim = edge2_mlat[(
                (center_time >= data_lim[sat_time_str].iloc[0])
                & (center_time <= data_lim[sat_time_str].iloc[-1]))]
            for xl, xr in zip(epi_llim, epi_rlim):
                line_ax.axhspan(xl, xr, color="#A7A6BA", alpha=0.25, zorder=0)


def gold_sat_plot(conjs, sat_data, gold_fold, center_time, edge1_mlat,
                  edge2_mlat, cmap=None, vmin=None, vmax=None,
                  mlat_lim=None, lon_lim=None, sat_time_str='time',
                  sat_lon_str='lon', sat_lat_str='lat', sat_mlat_str='mlat',
                  sat_n_str='ne', sat_name='Satellite', fs=12):
    """GOLD plus satellite conjunction 4 panel plot.

    Parameters
    ----------
    conjs : pd.DataFrame
        conjunction starting time (sat_t1) and ending time (sat_t2), satellite
        starting (sat_lon1) and ending (sat_lon2) longitudes,
        gold time (gold_time), and gold starting (gold_lon1) and ending
        (gold_lon2) longitudes
    sat_data : pd.DataFrame
        satellite data used for creating conjs
    gold_fold : string
        folder where gold data is
    center_time : array-like
        EPI center Time
    edge1_mlat : array-like
        EPI edge 1 magnetic latitude
    edge2_mlat : array-like
        EPI edge 2 magnetic latitude
    cmap : color bar type or string or NoneType
        color bar for heatmap
        if None (default), cmap = cmocean.cm.algae_r
    vmin, vmax : float or NoneType
        color range
        if None (default), vmin = 5 and vmax = maximum of GOLD OI
    mlat_lim, lon_lim : array-like or NoneType
        magnetic latitude limit for plot default is None [-40, 40]
        longiutde limit, default is derived from gdc['lon']
    sat_{}_str : strings
        satellite parameter column string names for time, lon, lat, mlat, and n
    sat_name : string
        Satellite name for labeling

    Returns
    -------
    fig : figure handle
        figure
    ax_list : axis handles
        for adding alphabet labels later
    """

    # Create Figure ----
    fig = plt.figure(figsize=(16, 16))
    plt.rcParams.update({'font.size': fs})
    # Master grid--------------------------------------------------------------
    gs_master = fig.add_gridspec(2, 2)
    ax_list = []
    for j in range(len(conjs)):
        # Open GOLD FILE
        tc = conjs['sat_t1'].iloc[j]
        st = dt.datetime(tc.year, tc.month, tc.day, tc.hour, tc.minute)
        gdc = gold_load.open_GOLD_NObarrel(st, fdir=gold_fold)

        # limit satelltie data
        data_lim = sat_data.copy()
        data_lim = data_lim[(
            (data_lim[sat_time_str] >= conjs['sat_t1'].iloc[j])
            & (data_lim[sat_time_str] <= conjs['sat_t2'].iloc[j]))]

        # Subfigure -----------------------------------------------------------
        gs_sub = gs_master[j].subgridspec(
            nrows=2, ncols=2, height_ratios=[0.1, 1], hspace=0.08, wspace=0.08)
        gs_gold = gs_sub[1, 0]
        cs_gold = gs_sub[0, 0]
        gs_sat = gs_sub[1, 1]

        # GOLD axis -----------------------------------------------------------
        gx = fig.add_subplot(gs_gold, projection=ccrs.PlateCarree())
        cax = fig.add_subplot(cs_gold)

        ax_list.append(gx)
        # Sat Panel -----------------------------------------------------------
        ax = fig.add_subplot(gs_sat)
        ax_list.append(ax)
        # plot the data
        gold_conj_pans(gx, cax, ax, gdc, data_lim, center_time, edge1_mlat,
                       edge2_mlat, cmap=cmap, vmin=vmin, vmax=vmax,
                       mlat_lim=mlat_lim, lon_lim=lon_lim,
                       sat_time_str=sat_time_str, sat_lon_str=sat_lon_str,
                       sat_lat_str=sat_lat_str, sat_mlat_str=sat_mlat_str,
                       sat_n_str=sat_n_str, sat_name=sat_name, fs=fs)

    return fig, ax_list


def gold_sat_dailies(sat_data, gold_fold, center_time, edge1_mlat, edge2_mlat,
                     fig_dir, max_mindif=60, fig_all=True, cmap=None,
                     vmin=None, vmax=None, mlat_lim=None, lon_lim=None,
                     sat_time_str='time', sat_lon_str='lon', sat_lat_str='lat',
                     sat_mlat_str='mlat', sat_n_str='ne',
                     sat_name='Satellite', fs=12):
    """GOLD plus satellite conjunction daily plots.

    Parameters
    ----------
    sat_data : pd.DataFrame
        satellite data used for creating conjs
    gold_fold : string
        folder where gold data is
    center_time : array-like
        EPI center Time
    edge1_mlat : array-like
        EPI edge 1 magnetic latitude
    edge2_mlat : array-like
        EPI edge 2 magnetic latitude
    fig_dir : string
        figure directory
    max_mindif : float
        maximum minute difference between satellite time and GOLD time
        default 60
    fig_all: boolean
        if True (default), figures will include all possible conjunctions
        if False, figures will only include non-repeated latest and fullest
        in time conjunctions
    cmap : color bar type or string or NoneType
        color bar for heatmap
        if None (default), cmap = cmocean.cm.algae_r
    vmin, vmax : float or NoneType
        color range
        if None (default), vmin = 5 and vmax = maximum of GOLD OI
    mlat_lim, lon_lim : array-like or NoneType
        magnetic latitude limit for plot default is None [-40, 40]
        longiutde limit, default is derived from gdc['lon']
    sat_{}_str : strings
        satellite parameter column string names for time, lon, lat, mlat, and n
    sat_name : string
        Satellite name for labeling
        default "Satellite"

    Returns
    -------
    fig : figure handle
        figure
    ax_list : axis handles
        for adding alphabet labels later
    """

    # saving strings
    sname = sat_name.replace(" ", "")
    stime = sat_data[sat_time_str].iloc[0]
    sday = stime.strftime('%d%B%Y')
    date_dir = os.path.join(fig_dir, stime.strftime('%Y'),
                            stime.strftime('%m'), stime.strftime('%Y%m%d'))

    # create GOLD files
    gold_load.GOLD_no_barrel(stime, gold_fold, lon_range=None,
                             closest=False)
    etime = stime + dt.timedelta(days=1)
    gold_load.GOLD_no_barrel(etime, gold_fold, lon_range=None,
                             closest=False)

    # get conjunctions from sat_data
    condf = conjunctions.goldcon(
        sat_data[sat_time_str], sat_data[sat_mlat_str], sat_data[sat_lon_str],
        max_mindif=max_mindif)

    if len(condf) > 0:
        Path(date_dir).mkdir(parents=True, exist_ok=True)

    # create figures for all conjunctions or only non-repeats based on fig_all
    if fig_all:
        lc = len(condf)
        num_figs = np.ceil(lc / 4)
        for i in range(int(num_figs)):
            i1 = i * 4
            i2 = (i + 1) * 4
            if i2 > lc:
                i2 = lc
            condf_plot = condf.iloc[i1:i2]

            fig, ax_list = gold_sat_plot(
                condf_plot, sat_data, gold_fold, center_time, edge1_mlat,
                edge2_mlat, cmap=cmap, vmin=vmin, vmax=vmax, mlat_lim=mlat_lim,
                lon_lim=lon_lim, sat_time_str=sat_time_str,
                sat_lon_str=sat_lon_str, sat_lat_str=sat_lat_str,
                sat_mlat_str=sat_mlat_str, sat_n_str=sat_n_str,
                sat_name=sat_name, fs=fs)

            # save figs
            fig_name = f'GOLD_{sname}_conjunction_{sday}_{i}.png'
            fig_save_dir = os.path.join(date_dir, fig_name)
            fig.savefig(fig_save_dir, bbox_inches='tight', dpi=300)
            plt.close(fig)
    else:
        condf_lim = conjunctions.reduce_conj(condf)

        lc = len(condf_lim)
        num_figs = np.ceil(lc / 4)
        for i in range(int(num_figs)):
            i1 = i * 4
            i2 = (i + 1) * 4
            if i2 > lc:
                i2 = lc
            condf_plot = condf_lim.iloc[i1:i2]

            fig, ax_list = gold_sat_plot(
                condf_plot, sat_data, gold_fold, center_time, edge1_mlat,
                edge2_mlat, cmap=cmap, vmin=vmin, vmax=vmax, mlat_lim=mlat_lim,
                lon_lim=lon_lim, sat_time_str=sat_time_str,
                sat_lon_str=sat_lon_str, sat_lat_str=sat_lat_str,
                sat_mlat_str=sat_mlat_str, sat_n_str=sat_n_str,
                sat_name=sat_name, fs=fs)

            fig_name = f'GOLD_{sname}_reduced_conjunction_{sday}_{i}.png'
            fig_save_dir = os.path.join(date_dir, fig_name)
            fig.savefig(fig_save_dir, bbox_inches='tight', dpi=300)
            plt.close(fig)
