"""Midlat trough plotting.
Updated 21 February 2026
alanahco

List Functions
--------------
trough_plot
"""
import numpy as np
from matplotlib import patches
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from pySpyEPI.plotting import plotL1
from pySpyEPI.detection.midlat_trough import zero_crossings, find_trough_lat
from pySpyEPI.utils.calc import moving_average_same_length


def trough_plot(data_df, Sat, equator_bound=35, auroral_bound=75, set_lat=55,
                trend_ylim=None, lat_textloc=1.5 * 10**5, fs=15):
    """Plot the mechanics of the find_trough_lat function (3 panels).

    Parameters
    ----------
    data_df : DataFrame
        IRR data inclduing magnetic latitude and Ne
    Sat : str
        swarm satellite A , B, or C
    equator_bound : float kwarg (35)
        equatorward magnetic latitude to search down to
        set as a positive number, will be used as +/-equator_bound
    auroral_bound : float kwarg (75)
        autroalward magnetic latitude to search up to
        set as a positive number, will be used as +/-auroral_bound
    set_lat : float kwarg (55)
        if a trough cannot be found, set_lat will be used in its place
        set as a positive number, will be used as +/-set_lat
    trend_ylim : list
        y axis limit for trend axes
        if None, no limit is set
    lat_textloc : float
        y axis location of set lat text on trend plot
        default is 1.5 * 10**5
    fs : int
        fontsize
        Default: 15

    Returns
    -------
    fig : figure handle
        figure with three panels, 2 on top and 1 on the bottom
        Top left panel: southern hemisphere filtered Ne
        Top right panel: northen hemipshere filtered Ne
        Bottom panel: Ne ranging +/- auroral_bound
    """

    # find trough lats without worrying about set_lat
    data_lat_filt = data_df[(
        (data_df["mlat"] >= -auroral_bound)
        & (data_df["mlat"] <= auroral_bound))]

    # Find trough lats
    trough_lats_any = find_trough_lat(
        data_lat_filt, equator_bound=equator_bound,
        auroral_bound=auroral_bound, set_lat=auroral_bound)

    ax_list = []
    # Create figure 3 panels 2 on top and one on the bottom
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], wspace=0.1, hspace=0.3)

    # font size
    plt.rcParams.update({'font.size': fs})

    # Southern Hemsiphere -----------------------------------------------------
    axaa = fig.add_subplot(gs[0, 0])
    ax_list.append(axaa)
    # filter data for southernn hemisphere
    data_lat_filt = data_df[(
        (data_df["mlat"] >= -auroral_bound)
        & (data_df["mlat"] <= -equator_bound))]
    param = "ne"
    param2 = "mlat"

    # filter data ------------
    window_size = 320
    moving_averages = moving_average_same_length(data_lat_filt[param].values,
                                                 window_size)
    filt10 = moving_average_same_length(data_lat_filt[param].values, 10)
    detrended = filt10 - moving_averages
    zs, ms, zi = zero_crossings(detrended, data_lat_filt["mlat"].values)

    # plot the zero crossings -------------
    axaa.plot(data_lat_filt[param2], detrended, "-", color="#A7A6BA",
              linewidth=2)
    axaa.set_ylabel(r'$N_{10s} - N_{320s}$ [$cm^{-3}$]', fontsize=fs)
    axaa.scatter(zs, np.zeros(len(zs)), marker="o",
                 color='#2B547E', label="Zero Crossings", zorder=2, s=40)
    axaa.axvline(x=-set_lat, color='k', linestyle=':',
                 linewidth=2)

    deg_south = r"$^{\circ}S$"
    axaa.text(-set_lat - 6, lat_textloc, f"{set_lat}{deg_south}", color='k')

    axaa.set_xlim(-auroral_bound, -equator_bound)

    # set legend and tick lengths
    axaa.xaxis.set_major_locator(mticker.MultipleLocator(10))
    axaa.xaxis.set_minor_locator(mticker.LinearLocator(8))
    axaa.tick_params(axis='both', which='major', labelsize=fs, width=1,
                     length=10, color='gray')
    axaa.tick_params(axis='both', which='minor', labelsize=fs, width=1,
                     length=5, color='gray')
    axaa.spines['top'].set_visible(False)
    axaa.spines['right'].set_visible(False)
    axaa.spines['bottom'].set_color('gray')
    axaa.spines['left'].set_color('gray')
    axaa.set_ylim(trend_ylim)
    plotL1.format_latitude_labels(axaa, xy='x', dec=0)
    axaa.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axaa.yaxis.grid(True)
    axaa.tick_params(right=False, labelright=False)

    # add square --------]=
    # 1. Define the point and the square size
    point_x = np.min(trough_lats_any)
    point_y = 0
    square_x = 2.0  # Width and height of the square
    square_y = np.nanmax(detrended) / 7

    # 2. Calculate the bottom-left corner of the square
    bottom_left_x = point_x - (square_x / 2)
    bottom_left_y = point_y - (square_y / 2)
    square = patches.Rectangle((bottom_left_x, bottom_left_y),
                               width=square_x, height=square_y,
                               edgecolor='#A81B45', facecolor='none',
                               linewidth=2, zorder=2)
    axaa.add_patch(square)

    # Northern Hemsiphere -----------------------------------------------------
    axaa = fig.add_subplot(gs[0, 1])
    ax_list.append(axaa)
    # filter data for northern hemisphere
    data_lat_filt = data_df[(
        (data_df["mlat"] >= equator_bound)
        & (data_df["mlat"] <= auroral_bound))]

    # filter data ----
    moving_averages = moving_average_same_length(data_lat_filt[param],
                                                 window_size)
    filt10 = moving_average_same_length(data_lat_filt[param], 10)
    detrended = filt10 - moving_averages
    zs, ms, zi = zero_crossings(detrended, data_lat_filt["mlat"].values)

    # plot the zero crossings
    axaa.plot(data_lat_filt[param2], detrended, "-",
              color="#A7A6BA", linewidth=2)  # (10s $-$ 320s) #14944f
    axaa.axvline(x=set_lat, color='k', linestyle=':',
                 linewidth=2)
    axaa.scatter(zs, np.zeros(len(zs)), marker="o",
                 color='#2B547E', label="Zero Crossings", zorder=2, s=40)
    axaa.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axaa.tick_params(left=True, labelleft=False)

    deg_north = r"$^{\circ}N$"
    axaa.text(set_lat - 6, lat_textloc, f"{set_lat}{deg_north}", color='k')

    # set plot params -------
    axaa.set_xlim(equator_bound, auroral_bound)
    axaa.legend()
    axaa.xaxis.set_major_locator(mticker.MultipleLocator(10))
    axaa.xaxis.set_minor_locator(mticker.LinearLocator(8))
    axaa.tick_params(axis='both', which='major', labelsize=fs, width=1,
                     length=10, color='gray')
    axaa.tick_params(axis='both', which='minor', labelsize=fs, width=1,
                     length=5, color='gray')
    axaa.spines['top'].set_visible(False)
    axaa.spines['right'].set_visible(False)
    axaa.spines['bottom'].set_color('gray')
    axaa.spines['left'].set_color('gray')
    axaa.set_ylim(trend_ylim)
    axaa.yaxis.grid(True)
    plotL1.format_latitude_labels(axaa, xy='x', dec=0)

    # add square --------
    # 1. Define the point and the square size
    point_x = np.max(trough_lats_any)
    point_y = 0
    square_x = 2.0  # Width and height of the square

    # 2. Calculate the bottom-left corner of the square
    bottom_left_x = point_x - (square_x / 2)
    bottom_left_y = point_y - (square_y / 2)

    square = patches.Rectangle((bottom_left_x, bottom_left_y),
                               width=square_x, height=square_y,
                               edgecolor='#A81B45', facecolor='none',
                               linewidth=2, zorder=2)
    axaa.add_patch(square)

    # All latitudes -----------------------------------------------------------
    axaa = fig.add_subplot(gs[1, :])
    ax_list.append(axaa)
    data_lat_filt = data_df[(
        (data_df["mlat"] >= -auroral_bound)
        & (data_df["mlat"] <= auroral_bound))]

    # Find trough lats
    trough_lats = find_trough_lat(data_lat_filt, equator_bound=equator_bound,
                                  auroral_bound=auroral_bound, set_lat=set_lat)

    # plot all
    axaa.plot(data_lat_filt[param2], data_lat_filt[param], "-",
              color="#F48FB1", linewidth=2)
    axaa.axvline(x=trough_lats[0], color='#A81B45',
                 linestyle='--', linewidth=2)
    axaa.axvline(x=trough_lats[1], color='#A81B45', linestyle='--',
                 linewidth=2)
    axaa.set_ylabel(param, fontsize=fs)
    axaa.set_xlim(-auroral_bound, auroral_bound)

    axaa.set_ylabel(r'$N_p$ [$cm^{-3}$]', fontsize=fs)
    axaa.set_xlabel('Magnetic Latitude', fontsize=fs)

    axaa.xaxis.set_major_locator(mticker.MultipleLocator(15))
    axaa.xaxis.set_minor_locator(mticker.LinearLocator(29))
    axaa.tick_params(axis='both', which='major', labelsize=fs, width=1,
                     length=10, color='gray')
    axaa.tick_params(axis='both', which='minor', labelsize=fs, width=1,
                     length=5, color='gray')
    axaa.spines['top'].set_visible(False)
    axaa.spines['right'].set_visible(False)
    axaa.spines['bottom'].set_color('gray')
    axaa.spines['left'].set_color('gray')
    axaa.xaxis.grid(True)
    plotL1.format_latitude_labels(axaa, xy='x', dec=0)
    axaa.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    # plot title
    dash = '$-$'
    f_time_1 = data_df["time"].iloc[0].strftime('%H:%M')
    f_time_2 = data_df["time"].iloc[-1].strftime('%H:%M')
    title_date = data_df['time'].iloc[0].strftime('%d %b %Y')
    titl = f'{title_date} {f_time_1} UT {dash} {f_time_2}'
    plt.suptitle(f'Swarm {Sat} {titl} UT',
                 y=0.92, fontsize=fs)
    return fig, ax_list
