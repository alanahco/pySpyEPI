""" Level 1 plotting functions
Updated 20 February 2026
List of Functions:
    - format_latitude_labels
    - plot_maglat_panel
    - add_lat_lines
    - adjust_axes_time
    - make_legend
    - add_peak_shading
    - add_alphalabels
    - get_terminator
    - color_gradient
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pandas as pd
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import cmocean
import pydarn
import matplotlib.colors as mcolors
import colorsys


def format_latitude_labels(ax, xy='x', dec=0):
    """
    Formats the latitude axis labels with degree symbols and N/S suffixes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The Matplotlib axes object
    xy : str kwarg
        'x' (defualt) or 'y' depending on which axis you want to have
        degree symbol N/S formatting
    dec : float
        number of decimal places, if 0, dec will be .0f, if anything else,
        deccimals will be .1f
    """

    def latitude_formatter(latitude, pos):
        if latitude > 0:
            if dec == 0:
                return f"{latitude:.0f}°N"
            else:
                return f"{latitude:.1f}°N"
        elif latitude < 0:
            if dec == 0:
                return f"{abs(latitude):.0f}°S"
            else:
                return f"{abs(latitude):.0f}°S"
        else:
            return "0°"
    if xy == 'x':
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(latitude_formatter))
    elif xy == 'y':
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(latitude_formatter))


def format_longitude_labels(ax, xy='x', dec=0):
    """Format the longitude axis labels with degree symbols and E/W suffixes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The Matplotlib axes object
    xy : str kwarg
        'x' (defualt) or 'y' depending on which axis you want to have
        degree symbol E/W formatting
    dec : int
        = 0, no decimal, = 1, 1 decimal
    """

    def longitude_formatter(longitude, pos):
        if longitude > 0:
            if dec == 0:
                return f"{longitude:.0f}°E"
            else:
                return f"{longitude:.1f}°E"
        elif longitude < 0:
            if dec == 0:
                return f"{abs(longitude):.0f}°W"
            else:
                return f"{abs(longitude):.1f}°W"
        else:
            return "0°"
    if xy == 'x':
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(longitude_formatter))
    elif xy == 'y':
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(longitude_formatter))


def plot_maglat_panel(data_df, ax=None, cax=None, fs=15,
                      cmap_name=cmocean.cm.phase, glat_bool=True,
                      time_str='time', lon_str='lon',
                      lt_hr_str='lt_hour', mlat_str='mlat', lt_str='lt',
                      lat_str='lat', lat_lines=True, maj_num=11):
    """Plots Mag_Lat vs Longitude.

    Parameters
    ----------
    data_df : pd.DataFrame
        Must contain magnetic latitude, geo longitude, local time, geo lat, and
        time
    ax : matplotlib.axes.Axes or None
        If provided, plots into this axis. Otherwise creates a new one.
    cax : matplotlib.axes.Axes or None
        If provided, places colorbar in this axis.
    fs : int
        Font size for labels.
    cmap_name : color bar
        Name of colormap (default 'viridis')
    glat_bool : boolean kwrg True
        set to True will plot top x axis as geographic latitude
    time_str, lon_str, lt_hr_str, mlat_str, lt_str, lat_str : strings
        variable names
    lat_lines : bool
        if True, lat lines are plotted (default)
        if False, not plotted
    maj_num : int
        number of major ticks
        default 11
    Returns
    -------
    ax : matplotlib.axes.Axes
        The matplotlib axis containing the plot.

    Notes
    -----
        - LT_hour as line gradient
        - Magnetic Lat on bottom axis
        - Geographic Lat on top axis
        - Red line marking Geo_Lat = 0°
    """

    # Compute LT_hour if not already present
    if lt_hr_str not in data_df.columns:
        data_df[lt_hr_str] = data_df[lt_str].apply(
            lambda tim: tim.hour + tim.minute / 60 + tim.second / 3600)

    # Use time for x-spacing
    x = data_df[time_str]
    y = data_df[lon_str]
    lt = data_df[lt_hr_str]

    # Prepare figure
    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 5))

    # Scatter plot
    ax.scatter(x, y, c=lt, cmap=cmap_name,
               norm=Normalize(vmin=lt.min(), vmax=lt.max()),
               marker='.', lw=0.05)

    ax.set_ylabel("Geographic Longitude", fontsize=fs)

    # format longitude labels
    if abs(min(y)) - abs(max(y)) < 4:
        dec = 1
    else:
        dec = 0
    format_longitude_labels(ax, xy='y', dec=dec)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_xlim([data_df[time_str].iloc[0], data_df[time_str].iloc[-1]])

    # Set x-ticks using time, but label with magnetic latitude
    xtick_locs = np.round(np.linspace(0, len(x) - 1, maj_num)).astype(int)
    xtimes = x.iloc[xtick_locs]
    xlabels = data_df[mlat_str].iloc[xtick_locs].round(2)
    ax.set_xticks(xtimes)

    xlab_NS = []
    for lab_c in xlabels:
        if lab_c > 0:
            xlab_NS.append(f"{lab_c:.1f}°N")
        elif lab_c < 0:
            xlab_NS.append(f"{abs(lab_c):.1f}°S")
        else:
            xlab_NS.append("0°")

    ax.set_xticklabels(xlab_NS)

    ax.set_xlabel("Magnetic Latitude", fontsize=fs)
    ax.tick_params(axis='x', labelsize=fs)

    # Optional: add secondary top axis for geographic latitude
    if glat_bool:
        secax = ax.twiny()
        secax.set_xlim(ax.get_xlim())
        secax.set_xticks(xtimes)

        xgeo_NS = []
        for lab_c in data_df[lat_str].iloc[xtick_locs].round(4):
            if lab_c > 0:
                xgeo_NS.append(f"{lab_c:.1f}°N")
            elif lab_c < 0:
                xgeo_NS.append(f"{abs(lab_c):.1f}°S")
            else:
                xgeo_NS.append("0°")
        secax.set_xticklabels(xgeo_NS)
        secax.set_xlabel("Geographic Latitude", fontsize=fs)
        secax.tick_params(axis='x', labelsize=fs)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_name,
                               norm=Normalize(vmin=lt.min(), vmax=lt.max()))
    sm.set_array([])
    if cax is not None:
        cb = ax.figure.colorbar(sm, cax=cax)
    else:
        cb = ax.figure.colorbar(sm, ax=ax)

    # set tick labels for color bar
    tick_locator = mticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()

    # add label
    cb.set_label("Local Time [h]", fontsize=fs, rotation=270, labelpad=20)
    cb.ax.tick_params(labelsize=fs)

    ax.tick_params(axis='both', which='major', labelsize=fs, width=1.5,
                   length=10)
    ax.xaxis.grid()

    if lat_lines:
        # Add vertical line at Geo_Lat and Mag_Lat 0°
        lat_val = 0
        add_lat_lines(data_df, ax, lat_val, lat_type='mlat')
        add_lat_lines(data_df, ax, lat_val, lat_type='lat')

        # add magnetic latitude lines
        lat_val = 25
        add_lat_lines(data_df, ax, lat_val, lat_type='mlat')
        ax.legend(loc='best', fontsize=fs - 2)

    return ax


def add_lat_lines(data_df, ax, lat_val, lat_type='mlat', col=None,
                  dir='NS', style=None, t_str='time'):
    """ Add latitude lines to a plot
    Parameters
    ----------
    data_df : dataframe
        dataframe that is used for axis data
    ax : axis
    lat_val : double
        desired latitude location
    lat_type : str kwarg
        default is mlat for magnetic latitude
        other option is lat for geographic latitude
    col : NoneType or string
        color for lat lines
        if None, magenta for maglat and red for geo lat
    dir : string
        NS, lat_val in both hemispheres, N or S in one hemisphere
    style : str
        line style for 0 degrees
        if None (default), -- is used for mlat and : is used for lat
    Retruns
    -------
    updated axis
    """

    # update params based on what is given
    if lat_type == 'mlat':
        lab = 'MLat'
        if col is None:
            col = 'magenta'
    elif lat_type == 'lat':
        lab = 'GLat'
        if col is None:
            col = '#008080'

    if lat_val != 0:
        if lat_type == 'mlat':
            if style is None:
                style = '-'
        elif lat_type == 'lat':
            lab = 'GLat'
            if style is None:
                style = '--'
        # Find closest location closest to + latitude input
        if dir == 'NS':
            h = np.where(
                min(data_df[lat_type], key=lambda x: abs(x + lat_val))
                == data_df[lat_type])[0][0]

            if abs(min(data_df[lat_type],
                       key=lambda x: abs(x + lat_val)) + lat_val) < 0.5:
                ax.axvline(x=data_df[t_str].iloc[h], color=col,
                           linestyle=style, linewidth=2.5, alpha=0.5)

        # Find closest location closest to - latitude input
        g = np.where(
            min(data_df[lat_type], key=lambda x: abs(x - lat_val))
            == data_df[lat_type])[0][0]
        pm = '\u00B1'
        if abs(min(data_df[lat_type],
                   key=lambda x: abs(x - lat_val)) - lat_val) < 0.5:
            ax.axvline(x=data_df[t_str].iloc[g], color=col,
                       linestyle=style, linewidth=2.5, alpha=0.5,
                       label=f'{pm}{abs(lat_val)}° {lab}')
    else:
        if lat_type == 'mlat':
            if style is None:
                style = '--'
        elif lat_type == 'lat':
            if style is None:
                style = ':'
        # Add vertical line at 0°
        zero_lat = data_df[lat_type].to_numpy()
        zero_idx = np.argmin(np.abs(zero_lat))

        time_at_zero = data_df[t_str].iloc[zero_idx]
        ax.axvline(time_at_zero, color=col, linestyle=style, lw=2.5,
                   alpha=0.5, label=f'0° {lab}')


def adjust_axes_time(x_val, ax, xtime, xlab='x', ylab='y',
                     time_bool=False, leg=True, fs=14, x_on=True, sci_y=False,
                     minor_ticks=True, xgrid=True, maj_num=11):
    """ Adjust axes based on TIME.
    Parameters
    ----------
    x_val : array-like
        tick markers for x axis
    ax : axis handle
    xtick_loc : array-like
        locations of x ticks
    xtime : array-like
        time of x array
    xlab : kwarg str
        x axis label default 'x'
    ylab : kwarg str
        y axis label default 'y'
    time_bool : boolean kwarg
        default False, if the x axis is time
    leg : boolean kwarg
        set legend if True (default)
    fs : int kwarg (14)
        fontsize
    x_on : bool kwarg (True)
        if True, xlabels will be plotted
    sci_y : bool
        y axis labels in scientific notation
        default False (pre-set notation)
    minor_ticks : bool
        adds minor ticks if True (default)
    xgrid : bool
        adds xgrid if True (default)
    maj_num : int
        number of major ticks
        default 11

    Returns
    -------
    adjusted axes
    """

    # set legend if specified
    if leg:
        ax.legend()

    # restric axes
    ax.set_xlim([xtime[0], xtime[-1]])
    xtick_loc = np.round(np.linspace(0, len(x_val) - 1, maj_num)).astype(int)

    # set_ticks
    ax.set_xticks(xtime[xtick_loc])

    if minor_ticks:
        # Get the major tick locations
        major_ticks = ax.get_xticks()

        # Compute the midpoints
        minor_ticks = (major_ticks[:-1] + major_ticks[1:]) / 2

        # Set the minor ticks
        ax.xaxis.set_minor_locator(mticker.FixedLocator(minor_ticks))

    # set y axis label
    ax.set_ylabel(ylab, fontsize=fs)

    if x_on:

        if time_bool:
            # set time ticks
            t_labels = pd.to_datetime(x_val[xtick_loc]).strftime("%H:%M")
            ax.set_xticklabels(t_labels)
        else:
            # round to two value places
            ax.set_xticklabels(np.round(x_val[xtick_loc], 2))

        ax.set_xlabel(xlab, fontsize=fs)
    else:
        ax.set_xticklabels([])

    # set tick lengths
    ax.tick_params(axis='both', which='major', labelsize=fs, width=1.5,
                   length=10)
    ax.tick_params(axis='both', which='minor', labelsize=fs, width=1.5,
                   length=5)
    if xgrid:
        ax.xaxis.grid()

    if sci_y:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(fs)


def make_legend(leg_ax, leg_labs, leg_cols, leg_styles, modes, ls_marker=None,
                leg_alphas=None, lwidths=None, **kwargs):
    """
    Create a custom legend on a given axis. Chat GPT aid

    Parameters
    ----------
    leg_ax : matplotlib axis
        Axis to place the legend on.
    leg_labs : list of str
        Labels for the legend entries.
    leg_cols : list of str
        Colors for the legend entries.
    leg_styles : list of str
        Marker styles (if scatter) or line styles (if line).
    modes : list of str
        Type of legend entry for each label.
        Options: "line", "scatter", "shading", "line+shading", and
        "line+scatter", "hollowmarker" (maybe)
    leg_alphas : None type or list-like
        if None (defualt), line alphas assumed to be 1
        if provided, alpha (transparency) will be used for legend
    ls_marker : None type or list-like
        if None (defualt), scatter mark for line+scatter is assumed to be 'o'
        if provided as a list of strings, different markers will be used for
        line+scatter
    lwidths : None type or list like
        if None, all linewidths are set to 1
        if provided as a list of floats used for linewidths
    kwargs : dict
        Additional keyword arguments passed to ax.legend().
    """
    handles = []

    # fill alphas with ones
    if leg_alphas is None:
        leg_alphas = np.ones(len(leg_labs))

    # if ls_marker is not provided assumed 'o'
    if ls_marker is None:
        ls_marker = ['o'] * len(leg_labs)

    if lwidths is None:
        lwidths = np.ones(len(leg_labs))

    # go through different modes and add info for legend
    for lab, col, style, mode, alph, ls_m, lw in zip(
            leg_labs, leg_cols, leg_styles, modes, leg_alphas, ls_marker,
            lwidths):
        if mode == "line":
            h = Line2D([], [], color=col, linestyle=style, marker=None,
                       label=lab, alpha=alph, linewidth=lw)
        elif mode == "scatter":
            h = Line2D([], [], marker=style, linestyle="None",
                       markerfacecolor=col, markeredgecolor=col,
                       markersize=8, label=lab, alpha=alph)
        elif mode == "shading":
            h = Patch(facecolor=col, edgecolor="none", alpha=alph, label=lab)
        elif mode == "line+shading":
            line = Line2D([], [], color=col, linestyle=style, linewidth=lw)
            patch = Patch(facecolor=col, edgecolor="none", alpha=alph)
            h = (line, patch)  # composite handle
        elif mode == "line+scatter":
            h = Line2D([], [], color=col, linestyle=style, marker=ls_m,
                       label=lab, alpha=alph, linewidth=lw)
        elif mode == "errorbar":
            line = Line2D([], [], marker=style, linestyle="None",
                          markerfacecolor=col, markeredgecolor=col,
                          markersize=8, label=lab, alpha=alph)

            # vertical error bar (no horizontal caps)
            err = Line2D([], [], color=col, linestyle='None',
                         marker='|', markersize=20, markeredgewidth=1.5,
                         alpha=alph)
            h = (line, err)
        else:
            raise ValueError(f"Unknown mode: {mode}")
        handles.append(h)

    leg_ax.legend(handles=handles, labels=leg_labs, **kwargs)
    leg_ax.axis("off")  # hide axes completely)


def add_peak_shading(ax_list, x_data, properties, col='yellow', alpha=0.1,
                     axis='x'):
    """ Add shading to each panel in ax_list
    Parameters
    ----------
    ax_list : list
        list of figure axes from fig.axes
    x_data : array-like
        dataset for x axis
    properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
        used for shading
    col : kwarg str ('yellow'
        color of shading
    axis : string
        horizontal ('y') or vertical shading ('x')
        defualt 'x'
    Example
    -------
    # barrel roll
    barrel_start=8
    det_filt=5
    prom_perc=10
    peak_width=3

    x_secdf = (IRR_data['Time'].dt.hour*3600 + IRR_data['Time'].dt.minute*60
               + IRR_data['Time'].dt.second)
    x_sec = x_secdf.values
    ne_copy = IRR_data['Ne'].copy()
    ne = ne_copy.values
    ne[IRR_data["Ne_flag"] > 2] = np.nan
    ne_df, peaks, properties = get_barrel_info(x_sec, ne,
                                               barrel_start=barrel_start,
                                               det_filt=det_filt,
                                               prom_perc=prom_perc,
                                               width=peak_width)

    # Get Peak flags
    peak_flag = process_peaks(IRR_data, peaks, properties)

    mask_keep = (peak_flag == 1)
    peaks_keep = peaks[mask_keep]
    prop_keep = (
        {key: value[mask_keep] for key, value in properties.items()
         if isinstance(value, np.ndarray)})

    # Get figure axes
    ax_list = fig.axes
    axs = [ax_list[0], ax_list[3], ax_list[5], ax_list[7]]
    add_peak_shading(axs, IRR_data['Time'].values, prop_keep, col='yellow')
    """
    # get the left and right sides to shade
    xl = np.ceil(properties["left_ips"]).astype(int)
    xr = np.trunc(properties["right_ips"]).astype(int)

    # iterate through axes list
    for ax in ax_list:

        # iterate through peaks edges
        for x in range(len(xl)):
            if axis == 'x':
                ax.axvspan(x_data[xl[x]], x_data[xr[x]], color=col,
                           alpha=alpha, zorder=0)
            elif axis == 'y':
                ax.axhspan(x_data[xl[x]], x_data[xr[x]], color=col,
                           alpha=alpha, zorder=0)


def add_alphalabels(ax_list, inout='in', color='black', yperc=5, xperc=1,
                    start_num=0, bold=False, fs=14):
    """Add alphabetical labels to a plot for publication purposes.

    Parameters
    ----------
    ax_list : list-like
        list of axes handles
    inout : string
        'in' or 'out' depending on if the labels should be on the inside or
        outside of the panel
        default 'in'
    color : string
        text color
        default 'black'
    yperc : float
        location  of label based on percent of y limit from max
        loc = ymax +/- yperc*ydif/100
        + or - based on in or out
        default 5%
    xperc : float
        location  of label based on percent of xlim form lower limit
        loc = xmin + xperc*xdif/100
        + or - based on in or out
        default 1%
    start_num : int
        if 0 (default) will start at (a). if 1, will start at (b) and so on
    bold : boolean
        if True, labels will be bold
        default False
    fs : float
        fontsize
        default 14
    """
    # Iterate through ax_list
    alph_list = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)',
                 '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)',
                 '(q)', '(r)', '(s)', '(t)', '(u)', '(v)', '(w)', '(x)',
                 '(y)', '(z)']

    for num, ax in enumerate(ax_list):
        # get ylims
        yls = ax.get_ylim()

        # get y limit dif
        ydif = max(yls) - min(yls)

        # get xlims
        xls = ax.get_xlim()
        xdif = max(xls) - min(xls)

        if inout == 'in':
            io = -1
        else:
            io = 1

        if ax.get_yscale() == 'log':
            yloc = yls[1] * (10 ** (io * 2.5 * yperc / 100))

        else:
            yloc = max(yls) + io * ydif * yperc / 100

        xloc = min(xls) + xdif * xperc / 100

        alpha_num = num + start_num
        if bold:
            ax.text(xloc, yloc, alph_list[alpha_num], color=color, fontsize=fs,
                    fontweight='bold')
        else:
            ax.text(xloc, yloc, alph_list[alpha_num], color=color, fontsize=fs)


def get_terminator(date_term, alt_km=300):
    """Get terminator lats and lons for plotting.

    Parameters
    ----------
    date_term : datetime
        Set the date and time for the terminator
    alt_km :float
        altitude of terminator
        default 300 km

    Returns
    -------
    lats : array-like
        latitudes of terminator
    lons : array-like
        longitudes of terminator
    """

    # Get antisolar position and the arc (terminator) at the given height
    antisolarpsn, arc, ang = pydarn.terminator(date_term, alt_km)

    # antisolarpsn contains the latitude and longitude of the antisolar point
    # arc represents the radius of the terminator arc

    # Now, you can directly use the geographic coordinates from antisolarpsn.
    lat_antisolar = antisolarpsn[1]
    lon_antisolar = antisolarpsn[0]

    # Get positions along the terminator arc in geographic coordinates
    lats = []
    lons = []

    for b in range(-180, 180, 1):  # Iterate over longitudes from -180 to 180
        lat, lon = pydarn.GeneralUtils.new_coordinate(
            lat_antisolar, lon_antisolar, arc, b, R=pydarn.Re)
        lats.append(lat)
        lons.append(lon)

    return lats, lons


def color_gradient(color, n, light_range=(0.25, 0.75)):
    """Generate n shades of a color.

    Parameters
    ----------
    color : str
        Any valid Matplotlib color.
    n : int
        Number of shades.
    light_range : tuple
        (darkest, lightest) lightness values in [0, 1].

    Returns
    -------
    color_list : list of str
        Hex color codes.
    """
    r, g, b = mcolors.to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)

    lightness = np.linspace(light_range[1], light_range[0], n)

    color_list = [
        mcolors.to_hex(colorsys.hls_to_rgb(h, l_new, s))
        for l_new in lightness
    ]

    return color_list
