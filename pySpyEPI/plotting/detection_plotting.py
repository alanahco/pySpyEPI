"""Detection Plotting tools and plots
alanahco
09 July 2026
List Functions
--------------
flag_markers
algorithm_example
"""
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import math
from pySpyEPI.barrel import barrel_roll
from pySpyEPI.plotting import barrel_plot
import warnings
from pySpyEPI.plotting import plotL1


def flag_markers(ax, peaks, properties, stats_df, barrel_df,
                 col_pass='C1', col_fail='#32CD32', col_mlat_flag='#E8BCF0',
                 use_barrel=True, fail_gradient=True, light_range=None):
    """Plot the flags in different styles

    Parameters
    ----------
    ax : axis handle
        axis for plotting on
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/reference/generated/
            scipy.signal.find_peaks.html
        for detailed explanation
    stats_df : array-like
        df containing flagged peaks
        from cand info
    barrel_df : pd.DataFrame
        Created by pySpyEPI.barrel.barrel_roll.triple_barrel
    col_pass, col_fail, col_mlat_flag,: strings
        color codes for flags
        Defaults in order #028A0F, #32CD32, #B2FBA5
    Returns
    -------
    f_marks : list-like
        markers used for legend
    f_labs : list-like
        labels for legend
    f_cols : list-like
        colors used for legend
    f_mode : list-like
        'scatter' for same length as f_marks for legend

    Notes
    -----
    green circle: keep
    red upside down triangle: > 30% other flags were frequency flags
    red triangle : delta Ne 20s frequency peak >= 0.2 Hz
    red star : bubble in bubble
    red plus sign : percent depth < 10%
    red x : quality flag
    red square : prominence < 10**4
    red diamond: delta Ne 40s flag

    """
    f_marks = []
    f_labs = []
    f_cols = []
    f_mode = []
    f_flag = []

    # if color fail should be a gradient
    if fail_gradient:
        if light_range is None:
            light_range = [0.05, 0.78]
        fail_colors = plotL1.color_gradient(col_fail, 8,
                                            light_range=light_range)
    else:
        fail_colors = [col_fail] * 9

    stats_df = stats_df.copy()
    # replace all flags less than 0 to one negative number so that
    # only one flag marker is made for f_marks for legend
    epi_fl = stats_df['epi_flag'].values
    mlat_fl = stats_df['mlat_flag'].values

    ml_bad = (mlat_fl >= 4) & (mlat_fl != 5.1) & (mlat_fl != 4.1)
    epi_fl[(epi_fl == 1) & (ml_bad)] = -20

    for mk in np.unique(epi_fl):

        # base the marker off of the flag
        if mk == -20:
            marker = "o"
            col = col_mlat_flag
            flab = "mlat flag"
        elif mk == 1:
            marker = "o"
            col = col_pass
            flab = "EPI"
        elif mk == 2:
            marker = "^"
            col = fail_colors[0]
            flab = "Freq"
        elif mk == 3:
            marker = "*"
            col = fail_colors[1]
            flab = "Nested"
        elif mk == 4:
            marker = "P"
            col = fail_colors[2]
            flab = '%$_{depth}$'
        elif mk == 7:
            marker = "x"
            col = fail_colors[3]
            flab = "Quality"
        elif mk == 8:
            marker = "s"
            col = fail_colors[4]
            flab = "Prominence"
        elif mk == 9:
            marker = "D"
            col = fail_colors[5]
            flab = '$dN_{WDR}$'
        elif (mk < 0) & (mk > -20):
            marker = "v"
            col = fail_colors[6]
            flab = "% Freq"
        elif mk == 6:
            marker = "_"
            col = fail_colors[7]
            flab = "user"

        mask_mark = (epi_fl == mk)
        peaks_mark = peaks[mask_mark]
        prop_mark = (
            {key: value[mask_mark] for key, value in properties.items()
             if isinstance(value, np.ndarray)})

        if ax is not None:
            barrel_plot.barrel_ax(barrel_df, peaks_mark, prop_mark, ax,
                                  detrended=True, marker=marker,
                                  mark_color=col)

        # get flag makers and labels
        f_marks.append(marker)
        f_labs.append(flab)
        f_cols.append(col)
        f_mode.append('scatter')
        f_flag.append(mk)

    f_marks = np.array(f_marks)
    f_labs = np.array(f_labs)
    f_cols = np.array(f_cols)
    f_mode = np.array(f_mode)
    f_flag = np.array(f_flag)

    return f_marks, f_labs, f_cols, f_mode, f_flag


def algorithm_example(data_df, barrel_df, peaks, properties, stats_df, sat,
                      barrel_radii=None, envelope_radius=8, e_up=0.1, e_lo=0.05,
                      leg_loc='upper left', leg_anchor=None,
                      x_type='time', x_str='time', t_str='time', n_str='ne',
                      scale_n=True, exp_inc=True,
                      ex_scale=2, day_include=True, xlab='Time (s)',
                      wdr_ylim=None, n_unit='$cm^{-3}$', fs=15):

    """Plot an example of the algorithm for EPI detection.

    Parameters
    ----------
    data_df : pd.DataFrame
        Data from in-situ satellite, used for detection
    barrel_df : pd.DataFrame
        DataFrame containing info about barrel parameters
        Created by pySpyEPI.barrel.barrel_roll.triple_barrel
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/reference/generated/
            scipy.signal.find_peaks.html
        for detailed explanation
    stats_df : pd.DataFrame
        DataFrame containing EPI candidate info
        Created by pySpyEPI.detection.stats.cand_info.epi_stats
    sat : string
        Satellite name
        e.g. 'Swarm A'
    barrel_radii : list-like or NoneType
        list of 2 radii for barrel example
        default : None, use [8, 80]
    envelope : float
        barrel radius from barrel_radii for envelope purposes
        default 8
    e_up, e_lo : float
        upper and lower scaling factor for envelope
        default : 0.1, 0.05
    leg_loc : string
        location of legend
        default : 'upper_left'
    leg_anchor : list-like
        anchor for legend toa adjust placement
        default : NoneType, use [0.045, 1]
    x_type : string
        'time' (default) or other
        type of barrel radius
    x_str : string
        column name of x_type
        default : 'time'
    t_str, n_str : string
        column names of time and density
        default : 'time', 'ne'
    scale_n : boolean
        if True (default), density will be scaled
            by x_span / (np.nanmax(ne) * 10**exp)
    exp_inc : boolean
        if True (default), ne scaled by x_span / (np.nanmax(ne) * 10**exp)
            where exp = int(math.log10(x_span)) - ex_scale
        if False, ne scaled by x_span / (np.nanmax(ne) * 10**0)
    ex_scale : float
        scaling param for the exponent.
        exp = int(math.log10(x_span)) - ex_scale
    day_include : boolean
        if True (default), day will be included in title
        if False, only times will be included in title
    xlab :string
        x label default 'Time (s)',
    wdr_ylim : list-like or NoneType
        y axis limit for dN_WDR panel
        default : None, no limit applied
    n_unit : string
        unit of density
        default $cm^{-3}$
    fs : float
        fontsize
        default : 15

    Return
    ------
    fig : figure handle
        panel 1 : scaled density with barrel examples
        panel 2 : density barrel trends
        panel 3 : dN_WDR
        panel 4 : dN_LSDD with EPI candidates
    ax_list : list-like
        list of axes for alphabet labeling
    """
    # Create Figure
    fig = plt.figure(figsize=(12, 12))
    gs = gridspec.GridSpec(nrows=4, ncols=1)
    plt.rcParams.update({'font.size': fs})
    ax_list = []

    # Establish NoneTypes
    if barrel_radii is None:
        barrel_radii = [80, 8]
    if leg_anchor is None:
        leg_anchor = [0.045, 1]
    if envelope_radius not in barrel_radii:
        warnings.warn("envelope_radius not in barrel_radii")

    # scale x and y axes ------------------------------------------------------
    if x_type.lower() == 'time':
        x_secdf = (data_df[x_str].dt.hour * 3600
                   + data_df[x_str].dt.minute * 60
                   + data_df[x_str].dt.second
                   + data_df[x_str].dt.microsecond / 10 ** 6)
        x_sec = x_secdf.values
        if x_sec[0] > x_sec[-1]:
            x_sc = np.linspace(0, len(x_sec) - 1, len(x_sec))
        else:
            x_sc = x_sec - x_sec.min()
    else:
        # if not time, don't scale
        x_sc = data_df[x_str].values

    # get span of x
    x_span = np.nanmax(x_sc) - np.nanmin(x_sc)

    # Scale the density ---
    ne = data_df[n_str].values
    if scale_n:
        x_span = np.nanmax(x_sc) - np.nanmin(x_sc)
        if exp_inc:
            exp = int(math.log10(x_span)) - ex_scale
        else:
            exp = 0
        scal_param = x_span / (np.nanmax(ne) * 10**exp)
        ne_sc = ne * scal_param
    else:
        ne_sc = ne
    xlim = [x_sc.min(), x_sc.max()]

    # Scaled Density Panel ----------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])

    # plot scaled density
    ax.plot(x_sc, ne_sc, color="#F48FB1", linestyle='-', label='$N_{sc}$')

    rcols = ['#702670', '#5a88d5']
    rstyles = ['--', '-.']

    # plot barrel radii
    for r_sc, rcol, rstyle in zip(barrel_radii, rcols, rstyles):
        contact_x, contact_y = barrel_roll.simple_barrel(x_sc, ne_sc, r_sc,
                                                         direction='forward')
        blab = fr'$B_{{{r_sc}s}}$'
        ax.plot(contact_x, contact_y, linestyle=rstyle, linewidth=2,
                color=rcol, label=blab)

        if r_sc == envelope_radius:
            elab = fr'$E_{{{r_sc}s}}$'
            BRC_upper = contact_y + e_up * max(contact_y)
            BRC_lower = contact_y - e_lo * min(contact_y)

            ax.fill_between(contact_x, BRC_lower, BRC_upper, color='gray',
                            alpha=0.2, label=elab)
    # adjust axes
    ax.set_xlim(xlim)
    ax.xaxis.grid(True)

    # adjust ticks
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('gray')
    ax.spines['left'].set_color('gray')
    ax.tick_params(axis='both', which='major', labelsize=fs, length=10,
                   color='gray')
    ax.tick_params(bottom=True, labelbottom=False, top=False, labeltop=False)
    ax_list.append(ax)

    # add y label
    if scale_n:
        ax.set_ylabel('Scaled Density', fontsize=fs)
    else:
        ax.set_ylabel(f'Density [{n_unit}]', fontsize=fs)

    # legend
    ax.legend(loc=leg_loc, bbox_to_anchor=(leg_anchor))

    # Trend Panel -------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])

    # Large Scale Trend
    ax.plot(x_sc, barrel_df['barrel_ne'], color='#702670', linestyle='--',
            label='$N_{LST}$', linewidth=2)

    # narrow barrels
    ax.plot(x_sc, barrel_df['narlo_barrel'], color='#C64B8C', linestyle='-',
            alpha=0.5)
    ax.plot(x_sc, barrel_df['narup_barrel'], color='#C64B8C', linestyle='-',
            alpha=0.5)
    ax.fill_between(x_sc, barrel_df['narlo_barrel'], barrel_df['narup_barrel'],
                    color='#C64B8C', alpha=0.2, label=r'$B_{a}$ & $B_{b}$')

    # weighted barrel
    ax.plot(x_sc, barrel_df['weight_barrel'], color='#5a88d5',
            linestyle='-.', zorder=2, label='$N_{WB}$', linewidth=2)

    # adjust axes
    ax.set_xlim(xlim)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.xaxis.grid(True)

    # adjust ticks
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('gray')
    ax.spines['left'].set_color('gray')
    ax.tick_params(axis='both', which='major', labelsize=fs, length=10,
                   color='gray')
    ax.tick_params(bottom=True, labelbottom=False, top=False, labeltop=False)
    ax.set_ylabel('Density [$cm^{-3}$]', fontsize=fs)
    ax_list.append(ax)

    # legend
    ax.legend(loc=leg_loc, bbox_to_anchor=(leg_anchor))

    # WDR Panel ---------------------------------------------------------------
    ax = fig.add_subplot(gs[2, 0])

    # wide upper barrel
    ax.plot(x_sc, barrel_df['dne_wdr'], color='#5a88d5', linestyle='-')

    # calculate 15th percentile of whole dne_wdr
    wdr_15 = np.nanpercentile(barrel_df['dne_wdr'], 15)
    ax.axhline(y=0.995, color='#3C4142', linestyle=':', linewidth=2, alpha=1,
               label='0.995')
    ax.axhline(y=0.95, color='#7D7F7C', linestyle='-.', linewidth=2, alpha=1,
               label='0.95')
    ax.axhline(y=wdr_15, color='#B6B6B4', linestyle='--', linewidth=2, alpha=1,
               label='$15^{th}$ %-ile')

    # adjust axes
    ax.set_xlim(xlim)
    ax.set_ylim(wdr_ylim)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.xaxis.grid(True)

    # adjust ticks
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('gray')
    ax.spines['left'].set_color('gray')
    ax.tick_params(axis='both', which='major', labelsize=fs, length=10,
                   color='gray')
    ax.tick_params(bottom=True, labelbottom=False, top=False, labeltop=False)
    ax.set_ylabel('$dN_{WDR}$', fontsize=fs)
    ax_list.append(ax)

    # Legend
    ax.legend(loc=leg_loc, bbox_to_anchor=(leg_anchor))

    # LSDD Panel --------------------------------------------------------------
    ax = fig.add_subplot(gs[3, 0])

    # wide upper barrel
    ax.plot(x_sc, barrel_df['dne_dd_filt'], color='#702670', linestyle='-')

    keep_mask = stats_df['epi_flag'] == 1

    ax.scatter(x_sc[peaks[keep_mask]],
               barrel_df['dne_dd_filt'].iloc[peaks[keep_mask]],
               marker='o', color='C1', zorder=2, s=35, label='Confirmed')

    f_marks, f_labs, f_cols, f_mode, f_flag = flag_markers(
        None, peaks, properties, stats_df, barrel_df, col_pass='C1',
        col_fail='#32CD32', col_mlat_flag='#B2FBA5')

    epi_fl = stats_df['epi_flag'].values
    mlat_fl = stats_df['mlat_flag'].values

    ml_bad = (mlat_fl >= 4) & (mlat_fl != 5.1) & (mlat_fl != 4.1)
    epi_fl[(epi_fl == 1) & (ml_bad)] = -20

    for ef in np.unique(epi_fl):
        if ef != 1:
            f_mask = (ef == epi_fl)
            e_mask = (ef == f_flag)
            ax.scatter(
                x_sc[peaks[f_mask]],
                barrel_df['dne_dd_filt'].iloc[peaks[f_mask]],
                marker=f_marks[e_mask][0], color=f_cols[e_mask][0], zorder=2,
                s=55, label=f_labs[e_mask][0])

    # plot prominences
    ymin = (barrel_df['dne_dd_filt'].iloc[peaks]
            - properties['prominences'])
    ymax = barrel_df['dne_dd_filt'].iloc[peaks]
    ax.vlines(x=x_sc[peaks], ymin=ymin, ymax=ymax,
              color="C0", zorder=0)
    xl = np.ceil(properties["left_ips"]).astype(int)
    xr = np.trunc(properties["right_ips"]).astype(int)
    ax.hlines(y=properties["width_heights"],
              xmin=x_sc[xl],
              xmax=x_sc[xr], color="C0", zorder=0)

    # Legend
    ax.legend(loc=leg_loc, bbox_to_anchor=(leg_anchor))

    lsdd_string = '$dN_{LSDD}$'
    # adjust axes
    ax.set_xlim(xlim)
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
    ax.xaxis.grid(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('gray')
    ax.spines['left'].set_color('gray')
    ax.tick_params(axis='both', which='major', labelsize=fs, length=10,
                   color='gray')
    ax.tick_params(top=False, labeltop=False)
    ax.set_ylabel(f'{lsdd_string} [{n_unit}]', fontsize=fs)
    ax_list.append(ax)
    ax.set_xlabel(xlab)

    # Set super title ---------------------------------------------------------
    f_time_1 = data_df[t_str].iloc[0].strftime('%H:%M')
    f_time_2 = data_df[t_str].iloc[-1].strftime('%H:%M')
    title_date = data_df[t_str].iloc[0].strftime('%d %b %Y')

    dash = '$-$'
    if day_include:
        titl = f'Swarm {sat} {title_date} {f_time_1} UT {dash} {f_time_2}'
    else:
        titl = f'{f_time_1} UT {dash} {f_time_2}'
    plt.suptitle(f'{titl} UT', y=0.91, x=0.5, fontsize=fs)

    return fig, ax_list


def detection_plot(data_df, x_str, barrel_df, stats_df, peaks, properties,
                   flag_mark=True, xtype='datetime', leg_ax=None, xlim=None,
                   xlab=None, n_type='p', n_unit='$cm^{-3}$', lon_lim=None,
                   time_str='time', lon_str='lon', lt_hr_str='lt_hr',
                   mlat_str='mlat', lt_str='lt', lat_str='lat', shade_epi=True,
                   shade_col="#A7A6BA", shade_alpha=0.2, sat_name=None, fs=13):
    """Create figure showing detection parameters.

    Parameters
    ----------
    data_df : pd.DataFrame
        data to plot
    x_str : string
        column name for x axis string
    barrel_df : pd.DataFrame
        DataFrame containing info about barrel parameters
        Created by pySpyEPI.barrel.barrel_roll.triple_barrel
    stats_df : pd.DataFrame
        DataFrame containing EPI candidate info
        Created by pySpyEPI.detection.stats.cand_info.epi_stats
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/reference/generated/
            scipy.signal.find_peaks.html
        for detailed explanation
    flag_mark : boolean
        if True (default), discarded candidates are marked by removal reason
    xtype : string
        type of data for x_str.
        Default 'datetime', x axis will be in format H:M
        if longitude or latitude, degrees direction added to labels
    leg_ax : axis_handle or NoneType
        if given, legend will be put on axis
        if None (default), legend will be put on ax in best location
    xlim : NoneType or array-like
        x axis limit
        if None (default), data_df[x_str] first and last values used to limit
    xlab : string
        x axis label
        if None, default is x_str
    n_type : string
        type of density
        will appear on figure as $N_{n_type}$
        to save space recommend 'e' for electron, 'i' for ion,
        and 'p' for plasma
        'p' is default
    n_unit : string
        unit for density
        default cm^-3
    time_str, lon_str, lt_hr_str, mlat_str, lt_str, lat_str : strings
        variable column names for time, longitude, local time hour,
        magnetic latitude, local time, latitutde
    shade_epi : boolean
        if True (default), shading will be added to confirmed EPI halfwidths
        if False, no shading will be added
    shade_col : string
        color of shading for EPI halfwidths
        default "#A7A6BA" (purplish gray)
    shade_alpha : float
        shade transparency
        default 0.2
    sat_name : string or NoneType
        satellite name for title
        default None
    fs : float
        fontsize
        default 13

    Returns
    -------
    fig : figure handle
        figure containing coordinate axis, dN_WDR, Np, and dN_LSDD
    ax_list : list of axes handles
        list of axes with data plotted (not legend axes)
    """

    # establish figure
    fig = plt.figure(figsize=(12, 12))
    ax_list = []
    gs = gridspec.GridSpec(nrows=4, ncols=2, wspace=0.1, hspace=0.4,
                           width_ratios=[1, 0.05])
    n = 0
    xval = data_df[x_str].values

    if xlab is None:
        xlab = x_str

    if xtype.lower() == 'datetime':
        time_bool = True

    # -------------- Coordinate Panel -----------------------------------------
    ax = fig.add_subplot(gs[n, 0])
    ax_list.append(ax)
    cax = fig.add_subplot(gs[n, 1])
    plotL1.plot_maglat_panel(
        data_df, ax=ax, cax=cax, fs=fs, glat_bool=True,
        time_str=time_str, lon_str=lon_str, lt_hr_str=lt_hr_str,
        mlat_str=mlat_str, lt_str=lt_str, lat_str=lat_str, lat_lines=True)
    n += 1
    ax.set_ylim(lon_lim)

    # -------------- WDR Panel ------------------------------------------------
    ax = fig.add_subplot(gs[n, 0])
    ax_list.append(ax)
    ax.plot(data_df[x_str], barrel_df['dne_wdr'], color='#5a88d5',
            zorder=1)

    plotL1.adjust_axes_time(
        xval, ax, xval, xlab=xlab, ylab='$dN_{WDR}$',
        time_bool=time_bool, leg=False, fs=fs, x_on=True, sci_y=False,
        minor_ticks=True, xgrid=True, maj_num=11)

    if xtype.lower() == 'latitude':
        plotL1.format_latitude_labels(ax, xy='x', dec=1)
    elif xtype.lower() == 'longitude':
        plotL1.format_longitude_labels(ax, xy='x', dec=1)
    n += 1

    # -------------- Density Panel --------------------------------------------
    ax = fig.add_subplot(gs[n, 0])
    ax_list.append(ax)
    lax = fig.add_subplot(gs[n, 1])
    ax.plot(data_df[x_str], barrel_df['ne'], color='#F48FB1')
    ax.plot(data_df[x_str], barrel_df['barrel_ne'], color='#702670',
            linestyle='--')
    ax.plot(data_df[x_str], barrel_df['weight_barrel'], color='#5a88d5',
            linestyle='-.')

    plotL1.adjust_axes_time(
        xval, ax, xval, xlab=xlab, ylab=f'Density {n_unit}',
        time_bool=time_bool, leg=False, fs=fs, x_on=True, sci_y=True,
        minor_ticks=True, xgrid=True)

    if xtype.lower() == 'latitude':
        plotL1.format_latitude_labels(ax, xy='x', dec=1)
    elif xtype.lower() == 'longitude':
        plotL1.format_longitude_labels(ax, xy='x', dec=1)

    leg_labs = [fr'$N_{{{n_type}}}$', '$N_{LST}$', '$N_{WB}$']
    leg_cols = ['#F48FB1', '#702670', '#5a88d5']
    modes = ['line', 'line', 'line']
    leg_styles = ['-', '--', '-.']
    lwidths = [2] * 3

    plotL1.make_legend(lax, leg_labs=leg_labs, leg_cols=leg_cols,
                       leg_styles=leg_styles, modes=modes, lwidths=lwidths,
                       loc='center left', fontsize=fs)
    n += 1

    # -------------- LSDD Panel -----------------------------------------------
    ax = fig.add_subplot(gs[n, 0])
    lax = fig.add_subplot(gs[n, 1])
    ax_list.append(ax)
    lsdd_ax(ax, data_df, x_str, barrel_df, stats_df, peaks, properties,
            flag_mark=flag_mark, xtype=xtype, leg_ax=lax, xlim=None,
            xlab=xlab, leg_loc='center_left')

    # add shading
    if shade_epi:
        keep_mask = stats_df['epi_flag'] == 1
        prop_kp = ({key: value[keep_mask] for key, value in properties.items()
                    if isinstance(value, np.ndarray)})
        plotL1.add_peak_shading(ax_list, data_df[x_str].values, prop_kp,
                                col=shade_col, alpha=shade_alpha, axis='x')

    title_date = data_df[time_str].iloc[0].strftime('%d %b %Y')
    ax.set_title(title_date, y=-0.45, x=0.95, fontsize=fs)

    # Set super title ---------------------------------------------------------
    if sat_name:
        plt.suptitle(f'{sat_name}', y=0.94, x=0.48, fontsize=fs)

    return fig, ax_list


def lsdd_ax(ax, data_df, x_str, barrel_df, stats_df, peaks, properties,
            flag_mark=True, xtype='datetime', leg_ax=None, xlim=None,
            xlab=None, n_unit='$cm^{-3}$', leg_loc='center left', fs=13):
    """Plot dN_LSDD panel.

    Parameters
    ----------
    ax : axis handle
        axis for plotting
    data_df : pd.DataFrame
        data to plot
    x_str : string
        column name for x axis string
    barrel_df : pd.DataFrame
        DataFrame containing info about barrel parameters
        Created by pySpyEPI.barrel.barrel_roll.triple_barrel
    stats_df : pd.DataFrame
        DataFrame containing EPI candidate info
        Created by pySpyEPI.detection.stats.cand_info.epi_stats
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/reference/generated/
            scipy.signal.find_peaks.html
        for detailed explanation
    flag_mark : boolean
        if True (default), discarded candidates are marked by removal reason
    xtype : string
        type of data for x_str.
        Default 'datetime', x axis will be in format H:M
        if longitude or latitude, degrees direction added to labels
    leg_ax : axis_handle or NoneType
        if given, legend will be put on axis
        if None (default), legend will be put on ax in best location
    xlim : NoneType or array-like
        x axis limit
        if None (default), data_df[x_str] first and last values used to limit
    xlab : string
        x axis label
        if None, default is x_str
    n_unit : string
        unit for y axis
        default cm^-3
    fs : float
        fontsize
        default 13

    """

    ax.plot(data_df[x_str], barrel_df['dne_dd_filt'], color='#702670',
            zorder=1)

    keep_mask = stats_df['epi_flag'] == 1
    disc_mask = stats_df['epi_flag'] != 1

    if len(peaks[keep_mask]) > 0:
        leg_labs = ['EPI']
        leg_cols = ['C1']
        leg_marks = ['o']
    else:
        leg_labs = []
        leg_cols = []
        leg_marks = []

    # plot EPIs
    ax.scatter(data_df[x_str].iloc[peaks[keep_mask]],
               barrel_df['dne_dd_filt'].iloc[peaks[keep_mask]],
               marker='o', color='C1', zorder=2, s=35, label='EPI')

    # add flag markers if specified
    if flag_mark:

        # get colors and markers
        f_marks, f_labs, f_cols, f_mode, f_flag = flag_markers(
            None, peaks, properties, stats_df, barrel_df, col_pass='C1',
            col_fail='#32CD32', col_mlat_flag='#E8BCF0')

        epi_fl = stats_df['epi_flag'].values
        mlat_fl = stats_df['mlat_flag'].values

        ml_bad = (mlat_fl >= 4) & (mlat_fl != 5.1) & (mlat_fl != 4.1)
        epi_fl[(epi_fl == 1) & (ml_bad)] = -20

        # plot discarded EPIs
        for ef in np.unique(epi_fl):
            if ef != 1:
                f_mask = (ef == epi_fl)
                e_mask = (ef == f_flag)
                ax.scatter(
                    data_df[x_str].iloc[peaks[f_mask]],
                    barrel_df['dne_dd_filt'].iloc[peaks[f_mask]],
                    marker=f_marks[e_mask][0], color=f_cols[e_mask][0],
                    zorder=2, s=55, label=f_labs[e_mask][0])

                # add flags to list for legend
                if f_labs[e_mask][0] not in leg_labs:
                    leg_labs.append(f_labs[e_mask][0])
                    leg_cols.append(f_cols[e_mask][0])
                    leg_marks.append(f_marks[e_mask][0])
    else:
        # if flag_mark is False, only plot the discarded as X's
        ax.scatter(
            data_df[x_str].iloc[peaks[disc_mask]],
            barrel_df['dne_dd_filt'].iloc[peaks[disc_mask]], marker='X',
            color='#32CD32', zorder=2, s=35, label='Discarded')

        leg_labs.append('Discarded')
        leg_cols.append('#32CD32')
        leg_marks.append('X')

    # plot prominences
    ymin = (barrel_df['dne_dd_filt'].iloc[peaks]
            - properties['prominences'])
    ymax = barrel_df['dne_dd_filt'].iloc[peaks]
    ax.vlines(x=data_df[x_str].iloc[peaks], ymin=ymin, ymax=ymax,
              color="C0", zorder=0)
    xl = np.ceil(properties["left_ips"]).astype(int)
    xr = np.trunc(properties["right_ips"]).astype(int)
    ax.hlines(y=properties["width_heights"],
              xmin=data_df[x_str].iloc[xl],
              xmax=data_df[x_str].iloc[xr], color="C0", zorder=0)

    # Legend
    if leg_ax:
        leg_modes = ['scatter'] * len(leg_labs)
        plotL1.make_legend(leg_ax, leg_labs=leg_labs, leg_cols=leg_cols,
                           leg_styles=leg_marks, modes=leg_modes,
                           loc=leg_loc, fontsize=fs)
    else:
        ax.legend()

    # set x limit
    if xlim:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim([data_df[x_str].iloc[0], data_df[x_str].iloc[-1]])

    if xlab is None:
        xlab = x_str

    # adjust axis by type
    if xtype.lower() == 'datetime':
        time_bool = True
    else:
        time_bool = False
    xval = data_df[x_str].values

    xval = data_df[x_str].values
    plotL1.adjust_axes_time(xval, ax, xval, xlab=xlab,
                            ylab='$dN_{LSDD}$ [$cm^{-3}$]',
                            time_bool=time_bool, leg=False,
                            fs=fs, x_on=True, sci_y=True, minor_ticks=True,
                            xgrid=True)
    if xtype.lower() == 'latitude':
        plotL1.format_latitude_labels(ax, xy='x', dec=1)
    elif xtype.lower() == 'longitude':
        plotL1.format_longitude_labels(ax, xy='x', dec=1)

    # adjust axes parameters
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.xaxis.grid(True)
    ax.spines['bottom'].set_color('gray')
    ax.spines['left'].set_color('gray')
    ax.tick_params(axis='both', which='major', labelsize=fs, length=10,
                   color='gray')
    ax.tick_params(top=False, labeltop=False)

    # add labels
    lsdd_string = '$dN_{LSDD}$'
    ax.set_ylabel(f'{lsdd_string} [{n_unit}]', fontsize=fs)
    ax.set_xlabel(xlab)
