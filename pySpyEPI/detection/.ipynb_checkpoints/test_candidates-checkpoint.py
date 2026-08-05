"""Barrel Rolling utility functions
Created 19 February 2026
alanahco
List Functions
--------------
test_peaks
no_doubbubs
"""
# Statistical Functions
import pandas as pd
import numpy as np
from pySpyEPI.utils import calc


def test_peaks(ne, peaks, properties, barrel_up, freq=1, nflag_mask=None,
               fft_test=False, nest_test=True, prom_test=True,
               percent_test=True, wdr_test=True, utest_mask=None,
               min_height=10**4, min_perc=10, wdr50_set=0.995, wdr20_set=0.95):
    """Flag the peaks.

    Parameters
    ----------
    ne : array
        array of densities
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    barrel_up : dataframe
        dataframe of barrel values including fitlered Ne and dne_dd_filt
    freq : int
        data frequency
        default 1 Hz
    nflag_mask : array-like or NoneType
        density flags mask
        where True values are good data
        False, bad data
        default None, no density flag assumed
    fft_test : bool
        default False
        fft of delta ne 20s for data will be calculated and erroneous
        data will be flagged (SWARM specific)
    nest_test : bool
        default True
        if True, test for nested bubbles
    prom_test : bool
        default True
        if True, test for a minimum prominence
    percent_test : bool
        if True, test for a minimum percent depth
    utest_mask : array-like or NoneType
        Additional test from user as a mask
        True values are assumed to be good and False values assumed to be for
        removal
        array must be same length as ne
    wdr_test : bool
        if True, test the weighted density ratio
    min_height: float
        minimum height for a prominence to be considered a prominence
        default 10**4
    min_perc : float
        minimum percent depth
        default 10%
    wdr50_set : float
        default 0.995
        minimum wdr for median of EPI WDR
    wdr20_set : float
        default 0.95
        if 25th percentile WDR of curve is below wdr20_set, use wdr20_set to
        test 20th percentile of EPI WDR
    Returns
    -------
    info_df : pd.DataFrame
        array of peak flags or DF of peak info that went into determining flags
        and flags

    Notes
    -----
    Peak flags and their meanings
        1 - peaks to keep
        2 - peak flagged because frequency is > 0.2 Hz at max amplitude
            or a > 2 Hz freq of amplitude > 0.5 * max(amplitude)
        3 - peak flagged because of bub in bub
        4 - peak flagged because prominence does not meet threshold
        5 - peak flagged because the peak and/ or edges are within +/- 0.5 mlat
            REMOVED (08 Oct 2025)
        6 - peak flagged from user specified test
        7 - peak flagged becuase during the halfwidth, Ne is flagged
        8 - peak flagged because the prominence is < min_height
        9 : peak flagged becuase of the WDR (weighted density ratio)
                l20 : 20th percentile of EPI left edge (edge 1) to center
                r20 : 20th percentile of EPI center to right edge (edge 2)
                med_l : median of left edge (edge 1) to center
                med_r : median of center to right edge (edge 2)
                lr_15 : 15th percentile of pass Np / W_B
                lr_set : 0.995 (0.5 %)
                Remove if any of the following are True
                (l20 > lr_15)
                (r20 > lr_15)
                (med_l > lr_set)
                (med_r > lr_set)
    The flags go in this order and will overwrite oen another:
        user specified (6)
        bubble in bubble (3)
        WDR (9)
        prom < 10**4 (8)
        Height Flag (4)
        frequency flag (2)
        flagged Np (7)

        * a negative flag means that > 30% of the flags are 2, so the other
        ones should not be used at all.

    """
    # Establish an array of 1s
    peak_flag = np.ones(len(peaks))

    # Left and Right indices
    xleft = np.trunc(properties["left_ips"]).astype(int)
    xright = np.ceil(properties["right_ips"]).astype(int)

    # Bubble in bubble-----------------------
    bub_in_bub = no_doubbubs(peaks, properties)

    if nest_test:
        peak_flag[bub_in_bub == 1] = 3

    # Calculate FFT to remove artificial signals ------------------------------
    if fft_test:
        data4fft = calc.calc_delta_ne(ne, t=20, freq=freq, divNe=True)
        spec, time_axis, freq_axis = calc.compute_dynamic_spectra(
            data4fft, f=1, win_size=20, step_size=10)

        # check if any spec above 0.2 Hz is greater than spec_check
        if len(spec) != 0:
            spec_check = np.max(spec) / 2
        else:
            spec_check = np.nan

        # get most prominent frequencies:
        high_freqs = []
        high_specs = []
        amp_list = []
        for hf in range(len(time_axis)):
            high_freqs.append(freq_axis[spec[hf, :].argmax()])
            amp_list.append(max(spec[hf, :]))
            if np.any(spec[hf, (freq_axis > 0.2)] > spec_check):
                high_specs.append(1)
            else:
                high_specs.append(0)

        high_freqs = np.array(high_freqs)
        high_specs = np.array(high_specs)
        amp_list = np.array(amp_list)

    # Keep track of following params:
    rlo_vec = []
    llo_vec = []
    rmed_vec = []
    hmed_vec = []
    lmed_vec = []
    hlo_vec = []
    depth_perc_vec = []
    max_freq_vec = []
    max_spec_vec = []
    flag_vec = []
    user_vec = []

    wdne = barrel_up['dne_wdr'].values
    lr_15 = np.nanpercentile(wdne, 15)
    lr_25 = np.nanpercentile(wdne, 25)
    lr_50 = np.nanpercentile(wdne, 50)

    # iterate through all identified peaks
    for i in range(len(peaks)):
        xr = xright[i]
        xl = xleft[i]
        xp = peaks[i]

        lval = (wdne[xl:xp + 1])
        rval = (wdne[xp:xr + 1])
        hval = (wdne[xl:xr + 1])
        med_l = (np.nanmedian(lval))
        med_r = (np.nanmedian(rval))
        med_h = (np.nanmedian(hval))

        l20 = np.nanpercentile(lval, 20)
        r20 = np.nanpercentile(rval, 20)
        h20 = np.nanpercentile(hval, 20)

        llo_vec.append(r20)
        rlo_vec.append(l20)
        hlo_vec.append(h20)

        lmed_vec.append(med_l)
        rmed_vec.append(med_r)
        hmed_vec.append(med_h)

        if wdr_test:

            if lr_15 < wdr20_set:
                t15 = wdr20_set
            else:
                t15 = lr_15

            if ((l20 > t15) | (r20 > t15) | (med_l > wdr50_set)
                    | (med_r > wdr50_set)):
                peak_flag[i] = 9

        if percent_test:
            # check the depth percentage
            ne_perc = (abs(barrel_up["dne_dd_filt"].iloc[xp]
                       / barrel_up["barrel_ne"].iloc[xp])) * 100

            # if it is less than minimum percent, remove
            if ne_perc < min_perc:
                peak_flag[i] = 4

        if prom_test:
            # if prom is < 10**4
            if barrel_up["dne_dd_filt"].iloc[xp] < min_height:
                peak_flag[i] = 8

        # save in vector
        depth_perc_vec.append(ne_perc)

        # If the barrel roll is < 10x the maximum of the whole,
        # set harder standards

        # For the FFT, the xp/xr/xl (indices) will align with time_axis since
        # both start at 0
        if fft_test:
            idx_after = -1
            idx_before = -1

            # Get indices of time axis closest to time_axis for xr and xl
            if np.any(time_axis > xr):
                idx_after = np.where(time_axis > xr)[0][0]
            if np.any(time_axis < xl):
                idx_before = np.where(time_axis < xl)[0][-1]

            # look for highest freq inside peak width
            # i believe idx_before or idx_after will always be defined
            if (idx_after == -1) & (idx_before != -1):
                max_freq = high_freqs[idx_before]
                max_spec = high_specs[idx_before]

            elif (idx_after != -1) & (idx_before == -1):
                max_freq = high_freqs[idx_after]
                max_spec = high_specs[idx_after]

            elif (idx_after != -1) & (idx_before != -1):
                # Look for freq associated with maximum amplitude
                amp_check = amp_list[idx_before:idx_after + 1]
                freq_check = high_freqs[idx_before:idx_after + 1]
                max_freq = freq_check[amp_check.argmax()]

                if np.any(high_specs[idx_before:idx_after + 1] > 0):
                    max_spec = 1
                else:
                    max_spec = 0

            if max_freq >= 0.2:
                peak_flag[i] = 2

            # save as vectors
            max_freq_vec.append(max_freq)
            max_spec_vec.append(max_spec)

        if nflag_mask is not None:
            # Density Flag-----------------------------------------------------
            # get ne flag for Halfwidth region representing bad Ne data
            H_flag = nflag_mask[xl:xr + 1]

            # if all values are True, don't discard EPI
            # if any are False, discard EPI
            if np.all(H_flag):
                flag_vec.append(0)
            else:
                peak_flag[i] = 7
                flag_vec.append(1)

        # user specified test -------------------------------------------------
        if utest_mask is not None:
            Huse = utest_mask[xl:xr + 1]

            # if all values are True, don't discard EPI
            # if any are False, discard EPI
            if np.all(Huse):
                user_vec.append(0)
            else:
                # update peak_flag if not already changed by a different flag
                if peak_flag[i] == 1:
                    peak_flag[i] = 6
                    user_vec.append(1)

    # if more than 50% of the flags are 2 (bad data),
    # set the peak flag of the remaining values to -1
    if fft_test:
        # if a dub bub has a bad frequency, remove the outer bub too #
        max_freq_vec = np.array(max_freq_vec)
        if np.any(max_freq_vec[bub_in_bub == 1] >= 0.2):

            # check outer bubble
            for pb in peaks[(bub_in_bub == 1) & (max_freq_vec >= 0.2)]:
                for xl, xr, xp in zip(xleft, xright, peaks):
                    if xp != pb:
                        if (xl < pb) & (xr > pb):
                            peak_flag[xl == xleft] = -1 * peak_flag[xl == xleft]
        if len(peak_flag) != 0:
            if len(peak_flag[peak_flag == 2]) / len(peak_flag) > 0.3:
                peak_flag[peak_flag != 2] = -1 * peak_flag[peak_flag != 2]

    info_df = pd.DataFrame()
    info_df['epi_flag'] = peak_flag
    info_df['p20_wdr_l'] = np.array(llo_vec)
    info_df['p20_wdr_r'] = np.array(rlo_vec)
    info_df['p20_wdr_h'] = np.array(hlo_vec)
    info_df['med_wdr_h'] = np.array(hmed_vec)
    info_df['med_wdr_r'] = np.array(rmed_vec)
    info_df['med_wdr_l'] = np.array(lmed_vec)
    info_df['percent_depth'] = np.array(depth_perc_vec)
    info_df['nested_bub'] = bub_in_bub

    info_df['pass_wdr_p15'] = [lr_15] * len(peak_flag)
    info_df['pass_wdr_p25'] = [lr_25] * len(peak_flag)
    info_df['pass_wdr_p25'] = [lr_50] * len(peak_flag)

    if nflag_mask is not None:
        info_df['flagged_np'] = np.array(flag_vec)

    if fft_test:
        info_df['max_freq_dnp20s'] = np.array(max_freq_vec)
        info_df['high_freq'] = np.array(max_spec_vec)

    if utest_mask is not None:
        info_df['user_test'] = np.array(user_vec)

    return info_df


def no_doubbubs(peaks, properties):
    """Return array of 0s and 1s where 1 indicates nested bubble.
    Parameters
    ----------
    peaks, properties : array-like
        list of peak indices and their properties
        see https://docs.scipy.org/doc/scipy/
        reference/generated/scipy.signal.find_peaks.html
        for detailed explanation
    Returns
    -------
    bub_in_bub: array-like
        array of same length as peaks
        if a peak is inside the width of another peak, then
        it is marked 1 as a bubble in bubble
    """
    # Get edge info
    xleft = np.trunc(properties["left_ips"]).astype(int)
    xright = np.ceil(properties["right_ips"]).astype(int)
    bub_in_bub = np.zeros(len(xleft))

    # Removing bubbles inside of bubbles
    # Go through bubble list
    for ip, xp in enumerate(peaks):
        lef = xleft[ip]
        rig = xright[ip]

        # get the sum of peaks that fall in this range
        dubs = peaks[(peaks > lef) & (peaks < rig)]
        dubl = xleft[(peaks > lef) & (peaks < rig)]
        dubr = xright[(peaks > lef) & (peaks < rig)]

        if len(dubs) > 1:
            sub = dubr - dubl
            # Keep maximum of sub and remove minimum
            minis = dubs[dubs != dubs[np.argmax(sub)]]

            # ger indices of all minimums
            for mm in minis:
                min_rem = np.where(peaks == mm)[0]

                # set bub_in_bub to 1
                bub_in_bub[min_rem] = 1

    # initiate array
    bub_in_bub = np.array(bub_in_bub)

    return bub_in_bub
