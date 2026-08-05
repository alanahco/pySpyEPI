"""Flag detections.
Created 19 February 2026
alanahco
List Functions
--------------
mlat_test
"""
import numpy as np


def mlat_flag(p_mlat, flags):
    """Test the magnetic latitude of the detections and flag accordingly.

    Parameters
    ----------
    p_mlat : array-like
        magnetic latitude of the peaks
    flags : array-like
        flags of each peak from process_peaks
    Returns
    -------
    mflag : array-like
        magnetic latitude flags desribed in notes
    Notes
    -----
    Depletions above 40 degrees need to be treated skeptically
    Additionally if the only depletions are "super", be skeptical
    0 : no flag
    5, 10, 15 : abs(peak mlat) >= 40 degrees
        5 : at least one 1 flag between +/- 25
        10 : no flags between +/- 25
        15 : flags present but no 1 flag between +/- 25
        Value can be X.0, X.1, or X.2
        X.0 means no flags between 25 and 40
        X.1 means a 1 flag between 25 -> 35 AND 35 -> 40
        X.2 means a 1 flag between 25 and 35 but not 35 to 40
        X.3 means a 1 flag between 35 and 40 but not 25 to 35
        X.4 means a non 1 flag between 25 and 40
    4, 8, 12 : abs(peak mlat) < 40 degrees and >= 35
        4 : at least one 1 flag between +/- 25
        8 : no flags between +/- 25
        12 : flags present but no 1 flag between +/- 25
        Value can be X.0, X.1, or X.2
        X.0 means no flags between 25 and 35
        X.1 means a 1 flag between 25 and 35
        X.2 means a non 1 flag between 25 and 35
    3, 6, 9 : abs(peak mlat) < 35 degrees and > 25
        3 : at least one 1 flag between +/- 25
        6 : no flags between +/- 25
        9 : flags present but no 1 flag between +/- 25
    Generally best to keep flags 0, 3, 4.1, and maybe 5.1
    """
    # initialize magnetic latitude flag
    mflag = []

    # north check
    # get flags between 0 and 25
    # if there are no flags, ch_n = 0
    # if there are flags and one of them is 1, ch_n = 1
    # if there are flags and none of them are 1, ch_n = 2
    f_north = flags[(p_mlat <= 25) & (p_mlat >= 0)]
    if len(f_north) == 0:
        ch_n = 0
    else:
        if np.any(f_north == 1):
            ch_n = 1
        else:
            ch_n = 2

    # check midlats
    f_north_mids = flags[(p_mlat > 25) & (p_mlat < 35)]

    if len(f_north_mids) == 0:
        ch_nm = 0
    else:
        if np.any(f_north_mids == 1):
            ch_nm = 1
        else:
            ch_nm = 2

    # check high midlats
    f_north_highs = flags[(p_mlat >= 35) & (p_mlat < 40)]

    if len(f_north_highs) == 0:
        ch_nh = 0
    else:
        if np.any(f_north_highs == 1):
            ch_nh = 1
        else:
            ch_nh = 2

    # south check
    # get flags between -25 and 0
    # if there are no flags, ch_n = 0
    # if there are flags and one of them is 1, ch_n = 1
    # if there are flags and none of them are 1, ch_n = 2
    f_south = flags[(p_mlat >= -25) & (p_mlat < 0)]
    if len(f_south) == 0:
        ch_s = 0
    else:
        if np.any(f_south == 1):
            ch_s = 1
        else:
            ch_s = 2

    # check midlats
    f_south_mids = flags[(p_mlat > -35) & (p_mlat < -25)]

    if len(f_south_mids) == 0:
        ch_sm = 0
    else:
        if np.any(f_south_mids == 1):
            ch_sm = 1
        else:
            ch_sm = 2

    # check high midlats
    f_south_highs = flags[(p_mlat >= 35) & (p_mlat < 40)]

    if len(f_south_highs) == 0:
        ch_sh = 0
    else:
        if np.any(f_south_highs == 1):
            ch_sh = 1
        else:
            ch_sh = 2

    # mflag if above 40 degrees (abs val)
    # mflag if above 25 degrees and below 40 (abs val)
    for i in range(len(flags)):

        # get mlat for peak
        m_check = p_mlat[i]

        # check if it is south or north and set ch accordingly
        if m_check < 0:
            ch = ch_s
            ch_m = ch_sm
            ch_h = ch_sh
        else:
            ch = ch_n
            ch_m = ch_nm
            ch_h = ch_nh

        # check 35 to 40 AND 25 to 35 for 40 up
        if (ch_m == 1) & (ch_h == 1):
            ch_add = 0.1
        elif (ch_m == 1) & (ch_h != 1):
            ch_add = 0.2
        elif (ch_m == 0) & (ch_h == 0):
            ch_add = 0.0
        elif (ch_m != 1) & (ch_h == 1):
            ch_add = 0.3
        elif (ch_m == 2) | (ch_h == 2):
            ch_add = 0.3

        # adjust ch midlat for 35 to 40 mlat
        ch_m_add = ch_m / 10

        # flag based on ch and lat location
        if abs(m_check) >= 40:
            if ch == 0:
                mflag.append(10 + ch_add)
            elif ch == 1:
                mflag.append(5 + ch_add)
            elif ch == 2:
                mflag.append(15 + ch_add)
        elif (abs(m_check) < 40) & (abs(m_check) >= 35):
            if ch == 0:
                mflag.append(8 + ch_m_add)
            elif ch == 1:
                mflag.append(4 + ch_m_add)
            elif ch == 2:
                mflag.append(12 + ch_m_add)
        elif (abs(m_check) < 35) & (abs(m_check) > 25):
            if ch == 0:
                mflag.append(6)
            elif ch == 1:
                mflag.append(3)
            elif ch == 2:
                mflag.append(9)
        else:
            mflag.append(0)

    # convert to array
    mflag = np.array(mflag)
    return mflag
