"""Calculation Functions
Created 19 February 2026
alanahco

List Functions
--------------
detect_outliers
moving_median
moving_std
calc_delta_ne
bandpass_filter
highpass_filter
naninterp
compute_dynamic_spectra
safe_nanmedian
safe_nanmean
safe_nanstd
moving_average_same_length
is_dark
"""
import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt
from numpy.fft import fft, fftfreq
import pydarn


def detect_outliers(arr):
    """ Detect outliers in an array
    Parameters:
    -----------
    arr: numpy array object
        set of numbers
    Returns:
    --------
    outlier_indices: numpy array object
        array of indices where arr has outliers
    Notes:
    ------
    Uses InterQuartile Range (IQR)
    IQR = q3-q1
    outlier > q3 + 1.5*IQR
    outlier < q1 - 1.5*IQR
    """
    arr = np.array(arr)
    if len(arr[np.isfinite(arr)]) > 0:
        q1 = np.percentile(arr[np.isfinite(arr)], 25)
        q3 = np.percentile(arr[np.isfinite(arr)], 75)
        IQR = q3 - q1
        upper_lim = q3 + 1.5 * IQR
        lower_lim = q1 - 1.5 * IQR
        outlier_indices = np.where((arr > upper_lim) | (arr < lower_lim))[0]
    else:
        outlier_indices = []

    return outlier_indices


def moving_median(arr, window):
    """ Calculate moving median with NaN handling and same-length output
    Parameters
    ----------
    arr : array-like
        array for a moving mean to be calculated over
    window : int
        window size, in number of points
    Return
    ------
    med_arr : array of same length as arr
        median filtered array
    """
    med_arr = pd.Series(arr).rolling(window, center=True,
                                     min_periods=1).median().to_numpy()
    return med_arr


def moving_std(arr, window=3):
    """ calculate a moving standard deviation
    Parameters
    ----------
    arr : array-like
        data set for std to be calculated on
    window : int kwarg default 3
        number of data points to include in the standard deviation
    Returns
    -------
    std_arr : array-like with len of arr
        array of standard deviations
    """
    std_arr = pd.Series(arr).rolling(window, center=True,
                                     min_periods=1).std().to_numpy()
    return std_arr


def calc_delta_ne(ne, t=40, freq=1, divNe=True):
    """calculate delta Ne at set t as defined by
    Parameters
    ----------
    ne : array-like
        electron density
    t : kwarg int
        default 40 seconds, length of median filter
    difNe : boolean
        if True, returned parameter will be delta ne / ne (default)
        if False, returned parameter will be delta ne
    Returns
    -------
    delta_ne/ne : array-like
        median filtered Ne divided by Ne for a unitless measurement
        with the same length as original input ne array
    Notes
    -----
    The following definitions come from the product definitions for IPIR
    dataset, the dataset given uses the absolute value and does not divide by
    ne again.
    Delta_Ne10s = Derived by subtracting Ne by its median filtered value
        in 10 seconds; indicates the electron density fluctuations
        smaller than 75 km
    Delta_Ne20s = Derived by subtracting Ne by its median filtered value
        in 20 seconds; indicates the electron density fluctuations
        smaller than 150 km
    Delta_Ne40s = Derived by subtracting Ne by its median filtered value
        in 40 seconds; indicates the electron density fluctuations
        smaller than 300 km
    """
    # adjust window by freq. e.g. 20 seconds, 2 Hz, 40 point window
    t_move = t * freq
    med_filt = moving_median(ne, t_move)
    delta_ne = ne - med_filt

    if divNe:
        return delta_ne / ne
    else:
        return delta_ne


def bandpass_filter(data, fs, order=4, lowcut=None, highcut=None, type='band'):
    """Apply a Butterworth bandpass filter.

    Parameters:
        data : array-like
            The input signal to filter.
        fs : float
            Sampling frequency in Hz.
        lowcut : float
            Low frequency cutoff in Hz.
        highcut : float
            High frequency cutoff in Hz.
        order : int
            Filter order.
        type : str kwarg 'band' default
            determines what type of filter
            including 'band', 'high', and
            'low
    Returns:
        Filtered data (same shape as input).
    """
    nyq = 0.5 * fs  # Nyquist Frequency
    if type == 'band':
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype=type)  # Bandpass
    elif type == 'high':
        high = highcut / nyq
        b, a = butter(order, [high], btype=type)  # Highpass
    elif type == 'low':
        low = lowcut / nyq
        b, a = butter(order, [low], btype=type)  # Lowpass
    y = filtfilt(b, a, data)
    return y


def highpass_filter(data, cutoff, fs, order):
    """ High Pass filter for magnetic field data
    Parameters
    ----------
    data : array-like
        dataset to be filtered
    cutoff : float
        cutoff frequency
    fs : int
        data cadence
    order : int
        order for butterworth filter
    Returns
    -------
    filtered data
    """
    nyq = 0.5 * fs
    norm_cutoff = cutoff / nyq
    b, a = butter(order, norm_cutoff, btype='high')
    return filtfilt(b, a, data)


def naninterp(data):
    """
    Linearly interpolate over NaNs in 1D or 2D NumPy arrays.

    Parameters
    ----------
    data: np.ndarray, can be 1D or 2D. In 2D, rows = time, cols = components.

    Returns
    -------
    np.ndarray of same shape as input, with NaNs linearly interpolated.

    Notes
    -----
    raises value error if the input array is not 1D or 2D
    """
    data = np.asarray(data)
    if data.ndim == 1:
        x = np.arange(data.shape[0])
        y = data.copy()
        nans = np.isnan(y)
        if not np.all(nans):  # Avoid all-NaN case
            y[nans] = np.interp(x[nans], x[~nans], y[~nans])
        return y
    elif data.ndim == 2:
        data_interp = data.copy()
        for i in range(data.shape[1]):
            x = np.arange(data.shape[0])
            y = data[:, i]
            nans = np.isnan(y)
            if np.all(nans):
                continue
            y[nans] = np.interp(x[nans], x[~nans], y[~nans])
            data_interp[:, i] = y
        return data_interp
    else:
        raise ValueError("Input array must be 1D or 2D.")


def compute_dynamic_spectra(signal, f=50, win_size=1024, step_size=500):
    """
    Compute the dynamic spectrum (sliding FFT) of a signal with proper
    amplitude normalization.

    Parameters
    ----------
    signal : 1D numpy array
    fs : int
        sampling frequency in Hz
    win_size : int
        number of samples per FFT window (default 1024)
    step_size : int
        number of samples to shift the window by (default 500)

    Returns
    -------
    spec : 2D array
        2D array of FFT magnitudes in signal unit (time x frequency)
    time_axis : 1D array
        array of center times (in seconds)
    freq_axis : 1D array
        array of frequency bins (Hz)
    Notes
    -----
    default: a 50Hz signal, with a window size of 1024 (~20 seconds)
    by stepsize 500 (~10 seconds)
    """
    n = len(signal)
    freq_axis = fftfreq(win_size, d=1 / f)[:win_size // 2]

    spec = []
    time_axis = []

    for start in range(0, n - win_size + 1, step_size):
        end = start + win_size
        window = signal[start:end]
        # Normalize FFT magnitude
        fft_vals = (2.0 / win_size) * np.abs(fft(window))[:win_size // 2]
        spec.append(fft_vals)
        center_time = (start + end) / 2 / f  # time in seconds
        time_axis.append(center_time)

    spec = np.array(spec)  # shape: (n_windows, n_freq_bins)
    time_axis = np.array(time_axis)

    return spec, time_axis, freq_axis


def safe_nanmedian(arr):
    """ Get mean of array if it is not all nan or length equal to 0
    """
    if arr.size == 0 or np.all(np.isnan(arr)):
        return np.nan
    return np.nanmedian(arr)


def safe_nanmean(arr):
    """ Get mean of array if it is not all nan or length equal to 0
    """
    if arr.size == 0 or np.all(np.isnan(arr)):
        return np.nan
    return np.nanmean(arr)


def safe_nanstd(arr):
    """ Get mean of array if it is not all nan or length equal to 0
    """
    if arr.size == 0 or np.all(np.isnan(arr)):
        return np.nan
    return np.nanstd(arr)


def moving_average_same_length(arr, window):
    """
    Computes the moving average with output length same as input.

    Parameters
    ----------
    arr : array-like (1D)
    window : int
        number of points to average over

    Returns
    -------
    result : np.ndarray
        Moving-averaged array with same length as input.
    """
    result = pd.Series(arr).rolling(window, center=True,
                                    min_periods=1).mean().to_numpy()

    return result


def is_dark(lat, lon, time, height_km=300):
    """
    Returns:
        1 -> dark
        0 -> sunlit
    """

    # antisolar point
    antisolarpsn, arc, ang = pydarn.terminator(time, height_km)

    lat0 = np.radians(antisolarpsn[1])
    lon0 = np.radians(antisolarpsn[0])

    lat1 = np.radians(lat)
    lon1 = np.radians(lon)

    # great-circle angular distance
    cos_d = (
        np.sin(lat0) * np.sin(lat1)
        + np.cos(lat0) * np.cos(lat1) * np.cos(lon1 - lon0)
    )

    # numerical safety
    cos_d = np.clip(cos_d, -1, 1)

    d = np.degrees(np.arccos(cos_d))

    # inside antisolar hemisphere = darkness
    return int(d < 90)
