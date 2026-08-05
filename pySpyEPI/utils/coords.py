"""Coordinate Transformations.
Created 19 February 2026
alanahco

List Functions
--------------
longitude_to_local_time
compute_qd_coords
earth_radius
compute_magnetic_coords
NEC_2_MFA
track_to_nec_rotation
coordinate_transform
sat2nec_errors
nec2mfa_errors
compute_declination_inclination
NEC_2_MFA_Rotation
ecef2ned
geodetic_to_ecef
"""
import pandas as pd
import numpy as np
import apexpy
from scipy.signal import butter, filtfilt
from pySpyEPI.utils import calc


def longitude_to_local_time(longitude, utc_time):
    """Convert Longiutde to local time.

    Parameters
    ----------
    longitude : array-like
        longitudes
    utc_time : array-like
        time in UT
    Returns
    -------
    local_times : array-like
        local times array
    """
    offset_sec = (3600 * np.array(longitude)) / 15
    offset = pd.to_timedelta(offset_sec, unit='s')
    local_times = pd.to_datetime(utc_time) + offset
    return local_times


def compute_qd_coords(lat, lon, rad_km, time):
    """ compute magnetic coordinates from geographic

    Parameters
    ----------
    lat : array-like
        latitudes
    lon : array-like
        longitudes
    time : array-like
        time

    Returns
    -------
    mlat : array-like
        magnetic latitude
    mlon : array-like
        magnetic longitude
    """
    apex = apexpy.Apex(date=time[0])
    mlat, mlon = apex.convert(lat, lon, 'geo', 'qd')

    # Get L shells
    aalt = apex.get_apex(mlat, rad_km)
    L_shells = 1.0 + aalt / apex.RE

    return mlat, mlon, L_shells


def earth_radius(lat, Re=6378137, Rp=6356752):
    """Convert altitude in radius to m using earth's radius at a given lat.

    Parameters
    ----------
    lat : array-like
        latitude array
    Re : float
        Radius of Earth's equator in meters (default=6378137)
    Rp : float
        Radius of Earth's poles in meters (default=6356752)

    Returns
    -------
    Rearth : array-like
        Earth's radius in m at given latitudes
    """

    Rearth = []

    # iterate through latitudes
    for i, l in enumerate(lat):
        # convert latitude to raidans
        lat_rad = l * (np.pi / 180)

        # caluclat earth's raidus at a specific altitude
        eq_top = (((Re ** 2 * np.cos(lat_rad)) ** 2)
                  + ((Rp ** 2 * np.sin(lat_rad)) ** 2))
        eq_bot = (((Re * np.cos(lat_rad)) ** 2)
                  + ((Rp * np.sin(lat_rad)) ** 2))

        # take square root and append
        Rearth.append((eq_top / eq_bot) ** 0.5)

    # convert to array
    Rearth = np.array(Rearth)

    return Rearth


def compute_magnetic_coords(lats, lons, rad_km, time):
    """Compute magnetic coordinates from geographic.

    Parameters
    ----------
    lat : array-like
        latitudes
    lon : array-like
        longitudes
    rad_km : int
        altitude of satellite in km
    time : array-like
        time

    Returns
    -------
    mlat : array-like
        magnetic latitude
    mlon : array-like
        magnetic longitude
    L_shells : array-like
        L shells for each point
    """
    # Get Apex
    apex_swarm = apexpy.Apex(date=time[0])

    # Calculate at certain altitude
    alat, alon = apex_swarm.geo2apex(lats, lons, rad_km)

    # Get L shells
    aalt = apex_swarm.get_apex(alat, rad_km)
    L_shells = 1.0 + aalt / apex_swarm.RE

    return alat, alon, L_shells


def NEC_2_MFA(B_NEC):
    """ converts magnetic field in NEC coordinates to MFA coordinates

    Parameters
    ----------
    B_NEC : 3D array
        Magnetic field in NEC coordinate system
    Returns
    -------
    B_MFA : 3D array
        Magnetic field in MFA coordinate system
    B_MFA_hp60_Z : array-like
        HP filt Z in 60 sec
    B_MFA_hp1_Y : array-like
         HP filt Y in 1 sec
    B_MFA_hp1_Z : array-like
        HP filt Z in 1 sec
    """
    fs = 50  # Hz
    B_NEC_interp = calc.naninterp(B_NEC)

    # --- Step 1: Smooth NEC field with 60 s low-pass Butterworth filter ---
    cutoff = 1 / 60  # Hz
    order = 4
    b, a = butter(order, cutoff / (fs / 2), btype='low')

    B_NEC_smooth_full = np.empty_like(B_NEC)
    for i in range(3):
        B_NEC_smooth_full[:, i] = filtfilt(b, a, B_NEC_interp[:, i])

    # Convert directly using inclination/declination-based rotation
    B_MFA_NO, R = NEC_2_MFA_Rotation(B_NEC_smooth_full)

    B_MFA = np.einsum('nij,nj->ni', R, B_NEC_interp)

    # High-pass filters
    B_MFA_hp60_par = calc.highpass_filter(B_MFA[:, 0], 1 / 60, fs, order=2)
    B_MFA_hp2_zon = calc.highpass_filter(B_MFA[:, 1], 0.5, fs, order=4)
    B_MFA_hp2_par = calc.highpass_filter(B_MFA[:, 0], 0.5, fs, order=4)

    return B_MFA, B_MFA_hp60_par, B_MFA_hp2_zon, B_MFA_hp2_par


def track_to_nec_rotation(VN, VE, VC, C_is_up=True):
    """ Calculate rotation matrix to go from
    satellite track coordinate system to NEC
    Parameters
    ----------
    VN : array-like
        Satellite Velocity in Northward Direction
    VE : array-like
        Satellite Velocity in Eastward direction
    VC : array-like
        Satellite Velocity in Center direction
    C_is_up : bool kwarg
        Default True indicating that input coordinates
        assume that positive is upward
    Returns
    -------
    R : matrix
        rotation matrix to convert from Track to NEC
    Notes
    -----
    This uses satellite velocity in NEC not ion drift velocity in NEC
    """
    # compute azimuth/elevation of velocity
    if not C_is_up:
        VC = -VC

    # Calculate angles
    psi = np.arctan2(VE, VN)
    gamma = np.arctan2(VC, np.hypot(VN, VE))

    # Calculate sin and cos
    cg, sg = np.cos(gamma), np.sin(gamma)
    cp, sp = np.cos(psi), np.sin(psi)
    cg = np.cos(gamma)  # shape = (N,)

    # Set up Rotation Matrix
    R = np.zeros((VN.size, 3, 3))  # shape = (N, 3, 3)

    # Fill Rotation Matrix
    R[:, 0, 0] = cg * cp
    R[:, 0, 1] = -sp
    R[:, 0, 2] = -sg * cp

    R[:, 1, 0] = cg * sp
    R[:, 1, 1] = cp
    R[:, 1, 2] = -sg * sp

    R[:, 2, 0] = sg
    R[:, 2, 1] = 0.0
    R[:, 2, 2] = cg

    return R  # NEC <- Track


def coordinate_transform(V_along, V_cross, V_radial, R):
    """ Transform Vector in coordinates from 1 system to
    another (satellite to NEC or NEC to MFA)
    Vector can be velocity, magnetic field,
    electric field etc.
    Parameters
    ----------
    V_along : array-like
        Vector in along track direction
    V_cross : array-like
        Vector in cross-track direction
    V_radial : array-like
        Vetor in radial direction
    R : matrix
        rotation matrix defined by track_to_nec_rotation
    Returns
    -------
    VNEC : array-like
        3D array of NEC component Velocity
    """
    # Stack V track vectors
    Vtrack = np.stack([V_along, V_cross, V_radial], axis=1)  # shape = (N,3)

    # Convert to V in NEC coordinates 3D array
    VNEC = np.einsum('nij,nj->ni', R, Vtrack)

    return VNEC


def sat2nec_errors(sigma_Vixh, sigma_Vixv, sigma_Viy, sigma_Viz, R):
    """ Transform the errors from satellite track to NEC
    Parameters
    ----------
    sigma_Vixh : array-like
        error for Vixh
    sigma_Vixv : array-like
        error for Vixv
    sigma_Viy : array-like
        error for Viy
    sigma_Viz : array-like
        error for Viz
    R : matrix
        rotation matrix defined by track_to_nec_rotation
    Returns
    -------
    sigma_NEC: 3D array
        errors for each component in a 3D array
    """
    N = sigma_Vixh.size
    C_track = np.zeros((N, 3, 3))
    C_track[:, 0, 0] = sigma_Vixh ** 2 + sigma_Vixv ** 2   # along-track
    C_track[:, 1, 1] = sigma_Viy ** 2                    # cross-track
    C_track[:, 2, 2] = sigma_Viz ** 2                    # radial

    # Propagate errors into NEC frame
    C_NEC = np.einsum('nij,njk,nlk->nil', R, C_track, R)  # R * C_track * R^T

    # Extract 1-sigma NEC errors
    sigma_N = np.sqrt(C_NEC[:, 0, 0])
    sigma_E = np.sqrt(C_NEC[:, 1, 1])
    sigma_C = np.sqrt(C_NEC[:, 2, 2])

    # Create 3D array
    sigma_NEC = np.stack([sigma_N, sigma_E, sigma_C], axis=1)

    return sigma_NEC


def nec2mfa_errors(sigma_N, sigma_E, sigma_C, R):
    """Transform the errors from satellite track to NEC.
    Parameters
    ----------
    sigma_Vixh : array-like
        error for Vixh
    sigma_Vixv : array-like
        error for Vixv
    sigma_Viy : array-like
        error for Viy
    sigma_Viz : array-like
        error for Viz
    R : matrix
        rotation matrix defined by track_to_nec_rotation

    Returns
    -------
    sigma_MFA : 3D array
        error array in MFA coor
    """
    N = sigma_N.size
    C_track = np.zeros((N, 3, 3))
    C_track[:, 0, 0] = sigma_N
    C_track[:, 1, 1] = sigma_E
    C_track[:, 2, 2] = sigma_C

    # Propagate errors into NEC frame
    C_MFA = np.einsum('nij,njk,nlk->nil', R, C_track, R)  # R * C_track * R^T

    # Extract 1-sigma NEC errors
    sigma_par = np.sqrt(C_MFA[:, 0, 0])
    sigma_zon = np.sqrt(C_MFA[:, 1, 1])
    sigma_mer = np.sqrt(C_MFA[:, 2, 2])

    sigma_MFA = np.stack([sigma_par, sigma_zon, sigma_mer], axis=1)

    return sigma_MFA


def compute_declination_inclination(B_N, B_E, B_C):
    """Compute declination and inclination angles from B field (NEC).

    Parameters
    ----------
    B_N : array-like
        north component magnetic field array
    B_E : array-like
        east component magnetic field array
    B_C : array-like
        center component magnetic field array

    Returns
    -------
    declination_deg : float
        declination angle in degrees
    inclination_deg : float
        inclination angle in degrees
    """
    # Horizontal field magnitude
    B_H = np.sqrt(B_N**2 + B_E**2)

    # Declination: angle between north and horizontal field
    declination_rad = np.arctan2(B_E, B_N)
    declination_deg = np.degrees(declination_rad)  # Convert to degrees

    # Inclination: angle between horizontal plane and total field
    inclination_rad = np.arctan2(B_C, B_H)
    inclination_deg = np.degrees(inclination_rad)  # Convert to degrees

    return declination_deg, inclination_deg


def NEC_2_MFA_Rotation(B_NEC):
    """
    Transform magnetic field vectors from NEC to MFA coordinates.

    Parameters
    ----------
    B_NEC : array of shape (n, 3)
        magnetic field vectors in NEC (North, East, Center)
    Returns
    -------
    B_MFA: array of shape (n, 3)
        magnetic field vectors in MFA (parallel, zonal, meridional)
    """
    D_deg, I_deg = compute_declination_inclination(B_NEC[:, 0], B_NEC[:, 1],
                                                   B_NEC[:, 2])
    D_rad = np.radians(D_deg)
    I_rad = np.radians(I_deg)
    n = len(D_rad)

    # Preallocate rotation matrices
    R = np.zeros((n, 3, 3))

    # Fill in the rotation matrices
    R[:, 0, 0] = np.cos(I_rad) * np.cos(D_rad)
    R[:, 0, 1] = -np.sin(D_rad)
    R[:, 0, 2] = np.sin(I_rad) * np.cos(D_rad)

    R[:, 1, 0] = np.cos(I_rad) * np.sin(D_rad)
    R[:, 1, 1] = np.cos(D_rad)
    R[:, 1, 2] = np.sin(I_rad) * np.sin(D_rad)

    R[:, 2, 0] = np.sin(I_rad)
    R[:, 2, 1] = 0
    R[:, 2, 2] = -np.cos(I_rad)

    # Transform B_NEC -> B_MFA
    B_MFA = np.einsum('nij,nj->ni', R, B_NEC)

    return B_MFA, R


def ecef2ned(time, lat, lon, alt):
    """Compute NED positions from ECEF

    Parameters
    ----------
    time : array-like
        time array
    lat : array-like
        latitude array
    lon : array-like
        longitude array
    alt : array-like
        altitude array

    Returns
    -------
    v_north, v_east, v_down : array-like
        velocity arrays in NED (NEC) coordinates
    """

    # Compute ECEF positions
    ecef = geodetic_to_ecef(lat, lon, alt)  # shape (N, 3)

    # Time differences (seconds)
    dt = np.diff(time).astype('timedelta64[s]').astype(float)

    # Compute finite differences for ECEF velocities (m/s)
    v_ecef = np.full_like(ecef, np.nan)
    v_ecef[1:] = np.diff(ecef, axis=0) / dt[:, None]

    # Transform ECEF velocity to NED
    v_ned = np.full_like(v_ecef, np.nan)
    for i in range(1, len(lat)):
        R = ecef_to_ned_matrix(lat[i], lon[i])
        v_ned[i] = R @ v_ecef[i]

    # Split into components
    v_north, v_east, v_down = v_ned[:, 0], v_ned[:, 1], v_ned[:, 2]

    return v_north, v_east, v_down


def geodetic_to_ecef(lat, lon, alt):
    """Convert geodetic coordinates to ECEF (x, y, z) in meters.

    Parameters
    ----------
    lat : array-like
        latitudes
    lon : array-like
        longitudes
    alt : array-like
        altitudes

    Returns
    -------
    ECEF coordinates
    """

    a = 6378137.0          # Semi-major axis (m)
    e2 = 6.69437999014e-3  # First eccentricity squared

    lat, lon = np.radians(lat), np.radians(lon)
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
    x = (N + alt) * np.cos(lat) * np.cos(lon)
    y = (N + alt) * np.cos(lat) * np.sin(lon)
    z = (N * (1 - e2) + alt * (1 - e2)) * np.sin(lat)

    return np.column_stack((x, y, z))


def ecef_to_ned_matrix(lat, lon):
    """Rotation matrix from ECEF to NED coordinates.

    Parameters
    ----------
    lat : array-like
        latitudes
    lon : array-like
        longitudes

    Returns
    -------
    Rotation Matrix
    """
    lat, lon = np.radians(lat), np.radians(lon)
    return np.array([
        [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)],
        [-np.sin(lon), np.cos(lon), 0],
        [-np.cos(lat) * np.cos(lon), -np.cos(lat) * np.sin(lon), -np.sin(lat)]
    ])
