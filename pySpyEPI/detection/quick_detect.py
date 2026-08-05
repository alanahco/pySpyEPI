"""Quick Detection with standard parameters.
alanahco
Created 14 April 2026

List Functions
--------------
detect_swarm
"""
from pySpyEPI.io.load import swarm_load
from pySpyEPI.barrel import barrel_roll
from pySpyEPI.detection.stats import cand_info
from pySpyEPI.detection import midlat_trough as ml_tr


def detect_swarm(st, ed, satellite, fdir, EFI_bool=False, trough=False):
    """Load and detect Swarm EPIs from given time.

    Parameters
    ----------
    st : datetime
        starting time of swarm pass
    ed : datetime
        ending time of swarm pass
    satelltie : string
        Swarm Satellite letter 'A', 'B', or 'C'
    fdir : string
        data directory
    EFI_bool : boolean
        if True, EFI data will be returned
        if False (default), None will be returned
    trough : boolean
        if True, trough lats will be calculated
        if False, +'- 45 degrees will be used
    """

    # Load data bsed on st and ed ---------------------------------------------
    data, IRR_bool, EFI_bool2 = swarm_load.load_all(
        st, ed, satellite, fdir, EFI=False)

    if IRR_bool:
        irr = data['IRR']

    # find trough or limit by magnetic latitude -------------------------------
    if trough:
        tr_lats = ml_tr.find_trough_lat(
            irr, equator_bound=35, auroral_bound=70, set_lat=45, ne_str='ne')
        lat_mask = ((irr['mlat'] >= tr_lats[1]) & (irr['mlat'] <= tr_lats[0]))
        irr = irr.copy()
        irr = irr[lat_mask]
    else:
        irr = irr.copy()
        irr = irr[abs(irr['mlat']) <= 45]

    # Get and limit EFI data --------------------------------------------------
    if EFI_bool:
        efi = swarm_load.load_EFI(irr['time'].iloc[0], irr['time'].iloc[-1],
                                  satellite, fdir)
    else:
        efi = None

    # Barrel Roll -------------------------------------------------------------
    barrel_df, peaks, properties = barrel_roll.triple_barrel(
        irr, barrel_start=8, det_filt=5, peak_width=6, num_barrels=5,
        prom_height=5000, x_type='time', x_str='time', ne_str='ne',
        e_up_large=0.1, e_lo_large=0.05, e_up_small=0.1,
        e_lo_small=0.05, upper_weight=2, lower_weight=1)

    # get candidate info ------------------------------------------------------
    time2 = irr['time'].values
    time2 = time2.astype('datetime64[us]')
    lat2 = irr['lat'].values
    lon2 = irr['lon'].values
    mlat2 = irr['mlat'].values
    mlon2 = irr['mlon'].values
    lt2 = irr['lt'].values
    lt2 = lt2.astype('datetime64[us]')
    lshell2 = irr['l_shell'].values
    alt2 = irr['altitude'].values

    stats_df = cand_info.epi_stats(time2, lat2, lon2, mlat2, mlon2, lt2,
                                   lshell2, alt2, peaks, properties, barrel_df,
                                   satellite, fft_test=True)

    return irr, efi, barrel_df, peaks, properties, stats_df
