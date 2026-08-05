"""load DMSP data from cda and files from madrigal.
Updated 21 February 2026
alanahco

Note
----
from IPython.core.display import HTML was original
from IPython.display import display, HTML recommended but untested

List Functions
--------------
cda_dmsp
load_dmsp_file
"""
import os
import glob
import h5py
import pandas as pd
import numpy as np
import datetime as dt
from IPython.display import display, HTML
from cdasws import CdasWs, TimeInterval
from cdasws.datarepresentation import DataRepresentation as dr
from pySpyEPI.utils import coords
from pySpyEPI.utils.sats import dmsp_utils as ddu


def cda_dmsp(st, ed, sat_num):
    """Donwload and read DMSP data from cdaWeb
    Pramaters:
    st : datetime
        time of desired data start
    ed : datetime
        time of desired data end
    sat_num : integer
        satetllite number e.g. 15, 16, 17, 18 etc.
    """
    # convert datetime to datestring format
    t1 = st.strftime('%Y-%m-%dT%H:%M:%SZ')
    t2 = ed.strftime('%Y-%m-%dT%H:%M:%SZ')

    # convert satellite to string
    sat = f'F{sat_num}'

    # Access CDA web
    cdas = CdasWs()

    # Get date for data access
    today = dt.datetime.now().strftime('%Y-%B-%d')

    # display link information
    def display_link(url: str, text: str) -> None:
        display(HTML(f'<a href="{url}" target="_blank">{text}</a>'))

    # Get dataset
    datasets = cdas.get_datasets(observatoryGroup='DMSP')

    for index, dataset in enumerate(datasets):
        dataset_id = dataset["Id"]
        # dataset_label = dataset["Label"]

        # select this dataset for use below change if want different dataset
        if 'DMSP-' + sat + '_SSIES-3_THERMAL-PLASMA' in dataset_id:
            # Stop searching if dataset is found
            break

    # get dataset info
    doi = dataset['Doi']

    # display info about dataset that was selected
    print(f'\nSelected dataset: {dataset_id} (doi:{doi})')
    display_link(dataset['Notes'], 'Notes')
    display(HTML(f'<b>Data Citation</b>: {cdas.get_citation(doi)}. {today}.'))

    # set up variables
    variables = cdas.get_variables(doi)
    var_names = []
    for index, variable in enumerate(variables):
        name = variable['Name']
        var_names.append(name)
    time_interval = TimeInterval(t1, t2)
    _, data = cdas.get_data(doi, var_names, time_interval,
                            dataRepresentation=dr.XARRAY)

    # convert x array to dataframe
    DMSP_data = data[var_names].to_dataframe()

    # convert longitude to -180 to 180 scale
    DMSP_data.loc[DMSP_data["glon"] > 180, "glon"] -= 360

    # get local time
    DMSP_data["lt"] = coords.longitude_to_local_time(DMSP_data["glon"],
                                                     DMSP_data.index)
    return DMSP_data


def load_dmsp_file(st, ed, sat, fdir):
    """Load DMSP files downloaded from Madrigal.

    Parameters
    ----------
    st : datetime
        starting time
    ed : dateteim
        ending time
    sat : int
        16, 17, or 18 (maybe others depending on the year)
    fdir : str
        File directory

    Returns
    -------
    dmsp_data : dictionary
        contains uf, ut, and spec data
    uf_bool : bool
        True if uf data exists
    ut_bool : bool
        True if ut data exists
    spec_bool : bool
        True if spec data exists

    Notes
    -----
    There are often 3 dmsp files (ut, e, and s1). The e file provides
    the spectrum. The ut file has the most updated data and does not need
    additional processing. The s1 file provides magnetic field data and can
    be used as a more raw subsitute for the ut file if it does not exist. Most
    later dates will have a ut file.
    Adapted from code created by Grace Kwon.
    """
    # load files
    date_yr = st.strftime('%Y')
    date_str = st.strftime('%Y%m%d')

    fname = f"*{date_str}_{sat}*.hdf5"
    filename = os.path.join(fdir, date_yr, fname)
    fnum = glob.glob(filename, recursive=True)

    # check which files exist
    ut_exists = any('ut' in fn for fn in fnum)

    # intialize params
    e_dict = {}
    dmsp_ut = pd.DataFrame()
    dmsp_s1 = pd.DataFrame()
    df = pd.DataFrame()
    mdf = pd.DataFrame()
    # iterate through files
    for i in range(len(fnum)):

        # get just file name so not confused with path
        base_name = os.path.basename(fnum[i])

        #  Pull n, t, v data from  ut file ------------------------------------
        if 'ut' in base_name:
            with h5py.File(fnum[i], 'r') as file:
                d = file['/Data/Table Layout'][:]

            ut_time = pd.to_datetime(
                [f"{y}-{m:02d}-{d:02d} {h:02d}:{min_:02d}:{s:02d}"
                 for y, m, d, h, min_, s in zip(d['year'], d['month'],
                                                d['day'], d['hour'], d['min'],
                                                d['sec'])])
            data = {
                'time': ut_time,
                'mlt': d['mlt'],
                'mlat': d['mlat'],
                'lat': d['gdlat'],
                'lon': d['glon'],
                'alt': d['gdalt'],
                'ni': ddu.rem_flags(d['ni'], d['ni_rpa_flag']),
                'ti': ddu.rem_flags(d['ti'], d['ti_rpa_flag']),
                'te': ddu.rem_flags(d['te'], d['ti_rpa_flag']),
                'alongV': ddu.rem_flags(d['ion_v_sat_for'],
                                        d['ion_v_for_flag']),
                'crossV': -ddu.rem_flags(d['ion_v_sat_left'],
                                         d['ion_v_left_flag']),
                'vertV': ddu.rem_flags(d['vert_ion_v'], d['ion_v_up_flag'])
            }

            data_uf = {
                'time': ut_time,
                'mlt': d['mlt'],
                'mlat': d['mlat'],
                'lat': d['gdlat'],
                'lon': d['glon'],
                'alt': d['gdalt'],
                'ni': d['ni'],
                'ti': d['ti'],
                'te': d['te'],
                'alongV': d['ion_v_sat_for'],
                'crossV': -d['ion_v_sat_left'],
                'vertV': d['vert_ion_v']
            }

            # Create a DataFrame
            df = pd.DataFrame(data)
            df_uf = pd.DataFrame(data_uf)
            dmsp_ut = ddu.join_magvel(df)

        # Pull magnetic field data from bc file -------------------------------
        elif 's1' in base_name:
            with h5py.File(fnum[i], 'r') as file:
                d = file['/Data/Table Layout'][:]

            mag_time = pd.to_datetime(
                [f"{y}-{m:02d}-{d:02d} {h:02d}:{min_:02d}:{s:02d}" for y, m, d,
                 h, min_, s in zip(d['year'], d['month'], d['day'], d['hour'],
                                   d['min'], d['sec'])])

            magdata = {
                'magtime': mag_time,
                'alongB': d["b_forward"],
                'crossB': d["b_perp"],
                'vertB': d["bd"]
            }
            mdf = pd.DataFrame(magdata)

            if not ut_exists:
                ut_time = mag_time
                data_uf = {
                    'time': ut_time,
                    'mlt': d['mlt'],
                    'mlat': d['mlat'],
                    'lat': d['gdlat'],
                    'lon': d['glon'],
                    'alt': d['gdalt'],
                    'ni': d['ne'],
                    'ti': d['ti'],
                    'te': d['te'],
                    'alongV': d['ion_v_sat_for'],
                    'crossV': -d['ion_v_sat_left'],
                    'vertV': d['vert_ion_v']
                }
                df_uf = pd.DataFrame(data_uf)
                dmsp_s1 = ddu.join_magvel(df_uf)

        # Pull spectrum data from e file --------------------------------------
        elif 'e' in base_name:
            with h5py.File(fnum[i], 'r') as file:
                d = file['/Data/Table Layout'][:]

            e_gdlat = d['gdlat'][::19]
            e_glon = d['glon'][::19]
            e_gdalt = d['gdalt'][::19]

            e_time = pd.to_datetime(
                [f"{y}-{m:02d}-{d:02d} {h:02d}:{min_:02d}:{s:02d}"
                 for y, m, d, h, min_, s in zip(d['year'], d['month'],
                                                d['day'], d['hour'], d['min'],
                                                d['sec'])])

            e_time = np.array(e_time)
            timelen = int(len(e_time) / 19)
            e_time = e_time.reshape(timelen, 19)[:, 0]

            e_ch_energy = d['ch_energy']
            e_ch_energy = np.rot90(e_ch_energy.reshape(timelen, 19), 2)[0]

            # make sure log10 not computed where value is 0 or negative
            e_ion_d_energy = np.full_like(d['ion_d_ener'], np.nan, dtype=float)

            np.log10(
                d['ion_d_ener'],
                where=d['ion_d_ener'] > 0,
                out=e_ion_d_energy
            )

            e_el_d_energy = np.full_like(d['el_d_ener'], np.nan, dtype=float)

            np.log10(
                d['el_d_ener'],
                where=d['el_d_ener'] > 0,
                out=e_el_d_energy
            )

            e_el_d_energy = np.rot90(e_el_d_energy.reshape(timelen, 19), 1)
            e_ion_d_energy = np.rot90(e_ion_d_energy.reshape(timelen, 19), 1)

            e_dict = {
                "time": pd.to_datetime(e_time).to_pydatetime(),
                "gdlat": e_gdlat,
                "glon": e_glon,
                "gdalt": e_gdalt,
                "ch_energy": e_ch_energy,
                "el_d_energy": e_el_d_energy,
                "ion_d_energy": e_ion_d_energy,
            }

    # spec
    if not e_dict:
        spec_bool = False
    else:
        spec_bool = True

    # ut
    if len(dmsp_ut) == 0:
        ut_bool = False
    else:
        ut_bool = True

    # uf
    if len(dmsp_s1) == 0:
        uf_bool = False
    else:
        uf_bool = True

    # mag
    if len(mdf) == 0:
        mag_bool = False
    else:
        mag_bool = True

    # Mask data by time if exists
    # add local time
    if spec_bool:
        mask = (e_dict["time"] >= st) & (e_dict["time"] <= ed)
        filt_dc = {}

        for key, value in e_dict.items():
            if not isinstance(value, np.ndarray):
                continue

            # 1D time-based arrays
            if value.ndim == 1 and value.shape[0] == mask.size:
                filt_dc[key] = value[mask]

            # 2D (energy, time) arrays
            elif value.ndim == 2 and value.shape[1] == mask.size:
                filt_dc[key] = value[:, mask]

            # non-time arrays (e.g. energy bins)
            else:
                filt_dc[key] = value
        e_dict = filt_dc.copy()

    if uf_bool:
        dmsp_s1 = dmsp_s1[(dmsp_s1['time'] > st) & (dmsp_s1['time'] < ed)]
        lt = coords.longitude_to_local_time(dmsp_s1['lon'], dmsp_s1['time'])
        dmsp_s1 = dmsp_s1.copy()
        dmsp_s1['lt'] = lt

        dmsp_s1['lt_hr'] = dmsp_s1['lt'].apply(
            lambda tim: tim.hour + tim.minute / 60 + tim.second / 3600)

    if ut_bool:
        dmsp_ut = dmsp_ut[(dmsp_ut['time'] > st) & (dmsp_ut['time'] < ed)]
        lt = coords.longitude_to_local_time(dmsp_ut['lon'], dmsp_ut['time'])
        dmsp_ut = dmsp_ut.copy()
        dmsp_ut['lt'] = lt
        dmsp_ut['lt_hr'] = dmsp_ut['lt'].apply(
            lambda tim: tim.hour + tim.minute / 60 + tim.second / 3600)

    if mag_bool:
        mag_df = mdf[(mdf['magtime'] > st) & (mdf['magtime'] < ed)]
    else:
        mag_df = pd.DataFrame()

    dmsp_data = {}
    dmsp_data['s1'] = dmsp_s1
    dmsp_data['ut'] = dmsp_ut
    dmsp_data['spec'] = e_dict
    dmsp_data['mag'] = mag_df

    if mag_bool & ut_bool:
        time = dmsp_ut['time'].values
        gdlat = dmsp_ut['lat'].values
        glon = dmsp_ut['lon'].values
        gdalt = dmsp_ut['alt'].values

        alongV = dmsp_ut['alongV'].values
        crossV = dmsp_ut['crossV'].values
        vertV = dmsp_ut['vertV'].values
        alongB = mag_df['alongB'].values
        crossB = mag_df['crossB'].values
        vertB = mag_df['vertB'].values

        try:
            Vpar, Vzon, Vmer = ddu.rot_sat2field(time, gdlat, glon, gdalt,
                                                 alongV, crossV, vertV, alongB,
                                                 crossB, vertB)
        except Exception:
            mag_bool = False

        if mag_bool:
            dmsp_data['ut']['Vparallel'] = Vpar
            dmsp_data['ut']['Vzonal'] = Vzon
            dmsp_data['ut']['Vmeridional'] = Vmer

    return dmsp_data, uf_bool, ut_bool, spec_bool, mag_bool
