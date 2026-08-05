"""Download Swarm Files.
Updated 19 February 2026
alanahco
List Functions
--------------
download_and_unzip
download_advanced_files
downloadVIRES
"""
import datetime as dt
from dateutil.relativedelta import relativedelta
import glob
import requests
import os
import zipfile
from viresclient import SwarmRequest
import tqdm


def download_and_unzip(
        ymd, satellite, out_dir,
        s_url="https://swarm-diss.eo.esa.int/?do=download&file=swarm%2FLevel",
        level='1b', baseline='Latest_baselines',
        instrument1='EFI', instrument2='LP',
        f_end='0602', T1='000000', T2='235959', num_days=0, pr=False):
    """Download Swarm daily File.

        Parameters
        ----------
        ymd: datetime object
            year month day of desired swarm file
        satellite : str
            satellite string 'A', 'B', or 'C'
        out_dir : str
            directory string for file output
        bse_url : str
            Base URL where data can be found before Level specification
        level : str kwarg
            This works for Level1b, not tested n Level2daily
        baseline : str kwarg
            desired baseline 'Latest_baselines' (defaut) or
            'Entire_mission_data' (not tested)
        instrument1 : str kwarg
            desired insturment default is 'EFI'
            for Electric Field Instrument
        instrument2 : str kwarg
            desired dataset from instrument1 default is 'LP'
            for Langmuir Probe
        f_end : str kwarg
            For different data products there are different numbers at the end
            The most common for EFIxLP is (Default) '0602_MDR_EFI_LP' where
             0602 represents the file version
             MDR_EFI_LP represents the Record Type
        T1 : str kwarg
            starting Time string format "HHMMSS"
            Most files contain "000000" to start, but if the file is not the
            whole day it will be something else
            Check website if download fails
        T2 : str kwarg
            ending Time string format "HHMMSS"
            Most files contain "235959" to end, but if the file is not the
            whole day it will be something else
            Check website if download fails
        num_days : int kwarg
            number of days from starting date to be downloaded after initial
            file (default is 0)
        pr : boolean
            if True, already existing filenames will be printed
            default False
        Returns
        -------
            No returns, just file downloaded to out_dir

        Notes
        -----
            Default is an EFI file
            Options found at https://swarm-diss.eo.esa.int/#
            File format found at
                https://swarmhandbook.earth.esa.int/article/product
    """
    # Adjsut the name based on if it is level 1b or level 2daily
    if level == '1b':
        full_url = (s_url + level + "%2F" + baseline + "%2F" + instrument1
                    + 'x_' + instrument2)
    elif level == '2daily':
        full_url = (s_url + level + "%2F" + baseline + "%2F" + instrument1
                    + "%2F" + instrument2)
    # Out Folder
    yer = ymd.year
    mnth = ymd.month
    dy = ymd.day
    out_folder = f'{out_dir}/{instrument1}/Sat_{satellite}/{yer}/'

    # Make the path if it does not exist
    if not os.path.exists(out_folder):
        print(f'Making path {out_folder}')
        os.makedirs(out_folder)

    # Start at first day and go for num_days
    start_date = dt.datetime(yer, mnth, dy)
    end_date = start_date + relativedelta(days=num_days)

    # Start with start date and go until end date is reached
    current_date = start_date
    while current_date <= end_date:
        date_str = current_date.strftime("%Y%m%d")
        f_bse = "SW_OPER_"
        d_str = date_str + "T" + T1 + "_" + date_str + "T" + T2 + "_" + f_end
        d_check = date_str + "T" + T1 + "_" + date_str + "T" + T2 + "_"

        if level == '1b':
            filename = (f_bse + instrument1 + satellite + "_" + instrument2
                        + "_1B_" + d_str + ".CDF.ZIP")
        elif level == '2daily':
            filename = (f_bse + instrument1 + satellite + instrument2 + "_2F_"
                        + d_str + ".ZIP")

        # Set Full File URL
        file_url = full_url + "%2FSat_" + satellite + "%2F" + filename
        zip_path = os.path.join(out_folder, filename)
        current_date = current_date + dt.timedelta(days=1)
        extract_folder = os.path.join(out_folder, date_str)

        # Find file if it already exists or another version (e.g. 0602 vs 0701)
        if level == '1b':
            efile = (f_bse + instrument1 + satellite + "_" + instrument2
                     + "_1B_" + d_check + "*.cdf")
        elif level == '2daily':
            efile = (f_bse + instrument1 + satellite + instrument2 + "_2F_"
                     + "*.cdf")

        extracted_files = os.path.join(extract_folder, efile)
        found_file = extracted_files

        if len(glob.glob(extracted_files)) > 0:
            found_file = glob.glob(extracted_files)[0]
        if os.path.exists(found_file):
            if pr:
                print(f"File already exists: {found_file}.Skipping download.")
        else:
            # Download file
            response = requests.get(file_url)
            if response.status_code == 200:
                with open(zip_path, 'wb') as f:
                    f.write(response.content)
                print("Downloading: " + filename)

                # Unzip file into date folder
                extract_folder = os.path.join(out_folder, date_str)
                os.makedirs(extract_folder, exist_ok=True)

                try:
                    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(extract_folder)
                    print("Extracted to:" + extract_folder)
                    os.remove(zip_path)
                except zipfile.BadZipFile:
                    print(f"Failed filename {filename} does not exist")


def download_advanced_files(ymd, satellite, out_dir, s_url=None,
                            data_type='Plasma_Data',
                            dataset='16Hz_TII_Cross-track_Dataset',
                            f_end='0401', T1='000000', T2='235959', num_days=0,
                            pr=False):
    """ Download Swarm daily File
        Parameters
        ----------
        ymd: datetime object
            year month day of desired swarm file
        satellite : str
            satellite string 'A', 'B', or 'C'
        out_dir : str
            directory string for file output
        s_url : str or NoneType
            if None, s_url is set to advanced files
            Base URL where data can be found before Level specification
        data_type : str kwarg
            default Plasma_Data
            Probably will break for other data types
        data_set : str kwarg
            dataset type
            defualt 16Hz_TII_Cross-track_Dataset
            for TII velocity dataset
        f_end : str kwarg
            For different data products there are different numbers at the end
            The most common for EFIxLP is (Default) '0602_MDR_EFI_LP' where
             0602 represents the file version
             MDR_EFI_LP represents the Record Type
        T1 : str kwarg
            starting Time string format "HHMMSS"
            Most files contain "000000" to start, but if the file is not the
            whole day it will be something else
            Check website if download fails
        T2 : str kwarg
            ending Time string format "HHMMSS"
            Most files contain "235959" to end, but if the file is not the
            whole day it will be something else
            Check website if download fails
        num_days : int kwarg
            number of days from starting date to be downloaded after initial
            file (default is 0)
        pr : boolean
            if True, already existing filenames will be printed
            default False
        Returns
        -------
            No returns, just file downloaded to out_dir
        Notes
        -----
            Default is an EFI file
            Options found at https://swarm-diss.eo.esa.int/#
            File format found at
                https://swarmhandbook.earth.esa.int/article/product
    """
    if s_url is None:
        url1 = "https://swarm-diss.eo.esa.int/"
        url2 = "?do=download&file=swarm%2FAdvanced"
        s_url = url1 + url2
    # Adjsut the name based on if it is level 1b or level 2daily
    full_url = (s_url + "%2F" + data_type + "%2F" + dataset)

    # Out Folder
    yer = ymd.year
    mnth = ymd.month
    dy = ymd.day
    out_folder = os.path.join(out_dir, 'Sat_' + satellite, str(yer))

    # Make the path if it does not exist
    if not os.path.exists(out_folder):
        print(f'Making path {out_folder}')
        os.makedirs(out_folder)

    # Start at first day and go for num_days
    start_date = dt.datetime(yer, mnth, dy)
    end_date = start_date + relativedelta(days=num_days)

    # Start with start date and go until end date is reached
    current_date = start_date
    while current_date <= end_date:
        date_str = current_date.strftime("%Y%m%d")

        # Exapmle:
        # SW_EXPT_EFIA_TCT16_20131212T000000_20131212T235959_0401.ZIP
        f_bse = f'SW_EXPT_EFI{satellite}_TCT16_'
        d_str = date_str + "T" + T1 + "_" + date_str + "T" + T2 + "_" + f_end

        filename = (f_bse + d_str + ".ZIP")

        # Set Full File URL
        file_url = full_url + "%2F%2FSat_" + satellite + "%2F" + filename
        zip_path = os.path.join(out_folder, filename)
        current_date = current_date + dt.timedelta(days=1)
        extract_folder = os.path.join(out_folder, date_str)

        # Find file if it already exists
        efile = (f_bse + d_str + "*.cdf")

        extracted_files = os.path.join(extract_folder, efile)
        found_file = extracted_files

        if len(glob.glob(extracted_files)) > 0:
            found_file = glob.glob(extracted_files)[0]
        if os.path.exists(found_file):
            if pr:
                print(f"File already exists: {found_file}.Skipping download.")
        else:
            # Download file
            response = requests.get(file_url)
            if response.status_code == 200:
                with open(zip_path, 'wb') as f:
                    f.write(response.content)
                print("Downloading: " + filename)

                # Unzip file into date folder
                extract_folder = os.path.join(out_folder, date_str)
                os.makedirs(extract_folder, exist_ok=True)

                try:
                    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(extract_folder)
                    print("Extracted to:" + extract_folder)
                    os.remove(zip_path)
                except zipfile.BadZipFile:
                    print(f"Failed filename {filename} does not exist")


def downloadVIRES(ymd, satellite, base_dir, collection, measurements,
                  f_end, num_days=0):
    """
    Download data from Swarm VIRES
    Parameters
    ----------
    ymd: datetime object
            year month day of desired swarm file
        satellite : str
            satellite string 'A', 'B', or 'C'
        base_dir : str
            directory string for file output
        collection : str
            VIRES collection
            from viresclient import SwarmRequest
            request = SwarmRequest()
            request.available_collections()
        measurements : list of strings
            variables to download
        f_end : str kwarg
            ending of the filename such as .nc
        num_days : int kwarg
            number of days from starting date to be downloaded after initial
            file (default is 0)
        Returns
        -------
            No returns, just file downloaded to out_dir
    Notes
    ----
    For SW_OPER_MAG{satellite}_HR_1B,
        base_dir = "/Users/alanahco/University of Michigan Dropbox
            /ENGIN-Zou-Research-Group/SWARM/MAG_VIRES"
        collection = SW_OPER_MAG{satellite}_HR_1B
        measurments =
            ["B_VFM", "B_NEC", "Flags_B", "Flags_q","B_error"]
        f_end = "MDR_MAG_HR.nc"
    For SW_EXPT_EFI{satellite}_TCT16,
        base_dir = "/Users/alanahco/University of Michigan Dropbox
            /ENGIN-Zou-Research-Group/SWARM/EFI_TII_VIRES"
        collection = f'SW_EXPT_EFI{Sat}_TCT16'
        measurements = ['VsatC', 'VsatE', 'VsatN', 'Bx', 'By', 'Bz', 'Ehx',
                'Ehy', 'Ehz', 'Evx', 'Evy', 'Evz', 'Vicrx', 'Vicry', 'Vicrz',
                'Vixv', 'Vixh', 'Viy', 'Viz', 'Vixv_error', 'Vixh_error',
                'Viy_error', 'Viz_error', 'Latitude_QD', 'MLT_QD',
                'Calibration_flags', 'Quality_flags']
        f_end = "0401.nc"
    """
    # Set up Dates
    yer = ymd.year
    mnth = ymd.month
    dy = ymd.day

    # start and end dates depending on ymd and num_days
    start_date = dt.datetime(yer, mnth, dy)
    end_date = start_date + relativedelta(days=num_days)

    # Start with start date and go until end date is reached
    current_date = start_date

    while current_date <= end_date:

        # Get end of day
        eod = current_date + relativedelta(days=1) - relativedelta(seconds=1)

        # date strings
        start_date_string = current_date.strftime("%Y-%m-%dT%H:%M:%SZ")
        end_date_string = eod.strftime("%Y-%m-%dT%H:%M:%SZ")

        # Build output path
        out_dir = os.path.join(base_dir, 'Sat_' + satellite, str(yer),
                               f"{current_date.strftime('%Y%m%d')}")

        # Make the path if it does not exist
        if not os.path.exists(out_dir):
            print(f'Making path {out_dir}')
            os.makedirs(out_dir)

        filename = (
            f"{collection}"
            f"{current_date.strftime('%Y%m%dT%H%M%S')}_"
            f"{(eod).strftime('%Y%m%dT%H%M%S')}_{f_end}"
        )

        full_path = os.path.join(out_dir, filename)

        if os.path.exists(full_path):
            print(f"Already exists: {full_path} — skipping.")
            current_date = current_date + dt.timedelta(days=1)
            continue

        mnth = current_date.month
        dy = current_date.day
        print(f"\nDownloading for {yer}-{mnth:02}-{dy:02}")

        request = SwarmRequest()

        request.set_collection(collection)

        request.set_products(
            measurements=measurements,
        )

        tqdm.tqdm.disable = True
        # Get the data
        data = request.get_between(start_date_string, end_date_string)
        ds = data.as_xarray()

        # If the dtype is category then it will not save as a netcdf
        for var in ds.variables:
            if str(ds[var].dtype).startswith('category'):
                ds[var] = ds[var].astype(str)

        # Save
        ds.to_netcdf(full_path)
        print(f"Saved: {full_path}")

        # Update Current Date
        current_date = current_date + dt.timedelta(days=1)
