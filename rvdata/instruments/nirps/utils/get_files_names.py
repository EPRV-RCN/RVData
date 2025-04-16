'''
RVData/rvdata/instruments/nirps/utils/get_files_names.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Wen Mar 07 2025
Last Modified: Fri Mar 21 2025
Version: 1.0.0
Description:
Retrieves related FITS file names based on a given raw file path.
Handles DRIFT file validation and adjusts BLAZE file paths based on timestamps.

---------------------
Libraries
---------------------
'''
from astropy.io import fits
from datetime import datetime, timedelta
import os

import rvdata.instruments.nirps.config.config as config


def get_files_names(full_path: str, directory_structure: str) -> dict:
    """
    This function retrieves the names of related FITS files based on a given file's path
    and constructs a dictionary containing paths to these files.

    Args:
        full_path (str): The full path to the raw file.
        directory_structure (str): Type of database architecture that stores
            resources. Must be either 'dace' or 'standard'.

    Returns:
        dict: A dictionary with keys representing the file types and values
            containing their respective paths.
    """

    # Get the directory path and base file name from the full path
    repo_path = os.path.dirname(full_path)
    base_file_name = os.path.basename(full_path)

    if directory_structure == 'dace':
        repo_path = repo_path.replace(
            "NIRPSRAW/raw", f"NIRPSDRS/{config.DRS_VERSION}/reduced"
        )

    # Construct paths for the S2D and BLAZE FITS files (both A and B versions)
    s2d_blaze_file_A = os.path.join(
        repo_path, 'r.'+base_file_name[:-5]+'_S2D_BLAZE_A.fits'
    )
    s2d_blaze_file_B = os.path.join(
        repo_path, 'r.'+base_file_name[:-5]+'_S2D_BLAZE_B.fits'
    )
    drift_file_B = os.path.join(
        repo_path, 'r.'+base_file_name[:-5]+'_DRIFT_MATRIX_B.fits'
    )
    telluric_file_A = os.path.join(
        repo_path, 'r.'+base_file_name[:-5]+'_S2D_BLAZE_TELL_CORR_A.fits'
    )
    skysub_file_A = os.path.join(
        repo_path, 'r.'+base_file_name[:-5]+'_S2D_BLAZE_SKYSUB_A.fits'
    )

    if not os.path.isfile(drift_file_B):
        with fits.open(full_path) as hdu_raw:
            dpr_type = hdu_raw['PRIMARY'].header['HIERARCH ESO DPR TYPE']
            if dpr_type.split(",")[1] == 'SKY':
                print(
                    'SKY type doesn\'t have any DRIFT correction, '
                    'DRIFT extension will be generated with zeros'
                )
                drift_file_B = None
            else:
                print(
                    'No DRIFT_MATRIX_B file found, '
                    'DRIFT extension will be generated with zeros'
                )

    # Open the S2D BLAZE FITS file (_A version) to retrieve the BLAZE file names
    # These names are stored in specific header fields: HIERARCH ESO PRO REC1 CALn NAME
    # On Windows, ":" in file names is replaced with "_" as ":" is not allowed in
    # Windows file systems.
    with fits.open(s2d_blaze_file_A) as hdul:
        for i in hdul['PRIMARY'].header['ESO PRO REC1 CAL* CATG']:
            if 'BLAZE_A' == hdul['PRIMARY'].header[i]:
                if (os.name == 'nt'):
                    # For Windows: Replace ":" with "_" in file names
                    blaze_file_A = adjust_repo_path(
                        repo_path,
                        hdul["PRIMARY"].header[i[:-4]+'NAME'].replace(":", "_"),
                        directory_structure
                    )
                else:
                    # For non-Windows systems: Use file names as they are
                    blaze_file_A = adjust_repo_path(
                        repo_path,
                        hdul["PRIMARY"].header[i[:-4]+'NAME'],
                        directory_structure
                    )

            if 'BLAZE_B' == hdul['PRIMARY'].header[i]:
                if (os.name == 'nt'):
                    # For Windows: Replace ":" with "_" in file names
                    blaze_file_B = adjust_repo_path(
                        repo_path,
                        hdul["PRIMARY"].header[i[:-4]+'NAME'].replace(":", "_"),
                        directory_structure
                    )
                else:
                    # For non-Windows systems: Use file names as they are
                    blaze_file_B = adjust_repo_path(
                        repo_path,
                        hdul["PRIMARY"].header[i[:-4]+'NAME'],
                        directory_structure
                    )

    # Construct a dictionary of all the file paths
    # `full_path` corresponds to the raw file path (the input to this function)
    names = {
        "raw_file": full_path,
        "s2d_blaze_file_A": s2d_blaze_file_A,
        "s2d_blaze_file_B": s2d_blaze_file_B,
        "blaze_file_A": blaze_file_A,
        "blaze_file_B": blaze_file_B,
        "drift_file_B": drift_file_B,
        "telluric_file_A": telluric_file_A,
        "skysub_file_A": skysub_file_A
    }
    return names


def adjust_repo_path(
        repo_path: str, blaze_filename: str, directory_structure: str
) -> str:
    """
    Adjusts the repository path based on the timestamp in the BLAZE file name.
    If the file's timestamp is before noon, it should be placed in the previous day's
    directory.

    Args:
        repo_path (str): The original repository path.
        blaze_filename (str): The BLAZE file name containing the timestamp.
        directory_structure (str): Type of database architecture that stores
            resources. Must be either 'dace' or 'standard'.

    Returns:
        str: The corrected file path.
    """
    try:
        if (directory_structure == 'dace'):
            # Extract the date and time from the filename
            filename_parts = blaze_filename.split("_BLAZE")[0]

            timestamp_str = filename_parts.split("r.NIRPS.")[1]

            # Convert timestamp to datetime object
            if ':' in timestamp_str:
                file_datetime = datetime.strptime(timestamp_str, "%Y-%m-%dT%H:%M:%S.%f")
            else:
                file_datetime = datetime.strptime(timestamp_str, "%Y-%m-%dT%H_%M_%S.%f")

            # Extract the current repo date
            repo_date = file_datetime.date()

            # If file's time is before noon, shift repo_path to the previous day
            if file_datetime.hour < 12:
                new_repo_date = repo_date - timedelta(days=1)
                new_repo_path = os.path.join(
                    os.path.dirname(repo_path), new_repo_date.strftime("%Y-%m-%d")
                )
            else:
                # No change if file was created after noon
                new_repo_path = os.path.join(
                    os.path.dirname(repo_path), repo_date.strftime("%Y-%m-%d")
                )
            return os.path.join(new_repo_path, blaze_filename)
        return os.path.join(repo_path, blaze_filename)

    except Exception as e:
        print(f"Error processing repo path for BLAZE file '{blaze_filename}': {e}")

        # Return original path if something goes wrong
        return os.path.join(repo_path, blaze_filename)
