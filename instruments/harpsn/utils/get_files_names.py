'''
RVData/instruments/harpsn/utils/get_files_names.py

UNIGE-ESO - EPRV
Author: Loris JACQUES
Created: Wed Feb 26 2025
Last Modified: Wed Feb 26 2025
Version: 1.0.0
Description: 
'''

'''
---------------------
external libraries
---------------------
'''
from astropy.io import fits
from datetime import datetime, timedelta
import os

'''
---------------------
internal libraries
---------------------
'''


def get_files_names(full_path:str) -> dict:
    """
    This function retrieves the names of related FITS files based on a given file's path
    and constructs a dictionary containing paths to these files.

    Args:
        full_path (str): The full path to the raw file.

    Returns:
        dict: A dictionary with keys representing the file types and values containing their respective paths.
    """
    # Get the directory path and base file name from the full path
    repo_path = os.path.dirname(full_path)
    base_file_name = os.path.basename(full_path)
    
    # Construct paths for the S2D and BLAZE FITS files (both A and B versions)
    s2d_blaze_file_A = os.path.join(repo_path, 'r.'+base_file_name[:-5]+'_S2D_BLAZE_A.fits')
    s2d_blaze_file_B = os.path.join(repo_path, 'r.'+base_file_name[:-5]+'_S2D_BLAZE_B.fits')
    drift_file_B = os.path.join(repo_path, 'r.'+base_file_name[:-5]+'_DRIFT_MATRIX_B.fits')


    if not os.path.isfile(drift_file_B):
        with fits.open(full_path) as hdu_raw:
            dpr_type = hdu_raw['PRIMARY'].header['HIERARCH TNG DPR TYPE']
            if dpr_type.split(",")[1] == 'SKY':
                print('SKY type doesn\'t have any DRIFT correction')
                drift_file_B = None
            elif dpr_type.split(",")[1] == 'DARK':
                print('DARK type doesn\'t have any DRIFT correction')
                drift_file_B = None
            else:
                print('ERROR: NO DRIFT FILE FOUND')
                return

    # Open the S2D BLAZE FITS file (A version) to retrieve the BLAZE file names
    # These names are stored in specific header fields: HIERARCH ESO PRO REC1 CALn NAME
    with fits.open(s2d_blaze_file_A) as hdul:
        for i in hdul['PRIMARY'].header['ESO PRO REC1 CAL* CATG']:
            if 'BLAZE_A' == hdul['PRIMARY'].header[i]:
                blaze_file_A = adjust_repo_path(
                    repo_path, 
                    hdul["PRIMARY"].header[i[:-4]+'NAME']
                )
            if 'BLAZE_B' == hdul['PRIMARY'].header[i]:
                blaze_file_B = adjust_repo_path(
                    repo_path, 
                    hdul["PRIMARY"].header[i[:-4]+'NAME']
                )
    
    # Construct a dictionary of all the file paths
    # `full_path` corresponds to the raw file path (the input to this function)
    names = {
        "raw_file": full_path,
        "s2d_blaze_file_A": s2d_blaze_file_A, 
        "s2d_blaze_file_B": s2d_blaze_file_B,
        "blaze_file_A": blaze_file_A,
        "blaze_file_B": blaze_file_B,
        "drift_file_B": drift_file_B
    }
    return names


def adjust_repo_path(repo_path: str, blaze_filename: str) -> str:
    """
    Adjusts the repository path based on the timestamp in the BLAZE file name.
    If the file's timestamp is before noon, it should be placed in the previous day's directory.

    Args:
        repo_path (str): The original repository path.
        blaze_filename (str): The BLAZE file name containing the timestamp.

    Returns:
        str: The corrected file path.
    """
    try:
        # Extract the date and time from the filename
        filename_parts = blaze_filename.split("_BLAZE")[0]

        timestamp_str = filename_parts.split("r.HARPN.")[1]  
   
        # Convert timestamp to datetime object
        file_datetime = datetime.strptime(timestamp_str, "%Y-%m-%dT%H-%M-%S.%f")

        # Extract the current repo date
        repo_date = file_datetime.date()

        # If file's time is before noon, shift repo_path to the previous day
        if file_datetime.hour < 12:
            new_repo_date = repo_date - timedelta(days=1)
            new_repo_path = os.path.join(os.path.dirname(repo_path), new_repo_date.strftime("%Y-%m-%d"))
            return os.path.join(new_repo_path, blaze_filename)

        return os.path.join(repo_path, blaze_filename)  # No change if file was created after noon

    except Exception as e:
        print(f"Error processing repo path for BLAZE file '{blaze_filename}': {e}")
        return os.path.join(repo_path, blaze_filename)  # Return original path if something goes wrong

