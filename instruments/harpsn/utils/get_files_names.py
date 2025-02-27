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
import os

'''
---------------------
internal libraries
---------------------
'''


# TODO get the file adaptativ to every configuration of fibers number
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
                blaze_file_A = os.path.join(
                    repo_path, 
                    hdul["PRIMARY"].header[i[:-4]+'NAME']
                )
            if 'BLAZE_B' == hdul['PRIMARY'].header[i]:
                blaze_file_B = os.path.join(
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