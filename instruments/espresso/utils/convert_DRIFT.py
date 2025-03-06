'''
RVData/instruments/espresso/utils/convert_DRIFT.py

UNIGE-ESO - EPRV
Author: Loris JACQUES
Created: Mon Mar 03 2025
Last Modified: Mon Mar 03 2025
Version: 1.0.0
Description: 
'''

'''
---------------------
external libraries
---------------------
'''
from astropy.io import fits
import numpy as np

'''
---------------------
internal libraries
---------------------
'''
from core.models.level2 import RV2
import instruments.espresso.config.config as config

def convert_DRIFT(RV2: RV2, file_path: str) -> None:
    """
    Processes a FITS file and converts its data into a 'DRIFT' extension, 
    which is then added to or updated in the provided RV2 object.

    Args:
        RV2 (RV2): The target object where the 'DRIFT' extension will be added or updated.
        file_path (str): The path to the FITS file containing the drift data.

    Returns:
        None: The function modifies the RV2 object in place by adding or updating the 'DRIFT' extension.
    """

    if(file_path != None):
        with fits.open(file_path) as hdul:
            # Extract drift data from the FITS file (2nd HDU)
            drift_hdu = fits.ImageHDU(
                data = hdul[1].data,
                header = hdul[1].header
            )
    else:
        # If no file is provided, create an empty ImageHDU with default dimensions
        # This case occurs when Fiber B is SKY or DARK.
        drift_hdu = fits.ImageHDU(
            data = np.zeros((config.NUMORDER, config.num_pixel))
        )

    # Update the header with relevant metadata
    drift_hdu.header['EXTNAME'] = 'DRIFT'
    drift_hdu.header['CTYPE1'] = ('Pixels', 'Name of axis 1')
    drift_hdu.header['CTYPE2'] = ('Order-N', 'Name of axis 2')

    # Check if the extension already exists in the RV2 object
    if(drift_hdu.header['EXTNAME'] not in RV2.extensions):
        # If the extension does not exist, create it
        RV2.create_extension(
            ext_name = drift_hdu.header['EXTNAME'], 
            ext_type = 'ImageHDU', 
            header = drift_hdu.header, 
            data = drift_hdu.data
        )
    else:
        # If the extension exists, update its data and header
        RV2.set_header(drift_hdu.header['EXTNAME'], drift_hdu.header)
        RV2.set_data(drift_hdu.header['EXTNAME'], drift_hdu.data)


def add_nan_row(matrix: np.ndarray, row_index: int) -> np.ndarray:
    """
    Inserts a row of NaN values at a specified index in a 2D NumPy array.

    Parameters:
        matrix (np.ndarray): The original 2D array.
        row_index (int): The index at which the NaN row should be inserted.

    Returns:
        matrix_updated (np.ndarray): A new array with the NaN row inserted.
    """

    # Force the array to be of type float to avoid insertion issues
    matrix = matrix.astype(np.float64)

    # Create a row filled with NaN values
    nan_row = np.full((1, matrix.shape[1]), np.nan)

    # Insert the NaN row into the array
    matrix_updated = np.insert(matrix, row_index, nan_row, axis=0)
    return matrix_updated
