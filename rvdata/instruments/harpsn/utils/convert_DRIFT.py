'''
RVData/rvdata/instruments/harpsn/utils/convert_DRIFT.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Wed Feb 26 2025
Last Modified: Wed Feb 26 2025
Version: 1.0.0
Description:

---------------------
Libraries
---------------------
'''
from astropy.io import fits
import numpy as np

from rvdata.core.models.level2 import RV2
import rvdata.instruments.harpsn.config.config as config


def convert_DRIFT(RV2: RV2, file_path: str) -> None:
    """
    Processes a FITS file and converts its data into a 'DRIFT' extension,
    which is then added to or updated in the provided RV2 object.

    Args:
        RV2 (RV2): The target object where the 'DRIFT' extension will be added
            or updated.
        file_path (str): The path to the FITS file containing the drift data.

    Returns:
        None: The function modifies the RV2 object in place by adding or
            updating the 'DRIFT' extension.
    """
    if (file_path is not None):
        with fits.open(file_path) as hdul:
            # Extract drift data from the FITS file (2nd HDU)
            drift_data = hdul[1].data

            drift_hdu = fits.ImageHDU(
                data=drift_data,
                header=hdul[1].header
            )
    else:
        # If no file is provided, create an empty ImageHDU with default
        # dimensions. This case occurs when Fiber B is SKY or DARK.
        drift_hdu = fits.ImageHDU(
            data=np.zeros((config.NUMORDER, config.num_pixel))
        )

    # Update the header with relevant metadata
    drift_hdu.header['EXTNAME'] = 'DRIFT'
    drift_hdu.header['CTYPE1'] = ('Pixels', 'Name of axis 1')
    drift_hdu.header['CTYPE2'] = ('Order-N', 'Name of axis 2')

    # Check if the extension already exists in the RV2 object
    if (drift_hdu.header['EXTNAME'] not in RV2.extensions):
        # If the extension does not exist, create it
        RV2.create_extension(
            ext_name=drift_hdu.header['EXTNAME'],
            ext_type='ImageHDU',
            header=drift_hdu.header,
            data=drift_hdu.data
        )
    else:
        # If the extension exists, update its data and header
        RV2.set_header(drift_hdu.header['EXTNAME'], drift_hdu.header)
        RV2.set_data(drift_hdu.header['EXTNAME'], drift_hdu.data)
