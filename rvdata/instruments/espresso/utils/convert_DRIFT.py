'''
RVData/rvdata/instruments/espresso/utils/convert_DRIFT.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Mon Mar 03 2025
Last Modified: Mon Mar 03 2025
Version: 1.0.0
Description:
Extracts 'DRIFT' calibration data from a FITS file and stores it in an `RV2`
object. If no file is provided, an empty DRIFT extension is created.

---------------------
Libraries
---------------------
'''
from astropy.io import fits
import numpy as np

from rvdata.core.models.level2 import RV2
import rvdata.instruments.espresso.config.config as config


def convert_DRIFT(
        RV2: RV2, file_path: str,
        trace_ind_start: int, slice_nb: int
) -> None:
    """
    Processes a FITS file and converts its data into a 'DRIFT' extension,
    which is then added to or updated in the provided RV2 object.

    Args:
        RV2 (RV2): The target object where the 'DRIFT' extension will be added
            or updated.
        file_path (str): The path to the FITS file containing the drift data.
        trace_ind_start (int): The starting index for naming the extensions
            (e.g., TRACE<trace_ind_start>_DRIFT).
        slice_nb (int): The number of slices to extract from the FITS file,
            each corresponding to one extension.
    Returns:
        None: The function modifies the RV2 object in place by adding or
            updating the 'DRIFT' extension.
    """
    # Loop through each slice from 1 to slice_nb
    for slice in range(1, slice_nb+1):
        if (file_path is not None):
            with fits.open(file_path) as hdul:
                # Extract drift data from the FITS file (2nd HDU)
                drift_hdu = fits.ImageHDU(
                    data=hdul[1].data[slice-1::slice_nb, :],
                    header=hdul[1].header
                )
        else:
            # If no file is provided, create an empty ImageHDU with default
            # dimensions. This case occurs when Fiber B is SKY or DARK.
            drift_hdu = fits.ImageHDU(
                data=np.zeros((config.NUMORDER, config.num_pixel), dtype=np.float32)
            )

        # Update the header with relevant metadata
        drift_hdu.header['EXTNAME'] = (
            'TRACE'+str(trace_ind_start+slice-1)+'_DRIFT'
        )
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
