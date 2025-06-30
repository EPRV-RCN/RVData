'''
RVData/rvdata/instruments/nirps/utils/convert_SKYSUB.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Fri Mar 21 2025
Last Modified: Fri Mar 21 2025
Version: 1.0.0
Description:
This script extracts 'SKYSUB' datas from a FITS file and stores them in
an `RV2` object as extensions (e.g., 'TRACE<X>_SKYSUB_SCI',
'TRACE<X>_SKYSUB_ERR', 'TRACE<X>_SKYSUB_QUALDATA'). It ensures proper
naming, metadata handling, and updates existing extensions if necessary.

---------------------
Libraries
---------------------
'''
from astropy.io import fits

from rvdata.core.models.level2 import RV2
import rvdata.instruments.nirps.config.config as config


def convert_SKYSUB(
        RV2: RV2, file_path: str,
        trace_ind_start: int, slice_nb: int
) -> None:
    """
    This function processes a FITS file and converts its data into multiple 'SKYSUB'
    extensions, which are then added to or updated in the provided RV2 object.

    Args:
        RV2 (RV2): The target object where the 'SKYSUB' extensions will be added or
            updated.
        file_path (str): The path to the FITS file containing the data to process.
        trace_ind_start (int): The starting index for naming the extensions
            (e.g., TRACE<trace_ind_start>_SKYSUB_SCI).
        slice_nb (int): The number of slices to extract from the FITS file, each
            corresponding to one extension.

    Returns:
        None: Modifies the RV2 object in place by adding or updating the 'SKYSUB'
            extensions.
    """

    with fits.open(file_path) as hdul:
        for field in config.extnames_skysub.keys():
            # Loop through each slice from 1 to slice_nb
            for slice in range(1, slice_nb+1):
                # Extract the corresponding data for this slice
                # Each slice takes every slice_nb-th row starting from (slice-1)
                skysub_hdu = fits.ImageHDU(
                    data=hdul[field].data[slice-1::slice_nb, :],
                    header=hdul[1].header
                )

                # Update the header of the new HDU with relevant metadata
                skysub_hdu.header['EXTNAME'] = (
                    'TRACE'
                    + str(trace_ind_start+slice-1)
                    + config.extnames_skysub[field]
                )
                skysub_hdu.header['CTYPE1'] = ('Pixels', 'Name of axis 1')
                skysub_hdu.header['CTYPE2'] = ('Order-N', 'Name of axis 2')

                # Check if the extension already exists in the RV2 object
                if (skysub_hdu.header['EXTNAME'] not in RV2.extensions):
                    # If the extension does not exist, create it
                    RV2.create_extension(
                        ext_name=skysub_hdu.header['EXTNAME'],
                        ext_type='ImageHDU',
                        header=skysub_hdu.header,
                        data=skysub_hdu.data
                    )
                else:
                    # If the extension exists, update its data and header
                    RV2.set_header(skysub_hdu.header['EXTNAME'], skysub_hdu.header)
                    RV2.set_data(skysub_hdu.header['EXTNAME'], skysub_hdu.data)
