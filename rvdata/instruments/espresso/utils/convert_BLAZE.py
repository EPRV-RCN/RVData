"""
RVData/rvdata/instruments/espresso/utils/convert_S2D_BLAZE.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Mon Mar 03 2025
Last Modified: Mon Mar 03 2025
Version: 1.0.0
Description:
This script extracts 'BLAZE' calibration data from a FITS file and stores it in
an `RV2` object as extensions (e.g., 'TRACE<X>_BLAZE'). It ensures proper
naming, metadata handling, and updates existing extensions if necessary.

---------------------
Libraries
---------------------
"""

from astropy.io import fits

from rvdata.core.models.level2 import RV2


def convert_BLAZE(
    RV2: RV2, file_path: str, trace_ind_start: int, slice_nb: int
) -> None:
    """
    This function processes a FITS file and converts its data into multiple
    'BLAZE' extensions, which are then added to or updated in the provided RV2
    object.

    Args:
        RV2 (RV2): The target object where the 'BLAZE' extensions will be
            added or updated.
        file_path (str): The path to the FITS file containing the data to
            process.
        trace_ind_start (int): The starting index for naming the extensions
            (e.g., TRACE<trace_ind_start>_BLAZE).
        slice_nb (int): The number of slices to extract from the FITS file,
            each corresponding to one extension.

    Returns:
        None: Modifies the RV2 object in place by adding or updating the
            'BLAZE' extensions.
    """

    with fits.open(file_path) as hdul:
        # Loop through each slice from 1 to slice_nb
        for slice in range(1, slice_nb + 1):
            # Extract the corresponding data for this slice
            # Each slice takes every slice_nb-th row starting from (slice-1)
            single_cam_values = hdul[1].data[slice - 1 :: slice_nb, :]

            # Update the header of the new HDU with relevant metadata
            blaze_hdu = fits.ImageHDU(data=single_cam_values, header=hdul[1].header)

            # Update the header of the new HDU with relevant metadata
            blaze_hdu.header["EXTNAME"] = (
                "TRACE" + str(trace_ind_start + slice - 1) + "_BLAZE"
            )
            blaze_hdu.header["CTYPE1"] = ("Pixels", "Name of axis 1")
            blaze_hdu.header["CTYPE2"] = ("Order-N", "Name of axis 2")

            # Check if the extension already exists in the RV2 object
            if blaze_hdu.header["EXTNAME"] not in RV2.extensions:
                # If the extension does not exist, create it
                RV2.create_extension(
                    ext_name=blaze_hdu.header["EXTNAME"],
                    ext_type="ImageHDU",
                    header=blaze_hdu.header,
                    data=blaze_hdu.data,
                )
            else:
                # If the extension exists, update its data and header
                RV2.set_header(blaze_hdu.header["EXTNAME"], blaze_hdu.header)
                RV2.set_data(blaze_hdu.header["EXTNAME"], blaze_hdu.data)
