"""
RVData/rvdata/instruments/harps/utils/convert_S2D_BLAZE.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Wed Feb 26 2025
Last Modified: Wed Feb 26 2025
Version: 1.0.0
Description:
Extracts and processes data from an S2D_BLAZE FITS file. Stores key
calibration values in an `RV2` object, applies a Doppler shift correction
if needed, and organizes data into FITS extensions. Inserts a NaN row at
a specified index and updates existing data if necessary.

---------------------
Libraries
---------------------
"""

from astropy.io import fits
from astropy.constants import c
import numpy as np


from rvdata.core.models.level2 import RV2
import rvdata.instruments.harps.config.config as config


def convert_S2D_BLAZE(
    RV2: RV2, file_path: str, trace_ind_start: int, slice_nb: int
) -> None:
    """
    Extracts and processes relevant data from an S2D_BLAZE FITS file.

    This function reads the specified FITS file, extracts relevant
    information, and stores the processed data in the RV2 object. The function
    ensures that certain header values are only set once during the first
    iteration. It also applies a Doppler shift correction when necessary and
    organizes the extracted data into new FITS extensions.

    Args:
        RV2 (RV2): The RV2 object where extracted data will be stored.
        file_path (str): The path to the FITS file.
        trace_ind_start (int): The starting index for trace processing.
        slice_nb (int): The number of slices to process.

    Returns:
        None
    """

    with fits.open(file_path) as hdul:
        #  Execute only on the first iteration
        if trace_ind_start == 1:
            # Set the instrument header from the primary header
            RV2.set_header("INSTRUMENT_HEADER", hdul["PRIMARY"].header)

            # Retrieve barycentric correction and bjd values
            barycorr_kms_data = hdul["PRIMARY"].header["HIERARCH ESO QC BERV"]
            bjd_tdb_data = hdul["PRIMARY"].header["HIERARCH ESO QC BJD"]

            barycorr_z_data = (barycorr_kms_data / (c.to("km/s"))).value

            # Add the extensions
            add_fits_extension(RV2, "BARYCORR_KMS", barycorr_kms_data)
            add_fits_extension(RV2, "BARYCORR_Z", barycorr_z_data)
            add_fits_extension(RV2, "BJD_TDB", bjd_tdb_data)

        # Loop through configured fields and process data
        for field in config.extnames.keys():
            for slice in range(1, slice_nb + 1):
                # We extract the values of the specific slice
                single_cam_values = hdul[field].data[slice - 1 :: slice_nb, :]

                # Remove barycentric correction if applicable
                if "BARY" in field:
                    single_cam_values = doppler_shift(
                        single_cam_values, RV2.data["BARYCORR_KMS"][0]
                    )

                # Insert a row of NaN at the specified index if len don't match
                if len(single_cam_values) != config.NUMORDER:
                    single_cam_values = add_nan_row(
                        single_cam_values, config.empty_raw_order
                    )

                # We create a new HDU with the values of the slice and the
                # original header
                hdu_l2 = fits.ImageHDU(
                    data=single_cam_values.copy(), header=hdul[field].header
                )

                # Update header fields
                hdu_l2.header["EXTNAME"] = (
                    "TRACE" + str(trace_ind_start + slice - 1) + config.extnames[field]
                )
                hdu_l2.header["CTYPE1"] = (config.extnames[field][1:], "Name of axis 1")
                hdu_l2.header["CTYPE2"] = ("Order-N", "Name of axis 2")

                # Add or update the extension in RV2
                if hdu_l2.header["EXTNAME"] not in RV2.extensions:
                    RV2.create_extension(
                        ext_name=hdu_l2.header["EXTNAME"],
                        ext_type="ImageHDU",
                        header=hdu_l2.header,
                        data=hdu_l2.data,
                    )
                else:
                    RV2.set_data(ext_name=hdu_l2.header["EXTNAME"], data=hdu_l2.data)
                    RV2.set_header(
                        ext_name=hdu_l2.header["EXTNAME"], header=hdu_l2.header
                    )


def doppler_shift(wave: np.ndarray, rv: float) -> np.ndarray:
    """
    Performs the doppler shift on the wavelength values.

    Args:
        wave (np.ndarray): the original wavelength values
        rv (float): the radial velocity of the object in km/s

    Returns:
        wave_shifted (np.ndarray): the doppler shifted wavelength values
    """

    wave_shifted = wave + wave * rv / (c / 1e3).value
    return wave_shifted


def add_fits_extension(rv2_obj: RV2, name: str, value: float) -> None:
    """
    Creates and adds a new FITS extension to an RV2 object.

    This function generates a FITS extension containing a single value and
    associates it with the given name. The extension is then added to the RV2
    object for further processing.

    Args:
        rv2_obj (RV2): The RV2 object where the extension will be stored.
        name (str): The name of the extension.
        value (float): The value to store in the FITS extension.

    Returns:
        None
    """

    # Create an ImageHDU with the specified value
    hdu = fits.ImageHDU(data=np.ones(1) * value)
    hdu.header["EXTNAME"] = name
    hdu.header["CTYPE1"] = (name, "Name of axis 1")

    # Add the header and data to the RV2 object
    rv2_obj.set_header(name, hdu.header)
    rv2_obj.set_data(name, hdu.data)


def add_nan_row(matrix: np.ndarray, row_index: int) -> np.ndarray:
    """
    Inserts a row of NaN values at a specified index in a 2D NumPy array.

    Args:
        matrix (np.ndarray): The original 2D array.
        row_index (int): The index at which the NaN row should be inserted.

    Returns:
        matrix_updated (np.ndarray): A new array with the NaN row inserted.
    """

    # Get the original dtype
    dtype = matrix.dtype

    if dtype.itemsize <= np.dtype(np.int16).itemsize:
        if dtype == np.uint16:
            # Create a row filled with '16384' values (same as bad pixel)
            nan_row = np.full((1, matrix.shape[1]), 16384, dtype=dtype)
        else:
            raise ValueError(
                "The extension data type does not support the addition of a row"
                " for the missing order. The current data type of the extension "
                "is '{0}', which is incompatible with the intended operation. "
                "Please ensure the extension's data type allows insertion of NaN"
                " values or consider using a compatible type.".format(dtype)
            )
    else:
        # Create a row filled with NaN values
        nan_row = np.full((1, matrix.shape[1]), np.nan, dtype=dtype)

    # Insert the NaN row into the array
    matrix_updated = np.insert(matrix, row_index, nan_row, axis=0)
    return matrix_updated
