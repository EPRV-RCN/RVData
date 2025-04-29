'''
RVData/rvdata/instruments/nirps/utils/validate_fits_file.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Wen Mar 07 2025
Last Modified: Wen Mar 07 2025
Version: 1.0.0
Description:
Validates a FITS file to ensure it meets the necessary criteria for conversion.

---------------------
Libraries
---------------------
'''
from astropy.io import fits

import rvdata.instruments.nirps.config.config as config


def validate_fits_file(path: str) -> None:
    """
    Validates whether a FITS file meets the required conditions for conversion.

    Args:
        path (str): The path to the FITS file.

    Raises:
        ValueError: If the FITS file does not meet the conversion criteria.
    """

    with fits.open(path) as hdu_raw:
        # Check required DPR category
        dpr_catg = hdu_raw['PRIMARY'].header['HIERARCH ESO DPR CATG']
        if dpr_catg != config.DPR_CATG_REQUIRED:
            print("Not translatable")
            raise ValueError(
                f"Error: File {path} is '{dpr_catg}' instead of 'SCIENCE'."
                " Conversion not possible."
            )

        # Check excluded objects
        object_name = hdu_raw['PRIMARY'].header['OBJECT']
        if object_name in config.EXCLUDE_OBJECTS:
            print("Not translatable")
            raise ValueError(
                f"Error: File {path} corresponds to an observation of {object_name}."
                " Conversion not possible."
            )

        # Check excluded DPR types
        dpr_type = hdu_raw['PRIMARY'].header['HIERARCH ESO DPR TYPE'].split(",")[1]
        if dpr_type in config.EXCLUDE_DPR_TYPES:
            print("Not translatable")
            raise ValueError(
                f"Error: File {path} corresponds to a '{dpr_type}' observation."
                " Conversion not possible."
            )

        # Check excluded PROGRAM
        program = hdu_raw['PRIMARY'].header['HIERARCH ESO OBS PROG ID']
        if program in config.EXCLUDE_PROGRAMS:
            print("Not translatable")
            raise ValueError(
                f"Error: File {path} corresponds to a specific program'{program}'"
                ". Conversion not possible."
            )

        print(dpr_catg, dpr_type, object_name)
