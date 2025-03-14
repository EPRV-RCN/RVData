'''
RVData/rvdata/instruments/harpsn/utils/validate_fits_file.py

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

import rvdata.instruments.harpsn.config.config as config


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
        dpr_catg = hdu_raw['PRIMARY'].header['HIERARCH TNG DPR CATG']
        if dpr_catg != config.DPR_CATG_REQUIRED:
            print("Not translatable")
            raise ValueError(
                f"Error: File {path} is '{dpr_catg}' instead of 'SCIENCE'."
                " Conversion not possible."
            )

        # Check excluded objects
        object_name = hdu_raw['PRIMARY'].header['TNG OBS TARG NAME']
        if object_name in config.EXCLUDE_OBJECTS:
            print("Not translatable")
            raise ValueError(
                f"Error: File {path} corresponds to an observation of "
                "{object_name}. Conversion not possible."
            )

        # Check excluded DPR types
        dpr_type = (
            hdu_raw['PRIMARY'].header['HIERARCH TNG DPR TYPE'].split(",")[1]
            )
        if dpr_type in config.EXCLUDE_DPR_TYPES:
            print("Not translatable")
            raise ValueError(
                f"Error: File {path} corresponds to a '{dpr_type}' "
                "observation. Conversion not possible."
            )

        print(dpr_catg, dpr_type, object_name)
