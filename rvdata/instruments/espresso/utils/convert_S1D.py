"""
RVData/rvdata/instruments/espresso/utils/convert_S2D_BLAZE.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Mon Mar 03 2025
Last Modified: Mon Mar 03 2025
Version: 1.0.0
Description:
Extracts and processes data from an S2D_BLAZE FITS file. Stores key
calibration values in an `RV3` object, applies a Doppler shift correction
if needed, and organizes data into FITS extensions.

---------------------
Libraries
---------------------
"""

from astropy.io import fits
from astropy.table import Table as AstropyTable, vstack
from rvdata.core.models.level3 import RV3
import rvdata.instruments.espresso.config.config as config
from rvdata.core.models.definitions import LEVEL3_EXTENSIONS
import os


ext_descript = {
        'STITCHED_TELLCORR_SCI_FLUX' : 'Order stitched blaze and telluric corrected flux co-added across all science traces',
        'STITCHED_TELLCORR_SCI_WAVE': 'Order stitched barycentric- and drift-corrected wavelength solution for STITCHED_TELLCORR_SCI_FLUX',
        'STITCHED_TELLCORR_SCI_VAR' : 'Order stitched variance for STITCHED_TELLCORR_SCI_FLUX',
        'STITCHED_SKYSUB_SCI_FLUX' : 'Order stitched blaze-corrected and sky-subtracted flux co-added across all science traces',
        'STITCHED_SKYSUB_SCI_WAVE': 'Order stitched barycentric- and drift-corrected wavelength solution for STITCHED_SKYSUB_SCI_FLUX',
        'STITCHED_SKYSUB_SCI_VAR' : 'Order stitched variance for STITCHED_SKYSUB_SCI_FLUX'
    }


def add_to_ext_descript(RV4, ext_name, description):
    """Add a row to the EXT_DESCRIPT table for the given extension."""
    row = AstropyTable({"Name": [ext_name], "Description": [description]})
    RV4.data["EXT_DESCRIPT"] = vstack([RV4.data["EXT_DESCRIPT"], row])


def convert_S1D(
    RV3: RV3, file_names: dict,
) -> None:
    """
    Extracts and processes relevant data from an ESPRESSO S1D FITS file.

    This function reads the specified FITS file, extracts relevant
    information, and stores the processed data in the RV3 object. The function
    ensures that certain header values are only set once during the first
    iteration. It also applies a Doppler shift correction when necessary and
    organizes the extracted data into new FITS extensions.

    Args:
        RV3 (RV3): The RV3 object where extracted data will be stored.
        file_names (dict): A dictionary containing paths to the relevant FITS files.

    Returns:
        None
    """
    
    RV3.headers['PRIMARY']["BLZCORR"] = (True, 'Has blaze been removed?')
    RV3.headers['PRIMARY']["LMPCORR"] = (True, 'Has lamp SED been removed?')
    RV3.headers['PRIMARY']["SEDCORR"] = (False, 'Has SED been removed?')
    for file in ['s1d_A', 's1d_B', 's1d_skysub_A', 's1d_tell_corr_A']:
        if not file_names[file] or not os.path.exists(file_names[file]):
            print(f"File {file_names[file]} does not exist.")
            continue
        with fits.open(file_names[file]) as hdul:
            data = hdul[1].data
            for col in config.extnames_s1d.keys():
                for extname in config.extnames_s1d[col]:
                    ext_type = LEVEL3_EXTENSIONS.query(f"Name == '{extname}'")["DataType"].values[0]
                    if ("_B" in file):
                        continue
                    corr = 'SKYSUBCORR' if 'skysub' in file else ('TELLCORR' if 'tell_corr' in file else 'CORR')
                    extname = extname.replace('CORR', corr)
                    if extname not in RV3.extensions:
                        # If the extension does not exist, create it
                        RV3.create_extension(
                            ext_name=extname,
                            ext_type=ext_type,
                            data=data[col],
                        )
                    else:
                        # If the extension exists, update its data and header
                        RV3.set_data(extname, data[col])
                    if extname not in RV3.data['EXT_DESCRIPT']['Name']:
                        add_to_ext_descript(RV3, extname, ext_descript[extname])
