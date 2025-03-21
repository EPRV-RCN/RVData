'''
RVData/rvdata/instruments/nirps/utils/convert_RAW.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Wen Mar 07 2025
Last Modified: Wen Mar 07 2025
Version: 1.0.0
Description:
Converts raw FITS data into an RV2 object, extracting relevant data
from specified extensions (image and binary tables, etc) and updating
the object with the necessary metadata for the EPRV project.

---------------------
Libraries
---------------------
'''
from astropy.io import fits

from rvdata.core.models.level2 import RV2
import rvdata.instruments.nirps.config.config as config


def convert_RAW(RV2: RV2, file_path: str) -> None:
    """
    Convert and integrate raw FITS data into an RV2 object.

    This function reads a FITS file, extracts the relevant data from configured
    extensions, and adds or updates them in the given RV2 object. It supports both
    image (IMAGEHDU) and binary table (BINTABLEHDU) extensions.

    Parameters:
        RV2 (RV2): The RV2 object where the converted data will be stored.
        file_path (str): Path to the FITS file to be processed.
    Returns:
        None : The function modifies the RV2 object in place by adding or
            updating the extnames_raw extensions.
    """

    with fits.open(file_path) as hdul:
        # Loop through all predefined extensions in config.extnames_raw
        for field in config.extnames_raw.keys():
            field_info = config.extnames_raw.get(field, {})  # Retrieve field details
            field_type = field_info.get('type')  # Get the extension type

            try:
                # Extract the data from the current field in the FITS file
                raw_data = hdul[field].data

                # Create the appropriate HDU type based on field type
                if field_type == 'ImageHDU':
                    raw_hdu = fits.ImageHDU(
                        data=raw_data,
                        header=hdul[field].header
                    )
                elif field_type == 'BinTableHDU':
                    if field_info.get('name') == 'EXPMETER':
                        raw_hdu = fix_tunit_keywords(hdul[field])
                    else:
                        raw_hdu = fits.BinTableHDU(
                            data=raw_data,
                            header=hdul[field].header
                        )
                else:
                    # Skip unknown types
                    continue

                # Update the header with relevant metadata
                raw_hdu.header['EXTNAME'] = field_info.get('name')
                # Check if the extension already exists in the RV2 object
                if (raw_hdu.header['EXTNAME'] not in RV2.extensions):
                    # If the extension does not exist, create it
                    RV2.create_extension(
                        ext_name=raw_hdu.header['EXTNAME'],
                        ext_type=field_type,
                        header=raw_hdu.header,
                        data=raw_hdu.data
                    )
                else:
                    # If the extension exists, update its data and header
                    RV2.set_header(raw_hdu.header['EXTNAME'], raw_hdu.header)
                    RV2.set_data(raw_hdu.header['EXTNAME'], raw_hdu.data)
            except Exception:
                # Skip if the extension is not find
                raise Exception
                continue


def fix_tunit_keywords(hdu) -> None:
    """
    Corrects non-standard TUNIT keywords in a FITS HDU header using values
    from config.py.

    Parameters:
        hdu (fits.BinTableHDU or fits.ImageHDU): The FITS HDU to process.

    Returns:
        None: The function modifies the HDU header in place.
    """
    # Retrieve the existing columns as ColDefs
    col_defs = hdu.columns

    # Create a list of new columns with corrected units
    new_cols = []
    for col in col_defs:
        new_unit = config.TUNIT_FIXES.get(col.unit, col.unit)
        new_col = fits.Column(
            name=col.name, format=col.format,
            unit=new_unit, array=hdu.data[col.name]
        )
        new_cols.append(new_col)

    # Recreate a BinTableHDU with the corrected columns
    new_hdu = fits.BinTableHDU.from_columns(new_cols, header=hdu.header)
    return new_hdu
