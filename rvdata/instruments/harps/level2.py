"""
RVData/rvdata/instruments/harps/level2.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Tue Jan 07 2025
Last Modified: Tue Jan 07 2025
Version: 1.0.0

---------------------
Libraries
---------------------
"""

import os
from astropy.io import fits
import pandas as pd
import numpy as np
import rvdata.instruments.harps.config.config as config
from rvdata.instruments.harps.utils import (
    convert_S2D_BLAZE,
    convert_BLAZE,
    convert_DRIFT,
    get_files_names,
    create_PRIMARY,
    validate_fits_file,
)
from rvdata.core.models.level2 import RV2


# HARPS Level2 Reader


class HARPSRV2(RV2):
    """
    Read HARPS Level 1 and Level 2 files and convert them into the EPRV
    standard format.

    This class extends the `RV2` base class to handle the reading of HARPS
    (High Accuracy Radial velocity Planet Searcher) Level 1 and Level 2 files,
    combining information from both sources to produce a standardized EPRV
    output. It processes various FITS extensions and organizes flux,
    wavelength, variance, and metadata into a structured Python object.

    Methods
    -------
    do_conversion(hdul: fits.HDUList) -> None
        Reads the input FITS HDU list, extracts specific extensions related to
        the science data for different chips and fibers, and stores them in a
        standardized format.

        - The method validates the FITS file before conversion to ensure it
          meets the required criteria.
        - Retrieves necessary file paths for additional files required in the
          processing.
        - Converts the spectral blaze functions (``S2D_BLAZE_A``,
          ``S2D_BLAZE_B``, ``BLAZE_A``, ``BLAZE_B``) for the different fibers.
        - Processes and converts the drift file for instrumental calibration.
        - Creates the ``PRIMARY`` header and necessary metadata.
        - For now: Removes unused or redundant extensions such as ``RECEIPT``
          and ``DRP_CONFIG``.

    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions, where the keys are
        extension names, and the values are the respective data arrays.

    header : dict
        A dictionary containing metadata headers from the FITS files, with
        each extension's metadata stored under its respective key.

    Notes
    -----
    - The ``do_conversion`` method processes and extracts science and
      calibration data.
    - The method ensures the FITS file meets the required criteria before
      conversion.
    - Blaze correction functions are processed and stored for each fiber.
    - The drift file is processed separately for calibration.
    - Unused extensions (like ``RECEIPT`` and ``DRP_CONFIG``) are removed from the
      final output.

    Example
    -------
    >>> from core.models.level2 import RV2
    >>> rv2_obj = HARPSRV2.from_fits("harps_level1_file.fits")
    >>> rv2_obj.to_fits("standard_level2.fits")
    """

    def do_conversion(
        self, hdul: fits.HDUList, directory_structure: str = "standard"
    ) -> None:
        """
        Converts FITS files based on certain conditions and configurations.

        This method performs several processing steps:

        1. Validates the FITS file structure before conversion.
        2. Retrieves paths for required additional files (e.g., blaze
           functions, drift corrections).
        3. Converts and stores spectral blaze functions for different fibers.
        4. Converts the drift correction data.
        5. Creates the `PRIMARY` header and integrates necessary metadata.
        6. Cleans up unused extensions like `RECEIPT` and `DRP_CONFIG`.

        Args:
            hdul (fits.HDUList): The FITS HDU list to be processed.
            directory_structure (str): Type of database architecture that
            stores resources. Must be either 'dace' or 'standard'.

        Raises:
            ValueError: If the FITS file is invalid and does not meet the
                required criteria for conversion.

        """

        path = os.path.join(self.dirname, self.filename)

        # Validate the FITS file before conversion.
        # If it does not meet the criteria, raise an error
        try:
            validate_fits_file(path)
            print("File is valid for conversion!")
        except ValueError as e:
            raise ValueError(e)

        # Retrieve the paths for the necessary files
        names = get_files_names(path, directory_structure)

        # Convert S2D_BLAZE_A, S2D_BLAZE_B, BLAZE_A, and BLAZE_B files
        trace_ind_start = 1

        with fits.open(path) as hdu_raw:
            dpr_type = hdu_raw["PRIMARY"].header["HIERARCH ESO DPR TYPE"].split(",")[1]
        fibers = config.fiber.get(dpr_type, {})

        for fiber in fibers:
            convert_S2D_BLAZE(
                self, names["s2d_blaze_file_" + fiber], trace_ind_start, config.slice_nb
            )
            convert_BLAZE(
                self, names["blaze_file_" + fiber], trace_ind_start, config.slice_nb
            )
            trace_ind_start += config.slice_nb

        # Convert the Drift file
        convert_DRIFT(self, names["drift_file_B"])

        # Create the PRIMARY header
        nb_fiber = len(fibers)
        nb_trace = nb_fiber * config.slice_nb
        create_PRIMARY(self, names, nb_trace, config.slice_nb)
        # Filling the EXT_DESCRIPT and ORDER_TABLE extensions
        try:
            # Get the parent directory of the "utils" folder
            base_dir = os.path.dirname(os.path.realpath(__file__))

            # Properly construct the file path
            ext_descript_path = os.path.join(
                base_dir, "config", "ext_descript.csv"
            )
            ext_descript_df = pd.read_csv(ext_descript_path)
            self.set_data('EXT_DESCRIPT', ext_descript_df)
        except Exception as e:
            print('Error while setting EXT_DESCRIPT data:', e)

        try:
            base_dir = os.path.dirname(os.path.realpath(__file__))

            # Properly construct the file path
            table_order_path = os.path.join(base_dir, "config", "table_order.csv")
            table_order_df = pd.read_csv(table_order_path, sep='\t')[['physical_order ', 'start_wav(nm) ', 'end_wav(nm) ']]
            table_order_df['index_order'] = np.linspace(0, 71, 72, dtype='int')[::-1]
            self.set_data('ORDER_TABLE', table_order_df)
        except Exception as e:
            print('Error while setting ORDER_TABLE data:', e)
        # Remove empty extensions
        rm_list = []
        for key, value in self.headers.items():
            if len(value) == 0:
                print(f"Removing empty extension: {key}")
                rm_list.append(key)

        # for key in rm_list:
        #    self.del_extension(key)
