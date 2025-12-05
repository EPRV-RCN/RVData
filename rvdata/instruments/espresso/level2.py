"""
RVData/rvdata/instruments/espresso/level2.py

UNIGE-ESO - EPRV
Author: Loris JACQUES & Emile FONTANET
Created: Mon Mar 03 2025
Last Modified: Mon Mar 03 2025
Version: 1.0.0

---------------------
Libraries
---------------------
"""

from astropy.io import fits
import os
import numpy as np
import pandas as pd
from rvdata.core.models.level2 import RV2
import rvdata.instruments.espresso.config.config as config
from rvdata.instruments.espresso.utils import (
    convert_S2D_BLAZE,
    convert_BLAZE,
    convert_DRIFT,
    get_files_names,
    create_PRIMARY,
    validate_fits_file,
    convert_RAW,
    convert_TELLURIC,
    convert_SKYSUB,
)

# ESPRESSO Level2 Reader


class ESPRESSORV2(RV2):
    """
    Read ESPRESSO Level 1 and Level 2 files and convert them into the EPRV
    standard format.

    This class extends the `RV2` base class to handle the reading of ESPRESSO
    (Echelle SPectrograph for Rocky Exoplanets and Stable Spectroscopic
    Observations) Level 1 and Level 2 files, combining information from both
    sources to produce a standardized EPRV output. It processes various FITS
    extensions and organizes flux, wavelength, variance, and metadata into a
    structured Python object.

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
        >>> rv2_obj = ESPRESSORV2.from_fits("espresso_level1_file.fits")
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
        5. Creates the ``PRIMARY`` header and integrates necessary metadata.
        6. Cleans up unused extensions like ``RECEIPT`` and ``DRP_CONFIG``.

        Parameters
        ----------
        hdul : fits.HDUList
            The FITS HDU list to be processed.
        directory_structure : str
            Type of database architecture that stores resources. Must be either
            'dace' or 'standard'.

        Raises
        ------
        ValueError
            If the FITS file is invalid and does not meet the required criteria
            for conversion.

        :noindex:
        """

        path = os.path.join(self.dirname, self.filename)
        # Validate the FITS file before conversion. If it does not meet the
        # criteria, raise an error
        try:
            validate_fits_file(path)
            print("File is valid for conversion!")
        except ValueError as e:
            raise ValueError(e)

        # Retrieve the paths for the necessary files
        names = get_files_names(path, directory_structure)

        # Convert RAW, S2D_BLAZE_A, S2D_BLAZE_B, BLAZE_A, BLAZE_B,
        # DRIFT_B,TELLURIC and SKYSUB files
        trace_ind_start = 1

        with fits.open(path) as hdu_raw:
            dpr_type = hdu_raw["PRIMARY"].header["HIERARCH ESO DPR TYPE"].split(",")[1]
            slice_nb = config.slice_nb[
                hdu_raw["PRIMARY"].header["HIERARCH ESO INS MODE"]
            ]

        fibers = config.fiber.get(dpr_type, {})
        convert_RAW(self, path)

        for fiber in fibers:
            convert_S2D_BLAZE(
                self, names["s2d_blaze_file_" + fiber], trace_ind_start, slice_nb
            )
            convert_BLAZE(self, names["blaze_file_" + fiber], trace_ind_start, slice_nb)
            if fiber == "A":
                try:
                    convert_TELLURIC(
                        self,
                        names["telluric_file_" + fiber],
                        trace_ind_start,
                        slice_nb,
                    )
                    print("TRACEi_TELLURIC_x extensions " "have been generated.")
                except Exception:
                    print(
                        "No TELLURIC file found, TRACEi_TELLURIC_x extensions "
                        "will not be generated."
                    )

                try:
                    convert_SKYSUB(
                        self,
                        names["skysub_file_" + fiber],
                        trace_ind_start,
                        slice_nb,
                    )
                    print("TRACEi_SKYSUB_x extensions " "have been generated.")
                except Exception:
                    print(
                        "No SKYSUB file found, TRACEi_SKYSUB_x extensions "
                        "will not be generated."
                    )

            if fiber == "B":
                convert_DRIFT(
                    self, names["drift_file_" + fiber], trace_ind_start, slice_nb
                )

            trace_ind_start += slice_nb

        # Create the PRIMARY header
        nb_fiber = len(fibers)
        nb_trace = nb_fiber * slice_nb
        create_PRIMARY(self, names, nb_trace, slice_nb)

        # Filling the EXT_DESCRIPT and ORDER_TABLE extensions
        try:
            # Get the parent directory of the "utils" folder
            base_dir = os.path.dirname(os.path.realpath(__file__))

            # Properly construct the file path
            ext_descript_path = os.path.join(base_dir, "config", "ext_descript.csv")
            ext_descript_df = pd.read_csv(ext_descript_path)
            self.set_data("EXT_DESCRIPT", ext_descript_df)
        except Exception as e:
            print("Error while setting EXT_DESCRIPT data:", e)
        try:
            base_dir = os.path.dirname(os.path.realpath(__file__))

            # Properly construct the file path
            table_order_path = os.path.join(base_dir, "config", "table_order.csv")
            self.set_data("ORDER_TABLE", pd.read_csv(table_order_path))
        except Exception as e:
            print("Error while setting ORDER_TABLE data:", e)
        # Remove empty extensions
        rm_list = []
        for key, value in self.headers.items():
            if len(value) == 0:
                rm_list.append(key)
        # We do not remove the extensions for now, as it is not needed
        # for key in rm_list:
        #    self.del_extension(key)

        # Set extension description table - get rid of extensions not present
        ext_name_list = np.array(list(self.extensions))
        _, x_inds, _ = np.intersect1d(
            ext_descript_df["Name"].values, ext_name_list, return_indices=True
        )
        i_to_drop = np.setdiff1d(np.arange(ext_descript_df.shape[0]), x_inds)

        ext_descript_df.drop(i_to_drop, inplace=True)
        ext_descript_df.reset_index(inplace=True, drop=True)

        # Sort the extension description table to match the data object
        i_for_sort = []
        for name in ext_name_list:
            i_for_sort.append(np.where(name == ext_descript_df["Name"].values)[0][0])
        ext_descript_df = ext_descript_df.iloc[i_for_sort]

        self.set_data("EXT_DESCRIPT", ext_descript_df)
