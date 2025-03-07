'''
RVData/instruments/harps/level2.py

UNIGE-ESO - EPRV
Author: Loris JACQUES
Created: Tue Jan 07 2025
Last Modified: Tue Jan 07 2025
Version: 1.0.0
Description: 
'''

'''
---------------------
external libraries
---------------------
'''
from astropy.io import fits
import os

'''
---------------------
internal libraries
---------------------
'''
from core.models.level2 import RV2
import instruments.harps.config.config as config
from instruments.harps.utils import convert_S2D_BLAZE, convert_BLAZE, convert_DRIFT, get_files_names, create_PRIMARY, validate_fits_file

# HARPS Level2 Reader
class HARPSRV2(RV2):
    """
    Read HARPS Level 1 and Level 2 files and convert them into the EPRV standard format.

    This class extends the `RV2` base class to handle the reading of HARPS (High Accuracy 
    Radial velocity Planet Searcher) Level 1 and Level 2 files, combining information from 
    both sources to produce a standardized EPRV output. It processes various FITS extensions 
    and organizes flux, wavelength, variance, and metadata into a structured Python object.

    Methods
    -------
    do_conversion(hdul: fits.HDUList) -> None:
        Reads the input FITS HDU list, extracts specific extensions related to the science
        data for different chips and fibers, and stores them in a standardized format.

        - The method validates the FITS file before conversion to ensure it meets the 
          required criteria.
        - Retrieves necessary file paths for additional files required in the processing.
        - Converts the spectral blaze functions (`S2D_BLAZE_A`, `S2D_BLAZE_B`, `BLAZE_A`, `BLAZE_B`)
          for the different fibers.
        - Processes and converts the drift file for instrumental calibration.
        - Creates the `PRIMARY` header and necessary metadata.
        - For now : Removes unused or redundant extensions such as `RECEIPT` and `DRP_CONFIG`.

    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions, where the keys are extension 
        names, and the values are the respective data arrays.

    header : dict
        A dictionary containing metadata headers from the FITS files, with each extension's 
        metadata stored under its respective key.

    Notes
    -----
    - The `do_conversion` method processes and extracts science and calibration data.
    - The method ensures the FITS file meets the required criteria before conversion.
    - Blaze correction functions are processed and stored for each fiber.
    - The drift file is processed separately for calibration.
    - Unused extensions (like `RECEIPT` and `DRP_CONFIG`) are removed from the final output.

    Example
    -------
    >>> from core.models.level2 import RV2
    >>> rv2_obj = HARPSRV2.from_fits("harps_level1_file.fits")
    >>> rv2_obj.to_fits("standard_level2.fits")
    """
    

    def do_conversion(self, hdul: fits.HDUList) -> None:
        """
        Converts FITS files based on certain conditions and configurations.

        This method performs several processing steps:
        1. Validates the FITS file structure before conversion.
        2. Retrieves paths for required additional files (e.g., blaze functions, drift corrections).
        3. Converts and stores spectral blaze functions for different fibers.
        4. Converts the drift correction data.
        5. Creates the `PRIMARY` header and integrates necessary metadata.
        6. Cleans up unused extensions like `RECEIPT` and `DRP_CONFIG`.

        Args:
            hdul (fits.HDUList): The FITS HDU list to be processed.
        
        Raises:
            ValueError: If the FITS file is invalid and does not meet the required criteria for conversion.
        """

        path = os.path.join(self.dirname, self.filename)

        # Validate the FITS file before conversion. If it does not meet the criteria, raise an error
        try :
            validate_fits_file(path)
            print("File is valid for conversion!")
        except ValueError as e:
            raise ValueError(e)

        # Retrieve the paths for the necessary files
        names = get_files_names(path)

        # Convert S2D_BLAZE_A, S2D_BLAZE_B, BLAZE_A, and BLAZE_B files
        trace_ind_start = 1

        with fits.open(path) as hdu_raw:
            dpr_type = hdu_raw['PRIMARY'].header['HIERARCH ESO DPR TYPE'].split(",")[1]
        fibers = config.fiber.get(dpr_type, {})

        for fiber in fibers:
            convert_S2D_BLAZE(self, names["s2d_blaze_file_"+fiber], trace_ind_start, config.slice_nb)
            convert_BLAZE(self, names["blaze_file_"+fiber], trace_ind_start, config.slice_nb)
            trace_ind_start+=config.slice_nb

        # Convert the Drift file
        convert_DRIFT(self, names["drift_file_B"])

        # Create the PRIMARY header
        nb_fiber = len(fibers)
        nb_trace = nb_fiber * config.slice_nb
        create_PRIMARY(self, names, nb_trace, config.slice_nb)

        # Remove empty extensions
        self.del_extension('RECEIPT')
        self.del_extension('DRP_CONFIG')
        print('end')

