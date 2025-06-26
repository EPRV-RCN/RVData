from astropy.io import fits
# from astropy.table import Table
# import numpy as np
# import pandas as pd
# import os
# from collections import OrderedDict

# # import base class
from rvdata.core.models.level4 import RV4
# from rvdata.core.models.definitions import LEVEL4_EXTENSIONS


# KPF Level2 Reader
class NEIDRV4(RV4):
    """
    Data model and reader for RVData Level 4 (RV) data constructed from
    NEID Level 2 pipeline products.

    This class extends the `RV4` base class to handle NEID data. It 
    reads the relevant extensions from a NEID Level 2 FITS file and 
    organizes them into a standardized format.

    Parameters
    ----------
    Inherits all parameters from :class:`RV4`.

    Attributes
    ----------
    extensions : dict
        Dictionary of all created extensions (e.g., 'RV1', 'CCF1', etc.),
        mapping extension names to their data arrays.
    headers : dict
        Dictionary of headers for each extension, mapping extension names to
        their FITS headers.
    data : dict
        Dictionary of data arrays for each extension.

    Notes
    -----
    To construct an RVData Level 4 object, a NEID Level 2 FITS file is required.
    The classmethod `from_fits` should be used to instantiate the object 
    from these files. The `_read` method is not intended to be called 
    directly by users.

    Example
    -------
    >>> from rvdata.instruments.neid.level4 import NEIDRV4
    >>> neid_rv4_obj = NEIDRV4.from_fits("neidL2_YYYYMMDDTHHMMSS.fits", instrument="NEID")
    >>> neid_rv4_obj.to_fits("neid_L4_standard.fits")
    """

    def _read(self, hdul2: fits.HDUList, **kwargs) -> None:

        # Set up the primary header (take code from L2 translator)

        # RV1 - turn the CCFS extension header into a table

        # CCF1 - take the CCFS data extension, unsure about header

        # Diagnostics1 - take the activity extension table? also some stuff from the CCFS header
        