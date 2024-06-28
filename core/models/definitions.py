
from astropy.io import fits

import numpy as np
import os
import pandas as pd
from collections import OrderedDict

## Header keywords required by all levels of data
# defined in a series of CSV files
LEVEL2_HEADER_FILE = os.path.abspath(os.path.dirname(__file__)) + '/headers/L2.csv'

# Minimum level 1 extensions should be defined here
# as a dictionary with the name of the extensions as keys
# and the fits data type as the values
LEVEL2_EXTENSIONS = {'PRIMARY': fits.PrimaryHDU,
                     'RECEIPT': fits.BinTableHDU,
                     'CONFIG': fits.BinTableHDU,
                     'TELEMETRY': fits.BinTableHDU,

                     'SCI_FLUX': fits.ImageHDU,
                     'SKY_FLUX': fits.ImageHDU,
                     'CAL_FLUX': fits.ImageHDU,
                     'SCI_VAR': fits.ImageHDU,
                     'SKY_VAR': fits.ImageHDU,
                     'CAL_VAR': fits.ImageHDU,
                     'SCI_WAVE': fits.ImageHDU,
                     'SKY_WAVE': fits.ImageHDU,
                     'CAL_WAVE': fits.ImageHDU,
                     'TELLURIC': fits.BinTableHDU,
                     
                     'BARY_CORR': fits.BinTableHDU
                    }

# mapping between fits extension data types and Python object data types
FITS_TYPE_MAP = {fits.PrimaryHDU: OrderedDict,
                 fits.ImageHDU: np.array,
                 fits.BinTableHDU: pd.DataFrame}
