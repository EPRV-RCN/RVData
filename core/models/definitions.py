
from astropy.io import fits
from specutils import SpectrumCollection

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

                     'SCI1': fits.ImageHDU,
                     'SKY1': fits.ImageHDU,
                     'CAL1': fits.ImageHDU,
                     
                     'BARY_CORR': fits.BinTableHDU
                    }

# mapping between fits extension data types and Python object data types
FITS_TYPE_MAP = {fits.PrimaryHDU: OrderedDict,
                 fits.ImageHDU: SpectrumCollection,
                 fits.BinTableHDU: pd.DataFrame}

INSTRUMENT_READERS = {'KPF': 
                      {'module': 'instruments.kpf.level2', 'class': 'KPFRV2', 'method': '_read'}
                     }
