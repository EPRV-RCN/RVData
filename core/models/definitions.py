from astropy.io import fits
from specutils import SpectrumCollection

import os
import pandas as pd
from collections import OrderedDict

# Header keywords required by all levels of data
# defined in a series of CSV files
LEVEL2_HEADER_FILE = os.path.abspath(os.path.dirname(__file__)) + "/headers/L2.csv"

# Base extensions for all data levels
# extensions should be defined here
# as a dictionary with the name of the extensions as keys
# and the fits data type as the values
BASE_EXTENSIONS = {
    "PRIMARY": fits.PrimaryHDU,
    "INSTRUMENT_HEADER": fits.BinTableHDU,
    "RECEIPT": fits.BinTableHDU,
    "CONFIG": fits.BinTableHDU,
}

# Minimum level 1 extensions, used in addition to BASE_EXTENSIONS
LEVEL1_EXTENSIONS = {}

# Minimum level 2 extensions, used in addtion to BASE_EXTENSIONS
LEVEL2_EXTENSIONS = {
    "SCI1": fits.ImageHDU,
    "SKY1": fits.ImageHDU,
    "CAL1": fits.ImageHDU,
    "BARY_KMS": fits.ImageHDU,
    "BARY_Z": fits.ImageHDU,
    "BJD": fits.ImageHDU,
}

# mapping between fits extension data types and Python object data types
FITS_TYPE_MAP = {
    fits.PrimaryHDU: OrderedDict,
    fits.ImageHDU: SpectrumCollection,
    fits.BinTableHDU: pd.DataFrame,
}

INSTRUMENT_READERS = {
    "KPF": {"module": "instruments.kpf.level2", "class": "KPFRV2", "method": "_read"}
}
