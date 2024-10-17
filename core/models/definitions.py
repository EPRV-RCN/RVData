from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

# Extensions (FITS HDUs) should be defined here.
# Definition is in the form of a list of dicts.
# Columns headers for tables are a list of strings.

# Base extensions for all data levels
BASE_EXTENSIONS = [
    {
        "name": "PRIMARY",
        "fits_type": fits.PrimaryHDU,
        "py_type": OrderedDict,
        "columns": [],
    },
    {
        "name": "INSTRUMENT_HEADER",
        "fits_type": fits.BinTableHDU,
        "py_type": pd.DataFrame,
        "columns": [],
    },
    {
        "name": "RECEIPT",
        "fits_type": fits.BinTableHDU,
        "py_type": pd.DataFrame,
        "columns": [],
    },
    {
        "name": "DRP_CONFIG",
        "fits_type": fits.BinTableHDU,
        "py_type": pd.DataFrame,
        "columns": [],
    },
]

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
    "PrimaryHDU": OrderedDict,
    "ImageHDU": np.array,
    "BinTableHDU": pd.DataFrame,
}

# Header keywords required by all levels of data
# defined in a series of CSV files
config_path = Path("core/models/config")

LEVEL2_EXTENSIONS = pd.read_csv(config_path / "L2-extensions.csv")

# LEVEL2_HEADER_FILE = os.path.abspath(os.path.dirname(__file__)) + "/headers/L2.csv"

INSTRUMENT_READERS = {
    "KPF": {"module": "instruments.kpf.level2", "class": "KPFRV2", "method": "_read"}
}


def read_csv(filename):
    dataframe = pd.read_csv(filename)
