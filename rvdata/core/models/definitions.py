from collections import OrderedDict
import importlib

import numpy as np
import pandas as pd

# mapping between fits extension data types and Python object data types
FITS_TYPE_MAP = {
    "PrimaryHDU": OrderedDict,
    "ImageHDU": np.array,
    "BinTableHDU": pd.DataFrame,
}

# Header keywords required by all levels of data are defined in a series
# of CSV files in this directory.
config_path = importlib.resources.files(__package__) / "config"

LEVEL2_EXTENSIONS = pd.read_csv(config_path / "L2-extensions.csv")
LEVEL2_PRIMARY_KEYWORDS = pd.read_csv(config_path / "L2-PRIMARY-keywords.csv")

# Dictionary of instrument readers
INSTRUMENT_READERS = {
    "KPF": {"module": "rvdata.instruments.kpf.level2", "class": "KPFRV2", "method": "_read"},
    "ESPRESSO": {
        "module": "rvdata.instruments.espresso.level2",
        "class": "ESPRESSORV2",
        "method": "do_conversion",
    },
    "HARPS": {
        "module": "rvdata.instruments.harps.level2",
        "class": "HARPSRV2",
        "method": "do_conversion",
    },
    "HARPSN": {
        "module": "rvdata.instruments.harpsn.level2",
        "class": "HARPSNRV2",
        "method": "do_conversion",
    },
    "NEID": {
        "module": "rvdata.instruments.neid.level2",
        "class": "NEIDRV2",
        "method": "_read",
    },
}
