from collections import OrderedDict
import importlib.resources
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
try:
    config_path = importlib.resources.files("rvdata.core.models.config")
except Exception:
    config_path = None  # fallback for Sphinx/doc builds

if config_path is not None:
    LEVEL2_EXTENSIONS = pd.read_csv(config_path / "L2-extensions.csv")
    LEVEL2_PRIMARY_KEYWORDS = pd.read_csv(config_path / "L2-PRIMARY-keywords.csv")
    LEVEL3_EXTENSIONS = pd.read_csv(config_path / "L3-extensions.csv")
    LEVEL3_PRIMARY_KEYWORDS = pd.read_csv(config_path / "L3-PRIMARY-keywords.csv")
    LEVEL4_EXTENSIONS = pd.read_csv(config_path / "L4-extensions.csv")
    LEVEL4_PRIMARY_KEYWORDS = pd.read_csv(config_path / "L4-PRIMARY-keywords.csv")
else:
    LEVEL2_EXTENSIONS = None
    LEVEL2_PRIMARY_KEYWORDS = None
    LEVEL3_EXTENSIONS = None
    LEVEL3_PRIMARY_KEYWORDS = None
    LEVEL4_EXTENSIONS = None
    LEVEL4_PRIMARY_KEYWORDS = None

# Dictionary of instrument readers
INSTRUMENT_READERS = {
    "KPF": {
        "module": "rvdata.instruments.kpf.level2",
        "class": "KPFRV2",
        "method": "_read",
    },
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
    "EXPRES": {
        "module": "rvdata.instruments.expres.level2",
        "class": "EXPRESRV2",
        "method": "_read",
    },
}
