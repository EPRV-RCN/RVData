from collections import OrderedDict
from pathlib import Path

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
config_path = Path("core/models/config")

LEVEL2_EXTENSIONS = pd.read_csv(config_path / "L2-extensions.csv")
LEVEL2_PRIMARY_KEYWORDS = pd.read_csv(config_path / "L2-PRIMARY-keywords.csv")

# Dictionary of instrument readers
INSTRUMENT_READERS = {
    "KPF": {"module": "instruments.kpf.level2", "class": "KPFRV2", "method": "_read"}
}
