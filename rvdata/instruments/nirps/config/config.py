# Required category for conversion
DPR_CATG_REQUIRED = "SCIENCE"

# Objects that should be excluded from conversion
EXCLUDE_OBJECTS = {"SUN", "solar_spectrum", "Sun"}

# DPR types that should be excluded from conversion
EXCLUDE_DPR_TYPES = {"CIRPOL"}

# PROGRAM that should be excluded from conversion
EXCLUDE_PROGRAMS = {}

fiber = {
    "FP": {'A': 'SCI', 'B': 'FP'},
    "SKY": {'A': 'SCI', 'B': 'SKY'}
}

# Number of slices in a fiber
slice_nb = 1

# Number of pixels per order
num_pixel = 4084

# Number of orders in a fiber
NUMORDER = 71

data_format = "L2"  # Can either be original or L2

# DRS version for proper file search functionality.
DRS_VERSION = "DRS-3.2.0"

# Allows the conversion of the RAW file
extnames_raw = {
    'posemeter': {'name': 'EXPMETER', 'type': 'BinTableHDU'},
    'GUIDING': {'name': 'GUIDINGIMAGE', 'type': 'ImageHDU'}
}

extnames_telluric = {
    'SCIDATA': '_TELLURIC_FLUX',
    'ERRDATA': '_TELLURIC_VAR',
    'QUALDATA': '_TELLURIC_QUALDATA'
}

extnames_skysub = {
    'SCIDATA': '_SKYSUB_FLUX',
    'ERRDATA': '_SKYSUB_VAR',
    'QUALDATA': '_SKYSUB_QUALDATA'
}

# Allows the conversion of S2D_BLAZE files
extnames = {
    'SCIDATA': '_FLUX',
    'ERRDATA': '_VAR',
    'WAVEDATA_VAC_BARY': '_WAVE',
    'QUALDATA': '_QUALDATA',
    'DLLDATA_VAC_BARY': '_DISP'
}

# Allows the correction of TUNIT Keyword
TUNIT_FIXES = {
    "sec": "s",
    "counts": "count",
    "days": "d",
    "ADU": "adu"
}

# Define the time ranges of instrument versions
INSTRUMENT_VERSIONS = [
    {"version": 'NIRPS',
     "start_date": "2020-08-11",
     "end_date": None}  # None = Until now
]

# Parameters for Simbad queries
# Allows setting the timeout after which a query is considered failed without a response
timeout = 30
