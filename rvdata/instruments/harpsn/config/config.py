# Required category for conversion
DPR_CATG_REQUIRED = "SCIENCE"

# Objects that should be excluded from conversion
EXCLUDE_OBJECTS = {"SUN", "solar_spectrum", "Sun"}

# DPR types that should be excluded from conversion
EXCLUDE_DPR_TYPES = {"CIRPOL"}

# PROGRAM that should be excluded from conversion
EXCLUDE_PROGRAMS = {}

fiber = {
    "WAVE": {'A': 'SCI', 'B': 'FP'},
    "SKY": {'A': 'SCI', 'B': 'SKY'},
    "DARK": {'A': 'SCI'}
}

data_format = "L2"  # Can either be original or L2

# Number of slices in a fiber
slice_nb = 1

# Number of pixels per order
num_pixel = 4096

# Number of orders in a fiber
NUMORDER = 69

# Allows the conversion of the S2D_BLAZE files
extnames = {
    'SCIDATA': '_FLUX',
    'ERRDATA': '_VAR',
    'WAVEDATA_VAC_BARY': '_WAVE',
    'QUALDATA': '_QUALDATA',
    'DLLDATA_VAC_BARY': '_DISP'
}

# Define the time ranges of instrument versions
INSTRUMENT_VERSIONS = [
    {
        "version": 'HARPN',
        "start_date": "2012-11-10",
        "end_date": None  # None = Until now
    }
]

# Parameters for Simbad queries. Allows setting the timeout after which a
# query is considered failed without a response
timeout = 30
