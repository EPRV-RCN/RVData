# Required category for conversion
DPR_CATG_REQUIRED = "SCIENCE"

# Objects that should be excluded from conversion
EXCLUDE_OBJECTS = {"SUN", "solar_spectrum", "Sun"}

# DPR types that should be excluded from conversion
EXCLUDE_DPR_TYPES = {"CIRPOL"}

# PROGRAM that should be excluded from conversion
EXCLUDE_PROGRAMS = {}

fiber = {"FP": {"A": "SCI", "B": "FP"}, "SKY": {"A": "SCI", "B": "SKY"}}

# Number of slices in a fiber
slice_nb = 2

# Number of pixels per order
num_pixel = 9211

# Number of orders in a fiber
NUMORDER = 85

data_format = "L2"  # Can either be original or L2

# DRS version for proper file search functionality.
DRS_VERSION = "DRS-3.3.10"

# Allows the conversion of the RAW file
extnames_raw = {
    "Exp Meter bin table": {"name": "EXPMETER", "type": "BinTableHDU"},
    "FS1INT": {"name": "GUIDINGIMAGE", "type": "ImageHDU"},
    "PS1": {"name": "PUPILIMAGE", "type": "ImageHDU"},
}

# Allows the conversion of TELLURIC files
extnames_telluric = {
    "SCIDATA": "_TELLURIC_FLUX",
    "ERRDATA": "_TELLURIC_VAR",
    "QUALDATA": "_TELLURIC_QUALDATA",
}

# Allows the conversion of SKYSUB files
extnames_skysub = {
    "SCIDATA": "_SKYSUB_FLUX",
    "ERRDATA": "_SKYSUB_VAR",
    "QUALDATA": "_SKYSUB_QUALDATA",
}

# Allows the conversion of S2D_BLAZE files
extnames = {
    "SCIDATA": "_FLUX",
    "ERRDATA": "_VAR",
    "WAVEDATA_VAC_BARY": "_WAVE",
    "QUALDATA": "_QUALDATA",
    "DLLDATA_VAC_BARY": "_DISP",
}

# Allows the correction of TUNIT Keyword
TUNIT_FIXES = {"sec": "s", "counts": "count", "ADU": "adu", "days": "d"}

# Define the time ranges of instrument versions
INSTRUMENT_VERSIONS = [
    {"version": "ESPRESSO18", "start_date": "2017-11-27", "end_date": "2019-06-14"},
    {
        "version": "ESPRESSO19",
        "start_date": "2019-06-24",
        "end_date": None,  # None = Until now
    },
]

# Parameters for Simbad queries. Allows setting the timeout after which a
# query is considered failed without a response
timeout = 30
