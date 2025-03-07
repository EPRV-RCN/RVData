# Required category for conversion
DPR_CATG_REQUIRED = "SCIENCE"

# Objects that should be excluded from conversion
EXCLUDE_OBJECTS = {"SUN", "solar_spectrum", "Sun"}

# DPR types that should be excluded from conversion
EXCLUDE_DPR_TYPES = {"CIRPOL"}

fiber = {
    "FP": {'A': 'SCI', 'B': 'FP'},
    "SKY": {'A': 'SCI', 'B': 'SKY'}
}

# Nombre de slice dans une fibre 
slice_nb = 2

# Nombre de pixels par ordre
num_pixel = 9211

# Nombre d'ordre' dans une fibre 
NUMORDER = 85

data_format = "L2"#Can either be original or L2
slices= [0,1]

# Permet la convertion du fichier RAW
extnames_raw = {
    'Exp Meter bin table': {'name': 'EXPMETER', 'type': 'BinTableHDU'},
    'FS1INT': {'name': 'PUPILIMAGE', 'type': 'ImageHDU'},
    'PS1': {'name': 'GUIDINGIMAGE', 'type': 'ImageHDU'}
}

# Permet la convertion des fichiers S2D_BLAZE
extnames = {
    'SCIDATA' : '_FLUX',
    'ERRDATA' : '_VAR',
    'WAVEDATA_VAC_BARY' : '_WAVE',
    'QUALDATA' : '_QUALDATA',
    'DLLDATA_VAC_BARY' : '_DISP'
}

# Définir les plages temporelles des versions d'instruments
INSTRUMENT_VERSIONS = [
    {"version": 'ESPRESSO18', "start_date": "2017-11-27", "end_date": "2019-06-12"},
    {"version": 'ESPRESSO19', "start_date": "2019-06-24", "end_date": None}  # None = Jusqu'à maintenant
]

# Parametres des requetes Simbad
# Permet de régler le temps avec lequel une requete est considérée comme échoué sans réponse.
timeout = 30