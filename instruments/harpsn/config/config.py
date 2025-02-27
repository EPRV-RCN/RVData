fiber = {
    "WAVE": {'A': 'SCI', 'B': 'FP'},
    "SKY": {'A': 'SCI', 'B': 'SKY'},
    "DARK": {'A': 'SCI'}
}


data_format = "L2"#Can either be original or L2

# Nombre de slice dans une fibre 
slice_nb = 1

# Nombre de pixels par ordre
num_pixel = 4096

# Nombre d'ordre' dans une fibre 
NUMORDER = 69

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
    {"version": 'HARPN', "start_date": "2012-11-10", "end_date": None}  # None = Jusqu'à maintenant
]

# Parametres des requetes Simbad
# Permet de régler le temps avec lequel une requete est considérée comme échoué sans réponse.
timeout = 30