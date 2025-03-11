import pandas as pd
from astropy.io import fits

def _test_l2_extensions(inpfile):
    reference_extensions = pd.read_csv('core/models/config/L2-extensions.csv')
    extdf = reference_extensions
    hdul = fits.open(inpfile)
    for i, row in extdf.iterrows():
        ext = row['Name']
        req = row['Required']
        if req:
            assert ext in hdul, f"Extension {ext} not found in {inpfile}"
    
    hdul.close()

if __name__ == '__main__':
    _test_l2_extensions('rvstandard.fits')

