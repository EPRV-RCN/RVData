from astropy.io import fits
import astropy.units as u
from astropy.nddata import VarianceUncertainty
from specutils import SpectrumCollection
from specutils.utils.wcs_utils import gwcs_from_array
import numpy as np

# import base class
from core.models.level2 import RV2

# NEID Level2 Reader
class NEIDRV2(RV2):
    """
    Read a NEID level 1 file and convert it to the EPRV standard format Python object.

    This class extends the `RV2` base class to handle the reading of NEID Level 2 files and 
    converts them into a standardized EPRV level 2 data format. Each extension from the FITS file 
    is read, and relevant data, including flux, wavelength, variance, and metadata, are stored as 
    attributes of the resulting Python object.

    Methods
    -------
    _read(hdul: fits.HDUList) -> None:
        Reads the input FITS HDU list, extracts specific extensions related to the science 
        data for different chips and fibers, and stores them in a standardized format.
        - The method processes data from different fibers depending on NEID OBS-MODE:
          SCI/SKY/CAL for HR mode and SCI/SKY for HE mode. (note about CAL extension in HE)
        - For each fiber, the flux, wavelength, variance, and metadata are extracted and stored as 
          a `SpectrumCollection` object.

    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions (e.g., `C1_SCI1`, `C1_SKY1`, `C1_CAL1`) 
        where the keys are the extension names and the values are `SpectrumCollection` objects 
        for each respective dataset.
    
    header : dict
        A dictionary containing metadata headers from the FITS file, with each extension's 
        metadata stored under its respective key.

    Notes
    -----
    - The `_read` method processes and organizes science and calibration data from all SCI, SKY, 
      and CAL fibers.
    - The method converts the flux, wavelength, and variance for each extension into 
      `SpectrumCollection` objects.
    
    Example
    -------
    >>> from astropy.io import fits
    >>> hdul = fits.open('neidL1_YYYYMMDDTHHMMSS.fits')
    >>> rv2_obj = NEIDRV2()
    >>> rv2_obj._read(hdul)
    """

    def _read(self, hdul: fits.HDUList) -> None:
        for fiber in ['SCI', 'SKY', 'CAL']:
            flux_ext = f'{fiber}FLUX'
            wave_ext = f'{fiber}WAVE'
            var_ext = f'{fiber}VAR'
            out_ext = f'{fiber}1'

            if fiber == 'CAL' and hdul[0].header['OBS-MODE'] == 'HE':
              flux = u.Quantity(np.array([]), unit=u.electron)
              wave = u.Quantity(np.array([]), unit='AA')
              wcs = np.array([gwcs_from_array(x) for x in wave])
              meta = hdul[flux_ext].header

              spec = SpectrumCollection(flux=flux, 
                                        spectral_axis=wave, 
                                        wcs=wcs, meta=meta)
            else:
              flux = u.Quantity(hdul[flux_ext].data, unit=u.electron)
              wave = u.Quantity(hdul[wave_ext].data, unit='AA')
              wcs = np.array([gwcs_from_array(x) for x in wave])
              var = VarianceUncertainty(hdul[var_ext].data, unit=u.electron)
              meta = hdul[flux_ext].header

              spec = SpectrumCollection(flux=flux, 
                                        spectral_axis=wave,
                                        uncertainty=var,
                                        wcs=wcs, meta=meta)
            
            if out_ext not in self.extensions.keys():
                self.create_extension(out_ext, SpectrumCollection)
            setattr(self, out_ext, spec)
            self.header[out_ext] = meta
