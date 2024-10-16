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
        # Output original primary header to own extension (will clean up after definitions)
        instrument_header_ext = 'INSTRUMENT_HEADER'
        if instrument_header_ext not in self.extensions.keys():
            self.create_extension(instrument_header_ext, SpectrumCollection)
        self.header[instrument_header_ext] = hdul[0].header
        
        # Check observation mode to set fiber list
        if hdul[0].header['OBS-MODE'] == 'HR':
            fiber_list = ['SCI', 'SKY', 'CAL']
        elif hdul[0].header['OBS-MODE'] == 'HE':
            fiber_list = ['SCI', 'SKY']

        for fiber in fiber_list:
            flux_ext = f'{fiber}FLUX'
            wave_ext = f'{fiber}WAVE'
            var_ext = f'{fiber}VAR'
            out_ext = f'{fiber}1'

            # Extracted flux
            flux = u.Quantity(hdul[flux_ext].data, unit=u.electron)

            # Wavelength array and associated WCS
            wave = u.Quantity(hdul[wave_ext].data, unit='AA')
            wcs = np.array([gwcs_from_array(x) for x in wave])

            # Variance array
            var = VarianceUncertainty(hdul[var_ext].data, unit=u.electron)

            # Header
            meta = hdul[flux_ext].header

            # Construct spectrum collection and output to the data object
            spec = SpectrumCollection(flux=flux, 
                                      spectral_axis=wave,
                                      uncertainty=var,
                                      wcs=wcs, 
                                      meta=meta)
            
            if out_ext not in self.extensions.keys():
                self.create_extension(out_ext, SpectrumCollection)
            setattr(self, out_ext, spec)
            self.header[out_ext] = meta

            ### For dealing with blaze, will need NEID L2 file.
            
            # blaze_ext = f'{fiber}BLAZE'
            # out_ext = f'{fiber}1_BLAZE'

            # # Blaze flux from extension, wavelength from target flux extension
            # blaze_flux = u.Quantity(hdul2[blaze_ext].data, unit=u.electron)
            # wave = u.Quantity(hdul[wave_ext].data, unit='AA')
            # wcs = np.array([gwcs_from_array(x) for x in wave])

            # # Header
            # meta = hdul2[blaze_ext].header

            # # Construct spectrum collection and output to the data object
            # spec = SpectrumCollection(flux=flux, 
            #                           spectral_axis=wave,
            #                           wcs=wcs, 
            #                           meta=meta)
            
            # if out_ext not in self.extensions.keys():
            #     self.create_extension(out_ext, SpectrumCollection)
            # setattr(self, out_ext, spec)
            # self.header[out_ext] = meta

        # Add BJD and barycentric correction extensions
        bary_kms = np.array([hdul[0].header[f'SSBRV{173-order:03d}'] for order in range(122)])
        bary_z = np.array([hdul[0].header[f'SSBZ{173-order:03d}'] for order in range(122)])
        bjd = np.array([hdul[0].header[f'SSBJD{173-order:03d}'] for order in range(122)])

        setattr(self, "BARY_KMS", bary_kms)
        setattr(self, "BARY_Z", bary_z)
        setattr(self, "BJD", bjd)

