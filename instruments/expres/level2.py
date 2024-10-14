from astropy.io import fits
from astropy.constants import c
import astropy.units as u
from astropy.nddata import VarianceUncertainty
from specutils import SpectrumCollection
from specutils.utils.wcs_utils import gwcs_from_array
import numpy as np

# import base class
from core.models.level2 import RV2

# EXPRES Level2 Reader
class EXPRESRV2(RV2):
    """
    Read EXPRES extracted file and convert it to the EPRV standard format Python object.

    This class extends the `RV2` base class to handle the reading of EXPRES
    (EXtreme PREcision Spectrograph) files and converts them into a standardized EPRV 
    format. Each extension from the FITS file is read, and relevant data, including flux, 
    wavelength, variance, and metadata, are stored as attributes of the resulting Python object.

    Methods
    -------
    _read(hdul: fits.HDUList) -> None:
        Reads the input FITS HDU list, extracts specific extensions related to the science 

    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions (e.g., `C1_SCI1`) 
        where the keys are the extension names and the values are `SpectrumCollection` objects 
        for each respective dataset.
    
    header : dict
        A dictionary containing metadata headers from the FITS file, with each extension's 
        metadata stored under its respective key.

    Notes
    -----
    - The `_read` method processes science and calibration data, 
      and it extracts and organizes the science data.
    - The method converts the flux, wavelength, and variance for each extension into 
      `SpectrumCollection` objects.
    - Unused extensions are removed from the object.
    
    Example
    -------
    >>> from astropy.io import fits
    >>> hdul = fits.open('expres_level2_file.fits')
    >>> rv2_obj = EXPRESRV2()
    >>> rv2_obj._read(hdul)
    """

    def _read(self, hdul: fits.HDUList) -> None:
        data = hdul[1].data.copy()
        head1 = hdus[0].header
        head2 = hdus[1].header

        # Wavelength
        wave_arr = data['excalibur'] #data['wavelength']
        wave = u.Quantity(wave_arr, unit='AA')
        wcs = np.array([gwcs_from_array(x) for x in wave])

        bary_arr = data['bary_excalibur'] # data['bary_wavelength']
        bary = u.Quantity(bary_arr, unit='AA')
        berv_z = u.Quantity(1-bary_arr/wave_arr,unit='z')
        berv_kms = u.Quantity(1-bary_arr/wave_arr*c.to('km/s'),unit=u.km/u.s)

        # Intensity
        blaze = data['blaze']
        flux = u.Quantity(data['spectrum']*blaze, unit=u.electron)
        flat = u.Quantity(blaze, unit=u.electron)
        #cont = u.Quantity(data['continuum'])

        # Uncertainty
        unct = VarianceUncertainty(data['uncertainty']**2, unit=u.electron)
        flat_unct = VarianceUncertainty(blaze/30, unit=u.electron) # maybe?

        # Spectrum Collection Objects
        spec = SpectrumCollection(flux=flux,
                                  spectral_axis=wave,
                                  uncertainty=unct,
                                  wcs=wcs, meta=head1)
        self.create_extension('SCI1', SpectrumCollection)
        setattr(self, 'SCI1', spec)

        blaz = SpectrumCollection(flux=flat,
                                  spectral_axis=wave,
                                  uncertainty=flat_unct,
                                  wcs=wcs, meta=head1)
        self.create_extension('BLZ1', SpectrumCollection)
        setattr(self, 'BLZ1', blaz)

        # Delete Unused Extensions
        #self.del_extension('tellurics')