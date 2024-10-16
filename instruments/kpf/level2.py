from astropy.io import fits
import astropy.units as u
from astropy.nddata import VarianceUncertainty
from specutils import SpectrumCollection
from specutils.utils.wcs_utils import gwcs_from_array
import numpy as np

# import base class
from core.models.level2 import RV2

# KPF Level2 Reader
class KPFRV2(RV2):
    """
    Read a KPF level 2 file and convert it to the EPRV standard format Python object.

    This class extends the `RV2` base class to handle the reading of KPF (Keck Planet Finder) 
    Level 2 files and converts them into a standardized EPRV 
    format. Each extension from the FITS file is read, and relevant data, including flux, 
    wavelength, variance, and metadata, are stored as attributes of the resulting Python object.

    Methods
    -------
    _read(hdul: fits.HDUList) -> None:
        Reads the input FITS HDU list, extracts specific extensions related to the science 
        data for different chips and fibers, and stores them in a standardized format.

        - The method processes science data (`SCI_FLUX`, `SCI_WAVE`, `SCI_VAR`) from both 
          the GREEN and RED chips and different fibers (`SKY`, `CAL`).
        - For each chip and fiber, the flux, wavelength, variance, and metadata are extracted 
          and stored as a `SpectrumCollection` object.
        - Deletes unused extensions such as `RED_TELLURIC`, `GREEN_TELLURIC`, and `TELEMETRY`.

    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions (e.g., `C1_SCI1`, `C1_SKY1`, `C2_CAL1`) 
        where the keys are the extension names and the values are `SpectrumCollection` objects 
        for each respective dataset.
    
    header : dict
        A dictionary containing metadata headers from the FITS file, with each extension's 
        metadata stored under its respective key.

    Notes
    -----
    - The `_read` method processes science and calibration data from the GREEN and RED chips, 
      and it extracts and organizes data for both the SCI, SKY, and CAL fibers.
    - The method converts the flux, wavelength, and variance for each extension into 
      `SpectrumCollection` objects.
    - Unused extensions (like `RED_TELLURIC`, `GREEN_TELLURIC`, and `TELEMETRY`) are removed 
      from the object.
    
    Example
    -------
    >>> from astropy.io import fits
    >>> hdul = fits.open('kpf_level2_file.fits')
    >>> rv2_obj = KPFRV2()
    >>> rv2_obj._read(hdul)
    """

    def _read(self, hdul: fits.HDUList) -> None:
        for c, chip in enumerate(['GREEN', 'RED']):
            for i in range(1,4):
                flux_ext = f'{chip}_SCI_FLUX{i}'
                wave_ext = f'{chip}_SCI_WAVE{i}'
                var_ext = f'{chip}_SCI_VAR{i}'
                out_ext = f'C{str(c+1)}_SCI{i}'

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
            
            for fiber in ['SKY', 'CAL']:
                flux_ext = f'{chip}_{fiber}_FLUX'
                wave_ext = f'{chip}_{fiber}_WAVE'
                var_ext = f'{chip}_{fiber}_VAR'
                out_ext = f'C{str(c+1)}_{fiber}1'

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
    
        self.del_extension('RED_TELLURIC')
        self.del_extension('GREEN_TELLURIC')
        self.del_extension('TELEMETRY')