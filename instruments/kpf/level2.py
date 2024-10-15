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
    Read a KPF level 1 file and convert it to the EPRV standard format Python object.

    This class extends the `RV2` base class to handle the reading of KPF (Keck Planet Finder) 
    Level 1 files and converts them into a standardized EPRV 
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
    >>> from core.models.level2 import RV2
    >>> rv2_obj = RV2.from_fits("kpf_level1_file.fits")
    >>> rv2_obj.to_fits("standard_level2.fits")
    """

    def _read(self, hdul: fits.HDUList) -> None:
        for i in range(1,4):
            flux_array = None
            wave_array = None
            var_array = None
            for c, chip in enumerate(['GREEN', 'RED']):
                flux_ext = f'{chip}_SCI_FLUX{i}'
                wave_ext = f'{chip}_SCI_WAVE{i}'
                var_ext = f'{chip}_SCI_VAR{i}'

                if flux_array is None:
                    flux_array = hdul[flux_ext].data
                    meta = hdul[flux_ext].header
                else:
                    flux_array = np.concatenate((flux_array, hdul[flux_ext].data), axis=0)

                if wave_array is None:
                    wave_array = hdul[wave_ext].data
                else:
                    wave_array = np.concatenate((wave_array, hdul[wave_ext].data), axis=0)
                  
                if var_array is None:
                    var_array = hdul[var_ext].data
                else:
                    var_array = np.concatenate((var_array, hdul[var_ext].data), axis=0)
        
        out_ext = f'SCI1'

        flux = u.Quantity(flux_array, unit=u.electron)
        wave = u.Quantity(wave_array, unit='AA')
        wcs = np.array([gwcs_from_array(x) for x in wave])
        var = VarianceUncertainty(var_array, unit=u.electron)

        spec = SpectrumCollection(flux=flux, 
                                  spectral_axis=wave,
                                  uncertainty=var,
                                  wcs=wcs, meta=meta)

        if out_ext not in self.extensions.keys():
            self.create_extension(out_ext, SpectrumCollection)
        setattr(self, out_ext, spec)
        self.header[out_ext] = meta        
        
        for fiber in ['SKY', 'CAL']:
            flux_array = None
            wave_array = None
            var_array = None
            out_ext = f'{fiber}1'

            for c, chip in enumerate(['GREEN', 'RED']):
                flux_ext = f'{chip}_{fiber}_FLUX'
                wave_ext = f'{chip}_{fiber}_WAVE'
                var_ext = f'{chip}_{fiber}_VAR'

                if flux_array is None:
                    flux_array = hdul[flux_ext].data
                    meta = hdul[flux_ext].header
                else:
                    flux_array = np.concatenate((flux_array, hdul[flux_ext].data), axis=0)

                if wave_array is None:
                    wave_array = hdul[wave_ext].data
                else:
                    wave_array = np.concatenate((wave_array, hdul[wave_ext].data), axis=0)
                  
                if var_array is None:
                    var_array = hdul[var_ext].data
                else:
                    var_array = np.concatenate((var_array, hdul[var_ext].data), axis=0)

            flux = u.Quantity(flux_array, unit=u.electron)
            wave = u.Quantity(wave_array, unit='AA')
            wcs = np.array([gwcs_from_array(x) for x in wave])
            var = VarianceUncertainty(var_array, unit=u.electron)
            meta = hdul[flux_ext].header

            spec = SpectrumCollection(flux=flux, 
                                      spectral_axis=wave,
                                      uncertainty=var,
                                      wcs=wcs, meta=meta)
            
            if out_ext not in self.extensions.keys():
                self.create_extension(out_ext, SpectrumCollection)
            setattr(self, out_ext, spec)
            self.header[out_ext] = meta
    
        bary = hdul['BARY_CORR'].data
        bary_kms = bary['BARYVEL'] / 1000.
        setattr(self, 'BARY_KMS', bary_kms)
        setattr(self, 'BARY_Z', bary_kms / 3e5)  # aproximate!!!
        setattr(self, 'BJD', bary['PHOTON_BJD'])

        extra_headers = []
        for key in self.header.keys():
            if key not in self.extensions.keys():
                extra_headers.append(key)

        for key in extra_headers:
            del self.header[key]
