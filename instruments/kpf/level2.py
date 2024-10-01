from astropy.io import fits
import astropy.units as u
from astropy.nddata import VarianceUncertainty
from specutils import SpectrumCollection
from specutils.utils.wcs_utils import gwcs_from_array
import numpy as np

from core.models.level2 import RV2

class KPFRV2(RV2):

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