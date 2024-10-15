"""
Level 2 Data Model for RV spectral data
"""
# Standard dependencies
import copy
import warnings

# External dependencies
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from specutils import SpectrumCollection
from astropy.nddata import VarianceUncertainty
from specutils.utils.wcs_utils import gwcs_from_array
import numpy as np
import pandas as pd

import core.models.base
from core.models import definitions
from core.tools.headers import to_ascii_safe

class RV2(core.models.base.RVDataModel):
    """
    The level 2 RV data. Initialized with empty fields.
    Attributes inherited from RVDataModel, additional attributes below.

    """

    def __init__(self):
        super().__init__()
        self.level = 2
        extensions = copy.copy(definitions.LEVEL2_EXTENSIONS)
        python_types = copy.copy(definitions.FITS_TYPE_MAP)
        # add empty level2 extensions and empty headers for each extension
        for key, value in extensions.items():
            if key not in ['PRIMARY', 'RECEIPT', 'CONFIG']:
                if python_types[value] == SpectrumCollection:
                    atr = np.array([])
                else:    
                    atr = python_types[value]([])
                self.header[key] = fits.Header()
            else:
                continue
            self.create_extension(key, python_types[value])
            setattr(self, key, atr)

        # add level2 header keywords for PRIMARY header
        self.header_definitions = pd.read_csv(definitions.LEVEL2_HEADER_FILE)
        for i, row in self.header_definitions.iterrows():
            # ext_name = row['Ext']
            ext_name = 'PRIMARY'
            if ext_name not in self.header.keys():
                continue
            key = row['Keyword']
            if key is np.nan:
                continue
            val = to_ascii_safe(str(row['Value']))
            desc = to_ascii_safe(str(row['Description']))

            if val is np.nan:
                val = None
            if desc is np.nan:
                desc = None

            self.header[ext_name][key] = (val, desc)

    def _read(self, hdul: fits.HDUList) -> None:
        extension_names = [hdu.name for hdu in hdul]

        for i in range(1,10):  # limit 10 science tracesa
            flux_ext = f'SCI{i}_FLUX'
            wave_ext = f'SCI{i}_WAVE'
            var_ext = f'SCI{i}_VAR'
            out_ext = f'SCI{i}'

            if flux_ext not in extension_names:
                continue

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
            for i in range(1,10):  # limit 10 sky or cal fibers
                flux_ext = f'{fiber}{i}_FLUX'
                wave_ext = f'{fiber}{i}_WAVE'
                var_ext = f'{fiber}{i}_VAR'
                out_ext = f'{fiber}1'

                if flux_ext not in extension_names:
                    continue

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

        setattr(self, 'BARY_KMS', hdul['BARY_KMS'].data)
        self.header['BARY_KMS'] = hdul['BARY_KMS'].header
        setattr(self, 'BARY_Z', hdul['BARY_Z'].data)
        self.header['BARY_Z'] = hdul['BARY_Z'].header
        setattr(self, 'BJD', hdul['BJD'].data)
        self.header['BJD'] = hdul['BJD'].header
    
    def info(self):
        '''
        Pretty print information about this data to stdout 
        '''
        if self.filename is not None:
            print('File name: {}'.format(self.filename))
        else: 
            print('Empty {:s} Data product'.format(self.__class__.__name__))
        # a typical command window is 80 in length
        head_key = '|{:20s} |{:20s} \n{:40}'.format(
            'Header Name', '# Cards',
            '='*80 + '\n'
        )

        for key, value in self.header.items():
            row = '|{:20s} |{:20} \n'.format(key, len(value))
            head_key += row
        print(head_key)
        head = '|{:20s} |{:20s} |{:20s} \n{:40}'.format(
            'Extension Name', 'Data Type', 'Data Dimension',
            '='*80 + '\n'
        )

        for name in self.extensions.keys():
            if name == 'PRIMARY':
                continue
            
            ext = getattr(self, name)
            if isinstance(ext, SpectrumCollection):
                row = '|{:20s} |{:20s} |{:20s}\n'.format(name, 'spectrum',
                                                        str(ext.spectral_axis.shape))
                head += row
            elif isinstance(ext, pd.DataFrame):
                row = '|{:20s} |{:20s} |{:20s}\n'.format(name, 'table',
                                                        str(len(ext)))
                head += row
        print(head)
        
    def _create_hdul(self):
        '''
        Create an hdul in FITS format. 
        This is used by the base model for writing data context to file
        '''
        hdu_list = []
        hdu_definitions = self.extensions.items()
        for key, value in hdu_definitions:
            hduname = key
            if value == fits.PrimaryHDU:
                head = self.header[key]
                hdu = fits.PrimaryHDU(header=head)
                hdu_list.insert(0, hdu)
            elif value == fits.ImageHDU:
                data = getattr(self, key)
                if isinstance(data, SpectrumCollection):
                    flux = np.array(getattr(self, key).flux)
                    wave = np.array(getattr(self, key).spectral_axis)
                    var = getattr(self, key).uncertainty.array

                    for name, data in zip(['FLUX', 'WAVE', 'VAR'], [flux, wave, var]):
                        ndim = len(data.shape)
                        self.header[key]['NAXIS'] = ndim
                        if ndim == 0:
                            self.header[key]['NAXIS1'] = 0
                        else:
                            for d in range(ndim):
                                self.header[key]['NAXIS{}'.format(d+1)] = data.shape[d]
                        head = self.header[key]
                        hdu = fits.ImageHDU(data=data, header=head)
                        hduname = key+'_'+name
                        hdu.name = hduname
                        hdu_list.append(hdu)
                else:
                    if data is None:
                        ndim = 0
                    else:
                        ndim = len(data.shape)
                    self.header[key]['NAXIS'] = ndim
                    if ndim == 0:
                        self.header[key]['NAXIS1'] = 0
                    else:
                        for d in range(ndim):
                            self.header[key]['NAXIS{}'.format(d+1)] = data.shape[d]
                    head = self.header[key]
                    try:
                        hdu = fits.ImageHDU(data=data, header=head)
                        hdu.name = hduname
                        hdu_list.append(hdu)
                    except KeyError as ke:
                        print("KeyError exception raised: -->ke=" + str(ke))
                        print("Attempting to handle it...")
                        if str(ke) == '\'bool\'':
                            data = data.astype(float)
                            print("------>SHAPE=" + str(data.shape))
                            hdu = fits.ImageHDU(data=data, header=head)
                            hdu_list.append(hdu)
                        else:
                            raise KeyError("A different error...")
            elif value == fits.BinTableHDU:
                table = Table.from_pandas(getattr(self, key))
                self.header[key]['NAXIS1'] = len(table)
                head = self.header[key]
                hdu = fits.BinTableHDU(data=table, header=head)
                hdu.name = hduname
                hdu_list.append(hdu)
            else:
                print("Can't translate {} into a valid FITS format."\
                      .format(type(getattr(self, key))))
                continue

        return hdu_list
