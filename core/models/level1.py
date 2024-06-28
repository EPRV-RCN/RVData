"""
Level 1 Data Model
"""
# Standard dependencies
import copy
import warnings

# External dependencies
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd

from core.models.base import RVDataModel
from core.models import definitions

class RV1(RVDataModel):
    """
    The level 1 RV data. Initialized with empty fields.
    Attributes inherited from RVDataModel, additional attributes below.

    """

    def __init__(self):
        super().__init__()
        self.level = 1
        extensions = copy.copy(definitions.LEVEL1_EXTENSIONS)
        python_types = copy.copy(definitions.FITS_TYPE_MAP)
        # add empty level1 extensions and empty headers for each extension
        for key, value in extensions.items():
            if key not in ['PRIMARY', 'RECEIPT', 'CONFIG']:
                if python_types[value] == np.ndarray:
                    atr = np.array([])
                else:    
                    atr = python_types[value]([])
                self.header[key] = fits.Header()
            else:
                continue
            self.create_extension(key, python_types[value])
            setattr(self, key, atr)

        # add level1 header keywords for PRIMARY header
        self.header_definitions = pd.read_csv(definitions.LEVEL1_HEADER_FILE)
        for i, row in self.header_definitions.iterrows():
            ext_name = row['Ext']
            if ext_name not in self.header.keys():
                continue
            key = row['Keyword']
            val = row['Value']
            desc = row['Description']
            if val is np.nan:
                val = None
            if desc is np.nan:
                desc = None
            self.header[ext_name][key] = (val, desc)

    def _read(self, hdul: fits.HDUList) -> None:
        '''
        Parse the HDUL based on RV standard

        Args:
            hdul (fits.HDUList): List of HDUs parsed with astropy.

        '''
        for hdu in hdul:
            if isinstance(hdu, fits.ImageHDU):
                if hdu.name not in self.extensions:
                    self.create_extension(hdu.name, np.ndarray)
                setattr(self, hdu.name, hdu.data)
            elif isinstance(hdu, fits.BinTableHDU):
                if hdu.name not in self.extensions:
                    self.create_extension(hdu.name, pd.DataFrame)
                table = Table(hdu.data).to_pandas()
                setattr(self, hdu.name, table)
            elif hdu.name != 'PRIMARY' and hdu.name != 'RECEIPT':
                warnings.warn("Unrecognized extension {} of type {}".format(hdu.name, type(hdu)))
                continue
            
            self.header[hdu.name] = hdu.header
    
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
            if isinstance(ext, (np.ndarray, np.generic)):
                row = '|{:20s} |{:20s} |{:20s}\n'.format(name, 'image',
                                                        str(ext.shape))
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
            if value == fits.PrimaryHDU:
                head = self.header[key]
                hdu = fits.PrimaryHDU(header=head)
            elif value == fits.ImageHDU:
                data = getattr(self, key)
                if data is None:
                    ndim = 0
                    # data = np.array([])
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
                except KeyError as ke:
                    print("KeyError exception raised: -->ke=" + str(ke))
                    print("Attempting to handle it...")
                    if str(ke) == '\'bool\'':
                        data = data.astype(float)
                        print("------>SHAPE=" + str(data.shape))
                        hdu = fits.ImageHDU(data=data, header=head)
                    else:
                        raise KeyError("A different error...")
            elif value == fits.BinTableHDU:
                table = Table.from_pandas(getattr(self, key))
                self.header[key]['NAXIS1'] = len(table)
                head = self.header[key]
                hdu = fits.BinTableHDU(data=table, header=head)
            else:
                print("Can't translate {} into a valid FITS format."\
                      .format(type(getattr(self, key))))
                continue
            hdu.name = key
            if hdu.name == 'PRIMARY':
                hdu_list.insert(0, hdu)
            else:
                hdu_list.append(hdu)

        return hdu_list
    
if __name__ == "__main__":
    pass
