"""
Level 2 Data Model for RV spectral data
"""

# External dependencies
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd

import core.models.base
from core.models.definitions import LEVEL2_EXTENSIONS


class RV2(core.models.base.RVDataModel):
    """
    The level 2 RV data. Initialized with empty fields.
    Attributes inherited from RVDataModel, additional attributes below.

    """

    def __init__(self):
        super().__init__()
        self.level = 2

        # TODO: initialize header keywords for each extension

        for i, row in LEVEL2_EXTENSIONS.iterrows():
            if row["Required"]:
                # TODO: set description and comment
                self.create_extension(row["Name"], row["DataType"])

    def _read(self, hdul: fits.HDUList) -> None:
        l2_ext = LEVEL2_EXTENSIONS.set_index("Name")

        for hdu in hdul:
            if 'TRACE' in hdu.name:
                t1 = 'TRACE1_' + hdu.name.split('_')[1]
                fits_type = l2_ext.loc[t1]["DataType"]
            else:
                fits_type = l2_ext.loc[hdu.name]["DataType"]
            if hdu.name not in self.extensions.keys():
                self.create_extension(hdu.name, fits_type)

            if hdu.name == "PRIMARY":
                pass
            elif fits_type == "ImageHDU":
                data = np.array(hdu.data)
                self.set_data(hdu.name, data)
            elif fits_type == "BinTableHDU":
                data = Table(hdu.data).to_pandas()
                self.set_data(hdu.name, data)

            self.set_header(hdu.name, hdu.header)

    def info(self):
        """
        Pretty print information about this data to stdout
        """
        if self.filename is not None:
            print("File name: {}".format(self.filename))
        else:
            print("Empty {:s} Data product".format(self.__class__.__name__))
        # a typical command window is 80 in length
        head_key = "|{:20s} |{:20s} \n{:40}".format(
            "Header Name", "# Cards", "=" * 80 + "\n"
        )

        for key, value in self.headers.items():
            if value is None:
                length = 0
            else:
                length = len(value)
            row = "|{:20s} |{:20} \n".format(key, length)
            head_key += row
        print(head_key)
        head = "|{:20s} |{:20s} |{:20s} \n{:40}".format(
            "Extension Name", "Data Type", "Data Dimension", "=" * 80 + "\n"
        )

        for name in self.extensions.keys():
            if name == "PRIMARY":
                continue

            ext = self.data[name]
            if isinstance(ext, np.ndarray):
                row = "|{:20s} |{:20s} |{:20s}\n".format(name, "array", str(ext.shape))
                head += row
            elif isinstance(ext, pd.DataFrame):
                row = "|{:20s} |{:20s} |{:20s}\n".format(name, "table", str(len(ext)))
                head += row
        print(head)

    def _create_hdul(self):
        """
        Create an hdul in FITS format.
        This is used by the base model for writing data context to file
        """
        hdu_list = []
        hdu_definitions = self.extensions.items()
        for key, value in hdu_definitions:
            hduname = key
            if value == "PrimaryHDU":
                head = fits.Header(self.headers[key])
                hdu = fits.PrimaryHDU(header=head)
                hdu_list.insert(0, hdu)
            elif value == "ImageHDU":
                data = self.data[key]
                if data is None:
                    ndim = 0
                else:
                    ndim = len(data.shape)
                self.headers[key]["NAXIS"] = ndim
                if ndim == 0:
                    self.headers[key]["NAXIS1"] = 0
                else:
                    for d in range(ndim):
                        self.headers[key]["NAXIS{}".format(d + 1)] = data.shape[d]
                head = fits.Header(self.headers[key])
                try:
                    hdu = fits.ImageHDU(data=data, header=head)
                    hdu.name = hduname
                    hdu_list.append(hdu)
                except KeyError as ke:
                    print("KeyError exception raised: -->ke=" + str(ke))
                    print("Attempting to handle it...")
                    if str(ke) == "'bool'":
                        data = data.astype(float)
                        print("------>SHAPE=" + str(data.shape))
                        hdu = fits.ImageHDU(data=data, header=head)
                        hdu_list.append(hdu)
                    else:
                        raise KeyError("A different error...")
            elif value == "BinTableHDU":
                table = Table.from_pandas(self.data[key])
                self.headers[key]["NAXIS1"] = len(table)
                head = fits.Header(self.headers[key])
                hdu = fits.BinTableHDU(data=table, header=head)
                hdu.name = hduname
                hdu_list.append(hdu)
            else:
                print(
                    "Can't translate {} into a valid FITS format.".format(
                        type(self.data[key])
                    )
                )
                continue

        return hdu_list
