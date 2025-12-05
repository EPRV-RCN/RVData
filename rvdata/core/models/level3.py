"""
Level 3 Data Model for RV spectral data
"""

# External dependencies
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd

import rvdata.core.models.base
from rvdata.core.models.definitions import LEVEL3_EXTENSIONS, LEVEL3_PRIMARY_KEYWORDS
from rvdata.core.tools import stitch_spectrum
from rvdata.core.tools.utils import create_configdict_from_file


class RV3(rvdata.core.models.base.RVDataModel):
    """
    The level 3 RV data. Initialized with empty fields.
    Attributes inherited from RVDataModel, additional attributes below.

    """

    def __init__(self):
        super().__init__()
        self.level = 3

        for i, row in LEVEL3_EXTENSIONS.iterrows():
            if row["Required"]:
                # TODO: set description and comment
                if row["Name"] not in self.extensions.keys():
                    self.create_extension(row["Name"], row["DataType"])

        # Add EXT_DESCRIPT as a DataFrame, dropping the Comments column
        ext_descript = (
            LEVEL3_EXTENSIONS.copy().query("Required == True").reset_index(drop=True)
        )
        if "Comments" in ext_descript.columns:
            ext_descript = ext_descript.drop(columns=["Comments"])
        self.set_data("EXT_DESCRIPT", ext_descript)

        # TODO: initialize the LEVEL3_PRIMARY_KEYWORDS in the Level 3 PRIMARY header
        # Add L3 specific entries to the primary header

        for i, row in LEVEL3_PRIMARY_KEYWORDS.iterrows():
            key = row["Keyword"]
            value = row["Default"]
            try:
                if row["Data type"].lower() == "uint":
                    self.headers["PRIMARY"][key] = int(value)
                elif row["Data type"].lower() == "float":
                    self.headers["PRIMARY"][key] = float(value)
                elif row["Data type"].lower() == "string":
                    self.headers["PRIMARY"][key] = str(value)
                elif row["Data type"].lower() == "double":
                    self.headers["PRIMARY"][key] = np.float64(value)
                elif row["Data type"].lower() == "boolean":
                    self.headers["PRIMARY"][key] = value.upper() == "TRUE"
                else:
                    print(f"Unknown type {row['Data type']} for keyword {key}")
            except (TypeError, AttributeError, ValueError):
                print(
                    f"Cannot convert value {value} for keyword {key} to type {row['Data type']}"
                )

    def _read(self, hdul: fits.HDUList) -> None:
        l3_ext = LEVEL3_EXTENSIONS.set_index("Name")
        for hdu in hdul:
            fits_type = l3_ext.loc[hdu.name]["DataType"]
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

    @staticmethod
    def convert_level2_to_level3(l2obj) -> None:
        """
        Read data from a Level 2 RVDataModel object and populate Level 3 fields
        """

        l3obj = RV3()

        # Set up the primary header
        l3prihdr = l2obj.headers["PRIMARY"]
        l3prihdr["DATALVL"] = 3

        # TODO iterate over traces
        # read the wavelength, flux, and blaze data

        sci_flx = l2obj.data["TRACE1_FLUX"].astype(np.float64)
        sci_wav = l2obj.data["TRACE1_WAVE"].astype(np.float64)
        sci_blz = l2obj.data["TRACE1_BLAZE"].astype(np.float64)

        # get instrument stitching config
        inst = l3prihdr["INSTRUME"].lower()
        stitch_config = create_configdict_from_file(
            f"rvdata/instruments/{inst}/config/{inst}_level3.config"
        )

        # stitch the orders
        try:
            st_wave, st_flux = stitch_spectrum.stitch_orders(
                sci_wav, sci_flx, sci_blz, inst_stitch_config=stitch_config
            )
            # save the stitched spectrum
            l3prihdr["BLZCORR"] = True
            l3prihdr["LMPCORR"] = True
            l3prihdr["SEDCORR"] = False
            l3prihdr["INTERPMD"] = "BINDENSITY"
            l3prihdr["FLXNRMMD"] = "None"
            l3prihdr["DISPCORR"] = True

        except Exception as e:
            print(f"Error stitching orders: {e}")
            l3prihdr["BLZCORR"] = False
            l3prihdr["LMPCORR"] = False
            l3prihdr["SEDCORR"] = False
            l3prihdr["INTERPMD"] = "None"
            l3prihdr["FLXNRMMD"] = "None"
            l3prihdr["DISPCORR"] = False

        l3obj.set_header("PRIMARY", l3prihdr)
        l3obj.set_header("INSTRUMENT_HEADER", l2obj.headers["INSTRUMENT_HEADER"])
        l3obj.set_data("STITCHED_CORR_SCI_WAVE", st_wave)
        l3obj.set_data("STITCHED_CORR_SCI_FLUX", st_flux)

        return l3obj

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
