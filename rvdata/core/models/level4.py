"""
Level 2 Data Model for RV spectral data
"""

# External dependencies
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd

import rvdata.core.models.base
from rvdata.core.models.definitions import (
    BASE_DRP_CONFIG_COLUMNS,
    BASE_RECEIPT_COLUMNS,
    LEVEL2_PRIMARY_KEYWORDS,
    LEVEL4_EXTENSIONS,
    LEVEL4_PRIMARY_KEYWORDS,
    LEVEL4_RV_TABLE_COLUMNS,
)
from rvdata.core.tools.headers import parse_value_to_datatype


class RV4(rvdata.core.models.base.RVDataModel):
    """
    The level 4 RV data. Initialized with empty fields.
    Attributes inherited from RVDataModel, additional attributes below.

    """

    def __init__(self):
        super().__init__()
        self.level = 4

        for _, row in LEVEL4_EXTENSIONS.iterrows():
            if row["Required"] and row["Name"] not in self.extensions.keys():
                self.create_extension(row["Name"], row["DataType"])

        # initialize PRIMARY header keywords to defaults with units and descriptions
        for _, row in pd.concat(
            [LEVEL2_PRIMARY_KEYWORDS, LEVEL4_PRIMARY_KEYWORDS]
        ).iterrows():
            if row["Required"]:
                keyword = row["Keyword"].split()[0]
                datatype = row["DataType"]
                default = row["Default"]
                description = row["Description"]
                units = row["Units"]
                if pd.isna(units) or units == "" or units.lower() == "N/A".lower():
                    unitstr = ""
                else:
                    unitstr = f"[{units}] "
                self.headers["PRIMARY"][keyword] = parse_value_to_datatype(
                    keyword, datatype, (default, f"{unitstr}{description}")
                )

        # Add EXT_DESCRIPT as a DataFrame
        # Only use the Name and Description columns
        ext_descript = (
            LEVEL4_EXTENSIONS.copy()
            .query("Required == True")[["Name", "Description"]]
            .reset_index(drop=True)
        )
        self.set_data("EXT_DESCRIPT", ext_descript)

        # Initialize INSTRUMENT_HEADER with a dummy zero image
        self.set_data("INSTRUMENT_HEADER", np.zeros((1,), dtype=np.float32))

        # Initialize RECEIPT with receipt columns
        receipt_columns = BASE_RECEIPT_COLUMNS["Name"].tolist()
        self.set_data("RECEIPT", pd.DataFrame(columns=receipt_columns))

        # Initialize DRP_CONFIG with columns from definition
        drp_config_columns = BASE_DRP_CONFIG_COLUMNS["Name"].tolist()
        self.set_data("DRP_CONFIG", pd.DataFrame(columns=drp_config_columns))

        # Initialize RV1 with columns from definition
        order_table_columns = LEVEL4_RV_TABLE_COLUMNS["Name"].tolist()
        self.set_data("RV1", pd.DataFrame(columns=order_table_columns))

    def _read(self, hdul: fits.HDUList) -> None:
        l4_ext = LEVEL4_EXTENSIONS.set_index("Name")
        for hdu in hdul:
            if "RV" in hdu.name:
                fits_type = "BinTableHDU"
            elif "CCF" in hdu.name:
                fits_type = "ImageHDU"
            else:
                fits_type = l4_ext.loc[hdu.name]["DataType"]
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
                head = fits.Header()
                for keyword, content in self.headers[key].items():
                    head[keyword] = content
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
