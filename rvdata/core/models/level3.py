"""
Level 3 Data Model for RV spectral data
"""

# External dependencies
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd

import rvdata.core.models.base
from rvdata.core.models.definitions import (
    BASE_DRP_CONFIG_COLUMNS,
    BASE_ORDER_TABLE_COLUMNS,
    BASE_RECEIPT_COLUMNS,
    LEVEL2_PRIMARY_KEYWORDS,
    LEVEL3_EXTENSIONS,
    LEVEL3_PRIMARY_KEYWORDS,
)
from rvdata.core.tools.headers import parse_value_to_datatype
from rvdata.core.tools.utils import create_configdict_from_file
import rvdata.core.tools.stitch_spectrum as stitch_spectrum


class RV3(rvdata.core.models.base.RVDataModel):
    """
    The level 3 RV data. Initialized with empty fields.
    Attributes inherited from RVDataModel, additional attributes below.

    """

    def __init__(self):
        super().__init__()
        self.level = 3

        for _, row in LEVEL3_EXTENSIONS.iterrows():
            if row["Required"] and row["Name"] not in self.extensions.keys():
                self.create_extension(row["Name"], row["DataType"])

        # initialize PRIMARY header keywords to defaults with units and descriptions
        for _, row in pd.concat(
            [LEVEL2_PRIMARY_KEYWORDS, LEVEL3_PRIMARY_KEYWORDS]
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
            LEVEL3_EXTENSIONS.copy()
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

        # Initialize ORDER_TABLE with columns from definition
        order_table_columns = BASE_ORDER_TABLE_COLUMNS["Name"].tolist()
        self.set_data("ORDER_TABLE", pd.DataFrame(columns=order_table_columns))

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

    def convert_level2_to_level3(self, l2obj) -> None:
        """
        Read data from a Level 2 RVDataModel object and populate Level 3 fields

        Parameters
        ----------
        l2obj : RV2, RVDataModel
            A Level 2 RVDataModel object containing spectral data and headers to convert.

        Returns
        -------
        None
            This method modifies the current RV3 object in place.
        """

        # Set up the primary header
        l3prihdr = l2obj.headers["PRIMARY"]
        l3prihdr["DATALVL"] = "L3"

        # get instrument stitching config
        inst = l3prihdr["INSTRUME"].lower()
        inst_stitch_config = create_configdict_from_file(
            f"rvdata/instruments/{inst}/config/{inst}_level3.config"
        )

        traces = np.arange(1, l3prihdr["NUMTRACE"] + 1)
        traces2stitch = []
        for trace_num in traces:
            if l3prihdr[f"CLSRC{trace_num}"] is None:
                traces2stitch.append(trace_num)
            else:
                continue

        st_wav: dict[int, np.ndarray] = {}
        st_flx: dict[int, np.ndarray] = {}
        st_var: dict[int, np.ndarray] = {}
        # stitch the orders for each trace
        try:
            for trace_num in traces2stitch:
                # read the wavelength, flux, and blaze data
                sci_flx = l2obj.data[f"TRACE{trace_num}_FLUX"].astype(np.float64)
                sci_wav = l2obj.data[f"TRACE{trace_num}_WAVE"].astype(np.float64)
                sci_blz = l2obj.data[f"TRACE{trace_num}_BLAZE"].astype(np.float64)
                sci_var = l2obj.data[f"TRACE{trace_num}_VAR"].astype(np.float64)
                order_table = l2obj.data["ORDER_TABLE"]

                # Masking invalid data
                sci_flxm, sci_wavm, sci_blzm, sci_varm, spec_mask = (
                    stitch_spectrum.mask_bad_spectrum_data(
                        sci_flx, sci_wav, sci_blz, sci_var
                    )
                )

                # Calculate the normalized blaze function
                sci_blzme = stitch_spectrum.calculate_normalized_blaze_function(
                    sci_wavm, sci_blzm, inst_stitch_config, order_table
                )

                # Deblazed science flux and deblazed science variance
                sci_dflxm = np.divide(sci_flxm, sci_blzme)
                sci_dvarm = np.divide(sci_varm, sci_blzme**2)
                sci_dcovm = sci_dvarm[None, :]

                # Define a common wavelength grid with constant velocity spacing
                wavegrid = stitch_spectrum.get_wavelength_grid_with_constant_velocity(
                    inst_stitch_config["wavegrid_start"],
                    inst_stitch_config["wavegrid_end"],
                    inst_stitch_config["velpix"],
                )

                # Stitch the deblazed spectrum into wavelength grid using inverse-variance weighting
                st_wav[trace_num], st_flx[trace_num], st_var[trace_num] = (
                    stitch_spectrum.stitch_deblazed_spectrum(
                        wavegrid, sci_wavm, sci_dflxm, sci_dcovm
                    )
                )

        except Exception as e:
            print(f"Error stitching orders: {e}")
            l3prihdr["BLZCORR"] = False
            l3prihdr["LMPCORR"] = False
            l3prihdr["SEDCORR"] = False
            l3prihdr["INTERPMD"] = "None"
            l3prihdr["FLXNRMMD"] = "None"
            l3prihdr["DISPCORR"] = False
        else:
            # all good, set the header keywords
            l3prihdr["BLZCORR"] = True
            l3prihdr["LMPCORR"] = True
            l3prihdr["SEDCORR"] = False
            l3prihdr["INTERPMD"] = "BINDENSITY"
            l3prihdr["FLXNRMMD"] = "None"
            l3prihdr["DISPCORR"] = True

        # TODO: if there are multiple traces, co-add all traces to produce the "SCI" extensions
        #       For now, just copy the first trace
        self.set_data("STITCHED_CORR_SCI_WAVE", st_wav[1])
        self.set_data("STITCHED_CORR_SCI_FLUX", st_flx[1])
        self.set_data("STITCHED_CORR_SCI_VAR", st_var[1])
        # TODO set data for STICHED_CORR_TRACE{n}_WAVE/FLUX/VAR if multiple traces

        self.set_header("PRIMARY", l3prihdr)
        self.set_header("INSTRUMENT_HEADER", l2obj.headers["INSTRUMENT_HEADER"])
