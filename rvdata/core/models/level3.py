"""
Level 3 Data Model for RV spectral data
"""

import re

import numpy as np
import pandas as pd
import importlib.resources
from astropy.io import fits
from astropy.table import Table, vstack

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

    def _get_min_bit_depth(self, ext_name):
        """Look up MinBitDepth for an ImageHDU extension from the L3 config."""
        # Handle multiplicity:
        #   STITCHED_CORR_TRACE2_WAVE -> STITCHED_CORR_TRACE1_WAVE
        #   STITCHED_CUSTOMCORR2_TRACE2_WAVE -> STITCHED_CUSTOMCORR1_TRACE1_WAVE
        canonical = re.sub(r'(?<=CUSTOMCORR)\d+', '1', ext_name)
        canonical = re.sub(r'(?<=TRACE)\d+', '1', canonical)
        row = LEVEL3_EXTENSIONS[LEVEL3_EXTENSIONS["Name"] == canonical]
        if row.empty:
            row = LEVEL3_EXTENSIONS[LEVEL3_EXTENSIONS["Name"] == ext_name]
        if row.empty:
            return None
        val = row.iloc[0]["MinBitDepth"]
        if pd.isna(val):
            return None
        return int(val)

    def _read(self, hdul: fits.HDUList) -> None:
        import warnings

        l3_ext = LEVEL3_EXTENSIONS.set_index("Name")
        for hdu in hdul:
            # Get DataType from predefined extensions, or infer from HDU type
            if hdu.name in l3_ext.index:
                fits_type = l3_ext.loc[hdu.name]["DataType"]
            elif isinstance(hdu, fits.PrimaryHDU):
                fits_type = "PrimaryHDU"
            elif isinstance(hdu, fits.ImageHDU):
                fits_type = "ImageHDU"
                # Warn about non-standard extensions (except dynamic TRACE extensions)
                if not hdu.name.startswith("STITCHED_CORR_TRACE"):
                    warnings.warn(
                        f"Non-standard extension '{hdu.name}' found in L3 file. "
                        "This may indicate a malformed file.",
                        UserWarning,
                    )
            elif isinstance(hdu, fits.BinTableHDU):
                fits_type = "BinTableHDU"
                if hdu.name not in ["EXT_DESCRIPT", "RECEIPT", "DRP_CONFIG", "ORDER_TABLE"]:
                    warnings.warn(
                        f"Non-standard table extension '{hdu.name}' found in L3 file. "
                        "This may indicate a malformed file.",
                        UserWarning,
                    )
            else:
                continue  # Skip unknown HDU types

            if hdu.name not in self.extensions.keys():
                self.create_extension(hdu.name, fits_type)

            if hdu.name == "PRIMARY":
                pass
            elif fits_type == "ImageHDU":
                data = np.array(hdu.data)
                self.set_data(hdu.name, data)
            elif fits_type == "BinTableHDU":
                data = Table.read(hdu)
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
            elif isinstance(ext, Table):
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

        # get the order table from level 2 data
        order_table = l2obj.data["ORDER_TABLE"]

        # get instrument stitching config
        inst = l3prihdr["INSTRUME"].lower()
        inst_stitch_config = create_configdict_from_file(
            importlib.resources.files("rvdata").joinpath(
                f"instruments/{inst}/config/{inst}_level3.config"
            )
        )

        # determine which traces need to be stitched
        traces = np.arange(1, l3prihdr["NUMTRACE"] + 1)
        traces2stitch = []
        for trace_num in traces:
            this_trace = str(l3prihdr[f"TRACE{trace_num}"]).strip()
            clsrc = l3prihdr[f"CLSRC{trace_num}"]
            # Check for None (Python None or string "None" from FITS header)
            clsrc_is_none = clsrc is None or (isinstance(clsrc, str) and clsrc.strip().lower() == "none")
            if clsrc_is_none and (this_trace.lower() != "sky"):
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

        # Set the STITCHED_CORR_SCI_WAVE/FLUX/VAR extensions
        # if only one trace to stitch, set SCI extension as that trace
        if len(traces2stitch) == 1:
            sci_trace_num = traces2stitch[0]
            try:
                self.set_data("STITCHED_CORR_SCI_WAVE", st_wav[sci_trace_num])
                self.set_data("STITCHED_CORR_SCI_FLUX", st_flx[sci_trace_num])
                self.set_data("STITCHED_CORR_SCI_VAR", st_var[sci_trace_num])
            except KeyError:
                pass
        elif len(traces2stitch) > 1:
            # set STITCHED_CORR_TRACE{n}_WAVE/FLUX/VAR for each trace
            # Note: The L3 standard defines STITCHED_CORR_TRACE1_* as optional extensions.
            # For multi-trace instruments (e.g., KPF with 3 science fibers), we dynamically
            # create TRACE{n} extensions for each science trace. This is an intentional
            # extension of the standard to support instruments with multiple science fibers.
            for trace_num in traces2stitch:
                try:
                    # Create extensions if they don't exist (they're optional)
                    for suffix in ["WAVE", "FLUX", "VAR"]:
                        ext_name = f"STITCHED_CORR_TRACE{trace_num}_{suffix}"
                        if ext_name not in self.extensions:
                            self.create_extension(ext_name, "ImageHDU")
                    self.set_data(
                        f"STITCHED_CORR_TRACE{trace_num}_WAVE", st_wav[trace_num]
                    )
                    self.set_data(
                        f"STITCHED_CORR_TRACE{trace_num}_FLUX", st_flx[trace_num]
                    )
                    self.set_data(
                        f"STITCHED_CORR_TRACE{trace_num}_VAR", st_var[trace_num]
                    )
                except KeyError:
                    pass
            # TODO: if there are multiple traces, co-add all traces to produce the "SCI" extensions

            # Update EXT_DESCRIPT table with new trace extensions
            ext_descript = self.data["EXT_DESCRIPT"]
            for trace_num in traces2stitch:
                for suffix in ["WAVE", "FLUX", "VAR"]:
                    ext_name = f"STITCHED_CORR_TRACE{trace_num}_{suffix}"
                    if ext_name not in ext_descript["Name"]:
                        new_row = Table(
                            {"Name": [ext_name], "Description": [f"Stitched trace {trace_num} {suffix.lower()}"]}
                        )
                        ext_descript = vstack([ext_descript, new_row])
            self.set_data("EXT_DESCRIPT", ext_descript)

        # set the order table from level 2 data into level 3 object
        self.set_data("ORDER_TABLE", order_table)
        # set the headers into the level 3 object
        self.set_header("PRIMARY", l3prihdr)
        self.set_header("INSTRUMENT_HEADER", l2obj.headers["INSTRUMENT_HEADER"])
        self.set_header("ORDER_TABLE", l2obj.headers["ORDER_TABLE"])

        # Inherit receipt from L2 object and add conversion entry
        self.receipt = l2obj.receipt.copy()
        self.receipt_add_entry("convert_level2_to_level3", "PASS")
