from astropy.io import fits

# from astropy.table import Table
import numpy as np
import pandas as pd
import os
from collections import OrderedDict

from rvdata.core.models.level4 import RV4

# NEID specific utility functions
from rvdata.instruments.neid.utils import make_neid_primary_header


# KPF Level2 Reader
class NEIDRV4(RV4):
    """
    Data model and reader for RVData Level 4 (RV) data constructed from
    NEID Level 2 pipeline products.

    This class extends the `RV4` base class to handle NEID data. It
    reads the relevant extensions from a NEID Level 2 FITS file and
    organizes them into a standardized format.

    Parameters
    ----------
    Inherits all parameters from :class:`RV4`.

    Attributes
    ----------
    extensions : dict
        Dictionary of all created extensions (e.g., 'RV1', 'CCF1', etc.),
        mapping extension names to their data arrays.
    headers : dict
        Dictionary of headers for each extension, mapping extension names to
        their FITS headers.
    data : dict
        Dictionary of data arrays for each extension.

    Notes
    -----
    To construct an RVData Level 4 object, a NEID Level 2 FITS file is required.
    The classmethod `from_fits` should be used to instantiate the object
    from these files. The `_read` method is not intended to be called
    directly by users.

    Example
    -------
    >>> from rvdata.instruments.neid.level4 import NEIDRV4
    >>> neid_rv4_obj = NEIDRV4.from_fits("neidL2_YYYYMMDDTHHMMSS.fits", instrument="NEID")
    >>> neid_rv4_obj.to_fits("neid_L4_standard.fits")
    """

    def _read(self, hdul: fits.HDUList, **kwargs) -> None:

        # Set up extension description table
        ext_table = {
            "extension_name": [],
            "description": [],
        }

        # Set up the primary header
        phead = make_neid_primary_header.make_base_primary_header(hdul[0].header)
        phead["DATALVL"] = "L4"

        # Add RV specific entries to the primary header
        phead["BJDTDB"] = hdul["CCFS"].header["CCFJDMOD"]
        phead["RV"] = hdul["CCFS"].header["CCFRVMOD"]
        phead["RVERR"] = hdul["CCFS"].header["DVRMSMOD"]
        phead["RVMETHOD"] = "CCF"
        phead["SYSVEL"] = hdul["PRIMARY"].header["QRV"]

        ext_table["extension_name"].append("PRIMARY")
        ext_table["description"].append("EPRV Standard FITS HEADER (no data)")

        # Instrument header
        self.set_header("INSTRUMENT_HEADER", hdul["PRIMARY"].header)
        ext_table["extension_name"].append("INSTRUMENT_HEADER")
        ext_table["description"].append("Primary header of native instrument file")

        # RV1 - turn the CCFS extension header into a table

        neid_fsr = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "config/neid_fsr.csv")
        )

        rv_table_data = OrderedDict(
            {
                "BJD_TDB": np.array(
                    [hdul[0].header[f"SSBJD{173-order:03d}"] for order in range(122)]
                ),
                "RV": np.array(
                    [
                        (
                            hdul["CCFS"].header[f"CCFRV{173-order:03d}"]
                            if not (hdul["CCFS"].data[order] == 0).all()
                            else np.nan
                        )
                        for order in range(122)
                    ]
                ),
                "RV_error": np.full(122, np.nan),
                "BC_vel": np.array(
                    [hdul[0].header[f"SSBRV{173-order:03d}"] for order in range(122)]
                ),
                "wave_start": np.full(122, np.nan),
                "wave_end": np.full(122, np.nan),
                "pixel_start": np.full(122, np.nan),
                "pixel_end": np.full(122, np.nan),
                "order_index": np.arange(122),
                "echelle_order": 173 - np.arange(122),
                "weight": np.array(
                    [
                        (
                            hdul["CCFS"].header[f"CCFWT{173-order:03d}"]
                            if hdul["CCFS"].header[f"CCFWT{173-order:03d}"] is not None
                            else np.nan
                        )
                        for order in range(122)
                    ]
                ),
            }
        )

        # Add BERV to the primary header and write primary header
        order_bjd_sort = rv_table_data["BJD_TDB"][np.argsort(rv_table_data["BJD_TDB"])]
        order_berv_sort = rv_table_data["BC_vel"][np.argsort(rv_table_data["BJD_TDB"])]
        phead["BERV"] = np.interp(phead["BJDTDB"], order_bjd_sort, order_berv_sort)

        self.set_header("PRIMARY", phead)

        # Add information about the wavelength/pixel extents of the RV computation per order
        for order in range(122):
            if (
                np.isfinite(neid_fsr["fsr_start"].values[order])
                and np.isfinite(rv_table_data["RV"][order])
            ):
                fsr_pixel_start = int(neid_fsr["fsr_start"].values[order])
                fsr_pixel_end = int(neid_fsr["fsr_end"].values[order])

                rv_table_data["wave_start"][order] = hdul["SCIWAVE"].data[
                    order, fsr_pixel_start
                ]
                rv_table_data["wave_end"][order] = hdul["SCIWAVE"].data[
                    order, fsr_pixel_end
                ]

                rv_table_data["pixel_start"][order] = fsr_pixel_start
                rv_table_data["pixel_end"][order] = fsr_pixel_end

        self.set_data("RV1", pd.DataFrame(rv_table_data))
        ext_table["extension_name"].append("RV1")
        ext_table["description"].append(
            "Order-wise RV measurement table for NEID Science fiber trace"
        )

        # CCF extension

        # Translate some NEID CCFS extension header keys to standard keys
        ccf_header_map = {
            "VELSTART": "CCFSTART",
            "VELSTEP": "CCFSTEP",
            "CCFMASK": "CCFMASK",
        }
        ccf_header = OrderedDict()

        for ccf_skey, ccf_ikey in ccf_header_map.items():
            ccf_header[ccf_skey] = hdul["CCFS"].header[ccf_ikey]
        ccf_header["VELNSTEP"] = hdul["CCFS"].data.shape[1]

        self.create_extension(
            "CCF1", "ImageHDU", data=hdul["CCFS"].data, header=ccf_header
        )
        ext_table["extension_name"].append("CCF1")
        ext_table["description"].append("Order-wise CCFs for NEID Science fiber trace")

        # Diagnostics extension - for now just put in the activity extension directly

        diagnostics_table_data = {"metric_name": [], "value": [], "uncertainty": []}

        # Activity indicators in NEID extension to include in diagnostics table
        neid_activity_use = {
            "CaIIHK": "CaIIHK",
            "HeI_1": "HeI",
            "NaI": "NaI",
            "Ha06_1": "Halpha06",
            "Ha16_1": "Halpha16",
            "CaI_1": "CaI",
            "CaIRT1": "CaIRT1",
            "CaIRT2": "CaIRT2",
            "CaIRT3": "CaIRT3",
            "NaINIR": "NaINIR",
            "PaDelta": "PaDelta",
            "Mn539": "Mn539",
        }

        for index_name_neid, index_name_std in neid_activity_use.items():
            # Location in the NEID table with that index
            table_loc = np.where(hdul["ACTIVITY"].data["INDEX"] == index_name_neid)[0]

            diagnostics_table_data["metric_name"].append(index_name_std)
            diagnostics_table_data["value"].append(
                hdul["ACTIVITY"].data["VALUE"][table_loc][0]
            )
            diagnostics_table_data["uncertainty"].append(
                hdul["ACTIVITY"].data["UNCERTAINTY"][table_loc][0]
            )

            # Also output the telluric corrected version
            table_loc = np.where(
                hdul["ACTIVITY"].data["INDEX"] == (index_name_neid + "_tellcorr")
            )[0]

            diagnostics_table_data["metric_name"].append(index_name_std + "_tellcorr")
            diagnostics_table_data["value"].append(
                hdul["ACTIVITY"].data["VALUE"][table_loc][0]
            )
            diagnostics_table_data["uncertainty"].append(
                hdul["ACTIVITY"].data["UNCERTAINTY"][table_loc][0]
            )

        # Add CCF activity indicators as well
        diagnostics_table_data["metric_name"].append("CCF_FWHM")
        diagnostics_table_data["value"].append(hdul["CCFS"].header["FWHMMOD"])
        diagnostics_table_data["uncertainty"].append(np.nan)

        diagnostics_table_data["metric_name"].append("CCF_BIS")
        diagnostics_table_data["value"].append(hdul["CCFS"].header["BISMOD"])
        diagnostics_table_data["uncertainty"].append(hdul["CCFS"].header["EBISMOD"])

        # Output the diagnostics table
        self.create_extension(
            "DIAGNOSTICS1",
            "BinTableHDU",
            data=pd.DataFrame(diagnostics_table_data),
        )
        ext_table["extension_name"].append("DIAGNOSTICS1")
        ext_table["description"].append(
            "Table of activity diagnostics for NEID science fiber trace"
        )

        # Set extension description table
        self.set_data("EXT_DESCRIPT", pd.DataFrame(ext_table))
