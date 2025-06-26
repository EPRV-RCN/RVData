from astropy.io import fits

# from astropy.table import Table
import numpy as np
import pandas as pd
import os
from collections import OrderedDict

# # import base class
from rvdata.core.models.level4 import RV4

# from rvdata.core.models.definitions import LEVEL4_EXTENSIONS


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

        # Instrument header
        self.set_header("INSTRUMENT_HEADER", hdul["PRIMARY"].header)

        # Set up the primary header - code from L2 translator to handle OBSMODE dependent entries

        # Set up for obs-mode dependent primary header entries
        mode_dep_phead = {}
        catalogue_map = {
            "CID": "QOBJECT",
            "CRA": "QRA",
            "CDEC": "QDEC",
            "CEQNX": "QEQNX",
            "CEPCH": "QEPOCH",
            "CPLX": "QPLX",
            "CPMR": "QPMRA",
            "CPMD": "QPMDEC",
            "CRV": "QRV",
            "CZ": "QZ",
        }

        # Check observation mode to set number of traces
        if hdul[0].header["OBS-MODE"] == "HR":
            fiber_list = ["SCI", "SKY", "CAL"]
            mode_dep_phead["CLSRC3"] = hdul[0].header["CAL-OBJ"]
        elif hdul[0].header["OBS-MODE"] == "HE":
            fiber_list = ["SCI", "SKY"]
        mode_dep_phead["NUMTRACE"] = len(fiber_list)

        for i_fiber, fiber in enumerate(fiber_list):
            mode_dep_phead[f"TRACE{i_fiber+1}"] = hdul[0].header[f"{fiber}-OBJ"]

            if hdul[0].header["OBSTYPE"] == "Cal":
                mode_dep_phead[f"CLSRC{i_fiber+1}"] = hdul[0].header[f"{fiber}-OBJ"]

            if hdul[0].header[f"{fiber}-OBJ"] == hdul[0].header["QOBJECT"]:
                for pkey, ikey in catalogue_map.items():
                    mode_dep_phead[f"{pkey}{i_fiber+1}"] = hdul[0].header[ikey]
                mode_dep_phead[f"CSRC{i_fiber+1}"] = "GAIADR2"

        # Set up data standard primary header
        hmap_path = os.path.join(os.path.dirname(__file__), "config/header_map.csv")
        headmap = pd.read_csv(hmap_path, header=0)

        phead = fits.PrimaryHDU().header
        ihead = self.headers["INSTRUMENT_HEADER"]
        for i, row in headmap.iterrows():
            skey = row["STANDARD"]
            instkey = row["INSTRUMENT"]
            if row["MODE_DEP"] != "Y":
                if pd.notnull(instkey):
                    instval = ihead[instkey]
                else:
                    instval = row["DEFAULT"]
                if pd.notnull(instval):
                    phead[skey] = instval
                else:
                    phead[skey] = None
            else:
                if skey in mode_dep_phead.keys():
                    phead[skey] = mode_dep_phead[skey]
                else:
                    continue

        # Add instrument era
        eramap = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "config/neid_inst_eras.csv")
        )
        era_time_diffs = phead["JD_UTC"] - eramap["startdate"].values
        era = eramap["era"].values[
            np.argmin(era_time_diffs[np.where(era_time_diffs >= 0)[0]])
        ]
        phead["INSTERA"] = era

        self.set_header("PRIMARY", phead)

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
                        hdul["CCFS"].header[f"CCFRV{173-order:03d}"]
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
                        hdul["CCFS"].header[f"CCFWT{173-order:03d}"]
                        for order in range(122)
                    ]
                ),
            }
        )

        for order in range(122):
            if (
                np.isfinite(neid_fsr["fsr_start"].values[order])
                and (rv_table_data["RV"][order] != None)
                and (rv_table_data["RV"][order] != 0)
            ):
                fsr_pixel_start = int(neid_fsr["fsr_start"].values[order])
                fsr_pixel_end = int(neid_fsr["fsr_end"].values[order])

                rv_table_data["wave_start"][order] = hdul["sciwave"].data[
                    order, fsr_pixel_start
                ]
                rv_table_data["wave_end"][order] = hdul["sciwave"].data[
                    order, fsr_pixel_end
                ]

                rv_table_data["pixel_start"][order] = fsr_pixel_start
                rv_table_data["pixel_end"][order] = fsr_pixel_end

        self.set_data("RV1", pd.DataFrame(rv_table_data))

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

        # Diagnostics extension - for now just put in the activity extension directly

        self.create_extension(
            "DIAGNOSTICS1",
            "BinTableHDU",
            data=hdul["ACTIVITY"].data,
            header=hdul["ACTIVITY"].header,
        )
