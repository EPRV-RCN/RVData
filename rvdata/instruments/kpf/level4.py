from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd
import os
from collections import OrderedDict

# import base class
from rvdata.core.models.level4 import RV4
from rvdata.core.models.definitions import LEVEL4_EXTENSIONS


# KPF Level2 Reader
class KPFRV4(RV4):
    """
    Data model and reader for RVData Level 4 (RV) data constructed from
    KPF Level 2 pipeline products.

    This class extends the `RV4` base class to handle Keck Planet Finder (KPF)
    data. It reads the relevant science and calibration extensions
    from a KPF Level 2 file and organizes them into a
    standardized format.

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
    To construct an RVData Level 4 object, a KPF Level 2 FITS file is required.
    The classmethod `from_fits` should be used to
    instantiate the object from these files. The `_read` method is not intended
    to be called directly by users.

    Example
    -------
    >>> from rvdata.instruments.kpf.level4 import KPFRV4
    >>> obj = KPFRV4.from_fits("kpf_L2.fits"")
    >>> obj.to_fits("kpf_L4_standard.fits")
    """

    def _read(self, hdul2: fits.HDUList, **kwargs) -> None:

        # construct the extensions
        ccf_array = None
        for c, chip in enumerate(["GREEN", "RED"]):
            ccf_ext = f"{chip}_CCF_RW"

            if ccf_array is None:
                ccf_array = hdul2[ccf_ext].data
                ccf_meta = OrderedDict(hdul2[ccf_ext].header)
            else:
                ccf_array = np.concatenate((ccf_array, hdul2[ccf_ext].data), axis=1)

        self.create_extension("CCF1", "ImageHDU", data=ccf_array, header=ccf_meta)

        # set the primary header
        hmap_path = os.path.join(os.path.dirname(__file__), "config/header_map.csv")
        headmap = pd.read_csv(hmap_path, header=0)
        phead = fits.PrimaryHDU().header
        ihead = hdul2["PRIMARY"].header
        for i, row in headmap.iterrows():
            skey = row["STANDARD"]
            kpfkey = row["INSTRUMENT"]
            if pd.notnull(kpfkey):
                kpfval = ihead[kpfkey]
            else:
                kpfval = row["DEFAULT"]
            if pd.notnull(kpfval):
                phead[skey] = kpfval
            else:
                phead[skey] = None

        self.set_header("PRIMARY", phead)
        self.set_header("INSTRUMENT_HEADER", ihead)

        rv_ext = "RV"
        # column mapping between KPF and data standard, keys are standard names
        # and values are KPF names
        colmap = {
            "BJD_TDB": "CCFBJD",
            "RV_TRACE2": "orderlet1",
            "RV_TRACE3": "orderlet2",
            "RV_TRACE4": "orderlet3",
            "COMBINED_RV": "RV",
            "COMBINED_RV_ERR": "RV error",
            "BCVEL": "Bary_RVC",
            "wave_start": "s_wavelength",
            "wave_end": "e_wavelength",
            "order_index": "order no.",
            "weight": "CCF Weights",
        }

        sysvel = hdul2["PRIMARY"].header.get("TARGRADV", 0.0)
        arr = Table(hdul2[rv_ext].data).to_pandas()
        rvdata = pd.DataFrame(arr)
        for std_col, kpf_col in colmap.items():
            if kpf_col in rvdata.columns:
                rvdata[std_col] = rvdata[kpf_col]
                rvdata.drop(columns=[kpf_col], inplace=True)
            else:
                print(f"Warning: {kpf_col} not found in KPF data, skipping.")

        for col in rvdata.columns:
            if col not in colmap.keys():
                rvdata.drop(columns=[col], inplace=True)

        rvdata["echelle_order"] = 137 - rvdata["order_index"]
        rvdata["COMBINED_RV"] = rvdata["COMBINED_RV"] - sysvel
        for i in [2, 3, 4]:
            rvdata["RV_TRACE{}".format(i)] = rvdata["RV_TRACE{}".format(i)] - sysvel

        rv_header = OrderedDict(hdul2[rv_ext].header)
        rv_header["RVSTART"] = ccf_meta["STARTV"]
        rv_header["RVSTEP"] = ccf_meta["STEPV"]
        rv_header["MASK"] = ccf_meta["SCI_MASK"]

        self.set_data("RV1", rvdata)
        self.set_header("RV1", rv_header)

        self.set_header("DRP_CONFIG", OrderedDict(hdul2["CONFIG"].header))
        self.set_data("DRP_CONFIG", Table(hdul2["CONFIG"].data).to_pandas())

        self.set_header("RECEIPT", OrderedDict(hdul2["RECEIPT"].header))
        self.set_data("RECEIPT", Table(hdul2["RECEIPT"].data).to_pandas())

        all_exts = list(self.extensions.keys())
        for ext_name in all_exts:
            if ext_name not in LEVEL4_EXTENSIONS["Name"].values:
                self.del_extension(ext_name)
