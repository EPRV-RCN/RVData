#!/usr/bin/env python
# coding: utf-8

from astropy.io import fits
import numpy as np
import pandas as pd
from collections import OrderedDict
import pickle
from rvdata.core.models.level4 import RV4
# from rvdata.core.models.level2 import RV2


class MAROONXRV4(RV4):
    """
    Data model and reader for RVData Level 4 (RV) data constructed from
    MAROONX Level 2 standard product.

    This class extends the `RV4` base class to handle MAROONX data. It
    reads the relevant extensions from a MAROONX Level 2 FITS file and
    the MAROONX RV data product pickle file, and
    organizes them into a standardized format.

    Parameters
    ----------
    Inherits all parameters from :class:`RV4`.

    Methods
    -------
    createL4(self, hdul: fits.HDUList, RVfile: str, channel: str) -> None:
        Reads the standard L2 fits file and SERVAL RV pickle file, extracts
        orderwise RV, activity indicators, and stores them
        in a standardized format.

    Attributes
    ----------
    extensions : dict
        Dictionary of all created extensions (e.g., 'RV1', 'DIAGNOSTICS1',etc.)
        mapping extension names to their data arrays.
    headers : dict
        Dictionary of headers for each extension, mapping extension names to
        their FITS headers.
    data : dict
        Dictionary of data arrays for each extension.

    Notes
    -----
    To construct an RVData Level 4 object, a MAROONX standard Level 2 FITS file
    and pickle file containing MAROONX RV data product generated
    with SERVAL, is required.
    The classmethod 'from_fits' should be used to instantiate the object
    from standard Level 2 FITS. 'createL4' should be
    used to construct Level 4 object.
    The '_read' method is not intended to be called directly by users.

    Example
    -------
    >>> from core.models.level2 import RV2
    >>> from core.models.level4 import RV4
    >>> L2_obj = RV2.from_fits("MAROONXRED_SL2_YYYYMMDDTHHMMSS.fits")
    >>> MX_rv4_obj.createL4(fits.open(L2_obj.filename), RVfile, channel)
    >>> MX_rv4_obj.to_fits(f'MAROONXRED_SL4_YYYYMMDDTHHMMSS.fits')
    """
    def _read(self, hdul: fits.HDUList) -> None:
        raise RuntimeError(
            "Incorrect usage: MAROONXRV4 does not support RV4.from_fits()\n"
            "Use createL4(hdul, RVfile, channel) instead.\n"
            "Example:\n"
            "    rv4 = MAROONXRV4()\n"
            "    L2 = RV2.from_fits('MAROONXRED_SL2.fits')\n"
            "    rv4.createL4(fits.open(L2.filename), 'serval.pkl', 'RED')\n"
            "    rv4.to_fits('MAROONXRED_SL4.fits')"
        )

    def createL4(self, hdul: fits.HDUList, RVfile: str, channel: str) -> None:
        """
        Reads a standard L2 fits for MAROONX and creates
        a standard L4 fits data product
        """
        ext_table = {
            "Name": [],
            "Description": [],
        }

        if channel.upper() == "RED":
            serval_min = 69
            serval_max = 94
        elif channel.upper() == "BLUE":
            serval_min = 93
            serval_max = 122
        else:
            raise ValueError("Invalid channel specified. Use 'RED' or 'BLUE'.")

        with open(RVfile, "rb") as f:
            data = pickle.load(f)
        # Extract relevant data from the pickle file
        bjds_pk = np.array(data["bjd"])
        rvs_comb = np.array(data["RV"])
        erv_comb = np.array(data["e_RV"])
        rvs_ord = np.array(data["rv"])  # orderwise
        erv_ord = np.array(data["e_rv"])  # orderwise

        bjd_fits = float(hdul["BJD_TDB"].data[0])
        order_tab = hdul["ORDER_TABLE"].data         
        phdr = hdul["PRIMARY"].header
        ihdr = hdul[1].header
        berv = float(ihdr['HIERARCH BERV_FLUXWEIGHTED_FRD'])/1000.0
        orders = order_tab["echelle_order"].astype(int)

        # Find the index in the pickle data that matches the BJD from FITS
        tol = 1e-4
        diff = np.abs(bjds_pk - bjd_fits)
        idx = np.argmin(diff)
        if diff[idx] > tol:
            raise ValueError("No matching BJD found within tolerance!")
        print(f"Matched BJD to SERVAL index: {idx}")

        rv_serval = rvs_ord[idx][::-1]
        erv_serval = erv_ord[idx][::-1]
        rv_sel = float(rvs_comb[idx])
        erv_sel = float(erv_comb[idx])
        bjd_tdb = float(bjds_pk[idx])
        rv_full = np.full(len(orders), np.nan)
        erv_full = np.full(len(orders), np.nan)
        valid = (orders >= serval_min) & (orders <= serval_max)
        rv_full[valid] = rv_serval
        erv_full[valid] = erv_serval

        # Set up the primary header
        phead = phdr
        phead["DATALVL"] = 'L4'
        # Add RV specific entries to the primary header
        phead["BJDTDB"] = bjd_tdb
        phead["RV"] = rv_sel/1000
        phead["RVERR"] = erv_sel/1000
        phead["RVMETHOD"] = "SERVAL"
        phead["BERV"] = berv
        phead["SYSVEL"] = float(phdr['CRV2'])
        self.set_header("PRIMARY", phead)
        ext_table["Name"].append("PRIMARY")
        ext_table["Description"].append("EPRV Standard FITS HEADER")

        # Instrument header
        self.set_header("INSTRUMENT_HEADER", ihdr)
        ext_table["Name"].append("INSTRUMENT_HEADER")
        ext_table["Description"].append(
            "Primary header of native instrument file")
        ext_table["Name"].append("RECEIPT")
        ext_table["Description"].append(
            "The list of operations performed on the file")
        ext_table["Name"].append("DRP_CONFIG")
        ext_table["Description"].append(
            "Pipeline details to go from the raw file to this file")
        # Order-wise RV data
        rv_table_data = OrderedDict([
            ("BJD_TDB", np.full(len(orders), bjd_tdb)),
            ("echelle_order", orders),
            ("order_index", order_tab["order_index"]),
            ("wave_start", order_tab["wave_start"]),
            ("wave_end", order_tab["wave_end"]),
            ("rv", rv_full),
            ("rv_err", erv_full),
            ("BERV", berv)
        ])
        self.set_data("RV1", pd.DataFrame(rv_table_data))
        ext_table["Name"].append("RV1")
        ext_table["Description"].append(
            "Order-wise RV table for MAROONX from SERVAL"
        )
        # Activity data
        activity_keys = [
            'dLW', 'e_dLW',
            'crx', 'e_crx',
            'irt_ind1_v', 'irt_ind1_e',
            'irt_ind2_v', 'irt_ind2_e',
            'irt_ind3_v', 'irt_ind3_e',
            'halpha_ind1_v', 'halpha_ind1_e',
            'nad_ind1_v', 'nad_ind1_e',
            'nad_ind2_v', 'nad_ind2_e',
        ]
        BJD_column = []
        metric_names = []
        values = []
        errors = []
        for key in activity_keys:
            if key.startswith("e_") or key.endswith("_e"):
                continue
            if "ind" in key:
                err_key = key.replace("_v", "_e")
                metric_name = key.replace("_v", "")
            else:
                err_key = f"e_{key}"
                metric_name = key
            value = data[key][idx]
            error = data[err_key][idx]
            BJD_column.append(bjd_tdb)
            metric_names.append(metric_name)
            values.append(value)
            errors.append(error)
        activity_dict = OrderedDict([
            ("BJD_TDB",     BJD_column),
            ("metric_name", metric_names),
            ("Value",       values),
            ("Error",       errors),
        ])

        # Output the diagnostics table
        self.create_extension(
            "DIAGNOSTICS1",
            "BinTableHDU",
            data=pd.DataFrame(activity_dict),
        )
        ext_table["Name"].append("DIAGNOSTICS1")
        ext_table["Description"].append(
            "Table of activity results for MAROONX from SERVAL"
        )
        # Set extension Description table
        self.set_data("EXT_DESCRIPT", pd.DataFrame(ext_table))
