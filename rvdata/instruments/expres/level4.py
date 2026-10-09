"""
RVData Level 4 reader for EXPRES.

The native input is an EXPRES "ccf" file (one per exposure, same basename as
the fitspec file, conventionally in a sibling ``ccf/`` directory next to
``fitspec/``):

* HDU 0: CCF header (pipeline version, mask, combined velocity ``V`` and
  error ``E_V`` in cm/s, barycentric MJD, quality flags).
* HDU 1 (no EXTNAME): combined CCF, columns ``V_grid`` (cm/s), ``ccf``, ``e_ccf``.
* HDU 2 (no EXTNAME): per-order CCFs, columns ``orders`` (absolute echelle
  order), ``ccfs``, ``errs``, ``v`` and ``e_v`` (cm/s).

The standard PRIMARY and INSTRUMENT_HEADER are inherited from the Level 2
translation of the matching fitspec file, which also supplies the
wavelength extent of each order and the activity indicators.
"""

import os
from collections import OrderedDict

import numpy as np
import pandas as pd
from astropy import constants as const
from astropy.io import fits

from rvdata.core.models.level4 import RV4
from rvdata.instruments.expres.level2 import EXPRESRV2, _is_structural

CM_PER_KM = 1e5

# fitspec HDU 1 header card -> (metric name, uncertainty card or None)
ACTIVITY_CARDS = OrderedDict(
    [
        ("S-VALUE", ("S-VALUE", None)),
        ("HALPHA", ("HALPHA", None)),
        ("HWIDTH", ("HWIDTH", None)),
        ("CCFFWHM", ("CCFFWHM", "CCFFWHME")),
        ("BIS", ("BIS", None)),
        ("CCFVSPAN", ("CCFVSPAN", None)),
        ("BIGAUSS", ("BIGAUSS", None)),
        ("SKEWNORM", ("SKEWNORM", None)),
        ("QUALITY", ("QUALITY", None)),
        ("WATRCOL", ("WATRCOL", None)),
    ]
)
# ccf HDU 0 header card -> metric name
CCF_CARDS = OrderedDict([("SNR", "SNR"), ("CHI2", "CHI2"), ("EXPCOUNT", "EXPCOUNT")])


def find_fitspec(ccf_path, l1_file=None):
    """Locate the fitspec file that belongs to a ccf file.

    Tries ``l1_file`` if given, then ``<ccf dir>/../fitspec/<basename>``.
    Raises FileNotFoundError listing every path tried.
    """
    tried = []
    if l1_file is not None:
        tried.append(os.fspath(l1_file))
    if ccf_path is not None:
        ccf_dir, base = os.path.split(os.path.abspath(ccf_path))
        tried.append(os.path.normpath(os.path.join(ccf_dir, os.pardir, "fitspec", base)))
    for path in tried:
        if os.path.exists(path):
            return path
    raise FileNotFoundError(
        "EXPRES Level 4 needs the fitspec file that matches the ccf file. "
        "Tried: " + ", ".join(tried) + ". Pass l1_file=<path to fitspec file> "
        "or place it at <ccf dir>/../fitspec/<same basename>."
    )


class EXPRESRV4(RV4):
    """
    Data model and reader for RVData Level 4 (RV) data constructed from an
    EXPRES ccf file plus its matching fitspec file.

    Example
    -------
    >>> from rvdata.core.models.level4 import RV4
    >>> l4 = RV4.from_fits("ccf/Sun_20240126.5073.fits", instrument="EXPRES")
    >>> l4 = RV4.from_fits("Sun_20240126.5073.ccf.fits", instrument="EXPRES",
    ...                    l1_file="Sun_20240126.5073.fitspec.fits")
    >>> l4.to_fits()
    """

    def _read(self, hdul: fits.HDUList, l1_file=None, **kwargs) -> None:
        ccf_head = hdul[0].header
        combined = hdul[1].data
        per_order = hdul[2].data

        fitspec_path = find_fitspec(hdul.filename(), l1_file)
        l2 = EXPRESRV2()
        with fits.open(fitspec_path, memmap=False) as l1hdul:
            l2.read(l1hdul, instrument="EXPRES")
            l1_head1 = l1hdul[1].header.copy()
            l1_head2 = l1hdul[2].header.copy()

        ext_table = {"Name": [], "Description": []}

        def describe(name, description):
            ext_table["Name"].append(name)
            ext_table["Description"].append(description)

        c_kms = const.c.to("km/s").value
        bjd_tdb = float(l2.data["BJD_TDB"][0])
        berv_kms = float(l1_head2["wtd_mdpt_bc"]) * c_kms

        phead = l2.headers["PRIMARY"].copy()
        phead["DATALVL"] = "L4"
        phead["BJDTDB"] = bjd_tdb
        phead["RV"] = float(ccf_head["V"]) / CM_PER_KM
        phead["RVERR"] = float(ccf_head["E_V"]) / CM_PER_KM
        phead["RVMETHOD"] = "CCF"
        phead["BERV"] = berv_kms
        # The pipeline subtracts no systemic velocity; for the Sun the
        # systemic velocity is zero by definition, otherwise it is unknown.
        phead["SYSVEL"] = 0.0 if phead["ISSOLAR"] else None
        self.set_header("PRIMARY", phead)
        describe("PRIMARY", "EPRV Standard FITS HEADER (no data)")

        self.set_header("INSTRUMENT_HEADER", l2.headers["INSTRUMENT_HEADER"])
        describe("INSTRUMENT_HEADER", "Inherited EXPRES fitspec primary header (no data)")
        describe("RECEIPT", "Table of operations that have been performed on this file")

        entries = list(l2.data["DRP_CONFIG"]["ENTRY"].value)
        for card in ccf_head.cards:
            if not card.keyword or _is_structural(card.keyword):
                continue
            entries.append(f"ccf:{card.keyword} = {card.value}")
        self.set_data("DRP_CONFIG", pd.DataFrame({"ENTRY": entries}))
        describe("DRP_CONFIG", "Pipeline details (settings etc) to go from native data to L4")
        describe("EXT_DESCRIPT", "Table describing contents of each extension")

        # Per-order RVs, matched to the fitspec order table by echelle order.
        order_table = l2.data["ORDER_TABLE"]
        l1_orders = order_table["ECHELLE_ORDER"].value.astype(int)
        ccf_orders = per_order["orders"].astype(int)
        rows = np.array([int(np.where(l1_orders == o)[0][0]) for o in ccf_orders])
        n = len(ccf_orders)
        rv1 = OrderedDict(
            {
                "BJD_TDB": np.full(n, bjd_tdb, dtype=np.float64),
                "RV": per_order["v"].astype(np.float64) / CM_PER_KM,
                "RV_ERR": per_order["e_v"].astype(np.float64) / CM_PER_KM,
                "BERV": np.full(n, berv_kms, dtype=np.float64),
                "WAVE_START": order_table["WAVE_START"].value[rows],
                "WAVE_END": order_table["WAVE_END"].value[rows],
                "ORDER_INDEX": rows,
                "ECHELLE_ORDER": ccf_orders,
            }
        )
        self.set_data("RV1", pd.DataFrame(rv1))
        self.set_header(
            "RV1",
            fits.Header(
                {
                    "RVMETHOD": "CCF",
                    "SKYRMVD": False,
                    "TELLRMVD": bool(ccf_head.get("DIV_TELL", False)),
                }
            ),
        )
        describe("RV1", "Order-wise CCF RVs for the EXPRES science trace (km/s)")

        vgrid = combined["V_grid"].astype(np.float64) / CM_PER_KM
        ccf_header = OrderedDict(
            {
                "CCFSTART": float(vgrid[0]),
                "CCFSTEP": float(vgrid[1] - vgrid[0]),
                "VELNSTEP": int(len(vgrid)),
                "CCFMASK": str(ccf_head["MASK"]),
            }
        )
        self.create_extension(
            "CCF1",
            "ImageHDU",
            data=per_order["ccfs"].astype(np.float64),
            header=ccf_header,
        )
        describe("CCF1", "Order-wise CCFs for the EXPRES science trace, same row order as RV1")

        self.create_extension(
            "CUSTOM_CCF1",
            "ImageHDU",
            data=combined["ccf"].astype(np.float64)[np.newaxis, :],
            header=OrderedDict(ccf_header),
        )
        describe("CUSTOM_CCF1", "Combined (order-summed) CCF from which the PRIMARY RV was derived")
        custom_rv = OrderedDict(
            {
                "BJD_TDB": np.array([bjd_tdb], dtype=np.float64),
                "RV": np.array([phead["RV"]], dtype=np.float64),
                "RV_ERR": np.array([phead["RVERR"]], dtype=np.float64),
                "BERV": np.array([berv_kms], dtype=np.float64),
                "WAVE_START": np.array([np.min(rv1["WAVE_START"])]),
                "WAVE_END": np.array([np.max(rv1["WAVE_END"])]),
            }
        )
        self.create_extension(
            "CUSTOM_RV1",
            "BinTableHDU",
            data=pd.DataFrame(custom_rv),
            header=fits.Header(
                {"RVMETHOD": "CCF", "SKYRMVD": False, "TELLRMVD": bool(ccf_head.get("DIV_TELL", False))}
            ),
        )
        describe("CUSTOM_RV1", "RV from the combined CCF (same value as PRIMARY RV/RVERR)")

        diag = {"metric_name": [], "value": [], "uncertainty": []}
        for card, (name, err_card) in ACTIVITY_CARDS.items():
            if card not in l1_head1:
                continue
            diag["metric_name"].append(name)
            diag["value"].append(float(l1_head1[card]))
            diag["uncertainty"].append(
                float(l1_head1[err_card]) if err_card and err_card in l1_head1 else np.nan
            )
        for card, name in CCF_CARDS.items():
            if card not in ccf_head:
                continue
            diag["metric_name"].append(name)
            diag["value"].append(float(ccf_head[card]))
            diag["uncertainty"].append(np.nan)
        self.create_extension("DIAGNOSTICS1", "BinTableHDU", data=pd.DataFrame(diag))
        describe("DIAGNOSTICS1", "Activity indicators and CCF quality metrics from the EXPRES pipeline (native units)")

        self.set_data("EXT_DESCRIPT", pd.DataFrame(ext_table))
