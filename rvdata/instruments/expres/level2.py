"""
RVData Level 2 reader for EXPRES.

The native input is an EXPRES "fitspec" file (optimally extracted spectrum,
one FITS file per exposure, produced by the EXPRES pipeline described in
Petersburg et al. 2020, AJ 159, 187):

* HDU 0 ``PRIMARY``: observation header (no data).
* HDU 1 ``optimal``: BinTable with one row per echelle order and the columns
  ``spectrum``, ``uncertainty``, ``blaze``, ``wavelength`` (vacuum Angstrom),
  ``bary_wavelength``, ``pixel_mask``, ``tellurics``, ``orders`` (absolute
  echelle order) and more. Its header carries the pipeline version, the
  barycentric MJD (``BARYMJD``, TDB, at the Sun for solar data) and the
  activity indicators.
* HDU 2 ``EXPOSURE METER + BARY_CORR``: chromatic exposure meter time series
  and the barycentric correction inputs/outputs.

Every standard keyword is filled from a native value or a verified constant;
anything the native file does not carry is left undefined rather than guessed.
"""

import os
import warnings
from collections import OrderedDict
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from astropy import constants as const
from astropy.io import fits
from astropy.time import Time

import rvdata
from rvdata.core.models.level2 import RV2

_CONFIG_DIR = os.path.join(os.path.dirname(__file__), "config")

# Instrument eras: version tag and UT start date of each permanent change.
expres_epochs, epoch_start_isot = np.loadtxt(
    os.path.join(_CONFIG_DIR, "expres_epochs.csv"),
    delimiter=",",
    skiprows=1,
    dtype=str,
    encoding="utf-8-sig",
).T
epoch_start_mjd = Time(epoch_start_isot).mjd

# standard keyword -> native keyword (blank = filled in code or undefined)
header_map = (
    pd.read_csv(os.path.join(_CONFIG_DIR, "expres_header_map.csv"))
    .fillna("")
    .set_index("standard")
)

# Native OBSTYPE -> standard OBSTYPE
obstype_map = {
    "Science": "Sci",
    "Solar": "Sci",
    "Calibration": "Cal",
    "Dark": "Cal",
    "ThAr": "Cal",
    "Quartz": "Cal",
    "LFC": "Cal",
}

# FITS structural cards that are not pipeline configuration
_STRUCTURAL_PREFIXES = (
    "XTENSION", "BITPIX", "NAXIS", "PCOUNT", "GCOUNT", "TFIELDS",
    "TTYPE", "TFORM", "TDIM", "TUNIT", "EXTNAME", "SIMPLE", "EXTEND",
    "COMMENT", "HISTORY",
)


def _is_structural(key):
    return key == "" or any(key.startswith(p) for p in _STRUCTURAL_PREFIXES)


def drp_flag(hdul):
    """Pass if every column the translator needs is present in HDU 1."""
    needed = {
        "spectrum", "uncertainty", "blaze", "wavelength", "bary_wavelength",
        "pixel_mask", "tellurics", "orders",
    }
    return "Pass" if needed.issubset(set(hdul[1].columns.names)) else "Fail"


def instrument_era(mjd):
    """INSTERA tag for an observation at the given MJD, or None if the MJD
    predates the first recorded era (UNDEFINED rather than the newest era)."""
    if mjd < epoch_start_mjd[0]:
        return None
    return str(expres_epochs[np.sum(mjd >= epoch_start_mjd) - 1])


class EXPRESRV2(RV2):
    """
    Data model and reader for RVData Level 2 data constructed from an EXPRES
    fitspec file.

    Parameters
    ----------
    Inherits all parameters from :class:`RV2`.

    Notes
    -----
    Use the classmethod ``from_fits``:

    >>> from rvdata.core.models.level2 import RV2
    >>> l2 = RV2.from_fits("fitspec/Sun_20240126.5073.fits", instrument="EXPRES")
    >>> l2.to_fits()
    """

    def _read(self, hdul: fits.HDUList, **kwargs) -> None:
        head0 = hdul[0].header
        head1 = hdul[1].header
        head2 = hdul[2].header
        data = hdul[1].data

        raw_obstype = str(head0["OBSTYPE"]).strip()
        is_solar = str(head0["OBJECT"]).strip() == "Sun" or raw_obstype == "Solar"
        obstype = obstype_map.get(raw_obstype)
        if obstype is None:
            warnings.warn(
                f"EXPRES OBSTYPE {raw_obstype!r} not recognised; "
                "OBSTYPE left undefined"
            )
        jd_utc = Time(head0["DATE-SHT"], format="isot", scale="utc").jd
        mid_mjd = Time(head0["MIDPOINT"], format="isot", scale="utc").mjd
        instflag = "Pass" if bool(head0.get("EXPMTR", False)) else "Fail"
        drpflag = drp_flag(hdul)
        version = rvdata.__version__ or ""

        # Keywords that are computed rather than copied.
        computed = {
            "ORGANIZA": "Yale",
            "DATALVL": "L2",
            "OBSTYPE": obstype,
            "BINNING": str(head0["CCDBIN"]).replace(" ", "")[1:-1].replace(",", "x"),
            "NUMTRACE": 1,
            "NUMORDER": int(head1["NAXIS2"]),
            "FILENAME": "",  # overwritten by to_fits with the standard name
            "DATE": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3],
            "JD_UTC": jd_utc,
            "TRACE1": "SCI" if obstype == "Sci" else "CAL",
            "CLSRC1": None if obstype == "Sci" else str(head0["OBJECT"]).strip(),
            "CSRC1": "SOLAR SYSTEM" if is_solar else None,
            "CID1": "Sun" if is_solar else None,
            # REQRA/REQDEC for the Sun are the apparent pointing at the
            # observation, so the catalog epoch is the observation epoch.
            "CEPCH1": Time(head0["MIDPOINT"], format="isot").decimalyear if is_solar else None,
            "CRV1": 0.0 if is_solar else None,
            "OBSERVAT": "Lowell Observatory",
            "NUMTEL": 1,
            # barycorrpy (HDU 2 kwargs LAT/LONGI/ALT) builds the site with
            # astropy EarthLocation.from_geodetic, whose ellipsoid is WGS84.
            "GEOSYS": "WGS84",
            "OBSLON": float(head2.get("LONGI", head0["SITELONG"])),
            "OBSLAT": float(head2.get("LAT", head0["SITELAT"])),
            "OBSALT": float(head2.get("ALT", head0["SITEELEV"])),
            "ISSOLAR": is_solar,
            "DRPTAG": str(head1["VERSION"]),
            "EPRVTAG": f"v{version}" if version else None,
            "VOCLASS": f"EPRVSTANDARDv{version}" if version else None,
            "INSTERA": instrument_era(mid_mjd),
            "EXTRACT": str(head1.get("EXTNAME", "")) or None,
            "FULLCOMP": "Yes",
            "INSTFLAG": instflag,
            "DRPFLAG": drpflag,
            "SUMMFLAG": "Pass" if (instflag == "Pass" and drpflag == "Pass") else "Fail",
            "DQLVL0": 0,
            "DQLVL1": 0,
            "DQLVL2": 0,
        }

        # A real fits.Header (not a plain dict), following the NEID reader's
        # style: RVDataModel.read() recasts PRIMARY header values by keyword
        # after _read() returns, assigning (value, comment) tuples back into
        # self.headers["PRIMARY"]. A fits.Header interprets that assignment
        # as value+comment and still returns a plain scalar on lookup; a
        # plain dict would instead store the literal tuple.
        standard_head = fits.PrimaryHDU().header
        for key in header_map.index:
            native_key = header_map.loc[key, "expres"]
            required = header_map.loc[key, "required"] == "Y"
            if key in computed:
                value = computed[key]
            elif native_key and native_key in head0:
                value = head0[native_key]
            else:
                value = None
            if value is None and not required:
                continue
            standard_head[key] = value
        self.set_header("PRIMARY", standard_head)

        ext_table = {"Name": [], "Description": []}

        def describe(name, description):
            ext_table["Name"].append(name)
            ext_table["Description"].append(description)

        describe("PRIMARY", "EPRV Standard FITS HEADER (no data)")

        self.set_header("INSTRUMENT_HEADER", head0)
        describe("INSTRUMENT_HEADER", "Inherited EXPRES fitspec primary header (no data)")
        describe("RECEIPT", "Table of operations that have been performed on this file")

        # DRP_CONFIG: every non-structural card of the extraction and exposure
        # meter headers, prefixed by the native extension it came from.
        entries = []
        for prefix, header in (("optimal", head1), ("expmeter", head2)):
            for card in header.cards:
                if not card.keyword or _is_structural(card.keyword):
                    continue
                entries.append(f"{prefix}:{card.keyword} = {card.value}")
        self.set_data("DRP_CONFIG", pd.DataFrame({"ENTRY": entries}))
        describe("DRP_CONFIG", "Pipeline details (settings etc) to go from native data to L2")
        describe("EXT_DESCRIPT", "Table describing contents of each extension")

        wave = data["wavelength"].astype(np.float64)
        spectrum = data["spectrum"].astype(np.float64)
        uncertainty = data["uncertainty"].astype(np.float64)
        blaze = data["blaze"].astype(np.float64)
        bary_wave = data["bary_wavelength"].astype(np.float64)
        pixel_mask = data["pixel_mask"].astype(bool)
        orders = data["orders"].astype(int)

        self.set_data(
            "ORDER_TABLE",
            pd.DataFrame(
                {
                    "ECHELLE_ORDER": orders,
                    "ORDER_INDEX": np.arange(len(orders)),
                    "WAVE_START": np.nanmin(wave, axis=1),
                    "WAVE_END": np.nanmax(wave, axis=1),
                }
            ),
        )
        describe("ORDER_TABLE", "Table capturing the wavelength extent of each order in Trace 1")

        # Native "spectrum" is the extracted flux divided by the blaze, and
        # "uncertainty" is its 1-sigma error on the same scale.
        self.set_data("TRACE1_FLUX", spectrum * blaze)
        describe("TRACE1_FLUX", "Extracted flux in trace 1 (native spectrum x blaze)")
        self.set_data("TRACE1_WAVE", wave)
        describe("TRACE1_WAVE", "Vacuum wavelength solution for trace 1 (Angstrom)")
        self.set_data("TRACE1_VAR", (uncertainty * blaze) ** 2)
        describe("TRACE1_VAR", "Variance of TRACE1_FLUX")
        self.set_data("TRACE1_BLAZE", blaze)
        blaze_head = fits.Header()
        blaze_head["BLZNORM"] = (False, "EXPRES blaze is unnormalized e- counts")
        self.set_header("TRACE1_BLAZE", blaze_head)
        describe("TRACE1_BLAZE", "Blaze function for trace 1")

        self.create_extension(
            "TRACE1_QUALITY", "ImageHDU", data=(~pixel_mask).astype(np.uint8)
        )
        describe("TRACE1_QUALITY", "Pixel quality for trace 1: 0 = good, 1 = masked by the pipeline")

        # Barycentric correction as applied by the pipeline:
        # bary_wavelength = wavelength * (1 + z)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            bary_z = bary_wave / wave - 1.0
        self.set_data("BARYCORR_KMS", bary_z * const.c.to("km/s").value)
        describe(
            "BARYCORR_KMS",
            "Barycentric correction per pixel in km/s (lambda_bary = lambda * "
            "(1 + v/c)); for solar data this includes the solar gravitational "
            "redshift, not a purely kinematic velocity",
        )
        self.set_data("BARYCORR_Z", bary_z)
        describe("BARYCORR_Z", "Barycentric correction per pixel as redshift z")

        # BARYMJD is the barycentric (TDB) photon-weighted midpoint; for solar
        # data it is the emission time at the Sun (light-travel corrected).
        self.set_data(
            "BJD_TDB", np.array([head1["BARYMJD"] + 2400000.5], dtype=np.float64)
        )
        describe("BJD_TDB", "Photon-weighted midpoint, BJD_TDB (at the Sun for solar data)")

        # Exposure meter: one row per time sample, one column per wavelength.
        expm = hdul[2].data
        expm_counts = np.array([row["expm_specs"] for row in expm], dtype=np.float64).T
        expm_times = expm["midpoints"].astype(np.float64)
        expm_waves = np.array(expm["wavelengths"][0], dtype=np.float64)
        expm_table = OrderedDict({"TIME": expm_times})
        for wl, counts in zip(expm_waves, expm_counts):
            expm_table[f"{wl:.3f}"] = counts
        # Strip the native HDU2 structural cards (TFIELDS, TTYPE1, EXTNAME,
        # etc.) so the in-memory header does not misdescribe the EXPMETER
        # table, which has different columns than the native extension.
        expm_head = fits.Header()
        for card in head2.cards:
            if card.keyword and not _is_structural(card.keyword):
                expm_head.append(card)
        self.create_extension(
            "EXPMETER", "BinTableHDU", header=expm_head, data=pd.DataFrame(expm_table)
        )
        describe("EXPMETER", "Chromatic exposure meter counts; TIME in seconds from exposure start, one column per wavelength (nm)")

        self.create_extension(
            "TRACE1_TELLURIC", "ImageHDU", data=data["tellurics"].astype(np.float64)
        )
        describe("TRACE1_TELLURIC", "SELENITE telluric model for trace 1 (unphysical where TRACE1_QUALITY = 1)")

        self.set_data("EXT_DESCRIPT", pd.DataFrame(ext_table))
