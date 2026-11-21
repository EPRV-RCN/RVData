"""HARPS and HARPS-N level 3 and level 4 translators, on small synthetic DRS files.

Both instruments run the ESPRESSO-derived DRS, so their S1D and CCF products
have ESPRESSO's layout. HARPS keeps the "HIERARCH ESO" keywords; HARPS-N writes
the same QC and RV keywords as "HIERARCH TNG" (checked on a DRS 3.0.1 header).
"""

import importlib
import importlib.resources

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from rvdata.core.models.level3 import RV3
from rvdata.core.models.level4 import RV4

INSTRUMENTS = {"harps": "ESO", "harpsn": "TNG"}
# Real ESO archive names carry colons (HARPS.2023-01-01T00:00:00.000.fits),
# but Windows forbids ':' in filenames and these tests write the derived
# products to disk. The L3/L4 name derivation is pure string concatenation,
# so a colon-free stand-in exercises the same code.
RAW = "HARPS.2023-01-01T00-00-00.000.fits"
NPIX, NVEL, NORDER = 6, 5, 3


def _utils(inst):
    return importlib.import_module(f"rvdata.instruments.{inst}.utils")


def _write_s1d(path, scale):
    flux = scale * np.arange(100.0, 100.0 + NPIX)
    table = Table(
        {
            "wavelength": np.linspace(5000.0, 5001.0, NPIX),
            "flux": flux,
            "error": np.sqrt(flux),  # the DRS writes 1-sigma errors
        }
    )
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(table)]).writeto(path)
    return np.sqrt(flux)


def _ccf_header(kw, rv=-10.5):
    h = fits.Header()
    values = {
        "QC BJD": 2459945.5,
        "QC BERV": 12.3,
        "QC CCF RV": rv,
        "QC CCF RV ERROR": 0.0012,
        "RV START": -20.0,
        "RV STEP": 0.5,
        "QC CCF MASK": "G2",
        "QC CCF FWHM": 7.1,
        "QC CCF FWHM ERROR": 0.01,
        "QC CCF CONTRAST": 45.0,
        "QC CCF CONTRAST ERROR": 0.1,
        "QC CCF BIS SPAN": 0.02,
        "QC CCF BIS SPAN ERROR": 0.003,
    }
    for key, value in values.items():
        h[f"HIERARCH {kw} {key}"] = (value, key.lower())
    return h


def _write_ccf(path, kw, rv=-10.5):
    data = np.ones((NORDER + 1, NVEL))
    fits.HDUList(
        [
            fits.PrimaryHDU(header=_ccf_header(kw, rv)),
            fits.ImageHDU(data, name="SCIDATA"),
            fits.ImageHDU(data / 10, name="ERRDATA"),
        ]
    ).writeto(path)


def _product(tmp_path, suffix):
    return str(tmp_path / ("r." + RAW[:-5] + f"_{suffix}.fits"))


@pytest.mark.parametrize("inst", INSTRUMENTS)
@pytest.mark.parametrize(
    "level, keys",
    [
        (3, {"s1d_A": "S1D_A", "s1d_B": "S1D_B",
             "s1d_skysub_A": "S1D_SKYSUB_A", "s1d_tell_corr_A": "S1D_TELL_CORR_A"}),
        (4, {"ccf_file": "CCF_A", "ccf_file_B": "CCF_B",
             "ccf_tel_corr_file": "CCF_TELL_CORR_A", "ccf_sky_file": "CCF_SKYSUB_A"}),
    ],
)  # fmt: skip
def test_get_files_names_finds_level_products(tmp_path, inst, level, keys):
    raw = str(tmp_path / RAW)
    names = _utils(inst).get_files_names(raw, "standard", level=level)
    assert names["raw_file"] == raw
    for key, suffix in keys.items():
        assert names[key] == _product(tmp_path, suffix)


def _l3_required_keywords():
    config = importlib.resources.files("rvdata.core.models.config")
    table = pd.read_csv(config / "L3-PRIMARY-keywords.csv")
    return table.loc[table["Required"].astype(bool), "Keyword"].tolist()


@pytest.mark.parametrize("inst", INSTRUMENTS)
def test_convert_s1d(tmp_path, inst):
    names, errors = {}, {}
    for key, suffix, scale in [
        ("s1d_A", "S1D_A", 1.0),
        ("s1d_B", "S1D_B", 2.0),
        ("s1d_tell_corr_A", "S1D_TELL_CORR_A", 4.0),
    ]:
        names[key] = _product(tmp_path, suffix)
        errors[key] = _write_s1d(names[key], scale)
    names["s1d_skysub_A"] = _product(tmp_path, "S1D_SKYSUB_A")  # not written

    rv3 = RV3()
    # create_PRIMARY runs first and replaces PRIMARY with the level 2 header.
    rv3.set_header("PRIMARY", fits.Header({"DATALVL": "L3"}))
    _utils(inst).convert_S1D(rv3, names)

    primary = rv3.headers["PRIMARY"]
    assert [k for k in _l3_required_keywords() if k not in primary] == []
    assert primary["FLXNRMMD"] == "None"
    np.testing.assert_allclose(rv3.data["STITCHED_CORR_SCI_VAR"], errors["s1d_A"] ** 2)
    np.testing.assert_allclose(
        rv3.data["STITCHED_TELLCORR_SCI_VAR"], errors["s1d_tell_corr_A"] ** 2
    )
    assert "STITCHED_SKYSUB_SCI_FLUX" not in rv3.extensions


@pytest.mark.parametrize("inst, kw", INSTRUMENTS.items())
def test_convert_ccf_reads_instrument_keywords(tmp_path, inst, kw):
    names = {
        "ccf_file": _product(tmp_path, "CCF_A"),
        "ccf_file_B": _product(tmp_path, "CCF_B"),  # not written
        "ccf_tel_corr_file": _product(tmp_path, "CCF_TELL_CORR_A"),
        "ccf_sky_file": _product(tmp_path, "CCF_SKYSUB_A"),  # not written
    }
    _write_ccf(names["ccf_file"], kw)
    _write_ccf(names["ccf_tel_corr_file"], kw, rv=-10.4)

    rv4 = RV4()
    rv4.set_header("PRIMARY", fits.Header({"DATALVL": "L4"}))  # as create_PRIMARY
    _utils(inst).convert_CCF(rv4, names)

    primary = rv4.headers["PRIMARY"]
    assert (primary["RV"], primary["RVERR"], primary["BERV"]) == (-10.5, 0.0012, 12.3)
    assert primary["BJDTDB"] == 2459945.5

    rv1 = rv4.data["RV1"]
    assert list(rv1["RV"]) == [-10.5]
    assert list(rv4.data["RV_TEL_CORR"]["RV"]) == [-10.4]
    assert "RV2" not in rv4.extensions

    ccf = rv4.headers["CCF1"]
    assert (ccf["VELSTART"], ccf["VELSTEP"], ccf["VELNSTEP"]) == (-20.0, 0.5, NVEL)
    assert ccf["CCFMASK"] == "G2"

    diag = rv4.data["DIAGNOSTICS1"]
    assert dict(zip(diag["metric_name"], diag["value"])) == {
        "FWHM": 7.1, "CONTRAST": 45.0, "BIS SPAN": 0.02,
    }  # fmt: skip
