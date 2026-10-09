"""Regression tests for the EXPRES translator (L2, L3, L4).

The fixture is one EXPRES solar exposure (2024-01-26, exposure 5073) as the
native "fitspec" (extracted spectrum) and "ccf" (RV) products. It is
downloaded from the RVData fixture host, or read from a local directory
named by the EXPRES_FIXTURE_DIR environment variable with the same layout:

    $EXPRES_FIXTURE_DIR/fitspec/Sun_20240126.5073.fits
    $EXPRES_FIXTURE_DIR/ccf/Sun_20240126.5073.fits
"""
import os

import numpy as np
import pytest
import requests
from astropy import constants as const
from astropy.io import fits

from rvdata.core.models.base import RVDataModel
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level3 import RV3
from rvdata.core.models.level4 import RV4
from rvdata.instruments.expres.level2 import EXPRESRV2
from rvdata.tests.regression.compliance import (
    check_l2_compliance,
    check_l3_compliance,
    check_l4_compliance,
)

FIXTURE_BASENAME = "Sun_20240126.5073.fits"
FIXTURE_URL_ROOT = "http://grinnell.as.arizona.edu/~rvdata/expres"
LOCAL_FIXTURE_DIR = "expres_fixtures"


def download_file(url, filename):
    response = requests.get(url)
    response.raise_for_status()
    with open(filename, "wb") as file:
        file.write(response.content)


def expres_fixture_files():
    """Return (fitspec_path, ccf_path), downloading if needed."""
    root = os.environ.get("EXPRES_FIXTURE_DIR", LOCAL_FIXTURE_DIR)
    paths = {}
    for kind in ("fitspec", "ccf"):
        path = os.path.join(root, kind, FIXTURE_BASENAME)
        if not os.path.exists(path):
            os.makedirs(os.path.dirname(path), exist_ok=True)
            download_file(f"{FIXTURE_URL_ROOT}/{kind}/{FIXTURE_BASENAME}", path)
        paths[kind] = path
    return paths["fitspec"], paths["ccf"]


# ---------------------------------------------------------------- Level 2


def test_expres_l2_compliance():
    fitspec, _ = expres_fixture_files()
    l2 = RV2.from_fits(fitspec, instrument="EXPRES")
    out = l2.to_fits()
    base = os.path.basename(out)
    assert RVDataModel.FILENAME_PATTERN.match(base), base
    assert base.startswith("expres_SL2_"), base
    check_l2_compliance(out)


def test_expres_l2_primary_header_values():
    fitspec, _ = expres_fixture_files()
    l2 = RV2.from_fits(fitspec, instrument="EXPRES")
    h = l2.headers["PRIMARY"]
    with fits.open(fitspec) as hdul:
        h0 = hdul[0].header
        h2 = hdul[2].header
    assert h["INSTRUME"] == "EXPRES"
    assert h["DATALVL"] == "L2"
    assert h["OBSTYPE"] == "Sci"
    assert h["ISSOLAR"] is True
    assert h["NUMTRACE"] == 1
    assert h["TRACE1"] == "SCI"
    assert h["CLSRC1"] is None
    assert h["CSRC1"] == "SOLAR SYSTEM"
    assert h["CID1"] == "Sun"
    assert h["CRV1"] == 0.0
    assert h["NUMTEL"] == 1
    assert h["TELESCOP"] == h0["TELESCP"]
    assert h["TELEID1"] == h0["TELESCP"]
    assert h["GEOSYS"] == "WGS84"
    assert h["OBSLON"] == h2["LONGI"]
    assert h["OBSLAT"] == h2["LAT"]
    assert h["OBSALT"] == h2["ALT"]
    assert h["BINNING"] == "1x1"
    assert h["NUMORDER"] == 86
    assert h["DATE-OBS"] == h0["DATE-SHT"]
    assert h["EXPTIME"] == float(h0["AEXPTIME"])
    assert h["INSTERA"] == "5.10.0"
    assert h["DRPTAG"] == "0.4.1"
    assert h["EXTRACT"] == "optimal"
    assert h["SUMMFLAG"] in ("Pass", "Fail")
    assert h["DQLVL0"] == 0 and h["DQLVL1"] == 0 and h["DQLVL2"] == 0
    assert h["FULLCOMP"] == "Yes"


def test_expres_l2_flux_var_wave_blaze():
    fitspec, _ = expres_fixture_files()
    l2 = RV2.from_fits(fitspec, instrument="EXPRES")
    with fits.open(fitspec) as hdul:
        d = hdul[1].data
        spectrum = d["spectrum"].astype(np.float64)
        unc = d["uncertainty"].astype(np.float64)
        blaze = d["blaze"].astype(np.float64)
        wave = d["wavelength"].astype(np.float64)
    flux = l2.data["TRACE1_FLUX"]
    var = l2.data["TRACE1_VAR"]
    assert flux.shape == (86, 7920)
    np.testing.assert_allclose(flux, spectrum * blaze, equal_nan=True)
    np.testing.assert_allclose(var, (unc * blaze) ** 2, equal_nan=True)
    np.testing.assert_array_equal(l2.data["TRACE1_BLAZE"], blaze)
    assert l2.data["TRACE1_WAVE"].dtype == np.float64
    np.testing.assert_array_equal(l2.data["TRACE1_WAVE"], wave)


def test_expres_l2_quality_and_order_table():
    fitspec, _ = expres_fixture_files()
    l2 = RV2.from_fits(fitspec, instrument="EXPRES")
    with fits.open(fitspec) as hdul:
        orders = hdul[1].data["orders"].astype(int)
        pixel_mask = hdul[1].data["pixel_mask"].astype(bool)
        wave = hdul[1].data["wavelength"].astype(np.float64)
    ot = l2.data["ORDER_TABLE"]
    np.testing.assert_array_equal(ot["ECHELLE_ORDER"].value, orders)
    assert ot["ECHELLE_ORDER"].value[0] == 160
    assert ot["ECHELLE_ORDER"].value[-1] == 75
    np.testing.assert_array_equal(ot["ORDER_INDEX"].value, np.arange(86))
    assert np.all(np.isfinite(ot["WAVE_START"].value))
    assert np.all(np.isfinite(ot["WAVE_END"].value))
    np.testing.assert_allclose(ot["WAVE_START"].value, wave.min(axis=1))
    np.testing.assert_allclose(ot["WAVE_END"].value, wave.max(axis=1))
    quality = l2.data["TRACE1_QUALITY"]
    assert quality.dtype == np.uint8
    np.testing.assert_array_equal(quality, (~pixel_mask).astype(np.uint8))
    # order 75 is entirely masked on this exposure
    assert quality[-1].all()


def test_expres_l2_time_and_barycentric():
    fitspec, _ = expres_fixture_files()
    l2 = RV2.from_fits(fitspec, instrument="EXPRES")
    with fits.open(fitspec) as hdul:
        barymjd = hdul[1].header["BARYMJD"]
        wtd_mdpt = hdul[2].header["wtd_mdpt"]
        z_hdr = hdul[2].header["wtd_mdpt_bc"]
        wave = hdul[1].data["wavelength"].astype(np.float64)
        bwave = hdul[1].data["bary_wavelength"].astype(np.float64)
    bjd = l2.data["BJD_TDB"]
    assert bjd.dtype == np.float64
    assert bjd.shape == (1,)
    assert abs(bjd[0] - (barymjd + 2400000.5)) < 1e-9
    # BJD_TDB at the Sun is ~422 s before the UTC midpoint on this date
    # (TDB-UTC = +69.2 s, light travel Sun->Earth = -491.3 s)
    assert -430 < (bjd[0] - wtd_mdpt) * 86400 < -415
    z = l2.data["BARYCORR_Z"]
    v = l2.data["BARYCORR_KMS"]
    assert z.shape == wave.shape and v.shape == wave.shape
    np.testing.assert_allclose(z, bwave / wave - 1, equal_nan=True, rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        v, z * const.c.to("km/s").value, equal_nan=True, rtol=1e-12
    )
    # sign convention: bary_wavelength = wavelength * (1 + z); on this
    # exposure z is negative (about -0.59 km/s) and matches the header
    assert abs(np.nanmedian(z) - z_hdr) < 2e-8
    assert np.nanmedian(v) < 0
    assert h_float(l2.headers["PRIMARY"]["JD_UTC"]) == pytest.approx(
        _jd(fitspec), abs=1e-9
    )


def h_float(x):
    return float(x[0] if isinstance(x, tuple) else x)


def _jd(fitspec):
    from astropy.time import Time

    with fits.open(fitspec) as hdul:
        return Time(hdul[0].header["DATE-SHT"], format="isot", scale="utc").jd


def test_expres_l2_ext_descript_and_drp_config():
    fitspec, _ = expres_fixture_files()
    l2 = RV2.from_fits(fitspec, instrument="EXPRES")
    ext = l2.data["EXT_DESCRIPT"]
    assert list(ext.columns) == ["Name", "Description"]
    for name in (
        "PRIMARY", "INSTRUMENT_HEADER", "RECEIPT", "DRP_CONFIG", "EXT_DESCRIPT",
        "ORDER_TABLE", "TRACE1_FLUX", "TRACE1_WAVE", "TRACE1_VAR", "TRACE1_BLAZE",
        "TRACE1_QUALITY", "BARYCORR_KMS", "BARYCORR_Z", "BJD_TDB", "EXPMETER",
        "TRACE1_TELLURIC",
    ):
        assert name in ext["Name"].value, name
    drp = l2.data["DRP_CONFIG"]
    assert list(drp.columns) == ["ENTRY"]
    entries = list(drp["ENTRY"].value)
    assert any(e.startswith("optimal:VERSION = 0.4.1") for e in entries)
    assert any(e.startswith("optimal:BARYMJD = ") for e in entries)
    assert any(e.startswith("expmeter:SolSystemTarget = Sun") for e in entries)
    assert not any(e.startswith("optimal:TTYPE") for e in entries)


def _l2_from_patched_header(fitspec, **cards):
    """Read the fixture with PRIMARY cards replaced, without touching disk."""
    with fits.open(fitspec, memmap=False) as hdul:
        for key, value in cards.items():
            hdul[0].header[key] = value
        l2 = EXPRESRV2()
        l2.read(hdul, instrument="EXPRES")
    return l2


def test_expres_l2_stellar_branch():
    fitspec, _ = expres_fixture_files()
    l2 = _l2_from_patched_header(fitspec, OBJECT="HD 10700", OBSTYPE="Science")
    h = l2.headers["PRIMARY"]
    assert h["ISSOLAR"] is False
    assert h["OBSTYPE"] == "Sci"
    assert h["TRACE1"] == "SCI"
    assert h["CLSRC1"] is None
    assert h["OBJECT"] == "HD 10700"
    assert h["CSRC1"] is None
    assert h["CID1"] is None
    assert h["CEPCH1"] is None
    assert h["CRV1"] is None


def test_expres_l2_cal_branch():
    fitspec, _ = expres_fixture_files()
    l2 = _l2_from_patched_header(fitspec, OBJECT="ThAr", OBSTYPE="ThAr")
    h = l2.headers["PRIMARY"]
    assert h["OBSTYPE"] == "Cal"
    assert h["TRACE1"] == "CAL"
    assert h["CLSRC1"] == "ThAr"
    assert h["ISSOLAR"] is False


# ---------------------------------------------------------------- Level 3


def test_expres_l3_compliance():
    fitspec, _ = expres_fixture_files()
    l3 = RV3.from_fits(fitspec, instrument="EXPRES")
    out = l3.to_fits()
    base = os.path.basename(out)
    assert RVDataModel.FILENAME_PATTERN.match(base), base
    assert base.startswith("expres_SL3_"), base
    check_l3_compliance(out)


def test_expres_l3_stitched_range():
    fitspec, _ = expres_fixture_files()
    l3 = RV3.from_fits(fitspec, instrument="EXPRES")
    wave = l3.data["STITCHED_CORR_SCI_WAVE"]
    flux = l3.data["STITCHED_CORR_SCI_FLUX"]
    var = l3.data["STITCHED_CORR_SCI_VAR"]
    assert wave.dtype == np.float64
    assert wave.shape == flux.shape == var.shape
    assert wave.min() >= 3800.0 and wave.max() <= 8100.0
    assert np.all(np.diff(wave) > 0)
    finite = np.isfinite(flux)
    assert finite.mean() > 0.5
    # the red end (order 76, ~8000-8100 A) stitches without raising even
    # though order 75 has no finite flux at all
    assert np.isfinite(flux[(wave > 7900) & (wave < 8000)]).mean() > 0.3
    assert l3.headers["PRIMARY"]["DATALVL"] == "L3"
    assert l3.headers["PRIMARY"]["INSTRUME"] == "EXPRES"


def test_expres_l3_does_not_mutate_l2():
    fitspec, _ = expres_fixture_files()
    with fits.open(fitspec, memmap=False) as hdul:
        l2 = EXPRESRV2()
        l2.read(hdul, instrument="EXPRES")
    keys_before = set(l2.headers["PRIMARY"].keys())
    l3 = RV3()
    l3.convert_level2_to_level3(l2)
    assert l2.headers["PRIMARY"]["DATALVL"] == "L2"
    assert set(l2.headers["PRIMARY"].keys()) == keys_before
    assert l3.headers["PRIMARY"] is not l2.headers["PRIMARY"]
