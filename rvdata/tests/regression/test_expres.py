"""Regression tests for the EXPRES translator (L2, L3, L4).

The fixture is one EXPRES solar exposure (2024-01-26, exposure 5073) as the
native "fitspec" (extracted spectrum) and "ccf" (RV) products. It is
downloaded from the RVData fixture host, or read from a local directory
named by the EXPRES_FIXTURE_DIR environment variable with the same layout:

    $EXPRES_FIXTURE_DIR/fitspec/Sun_20240126.5073.fits
    $EXPRES_FIXTURE_DIR/ccf/Sun_20240126.5073.fits
"""
import os
import re

import numpy as np
import pytest
import requests
from astropy import constants as const
from astropy.io import fits

from rvdata.core.models.base import RVDataModel
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level3 import RV3
from rvdata.core.models.level4 import RV4
from rvdata.instruments.expres.level2 import EXPRESRV2, instrument_era
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
    assert instrument_era(50000.0) is None
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
    assert l2.headers["TRACE1_BLAZE"]["BLZNORM"] is False
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
    assert np.isfinite(flux[(wave > 7900) & (wave < 8000)]).mean() > 0.6
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


# ---------------------------------------------------------------- Level 4


def test_expres_l4_compliance():
    _, ccf = expres_fixture_files()
    l4 = RV4.from_fits(ccf, instrument="EXPRES")
    out = l4.to_fits()
    base = os.path.basename(out)
    assert RVDataModel.FILENAME_PATTERN.match(base), base
    assert base.startswith("expres_SL4_"), base
    check_l4_compliance(out)


def test_expres_l4_primary_and_rv1():
    fitspec, ccf = expres_fixture_files()
    l4 = RV4.from_fits(ccf, instrument="EXPRES")
    with fits.open(ccf) as c:
        ch = c[0].header
        per_order = c[2].data
        vgrid = c[1].data["V_grid"].astype(np.float64)
        combined = c[1].data["ccf"].astype(np.float64)
    with fits.open(fitspec) as f:
        barymjd = f[1].header["BARYMJD"]
        z_hdr = f[2].header["wtd_mdpt_bc"]
        orders_l1 = f[1].data["orders"].astype(int)
        wave_l1 = f[1].data["wavelength"].astype(np.float64)
    h = l4.headers["PRIMARY"]
    assert h["DATALVL"] == "L4"
    assert h["INSTRUME"] == "EXPRES"
    assert h["ISSOLAR"] is True
    assert h["RVMETHOD"] == "CCF"
    assert h["RV"] == pytest.approx(ch["V"] / 1e5)
    assert h["RVERR"] == pytest.approx(ch["E_V"] / 1e5)
    assert h["BJDTDB"] == pytest.approx(barymjd + 2400000.5, abs=1e-9)
    assert h["BERV"] == pytest.approx(z_hdr * const.c.to("km/s").value, rel=1e-9)
    assert h["SYSVEL"] == 0.0
    # DRPTAG at L4 is the ccf file's pipeline version, not the fitspec's
    assert h["DRPTAG"] == ch["VERSION"]

    rv1 = l4.data["RV1"]
    assert len(rv1) == len(per_order) == 75
    np.testing.assert_array_equal(rv1["ECHELLE_ORDER"].value, per_order["orders"].astype(int))
    np.testing.assert_allclose(rv1["RV"].value, per_order["v"].astype(np.float64) / 1e5)
    np.testing.assert_allclose(rv1["RV_ERR"].value, per_order["e_v"].astype(np.float64) / 1e5)
    assert np.all(rv1["BJD_TDB"].value == h["BJDTDB"])
    assert np.all(rv1["BERV"].value == h["BERV"])
    # ORDER_INDEX and wavelength extent come from the matching fitspec row
    for i, order in enumerate(rv1["ECHELLE_ORDER"].value[:5]):
        row = int(np.where(orders_l1 == order)[0][0])
        assert rv1["ORDER_INDEX"].value[i] == row
        assert rv1["WAVE_START"].value[i] == pytest.approx(wave_l1[row].min())
        assert rv1["WAVE_END"].value[i] == pytest.approx(wave_l1[row].max())
    rv1h = l4.headers["RV1"]
    assert rv1h["RVMETHOD"] == "CCF"
    assert rv1h["SKYRMVD"] is False
    assert rv1h["TELLRMVD"] == bool(ch["DIV_TELL"])
    assert rv1h["RED_ORD"] == ch["RED_ORD"]
    assert rv1h["BLUE_ORD"] == ch["BLUE_ORD"]

    ccf1 = l4.data["CCF1"]
    assert ccf1.shape == (75, 1001)
    np.testing.assert_array_equal(ccf1, per_order["ccfs"].astype(np.float64))
    c1h = l4.headers["CCF1"]
    assert c1h["VELSTART"] == pytest.approx(vgrid[0] / 1e5)
    assert c1h["VELSTEP"] == pytest.approx((vgrid[1] - vgrid[0]) / 1e5)
    assert c1h["VELNSTEP"] == 1001
    assert c1h["CCFMASK"] == ch["MASK"]
    assert c1h["RED_ORD"] == ch["RED_ORD"]
    assert c1h["BLUE_ORD"] == ch["BLUE_ORD"]

    cc = l4.data["CUSTOM_CCF1"]
    assert cc.shape == (1, 1001)
    np.testing.assert_array_equal(cc[0], combined)
    crv = l4.data["CUSTOM_RV1"]
    assert len(crv) == 1
    assert crv["RV"].value[0] == pytest.approx(ch["V"] / 1e5)
    assert crv["RV_ERR"].value[0] == pytest.approx(ch["E_V"] / 1e5)

    diag = l4.data["DIAGNOSTICS1"]
    assert list(diag.colnames) == ["metric_name", "value", "uncertainty"]
    names = list(diag["metric_name"].value)
    for n in ("HALPHA", "HWIDTH", "CCFFWHM", "BIS", "CCFVSPAN", "BIGAUSS",
              "SKEWNORM", "QUALITY", "SNR", "CHI2", "EXPCOUNT", "WATRCOL"):
        assert n in names, n
    fwhm = diag[diag["metric_name"] == "CCFFWHM"][0]
    assert np.isfinite(fwhm["value"]) and np.isfinite(fwhm["uncertainty"])
    ext = l4.data["EXT_DESCRIPT"]
    assert list(ext.colnames) == ["Name", "Description"]
    for n in ("PRIMARY", "INSTRUMENT_HEADER", "RECEIPT", "DRP_CONFIG",
              "EXT_DESCRIPT", "RV1", "CCF1", "CUSTOM_CCF1", "CUSTOM_RV1",
              "DIAGNOSTICS1"):
        assert n in ext["Name"].value, n
    drp = l4.data["DRP_CONFIG"]
    assert any(e.startswith("ccf:MASK = ESPRESSO_G2.fits") for e in drp["ENTRY"].value)


def test_expres_l4_explicit_l1_file(tmp_path):
    fitspec, ccf = expres_fixture_files()
    alone = tmp_path / "ccf_only" / FIXTURE_BASENAME
    alone.parent.mkdir()
    alone.write_bytes(open(ccf, "rb").read())
    l4 = RV4.from_fits(str(alone), instrument="EXPRES", l1_file=fitspec)
    assert l4.headers["PRIMARY"]["DATALVL"] == "L4"
    assert len(l4.data["RV1"]) == 75


def test_expres_l4_missing_fitspec(tmp_path):
    _, ccf = expres_fixture_files()
    alone = tmp_path / "ccf" / FIXTURE_BASENAME
    alone.parent.mkdir()
    alone.write_bytes(open(ccf, "rb").read())
    with pytest.raises(FileNotFoundError) as err:
        RV4.from_fits(str(alone), instrument="EXPRES")
    msg = str(err.value)
    assert str(tmp_path / "fitspec" / FIXTURE_BASENAME) in msg
    assert "l1_file" in msg


def test_expres_l4_diagnostics_optional_cards(tmp_path):
    """An older fitspec without activity cards still translates."""
    fitspec, ccf = expres_fixture_files()
    root = tmp_path
    (root / "fitspec").mkdir()
    (root / "ccf").mkdir()
    with fits.open(fitspec, memmap=False) as hdul:
        for key in ("S-VALUE", "HALPHA", "HWIDTH", "CCFFWHM", "CCFFWHME", "BIS",
                    "CCFVSPAN", "BIGAUSS", "SKEWNORM", "QUALITY", "WATRCOL"):
            if key in hdul[1].header:
                del hdul[1].header[key]
        hdul.writeto(root / "fitspec" / FIXTURE_BASENAME)
    (root / "ccf" / FIXTURE_BASENAME).write_bytes(open(ccf, "rb").read())
    l4 = RV4.from_fits(str(root / "ccf" / FIXTURE_BASENAME), instrument="EXPRES")
    names = list(l4.data["DIAGNOSTICS1"]["metric_name"].value)
    assert "HALPHA" not in names
    assert "SNR" in names and "CHI2" in names


def test_expres_l4_rv1_weight():
    _, ccf = expres_fixture_files()
    l4 = RV4.from_fits(ccf, instrument="EXPRES")
    with fits.open(ccf) as c:
        ch = c[0].header
        per_order = c[2].data
        orders = per_order["orders"].astype(int)
        e_v = per_order["e_v"].astype(np.float64)
    rv1 = l4.data["RV1"]
    assert "WEIGHT" in rv1.colnames
    weight = rv1["WEIGHT"].value
    red_ord, blue_ord = ch["RED_ORD"], ch["BLUE_ORD"]
    in_window = (orders >= red_ord) & (orders <= blue_ord)
    finite = np.isfinite(e_v)
    # fixture has four non-finite RV_ERR rows, at these echelle orders
    non_finite_orders = set(orders[~finite].tolist())
    assert non_finite_orders == {109, 106, 97, 81}
    assert np.all(weight[~finite] == 0.0)
    assert np.all(weight[~in_window] == 0.0)
    expected_ones = in_window & finite
    assert np.all(weight[expected_ones] == 1.0)
    assert np.all(weight[~expected_ones] == 0.0)
    assert (weight == 1.0).sum() == expected_ones.sum()


def test_expres_l4_order_mismatch(tmp_path):
    """A ccf echelle order missing from the fitspec ORDER_TABLE raises a
    clear error instead of a bare numpy IndexError."""
    fitspec, ccf = expres_fixture_files()
    root = tmp_path
    (root / "fitspec").mkdir()
    (root / "ccf").mkdir()
    with fits.open(fitspec, memmap=False) as hdul:
        orders = hdul[1].data["orders"].astype(int)
        assert orders[0] == 160
        keep = np.arange(len(orders)) != 0
        hdul[1].data = hdul[1].data[keep]
        hdul.writeto(root / "fitspec" / FIXTURE_BASENAME)
    (root / "ccf" / FIXTURE_BASENAME).write_bytes(open(ccf, "rb").read())
    with pytest.raises(ValueError, match="order 160"):
        RV4.from_fits(str(root / "ccf" / FIXTURE_BASENAME), instrument="EXPRES")


def test_expres_l4_stellar_branch(tmp_path):
    """A non-solar target leaves SYSVEL undefined at L4."""
    fitspec, ccf = expres_fixture_files()
    root = tmp_path
    (root / "fitspec").mkdir()
    (root / "ccf").mkdir()
    with fits.open(fitspec, memmap=False) as hdul:
        hdul[0].header["OBJECT"] = "HD 10700"
        hdul[0].header["OBSTYPE"] = "Science"
        hdul.writeto(root / "fitspec" / FIXTURE_BASENAME)
    (root / "ccf" / FIXTURE_BASENAME).write_bytes(open(ccf, "rb").read())
    l4 = RV4.from_fits(str(root / "ccf" / FIXTURE_BASENAME), instrument="EXPRES")
    h = l4.headers["PRIMARY"]
    assert h["ISSOLAR"] is False
    assert h["SYSVEL"] is None


def test_expres_l2_snr_consistency():
    """The pipeline's per-order SNR card is an independent check on the
    flux/variance scaling: max(flux / sqrt(var)) in that order should be
    close to the pipeline's own per-pixel SNR number."""
    fitspec, ccf = expres_fixture_files()
    l2 = RV2.from_fits(fitspec, instrument="EXPRES")
    with fits.open(ccf) as c:
        snr_card = c[0].header["SNR"]
        snr_comment = c[0].header.comments["SNR"]
    order = int(re.search(r"order (\d+)", snr_comment).group(1))
    ot = l2.data["ORDER_TABLE"]
    row = int(np.where(ot["ECHELLE_ORDER"].value == order)[0][0])
    flux = l2.data["TRACE1_FLUX"][row]
    var = l2.data["TRACE1_VAR"][row]
    quality = l2.data["TRACE1_QUALITY"][row]
    good = (quality == 0) & np.isfinite(flux) & np.isfinite(var) & (var > 0)
    snr = np.max(flux[good] / np.sqrt(var[good]))
    assert snr == pytest.approx(snr_card, rel=0.15)
