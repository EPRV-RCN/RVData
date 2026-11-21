"""ESPRESSO S1D to L3 conversion, on small synthetic S1D files."""

import importlib.resources

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from rvdata.core.models.level3 import RV3
from rvdata.instruments.espresso.utils import convert_S1D

NPIX = 6
FILES = {
    "s1d_A": ("CORR", 1.0),
    "s1d_B": (None, 2.0),  # fiber B is skipped
    "s1d_skysub_A": ("SKYSUB", 3.0),
    "s1d_tell_corr_A": ("TELLCORR", 4.0),
}


def _write_s1d(path, scale):
    wave = np.linspace(5000.0, 5001.0, NPIX)
    flux = scale * np.arange(100.0, 100.0 + NPIX)
    error = np.sqrt(flux)  # the DRS writes 1-sigma errors, in e-
    table = Table({"wavelength": wave, "flux": flux, "error": error})
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(table)]).writeto(path)
    return flux, error


@pytest.fixture
def converted(tmp_path):
    names, errors = {}, {}
    for key, (corr, scale) in FILES.items():
        path = tmp_path / f"{key}.fits"
        _, errors[key] = _write_s1d(path, scale)
        names[key] = str(path)
    rv3 = RV3()
    # ESPRESSORV3 calls create_PRIMARY first, which replaces the whole PRIMARY
    # header (and with it RV3's level 3 defaults) by the level 2 one.
    rv3.set_header("PRIMARY", fits.Header({"DATALVL": "L3"}))
    convert_S1D(rv3, names)
    return rv3, errors


def _l3_required_keywords():
    config = importlib.resources.files("rvdata.core.models.config")
    table = pd.read_csv(config / "L3-PRIMARY-keywords.csv")
    return table.loc[table["Required"].astype(bool), "Keyword"].tolist()


def test_writes_every_required_l3_keyword(converted):
    rv3, _ = converted
    missing = [k for k in _l3_required_keywords() if k not in rv3.headers["PRIMARY"]]
    assert missing == []


def test_flux_is_not_normalized(converted):
    rv3, _ = converted
    assert rv3.headers["PRIMARY"]["FLXNRMMD"] == "None"


@pytest.mark.parametrize("key", ["s1d_A", "s1d_skysub_A", "s1d_tell_corr_A"])
def test_variance_is_error_squared(converted, key):
    rv3, errors = converted
    extname = f"STITCHED_{FILES[key][0]}_SCI_VAR"
    np.testing.assert_allclose(rv3.data[extname], errors[key] ** 2)


def test_flux_and_wave_copied(converted):
    rv3, _ = converted
    np.testing.assert_allclose(
        rv3.data["STITCHED_CORR_SCI_FLUX"], np.arange(100.0, 100.0 + NPIX)
    )
    np.testing.assert_allclose(
        rv3.data["STITCHED_CORR_SCI_WAVE"], np.linspace(5000.0, 5001.0, NPIX)
    )
