"""Unit tests for the compliance checker itself (rvdata/tests/regression)."""

import importlib

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from rvdata.core.models.level2 import RV2
from rvdata.tests.regression import compliance


def _required_keys(*csv_names):
    config = importlib.resources.files("rvdata.core.models.config")
    frames = [pd.read_csv(config / name) for name in csv_names]
    table = pd.concat(frames, ignore_index=True)
    keys = table.loc[table["Required"].map(compliance._is_required), "Keyword"]
    return [key.split("...")[0].strip() for key in keys]


@pytest.mark.parametrize(
    "value, expected",
    [
        (True, True),
        (False, False),
        (np.bool_(True), True),
        ("True", True),
        ("False", False),
        ("Yes", True),
        ("No", False),
        (" y ", True),
        (float("nan"), False),
    ],
)
def test_is_required(value, expected):
    assert compliance._is_required(value) is expected


def test_required_column_parses_as_bool():
    # The reason `row["Required"] == "Y"` never matched: pandas reads the
    # True/False column as bool.
    config = importlib.resources.files("rvdata.core.models.config")
    for name in ("L3-PRIMARY-keywords.csv", "L4-PRIMARY-keywords.csv"):
        assert pd.read_csv(config / name)["Required"].dtype == bool


@pytest.mark.parametrize(
    "check", [compliance.check_l3_header, compliance.check_l4_header]
)
def test_header_check_rejects_empty_header(check):
    with pytest.raises(AssertionError, match="not found in header"):
        check(fits.Header())


@pytest.mark.parametrize(
    "check, level_csv",
    [
        (compliance.check_l3_header, "L3-PRIMARY-keywords.csv"),
        (compliance.check_l4_header, "L4-PRIMARY-keywords.csv"),
    ],
)
def test_header_check_accepts_all_required_keywords(check, level_csv):
    header = fits.Header()
    for key in _required_keys("L2-PRIMARY-keywords.csv", level_csv):
        header[key] = "x"
    check(header)


@pytest.mark.parametrize(
    "check, level_csv, missing",
    [
        (compliance.check_l3_header, "L3-PRIMARY-keywords.csv", "BLZCORR"),
        (compliance.check_l4_header, "L4-PRIMARY-keywords.csv", "RVMETHOD"),
    ],
)
def test_header_check_names_missing_level_keyword(check, level_csv, missing):
    header = fits.Header()
    for key in _required_keys("L2-PRIMARY-keywords.csv", level_csv):
        if key != missing:
            header[key] = "x"
    with pytest.raises(AssertionError, match=missing):
        check(header)


def test_l2_extension_check_opens_file_once(tmp_path, monkeypatch):
    # check_l2_extensions used to call fits.open a second time inside its
    # own `with fits.open(...)` block and never close that handle.
    path = tmp_path / "rv2.fits"
    RV2().to_fits(str(path))

    opened = []
    real_open = fits.open

    def counting_open(*args, **kwargs):
        hdul = real_open(*args, **kwargs)
        opened.append(hdul)
        return hdul

    monkeypatch.setattr(compliance.fits, "open", counting_open)
    try:
        compliance.check_l2_extensions(str(path))
    except AssertionError:
        pass  # this test is about file handles, not the content of a bare RV2
    assert len(opened) == 1
