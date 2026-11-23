"""NEID translator regressions, on small synthetic NEID L2 files.

The headers carry the keys the translator reads, with the values that broke it:
solar frames (blank QOBJECT), calibration frames, and a missing DRIFTRV0.
"""

import csv
import importlib.resources

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from rvdata.instruments.neid.level2 import NEIDRV2
from rvdata.instruments.neid.level4 import NEIDRV4
from rvdata.instruments.neid.utils.make_neid_primary_header import (
    make_base_primary_header,
)
from rvdata.tests.regression import compliance

NORDER, NPIX = 122, 8
ORDERS = range(173, 51, -1)  # NEID echelle orders, index 0 first
ACTIVITY_INDICES = [
    "CaIIHK", "HeI_1", "NaI", "Ha06_1", "Ha16_1", "CaI_1",
    "CaIRT1", "CaIRT2", "CaIRT3", "NaINIR", "PaDelta", "Mn539",
]  # fmt: skip


def _header_map_instrument_keys():
    config = importlib.resources.files("rvdata.instruments.neid") / "config"
    with (config / "header_map.csv").open() as f:
        return {row["INSTRUMENT"] for row in csv.DictReader(f) if row["INSTRUMENT"]}


def _primary(**overrides):
    """A NEID primary header for a stellar HR science frame, with overrides."""
    h = fits.Header()
    for key in sorted(_header_map_instrument_keys()):
        h[key] = "x"
    h.update(
        {
            "OBSJD": 2460000.5,
            "EXPTIME": 300.0,
            "DQLEVEL0": 0,
            "DQLEVEL1": 0,
            "DQLEVEL2": 0,
            "OBS-MODE": "HR",
            "OBSTYPE": "Sci",
            "OBJECT": "HD 4628",
            "SCI-OBJ": "HD 4628",
            "SKY-OBJ": "Sky",
            "CAL-OBJ": "Etalon",
            "QOBJECT": "HD 4628",
            "QRA": "04:52:37.12",
            "QDEC": "+05:56:46.5",
            "QEQNX": 2000.0,
            "QEPOCH": 2000.0,
            "QPLX": 134.5,
            "QPMRA": 630.0,
            "QPMDEC": -1104.0,
            "QRV": -10.2,
            "QZ": 0.0,
            "DRIFTRV0": 1.5,
        }
    )
    for order in ORDERS:
        h[f"SSBRV{order:03d}"] = 10.0
        h[f"SSBZ{order:03d}"] = 3.3e-5
        h[f"SSBJD{order:03d}"] = 2460000.6
    h.update(overrides)
    return h


# A solar-feed frame as the archive writes it (neidL2_20240617T172049): the target
# fields are blank or zero because the Sun is not a catalogue target.
SOLAR = dict(
    OBSTYPE="Eng", OBJECT="Sun", **{"SCI-OBJ": "Sun", "SKY-OBJ": "Dark"},
    QOBJECT="", QRA="0", QDEC="0", QEQNX="", QEPOCH="2000.0", QRV="0", QZ="",
)  # fmt: skip


def _ccfs_header():
    h = fits.Header(
        {
            "CCFJDMOD": 2460000.6,
            "CCFRVMOD": -10.1,
            "DVRMSMOD": 0.0005,
            "CCFSTART": -100.0,
            "CCFSTEP": 0.5,
            "CCFMASK": "G2",
            "FWHMMOD": 6.5,
            "BISMOD": 0.01,
            "EBISMOD": 0.002,
        }
    )
    for order in ORDERS:
        h[f"CCFRV{order:03d}"] = -10.1
        h[f"CCFWT{order:03d}"] = 1.0
    return h


def _activity(drop=()):
    names = [n for i in ACTIVITY_INDICES for n in (i, i + "_tellcorr")]
    names = [n for n in names if n not in drop]
    values = np.arange(len(names), dtype=float)
    return Table({"INDEX": names, "VALUE": values, "UNCERTAINTY": values / 10})


def _neid_l2(primary=None, activity=None):
    """A minimal NEID L2 HDUList: 122 orders x 8 pixels, every extension NEIDRV2/4 read."""
    image = np.ones((NORDER, NPIX))
    wave = 4000.0 + np.arange(NORDER * NPIX, dtype=float).reshape(NORDER, NPIX)
    hdus = [fits.PrimaryHDU(header=primary if primary is not None else _primary())]
    for fiber in ("SCI", "SKY", "CAL"):
        hdus += [
            fits.ImageHDU(image, name=f"{fiber}FLUX"),
            fits.ImageHDU(wave.copy(), name=f"{fiber}WAVE"),
            fits.ImageHDU(image, name=f"{fiber}VAR"),
            fits.ImageHDU(image, name=f"{fiber}BLAZE"),
        ]
    hdus += [
        fits.ImageHDU(np.ones((5, 3, 4)), name="EXPMETER"),
        fits.ImageHDU(np.ones((NORDER, NPIX, 2)), name="TELLURIC"),
        # All-zero CCFs mark every order as having no RV, so NEIDRV4 never
        # indexes SCIWAVE at the real pixel ranges in neid_fsr.csv.
        fits.ImageHDU(np.zeros((NORDER, 5)), header=_ccfs_header(), name="CCFS"),
        fits.BinTableHDU(activity if activity is not None else _activity(), name="ACTIVITY"),
    ]
    return fits.HDUList(hdus)


def test_science_frame_names_no_calibration_source():
    # header_map.csv gives CLSRC1/CLSRC2 the default "None", which pandas read as
    # NaN, so the keyword was written with no value.
    phead = make_base_primary_header(_primary())
    assert phead["CLSRC1"] == "None"
    assert phead["CLSRC2"] == "None"


def test_calibration_frame_names_each_fiber_source():
    phead = make_base_primary_header(
        _primary(OBSTYPE="Cal", OBJECT="LFC", **{"SCI-OBJ": "LFC", "SKY-OBJ": "LFC"})
    )
    assert (phead["CLSRC1"], phead["CLSRC2"], phead["CLSRC3"]) == ("LFC", "LFC", "Etalon")


def test_stellar_frame_copies_catalogue_target():
    phead = make_base_primary_header(_primary())
    assert phead["CID1"] == "HD 4628"
    assert phead["CSRC1"] == "GAIADR2"
    assert phead["CRV1"] == pytest.approx(-10.2)


@pytest.mark.parametrize("key", ["CSRC1", "CID1", "CRA1", "CDEC1", "CEQNX1", "CEPCH1", "CRV1"])
def test_solar_frame_keeps_required_catalogue_keywords(key):
    # SCI-OBJ "Sun" never matched the blank QOBJECT, so these were left out.
    phead = make_base_primary_header(_primary(**SOLAR))
    assert key in phead


def test_solar_frame_marks_catalogue_values_unknown():
    phead = make_base_primary_header(_primary(**SOLAR))
    assert phead["CID1"] == "UNKNOWN"
    assert phead["CRV1"] is None  # Float keyword: undefined value, not a string


def test_solar_primary_header_passes_l2_compliance():
    compliance.check_l2_header(make_base_primary_header(_primary(**SOLAR)))


def test_l2_reads_frame_without_drift():
    # float(None) raised TypeError on frames whose DRIFTRV0 is blank.
    rv2 = NEIDRV2()
    rv2._read(_neid_l2(_primary(DRIFTRV0=None)))
    assert np.isnan(rv2.data["TRACE1_DRIFT"]).all()


def test_l2_reads_drift_in_kms():
    rv2 = NEIDRV2()
    rv2._read(_neid_l2(_primary(DRIFTRV0=1.5)))
    np.testing.assert_allclose(rv2.data["TRACE1_DRIFT"], [0.0015])


def test_l2_blaze_extensions_record_blznorm(tmp_path):
    rv2 = NEIDRV2()
    rv2._read(_neid_l2())
    path = tmp_path / "neid_SL2_20230101T000000.fits"
    rv2.to_fits(str(path))
    with fits.open(path) as hdul:
        for trace in (1, 2, 3):
            assert hdul[f"TRACE{trace}_BLAZE"].header["BLZNORM"] is False


def test_l4_reads_frame_missing_an_activity_index():
    # Indexing an empty match raised IndexError when the pipeline skipped an index.
    rv4 = NEIDRV4()
    rv4._read(_neid_l2(activity=_activity(drop=("PaDelta", "PaDelta_tellcorr"))))
    table = rv4.data["DIAGNOSTICS1"]
    diag = dict(zip(table["metric_name"], zip(table["value"], table["uncertainty"])))
    # The table masks NaN (as it already did for a missing EFWHMMOD).
    assert np.ma.is_masked(diag["PaDelta"][0]) or np.isnan(diag["PaDelta"][0])
    assert np.ma.is_masked(diag["PaDelta_tellcorr"][1]) or np.isnan(diag["PaDelta_tellcorr"][1])
    assert diag["CaIIHK"] == (0.0, 0.0)
