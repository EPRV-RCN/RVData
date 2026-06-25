# Tests for RVDataModel header serialization (core.models.base)
from collections import OrderedDict

from astropy.io import fits

from rvdata.core.models.level2 import RV2


def _write_minimal_l2(path):
    """Write a minimal standard-format L2 file whose PRIMARY header carries
    comments on both a standard keyword (recast on read) and a non-standard one."""
    primary = fits.PrimaryHDU()
    primary.header["INSTRUME"] = ("KPF", "Instrument name")  # standard keyword
    primary.header["EXPTIME"] = (60.0, "[s] Exposure time")  # standard keyword
    primary.header["MYNOTE"] = ("hi", "a non-standard note")  # not recast
    fits.HDUList([primary]).writeto(path, overwrite=True)


def test_roundtrip_preserves_primary_comments(tmp_path):
    """from_fits -> to_fits must not drop PRIMARY header comments.

    Regression: from_fits stores PRIMARY as an astropy fits.Header (comments in
    a parallel store). Two sites silently dropped those comments on round-trip:
    the read-path keyword recast overwrote them with "", and _create_hdul
    rebuilt PRIMARY from bare scalars. Comments like the [s]/[km/s] unit
    annotations vanished from the written file.
    """
    src = tmp_path / "min_L2.fits"
    _write_minimal_l2(src)

    model = RV2.from_fits(str(src))
    model.to_fits(out_filedir=str(tmp_path), out_filename="roundtrip.fits")

    written = fits.open(tmp_path / "roundtrip.fits")[0].header
    assert written.comments["INSTRUME"] == "Instrument name"
    assert written.comments["EXPTIME"] == "[s] Exposure time"
    assert written.comments["MYNOTE"] == "a non-standard note"


def test_create_hdul_preserves_primary_comments_from_ordereddict():
    """The translator-built form (OrderedDict of (value, comment) tuples) must
    keep working and retain comments through _create_hdul."""
    model = RV2()
    primary = OrderedDict()
    primary["INSTRUME"] = ("KPF", "Instrument name")
    primary["EXPTIME"] = (60.0, "[s] Exposure time")
    model.headers["PRIMARY"] = primary

    written = next(h for h in model._create_hdul() if h.name == "PRIMARY").header
    assert written.comments["INSTRUME"] == "Instrument name"
    assert written.comments["EXPTIME"] == "[s] Exposure time"
