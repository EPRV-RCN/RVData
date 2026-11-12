"""Reader and import regressions: names the readers must accept, and
modules that must import without optional system pieces."""

import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from rvdata.core.models.level2 import RV2
from rvdata.core.models.level4 import RV4


def _image(name, shape=(3, 5), value=1.0):
    return fits.ImageHDU(np.full(shape, value), name=name)


@pytest.mark.parametrize(
    "extname",
    ["CUSTOM1_TRACE1_FLUX", "CUSTOM1_TRACE1_WAVE", "CUSTOM2_TRACE3_VAR"],
)
def test_rv2_reads_custom_trace_extensions(extname):
    # "CUSTOM1_TRACE1_FLUX".split("_")[1] is "TRACE1", so the old reader
    # looked up "TRACE1_TRACE1" and raised KeyError.
    hdul = fits.HDUList(
        [fits.PrimaryHDU(), _image("TRACE1_FLUX"), _image(extname, value=7.0)]
    )
    rv2 = RV2()
    rv2._read(hdul)
    assert extname in rv2.extensions
    np.testing.assert_array_equal(rv2.data[extname], np.full((3, 5), 7.0))


def test_rv2_reads_extra_traces():
    hdul = fits.HDUList(
        [fits.PrimaryHDU(), _image("TRACE2_FLUX", value=2.0), _image("TRACE2_WAVE")]
    )
    rv2 = RV2()
    rv2._read(hdul)
    np.testing.assert_array_equal(rv2.data["TRACE2_FLUX"], np.full((3, 5), 2.0))


def test_rv2_reads_plain_drift_extension():
    # The HARPS and HARPS-N translators write a "DRIFT" image, which is not a
    # row in L2-extensions.csv (only TRACE1_DRIFT is), so reading back raised
    # KeyError.
    hdul = fits.HDUList([fits.PrimaryHDU(), _image("DRIFT", value=3.0)])
    rv2 = RV2()
    rv2._read(hdul)
    np.testing.assert_array_equal(rv2.data["DRIFT"], np.full((3, 5), 3.0))


def test_rv2_unknown_extension_warns_instead_of_crashing():
    hdul = fits.HDUList([fits.PrimaryHDU(), _image("SOMETHING_ELSE")])
    rv2 = RV2()
    with pytest.warns(UserWarning, match="SOMETHING_ELSE"):
        rv2._read(hdul)
    assert "SOMETHING_ELSE" in rv2.extensions


def test_rv4_reads_second_diagnostics_table():
    # Only DIAGNOSTICS1 is in L4-extensions.csv, so DIAGNOSTICS2 raised KeyError.
    diag = Table({"metric_name": ["SHK"], "value": [0.17]})
    hdul = fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.BinTableHDU(diag, name="DIAGNOSTICS1"),
            fits.BinTableHDU(diag, name="DIAGNOSTICS2"),
        ]
    )
    rv4 = RV4()
    rv4._read(hdul)
    assert "DIAGNOSTICS2" in rv4.extensions
    assert list(rv4.data["DIAGNOSTICS2"]["value"]) == pytest.approx([0.17])


def _python(code, **env_overrides):
    env = dict(os.environ)
    env.update(env_overrides)
    return subprocess.run(
        [sys.executable, "-c", code], env=env, capture_output=True, text=True
    )


def test_import_without_git_executable():
    # GitPython raises ImportError at import time when it cannot find git,
    # which made every rvdata import fail on minimal compute nodes.
    result = _python(
        "import rvdata.core.models.level2",
        GIT_PYTHON_GIT_EXECUTABLE="/nonexistent/git",
        GIT_PYTHON_REFRESH="raise",
    )
    assert result.returncode == 0, result.stderr


def test_receipt_without_git_executable():
    result = _python(
        "from rvdata.core.models.level2 import RV2\n"
        "rv2 = RV2()\n"
        "rv2.receipt_add_entry('test', '', 'PASS')\n"
        "print(len(rv2.data['RECEIPT']))",
        GIT_PYTHON_GIT_EXECUTABLE="/nonexistent/git",
        GIT_PYTHON_REFRESH="raise",
    )
    assert result.returncode == 0, result.stderr


def test_level3_import_without_bindensity():
    # Reading and writing L3 files needs no resampling, so a missing
    # bindensity build (no wheel for the interpreter) must not block it.
    result = _python(
        "import sys\n"
        "sys.modules['bindensity'] = None\n"
        "import rvdata.core.models.level3\n"
        "import rvdata.core.tools.stitch_spectrum"
    )
    assert result.returncode == 0, result.stderr
