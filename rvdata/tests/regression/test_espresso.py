import os
import re
from pathlib import Path

import requests

from rvdata.core.models.base import RVDataModel
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level3 import RV3
from rvdata.core.models.level4 import RV4
from rvdata.tests.regression.compliance import (
    check_l2_compliance,
    check_l3_compliance,
    check_l4_compliance,
)

FILE_URLS = {
    "ESPRESSO": {
        "raw": "http://grinnell.as.arizona.edu/~rvdata/espresso/ESPRE.2017-12-03T02-09-40.348.fits",
        "S2D_BLAZE_A": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T02-09-40.348_S2D_BLAZE_A.fits",
        "S2D_BLAZE_B": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T02-09-40.348_S2D_BLAZE_B.fits",
        "BLAZE_A": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T10-43-59.835_BLAZE_A.fits",
        "BLAZE_B": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T10-43-59.835_BLAZE_B.fits",
        "S1D_A": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T02-09-40.348_S1D_A.fits",
        "S1D_B": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T02-09-40.348_S1D_B.fits",
        "S1D_TELL_CORR_A": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T02-09-40.348_S1D_TELL_CORR_A.fits",
        "DRIFT_MATRIX_B": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T02-09-40.348_DRIFT_MATRIX_B.fits",
        "CCF_A": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T02-09-40.348_CCF_A.fits",
        "CCF_TELL_CORR_A": "http://grinnell.as.arizona.edu/~rvdata/espresso/r.ESPRE.2017-12-03T02-09-40.348_CCF_TELL_CORR_A.fits",
    }
}


def download_instrument_files(instrument: str = "ESPRESSO") -> dict[str, Path]:
    """Download all files for the specified instrument."""

    # CA bundle path relative to this test file
    # CA_BUNDLE_PATH = Path(__file__).parent / "fixtures" / "dace-unige-ch-chain.pem"

    local_files = {}

    for key, url in FILE_URLS[instrument].items():
        # get the file name from the URL
        filename = url.rsplit("/", 1)[-1]
        if os.name == "nt":
            # For Windows: colons are not allowed in filenames, use underscores
            # to match what get_files_names.py produces via replace(":", "_")
            filename = re.sub(r"T(\d{2})-(\d{2})-(\d{2})", r"T\1_\2_\3", filename)
        else:
            # Restore colons in ISO timestamp time portions (HH-MM-SS -> HH:MM:SS)
            filename = re.sub(r"T(\d{2})-(\d{2})-(\d{2})", r"T\1:\2:\3", filename)
        filepath = Path(filename)

        if not filepath.exists():
            # response = requests.get(url, verify=str(CA_BUNDLE_PATH), timeout=30)
            response = requests.get(url)

            response.raise_for_status()
            filepath.write_bytes(response.content)

        local_files[key] = filepath

    return local_files


def test_espresso():
    """Test ESPRESSO data processing for L2 and L4 standards."""
    files = download_instrument_files("ESPRESSO")
    raw_file = files["raw"]

    # Test Level 2 - use auto-generated filename
    espr2 = RV2.from_fits(str(raw_file), instrument="ESPRESSO")
    l2_standard = espr2.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l2_standard)), \
        f"L2 filename '{os.path.basename(l2_standard)}' does not match EPRV convention"
    assert os.path.basename(l2_standard).startswith("espresso_SL2_"), \
        f"L2 filename should start with 'espresso_SL2_', got '{os.path.basename(l2_standard)}'"
    check_l2_compliance(l2_standard)

    # Test Level 3 - use auto-generated filename
    espr3 = RV3.from_fits(str(raw_file), instrument="ESPRESSO")
    l3_standard = espr3.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l3_standard)), \
        f"L3 filename '{os.path.basename(l3_standard)}' does not match EPRV convention"
    assert os.path.basename(l3_standard).startswith("espresso_SL3_"), \
        f"L3 filename should start with 'espresso_SL3_', got '{os.path.basename(l3_standard)}'"
    check_l3_compliance(l3_standard)

    # Test Level 4 - use auto-generated filename
    espr4 = RV4.from_fits(str(raw_file), instrument="ESPRESSO")
    l4_standard = espr4.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l4_standard)), \
        f"L4 filename '{os.path.basename(l4_standard)}' does not match EPRV convention"
    assert os.path.basename(l4_standard).startswith("espresso_SL4_"), \
        f"L4 filename should start with 'espresso_SL4_', got '{os.path.basename(l4_standard)}'"
    check_l4_compliance(l4_standard)


def test_espresso_benchmark(benchmark):
    """Benchmark wrapper for test_espresso."""
    download_instrument_files("ESPRESSO")
    benchmark(test_espresso)


if __name__ == "__main__":
    test_espresso()
