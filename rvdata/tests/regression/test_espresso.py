import requests
from pathlib import Path
from typing import Dict
from urllib.parse import urlparse

from rvdata.core.models.level2 import RV2
from rvdata.core.models.level4 import RV4
from rvdata.tests.regression.compliance import (
    check_l2_extensions,
    check_l2_header,
    check_l4_extensions,
    check_l4_header,
)

FILE_URLS = {
    "ESPRESSO": {
        "raw": "https://dace.unige.ch/downloads/EPRV_standard_data/ESPRE.2017-12-03T02:09:40.348.fits",
        "S2D_BLAZE_A": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S2D_BLAZE_A.fits",
        "S2D_BLAZE_B": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S2D_BLAZE_B.fits",
        "BLAZE_A": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T10:43:59.835_BLAZE_A.fits",
        "BLAZE_B": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T10:43:59.835_BLAZE_B.fits",
        "S1D_A": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S1D_A.fits",
        "S1D_B": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S1D_B.fits",
        "S1D_TELL_CORR_A": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_S1D_TELL_CORR_A.fits",
        "DRIFT_MATRIX_B": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_DRIFT_MATRIX_B.fits",
        "CCF_A": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_CCF_A.fits",
        "CCF_TELL_CORR_A": "https://dace.unige.ch/downloads/EPRV_standard_data/r.ESPRE.2017-12-03T02:09:40.348_CCF_TELL_CORR_A.fits",
    }
}


def download_instrument_files(instrument: str = "ESPRESSO") -> Dict[str, Path]:
    """Download all files for the specified instrument."""
    local_files = {}

    for key, url in FILE_URLS[instrument].items():
        # Extract filename from URL and sanitize for Windows
        filename = Path(urlparse(url).path).name
        filename = filename.replace(":", "_")  # Windows-safe
        filepath = Path(filename)

        if not filepath.exists():
            response = requests.get(url, verify=False)
            response.raise_for_status()
            filepath.write_bytes(response.content)

        local_files[key] = filepath

    return local_files


def test_espresso():
    """Test ESPRESSO data processing for L2 and L4 standards."""
    files = download_instrument_files("ESPRESSO")
    raw_file = files["raw"]

    # Test Level 2
    espr2 = RV2.from_fits(str(raw_file), instrument="ESPRESSO")
    l2_standard = Path("esp_L2_standard.fits")
    espr2.to_fits(str(l2_standard))
    l2_obj = RV2.from_fits(str(l2_standard))
    check_l2_extensions(str(l2_standard))
    check_l2_header(l2_obj.headers["PRIMARY"])

    # Test Level 4
    espr4 = RV4.from_fits(str(raw_file), instrument="ESPRESSO")
    l4_standard = Path("espr_L4_standard.fits")
    espr4.to_fits(str(l4_standard))
    l4_obj = RV4.from_fits(str(l4_standard))
    check_l4_extensions(str(l4_standard))
    check_l4_header(l4_obj.headers["PRIMARY"])


def test_espresso_benchmark(benchmark):
    """Benchmark wrapper for test_espresso."""
    download_instrument_files("ESPRESSO")
    benchmark(test_espresso)


if __name__ == "__main__":
    test_espresso()
