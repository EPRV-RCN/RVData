import requests
import os
from astropy.io import fits
from astropy.time import Time
from rvdata.core.models.base import RVDataModel
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level3 import RV3
from rvdata.core.models.level4 import RV4

from rvdata.tests.regression.compliance import (
    check_l2_compliance, check_l3_compliance, check_l4_compliance,
)


file_urls = {
    "KPF": [
        "http://grinnell.as.arizona.edu/~rvdata/kpf/KP.20250208.17485.59.fits",
        "http://grinnell.as.arizona.edu/~rvdata/kpf/KP.20250208.17485.59_L1.fits",
        "http://grinnell.as.arizona.edu/~rvdata/kpf/KP.20241022.41656.30_L2.fits",
    ]
}


def download_file(url, filename):
    response = requests.get(url)
    response.raise_for_status()  # Check if the request was successful
    with open(filename, "wb") as file:
        file.write(response.content)


def download_files():
    l0file = "kpf_L0.fits"
    l1file = "kpf_L1.fits"
    l2file = "kpf_L2.fits"
    if not os.path.exists(l0file):
        download_file(file_urls["KPF"][0], l0file)
    if not os.path.exists(l1file):
        download_file(file_urls["KPF"][1], l1file)
    if not os.path.exists(l2file):
        download_file(file_urls["KPF"][2], l2file)

    return l0file, l1file, l2file


def assert_full_dateobs(standard_file):
    # Native KPF headers keep only the date in DATE-OBS (the time lives in
    # DATE-BEG); the standard requires a full ISO datetime, and a date-only
    # value silently becomes T000000 in generated filenames.
    dateobs = str(fits.getheader(standard_file)["DATE-OBS"])
    assert "T" in dateobs, \
        f"DATE-OBS '{dateobs}' in {standard_file} lacks a time component"


def assert_jd_utc_matches_dateobs(standard_file):
    # The standard defines JD_UTC as the full Julian Date of DATE-OBS. Native
    # KPF only carries MJD-OBS, which was once copied across unconverted,
    # leaving ~60605 where ~2460605 belongs (issue #186).
    header = fits.getheader(standard_file)
    jd_utc = header["JD_UTC"]
    expected = Time(header["DATE-OBS"], scale="utc").jd
    assert abs(jd_utc - expected) < 1e-5, \
        f"JD_UTC {jd_utc} in {standard_file} is not the JD of DATE-OBS " \
        f"{header['DATE-OBS']} ({expected})"


def test_kpf():
    l0file, l1file, l2file = download_files()

    # Check L2 - use auto-generated filename
    kpf2 = RV2.from_fits(l1file, l0file=l0file, instrument="KPF")
    l2_standard = kpf2.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l2_standard)), \
        f"L2 filename '{os.path.basename(l2_standard)}' does not match EPRV convention"
    assert os.path.basename(l2_standard).startswith("kpf_SL2_"), \
        f"L2 filename should start with 'kpf_SL2_', got '{os.path.basename(l2_standard)}'"
    assert_full_dateobs(l2_standard)
    assert_jd_utc_matches_dateobs(l2_standard)
    check_l2_compliance(l2_standard)

    # Check L3
    # Note: L3 creation requires an RVData-standard L2 file (created above).
    # Native KPF L2 files must first be converted using KPFRV2.from_fits().
    kpf3 = RV3.from_fits(l2_standard, instrument="KPF")
    l3_standard = kpf3.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l3_standard)), \
        f"L3 filename '{os.path.basename(l3_standard)}' does not match EPRV convention"
    assert os.path.basename(l3_standard).startswith("kpf_SL3_"), \
        f"L3 filename should start with 'kpf_SL3_', got '{os.path.basename(l3_standard)}'"
    check_l3_compliance(l3_standard)

    # Check L4 - use auto-generated filename
    kpf4 = RV4.from_fits(l2file, instrument="KPF")
    l4_standard = kpf4.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l4_standard)), \
        f"L4 filename '{os.path.basename(l4_standard)}' does not match EPRV convention"
    assert os.path.basename(l4_standard).startswith("kpf_SL4_"), \
        f"L4 filename should start with 'kpf_SL4_', got '{os.path.basename(l4_standard)}'"
    assert_full_dateobs(l4_standard)
    assert_jd_utc_matches_dateobs(l4_standard)
    check_l4_compliance(l4_standard)


def test_kpf_benchmark(benchmark):
    # run test_kpf() once to download the files
    _ = download_files()
    # now run it again with benchmark
    benchmark(test_kpf)


if __name__ == "__main__":
    test_kpf()
