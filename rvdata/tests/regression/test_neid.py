import requests
import os
from rvdata.core.models.base import RVDataModel
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level3 import RV3
from rvdata.core.models.level4 import RV4

from rvdata.tests.regression.compliance import (
    check_l2_extensions, check_l2_header,
    check_l3_extensions, check_l3_header,
    check_l4_extensions, check_l4_header,
    check_l4_rv_columns,
    check_order_table_columns, check_receipt_columns,
    check_drp_config_columns, check_telemetry_columns,
)

file_urls = {
    "NEID": {
        "NATIVE_L2": "http://grinnell.as.arizona.edu/~rvdata/neid/neidL2_20231010T020006.fits",
    }
}


def download_file(url, filename):
    response = requests.get(url)
    response.raise_for_status()  # Check if the request was successful
    with open(filename, "wb") as file:
        file.write(response.content)


def download_files():
    native_l2_file = "neid_L2.fits"
    if not os.path.exists(native_l2_file):
        download_file(file_urls["NEID"]["NATIVE_L2"], native_l2_file)

    return native_l2_file


def test_neid():
    native_l2_file = download_files()

    # Check L2 - use auto-generated filename
    neidl2 = RV2.from_fits(native_l2_file, instrument="NEID")
    standard_l2_file = neidl2.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(standard_l2_file)), \
        f"L2 filename '{standard_l2_file}' does not match EPRV convention"
    assert standard_l2_file.startswith("neid_SL2_"), \
        f"L2 filename should start with 'neid_SL2_', got '{standard_l2_file}'"
    neidl2_obj = RV2.from_fits(standard_l2_file)

    check_l2_extensions(standard_l2_file)
    check_l2_header(neidl2_obj.headers["PRIMARY"])
    check_order_table_columns(standard_l2_file)
    check_receipt_columns(standard_l2_file)
    check_drp_config_columns(standard_l2_file)
    check_telemetry_columns(standard_l2_file)

    # Check L3 - use auto-generated filename
    neidl3 = RV3.from_fits(native_l2_file, instrument="NEID")
    standard_l3_file = neidl3.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(standard_l3_file)), \
        f"L3 filename '{standard_l3_file}' does not match EPRV convention"
    assert standard_l3_file.startswith("neid_SL3_"), \
        f"L3 filename should start with 'neid_SL3_', got '{standard_l3_file}'"
    neidl3_obj = RV3.from_fits(standard_l3_file)

    check_l3_extensions(standard_l3_file)
    check_l3_header(neidl3_obj.headers["PRIMARY"])
    check_order_table_columns(standard_l3_file)
    check_receipt_columns(standard_l3_file)
    check_drp_config_columns(standard_l3_file)

    # Check L4 - use auto-generated filename
    neidl4 = RV4.from_fits(native_l2_file, instrument="NEID")
    standard_l4_file = neidl4.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(standard_l4_file)), \
        f"L4 filename '{standard_l4_file}' does not match EPRV convention"
    assert standard_l4_file.startswith("neid_SL4_"), \
        f"L4 filename should start with 'neid_SL4_', got '{standard_l4_file}'"
    neidl4_obj = RV4.from_fits(standard_l4_file)

    check_l4_extensions(standard_l4_file)
    check_l4_header(neidl4_obj.headers["PRIMARY"])
    check_l4_rv_columns(standard_l4_file)
    check_receipt_columns(standard_l4_file)
    check_drp_config_columns(standard_l4_file)


def test_neid_benchmark(benchmark):
    # run test_kpf() once to download the files
    _ = download_files()
    # now run it again with benchmark
    benchmark(test_neid)


if __name__ == "__main__":
    test_neid()
