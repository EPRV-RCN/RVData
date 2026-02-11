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


def test_kpf():
    l0file, l1file, l2file = download_files()

    # Check L2 - use auto-generated filename
    kpf2 = RV2.from_fits(l1file, l0file=l0file, instrument="KPF")
    l2_standard = kpf2.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l2_standard)), \
        f"L2 filename '{l2_standard}' does not match EPRV convention"
    assert l2_standard.startswith("kpf_SL2_"), \
        f"L2 filename should start with 'kpf_SL2_', got '{l2_standard}'"
    l2_obj = RV2.from_fits(l2_standard)

    check_l2_extensions(l2_standard)
    check_l2_header(l2_obj.headers['PRIMARY'])
    check_order_table_columns(l2_standard)
    check_receipt_columns(l2_standard)
    check_drp_config_columns(l2_standard)
    check_telemetry_columns(l2_standard)

    # Check L3
    # Note: L3 creation requires an RVData-standard L2 file (created above).
    # Native KPF L2 files must first be converted using KPFRV2.from_fits().
    kpf3 = RV3.from_fits(l2_standard, instrument="KPF")
    l3_standard = kpf3.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l3_standard)), \
        f"L3 filename '{l3_standard}' does not match EPRV convention"
    assert l3_standard.startswith("kpf_SL3_"), \
        f"L3 filename should start with 'kpf_SL3_', got '{l3_standard}'"
    l3_obj = RV3.from_fits(l3_standard)

    check_l3_extensions(l3_standard)
    check_l3_header(l3_obj.headers['PRIMARY'])
    check_order_table_columns(l3_standard)
    check_receipt_columns(l3_standard)
    check_drp_config_columns(l3_standard)

    # Check L4 - use auto-generated filename
    kpf4 = RV4.from_fits(l2file, instrument="KPF")
    l4_standard = kpf4.to_fits()  # Auto-generate filename
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l4_standard)), \
        f"L4 filename '{l4_standard}' does not match EPRV convention"
    assert l4_standard.startswith("kpf_SL4_"), \
        f"L4 filename should start with 'kpf_SL4_', got '{l4_standard}'"
    l4_obj = RV4.from_fits(l4_standard)

    check_l4_extensions(l4_standard)
    check_l4_header(l4_obj.headers['PRIMARY'])
    check_l4_rv_columns(l4_standard)
    check_receipt_columns(l4_standard)
    check_drp_config_columns(l4_standard)


def test_kpf_benchmark(benchmark):
    # run test_kpf() once to download the files
    _ = download_files()
    # now run it again with benchmark
    benchmark(test_kpf)


if __name__ == "__main__":
    test_kpf()
