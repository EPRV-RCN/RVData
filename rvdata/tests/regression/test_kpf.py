import requests
import os
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level4 import RV4

# from rvdata.instruments.kpf.level2 import KPFRV2
from rvdata.tests.regression.compliance import check_l2_extensions, check_l2_header
from rvdata.tests.regression.compliance import check_l4_extensions, check_l4_header


file_urls = {
    "KPF": [
        "https://zenodo.org/records/15047493/files/KP.20250208.17485.59.fits?download=1",
        "https://zenodo.org/records/15047493/files/KP.20250208.17485.59_L1.fits?download=1",
        "https://zenodo.org/records/15047493/files/KP.20241022.41656.30_L2.fits?download=1",
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
    kpf2 = RV2.from_fits(l1file, l0file=l0file, instrument="KPF")
    l2_standard = "./kpf_L2_standard.fits"
    kpf2.to_fits(l2_standard)
    l2_obj = RV2.from_fits(l2_standard)

    check_l2_extensions(l2_standard)
    check_l2_header(l2_obj.headers['PRIMARY'])

    kpf4 = RV4.from_fits(l2file, instrument="KPF")
    l4_standard = "./kpf_L4_standard.fits"
    kpf4.to_fits(l4_standard)
    l4_obj = RV4.from_fits(l4_standard)

    check_l4_extensions(l4_standard)
    check_l4_header(l4_obj.headers['PRIMARY'])


def test_kpf_benchmark(benchmark):
    # run test_kpf() once to download the files
    _ = download_files
    # now run it again with benchmark
    benchmark(test_kpf)


if __name__ == "__main__":
    test_kpf()
