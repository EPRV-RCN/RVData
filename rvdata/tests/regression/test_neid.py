import requests
import os
from rvdata.core.models.level2 import RV2
from rvdata.tests.regression.compliance import check_l2_extensions  # , check_l2_header


file_urls = {
    "NEID": [
        "https://zenodo.org/records/15009170/files/KP.20240715.46231.07_L1.fits?download=1", #just needs one NEID L2
    ]
}


def download_file(url, filename):
    response = requests.get(url)
    response.raise_for_status()  # Check if the request was successful
    with open(filename, "wb") as file:
        file.write(response.content)


def test_neid():
    l2file = "neid_L2.fits"
    if not os.path.exists(l2file):
        download_file(file_urls["NEID"][0], l2file)

    neid_standardl2 = RV2.from_fits(l2file, instrument="NEID")
    standard_out = "./neid_L2_standard.fits"
    neid_standardl2.to_fits(standard_out)
    # l2_standard = RV2.from_fits(standard_out)

    check_l2_extensions(standard_out)
    # check_l2_header(l2_standard.headers['PRIMARY'])


def test_neid_benchmark(benchmark):
    # run test_kpf() once to download the files
    test_neid()
    # now run it again with benchmark
    benchmark(test_neid)


if __name__ == "__main__":
    test_neid()
