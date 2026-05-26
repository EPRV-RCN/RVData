#!/usr/bin/env python
# coding: utf-8

import requests
import os
from astropy.io import fits
from rvdata.core.models.base import RVDataModel
from rvdata.core.models.level3 import RV3
from rvdata.instruments.maroonx.level2 import MAROONXRV2
from rvdata.instruments.maroonx.level4 import MAROONXRV4

from rvdata.tests.regression.compliance import (
    check_l2_compliance, check_l3_compliance, check_l4_compliance,
)

pathname = "http://grinnell.as.arizona.edu/~rvdata/maroonx/"
file_urls = {
    "MAROONX": {
        "MX_OP": pathname + "20240923T092056Z_SOOOE_x_0040.hd5",
        "FLAT": pathname + "20240916T19_masterflat__FFFFF_x__blaze.hd5",
        "SERVAL_BLUE": pathname + "tauCet_activity_results_blue.pkl",
        "SERVAL_RED": pathname + "tauCet_activity_results_red.pkl"
    }
}


def download_file(url, filename):
    response = requests.get(url)
    response.raise_for_status()  # Check if the request was successful
    with open(filename, "wb") as file:
        file.write(response.content)


def download_files():
    file_mx = "20240923T092056Z_SOOOE_x_0040.hd5"
    flat_file = "20240916T19_masterflat__FFFFF_x__blaze.hd5"
    rv_blue = "tauCet_activity_results_blue.pkl"
    rv_red = "tauCet_activity_results_red.pkl"

    if not os.path.exists(file_mx):
        download_file(file_urls["MAROONX"]['MX_OP'], file_mx)
    if not os.path.exists(flat_file):
        download_file(file_urls["MAROONX"]['FLAT'], flat_file)
    if not os.path.exists(rv_blue):
        download_file(file_urls["MAROONX"]['SERVAL_BLUE'], rv_blue)
    if not os.path.exists(rv_red):
        download_file(file_urls["MAROONX"]['SERVAL_RED'], rv_red)

    return file_mx, flat_file, rv_blue, rv_red


def test_maroonx():
    file_mx, flat_file, rv_blue, rv_red = download_files()

    # Check L2
    mx2 = MAROONXRV2()
    mx2.createL2(file_mx, flat_file)
    l2_standardB, l2_standardR = mx2.write_camera_fits(file_mx, out_dir='.')
    # Check red camera
    match_r = RVDataModel.FILENAME_PATTERN.match(
        os.path.basename(l2_standardR))
    assert match_r, f"L2 filename '{l2_standardR}' does \
          not match EPRV convention"
    assert l2_standardR.startswith("maroonxred_SL2_"), \
        f"L2 filename should start with 'maroonxred_SL2_', \
            got '{l2_standardR}'"
    check_l2_compliance(l2_standardR)
    # Check blue camera
    match_b = RVDataModel.FILENAME_PATTERN.match(
        os.path.basename(l2_standardB))
    assert match_b, f"L2 filename '{l2_standardB}' does \
          not match EPRV convention"
    assert l2_standardB.startswith("maroonxblue_SL2_"), \
        f"L2 filename should start with 'maroonxblue_SL2_', \
            got '{l2_standardB}'"
    check_l2_compliance(l2_standardB)

    # Check L3
    # Note: L3 creation requires an RVData-standard L2 file (created above).
    mx3 = RV3.from_fits(l2_standardR, instrument="MAROONX")
    l3_standard = mx3.to_fits()
    assert RVDataModel.FILENAME_PATTERN.match(os.path.basename(l3_standard)), \
        f"L3 filename '{l3_standard}' does not match EPRV convention"
    assert l3_standard.startswith("maroonxred_SL3_"), \
        f"L3 filename should start with 'maroonxred_SL3_', got '{l3_standard}'"
    check_l3_compliance(l3_standard)

    # Check L4
    mx4 = MAROONXRV4()
    mx4.createL4(fits.open(l2_standardR), rv_red, channel="RED")
    l4_standard = mx4.to_fits()
    match4 = RVDataModel.FILENAME_PATTERN.match(os.path.basename(l4_standard))
    assert match4, (
        f"L4 filename '{l4_standard}' does not match EPRV convention"
    )
    assert l4_standard.startswith("maroonxred_SL4_"), \
        f"L4 filename should start with 'maroonxred_SL4_', got '{l4_standard}'"
    check_l4_compliance(l4_standard)


def test_MX_benchmark(benchmark):
    # run test_maroonx() once to download the files
    _ = download_files()
    # now run it again with benchmark
    benchmark(test_maroonx)


if __name__ == "__main__":
    test_maroonx()
