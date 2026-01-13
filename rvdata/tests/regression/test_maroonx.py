#!/usr/bin/env python
# coding: utf-8

import requests
import os
from astropy.io import fits
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level4 import RV4
from rvdata.instruments.maroonx.level2 import MAROONXRV2
from rvdata.instruments.maroonx.level4 import MAROONXRV4

from rvdata.tests.regression.compliance import check_l2_extensions, check_l2_header
from rvdata.tests.regression.compliance import check_l4_extensions, check_l4_header

file_urls = {
    "MAROONX": {
        "MX_OP": "http://grinnell.as.arizona.edu/~rvdata/maroonx/20240923T092056Z_SOOOE_x_0040.hd5",
        "FLAT": "http://grinnell.as.arizona.edu/~rvdata/maroonx/20240916T19_masterflat__FFFFF_x__blaze.hd5",
        "SERVAL_BLUE": "http://grinnell.as.arizona.edu/~rvdata/maroonx/tauCet_activity_results_blue.pkl",
        "SERVAL_RED": "http://grinnell.as.arizona.edu/~rvdata/maroonx/tauCet_activity_results_red.pkl"
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
    timestamp = os.path.basename(file_mx).split("_")[0].replace("Z", "")

    rv2_obj = MAROONXRV2()
    rv2_obj.createL2(file_mx, flat_file)
    blue_SL2_fits, red_SL2_fits = rv2_obj.write_camera_fits(file_mx,
                                                            out_dir='.')
    l2_obj_blue = RV2.from_fits(blue_SL2_fits)
    l2_obj_red = RV2.from_fits(red_SL2_fits)

    # Check blue camera
    check_l2_extensions(blue_SL2_fits)
    check_l2_header(l2_obj_blue.headers['PRIMARY'])

    # Check red camera
    check_l2_extensions(red_SL2_fits)
    check_l2_header(l2_obj_red.headers['PRIMARY'])

    rv4_obj = MAROONXRV4()
    # Check blue camera
    rv4_obj.createL4(fits.open(l2_obj_blue.filename), rv_blue, channel="BLUE")
    l4_standard_blue = f"./MAROONXBLUE_SL4_{timestamp}.fits"
    rv4_obj.to_fits(l4_standard_blue)
    l4_obj_blue = RV4.from_fits(l4_standard_blue)

    check_l4_extensions(l4_standard_blue)
    check_l4_header(l4_obj_blue.headers['PRIMARY'])

    # Check red camera
    rv4_obj.createL4(fits.open(l2_obj_red.filename), rv_red, channel="RED")
    l4_standard_red = f"./MAROONXRED_SL4_{timestamp}.fits"
    rv4_obj.to_fits(l4_standard_red)
    l4_obj_red = RV4.from_fits(l4_standard_red)

    check_l4_extensions(l4_standard_red)
    check_l4_header(l4_obj_red.headers['PRIMARY'])


def test_MX_benchmark(benchmark):
    # run test_maroonx() once to download the files
    _ = download_files()
    # now run it again with benchmark
    benchmark(test_maroonx)


if __name__ == "__main__":
    test_maroonx()
