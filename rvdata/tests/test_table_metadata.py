"""
Tests that BinTableHDU column metadata (TUNIT, TFORM, TDISP) survives
a read/write round-trip when stored as astropy Table internally.
"""

import os
import tempfile

import numpy as np
from astropy.io import fits
from astropy.table import Table

from rvdata.core.models.level4 import RV4


def _make_fits_with_metadata(path):
    """Create a minimal FITS file with a BinTableHDU that has TUNIT/TFORM/TDISP."""
    # Primary HDU with required keywords
    phdu = fits.PrimaryHDU()
    phdu.header["DATALVL"] = "L4"
    phdu.header["INSTRUME"] = "TEST"
    phdu.header["DATE-OBS"] = "2024-01-01T00:00:00"

    # Build a BinTableHDU with explicit column metadata
    col1 = fits.Column(
        name="BJD_TDB", format="D", unit="d", disp="F15.7",
        array=np.array([2460000.5, 2460001.5]))
    col2 = fits.Column(
        name="RV", format="E", unit="m/s", disp="F10.3",
        array=np.array([1.23, 4.56]))
    col3 = fits.Column(
        name="RV_ERR", format="E", unit="m/s", disp="F10.3",
        array=np.array([0.01, 0.02]))

    rv_hdu = fits.BinTableHDU.from_columns([col1, col2, col3], name="RV1")

    # EXT_DESCRIPT table (required by RV4)
    ext_col1 = fits.Column(
        name="Name", format="20A",
        array=np.array(["RV1", "EXT_DESCRIPT", "INSTRUMENT_HEADER",
                        "RECEIPT", "DRP_CONFIG"]))
    ext_col2 = fits.Column(
        name="Description", format="60A",
        array=np.array(["RV table", "Extension descriptions",
                        "Instrument header", "Receipt", "DRP config"]))
    ext_hdu = fits.BinTableHDU.from_columns([ext_col1, ext_col2], name="EXT_DESCRIPT")

    # INSTRUMENT_HEADER as ImageHDU
    inst_hdu = fits.ImageHDU(data=np.zeros((1,), dtype=np.float32), name="INSTRUMENT_HEADER")

    # RECEIPT table
    receipt_col = fits.Column(name="Time", format="30A", array=np.array([""]))
    receipt_hdu = fits.BinTableHDU.from_columns([receipt_col], name="RECEIPT")

    # DRP_CONFIG table
    drp_col = fits.Column(name="Keyword", format="20A", array=np.array([""]))
    drp_hdu = fits.BinTableHDU.from_columns([drp_col], name="DRP_CONFIG")

    hdul = fits.HDUList([phdu, ext_hdu, inst_hdu, rv_hdu, receipt_hdu, drp_hdu])
    hdul.writeto(path, overwrite=True)
    hdul.close()


def test_bintable_metadata_roundtrip():
    """TUNIT, TFORM, and TDISP survive a read/write cycle."""
    with tempfile.TemporaryDirectory() as tmpdir:
        src = os.path.join(tmpdir, "test_SL4_20240101T000000.fits")
        dst = os.path.join(tmpdir, "out_SL4_20240101T000000.fits")

        _make_fits_with_metadata(src)

        # Read original metadata
        with fits.open(src) as hdul_orig:
            orig_header = hdul_orig["RV1"].header
            orig_tunits = {i: orig_header.get(f"TUNIT{i}") for i in range(1, 4)}
            orig_tforms = {i: orig_header.get(f"TFORM{i}") for i in range(1, 4)}
            orig_tdisps = {i: orig_header.get(f"TDISP{i}") for i in range(1, 4)}

        # Round-trip through RV4
        obj = RV4.from_fits(src)
        obj.to_fits(dst)

        # Read back and compare
        with fits.open(dst) as hdul_out:
            out_header = hdul_out["RV1"].header
            for i in range(1, 4):
                assert out_header.get(f"TUNIT{i}") == orig_tunits[i], \
                    f"TUNIT{i} mismatch: {out_header.get(f'TUNIT{i}')} != {orig_tunits[i]}"
                assert out_header.get(f"TFORM{i}") == orig_tforms[i], \
                    f"TFORM{i} mismatch: {out_header.get(f'TFORM{i}')} != {orig_tforms[i]}"
                assert out_header.get(f"TDISP{i}") == orig_tdisps[i], \
                    f"TDISP{i} mismatch: {out_header.get(f'TDISP{i}')} != {orig_tdisps[i]}"


def test_internal_data_is_astropy_table():
    """Verify that BinTableHDU data is stored as astropy Table, not DataFrame."""
    with tempfile.TemporaryDirectory() as tmpdir:
        src = os.path.join(tmpdir, "test_SL4_20240101T000000.fits")
        _make_fits_with_metadata(src)

        obj = RV4.from_fits(src)

        for ext_name, ext_type in obj.extensions.items():
            if ext_type == "BinTableHDU":
                assert isinstance(obj.data[ext_name], Table), \
                    f"Extension {ext_name} should be astropy Table, got {type(obj.data[ext_name])}"
