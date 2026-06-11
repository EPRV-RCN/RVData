"""Tests for MinBitDepth enforcement in set_data() and compliance checks."""

import os
import tempfile
import warnings

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from rvdata.core.models.level2 import RV2
from rvdata.tests.regression.compliance import (
    _check_image_bitdepth,
    _check_table_column_bitdepth,
)


class TestSetDataUpcast:
    """Test that set_data() upcasts under-precision ImageHDU data."""

    def test_float32_wave_upcasted_to_float64(self):
        """TRACE1_WAVE requires MinBitDepth=64; float32 should be upcasted."""
        rv2 = RV2()
        data = np.ones((10, 100), dtype=np.float32)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rv2.set_data("TRACE1_WAVE", data)
            assert len(w) == 1
            assert "MinBitDepth=64" in str(w[0].message)
            assert "Upcasting" in str(w[0].message)
        assert rv2.data["TRACE1_WAVE"].dtype == np.float64

    def test_float64_wave_no_warning(self):
        """float64 data should not trigger a warning."""
        rv2 = RV2()
        data = np.ones((10, 100), dtype=np.float64)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rv2.set_data("TRACE1_WAVE", data)
            upcast_warnings = [x for x in w if "MinBitDepth" in str(x.message)]
            assert len(upcast_warnings) == 0
        assert rv2.data["TRACE1_WAVE"].dtype == np.float64

    def test_flux_no_min_bit_depth(self):
        """TRACE1_FLUX has no MinBitDepth; float32 should be accepted as-is."""
        rv2 = RV2()
        data = np.ones((10, 100), dtype=np.float32)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rv2.set_data("TRACE1_FLUX", data)
            upcast_warnings = [x for x in w if "MinBitDepth" in str(x.message)]
            assert len(upcast_warnings) == 0
        assert rv2.data["TRACE1_FLUX"].dtype == np.float32

    def test_quality_uint8_satisfies_min8(self):
        """TRACE1_QUALITY requires MinBitDepth=8; uint8 should pass."""
        rv2 = RV2()
        rv2.create_extension("TRACE1_QUALITY", "ImageHDU")
        data = np.zeros((10, 100), dtype=np.uint8)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rv2.set_data("TRACE1_QUALITY", data)
            upcast_warnings = [x for x in w if "MinBitDepth" in str(x.message)]
            assert len(upcast_warnings) == 0
        assert rv2.data["TRACE1_QUALITY"].dtype == np.uint8

    def test_empty_array_no_warning(self):
        """Empty arrays should not trigger a warning."""
        rv2 = RV2()
        data = np.array([], dtype=np.float32)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rv2.set_data("TRACE1_WAVE", data)
            upcast_warnings = [x for x in w if "MinBitDepth" in str(x.message)]
            assert len(upcast_warnings) == 0

    def test_uint8_upcasted_preserves_signedness(self):
        """Unsigned integer data should upcast to unsigned target dtype."""
        rv2 = RV2()
        rv2.create_extension("TRACE1_QUALITY", "ImageHDU")
        # uint8 has 8 bits, but if MinBitDepth were higher it should stay unsigned.
        # Use BJD_TDB (MinBitDepth=64) with uint8 data to test signedness preservation.
        data = np.ones((10, 100), dtype=np.uint8)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # BJD_TDB requires MinBitDepth=64, uint8 should upcast to uint64
            rv2.set_data("BJD_TDB", data)
            upcast_warnings = [x for x in w if "MinBitDepth" in str(x.message)]
            assert len(upcast_warnings) == 1
        assert rv2.data["BJD_TDB"].dtype == np.uint64

    def test_int8_upcasted_preserves_signedness(self):
        """Signed integer data should upcast to signed target dtype."""
        rv2 = RV2()
        data = np.ones((10, 100), dtype=np.int8)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rv2.set_data("BJD_TDB", data)
            upcast_warnings = [x for x in w if "MinBitDepth" in str(x.message)]
            assert len(upcast_warnings) == 1
        assert rv2.data["BJD_TDB"].dtype == np.int64

    def test_multiplicity_trace2_wave_upcasted(self):
        """TRACE2_WAVE should inherit TRACE1_WAVE's MinBitDepth=64."""
        rv2 = RV2()
        rv2.create_extension("TRACE2_WAVE", "ImageHDU")
        data = np.ones((10, 100), dtype=np.float32)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rv2.set_data("TRACE2_WAVE", data)
            assert any("MinBitDepth=64" in str(x.message) for x in w)
        assert rv2.data["TRACE2_WAVE"].dtype == np.float64


class TestComplianceImageBitdepth:
    """Test _check_image_bitdepth compliance helper."""

    def _write_l2_with_wave_dtype(self, filepath, dtype):
        """Write a minimal L2-like FITS file with TRACE1_WAVE at given dtype."""
        rv2 = RV2()
        wave_data = np.linspace(5000, 6000, 100, dtype=np.float64).reshape(1, 100)
        flux_data = np.ones((1, 100), dtype=np.float32)
        rv2.set_data("TRACE1_WAVE", wave_data)
        rv2.set_data("TRACE1_FLUX", flux_data)
        rv2.set_data("TRACE1_VAR", flux_data)
        rv2.set_data("TRACE1_BLAZE", flux_data)
        rv2.set_data("BARYCORR_KMS", np.array([1.0]))
        rv2.set_data("BARYCORR_Z", np.array([1e-5]))
        rv2.set_data("BJD_TDB", np.array([2460000.0], dtype=np.float64))
        rv2.to_fits(out_filename=filepath)

        # If we need to downcast for testing, rewrite just that extension
        if dtype != np.float64:
            with fits.open(filepath, mode='update') as hdul:
                hdul["TRACE1_WAVE"].data = wave_data.astype(dtype)

    def test_compliant_file_passes(self):
        """A file with float64 TRACE1_WAVE should pass."""
        import importlib
        import pandas as pd

        extdf = pd.read_csv(
            importlib.resources.files("rvdata.core.models.config")
            / "L2-extensions.csv"
        )
        with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as f:
            filepath = f.name
        try:
            self._write_l2_with_wave_dtype(filepath, np.float64)
            with fits.open(filepath) as hdul:
                _check_image_bitdepth(hdul, extdf)  # should not raise
        finally:
            os.unlink(filepath)

    def test_noncompliant_file_fails(self):
        """A file with float32 TRACE1_WAVE should fail the check."""
        import importlib
        import pandas as pd

        extdf = pd.read_csv(
            importlib.resources.files("rvdata.core.models.config")
            / "L2-extensions.csv"
        )
        with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as f:
            filepath = f.name
        try:
            self._write_l2_with_wave_dtype(filepath, np.float32)
            with fits.open(filepath) as hdul:
                with pytest.raises(AssertionError, match="MinBitDepth=64"):
                    _check_image_bitdepth(hdul, extdf)
        finally:
            os.unlink(filepath)


class TestComplianceTableColumnBitdepth:
    """Test _check_table_column_bitdepth compliance helper."""

    def test_float64_columns_pass(self):
        """RV1 table with float64 BJD_TDB should pass."""
        import importlib
        import pandas as pd

        coldf = pd.read_csv(
            importlib.resources.files("rvdata.core.models.config")
            / "L4-RV_TABLE-columns.csv"
        )
        table = Table({
            "BJD_TDB": np.array([2460000.0], dtype=np.float64),
            "RV": np.array([1.0], dtype=np.float32),
            "RV_ERR": np.array([0.1], dtype=np.float32),
            "BERV": np.array([10.0], dtype=np.float32),
            "WAVE_START": np.array([5000.0], dtype=np.float64),
            "WAVE_END": np.array([6000.0], dtype=np.float64),
        })
        hdu = fits.BinTableHDU(data=table, name="RV1")
        hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
        _check_table_column_bitdepth(hdul, coldf, "RV1")  # should not raise

    def test_float32_bjd_fails(self):
        """RV1 table with float32 BJD_TDB should fail."""
        import importlib
        import pandas as pd

        coldf = pd.read_csv(
            importlib.resources.files("rvdata.core.models.config")
            / "L4-RV_TABLE-columns.csv"
        )
        table = Table({
            "BJD_TDB": np.array([2460000.0], dtype=np.float32),
            "RV": np.array([1.0], dtype=np.float32),
            "RV_ERR": np.array([0.1], dtype=np.float32),
            "BERV": np.array([10.0], dtype=np.float32),
            "WAVE_START": np.array([5000.0], dtype=np.float64),
            "WAVE_END": np.array([6000.0], dtype=np.float64),
        })
        hdu = fits.BinTableHDU(data=table, name="RV1")
        hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
        with pytest.raises(AssertionError, match="MinBitDepth=64"):
            _check_table_column_bitdepth(hdul, coldf, "RV1")
