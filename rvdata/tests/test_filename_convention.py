# Tests for filename convention enforcement
import warnings
from collections import OrderedDict

import pytest

from rvdata.core.models.base import RVDataModel
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level3 import RV3
from rvdata.core.models.level4 import RV4


class TestFilenamePattern:
    """Test the filename pattern matching."""

    def test_valid_filenames(self):
        """Test that valid filenames match the pattern."""
        valid_filenames = [
            "kpf_SL2_20250208T045125.fits",
            "neid_SL3_20241022T123456.fits",
            "espresso_SL4_20220824T033756.fits",
            "HARPS_SL2_20200101T000000.fits",
        ]
        for filename in valid_filenames:
            assert RVDataModel.FILENAME_PATTERN.match(filename), f"{filename} should be valid"

    def test_invalid_filenames(self):
        """Test that invalid filenames don't match the pattern."""
        invalid_filenames = [
            "kpf_L2_20250208T045125.fits",  # Missing 'S'
            "kpf_SL2_standard.fits",  # Wrong datetime format
            "kpf_SL5_20250208T045125.fits",  # Invalid level (5)
            "kpf_SL2_2025-02-08T04:51:25.fits",  # Wrong datetime format (with dashes/colons)
            "kpf_SL2_20250208T045125.fit",  # Wrong extension
            "_SL2_20250208T045125.fits",  # Missing instrument
        ]
        for filename in invalid_filenames:
            assert not RVDataModel.FILENAME_PATTERN.match(filename), f"{filename} should be invalid"


class TestGenerateStandardFilename:
    """Test the generate_standard_filename method."""

    def test_generate_filename_level2(self):
        """Test filename generation for level 2 data."""
        obj = RV2()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = "KPF"
        obj.headers["PRIMARY"]["DATE-OBS"] = "2025-02-08T04:51:25.123"

        filename = obj.generate_standard_filename()
        assert filename == "kpf_SL2_20250208T045125.fits"

    def test_generate_filename_level3(self):
        """Test filename generation for level 3 data."""
        obj = RV3()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = "NEID"
        obj.headers["PRIMARY"]["DATE-OBS"] = "2024-10-22T12:34:56"

        filename = obj.generate_standard_filename()
        assert filename == "neid_SL3_20241022T123456.fits"

    def test_generate_filename_level4(self):
        """Test filename generation for level 4 data."""
        obj = RV4()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = "ESPRESSO"
        obj.headers["PRIMARY"]["DATE-OBS"] = "2022-08-24T03:37:56.276"

        filename = obj.generate_standard_filename()
        assert filename == "espresso_SL4_20220824T033756.fits"

    def test_generate_filename_with_tuple_headers(self):
        """Test filename generation when headers are (value, comment) tuples."""
        obj = RV2()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = ("KPF", "Instrument name")
        obj.headers["PRIMARY"]["DATE-OBS"] = ("2025-02-08T04:51:25", "Observation date")

        filename = obj.generate_standard_filename()
        assert filename == "kpf_SL2_20250208T045125.fits"

    def test_generate_filename_missing_instrume(self):
        """Test that missing INSTRUME raises ValueError."""
        obj = RV2()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["DATE-OBS"] = "2025-02-08T04:51:25"

        with pytest.raises(ValueError, match="INSTRUME"):
            obj.generate_standard_filename()

    def test_generate_filename_missing_date_obs(self):
        """Test that missing DATE-OBS raises ValueError."""
        obj = RV2()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = "KPF"

        with pytest.raises(ValueError, match="DATE-OBS"):
            obj.generate_standard_filename()


class TestValidateFilename:
    """Test the validate_filename method."""

    def test_validate_valid_filename(self):
        """Test validation of a valid filename."""
        obj = RV2()
        assert obj.validate_filename("kpf_SL2_20250208T045125.fits")

    def test_validate_invalid_filename(self):
        """Test validation of an invalid filename."""
        obj = RV2()
        assert not obj.validate_filename("kpf_L2_standard.fits")

    def test_validate_filename_with_path(self):
        """Test validation strips directory path."""
        obj = RV2()
        assert obj.validate_filename("/some/path/kpf_SL2_20250208T045125.fits")
        assert not obj.validate_filename("/some/path/kpf_L2_standard.fits")


class TestCheckFilenameConvention:
    """Test the check_filename_convention method."""

    def test_check_valid_filename_no_warning(self):
        """Test that valid filenames don't produce warnings."""
        obj = RV2()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = "KPF"
        obj.headers["PRIMARY"]["DATE-OBS"] = "2025-02-08T04:51:25"

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = obj.check_filename_convention("kpf_SL2_20250208T045125.fits")
            assert result is True
            assert len(w) == 0

    def test_check_invalid_filename_warning(self):
        """Test that invalid filenames produce warnings with suggestion."""
        obj = RV2()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = "KPF"
        obj.headers["PRIMARY"]["DATE-OBS"] = "2025-02-08T04:51:25"

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = obj.check_filename_convention("kpf_L2_standard.fits")
            assert result is False
            assert len(w) == 1
            assert "does not follow the EPRV naming convention" in str(w[0].message)
            assert "kpf_SL2_20250208T045125.fits" in str(w[0].message)

    def test_check_invalid_filename_no_headers(self):
        """Test warning when headers are missing for suggestion."""
        obj = RV2()
        # No PRIMARY header set

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = obj.check_filename_convention("bad_filename.fits")
            assert result is False
            assert len(w) == 1
            assert "does not follow the EPRV naming convention" in str(w[0].message)


class TestToFitsAutoFilename:
    """Test that to_fits() auto-generates correct filenames."""

    def test_to_fits_auto_filename_generation(self):
        """Test that to_fits() with no args uses generate_standard_filename()."""
        obj = RV2()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = "TEST"
        obj.headers["PRIMARY"]["DATE-OBS"] = "2025-02-08T04:51:25"

        # Verify generate_standard_filename() produces the expected result
        expected_filename = "test_SL2_20250208T045125.fits"
        assert obj.generate_standard_filename() == expected_filename

    def test_to_fits_returns_filename(self):
        """Test that to_fits() returns the filename it wrote to."""
        import tempfile
        import os

        obj = RV2()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = "TEST"
        obj.headers["PRIMARY"]["DATE-OBS"] = "2025-02-08T04:51:25"
        obj.extensions["PRIMARY"] = None

        with tempfile.TemporaryDirectory() as tmpdir:
            # Use explicit path to avoid changing directory (which breaks git repo lookup)
            expected_filename = "test_SL2_20250208T045125.fits"
            full_path = os.path.join(tmpdir, expected_filename)
            filename = obj.to_fits(full_path)
            assert filename == full_path
            assert os.path.exists(filename)
            assert os.path.basename(filename) == expected_filename

    def test_to_fits_with_explicit_filename_returns_that_filename(self):
        """Test that to_fits() with explicit filename returns that filename."""
        import tempfile
        import os

        obj = RV2()
        obj.headers["PRIMARY"] = OrderedDict()
        obj.headers["PRIMARY"]["INSTRUME"] = "TEST"
        obj.headers["PRIMARY"]["DATE-OBS"] = "2025-02-08T04:51:25"
        obj.extensions["PRIMARY"] = None

        with tempfile.TemporaryDirectory() as tmpdir:
            explicit_name = os.path.join(tmpdir, "my_custom_name.fits")
            filename = obj.to_fits(explicit_name)
            assert filename == explicit_name
            assert os.path.exists(filename)
