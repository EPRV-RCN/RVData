# Tests for core.tools.headers module
import warnings

from astropy.io.fits.card import Undefined

from rvdata.core.tools.headers import parse_value_to_datatype, to_ascii_safe


def test_to_ascii_safe_clean():
    assert to_ascii_safe("abc") == "abc"


def test_to_ascii_safe_unicode():
    assert to_ascii_safe("äbc") == "abc"


def test_parse_value_to_datatype_none():
    """Test that None values are handled correctly."""
    result = parse_value_to_datatype("TEST", "uint", None)
    assert result == (None, "")


def test_parse_value_to_datatype_undefined_object():
    """Test that astropy Undefined objects are handled correctly."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # Turn warnings into errors
        result = parse_value_to_datatype("TEST", "uint", Undefined())
    assert result == (None, "")


def test_parse_value_to_datatype_undefined_string():
    """Test that 'UNDEFINED' string is handled correctly."""
    result = parse_value_to_datatype("TEST", "uint", "UNDEFINED")
    assert result == (None, "")
    result = parse_value_to_datatype("TEST", "uint", "undefined")
    assert result == (None, "")


def test_parse_value_to_datatype_empty_string():
    """Test that empty and whitespace strings are handled correctly."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # Turn warnings into errors
        result = parse_value_to_datatype("TEST", "uint", "")
        assert result == (None, "")
        result = parse_value_to_datatype("TEST", "uint", " ")
        assert result == (None, "")
        result = parse_value_to_datatype("TEST", "uint", "   ")
        assert result == (None, "")


def test_parse_value_to_datatype_valid_values():
    """Test that valid values are still converted correctly."""
    assert parse_value_to_datatype("TEST", "uint", 42) == (42, "")
    assert parse_value_to_datatype("TEST", "uint", "42") == (42, "")
    assert parse_value_to_datatype("TEST", "float", 3.14) == (3.14, "")
    assert parse_value_to_datatype("TEST", "string", "hello") == ("hello", "")
    assert parse_value_to_datatype("TEST", "boolean", True) == (True, "")
    assert parse_value_to_datatype("TEST", "boolean", "True") == (True, "")
