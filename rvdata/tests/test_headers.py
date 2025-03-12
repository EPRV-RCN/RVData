# Tests for core.tools.headers module
from rvdata.core.tools.headers import to_ascii_safe


def test_to_ascii_safe_clean():
    assert to_ascii_safe("abc") == "abc"


def test_to_ascii_safe_unicode():
    assert to_ascii_safe("äbc") == "abc"
