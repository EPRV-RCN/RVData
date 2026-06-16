"""Tests for RECEIPT extension population and round-tripping.

Covers the three issues raised in EPRV-RCN/RVData#168:
  (a) row keys written by receipt_add_entry match the BASE-RECEIPT-columns.csv
      schema (UPPERCASE: TIME, CODE_RELEASE, BRANCH_NAME, COMMIT_HASH,
      FUNCTION, ARGS, STATUS),
  (b) self.receipt entries actually land in the RECEIPT BinTableHDU on disk
      (the sync is done in to_fits via _sync_receipt_to_extension),
  (c) RECEIPT extension columns are string-typed even when no entries have
      been added (empty-init uses object dtype).

Also covers the ``@receipt_logged`` decorator introduced in PR #170:
  - successful calls log STATUS=PASS with a key=value ARGS string,
  - exceptions log STATUS=FAIL and are re-raised unchanged,
  - convert_level2_to_level3 is auto-instrumented via the decorator.
"""

import os
import tempfile

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from rvdata.core.models.base import receipt_logged
from rvdata.core.models.definitions import BASE_RECEIPT_COLUMNS
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level3 import RV3

EXPECTED_COLUMNS = BASE_RECEIPT_COLUMNS["Name"].tolist()


def _byte_string_dtype(col):
    """Astropy serializes string columns as fixed-width bytes (S<n>)."""
    return col.dtype.kind in ("S", "U", "O")


def test_receipt_add_entry_uses_csv_column_names():
    rv2 = RV2()
    rv2.receipt_add_entry("foo", "bar=1", "PASS")
    assert list(rv2.receipt.columns) == EXPECTED_COLUMNS
    row = rv2.receipt.iloc[-1]
    assert row["FUNCTION"] == "foo"
    assert row["ARGS"] == "bar=1"
    assert row["STATUS"] == "PASS"


def test_receipt_extension_empty_columns_are_string_typed():
    """RECEIPT is initialised string-typed even before any entries are added."""
    rv2 = RV2()
    receipt_ext = rv2.data["RECEIPT"]
    assert isinstance(receipt_ext, Table)
    assert list(receipt_ext.colnames) == EXPECTED_COLUMNS
    for col in EXPECTED_COLUMNS:
        assert not np.issubdtype(receipt_ext[col].dtype, np.floating), (
            f"RECEIPT column {col!r} initialised as floating dtype "
            f"{receipt_ext[col].dtype}; expected string/object dtype"
        )


def test_receipt_roundtrips_to_fits():
    rv2 = RV2()
    rv2.set_header("PRIMARY", {"INSTRUME": "TEST", "DATE-OBS": "2026-01-01T00:00:00"})
    rv2.receipt_add_entry("user_step", "x=1", "PASS")

    with tempfile.TemporaryDirectory() as tmp:
        fn = os.path.join(tmp, "test_SL2_20260101T000000.fits")
        rv2.to_fits(out_filename=fn)
        with fits.open(fn) as hdul:
            assert "RECEIPT" in [hdu.name for hdu in hdul]
            receipt = Table.read(hdul["RECEIPT"])

    assert list(receipt.colnames) == EXPECTED_COLUMNS
    # Two rows: the explicit user_step + the implicit to_fits entry.
    assert len(receipt) == 2
    functions = [str(v) for v in receipt["FUNCTION"]]
    statuses = [str(v) for v in receipt["STATUS"]]
    args = [str(v) for v in receipt["ARGS"]]
    assert "user_step" in functions
    assert "to_fits" in functions
    assert set(statuses) == {"PASS"}
    # to_fits should log its output path as a key=value pair (PR #170 fix).
    to_fits_args = args[functions.index("to_fits")]
    assert to_fits_args.startswith("out_filepath="), (
        f"to_fits ARGS should be 'out_filepath=<path>', got {to_fits_args!r}"
    )
    # Every column should be string-typed (bytes/str/object — never float).
    for col in EXPECTED_COLUMNS:
        assert _byte_string_dtype(receipt[col]), (
            f"RECEIPT column {col!r} has non-string dtype {receipt[col].dtype}"
        )


def test_receipt_logged_decorator_logs_pass_with_keyvalue_args():
    """Decorated methods record a PASS entry with key=value ARGS formatting."""

    class _Stub(RV2):
        @receipt_logged
        def custom_step(self, alpha, beta=7):
            return alpha + beta

    stub = _Stub()
    result = stub.custom_step(3, beta=4)
    assert result == 7
    row = stub.receipt.iloc[-1]
    assert row["FUNCTION"] == "custom_step"
    assert row["STATUS"] == "PASS"
    assert "alpha=3" in row["ARGS"]
    assert "beta=4" in row["ARGS"]


def test_receipt_logged_decorator_logs_fail_and_reraises():
    """An exception in a decorated method records FAIL and propagates."""

    class _Boom(RV2):
        @receipt_logged
        def kaboom(self, msg):
            raise ValueError(msg)

    boom = _Boom()
    with pytest.raises(ValueError, match="nope"):
        boom.kaboom("nope")
    row = boom.receipt.iloc[-1]
    assert row["FUNCTION"] == "kaboom"
    assert row["STATUS"] == "FAIL"
    assert "msg=nope" in row["ARGS"]


def test_receipt_logged_skips_self_and_cls():
    """The ARGS string never includes 'self' or 'cls'."""

    class _NoSelf(RV2):
        @receipt_logged
        def echo(self, x):
            return x

    obj = _NoSelf()
    obj.echo("hello")
    args = obj.receipt.iloc[-1]["ARGS"]
    assert "self" not in args
    assert args == "x=hello"


def test_convert_level2_to_level3_logs_receipt_entry():
    """RV3.convert_level2_to_level3 is decorated; entry appears on success."""
    rv2 = RV2()
    rv2.set_header("PRIMARY", {
        "INSTRUME": "TEST",
        "DATE-OBS": "2026-01-01T00:00:00",
        "DATALVL": "L2",
        "NUMTRACE": 0,
    })
    rv2.receipt_add_entry("prior_l2_step", "k=v", "PASS")

    rv3 = RV3()
    # convert_level2_to_level3 hits the stitching path which needs an
    # INSTRUME-specific config we don't have here; we only care that the
    # decorator runs around it. Catch the inevitable error from the
    # missing config and verify the FAIL entry was logged.
    try:
        rv3.convert_level2_to_level3(rv2)
    except Exception:
        pass

    functions = list(rv3.receipt["FUNCTION"])
    assert "convert_level2_to_level3" in functions, functions
    row = rv3.receipt[rv3.receipt["FUNCTION"] == "convert_level2_to_level3"].iloc[-1]
    # l2obj should render as <RV2> (default object repr is collapsed).
    assert "l2obj=<RV2>" in row["ARGS"], row["ARGS"]


def test_receipt_info_does_not_raise():
    """receipt_info historically referenced stale Title Case columns;
    after the PR #170 fix it should print without raising."""
    rv2 = RV2()
    rv2.receipt_add_entry("step_a", "x=1", "PASS")
    # Should not raise (previously KeyError on 'Time'/'Module_Name'/'Status').
    rv2.receipt_info()


def test_from_fits_syncs_receipt_extension():
    """RECEIPT extension is in sync with self.receipt right after from_fits,
    without needing a write/read roundtrip (raised in review of #170).

    read() (decorated) and from_fits both append rows in memory; without an
    explicit sync, self.data["RECEIPT"] lags self.receipt until the next write.
    """
    rv2 = RV2()
    rv2.set_header("PRIMARY", {"INSTRUME": "TEST", "DATE-OBS": "2026-01-01T00:00:00"})
    rv2.receipt_add_entry("user_step", "x=1", "PASS")

    with tempfile.TemporaryDirectory() as tmp:
        fn = os.path.join(tmp, "test_SL2_20260101T000000.fits")
        rv2.to_fits(out_filename=fn)
        reloaded = RV2.from_fits(fn)

    receipt_ext = reloaded.data["RECEIPT"]
    assert isinstance(receipt_ext, Table)
    # The extension must reflect every row in the live DataFrame, including the
    # read and from_fits entries added in memory during the read.
    assert len(receipt_ext) == len(reloaded.receipt)
    functions = [str(v) for v in receipt_ext["FUNCTION"]]
    assert "from_fits" in functions
