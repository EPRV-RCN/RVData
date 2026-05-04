"""Regression tests for echelle_order_to_order_index table-type compatibility."""

import pandas as pd
import pytest
from astropy.table import Table

from rvdata.core.tools.stitch_spectrum import echelle_order_to_order_index


@pytest.fixture
def rows():
    return {
        "ECHELLE_ORDER": [92, 93, 94, 95, 96],
        "ORDER_INDEX": [0, 1, 2, 3, 4],
    }


@pytest.fixture(params=["dataframe", "astropy_table"])
def order_table(request, rows):
    if request.param == "dataframe":
        return pd.DataFrame(rows)
    return Table(rows)


def test_returns_matching_order_index(order_table):
    assert echelle_order_to_order_index(94, order_table) == 2


def test_returns_int(order_table):
    assert isinstance(echelle_order_to_order_index(92, order_table), int)


def test_missing_echelle_order_returns_none(order_table):
    assert echelle_order_to_order_index(999, order_table) is None


def test_missing_column_returns_none(rows):
    table = Table({"ECHELLE_ORDER": rows["ECHELLE_ORDER"]})
    assert echelle_order_to_order_index(94, table) is None


def test_astropy_table_without_indices_matches_dataframe(rows):
    df = pd.DataFrame(rows)
    table = Table(rows)
    assert table.indices == []
    for order in rows["ECHELLE_ORDER"]:
        assert echelle_order_to_order_index(
            order, table
        ) == echelle_order_to_order_index(order, df)
