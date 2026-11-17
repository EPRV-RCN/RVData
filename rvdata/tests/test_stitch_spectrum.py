"""stitch_deblazed_spectrum with the NaN wavelengths that mask_bad_spectrum_data writes."""

import numpy as np
import pytest

from rvdata.core.tools.stitch_spectrum import (
    get_wavelength_grid_with_constant_velocity,
    mask_bad_spectrum_data,
    stitch_deblazed_spectrum,
)

NPIX = 400


@pytest.fixture
def orders():
    """Two overlapping orders, with bad pixels inside and at both ends of each."""
    wav = np.array([np.linspace(5000.0, 5010.0, NPIX), np.linspace(5008.0, 5018.0, NPIX)])
    flx = 1000.0 + 50.0 * np.sin(3.0 * wav)
    flx[0, 100:110] = np.nan
    flx[1, 250:253] = 0.0
    flx[:, :4] = np.nan
    flx[:, -4:] = np.nan
    return wav, flx, flx.copy(), np.ones_like(wav)


def _stitch(wav, flx, var):
    grid = get_wavelength_grid_with_constant_velocity(5001.0, 5017.0, 1000.0)
    return stitch_deblazed_spectrum(grid, wav, flx, var[None, :])


def test_masked_wavelengths_give_the_same_stitch(orders):
    wav, flx, var, blz = orders
    flxm, wavm, _, varm, _ = mask_bad_spectrum_data(flx, wav, blz, var)
    assert np.isnan(wavm).any()  # the masking NaNs the wavelength too

    expected = _stitch(wav, flxm, varm)  # edges intact, bad flux NaN
    got = _stitch(wavm, flxm, varm)  # what convert_level2_to_level3 passes

    for e, g in zip(expected, got):
        np.testing.assert_allclose(g, e, rtol=1e-10)  # NaNs must match too
    assert np.isfinite(expected[1]).mean() > 0.9


def test_fully_masked_order_is_skipped(orders):
    wav, flx, var, blz = orders
    flxm, wavm, _, varm, _ = mask_bad_spectrum_data(flx, wav, blz, var)
    wavm[1] = np.nan
    flxm[1] = np.nan
    varm[1] = np.nan

    st_wav, st_flx, st_var = _stitch(wavm, flxm, varm)
    alone = _stitch(wavm[:1], flxm[:1], varm[:1])

    np.testing.assert_array_equal(st_flx, alone[1])
    np.testing.assert_array_equal(st_var, alone[2])
    assert np.isfinite(st_flx[st_wav < 5009.5]).mean() > 0.9
