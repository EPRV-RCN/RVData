from scipy.constants import codata


LIGHT_SPEED_MPS = codata.value("speed of light in vacuum")


def corr_waves_RV(waves, rv_mps):
    """Correct wavelengths for an RV shift in m/s."""

    if rv_mps is None:
        return waves
    return waves / (1.0 + float(rv_mps) / LIGHT_SPEED_MPS)
