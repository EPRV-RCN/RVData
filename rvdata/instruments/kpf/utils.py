from astropy.time import Time

MJD_TO_JD_OFFSET = 2400000.5


def jd_utc(date_obs, native_header):
    """
    Full Julian Date (UTC) of DATE-OBS, as the RVData standard defines JD_UTC.

    Native KPF headers carry no full-JD card, only ``MJD-OBS``. The JD is
    computed from the standard DATE-OBS (mapped from native ``DATE-BEG``) so
    the two cards agree exactly; ``MJD-OBS`` is the fallback when DATE-OBS is
    missing or date-only.

    Parameters
    ----------
    date_obs : str or None
        The standard DATE-OBS value (ISO datetime).
    native_header : Mapping
        The native KPF primary header.

    Returns
    -------
    float or None
        The full JD, or None if neither source is available.
    """
    if date_obs and "T" in str(date_obs):
        return Time(str(date_obs), scale="utc").jd
    mjd = native_header.get("MJD-OBS")
    if mjd is not None:
        return float(mjd) + MJD_TO_JD_OFFSET
    return None
