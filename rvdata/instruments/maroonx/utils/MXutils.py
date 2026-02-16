#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from astropy.coordinates import Angle, SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
from collections import OrderedDict
import os
from astroquery.simbad import Simbad


def pad_array(arr, target_shape, pad_value=np.nan):
    """
    Pad a 1D array with NaNs to match the target length.
    If the input array is shorter than the specified target length,
    NaNs are appended to its end.

    Parameters
    ----------
    arr: np.ndarray
        Input 1D array to pad.
    target_shape: tuple
        Desired shape of the output array.
    pad_value: float
        Value to pad the array with. Default is np.nan.

    Returns
    -------
    padded: np.ndarray
        Padded array of length 'target_shape'.
    """
    arr = np.atleast_2d(arr)
    current_shape = arr.shape
    padded = np.full(target_shape, pad_value, dtype=arr.dtype)
    min_rows = min(current_shape[0], target_shape[0])
    min_cols = min(current_shape[1], target_shape[1])
    padded[:min_rows, :min_cols] = arr[:min_rows, :min_cols]
    return padded


def clean_key(header):
    """
    Standardize a metadata key to ensure compatibility with FITS standards.

    Parameters
    ----------
    header: dict
        Original header dictionary.

    Returns
    -------
    cleanheader: dict
        Cleaned header with standardized keys.
    """
    cleanheader = {}
    for key, value in header.items():
        k = key.replace("MAROONX", "")
        cleanheader[f"HIERARCH {k}" if len(k) > 8 else k] = value
    return cleanheader


def get_catalogInfo(object_name):
    """
    Fetches and returns the Gaia catalog information
    of a given object from Simbad.

    Parameters
    ----------
    object_name: str
        The name of the object to be searched in Simbad.

    Returns:
    --------
    cat: dict
        Contains Gaia DR3/DR2 source identifier, systemic
        radial velocity (km/s), catalog parallax in mas,
        catalog proper motion RA and dec in arcsec/yr.
    """
    cat = {}
    simbad = Simbad()
    simbad.add_votable_fields("ids", "rvz_radvel", "plx_value",
                              "pmdec", "pmra")
    res = simbad.query_object(object_name)
    if res is None:
        raise RuntimeError(f"SIMBAD returned no data for '{object_name}'")

    for row in res["ids"][0].split("|"):
        if "gaia dr3" in row.lower():
            cat["csrc"] = row[:8]
            cat["cid"] = row[5:]
            break
        elif "gaia dr2" in row.lower():
            cat["csrc"] = row[:8]
            cat["cid"] = row[5:]
    cat["crv"] = float(res["rvz_radvel"][0])
    cat["cplx"] = float(res["plx_value"][0])
    cat["cpmr"] = float(res["pmra"][0])/1000
    cat["cpmd"] = float(res["pmdec"][0])/1000

    return cat


def fetch_instEra(jdval):
    """
    Reads the instrument era table and returns the era within which a
    given observation's JD_UTC_START falls.

    Parameters
    ----------
    jdval : float
        JD_UTC_START of the observation.

    Returns
    -------
    era : str
        The corresponding era number, or raises an error if
        era not found.
    """

    era_path = os.path.join(os.getcwd(), 'config', 'MX_inst_era.csv')
    instera = pd.read_csv(era_path)
    current_jd = Time.now().jd
    instera["End_JD"] = instera["End_JD"].fillna(current_jd)

    for _, row in instera.iterrows():
        start_jd = row["Start_JD"]
        end_jd = row["End_JD"]

        if start_jd <= jdval <= end_jd:
            return row["Era"]

    raise ValueError(
        f"JD {jdval} does not fall in any defined instrument era.")


def standardizeMXHeader(prime, file, obstype, channel, datalvl='L2'):
    """
    Convert MAROON-X header keywords into the standardized EPRV header format.

    Parameters
    ----------
    prime: dict
        Original MAROON-X header dictionary.
    file: str
        Filename of the HDF5 file being processed.
    obstype: str
        Observation type ('SCI' for science, 'CAL' for calibration).
    channel: str
        Flag to specify blue or red camera.
    datalvl: str
        The Data product standard level.

    Returns
    -------
    phead: dict
        Dictionary containing standardized header keywords.
    """

    # Build path to header_map.csv relative to this module, not the current
    # working directory, so it works regardless of where the process is started.
    module_dir = os.path.dirname(__file__)
    maroonx_root = os.path.dirname(module_dir)
    hmap_path = os.path.join(maroonx_root, 'config', 'header_map.csv')
    headmap = pd.read_csv(hmap_path, header=0)

    name = file
    target = prime['MAROONX TELESCOPE TARGETNAME']
    ra_deg = prime['MAROONX TELESCOPE TARGETRA']
    dec_deg = prime['MAROONX TELESCOPE TARGETDEC']
    ra_hms = Angle(ra_deg, unit=u.deg).to_string(
        unit=u.hour, sep=':', precision=3, pad=True)
    dec_dms = Angle(dec_deg, unit=u.deg).to_string(
        unit=u.deg, sep=':', precision=3, pad=True, alwayssign=True)
    if channel.upper() == "RED":
        numorder = 28
    else:
        numorder = 34

    catg = get_catalogInfo(target)
    era = fetch_instEra(float(prime['JD_UTC_START']))
    header_conv = {
        'OBSTYPE': obstype,
        'NUMORDER': numorder,
        'FILENAME': name,
        'DATALVL': datalvl,
        'CSRC2': catg['csrc'], 'CSRC3': catg['csrc'], 'CSRC4': catg['csrc'],
        'CID2': catg['cid'], 'CID3': catg['cid'], 'CID4': catg['cid'],
        'CRA2': ra_hms, 'CRA3': ra_hms, 'CRA4': ra_hms,
        'CDEC2': dec_dms, 'CDEC3': dec_dms, 'CDEC4': dec_dms,
        'CEQNX2': 2000, 'CEQNX3': 2000, 'CEQNX4': 2000,
        'CEPCH2': 2016.0, 'CEPCH3': 2016.0, 'CEPCH4': 2016.0,
        'CRV2': catg['crv'], 'CRV3': catg['crv'], 'CRV4': catg['crv'],
        'CPLX2': catg['cplx'], 'CPLX3': catg['cplx'], 'CPLX4': catg['cplx'],
        'CPMR2': catg['cpmr'], 'CPMR3': catg['cpmr'], 'CPMR4': catg['cpmr'],
        'CPMD2': catg['cpmd'], 'CPMD3': catg['cpmd'], 'CPMD4': catg['cpmd'],
        'INSTERA': era,
    }

    phead = OrderedDict()
    ihead = prime
    for i, row in headmap.iterrows():
        skey = row['STANDARD']
        if skey in header_conv:
            phead[skey] = header_conv[skey]
        else:
            mxkey = row['INSTRUMENT']
            if pd.notnull(mxkey):
                mxval = ihead.get(mxkey, None)
            else:
                mxval = row['DEFAULT']
            phead[skey] = mxval if pd.notnull(mxval) else None

    return phead


def compute_bjd_from_header(hdr):
    """
    Calculates the BJD_TDB value.

    Parameters
    ----------
    hdr: dict
        Original header dictionary.

    Returns
    -------
    t_bjd: float
        Calculated BJD_TDB value.
    """
    lon = -155.469047
    lat = 19.823801
    alt = 4213.0
    jd_utc = float(hdr["JD_UTC_FLUXWEIGHTED_FRD"])
    t = Time(jd_utc, format='jd', scale='utc')

    location = EarthLocation.from_geodetic(
        lat=lat*u.deg, lon=lon*u.deg, height=alt*u.m)

    ra = hdr["MAROONX TELESCOPE TARGETRA"]
    dec = hdr["MAROONX TELESCOPE TARGETDEC"]

    target = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")

    ltt_bary = t.light_travel_time(target, location=location)
    t_bjd = t.tdb + ltt_bary

    return t_bjd.value


def create_wavehead(spec, fiber, channel):
    """
    Construct WCS-compliant wavelength header for a given
    fiber for blue and red camera.

    Parameters
    ----------
    spec: pd.DataFrame
        Spectrum data, can be blue or red.
    fiber: int
        Fiber number for which to generate the wavelength metadata.
    channel: str
        Flag to specify blue or red camera.

    Returns
    -------
    wave_meta: dict
        Dictionary of FITS-compliant wavelength metadata for the fiber.
    """

    if channel.upper() == "RED":
        CRPIX = 2018
        PS_2 = PS_4 = '[1, 4036]'
        naxis1 = 4036
        orders = spec.index.levels[1]   # Red orders 67–94
    else:
        CRPIX = 1977
        PS_2 = PS_4 = '[1, 3954]'
        naxis1 = 3954
        orders = spec.index.levels[1]   # Blue orders 91–124

    wave_meta = {}
    wave_meta["XTENSION"] = "IMAGE"
    wave_meta["BITPIX"] = 16
    wave_meta["NAXIS"] = 2
    wave_meta["NAXIS1"] = naxis1
    wave_meta["NAXIS2"] = len(orders)
    wave_meta["WCSAXES"] = len(orders)
    wave_meta["PCOUNT"] = 0
    wave_meta["GCOUNT"] = 0
    wave_meta["EXTNAME"] = f"{channel}_TRACE{fiber}_WAVE"

    i = 1
    for order in sorted(orders):
        wave = spec['wavelengths'][fiber][order]
        if wave is None:
            continue
        CRVAL = float(wave[CRPIX - 1])
        CDELT = float(wave[CRPIX] - wave[CRPIX - 1])
        wave_meta[f"CTYPE{i}"] = "WAVE"
        wave_meta[f"CPDIS{i}"] = "Non-Parametric"
        wave_meta[f"CUNIT{i}"] = "nanometer"
        wave_meta[f"CRPIX{i}"] = CRPIX
        wave_meta[f"CRVAL{i}"] = CRVAL
        wave_meta[f"CDELT{i}"] = CDELT
        wave_meta[f"PS{i}_0"] = "Spline"
        wave_meta[f"PS{i}_1"] = 0
        wave_meta[f"PS{i}_2"] = PS_2
        wave_meta[f"PS{i}_3"] = order
        wave_meta[f"PS{i}_4"] = PS_4
        i += 1

    return wave_meta


def make_orderTable(spec):
    """
    Generate a combined echelle order table from spectra with
    the starting/ending wavelengths for each order.

    Parameters
    ----------
    spec: pd.DataFrame
        Spectrum data, can be blue or red.

    Returns
    -------
    order_table: pd.DataFrame
        A DataFrame containing echelle order number, internal index, and
        wavelength start/end values for each order.
    """
    orders = spec.index.levels[1]

    wave_start = []
    wave_end = []

    for order in orders:
        w = spec['wavelengths'][3][order]  # fiber 3 is reference
        wave_start.append(np.nanmin(w))
        wave_end.append(np.nanmax(w))

    order_table = pd.DataFrame({
        "echelle_order": np.array(orders, dtype=np.float32),
        "order_index":   np.arange(len(orders), dtype=np.float32),
        "wave_start":    np.array(wave_start, dtype=np.float32),
        "wave_end":      np.array(wave_end, dtype=np.float32),
    })

    return order_table
