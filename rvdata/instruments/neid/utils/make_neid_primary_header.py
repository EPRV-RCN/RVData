from astropy.io import fits
import numpy as np
import os
import pandas as pd


# Required catalogue keywords for a science trace whose object has no catalogue entry
# (e.g. the Sun on solar-feed frames, where QOBJECT is blank). Values follow the
# Default column of L2-PRIMARY-keywords.csv; Float keywords get an undefined value
# rather than the string "UNDEFINED".
UNKNOWN_CATALOGUE = {
    "CSRC": "UNKNOWN",
    "CID": "UNKNOWN",
    "CRA": "UNKNOWN",
    "CDEC": "UNKNOWN",
    "CEQNX": None,
    "CEPCH": None,
    "CRV": None,
}


def get_neid_instrument_era(obs_jd):
    """Function to determine which NEID instrument era a given
       observation falls under.

    Parameters
    ----------
    obs_jd : float
        The JD of the observation.

    Returns
    -------
    era : int
        The observation's NEID instrument RV era.
    """

    # Read in the map of dates to NEID instrument RV eras
    eramap = pd.read_csv(
        os.path.join(
            os.path.dirname(os.path.dirname(__file__)), "config/neid_inst_eras.csv"
        )
    )

    # Get the era based on the given observation time
    era_time_diffs = obs_jd - eramap["startdate"].values
    era = eramap["era"].values[
        np.argmin(era_time_diffs[np.where(era_time_diffs >= 0)[0]])
    ]

    return era


def make_base_primary_header(inst_pri_hdr):
    """Function to convert the primary header from a native NEID data
       format FITS file into the standardized data format primary header

    Parameters
    ----------
    inst_pri_hdr : astropy.io.fits.Header
        The NEID data file primary header

    Returns
    -------
    phead : astropy.io.fits.Header
        The standardized data format primary header
    """

    # Header key fix
    for in_hdrkeynam in ("DQLEVEL0", "DQLEVEL1", "DQLEVEL2"):
        if in_hdrkeynam not in inst_pri_hdr.keys():
            continue
        in_hdrkeyval = inst_pri_hdr[in_hdrkeynam]
        if isinstance(in_hdrkeyval, str):
            try:
                out_hdrkeyval = int(in_hdrkeyval.strip())
            except ValueError:
                out_hdrkeyval = None
        elif isinstance(in_hdrkeyval, int):
            if in_hdrkeyval < 0:
                out_hdrkeyval = None
            else:
                out_hdrkeyval = int(in_hdrkeyval)
        else:
            out_hdrkeyval = None
        inst_pri_hdr[in_hdrkeynam] = out_hdrkeyval

    # Set up for obs-mode dependent primary header entries
    mode_dep_phead = {}

    # Create map of keys for target catalogue information
    catalogue_map = {
        "CID": "QOBJECT",
        "CRA": "QRA",
        "CDEC": "QDEC",
        "CEQNX": "QEQNX",
        "CEPCH": "QEPOCH",
        "CPLX": "QPLX",
        "CPMR": "QPMRA",
        "CPMD": "QPMDEC",
        "CRV": "QRV",
        "CZ": "QZ",
    }

    # Check observation mode to set number of traces
    if inst_pri_hdr["OBS-MODE"] == "HR":
        fiber_list = ["SCI", "SKY", "CAL"]
        mode_dep_phead["CLSRC3"] = inst_pri_hdr["CAL-OBJ"]
    elif inst_pri_hdr["OBS-MODE"] == "HE":
        fiber_list = ["SCI", "SKY"]
    mode_dep_phead["NUMTRACE"] = len(fiber_list)

    # Read in the trace information for each fiber
    for i_fiber, fiber in enumerate(fiber_list):
        # TRACE<n> records the *type* of object in the trace (SCI/CAL/SKY), not
        # the object itself -- the fiber role is already exactly that value.
        # The object in each fiber is recorded separately via CLSRC<n>/CID<n>.
        mode_dep_phead[f"TRACE{i_fiber + 1}"] = fiber

        if inst_pri_hdr["OBSTYPE"] == "Cal":
            mode_dep_phead[f"CLSRC{i_fiber + 1}"] = inst_pri_hdr[f"{fiber}-OBJ"]

        if inst_pri_hdr[f"{fiber}-OBJ"] == inst_pri_hdr["QOBJECT"]:
            for pkey, ikey in catalogue_map.items():
                mode_dep_phead[f"{pkey}{i_fiber + 1}"] = inst_pri_hdr[ikey]
            mode_dep_phead[f"CSRC{i_fiber + 1}"] = "GAIADR2"
        elif fiber == "SCI":
            for pkey, value in UNKNOWN_CATALOGUE.items():
                mode_dep_phead[f"{pkey}{i_fiber + 1}"] = value

    # Set up data standard primary header
    hmap_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)), "config/header_map.csv"
    )
    # Only blank cells are missing: the DEFAULT column holds the literal string "None"
    # (CLSRC1/CLSRC2), which pandas would otherwise read as NaN.
    headmap = pd.read_csv(hmap_path, header=0, keep_default_na=False, na_values=[""])

    phead = fits.PrimaryHDU().header
    for i, row in headmap.iterrows():
        skey = row["STANDARD"]
        instkey = row["INSTRUMENT"]
        if skey in mode_dep_phead.keys():
            # Values set for this observation (e.g. CLSRC on calibration frames) win
            # over the static map, including for keys that have a DEFAULT.
            phead[skey] = mode_dep_phead[skey]
        elif row["MODE_DEP"] != "Y":
            if pd.notnull(instkey):
                instval = inst_pri_hdr[instkey]
            else:
                instval = row["DEFAULT"]
            if pd.notnull(instval):
                phead[skey] = instval
            else:
                phead[skey] = None
        else:
            continue

    # Add instrument era
    phead["INSTERA"] = get_neid_instrument_era(phead["JD_UTC"])

    # Check to see if it is a solar observation and change relevant entries
    if inst_pri_hdr["OBJECT"] == "Sun":
        phead["TELESCOP"] = "NEID Solar Feed"
        phead["TELEID1"] = "NEID Solar Feed"
        phead["ISSOLAR"] = True
    else:
        phead["ISSOLAR"] = False

    return phead
