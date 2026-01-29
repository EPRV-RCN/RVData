from astropy.io import fits
import numpy as np
import os
import pandas as pd


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
        mode_dep_phead[f"TRACE{i_fiber + 1}"] = inst_pri_hdr[f"{fiber}-OBJ"]

        if inst_pri_hdr["OBSTYPE"] == "Cal":
            mode_dep_phead[f"CLSRC{i_fiber + 1}"] = inst_pri_hdr[f"{fiber}-OBJ"]

        if inst_pri_hdr[f"{fiber}-OBJ"] == inst_pri_hdr["QOBJECT"]:
            for pkey, ikey in catalogue_map.items():
                mode_dep_phead[f"{pkey}{i_fiber + 1}"] = inst_pri_hdr[ikey]
            mode_dep_phead[f"CSRC{i_fiber + 1}"] = "GAIADR2"

    # Set up data standard primary header
    hmap_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)), "config/header_map.csv"
    )
    headmap = pd.read_csv(hmap_path, header=0)

    phead = fits.PrimaryHDU().header
    for i, row in headmap.iterrows():
        skey = row["STANDARD"]
        instkey = row["INSTRUMENT"]
        if row["MODE_DEP"] != "Y":
            if pd.notnull(instkey):
                instval = inst_pri_hdr[instkey]
            else:
                instval = row["DEFAULT"]
            if pd.notnull(instval):
                phead[skey] = instval
            else:
                phead[skey] = None
        else:
            if skey in mode_dep_phead.keys():
                phead[skey] = mode_dep_phead[skey]
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
