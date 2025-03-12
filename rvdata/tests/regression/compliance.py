import importlib

import pandas as pd
from astropy.io import fits


def check_l2_extensions(inpfile):
    l2_csv = importlib.resources.files("rvdata.core.models.config") / "L2-extensions.csv"
    reference_extensions = pd.read_csv(l2_csv)
    extdf = reference_extensions
    hdul = fits.open(inpfile)
    for i, row in extdf.iterrows():
        ext = row["Name"]
        req = row["Required"]
        if req:
            assert ext in hdul, f"Extension {ext} not found in {inpfile}"

    hdul.close()


def check_l2_header(header):
    ref_csv = importlib.resources.files("rvdata.core.models.config") / "L2-PRIMARY-keywords.csv"
    reference_header = pd.read_csv(ref_csv)
    for i, row in reference_header.iterrows():
        key = row["Keyword"]
        req = -row["Optional"]
        if req:
            assert key in header, f"Keyword {key} not found in header"


if __name__ == "__main__":
    check_l2_extensions("rvstandard.fits")
