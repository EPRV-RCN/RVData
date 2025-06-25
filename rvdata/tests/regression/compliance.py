import importlib

import pandas as pd
from astropy.io import fits


def check_l2_extensions(inpfile):
    l2_csv = (
        importlib.resources.files("rvdata.core.models.config") / "L2-extensions.csv"
    )
    reference_extensions = pd.read_csv(l2_csv)
    extdf = reference_extensions
    with fits.open(inpfile) as hdul:
        hdul = fits.open(inpfile)
        for i, row in extdf.iterrows():
            ext = row["Name"]
            req = row["Required"]
            if req:
                assert ext in hdul, f"Required extension {ext} not found in {inpfile}"

        ext_table = pd.DataFrame(hdul["EXT_DESCRIPT"].data)
        for i, row in ext_table.iterrows():
            extname = row["Name"]
            req = row["Required"]
            if req:
                assert (
                    extname in hdul
                ), f"Extension {extname} not found in data but present in EXT_DESCRIPT table."

        # Check every extension in the data (except PRIMARY) has an entry in EXT_DESCRIPT
        ext_names_in_table = ext_table["Name"].tolist()
        for hdu in hdul:
            if hdu.name == "PRIMARY":
                continue
            assert (
                hdu.name in ext_names_in_table
            ), f"Extension {hdu.name} present in data but missing from EXT_DESCRIPT table."


def check_l2_header(header):
    ref_csv = (
        importlib.resources.files("rvdata.core.models.config")
        / "L2-PRIMARY-keywords.csv"
    )
    reference_header = pd.read_csv(ref_csv)
    for i, row in reference_header.iterrows():
        key = row["Keyword"]
        req = row["Required"] == "Y"
        if "#" in key or "..." in key:
            # print(f"Stripping multi keyword: {key}")
            key = key.split("...")[0].strip()
        if req and (key not in header):
            print(key, req)
            assert key in header, f"Keyword {key} not found in header"
        elif req and (key in header):
            value = header[key]
            print(f"{key} = {value} ✓")
        else:
            continue


def check_l4_extensions(inpfile):
    l2_csv = (
        importlib.resources.files("rvdata.core.models.config") / "L4-extensions.csv"
    )
    reference_extensions = pd.read_csv(l2_csv)
    extdf = reference_extensions
    with fits.open(inpfile) as hdul:
        for i, row in extdf.iterrows():
            ext = row["Name"]
            req = row["Required"]
            if req:
                assert ext in hdul, f"Required extension {ext} not found in {inpfile}"

    # commented out because EXT_DESCRIPT is not used in L4 but I suspect we might add it later
    # ext_table = pd.DataFrame(hdul["EXT_DESCRIPT"].data)
    # for i, row in ext_table.iterrows():
    #     extname = row["Name"]
    #     req = row["Required"]
    #     if req:
    #         assert (
    #             extname in hdul
    #         ), f"Extension {extname} not found in data but present in EXT_DESCRIPT table."

    # Check every extension in the data (except PRIMARY) has an entry in EXT_DESCRIPT
    # ext_names_in_table = ext_table["Name"].tolist()
    # for hdu in hdul:
    #     if hdu.name == "PRIMARY":
    #         continue
    #     assert (
    #         hdu.name in ext_names_in_table
    #     ), f"Extension {hdu.name} present in data but missing from EXT_DESCRIPT table."

    hdul.close()


def check_l4_header(header):
    ref_csv = (
        importlib.resources.files("rvdata.core.models.config")
        / "L4-PRIMARY-keywords.csv"
    )
    reference_header = pd.read_csv(ref_csv)
    for i, row in reference_header.iterrows():
        key = row["Keyword"]
        req = row["Required"] == "Y"
        if "#" in key or "..." in key:
            # print(f"Stripping multi keyword: {key}")
            key = key.split("...")[0].strip()
        if req and (key not in header):
            print(key, req)
            assert key in header, f"Keyword {key} not found in header"
        elif req and (key in header):
            value = header[key]
            print(f"{key} = {value} ✓")
        else:
            continue


if __name__ == "__main__":
    check_l2_extensions("rvstandard.fits")
