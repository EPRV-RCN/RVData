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
        # check that all required extensions exist
        hdul = fits.open(inpfile)
        for i, row in extdf.iterrows():
            ext = row["Name"]
            req = row["Required"]
            if req:
                assert ext in hdul, f"Required extension {ext} not found in {inpfile}"

        # check that all extensions from the EXT_DESCRIPT table exist
        ext_table = pd.DataFrame(hdul["EXT_DESCRIPT"].data)
        for i, row in ext_table.iterrows():
            extname = row["Name"]
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
        key = row["Keyword"].split()[0]
        req = row["Required"]
        if req and (key not in header):
            print(key, req)
            assert key in header, f"Keyword {key} not found in header"
        elif req and (key in header):
            value = header[key]
            print(f"{key} = {value} ✓")
        else:
            continue


def check_l3_extensions(inpfile):
    l3_csv = (
        importlib.resources.files("rvdata.core.models.config") / "L3-extensions.csv"
    )
    reference_extensions = pd.read_csv(l3_csv)
    extdf = reference_extensions
    with fits.open(inpfile) as hdul:
        # check that all required extensions exist
        for i, row in extdf.iterrows():
            ext = row["Name"]
            req = row["Required"]
            if req:
                assert ext in hdul, f"Required extension {ext} not found in {inpfile}"

        # check that all extensions from the EXT_DESCRIPT table exist
        ext_table = pd.DataFrame(hdul["EXT_DESCRIPT"].data)
        for i, row in ext_table.iterrows():
            extname = row["Name"]
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

    hdul.close()


def check_l3_header(header):
    ref_csv_l2 = (
        importlib.resources.files("rvdata.core.models.config")
        / "L2-PRIMARY-keywords.csv"
    )
    ref_csv_l3 = (
        importlib.resources.files("rvdata.core.models.config")
        / "L3-PRIMARY-keywords.csv"
    )
    reference_header = pd.concat(
        [pd.read_csv(ref_csv_l2), pd.read_csv(ref_csv_l3)], ignore_index=True
    )
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
    l4_csv = (
        importlib.resources.files("rvdata.core.models.config") / "L4-extensions.csv"
    )
    reference_extensions = pd.read_csv(l4_csv)
    extdf = reference_extensions
    with fits.open(inpfile) as hdul:
        # check that all required extensions exist
        for i, row in extdf.iterrows():
            ext = row["Name"]
            req = row["Required"]
            if req:
                assert ext in hdul, f"Required extension {ext} not found in {inpfile}"

        # check that all extensions from the EXT_DESCRIPT table exist
        ext_table = pd.DataFrame(hdul["EXT_DESCRIPT"].data)
        for i, row in ext_table.iterrows():
            extname = row["Name"]
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

    hdul.close()


def check_l4_header(header):
    ref_csv_l2 = (
        importlib.resources.files("rvdata.core.models.config")
        / "L2-PRIMARY-keywords.csv"
    )
    ref_csv_l4 = (
        importlib.resources.files("rvdata.core.models.config")
        / "L4-PRIMARY-keywords.csv"
    )
    reference_header = pd.concat(
        [pd.read_csv(ref_csv_l2), pd.read_csv(ref_csv_l4)], ignore_index=True
    )
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


def _check_table_columns(
    inpfile, ext_name, csv_filename,
    allowed_extra_prefixes=(), strict=False,
):
    """Check that columns in a BinTableHDU match the standard names from a CSV spec.

    Parameters
    ----------
    inpfile : str
        Path to FITS file.
    ext_name : str
        Name of the extension to check.
    csv_filename : str
        Filename of the CSV spec (in rvdata.core.models.config).
    allowed_extra_prefixes : tuple of str
        Column name prefixes that are allowed beyond the standard set
        (only used when strict=True).
    strict : bool
        If True, reject any non-standard columns (unless they match
        allowed_extra_prefixes). If False, only check that required
        columns are present.
    """
    from astropy.table import Table

    csv_path = (
        importlib.resources.files("rvdata.core.models.config") / csv_filename
    )
    reference_columns = pd.read_csv(csv_path)
    all_standard_columns = reference_columns["Name"].tolist()

    # Determine required columns: if CSV has a Required column use it,
    # otherwise all columns are required
    if "Required" in reference_columns.columns:
        required_columns = reference_columns.loc[
            reference_columns["Required"] == "Yes", "Name"
        ].tolist()
    else:
        required_columns = all_standard_columns

    with fits.open(inpfile) as hdul:
        assert ext_name in hdul, (
            f"{ext_name} extension not found in {inpfile}"
        )

        table = Table(hdul[ext_name].data).to_pandas()
        actual_columns = table.columns.tolist()

        # Check all required columns are present
        for col in required_columns:
            assert col in actual_columns, (
                f"Required {ext_name} column '{col}' not found in {inpfile}. "
                f"Actual columns: {actual_columns}"
            )

        # Optionally check no non-standard column names
        if strict:
            for col in actual_columns:
                if col in all_standard_columns:
                    continue
                if any(col.startswith(prefix) for prefix in allowed_extra_prefixes):
                    continue
                assert False, (
                    f"Non-standard {ext_name} column '{col}' found in {inpfile}. "
                    f"Standard columns are: {all_standard_columns}"
                )


def check_l4_rv_columns(inpfile):
    """Check that RV1 columns use standard names from the CSV spec."""
    _check_table_columns(
        inpfile, "RV1", "L4-RV_TABLE-columns.csv",
        allowed_extra_prefixes=("RV_TRACE",),
        strict=True,
    )


def check_order_table_columns(inpfile):
    """Check that ORDER_TABLE columns use standard names from the CSV spec."""
    _check_table_columns(inpfile, "ORDER_TABLE", "BASE-ORDER_TABLE-columns.csv")


def check_receipt_columns(inpfile):
    """Check that RECEIPT columns use standard names from the CSV spec."""
    _check_table_columns(inpfile, "RECEIPT", "BASE-RECEIPT-columns.csv")


def check_drp_config_columns(inpfile):
    """Check that DRP_CONFIG columns use standard names from the CSV spec."""
    _check_table_columns(inpfile, "DRP_CONFIG", "BASE-DRP_CONFIG-columns.csv")


def check_telemetry_columns(inpfile):
    """Check that TELEMETRY columns use standard names from the CSV spec.

    Only checks if the TELEMETRY extension is present (it is optional).
    """
    with fits.open(inpfile) as hdul:
        if "TELEMETRY" not in hdul:
            return
    _check_table_columns(inpfile, "TELEMETRY", "L2-TELEMETRY-columns.csv")


if __name__ == "__main__":
    check_l2_extensions("rvstandard.fits")
