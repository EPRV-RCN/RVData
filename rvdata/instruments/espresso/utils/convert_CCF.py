from astropy.io import fits
from rvdata.core.models.level4 import RV4
import pandas as pd
import numpy as np
import os


def add_to_ext_descript(RV4, ext_name, description):
    """Add a row to the EXT_DESCRIPT table for the given extension."""
    row = pd.DataFrame({"Name": [ext_name], "Description": [description]})
    RV4.data["EXT_DESCRIPT"] = pd.concat(
        [RV4.data["EXT_DESCRIPT"], row], ignore_index=True
    )


def convert_CCF(RV4: RV4, file_path: str) -> None:
    # Adding the primary header entries
    with fits.open(file_path["ccf_file"]) as hdul:
        RV4.headers["PRIMARY"]["BERV"] = (
            hdul["primary"].header["HIERARCH ESO QC BERV"],
            "[km/s] Barycentric RV",
        )
        RV4.headers["PRIMARY"]["RV"] = (
            hdul["primary"].header["HIERARCH ESO QC CCF RV"],
            "[km/s] Radial velocity value ",
        )
        RV4.headers["PRIMARY"]["RVERR"] = (
            hdul["primary"].header["HIERARCH ESO QC CCF RV ERROR"],
            "[km/s] Error on the radial velocity value ",
        )
        RV4.headers["PRIMARY"]["RVMETHOD"] = ("CCF", "RV derivation method")
        RV4.headers["PRIMARY"]["SYSVEL"] = (0, "Systemic velocity subtracted from RV")
        RV4.headers["PRIMARY"]["BJDTDB"] = (
            hdul["primary"].header["HIERARCH ESO QC BJD"],
            "[JD] Barycentric Julian Day in Terrestrial Dynamical time",
        )
    # Creating the RV extension
    convert_RV(RV4, file_path["ccf_file"], trace_nb=1)
    # Creating the CCF extension
    add_CCFs(RV4, file_path["ccf_file"])
    # Creating the diagnostic extension
    convert_DIAGNOSTICS(RV4, file_path["ccf_file"])
    if os.path.exists(file_path["ccf_file_B"]):
        convert_RV(RV4, file_path["ccf_file_B"], trace_nb=2)
        add_CCFs(RV4, file_path["ccf_file_B"], trace_nb=2)
        convert_DIAGNOSTICS(RV4, file_path["ccf_file_B"], trace_nb=2)

    convert_custom_files(RV4, file_path)


def convert_RV(RV4, file_path, trace_nb=1):
    with fits.open(file_path) as hdul:
        header = hdul["PRIMARY"].header
        header_ = fits.Header()

        header_["RVMETHOD"] = ("CCF", "RV extraction method")
        header_["SKYRMVD"] = (False,)
        header_["TELLRMVD"] = (False,)
        tab = pd.DataFrame(
            {
                "BJD_TDB": [header["HIERARCH ESO QC BJD"]],
                "RV": [header["HIERARCH ESO QC CCF RV"]],
                "RV_error": [header["HIERARCH ESO QC CCF RV ERROR"]],
                "BC_vel": [header["HIERARCH ESO QC BERV"]],
                "wave_start": [np.nan],
                "wave_end": [np.nan],
            }
        )
        if f"RV{trace_nb}" not in RV4.extensions:
            RV4.create_extension(
                ext_name=f"RV{trace_nb}",
                ext_type="BinTableHDU",
                header=header_,
                data=tab,
            )
            add_to_ext_descript(
                RV4, f"RV{trace_nb}", f"Radial velocity table for trace {trace_nb}"
            )
        else:
            RV4.set_data(f"RV{trace_nb}", tab)
            RV4.set_header(f"RV{trace_nb}", header_)


def add_CCFs(RV4, file_path, trace_nb=1):
    with fits.open(file_path) as hdul:
        primary_ = hdul["PRIMARY"].header
        ccf_data = hdul["scidata"].data
        ccf_header = hdul["scidata"].header
        ccf_header["VELSTART"] = (
            primary_["HIERARCH ESO RV START"],
            primary_.comments["HIERARCH ESO RV START"],
        )
        ccf_header["VELSTEP"] = (
            primary_["HIERARCH ESO RV STEP"],
            primary_.comments["HIERARCH ESO RV STEP"],
        )
        ccf_header["CCFMASK"] = (
            primary_["HIERARCH ESO QC CCF MASK"],
            primary_.comments["HIERARCH ESO QC CCF MASK"],
        )
        ccf_header["VELNSTEP"] = (ccf_data.shape[1], "Number of Steps in Velocity Grid")
        if f"CCF{trace_nb}" in RV4.extensions:
            RV4.set_data(f"CCF{trace_nb}", ccf_data)
            RV4.set_header(f"CCF{trace_nb}", ccf_header)
        else:
            RV4.create_extension(
                ext_name=f"CCF{trace_nb}",
                ext_type="ImageHDU",
                data=ccf_data,
                header=ccf_header,
            )
            add_to_ext_descript(
                RV4,
                f"CCF{trace_nb}",
                f"Cross-correlation function for trace {trace_nb}",
            )


def convert_DIAGNOSTICS(RV4, file_path, trace_nb=1):

    with fits.open(file_path) as hdul:
        header = hdul["primary"].header
        diag = []
        for diag_ in ["FWHM", "CONTRAST", "BIS SPAN"]:
            row = {
                "BJD_TDB": header["HIERARCH ESO QC BJD"],
                "metric_name": diag_,
                "value": header[f"HIERARCH ESO QC CCF {diag_}"],
                "error": header[f"HIERARCH ESO QC CCF {diag_} ERROR"],
            }
            diag.append(row)
        diag = pd.DataFrame(diag)
        RV4.create_extension(
            ext_name=f"DIAGNOSTICS{trace_nb}",
            ext_type="BinTableHDU",
            header=None,
            data=diag,
        )
        add_to_ext_descript(
            RV4, f"DIAGNOSTICS{trace_nb}", f"CCF diagnostics for trace {trace_nb}"
        )


def convert_custom_files(RV4, file_path):
    ext_names = {"ccf_sky_file": "SKYSUB", "ccf_tel_corr_file": "TEL_CORR"}
    for file in ["ccf_sky_file", "ccf_tel_corr_file"]:
        if os.path.isfile(file_path[file]):
            with fits.open(file_path[file]) as hdul:
                data = hdul["scidata"].data
                header = hdul["primary"].header
                header_RV = fits.Header()
                header_CCF = fits.Header()
                header_RV["RVMETHOD"] = ("CCF", "RV extraction method")
                header_RV["SKYRMVD"] = ("sky" in file, "Sky removed?")
                header_RV["TELLRMVD"] = ("tel" in file, "Telluric removed?")
                tab = pd.DataFrame(
                    {
                        "BJD_TDB": [header["HIERARCH ESO QC BJD"]],
                        "RV": [header["HIERARCH ESO QC CCF RV"]],
                        "RV_error": [header["HIERARCH ESO QC CCF RV ERROR"]],
                        "BC_vel": [header["HIERARCH ESO QC BERV"]],
                        "wave_start": [np.nan],
                        "wave_end": [np.nan],
                    }
                )
                RV4.create_extension(
                    ext_name=f"RV_{ext_names[file]}",
                    ext_type="BinTableHDU",
                    data=tab,
                    header=header_RV,
                )
                add_to_ext_descript(
                    RV4,
                    f"RV_{ext_names[file]}",
                    f"Radial velocity from {ext_names[file]} CCF",
                )

                header_CCF["VELSTART"] = (
                    header["HIERARCH ESO RV START"],
                    header.comments["HIERARCH ESO RV START"],
                )
                header_CCF["VELSTEP"] = (
                    header["HIERARCH ESO RV STEP"],
                    header.comments["HIERARCH ESO RV STEP"],
                )
                header_CCF["CCFMASK"] = (
                    header["HIERARCH ESO QC CCF MASK"],
                    header.comments["HIERARCH ESO QC CCF MASK"],
                )
                header_CCF["VELNSTEP"] = (
                    data.shape[1],
                    "Number of Steps in Velocity Grid",
                )

                RV4.create_extension(
                    ext_name=f"CCF_{ext_names[file]}",
                    ext_type="ImageHDU",
                    data=data,
                    header=header_CCF,
                )
                add_to_ext_descript(
                    RV4,
                    f"CCF_{ext_names[file]}",
                    f"Cross-correlation function from {ext_names[file]} CCF",
                )

        else:
            print(f"File {file_path[file]} not found. Skipping...")
