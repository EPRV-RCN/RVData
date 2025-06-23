import numpy as np
import os
from astropy.io import fits
from astropy.time import Time

# import astropy.units as u
from astropy.constants import c
from collections import OrderedDict
import pandas as pd

# import base class
from rvdata.core.models.level2 import RV2

expres_epochs, epoch_start_isot = np.loadtxt(
    os.path.join(os.path.join(os.path.dirname(__file__), "config/expres_epochs.csv")),
    delimiter=",",
    skiprows=1,
    dtype="str",
).T
epoch_start_mjd = Time(epoch_start_isot).mjd

header_map = pd.read_csv(
    os.path.join(os.path.dirname(__file__), "config/expres_header_map.csv")
).set_index("standard")
header_map.fillna("", inplace=True)

static_headers = {"ORGANIZA": "Yale", "DATALVL": "L2", "NUMTRACE": 1}
obstype_map = {
    "Science": "Sci",
    "Solar": "Sci",
    "Calibration": "Cal",
    "Dark": "Cal",
    "ThAr": "Cal",
    "Quartz": "Cal",
    "LFC": "Cal",
}


# EXPRES Level2 Reader
class EXPRESRV2(RV2):
    """
    Read EXPRES extracted file and convert it to the EPRV standard format Python object.

    This class extends the `RV2` base class to handle the reading of EXPRES
    (EXtreme PREcision Spectrograph) files and converts them into a standardized EPRV
    format. Each extension from the FITS file is read, and relevant data, including flux,
    wavelength, variance, and metadata, are stored as attributes of the resulting Python object.

    Methods
    -------
    _read(hdul: fits.HDUList) -> None:
        Reads the input FITS HDU list, extracts specific extensions related to the science

    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions (e.g., `C1_SCI1`)
        where the keys are the extension names and the values are `SpectrumCollection` objects
        for each respective dataset.

    header : dict
        A dictionary containing metadata headers from the FITS file, with each extension's
        metadata stored under its respective key.

    Notes
    -----
    - The `_read` method processes science and calibration data,
      and it extracts and organizes the science data.
    - The method converts the flux, wavelength, and variance for each extension into
      `SpectrumCollection` objects.
    - Unused extensions are removed from the object.

    Example
    -------
    >>> from astropy.io import fits
    >>> hdul = fits.open('expres_level2_file.fits')
    >>> rv2_obj = EXPRESRV2()
    >>> rv2_obj._read(hdul)
    """

    def _read(self, hdul: fits.HDUList) -> None:
        # Set up extension description table
        data = hdul[1].data.copy()
        ext_table = {
            "extension_name": [],
            "description": [],
        }
        head0 = hdul[0].header
        primary_header = hdul[0].header
        #        fitspec_header = hdul[1].header
        expmeter_header = hdul[2].header
        expmeter_data = hdul[2].data.copy()
        itrace = 1

        self.header_funcs = {
            "OBSTYPE": (lambda hdul: obstype_map[hdul[0].header["OBSTYPE"]]),
            "BINNING": (
                lambda hdul: hdul[0]
                .header["CCDBIN"]
                .replace(" ", "")[1:-1]
                .replace(",", "x")
            ),
            "NUMORDER": (lambda hdul: hdul[1].header["NAXIS2"]),
            "FILENAME": (
                lambda hdul: f"EXPRESL2_{hdul[0].header['MIDPOINT'][2:-1].replace('-', '')}.fits"
            ),
            "JD_UTC": (lambda hdul: Time(hdul[0].header["DATE-SHT"]).jd),
            "INSTERA": (
                lambda hdul: expres_epochs[
                    np.sum(Time(hdul[0].header["MIDPOINT"]).mjd >= epoch_start_mjd) - 1
                ]
            ),
            "INSTFLAG": (
                lambda hdul: (
                    "Pass"
                    if (bool(hdul[0].header["EXPMTR"]) & bool(hdul[0].header["EXPMTR"]))
                    else "Fail"
                )
            ),
            "DRPFLAG": (lambda hdul: drpFlag(hdul)),
            "DRPTAG": (lambda hdul: hdul[1].header["VERSION"]),
            "VERSION": (lambda hdul: hdul[1].header["VERSION"]),
            "EXTRACT": (lambda hdul: hdul[1].header["EXTNAME"]),
        }

        # 0: Primary with just EPRV Standard FITS Headers
        standard_head = OrderedDict()
        for key in header_map.index:
            if key in static_headers.keys():  # Keywords that never change
                standard_head[key] = static_headers[key]
            elif key in self.header_funcs.keys():  # Keywords that require processing
                standard_head[key] = self.header_funcs[key](hdul)
            else:
                _ = header_map.loc[key, "expres"]
                if not _:
                    expres_val = ""

                else:
                    expres_val = head0[_]
                if header_map.loc[key, "required"] == "N" and not expres_val:
                    continue
                standard_head[key] = expres_val
        primary_header = standard_head
        self.set_header("PRIMARY", primary_header)
        ext_table["extension_name"].append("PRIMARY")
        ext_table["description"].append("EPRV Standard Header")

        # 1: Instrument Header
        self.set_header("INSTRUMENT_HEADER", head0)
        ext_table["extension_name"].append("INSTRUMENT_HEADER")
        ext_table["description"].append("Primary header of native instrument file")

        # 2: Receipt
        ext_table["extension_name"].append("RECEIPT")
        ext_table["description"].append("Receipt")

        # 3: DRP_CONFIG
        ext_table["extension_name"].append("DRP_CONFIG")
        ext_table["description"].append("drp configuration information")

        # 4: EXT_DESCRIPT
        ext_table["extension_name"].append("EXT_DESCRIPT")
        ext_table["description"].append("extension config information")

        # 5: ORDER_TABLE
        order_table_data = pd.DataFrame(
            {
                "echelle_order": 160 - np.arange(hdul[1].data["wavelength"].shape[0]),
                "order_index": np.arange(hdul[1].data["wavelength"].shape[0]),
                "wave_start": np.nanmin(hdul[1].data["wavelength"].data, axis=1),
                "wave_end": np.nanmax(hdul[1].data["wavelength"].data, axis=1),
            }
        )
        self.set_data("ORDER_TABLE", order_table_data)
        ext_table["extension_name"].append("ORDER_TABLE")
        ext_table["description"].append("Table of echelle order information")

        # Spectrum data

        itrace = 1
        blaze = data["blaze"]

        # 6: TRACE1_FLUX
        spec = data["spectrum"] * blaze
        self.set_data(f"TRACE{itrace}_FLUX", spec)
        ext_table["extension_name"].append(f"TRACE{itrace}_FLUX")
        ext_table["description"].append("Flux")

        # 7: TRACE1_WAVE
        wave = data["wavelength"]
        self.set_data(f"TRACE{itrace}_WAVE", wave)
        ext_table["extension_name"].append(f"TRACE{itrace}_WAVE")
        ext_table["description"].append("Wavelength solution")

        # 8: TRACE1_VAR
        variance = data["uncertainty"] ** 2.0
        self.set_data(f"TRACE{itrace}_VAR", variance)
        ext_table["extension_name"].append(f"TRACE{itrace}_VAR")
        ext_table["description"].append("Variance")

        # 9: TRACE1_BLAZE
        self.set_data(f"TRACE{itrace}_BLAZE", blaze)
        ext_table["extension_name"].append(f"TRACE{itrace}_BLAZE")
        ext_table["description"].append("Blaze function")

        # # 10+11: BARYCORR_KMS + BARYCORR_Z
        bary_arr = data["bary_wavelength"]  # data['bary_wavelength']
        berv_kms = ((1 - bary_arr / wave) * c.to("km/s")).value
        berv_z = 1 - bary_arr / wave
        self.set_data("BARYCORR_KMS", berv_kms)
        ext_table["extension_name"].append("BARYCORR_KMS")
        ext_table["description"].append(
            "Barycentric correction velocity per order in km/s"
        )

        self.set_data("BARYCORR_Z", berv_z)
        ext_table["extension_name"].append("BARYCORR_Z")
        ext_table["description"].append(
            "Barycentric correction velocity per order in redshift (z)"
        )

        # # 12: BJD_TDB - need to convert this to by pixel
        self.set_data(
            "BJD_TDB", np.array([hdul[2].header["wtd_mdpt"]]).astype(np.float64)
        )
        ext_table["extension_name"].append("BJD_TDB")
        ext_table["description"].append(
            "Photon weighted midpoint, barycentric dynamical time (JD)"
        )

        # # 12: Exposure Meter
        expmeter_header = hdul[2].header
        expmeter_data = hdul[2].data.copy()
        expmeter_array = np.array([row[0] for row in expmeter_data]).T
        expmeter_times = hdul[2].data["midpoints"].astype(np.float64)
        expmeter_wavelengths = hdul[2].data["wavelengths"][0].astype(np.float64)

        expmeter_extension_data = {"time": expmeter_times}
        for i_wave, col_wavelength in enumerate(expmeter_wavelengths):
            expmeter_extension_data[str(col_wavelength)] = expmeter_array[i_wave]

        self.create_extension(
            "EXPMETER",
            "BinTableHDU",
            header=expmeter_header,
            data=expmeter_extension_data,
        )
        ext_table["extension_name"].append("EXPMETER")
        ext_table["description"].append("Chromatic exposure meter")

        # # 13: Telluric Model
        telluric = data["tellurics"]
        self.create_extension(f"TRACE{itrace}_TELLURIC", "ImageHDU", data=telluric)
        ext_table["extension_name"].append("TRACE1_TELLURIC")
        ext_table["description"].append("Telluric line and continuum absorption model")

        # Set extension description table
        self.set_data("EXT_DESCRIPT", pd.DataFrame(ext_table))

    # =============================================================================
    # Methods for standardizing header keywords

    def standardizeExpresHeader(self, hdul):
        head0 = hdul[0].header
        standard_head = OrderedDict()
        for key in header_map.index:
            if key in static_headers.keys():  # Keywords that never change
                standard_head[key] = static_headers[key]
            elif key in self.header_funcs.keys():  # Keywords that require processing
                standard_head[key] = self.header_funcs[key](hdul)
            else:
                _ = header_map.loc[key, "expres"]
                if not _:
                    expres_val = ""

                else:
                    expres_val = head0[_]
                if header_map.loc[key, "required"] == "N" and not expres_val:
                    continue
                standard_head[key] = expres_val
        return standard_head


def drpFlag(hdul):
    extensions_to_check = [
        "spectrum",
        "blaze",
        "wavelength",
        "bary_wavelength",
        "excalibur",
        "bary_excalibur",  # We don't actually need this here
        "continuum",
        "tellurics",
    ]
    extension_list = hdul[1].data.dtype.names
    percent_there = np.sum(
        [extn in extension_list for extn in extensions_to_check]
    ) / len(extensions_to_check)
    if percent_there == 1:
        return "Pass"
    return "Fail" if percent_there == 0 else "Warn"
