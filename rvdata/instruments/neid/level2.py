# Standard library imports
import os

# Third party library imports
from astropy.io import fits
import numpy as np
import pandas as pd

# import base class
from rvdata.core.models.level2 import RV2


# NEID Level2 Reader
class NEIDRV2(RV2):
    """
    Read a NEID Level 2 file and convert it to the EPRV L2 standard format Python object.

    This class extends the `RV2` base class to handle the reading of NEID Level 2 files and
    converts them into a standardized EPRV Level 2 data format. Relevant extensions are taken from
    the NEID FITS file and stored as attributes of the data standard Python object.

    Methods
    -------
    _read(hdul: fits.HDUList) -> None
        Reads the input FITS HDU list, extracts specific extensions related to the science data
        for different chips and fibers, and stores them in a standardized format.

        The method processes data from different fibers depending on the NEID observation mode.
        There are three traces (SCI/SKY/CAL) for HR mode and two traces (SCI/SKY) for HE mode.

    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions (e.g., `INSTRUMENT_HEADER`,
        `TRACE1_FLUX`, `EXPMETER`) where the keys are the extension names and the values are the
        data type (e.g., ImageHDU, BinTableHDU).

    headers : dict
        A dictionary containing metadata headers relevant for each extension. These are largely
        taken from the input NEID FITS file, but some metadata is added in the translator (e.g.,
        the instrument RV era).

    data : dict
        A dictionary containing the data entries for each extension. These are largely directly
        taken from the input NEID FITS file, but some are reorganized (e.g., exposure meter data).

    Notes
    -----
    To construct an RVData Level 2 object, only a NEID Level 2 FITS file is required. The
    classmethod `from_fits` should be used to instantiate the object from these files. The `_read`
    method is not intended to be called directly by users.

    Example
    -------
    >>> from rvdata.instruments.neid.level2 import NEIDRV2
    >>> hdul = fits.open('neidL2_YYYYMMDDTHHMMSS.fits')
    >>> neid_rv2_obj = NEIDRV2.from_fits("neidL2_YYYYMMDDTHHMMSS.fits", instrument="NEID")
    >>> neid_rv2_obj.to_fits("neid_L2_standard.fits")
    """

    def _read(self, hdul: fits.HDUList) -> None:
        """
        Read and process NEID Level 1 FITS file, extracting all relevant extensions.

        Parameters
        ----------
        hdul : fits.HDUList
            The FITS HDU list to be processed.
        """

        # Set up extension description table
        ext_table = {
            "extension_name": [],
            "description": [],
        }

        # Instrument header
        self.set_header("INSTRUMENT_HEADER", hdul["PRIMARY"].header)
        ext_table["extension_name"].append("INSTRUMENT_HEADER")
        ext_table["description"].append("Primary header of native instrument file")

        # Set up for obs-mode dependent primary header entries
        mode_dep_phead = {}
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

        # Order Table
        order_table_data = pd.DataFrame(
            {
                "echelle_order": 173 - np.arange(hdul["SCIWAVE"].data.shape[0]),
                "order_index": np.arange(hdul["SCIWAVE"].data.shape[0]),
                "wave_start": np.nanmin(hdul["SCIWAVE"].data, axis=1),
                "wave_end": np.nanmax(hdul["SCIWAVE"].data, axis=1),
            }
        )
        self.set_data("ORDER_TABLE", order_table_data)
        ext_table["extension_name"].append("ORDER_TABLE")
        ext_table["description"].append("Table of echelle order information")

        # Prepare fiber-related extensions

        # Check observation mode to set fiber list
        if hdul[0].header["OBS-MODE"] == "HR":
            fiber_list = ["SCI", "SKY", "CAL"]
            expmeter_index = 4
            mode_dep_phead["CLSRC3"] = hdul[0].header["CAL-OBJ"]
        elif hdul[0].header["OBS-MODE"] == "HE":
            fiber_list = ["SCI", "SKY"]
            expmeter_index = 3
        mode_dep_phead["NUMTRACE"] = len(fiber_list)

        for i_fiber, fiber in enumerate(fiber_list):
            mode_dep_phead[f"TRACE{i_fiber+1}"] = hdul[0].header[f"{fiber}-OBJ"]

            if hdul[0].header["OBSTYPE"] == "Cal":
                mode_dep_phead[f"CLSRC{i_fiber+1}"] = hdul[0].header[f"{fiber}-OBJ"]

            if hdul[0].header[f"{fiber}-OBJ"] == hdul[0].header["QOBJECT"]:
                for pkey, ikey in catalogue_map.items():
                    mode_dep_phead[f"{pkey}{i_fiber+1}"] = hdul[0].header[ikey]
                mode_dep_phead[f"CSRC{i_fiber+1}"] = "GAIADR2"

            # Set the input extension names for this fiber
            flux_ext = f"{fiber}FLUX"
            wave_ext = f"{fiber}WAVE"
            var_ext = f"{fiber}VAR"
            blaze_ext = f"{fiber}BLAZE"

            # Set the output extension name prefix for this fiber (1-indexed)
            out_prefix = f"TRACE{i_fiber+1}_"

            # Flux
            flux_array = hdul[flux_ext].data
            flux_meta = hdul[flux_ext].header

            # Wavelength
            wave_array = hdul[wave_ext].data
            wave_meta = hdul[wave_ext].header

            # Variance
            var_array = hdul[var_ext].data
            var_meta = hdul[var_ext].header

            # Blaze
            blaze_array = hdul[blaze_ext].data
            blaze_meta = hdul[blaze_ext].header

            # Output extensions into base model. If first fiber, extension already exists in object
            if i_fiber == 0:
                self.set_header(out_prefix + "FLUX", flux_meta)
                self.set_data(out_prefix + "FLUX", flux_array)

                self.set_header(out_prefix + "WAVE", wave_meta)
                self.set_data(out_prefix + "WAVE", wave_array)

                self.set_header(out_prefix + "VAR", var_meta)
                self.set_data(out_prefix + "VAR", var_array)

                self.set_header(out_prefix + "BLAZE", blaze_meta)
                self.set_data(out_prefix + "BLAZE", blaze_array)
            else:
                self.create_extension(
                    out_prefix + "FLUX", "ImageHDU", data=flux_array, header=flux_meta
                )

                self.create_extension(
                    out_prefix + "WAVE", "ImageHDU", data=wave_array, header=wave_meta
                )

                self.create_extension(
                    out_prefix + "VAR", "ImageHDU", data=var_array, header=var_meta
                )

                self.create_extension(
                    out_prefix + "BLAZE",
                    "ImageHDU",
                    data=blaze_array,
                    header=blaze_meta,
                )

            # Extension description table entries
            ext_table["extension_name"].append(out_prefix + "FLUX")
            ext_table["description"].append(
                f"Flux in NEID {hdul[0].header["OBS-MODE"]} {fiber} fiber"
            )

            ext_table["extension_name"].append(out_prefix + "WAVE")
            ext_table["description"].append(
                f"Wavelength solution for NEID {hdul[0].header["OBS-MODE"]} {fiber} fiber"
            )

            ext_table["extension_name"].append(out_prefix + "VAR")
            ext_table["description"].append(
                f"Flux variance in NEID {hdul[0].header["OBS-MODE"]} {fiber} fiber"
            )

            ext_table["extension_name"].append(out_prefix + "BLAZE")
            ext_table["description"].append(
                f"Blaze for NEID {hdul[0].header["OBS-MODE"]} {fiber} fiber"
            )

        # Barycentric correction and timing related extensions

        # Extract barycentric velocities, redshifts, and JDs from NEID primary header
        bary_kms = np.array(
            [hdul[0].header[f"SSBRV{173-order:03d}"] for order in range(122)]
        )
        bary_z = np.array(
            [hdul[0].header[f"SSBZ{173-order:03d}"] for order in range(122)]
        )
        bjd = np.array(
            [hdul[0].header[f"SSBJD{173-order:03d}"] for order in range(122)]
        )

        # Output (these currently do not have headers to inherit from the NEID data format)
        self.set_data("BARYCORR_KMS", bary_kms)
        ext_table["extension_name"].append("BARYCORR_KMS")
        ext_table["description"].append(
            "Barycentric correction velocity per order in km/s"
        )

        self.set_data("BARYCORR_Z", bary_z)
        ext_table["extension_name"].append("BARYCORR_Z")
        ext_table["description"].append(
            "Barycentric correction velocity per order in redshift (z)"
        )

        self.set_data("BJD_TDB", bjd)
        ext_table["extension_name"].append("BJD_TDB")
        ext_table["description"].append(
            "Photon weighted midpoint per order, barycentric dynamical time (JD)"
        )

        # Drift

        # Just set the value of driftrv0 from the header in km/s
        drift_data = np.array([hdul[0].header["driftrv0"] / 1e3])
        drift_meta = fits.Header(
            {"COMMENT": "NEID drift relative to start of observing session"}
        )
        self.create_extension("DRIFT", "ImageHDU", header=drift_meta, data=drift_data)
        ext_table["extension_name"].append("DRIFT")
        ext_table["description"].append("Instrument drift velocity in km/s")

        # Expmeter (316 time stamps, 122 wavelengths)
        expmeter_data = hdul["EXPMETER"].data[expmeter_index].astype(np.float64)

        # The array of exposure meter time stamps and wavelengths
        expmeter_times = hdul["EXPMETER"].data[0, 0].astype(np.float64)
        expmeter_wavelengths = hdul["EXPMETER"].data[1, :, 0].astype(np.float64)

        # Turn the exposure meter information into a dictionary, read in as a DataFrame
        expmeter_extension_data = {"time": expmeter_times}
        for i_wave, col_wavelength in enumerate(expmeter_wavelengths):
            expmeter_extension_data[str(col_wavelength)] = expmeter_data[i_wave]
        expmeter_extension_data = pd.DataFrame(expmeter_extension_data)

        self.create_extension("EXPMETER", "BinTableHDU", data=expmeter_extension_data)
        ext_table["extension_name"].append("EXPMETER")
        ext_table["description"].append("Chromatic exposure meter for science fiber")

        # Telemetry - Nothing for now

        # Telluric model (from NEID L2 extension - combine line and continuum models)
        self.create_extension(
            "TRACE1_TELLURIC",
            "ImageHDU",
            header=hdul["TELLURIC"].header,
            data=hdul["TELLURIC"].data[:, :, 0] * hdul["TELLURIC"].data[:, :, 1],
        )
        ext_table["extension_name"].append("TRACE1_TELLURIC")
        ext_table["description"].append(
            "Telluric model for science fiber (line and continuum absorption combined)"
        )

        # Sky model - Nothing for now

        # Ancillary spectra - Nothing for now

        # Images - Nothing for now

        # Standardized primary header

        hmap_path = os.path.join(os.path.dirname(__file__), "config/header_map.csv")
        headmap = pd.read_csv(hmap_path, header=0)

        phead = fits.PrimaryHDU().header
        ihead = self.headers["INSTRUMENT_HEADER"]
        for i, row in headmap.iterrows():
            skey = row["STANDARD"]
            instkey = row["INSTRUMENT"]
            if row["MODE_DEP"] != "Y":
                if pd.notnull(instkey):
                    instval = ihead[instkey]
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
        eramap = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "config/neid_inst_eras.csv")
        )
        era_time_diffs = phead["JD_UTC"] - eramap["startdate"].values
        era = eramap["era"].values[
            np.argmin(era_time_diffs[np.where(era_time_diffs >= 0)[0]])
        ]
        phead["INSTERA"] = era

        self.set_header("PRIMARY", phead)

        # Set extension description table
        self.set_data("EXT_DESCRIPT", pd.DataFrame(ext_table))
