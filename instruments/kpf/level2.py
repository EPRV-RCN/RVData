from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd
import os
from collections import OrderedDict

# import base class
from core.models.level2 import RV2


# KPF Level2 Reader
class KPFRV2(RV2):
    """
    Read a KPF level 1 file and convert it to the EPRV standard format Python object.

    This class extends the `RV2` base class to handle the reading of KPF (Keck Planet Finder)
    Level 1 files and converts them into a standardized EPRV
    format. Each extension from the FITS file is read, and relevant data, including flux,
    wavelength, variance, and metadata, are stored as attributes of the resulting Python object.

    Methods
    -------
    _read(hdul0: fits.HDUList, hdul1: fits.HDUlist) -> None:
        Reads the input FITS HDU list, extracts specific extensions related to the science
        data for different chips and fibers, and stores them in a standardized format.

        - The method processes science data (`SCI_FLUX`, `SCI_WAVE`, `SCI_VAR`) from both
          the GREEN and RED chips and different fibers (`SKY`, `CAL`).
        - For each chip and fiber, the flux, wavelength, variance, and metadata are extracted
          and stored as a `SpectrumCollection` object.
        - Deletes unused extensions such as `RED_TELLURIC`, `GREEN_TELLURIC`, and `TELEMETRY`.

    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions (e.g., `C1_SCI1`, `C1_SKY1`, `C2_CAL1`)
        where the keys are the extension names and the values are `SpectrumCollection` objects
        for each respective dataset.

    header : dict
        A dictionary containing metadata headers from the FITS file, with each extension's
        metadata stored under its respective key.

    Notes
    -----
    - The `_read` method processes science and calibration data from the GREEN and RED chips,
      and it extracts and organizes data for both the SCI, SKY, and CAL fibers.
    - The method converts the flux, wavelength, and variance for each extension into
      `SpectrumCollection` objects.
    - Unused extensions (like `RED_TELLURIC`, `GREEN_TELLURIC`, and `TELEMETRY`) are removed
      from the object.

    Example
    -------
    >>> from core.models.level2 import RV2
    >>> rv2_obj = RV2.from_fits("kpf_level0_file.fits", "kpf_level1_file.fits")
    >>> rv2_obj.to_fits("standard_level2.fits")
    """

    def _read(self, hdul1: fits.HDUList, **kwargs) -> None:
        hdul0 = fits.open(kwargs["l0file"])

        for i in range(1, 4):
            flux_array = None
            wave_array = None
            var_array = None
            out_prefix = f"TRACE{i+1}_"

            for c, chip in enumerate(["GREEN", "RED"]):
                flux_ext = f"{chip}_SCI_FLUX{i}"
                wave_ext = f"{chip}_SCI_WAVE{i}"
                var_ext = f"{chip}_SCI_VAR{i}"

                if flux_array is None:
                    flux_array = hdul1[flux_ext].data
                    flux_meta = OrderedDict(hdul1[flux_ext].header)
                else:
                    flux_array = np.concatenate(
                        (flux_array, hdul1[flux_ext].data), axis=0
                    )

                if wave_array is None:
                    wave_array = hdul1[wave_ext].data
                    wave_meta = OrderedDict(hdul1[wave_ext].header)
                else:
                    wave_array = np.concatenate(
                        (wave_array, hdul1[wave_ext].data), axis=0
                    )

                if var_array is None:
                    var_array = hdul1[var_ext].data
                    var_meta = OrderedDict(hdul1[var_ext].header)
                else:
                    var_array = np.concatenate((var_array, hdul1[var_ext].data), axis=0)

            self.create_extension(
                out_prefix + "FLUX", "ImageHDU", data=flux_array, header=flux_meta
            )
            self.create_extension(
                out_prefix + "WAVE", "ImageHDU", data=wave_array, header=wave_meta
            )
            self.create_extension(
                out_prefix + "VAR", "ImageHDU", data=var_array, header=var_meta
            )
            blaze = flux_array * 0.0 + 1.0
            self.create_extension(
                out_prefix + "BLAZE", "ImageHDU", data=blaze, header=flux_meta
            )

        for i, fiber in zip([1, 5], ["CAL", "SKY"]):
            flux_array = None
            wave_array = None
            var_array = None
            out_prefix = f"TRACE{i}_"

            for c, chip in enumerate(["GREEN", "RED"]):
                flux_ext = f"{chip}_{fiber}_FLUX"
                wave_ext = f"{chip}_{fiber}_WAVE"
                var_ext = f"{chip}_{fiber}_VAR"

                if flux_array is None:
                    flux_array = hdul1[flux_ext].data
                    flux_meta = OrderedDict(hdul1[flux_ext].header)
                else:
                    flux_array = np.concatenate(
                        (flux_array, hdul1[flux_ext].data), axis=0
                    )

                if wave_array is None:
                    wave_array = hdul1[wave_ext].data
                    wave_meta = OrderedDict(hdul1[wave_ext].header)
                else:
                    wave_array = np.concatenate(
                        (wave_array, hdul1[wave_ext].data), axis=0
                    )

                if var_array is None:
                    var_array = hdul1[var_ext].data
                    var_meta = OrderedDict(hdul1[var_ext].header)
                else:
                    var_array = np.concatenate((var_array, hdul1[var_ext].data), axis=0)

            if i == 1:
                self.set_header(out_prefix + "FLUX", flux_meta)
                self.set_data(out_prefix + "FLUX", flux_array)

                self.set_header(out_prefix + "WAVE", wave_meta)
                self.set_data(out_prefix + "WAVE", wave_array)

                self.set_header(out_prefix + "VAR", var_meta)
                self.set_data(out_prefix + "VAR", var_array)

                self.set_header(out_prefix + "BLAZE", flux_meta)
                self.set_data(out_prefix + "BLAZE", flux_array * 0.0 + 1.0)
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
                blaze = flux_array * 0.0 + 1.0
                self.create_extension(
                    out_prefix + "BLAZE", "ImageHDU", data=blaze, header=flux_meta
                )

        bary = hdul1["BARY_CORR"].data
        bary_kms = bary["BARYVEL"] / 1000.0

        self.set_header("DRIFT", wave_meta)
        self.set_data("DRIFT", np.zeros_like(flux_array))

        self.set_header("BARYCORR_KMS", OrderedDict(hdul1["BARY_CORR"].header))
        self.set_header("BARYCORR_Z", OrderedDict(hdul1["BARY_CORR"].header))
        self.set_data("BARYCORR_KMS", bary_kms)
        self.set_data("BARYCORR_Z", bary_kms / 3e5)  # aproximate!!!

        self.create_extension('EXPMETER', "BinTableHDU", data=hdul0["EXPMETER_SCI"].data, header=OrderedDict(hdul0["EXPMETER_SCI"].header))
        self.create_extension('TELEMETRY', "BinTableHDU", data=hdul1["TELEMETRY"].data, header=OrderedDict(hdul1["TELEMETRY"].header))

        self.set_header("BJD_TDB", OrderedDict(hdul1["BARY_CORR"].header))
        self.set_data("BJD_TDB", bary["PHOTON_BJD"])

        self.set_header("INSTRUMENT_HEADER", hdul1["PRIMARY"].header)

        self.set_header("DRP_CONFIG", OrderedDict(hdul1["CONFIG"].header))
        self.set_data("DRP_CONFIG", Table(hdul1["CONFIG"].data).to_pandas())

        self.set_header("RECEIPT", OrderedDict(hdul1["RECEIPT"].header))
        self.set_data("RECEIPT", Table(hdul1["RECEIPT"].data).to_pandas())

        hmap_path = os.path.join(os.path.dirname(__file__), 'config/header_map.csv')
        headmap = pd.read_csv(hmap_path, header=0)

        phead = fits.PrimaryHDU().header
        ihead = self.headers['INSTRUMENT_HEADER']
        for i, row in headmap.iterrows():
            skey = row['STANDARD']
            kpfkey = row['INSTRUMENT']
            if pd.notnull(kpfkey):
                kpfval = ihead[kpfkey]
            else:
                kpfval = row['DEFAULT']
            if pd.notnull(kpfval):
                phead[skey] = kpfval
            else:
                phead[skey] = None

        print(self.extensions)

        self.set_header("PRIMARY", phead)
