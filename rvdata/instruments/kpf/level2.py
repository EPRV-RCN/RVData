from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd
import os
import warnings
from collections import OrderedDict

# import base class
from rvdata.core.models.level2 import RV2


# KPF Level2 Reader
class KPFRV2(RV2):
    """
    Data model and reader for RVData Level 2 (RV) data constructed from
    KPF Level 0 and KPF Level 1 pipeline products.

    This class extends the `RV2` base class to handle Keck Planet Finder (KPF)
    data. It reads the relevant science and calibration extensions
    from both a KPF Level 0 and a KPF Level 1 FITS file, organizes them into a
    standardized format, and provides convenient access to flux, wavelength,
    variance, blaze, and metadata for each fiber and chip.

    Parameters
    ----------
    Inherits all parameters from :class:`RV2`.

    Attributes
    ----------
    extensions : dict
        Dictionary of all created extensions (e.g., 'TRACE2_FLUX',
        'TRACE2_WAVE', etc.), mapping extension names to their data arrays.
    headers : dict
        Dictionary of headers for each extension, mapping extension names to
        their FITS headers.
    data : dict
        Dictionary of data arrays for each extension.

    Notes
    -----
    To construct an RVData Level 2 object, both a KPF Level 0 and a KPF Level
    1 FITS file are required. The classmethod `from_fits` should be used to
    instantiate the object from these files. The `_read` method is not intended
    to be called directly by users.

    Example
    -------
    >>> from rvdata.instruments.kpf.level2 import KPFRV2
    >>> obj = KPFRV2.from_fits("kpf_L1.fits", l0file="kpf_L0.fits")
    >>> obj.to_fits("kpf_L2_standard.fits")
    """

    def _read(self, hdul1: fits.HDUList, **kwargs) -> None:
        hdul0 = fits.open(kwargs["l0file"])
        dateobs = hdul1["PRIMARY"].header["DATE-OBS"]

        blazedf = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "config/smooth_lamp_pattern.csv"),
            header=0,
        )
        for i, row in blazedf.iterrows():
            if dateobs >= row["UT_start_date"] and dateobs <= row["UT_end_date"]:
                blazefile = row["CALPATH"]
                blazepath = os.path.join(
                    os.path.dirname(__file__), "reference_fits", blazefile
                )
                blazeHDU = fits.open(blazepath)
                break

        for i in range(1, 4):
            flux_array = None
            wave_array = None
            var_array = None
            blaze_array = None
            out_prefix = f"TRACE{i+1}_"

            for c, chip in enumerate(["GREEN", "RED"]):
                flux_ext = f"{chip}_SCI_FLUX{i}"
                wave_ext = f"{chip}_SCI_WAVE{i}"
                var_ext = f"{chip}_SCI_VAR{i}"
                blaze_ext = f"{chip}_SCI_BLAZE{i}"

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

                if blaze_ext in hdul1:
                    blaze_data = hdul1[blaze_ext].data
                    blaze_meta = OrderedDict(hdul1[blaze_ext].header)
                else:
                    warnings.warn(
                        "Blaze extensions not found in KPF L1 file, using default."
                    )
                    blaze_data = blazeHDU[flux_ext].data
                    blaze_meta = OrderedDict(blazeHDU[flux_ext].header)

                if blaze_array is None:
                    blaze_array = blaze_data
                else:
                    blaze_array = np.concatenate((blaze_array, blaze_data), axis=0)

            self.create_extension(
                out_prefix + "FLUX", "ImageHDU", data=flux_array, header=flux_meta
            )
            self.create_extension(
                out_prefix + "WAVE", "ImageHDU", data=wave_array, header=wave_meta
            )
            self.create_extension(
                out_prefix + "VAR", "ImageHDU", data=var_array, header=var_meta
            )

            # normalize blaze for each order
            for i in range(blaze_array.shape[0]):
                blaze_array[i, :] = blaze_array[i, :] / np.nanpercentile(
                    blaze_array[i, :], 99
                )

            self.create_extension(
                out_prefix + "BLAZE", "ImageHDU", data=blaze_array, header=blaze_meta
            )

        for i, fiber in zip([1, 5], ["CAL", "SKY"]):
            flux_array = None
            wave_array = None
            var_array = None
            blaze_array = None
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

                if blaze_array is None:
                    blaze_array = blazeHDU[flux_ext].data
                    blaze_meta = OrderedDict(blazeHDU[flux_ext].header)
                else:
                    blaze_array = np.concatenate(
                        (blaze_array, blazeHDU[flux_ext].data), axis=0
                    )

            if i == 1:
                self.set_header(out_prefix + "FLUX", flux_meta)
                self.set_data(out_prefix + "FLUX", flux_array)

                self.set_header(out_prefix + "WAVE", wave_meta)
                self.set_data(out_prefix + "WAVE", wave_array)

                self.set_header(out_prefix + "VAR", var_meta)
                self.set_data(out_prefix + "VAR", var_array)

                self.set_header(out_prefix + "BLAZE", flux_meta)
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

        bary = hdul1["BARY_CORR"].data
        bary_kms = bary["BARYVEL"] / 1000.0

        self.set_header("BARYCORR_KMS", OrderedDict(hdul1["BARY_CORR"].header))
        self.set_header("BARYCORR_Z", OrderedDict(hdul1["BARY_CORR"].header))
        self.set_data("BARYCORR_KMS", bary_kms)
        self.set_data("BARYCORR_Z", bary_kms / 3e5)  # aproximate!!!

        self.create_extension(
            "EXPMETER",
            "BinTableHDU",
            data=hdul0["EXPMETER_SCI"].data,
            header=OrderedDict(hdul0["EXPMETER_SCI"].header),
        )
        self.create_extension(
            "TELEMETRY",
            "BinTableHDU",
            data=hdul1["TELEMETRY"].data,
            header=OrderedDict(hdul1["TELEMETRY"].header),
        )

        self.set_header("BJD_TDB", OrderedDict(hdul1["BARY_CORR"].header))
        self.set_data("BJD_TDB", bary["PHOTON_BJD"])

        self.set_header("INSTRUMENT_HEADER", hdul1["PRIMARY"].header)

        self.set_header("DRP_CONFIG", OrderedDict(hdul1["CONFIG"].header))
        self.set_data("DRP_CONFIG", Table(hdul1["CONFIG"].data).to_pandas())

        self.set_header("RECEIPT", OrderedDict(hdul1["RECEIPT"].header))
        self.set_data("RECEIPT", Table(hdul1["RECEIPT"].data).to_pandas())

        wavelengths = self.data["TRACE2_WAVE"]
        order_table_data = pd.DataFrame(
            {
                "echelle_order": 137 - np.arange(wavelengths.shape[0]),
                "order_index": np.arange(wavelengths.shape[0]),
                "wave_start": np.nanmin(wavelengths.data, axis=1),
                "wave_end": np.nanmax(wavelengths.data, axis=1),
            }
        )
        self.set_data("ORDER_TABLE", order_table_data)

        hmap_path = os.path.join(os.path.dirname(__file__), "config/header_map.csv")
        headmap = pd.read_csv(hmap_path, header=0)

        phead = RV2().headers["PRIMARY"]
        ihead = self.headers["INSTRUMENT_HEADER"]
        for i, row in headmap.iterrows():
            skey = row["STANDARD"]
            kpfkey = row["INSTRUMENT"]
            content = phead.get(skey, '')
            if len(content) == 2:
                description = content[1]
            else:
                description = ''
            if pd.notnull(kpfkey) and kpfkey in ihead.keys():
                kpfval = ihead[kpfkey]
            else:
                kpfval = row["DEFAULT"]
            if pd.notnull(kpfval):
                phead[skey] = (kpfval, description)
            else:
                phead[skey] = (None, description)

        phead["ISSOLAR"] = (ihead["OBJECT"].lower() == "socal", "Is this the Sun?")
        self.set_header("PRIMARY", phead)

        # overwrite EXT_DESCRIPT as a DataFrame, dropping the Comments column
        ext_file = os.path.join(
            os.path.dirname(__file__), "config", "L2-extensions.csv"
        )
        ext_descript = pd.read_csv(ext_file, header=0)
        if "Comments" in ext_descript.columns:
            ext_descript = ext_descript.drop(columns=["Comments"])
        self.set_data("EXT_DESCRIPT", ext_descript)
