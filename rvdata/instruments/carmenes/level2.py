"""
CARMENES Level 2 reader.

This module contains the CARMENES-to-RVData Level 2 translator scaffold.  The
base model already creates the required RVData L2 extensions; the methods below
are organized around filling those extensions from native CARMENES products.
"""

import os
import re
import warnings
from collections import OrderedDict

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy import constants
from astropy.time import Time
from astropy.coordinates import Angle
import astropy.units as u
from rvdata.instruments.carmenes.utils import (
    corr_waves_RV,
)

from rvdata.core.models.level2 import RV2


class CARMENESRV2(RV2):
    """
    Read CARMENES data products and convert them into RVData Level 2 format.

    The public entry point for native CARMENES files is
    ``CARMENESRV2.from_fits(..., instrument="CARMENES")``.  Do not call
    ``RV2.from_fits(..., instrument="CARMENES")`` for this translator: the
    base-class dispatcher currently creates a plain ``RV2`` instance, while the
    CARMENES translator relies on subclass-specific helper methods such as
    ``_validate_input`` and ``_populate_primary_header``.

    Expected implementation areas
    -----------------------------
    - ``_trace_spec``: define native flux/wavelength/variance/blaze HDU names
      for the science trace.
    - ``_populate_barycentric_extensions``: fill ``BARYCORR_KMS``,
      ``BARYCORR_Z``, and ``BJD_TDB``.
    - ``_populate_primary_header``: map the CARMENES primary header to the
      standardized RVData primary header.
    - Optional metadata methods for exposure meter, telemetry, DRP config,
      receipt, drift, tellurics, sky models, and extension descriptions.

    Example
    -------
    >>> from rvdata.instruments.carmenes.level2 import CARMENESRV2
    >>> rv2 = CARMENESRV2.from_fits("carmenes_file.fits", instrument="CARMENES")
    >>> rv2.to_fits("carmenes_L2_standard.fits")
    """

    instrument_name = "CARMENES"

    def generate_standard_filename(self) -> str:
        """Generate a CARMENES standard filename including the VIS/NIR channel."""

        if not hasattr(self, "channel") or self.channel not in ("vis", "nir"):
            raise ValueError(
                "CARMENES channel is not set; read a CARMENES file before "
                "generating a standard filename."
            )

        filename = super().generate_standard_filename()
        prefix, suffix = filename.split("_SL", 1)
        return f"{prefix}{self.channel}_SL{suffix}"

    def _read(self, hdul1: fits.HDUList, **kwargs) -> None:
        """
        Populate this RVData Level 2 object from a native CARMENES FITS file.

        Parameters
        ----------
        hdul1 : fits.HDUList
            Open CARMENES FITS HDU list passed in by ``RVDataModel.read``.
        **kwargs
            Reserved for auxiliary files or conversion options.
        """

        self._validate_input(hdul1, **kwargs)
        self._populate_instrument_header(hdul1)
        self._populate_trace_extensions(hdul1, **kwargs)
        self._populate_order_table(hdul1)
        self._populate_barycentric_extensions(hdul1, **kwargs)

        self._populate_optional_extensions(hdul1, **kwargs)

        self._populate_primary_header(hdul1, **kwargs)
        self._populate_extension_descriptions()

    # ------------------------------------------------------------------
    # High-level translation steps

    def _validate_input(self, hdul1: fits.HDUList, **kwargs) -> None:
        """Validate that the input looks like a CARMENES product, and find out, whether it is a VIS or a NIR file."""

        if "PRIMARY" not in hdul1:
            raise ValueError("CARMENES input must contain a PRIMARY HDU.")

        self.channel = str(hdul1["PRIMARY"].header.get("SUBSYS", "")).lower()
        if self.channel not in ("vis", "nir"):
            raise ValueError(
                "CARMENES channel must be either 'vis' or 'nir'; "
                f"got {self.channel!r}.")

        # One should probably provide here the info on, which fiber is being provided (A or B, sci or cal)
        # for now fiber A (sci)
        self.trace_type = 'sci'

    def _populate_instrument_header(self, hdul1: fits.HDUList) -> None:
        """Store the native primary header as ``INSTRUMENT_HEADER``."""

        # self.set_header("INSTRUMENT_HEADER", OrderedDict(hdul1["PRIMARY"].header))
        self.set_header("INSTRUMENT_HEADER", hdul1["PRIMARY"].header.copy())

    def _populate_trace_extensions(self, hdul1: fits.HDUList, **kwargs) -> None:
        """Populate ``TRACE1_FLUX/WAVE/VAR/BLAZE`` image extensions."""

        trace_spec = self._trace_spec(hdul1, **kwargs)
        out_prefix = "TRACE1_"

        flux_data, flux_header = self._read_image_hdu(hdul1, trace_spec["flux"])
        wave_data, wave_header = self._read_image_hdu(hdul1, trace_spec["wave"])
        fp_drift = hdul1["PRIMARY"].header.get("HIERARCH CARACAL SERVAL FP RV")
        wave_data = corr_waves_RV(wave_data, fp_drift)
        var_data, var_header = self._read_image_hdu(hdul1, trace_spec["var"])

        blaze_ext = trace_spec.get("blaze")
        if blaze_ext is None:
            blaze_data = np.ones_like(flux_data, dtype=float)
            blaze_header = fits.Header()
            blaze_header["BLZNORM"] = (True, "Blaze is normalized")
            blaze_header["BLAZESRC"] = (
                "NONE",
                "No native blaze provided; array set to unity",
            )
        else:
            blaze_data, blaze_header = self._read_image_hdu(hdul1, blaze_ext)

        self._set_or_create_image(out_prefix + "FLUX", flux_data, flux_header)
        self._set_or_create_image(out_prefix + "WAVE", wave_data, wave_header)
        self._set_or_create_image(out_prefix + "VAR", var_data, var_header)
        self._set_or_create_image(out_prefix + "BLAZE", blaze_data, blaze_header)

    def _populate_order_table(
        self, hdul1: fits.HDUList, wave_ext: str = "TRACE1_WAVE"
    ) -> None:
        """Build ``ORDER_TABLE`` from a populated wavelength extension."""

        if wave_ext not in self.data or self.data[wave_ext].size == 0:
            raise NotImplementedError(
                "Populate trace wavelength data before building ORDER_TABLE, "
                "or override _populate_order_table for CARMENES."
            )

        wavelengths = np.asarray(self.data[wave_ext])
        order_table = pd.DataFrame(
            {
                "ECHELLE_ORDER": self._echelle_orders(wavelengths),
                "ORDER_INDEX": np.arange(wavelengths.shape[0]),
                "WAVE_START": np.nanmin(wavelengths, axis=1),
                "WAVE_END": np.nanmax(wavelengths, axis=1),
            }
        )
        for column_name, values in self._order_table_extra_columns(
            hdul1, wavelengths
        ).items():
            values = np.asarray(values)
            if values.shape[0] != wavelengths.shape[0]:
                raise ValueError(
                    f"ORDER_TABLE column {column_name!r} has {values.shape[0]} "
                    f"values, expected {wavelengths.shape[0]}."
                )
            order_table[column_name] = values

        self.set_data("ORDER_TABLE", order_table)

    def _order_table_extra_columns(
        self, hdul1: fits.HDUList, wavelengths: np.ndarray
    ) -> dict[str, np.ndarray]:
        """
        Return optional per-order columns to append to ``ORDER_TABLE``.

        CARMENES stores one SNR and one reduced-chi value per order in the
        primary header using ``HIERARCH CARACAL FOX SNR n`` and
        ``HIERARCH CARACAL FOX RCHI n``, where ``n`` is the order index.
        """

        header = hdul1["PRIMARY"].header
        n_orders = wavelengths.shape[0]
        if self.channel == "vis":
            rchi = np.array(
                [
                    header[f"HIERARCH CARACAL FOX RCHI {order_index}"]
                    for order_index in range(n_orders)
                ],
                dtype=float,
            )
            snr = np.array(
                [
                    header[f"HIERARCH CARACAL FOX SNR {order_index}"]
                    for order_index in range(n_orders)
                ],
                dtype=float,
            )

            return {
                "SNR_PER_PIXEL": snr,
                "SQRT_REDUCED_CHI2": rchi,  # maybe np.sqrt(rchi)---ask Mathias!!!
            }
        elif self.channel == "nir":
            rchi_l = np.array(
                [
                    header[f"HIERARCH CARACAL FOX RCHI {2 * order_index}"]
                    for order_index in range(n_orders)
                ],
                dtype=float,
            )
            rchi_r = np.array(
                [
                    header[f"HIERARCH CARACAL FOX RCHI {2 * order_index + 1}"]
                    for order_index in range(n_orders)
                ],
                dtype=float,
            )
            snr_l = np.array(
                [
                    header[f"HIERARCH CARACAL FOX SNR {2 * order_index}"]
                    for order_index in range(n_orders)
                ],
                dtype=float,
            )
            snr_r = np.array(
                [
                    header[f"HIERARCH CARACAL FOX SNR {2 * order_index + 1}"]
                    for order_index in range(n_orders)
                ],
                dtype=float,
            )

            return {
                "SNR_PER_PIXEL_LEFT": snr_l,
                "SNR_PER_PIXEL_RIGHT": snr_r,
                "SQRT_REDUCED_CHI2_LEFT": rchi_l,  # maybe np.sqrt(rchi)---ask Mathias!!!
                "SQRT_REDUCED_CHI2_RIGHT": rchi_r,  # maybe np.sqrt(rchi)---ask Mathias!!!
            }

    def _populate_barycentric_extensions(
        self, hdul1: fits.HDUList, **kwargs
    ) -> None:
        """Populate barycentric correction and BJD extensions."""

        berv_kms = hdul1["PRIMARY"].header["HIERARCH CARACAL BERV"]
        bjd_tdb = hdul1["PRIMARY"].header["HIERARCH CARACAL BJD"] + 2400000.  # or whatever the real key is

        self.set_data("BARYCORR_KMS", np.array([berv_kms], dtype=float))
        self.set_data(
            "BARYCORR_Z",
            np.array([berv_kms / constants.c.to_value("km/s")], dtype=float),
        )
        self.set_data("BJD_TDB", np.array([bjd_tdb], dtype=float))

    def _populate_optional_extensions(self, hdul1: fits.HDUList, **kwargs) -> None:
        """
        Populate optional RVData L2 extensions when CARMENES products provide them.

        Optional files are provided as keyword arguments:

        - ``fiber_b_file``: CARMENES fiber B product.
        - ``raw_image_file``: raw detector image product.
        - ``exposure_meter_file``: exposure meter product.
        """

        fiber_b_file = kwargs.get("fiber_b_file")
        if fiber_b_file is not None:
            with fits.open(fiber_b_file, memmap=False) as hdul_fiber_b:
                self._populate_fiber_b_extensions(
                    hdul1, hdul_fiber_b=hdul_fiber_b, **kwargs
                )

        raw_image_file = kwargs.get("raw_image_file")
        if raw_image_file is not None:
            with fits.open(raw_image_file, memmap=False) as hdul_raw_image:
                self._populate_raw_image_extensions(
                    hdul1, hdul_raw_image=hdul_raw_image, **kwargs
                )

        exposure_meter_file = kwargs.get("exposure_meter_file")
        if exposure_meter_file is not None:
            with fits.open(exposure_meter_file, memmap=False) as hdul_exposure_meter:
                self._populate_exposure_meter_extensions(
                    hdul1, hdul_exposure_meter=hdul_exposure_meter, **kwargs
                )

    def _populate_fiber_b_extensions(
        self, hdul1: fits.HDUList, hdul_fiber_b: fits.HDUList, **kwargs
    ) -> None:
        """Populate extensions derived from an optional fiber B file."""

        trace_spec = self._trace_spec(hdul_fiber_b, **kwargs)
        out_prefix = "fib_B_"

        flux_data, flux_header = self._read_image_hdu(
            hdul_fiber_b, trace_spec["flux"]
        )
        wave_data, wave_header = self._read_image_hdu(
            hdul_fiber_b, trace_spec["wave"]
        )
        fp_drift = hdul1["PRIMARY"].header.get("HIERARCH CARACAL SERVAL FP RV")
        wave_data = corr_waves_RV(wave_data, fp_drift)
        var_data, var_header = self._read_image_hdu(hdul_fiber_b, trace_spec["var"])

        blaze_ext = trace_spec.get("blaze")
        if blaze_ext is None:
            blaze_data = np.ones_like(flux_data, dtype=float)
            blaze_header = fits.Header()
            blaze_header["BLZNORM"] = (True, "Blaze is normalized")
            blaze_header["BLAZESRC"] = (
                "NONE",
                "No native blaze provided; array set to unity",
            )
        else:
            blaze_data, blaze_header = self._read_image_hdu(hdul_fiber_b, blaze_ext)

        self._set_or_create_image(out_prefix + "FLUX", flux_data, flux_header)
        self._set_or_create_image(out_prefix + "WAVE", wave_data, wave_header)
        self._set_or_create_image(out_prefix + "VAR", var_data, var_header)
        self._set_or_create_image(out_prefix + "BLAZE", blaze_data, blaze_header)

    def _populate_raw_image_extensions(
        self, hdul1: fits.HDUList, hdul_raw_image: fits.HDUList, **kwargs
    ) -> None:
        """Populate extensions derived from an optional raw image file."""

        if self.channel == "vis":
            image_hdus = [
                hdu
                for hdu in hdul_raw_image
                if isinstance(hdu, (fits.PrimaryHDU, fits.ImageHDU))
                and hdu.data is not None
            ]
            if len(image_hdus) != 1:
                raise ValueError(
                    "CARMENES VIS raw image file must contain exactly one "
                    f"image extension with data; found {len(image_hdus)}."
                )
            image_hdu = image_hdus[0]
            self._set_or_create_image(
                "RAW_IMAGE", np.asarray(image_hdu.data), image_hdu.header
            )
        elif self.channel == "nir":
            for sca_ext in ("SCA2", "SCA1"):
                image_data, image_header = self._read_image_hdu(
                    hdul_raw_image, sca_ext
                )
                self._set_or_create_image(
                    f"RAW_IMAGE_{sca_ext}", image_data, image_header
                )
        else:
            raise ValueError(
                f"CARMENES channel must be either 'vis' or 'nir'; "
                f"got {self.channel!r}."
            )

    def _populate_exposure_meter_extensions(
        self, hdul1: fits.HDUList, hdul_exposure_meter: fits.HDUList, **kwargs
    ) -> None:
        """Populate extensions derived from an optional exposure meter file."""

        pass

    def _populate_primary_header(self, hdul1: fits.HDUList, **kwargs) -> None:
        """Populate the standardized RVData primary header."""

        hmap_path = os.path.join(
            os.path.dirname(__file__), "config", "header_map_carm.csv"
        )
        headmap = pd.read_csv(hmap_path, header=0)

        phead = RV2().headers["PRIMARY"]
        ihead = self.headers["INSTRUMENT_HEADER"]

        for _, row in headmap.iterrows():
            skey = row["STANDARD"]
            carmenes_key = row["INSTRUMENT"]
            content = phead.get(skey, "")
            description = content[1] if len(content) == 2 else ""
            if "DESCRIPTION" in headmap.columns and pd.notnull(row["DESCRIPTION"]):
                description = row["DESCRIPTION"]

            if pd.notnull(carmenes_key) and carmenes_key in ihead:
                value = ihead[carmenes_key]
            else:
                value = row["DEFAULT"]

            phead[skey] = (value if pd.notnull(value) else None, description)

        # set header keywords here

        utc_from_jd = Time(ihead["HIERARCH CARACAL UTC"], format="jd", scale="utc")
        self._set_primary_value(phead, "DATE", utc_from_jd.isot)
        jd_start = ihead["MJD-OBS"] + 2400000.5
        self._set_primary_value(phead, "JD_UTC", jd_start)

        # INSTERA can be used to track changes to the instrument (maybe in NIR useful?)
        # FULLCOMP could be set to "No", as long not compatible to EPRV standard

        self._set_primary_value(phead, "NUMTRACE", 1)
        self._set_primary_value(phead, "NUMORDER", self.data["TRACE1_WAVE"].shape[0])
        trace1 = str(ihead["HIERARCH CARACAL CATG"]).split(",", 1)[0].strip()
        self._set_primary_value(phead, "TRACE1", trace1)

        ra_deg = ihead["RA"]
        dec_deg = ihead["DEC"]

        cra1 = Angle(ra_deg, unit=u.deg).to_string(
            unit=u.hour,
            sep=":",
            precision=3,
            pad=True,
            )

        cdec1 = Angle(dec_deg, unit=u.deg).to_string(
            unit=u.deg,
            sep=":",
            precision=3,
            alwayssign=True,
            pad=True,
        )

        self._set_primary_value(phead, "CRA1", cra1)
        self._set_primary_value(phead, "CDEC1", cdec1)

        karmn_id = self._plain_value(phead.get("CID1", ""))
        if karmn_id in ("", None, "UNKNOWN"):
            karmn_id = ihead.get("OBJECT", "")

        self.catalog_data = None
        carmencita_data = None
        carmencita_query = kwargs.get("carmencita_query", True)
        if carmencita_query and karmn_id not in ("", None, "UNKNOWN"):
            carmencita_data = self.carmencita_queryID(karmn_id)
            self.catalog_data = carmencita_data

        query_catalog = kwargs.get("query_catalog", False)
        if query_catalog and karmn_id not in ("", None, "UNKNOWN"):
            catalog_data = self.simbad_queryID(karmn_id)
            if catalog_data is not None:
                self.catalog_data = catalog_data

        if self.catalog_data is not None:
            catalog_comments = self.catalog_data.get("comment_suffixes", {})
            self._set_primary_value(
                phead, "CSRC1", self.catalog_data["catalog_source"]
            )
            self._set_primary_value(
                phead, "CID1", self.catalog_data["catalog_identifier"]
            )
            self._set_primary_value(
                phead,
                "CRA1",
                self.catalog_data["ra_sexagesimal"],
                self._comment_with_suffix(
                    phead, "CRA1", catalog_comments.get("CRA1")
                ),
            )
            self._set_primary_value(
                phead,
                "CDEC1",
                self.catalog_data["dec_sexagesimal"],
                self._comment_with_suffix(
                    phead, "CDEC1", catalog_comments.get("CDEC1")
                ),
            )
            self._set_primary_value(phead, "CEQNX1", self.catalog_data["equinox"])
            self._set_primary_value(phead, "CEPCH1", self.catalog_data["epoch"])
            self._set_primary_value(
                phead,
                "CRV1",
                self.catalog_data["systemic_rv_kms"],
                self._comment_with_suffix(
                    phead, "CRV1", catalog_comments.get("CRV1")
                ),
            )
            self._set_primary_value(
                phead,
                "CPLX1",
                self.catalog_data["parallax_mas"],
                self._comment_with_suffix(
                    phead, "CPLX1", catalog_comments.get("CPLX1")
                ),
            )
            self._set_primary_value(
                phead,
                "CPMR1",
                self.catalog_data["pmra_arcsec_per_yr"],
                self._comment_with_suffix(
                    phead, "CPMR1", catalog_comments.get("CPMR1")
                ),
            )
            self._set_primary_value(
                phead,
                "CPMD1",
                self.catalog_data["pmdec_arcsec_per_yr"],
                self._comment_with_suffix(
                    phead, "CPMD1", catalog_comments.get("CPMD1")
                ),
            )
            self._set_primary_value(
                phead,
                "CZ1",
                self.catalog_data["catalog_z"],
                self._comment_with_suffix(phead, "CZ1", catalog_comments.get("CZ1")),
            )
            self._set_primary_value(phead, "CCLRN1", self.catalog_data["color_name"])
            self._set_primary_value(phead, "CCLR1", self.catalog_data["color_value"])

        if carmencita_data is not None:
            aliases = carmencita_data.get("aliases")
            if aliases not in ("", None, "UNKNOWN"):
                self._set_primary_value(phead, "ALIASES", aliases)

        dq_keys = ("DQLVL0", "DQLVL1", "DQLVL2")

        try:
            all_ok = all(int(self._plain_value(phead[key])) == 0 for key in dq_keys)
        except (KeyError, TypeError, ValueError):
            all_ok = False

        self._set_primary_value(phead, "SUMMFLAG", "Pass" if all_ok else "Fail")

        self._set_primary_value(phead, "CHANNEL", self.channel, "CARMENES channel")
        # now populate optional keywords and CARMENES own keywords

        if "LST" in ihead:
            self._set_primary_value(phead, "TLST1", self._hours_to_sexagesimal(ihead["LST"]))

        if "HIERARCH CAHA TEL POS SET RA" in ihead:
            self._set_primary_value(
                phead,
                "TRA1",
                self._ra_deg_to_sexagesimal(ihead["HIERARCH CAHA TEL POS SET RA"]),
            )
        if "HIERARCH CAHA TEL POS SET DEC" in ihead:
            self._set_primary_value(
                phead,
                "TDEC1",
                self._dec_deg_to_sexagesimal(ihead["HIERARCH CAHA TEL POS SET DEC"]),
            )

        # if fiber B is extracted, the extracted flux will be blazed; then we need to change it here

        self.set_header("PRIMARY", phead)

    @staticmethod
    def _plain_value(value):
        """Return the scalar value from a FITS-style (value, comment) tuple."""

        return value[0] if isinstance(value, tuple) else value

    @staticmethod
    def _hours_to_sexagesimal(hours) -> str:
        """Convert decimal hours to a sexagesimal HH:MM:SS.sss string."""

        hours = float(hours) % 24.0
        return Angle(hours, unit=u.hourangle).to_string(
            unit=u.hourangle,
            sep=":",
            precision=3,
            pad=True,
        )

    @staticmethod
    def _set_primary_value(phead: OrderedDict, key: str, value, comment=None) -> None:
        """Set a primary header value while preserving the base comment."""

        content = phead.get(key, "")
        if comment is None:
            comment = content[1] if len(content) == 2 else ""
        phead[key] = (value, comment)

    @staticmethod
    def _comment_with_suffix(phead: OrderedDict, key: str, suffix) -> str | None:
        """Return the current header comment with a catalog reference appended."""

        if suffix is None or pd.isna(suffix):
            return None

        suffix = str(suffix).strip()
        if suffix in ("", "UNKNOWN"):
            return None

        content = phead.get(key, "")
        comment = content[1] if len(content) == 2 else ""
        comment = "" if pd.isna(comment) else str(comment).strip()
        if not comment:
            return suffix
        if suffix in comment.split():
            return comment
        return f"{comment} {suffix}"

    @staticmethod
    def _is_carmenes_id(object_id: str) -> bool:
        """Return True when an identifier looks like a CARMENES J-name."""

        return re.match(r"^J\S+", str(object_id).strip()) is not None

    def carmencita_queryID(self, object_id: str) -> dict | None:
        """
        Resolve an identifier against the local Carmencita CSV subset.

        The local table is keyed primarily by Karmn identifier.  If the input is
        not a Karmn identifier, the object name column is tried as a fallback.
        """

        catalog_path = os.path.join(
            os.path.dirname(__file__), "catalogs", "carmencita_108_subset.csv"
        )
        if not os.path.exists(catalog_path):
            warnings.warn(
                f"Carmencita catalog file not found: {catalog_path}",
                RuntimeWarning,
                stacklevel=2,
            )
            return None

        catalog = pd.read_csv(catalog_path)
        lookup_id = str(object_id).strip()
        if lookup_id.lower().startswith("karmn "):
            lookup_id = lookup_id[6:].strip()

        karmn_values = catalog["Karmn"].astype(str).str.strip()
        match = catalog[karmn_values == lookup_id]
        if match.empty:
            name_values = catalog["Name"].astype(str).str.strip()
            match = catalog[name_values == lookup_id]
        if match.empty:
            return None

        row = match.iloc[0]
        karmn = self._table_value(row, "Karmn")
        name = self._table_value(row, "Name")
        gaia_id = self._table_value(row, "Gaia_id")
        gaia_dr3_id = None
        if gaia_id is not None:
            gaia_dr3_id = f"Gaia DR3 {int(gaia_id)}"

        systemic_rv = self._table_value(row, "Vr_kms-1")
        catalog_z = None
        if systemic_rv is not None:
            catalog_z = float(systemic_rv) / constants.c.to("km/s").value

        bp_mag = self._table_value(row, "BP_mag")
        rp_mag = self._table_value(row, "RP_mag")
        color_value = None
        if bp_mag is not None and rp_mag is not None:
            color_value = float(bp_mag) - float(rp_mag)

        pmra = self._table_value(row, "muRA_masa-1")
        pmdec = self._table_value(row, "muDE_masa-1")
        parallax = self._table_value(row, "pi_mas")
        ref01 = self._table_value(row, "Ref01")
        ref07 = self._table_value(row, "Ref07")
        ref08 = self._table_value(row, "Ref08")
        ref10 = self._table_value(row, "Ref10")

        return {
            "catalog_source": "Carmencita",
            "catalog_identifier": str(karmn) if karmn is not None else lookup_id,
            "aliases": str(name) if name is not None else None,
            "comment_suffixes": {
                "CRA1": ref01,
                "CDEC1": ref01,
                "CPMR1": ref07,
                "CPMD1": ref07,
                "CPLX1": ref08,
                "CRV1": ref10,
                "CZ1": ref10,
            },
            "gaia_dr3_id": gaia_dr3_id,
            "ra_sexagesimal": self._table_value(row, "RA_J2000"),
            "dec_sexagesimal": self._table_value(row, "DE_J2000"),
            "equinox": 2000.0,
            "epoch": 2000.0,
            "systemic_rv_kms": (
                float(systemic_rv) if systemic_rv is not None else None
            ),
            "parallax_mas": float(parallax) if parallax is not None else None,
            "pmra_arcsec_per_yr": (
                float(pmra) / 1000.0 if pmra is not None else None
            ),
            "pmdec_arcsec_per_yr": (
                float(pmdec) / 1000.0 if pmdec is not None else None
            ),
            "catalog_z": catalog_z,
            "color_name": "Gaia BP-RP",
            "color_value": color_value,
        }

    def simbad_queryID(self, object_id: str) -> dict | None:
        """
        Resolve a CARMENES identifier through SIMBAD and Gaia DR3.

        Returns
        -------
        dict or None
            Gaia DR3 identifier, RA/Dec in sexagesimal, equinox, epoch,
            systemic radial velocity in km/s, parallax, proper motion, and
            catalog redshift. ``None`` is returned when the identifier cannot
            be resolved to a Gaia DR3 source.
        """

        try:
            from astroquery.simbad import Simbad

            simbad_object_id = str(object_id).strip()
            if (
                self._is_carmenes_id(simbad_object_id)
                and not simbad_object_id.lower().startswith("karmn ")
            ):
                simbad_object_id = f"Karmn {simbad_object_id}"

            simbad = Simbad()
            simbad.add_votable_fields("ids", "rvz_radvel")
            result = simbad.query_object(simbad_object_id)
            if result is None or len(result) == 0:
                return None

            row = result[0]
            gaia_dr3_id = self._gaia_dr3_id_from_simbad_ids(row["ids"])
            if gaia_dr3_id is None:
                return None

            source_id = gaia_dr3_id.removeprefix("Gaia DR3").strip()
            gaia_row = self._query_gaia_dr3_source(source_id)
            if gaia_row is None:
                return None

            ra_deg = self._table_value(gaia_row, "ra")
            dec_deg = self._table_value(gaia_row, "dec")
            epoch = self._table_value(gaia_row, "ref_epoch")
            systemic_rv = self._table_value(gaia_row, "radial_velocity")
            parallax = self._table_value(gaia_row, "parallax")
            pmra = self._table_value(gaia_row, "pmra")
            pmdec = self._table_value(gaia_row, "pmdec")
            bp_mag = self._table_value(gaia_row, "phot_bp_mean_mag")
            rp_mag = self._table_value(gaia_row, "phot_rp_mean_mag")
            if systemic_rv is None:
                systemic_rv = self._table_value(row, "rvz_radvel")
            catalog_z = None
            if systemic_rv is not None:
                catalog_z = float(systemic_rv) / constants.c.to("km/s").value
            color_value = None
            if bp_mag is not None and rp_mag is not None:
                color_value = float(bp_mag) - float(rp_mag)

            return {
                "catalog_source": "Gaia DR3",
                "catalog_identifier": gaia_dr3_id,
                "gaia_dr3_id": gaia_dr3_id,
                "ra_sexagesimal": self._ra_deg_to_sexagesimal(ra_deg),
                "dec_sexagesimal": self._dec_deg_to_sexagesimal(dec_deg),
                "equinox": 2000.0,
                "epoch": float(epoch) if epoch is not None else None,
                "systemic_rv_kms": (
                    float(systemic_rv) if systemic_rv is not None else None
                ),
                "parallax_mas": float(parallax) if parallax is not None else None,
                "pmra_arcsec_per_yr": (
                    float(pmra) / 1000.0 if pmra is not None else None
                ),
                "pmdec_arcsec_per_yr": (
                    float(pmdec) / 1000.0 if pmdec is not None else None
                ),
                "catalog_z": catalog_z,
                "color_name": "Gaia BP-RP",
                "color_value": color_value,
            }
        except Exception as exc:
            warnings.warn(
                f"SIMBAD/Gaia catalog lookup failed for {object_id!r}: {exc}",
                RuntimeWarning,
                stacklevel=2,
            )
            return None

    @staticmethod
    def _gaia_dr3_id_from_simbad_ids(ids_value) -> str | None:
        """Extract the Gaia DR3 identifier from a SIMBAD ids field."""

        for name in str(ids_value).split("|"):
            name = name.strip()
            if name.lower().startswith("gaia dr3 "):
                return name
        return None

    @staticmethod
    def _query_gaia_dr3_source(source_id: str):
        """Return one Gaia DR3 source row for a numeric Gaia source id."""

        if not re.fullmatch(r"\d+", source_id):
            return None

        from astroquery.gaia import Gaia

        query = f"""
            SELECT TOP 1
                source_id, ra, dec, ref_epoch, radial_velocity,
                parallax, pmra, pmdec, phot_bp_mean_mag, phot_rp_mean_mag
            FROM gaiadr3.gaia_source
            WHERE source_id = {source_id}
        """
        result = Gaia.launch_job(query).get_results()
        if result is None or len(result) == 0:
            return None
        return result[0]

    @staticmethod
    def _table_value(row, key):
        """Return a scalar table-row value, converting masked values to None."""

        value = row[key]
        if np.ma.is_masked(value):
            return None
        if hasattr(value, "item"):
            value = value.item()
        if pd.isna(value):
            return None
        return value

    @staticmethod
    def _ra_deg_to_sexagesimal(ra_deg) -> str:
        """Convert right ascension in degrees to HH:MM:SS.sss."""

        return Angle(float(ra_deg), unit=u.deg).to_string(
            unit=u.hour,
            sep=":",
            precision=3,
            pad=True,
        )

    @staticmethod
    def _dec_deg_to_sexagesimal(dec_deg) -> str:
        """Convert declination in degrees to signed DD:MM:SS.sss."""

        return Angle(float(dec_deg), unit=u.deg).to_string(
            unit=u.deg,
            sep=":",
            precision=3,
            alwayssign=True,
            pad=True,
        )

    def _populate_extension_descriptions(self) -> None:
        """
        Populate ``EXT_DESCRIPT``.

        If ``config/ext_descript_carm.csv`` exists next to this module, it is used.
        Otherwise a compact table is generated from currently present
        extensions so interim products remain inspectable while the translator
        is under development.
        """

        ext_path = os.path.join(os.path.dirname(__file__), "config", "ext_descript_carm.csv")
        if os.path.exists(ext_path):
            ext_descript = pd.read_csv(ext_path)
            if "Comments" in ext_descript.columns:
                ext_descript = ext_descript.drop(columns=["Comments"])
            ext_descript = ext_descript[
                ext_descript["Name"].isin(self.extensions.keys())
            ].copy()
        else:
            ext_descript = pd.DataFrame(
                {
                    "Name": list(self.extensions.keys()),
                    "Description": [
                        "CARMENES translator scaffold placeholder"
                        for _ in self.extensions
                    ],
                }
            )

        self.set_data("EXT_DESCRIPT", ext_descript.reset_index(drop=True))

    # ------------------------------------------------------------------
    # Instrument-specific hooks to fill in

    def _trace_spec(self, hdul1: fits.HDUList, **kwargs) -> dict[str, object]:
        """
        Return native HDU mappings for the science trace.

        The returned dictionary should contain:

        - ``flux``: native flux HDU name.
        - ``wave``: native wavelength HDU name.
        - ``var``: native variance HDU name.
        - ``blaze``: optional native blaze HDU name.  Use ``None`` to create a
          temporary all-ones blaze array.

        Example
        -------
        return {
            "flux": "SCI_FLUX",
            "wave": "SCI_WAVE",
            "var": "SCI_VAR",
            "blaze": "SCI_BLAZE",
        }
        """

        return {
            "flux": "SPEC",
            "wave": "WAVE",
            "var": "SIG",
            "blaze": None,
        }

    def _echelle_orders(self, wavelengths: np.ndarray) -> np.ndarray:
        """Return physical echelle orders for the wavelength array rows."""

        if self.channel == "vis":
            return 118 - np.arange(wavelengths.shape[0])
        elif self.channel == "nir":
            return 63 - np.arange(wavelengths.shape[0])

    # ------------------------------------------------------------------
    # Small helpers

    @staticmethod
    def _require_hdu(hdul1: fits.HDUList, ext_name: str):
        """Return an HDU or raise a helpful error if it is missing."""

        if ext_name not in hdul1:
            raise KeyError(f"Required CARMENES HDU '{ext_name}' not found.")
        return hdul1[ext_name]

    def _read_image_hdu(
        self, hdul1: fits.HDUList, ext_name: str
    ) -> tuple[np.ndarray, fits.Header]:
        """Read an image HDU's data and header."""

        hdu = self._require_hdu(hdul1, ext_name)
        if hdu.data is None:
            raise ValueError(f"CARMENES HDU '{ext_name}' has no data.")
        return np.asarray(hdu.data), hdu.header

    @staticmethod
    def _copy_header(header: fits.Header | OrderedDict | None):
        """Copy FITS headers without losing COMMENT/HISTORY/HIERARCH behavior."""

        if header is None:
            return OrderedDict()
        if isinstance(header, fits.Header):
            return header.copy()
        return OrderedDict(header)

    def _set_or_create_image(
        self, ext_name: str, data: np.ndarray, header: fits.Header | OrderedDict
    ) -> None:
        """Set an existing image extension or create it if this trace is new."""

        ext_name = ext_name.upper()
        header = self._copy_header(header)
        if ext_name in self.extensions:
            self.set_header(ext_name, header)
            self.set_data(ext_name, data)
        else:
            self.create_extension(ext_name, "ImageHDU", data=data, header=header)

    def _set_or_create_table(
        self,
        ext_name: str,
        data: Table | pd.DataFrame,
        header: fits.Header | OrderedDict | None = None,
    ) -> None:
        """Set an existing table extension or create it if absent."""

        ext_name = ext_name.upper()
        if ext_name in self.extensions:
            if header is not None:
                self.set_header(ext_name, self._copy_header(header))
            self.set_data(ext_name, data)
        else:
            self.create_extension(
                ext_name,
                "BinTableHDU",
                data=data,
                header=self._copy_header(header),
            )
