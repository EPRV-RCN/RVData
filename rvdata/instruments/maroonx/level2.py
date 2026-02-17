#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from astropy.io import fits
import os

# import base class
from rvdata.core.models.level2 import RV2
from rvdata.instruments.maroonx.utils import MXutils


# MAROONX Level2 Reader
class MAROONXRV2(RV2):
    """
    Read a MAROONX hd5 data product and convert it to the EPRV standard
    format Python object.

    This class extends the 'RV2' base class to handle the reading of MAROONX
    hdf5 files and converts them into a standardized EPRV format.
    Each extension from the hdf5 file is read, and relevant data,
    including flux, wavelength, variance, and metadata, are stored
    as attributes of the resulting Python object.

    Methods
    -------
    createL2(self, file: str, flat: str) -> None:
        Reads the input HDF5 object file and flat file,
        extracts specific extensions related to the science
        data for different channels and fibers, and stores them
        in a standardized format.

        - The method processes data in different fibers (2, 3, 4, 5, 6)
          from both the BLUE and RED channel.
        - For each channel and fiber, the flux, wavelength, variance,
          and metadata are extracted and stored and separate
          BLUE and RED objects are created.

    Attributes
    ----------
    extensions : dict
        Dictionary of all created extensions (e.g., 'TRACE1_FLUX',
        'TRACE1_WAVE', etc.), mapping extension names to their data arrays.
    headers : dict
        Dictionary of headers for each extension, mapping extension names to
        their FITS headers. The wavelength header metadata is added using MX
        utility function in MXUtils.
    data : dict
        Dictionary of data arrays for each extension.

    Notes
    -----
    - The 'createL2' method processes science and calibration data
    from the BLUE and RED channels, and it extracts and organizes
    data for all the MAROON-X fibers.
    - To construct RVData Level 2 object, a MAROONX hdf5 data product
      and a MAROONX flat files are required.
      The classmethod 'from_fits' cannot be used to read the hd5 files,
    hence, 'createL2' method should be used as per the example below.
    The '_read' method is not intended to be called directly by users.

    Example
    -------
    >>> from core.models.level2 import RV2
    >>> file_MX = 'YYYYMMDDTHHMMSSZ_FFFFF_x_TTTT.hd5'
    >>> flat = 'YYYYMMDDTNN_masterflat__FFFFF_x__blaze.hd5'
    >>> MX_rv2_obj = MAROONXRV2()
    >>> MX_rv2_obj.createL2(file_MX, flat)
    >>> blue_fits, red_fits = MX_rv2_obj.write_camera_fits(file_MX,
                                                        out_dir='.')
    """
    def _read(self, hdul: fits.HDUList) -> None:
        raise RuntimeError(
            "Incorrect usage: MAROONXRV2 does not support RV2.from_fits()"
            "or FITS input.\n"
            "Use createL2('MX.hd5', 'flat.hd5') to build Level 2 products.\n"
            "Example:\n"
            "    rv2 = MAROONXRV2()\n"
            "    rv2.createL2('MXfile.hd5', 'flatfile.hd5')"
        )

    def createL2(self, file: str, flat: str) -> None:
        """
        Reads the HDF5 file and populate two RV2 objects for red and blue:
          - self.blue_product
          - self.red_product
        """
        # Opening HDF5 file and extracting datasets stored in the file
        try:
            store = pd.HDFStore(file, 'r')
        except ImportError as e:
            raise ImportError(
                "The 'tables' (pytables) package is required for reading "
                "MAROON-X HDF5 files. Install it with: "
                "pip install 'rv-data-standard[maroonx]'"
            ) from e
        spec_blue = store['spec_blue']
        header_blue = store['header_blue']
        spec_red = store['spec_red']
        header_red = store['header_red']
        store.close()

        # read blaze file
        bzfile = os.path.join(os.getcwd(), flat)
        store2 = pd.HDFStore(bzfile, 'r+')
        blaze_blue = store2['blaze_blue']
        blaze_red = store2['blaze_red']
        store2.close()

        obstype, flux_key, fiber_range = (
            ('CAL', 'box_extraction', range(1, 6))
            if any(x in file for x in ['EEEE_', 'LLLL_', 'LLLE_'])
            else ('SCI', 'optimal_extraction', range(1, 7))
        )

        # instantiate two RV2 objects for blue and red
        self.blue_product = RV2()
        self.red_product = RV2()

        # Primary header according to EPRV Standard FITS format
        std_primary_blue = MXutils.standardizeMXHeader(
            header_blue, file, obstype, channel="BLUE", datalvl='L2'
        )
        std_primary_red = MXutils.standardizeMXHeader(
            header_red, file, obstype, channel="RED", datalvl='L2'
        )
        self.blue_product.set_header("PRIMARY", std_primary_blue)
        self.red_product.set_header("PRIMARY", std_primary_red)

        # Original Instrument Header
        self.blue_product.set_header(
            "INSTRUMENT_HEADER", MXutils.clean_key(header_blue)
        )
        self.red_product.set_header(
            "INSTRUMENT_HEADER", MXutils.clean_key(header_red)
        )

        # Creating Order tables (separate for blue and red)
        order_table_blue = MXutils.make_orderTable(spec_blue)
        order_table_red = MXutils.make_orderTable(spec_red)
        self.blue_product.set_data("ORDER_TABLE", data=order_table_blue)
        self.red_product.set_data("ORDER_TABLE", data=order_table_red)

        # Extracting flux/wavelength/variance/blaze for blue and red
        blue_data, red_data = {}, {}

        for fiber in fiber_range:
            ref_shape_blue = (34, 3954)
            ref_shape_red = (28, 4036)

            if fiber == 1:
                # Set NaN arrays for fiber 1 as sky fiber is not currently used
                blue_flux = np.full(ref_shape_blue, np.nan)
                red_flux = np.full(ref_shape_red, np.nan)
                blue_wav = np.full(ref_shape_blue, np.nan)
                red_wav = np.full(ref_shape_red, np.nan)
                blue_var = np.full(ref_shape_blue, np.nan)
                red_var = np.full(ref_shape_red, np.nan)
                blue_blz = np.ones(ref_shape_blue)
                red_blz = np.ones(ref_shape_red)
            else:
                blue_flux = (
                    spec_blue['box_extraction'][fiber][:]
                    if fiber == 5
                    else spec_blue[flux_key][fiber][:]
                )
                red_flux = (
                    spec_red['box_extraction'][fiber][:]
                    if fiber == 5
                    else spec_red[flux_key][fiber][:]
                )
                blue_wav = spec_blue['wavelengths'][fiber][:]
                red_wav = spec_red['wavelengths'][fiber][:]
                blue_var = (
                    np.full(ref_shape_blue, np.nan)
                    if obstype == 'CAL'
                    else spec_blue['optimal_var'][fiber][:]
                )
                red_var = (
                    np.full(ref_shape_red, np.nan)
                    if obstype == 'CAL'
                    else spec_red['optimal_var'][fiber][:]
                )
                blue_blz = blaze_blue['blaze'][fiber][:]
                red_blz = blaze_red['blaze'][fiber][:]

            blue_data[fiber] = {
                'FLUX': np.vstack(blue_flux),
                'WAVE': np.vstack(blue_wav),
                'VAR': np.vstack(blue_var),
                'BLAZE': np.vstack(blue_blz),
            }
            red_data[fiber] = {
                'FLUX': np.vstack(red_flux),
                'WAVE': np.vstack(red_wav),
                'VAR': np.vstack(red_var),
                'BLAZE': np.vstack(red_blz),
            }

        for fiber in fiber_range:
            out_prefix = f"TRACE{fiber}_"
            blue = blue_data[fiber]
            red = red_data[fiber]

            for key in ['FLUX', 'WAVE', 'VAR', 'BLAZE']:
                if fiber == 1:
                    self.blue_product.set_data(out_prefix + key, blue[key])
                    self.red_product.set_data(out_prefix + key, red[key])
                else:
                    if key == 'WAVE':
                        # blue wave header:
                        wave_md_b = MXutils.create_wavehead(
                            spec_blue, fiber, channel="BLUE"
                        )
                        wave_hdr_b = fits.Header(wave_md_b)
                        self.blue_product.create_extension(
                            out_prefix + key, "ImageHDU",
                            data=blue[key], header=wave_hdr_b
                        )

                        # red wave header:
                        wave_md_r = MXutils.create_wavehead(
                            spec_red, fiber, channel="RED"
                        )
                        wave_hdr_r = fits.Header(wave_md_r)
                        self.red_product.create_extension(
                            out_prefix + key, "ImageHDU",
                            data=red[key], header=wave_hdr_r
                        )
                    else:
                        self.blue_product.create_extension(
                            out_prefix + key, "ImageHDU", data=blue[key]
                        )
                        self.red_product.create_extension(
                            out_prefix + key, "ImageHDU", data=red[key]
                        )

        berv_kms_b = float(header_blue['BERV_FLUXWEIGHTED_FRD'])/1000.0
        berv_kms_r = float(header_red['BERV_FLUXWEIGHTED_FRD'])/1000.0
        berv_z_b = berv_kms_b / 3e5
        berv_z_r = berv_kms_r / 3e5
        bjd_tdb_b = MXutils.compute_bjd_from_header(header_blue)
        bjd_tdb_r = MXutils.compute_bjd_from_header(header_red)

        drift_b = float(header_blue['Relative_Drift'].split()[0])/1000
        drift_r = float(header_red['Relative_Drift'].split()[0])/1000

        self.blue_product.set_data("BARYCORR_KMS", np.array([berv_kms_b]))
        self.red_product.set_data("BARYCORR_KMS", np.array([berv_kms_r]))
        self.blue_product.set_data("BARYCORR_Z", np.array([berv_z_b]))
        self.red_product.set_data("BARYCORR_Z", np.array([berv_z_r]))
        self.blue_product.set_data("BJD_TDB", np.array([bjd_tdb_b]))
        self.red_product.set_data("BJD_TDB", np.array([bjd_tdb_r]))
        self.blue_product.create_extension(
            "TRACE3_DRIFT", "ImageHDU", data=np.array([drift_b])
        )
        self.red_product.create_extension(
            "TRACE3_DRIFT", "ImageHDU", data=np.array([drift_r])
        )

        # Set up extension Description table
        module_dir = os.path.dirname(__file__)
        exd_path = os.path.join(module_dir, 'config', 'MX_ext_Descript.csv')
        ext_table = pd.read_csv(exd_path)

        if "Required" in ext_table.columns:
            ext_table = ext_table.drop(columns=["Required"])

        created_exts_blue = set(self.blue_product.extensions.keys())
        created_exts_red = set(self.red_product.extensions.keys())

        valid_blue = ext_table[
            ext_table["Name"].isin(created_exts_blue)].copy()
        valid_red = ext_table[ext_table["Name"].isin(created_exts_red)].copy()

        self.blue_product.set_data("EXT_DESCRIPT", data=valid_blue)
        self.red_product.set_data("EXT_DESCRIPT", data=valid_red)

    def write_camera_fits(self, input_filename, out_dir='.'):
        """
        Write the two channel-specific FITS files using a timestamp extracted
        from input_filename.

        Example filenames produced:
            MAROONXBLUE_SL2_YYYYMMDDTHHMMSS.fits
            MAROONXRED_SL2_YYYYMMDDTHHMMSS.fits
        """
        base = os.path.basename(input_filename)
        timestamp = base.split("_")[0].replace("Z", "")
        blue_name = os.path.join(out_dir, f"MAROONXBLUE_SL2_{timestamp}.fits")
        red_name = os.path.join(out_dir, f"MAROONXRED_SL2_{timestamp}.fits")

        self.blue_product.to_fits(blue_name)
        self.red_product.to_fits(red_name)

        return blue_name, red_name
