#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import Angle
import astropy.units as u
from collections import OrderedDict
import os

# import base class
from core.models.level2 import RV2

# MAROONX Level2 Reader
class MAROONXRV2(RV2):
    """
    Read a MAROONX level 1 file and convert it to the EPRV standard format Python object.

    This class extends the `RV2` base class to handle the reading of MAROONX 
    Level 1 files and converts them into a standardized EPRV
    format. Each extension from the FITS file is read, and relevant data, including flux,
    wavelength, variance, and metadata, are stored as attributes of the resulting Python object.

    Methods
    -------
    _read(self, file: str) -> None:
        Reads the input HDF5 file, extracts specific extensions related to the science
        data for different channels and fibers, and stores them in a standardized format.

        - The method processes data in different fibers (2, 3, 4, 5, 6) from both
          the BLUE and RED channel.
        - For each channel and fiber, the flux, wavelength, variance, and metadata are extracted
          and stored as a `SpectrumCollection` object.

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
    - The `_read` method processes science and calibration data from the BLUE and RED channels,
      and it extracts and organizes data for all the MAROON-X fibers.
    - The method converts the flux, wavelength, and variance for each extension into
      `SpectrumCollection` objects.
    - Unused extensions are removed from the object.

    Example
    -------
    >>> from core.models.level2 import RV2
    >>> file_MAROONX = 'YYYYMMDDTHHMMSSZ_FFFFF_x_TTTT.hd5'
    >>> rv2_obj = MAROONXRV2()
    >>> rv2_obj._read(file_MAROONX)
    """

    @staticmethod
    def pad_array(arr, target_shape, pad_value=np.nan):
        """
        Pad a 1D array with NaNs to match the target length.
        If the input array is shorter than the specified target length, NaNs are appended to its end. 
        """
        arr = np.atleast_2d(arr)
        current_shape = arr.shape
    
        padded_arr = np.full(target_shape, pad_value, dtype=np.float64)
    
        min_rows = min(current_shape[0], target_shape[0])
        min_cols = min(current_shape[1], target_shape[1])
        padded_arr[:min_rows, :min_cols] = arr[:min_rows, :min_cols]
    
        return padded_arr
        
    @staticmethod
    def clean_key(header):
        """
        Standardize a metadata key to ensure compatibility with FITS standards.
        """
        cleanheader = {}
        for key, value in header.items():
            key = key.replace("MAROONX", "")
            cleanheader[f"HIERARCH {key}" if len(key) > 8 else key] = value
        return cleanheader
        
    def standardizeMXHeader(self, prime, file, obstype):
        """
        Convert MAROON-X header keywords into the standardized EPRV header format.
        """
        hmap_path = os.path.join(os.getcwd(), 'config', 'header_map.csv')
        headmap = pd.read_csv(hmap_path, header=0)

        name = file
        ra_deg = prime['MAROONX TELESCOPE TARGETRA']  # Right Ascension in degrees
        dec_deg = prime['MAROONX TELESCOPE TARGETDEC']  # Declination in degrees
        ra_hms = Angle(ra_deg, unit=u.deg).to_string(unit=u.hour, sep=':', precision=3, pad=True)
        dec_dms = Angle(dec_deg, unit=u.deg).to_string(unit=u.degree, sep=':', precision=3, pad=True, alwayssign=True)
    
        header_conv = {
        'OBSTYPE': obstype,
        'FILENAME': name,
        'CRA2': ra_hms,
        'CRA3': ra_hms,
        'CRA4': ra_hms,
        'CDEC2': dec_dms,
        'CDEC3': dec_dms,
        'CDEC4': dec_dms,
        }
    
        phead = OrderedDict()
        ihead = prime
        for i, row in headmap.iterrows():
            skey = row['STANDARD']
            if skey in header_conv: # Converted
                phead[skey] = header_conv[skey]
            else:
                mxkey = row['INSTRUMENT']
                if pd.notnull(mxkey):
                    mxval = ihead[mxkey]
                else:
                    mxval = row['DEFAULT']
                if pd.notnull(mxval):
                    phead[skey] = mxval
                else:
                    phead[skey] = None

        return phead

    def create_wavehead(self, spec_blue, spec_red, fiber):
        """
        Construct WCS-compliant wavelength header for a given fiber for both blue and red channel combined.
        """
        wave_meta = {}
        orders_blue = spec_blue.index.levels[1]  # Blue orders: 91-124
        orders_red = spec_red.index.levels[1]    # Red orders: 67-94
    
        wave_meta["XTENSION"] = "IMAGE"
        wave_meta["BITPIX"] = 16
        wave_meta["NAXIS"] = 2
        wave_meta["NAXIS1"] = 4036
        wave_meta["NAXIS2"] = 62
        wave_meta["WCSAXES"] = 62
        wave_meta["PCOUNT"] = 0
        wave_meta["GCOUNT"] = 0
        wave_meta["EXTNAME"] = f"TRACE{fiber}_WAVE"
    
        i = 1  
    
        for order in sorted(set(orders_blue).union(orders_red)):
            if order in orders_red:
                CRPIX = 2018
                PS_2 = PS_4 = '[1, 4036]'
                wave = spec_red['wavelengths'][fiber][order]
                if wave is not None:
                    CRVAL = float(wave[CRPIX - 1])
                    CDELT = float(wave[CRPIX] - wave[CRPIX - 1])
    
                    wave_meta[f"CTYPE{i}"] = "WAVE"
                    wave_meta[f"CPDIS{i}"] = "Non-Parametric"
                    wave_meta[f"CUNIT{i}"] = "nanometer"
                    wave_meta[f"CRPIX{i}"] = CRPIX
                    wave_meta[f"CRVAL{i}"] = CRVAL
                    wave_meta[f"CDELT{i}"] = CDELT
                    wave_meta[f"PS{i}_0"] = "Spline"
                    wave_meta[f"PS{i}_1"] = 0
                    wave_meta[f"PS{i}_2"] = PS_2
                    wave_meta[f"PS{i}_3"] = order
                    wave_meta[f"PS{i}_4"] = PS_4
                    i += 1  
    
            if order in orders_blue:
                CRPIX = 1977
                PS_2 = PS_4 = '[1, 3954]'
                wave = spec_blue['wavelengths'][fiber][order]
                if wave is not None:
                    CRVAL = float(wave[CRPIX - 1])
                    CDELT = float(wave[CRPIX] - wave[CRPIX - 1])
    
                    wave_meta[f"CTYPE{i}"] = "WAVE"
                    wave_meta[f"CPDIS{i}"] = "Non-Parametric"
                    wave_meta[f"CUNIT{i}"] = "nanometer"
                    wave_meta[f"CRPIX{i}"] = CRPIX
                    wave_meta[f"CRVAL{i}"] = CRVAL
                    wave_meta[f"CDELT{i}"] = CDELT
                    wave_meta[f"PS{i}_0"] = "Spline"
                    wave_meta[f"PS{i}_1"] = 0
                    wave_meta[f"PS{i}_2"] = PS_2
                    wave_meta[f"PS{i}_3"] = order
                    wave_meta[f"PS{i}_4"] = PS_4
                    i += 1  
        return wave_meta

    def make_orderTable(self, spec_blue, spec_red):
        """
        Generate a combined echelle order table from blue and red channel spectra with 
        the starting/ending wavelengths for each order.
        """
        ordB = spec_blue.index.levels[1]
        ordR = spec_red.index.levels[1]
        e_order = np.hstack([ordR, ordB])
        i_order = np.arange(len(e_order), dtype=int)
        wave_start = []
        wave_end = []
        
        for orders, wavelengths in [(ordR, spec_red['wavelengths'][3]), (ordB, spec_blue['wavelengths'][3])]:
            for o in orders:
                wave = wavelengths[o]
                wave_start.append(np.nanmin(wave))
                wave_end.append(np.nanmax(wave))
                
        order_table = pd.DataFrame({
            "echelle_order": np.array(e_order, dtype=np.float32),
            "order_index": np.array(i_order, dtype=np.float32),
            "wave_start": np.array(wave_start, dtype=np.float32),
            "wave_end": np.array(wave_end, dtype=np.float32),
        })
        return order_table
    
    def _read(self, file: str) -> None:
        """
        Read data from an HDF5 file instead of a FITS file.
        """
        # Opening HDF5 file
        store = pd.HDFStore(file, 'r')
    
        # Extracting datasets stored in the hdf5 file
        self.spec_blue = store['spec_blue']
        self.header_blue = store['header_blue']
        self.spec_red = store['spec_red']
        self.header_red = store['header_red']
        store.close() 

        # Reading blaze file
        bzfile = os.path.join(os.getcwd(), '20231012T05_masterflat_backgroundsubtracted_FFFFF_x_blaze.hd5')
        store2 = pd.HDFStore(bzfile, 'r+')
        blaze_blue   = store2['blaze_blue']
        blaze_red    = store2['blaze_red']   
        store2.close()

        obstype, flux_key, fiber_range = ('CAL', 'box_extraction', range(1, 6)) if any(x in file for x in ['EEEE_', 'LLLL_', 'LLLE_']) else ('SCI', 'optimal_extraction', range(1, 7))    

         # Set up extension description table
        ext_table = {
            "extension_name": [],
            "description": [],
        }

        # Primary according to EPRV Standard FITS format
        self.set_header("PRIMARY", self.standardizeMXHeader(self.header_blue, file, obstype))
        ext_table['extension_name'].append("PRIMARY")
        ext_table['description'].append("EPRV Standard Header")
        
        # Original Instrument Header
        self.set_header("INSTRUMENT_HEADER", self.clean_key(self.header_blue))
        ext_table["extension_name"].append("INSTRUMENT_HEADER")
        ext_table["description"].append("Header inherited from native instrument file for Blue camera")
        self.create_extension("INSTRUMENT_HEADER_RED", "ImageHDU", header=self.clean_key(self.header_red))
        ext_table["extension_name"].append("INSTRUMENT_HEADER_RED")
        ext_table["description"].append("Header inherited from native instrument file for Red camera")       

        #Create order table
        order_table_data = self.make_orderTable(self.spec_blue, self.spec_red)
        self.create_extension("ORDER_TABLE", "BinTableHDU", data=order_table_data)
        ext_table["extension_name"].append("ORDER_TABLE")
        ext_table["description"].append("Table of echelle order information")

        # Extracting flux/wavelength/variance/blaze
        blue_data, red_data = {}, {}
        for fiber in fiber_range:
            ref_shape_blue, ref_shape_red = (34, 3954), (28, 4036)
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
                blue_flux = self.spec_blue['box_extraction'][fiber][:] if fiber == 5 else self.spec_blue[flux_key][fiber][:]
                red_flux = self.spec_red['box_extraction'][fiber][:] if fiber == 5 else self.spec_red[flux_key][fiber][:]

                blue_wav = self.spec_blue['wavelengths'][fiber][:]
                red_wav = self.spec_red['wavelengths'][fiber][:]
        
                blue_var = np.full(ref_shape_blue, np.nan) if obstype == 'CAL' else self.spec_blue['optimal_var'][fiber][:]
                red_var = np.full(ref_shape_red, np.nan) if obstype == 'CAL' else self.spec_red['optimal_var'][fiber][:]
                
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
       
        # Padding and Combining  
        for fiber in fiber_range:
            out_prefix = f"TRACE{fiber}_"
            blue, red = blue_data[fiber], red_data[fiber]
            max_cols = max(blue['FLUX'].shape[1], red['FLUX'].shape[1])
            
            for key in ['FLUX', 'WAVE', 'VAR', 'BLAZE']:
                blue[key] = self.pad_array(blue[key], (blue[key].shape[0], max_cols))
                red[key] = self.pad_array(red[key], (red[key].shape[0], max_cols))
        
            combined = {key: np.concatenate((blue[key], red[key]), axis=0) for key in blue.keys()}
                
            for key, dtype in zip(['FLUX', 'WAVE', 'VAR', 'BLAZE'], [np.float32, np.float64, np.float32, np.float32]):
                if fiber == 1:
                    self.set_data(out_prefix + key, combined[key])
                else:
                    if key == 'WAVE':
                        wave_md = self.create_wavehead(self.spec_blue, self.spec_red, fiber)
                        wave_meta = fits.Header(wave_md)
                        self.create_extension(out_prefix + key, "ImageHDU", data=combined[key], header=wave_meta)
                    else:
                        self.create_extension(out_prefix + key, "ImageHDU", data=combined[key])
                ext_table["extension_name"].append(out_prefix + key)
                ext_table["description"].append(f"{key} in fiber {fiber}")


        berv_kms = float(self.header_blue['BERV_FLUXWEIGHTED_FRD'])/1000.0
        ext_table["extension_name"].append("BARYCORR_KMS")
        ext_table["description"].append("Barycentric correction velocity in km/s")
        
        berv_z = berv_kms/3e5   #approx
        ext_table["extension_name"].append("BARYCORR_Z")
        ext_table["description"].append("Barycentric correction velocity in redshift (z)")
        
        bjd_tdb = float(self.header_blue['BERV_MIDPOINT']) 
        ext_table["extension_name"].append("BJD_TDB")
        ext_table["description"].append("Photon weighted midpoint, barycentric dynamical time (JD)")
        
        drift_b = float(self.header_blue['Relative_Drift'].split()[0])
        drift_r = float(self.header_red['Relative_Drift'].split()[0])
        drift = np.vstack([drift_b, drift_r])
        ext_table["extension_name"].append("DRIFT")
        ext_table["description"].append("Instrument drift velocity measured in blue and red camera in km/s")

        # OPTIONAL
        #exposuremeter: N/A at the moment
        #telemetry: N/A at the moment
        #telluric model: N/A at the moment
        #sky model: N/A at the moment
        #ancillary spectrum: N/A at the moment
        #image: N/A at the moment

        self.set_data("BARYCORR_KMS", np.array([berv_kms]))
        self.set_data("BARYCORR_Z", np.array([berv_z]))
        self.set_data("BJD_TDB", np.array([bjd_tdb]))
        self.set_data("DRIFT", drift)

        # create extension description table
        self.create_extension("EXT_DESCRIPT", "BinTableHDU", data=pd.DataFrame(ext_table))