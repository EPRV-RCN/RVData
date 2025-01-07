import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.constants import c
from collections import OrderedDict
import pandas as pd

# import base class
from core.models.level2 import RV2

expres_epochs, epoch_start_isot = np.loadtxt(
    "./config/expres_epochs.csv", delimiter=",", skiprows=1, dtype="str"
).T
epoch_start_mjd = Time(epoch_start_isot).mjd

header_map = pd.read_csv("./config/expres_header_map.csv").set_index("standard")
static_headers = {
    "ORGANIZA": "Yale",
    "DATALVL": "L2",
}
obstype_map = {
    'Science',: "Sci",
    'Calibration',: "Cal",
    'Dark',: "Cal",
    'ThAr',: "Cal",
    'Quartz',: "Cal",
    'LFC': "Cal",
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
        data = hdul[1].data.copy()
        primary_header = hdul[0].header
        fitspec_header = hdul[1].header
        expmeter_header = hdul[2].header
        expmeter_data   = hdul[2].data.copy()
        itrace = 1

        # Primary with just EPRV Standard FITS Headers
        self.create_extension("PRIMARY", self.standardizeExpresHeader(hdul))

        # Original Instrument Header
        self.create_extension("INSTRUMENT_HEADER", primary_header)

        # Receipt?

        # DRP Config?

        # Flux
        blaze = data["blaze"]
        spec = data["spectrum"] * blaze
        self.create_extension(f"TRACE{itrace}_FLUX", "ImageHDU", data=spec)
        # cont = data['continuum']

        # Wavelength
        wave = data["excalibur"]  # data['wavelength']
        self.create_extension(f"TRACE{itrace}_WAVE", "ImageHDU", data=wave)

        # Variance
        variance = data["uncertainty"] ** 2
        self.create_extension(f"TRACE{itrace}_VAR", "ImageHDU", data=variance)

        # Blaze
        self.create_extension(f"TRACE{itrace}_BLAZE", "ImageHDU", data=blaze)

        # Barycentric Correction
        bary_arr = data["bary_excalibur"]  # data['bary_wavelength']
        berv_kms = 1 - bary_arr / wave * c.to("km/s")
        berv_z = 1 - bary_arr / wave
        self.create_extension("BARYCORR_KMS", "ImageHDU", data=berv_kms)
        self.create_extension("BARYCORR_Z", "ImageHDU", data=berv_z)

        # Photon Weighted Midpoint
        #   (FORMAT DOES NOT MATCH BARYCORR_KMS)
        self.create_extension("BJD_TDB", "ImageHDU", data=expmeter_header["HIERARCH wtd_mdpt"])

        # Instrument Drift Map (REQUIRED)
        # EXPRES doesn't really calculate this...
        # I could make this a polynomial v. excalibur wavelength thing, but that's not really honest
        self.create_extension("DRIFT", "ImageHDU", data=np.zeros_like(wave)) # setting to zero assuming shifts expected

        # Exposure Meter
        self.create_extension("EXPMETER", "ImageHDU", 
                              header=expmeter_header, data=expmeter_data)

        # Telemetry (Optional)
        # Might not have this either tbh
        # self.create_extension("TELEMETRY",)

        # Telluric Model
        telluric = data['tellurics']
        self.create_extension(f"TRACE{itrace}_TELLURIC", "ImageHDU", data=telluric)

        # No sky model, ancillary spectrum, or image

    # =============================================================================
    # Methods for standardizing header keywords

    def drpFlag(hdul):
        extensions_to_check = ['spectrum','blaze',
                               'wavelength','bary_wavelength',
                               'excalibur','bary_excalibur',
                               'continuum','tellurics']
        extension_list = hdul[1].data.dtype.names
        percent_there = np.sum([extn in extension_list for extn in extensions_to_check])/len(extensions_to_check)
        if percent_there==1:
            return 'Pass'
        else:
            return 'Fail' if percent_there==0 else 'Warn'

    header_funcs = {
        'OBSTYPE': (lambda hdul: obstype_map[hdul[0].header['OBSTYPE']]),
        'BINNING': (lambda hdul: hdul[0].header["CCDBIN"].replace(' ','')[1:-1].replace(",",'x')),
        'NUMORDER': (lambda hdul: hdul[1].header['NAXIS2']),
        'FILENAME': (lambda hdul: f"EXPRESL2_{hdul[0].header['MIDPOINT'][2:-1].replace('-','')}.fits"),
        'JD_UTC': (lambda hdul: Time(hdul[0].header['DATE-SHT']).jd),
        'INSTERA': (lambda hdul: expres_epochs[np.sum(Time(hdul[0].header['MIDPOINT']).mjd >= epoch_start_mjd) - 1]),
        'INSTFLAG': (lambda hdul: 'Pass' if (bool(hdul[0].header['EXPMTR']) & bool(hdul[0].header['EXPMTR'])) else 'Fail'),
        'DRPFLAG': drpFlag,
    }

    def standardizeExpresHeader(hdul):
        head0 = hdul[0].header
        standard_head = OrderedDict()
        for key in header_map.index:
            if key in static_headers.keys(): # Keywords that never change
                standard_head[key] = static_headers[key]
            elif key in header_funcs.keys(): # Keywords that require processing
                standard_head[key] = header_funcs[key](head0)
            else:
                expres_val = head0[header_map.loc[key,'expres']]
                if header_map.loc[key,'required']=='N' and expres_val=='':
                    continue
                standard_head[key] = expres_val
        return standard_head
