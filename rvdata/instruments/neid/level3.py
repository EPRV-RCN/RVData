from astropy.io import fits
# from astropy.table import Table
# import numpy as np
import pandas as pd
import os
# from collections import OrderedDict

# import base class
from rvdata.core.models.level3 import RV3
# from rvdata.core.models.definitions import LEVEL3_EXTENSIONS
from rvdata.core.tools import stitch_spectrum


# NEID Level3 Reader
class NEIDRV3(RV3):
    """
    Data model and reader for RVData Level 3 (RV) data constructed from
    NEID Level 2 data.

    This class extends the `RV3` base class to handle NEID
    data. It reads the relevant science and calibration extensions
    from a NEID Level 2 file and organizes them into a
    standardized format.

    Parameters
    ----------
    Inherits all parameters from :class:`RV3`.

    Attributes
    ----------
    extensions : dict
        Dictionary of all created extensions
        (e.g., 'STITCHED_CORR_TRACE1_FLUX', 'STITCHED_CORR_TRACE1_WAVE', etc.),
        mapping extension names to their data arrays.
    headers : dict
        Dictionary of headers for each extension, mapping extension names to
        their FITS headers.
    data : dict
        Dictionary of data arrays for each extension.

    Notes
    -----
    To construct an RVData Level 3 object, a NEID Level 2 FITS file is required.
    The classmethod `from_fits` should be used to
    instantiate the object from these files. The `_read` method is not intended
    to be called directly by users.

    Example
    -------
    >>> from rvdata.instruments.neid.level3 import NEIDRV3
    >>> obj = NEIDRV3.from_fits("neid_L2.fits"")
    >>> obj.to_fits("neid_L3_standard.fits")
    """

    def _read(self, hdul2: fits.HDUList, **kwargs) -> None:

        # read the wavelength, flux, and blaze data
        sci_flx = hdul2["SCIFLUX"].data  # 4-116 order in NEID out of 122
        sci_wav = hdul2["SCIWAVE"].data
        sci_blz = hdul2["SCIBLAZE"].data

        # stitch the orders
        st_wave, st_flux = stitch_spectrum.stitch_orders(
            sci_wav, sci_flx, sci_blz, inst_stitch_config_sel="NEID"
        )

        # save the stitched spectrum
        self.set_data("STITCHED_CORR_TRACE1_FLUX", st_flux)
        self.set_data("STITCHED_CORR_TRACE1_WAVE", st_wave)

        # set the primary header
        hmap_path = os.path.join(os.path.dirname(__file__), "config/header_map.csv")
        headmap = pd.read_csv(hmap_path, header=0)
        phead = fits.PrimaryHDU().header
        ihead = hdul2["PRIMARY"].header
        for i, row in headmap.iterrows():
            skey = row["STANDARD"]
            neidkey = row["INSTRUMENT"]
            if pd.notnull(neidkey):
                neidval = ihead[neidkey]
            else:
                neidval = row["DEFAULT"]
            if pd.notnull(neidval):
                phead[skey] = neidval
            else:
                phead[skey] = None

        self.set_header("PRIMARY", phead)
        self.set_header("INSTRUMENT_HEADER", ihead)

        # self.set_header("DRP_CONFIG", OrderedDict(hdul2["CONFIG"].header))
        # self.set_data("DRP_CONFIG", Table(hdul2["CONFIG"].data).to_pandas())

        # self.set_header("RECEIPT", OrderedDict(hdul2["RECEIPT"].header))
        # self.set_data("RECEIPT", Table(hdul2["RECEIPT"].data).to_pandas())

        # all_exts = list(self.extensions.keys())
        # for ext_name in all_exts:
        #     if ext_name not in LEVEL3_EXTENSIONS["Name"].values:
        #         self.del_extension(ext_name)
