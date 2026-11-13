"""
RVData/rvdata/instruments/harps/level3.py

Level 3 translator for HARPS, mirroring the ESPRESSO one: the DRS S1D
products become the stitched level 3 spectra.
"""

from astropy.io import fits
import os
from rvdata.core.models.level3 import RV3
import rvdata.instruments.harps.config.config as config
from rvdata.instruments.harps.utils import (
    get_files_names,
    create_PRIMARY,
    validate_fits_file,
    convert_S1D,
)


class HARPSRV3(RV3):
    """
    Read a HARPS raw file and its DRS S1D products and convert them into the
    EPRV standard level 3 format.

    Example
    -------
        >>> from rvdata.instruments.harps.level3 import HARPSRV3
        >>> rv3_obj = HARPSRV3.from_fits("harps_raw_file.fits")
        >>> rv3_obj.to_fits("standard_level3.fits")
    """

    def do_conversion(
        self, hdul: fits.HDUList, directory_structure: str = "standard"
    ) -> None:
        """
        Validate the raw file, find its S1D products, build the ``PRIMARY``
        header and store the stitched spectra.

        Parameters
        ----------
        hdul : fits.HDUList
        directory_structure : str
            Type of database architecture that stores resources. Must be either
            'dace' or 'standard'.

        :noindex:
        """
        path = os.path.join(self.dirname, self.filename)
        validate_fits_file(path)

        names = get_files_names(path, directory_structure, level=3)
        with fits.open(names["s1d_A"]) as hdul_s1d:
            self.set_header("INSTRUMENT_HEADER", hdul_s1d["PRIMARY"].header)

        nb_trace = 2
        create_PRIMARY(self, names, nb_trace, config.slice_nb, level=3)
        convert_S1D(self, names)
