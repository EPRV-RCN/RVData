"""
RVData/rvdata/instruments/harps/level4.py

Level 4 translator for HARPS, mirroring the ESPRESSO one: the DRS CCF
products become the level 4 RV, CCF and diagnostics tables.
"""

from astropy.io import fits
import os
from rvdata.core.models.level4 import RV4
import rvdata.instruments.harps.config.config as config
from rvdata.instruments.harps.utils import (
    get_files_names,
    create_PRIMARY,
    validate_fits_file,
    convert_CCF,
)


class HARPSRV4(RV4):
    """
    Read a HARPS raw file and its DRS CCF products and convert them into the
    EPRV standard level 4 format.

    Example
    -------
        >>> from rvdata.instruments.harps.level4 import HARPSRV4
        >>> rv4_obj = HARPSRV4.from_fits("harps_raw_file.fits")
        >>> rv4_obj.to_fits("standard_level4.fits")
    """

    def do_conversion(
        self, hdul: fits.HDUList, directory_structure: str = "standard"
    ) -> None:
        """
        Validate the raw file, find its CCF products, build the ``PRIMARY``
        header and store the RV, CCF and diagnostics tables.

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

        names = get_files_names(path, directory_structure, level=4)
        with fits.open(names["ccf_file"]) as hdul_ccf:
            self.set_header("INSTRUMENT_HEADER", hdul_ccf["PRIMARY"].header)

        nb_trace = 2
        create_PRIMARY(self, names, nb_trace, config.slice_nb, level=4)
        convert_CCF(self, names)
