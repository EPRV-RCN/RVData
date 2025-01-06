from astropy.io import fits
# import base class
from core.models.level2 import RV2
import instruments.espresso.utils as utils
# KPF Level2 Reader
class ESPRESSORV2(RV2):
    """
    Read an ESPRESSO raw file and convert it to the EPRV standard format Python object. In order to run fully, this
    translator also needs to open other data products from the ESPRESSO pipeline, that should be stored in the same
    directory. The ESPRESSO pipeline products are: 2D spectrum on fiber A and B and blaze function for each fiber.
    This translator assumes that the files are named according to the ESPRESSO pipeline convention, namely:
    raw_file.fits, r.raw_file_S2D_BLAZE_A.fits, r.raw_file_S2D_BLAZE_B.fits. The names of the BLAZE files are
    obtained from the primary header.


    Methods
    -------
    _read(hdul: fits.HDUList) -> None:
        This function calls the do_conversion method from the utils.py file, which handles the conversion
    Attributes
    ----------
    extensions : dict
        A dictionary containing all the created extensions (e.g., `C1_SCI1`, `C1_SKY1`, `C2_CAL1`)
        where the keys are the extension names and the values are `SpectrumCollection` objects
        for each respective dataset.

    header : dict
        A dictionary containing metadata headers from the FITS file, with each extension's
        metadata stored under its respective key.

 
        

    Example
    -------
    >>> from core.models.level2 import RV2
    >>> rv2_obj = RV2.from_fits("ESPRE.2000-00-00T00:00:00.fits")
    >>> rv2_obj.to_fits("standard_level2.fits")
    """

    def _read(self, hdul: fits.HDUList) -> None:
        utils.do_conversion(self)
        return
        