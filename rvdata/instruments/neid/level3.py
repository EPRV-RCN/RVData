from astropy.io import fits

# import base class
from rvdata.core.models.level3 import RV3

# import NEID Level2 reader
from rvdata.instruments.neid.level2 import NEIDRV2


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
        (e.g., 'STITCHED_CORR_SCI_FLUX', 'STITCHED_CORR_SCI_WAVE', etc.),
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
        l2obj = NEIDRV2()
        l2obj.read(hdul2, instrument="NEID")
        self.convert_level2_to_level3(l2obj)
