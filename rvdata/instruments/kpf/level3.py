from astropy.io import fits

# import base class
from rvdata.core.models.level3 import RV3
from rvdata.core.models.level2 import RV2


# KPF Level3 Reader
class KPFRV3(RV3):
    """
    Data model and reader for RVData Level 3 (stitched spectrum) data
    constructed from KPF Level 2 data.

    This class extends the `RV3` base class to handle Keck Planet Finder (KPF)
    data. It reads the relevant science and calibration extensions
    from a KPF Level 2 file and produces a stitched 1D spectrum.

    KPF has 5 traces:
    - TRACE1: Calibration fiber (CAL)
    - TRACE2, TRACE3, TRACE4: Science fibers (SCI1, SCI2, SCI3)
    - TRACE5: Sky fiber (SKY)

    The three science traces (TRACE2, TRACE3, TRACE4) are stitched individually.
    When multiple traces are present, they are stored in STITCHED_CORR_TRACE{n}_*
    extensions (e.g., STITCHED_CORR_TRACE2_FLUX, STITCHED_CORR_TRACE3_FLUX, etc.).

    Note: The STITCHED_CORR_SCI_* extensions are currently not populated when
    multiple traces are present. A future enhancement will co-add all science
    traces to produce these combined "SCI" extensions.

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
    To construct an RVData Level 3 object, an RVData-standard Level 2 FITS file
    is required. Native KPF Level 2 files are not directly supported; they must
    first be converted to RVData-standard format using KPFRV2.from_fits().

    The classmethod `from_fits` should be used to instantiate the object from
    these files. The `_read` method is not intended to be called directly by users.

    Example
    -------
    >>> from rvdata.instruments.kpf.level3 import KPFRV3
    >>> obj = KPFRV3.from_fits("kpf_L2_standard.fits")
    >>> obj.to_fits("kpf_L3_standard.fits")
    """

    def _read(self, hdul2: fits.HDUList, **kwargs) -> None:
        # Read the RVData L2 standard file using base RV2 class
        # (KPFRV2 is for native KPF L0+L1 files, not standard L2 files)
        l2obj = RV2()
        l2obj.read(hdul2, instrument=None)
        self.convert_level2_to_level3(l2obj)
