#!/usr/bin/env python
# coding: utf-8

from astropy.io import fits
from rvdata.core.models.level2 import RV2
from rvdata.core.models.level3 import RV3


class MAROONXRV3(RV3):
    """
    Data model and reader for RVData Level 3 (stitched spectrum) data
    constructed from MAROON-X Level 2 data.

    This class extends the `RV3` base class to handle MAROON-X data.
    It reads a single camera's RVData standard Level 2 FITS file (blue or red)
    and produces a stitched 1D spectrum for the science fibers of that camera.

    MAROON-X has 5 traces per camera and a virtual trace:
    - TRACE1: Sky fiber (SKY)
    - TRACE2, TRACE3, TRACE4: Science fibers (SCI)
    - TRACE5: Simultaneous calibration fiber (CAL)
    - TRACE6: Flux weighted combination of the 3 science fibers (VIRTUAL)

    The three science traces (TRACE2, TRACE3, TRACE4) and virtual trace
    (TRACE6) are stitched individually and stored in STITCHED_CORR_TRACE{n}_*
    extensions (e.g., STITCHED_CORR_TRACE2_FLUX, STITCHED_CORR_TRACE2_WAVE).

    Parameters
    ----------
    Inherits all parameters from :class:`RV3`.

    Notes
    -----
    This class requires an RVData standard Level 2 FITS file as input.
    Blue and red camera data must be passed separately, one at a time.
    The classmethod `from_fits` should be used to instantiate the object.
    The `_read` method is not intended to be called directly by users.

    Example
    -------
    >>> from rvdata.core.models.level3 import RV3
    >>> blue_l3 = RV3.from_fits("maroonxblue_SL2_YYYYMMDDTHHMMSS.fits",
    ...                          instrument="MAROONX")
    >>> blue_l3.to_fits("maroonxblue_SL3_YYYYMMDDTHHMMSS.fits")
    >>> red_l3 = RV3.from_fits("maroonxred_SL2_YYYYMMDDTHHMMSS.fits",
    ...                         instrument="MAROONX")
    >>> red_l3.to_fits("maroonxred_SL3_YYYYMMDDTHHMMSS.fits")
    """

    def _read(self, hdul2: fits.HDUList, **kwargs) -> None:
        l2obj = RV2()
        l2obj.read(hdul2, instrument=None)
        self.convert_level2_to_level3(l2obj)
