"""RVData Level 3 reader for EXPRES: stitched 1D spectrum built from the L2 translation."""

from astropy.io import fits

from rvdata.core.models.level3 import RV3
from rvdata.instruments.expres.level2 import EXPRESRV2


class EXPRESRV3(RV3):
    """
    Data model and reader for RVData Level 3 (stitched spectrum) data
    constructed from an EXPRES fitspec file.

    The fitspec file is first translated to Level 2 with
    :class:`~rvdata.instruments.expres.level2.EXPRESRV2`, then the science
    trace is deblazed and stitched onto a common log-wavelength grid by
    :meth:`RV3.convert_level2_to_level3` using the stitching parameters in
    ``config/expres_level3.config``.

    Example
    -------
    >>> from rvdata.core.models.level3 import RV3
    >>> l3 = RV3.from_fits("fitspec/Sun_20240126.5073.fits", instrument="EXPRES")
    >>> l3.to_fits()
    """

    def _read(self, hdul: fits.HDUList, **kwargs) -> None:
        l2obj = EXPRESRV2()
        l2obj.read(hdul, instrument="EXPRES")
        self.convert_level2_to_level3(l2obj)
