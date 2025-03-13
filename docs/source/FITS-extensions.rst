

.. |missing| replace:: **TBD**

Overview of EPRV FITS Standard Extensions
*****************************************

General information
========================
#. The PRIMARY extension in any L2, L3, or L4 file generated as part of the EPRV Data Standard should contain the EPRV Data Standard FITS Header (Add link). 
#. The Multiplicity keyword indicates whether or not an extension can be duplicated to have multiple HDUs based on, e.g., the number of traces. 
#. The Required keyword indicates whether or not an extension must be present and meaningfully populated in the FITS file in order for the file to be compliant with the EPRV Data Standard. 


Level 2 FITS Extensions
=======================
A Level 2 data file will contain the order-by-order flux, wavelength, variance, and blaze which will be split into individual extensions. 
If the spectrograph produces multiple traces (e.g. a calibration trace, a sky trace, and one or more science traces) then the data 
will also be split into one flux extension per trace, one wavelength extension per trace, etc.

At Level 2 the only correction that has already been applied to the required wavelength solution (``TRACE1_WAVE``) is a drift correction. 
The information necessary to execute common EPRV corrections to the required wavelength and flux extensions (e.g. blaze corrections 
and barycentric corrections) is provided in separate extensions.

Optional extensions within the Level 2 format allow for users to apply additional, custom corrections to either the flux or 
wavelength solution (see the CUSTOM1_TRACE1 extensions below).

This data product is intended primarily for users who are interested in the improving EPRV data reduction and post processing techinques.

.. csv-table::
    :header-rows: 1
    :file: ../../rvdata/core/models/config/L2-extensions.csv


Level 3 FITS Extensions
========================
A Level 3 file will contain 1-D spectra that stitch together all extracted orders. The flux array will have been blaze 
corrected and the wavelength solution will have been both drift and barycentric corrected. The variance array will be 
calculated using the blaze-corrected flux. 

This stitched flux extension and its corresponding wavelength and variance extensions will exist for each of the individual traces 
(``STITCHED_CORR_TRACE1_FLUX``, ``STITCHED_CORR_TRACE1_WAVE``, and ``STITCHED_CORR_TRACE1_VAR``). An additional set of extensions 
will contain the co-added flux from all science fibers, the corresponding variance based upon that co-added flux, and the appropriate 
wavelength solution (``COADD_STITCHED_CORR_FLUX``, ``COADD_STITCHED_CORR_WAVE``, and ``COADD_STITCHED_CORR_VAR``). 

If your instrument pipeline does not already produce co-added spectra, they can be computed following Bourrier, Delisle et al., 2024. 
which makes use of the Bindensity package available at https://gitlab.unige.ch/jean-baptiste.delisle/bindensity.  
But documentation within each instrument's page on this ReadtheDocs should note exactly what co-adding process is used to combine 
the individual flux extensions.

Optional extensions within the Level 3 format allow for users to apply additional, custom corrections to either the flux or 
wavelength arrays (see the ``CUSTOMCORR`` extensions below). 

Level 3 inherits the PRIMARY extension from Level 2, along with the INSTRUMENT_HEADER, RECEIPT, and DRP_CONFIG extensions. 

This data product is intended primarily for users who are interested in other types of science that can be accomplished with high 
resolution visible / NIR spectra. Potential examples include science related to exoplanet atmospheres or stellar astrophysics.

.. csv-table::
    :header-rows: 1
    :file: ../../rvdata/core/models/config/L3-extensions.csv


Level 4 FITS Extensions
========================
A Level 4 file will contain derived data products such as RV measurements, cross correlation functions (CCFs), CCF metrics 
(e.g., FWHM, BIS), and stellar activity indicators. Only the RV measurements are required, to support non-CCF RV reduction methods.

L4 data products are expected to be derived from an L2 file, and inherit the PRIMARY, INSTRUMENT_HEADER, RECEIPT, and DRP_CONFIG 
extensions from that L2.

The RV extension is a Binary Table that contains the following columns: BJD_TDB, wave_start, wave_end, pixel_start, pixel_end, RV_Trace, RV_error, BC_vel, order index, echelle order 
Optional, additional columns include: RV_weight
The wave_start and wave_end values can be set to same number if reporting central wavelength of the segment.

The header for this extension should contain RVMETHOD, RVSTART, RVSTEP, MASK keywords.

.. csv-table::
    :header-rows: 1
    :file: ../../rvdata/core/models/config/L4-extensions.csv