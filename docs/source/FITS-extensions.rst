

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
A Level 2 data file will contain the order-by-order flux, wavelength, variance, and blaze which will be split into individual HDUs. 
If the spectrograph produces multiple traces (e.g. a calibration trace, a sky trace, and one or more science traces) then the data 
will also be split between traces. 

At Level 2 the only correction that has already been applied to the data is a drift correction. 
The information necessary to execute common EPRV corrections (e.g. blaze corrections and barycentric corrections) is provided.

This data product is intended primarily for users who are interested in the improving EPRV data reduction and post processing techinques.

.. csv-table::
    :header-rows: 1
    :file: ../../core/models/config/L2-extensions.csv



Level 3a FITS Extensions
========================
A Level 3a file is structured in the format as Level 2, but here the flux array will have been blaze corrected and the wavelength 
solution will have been both drift and barycentric corrected. The variance array will be calculated using the blaze-corrected flux. 

Optional extensions within the Level 3a format allow for users to apply additional, custom corrections to either the flux or 
wavelength arrays (see the CUSTOMCORR extensions below).

Level 3a inherits the PRIMARY extension from Level 2, along with the INSTRUMENT_HEADER, RECEIPT, and DRP_CONFIG extensions. 

This data product is intended primarily for users who are interested in improving EPRV measurement techniques.

.. csv-table::
    :header-rows: 1
    :file: ../../core/models/config/L3a-extensions.csv



Level 3b FITS Extensions
========================
A Level 3b file stitches together the individual orders from Level 3a into a single 1-D spectra that spans the entirety of the 
instrument's spectral range. As with Level 3a the flux array will have been blaze corrected and the wavelength 
solution will have been both drift and barycentric corrected. The variance array will again be calculated using the blaze-corrected flux. 

Optional extensions within the Level 3b format allow for users to apply additional, custom corrections to either the flux or 
wavelength arrays (see the CUSTOMCORR extensions below).

Level 3b inherits the PRIMARY extension from Level 2, along with the INSTRUMENT_HEADER, RECEIPT, and DRP_CONFIG extensions. 

This data product is intended primarily for users who are interested in other types of science that can be accomplished with high 
resolution visible / NIR spectra. Potential examples include science related to exoplanet atmospheres or stellar astrophysics.

.. csv-table::
    :header-rows: 1
    :file: ../../core/models/config/L3b-extensions.csv