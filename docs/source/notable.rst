

.. |missing| replace:: **TBD**

Noteable Extensions & Header Keywords
**************************************

The Receipt Extension
=====================

The EPRV data standard makes use of a receipt extension, inspired by the KPF data reduction pipeline. 
The receipt is structured as a pandas.Dateframe table and keeps track of the data process history, so 
that the information stored by this instance can be reproduced from the original data. Each row of the 
receipt should capture an individual data processing event.

Translators that generate/modify the content of an EPRV standard data product are expected to also write to the receipt. 
The receipt is inherited from lower order data products (e.g., L2) to higher order data products (e.g., L3) to preserve 
the data processing history.

Three string inputs are required: name, any relevant parameters, and a status. 
The receipt will also automatically fill in additional information, such as the time of execution, 
code release version, current branch, ect.

The Config Extension
====================

This extension contains a two-column dataframe that stores each line of the input configuration file.

This file is based on the configuration language that the Python ConfigPaser class implements for Python programs. 
The syntax of the configuration language primarily consists of ``key = value`` statements and ``# comments`` (or ``; comments``).

The creation of the key/value pairs, in general, follows the convention as those provided by the ConfigParser classs. 
As used in the KPF and NEID Pipelines, the keys consists of alphanumerics and underscores (e.g. KPFPIPE_TEST_DATA) while 
the value is in the form of a string, a number, a boolean, or a list.



PRIMARY Header Keywords
=======================

**MULTIPLICITY** : Keyword in the PRIMARY extension that indicates whether or not an extension can be duplicated to have
multiple HDUs based on, e.g., the number of traces.

**REQUIRED** : Keyword in the PRIMARY extension that indicates whether or not an extension must be present and 
meaningfully populated in the FITS file in order for the file to be compliant with the EPRV Data Standard.

**INSTERA** : Tag for the "instrument era" used to track permanent changes to instrument (e.g., 1.1.2). Major changes 
that result in RV offset are recorded in the first numeral space. Small but permanent changes, (changes in the calibration 
light sources, etc) that don’t necessarily introduce RV offsets but do represent an ongoing status change in the instrument, 
are recorded in the second numeral. Anything that could impact the calibration or science outputs goes here. Even smaller 
changes, or things that you only recognize retrospectively, are recorded in the third numeral. 

Note: Each team will place a csv file into the config directory with 1 line per era that has the start and end dates of 
each timespan. If this file is already maintained elsewhere by the instrument team then link to that instead so the 
file doesn’t need to be maintained in multiple places.
