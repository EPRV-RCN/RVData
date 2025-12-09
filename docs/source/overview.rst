.. |eaccute| unicode:: U+00E9
.. |Eaccute| unicode:: U+00C9
		       


Overview of the EPRV FITS Data Standard
***************************************

Purpose
=======

A new generation of Extreme Precision Radial Velocity (EPRV)
cross-dispersed echelle spectrographs have come on-line in the past
decade and are now producing large quantities of data in support of
stellar, exoplanet, and galactic astrophysics. The EPRV community is
pushing towards a goal of achieving measurement precision at the 10
cm/s level, sufficient for detecting an Earth mass planet in an
Earth-like orbit around a Sun-like star. The data from these
instruments is complicated, and must be handled carefully in each
instrument's Data Reduction Pipeline (DRP) to avoid degrading the
intrinsic data qualtiy. Numerous derived products, including radial
velocities, stellar activity indicators, line shape assessments, and
other quantities, are routinely produced by DRPs for use in scientific
analyses. Sufficient data now exists for the same targets across
multiple instruments to investigate effects from different hardware
and software choices, which can impact how current instruments are
used and also how future instruments are designed.

These EPRV instruments generate data in native formats that were
defined by individual instrument teams at the time of construction,
rather than utilizing a common standard, which has introduced
unecessary complication in processes such as building and maintaining
data archives, carrying out scientific investigations that utilize
data from multiple instruments, and making cross-comparisons of
intrument internal performances. To address these issues, the **EPRV
Data Standardization Project** has brought together an international
collaboration of instrumentalists and software developers from existing
and future echelle spectrograph projects to *define a standard data format
for multiple levels of processed echelle data*. The project has also
generated *translator* software to convert files from their
native instrument formats into the standard format.  This project
was started in 2024, and expects to release the first versions for the
*standard definition* and the *translator software* in the first half
of 2026.

New echelle spectrographs are encouraged to use the *standard* data
format for their native data. Existing instrument teams who do not
have a *translator* are encouraged to engage with the project to
generate code for their instruments.

The **EPRV Data Standardization Project** was funded by `NASA EPRV
Foundational Science Grant `22-EPRV22_2-0018 <https://nspires.nasaprs.com/external/solicitations/summary.do?solId=%7bE52C5EC5-E0FC-403E-1071-4802DB562F0C%7d&path=&method=init>`_.

Data Standard Overview
======================

The EPRV Data Standard stores data in FITS files and is compliant with
`FITS v4.0 <https://fits.gsfc.nasa.gov/fits_standard.html>`_. The
standard defines data products at different *Levels* of processing,
and places requirements on filenames, structure, and contents.

#. All compliant files will be named using the following scheme: **<inst>_SL<#>_<timestamp>.fits**

   * **<inst>** : Instrument name or abbreviation, in lower-case. (e.g., neid, kpf, expres). Different cameras are treated as separate instruments (e.g., for Maroon-X: mxred, mxblue)
   * **<SL#>** : Standard Level #, where # is typically 2, 3, 4, etc.
   * **<timestamp>** : Start of Exposure, corresponding to header DATE-OBS keyword. Limited to *YYYYMMDDTHHMMSS*. In cases where fractional seconds are required to distinguish files, add milliseconds as *sss*.

#. Standard data levels are defined as follows:

   *  Level 0 is not defined by the standard, but is envisioned to be
      the native raw data format returned by the instrument.
   *  Level 1 is not definied by the standard, but is envisioned to be
      a minimally processed image data that is assembled to resemble
      the detector focal plane, and stored in an image extension. If a
      detector has multiple amplifiers, Level 1 contains a complete
      image in a physical pixel reference frame. This data level can
      have overscan and bias corrections applied.  Details of headers
      and additional extensions are not defined by the data
      standard. It is allowable that Level 1 is identical to Level 0,
      if the lowest level product returned by the instrument is a
      fully assembled image.
   *  Level 2 will contain the 1D extracted order-by-order flux,
      wavelength, variance, and blaze in individual extensions. If the
      spectrograph produces N spectral traces, which corresponding to
      N object inputs in the detector optical input, (e.g., a
      calibration trace, a sky trace, one or more science traces),
      then each trace will be preserved in a separate set of
      extensions. The barycentric motion correction will be calculated
      and included in a separate extension, but not applied to the
      wavelength solution.
   *  Level 3 will contain a 1D spectra with all orders stitched
      together. The flux array will be blaze corrected and the
      wavelength solution will be corrected for instrument drift and
      barycentric motion. The variance array will be calculated using
      the blaze-corrected flux.
   *  Level 4 will contain derived data products such as RV
      measurements, cross-correlation functions and metrics, and
      stellar activity indicators.


Data Translator Overview
========================

Instruments writing data in a native format can translate their data
into the Standard format by utilizing code in
https://github.com/EPRV-RCN/RVData .

This repository contains RV2(), RV3(), and RV4() base classes to hold
data required for generating the standard data files. These classes
have a *write* method named *.tofits()* which outputs the data into
the standard format.

To integrate a new instrument into the translator, the instrument team
should provide the corresponding *read* method that loads their native
files into the base classes.  More details are provided in :ref:`new-reader` .

      
Team Members
============

The **EPRV Data Standardization Project** team is listed below.  Core team members who have made significant contributions to standard or software development are indicated in *italics*.

  * *Jennifer Burt, Jet Propulsion Laboratory, Prinicpal Investigator*

  * *Megan Bedell, Flatiron Institute, Co-I (EXPRES)*
  * *Chad Bender, Univeristy of Arizona, Co-I (NEID, HPF)*
  * Sean Carey, NASA Exoplanet Science Institute, Co-I
  * *Jonathan Crass, Ohio State University, Co-I (iLocater)*
  * *BJ. Fulton, NASA Exoplanet Science Institute, Co-I (NEID, KPF)*
  * *Samuel Halverson, Jet Propulsion Laboratory, Co-I (KPF, NEID, HPF)*
  * Andrew Howard, Caltech, Co-I (KPF)
  * *Daniel Krolikowski, University of Arizona, Co-I (NEID, HPF)*
  * *Sarah Logsdon, NSF NOIRLab, Co-I (NEID)*
  * *Lily Zhao, University of Chicago, Co-I (EXPRES)*

  * Vanessa Bailey, Jet Propulsion Laboratory, Collaborator
  * Charles Beichman, NASA Exoplanet Science Institute, Collaborator (PARVI)
  * Lars Buchave, Technical University of Denmark, Collaborator (2ES)
  * Rose Gibson, University of California Los Angeles, Collaborator (PARVI)
  * *Andreas Quirrenback, University of Heidelberg, Collaborator (CARMENES)*
  * Andrew Szentgyorgyi, Smithsonian Astrophysical Observatory, Collaborator (G-CLEF)

  * * |Eaccute| tienne Artigau, Universit |eaccute| \ de Montr |eaccute| al (Spirou, NIRPS)*
  * Neil Cook, Universit |eaccute| de Montr |eaccute| al (Spirou, NIRPS)
  * *Tanya Das, University of Chicago (MAROON-X)*
  * Joris De Ridder, KU Leuven (MARVEL)
  * *Xavier Dumusque, University of Geneva (HARPS, HARPS-N, EXPRESSO)*
  * *Michelle Edwards, NSF NOIRLab (KPNO)*
  * *Emile Fontanet, University of Geneva (HARPS, HARPS-N, ESPRESSO)*
  * Daniel Holdsworth, South Africian Large Telescope, (SALT-HRS)
  * *Marshall Johnson, Ohio State University (iLocater)*
  * *Gaspare lo Curto, European Southern Observatory (HARPS, ESPRESSO)*
  * Joe Ninan, TIFR (HPF, NEID)
  * *Cem Onyuksel, Smithsonian Astrophysical Observatory (G-CLEF)*
  * *Leonardo Paredes, University of Arizona, Collaborator (NEID, HPF)*
  * Saskia Prins, KU Leuven (MARVEL)
  * *Timothy Pickering, MMT (Astropy, Specutils, Specreduce, SALT-HRS)*
  * *Martino Romaniello, European Southern Observatory (ESO, ESPRESSO, HARPS, NIRPS, CRIRES)*
  * *Damien S |eaccute| gransan, University of Geneva, (HARPS, HARPS-N, NIRPS, EXPRESSO, HARPS-3)*
  * Petr Skoda, Czech Academy of Sciences (PLATOspec)
  * Julian Sturmer, University of Heidelberg (2ES)
