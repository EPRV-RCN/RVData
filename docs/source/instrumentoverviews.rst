

.. |missing| replace:: **TBD**

Overviews of EPRV Instruments
******************************

The standardized EPRV data format allows the flexibility to incorporate instrument-specific content (e.g. single or multiple spectral traces, 
ancillary images), pipeline specific corrections (e.g. telluric corrections, barycentric corrections), and/or RV measurements (e.g. different 
CCF masks, line-by-line, or other custom RV calculations) at each data level. In order to capture key differences, each Instrument Team is 
encouraged to provide a ReadtheDocs file following the template below. These details are not meant to be exhaustive or to supersede individual 
instrument documentation, but rather to provide end users with a quick introduction to key instrument design and pipeline algorithm decisions.


.. ESPRESSO ENTRY ..

ESPRESSO
=================

**Instrument Details**

* Instrument Name: ESPRESSO
* Telescope/Observatory: VLT-Paranal
* Wavelength Range: 380 - 788 nm
* Resolution: 70,000 - 190,000 (depending on the mode)
* Fiber vs. slit: Fiber

**Documentation & Data Locations**

* Link to instrument manual: https://www.eso.org/sci/facilities/paranal/instruments/espresso/ESPRESSO_User_Manual_P111_v2.pdf
* Link to data reduction pipeline manual: https://www.eso.org/sci/facilities/paranal/instruments/espresso/doc.html
* Link to data archive: https://archive.eso.org/cms.html
* Person/email for translator maintenance & bug reports: Emile Fontanet [emile.fontanet@unige.ch]
* Naming convention for all data standard products:
   * L2 : ESPRE_L2_2022-08-24T03:37:56.276.fits
   * L3 : ESPRE_L3_2022-08-24T03:37:56.276.fits
   * L4 : ESPRE_L4_2022-08-24T03:37:56.276.fits
* Date of last translator update and current version number: June 2025 

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - ESPRESSO18
     - 
     - 
     - 
   * - ESPRESSO19
     - 
     - 
     - 

**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1: Science fiber, slice #1
   * Trace 2: Science fiber, slice #2
   * Trace 3: Calibration fiber, slice #1
   * Trace 4: Calibration fiber, slice #2
   * Wavelength Solution: 
   * Extraction Type: 
* Corrections and Models
   * Barycentric Correction: applied to all orders based on the same BERV value
* Ancillary Data
   * Exposure Meter: 
   * Ancillary Spectrum: 

**Useful details or tutorials for working with 1-D data**

TBD


----------


..  EXPRES Entry
EXPRES
=================

**Instrument Details**

* Instrument Name:  EXPRES (Extreme Precision Spectrograph)
* Telescope/Observatory: Lowell Discovery Telescope, Lowell Observatory
* Wavelength Range: 390-780 nm
* Resolution: 137,000
* Fiber vs. slit: fiber

**Documentation & Data Locations**

* Link to instrument manual: N/A
* Link to data reduction pipeline manual: https://ui.adsabs.harvard.edu/abs/2020AJ....159..187P/abstract
* Link to data archive: 
* Person/email for translator maintenance & bug reports: Lily Ling Zhao ()
* Naming convention for all data standard products:
   * L2 : 
   * L3 : 
   * L4 : 
* Date of last translator update and current version number: 

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - 1.1.0
     - April, 15, 2018
     - May, 5, 2018
     - 
   * - 1.2.0 
     - May 23, 2018
     - June 15, 2018
     - 
   * - 2.1.0
     - June 16, 2018 
     - June 28, 2018
     - 
   * - 2.2.0
     - June 29, 2018 
     - July 15, 2018
     - 
   * - 3.1.0 
     - July 16, 2018 
     - November 18, 2018
     - 
   * - 3.2.0 
     - December 18, 2018 
     - December 31, 2018
     - 
   * - 4.1.0 
     - February 5, 2019 
     - June 4, 20191
     - 
   * - 4.2.0 
     - June 5, 2019 
     - July 1, 2019
     - 
   * - 5.1.0 
     - August 10, 2019 
     - October 4, 2019
     - 
   * - 5.2.0 
     - October 5, 2019 
     - December 18, 2019
     - 
   * - 5.3.0 
     - December 23, 2019 
     - June 24, 2020
     - 
   * - 5.4.0
     - June 25, 2020 
     - August 14, 2020
     - 
   * - 5.5.0 
     - August 15, 2020 
     - December 14, 2020
     - 
   * - 5.6.0 
     - December 15, 2020 
     - April 27, 2021
     - 
   * - 5.7.0 
     - April 28, 2021 
     - March 2, 2022
     - 
   * - 5.8.0 
     - March 2, 2022 
     - June 30, 2022
     - 
   * - 5.9.0 
     - June 30, 2022 
     - July 15, 2023
     - 
   * - 5.10.0 
     - July 15, 2023 
     - July 30, 2024
     - 
   * - 5.11.0
     - July 30, 2024 
     - Present
     - 

**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1: Science fiber
   * Custom1_Trace1_Wave: excalibur wavelengths (https://ui.adsabs.harvard.edu/abs/2021AJ....161...80Z/abstract)
   * Wavelength Solution: polynomial
   * Extraction Type: flat-relative
* Corrections and Models
   * Barycentric Correction: Per-pixel correction derived with barycorrpy applied to chromatic exposure meter
   * Drift Correction: None
   * Telluric Model: SELENITE (https://ui.adsabs.harvard.edu/abs/2019AJ....157..187L/abstract)
   * Sky Model: None
* Ancillary Data
   * Exposure Meter: R~100 spectra (450-710 nm) taken every second (data shape is number of seconds by number of wavelength bins)

**Useful details or tutorials for working with 1-D data**

TBD


----------


G-CLEF
=================

**Instrument Details**

* Instrument Name: G-CLEF (Giant Magellan Telescope Consortium Large Earth Finder)
* Telescope/Observatory: Magellan Clay Telescope (2026+) then move to Giant Magellan Telescope (2032+)
* Wavelength Range (A): 350 - 950 nm
* Resolution: 108,000 in EPRV mode
* Fiber vs. slit: multiple fiber traces

**Documentation & Data Locations**

* Link to instrument manual: TBD
* Link to data reduction pipeline manual: TBD
* Link to data archive: TBD
* Person/email for translator maintenance & bug reports: Cem Onyuksel (cem.onyuksel@cfa.harvard.edu)


----------


HARPS
=================

**Instrument Details**

* Instrument Name: HARPS
* Telescope/Observatory: 3.6m Telescope-La Silla
* Wavelength Range: 378 - 691 nm
* Resolution: 115k
* Fiber vs. slit: Fiber

**Documentation & Data Locations**

* Link to instrument manual: https://www.eso.org/sci/facilities/lasilla/instruments/harps/doc.html
* Link to data reduction pipeline manual: https://www.eso.org/sci/facilities/lasilla/instruments/harps/doc.html
* Link to data archive: https://archive.eso.org/cms.html
* Person/email for translator maintenance & bug reports: Emile Fontanet (emile.fontanet@unige.ch)
* Naming convention for all data standard products:
   * L2 : HARPS_L2_2022-08-24T03:37:56.276.fits
   * L3 : HARPS_L3_2022-08-24T03:37:56.276.fits
   * L4 : HARPS_L4_2022-08-24T03:37:56.276.fits
* Date of last translator update and current version number: June 2025

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - HARPS03
     - February 2003
     - June 3, 2015
     - First light to octagonal fiber upgrade
   * - HARPS15
     - June 3, 2015
     - Present
     - After the octagonal fiber upgrade


**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1: Science Fiber
   * Trace 2: Calibration Fiber
   * Wavelength Solution: 
   * Extraction Type: 
* Corrections and Models
   * Barycentric Correction: Correction applied to all orders based on the same BERV value.
* Ancillary Data
   * 

**Useful details or tutorials for working with 1-D data**

TBD


----------


HARPS-N
=================

**Instrument Details**

* Instrument Name: HARPS-N
* Telescope/Observatory: Telescopio Nazionale Galileo, Roque de los Muchachos Observatory (La Palma)
* Wavelength Range: 383 - 690 nm
* Resolution: 115k
* Fiber vs. slit: Fiber

**Documentation & Data Locations**

* Link to instrument manual: https://www.tng.iac.es/instruments/harps/data/usermanv3.1.pdf
* Link to data reduction pipeline manual: http://www.tng.iac.es/instruments/harps/data/HARPS-N_DRSUserManual_1.1.pdf
* Link to data archive: http://archives.ia2.inaf.it/tng/
* Person/email for translator maintenance & bug reports: Emile Fontanet (emile.fontanet@unige.ch)
* Naming convention for all data standard products:
   * L2 : HARPN_L2_2022-08-24T03:37:56.276.fits
   * L3 : HARPN_L3_2022-08-24T03:37:56.276.fits
   * L4 : HARPN_L4_2022-08-24T03:37:56.276.fits
* Date of last translator update and current version number: June 2025

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - 
     - 
     - 
     - 


**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1: Science Fiber
   * Trace 2: Calibration Fiber
   * Wavelength Solution: 
   * Extraction Type: 
* Corrections and Models
   * Barycentric Correction: Correction applied to all orders based on the same BERV value
* Ancillary Data
   * 
**Useful details or tutorials for working with 1-D data**

TBD


----------


HERMES
=================

**Instrument Details**

* Instrument Name: HERMES (High Efficiency and Resolution Mercator Echelle Spectrograph)
* Telescope/Observatory: Mercator Telescope
* Wavelength Range: 377-900 nm
* Resolution: 85,000
* Fiber vs. slit: Fiber 

**Documentation & Data Locations**

* Link to instrument manual: https://www.mercator.iac.es/static/doc/hermes_aa.pdf
* Link to data reduction pipeline manual: https://www.mercator.iac.es/instruments/hermes/drs/
* Link to data archive: N/A
* Person/email for translator maintenance & bug reports: saskia.prins@kuleuven.be
* Naming convention for all data standard products:
   * L2 : 
   * L3 : 
   * L4 : 
* Date of last translator update and current version number: 

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - 1.0 
     - 2009
     - July 2018
     - First generation fiber link
   * - 2.0
     - July 2018
     - Present
     - Second generation fiber link

**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1:
   * Trace 2:
   * Trace 3:
   * Wavelength Solution: 
   * Extraction Type: 
* Corrections and Models
   * Barycentric Correction: 
* Ancillary Data
   * 

**Useful details or tutorials for working with 1-D data**

TBD


----------


iLocater
=================

**Instrument Details**

* Instrument Name: iLocater
* Telescope/Observatory: Large Binocular Telescope, Mt. Graham International Observatory
* Wavelength Range: 0.97-1.31 microns
* Resolution: 190,000
* Fiber vs. slit: Fiber

**Documentation & Data Locations**

* Link to instrument manual: TBD, but a description of the instrument can be found at https://ilocater.nd.edu/design/
* Link to data reduction pipeline manual: We will be modifying the APERO pipeline, there is not yet iLocater specific 
   documentation but the existing documentation can be found at http://apero.exoplanets.ca/main/spirou/main_spirou.html
* Link to data archive: TBD
* Person/email for translator maintenance & bug reports: Marshall Johnson (johnson.7240@osu.edu)


----------


KPF
=================

**Instrument Details**

* Instrument Name: KPF
* Telescope/Observatory: Keck
* Wavelength Range: 445 - 870 nm
* Resolution: 100k
* Fiber vs. slit: fiber

**Documentation & Data Locations**

* Link to instrument manual: https://kpf-pipeline.readthedocs.io/en/latest/
* Link to data reduction pipeline manual: https://kpf-pipeline.readthedocs.io/en/latest/
* Link to data archive: https://koa.ipac.caltech.edu/
* Person/email for translator maintenance & bug reports: BJ Fulton [bjfulton@ipac.caltech.edu]
* Naming convention for all data standard products:
   * L2 : KPFL2_20250208T045125.fits
   * L3 : KPFL3_20250208T045125.fits
   * L4 : KPFL4_20250208T045125.fits
* Date of last translator update and current version number: 2025-06-17 v0.2.0

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - 0.5
     - 2000-01-01 0:00:00
     - 2022-11-09 0:00:00
     - Before First Light (engineering)
   * - 1
     - 2022-11-09 0:00:01
     - 2024-02-03 0:00:00
     - First science era: First Light to Service Mission #1
   * - 1.5
     - 2024-02-03 0:00:01
     - 2024-02-11 0:00:00
     - Service Mission #1 (engineering)
   * - 2
     - 2024-02-23 12:00:00
     - 2024-11-01 0:00:00
     - Second science era: After Service Mission #1
   * - 2.5
     - 2024-11-01 12:00:00 
     - 2025-01-01 0:00:00
     - Service Mission #2 (engineering) which extended into on-sky science observations because of quasi-fixed pattern noise
   * 
     - 2.6
     - 2025-01-01 12:00:00
     - 2025-03-28 23:00:00
     - After a full Green-side warmup (CCR pressure issue)
   *
     - 2.7
     - 2025-03-28 23:00:01
     - 2025-04-23 0:22:00
     - Service Mission #3 (engineering) (Red CCR installed; Green CCR misbehaving; quasi-fixed pattern noise)
   *
     - 3
     - 2025-04-23 0:22:01
     - 3000-01-01 0:00:00
     - After Service Mission #3 (quasi-fixed pattern noise reduced)

**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1: Calibration fiber
   * Trace 2: Science slice 1
   * Trace 3: Science slice 2
   * Trace 4: Science slice 3
   * Trace 5: Sky fiber
   * Wavelength Solution: 9th order Legendre polynomial fit to LFC data then interpolated between 2 nearest LFC wavelength solutions
   * Extraction Type: optimal extraction
* Corrections and Models
* Ancillary Data
   * Exposure Meter: Table of exposure meter measurements. One column per wavelength bin
   * Ancillary Spectrum: Separate spectrograph covering only the CaHK region. Stored in the CA_HK_SCI_WAVE and CA_HK_SCI_FLUX extensions

**Useful details or tutorials for working with 1-D data**

TBD


----------


..  MAROON-X Entry
MAROON-X
=================

**Instrument Details**

* Instrument Name: MAROON-X
* Telescope/Observatory: Gemini North
* Wavelength Range: 500 - 920 nm 
* Resolution: 85,000
* Fiber vs. slit: Fiber

**Documentation & Data Locations**

* Link to instrument manual: N/A
* Link to data reduction pipeline manual: N/A
* Link to data archive: https://archive.gemini.edu
* Person/email for translator maintenance & bug reports: Tanya Das (tanyadas@uchicago.edu)
* Naming convention for all data standard products:
   * L2 : 
   * L3 : 
   * L4 : 
* Date of last translator update and current version number: 06-20-2025; v1.0.0

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - 
     - 
     - 
     - 


**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1: Sky fiber 1 (currently empty)
   * Trace 2: Science fiber 2. Contains data from the blue and red channels stacked
   * Trace 3: Science fiber 3. Contains data from the blue and red channels stacked
   * Trace 4: Science fiber 4. Contains data from the blue and red channels stacked
   * Trace 5: Simultaneous calibration fiber 5 with etalon. Contains data from the blue and red channels stacked
   * Trace 6: Virtual fiber. Created by resampling fiber 2 and 4 onto the wavelength grid of fiber 3, then applying outlier rejection and proper weighting to combine the 3 fibers into a virtual one.
   * Wavelength Solution: Non-parametric
   * Extraction Type: Fibers 2-4 uses flat optimal extraction. Fiber 5 uses box (sum) extraction.
* Corrections and Models
   * Barycentric Correction: 
* Ancillary Data
   * 

**Useful details or tutorials for working with 1-D data**

TBD


----------


MARVEL
=================

**Instrument Details**

* Instrument Name: MARVEL (Mercator Array for Radial VELocity observations)
* Telescope/Observatory: 4-telescope MARVEL array at Mercator Observatory
* Wavelength Range: 380 - 950 nm
* Resolution: 90,000 or 135,000
* Fiber vs. slit: fiber (4 input fibers, 1 per telescope)

**Documentation & Data Locations**

* Link to instrument manual: https://fys.kuleuven.be/ster/instruments/marvel
* Link to data reduction pipeline manual: N/A
* Link to data archive: N/A
* Person/email for translator maintenance & bug reports: saskia.prins@kuleuven.be
* Naming convention for all data standard products:
   * L2 : 
   * L3 : 
   * L4 : 
* Date of last translator update and current version number: 

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - 
     - 
     - 
     - 


**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1:
   * Trace 2:
   * Trace 3:
   * Wavelength Solution: 
   * Extraction Type: 
* Corrections and Models
   * Barycentric Correction: 
* Ancillary Data
   * 

**Useful details or tutorials for working with 1-D data**

TBD


----------


NEID
=================

**Instrument Details**

* Instrument Name: NEID
* Telescope/Observatory: WIYN 3.5 m Telescope/Kitt Peak National Observatory
* Wavelength Range: 380-930 nm
* Resolution: R~110,000 (High Resolution Mode); R~70,000 (High Efficiency Mode)
* Fiber vs. slit: fiber-fed

**Documentation & Data Locations**

* Link to instrument manual: 
* Link to data reduction pipeline manual: https://neid.ipac.caltech.edu/docs/NEID-DRP/
* Link to data archive: https://neid.ipac.caltech.edu 
* Person/email for translator maintenance & bug reports: Chad Bender (cbender@arizona.edu) 
* Naming convention for all data standard products:
   * L2 : 
   * L3 : 
   * L4 : 
* Date of last translator update and current version number: 

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - Start Date
     - End Date
     - Comments
   * - 1.0
     - November 6, 2019
     - November 30, 2019
     - 
   * - 2.0
     - November 30, 2019 
     - February 27, 2020
     - 
   * - 3.0
     - February 27, 2020 
     - March 31, 2020
     - 
   * - 4.0
     - October 26, 2020 
     - June 16, 2022
     -  
   * - 5.0 
     - October 18, 2022 
     - August 19, 2024
     - 
   * - 6.0
     - August 24, 2024 
     - Present
     - 

Note: these are not in the default NEID data product headers, but are added in the translator

**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1: Science fiber (HR and HE modes)
   * Trace 2: Sky fiber (HR and HE modes)
   * Trace 3: Calibration fiber (only for HR mode)
   * Wavelength Solution: polynomial
   * Extraction Type: Flat relative extraction is used (optimal for continuum sources)
* Corrections and Models
   * Barycentric Correction: Barycorrpy is used to calculate a per-order barycentric correction at the per-order flux weighted midpoint as 
      measured by the chromatic exposure meter. The relevant extensions include array data of length 122, with a value for each echelle order.
   * Drift Correction: The provided drift extension is a single value in km/s that represents the drift relative to the calibration session at 
      the start of the observing window (evening calibrations for nighttime, morning calibrations for solar).
   * Telluric Model: a line and continuum absorption telluric model. The model is only applicable for zero-indexed orders 55 through 110. The rest 
      of the orders have the model values set to 1. The model includes water and molecular oxygen.
* Ancillary Data
   * Exposure Meter: This extension includes processed chromatic exposure meter data in table format. Rows are temporal sampling and columns are 
      wavelength sampling. The first column are the timestamps for each exposure meter frame. The rest of the columns are labeled with the wavelength 
      corresponding to that spectral bin. 

**Useful details or tutorials for working with 1-D data**

TBD


----------


NIRPS
=================

**Instrument Details**

* Instrument Name: NIRPS
* Telescope/Observatory:  La Silla 3.6m 
* Wavelength Range: 0.95 to 1.8 microns
* Resolution: HE: 90000, HA: 100000
* Fiber vs. slit: Fiber

**Documentation & Data Locations**

* Link to instrument manual: https://www.eso.org/public/teles-instr/lasilla/36/nirps/ 
* Link to data reduction pipeline manual: http://apero.exoplanets.ca 
* Link to data archive: https://archive.eso.org/cms.html 
* Person/email for translator maintenance & bug reports: 
* Naming convention for all data standard products:
   * L2 : 
   * L3 : 
   * L4 : 
* Date of last translator update and current version number: 

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - 
     - 
     - 
     - 

**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1:
   * Trace 2:
   * Trace 3:
   * Wavelength Solution: 
   * Extraction Type: 
* Corrections and Models
   * Barycentric Correction: 
* Ancillary Data
   * 

**Useful details or tutorials for working with 1-D data**

TBD


----------


SPIRou
=================

**Instrument Details**

* Instrument Name: SPIRou
* Telescope/Observatory: CFHT
* Wavelength Range: 0.9245 to 2.4 microns 
* Resolution: 80,000
* Fiber vs. slit: Fiber

**Documentation & Data Locations**

* Link to instrument manual: https://www.cfht.hawaii.edu/Instruments/SPIRou/
* Link to data reduction pipeline manual: http://apero.exoplanets.ca 
* Link to data archive: https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/search/?collection=CFHT&noexec=true 
* Person/email for translator maintenance & bug reports: 
* Naming convention for all data standard products:
   * L2 : 
   * L3 : 
   * L4 : 
* Date of last translator update and current version number: 

**Instrument era (INSTERA keyword) date ranges**

.. list-table::
   :widths: 25 25 25 50 
   :header-rows: 1

   * - INSTERA
     - UT_start_date
     - UT_end_date
     - Comments
   * - 
     - 
     - 
     - 


**Change log**

TBD

**Instrument-Specific Level 2 Header and Extension Details**

* Fundamental Parameters
   * Trace 1:
   * Trace 2:
   * Trace 3:
   * Wavelength Solution: 
   * Extraction Type: 
* Corrections and Models
   * Barycentric Correction: 
* Ancillary Data
   * 

**Useful details or tutorials for working with 1-D data**

TBD