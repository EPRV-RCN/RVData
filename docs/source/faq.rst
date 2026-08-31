Frequently Asked Questions (FAQ)
*****************************************

Using the Translator
========================

*How do I get started using the translator for my favorite spectrograph?*

  The translators ship in the ``rv-data-standard`` package, which you can install from PyPI:

  .. code-block:: bash

      pip install rv-data-standard

  Translators currently exist for KPF, NEID, ESPRESSO, MAROON-X, EXPRES, HARPS, and HARPS-N.
  To translate a native file, load it through the base class for the data level you want, passing
  the instrument name, and then write it out in the standard format:

  .. code-block:: python

      from rvdata.core.models.level2 import RV2

      l2 = RV2.from_fits("native_file.fits", instrument="NEID")
      l2.to_fits("neid_SL2_20231010T020006.fits")

  Worked, end-to-end examples for each instrument and data level are provided in the
  :doc:`tutorials`. If your instrument does not yet have a translator, see
  `Adding a New Translator to the Data Standard`_ below.

*I have an instrument-specific translator question -or- I found a bug in an instrument translator.  Who do I reach out to for help?*

  Please open an issue on the project's GitHub repository at
  https://github.com/EPRV-RCN/RVData/issues. This is the best way to ask a question or report a
  bug: it reaches the whole development team and creates a record that others can benefit from.
  For matters that need a direct contact, reach out to the :ref:`code-review-team`. When reporting
  a translator bug, please include the instrument, the data level, the version of
  ``rv-data-standard`` you are using, and, where possible, a minimal example file that reproduces
  the problem.

Translator Applications
========================

*Can I run a translator locally on my machine?*

  Yes. The translators run entirely locally — no server or network service is required. Once you
  have ``pip install rv-data-standard``\ ed the package, you can translate files on your own
  machine as shown above. Separately, some archives (e.g., NExScI and DACE) plan to host
  pre-translated data or offer in-line translation at download time, but that is a convenience and
  not a requirement for using the translators yourself.

*Can I use a translator for solar data? What about calibration data?*

  Solar observations are supported and are in fact an explicit motivating use case for the standard
  (Sun-as-a-Star science). Solar files are translated like any other observation, and the PRIMARY
  header ``ISSOLAR`` keyword flags whether the observation is of the Sun. For calibration light:
  the Level 2 format preserves every spectral trace from an exposure, so calibration and sky traces
  recorded alongside the science trace(s) are retained as separate trace extensions. Dedicated
  calibration-only frames (e.g., raw flats or arcs) correspond to Level 0/Level 1 data, which the
  standard intentionally does not define.

*How does the translator handle different RV algorithms? What if, for example, I want a file that is a stellar template for template matching?*

  Level 4 records the radial velocities produced by the instrument's pipeline regardless of how
  they were derived; the ``RVMETHOD`` PRIMARY keyword captures the method used (e.g., ``CCF``,
  ``SERVAL``, ``LBL``). When a CCF-based method is used, the optional CCF and CCF-metric extensions
  can carry the cross-correlation functions and their summary statistics. The standard does not
  currently define a dedicated stellar-template product. If a template — or any other new product —
  would be valuable for your science, that is exactly the kind of addition to raise with the Change
  Review Board as a feature request (see :ref:`change-review-board`).

.. _adding-a-new-translator:

Adding a New Translator to the Data Standard
============================================

*I want to develop a new instrument translator. How do I get started?*

  Start with the developer workflow in :doc:`interact` (fork the repo, clone it, and branch off
  ``develop``) and the step-by-step guide in :ref:`new-reader`. That guide walks through creating an
  instrument directory, writing a ``_read`` method that maps your native FITS file onto the standard
  data model, and registering your instrument in the ``INSTRUMENT_READERS`` dictionary in
  ``definitions.py``. The existing translators under ``rvdata/instruments/`` (KPF, NEID, ESPRESSO,
  and others) are the best reference implementations. Note that continuous integration requires
  native test data to validate against — contact the :ref:`code-review-team` to have fixture files
  hosted on the project's central servers.

*What does it mean to be compliant with the standard? E.g. Observatory required keywords vs Standard keywords*

  A file is compliant when it has the structure the standard requires for its data level: the
  required FITS extensions are present and meaningfully populated, and the PRIMARY header contains
  the full set of standard keywords. Every standard keyword must be present; those marked *Required*
  must hold meaningful values, while keywords that are not required may be set to ``Null``. The
  ``Required`` and ``Multiplicity`` keywords in the PRIMARY header indicate, respectively, whether an
  extension must be present and whether it may be duplicated (for example, once per trace). Keywords
  that your observatory or pipeline needs but that are not part of the standard do **not** make a
  file non-compliant: the complete native header is preserved verbatim in the ``INSTRUMENT_HEADER``
  extension, so no information is lost. See :doc:`primaryheader` for the full keyword table.

*What version of Python should I be using? Is there Python version control?*

  The package requires Python 3.12 or newer (``requires-python >= 3.12``), and continuous
  integration tests every change against Python 3.12, 3.13, and 3.14, so any of those is supported.
  Two versioning schemes are used in the project, and they are kept independent: the **translator
  code** uses `Semantic Versioning <https://semver.org/>`_ (MAJOR.MINOR.PATCH), while the **data
  standard** itself uses `Calendar Versioning <https://calver.org/>`_ (YYYY.MM). See :ref:`versioning`
  for details.

*The current data standard doesn't accommodate a header keyword that is important for my instrument. What do I do?*

  Two things. First, no information is lost in the meantime: your instrument's complete native header
  is preserved in the ``INSTRUMENT_HEADER`` extension, so the keyword remains available to anyone
  working with your files. Second, if the keyword is broadly useful it should be considered for
  addition to the standard — submit it as a Feature Request issue on the GitHub repository, labeled
  ``datastandard``. The Change Review Board reviews such suggestions twice a year; see
  :ref:`change-review-board` for the process and timing.

Change Review Board
========================

*What is the Change Review Board?*

  The Change Review Board (CRB) is the body responsible for evolving the data standard over time. It
  includes representatives from major EPRV-community stakeholders — NASA, ESO, KPNO, Gemini, Keck,
  DACE, and NExScI — and meets twice a year to review proposed updates and decide whether a new
  version of the standard should be released. See :ref:`change-review-board`.

*How do I submit change requests for the data standard?*

  Open a Feature Request issue on the project's GitHub repository
  (https://github.com/EPRV-RCN/RVData) and label it ``datastandard``. To be taken up at a given CRB
  meeting, submit your request at least one month beforehand; requests made closer to a meeting may
  be rolled to the following one.

*What is the process for assessing change requests?*

  The CRB meets in late April (supporting an August release) and late October (supporting a February
  release). The board reviews the suggestions submitted ahead of each meeting, and may reach out to
  the author for clarification or to iterate on the idea, before deciding whether to incorporate the
  change into a new version of the standard. Approved changes are announced via the GitHub README and
  documented on this Read the Docs site.

*How often are changes implemented?*

  At most twice per year. Updates to the standard are planned for February and August. When the CRB
  approves a change, the corresponding updates to the data standard and translator code are released
  roughly three months after the meeting — by February 1 and August 1 for the October and April
  meetings, respectively — giving the community time to update their translators as needed. See
  :ref:`versioning` and :ref:`change-review-board`.
