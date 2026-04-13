Governance and Version Control
******************************

.. _versioning:

Version Control and Release Schedule
====================================

Versioning of the Data Standard and versioning of the Translator Code
are maintained independently.

The Data Standard uses the `Calendar Versioning <https://calver.org/>`_
scheme of YYYY.MM.  Updates to the standard will happen no more than
twice per year and are planned for the months of February and
August.  The first release version is *2026.02*.

The Translator Code uses the `Semantic Versioning <https://semver.org/>`_
scheme of MAJOR.MINOR.PATCH. There is no
prescribed release schedule, although it is anticipated that minor
releases will be no more frequent than twice per year, also in
Feburary and August.



.. _change-review-board:

Change Review Board
===================

We expect this to be a living standard and have thus established a
Change Review Board (CRB) that includes representatives from some of
the major stakeholders in the EPRV community including NASA, ESO,
KPNO, Gemini, Keck, DACE & NExScI. The CRB will meet twice a year, in
late April (supporting an August release) and late October (supporting
a February release), to review suggested updates to the Data Standard
and determine if a new version should be released.


Suggestions for updates to the data standard should be submitted as a
Feature Request issue on the project’s github repo and labeled as
‘datastandard’. Suggestions made at least one month before a CRB
meeting will be discussed at that meeting and the board may reach out
to the suggesting author for additional input or to iterate on the
idea. Suggestions made less than one month of a CRB meeting may be
rolled to a future meeting.

If an update is made to the data standard then an announcement will be
posted via the github readme file and all changes will be detailed on
the project’s read the docs page. A corresponding update to the data
standard and translator code bases will be released ~3 months later
(by August 1 and February 1 for the April and October CRB meetings,
respectively) to allow the community time to update their translators
as necessary.

Changes to the underlying data standard code base and the data
translator code base will be reviewed by a separate Data Standard
Definition. However as changes to the data standard approved by the
board will likely result in changes to these code bases and so these
teams may work in tandem.


.. _code-review-team:

Code Review Team
================


The Code Review Team are the code developers responsible for
maintaining the translator codebase, and reviewing and implementing
proposed code changes.

This team is comprised of Jennifer Burt, BJ Fulton, and Chad Bender.
