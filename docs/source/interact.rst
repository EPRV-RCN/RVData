

.. |missing| replace:: **TBD**

Repository structure and how to interact with it
************************************************

This repo follows the Gitflow workflow: https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow

Branch Structure
=================
* *main*: Release branch from which tagged releases are generated.
* *develop*: Development branch where updates are aggregated between releases
* *feature_branch_name*: Feature branches should be forked off of develop, and should be named with a human readable intuitive name. Delete feature branches once merged into develop and work in them is complete.

Code Standards
==============

Code in this repository must be compatable with the `flake8 <https://flake8.pycqa.org/en/latest/>` standard and `Google style <https://flake8.pycqa.org/en/latest/>` docstrings.


Developer Overview
==================

We recommend that developers use the following workflow:

**First, setup your working environment**

#. Create a fork of the repository into your github account by clicking the Fork icon on the right top corner of the main github EPRV-RCN/RVdata repo page.

#. Clone your fork to your local computer:
   ``git clone git@github.com:YourGithubID/RVdata.git``

#. Setup a new remote named upstream that points to the project level repository:
   ``git remote add upstream git@github.com:EPRV-RCN/RVdata.git``


**To add a new feature or bugfix to the repository:**

#. Ensure that you are on the develop branch:
   ``git checkout develop``

#. Create a feature branch off of develop for your new work. Bugfix branches should prefix with bugfix.
   ``git checkout -b feature_branch_name``

#. Make your changes, commit them with a useful commit message, and push to your fork:
   ``git commit -m "Description of committed changes``
   ``git push origin feature_branch_name``

#. Setup a pull request from your feature branch into the project level develop branch. Use the Pull requests menu item on the top bar of github.com. Be sure that the pull request points into the **EPRV-RCN/RVdata develop** branch. This will trigger a review request. Two reviews from the core development team are required before merging.

#. Iterate with reviewers as needed, using the pull request to capture discussion comments. Once the reviewers approve the pull request, the original author is responsible for merging. **If you encounter a merge conflict, ask for help. Do not proceed!**


Testing & Integration
=====================

The repository uses continuous integration and GitHub actions to validate all code changes at the pull request stage. These changes must complete prior to merging.

Pull requests require two reviews from the reviewers listed on the :ref:`code-review-team`

CI tests require data files to validate against the data level
definitions. For instrument translators, these files should be native
data files for the instrument, which get translated and then validated
against the definition. For instruments which natively generate data
in the standard format, these files should be representative files at
each data level and will be validated against future base standard
code changes.  These files are stored on central servers - contact the
:ref:`code-review-team` to load new or updated files.
