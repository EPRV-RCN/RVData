

.. |missing| replace:: **TBD**

Repository structure and how to interact with it
************************************************

This repo follows the Gitflow workflow: https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow

Branch Structure
=================
* *main*: Release branch from which tagged releases are generated.
* *develop*: Development branch where updates are aggregated between releases
* *feature_branch_name*: Feature branches should be forked off of develop, and should be named with a human readable intuitive name. Delete feature branches once merged into develop and work in them is complete.


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

[In progress]


Noteable Extensions & Header Keywords
=====================================

**The Receipt Extension**

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

**The Config Extension**

This extension contains a two-column dataframe that stores each line of the input configuration file.

This file is based on the configuration language that the Python ConfigPaser class implements for Python programs. 
The syntax of the configuration language primarily consists of ``key = value`` statements and ``# comments`` (or ``; comments``).

The creation of the key/value pairs, in general, follows the convention as those provided by the ConfigParser classs. 
As used in the KPF and NEID Pipelines, the keys consists of alphanumerics and underscores (e.g. KPFPIPE_TEST_DATA) while 
the value is in the form of a string, a number, a boolean, or a list.



**Keywords**

*MULTIPLICITY* : Keyword in the PRIMARY extension that indicates whether or not an extension can be duplicated to have multiple HDUs based on, e.g., the number of traces.

*REQUIRED* : Keyword in the PRIMARY extension that indicates whether or not an extension must be present and meaningfully populated in the FITS file in order for the file to be compliant with the EPRV Data Standard.

*INSTERA* : Tag for the "instrument era" used to track permanent changes to instrument, e.g. 1.1.2
    Incrementing: Major changes that result in RV offset are recorded in the first numeral space. Small but permanent changes, (changes in the calibration light sources, etc) that don’t necessarily introduce RV offsets but do represent an ongoing status change in the instrument, are recorded in the second numeral. Anything that could impact the calibration or science outputs goes here
    Even smaller changes, or things that you only recognize retrospectively, are recorded in the third numeral
    Each team will place a csv file into the config directory with 1 line per era that has the start and end dates of each timespan. If this file is already maintained elsewhere by the instrument team then link to that instead so the file doesn’t need to be maintained in multiple places
