[![Documentation Status](https://readthedocs.org/projects/eprv-data-standard/badge/?version=latest)](https://eprv-data-standard.readthedocs.io/en/latest/?badge=latest)

# EPRV_DataStandard
Development of a community endorsed, standardized EPRV data format at the 1D extracted spectra level, and a set of translator tools

https://eprv-data-standard.readthedocs.io/en/latest/

## Objective

1) Establish a community recommendation for standardized RV data and telemetry formats based upon input from participating (E)PRV instrument teams and partners from the observatories and data archives that will host the data for each instrument

2) Develop prototype ‘translator’ tools that will re-write data from the native spectrograph outputs into this community standard format, and act as a template for other instruments to do the same.

## Current Data Standard Links

Note: These continue to be in development for the time being, but the links below capture the current data format as recommended by this group

Data Format Overview: [Google Slides Here](https://docs.google.com/presentation/d/1XLTaW4iWFUQiw2KesEsYpkgGDLyoTWRG5QBszSWgxRM/edit?usp=sharing)

Level 2 Header Keywords: [Google Doc Here](https://docs.google.com/spreadsheets/d/1jv40V6z0DQEPOsw4wZUHFUflI8z777dHpCWQLGGWUO4/edit?usp=sharing)


## Repository structure and how to interact with it

This repo follows the [Gitflow workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow)

### Branch Structure

- *main* - Release branch from which tagged releases are generated.
- *develop* - Development branch where updates are aggregated between releases
- *feature_branch_name* - Feature branches should be forked off of develop, and should be named with a human readable intuitive name.  Delete feature branches once merged into develop and work in them is complete.



### Developers
We recommend that developers use the following workflow:

**First, setup your working environment:**

1) Create a fork of the repository into your github account by clicking the *Fork* icon on the right top corner of the main github EPRV-RCN/RVdata repo page.

2) Clone your fork to your local computer:
```git clone git@github.com:YourGithubID/RVdata.git```

3) Setup a new *remote* named **upstream** that points to the project level repository:
```git remote add upstream git@github.com:EPRV-RCN/RVdata.git```

**To add a new feature or bugfix to the repository:**

1) Create a feature branch off of develop for your new work. Bugfix branches should prefix with *bugfix*.
```git checkout develop
git checkout -b feature_branch_name
```

2) Make your changes, commit them with a useful commit message, and push to your fork.
```git add new_or_updated_filenames.py
git commit -m "Description of committed changes"
git push origin feature_branch_name
```

3) Setup a pull request from your feature branch into the project level develop branch. Use the *Pull requests* menu item on the top bar of github.com. Be sure that the pull request points into the **EPRV-RCN/RVdata _develop_** branch. This will trigger a review request. Two reviews from the core development team are required before merging.

4) Iterate with reviewers as needed, using the pull request to capture discussion comments. Once the reviewers approve the pull request, the original author is responsible for merging. **If you encounter a merge conflict, ask for help. _Do not proceed!_**


