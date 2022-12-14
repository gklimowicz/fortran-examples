# fortran-examples
This repository contains submodule for many Fortran-related projects.

The goal is to have a wide variety of example Fortran codes, as found in the wild.

I plan to use them as inputs into analyses on Fortran source, such as the investigation for the proposed Fortran 202y preprocessor.

# So. Many. Directories.
The numerous top-level directories each contain the source code for a single "project".
For projects that come from a super-project (such as "GEOS-ESM"), all the directories begin with the super-project's name as a prefix (such as "GEOS-ESM-AeroApps).

# Metrics
File names that begin with "all-" contain aggregated lists and metrics for the sample projects.

| File                 | Contents                                            |
|:---------------------|:----------------------------------------------------|
| all-files            | Every file in every project (including non-Fortran) |
| all-fortran-files    | Every Fortran file in every project                 |
| all-fortran-files-lc | Get `wc -l` output for every Fortran file           |
| all-projects         | List of all projects                                |
| all-projects-lc      | Get `wc -l` for Fortran files in each project       |


# `provenance`
Originally, this was to contain the provenance of each of the repositories.
Since this is already saved in `.gitmodules` for the git-based projects,
this now only contains the sources for the non-git projects.

This file is managed manually for the time being.


# Makefile
The Makefile drives the collection of these repositories. Most of them come from the projects listed on https://fortran-lang.org.

## `make all`
Create project lists and metrics as described above.

## `make add-new-projects`
Search for new projects to add. Only looks at fortran-lang.org so far.

## `make update-existing-projects`
Does a `git pull` on each submodule to bring it up to date.

## `make stats`
Print a couple lines of interesting statistics about the repository as a whole.

<!--  LocalWords:  GEOS ESM AeroApps
 -->
