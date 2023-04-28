# fortran-examples
This repository contains submodules for many Fortran-related projects.

The goal is to have a wide variety of example Fortran codes, as found in the wild.

We plan to use them as inputs into analyses on Fortran source, such as the investigation for the proposed Fortran 202y preprocessor.


# Usage
To get all these example Fortran projects
```
git clone git@github.com:gklimowicz/fortran-examples.git
git submodule update --init
```

After what seems like forever, you should have some 12GB or so of Fortran example projects.

This works best on a case-sensitive file system.
Some projects seem to have a few files with filenames that clash on case-insensitive file systems.

To update after changes have been pushed to `fortran-examples`
```
git pull
git submodule update
```


# So. Many. Directories.
The numerous top-level directories each contain the source code for a single "project".
For projects that come from a super-project (such as "GEOS-ESM"), all the directories begin with the super-project's name as a prefix (such as "GEOS-ESM-AeroApps).

Although this can obscure the `Makefile` and `all-*` file lists in the top-level directory,
it felt silly to create a `src` directory and put all the projects underneath that.


# Metrics
File names that begin with "all-" contain aggregated lists and metrics for the sample projects.

| File                 | Contents                                            |
|:---------------------|:----------------------------------------------------|
| all-files            | Every file in every project (including non-Fortran) |
| all-fortran-files    | Every Fortran file in every project                 |
| all-fortran-files-lc | Get `wc -l` output for every Fortran file           |
| all-projects         | List of all projects                                |
| all-projects-lc      | Get `wc -l` for Fortran files in each project       |


# `origins`
Originally, this was to contain the provenance of each of the repositories.
Since this is already saved in `.gitmodules` for the git-based projects,
this now only contains the origin for sources for non-git projects.
(Note that this may not be a direct link to a `.tgz` or `.zip` file, but
may be the page where such a link can be found.)

This file is managed manually.


# Makefile
The Makefile drives the collection of these repositories. Most of them come from the projects listed on https://fortran-lang.org.

## `make all`
Create project lists and metrics as described above.

## `make add-new-projects`
Search for new projects to add.
Only looks at fortran-lang.org so far.

This is broken at the moment.
To be fixed soon, as fortran-lang.org is a great source for Fortran projects that are actually used.


## `make update-existing-projects`
Do a `git pull` on each project directory here to bring it up to date
with respect to its home repository.


## `make stats`
Print a couple lines of interesting statistics about the repository as a whole.

<!--  LocalWords:  GEOS ESM AeroApps
 -->
