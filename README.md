# fortran-examples

This repository cantains submodule for many Fortran-related projects.

The goal is to have a wide variety of example Fortran codes, as found in the wild.

I plan to use them as inputs into analyses on Fortran source, such as the investigation for the proposed Fortran 202y preprocessor.

# Makefile
The Makefile drives the collection of these repositories. Most of them come from the projects listed on https://fortran-lang.org.

## `make all`
This creates lists of
* all projects (just their directory names)
* all the Fortran files across all the projects
* a word count for each Fortran file

## `make add-new-projects`
Search for new projects to add. Only looks at fortran-lang.org so far.

## `make update-existing-projects`
Does a `git pull` on each submodule te bring it up to date.

## `make stats`
Print a couple lines of interesting statistics about the repository as a whole.
