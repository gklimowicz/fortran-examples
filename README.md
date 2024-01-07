# fortran-examples
This repository contains submodules for many Fortran-related projects.

The goal is to have a wide variety of example Fortran codes, as found in the wild.

We plan to use them to analyze Fortran source as found in nature to investigate proposed features for the Fortran 202y (e.g., the proposed preprocessor).


# Usage
To get all these example Fortran projects
```
git clone git@github.com:gklimowicz/fortran-examples.git
git submodule update --init
```

After what seems like forever, you should have some 45GB or so of Fortran example projects.

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

Although this can obscure the `Makefile` and `all-*` file lists in the top-level directory, it felt silly to create a `src` directory and put all the projects underneath that.


# `origins.txt`
Originally, this was to contain the provenance of each of the repositories.
Since this is already saved in `.gitmodules` for the git-based projects,
this now only contains the origin for sources for non-git projects.
(Note that this may not be a direct link to a `.tgz` or `.zip` file, but
may be the page where such a link can be found.)

This file is managed manually.


# `project-exceptions.txt`
There are projects mentioned on the `fortran-lang.org` and `Beliavsky` pages that don't actually contain Fortran source code that we are interested in (such as only contain `.fypp` files).
We ignore projects that are in the `project-exceptions.txt` list.


# `fortran-exceptions.txt`
There are files in projects that have names like Fortran files,
but don't actually contain Fortran source code (e.g., `Makefile.f90`).
We ignore files that are in the `file-exceptions.txt` list.


# `fortran-file-patterns.txt`
Of all the files found in the projects, files that match these `grep` patterns
are likely to be Fortran files.
(Those that aren't really Fortran should be in `fortran-exceptions.txt`.


# Makefile
The Makefile drives the collection of these repositories and the gathering of statistics.
Most of them come from the projects listed on [fortran-lang.org's projects pages](https://fortran-lang.org) and [Beliavsky's Fortran-code-on-GitHub list](https://github.com/Beliavsky/Fortran-code-on-GitHub).


## `make all`
Update Git submodules from remote repos, and create project lists and metrics as described above.


## `make update`
Update each Git project from its origin repository.
Runs `git submodule update --remote`.


## `make add-new-projects`
Search for new projects to add.
Only looks at [fortran-lang.org](https://fortran-lang.org/en/packages) and
[github.com/Beliavsky](https://github.com/Beliavsky/Fortran-code-on-GitHub) so far.


## `make update-existing-projects`
Do a `git pull` on each project directory here to bring it up to date
with respect to its home repository.


## `make stats.txt`
Print a couple lines of interesting statistics about the repository as a whole.

# Metrics
File names that begin with "all-" contain aggregated lists and metrics for the sample projects.

| File                        | Contents                                            |
|:----------------------------|:----------------------------------------------------|
| all-files.txt               | Every file in every project (including non-Fortran) |
| all-fortran-files.txt       | Every Fortran file in every project                 |
| all-fortran-files-attr.txt  | Attributes found in each Fortran file               |
| all-fortran-files-fixed.txt | Every fixed-form Fortran file in every project      |
| all-fortran-files-free.txt  | Every free-form Fortran file in every project       |
| all-fortran-files-lc.txt    | Line count for every Fortran file                   |
| all-projects.txt            | List of all projects                                |
| all-projects-lc.txt         | Line count for Fortran files in each project        |


## `all-fortran-files-attr.txt`
Each Fortran file is scanned to identify simple characteristics of the file.

| Attribute        | Description                                                         |
|:-----------------|:--------------------------------------------------------------------|
| form             | `fixed` if the file appears to be fixed-form Fortran                |
|                  | `free` if the file appears to be free-form Fortran                  |
| lines            | The number of lines in the file                                     |
| maxlinelength    | The maximum length of a line in the file                            |
| cpreprocessor    | The number of C preprocessor directives in the file                 |
| ccomment         | The number of fixed-form comments beginning with `C` in column 1    |
| dcomment         | The number of fixed-form comments beginning with `D` in column 1    |
| starcomment      | The number of fixed-form comments beginning with `*` in column 1    |
| fixedbang        | The number of times `!` comments appear in a fixed-form file        |
| continuations    | The total number of line continuations in the file                  |
| maxcontinuations | The longest sequence of continuation lines in the file              |
| text73           | The number of lines in fixed-form file with text in columns 73-80   |
| text133          | The number of lines in fixed-form file with text in columns 133-140 |
| ampcont          | The number of `&`-style continuations seen in a free-form file.     |
| ampampcont       | The number of `&`-style continuations with `&` on continued line.   |
| include          | The number of Fortran `INCLUDE` lines in the file.                  |
| openmpdir        | The number of OpenMP directives in the file.                        |
| openaccdir       | The number of OpenACC directives in the file.                       |
| otherdir         | The number of unidentifiable directives in the file.                |

## `duplicates-ok.txt`
This is a list of files (produced by `bin/find-duplicates`) that are not really duplicated projects. They are either completely different projects, modernizations, or some other variation on a theme.

<!--  LocalWords:  GEOS ESM AeroApps modernizations
 -->
