dnl configure.ac -- Main configuration script for WHIZARD
dnl
dnl Process this file with autoconf to produce a configure script.
dnl ************************************************************************
dnl configure.ac -- Main configuration script for WHIZARD
dnl configure.ac -- WHIZARD configuration
dnl
dnl Copyright (C) 1999-2022 by
dnl     Wolfgang Kilian <kilian@physik.uni-siegen.de>
dnl     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
dnl     Juergen Reuter <juergen.reuter@desy.de>
dnl     with contributions from
dnl     cf. main AUTHORS file
dnl
dnl WHIZARD is free software; you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.
dnl
dnl WHIZARD is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
dnl
dnl ***********************************************************************
dnl Environment variables that can be set by the user:
dnl   FC		Fortran compiler
dnl   FCFLAGS		Fortran compiler flags
dnl ***********************************************************************
dnl
dnl Start configuration
AC_INIT([WHIZARD],[3.1.0])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.12.2 color-tests parallel-tests])
AC_PREREQ([2.65])
AM_MAKE_INCLUDE

dnl Make make less verbose to improve signal/noise
AM_SILENT_RULES([yes])

########################################################################
### Package-specific initialization

AC_MSG_NOTICE([**************************************************************])
WO_CONFIGURE_SECTION([Start of package configuration])

### Further version information
PACKAGE_DATE="Dec 14 2022"
PACKAGE_STATUS="release"

AC_SUBST(PACKAGE_DATE)
AC_SUBST(PACKAGE_STATUS)

AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([Package name: AC_PACKAGE_NAME()])
AC_MSG_NOTICE([Version:      AC_PACKAGE_VERSION()])
AC_MSG_NOTICE([Date:         $PACKAGE_DATE])
AC_MSG_NOTICE([Status:       $PACKAGE_STATUS])
AC_MSG_NOTICE([**************************************************************])

########################################################################
###---------------------------------------------------------------------
### shared library versioning (not the same as the package version!)

LIBRARY_VERSION="-version-info 2:2:0"
AC_SUBST([LIBRARY_VERSION])

########################################################################
###---------------------------------------------------------------------
### Define the main package variables

### Source directory, for testing purposes
SRCDIR=`cd $srcdir && pwd`
AC_SUBST([SRCDIR])

### Build directory, for testing purposes
BUILDDIR=`pwd`
AC_SUBST([BUILDDIR])

### Location of installed libraries and such
eval BINDIR=$bindir
case $BINDIR in
NONE*) eval BINDIR=$prefix/bin ;;
esac
case $BINDIR in
NONE*) BINDIR="\${prefix}/bin" ;;
esac
AC_SUBST([BINDIR])

eval INCLUDEDIR=$includedir
case $INCLUDEDIR in
NONE*) eval INCLUDEDIR=$prefix/include ;;
esac
case $INCLUDEDIR in
NONE*) INCLUDEDIR="\${prefix}/include" ;;
esac
AC_SUBST([INCLUDEDIR])

eval LIBDIR=$libdir
case $LIBDIR in
NONE*) eval LIBDIR=$prefix/lib ;;
esac
case $LIBDIR in
NONE*) eval LIBDIR=$ac_default_prefix/lib ;;
esac
AC_SUBST([LIBDIR])

### Location of installed libraries and such
eval PKGLIBDIR=$libdir/$PACKAGE
case $PKGLIBDIR in
NONE*) eval PKGLIBDIR=$prefix/lib/$PACKAGE ;;
esac
case $PKGLIBDIR in
NONE*) PKGLIBDIR="\${prefix}/lib/$PACKAGE" ;;
esac
AC_SUBST([PKGLIBDIR])

### Location of installed system-independent data
eval PKGDATADIR=$datarootdir/$PACKAGE
case $PKGDATADIR in
NONE*) eval PKGDATADIR=$prefix/share/$PACKAGE ;;
esac
case $PKGDATADIR in
NONE*) PKGDATADIR="\${prefix}/share/$PACKAGE" ;;
esac
AC_SUBST([PKGDATADIR])

### Location of installed TeX files and such
eval PKGTEXDIR=$datarootdir/texmf/$PACKAGE
case $PKGTEXDIR in
NONE*) eval PKGTEXDIR=$prefix/share/texmf/$PACKAGE ;;
esac
case $PKGTEXDIR in
NONE*) PKGTEXDIR="\${prefix}/share/texmf/$PACKAGE" ;;
esac
AC_SUBST([PKGTEXDIR])

### Parent location of installed .mod files
### To be used in Fortran source
FMODDIR=$prefix/lib/mod
case $FMODDIR in
NONE*) FMODDIR="\${prefix}/lib/mod" ;;
esac
AC_SUBST([FMODDIR])

### To be used in Makefile.am
### Don't use ${libdir} since lib may be changed to lib64 by configure
fmoddir="\${prefix}/lib/mod"
AC_SUBST([fmoddir])

########################################################################
###---------------------------------------------------------------------
### Required programs and checks
### GNU Tools

WO_CONFIGURE_SECTION([Generic tools])

### Initialize LIBTOOL
LT_INIT(dlopen)

LT_PREREQ([2.4.1b])

AX_CHECK_GNU_MAKE()
AC_PROG_GREP()

AC_MSG_CHECKING([for the suffix of shared libraries])
case $host in
  *-darwin* | rhapsody*)
    SHRLIB_EXT="dylib"
    ;;
  cygwin* | mingw* | pw32* | cegcc* | os2*)
    SHRLIB_EXT="dll"
    ;;
  hpux9* | hpux10* | hpux11*)
    SHRLIB_EXT="sl"
    ;;
  *)
    SHRLIB_EXT="so"
    ;;
esac
if test "x$SHRLIB_EXT" != "x"; then
   SHRLIB_EXT=$SHRLIB_EXT
else
   SHRLIB_EXT="so"
fi
AC_MSG_RESULT([.$SHRLIB_EXT])
AC_SUBST(SHRLIB_EXT)

### Export whether the C compiler is GNU
AC_MSG_CHECKING([whether the C compiler is the GNU compiler])
if test "x$ac_cv_c_compiler_gnu" = "xyes"; then
   CC_IS_GNU=".true."
else
   CC_IS_GNU=".false."
fi
AC_MSG_RESULT([$ac_cv_c_compiler_gnu])
AC_SUBST([CC_IS_GNU])

AC_CHECK_HEADERS([quadmath.h])
if test "x$ac_cv_header_quadmath_h" = "xyes"; then
   CC_HAS_QUADMATH=".true."
else
   CC_HAS_QUADMATH=".false."
fi
AC_SUBST([CC_HAS_QUADMATH])

########################################################################
###---------------------------------------------------------------------
### Host system MAC OS X check for XCode

case $host in
     *-darwin*)
        WO_HLINE()
     	AC_MSG_NOTICE([Host is $host, checking for XCode])
        AC_PATH_PROG(XCODE_SELECT, xcode-select)
	# locate currently selected Xcode path
	if test "x$XCODE_SELECT" != "x"; then
	  AC_MSG_CHECKING(Xcode location)
	  DEVELOPER_DIR=`$XCODE_SELECT -print-path`
	  AC_MSG_RESULT([$DEVELOPER_DIR])
	else
	  DEVELOPER_DIR=/Developer
	fi
	AC_SUBST(DEVELOPER_DIR)

	XCODEPLIST=$DEVELOPER_DIR/Applications/Xcode.app/Contents/version.plist
	if test -r "$XCODEPLIST"; then
	  AC_MSG_CHECKING(Xcode version)
	  if test "x$DEFAULTS" != "x"; then
	    XCODE_VERSION=`$DEFAULTS read $DEVELOPER_DIR/Applications/Xcode.app/Contents/version CFBundleShortVersionString`
	  else
	    XCODE_VERSION=`tr -d '\r\n' < $XCODEPLIST | sed -e 's/.*<key>CFBundleShortVersionString<\/key>.<string>\([[0-9.]]*\)<\/string>.*/\1/'`
	  fi
	  AC_MSG_RESULT([$XCODE_VERSION])
	  AC_SUBST(XCODE_VERSION)
	fi

	AC_MSG_NOTICE([checking for Security Integrity Protocol (SIP)])
	AC_PATH_PROG(CSRUTIL, csrutil)
	if test "x$CSRUTIL" != "x"; then
	   SIP_CHECK=`$CSRUTIL status | $SED "s/System Integrity Protection status: //"`
	   if test "$SIP_CHECK" = "enabled."; then
	      SIP_ACTIVE="yes"
	   else
	      SIP_ACTIVE="no"
 	   fi
	else
	   SIP_ACTIVE="no"
	fi
	AC_MSG_CHECKING([Checking whether MAC OS X SIP is activated])
	AC_MSG_RESULT([$SIP_ACTIVE])
	AC_SUBST([SIP_ACTIVE])

	WO_HLINE()
	;;
     *)
        ;;
     esac

########################################################################
###---------------------------------------------------------------------
### Enable the distribution tools
### (default: disabled, to speed up compilation)

AC_ARG_ENABLE([distribution],
  [AS_HELP_STRING([--enable-distribution],
    [build the distribution incl. all docu (developers only) [[no]]])])

AC_CACHE_CHECK([whether we want to build the distribution],
[wo_cv_distribution],
[dnl
if test "$enable_distribution" = "yes"; then
  wo_cv_distribution=yes
else
  wo_cv_distribution=no
fi])

AM_CONDITIONAL([DISTRIBUTION],
  [test "$enable_distribution" = "yes"])

########################################################################
###---------------------------------------------------------------------

if test "$enable_shared" = no; then
  AC_MSG_ERROR([you've used --disable-shared which will not produce a working Whizard.])
fi

########################################################################
###---------------------------------------------------------------------
### We include the m4 macro tool here

AC_PATH_PROG(M4,m4,false)
if test "$M4" = false; then
   AM_CONDITIONAL([M4_AVAILABLE],[false])
else
   AM_CONDITIONAL([M4_AVAILABLE],[true])
fi

########################################################################
###---------------------------------------------------------------------
### Dynamic runtime linking

WO_CONFIGURE_SECTION([Dynamic runtime linking])

### Look for libdl (should provide 'dlopen' and friends)
AC_PROG_CC()
WO_PROG_DL()

### Define the conditional for static builds

if test "$enable_static" = yes; then
   AM_CONDITIONAL([STATIC_AVAILABLE],[true])
else
   AM_CONDITIONAL([STATIC_AVAILABLE],[false])
fi

########################################################################
###---------------------------------------------------------------------
### Noweb

WO_CONFIGURE_SECTION([Checks for 'noweb' system])

### Enable/disable noweb and determine locations of notangle, cpif, noweave
WO_PROG_NOWEB()

########################################################################
###---------------------------------------------------------------------
### LaTeX

WO_CONFIGURE_SECTION([Checks for 'LaTeX' system])

### Determine whether LaTeX is present
AC_PROG_LATEX()
AC_PROG_DVIPS()
AC_PROG_PDFLATEX()
AC_PROG_MAKEINDEX()
AC_PROG_PS2PDF()
AC_PROG_EPSPDF()
AC_PROG_EPSTOPDF()
if test "$EPSPDF" = "no" -a "$EPSTOPDF" = "no"; then
   AC_MSG_NOTICE([*********************************************************])
   AC_MSG_NOTICE([WARNING: eps(to)pdf n/a; O'Mega documentation will crash!])
   AC_MSG_NOTICE([WARNING:  this applies only to the svn developer version!])
   AC_MSG_NOTICE([*********************************************************])
fi
AC_PROG_SUPP_PDF()
AC_PROG_GZIP()
AC_PATH_PROG(ACROREAD,acroread,false)
AC_PATH_PROG(GHOSTVIEW,gv ghostview,false)
AC_PROG_DOT()

### Determine whether Metapost is present and whether event display is possible
AC_PROG_MPOST()
WO_CHECK_EVENT_ANALYSIS_METHODS()

### We put here the check for HEVEA components as well
WO_PROG_HEVEA()

########################################################################
###---------------------------------------------------------------------
### Fortran compiler

WO_CONFIGURE_SECTION([Fortran compiler checks])

### Determine default compiler to use
user_FCFLAGS="${FCFLAGS}"
AC_PROG_FC()

### Choose FC standard for PYTHIA6 F77 files
AC_PROG_F77([$FC])

### Determine compiler vendor and version
WO_FC_GET_VENDOR_AND_VERSION()

### Veto against old gfortran 4 versions
WO_FC_VETO_GFORTRAN_4()

### Veto against buggy gfortran 6.5.0 version
WO_FC_VETO_GFORTRAN_65()

### Veto against ifort 15/16/17/18
WO_FC_VETO_IFORT_15_18()

### Veto against ifort 19.0.0/1/2
WO_FC_VETO_IFORT_190012()

### Veto against ifort 21.0/1/2
WO_FC_VETO_IFORT_21012()

### Require extension '.f90' for all compiler checks
AC_FC_SRCEXT([f90])

### Determine flags and extensions
WO_FC_PARAMETERS()

### Determine flags for linking the Fortran runtime library
WO_FC_LIBRARY_LDFLAGS()

### Check for Fortran 95 features
WO_FC_CHECK_F95()

### Check for allocatable subobjects (TR15581)
WO_FC_CHECK_TR15581()

### Check for allocatable scalars
WO_FC_CHECK_ALLOCATABLE_SCALARS()

### Check for ISO C binding support
WO_FC_CHECK_C_BINDING()

### Check for procedures pointers and abstract interfaces
WO_FC_CHECK_PROCEDURE_POINTERS()

### Check for type extension and further OO features
WO_FC_CHECK_OO_FEATURES()

### Check for submodules (not yet used)
WO_FC_CHECK_TR19767()

### Check for F2003 command-line interface
WO_FC_CHECK_CMDLINE()

### Check for F2003-style access to environment variables
WO_FC_CHECK_ENVVAR()

### Check for the flush statement
WO_FC_CHECK_FLUSH()

### Check for iso_fortran_env
WO_FC_CHECK_ISO_FORTRAN_ENV()
WO_FC_CHECK_ISO_FORTRAN_ENV_2008()

### Turn on/off master switch for debugging features
WO_FC_SET_DEBUG()

### OpenMP threading activated upon request
AC_OPENMP()
WO_FC_SET_OPENMP()

### Profiling compilation enforced upon request
WO_FC_SET_PROFILING()

### Impure subroutines enforced upon request
WO_FC_SET_OMEGA_IMPURE()

### Find the extension of Fortran module files
WO_FC_MODULE_FILE([FC_MODULE_NAME], [FCMOD], [$FC], [f90])

###---------------------------------------------------------------------
### Check for the requested precision

WO_FC_CONFIGURE_KINDS([src/basics/kinds.f90])
AC_PROG_INSTALL()
${INSTALL} -d circe1/src
cp -a src/basics/kinds.f90 circe1/src
${INSTALL} -d circe2/src
cp -a src/basics/kinds.f90 circe2/src
${INSTALL} -d omega/src
cp -a src/basics/kinds.f90 omega/src
${INSTALL} -d vamp/src
cp -a src/basics/kinds.f90 vamp/src

########################################################################
# VAMP Fortran options for the configure script
########################################################################

WO_FC_SET_MPI()


########################################################################
###---------------------------------------------------------------------
### OCaml

WO_CONFIGURE_SECTION([Objective Caml checks])

### Check for ocamlc and its relatives
AC_PROG_OCAML()
if test "$enable_ocaml" != "no"; then
   AC_OCAML_VERSION_CHECK(405000)
   AC_PROG_OCAMLLEX()
   AC_PROG_OCAMLYACC()
   AC_PROG_OCAMLCP()
   AC_OCAML_BIGARRAY_MODULE()
   ### Ocamlweb is required to be newer than v0.9
   AC_PROG_OCAMLWEB(009000)
   AC_PROG_OCAML_LABLGTK()
   AC_PATH_PROGS([OCAMLDOT],[ocamldot],[no])
   AM_CONDITIONAL([OCAMLDOT_AVAILABLE],[test "$OCAMLDOT" != "no"])
   AC_PATH_PROGS([OCAMLDEP],[ocamldep],[no])
   AM_CONDITIONAL([OCAMLDEP_AVAILABLE],[test "$OCAMLDEP" != "no"])
   AC_PATH_PROGS([OCAMLDEFUN],[ocamldefun],[no])
else
   AC_MSG_NOTICE([WARNING: OCaml and O'Mega matrix elements disabled by request!])
   AM_CONDITIONAL([OCAMLWEB_AVAILABLE],[false])
   AM_CONDITIONAL([OCAMLDOT_AVAILABLE],[false])
   AM_CONDITIONAL([OCAMLDEP_AVAILABLE],[false])
fi

########################################################################
###---------------------------------------------------------------------
### C++

WO_CONFIGURE_SECTION([C++ compiler checks])
AC_PROG_CXX()
AC_CXX_LIBRARY_LDFLAGS()

########################################################################
###---------------------------------------------------------------------
### Checks for external interfaces

WO_CONFIGURE_SECTION([Checking for PYTHON / PYTHON API])

AX_PYTHON_DEVEL([>= '3.5'])
WO_PROG_PYTHON_API()

########################################################################
# O'Mega options for the configure script
########################################################################

########################################################################
###---------------------------------------------------------------------
### O'Mega UFO file paths

WO_CONFIGURE_SECTION([O'Mega UFO file paths])

AC_ARG_ENABLE([default-UFO-dir],
[  --enable-default-UFO-dir=directory
                          Read precomputed model tables from this directory,
                          which will be populated by an administrator at
                          install time [[default=$datadir/UFO, enabled]].],
[case "$enableval" in
   no) OMEGA_DEFAULT_UFO_DIR="."
       ;;
    *) OMEGA_DEFAULT_UFO_DIR="$enableval"
       ;;
esac],
[### use eval b/c $datadir defaults to unexpanded ${datarootdir}
case "$OMEGA_DEFAULT_UFO_DIR" in
   "") OMEGA_DEFAULT_UFO_DIR="${prefix}/omega/share/UFO"
      ;;
   *) eval OMEGA_DEFAULT_UFO_DIR="$datadir/UFO"
      ;;
esac])
AC_SUBST([OMEGA_DEFAULT_UFO_DIR])

case "$OMEGA_DEFAULT_UFO_DIR" in
  .|""|NONE*) OMEGA_DEFAULT_UFO_DIR="."
              ;;
           *) AC_MSG_NOTICE([Creating default UFO directory $OMEGA_DEFAULT_UFO_DIR])
              $MKDIR_P "$OMEGA_DEFAULT_UFO_DIR" 2>/dev/null
              chmod u+w "$OMEGA_DEFAULT_UFO_DIR" 2>/dev/null
              ;;
esac

###---------------------------------------------------------------------
### Recola

WO_CONFIGURE_SECTION([RECOLA])

WO_PROG_RECOLA()


########################################################################
###---------------------------------------------------------------------
### Libraries

###---------------------------------------------------------------------
### LHAPDF

WO_CONFIGURE_SECTION([LHAPDF])

WO_PROG_LHAPDF()

###---------------------------------------------------------------------
### ROOT

WO_CONFIGURE_SECTION([ROOT])

WO_ROOT_PATH(,[
        AC_DEFINE([HAVE_ROOT],,[Root library])
        AC_CHECK_LIB([dl],[dlopen],[],AC_MSG_WARN([Root libraries not linking properly]))
        ],AC_MSG_RESULT([The ROOT support of HepMC might not be working properly]))

###---------------------------------------------------------------------
### HepMC

WO_CONFIGURE_SECTION([HepMC])

WO_PROG_HEPMC()

###---------------------------------------------------------------------
### STDHEP

WO_CONFIGURE_SECTION([STDHEP])

WO_PROG_TIRPC()

AC_MSG_NOTICE([StdHEP v5.06.01 is included internally])

###---------------------------------------------------------------------
### LCIO

WO_CONFIGURE_SECTION([LCIO])

WO_PROG_LCIO()

###---------------------------------------------------------------------
### HDF5 (for events, grids etc.)

WO_CONFIGURE_SECTION([HDF5])

WO_PROG_HDF5(1.8.0,no)

###---------------------------------------------------------------------
### PYTHIA6, PYTHIA8 etc

WO_CONFIGURE_SECTION([SHOWERS PYTHIA6 PYTHIA8 MPI])

WO_PROG_QCD()

WO_PROG_PYTHIA8()

###---------------------------------------------------------------------
### HOPPET

WO_CONFIGURE_SECTION([HOPPET])

WO_PROG_HOPPET()

###---------------------------------------------------------------------
### FASTJET

WO_CONFIGURE_SECTION([FASTJET])

WO_PROG_FASTJET()

###---------------------------------------------------------------------
### GoSam

WO_CONFIGURE_SECTION([GOSAM])

WO_PROG_GOSAM()

###---------------------------------------------------------------------
### OpenLoops

WO_CONFIGURE_SECTION([OPENLOOPS])

WO_PROG_OPENLOOPS()

###---------------------------------------------------------------------
### LoopTools

WO_CONFIGURE_SECTION([LOOPTOOLS])

WO_PROG_LOOPTOOLS()


########################################################################
###---------------------------------------------------------------------
### Extra flags for helping the linker finding libraries

WO_CONFIGURE_SECTION([Handle linking with C++ libraries])

WO_PROG_STDCPP()

########################################################################
###---------------------------------------------------------------------
### Miscellaneous

WO_CONFIGURE_SECTION([Numerical checks])

### Disable irrelevant optimization for parameter files
### (default: disabled, to speed up compilation)

AC_ARG_ENABLE([optimization-for-parameter-files],
  [AS_HELP_STRING([--enable-optimization-for-parameter-files],
    [enable (useless) optimization for parameter files [[no]]])])

AC_CACHE_CHECK([whether we want optimization for parameter files],
[wo_cv_optimization_for_parfiles],
[dnl
if test "$enable_optimization_for_parameter_files" = "yes"; then
  wo_cv_optimization_for_parfiles=yes
else
  wo_cv_optimization_for_parfiles=no
fi])

AM_CONDITIONAL([OPTIMIZATION_FOR_PARFILES],
  [test "$enable_optimization_for_parameter_files" = "yes"])


########################################################################
###---------------------------------------------------------------------
### Wrapup

WO_CONFIGURE_SECTION([Finalize configuration])

### Main directory

AC_CONFIG_FILES([Makefile])

###---------------------------------------------------------------------
### Subdirectory src

AC_CONFIG_FILES([src/Makefile])

###---------------------------------------------------------------------
### Subdirectory python: WHIZARD's PYTHON/CYTHON interface

AC_CONFIG_FILES([python/Makefile])
AC_CONFIG_FILES([python/setup.py], [chmod u+x python/setup.py])
AC_CONFIG_LINKS([python/whizard_python.pyx:python/whizard_python.pyx])
AC_CONFIG_LINKS([python/cwhizard.pxd:python/cwhizard.pxd])

###---------------------------------------------------------------------
### Subdirectory src/hepmc

AC_CONFIG_FILES([src/hepmc/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/lcio

AC_CONFIG_FILES([src/lcio/Makefile])

###---------------------------------------------------------------------
### Subdirectory pythia6: WHIZARD's internal PYTHIA6 version

AC_CONFIG_FILES([pythia6/Makefile])

###---------------------------------------------------------------------
### Subdirectory tauola: WHIZARD's internal TAUOLA version

AC_CONFIG_FILES([tauola/Makefile])

###---------------------------------------------------------------------
### Subdirectory stdhep: WHIZARD's internal StdHep version

AC_CONFIG_FILES([mcfio/Makefile])
AC_CONFIG_FILES([stdhep/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/muli: multiple interactions

AC_CONFIG_FILES([src/muli/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/lhapdf5: dummy library as LHAPDF5 replacement

AC_CONFIG_FILES([src/lhapdf5/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/lhapdf: LHAPDF v6

AC_CONFIG_FILES([src/lhapdf/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/pdf_builtin: Builtin PDFs

AC_CONFIG_FILES([src/pdf_builtin/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/pdf_builtin: Electron PDFs

AC_CONFIG_FILES([src/qed_pdf/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/tauola

AC_CONFIG_FILES([src/tauola/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/xdr: XDR reader

AC_CONFIG_FILES([src/xdr/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/hoppet

AC_CONFIG_FILES([src/hoppet/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/fastjet

AC_CONFIG_FILES([src/fastjet/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/looptools

AC_CONFIG_FILES([src/looptools/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/shower: shower and all that

AC_CONFIG_FILES([src/pythia8/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/shower: shower and all that

AC_CONFIG_FILES([src/shower/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/noweb-frame: frame for whizard Noweb sources

AC_CONFIG_FILES([src/noweb-frame/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/basics: numeric kinds, strings

AC_CONFIG_FILES([src/basics/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/utilities: simple utilities

AC_CONFIG_FILES([src/utilities/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/testing: unit-test support

AC_CONFIG_FILES([src/testing/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/system: modules related to local setup and OS issues

AC_CONFIG_FILES([src/system/Makefile])

AC_CONFIG_FILES([src/system/system_dependencies.f90], [
# Be cautious: Special handling of brackets due to M4!
# 1. Configure lines containing @VARIABLE@ must be indented by exactly 7 spaces
#    to be split (^[[[:space:]]]\{7\}), the first non-space character has to
#    be an double-quote (\"). and the line has to have at least 110 characters where the
#    111th must not be a double quote, a whitespace, a slash or an ampersand ([^\" \/&]).
# 2. Appeand each 110-wide character block (\(.\{110\}[[^\"]]\), refer to a block as \1)
#    (without a trailing double-quote) with an ampersand, a literal newline character,
#    seven white-spaces and another ampersand. Repeat with remaining pattern space (be greedy).
# Note: The greedy options also allows us to parse the line beginning from each character
#       again with the search pattern.
   $SED "/^[[[:space:]]]\{7\}\".\{110,\}[[^\"]]/ s/\(.\{110\}[[^\" \/&]]\)/\1\&\\n       \&/g" \
      < src/system/system_dependencies.f90 \
      > src/system/system_dependencies.tmp
   mv -f src/system/system_dependencies.tmp src/system/system_dependencies.f90
])

AC_CONFIG_FILES([src/system/debug_master.f90])

###---------------------------------------------------------------------
### Subdirectory src/combinatorics: standard algorithms

AC_CONFIG_FILES([src/combinatorics/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/parsing: text-handling and parsing

AC_CONFIG_FILES([src/parsing/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/rng: random-number generation

AC_CONFIG_FILES([src/rng/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/expr_base: abstract expressions

AC_CONFIG_FILES([src/expr_base/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/physics: particle-physics related functions

AC_CONFIG_FILES([src/physics/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/qft: quantum (field) theory concepts as data types

AC_CONFIG_FILES([src/qft/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/types: HEP and other types for common use

AC_CONFIG_FILES([src/types/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/matrix_elements: process code and libraries

AC_CONFIG_FILES([src/matrix_elements/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/me_methods: specific process code and interface

AC_CONFIG_FILES([src/me_methods/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/particles: particle objects

AC_CONFIG_FILES([src/particles/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/beams: beams and beam structure

AC_CONFIG_FILES([src/beams/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/events: generic events and event I/O

AC_CONFIG_FILES([src/events/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/vegas: VEGAS Monte Carlo adaptive integration

AC_CONFIG_FILES([src/vegas/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/mci: multi-channel integration and event generation

AC_CONFIG_FILES([src/mci/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/phase_space: parameterization and evaluation

AC_CONFIG_FILES([src/phase_space/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/blha: BLHA support (NLO data record)

AC_CONFIG_FILES([src/blha/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/gosam: GoSAM support (NLO amplitudes)

AC_CONFIG_FILES([src/gosam/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/openloops: OpenLoops support (NLO amplitudes)

AC_CONFIG_FILES([src/openloops/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/recola: Recola support (NLO amplitudes)

AC_CONFIG_FILES([src/recola/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/fks: FKS subtraction algorithm

AC_CONFIG_FILES([src/fks/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/matching: matching algorithms

AC_CONFIG_FILES([src/matching/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/variables: Implementation of variable lists

AC_CONFIG_FILES([src/variables/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/model_features: Model access and methods

AC_CONFIG_FILES([src/model_features/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/models: Model-specific code

AC_CONFIG_FILES([src/models/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/threshold

AC_CONFIG_FILES([src/threshold/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/models/threeshl_bundle

AC_CONFIG_FILES([src/models/threeshl_bundle/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/process_integration

AC_CONFIG_FILES([src/process_integration/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/transforms: event transforms and event API

AC_CONFIG_FILES([src/transforms/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/whizard-core

AC_CONFIG_FILES([src/whizard-core/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/api

AC_CONFIG_FILES([src/api/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/main

AC_CONFIG_FILES([src/main/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/prebuilt

AC_CONFIG_FILES([src/prebuilt/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/feynmf

AC_CONFIG_FILES([src/feynmf/Makefile])

###---------------------------------------------------------------------
### Subdirectory src/gamelan: WHIZARD graphics package

AC_CONFIG_FILES([src/gamelan/Makefile])

AC_CONFIG_FILES([src/gamelan/whizard-gml], [chmod u+x src/gamelan/whizard-gml])
###---------------------------------------------------------------------
### Subdirectory share

AC_CONFIG_FILES([share/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/doc

AC_CONFIG_FILES([share/doc/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/models

AC_CONFIG_FILES([share/models/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/cuts

AC_CONFIG_FILES([share/cuts/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/beam-sim

AC_CONFIG_FILES([share/beam-sim/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/susy

AC_CONFIG_FILES([share/susy/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/examples

AC_CONFIG_FILES([share/examples/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/tests

AC_CONFIG_FILES([share/tests/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/muli

AC_CONFIG_FILES([share/muli/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/SM_tt_threshold_data

AC_CONFIG_FILES([share/SM_tt_threshold_data/Makefile])

###---------------------------------------------------------------------
### Subdirectory share/gui

AC_CONFIG_FILES([share/gui/Makefile])

###---------------------------------------------------------------------
### Subdirectory tests

AC_CONFIG_FILES([tests/Makefile])

AC_CONFIG_FILES([tests/models/Makefile])
AC_CONFIG_FILES([tests/models/UFO/Makefile])
AC_CONFIG_FILES([tests/models/UFO/SM/Makefile])
AC_CONFIG_FILES([tests/models/UFO/MSSM/Makefile])
##_CONFIG_FILES([tests/models/UFO/SMEFTsim_top_alphaScheme/Makefile])

AC_CONFIG_FILES([tests/unit_tests/Makefile])
AC_CONFIG_FILES([tests/functional_tests/Makefile])
AC_CONFIG_FILES([tests/ext_tests_mssm/Makefile])
AC_CONFIG_FILES([tests/ext_tests_nmssm/Makefile])
AC_CONFIG_FILES([tests/ext_tests_ilc/Makefile])
AC_CONFIG_FILES([tests/ext_tests_shower/Makefile])
AC_CONFIG_FILES([tests/ext_tests_nlo/Makefile])
AC_CONFIG_FILES([tests/ext_tests_nlo_add/Makefile])

AC_CONFIG_FILES([tests/unit_tests/run_whizard_ut.sh],
  [chmod u+x tests/unit_tests/run_whizard_ut.sh])
AC_CONFIG_FILES([tests/unit_tests/run_whizard_ut_c.sh],
  [chmod u+x tests/unit_tests/run_whizard_ut_c.sh])
AC_CONFIG_FILES([tests/unit_tests/run_whizard_ut_cc.sh],
  [chmod u+x tests/unit_tests/run_whizard_ut_cc.sh])
AC_CONFIG_FILES([tests/functional_tests/run_whizard.sh],
  [chmod u+x tests/functional_tests/run_whizard.sh])
AC_CONFIG_FILES([tests/ext_tests_mssm/run_whizard.sh],
  [chmod u+x tests/ext_tests_mssm/run_whizard.sh])
AC_CONFIG_FILES([tests/ext_tests_nmssm/run_whizard.sh],
  [chmod u+x tests/ext_tests_nmssm/run_whizard.sh])
AC_CONFIG_FILES([tests/ext_tests_ilc/run_whizard.sh],
  [chmod u+x tests/ext_tests_ilc/run_whizard.sh])
AC_CONFIG_FILES([tests/ext_tests_shower/run_whizard.sh],
  [chmod u+x tests/ext_tests_shower/run_whizard.sh])
AC_CONFIG_FILES([tests/ext_tests_nlo/run_whizard.sh],
  [chmod u+x tests/ext_tests_nlo/run_whizard.sh])
AC_CONFIG_FILES([tests/ext_tests_nlo_add/run_whizard.sh],
  [chmod u+x tests/ext_tests_nlo_add/run_whizard.sh])

###--------------------------------------------------------------------
### Subdirectory scripts

AC_CONFIG_FILES([scripts/Makefile])

AC_CONFIG_FILES([scripts/whizard-config], [chmod u+x scripts/whizard-config])
AC_CONFIG_FILES([scripts/whizard-setup.sh], [chmod u+x scripts/whizard-setup.sh])
AC_CONFIG_FILES([scripts/whizard-setup.csh], [chmod u+x scripts/whizard-setup.csh])


###--------------------------------------------------------------------
### CIRCE1 subdirectory files

AC_CONFIG_FILES([circe1/Makefile])
AC_CONFIG_FILES([circe1/src/Makefile])
AC_CONFIG_FILES([circe1/minuit/Makefile])
AC_CONFIG_FILES([circe1/tools/Makefile])
AC_CONFIG_FILES([circe1/share/Makefile])
AC_CONFIG_FILES([circe1/share/data/Makefile])
AC_CONFIG_FILES([circe1/share/doc/Makefile])


###--------------------------------------------------------------------
### CIRCE2 subdirectory files

AC_CONFIG_FILES([circe2/Makefile])
AC_CONFIG_FILES([circe2/src/Makefile])
AC_CONFIG_FILES([circe2/share/Makefile])
AC_CONFIG_FILES([circe2/share/doc/Makefile])
AC_CONFIG_FILES([circe2/share/examples/Makefile])
AC_CONFIG_FILES([circe2/share/data/Makefile])
AC_CONFIG_FILES([circe2/share/tests/Makefile])
AC_CONFIG_FILES([circe2/tests/Makefile])


AC_CONFIG_FILES([circe2/tests/test_wrapper.sh], [chmod u+x circe2/tests/test_wrapper.sh])
AC_CONFIG_FILES([circe2/tests/circe2_tool.sh], [chmod u+x circe2/tests/circe2_tool.sh])
AC_CONFIG_FILES([circe2/tests/generate.sh], [chmod u+x circe2/tests/generate.sh])


###--------------------------------------------------------------------
### OMEGA subdirectory files

AC_CONFIG_FILES([omega/Makefile])
AC_CONFIG_FILES([omega/bin/Makefile])
AC_CONFIG_FILES([omega/lib/Makefile])
AC_CONFIG_FILES([omega/models/Makefile])
AC_CONFIG_FILES([omega/src/Makefile])
AC_CONFIG_FILES([omega/share/Makefile])
AC_CONFIG_FILES([omega/share/doc/Makefile])
AC_CONFIG_FILES([omega/extensions/Makefile])
AC_CONFIG_FILES([omega/extensions/people/Makefile])
AC_CONFIG_FILES([omega/extensions/people/jr/Makefile])
AC_CONFIG_FILES([omega/extensions/people/tho/Makefile])
AC_CONFIG_FILES([omega/tests/Makefile])
AC_CONFIG_FILES([omega/tests/UFO/Makefile])
AC_CONFIG_FILES([omega/tests/UFO/SM/Makefile])
AC_CONFIG_FILES([omega/tests/UFO/MSSM/Makefile])
AC_CONFIG_FILES([omega/tests/UFO/SMEFTsim_top_alphaScheme/Makefile])
AC_CONFIG_FILES([omega/tools/Makefile])
AC_CONFIG_FILES([omega/scripts/Makefile])
AC_CONFIG_FILES([omega/scripts/omega-config], [chmod u+x omega/scripts/omega-config])
AC_CONFIG_FILES([omega/scripts/ufo-sanitizer], [chmod u+x omega/scripts/ufo-sanitizer])


# Copy config.mli to the build directory (otherwise ocamlc and/or
# ocamlopt would create one on their own).

###--------------------------------------------------------------------
AC_CONFIG_FILES([omega/src/config.ml])

case "$srcdir" in
 .) ;;
 *) $MKDIR_P ./omega/src
    rm -f ./omega/src/config.mli
    cp $srcdir/omega/src/config.mli ./omega/src/config.mli 1>/dev/null 2>&1;;
esac
###--------------------------------------------------------------------


###--------------------------------------------------------------------
### VAMP subdirectory files

AC_CONFIG_FILES([vamp/Makefile])
AC_CONFIG_FILES([vamp/src/Makefile])
AC_CONFIG_FILES([vamp/share/Makefile])
AC_CONFIG_FILES([vamp/share/doc/Makefile])
AC_CONFIG_FILES([vamp/tests/Makefile])


########################################################################
###---------------------------------------------------------------------
### Final output

AC_OUTPUT()

########################################################################
###---------------------------------------------------------------------
### Final output

WO_SUMMARY()

########################################################################
