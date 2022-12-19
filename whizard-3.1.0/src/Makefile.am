## Makefile.am -- Makefile for WHIZARD
##
## Process this file with automake to produce Makefile.in
#
# Copyright (C) 1999-2022 by
#     Wolfgang Kilian <kilian@physik.uni-siegen.de>
#     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
#     Juergen Reuter <juergen.reuter@desy.de>
#     with contributions from
#     cf. main AUTHORS file
#
# WHIZARD is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# WHIZARD is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
########################################################################

## Subdirectories to configure
## directories in one line do not depend on each other
SUBDIRS = \
  noweb-frame \
  basics \
  hepmc lcio lhapdf lhapdf5 looptools xdr \
  expr_base utilities \
  system \
  combinatorics pdf_builtin testing \
  hoppet \
  muli parsing physics qed_pdf \
  fastjet types qft threshold \
  models particles matrix_elements variables \
  events prebuilt rng \
  beams tauola \
  me_methods pythia8 shower \
  blha model_features \
  gosam openloops \
  recola phase_space \
  fks vegas mci \
  process_integration \
  matching \
  transforms \
  whizard-core \
  api \
  main

if MPOST_AVAILABLE
SUBDIRS += gamelan feynmf
endif

########################################################################
## Build
## (1) the toplevel WHIZARD library, which wraps various libraries
## built in the subdirectories
## (2) the WHIZARD test library that is linked to the main executable
## as the provider of internal unit tests
## (3) the wrapper O'Mega library, which contains omegalib and
## the W/O interface modules for the various models

lib_LTLIBRARIES = libwhizard.la libomega.la

libwhizard_la_SOURCES =
libwhizard_la_LDFLAGS = $(LIBRARY_VERSION)

## Collect the various partial libraries
libwhizard_la_LIBADD = \
  api/libapi.la \
  threshold/libthreshold.la \
  whizard-core/libwhizard_core.la \
  transforms/libtransforms.la \
  ../vamp/src/libvamp.la \
  ../circe1/src/libcirce1.la \
  ../circe2/src/libcirce2.la \
  shower/libshower.la \
  tauola/libtauola_interface.la \
  muli/libmuli.la \
  pdf_builtin/libpdf_builtin.la \
  model_features/libmodel_features.la \
  variables/libvariables.la \
  process_integration/libprocess_integration.la \
  matching/libmatching.la \
  fks/libfks.la \
  gosam/libgosam.la \
  openloops/liboloops.la \
  recola/libwo_recola.la \
  blha/libblha.la \
  vegas/libvegas.la \
  mci/libmci.la \
  phase_space/libphase_space.la \
  xdr/libWOStdHep.la \
  events/libevents.la \
  beams/libbeams.la \
  particles/libparticles.la \
  me_methods/libme_methods.la \
  matrix_elements/libmatrix_elements.la \
  types/libtypes.la \
  qft/libqft.la \
  physics/libphysics.la \
  qed_pdf/libqed_pdf.la \
  expr_base/libexpr_base.la \
  rng/librng.la \
  parsing/libparsing.la \
  combinatorics/libcombinatorics.la \
  system/libsystem.la \
  testing/libtesting.la \
  utilities/libutilities.la \
  basics/libbasics.la

if HEPMC_AVAILABLE
libwhizard_la_LIBADD += hepmc/libHepMCWrap.la
else
libwhizard_la_LIBADD += hepmc/libHepMCWrap_dummy.la
endif 

if LCIO_AVAILABLE
libwhizard_la_LIBADD += lcio/libLCIOWrap.la
else
libwhizard_la_LIBADD += lcio/libLCIOWrap_dummy.la
endif

## If (parts of) LHAPDF is not available, link in a dummy as replacements
if LHAPDF5_AVAILABLE
libwhizard_la_LIBADD += $(LDFLAGS_LHAPDF)
if LHAPDF5_HAS_PHOTON_DUMMY
libwhizard_la_LIBADD += lhapdf5/libLHAPDF5_dummy.la
endif
else
if !LHAPDF6_AVAILABLE
libwhizard_la_LIBADD += lhapdf5/libLHAPDF5_dummy.la
endif
endif

if LHAPDF6_AVAILABLE
libwhizard_la_LIBADD += $(LHAPDF_LIBS) lhapdf/libLHAPDFWrap.la
else
libwhizard_la_LIBADD += lhapdf/libLHAPDFWrap_dummy.la
endif

if HOPPET_AVAILABLE
libwhizard_la_LIBADD += $(LDFLAGS_HOPPET) hoppet/libhoppet.la
else
libwhizard_la_LIBADD += hoppet/libhoppet_dummy.la
endif

if FASTJET_AVAILABLE
libwhizard_la_LIBADD += $(FASTJET_LIBS) fastjet/libFastjetWrap.la
else
libwhizard_la_LIBADD += fastjet/libFastjetWrap_dummy.la
endif

if LOOPTOOLS_AVAILABLE
libwhizard_la_LIBADD += $(LDFLAGS_LOOPTOOLS) looptools/liblooptools.la
else
libwhizard_la_LIBADD += looptools/liblooptools_dummy.la
endif

libwhizard_la_LIBADD += $(PYTHIA8_LIBS) pythia8/libwo_pythia8.la

if IS_IFORT_DARWIN
libwhizard_la_LIBADD += $(FCLIBS)
endif

## -------------------------------------------------------------------
## WHIZARD unit-test library

check_LTLIBRARIES = libwhizard_ut.la

libwhizard_ut_la_SOURCES =

libwhizard_ut_la_LIBADD = \
  main/libwhizard_main_ut.la \
  api/libapi_ut.la \
  threshold/libthreshold_ut.la \
  whizard-core/libwhizard_core_ut.la \
  transforms/libtransforms_ut.la \
  matching/libmatching_ut.la \
  fks/libfks_ut.la \
  shower/libshower_ut.la \
  process_integration/libprocess_integration_ut.la \
  recola/libwo_recola_ut.la \
  blha/libblha_ut.la \
  model_features/libmodel_features_ut.la \
  vegas/libvegas_ut.la mci/libmci_ut.la \
  phase_space/libphase_space_ut.la \
  pythia8/libwo_pythia8_ut.la \
  events/libevents_ut.la \
  beams/libbeams_ut.la \
  particles/libparticles_ut.la \
  me_methods/libme_methods_ut.la \
  matrix_elements/libmatrix_elements_ut.la \
  types/libtypes_ut.la \
  qft/libqft_ut.la \
  qed_pdf/libqed_pdf_ut.la \
  physics/libphysics_ut.la \
  rng/librng_ut.la \
  parsing/libparsing_ut.la \
  combinatorics/libcombinatorics_ut.la \
  system/libsystem_ut.la \
  utilities/libutilities_ut.la

## -------------------------------------------------------------------
## WHIZARD C API test library

check_LTLIBRARIES += libwhizard_ut_c.la

libwhizard_ut_c_la_SOURCES =

libwhizard_ut_c_la_LIBADD = \
  api/libapi_ut_c.la

## -------------------------------------------------------------------
## WHIZARD C++ API test library

check_LTLIBRARIES += libwhizard_ut_cc.la

libwhizard_ut_cc_la_SOURCES =

libwhizard_ut_cc_la_LIBADD = \
  api/libapi_ut_cc.la

## -------------------------------------------------------------------
## O'Mega main library

libomega_la_SOURCES =

libomega_la_LIBADD = \
  ../omega/src/libomega_core.la \
  models/libmodels.la

if RECOLA_AVAILABLE
libwhizard_la_LIBADD += $(LDFLAGS_RECOLA)
libwhizard_ut_la_LIBADD += $(LDFLAGS_RECOLA)
endif

########################################################################
## Build a standalone program

bin_PROGRAMS = whizard

whizard_SOURCES =

## A dummy source tells libtool that the F90 compiler is used for linking
## Without dummy, libtool uses the C linker (default: ld)
nodist_EXTRA_whizard_SOURCES = dummy.f90

whizard_LDADD = main/libwhizard_main.la
whizard_LDADD += ./libwhizard.la
whizard_LDADD += prebuilt/libwhizard_prebuilt.la
whizard_LDADD += $(CXXLIBS)
whizard_LDADD += $(RPC_CFLAGS)
whizard_LDADD += $(LDFLAGS_LHAPDF)
whizard_LDADD += $(LDFLAGS_HEPMC)
whizard_LDADD += $(LDFLAGS_LCIO)
whizard_LDADD += $(LDFLAGS_HOPPET)
whizard_LDADD += $(FASTJET_LIBS)
whizard_LDADD += $(LDFLAGS_LOOPTOOLS)

## -------------------------------------------------------------------
## Build a standalone program for running the Fortran unit tests

check_PROGRAMS = whizard_ut

whizard_ut_SOURCES =

## A dummy source tells libtool that the F90 compiler is used for linking
## Without dummy, libtool uses the C linker (default: ld)
nodist_EXTRA_whizard_ut_SOURCES = dummy.f90

whizard_ut_LDADD = ./libwhizard_ut.la
whizard_ut_LDADD += ./libwhizard.la
whizard_ut_LDADD += prebuilt/libwhizard_prebuilt.la
whizard_ut_LDADD += $(CXXLIBS)
whizard_ut_LDADD += $(RPC_CFLAGS)
whizard_ut_LDADD += $(LDFLAGS_LHAPDF)
whizard_ut_LDADD += $(LDFLAGS_HEPMC)
whizard_ut_LDADD += $(LDFLAGS_LCIO)
whizard_ut_LDADD += $(LDFLAGS_HOPPET)
whizard_ut_LDADD += $(FASTJET_LIBS)
whizard_ut_LDADD += $(LDFLAGS_LOOPTOOLS)

## -------------------------------------------------------------------
## Build a standalone program for running the C interface test

check_PROGRAMS += whizard_ut_c

whizard_ut_c_SOURCES =

whizard_ut_c_LDADD = ./libwhizard_ut_c.la
whizard_ut_c_LDADD += ./libwhizard.la
whizard_ut_c_LDADD += prebuilt/libwhizard_prebuilt.la
whizard_ut_c_LDADD += $(CXXLIBS)
whizard_ut_c_LDADD += $(RPC_CFLAGS)
whizard_ut_c_LDADD += $(LDFLAGS_LHAPDF)
whizard_ut_c_LDADD += $(LDFLAGS_HEPMC)
whizard_ut_c_LDADD += $(LDFLAGS_LCIO)
whizard_ut_c_LDADD += $(LDFLAGS_HOPPET)
whizard_ut_c_LDADD += $(FASTJET_LIBS)
whizard_ut_c_LDADD += $(LDFLAGS_LOOPTOOLS)

## -------------------------------------------------------------------
## Build a standalone program for running the C++ interface test

check_PROGRAMS += whizard_ut_cc

## A dummy source tells libtool that the C++ compiler is used for linking
## Without dummy, libtool uses the C linker (default: ld)
whizard_ut_cc_SOURCES =
nodist_EXTRA_whizard_ut_cc_SOURCES = dummy_cc.cc

whizard_ut_cc_LDADD = ./libwhizard_ut_cc.la
whizard_ut_cc_LDADD += ./libwhizard.la
whizard_ut_cc_LDADD += prebuilt/libwhizard_prebuilt.la
whizard_ut_cc_LDADD += $(CXXLIBS)
whizard_ut_cc_LDADD += $(RPC_CFLAGS)
whizard_ut_cc_LDADD += $(LDFLAGS_LHAPDF)
whizard_ut_cc_LDADD += $(LDFLAGS_HEPMC)
whizard_ut_cc_LDADD += $(LDFLAGS_LCIO)
whizard_ut_cc_LDADD += $(LDFLAGS_HOPPET)
whizard_ut_cc_LDADD += $(FASTJET_LIBS)
whizard_ut_cc_LDADD += $(LDFLAGS_LOOPTOOLS)

########################################################################
## Default Fortran compiler options

AM_FCFLAGS =

## Profiling
if FC_USE_PROFILING
AM_FCFLAGS += $(FCFLAGS_PROFILING)
endif

## OpenMP
if FC_USE_OPENMP
AM_FCFLAGS += $(FCFLAGS_OPENMP)
endif

########################################################################
## Non-standard cleanup tasks

## Remove backup files
maintainer-clean-local:
	-rm -f *~

## Remove Module lists only after final subdir cleanup
distclean-local:
	for d in $(SUBDIRS); do \
	  rm -f $$d/Modules; \
	done
