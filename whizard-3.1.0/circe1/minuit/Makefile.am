# Makefile.am --
##
## Process this file with automake to produce Makefile.in
##
########################################################################
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

C_SOURCES = intrac.c

FC_SOURCES = \
   minuit.f  mncomd.f  mndxdi.f  mnexin.f  mnhess.f \
   mnlims.f  mnparm.f  mnprin.f  mnsave.f  mnstat.f \
   mnwarn.f  mnamin.f  mncont.f  mneig.f   mnfixp.f \
   mnimpr.f  mnline.f  mnpars.f  mnpsdf.f  mnscan.f \
   mnstin.f  mnwerr.f  mnbins.f  mncrck.f  mnemat.f \
   mnfree.f  mninex.f  mnmatu.f  mnpfit.f  mnrazz.f \
   mnseek.f  mntiny.f  stand.f   mncalf.f  mncros.f \
   mnerrs.f  mngrad.f  mninit.f  mnmigr.f  mnpint.f \
   mnread.f  mnset.f   mnunpt.f  mncler.f  mncuve.f \
   mneval.f  mnhelp.f  mninpu.f  mnmnos.f  mnplot.f \
   mnrn15.f  mnseti.f  mnvers.f  mncntr.f  mnderi.f \
   mnexcm.f  mnhes1.f  mnintr.f  mnmnot.f  mnpout.f \
   mnrset.f  mnsimp.f  mnvert.f

lib_LTLIBRARIES = libminuit.la
libminuit_la_SOURCES = $(FC_SOURCES) $(C_SOURCES)

EXTRA_DIST = $(FC_SOURCES) $(C_SOURCES) 

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:

AM_FFLAGS = 
AM_FCFLAGS = -I$(top_builddir)/src

########################################################################
## Default Fortran compiler options

## Profiling
if FC_USE_PROFILING
AM_FFLAGS += $(FCFLAGS_PROFILING)
AM_FCFLAGS += $(FCFLAGS_PROFILING)
endif

## OpenMP
if FC_USE_OPENMP
AM_FFLAGS += $(FCFLAGS_OPENMP)
AM_FCFLAGS += $(FCFLAGS_OPENMP)
endif

## Remove FC module files
clean-local: 
	-rm -f *.$(FCMOD)

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
maintainer-clean-local: maintainer-clean-backup

###           module="`basename $$src | sed 's/\.f[90][0358]//'`"; 
########################################################################
# MPI
########################################################################
###
### # The -mismatch_all is for mpi_send() etc.
### MPIFC = mpif90
### MPIFCFLAGS = # -mismatch_all

