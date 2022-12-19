## Makefile.am -- Makefile for executable WHIZARD test scripts
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

## Subdirectories to configure
SUBDIRS = \
  models \
  unit_tests \
  functional_tests

## Subdirs to be configured but not built
## In particular, not visited by 'make check'
DIST_SUBDIRS = \
  $(SUBDIRS) \
  ext_tests_mssm \
  ext_tests_nmssm \
  ext_tests_ilc \
  ext_tests_nlo \
  ext_tests_nlo_add \
  ext_tests_shower

## Invoke extended tests explicitly by 'make check-extended'
check-extended:
	$(MAKE) -C ext_tests_mssm check
	$(MAKE) -C ext_tests_nmssm check
	$(MAKE) -C ext_tests_ilc check
	$(MAKE) -C ext_tests_nlo check
	$(MAKE) -C ext_tests_nlo_add check
	$(MAKE) -C ext_tests_shower check

## make clean won't be executed automatically in the subdir.
## Note that after 'make distcheck', the Makefile is missing, so this fails,
## that failure should be ignored.

clean-local:
	-test -f ext_tests_mssm/Makefile && $(MAKE) -C ext_tests_mssm clean
	-test -f ext_tests_nmssm/Makefile && $(MAKE) -C ext_tests_nmssm clean
	-test -f ext_tests_ilc/Makefile && $(MAKE) -C ext_tests_ilc clean
	-test -f ext_tests_shower/Makefile && $(MAKE) -C ext_tests_shower clean
	-test -f ext_tests_nlo/Makefile && $(MAKE) -C ext_tests_nlo clean
	-test -f ext_tests_nlo_add/Makefile && $(MAKE) -C ext_tests_nlo_add clean
