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

## Invoke extended tests by 'make check' or 'make check-extended'
check-extended: check

########################################################################

WHIZARD_DRIVER = run_whizard.sh

TESTS_EXTENDED = \
    ilc_top_pair_360.run \
    ilc_top_pair_500.run \
    ilc_vbf_higgs_360.run \
    ilc_vbf_higgs_500.run \
    ilc_vbf_no_higgs_360.run \
    ilc_vbf_no_higgs_500.run \
    ilc_higgs_strahlung_360.run \
    ilc_higgs_strahlung_500.run \
    ilc_higgs_strahlung_background_360.run \
    ilc_higgs_strahlung_background_500.run \
    ilc_higgs_coupling_360.run \
    ilc_higgs_coupling_500.run \
    ilc_higgs_coupling_background_360.run \
    ilc_higgs_coupling_background_500.run

XFAIL_TESTS_EXTENDED =

TESTS_REQ_GAMELAN = #\
    $(TESTS_REQ_OCAML) \
    $(TESTS_EXTENDED)

TEST_DRIVERS_RUN = \
    $(TESTS_EXTENDED)

TEST_DRIVERS_SH = $(TEST_DRIVERS_RUN:.run=.sh)

########################################################################

TESTS =
XFAIL_TESTS =
TESTS_SRC =

if FC_PREC
else
TESTS += $(TESTS_EXTENDED)
XFAIL_TESTS += $(XFAIL_TESTS_EXTENDED)
endif

EXTRA_DIST = $(TEST_DRIVERS_SH) \
    $(TESTS_SRC)

########################################################################

VPATH = $(srcdir)

SUFFIXES = .sh .run

.sh.run:
	@rm -f $@
	@if test -f $(top_builddir)/share/tests/ext_tests_ilc/$*.sin; then \
	  $(SED) 's|@script@|$(top_builddir)/share/tests/ext_tests_ilc/$*|g' $< > $@; \
	elif test -f $(top_srcdir)/share/tests/ext_tests_ilc/$*.sin; then \
	  $(SED) 's|@script@|$(top_srcdir)/share/tests/ext_tests_ilc/$*|g' $< > $@; \
	fi
	@chmod +x $@

if MPOST_AVAILABLE
$(TESTS_REQ_GAMELAN): gamelan.sty
$(UNIT_TESTS_REQ_GAMELAN): gamelan.sty
gamelan.sty: $(top_builddir)/src/gamelan/gamelan.sty
	cp $< $@

$(top_builddir)/src/gamelan/gamelan.sty:
	$(MAKE) -C $(top_builddir)/src/gamelan gamelan.sty
endif

BUILT_SOURCES = \
    TESTFLAG  \
    HEPMC2_FLAG \
    HEPMC3_FLAG \
    LCIO_FLAG \
    FASTJET_FLAG \
    LHAPDF5_FLAG \
    LHAPDF6_FLAG \
    GAMELAN_FLAG \
    EVENT_ANALYSIS_FLAG \
    OCAML_FLAG \
    PYTHIA6_FLAG \
    STATIC_FLAG

ilc_top_pair_360.run: ilc_settings.sin
ilc_top_pair_500.run: ilc_settings.sin
ilc_vbf_higgs_360.run: ilc_settings.sin
ilc_vbf_higgs_500.run: ilc_settings.sin
ilc_higgs_strahlung_360.run: ilc_settings.sin
ilc_higgs_strahlung_500.run: ilc_settings.sin
ilc_higgs_strahlung_background_360.run: ilc_settings.sin
ilc_higgs_strahlung_background_500.run: ilc_settings.sin
ilc_higgs_coupling_360.run: ilc_settings.sin
ilc_higgs_coupling_500.run: ilc_settings.sin
ilc_higgs_coupling_background_360.run: ilc_settings.sin
ilc_higgs_coupling_background_500.run: ilc_settings.sin
ilc_settings.sin: $(top_builddir)/share/tests/ext_tests_ilc/ilc_settings.sin
	cp $< $@

# If this file is found in the working directory, WHIZARD
# will use the paths for the uninstalled version (source/build tree),
# otherwise it uses the installed version
TESTFLAG:
	touch $@

FASTJET_FLAG:
if FASTJET_AVAILABLE
	touch $@
endif

HEPMC2_FLAG:
if HEPMC2_AVAILABLE
	touch $@
endif

HEPMC3_FLAG:
if HEPMC3_AVAILABLE
	touch $@
endif

LCIO_FLAG:
if LCIO_AVAILABLE
	touch $@
endif

LHAPDF5_FLAG:
if LHAPDF5_AVAILABLE
	touch $@
endif

LHAPDF6_FLAG:
if LHAPDF6_AVAILABLE
	touch $@
endif

GAMELAN_FLAG:
if MPOST_AVAILABLE
	touch $@
endif

OCAML_FLAG:
if OCAML_AVAILABLE
	touch $@
endif

PYTHIA6_FLAG:
if PYTHIA6_AVAILABLE
	touch $@
endif

EVENT_ANALYSIS_FLAG:
if EVENT_ANALYSIS_AVAILABLE
	touch $@
endif

STATIC_FLAG:
if STATIC_AVAILABLE
	touch $@
endif

## installcheck runs the test scripts with the TESTFLAG removed.
installcheck-local: notestflag check-am
notestflag:
	rm -f TESTFLAG
.PHONY: notestflag

## Remove generated files
clean-local:
	rm -f gamelan.sty
	rm -f TESTFLAG GAMELAN_FLAG
	rm -f OCAML_FLAG FASTJET_FLAG HEPMC2_FLAG HEPMC3_FLAG LCIO_FLAG
	rm -f EVENT_ANALYSIS_FLAG PYTHIA6_FLAG LHAPDF5_FLAG
	rm -f LHAPDF6_FLAG STATIC_FLAG static_1.exe
	rm -f *.run *.log slha_test.out
	rm -f core* stdhep_rd
	rm -f *.f90 *.c *.$(FCMOD) *.o *.la
	rm -f *.makefile
	rm -f *.sin *.hbc
	rm -f *.phs *.vg *.vgb *.evt *.evx *.lhe *.hepmc *.dat *.debug *.vg2 *.vgx2 *.fds
	rm -f *.tmp *.hepevt *.hepevt.verb *.lha *.lha.verb *.slcio
	rm -f prc_omega_diags_1_p_i1_diags.out prc_omega_diags_1_p_i1_diags.toc
	rm -f *.hep *.up.hep *.[1-9] *.[1-9][0-9] *.[1-9][0-9][0-9]
	rm -f *.tex *.mp *.mpx *.t[1-9] *.t[1-9][0-9] *.t[1-9][0-9][0-9]
	rm -f *.ltp *.aux *.dvi *.ps *.pdf so_test.*
	rm -f *.tbl sps1ap_decays.slha bar structure_6[a-b].out
if FC_SUBMODULES
	rm -f *.smod
endif

## Remove backup files
maintainer-clean-local: maintainer-clean-fc
	-rm -f *~
.PHONY: maintainer-clean-local
