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

WHIZARD_UT_DRIVER = run_whizard_ut.sh
WHIZARD_C_TEST_DRIVER = run_whizard_c_test.sh
WHIZARD_CC_TEST_DRIVER = run_whizard_cc_test.sh

UNIT_TESTS = \
    analysis.run \
    commands.run \
    pdg_arrays.run \
    expressions.run \
    beams.run \
    su_algebra.run \
    bloch_vectors.run \
    polarizations.run \
    numeric_utils.run \
    binary_tree.run \
    array_list.run \
    iterator.run \
    md5.run \
    cputime.run \
    lexers.run \
    parser.run \
    color.run \
    os_interface.run \
    evaluators.run \
    formats.run \
    sorting.run \
    grids.run \
    solver.run \
    state_matrices.run \
    interactions.run \
    xml.run \
    lorentz.run \
    phs_points.run \
    sm_qcd.run \
    sm_qed.run \
    sm_physics.run \
    electron_pdfs.run \
    models.run \
    auto_components.run \
    radiation_generator.run \
    blha.run \
    particles.run \
    beam_structures.run \
    sf_aux.run \
    sf_mappings.run \
    sf_base.run \
    sf_pdf_builtin.run \
    sf_isr.run \
    sf_epa.run \
    sf_ewa.run \
    sf_circe1.run \
    sf_circe2.run \
    sf_gaussian.run \
    sf_beam_events.run \
    sf_escan.run \
    phs_base.run \
    phs_none.run \
    phs_single.run \
    phs_rambo.run \
    resonances.run \
    phs_trees.run \
    phs_forests.run \
    phs_wood.run \
    phs_fks_generator.run \
    fks_regions.run \
    real_subtraction.run \
    rng_base.run \
    rng_tao.run \
    rng_stream.run \
    selectors.run \
    vegas.run \
    vamp2.run \
    mci_base.run \
    mci_none.run \
    mci_midpoint.run \
    mci_vamp.run \
    mci_vamp2.run \
    integration_results.run \
    prclib_interfaces.run \
    particle_specifiers.run \
    process_libraries.run \
    prclib_stacks.run \
    slha_interface.run \
    prc_test.run \
    prc_template_me.run \
    parton_states.run \
    subevt_expr.run \
    processes.run \
    process_stacks.run \
    cascades.run \
    cascades2_lexer.run \
    cascades2.run \
    event_transforms.run \
    resonance_insertion.run \
    recoil_kinematics.run \
    isr_handler.run \
    epa_handler.run \
    decays.run \
    shower.run \
    shower_base.run \
    events.run \
    hep_events.run \
    whizard_lha.run \
    pythia8.run \
    eio_data.run \
    eio_base.run \
    eio_direct.run \
    eio_raw.run \
    eio_checkpoints.run \
    eio_lhef.run \
    eio_stdhep.run \
    eio_ascii.run \
    eio_weights.run \
    eio_dump.run \
    iterations.run \
    rt_data.run \
    dispatch.run \
    dispatch_rng.run \
    dispatch_mci.run \
    dispatch_phs.run \
    dispatch_transforms.run \
    process_configurations.run \
    event_streams.run \
    integrations.run \
    ttv_formfactors.run

XFAIL_UNIT_TESTS =

UNIT_TESTS_REQ_EV_ANA = \
		phs_wood_vis.run \
		prc_omega_diags.run \
		integrations_history.run

UNIT_TESTS_REQ_FASTJET = \
		jets.run

UNIT_TESTS_REQ_HEPMC2 = \
    hepmc2.run \
    eio_hepmc2.run \
    api_hepmc2.run

UNIT_TESTS_REQ_HEPMC3 = \
    hepmc3.run \
    eio_hepmc3.run \
    api_hepmc3.run

UNIT_TESTS_REQ_LCIO = \
    lcio.run \
    eio_lcio.run \
    api_lcio.run

UNIT_TESTS_REQ_OCAML = \
    prc_omega.run \
    compilations.run \
    compilations_static.run \
    restricted_subprocesses.run \
    simulations.run \
    api.run \
    api_c.run \
    api_cc.run

UNIT_TESTS_REQ_RECOLA = \
	 prc_recola.run

UNIT_TESTS_REQ_LHAPDF5 = \
		sf_lhapdf5.run
UNIT_TESTS_REQ_LHAPDF6 = \
		sf_lhapdf6.run

TEST_DRIVERS_RUN = \
    $(UNIT_TESTS) \
    $(UNIT_TESTS_REQ_HEPMC2) \
    $(UNIT_TESTS_REQ_HEPMC3) \
    $(UNIT_TESTS_REQ_LCIO) \
    $(UNIT_TESTS_REQ_FASTJET) \
    $(UNIT_TESTS_REQ_LHAPDF5) \
    $(UNIT_TESTS_REQ_LHAPDF6) \
    $(UNIT_TESTS_REQ_OCAML) \
    $(UNIT_TESTS_REQ_RECOLA)
TEST_DRIVERS_SH = $(TEST_DRIVERS_RUN:.run=.sh)

########################################################################

TESTS =
XFAIL_TESTS =
TESTS_SRC =

UNIT_TESTS += $(UNIT_TESTS_REQ_FASTJET)
UNIT_TESTS += $(UNIT_TESTS_REQ_HEPMC2)
UNIT_TESTS += $(UNIT_TESTS_REQ_HEPMC3)
UNIT_TESTS += $(UNIT_TESTS_REQ_LCIO)
UNIT_TESTS += $(UNIT_TESTS_REQ_LHAPDF5)
UNIT_TESTS += $(UNIT_TESTS_REQ_LHAPDF6)
UNIT_TESTS += $(UNIT_TESTS_REQ_OCAML)
UNIT_TESTS += $(UNIT_TESTS_REQ_EV_ANA)
UNIT_TESTS += $(UNIT_TESTS_REQ_RECOLA)

TESTS += $(UNIT_TESTS)
XFAIL_TESTS += $(XFAIL_UNIT_TESTS)


EXTRA_DIST = $(TEST_DRIVERS_SH)
		$(TESTS_SRC)

########################################################################
# Force building the whizard_ut executable in the main src directory.
# This depends on the unit-test libraries which will be built recursively.

WHIZARD_UT = ../../src/whizard_ut

$(TEST_DRIVERS_RUN): $(WHIZARD_UT)

$(WHIZARD_UT):
	$(MAKE) -C ../../src check

########################################################################
# Force building the whizard_c_test executable in the main src directory.
# This depends on the unit-test libraries which will be built recursively.

WHIZARD_C_TEST = ../../src/whizard_c_test

$(TEST_DRIVERS_RUN): $(WHIZARD_C_TEST)

$(WHIZARD_C_TEST): $(WHIZARD_UT)

########################################################################
# Force building the whizard_c_test executable in the main src directory.
# This depends on the unit-test libraries which will be built recursively.

WHIZARD_CC_TEST = ../../src/whizard_cc_test

$(TEST_DRIVERS_RUN): $(WHIZARD_CC_TEST)

$(WHIZARD_CC_TEST): $(WHIZARD_C_TEST)

########################################################################

VPATH = $(srcdir)

SUFFIXES = .sh .run

.sh.run:
	@rm -f $@
	@cp $< $@
	@chmod +x $@

sf_beam_events.run: test_beam_events.dat
test_beam_events.dat: $(top_builddir)/share/beam-sim/test_beam_events.dat
	cp $< $@

cascades2_lexer.run: cascades2_lexer_1.fds
cascades2_lexer_1.fds: $(top_srcdir)/share/tests/cascades2_lexer_1.fds
	cp $< $@
cascades2.run: cascades2_1.fds cascades2_2.fds
cascades2_1.fds: $(top_srcdir)/share/tests/cascades2_1.fds
	cp $< $@
cascades2_2.fds: $(top_srcdir)/share/tests/cascades2_2.fds
	cp $< $@

WT_OCAML_NATIVE_EXT=opt

if MPOST_AVAILABLE
commands.run: gamelan.sty sps1ap_decays.slha
gamelan.sty: $(top_builddir)/src/gamelan/gamelan.sty
	cp $< $@

$(top_builddir)/src/gamelan/gamelan.sty:
	$(MAKE) -C $(top_builddir)/src/gamelan gamelan.sty
else
commands.run: sps1ap_decays.slha
endif

sps1ap_decays.slha: $(top_builddir)/share/susy/sps1ap_decays.slha
	cp $< $@

if OCAML_AVAILABLE

UFO_TAG_FILE = __init__.py
UFO_MODELPATH = ../models/UFO

models.run: $(UFO_MODELPATH)/SM/$(UFO_TAG_FILE)

$(UFO_MODELPATH)/SM/$(UFO_TAG_FILE): $(top_srcdir)/omega/tests/UFO/SM/$(UFO_TAG_FILE)
	$(MAKE) -C $(UFO_MODELPATH)/SM all

endif

BUILT_SOURCES = \
    TESTFLAG  \
    HEPMC2_FLAG \
    HEPMC3_FLAG \
    LCIO_FLAG \
    FASTJET_FLAG \
    LHAPDF5_FLAG \
    LHAPDF6_FLAG \
    EVENT_ANALYSIS_FLAG \
    OCAML_FLAG \
    RECOLA_FLAG \
    PYTHIA6_FLAG \
    PYTHIA8_FLAG \
    STATIC_FLAG \
    ref-output \
    err-output

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

OCAML_FLAG:
if OCAML_AVAILABLE
	touch $@
endif

RECOLA_FLAG:
if RECOLA_AVAILABLE
	touch $@
endif

PYTHIA6_FLAG:
if PYTHIA6_AVAILABLE
	touch $@
endif

PYTHIA8_FLAG:
if PYTHIA8_AVAILABLE
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

# The reference output files are in the source directory.  Copy them here.
ref-output: $(top_srcdir)/share/tests/unit_tests/ref-output
	mkdir -p ref-output
	for f in $</*.ref; do cp $$f $@; done

# The output files from self tests go here.
err-output:
	mkdir -p err-output

## installcheck runs the test scripts with the TESTFLAG removed.
## NOTE: disabled for unit tests, because whizard_ut will not be installed.
installcheck-local:
#installcheck-local: notestflag check-am
notestflag:
	rm -f TESTFLAG
.PHONY: notestflag

### Remove DWARF debug information on MAC OS X
clean-macosx:
	-rm -rf compilations_static_1.dSYM
	-rm -rf compilations_static_2.dSYM
.PHONY: clean-macosx

## Remove generated files
clean-local: clean-macosx
	rm -f gamelan.sty
	rm -f TESTFLAG 
	rm -f OCAML_FLAG FASTJET_FLAG HEPMC2_FLAG HEPMC3_FLAG LCIO_FLAG 
	rm -f RECOLA_FLAG EVENT_ANALYSIS_FLAG PYTHIA6_FLAG PYTHIA8_FLAG 
	rm -f LHAPDF5_FLAG LHAPDF6_FLAG STATIC_FLAG static_1.exe
	rm -f *.run *.log slha_test.out
	rm -f core* 
	rm -f *.f90 *.c *.$(FCMOD) *.o *.la
	rm -f *.makefile
	rm -f *.grid output.rcl
	rm -rf err-output
	rm -rf ref-output
	rm -f *.sin *.hbc *_fks_regions.out
	rm -f *.phs *.vg *.vgb *.evt *.evx *.lhe *.hepmc *.dat *.debug *.mdl
	rm -f *.tmp *.hepevt *.hepevt.verb *.lha *.lha.verb *.slcio
	rm -f prc_omega_diags_1_p_i1_diags.out prc_omega_diags_1_p_i1_diags.toc
	rm -f *.hep *.up.hep *.hep.out *.[1-9] *.[1-9][0-9] *.[1-9][0-9][0-9]
	rm -f *.tex *.mp *.mpx *.t[1-9] *.t[1-9][0-9] *.t[1-9][0-9][0-9]
	rm -f *.ltp *.aux *.dvi *.ps *.pdf so_test.*
	rm -f *.tbl sps1ap_decays.slha bar structure_6[a-b].out
	rm -f *.fds api_*.out
	rm -f *.vg2 *.vegas *.grids grids_2_test
	rm -rf output_cll
	rm -rf *.dSYM
if FC_SUBMODULES
	rm -f *.smod
endif

## Remove backup files
maintainer-clean-local: maintainer-clean-fc
	-rm -f *~
.PHONY: maintainer-clean-local
