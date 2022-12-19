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

WHIZARD_DRIVER = run_whizard.sh

TESTS_DEFAULT = \
    empty.run \
    fatal.run \
    cmdline_1.run \
    structure_1.run \
    structure_2.run \
    structure_3.run \
    structure_4.run \
    structure_5.run \
    structure_6.run \
    structure_7.run \
    structure_8.run \
    vars.run \
    extpar.run \
    testproc_1.run \
    testproc_2.run \
    testproc_3.run \
    testproc_4.run \
    testproc_5.run \
    testproc_6.run \
    testproc_7.run \
    testproc_8.run \
    testproc_9.run \
    testproc_10.run \
    testproc_11.run \
    testproc_12.run \
    template_me_1.run \
    template_me_2.run \
    model_scheme_1.run \
    rebuild_1.run \
    rebuild_4.run \
    susyhit.run \
    helicity.run \
    libraries_4.run \
    job_id_1.run \
    pack_1.run

XFAIL_TESTS_DEFAULT =

TESTS_REQ_FASTJET = \
    analyze_4.run \
    bjet_cluster.run \
    nlo_7.run \
    nlo_8.run \
    openloops_12.run \
    openloops_13.run \
    openloops_14.run

TESTS_REQ_OCAML = \
    libraries_1.run \
    libraries_2.run \
    libraries_3.run \
    rebuild_2.run \
    rebuild_3.run \
    rebuild_5.run \
    defaultcuts.run \
    cuts.run \
    model_change_1.run \
    model_change_2.run \
    model_change_3.run \
    model_test.run \
    job_id_2.run \
    job_id_3.run \
    job_id_4.run \
    qedtest_1.run \
    qedtest_2.run \
    qedtest_3.run \
    qedtest_4.run \
    qedtest_5.run \
    qedtest_6.run \
    qedtest_7.run \
    qedtest_8.run \
    qedtest_9.run \
    qedtest_10.run \
    rambo_vamp_1.run \
    rambo_vamp_2.run \
    beam_setup_1.run \
    beam_setup_2.run \
    beam_setup_3.run \
    beam_setup_4.run \
    beam_setup_5.run \
    qcdtest_1.run \
    qcdtest_2.run \
    qcdtest_3.run \
    qcdtest_4.run \
    qcdtest_5.run \
    qcdtest_6.run \
    observables_1.run \
    observables_2.run \
    event_weights_1.run \
    event_weights_2.run \
    event_eff_1.run \
    event_eff_2.run \
    event_dump_1.run \
    event_dump_2.run \
    event_failed_1.run \
    reweight_1.run \
    reweight_2.run \
    reweight_3.run \
    reweight_4.run \
    reweight_5.run \
    reweight_6.run \
    reweight_7.run \
    reweight_8.run \
    reweight_9.run \
    reweight_10.run \
    analyze_1.run \
    analyze_2.run \
    analyze_5.run \
    analyze_6.run \
    colors.run \
    colors_2.run \
    colors_hgg.run \
    alphas.run \
    jets_xsec.run \
    lhef_1.run \
    lhef_2.run \
    lhef_3.run \
    lhef_4.run \
    lhef_5.run \
    lhef_6.run \
    lhef_7.run \
    lhef_8.run \
    lhef_9.run \
    lhef_10.run \
    lhef_11.run \
    stdhep_1.run \
    stdhep_2.run \
    stdhep_3.run \
    stdhep_4.run \
    stdhep_5.run \
    stdhep_6.run \
    select_1.run \
    select_2.run \
    fatal_beam_decay.run \
    smtest_1.run \
    smtest_2.run \
    smtest_3.run \
    smtest_4.run \
    smtest_5.run \
    smtest_6.run \
    smtest_7.run \
    smtest_8.run \
    smtest_9.run \
    smtest_10.run \
    smtest_11.run \
    smtest_12.run \
    smtest_13.run \
    smtest_14.run \
    smtest_15.run \
    smtest_16.run \
    smtest_17.run \
    photon_isolation_1.run \
    photon_isolation_2.run \
    resonances_1.run \
    resonances_2.run \
    resonances_3.run \
    resonances_4.run \
    resonances_5.run \
    resonances_6.run \
    resonances_7.run \
    resonances_8.run \
    resonances_9.run \
    resonances_10.run \
    resonances_11.run \
    resonances_12.run \
    resonances_13.run \
    mssmtest_1.run \
    mssmtest_2.run \
    mssmtest_3.run \
    sm_cms_1.run \
    ufo_1.run \
    ufo_2.run \
    ufo_3.run \
    ufo_4.run \
    ufo_5.run \
    ufo_6.run \
    nlo_1.run \
    nlo_2.run \
    nlo_3.run \
    nlo_4.run \
    nlo_5.run \
    nlo_6.run \
    nlo_9.run \
    nlo_10.run \
    nlo_decay_1.run \
    real_partition_1.run \
    fks_res_1.run \
    fks_res_2.run \
    fks_res_3.run \
    openloops_1.run \
    openloops_2.run \
    openloops_3.run \
    openloops_4.run \
    openloops_5.run \
    openloops_6.run \
    openloops_7.run \
    openloops_8.run \
    openloops_9.run \
    openloops_10.run \
    openloops_11.run \
    recola_1.run \
    recola_2.run \
    recola_3.run \
    recola_4.run \
    recola_5.run \
    recola_6.run \
    recola_7.run \
    recola_8.run \
    recola_9.run \
    powheg_1.run \
    powheg_2.run \
    spincor_1.run \
    show_1.run \
    show_2.run \
    show_3.run \
    show_4.run \
    show_5.run \
    method_ovm_1.run \
    multi_comp_1.run \
    multi_comp_2.run \
    multi_comp_3.run \
    multi_comp_4.run \
    flvsum_1.run \
    br_redef_1.run \
    decay_err_1.run \
    decay_err_2.run \
    decay_err_3.run \
    polarized_1.run \
    pdf_builtin.run \
    ep_1.run \
    ep_2.run \
    ep_3.run \
    circe1_1.run \
    circe1_2.run \
    circe1_3.run \
    circe1_4.run \
    circe1_5.run \
    circe1_6.run \
    circe1_7.run \
    circe1_8.run \
    circe1_9.run \
    circe1_10.run \
    circe1_photons_1.run \
    circe1_photons_2.run \
    circe1_photons_3.run \
    circe1_photons_4.run \
    circe1_photons_5.run \
    circe1_errors_1.run \
    circe2_1.run \
    circe2_2.run \
    circe2_3.run \
    ewa_1.run \
    ewa_2.run \
    ewa_3.run \
    ewa_4.run \
    isr_1.run \
    isr_2.run \
    isr_3.run \
    isr_4.run \
    isr_5.run \
    isr_6.run \
    epa_1.run \
    epa_2.run \
    epa_3.run \
    epa_4.run \
    isr_epa_1.run \
    ilc.run \
    gaussian_1.run \
    gaussian_2.run \
    beam_events_1.run \
    beam_events_2.run \
    beam_events_3.run \
    beam_events_4.run \
    energy_scan_1.run \
    restrictions.run \
    process_log.run \
    shower_err_1.run \
    parton_shower_1.run \
    parton_shower_2.run \
    hadronize_1.run \
    mlm_matching_fsr.run \
    user_prc_threshold_1.run \
    cascades2_phs_1.run \
    cascades2_phs_2.run \
    user_prc_threshold_2.run \
    vamp2_1.run \
    vamp2_2.run \
    vamp2_3.run

XFAIL_TESTS_REQ_OCAML = \
    colors_hgg.run \
    hadronize_1.run

TESTS_REQ_HEPMC = \
    hepmc_1.run \
    hepmc_2.run \
    hepmc_3.run \
    hepmc_4.run \
    hepmc_5.run \
    hepmc_6.run \
    hepmc_7.run \
    hepmc_8.run \
    hepmc_9.run \
    hepmc_10.run
XFAIL_TESTS_REQ_HEPMC =

TESTS_REQ_LCIO = \
    lcio_1.run \
    lcio_2.run \
    lcio_3.run \
    lcio_4.run \
    lcio_5.run \
    lcio_6.run \
    lcio_7.run \
    lcio_8.run \
    lcio_9.run \
    lcio_10.run \
    lcio_11.run \
    lcio_12.run \
    resonances_14.run \
    resonances_15.run
XFAIL_TESTS_REQ_LCIO =

TESTS_REQ_LHAPDF5 = \
    lhapdf5.run
TESTS_REQ_LHAPDF6 = \
    lhapdf6.run
XFAIL_TESTS_REQ_LHAPDF5 =
XFAIL_TESTS_REQ_LHAPDF6 =

TESTS_STATIC = \
    static_1.run \
    static_2.run
XFAIL_TESTS_STATIC =

TESTS_REQ_PYTHIA6 = \
    pythia6_1.run \
    pythia6_2.run \
    pythia6_3.run \
    pythia6_4.run \
    tauola_1.run \
    tauola_2.run \
    tauola_3.run \
    isr_5.run \
    mlm_pythia6_isr.run \
    mlm_matching_isr.run
XFAIL_TESTS_REQ_PYTHIA6 =

TESTS_REQ_PYTHIA8 = \
    pythia8_1.run \
    pythia8_2.run
XFAIL_TESTS_REQ_PYTHIA8 =

TESTS_REQ_EV_ANA = \
    analyze_3.run
XFAIL_TESTS_REQ_EV_ANA =

TESTS_REQ_GAMELAN = \
    analyze_3.run

TEST_DRIVERS_RUN = \
    $(TESTS_DEFAULT) \
    $(TESTS_REQ_OCAML) \
    $(TESTS_REQ_LHAPDF5) \
    $(TESTS_REQ_LHAPDF6) \
    $(TESTS_REQ_HEPMC) \
    $(TESTS_REQ_LCIO) \
    $(TESTS_REQ_FASTJET) \
    $(TESTS_REQ_PYTHIA6) \
    $(TESTS_REQ_PYTHIA8) \
    $(TESTS_REQ_EV_ANA) \
    $(TESTS_STATIC)

TEST_DRIVERS_SH = $(TEST_DRIVERS_RUN:.run=.sh)

########################################################################

TESTS =
XFAIL_TESTS =
TESTS_SRC =

TESTS += $(TESTS_DEFAULT)
XFAIL_TESTS += $(XFAIL_TESTS_DEFAULT)

TESTS += $(TESTS_REQ_OCAML)
XFAIL_TESTS += $(XFAIL_TESTS_REQ_OCAML)

TESTS += $(TESTS_REQ_HEPMC)
XFAIL_TESTS += $(XFAIL_TESTS_REQ_HEPMC)

TESTS += $(TESTS_REQ_LCIO)
XFAIL_TESTS += $(XFAIL_TESTS_REQ_LCIO)

TESTS += $(TESTS_REQ_FASTJET)
XFAIL_TESTS += $(XFAIL_TESTS_REQ_FASTJET)

TESTS += $(TESTS_REQ_LHAPDF5)
XFAIL_TESTS += $(XFAIL_TESTS_REQ_LHAPDF5)

TESTS += $(TESTS_REQ_LHAPDF6)
XFAIL_TESTS += $(XFAIL_TESTS_REQ_LHAPDF6)

TESTS += $(TESTS_REQ_PYTHIA6)
XFAIL_TESTS += $(XFAIL_TESTS_REQ_PYTHIA6)

TESTS += $(TESTS_REQ_PYTHIA8)
XFAIL_TESTS += $(XFAIL_TESTS_REQ_PYTHIA8)

TESTS += $(TESTS_REQ_EV_ANA)
XFAIL_TESTS += $(XFAIL_TESTS_REQ_EV_ANA)

TESTS += $(TESTS_STATIC)
XFAIL_TESTS += $(XFAIL_TESTS_STATIC)


EXTRA_DIST = $(TEST_DRIVERS_SH) \
    $(TESTS_SRC)

########################################################################

VPATH = $(srcdir)

SUFFIXES = .sh .run

.sh.run:
	@rm -f $@
	@if test -f $(top_builddir)/share/tests/functional_tests/$*.sin; then \
		$(SED) 's|@script@|$(top_builddir)/share/tests/functional_tests/$*|g' $< > $@; \
	elif test -f $(top_srcdir)/share/tests/functional_tests/$*.sin; then \
		$(SED) 's|@script@|$(top_srcdir)/share/tests/functional_tests/$*|g' $< > $@; \
	else \
		echo "$*.sin not found!" 1>&2; \
		exit 2; \
	fi
	@chmod +x $@

cmdline_1.run: cmdline_1_a.sin cmdline_1_b.sin
cmdline_1_a.sin: $(top_builddir)/share/tests/functional_tests/cmdline_1_a.sin
	cp $< $@
cmdline_1_b.sin: $(top_builddir)/share/tests/functional_tests/cmdline_1_b.sin
	cp $< $@

structure_2.run: structure_2_inc.sin
structure_2_inc.sin: $(top_builddir)/share/tests/functional_tests/structure_2_inc.sin
	cp $< $@

testproc_3.run: testproc_3.phs
testproc_3.phs: $(top_builddir)/share/tests/functional_tests/testproc_3.phs
	cp $< $@

static_1.run: static_1.exe.sin
static_1.exe.sin: $(top_builddir)/share/tests/functional_tests/static_1.exe.sin
	cp $< $@
static_2.run: static_2.exe.sin
static_2.exe.sin: $(top_builddir)/share/tests/functional_tests/static_2.exe.sin
	cp $< $@

susyhit.run: susyhit.in

model_test.run: tdefs.$(FCMOD) tglue.$(FCMOD) \
	threeshl.$(FCMOD) tscript.$(FCMOD)
tdefs.mod: $(top_builddir)/src/models/threeshl_bundle/tdefs.$(FCMOD)
	cp $< $@
tglue.mod: $(top_builddir)/src/models/threeshl_bundle/tglue.$(FCMOD)
	cp $< $@
tscript.mod: $(top_builddir)/src/models/threeshl_bundle/tscript.$(FCMOD)
	cp $< $@
threeshl.mod: $(top_builddir)/src/models/threeshl_bundle/threeshl.$(FCMOD)
	cp $< $@

WT_OCAML_NATIVE_EXT=opt

if OCAML_AVAILABLE
OMEGA_QED = $(top_builddir)/omega/bin/omega_QED.$(WT_OCAML_NATIVE_EXT)
OMEGA_QCD = $(top_builddir)/omega/bin/omega_QCD.$(WT_OCAML_NATIVE_EXT)
OMEGA_MSSM = $(top_builddir)/omega/bin/omega_MSSM.$(WT_OCAML_NATIVE_EXT)

omega_MSSM.$(WT_OMEGA_CACHE_SUFFIX): $(OMEGA_MSSM)
	$(OMEGA_MSSM) -initialize .

UFO_TAG_FILE = __init__.py
UFO_MODELPATH = ../models/UFO

ufo_1.run: ufo_1_SM/$(UFO_TAG_FILE)
ufo_2.run: ufo_2_SM/$(UFO_TAG_FILE)
ufo_3.run: ufo_3_models/ufo_3_SM/$(UFO_TAG_FILE)
ufo_4.run: ufo_4_models/ufo_4_SM/$(UFO_TAG_FILE)
ufo_5.run: ufo_5_SM/$(UFO_TAG_FILE)
ufo_6.run: ufo_6_MSSM/$(UFO_TAG_FILE)

ufo_1_SM/$(UFO_TAG_FILE): $(UFO_MODELPATH)/SM/$(UFO_TAG_FILE)
	mkdir -p ufo_1_SM
	cp $(UFO_MODELPATH)/SM/*.py ufo_1_SM
ufo_2_SM/$(UFO_TAG_FILE): $(UFO_MODELPATH)/SM/$(UFO_TAG_FILE)
	mkdir -p ufo_2_SM
	cp $(UFO_MODELPATH)/SM/*.py ufo_2_SM
ufo_3_models/ufo_3_SM/$(UFO_TAG_FILE): $(UFO_MODELPATH)/SM/$(UFO_TAG_FILE)
	mkdir -p ufo_3_models/ufo_3_SM
	cp $(UFO_MODELPATH)/SM/*.py ufo_3_models/ufo_3_SM
ufo_4_models/ufo_4_SM/$(UFO_TAG_FILE): $(UFO_MODELPATH)/SM/$(UFO_TAG_FILE)
	mkdir -p ufo_4_models/ufo_4_SM
	cp $(UFO_MODELPATH)/SM/*.py ufo_4_models/ufo_4_SM
ufo_5_SM/$(UFO_TAG_FILE): $(UFO_MODELPATH)/SM/$(UFO_TAG_FILE)
	mkdir -p ufo_5_SM
	cp $(UFO_MODELPATH)/SM/*.py ufo_5_SM
ufo_6_MSSM/$(UFO_TAG_FILE): $(UFO_MODELPATH)/MSSM/$(UFO_TAG_FILE)
	mkdir -p ufo_6_MSSM
	cp $(UFO_MODELPATH)/MSSM/*.py ufo_6_MSSM

ufo_5.run: ufo_5_test.slha
ufo_5_test.slha: $(top_builddir)/share/tests/functional_tests/ufo_5_test.slha
	cp $< $@

$(UFO_MODELPATH)/SM/$(UFO_TAG_FILE): $(top_srcdir)/omega/tests/UFO/SM/$(UFO_TAG_FILE)
	$(MAKE) -C $(UFO_MODELPATH)/SM all
$(UFO_MODELPATH)/MSSM/$(UFO_TAG_FILE): $(top_srcdir)/omega/tests/UFO/MSSM/$(UFO_TAG_FILE)
	$(MAKE) -C $(UFO_MODELPATH)/MSSM all

endif OCAML_AVAILABLE

if MPOST_AVAILABLE
$(TESTS_REQ_GAMELAN): gamelan.sty
gamelan.sty: $(top_builddir)/src/gamelan/gamelan.sty
	cp $< $@

$(top_builddir)/src/gamelan/gamelan.sty:
	$(MAKE) -C $(top_builddir)/src/gamelan gamelan.sty
endif

noinst_PROGRAMS =

if OCAML_AVAILABLE
noinst_PROGRAMS += resonances_1_count
resonances_1_count_SOURCES = resonances_1_count.f90
resonances_1.run: resonances_1_count

noinst_PROGRAMS += resonances_2_count
resonances_2_count_SOURCES = resonances_2_count.f90
resonances_2.run: resonances_2_count

noinst_PROGRAMS += resonances_3_count
resonances_3_count_SOURCES = resonances_3_count.f90
resonances_3.run: resonances_3_count

noinst_PROGRAMS += resonances_4_count
resonances_4_count_SOURCES = resonances_4_count.f90
resonances_4.run: resonances_4_count

noinst_PROGRAMS += resonances_9_count
resonances_9_count_SOURCES = resonances_9_count.f90
resonances_9.run: resonances_9_count

noinst_PROGRAMS += resonances_10_count
resonances_10_count_SOURCES = resonances_10_count.f90
resonances_10.run: resonances_10_count

noinst_PROGRAMS += resonances_11_count
resonances_11_count_SOURCES = resonances_11_count.f90
resonances_11.run: resonances_11_count

noinst_PROGRAMS += resonances_14_check
resonances_14_check_SOURCES = resonances_14_check.f90
resonances_14.run: resonances_14_check

noinst_PROGRAMS += resonances_15_check
resonances_15_check_SOURCES = resonances_15_check.f90
resonances_15.run: resonances_15_check

noinst_PROGRAMS += epa_2_count
epa_2_count_SOURCES = epa_2_count.f90
epa_2.run: epa_2_count

noinst_PROGRAMS += isr_epa_1_count
isr_epa_1_count_SOURCES = isr_epa_1_count.f90
isr_epa_1.run: isr_epa_1_count

noinst_PROGRAMS += isr_6_digest
isr_6_digest_SOURCES = isr_6_digest.f90
isr_6.run: isr_6_digest

noinst_PROGRAMS += analyze_6_check
analyze_6_check_SOURCES = analyze_6_check.f90
analyze_6.run: analyze_6_check

endif

if HEPMC_AVAILABLE
TESTS_SRC += $(hepmc_6_rd_SOURCES)
noinst_PROGRAMS += hepmc_6_rd
if HEPMC_IS_VERSION3
hepmc_6_rd_SOURCES = hepmc3_6_rd.cpp
else
hepmc_6_rd_SOURCES = hepmc2_6_rd.cpp
endif
hepmc_6_rd_CXXFLAGS = $(HEPMC_INCLUDES) $(AM_CXXFLAGS)
hepmc_6_rd_LDADD = $(LDFLAGS_HEPMC)
hepmc_6.run: hepmc_6_rd
endif

if LCIO_AVAILABLE
TESTS_SRC += $(lcio_rd_SOURCES)
noinst_PROGRAMS += lcio_rd
lcio_rd_SOURCES = lcio_rd.cpp
lcio_rd_CXXFLAGS = $(LCIO_INCLUDES) $(AM_CXXFLAGS)
lcio_rd_LDADD = $(LDFLAGS_LCIO)
lcio_1.run: lcio_rd
lcio_2.run: lcio_rd
lcio_3.run: lcio_rd
lcio_4.run: lcio_rd
lcio_5.run: lcio_rd
lcio_10.run: lcio_rd
lcio_11.run: lcio_rd
endif

stdhep_4.run: stdhep_rd
stdhep_5.run: stdhep_rd
stdhep_6.run: stdhep_rd
polarized_1.run: stdhep_rd
tauola_1.run: stdhep_rd
tauola_2.run: stdhep_rd
tauola_3.run: stdhep_rd
stdhep_rd: $(top_builddir)/src/xdr/stdhep_rd
	cp $< $@

susyhit.in: $(top_builddir)/share/tests/functional_tests/susyhit.in
	cp $< $@

BUILT_SOURCES = \
    TESTFLAG  \
    HEPMC2_FLAG \
    HEPMC3_FLAG \
    LCIO_FLAG \
    FASTJET_FLAG \
    LHAPDF5_FLAG \
    LHAPDF6_FLAG \
    GAMELAN_FLAG \
    MPI_FLAG \
    EVENT_ANALYSIS_FLAG \
    OCAML_FLAG \
    PYTHIA6_FLAG \
    PYTHIA8_FLAG \
    OPENLOOPS_FLAG \
    RECOLA_FLAG \
    GZIP_FLAG \
    STATIC_FLAG \
    ref-output

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

MPI_FLAG:
if FC_USE_MPI
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

PYTHIA8_FLAG:
if PYTHIA8_AVAILABLE
	touch $@
endif

OPENLOOPS_FLAG:
if OPENLOOPS_AVAILABLE
	touch $@
endif

RECOLA_FLAG:
if RECOLA_AVAILABLE
	touch $@
endif

EVENT_ANALYSIS_FLAG:
if EVENT_ANALYSIS_AVAILABLE
	touch $@
endif

GZIP_FLAG:
if GZIP_AVAILABLE
	touch $@
endif

STATIC_FLAG:
if STATIC_AVAILABLE
	touch $@
endif

# The reference output files are in the source directory.  Copy them here.
if FC_QUAD
ref-output: $(top_srcdir)/share/tests/functional_tests/ref-output
	mkdir -p ref-output
	for f in $</*.ref; do cp $$f $@; done
	for f in $</../ref-output-prec/*.ref; do cp $$f $@; done
	for f in $</../ref-output-quad/*.ref; do cp $$f $@; done
else
if FC_EXT
ref-output: $(top_srcdir)/share/tests/functional_tests/ref-output
	mkdir -p ref-output
	for f in $</*.ref; do cp $$f $@; done
	for f in $</../ref-output-prec/*.ref; do cp $$f $@; done
	for f in $</../ref-output-ext/*.ref; do cp $$f $@; done
else
ref-output: $(top_srcdir)/share/tests/functional_tests/ref-output
	mkdir -p ref-output
	for f in $</*.ref; do cp $$f $@; done
	for f in $</../ref-output-double/*.ref; do cp $$f $@; done
endif
endif

## installcheck runs the test scripts with the TESTFLAG removed.
installcheck-local: notestflag check-am
notestflag:
	rm -f TESTFLAG
.PHONY: notestflag

### Remove DWARF debug information on MAC OS X
clean-macosx:
	-rm -rf hepmc_6_rd.dSYM
	-rm -rf lcio_rd.dSYM
	-rm -rf static_1.exe.dSYM
	-rm -rf static_2.exe.dSYM
.PHONY: clean-macosx

## Remove generated files
clean-local: clean-macosx
	rm -f gamelan.sty
	rm -f TESTFLAG GAMELAN_FLAG MPI_FLAG RECOLA_FLAG
	rm -f OCAML_FLAG OPENLOOPS_FLAG FASTJET_FLAG HEPMC2_FLAG HEPMC3_FLAG
	rm -f LCIO_FLAG EVENT_ANALYSIS_FLAG PYTHIA6_FLAG PYTHIA8_FLAG
	rm -f LHAPDF5_FLAG LHAPDF6_FLAG
	rm -f GZIP_FLAG STATIC_FLAG static_1.exe static_2.exe
	rm -f *.run *.log slha_test.out ufo_5_test.slha
	rm -f core* stdhep_rd
	rm -f *.f90 *.f90.in *.c *.$(FCMOD) *.o *.la
	rm -f *_count.out *_check.out
	rm -f *_digest.out
	rm -f *.makefile
	rm -f *.grid
	rm -rf err-output
	rm -rf ref-output
	rm -f *.sin *.hbc *_fks_regions.out
	rm -f *.olc *.olp *.olp_parameters output.rcl
	rm -f *.phs *.vg *.vg2 *.vgx2 *.pg *.vgb *.evt *.evx *.lhe *.hepmc *.dat *.debug *.fds
	rm -f *.tmp *.hepevt *.hepevt.verb *.lha *.lha.verb *.slcio
	rm -f prc_omega_diags_1_p_i1_diags.out prc_omega_diags_1_p_i1_diags.toc
	rm -rf pack_1.1 pack_1.2 pack_1.3 pack_1.4
	rm -rf pack_1.1.tgz pack_1.2.tgz pack_1.3.tgz pack_1.4.tgz
	rm -f *.hep *.up.hep *.[1-9] *.[1-9][0-9] *.[1-9][0-9][0-9]
	rm -f *.tex *.mp *.mpx *.t[1-9] *.t[1-9][0-9] *.t[1-9][0-9][0-9]
	rm -f *.ltp *.aux *.dvi *.ps *.pdf so_test.*
	rm -f *.tbl sps1ap_decays.slha bar structure_6[a-b].out
	rm -f slhaspectrum.in suspect2.out suspect2_lha.out
	rm -f susyhit.in susyhit_slha.out suspect2_lha.in
	rm -rf job_id_3_x.8001 job_id_4_x.8001
	rm -rf ufo_1_SM ufo_2_SM ufo_3_models ufo_4_models ufo_5_SM ufo_6_MSSM
	rm -f ufo_1_SM.ufo.mdl ufo_2_SM.ufo.mdl
	rm -f ufo_3_SM.ufo.mdl ufo_4_SM.ufo.mdl
	rm -f ufo_5_SM.ufo.mdl ufo_6_MSSM.ufo.mdl
	rm -f *.$(WT_OMEGA_CACHE_SUFFIX)
	rm -rf output_cll
	rm -rf *.dSYM
if FC_SUBMODULES
	rm -f *.smod
endif

## Remove backup files
maintainer-clean-local: maintainer-clean-fc
	-rm -f *~
.PHONY: maintainer-clean-local
