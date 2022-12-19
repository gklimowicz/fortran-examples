# Makefile.am -- Makefile for O'Mega within and without WHIZARD
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

SUBDIRS = UFO
DIST_SUBDIRS = UFO

# OMEGA_SPLIT = -target:single_function
  OMEGA_SPLIT = -target:split_function 10
# OMEGA_SPLIT = -target:split_module 10
# OMEGA_SPLIT = -target:split_file 10

OMEGA_QED = $(top_builddir)/omega/bin/omega_QED$(OCAML_NATIVE_EXT)
OMEGA_QED_OPTS = $(OMEGA_SPLIT) -target:parameter_module parameters_QED

OMEGA_QCD = $(top_builddir)/omega/bin/omega_QCD$(OCAML_NATIVE_EXT)
OMEGA_QCD_OPTS = $(OMEGA_SPLIT) -target:parameter_module parameters_QCD

OMEGA_SYM = $(top_builddir)/omega/bin/omega_SYM$(OCAML_NATIVE_EXT)
OMEGA_SYM_OPTS = $(OMEGA_SPLIT) -target:parameter_module parameters_SYM

OMEGA_SM = $(top_builddir)/omega/bin/omega_SM$(OCAML_NATIVE_EXT)
OMEGA_SM_OPTS = $(OMEGA_SPLIT) -target:parameter_module parameters_SM

OMEGA_SM_CKM = $(top_builddir)/omega/bin/omega_SM_CKM$(OCAML_NATIVE_EXT)

OMEGA_SM_Higgs = $(top_builddir)/omega/bin/omega_SM_Higgs$(OCAML_NATIVE_EXT)

OMEGA_THDM = $(top_builddir)/omega/bin/omega_THDM$(OCAML_NATIVE_EXT)

OMEGA_THDM_CKM = $(top_builddir)/omega/bin/omega_THDM_CKM$(OCAML_NATIVE_EXT)

OMEGA_HSExt = $(top_builddir)/omega/bin/omega_HSExt$(OCAML_NATIVE_EXT)

OMEGA_Zprime = $(top_builddir)/omega/bin/omega_Zprime$(OCAML_NATIVE_EXT)

OMEGA_SM_top_anom = $(top_builddir)/omega/bin/omega_SM_top_anom$(OCAML_NATIVE_EXT)
OMEGA_SM_top_anom_OPTS = $(OMEGA_SPLIT) -target:parameter_module parameters_SM_top_anom

OMEGA_UFO = $(top_builddir)/omega/bin/omega_UFO$(OCAML_NATIVE_EXT)
OMEGA_UFO_MAJORANA = \
	$(top_builddir)/omega/bin/omega_UFO_Majorana$(OCAML_NATIVE_EXT)
OMEGA_UFO_OPTS = -target:parameter_module parameters_UFO
OMEGA_UFO_PATH = $(top_srcdir)/omega/tests/UFO

OMEGA_XXX = $(top_builddir)/omega/bin/omega_%%%$(OCAML_NATIVE_EXT)
OMEGA_XXX_OPTS = -target:parameter_module parameters_%%%
OMEGA_UFO_XXX_OPTS = \
	"-model:UFO_dir $(top_srcdir)/omega/tests/UFO/%%%/ -model:exec"
OMEGA_XXX_MAJORANA = \
	$(top_builddir)/omega/bin/omega_%%%_Majorana$(OCAML_NATIVE_EXT)
OMEGA_XXX_MAJORANA_LEGACY = \
	$(top_builddir)/omega/bin/omega_%%%_Majorana_legacy$(OCAML_NATIVE_EXT)

OMEGA_QED_VM = $(top_builddir)/omega/bin/omega_QED_VM$(OCAML_NATIVE_EXT)
OMEGA_QCD_VM = $(top_builddir)/omega/bin/omega_QCD_VM$(OCAML_NATIVE_EXT)
OMEGA_SM_VM = $(top_builddir)/omega/bin/omega_SM_VM$(OCAML_NATIVE_EXT)
OMEGA_SM_CKM_VM = $(top_builddir)/omega/bin/omega_SM_CKM_VM$(OCAML_NATIVE_EXT)
OMEGA_THDM_VM = $(top_builddir)/omega/bin/omega_THDM_VM$(OCAML_NATIVE_EXT)
OMEGA_THDM_CKM_VM = $(top_builddir)/omega/bin/omega_THDM_CKM_VM$(OCAML_NATIVE_EXT)
OMEGA_HSExt_VM = $(top_builddir)/omega/bin/omega_HSExt_VM$(OCAML_NATIVE_EXT)
OMEGA_Zprime_VM = $(top_builddir)/omega/bin/omega_Zprime_VM$(OCAML_NATIVE_EXT)
OMEGA_SM_Higgs_VM = $(top_builddir)/omega/bin/omega_SM_Higgs_VM$(OCAML_NATIVE_EXT)
OMEGA_XXX_VM = $(top_builddir)/omega/bin/omega_%%%_VM$(OCAML_NATIVE_EXT)
OMEGA_XXX_VM_PARAMS_OPTS = -params -target:parameter_module_external \
	parameters_%%% -target:wrapper_module %% -target:bytecode_file %

AM_FCFLAGS = -I$(top_builddir)/omega/src
AM_LDFLAGS =

########################################################################
## Default Fortran compiler options

## OpenMP
if FC_USE_OPENMP
AM_FCFLAGS += $(FCFLAGS_OPENMP)
AM_TESTS_ENVIRONMENT = \
	export OMP_NUM_THREADS=1;
endif

########################################################################

TESTS =
XFAIL_TESTS =
EXTRA_PROGRAMS =
EXTRA_DIST =

########################################################################

include $(top_srcdir)/omega/src/Makefile.ocaml

if OCAML_AVAILABLE

OCAMLFLAGS += -I $(top_builddir)/omega/src
OMEGA_CORE = $(top_builddir)/omega/src/omega_core.cmxa
OMEGA_MODELS = $(top_builddir)/omega/src/omega_models.cmxa

TESTS += omega_unit
EXTRA_PROGRAMS += omega_unit

omega_unit_SOURCES = omega_unit.ml

omega_unit: $(OMEGA_CORE) omega_unit.cmx
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) -o omega_unit \
		unix.cmxa $(OMEGA_CORE) omega_unit.cmx

omega_unit.cmx: omega_unit.ml

omega_unit.cmx: $(OMEGA_CORE)

endif

########################################################################

KINDS = $(top_builddir)/omega/src/kinds.lo

TESTS += test_omega95 test_omega95_bispinors
EXTRA_PROGRAMS += test_omega95 test_omega95_bispinors

test_omega95_SOURCES = test_omega95.f90 omega_testtools.f90
test_omega95_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

test_omega95_bispinors_SOURCES = test_omega95_bispinors.f90 omega_testtools.f90
test_omega95_bispinors_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

test_omega95.o test_omega95_bispinors.o: omega_testtools.o

if NOWEB_AVAILABLE

test_omega95.f90: $(top_srcdir)/omega/src/omegalib.nw
	$(NOTANGLE) -R[[$@]] $< | $(CPIF) $@
test_omega95_bispinors.f90: $(top_srcdir)/omega/src/omegalib.nw
	$(NOTANGLE) -R[[$@]] $< | $(CPIF) $@
omega_testtools.f90: $(top_srcdir)/omega/src/omegalib.nw
	$(NOTANGLE) -R[[$@]] $< | $(CPIF) $@

endif NOWEB_AVAILABLE

########################################################################

if OCAML_AVAILABLE

TESTS += test_qed_eemm
EXTRA_PROGRAMS += test_qed_eemm

test_qed_eemm_SOURCES = test_qed_eemm.f90 parameters_QED.f90
nodist_test_qed_eemm_SOURCES = amplitude_qed_eemm.f90
test_qed_eemm_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

amplitude_qed_eemm.f90: $(OMEGA_QED) Makefile
	$(OMEGA_QED) $(OMEGA_QED_OPTS) -target:module amplitude_qed_eemm \
	-scatter "e+ e- -> m+ m-" > $@

test_qed_eemm.o: amplitude_qed_eemm.o
test_qed_eemm.o: parameters_QED.o
amplitude_qed_eemm.o: parameters_QED.o

endif

########################################################################

EXTENDED_COLOR_TESTS = \
	$(srcdir)/fc_s.ects \
	$(srcdir)/fc_a.ects $(srcdir)/cf_a.ects $(srcdir)/fa_f.ects \
	$(srcdir)/ca_c.ects $(srcdir)/af_f.ects $(srcdir)/ac_c.ects \
	$(srcdir)/aa_a.ects \
	$(srcdir)/fc_fc.ects \
	$(srcdir)/aa_s.ects $(srcdir)/as_a.ects $(srcdir)/sa_a.ects

TESTS += ects
EXTRA_PROGRAMS += ects
EXTRA_DIST += ects_driver.sh $(EXTENDED_COLOR_TESTS)

# Explicitly state dependence on model files

ects.f90: $(OMEGA_QCD) $(OMEGA_SYM) $(OMEGA_SM)
ects.f90: ects_driver.sh $(EXTENDED_COLOR_TESTS)
	@if $(AM_V_P); then :; else echo "  ECTS_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/ects_driver.sh \
		$(OMEGA_XXX) $(EXTENDED_COLOR_TESTS) > $@

ects_SOURCES = color_test_lib.f90 \
	parameters_SM.f90 parameters_QED.f90 parameters_QCD.f90 parameters_SYM.f90
nodist_ects_SOURCES = ects.f90
ects_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

########################################################################

TESTS += cascade
# if there is some debugging output ...
# XFAIL_TESTS += cascade

CASCADE_TESTS = \
	bhabha-s-channel.cascade bhabha-t-channel.cascade bhabha-full.cascade \
	ww-onlycc.cascade ww-notgc.cascade \
	jjj-notgc.cascade \
	vbf-noh.cascade

cascade: cascade_driver.sh Makefile
	$(SED) -e 's|%%cascade_tests%%|$(CASCADE_TESTS)|' \
	  -e 's|%%srcdir%%|$(srcdir)|' \
	  -e 's|%%SED%%|$(SED)|' \
	  -e 's|%%top_builddir%%|$(top_builddir)|' \
	  -e 's|%%OCAML_NATIVE_EXT%%|$(OCAML_NATIVE_EXT)|' $< >$@
	chmod +x $@

EXTRA_DIST += cascade_driver.sh $(CASCADE_TESTS)

########################################################################

TESTS += phase_space

PHASE_SPACE_TESTS = eeee.phs qqggg.phs

phase_space: phase_space_driver.sh Makefile
	$(SED) -e 's|%%phase_space_tests%%|$(PHASE_SPACE_TESTS)|' \
	  -e 's|%%srcdir%%|$(srcdir)|' \
	  -e 's|%%SED%%|$(SED)|' \
	  -e 's|%%top_builddir%%|$(top_builddir)|' \
	  -e 's|%%OCAML_NATIVE_EXT%%|$(OCAML_NATIVE_EXT)|' $< >$@
	chmod +x $@

EXTRA_DIST += phase_space_driver.sh $(PHASE_SPACE_TESTS)

########################################################################

TESTS += fermi
# XFAIL_TESTS += fermi

EXTRA_PROGRAMS += fermi
EXTRA_DIST += fermi_driver.sh
EXTRA_DIST += fermi.list

FERMI_SUPPORT_F90 = \
	omega_interface.f90 omega_testtools.f90 tao_random_numbers.f90 \
	parameters_QED.f90 parameters_QCD.f90 parameters_SYM.f90 \
	parameters_SM.f90 parameters_MSSM.f90 parameters_SM_top_anom.f90
FERMI_SUPPORT_O = $(FERMI_SUPPORT_F90:.f90=.o)
fermi_lib.o: $(FERMI_SUPPORT_O)

FERMI_LIB_F90 = fermi_lib.f90 $(FERMI_SUPPORT_F90)
FERMI_LIB_O = $(FERMI_LIB_F90:.f90=.o)

run_fermi: fermi
	./fermi

fermi.f90: fermi_driver.sh $(OMEGA_QED) $(OMEGA_QCD) $(OMEGA_SYM)
fermi.f90: $(OMEGA_SM) $(OMEGA_SM_top_anom)
fermi.f90: fermi.list
	@if $(AM_V_P); then :; else echo "  FERMI_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/fermi_driver.sh \
		$(OMEGA_XXX) $(OMEGA_SPLIT) < $< > $@

fermi_SOURCES = $(FERMI_LIB_F90)
nodist_fermi_SOURCES = fermi.f90
fermi_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

fermi.o: $(FERMI_LIB_O)

########################################################################

TESTS += ward
EXTRA_PROGRAMS += ward
EXTRA_DIST += ward_driver.sh
EXTRA_DIST += ward_identities.list

WARD_SUPPORT_F90 = \
	omega_interface.f90 omega_testtools.f90 tao_random_numbers.f90 \
	parameters_QED.f90 parameters_QCD.f90 parameters_SYM.f90 \
	parameters_SM.f90 parameters_SM_top_anom.f90
WARD_SUPPORT_O = $(WARD_SUPPORT_F90:.f90=.o)
ward_lib.o: $(WARD_SUPPORT_O)

WARD_LIB_F90 = ward_lib.f90 $(WARD_SUPPORT_F90)
WARD_LIB_O = $(WARD_LIB_F90:.f90=.o)

run_ward: ward
	./ward

ward.f90: ward_driver.sh $(OMEGA_QED) $(OMEGA_QCD) $(OMEGA_SYM)
ward.f90: $(OMEGA_SM) $(OMEGA_SM_top_anom)
ward.f90: ward_identities.list
	@if $(AM_V_P); then :; else echo "  WARD_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/ward_driver.sh \
		$(OMEGA_XXX) $(OMEGA_SPLIT) < $< > $@

ward_SOURCES = $(WARD_LIB_F90)
nodist_ward_SOURCES = ward.f90
ward_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

ward.o: $(WARD_LIB_O)

########################################################################

EXTRA_PROGRAMS += ward_long
EXTRA_DIST += ward_identities_long.list

run_ward_long: ward_long
	./ward_long

ward_long.f90: ward_driver.sh
ward_long.f90: ward_identities_long.list
	@if $(AM_V_P); then :; else echo "  WARD_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/ward_driver.sh \
		$(OMEGA_XXX) $(OMEGA_SPLIT) < $< > $@

ward_long_SOURCES = $(WARD_LIB_F90)
nodist_ward_long_SOURCES = ward_long.f90
ward_long_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

# ward_long.o: ward_long.f90
# 	$(FCCOMPILE) -c -o $@ $(FCFLAGS_f90) -O0 $<

ward_long.o: $(WARD_LIB_O)

########################################################################

EXTRA_PROGRAMS += ward_fail
EXTRA_DIST += ward_identities_fail.list

run_ward_fail: ward_fail
	./ward_fail

ward_fail.f90: ward_driver.sh
ward_fail.f90: ward_identities_fail.list
	@if $(AM_V_P); then :; else echo "  WARD_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/ward_driver.sh \
		$(OMEGA_XXX) $(OMEGA_SPLIT) < $< > $@

ward_fail_SOURCES = $(WARD_LIB_F90)
nodist_ward_fail_SOURCES = ward_fail.f90
ward_fail_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

ward_fail.o: ward_fail.f90
	$(FCCOMPILE) -c -o $@ $(FCFLAGS_f90) -O0 $<

ward_fail.o: $(WARD_LIB_O)

########################################################################

TESTS += compare_split_function compare_split_module
EXTRA_PROGRAMS += compare_split_function compare_split_module
EXTRA_DIST += compare_driver.sh
EXTRA_DIST += comparisons.list

COMPARE_SUPPORT_F90 = $(WARD_SUPPORT_F90)
COMPARE_SUPPORT_O = $(WARD_SUPPORT_O)
compare_lib.o: $(COMPARE_SUPPORT_O)

COMPARE_LIB_F90 = compare_lib.f90 $(COMPARE_SUPPORT_F90)
COMPARE_LIB_O = $(COMPARE_LIB_F90:.f90=.o)

run_compare: compare_split_function compare_split_module
	./compare_split_function
	./compare_split_module

compare_split_function.f90: comparisons.list
	@if $(AM_V_P); then :; else echo "  COMPARE_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/compare_driver.sh SF \
	"$(OMEGA_XXX) -target:single_function" \
	"$(OMEGA_XXX) -target:split_function 10" < $< > $@

compare_split_module.f90: comparisons.list
	@if $(AM_V_P); then :; else echo "  COMPARE_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/compare_driver.sh SM \
	"$(OMEGA_XXX) -target:single_function" \
	"$(OMEGA_XXX) -target:split_module 10" < $< > $@

compare_split_function.f90 compare_split_module.f90: \
	compare_driver.sh $(OMEGA_QCD) $(OMEGA_SM)

compare_split_function_SOURCES = $(COMPARE_LIB_F90)
nodist_compare_split_function_SOURCES = compare_split_function.f90
compare_split_function_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

compare_split_module_SOURCES = $(COMPARE_LIB_F90)
nodist_compare_split_module_SOURCES = compare_split_module.f90
compare_split_module_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

compare_split_function.o compare_split_module.o: $(COMPARE_LIB_O)

########################################################################

if OCAML_AVAILABLE

TESTS += compare_majorana compare_majorana_legacy compare_majorana_UFO
# XFAIL_TESTS += compare_majorana_UFO
EXTRA_PROGRAMS += compare_majorana compare_majorana_legacy compare_majorana_UFO
EXTRA_DIST += compare_driver_majorana.sh compare_driver_majorana_UFO.sh
EXTRA_DIST += comparisons_majorana.list comparisons_majorana_legacy.list \
	comparisons_majorana_UFO.list

compare_majorana.f90: comparisons_majorana.list
	@if $(AM_V_P); then :; else echo "  COMPARE_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/compare_driver_majorana.sh Maj \
	"$(OMEGA_XXX)" "$(OMEGA_XXX_MAJORANA)" < $< > $@

compare_majorana_legacy.f90: comparisons_majorana_legacy.list
	@if $(AM_V_P); then :; else echo "  COMPARE_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/compare_driver_majorana.sh MajL \
	"$(OMEGA_XXX)" "$(OMEGA_XXX_MAJORANA_LEGACY)" < $< > $@

compare_majorana_UFO.f90: comparisons_majorana_UFO.list
	@if $(AM_V_P); then :; else echo "  COMPARE_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/compare_driver_majorana_UFO.sh MajU \
	"$(OMEGA_UFO)" "$(OMEGA_UFO_MAJORANA)" "$(OMEGA_UFO_PATH)" < $< > $@

compare_majorana.f90 compare_majorana_legacy.f90 compare_majorana_UFO.f90: \
	compare_driver_majorana.sh $(OMEGA_UFO) $(OMEGA_UFO_MAJORANA)

compare_majorana_SOURCES = $(COMPARE_LIB_F90)
nodist_compare_majorana_SOURCES = compare_majorana.f90
compare_majorana_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

compare_majorana_legacy_SOURCES = $(COMPARE_LIB_F90)
nodist_compare_majorana_legacy_SOURCES = compare_majorana_legacy.f90
compare_majorana_legacy_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

compare_majorana_UFO_SOURCES = $(COMPARE_LIB_F90) parameters_SM_UFO.f90
nodist_compare_majorana_UFO_SOURCES = compare_majorana_UFO.f90
compare_majorana_UFO_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

compare_majorana.o compare_majorana_legacy.o compare_majorana_UFO.o: $(COMPARE_LIB_O)
compare_majorana_UFO.o: parameters_SM_UFO.o

endif

########################################################################

if OCAML_AVAILABLE

# At quadruple or extended precision, these tests take waaaaaayyyy too long!
if FC_PREC
else

TESTS += compare_amplitude_UFO
# XFAIL_TESTS += compare_amplitude_UFO

EXTRA_PROGRAMS += compare_amplitude_UFO
EXTRA_DIST += compare_driver_UFO.sh
EXTRA_DIST += comparisons_UFO.list

compare_amplitude_UFO_SOURCES = \
	parameters_SM_from_UFO.f90 compare_lib.f90 \
	omega_interface.f90 omega_testtools.f90 tao_random_numbers.f90

compare_amplitude_UFO.f90: comparisons_UFO.list compare_driver_UFO.sh $(OMEGA_UFO)
	@if $(AM_V_P); then :; else echo "  COMPARE_DRIVER_UFO"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/compare_driver_UFO.sh UFO \
	  "$(OMEGA_XXX) -model:constant_width" \
	  "$(OMEGA_UFO) -model:UFO_dir $(top_srcdir)/omega/tests/UFO/%%%/ -model:exec" \
		< $< > $@
# -model:long_flavors

nodist_compare_amplitude_UFO_SOURCES = \
	compare_amplitude_UFO.f90 parameters_SM_UFO.f90
compare_amplitude_UFO_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

parameters_SM_from_UFO.o: parameters_SM_UFO.o
compare_amplitude_UFO.o: parameters_SM_UFO.o parameters_SM_from_UFO.o
compare_amplitude_UFO.o: $(COMPARE_LIB_O)

endif

parameters_SM_UFO.f90: $(OMEGA_UFO)
	$(OMEGA_UFO) \
	  -model:UFO_dir $(OMEGA_UFO_PATH)/SM/ -model:exec \
	  -target:parameter_module parameters_sm_ufo -params > $@
endif

########################################################################

if OCAML_AVAILABLE

# At quadruple or extended precision, these tests take waaaaaayyyy too long!
if FC_PREC
else

TESTS += fermi_UFO
# XFAIL_TESTS += fermi_UFO

# We need more work on the parameters to pass the tests
# at quadruple or extended precision.
if FC_PREC
XFAIL_TESTS += fermi_UFO
endif

EXTRA_PROGRAMS += fermi_UFO
EXTRA_DIST += fermi_driver_UFO.sh
EXTRA_DIST += fermi_UFO.list

FERMI_UFO_SUPPORT_F90 = \
	omega_interface.f90 omega_testtools.f90 tao_random_numbers.f90

FERMI_UFO_SUPPORT_O = $(FERMI_UFO_SUPPORT_F90:.f90=.o)
fermi_UFO_lib.o: $(FERMI_SUPPORT_O)

FERMI_UFO_LIB_F90 = fermi_lib.f90 $(FERMI_UFO_SUPPORT_F90)
FERMI_UFO_LIB_O = $(FERMI_UFO_LIB_F90:.f90=.o)

run_fermi_UFO: fermi_UFO
	./fermi_UFO

fermi_UFO.f90: fermi_UFO.list fermi_driver_UFO.sh $(OMEGA_UFO)
	@if $(AM_V_P); then :; else echo "  FERMI_UFO_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/fermi_driver_UFO.sh \
	  $(OMEGA_UFO) $(OMEGA_UFO_MAJORANA) $(OMEGA_UFO_PATH) \
	  $(OMEGA_SPLIT) < $< > $@

fermi_UFO_SOURCES = $(FERMI_UFO_LIB_F90)
nodist_fermi_UFO_SOURCES = fermi_UFO.f90 parameters_SM_UFO.f90
fermi_UFO_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

fermi_UFO.o: $(FERMI_UFO_LIB_O) parameters_SM_UFO.o

endif
endif

########################################################################

if OCAML_AVAILABLE

# At quadruple or extended precision, these tests take waaaaaayyyy too long!
if FC_PREC
else

TESTS += ward_UFO

# We need more work on the parameters to pass the tests
# at quadruple or extended precision.
if FC_PREC
XFAIL_TESTS += ward_UFO
endif

EXTRA_PROGRAMS += ward_UFO
EXTRA_DIST += ward_driver_UFO.sh
EXTRA_DIST += ward_identities_UFO.list

WARD_UFO_SUPPORT_F90 = \
	omega_interface.f90 omega_testtools.f90 tao_random_numbers.f90

WARD_UFO_SUPPORT_O = $(WARD_UFO_SUPPORT_F90:.f90=.o)
ward_UFO_lib.o: $(WARD_SUPPORT_O)

WARD_UFO_LIB_F90 = ward_lib.f90 $(WARD_UFO_SUPPORT_F90)
WARD_UFO_LIB_O = $(WARD_UFO_LIB_F90:.f90=.o)

run_ward_UFO: ward_UFO
	./ward_UFO

ward_UFO.f90: ward_identities_UFO.list ward_driver_UFO.sh $(OMEGA_UFO)
	@if $(AM_V_P); then :; else echo "  WARD_UFO_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/ward_driver_UFO.sh \
	  $(OMEGA_UFO) -model:UFO_dir $(top_srcdir)/omega/tests/UFO/SM/ \
	  $(OMEGA_SPLIT) < $< > $@

ward_UFO_SOURCES = $(WARD_UFO_LIB_F90)
nodist_ward_UFO_SOURCES = ward_UFO.f90 parameters_SM_UFO.f90
ward_UFO_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

ward_UFO.o: $(WARD_UFO_LIB_O) parameters_SM_UFO.o

endif

endif

########################################################################

TESTS += compare_amplitude_VM
EXTRA_PROGRAMS += compare_amplitude_VM
EXTRA_DIST += compare_driver_VM.sh compare_driver_VM_wrappers.sh
EXTRA_DIST += comparisons_VM.list

compare_amplitude_VM.f90: comparisons_VM.list comparisons_VM.wrappers.o
	@if $(AM_V_P); then :; else echo "  COMPARE_DRIVER_VM"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/compare_driver_VM.sh \
	"$(OMEGA_XXX) " "$(OMEGA_XXX_VM) " "$(OMEGA_XXX_VM_PARAMS_OPTS)" < $< > $@

comparisons_VM.wrappers.f90: comparisons_VM.list
	@if $(AM_V_P); then :; else echo "  COMPARE_DRIVER_VM_WRAPPERS"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/compare_driver_VM_wrappers.sh \
	"$(OMEGA_XXX) " "$(OMEGA_XXX_VM) " "$(OMEGA_XXX_VM_PARAMS_OPTS)" < $< > $@

# Explicitly state dependence on model files
compare_amplitude_VM.f90: compare_driver_VM.sh \
	$(OMEGA_QED)      $(OMEGA_QED_VM)      \
	$(OMEGA_QCD)      $(OMEGA_QCD_VM)      \
	$(OMEGA_SM)       $(OMEGA_SM_VM)       \
	$(OMEGA_SM_CKM)   $(OMEGA_SM_CKM_VM)   \
	$(OMEGA_SM_Higgs) $(OMEGA_SM_Higgs_VM) \
	$(OMEGA_THDM)     $(OMEGA_THDM_VM)     \
	$(OMEGA_THDM_CKM) $(OMEGA_THDM_CKM_VM) \
	$(OMEGA_HSExt)    $(OMEGA_HSExt_VM)    \
	$(OMEGA_Zprime)   $(OMEGA_Zprime_VM)   

COMPARE_EXTRA_MODELS = parameters_SM_CKM.f90 parameters_SM_Higgs.f90 \
	parameters_THDM.f90 parameters_THDM_CKM.f90 parameters_HSExt.f90 \
	parameters_Zprime.f90
compare_amplitude_VM_SOURCES = $(COMPARE_LIB_F90) $(COMPARE_EXTRA_MODELS)
nodist_compare_amplitude_VM_SOURCES = compare_amplitude_VM.f90 comparisons_VM.wrappers.f90
compare_amplitude_VM_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

compare_amplitude_VM.o: $(COMPARE_LIB_O)

########################################################################

if FC_USE_OPENMP

TESTS += test_openmp
EXTRA_PROGRAMS += test_openmp

TESTOPENMP_SUPPORT_F90 = $(WARD_SUPPORT_F90)
TESTOPENMP_SUPPORT_O = $(WARD_SUPPORT_O)

test_openmp_SOURCES = test_openmp.f90 $(TESTOPENMP_SUPPORT_F90)
nodist_test_openmp_SOURCES = amplitude_openmp.f90
test_openmp_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

amplitude_openmp.f90: $(OMEGA_QCD) Makefile
	$(OMEGA_QCD) $(OMEGA_QCD_OPTS) -target:openmp \
	-target:module amplitude_openmp -scatter "gl gl -> gl gl gl" > $@

test_openmp.o: amplitude_openmp.o
test_openmp.o: $(TESTOPENMP_SUPPORT_O)
amplitude_openmp.o: parameters_QCD.o

endif

########################################################################

EXTRA_PROGRAMS += benchmark_VM_vs_Fortran
EXTRA_DIST += benchmark_VM_vs_Fortran_driver.sh

BENCHMARK_LIB_F90 = benchmark_lib.f90 $(WARD_SUPPORT_F90)
BENCHMARK_LIB_O = $(BENCHMARK_LIB_F90:.f90=.o)

benchmark_VM_vs_Fortran.f90: benchmark_processes.list benchmark_processes.wrappers.o
	@if $(AM_V_P); then :; else echo "  BENCHMARK_VM_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/benchmark_VM_vs_Fortran_driver.sh \
	"$(OMEGA_XXX) " "$(OMEGA_XXX_VM) " "$(OMEGA_XXX_VM_PARAMS_OPTS)" < $< > $@

benchmark_processes.wrappers.f90: benchmark_processes.list
	@if $(AM_V_P); then :; else echo "  BENCHMARK_DRIVER_WRAPPERS"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/benchmark_driver_wrappers.sh \
	"$(OMEGA_XXX) " "$(OMEGA_XXX_VM) " "$(OMEGA_XXX_VM_PARAMS_OPTS)" < $< > $@

# Explicitly state dependence on model files
benchmark_VM_vs_Fortran.f90: benchmark_VM_vs_Fortran_driver.sh \
	$(OMEGA_QED) $(OMEGA_QED_VM) \
	$(OMEGA_QCD) $(OMEGA_QCD_VM) \
	$(OMEGA_SM)  $(OMEGA_SM_VM)

benchmark_VM_vs_Fortran_SOURCES = $(BENCHMARK_LIB_F90)
nodist_benchmark_VM_vs_Fortran_SOURCES = benchmark_VM_vs_Fortran.f90 benchmark_processes.wrappers.f90
benchmark_VM_vs_Fortran_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

benchmark_VM_vs_Fortran.o: $(BENCHMARK_LIB_O)

########################################################################

if FC_USE_OPENMP

EXTRA_PROGRAMS += benchmark_amp_parallel

benchmark_amp_parallel.f90: benchmark_processes.list benchmark_processes.wrappers.o
	@if $(AM_V_P); then :; else echo "  BENCHMARK_PARALLEL_DRIVER"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/benchmark_amp_parallel_driver.sh \
	"$(OMEGA_XXX) " "$(OMEGA_XXX_VM) " "$(OMEGA_XXX_VM_PARAMS_OPTS)" < $< > $@

# Explicitly state dependence on model files
benchmark_amp_parallel.f90: benchmark_amp_parallel_driver.sh \
	$(OMEGA_QED) $(OMEGA_QED_VM) \
	$(OMEGA_QCD) $(OMEGA_QCD_VM) \
	$(OMEGA_SM)  $(OMEGA_SM_VM)

benchmark_amp_parallel_SOURCES = $(BENCHMARK_LIB_F90)
nodist_benchmark_amp_parallel_SOURCES = benchmark_amp_parallel.f90 benchmark_processes.wrappers.f90
benchmark_amp_parallel_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

benchmark_amp_parallel.o: $(BENCHMARK_LIB_O)

endif

########################################################################

EXTRA_PROGRAMS += benchmark

run_benchmark: benchmark
	./benchmark

BENCHMARK_PROCESS = -scatter "gl gl -> gl gl gl"
BENCHMARK_SPLIT_SIZE = 10

benchmark_SOURCES = benchmark.f90 parameters_QCD.f90
nodist_benchmark_SOURCES = \
	amplitude_benchmark_v1.f90 amplitude_benchmark_v2.f90 \
	amplitude_benchmark_v3.f90 # amplitude_benchmark_v4.f90
benchmark_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

amplitude_benchmark_v1.f90: $(OMEGA_QCD) Makefile
	$(OMEGA_QCD) $(OMEGA_QCD_OPTS) -target:module amplitude_benchmark_v1 \
	$(BENCHMARK_PROCESS) -target:single_function > $@

amplitude_benchmark_v2.f90: $(OMEGA_QCD) Makefile
	$(OMEGA_QCD) $(OMEGA_QCD_OPTS) -target:module amplitude_benchmark_v2 \
	$(BENCHMARK_PROCESS) -target:split_function $(BENCHMARK_SPLIT_SIZE) > $@

amplitude_benchmark_v3.f90: $(OMEGA_QCD) Makefile
	$(OMEGA_QCD) $(OMEGA_QCD_OPTS) -target:module amplitude_benchmark_v3 \
	$(BENCHMARK_PROCESS) -target:split_module $(BENCHMARK_SPLIT_SIZE) > $@

amplitude_benchmark_v4.f90: $(OMEGA_QCD) Makefile
	$(OMEGA_QCD) $(OMEGA_QCD_OPTS) -target:module amplitude_benchmark_v4 \
	$(BENCHMARK_PROCESS) -target:split_file $(BENCHMARK_SPLIT_SIZE) > $@

benchmark.o: \
	amplitude_benchmark_v1.o amplitude_benchmark_v2.o \
	amplitude_benchmark_v3.o # amplitude_benchmark_v4.o
benchmark.o: parameters_QCD.o
amplitude_benchmark_v1.o amplitude_benchmark_v2.o \
	amplitude_benchmark_v3.o amplitude_benchmark_v4.o: parameters_QCD.o

########################################################################

EXTRA_PROGRAMS += benchmark_UFO_SM

run_benchmark_UFO_SM: benchmark_UFO_SM
	./benchmark_UFO_SM

# NB: This IS portable ...
UFO_SM = $(OMEGA_UFO_PATH)/SM/

BENCHMARK_UFO_SM_PROCESS = -scatter "e+ e- -> W+ W- Z Z"

benchmark_UFO_SM_SOURCES = \
	benchmark_UFO_SM.f90 parameters_SM_from_UFO.f90
nodist_benchmark_UFO_SM_SOURCES = \
	amplitude_benchmark_UFO_SM.f90 \
	amplitude_benchmark_UFO_SM_classic.f90 \
	parameters_SM_UFO.f90

benchmark_UFO_SM_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

amplitude_benchmark_UFO_SM_classic.f90: $(OMEGA_SM) Makefile
	$(OMEGA_SM) -target:module amplitude_benchmark_UFO_SM_classic \
	-target:parameter_module parameters_SM_from_UFO \
	$(BENCHMARK_UFO_SM_PROCESS) > $@

amplitude_benchmark_UFO_SM.f90: $(OMEGA_UFO) Makefile
	$(OMEGA_UFO) -model:UFO_dir $(UFO_SM) -model:exec \
	-target:module amplitude_benchmark_UFO_SM \
	-target:parameter_module parameters_SM_UFO \
	$(BENCHMARK_UFO_SM_PROCESS) > $@

benchmark_UFO_SM.o: \
	amplitude_benchmark_UFO_SM.o amplitude_benchmark_UFO_SM_classic.o

benchmark_UFO_SM.o: parameters_SM_UFO.o parameters_SM_from_UFO.o 
amplitude_benchmark_UFO_SM_classic.o: parameters_SM_from_UFO.o
amplitude_benchmark_UFO_SM.o: parameters_SM_UFO.o

########################################################################

EXTRA_PROGRAMS += benchmark_UFO_SMEFT

run_benchmark_UFO_SMEFT: benchmark_UFO_SMEFT
	./benchmark_UFO_SMEFT

# NB: This is NOT portable ...
UFO_SMEFT = /home/ohl/physics/SMEFT_mW_UFO/

BENCHMARK_UFO_SMEFT_PROCESS = -scatter "e+ e- -> W+ W- Z"

benchmark_UFO_SMEFT_SOURCES = benchmark_UFO_SMEFT.f90
nodist_benchmark_UFO_SMEFT_SOURCES = \
	amplitude_benchmark_UFO_SMEFT.f90 \
	amplitude_benchmark_UFO_SMEFT_opt.f90 \
	parameters_UFO_SMEFT.f90

benchmark_UFO_SMEFT_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

amplitude_benchmark_UFO_SMEFT.f90: $(OMEGA_UFO) Makefile
	$(OMEGA_UFO) -model:UFO_dir $(UFO_SMEFT) -model:exec \
	-target:module amplitude_benchmark_UFO_SMEFT \
	-target:parameter_module parameters_UFO_SMEFT \
	$(BENCHMARK_UFO_SMEFT_PROCESS) | $(SED) 's/g == 0/.false./' > $@

amplitude_benchmark_UFO_SMEFT_opt.f90: $(OMEGA_UFO) Makefile
	$(OMEGA_UFO) -model:UFO_dir $(UFO_SMEFT) -model:exec \
	-target:module amplitude_benchmark_UFO_SMEFT_opt \
	-target:parameter_module parameters_UFO_SMEFT \
	$(BENCHMARK_UFO_SMEFT_PROCESS) > $@

benchmark_UFO_SMEFT.o: \
	amplitude_benchmark_UFO_SMEFT.o amplitude_benchmark_UFO_SMEFT_opt.o

benchmark_UFO_SMEFT.o: parameters_UFO_SMEFT.o
amplitude_benchmark_UFO_SMEFT.o amplitude_benchmark_UFO_SMEFT_opt.o: \
	parameters_UFO_SMEFT.o

parameters_UFO_SMEFT.f90: $(OMEGA_UFO)
	$(OMEGA_UFO) -model:UFO_dir $(UFO_SMEFT) -model:exec \
	-target:parameter_module parameters_UFO_SMEFT -params > $@

########################################################################

if OCAML_AVAILABLE

TESTS += vertex_unit
EXTRA_PROGRAMS += vertex_unit
vertex_unit_SOURCES = vertex_unit.ml

vertex_unit: $(OMEGA_CORE) vertex_unit.cmx
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) -o vertex_unit \
		unix.cmxa $(OMEGA_CORE) $(OMEGA_MODELS) vertex_unit.cmx

vertex_unit.cmx: vertex_unit.ml

vertex_unit.cmx: $(OMEGA_CORE) $(OMEGA_MODELS)

endif

########################################################################

if OCAML_AVAILABLE

TESTS += ufo_unit
EXTRA_PROGRAMS += ufo_unit
ufo_unit_SOURCES = ufo_unit.ml

ufo_unit: $(OMEGA_CORE) ufo_unit.cmx
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) -o ufo_unit \
		unix.cmxa $(OMEGA_CORE) $(OMEGA_MODELS) ufo_unit.cmx

ufo_unit.cmx: ufo_unit.ml

ufo_unit.cmx: $(OMEGA_CORE) $(OMEGA_MODELS)

endif

########################################################################

if OCAML_AVAILABLE

TESTS += keystones_omegalib keystones_UFO
TESTS += keystones_omegalib_bispinors keystones_UFO_bispinors
# XFAIL_TESTS += keystones_UFO
# XFAIL_TESTS += keystones_UFO_bispinors

EXTRA_PROGRAMS += keystones_omegalib keystones_UFO
EXTRA_PROGRAMS += keystones_omegalib_bispinors keystones_UFO_bispinors

keystones_omegalib_SOURCES = omega_testtools.f90 keystones_tools.f90
nodist_keystones_omegalib_SOURCES = keystones_omegalib.f90
keystones_omegalib_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

keystones_UFO_SOURCES = omega_testtools.f90 keystones_tools.f90
nodist_keystones_UFO_SOURCES = keystones_UFO.f90
keystones_UFO_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

keystones_omegalib_bispinors_SOURCES = omega_testtools.f90 keystones_tools.f90
nodist_keystones_omegalib_bispinors_SOURCES = keystones_omegalib_bispinors.f90
keystones_omegalib_bispinors_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

keystones_UFO_bispinors_SOURCES = omega_testtools.f90 keystones_tools.f90
nodist_keystones_UFO_bispinors_SOURCES = keystones_UFO_bispinors.f90
keystones_UFO_bispinors_LDADD = $(KINDS) $(top_builddir)/omega/src/libomega_core.la

EXTRA_PROGRAMS += keystones_omegalib_generate keystones_UFO_generate
EXTRA_PROGRAMS += keystones_omegalib_bispinors_generate keystones_UFO_bispinors_generate
keystones_omegalib_generate_SOURCES = \
	keystones.ml keystones.mli keystones_omegalib_generate.ml
keystones_UFO_generate_SOURCES = \
	keystones.ml keystones.mli keystones_UFO_generate.ml
keystones_omegalib_bispinors_generate_SOURCES = \
	keystones.ml keystones.mli keystones_omegalib_bispinors_generate.ml
keystones_UFO_bispinors_generate_SOURCES = \
	keystones.ml keystones.mli keystones_UFO_bispinors_generate.ml

keystones_omegalib.f90: keystones_omegalib_generate
	./keystones_omegalib_generate -cat > $@

keystones_UFO.f90: keystones_UFO_generate
	./keystones_UFO_generate -cat > $@

keystones_omegalib_bispinors.f90: keystones_omegalib_bispinors_generate
	./keystones_omegalib_bispinors_generate -cat > $@

keystones_UFO_bispinors.f90: keystones_UFO_bispinors_generate
	./keystones_UFO_bispinors_generate -cat > $@

keystones_omegalib_generate: $(OMEGA_CORE) keystones_omegalib_generate.cmx
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) \
		-o keystones_omegalib_generate \
		unix.cmxa $(OMEGA_CORE) $(OMEGA_MODELS) \
		keystones.cmx keystones_omegalib_generate.cmx

keystones_UFO_generate: $(OMEGA_CORE) keystones_UFO_generate.cmx
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) \
		-o keystones_UFO_generate \
		unix.cmxa $(OMEGA_CORE) $(OMEGA_MODELS) \
		keystones.cmx keystones_UFO_generate.cmx

keystones_omegalib_bispinors_generate: $(OMEGA_CORE) keystones_omegalib_bispinors_generate.cmx
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) \
		-o keystones_omegalib_bispinors_generate \
		unix.cmxa $(OMEGA_CORE) $(OMEGA_MODELS) \
		keystones.cmx keystones_omegalib_bispinors_generate.cmx

keystones_UFO_bispinors_generate: $(OMEGA_CORE) keystones_UFO_bispinors_generate.cmx
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) \
		-o keystones_UFO_bispinors_generate \
		unix.cmxa $(OMEGA_CORE) $(OMEGA_MODELS) \
		keystones.cmx keystones_UFO_bispinors_generate.cmx

keystones_omegalib_generate.cmx: \
	keystones.cmi keystones.cmx keystones_omegalib_generate.ml
keystones_omegalib_generate.cmx: $(OMEGA_CORE) $(OMEGA_MODELS)

keystones_UFO_generate.cmx: \
	keystones.cmi keystones.cmx keystones_UFO_generate.ml
keystones_UFO_generate.cmx: $(OMEGA_CORE) $(OMEGA_MODELS)

keystones_omegalib_bispinors_generate.cmx: \
	keystones.cmi keystones.cmx keystones_omegalib_bispinors_generate.ml
keystones_omegalib_bispinors_generate.cmx: $(OMEGA_CORE) $(OMEGA_MODELS)

keystones_UFO_bispinors_generate.cmx: \
	keystones.cmi keystones.cmx keystones_UFO_bispinors_generate.ml
keystones_UFO_bispinors_generate.cmx: $(OMEGA_CORE) $(OMEGA_MODELS)

keystones.cmx: keystones.ml keystones.cmi
keystones.cmx: $(OMEGA_CORE) $(OMEGA_MODELS)
keystones.cmi: keystones.mli $(OMEGA_CORE)

endif

########################################################################

if RECOLA_AVAILABLE

TESTS += compare_amplitude_recola

# We need more work on the parameters to pass the tests
# at quadruple or extended precision
if FC_PREC
XFAIL_TESTS += compare_amplitude_recola
endif

EXTRA_PROGRAMS += compare_amplitude_recola
AM_FCFLAGS += $(RECOLA_INCLUDES)

compare_amplitude_recola_SOURCES = \
	parameters_SM_Higgs_recola.f90 \
	omega_interface.f90 compare_lib.f90 compare_lib_recola.f90 \
	omega_testtools.f90 tao_random_numbers.f90

nodist_compare_amplitude_recola_SOURCES = compare_amplitude_recola.f90

compare_amplitude_recola.f90: comparisons_recola.list compare_driver_recola.sh
	@if $(AM_V_P); then :; else echo "  COMPARE_DRIVER_RECOLA"; fi
	$(AM_V_at)$(SHELL) $(srcdir)/compare_driver_recola.sh \
	  "$(OMEGA_XXX) -model:constant_width" < $< > $@

compare_amplitude_recola.o: \
	omega_testtools.f90 compare_lib.o compare_lib_recola.o \
	tao_random_numbers.o \
	parameters_SM_Higgs_recola.o

compare_lib_recola.o:  \
	omega_testtools.f90 compare_lib.o tao_random_numbers.o \
	parameters_SM_Higgs_recola.o

compare_amplitude_recola_LDADD = \
	$(LDFLAGS_RECOLA) \
	$(KINDS) $(top_builddir)/omega/src/libomega_core.la

run_compare_recola: compare_amplitude_recola
	./compare_amplitude_recola

endif

########################################################################

installcheck-local:
	PATH=$(DESTDIR)$(bindir):$$PATH; export PATH; \
	LD_LIBRARY_PATH=$(DESTDIR)$(libdir):$(DESTDIR)$(pkglibdir):$$LD_LIBRARY_PATH; \
		export LD_LIBRARY_PATH; \
	omega_QED.opt $(OMEGA_QED_OPTS) -scatter "e+ e- -> m+ m-" \
		-target:module amplitude_qed_eemm > amplitude_qed_eemm.f90; \
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -I$(pkgincludedir) \
		-L$(DESTDIR)$(libdir) -L$(DESTDIR)$(pkglibdir) \
		$(srcdir)/parameters_QED.f90 amplitude_qed_eemm.f90 \
		$(srcdir)/test_qed_eemm.f90 -lomega_core; \
	./a.out

########################################################################

### Remove DWARF debug information on MAC OS X
clean-macosx:
	-rm -rf a.out.dSYM
	-rm -rf compare_amplitude_UFO.dSYM
	-rm -rf compare_amplitude_VM.dSYM
	-rm -rf compare_split_function.dSYM
	-rm -rf compare_split_module.dSYM
	-rm -rf ects.dSYM
	-rm -rf test_omega95.dSYM
	-rm -rf test_omega95_bispinors.dSYM
	-rm -rf test_qed_eemm.dSYM
	-rm -rf ward.dSYM
.PHONY: clean-macosx

clean-local: clean-macosx
	rm -f a.out gmon.out *.$(FCMOD) \
		*.o *.cmi *.cmo *.cmx amplitude_*.f90 \
		$(EXTRA_PROGRAMS) ects.f90 ward.f90 ward_UFO.f90 \
		fermi.f90 fermi_UFO.f90 compare_*.f90 \
		parameters_SM_UFO.f90 keystones_omegalib.f90 keystones_UFO.f90 \
		keystones_UFO_bispinors.f90 keystones_omegalib_bispinors.f90 \
		omega_testtools.f90 test_omega95*.f90 benchmark*.f90 \
		parameters_UFO_SMEFT.f90 \
		*.hbc *wrappers.f90 cascade phase_space \
		output.rcl recola.log
	rm -fr  output_cll

if FC_SUBMODULES
	-rm -f *.smod
endif

########################################################################
## The End.
########################################################################
