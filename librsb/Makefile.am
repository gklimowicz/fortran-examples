ACLOCAL_AMFLAGS=-I m4

AM_LDFLAGS=-no-undefined -version-info $(LIBRSB_ABI_VERSION)

SUBDIRS=
dist_doc_DATA=README AUTHORS README.md
bin_SCRIPTS=librsb-config
bin_PROGRAMS=rsbench$(EXEEXT)
noinst_PROGRAMS= # rsbenchxx # cg jacobi
EXTRAPROGRAMSFORTRAN=
EXTRAPROGRAMSC=ot$(EXEEXT)
LIBS+=$(FCLIBS)
# L1C 
if HAVE_LIBRSBPP
SUBDIRS+=@librsbpp@
endif
SUBDIRS+= .
SUBDIRS+=examples scripts bench blas_sparse doc m4
if HAVE_RSBLIB
SUBDIRS+=@rsblib@
endif
if HAVE_RSBTEST
SUBDIRS+=@rsbtest@
endif
if HAVE_GTEST
noinst_PROGRAMS+=rsb_gtests$(EXEEXT)
rsb_gtests_SOURCES=rsb_gtests.cpp
rsb_gtests_LIBTOOLFLAGS= --preserve-dup-deps
rsb_gtests_LDADD= $(default_ldadd) -lgtest_main -lgtest -lpthread -lstdc++fs
endif
if HAVE_SPARSE_BLAS_INTERFACE
EXTRAPROGRAMSC+=sbtc$(EXEEXT)
sbtc_SOURCES=sbtc.c
if HAVE_LIBRSBPP
sbtc_SOURCES+=rsb_dummy.cpp
endif
sbtc_CFLAGS=$(CFLAGS) -D SBTC_USE_RSB_H
sbtc_LDADD=$(default_ldadd)
endif
if HAVE_FORTRAN_EXAMPLES
if HAVE_SPARSE_BLAS_INTERFACE
EXTRAPROGRAMSFORTRAN+=sbtf$(EXEEXT)
endif
endif
EXTRA_PROGRAMS= $(EXTRAPROGRAMSC) $(EXTRAPROGRAMSFORTRAN)
EXTRA_LIBRARIES=
EXTRA_SCRIPTS=librsb-config

NODISTSRC=  rsb_mkl.h blas_sparse.h

RSBENCH_BUILT_SOURCES=	\
		rsb_test_matops.c rsb_test_matops.h	\
		rsb_mkl.c rsb_mkl.h

LIB_BUILT_SOURCES_M4=	\
		rsb_unroll.c rsb_unroll.h \
		rsb_krnl_vb.c rsb_krnl_vb.h \
		rsb_krnl_lb.c rsb_krnl_lb.h \
		rsb_krnl.c rsb_krnl.h \
		rsb_bench.c rsb_bench.h \
		rsb_types.h \
		rsb_mergesort.h rsb_mergesort.c \
		rsb_permute.h rsb_permute.c \
		rsb_krnl_bcss_l.h rsb_krnl_bcss_l.c \
		rsb_krnl_bcss_u.h rsb_krnl_bcss_u.c \
		rsb_krnl_bcss_spsv_u.h rsb_krnl_bcss_spsv_u.c \
		rsb_krnl_bcss_spmv_u.h rsb_krnl_bcss_spmv_u.c \
		rsb_krnl_bcss_misc_u.h rsb_krnl_bcss_misc_u.c \
		rsb_krnl_bcoo_spmv_u.h rsb_krnl_bcoo_spmv_u.c \
		rsb_krnl_bcss.h rsb_krnl_bcss.c \
		rsb_spmv.h rsb_spmv.c \
		rsb_merge.h rsb_merge.c \
		rsb_ompio.c rsb_ompio.h \
		rsb_util.h rsb_util.c \
		rsb_spgemm_csr.c rsb_spgemm_csr.h \
		rsb_spsum_misc.c rsb_spsum_misc.h \
		rsb_prec.h rsb_prec.c
if WANT_RSB_USE_MKL
LIB_BUILT_SOURCES_M4+= rsb_mkl.c
endif

LIB_BUILT_SOURCES_MAKEFILE=
if HAVE_SED
LIB_BUILT_SOURCES_MAKEFILE+= \
		rsb_stropts.c \
		rsb_strmif.c \
		rsb_strmif.h
endif
if HAVE_M4
LIB_BUILT_SOURCES_MAKEFILE+= rsb_types.h blas_sparse.h
endif

LIB_BUILT_SOURCES=	\
	$(LIB_BUILT_SOURCES_MAKEFILE) $(LIB_BUILT_SOURCES_M4)
	
# rsb_prec.h: $(NODISTSRC)

include_HEADERS=
nodist_include_HEADERS=rsb-config.h
EXTRAHEADERSFORTRAN=

# automake uses FCLD even if no F90 file is to be compiled (e.g. FC='') this is a workaround:
if !HAVE_FC
#FCLD=$(CC)
endif
# ...but this is a better one (no need to overwrite something set by automake):
librsb_la_LINK=$(LINK)
librsb_base_la_LINK=$(LINK)
librsb_nounroll_la_LINK=$(LINK)

EXTRASRCFILESSFORTRAN=
EXTRASRCPROGSSFORTRAN=
if HAVE_FC
if HAVE_SPARSE_BLAS_INTERFACE
EXTRASRCFILESSFORTRAN+=rsb_blas_sparse.F90
EXTRASRCPROGSSFORTRAN+=rsb_blas_sparse.F90
if WANT_BLAS_SPARSE_MOD_INSTALL
nodist_include_HEADERS+=blas_sparse.mod
endif
blas_sparse.mod: rsb_blas_sparse.$(OBJEXT)
endif
endif

if HAVE_FC
nodist_include_HEADERS+=rsb.mod
rsb.mod: rsb.$(OBJEXT)
# if WANT_BLAS_SPARSE_FI
# nodist_include_HEADERS+=blas_sparse.fi
# blas_sparse.fi: blas_sparse.m4
# 	$(M4) $(M4_FLAGS) -D RSB_M4_WANT_BLAS_SPARSE_INTERFACE=1 $< > $@
# endif
EXTRAHEADERSFORTRAN+=rsb.F90
if WANT_RSB_F90
rsb.F90: ch2icfb $(top_srcdir)/rsb.h $(top_srcdir)/scripts/rsb_h_to_rsb_fi.sh rsb_types.h
	if test -f $(srcdir)/rsb_types.h -a ! -f $(builddir)/rsb_types.h ; then cp -pv $(srcdir)/rsb_types.h $(builddir)/rsb_types.h ; fi # out-of-dir behaviour varies between installations
	if test -f ch2icfb ; then if SED=$(SED) GREP=$(GREP) $(BASH) $(top_srcdir)/scripts/rsb_h_to_rsb_fi.sh $(srcdir) $(abs_builddir) > $@ ; then true; else rm "$@"; echo "Failure generating $@!" ; exit ; fi ; else echo "Warning: Your system did not build ch2icfb for some reason --- skipping rebuild of "$@ ;true ; fi

noinst_PROGRAMS+=ch2icfb$(EXEEXT) 
ch2icfb_SOURCES=ch2icfb.c
ch2icfb_CFLAGS=$(CFLAGS)
ch2icfb_LDADD=
endif
else
EXTRASRCFILESSFORTRAN+=
EXTRASRCPROGSSFORTRAN+=
endif

if HAVE_FORTRAN_EXAMPLES
sbtf_LDADD=$(default_ldadd)
sbtf.$(OBJEXT): rsb.$(OBJEXT)
if HAVE_SPARSE_BLAS_INTERFACE
sbtf_SOURCES=sbtf.F90
EXTRASRCFILESSFORTRAN+=
EXTRASRCPROGSSFORTRAN+=sbtf.F90 psbtf.F90 psb_mvsv_tester.f90
else
sbtf_SOURCES=
EXTRASRCFILESSFORTRAN+=
EXTRASRCPROGSSFORTRAN+=
endif
else
EXTRASRCFILESSFORTRAN+=
EXTRASRCPROGSSFORTRAN+=
endif

BUILT_SOURCES=	$(LIB_BUILT_SOURCES) $(RSBENCH_BUILT_SOURCES) $(ot_SOURCES) $(LIB_SPBLAS_BUILT_SOURCES) $(EXTRASRCFILESSFORTRAN) rsb-config.h
#BUILT_SOURCES+=rsbenchxx.cpp

.PHONY: e
e:
	vim $(top_srcdir)/rsb.h

DIST_ARCHIVES_NOVERSION = librsb.tar.gz

.PHONY: tdist
tdist: dox
	$(MAKE) dist VERSION=trunk

.PHONY: rtdist
rtdist: tdist
	mv librsb-trunk.tar.gz $(DIST_ARCHIVES_NOVERSION) 

.PHONY: rdist
rdist: dox
	$(MAKE) dist

.PHONY: ddist
ddist: rdist
	gpg -sbv -u 0xe0e669c8ef1258b8 $(DIST_ARCHIVES)
	md5sum $(DIST_ARCHIVES) > $(DIST_ARCHIVES).md5
	gpg -sbav $(DIST_ARCHIVES)
	gpg --verify $(DIST_ARCHIVES).sig

.PHONY: rdist
rdistnv: dox
	$(MAKE) dist
	mv $(DIST_ARCHIVES) $(DIST_ARCHIVES_NOVERSION) 

.PHONY: ddist
ddistnv: rdist
	gpg -sbv -u 0xe0e669c8ef1258b8 $(DIST_ARCHIVES_NOVERSION)
	md5sum $(DIST_ARCHIVES_NOVERSION) > $(DIST_ARCHIVES_NOVERSION).md5
	gpg -sbav $(DIST_ARCHIVES_NOVERSION)
	gpg --verify $(DIST_ARCHIVES_NOVERSION).sig

.PHONY: bdist
bdist: dox $(BINDISTFILES) 
	rm -fR -- $(PACKAGE)-$(build)
	$(MKDIR_P) $(PACKAGE)-$(build)
	cp -fR $(BINDISTFILES) $(PACKAGE)-$(build)/
	tar cvzf $(PACKAGE)-$(build).tgz $(PACKAGE)-$(build)  --exclude .git  --exclude .deps
	rm -fR $(PACKAGE)-$(build)
	tar tvzf $(PACKAGE)-$(build).tgz

.PHONY: help
help:
	@echo -e "Alternatives (see the README for these):\n make clean\n make cleanall\n make realclean\n make all \n make qqtests \n make qtests \n make tests \n make install \n make itests \n make dist"

if HAVE_M4
LIB_CLEANALL_FILES=$(LIB_BUILT_SOURCES) $(LIB_SPBLAS_BUILT_SOURCES) $(RSBENCH_BUILT_SOURCES) \
	$(EXTRASRCFILESSFORTRAN)
else
LIB_CLEANALL_FILES=$(LIB_BUILT_SOURCES_MAKEFILE)
endif
if WANT_INTERNAL_HEADERS_INSTALL
LIB_CLEANALL_FILES+=rsb_stub.c rsb-librsb-internals.h
endif

.PHONY: cleanall-local
cleanall-local:
	rm -rf $(LIB_CLEANALL_FILES)

.PHONY: cleanall
cleanall: clean
	$(MAKE) cleanall-local
	$(MAKE) clean

.PHONY: gclean
gclean:
	rm -rf  *.gcov *.gcno

.PHONY: realclean-local
realclean-local:
	$(MAKE) gclean
	rm -f *.rpr

noinst_LTLIBRARIES=librsb_nounroll.la librsb_base.la librsb_spblas.la
lib_LTLIBRARIES=librsb.la

librsb_la_SOURCES=
if HAVE_FC
if WANT_RSB_F90
# after Cray Fortran 14:
librsb_la_SOURCES+=rsb.F90
endif
endif
#am_librsb_la_OBJECTS=$(am_librsb_base_la_OBJECTS) $(am_librsb_nounroll_la_OBJECTS) $(am_librsb_spblas_la_OBJECTS) # detached librsbpp.la
librsb_la_LIBADD=librsb_nounroll.la librsb_base.la librsb_spblas.la
if HAVE_LIBRSBPP
#am_librsb_la_OBJECTS+=$(abs_top_builddir)/librsbpp/librsbpp.la # detached librsbpp.la
librsb_la_LIBADD+=$(abs_top_builddir)/librsbpp/librsbpp.la $(RSB_RSBPP_LIBS)
endif

#librsb_la_LIBADD=-lgfortran
default_ldadd=-L$(top_builddir)  $(top_builddir)/librsb.la
if WANT_RSB_USE_MKL
default_ldadd+=${RSB_RSBENCH_LIBS}
endif

librsb_nounroll_la_CFLAGS=$(NOUNROLLCFLAGS) $(KERNELS_EXTRA_CFLAGS) $(OPENMP_CFLAGS)
librsb_nounroll_la_SOURCES=$(LIB_BUILT_SOURCES) $(NODISTSRC)
nodist_librsb_nounroll_la_SOURCES=$(NODISTSRC)
librsb_base_la_CFLAGS=$(BASELIB_EXTRA_CFLAGS) $(OPENMP_CFLAGS)
nodist_librsb_base_la_SOURCES=$(NODISTSRC)
librsb_base_la_SOURCES= 	\
			rsb_common.h \
			rsb_is.c  rsb_is.h \
			rsb_mio.c rsb_mio.h \
			rsb_op.c  rsb_op.h \
			rsb_bio.c rsb_bio.h \
			rsb_get.c rsb_get.h \
			rsb_set.c rsb_set.h \
			rsb_coo.c rsb_coo.h \
			rsb_csr.c rsb_csr.h \
			rsb_coo_check.c rsb_coo_check.h \
			rsb_coo_symm.c rsb_coo_symm.h \
			rsb_idx.c rsb_idx.h \
			rsb_srt.c rsb_srt.h \
			rsb_srtp.c rsb_srtp.h \
			rsb_src.c rsb_src.h \
			rsb_test_accuracy.c rsb_test_accuracy.h \
			rsb_clone.c rsb_clone.h \
			rsb_rec.h rsb_rec.c \
			rsb_render.c rsb_render.h \
			rsb_eps.c rsb_eps.h \
			rsb_msort_up.c rsb_msort_up.h \
			rsb_sys.c rsb_sys.h \
			rsb_blas_stuff.c rsb_blas_stuff.h \
			rsb_gen.c rsb_gen.h \
			rsb_perf.c rsb_perf.h \
			rsb_rsb.c rsb.h \
			rsb_err.c rsb_err.h \
			rsb_tune.c rsb_tune.h \
			rsb_struct.h \
			rsb_do.c rsb_do.h \
			rsb_internals.c rsb_internals.h \
			rsb_garbage.c rsb_garbage.h \
			rsb_mmio.c rsb_mmio.h \
			rsb_partition.c rsb_partition.h \
		       	rsb_lbl.h \
			rsb_mbw.c rsb_mbw.h \
			rsb_limiter.c rsb_limiter.h \
			rsb_fpb.c rsb_fpb.h \
			rsb_spgemm.c rsb_spgemm.h \
			rsb_spsum.c rsb_spsum.h \
			rsb_spsv.c rsb_spsv.h \
			rsb_lock.h rsb_lock.c \
			rsb_swt.h rsb_swt.c \
			rsb_init.h rsb_init.c \
			rsb_dump.h rsb_dump.c \
			rsb_cpmv.h rsb_cpmv.c \
			rsb_psblas.h \
			rsb_asm.h  rsb_asm.c  \
			rsb_user.c \
			rsb_coo2rec.c rsb_coo2rec.h \
			rsb_rec2coo.c rsb_rec2coo.h \
			rsb_rec2csr.c rsb_rec2csr.h \
			rsb_csr2coo.c rsb_csr2coo.h \
			$(EXTRASRCFILESSFORTRAN)

LIB_SPBLAS_BUILT_SOURCES= rsb_libspblas.c rsb_libspblas.h

librsb_spblas_la_SOURCES= \
		rsb_libspblas_handle.h	rsb_libspblas_handle.c	\
		$(LIB_SPBLAS_BUILT_SOURCES)

if HAVE_M4
rsb_libspblas.h: $(srcdir)/rsb_libspblas.m4 $(srcdir)/libspblas_macros.m4 blas_sparse/blas_enum.h
rsb_libspblas.c: $(srcdir)/rsb_libspblas.m4 $(srcdir)/libspblas_macros.m4
psb_rsb_mod.F90: $(srcdir)/psb_rsb_mod.m4 $(srcdir)/rsb_fortran_macros.m4
rsb_mod.F90: $(srcdir)/rsb_mod.m4 $(srcdir)/rsb_fortran_macros.m4
if HAVE_FC
rsb_blas_sparse.F90: $(srcdir)/rsb_blas_sparse.m4 $(srcdir)/rsb_fortran_macros.m4
endif
rsb_libspblas_handle.c: rsb_libspblas.h rsb_libspblas_handle.h
rsb_libspblas_tests.c: rsb_libspblas.h rsb_libspblas_tests.h
rsb_mkl.c: $(srcdir)/rsb_mkl.m4
rsb_mkl.h: $(srcdir)/rsb_mkl.m4
endif

if WANT_INTERNAL_HEADERS_INSTALL
rsb-librsb-internals.h: $(librsb_base_la_SOURCES) $(librsb_nounroll_la_SOURCES)
	( cat $(top_srcdir)/rsb_license_header.inc                                                                            ; \
       	echo '/*! Collated internal headers of librsb -- for inspection purposes only (not for usage). */'; echo ; \
	echo '/* @cond INNERDOC */' ; \
	echo '#ifndef RSB_LIBRSB_INTERNALS_H_INCLUDED'									; \
	echo '#define RSB_LIBRSB_INTERNALS_H_INCLUDED'									; \
	find $+ -iname '*.h' -exec 'cat' '{}' ';' | $(GREP) -v 'cond INNERDOC\|endcond' ; \
	echo '#endif /* RSB_LIBRSB_INTERNALS_H_INCLUDED */'								; \
	echo '/* @endcond */' ;  )> $@

rsb-incoming.h: $(librsb_base_la_SOURCES) $(librsb_nounroll_la_SOURCES) $(top_srcdir)/rsb-incoming.grep
	( cat $(top_srcdir)/rsb_license_header.inc                                                                            ; \
       	echo '/* Collated internal headers of librsb -- for experimental use only only (not for usage). */'; echo	; \
	echo '#ifndef RSB_LIBRSB_INCOMING_H_INCLUDED'									; \
	echo '#define RSB_LIBRSB_INCOMING_H_INCLUDED'									; \
	echo '#include <rsb.h>'									; \
	( cat `git ls-files | $(GREP) '\.h'$$` $(srcdir)/rsb_util.h ) | $(GREP) -f $(top_srcdir)/rsb-incoming.grep | $(SED) -f $(top_srcdir)/rsb-incoming.sed | tac ; \
	echo '#endif /* RSB_LIBRSB_INCOMING_H_INCLUDED */'								; )> $@

nodist_include_HEADERS+=rsb-librsb-internals.h rsb-incoming.h
endif

include_HEADERS+=rsb.h $(EXTRAHEADERSFORTRAN)
if HAVE_M4
include_HEADERS+= rsb_types.h blas_sparse.h
endif

rsb_test_matops.h: rsb_mkl.h 
BINDISTFILES=$(include_HEADERS) $(EXTRA_LIBRARIES) $(dist_doc_DATA) doc/html doc/man/man3 doc/man/*rsb* librsb-config librsb.pc
rsbench_SOURCES=  rsbench.c $(RSBENCH_BUILT_SOURCES) rsb_genmm.c rsb_mmls.c \
		rsb_pr.c     rsb_pr.h \
		rsb_pcnt.c rsb_pcnt.h \
		rsb_failure_tests.h 	rsb_failure_tests.c 	\
		rsb_libspblas_tests.h 	rsb_libspblas_tests.c 
if HAVE_LIBRSBPP
rsbench_SOURCES+=rsb_dummy.cpp
endif

# order-only prerequisites for make rsbench
am_rsbench_OBJECTS: | $(include_HEADERS) $(noinst_LTLIBRARIES) $(lib_LTLIBRARIES) $(nodist_include_HEADERS)

# rsbenchxx_SOURCES=rsbenchxx.cpp
rsbench_DEPENDENCIES=librsb.la
if WANT_BLAS_SPARSE_MOD_INSTALL
else
if HAVE_FC
if HAVE_SPARSE_BLAS_INTERFACE
rsbench_DEPENDENCIES+=blas_sparse.mod
endif
rsbench_DEPENDENCIES+=rsb.mod
endif
endif
# rsbenchxx_DEPENDENCIES= rsbenchxx.cpp $(top_builddir)/librsb_nounroll.a $(top_builddir)/librsb_base.a # uhm
# if  WANT_V
# endif
#L1C_SOURCES=L1C.c
rsbench_CFLAGS=${RSB_RSBENCH_CFLAGS}
if HAVE_LIBRSBPP
$(abs_top_builddir)/librsbpp/librsbpp.la:
	$(MAKE) -j -C $(abs_top_builddir)/librsbpp/
endif
rsbench_LDADD=librsb.la ${RSB_RSBENCH_LIBS}

if HAVE_LIBRSBPP
rsbench_LINK=$(CXXLINK)
ot_LINK=$(CXXLINK)
sbtc_LINK=$(CXXLINK)
sbtf_LINK=$(CXXLINK)
else
rsbench_LINK=$(LINK)
ot_LINK=$(LINK)
sbtc_LINK=$(LINK)
sbtf_LINK=$(LINK)
endif
# rsbenchxx_LDADD=$(default_ldadd)
#L1C_LDADD=$(default_ldadd)
ot_LDADD=$(default_ldadd)

MATRICES=pd.mtx A.mtx vf.mtx ch.mtx cs.mtx lt.mtx ut.mtx us.mtx pg.mtx pd.mtx.bin pd.mtx.bad.bin pd.mtx.one.bin
RSB_MFILES=psb_mvsv_tester.m psbtf.m ot.m sbtc.m sbtg.m sbtf.m sbtg-types.m

EXTRA_DIST=$(LIBSOURCES) \
	$(librsb_spblas_la_SOURCES) \
	$(ot_SOURCES) \
	$(MATRICES) \
	$(RSB_MFILES) \
	autogen.sh \
	scripts/test.sh \
	$(top_srcdir)/rsb_license_header.inc \
	rsb-config.h.hin \
	rsb-incoming.sed rsb-incoming.grep \
	rsb.F90 \
	librsb.pc.in \
	$(EXTRASRCPROGSSFORTRAN)

if HAVE_SPARSE_BLAS_INTERFACE
EXTRA_DIST+=$(sbtc_SOURCES)
endif

sources: $(EXTRA_DIST)
if HAVE_LIBRSBPP
	$(MAKE) sources -C $(abs_top_builddir)/librsbpp
endif

LIBSOURCES=$(LIB_BUILT_SOURCES) \
	rsb_krnl_bcss.m4   \
	rsb_krnl_bcss_u.m4 \
	rsb_krnl_bcss_l.m4 \
	rsb_krnl_bcss_spsv_u.m4 \
	rsb_krnl_bcss_spmv_u.m4 \
	rsb_krnl_bcss_misc_u.m4 \
	rsb_krnl_bcoo_spmv_u.m4 \
	rsb_krnl_bcoo_macros.m4 \
	rsb_merge.m4 \
	rsb_util.m4 \
	rsb_ompio.m4 \
	rsb_mkl.m4 \
	rsb_spgemm_csr.m4 \
	rsb_spsum_misc.m4 \
	rsb_prec.m4 \
	rsb_spmv.m4 \
	rsb_krnl.m4 \
	rsb_bench.m4 \
	rsb_types.m4 \
	rsb_unroll.m4 \
	rsb_mergesort.m4 \
	rsb_permute.m4 \
	rsb_krnl_vb_macros.m4 rsb_misc.m4 rsb_test_matops.m4 mergesort_macros.m4 \
	rsb_krnl_bcss_macros.m4 rsb_krnl_macros.m4 do_unroll.m4 wisdom.m4 \
	rsb_config.m4 \
	rsb_krnl_linked_lists.m4 rsb_krnl_lb.m4 \
	rsb_krnl_vb.m4 \
	L1C.m4 \
	libspblas_macros.m4 \
	rsb_fortran_macros.m4 \
	sbtg-types.m4 \
	rsb_mod.m4 rsb_blas_sparse.m4 \
	rsb_libspblas.m4

ot_SOURCES_BUILT=ot.c ot-spmv_uauz.c ot-infty_norm.c ot-scale.c ot-spmv_uaua.c ot-spmv_unua.c \
	ot-spmv_uxua.c ot-spmv_sasa.c ot-spsv_uxua.c ot-spmv_sxsa.c ot-spsv_sxsx.c \
	ot-spsv.c ot-spmv.c ot-rowssums.c
ot_SOURCES=$(ot_SOURCES_BUILT)
if HAVE_LIBRSBPP
ot_SOURCES+=rsb_dummy.cpp
endif

if WANT_OCTAVE_TESTING
psbtf.F90: psbtf.m sbtg.m sbtg-types.m
	( cd $(abs_srcdir); abs_builddir="$(abs_builddir)" $(OCTAVE) $(OCTAVE_FLAGS) psbtf.m ) > psbtf.F90

psb_mvsv_tester.f90: psb_mvsv_tester.m sbtg.m sbtg-types.m
	( cd $(abs_srcdir); abs_builddir="$(abs_builddir)" $(OCTAVE) $(OCTAVE_FLAGS) psb_mvsv_tester.m ) > psb_mvsv_tester.f90

if HAVE_FORTRAN_EXAMPLES
sbtf.F90: sbtf.m sbtg.m sbtg-types.m
	( cd $(abs_srcdir); abs_builddir="$(abs_builddir)" $(OCTAVE) $(OCTAVE_FLAGS) sbtf.m ) > sbtf.F90
endif

if HAVE_SPARSE_BLAS_INTERFACE
sbtc.c: sbtc.m sbtg.m sbtg-types.m
	( cd $(abs_srcdir); abs_builddir="$(abs_builddir)" $(OCTAVE) $(OCTAVE_FLAGS) sbtc.m ) > sbtc.c
endif
else
if HAVE_FORTRAN_EXAMPLES
sbtf.F90: sbtf.m
	{ echo "PROGRAM MAIN" ; echo "PRINT *,'sorry, you did not install octave, so the octave based tester is disabled'"; echo "END PROGRAM MAIN" ; } > $@
endif

psb_mvsv_tester.f90: psb_mvsv_tester.m
	{ echo "PROGRAM MAIN" ; echo "PRINT *,'sorry, you did not install octave, so the octave based tester is disabled'"; echo "END PROGRAM MAIN" ; } > $@

psbtf.F90: psbtf.m
	{ echo "PROGRAM MAIN" ; echo "PRINT *,'sorry, you did not install octave, so the octave based tester is disabled'"; echo "END PROGRAM MAIN" ; } > $@

if HAVE_SPARSE_BLAS_INTERFACE
sbtc.c: sbtc.m
	echo "#include <stdio.h>" > $@
	echo "int main(void){printf(\"sorry, you did not install octave, so the octave based tester is disabled\\n\");return 0;}" >> sbtc.c
endif
endif

if HAVE_SPARSE_BLAS_INTERFACE
else
sbtc.c: sbtc.m
	echo "#include <stdio.h>" > $@
	echo "int main(void){printf(\"sorry, you did not want Sparse BLAS.\\n\");return 0;}" >> $@
endif
if HAVE_FORTRAN_EXAMPLES
else
sbtf.F90: sbtc.m
	{ echo "PROGRAM MAIN" ; echo "PRINT *,'sorry, you did not want Sparse BLAS.'"; echo "END PROGRAM MAIN" ; } > $@
endif

if WANT_OCTAVE_TESTING_AND_INT
$(ot_SOURCES_BUILT): ot.m sbtg.m sbtg-types.m
	$(OCTAVE) $(OCTAVE_FLAGS) ot.m $(WANT_ROW_UNLOOP_FACTORS) $(WANT_COLUMN_UNLOOP_FACTORS) $(WANT_MATRIX_ALL_OPS),$(WANT_MATRIX_ALL_META_OPS) $(WANT_MATRIX_OPS),$(WANT_MATRIX_ALL_META_OPS) main > ot.c
	for o in `echo $(WANT_MATRIX_ALL_OPS),$(WANT_MATRIX_ALL_META_OPS) | $(SED) "s/,/ /g"`  ; do $(OCTAVE) $(OCTAVE_FLAGS) ot.m $(WANT_ROW_UNLOOP_FACTORS) $(WANT_COLUMN_UNLOOP_FACTORS) $(WANT_MATRIX_ALL_OPS) $(WANT_MATRIX_OPS) $$o > ot-$$o.c  ; done
else
$(ot_SOURCES_BUILT):
	echo "#include <stdio.h>" > ot.c
	echo "int main(void){printf(\"sorry, you did not install octave, so the octave based tester is disabled\\n\");return 0;}" >> ot.c
	for o in `echo $(WANT_MATRIX_ALL_OPS),$(WANT_MATRIX_ALL_META_OPS) | $(SED) "s/,/ /g"` ; do echo "int f_$$o(void){return 0;}" > ot-$$o.c  ; done
endif

if HAVE_SED
rsb_strmif.h:
	( cat $(top_srcdir)/rsb_license_header.inc  ; \
	echo '/* @cond INNERDOC */' ; \
	echo '/* This file was generated by the Makefile */' ; \
	echo '#ifndef RSB_STRMIF_H_INCLUDED'; \
	echo '#define RSB_STRMIF_H_INCLUDED'; \
	echo '#include "rsb.h"' ; \
	echo '#include "rsb_common.h"' ; \
	echo '#include "rsb_do.h"' ; \
	echo 'rsb_err_t rsb__do_get_matrix_info_from_string(const struct rsb_mtx_t *matrix, const rsb_char_t *mis, void* info, size_t buflen);'; \
	echo 'rsb_err_t rsb__dump_flags(const rsb_flags_t cflags, const rsb_char_t *pre, const rsb_char_t *intra, const rsb_char_t *post);'; \
	echo '#endif /* RSB_STRMIF_H_INCLUDED */'; \
	echo '/* @endcond */' ; \
	) > $@

rsb_strmif.c: $(top_srcdir)/rsb.h $(top_srcdir)/rsb_common.h
	( cat $(top_srcdir)/rsb_license_header.inc  ; \
	echo '/* @cond INNERDOC */' ; \
	echo '/* This file was generated by the Makefile */' ; \
	echo '#include "rsb.h"' ; \
	echo '#include "rsb_common.h"' ; \
	echo '#include "rsb_do.h"' ; \
	echo 'rsb_err_t rsb__do_get_matrix_info_from_string(const struct rsb_mtx_t *matrix, const rsb_char_t *mis, void* info, size_t buflen)' ; \
	echo '{ rsb_err_t errval=RSB_ERR_BADARGS; if(!matrix || !mis || !info)goto err;' ; \
	$(GREP) -E '^(.define_|.) RSB_MIF_' $(top_srcdir)/rsb.h | $(SED) 's/^. /#define /g;s/=0x/0x/g' | $(SED) 's/[[:space:]][[:space:]]*/ /g;s/\/.*(//g;s/).*\///g;s/\/.*(//g;s/).*\///g;' | $(CUT) -d ' ' -f 2,4 | $(SED) 's/^\([[:graph:]][[:graph:]]*\) \([[:graph:]][[:graph:]]*\)/if(0 == strcmp(mis,"\1")){ errval = rsb__do_get_matrix_info(matrix,\1,info,buflen); goto done;}/g;'; \
	echo 'done:';	\
	echo 'return errval;';	\
	echo 'err: return RSB_ERR_GENERIC_ERROR;';	\
	echo '}'; \
	echo ''; \
	echo 'rsb_err_t rsb__dump_flags(const rsb_flags_t cflags, const rsb_char_t *pre, const rsb_char_t *intra, const rsb_char_t *post)' ; \
	echo '{ rsb_err_t errval=RSB_ERR_BADARGS; int ef=0; rsb_flags_t flags=cflags;'; \
	$(GREP) -h '^.define RSB_FLAG_[[:graph:]]*[[:space:]]*0x' $(top_srcdir)/rsb.h $(top_srcdir)/rsb_common.h | $(SED) 's/^.define //g;s/[[:space:]].*$$//g;s/\(.*\)/if(\1\&\&RSB_DO_FLAG_HAS(flags,\1)){RSB_STDOUT("%s%s",ef++?intra:pre,"\1");RSB_DO_FLAG_DEL(flags,\1);}/g'; \
	echo 'if(flags){RSB_STDOUT("%s0x%x",ef++?intra:pre,flags);}'; \
	echo 'if(ef)RSB_STDOUT("%s",post);';	\
	echo 'return errval;';	\
	echo '}'; \
	echo '/* @endcond */' ; \
	) > $@

rsb_stropts.c: $(top_srcdir)/rsb.h
	( cat $(top_srcdir)/rsb_license_header.inc  ; \
	echo '/* @cond INNERDOC */' ; \
	echo '/* This file was generated by the Makefile */' ; \
	echo '#include "rsb.h"' ; \
	echo '#include "rsb_common.h"' ; \
	echo 'rsb_err_t rsb__stropts_set(const rsb_char_t *opn, const rsb_char_t *arg)' ; \
	echo '{ rsb_err_t errval=RSB_ERR_NO_ERROR; if(!opn || !arg)goto err;' ; \
	$(GREP) -E '^(.define|.) RSB_IO_WANT_' $(top_srcdir)/rsb.h | $(SED) 's/^. /#define /g;s/=0x/0x/g' | $(GREP) 'rsb_int_t\|rsb_char_t\|rsb_real_t'| $(SED) 's/[[:space:]][[:space:]]*/ /g;s/\/.*(//g;s/).*\///g;' | $(CUT) -d ' ' -f 2,4,5 | $(SED) 's/^\([[:graph:]][[:graph:]]*\) \(const \)*\([[:graph:]][[:graph:]]*\)/if(0 == strcmp(opn,"\1")){ \2\3 RSB_DO_REINIT_SINGLE_VALUE_SET(\1,\&val,errval); goto done;}/g; s/\(rsb_char_t\*\)/\1 val = arg;/g;s/\(rsb_int_t\)/\1 val = rsb__util_atoi(arg);/g; s/\(rsb_real_t\)/\1 val = rsb__util_atof(arg);/g'; \
	echo 'done:';	\
	echo 'return errval;';	\
	echo 'err: return RSB_ERR_GENERIC_ERROR;';	\
	echo '}'; \
	echo '/* @endcond */' ; \
	) > $@
endif

.PHONY: feedback
feedback: rsbench$(EXEEXT)
	./rsbench$(EXEEXT) -O r 

M4_FLAGS=	\
		-I $(srcdir) \
		-D WANT_SPSM_DIAG_CHECK=$(WANT_SPSM_DIAG_CHECK) \
		-D WANT_HALFWORD_INDICES=$(WANT_HALFWORD_INDICES) \
		-D WANT_ROW_UNLOOP_FACTORS=$(WANT_ROW_UNLOOP_FACTORS) \
		-D WANT_COLUMN_UNLOOP_FACTORS=$(WANT_COLUMN_UNLOOP_FACTORS) \
		-D WANT_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR='$(WANT_SIMPLE_LOOP_UNROLL_DEFAULT_FACTOR)' \
		-D WANT_COLUMN_UNLOOP_FACTORS=$(WANT_COLUMN_UNLOOP_FACTORS) \
		-D WANT_MATRIX_OPS='$(WANT_MATRIX_OPS)' \
		-D WANT_TYPES='$(WANT_TYPES)' \
		-D WANT_MATRIX_ALL_OPS='$(WANT_MATRIX_ALL_OPS)' \
		-D WANT_MATRIX_ALL_TYPES='$(WANT_MATRIX_ALL_TYPES)' \
		-D WANT_MATRIX_ALL_META_OPS='$(WANT_MATRIX_ALL_META_OPS)' \
		-D WANT_MATRIX_BCSS_STORAGE='$(WANT_MATRIX_BCSS_STORAGE)' \
		-D WANT_MATRIX_BCOO_STORAGE='$(WANT_MATRIX_BCOO_STORAGE)' \
		-D WANT_MATRIX_VB_STORAGE='$(WANT_MATRIX_VB_STORAGE)' \
		-D WANT_LOOPING_KERNELS='$(WANT_LOOPING_KERNELS)' \
		-D WANT_MATRIX_STORAGE='$(WANT_MATRIX_STORAGE)' \
		-D WANT_LONG_IDX='$(WANT_LONG_IDX)' \
		-D WANT_IHI='$(WANT_IHI)'

# M4 macro files included by other macro files (rsb_krnl.m4 rsb_mergesort.m4 rsb_bench.m4 rsb_permute.m4 rsb_krnl_bcss_spsv_u.m4 rsb_krnl_bcss_spmv_u.m4 rsb_krnl_bcss_misc_u.m4 rsb_krnl_vb.m4 rsb_krnl_lb.m4  are all not):
RSB_KERNELS_MACROS=rsb_krnl_macros.m4 rsb_krnl_linked_lists.m4 wisdom.m4 rsb_misc.m4 rsb_krnl_bcss_macros.m4 \
	 rsb_krnl_bcoo_macros.m4 mergesort_macros.m4 rsb_krnl_vb_macros.m4
rsb_mergesort.c rsb_permute.c rsb_krnl_vb.c rsb_krnl_vb.h rsb_unroll.h rsb_krnl_lb.h rsb_krnl_lb.c rsb_krnl.c rsb_bench.c rsb_krnl_bcss_spsv_u.c rsb_krnl_bcss_spmv_u.c rsb_krnl_bcss_misc_u.c rsb_krnl_bcss.c rsb_krnl_bcoo.c rsb_krnl_bcoo_spmv_u.c rsb_krnl.c: $(RSB_KERNELS_MACROS)

blas_sparse.h: $(srcdir)/rsb_libspblas.h
	@if test ! -f $< ; then echo 'It seems like your "make" program version is inadequate here! GNU version works best here.'; fi
	cp -p $< $@ 

#rsb-config.h: config.h
#	( cat $(top_srcdir)/rsb_license_header.inc                                                                            ; \
#	echo '/* This header file is not intended to be included librsb programs: it is only for inspection. */'; \
#	echo '#ifndef RSB_CONFIG_H_INCLUDED'									; \
#	echo '#define RSB_CONFIG_H_INCLUDED'									; \
#	cat  $< | $(SED) 's/^#define /#define RSB_/g;s/ RSB_RSB_/ RSB_/g'   					; \
#	echo '/* #endif RSB_CONFIG_H_INCLUDED */'								; \
#	echo '#endif'												)> $@

# if WANT_CXX_TEST_RSBENCH
# rsbenchxx.cpp: $(srcdir)/rsbench.c
# 	cp $(srcdir)/rsbench.c rsbenchxx.cpp
# else
# endif

if HAVE_M4
sbtg-types.m: sbtg-types.m4 rsb_types.h
	$(M4) $(M4_FLAGS) $< > $@

.m4.c:
	$(M4) $(M4_FLAGS) $< > $@
if WANT_BUGGY_GCC_BUG_103995_AWK_WORKAROUND
	if test $@ == rsb_krnl_bcss_spmv_u.c ; then \
		mv $@ rsb_krnl_bcss_spmv_u.vanilla.c ; \
		awk -f $(top_srcdir)/scripts/librsb_gcc11.awk < rsb_krnl_bcss_spmv_u.vanilla.c > $@ ; \
		rm -f rsb_krnl_bcss_spmv_u.vanilla.c ; \
	fi
endif
if WANT_BUGGY_GCC_BUG_103995_COCCI_WORKAROUND
	if test $@ == rsb_krnl_bcss_spmv_u.c ; then \
		mv $@ rsb_krnl_bcss_spmv_u.vanilla.c ; \
		$(SPATCH) --local-includes --sp-file $(top_srcdir)/scripts/librsb_gcc11.cocci rsb_krnl_bcss_spmv_u.vanilla.c > gcc11.patch; \
		patch < gcc11.patch; \
		mv rsb_krnl_bcss_spmv_u.vanilla.c $@ ; \
		rm gcc11.patch; \
	fi
endif

.m4.h:
	$(M4) $(M4_FLAGS) -D ONLY_WANT_HEADERS=1 $< > $@

.m4.F90:
	$(M4) $(M4_FLAGS) $< > $@
endif

# AM_LDFLAGS=$(LIBS)

.PHONY: devsplinttest
devsplinttest:
	for f in $(librsb_nounroll_la_SOURCES) $(librsb_base_la_SOURCES) ; do splint -I`gcc -print-search-dirs | $(GREP) install: | $(SED) s/install:.//g`  -preproc -DHAVE_CONFIG_H $$f > $$f-splint.txt ;done

.PHONY: devtests
devtests:
if WANT_BUILD_DOC
	$(MAKE) tests -C $(abs_builddir)/doc
endif
	abs_top_srcdir="$(abs_top_srcdir)" $(BASH) $(srcdir)/scripts/devtests.sh
	PATH="`pwd`:$(PATH)" $(BASH) $(srcdir)/scripts/benchmark.sh --nrhs 1,2,4,8 $(srcdir)/pd.mtx $(srcdir)/cs.mtx
	@echo "	[*] dev test terminated successfully !"

.PHONY: tests
test: tests

README.md: README
	cp $< $@.tmp;
	$(GREP) '^\ - [A-Z]' README | while read; do ID=$${REPLY/ - /}; ID=$${ID//,/}; ID=$${ID,,}; ID=$${ID// /}; $(GREP) "^$${REPLY/ - /}" README -C0 | $(SED) "s/\(^$${REPLY/ - /}.*\)/## <a id=\"$${ID}\"><\/a>\1 ##/g;s/$${REPLY}/ - [$${REPLY/ - /}]\(#$${ID}\)/g;" -i $@.tmp; done;
	mv $@.tmp $@;

scripts/readme-tests.sh: README
	echo 'if test x"$${srcdir}" = x ; then srcdir=. ; fi' > $@
	LANG=C $(GREP) '^	*\(make \)**\./\(rsbench\|sbtc\|sbtf\)\|\(^	*test\> -f\)' $< | $(SED) 's/#.*$$//g;s/$$/ || exit 255/g' | $(SED) 's/A.mtx/$${srcdir}\/A.mtx/g' >> $@

.PHONY: mtests
mtests: rsbench$(EXEEXT) $(srcdir)/scripts/readme-tests.sh librsb-config
	srcdir="$(abs_srcdir)" $(BASH) -ex $(srcdir)/scripts/readme-tests.sh
	srcdir="$(abs_builddir)" $(BASH) $(srcdir)/scripts/doc-tests.sh
	if ! ./librsb-config --help ; then echo "Problem executing the librsb-config script!"; false; fi;
if WANT_BUILD_DOC
	if which lynx; then for f in doc/html/*.html ; do if lynx -dump $$f | $(GREP) '\\\(see\|code\)' ; then echo "Bad Doxygen generated in $$f" ; else true; fi ; done; fi
endif
	if ./rsbench$(EXEEXT)  -C | $(GREP) 'type char codes.*:*[SDCZ]' ; then cd examples ; $(MAKE) tests ; fi

.PHONY: gtests
gtests:
if HAVE_GTEST
	srcdir="$(abs_srcdir)" ./rsb_gtests$(EXEEXT)
else
	@echo " [*] Skipping tests based on Google Test (not detected at configure time)"
endif

.PHONY: check
check:	qtests

.PHONY: qtests
qtests:	all
	@echo " [*] beginning quick test..."
	$(MAKE) gtests
if HAVE_RUNNABLE_TESTS
	$(MAKE) mtests -C .
	./rsbench$(EXEEXT) -Q 30.0Q
else
	@echo " [!] Skipping several tests due to unusual configuration."
endif
if HAVE_LIBRSBPP
	$(MAKE) qtests -C @librsbpp@
endif
if HAVE_RSBLIB
	$(MAKE) qtests -C @rsblib@
endif
if HAVE_RSBTEST
	$(MAKE) qtests -C @rsbtest@
endif
	@echo " [*] quick test terminated successfully !"

.PHONY: qqtests
qqtests:	all
	@echo " [*] beginning quick quiet test..."
	$(MAKE) gtests         > /dev/null 2> /dev/null
	$(MAKE) mtests -C . > /dev/null 2> /dev/null
	./rsbench$(EXEEXT) -Q 30.0Q
if HAVE_LIBRSBPP
	$(MAKE) qtests -C @librsbpp@ > /dev/null 2> /dev/null
endif
if HAVE_RSBLIB
	$(MAKE) qtests -C @rsblib@ > /dev/null 2> /dev/null
endif
if HAVE_RSBTEST
	$(MAKE) qtests -C @rsbtest@ > /dev/null 2> /dev/null
endif
	@echo " [*] quick test terminated successfully !"

.PHONY: tests
tests:	$(EXTRA_LIBRARIES) ot$(EXEEXT) $(EXTRAPROGRAMSC) rsbench$(EXEEXT)
	$(MAKE) qtests -C .
if HAVE_RUNNABLE_TESTS
	srcdir="$(srcdir)" $(BASH) $(srcdir)/scripts/test.sh
	./ot$(EXEEXT)
if HAVE_SPARSE_BLAS_INTERFACE
	./sbtc$(EXEEXT)
endif
endif
	@echo "	[*] full test terminated successfully !"

.PHONY: btests
btests: 
	$(MAKE) hinfo.log
	$(BASH) $(srcdir)/scripts/test_configure.sh
	$(BASH) $(srcdir)/scripts/configure_for_debug.sh
	$(BASH) $(srcdir)/scripts/dist.sh
	$(BASH) $(srcdir)/scripts/dist_serial.sh
	$(BASH) $(srcdir)/scripts/dist_long_indices.sh
	$(BASH) $(srcdir)/scripts/dist_nofortran.sh
	$(BASH) $(srcdir)/scripts/build_optimized.sh CC=$(CC)
if HAVE_GPROF
	$(BASH) $(srcdir)/scripts/gprof.sh
else
	@echo "	[*] no gprof detected -- skipping test !"
endif
if HAVE_GCOV
	$(BASH) $(srcdir)/scripts/gcov.sh
	$(BASH) $(srcdir)/scripts/gcov_run.sh
else
	@echo "	[*] no gcov detected -- skipping test !"
endif
if HAVE_EFENCE_TEST
	PATH="$(prefix)/bin:./:$(PATH)" $(BASH) $(srcdir)/scripts/dev_brute_tests.sh
	@echo "	[*] brute force package testing terminated successfully !"
else
	@echo "	[*] no efence library detected -- skipping test !"
endif

.PHONY: itests
itests:
	$(MAKE) itests -C examples || { echo " [!] post-installation test failed !"; false; }
	@echo " [*] post-installation test terminated successfully !"

.PHONY: wc
wc:
	wc *.c
	wc *.h
	wc *.m4
	wc *.m
	wc *.F90
	cat *.m *.m4 *.c *.h *.F90 |wc
	cat `git ls-files| $(GREP) .c$$` | wc
	cat `git ls-files| $(GREP) .h$$` | wc
	cat `git ls-files| $(GREP) .m4$$` | wc
	cat `git ls-files| $(GREP) .m$$` | wc
	cat `git ls-files| $(GREP) .F90$$` | wc
	cat `git ls-files| $(GREP) '\.\(m\|m4\|c\|h\|F90\)$$'` | wc
if HAVE_LIBRSBPP
	cd @librsbpp@ && cat `git ls-files| $(GREP) '\.\(cpp\|hpp\|c\|h\)$$'` | wc
endif
if HAVE_RSBLIB
	cd @rsblib@ && cat `git ls-files| $(GREP) '\.\(cpp\|hpp\|c\|h\)$$'` | wc
endif
if HAVE_RSBTEST
	cd @rsbtest@ && cat `git ls-files| $(GREP) '\.\(cpp\|hpp\|c\|h\)$$'` | wc
endif

.PHONY: examples/make.sh
examples/make.sh:
	$(MAKE) && cd examples && $(MAKE) 

.PHONY: doxonly
doxonly: examples/make.sh
	$(MAKE) makedox -C doc

.PHONY: dox
dox: doxonly

.PHONY: install-exec-hook
install-exec-hook:
	$(MKDIR_P) "$(DESTDIR)$(docdir)"
if HAVE_PKGCONFIG_INSTALL
	$(MKDIR_P) "$(DESTDIR)$(libdir)/pkgconfig"
	$(INSTALL_DATA) librsb.pc "$(DESTDIR)$(libdir)/pkgconfig/"
endif

.PHONY: uninstall-hook
uninstall-hook:
if HAVE_PKGCONFIG_INSTALL
	if test -f "$(DESTDIR)$(libdir)/pkgconfig/librsb.pc" ; then  rm "$(DESTDIR)$(libdir)/pkgconfig/librsb.pc" ; fi
	if test -d "$(DESTDIR)$(libdir)/pkgconfig" ; then rmdir "$(DESTDIR)$(libdir)/pkgconfig" || true ; fi
endif

hinfo.log: all
	$(BASH) $(top_srcdir)/scripts/hinfo.sh 2>&1 | cat >  hinfo.log

# NOTE: The following target is only for experimental purposes.
#shared: $(top_builddir)/librsb.so
#$(top_builddir)/librsb.so: $(am_librsb_la_OBJECTS)
#	$(CC) -o $@ -shared $(am_librsb_la_OBJECTS)

if WANT_INTERNAL_HEADERS_INSTALL
rsb_stub.c: rsb_rsb.o
	( cat $(top_srcdir)/rsb_license_header.inc; nm rsb_rsb.o  | $(GREP) 'T rsb' | sort | $(SED) 's/.*rsb_/rsb_/g;s/^/int /g;s/$$/(){return 0;}/g'; echo  'int rsb_do_get_nnz_element(){return 0;} /* this is for sparsersb */' ; ) > rsb_stub.c && $(CC) -c rsb_stub.c                                                    
endif

if HAVE_FC
CLEANFILES = *.mod
endif

distclean-local:
	$(MAKE) realclean-local
	$(MAKE) cleanall-local

.DELETE_ON_ERROR:
