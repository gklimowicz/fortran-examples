
subdir=examples

CLEANFILES=make.sh
noinst_SCRIPTS=make.sh bench.sh benchex.sh
EXPSOURCES_RSB=$(srcdir)/hello.c $(srcdir)/transpose.c $(srcdir)/power.c $(srcdir)/autotune.c $(srcdir)/backsolve.c
EXPSOURCES_RSB+=$(srcdir)/cplusplus.cpp
EXPSOURCES_RSB+=$(srcdir)/snippets.c
EXPSOURCES_NSB=$(srcdir)/hello-spblas.c $(srcdir)/io-spblas.c
EXPSOURCES_FEX=$(srcdir)/fortran.F90 $(srcdir)/fortran_rsb_fi.F90
EXPSOURCES=$(EXPSOURCES_RSB) $(EXPSOURCES_NSB) $(EXPSOURCES_FEX)
EXSTUFF=$(top_srcdir)/pd.mtx $(top_srcdir)/vf.mtx make.sh $(srcdir)/benchex.sh $(srcdir)/bench.sh $(EXPSOURCES) \
	   $(builddir)/Makefile_demo.am $(srcdir)/configure.ac
EXTRA_DIST=$(EXSTUFF) $(srcdir)/make.sh.in \
	   $(srcdir)/Makefile_demo.am.in
LIBRSB_LIB=$(abs_top_builddir)/librsb.la

LIBS+=$(FCLIBS)
if HAVE_LIBRSBPP
EXAMPLES_LINK=$(CXXLINK)
else
EXAMPLES_LINK=$(LINK)
endif

if HAVE_C_EXAMPLES
hello_DEPENDENCIES=$(LIBRSB_LIB)
autotune_DEPENDENCIES=$(LIBRSB_LIB)
backsolve_DEPENDENCIES=$(LIBRSB_LIB)
transpose_DEPENDENCIES=$(LIBRSB_LIB)
power_DEPENDENCIES=$(LIBRSB_LIB)
snippets_DEPENDENCIES=$(LIBRSB_LIB)
hello_LDADD=$(default_ldadd)
autotune_LDADD=$(default_ldadd)
backsolve_LDADD=$(default_ldadd)
transpose_LDADD=$(default_ldadd)
power_LDADD=$(default_ldadd)
snippets_LDADD=$(default_ldadd)
hello_LINK=$(EXAMPLES_LINK)
autotune_LINK=$(EXAMPLES_LINK)
backsolve_LINK=$(EXAMPLES_LINK)
transpose_LINK=$(EXAMPLES_LINK)
power_LINK=$(EXAMPLES_LINK)
snippets_LINK=$(EXAMPLES_LINK)
EXTRAPROGRAMSC=hello$(EXEEXT) snippets$(EXEEXT) transpose$(EXEEXT) power$(EXEEXT) autotune$(EXEEXT) backsolve$(EXEEXT)
if HAVE_SPARSE_BLAS_INTERFACE
io_spblas_DEPENDENCIES=$(LIBRSB_LIB)
io_spblas_LDADD=$(default_ldadd)
hello_spblas_LDADD=$(default_ldadd)
hello_spblas_DEPENDENCIES=$(LIBRSB_LIB)
io_spblas_LINK=$(EXAMPLES_LINK)
hello_spblas_LINK=$(EXAMPLES_LINK)
EXTRAPROGRAMSSPBLAS=hello-spblas$(EXEEXT) io-spblas$(EXEEXT)
else
EXTRAPROGRAMSSPBLAS=
endif
else
EXTRAPROGRAMSSPBLAS=
EXTRAPROGRAMSC=
endif

# NOTE: see main Makefile's default_ldadd 
default_ldadd=-L$(top_builddir) $(LIBRSB_LIB)
default_ldadd+=${RSB_RSBENCH_LIBS}
AM_CFLAGS= -I$(top_builddir) -I$(top_srcdir)
AM_FCFLAGS=-I$(top_builddir) -I$(top_srcdir)
tests:	$(noinst_PROGRAMS) all
	if test $(abs_top_builddir) != $(abs_top_srcdir) ; then cp $(abs_top_srcdir)/pd.mtx $(abs_top_srcdir)/vf.mtx $(abs_srcdir) ; fi
if HAVE_RUNNABLE_TESTS
	(                                                       PATH="$(abs_top_builddir):$(PATH)" $(BASH) $(abs_srcdir)/bench.sh;   )
endif
	for mf in pd.mtx vf.mtx ; do if test -f $(abs_builddir)/$$mf ; then true; else cp -p $(abs_top_srcdir)/$$mf $(abs_builddir)/$$mf ; fi; done
	for ii in $(noinst_PROGRAMS) ; do echo $(abs_builddir)/$$ii ; if  $(abs_builddir)/$$ii ; then true ; else exit -1 ;fi ; done
if HAVE_C_EXAMPLES
	./autotune$(EXEEXT)  $(abs_top_srcdir)/pd.mtx
endif

if HAVE_PKGCONFIG_INSTALL
PKG_CONFIG_PATH="$(DESTDIR)$(libdir)/pkgconfig"
endif

.PHONY: itests
itests:
	test -w $(DESTDIR)$(docdir)/examples/make.sh || { echo "Did you install librsb before running 'make itests' ?"; false; }
if HAVE_RUNNABLE_TESTS
	(                                                       PATH="$(bindir):$(PATH)" $(BASH) $(DESTDIR)$(docdir)/examples/bench.sh;   )
	( cd $(DESTDIR)$(docdir)/examples/ && test -f pd.mtx && PATH="$(bindir):$(PATH)" $(BASH) $(DESTDIR)$(docdir)/examples/benchex.sh; )
endif
if WANT_CPP_EXAMPLES
if HAVE_AUTORECONF
	( $(OPENBSD_AUTORECONF_EXPORTS) cd $(DESTDIR)$(docdir)/examples/ && touch NEWS README AUTHORS ChangeLog && autoreconf -i && PATH="$(bindir):$(PATH)" ./configure --disable-dependency-tracking --host="$(host)" CC="`$(bindir)/librsb-config --cc`" CXX="`$(bindir)/librsb-config --cxx`" FC="`$(bindir)/librsb-config --fc`" MAKE="$(MAKE)" && $(MAKE) clean && $(MAKE) && $(MAKE) test && $(MAKE) distclean; )
endif
endif
if WANT_CPP_EXAMPLES
if HAVE_CMP
if HAVE_RSBLIB
	cd ${srcdir}/../rsblib/examples; \
		for f in *.cpp; do $(CMP) $$f $(DESTDIR)$(docdir)/examples/$$f; done
endif
	for f in $(DIST_SOURCES); do $(CMP) $(abs_srcdir)/$$f $(DESTDIR)$(docdir)/examples/$$f; done
endif
endif
	cd $(DESTDIR)$(docdir)/examples && WANT_CLEANUP=true PKG_CONFIG_PATH="$(PKG_CONFIG_PATH)" srcdir='' builddir='' LIBRSB_CONFIG="$(bindir)/librsb-config" $(BASH) make.sh

if HAVE_FORTRAN_EXAMPLES
fortran_rsb_fi_DEPENDENCIES=$(LIBRSB_LIB)
fortran_rsb_fi_LDADD=$(default_ldadd)
fortran_rsb_fi_SOURCES=fortran_rsb_fi.F90
fortran_rsb_fi_LIBTOOLFLAGS= --preserve-dup-deps
if HAVE_SPARSE_BLAS_INTERFACE
fortran_DEPENDENCIES=$(LIBRSB_LIB)
fortran_LDADD=$(default_ldadd)
fortran_SOURCES=fortran.F90
fortran_LIBTOOLFLAGS= --preserve-dup-deps
EXTRAPROGRAMSFORTRAN=fortran$(EXEEXT)
else
EXTRAPROGRAMSFORTRAN=
endif
if WANT_FORTRAN_LINKER
fortran_rsb_fi_LINK=$(FCLINK)
fortran_LINK=$(FCLINK)
else
fortran_rsb_fi_LINK=$(EXAMPLES_LINK)
fortran_LINK=$(EXAMPLES_LINK)
endif
EXTRAPROGRAMSFORTRAN+=fortran_rsb_fi$(EXEEXT)
else
EXTRAPROGRAMSFORTRAN=
endif

noinst_PROGRAMS=$(EXTRAPROGRAMSC) $(EXTRAPROGRAMSSPBLAS) $(EXTRAPROGRAMSFORTRAN)
if HAVE_LIBRSBPP
if WANT_CPP_EXAMPLES
# Adding a C++ example implies CXXLINK to be set.
noinst_PROGRAMS+=cplusplus$(EXEEXT)
cplusplus_SOURCES=cplusplus.cpp
cplusplus_LDADD=$(default_ldadd)
endif
endif

all: $(EXTRA_DIST)

pd.mtx: $(top_srcdir)/pd.mtx
	cp $< $@ 

vf.mtx: $(top_srcdir)/vf.mtx
	cp $< $@ 

install-data-local:
	$(MKDIR_P) "$(DESTDIR)$(docdir)"
	$(MKDIR_P) "$(DESTDIR)$(docdir)/examples/"
	$(INSTALL_DATA) $(EXSTUFF) "$(DESTDIR)$(docdir)/examples/"
	chmod +x "$(DESTDIR)$(docdir)/examples/make.sh"
	mv "$(DESTDIR)$(docdir)/examples/Makefile_demo.am" $(DESTDIR)$(docdir)/examples/Makefile.am

uninstall-local:
if HAVE_BASENAME
	for f in $(EXSTUFF) ; do b="`$(BASENAME) $$f`"; if test -f "$(DESTDIR)$(docdir)/examples/"$$b ; then  rm "$(DESTDIR)$(docdir)/examples/"$$b ; fi ; done
else
	for f in $(EXSTUFF) ; do if test -f "$(DESTDIR)$(docdir)/examples/"$$f ; then  rm "$(DESTDIR)$(docdir)/examples/"$$f ; fi ; done # broken on out-of-dir
endif
	if test -d "$(DESTDIR)$(docdir)/examples" ; then rmdir "$(DESTDIR)$(docdir)/examples" || true ; fi
	if test -f "$(DESTDIR)$(docdir)/examples/Makefile.am" ; then rm "$(DESTDIR)$(docdir)/examples/Makefile.am" || true ; fi
	if test -d "$(DESTDIR)$(docdir)" ; then rmdir "$(DESTDIR)$(docdir)" || true ; fi

.PHONY: realclean
realclean:
	rm -f *.eps *.rpr *.tex *.txt *.mtx
	$(MAKE) clean
