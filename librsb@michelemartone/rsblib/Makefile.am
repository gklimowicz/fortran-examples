ACLOCAL_AMFLAGS=-I m4
WARN_CXXFLAGS=-Wall -pedantic-errors -Wno-sign-compare -Wextra
AM_CXXFLAGS=
AM_CXXFLAGS+= $(OPENMP_CXXFLAGS)
if WANT_DBG
AM_CXXFLAGS+= ${WARN_CXXFLAGS}
endif
SUBDIRS= .
if WANT_PROGRAMS
SUBDIRS+=examples
endif
EXTRA_DIST=A.mtx T.mtx Doxyfile
if WANT_PROGRAMS
noinst_PROGRAMS=rsb$(EXEEXT) rgt$(EXEEXT) rgt_ne$(EXEEXT) rgt_rv$(EXEEXT)
rsb_SOURCES=rsb.cpp
endif
include_HEADERS=rsb.hpp

.PHONY: tests gtests
if WANT_PROGRAMS
tests: rsb$(EXEEXT) rgt$(EXEEXT) examples
	$(MAKE) -C examples tests
if HAVE_RUNNABLE_TESTS
	RSBP_QUIET=1 ./$< $(abs_top_srcdir)/T.mtx
endif
	@echo " [*] tests terminated successfully !"

gtests: rgt$(EXEEXT) rgt_ne$(EXEEXT) rgt_rv$(EXEEXT)
if HAVE_GTEST
	./rgt$(EXEEXT)
	./rgt_ne$(EXEEXT)
	./rgt_rv$(EXEEXT)
endif
else
tests:
	@echo " [*] Skipping tests based on Google Test (not detected at configure time)"
gtests:
	@echo " [*] Skipping tests based on Google Test (not detected at configure time)"
endif

.PHONY: qtests
qtests: gtests tests

.PHONY: dox doc
if HAVE_DOXYGEN
$(rsblib_MANS) man: dox
rsblib_MANS=$(abs_builddir)/man/man3/R*3 \
	  $(abs_builddir)/man/man3/rsb.hpp.3*
if WANT_RSBLIB_AS_SUBDIR
man3_MANS=
else
man3_MANS=$(rsblib_MANS)
endif
$(abs_builddir)/man $(abs_builddir)/html: dox
EXTRA_DIST+=man html
dox: doc
$(abs_builddir)/Doxyfile:
	if test $(abs_builddir)/Doxyfile != $(abs_srcdir)/Doxyfile ; then \
		sed "s|^OUTPUT_DIRECTORY *=|OUTPUT_DIRECTORY=$(abs_builddir)|g" $(abs_srcdir)/Doxyfile > $@ ; fi
doc: $(abs_builddir)/Doxyfile $(abs_builddir)/html/index.html
$(abs_builddir)/html/index.html:
	cd $(abs_srcdir) && $(DOXYGEN) $(abs_builddir)/Doxyfile
	cd $(abs_builddir)/man/man3/ && $(MKDIR_P) discarded && mv *.3 discarded && mv discarded/rsb.hpp.3 discarded/RsbLib.3 discarded/RsbMatrix.3 .
	rm -fR $(abs_builddir)/man/man3/discarded/
endif

.PHONY: e
e:
	$(EDITOR) rsb.cpp

.PHONY: wc
wc:
	wc *.cpp

.PHONY: t x
t x: tests$(EXEEXT) gtests$(EXEEXT)
if WANT_PROGRAMS
rgt_LDADD=$(RGT_LDADD)
rgt_SOURCES=rgt.cpp
rgt_CXXFLAGS=${AM_CXXFLAGS}
rgt_DEPENDENCIES=rsb.hpp
rgt_ne_LDADD=$(RGT_LDADD)
rgt_ne_SOURCES=rgt.cpp
rgt_ne_CXXFLAGS=-DRSBP_NOTHROW=1 ${AM_CXXFLAGS}
rgt_ne_DEPENDENCIES=rsb.hpp
rgt_rv_LDADD=$(RGT_LDADD)
rgt_rv_SOURCES=rgt.cpp
rgt_rv_CXXFLAGS=-DRSBP_NOTHROW=1 -DRSBP_WANT_REV=1 ${AM_CXXFLAGS}
rgt_rv_DEPENDENCIES=rsb.hpp
endif

.PHONY: covrun
covrun:
	$(MAKE) x

.PHONY: coverage
coverage:
	$(MAKE) clean
	$(MAKE) CXXFLAGS='-std=c++17 -pipe --coverage -O0 -g -Wcoverage-mismatch -Wall -pedantic  ${WARN_CXXFLAGS}'
	lcov  --directory `pwd` --zerocounters # deletes *.gcda
	$(MAKE) covrun
#	for f in *.o ; do gcov -f $${f/.o/}  ; done # STL
	gcov -f *.o *.hpp *.cpp
	rm -f *.info
	lcov --capture --directory `pwd` --output-file coverage.info --no-external
	lcov -a coverage.info  -o coverage-total.info
	genhtml coverage-total.info --highlight --legend --no-branch-coverage --function-coverage --branch-coverage  --output-directory coverage-info-dir
if HAVE_WHICH
	if which lynx  ; then lynx  -dump coverage-info-dir/index.html; fi
	if which links ; then links -dump coverage-info-dir/index.html | head -n 30 ; fi
endif

install-data-local:
	$(MAKE) -C examples $@

uninstall-local:
	$(MAKE) -C examples $@

.PHONY: realclean-local
realclean-local: clean
	rm -f *.eps *.bin examples/*.bin

.PHONY: distclean-local
distclean-local: realclean-local
	true
