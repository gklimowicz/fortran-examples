ACLOCAL_AMFLAGS=-I m4

export EXEEXT

CXXFLAGS_SANITIZE=-fsanitize=null -fsanitize=bounds -fsanitize=return -fsanitize=bool -Wold-style-cast 
WARN_CXXFLAGS=-Wall -pedantic-errors -Wno-sign-compare -Wextra
ADD_COV_LDLIBS=
AM_CXXFLAGS= $(OPENMP_CXXFLAGS)
if WANT_DBG
AM_CXXFLAGS+= ${WARN_CXXFLAGS}
endif
if WANT_NO_BLOCKINGS
AM_CXXFLAGS+= -D RSBPP_PRESET_SPMM_BLOCKINGS=$(RSBPP_PRESET_SPMM_BLOCKINGS)
endif

include_HEADERS=
noinst_HEADERS=  rsbpp.hpp openmp_allocator.hpp  \
		 librsbpp.h \
		 rsbck.hpp rsbpp_coo.hpp rsbpp_csr.hpp
if WANT_PROGRAMS
noinst_PROGRAMS= rsbpp$(EXEEXT) rsbtt$(EXEEXT) rsbct$(EXEEXT)
else
noinst_PROGRAMS=
endif

LIBRSBPP_BUILT_SOURCE_FILES=rsbck.hpp rsbpp_coo.hpp rsbpp_csr.hpp
LIBRSBPP_BUILT_SOURCE_FILES_M4=rsbck.m4 rsbpp_coo.m4 rsbpp_csr.m4

EXTRA_DIST=rsbpp.m4 \
	$(LIBRSBPP_BUILT_SOURCE_FILES) \
	$(LIBRSBPP_BUILT_SOURCE_FILES_M4) \
	A.mtx G.mtx S.mtx test.sh

sources: $(LIBRSBPP_BUILT_SOURCE_FILES)

if HAVE_LIBRSB_M4_MACROS_DIR
rsbck.hpp: rsbck.m4
	$(M4) -I .. $< > $@

rsbpp_coo.hpp: rsbpp_coo.m4
	$(M4) -I .. $< > $@

rsbpp_csr.hpp: rsbpp_csr.m4
	$(M4) -I .. $< > $@
endif

cleanall:
if HAVE_LIBRSB_M4_MACROS_DIR
	rm -f $(LIBRSBPP_BUILT_SOURCE_FILES)
endif
if !HAVE_LIBRSB_M4_MACROS_DIR
	echo "Note: sources are not being cleaned up -- code generation is disabled."
endif

if WANT_PROGRAMS
if HAVE_GTEST
noinst_PROGRAMS+= rsbgt
rsbgt_LDADD= $(RPPGT_LDADD) librsbpp.la
rsbgt_DEPENDENCIES=librsbpp.la
rsbgt_SOURCES=rsbgt.cpp
endif

rsbpp_SOURCES=rsbpp.cpp
rsbpp_LDADD= librsbpp.la

rsbtt_SOURCES=rsbtt.cpp
rsbtt_LDADD= librsbpp.la

rsbct_SOURCES=rsbct.c
rsbct_LDADD= librsbpp.la
rsbct_LINK=$(CXXLINK)
endif

PHONIES = x g p4 p3 p s t b covrun coverage
.PHONY: all $(PHONIES) e wc
if WANT_PROGRAMS
all: rsbpp$(EXEEXT) rsbtt$(EXEEXT) rsbct$(EXEEXT)
else
all: librsbpp.la
endif

if WANT_CONVENIENCE_ONLY
noinst_LTLIBRARIES=librsbpp.la
else
lib_LTLIBRARIES=librsbpp.la
endif
librsbpp_la_SOURCES=rsbpp_coo.cpp rsbpp_csr.cpp
# librsbpp_la_LIBADD=-lstdc++

e:
	$(EDITOR) rsbpp.cpp

if WANT_PROGRAMS
x: rsbpp$(EXEEXT) g
	make s p p3 p4

if HAVE_GTEST
g: rsbgt
	./rsbgt
else
g:
	echo "Skipping tests based on Google Test (not detected at configure time)"
endif

p4: rsbpp$(EXEEXT)
	OMP_NUM_THREADS=4 ./rsbpp
p3: rsbpp$(EXEEXT)
	OMP_NUM_THREADS=3 ./rsbpp$(EXEEXT)
p: rsbpp$(EXEEXT)
	OMP_NUM_THREADS=2 ./rsbpp$(EXEEXT)
s: rsbpp$(EXEEXT)
	OMP_NUM_THREADS=1 ./rsbpp$(EXEEXT)
t: all
if HAVE_RUNNABLE_TESTS
	./rsbtt$(EXEEXT)
if !WANT_MINIMAL_LIBRSBPP
	if ! test -f G.mtx ; then cp -p $(abs_srcdir)/G.mtx . ; fi ; $(SHELL) $(abs_srcdir)/test.sh
endif
endif

b: rsbpp$(EXEEXT)
	OMP_NUM_THREADS=1       RSB_NUM_THREADS=1       ./rsbpp$(EXEEXT) C1000000m''4000''M''4000''I1r8,16tNorsF
	OMP_NUM_THREADS=`nproc` RSB_NUM_THREADS=`nproc` ./rsbpp$(EXEEXT) C1000000m''4000''M''4000''I1r8,16tNorsF

covrun:
	./rsbpp$(EXEEXT) **vqm1M10DLUFS
	make t
	make x

coverage:
	make clean
	make CXXFLAGS='-std=c++17 -pipe --coverage -O0 -g -Wcoverage-mismatch -Wall -pedantic  ${WARN_CXXFLAGS}' \
		ADD_COV_LDLIBS=-lgcov
	lcov  --directory `pwd` --zerocounters # deletes *.gcda
	make covrun
#	for f in *.o ; do gcov -f $${f/.o/}  ; done # STL
	gcov -f *.o
	rm -f *.info
	lcov --capture --directory `pwd` --output-file coverage.info
	lcov -a coverage.info  -o coverage-total.info
	genhtml coverage-total.info --highlight --legend --no-branch-coverage --function-coverage --branch-coverage  --output-directory coverage-info-dir
if HAVE_WHICH
	if which lynx  ; then lynx  -dump coverage-info-dir/index.html; fi
	if which links ; then links -dump coverage-info-dir/index.html | head -n 30 ; fi
endif
endif

if !WANT_PROGRAMS
$(PHONIES):
	@ echo "no program to test!"
endif

qtests: t g

wc:
	wc *.cpp *.hpp

.PHONY: realclean-local
realclean-local:
	true

.PHONY: distclean-local
distclean-local: realclean
	true

