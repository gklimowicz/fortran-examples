subdir=examples

EXTRA_DIST=bench.gpi
DEFAULT_INCLUDES= -I$(top_builddir)

noinst_PROGRAMS=
if WANT_EXAMPLES
noinst_PROGRAMS+= \
	assemble$(EXEEXT) \
	autotune$(EXEEXT) \
	bench$(EXEEXT) \
	build$(EXEEXT) \
	misc$(EXEEXT) \
	twonnz$(EXEEXT) \
	example$(EXEEXT) \
	mtx2bin$(EXEEXT) \
	render$(EXEEXT) \
	span$(EXEEXT)
endif

tests: $(noinst_PROGRAMS)
if HAVE_RUNNABLE_TESTS
	if test ! -f ../A.mtx ; then cp $(abs_top_srcdir)/A.mtx ..; fi
if WANT_EXAMPLES
if !WANT_CUSTOM_NUMERICAL_TYPES
	for f in $(noinst_PROGRAMS) ; do ./$$f || exit -1 ; done
	RSB_USER_SET_MEM_HIERARCHY_INFO='L1:8/128/2M' ./bench$(EXEEXT)
	./render$(EXEEXT) ../A.mtx
	rm -f A.mtx.bin && ./mtx2bin$(EXEEXT) ../A.mtx A.mtx.bin     && test -f A.mtx.bin
	rm -f A.mtx.bin && ./mtx2bin$(EXEEXT) ../A.mtx A.mtx.bin S   && test -f A.mtx.bin
	rm -f A.mtx.bin && ./mtx2bin$(EXEEXT) ../A.mtx A.mtx.bin D   && test -f A.mtx.bin
	rm -f A.mtx.bin && ./mtx2bin$(EXEEXT) ../A.mtx A.mtx.bin C   && test -f A.mtx.bin
	rm -f A.mtx.bin && ./mtx2bin$(EXEEXT) ../A.mtx A.mtx.bin Z   && test -f A.mtx.bin
endif
	./mtx2bin$(EXEEXT) ../A.mtx A.mtx.bin '?' ; test $$? != 0
	./mtx2bin$(EXEEXT) ../non-existent.mtx A.mtx.bin 'S' ; test $$? != 0
endif
endif

default_ldadd=-L$(top_builddir)
AM_CXXFLAGS= -I$(top_srcdir) $(OPENMP_CXXFLAGS)

if WANT_EXAMPLES
assemble_LDADD=
assemble_SOURCES=assemble.cpp
assemble_DEPENDENCIES=$(top_srcdir)/rsb.hpp

autotune_LDADD=
autotune_SOURCES=autotune.cpp
autotune_DEPENDENCIES=$(top_srcdir)/rsb.hpp

bench_LDADD=
bench_SOURCES=bench.cpp
bench_DEPENDENCIES=$(top_srcdir)/rsb.hpp

build_LDADD=
build_SOURCES=build.cpp
build_DEPENDENCIES=$(top_srcdir)/rsb.hpp

misc_LDADD=
misc_SOURCES=misc.cpp
misc_DEPENDENCIES=$(top_srcdir)/rsb.hpp

twonnz_LDADD=
twonnz_SOURCES=twonnz.cpp
twonnz_DEPENDENCIES=$(top_srcdir)/rsb.hpp

example_LDADD=
example_SOURCES=example.cpp
example_DEPENDENCIES=$(top_srcdir)/rsb.hpp

if HAVE_STDCPPFS
mtx2bin_LDADD= -lstdc++fs
else
mtx2bin_LDADD=
endif
mtx2bin_SOURCES=mtx2bin.cpp
mtx2bin_DEPENDENCIES=$(top_srcdir)/rsb.hpp

render_LDADD=
render_SOURCES=render.cpp
render_DEPENDENCIES=$(top_srcdir)/rsb.hpp

span_LDADD=
span_SOURCES=span.cpp
span_DEPENDENCIES=$(top_srcdir)/rsb.hpp
endif

if WANT_RSBLIB_AS_SUBDIR
RSBLIB_EXAMPLES_DEST_DIR="$(DESTDIR)$(docdir)/../librsb"
else
RSBLIB_EXAMPLES_DEST_DIR="$(DESTDIR)$(docdir)"
endif

install-data-local:
	$(MKDIR_P) "$(RSBLIB_EXAMPLES_DEST_DIR)"
	$(MKDIR_P) "$(RSBLIB_EXAMPLES_DEST_DIR)/examples/"
	for f in $(DIST_SOURCES) ; do $(INSTALL_DATA) $(srcdir)/$$f "$(RSBLIB_EXAMPLES_DEST_DIR)/examples/"$$f; done
if WANT_RSBLIB_AS_SUBDIR
	rmdir "$(DESTDIR)$(docdir)" || true; # workaround to avoid leaving unnecessary empty doc/rsblib dir
endif

uninstall-local:
if WANT_RSBLIB_AS_SUBDIR
	$(MKDIR_P) "$(RSBLIB_EXAMPLES_DEST_DIR)";
	ls -l "$(RSBLIB_EXAMPLES_DEST_DIR)";
endif
	for f in $(DIST_SOURCES) ; do if test -f "$(RSBLIB_EXAMPLES_DEST_DIR)/examples/"$$f ; then  rm "$(RSBLIB_EXAMPLES_DEST_DIR)/examples/"$$f ; fi ; done
	if test -d "$(RSBLIB_EXAMPLES_DEST_DIR)/examples" ; then rmdir "$(RSBLIB_EXAMPLES_DEST_DIR)/examples" || true ; fi
	if test -d "$(RSBLIB_EXAMPLES_DEST_DIR)" ; then rmdir "$(RSBLIB_EXAMPLES_DEST_DIR)" || true ; fi
if WANT_RSBLIB_AS_SUBDIR
	if test -d "$(DESTDIR)$(docdir)" ; then rmdir "$(DESTDIR)$(docdir)" || true ; fi;
endif
