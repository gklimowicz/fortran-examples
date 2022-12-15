dnl m4 -D M4_RSB_MOD_HDR=... -D M4_RSB_MOD_ROOTDIR=$(prefix) -D M4_RSB_MOD_EXTRA_DOC=... -D M4_RSB_MOD_EXTRA_CMD=...
changecom(`')dnl
#%Module 
M4_RSB_MOD_HDR

module-whatis "Libraries:optimized kernels (Sparse BLAS):Recursive Sparse Blocks library"
set WWWDoc "http://librsb.sourceforge.net/"

set ROOTDIR			"M4_RSB_MOD_ROOTDIR"
setenv LIBRSB_CONFIG		"$ROOTDIR/bin/librsb-config"
prepend-path PATH        	"$ROOTDIR/bin"
prepend-path MANPATH        	"$ROOTDIR/share/man"
prepend-path LD_LIBRARY_PATH 	"$ROOTDIR/lib/"
M4_RSB_MOD_EXTRA_CMD
proc ModulesHelp {} {
  set LIBRSB_CONFIG		{${LIBRSB_CONFIG}}
  global WWWDoc ROOTDIR
  set manfiles [exec ls -C "$ROOTDIR/share/man/man3"]
  set manpages [subst [regsub -all {.3} $manfiles ""]]

  puts stderr "
Recursive Sparse Blocks library, a Sparse BLAS library.

To test the examples:
	cp -fR $ROOTDIR/share/doc/librsb/examples/ ~/rsb-examples && cd ~/rsb-examples && ./make.sh

Local HTML documentation installed in:
	$ROOTDIR/share/doc/librsb/html/index.html

Local Man pages:
        $manpages

For compilation and linkage flags, use the LIBRSB_CONFIG environment variable:
	CFLAGS+=\\ \$( ${LIBRSB_CONFIG} --I_opts)
	CXXFLAGS+=\\ \$( ${LIBRSB_CONFIG} --I_opts)
	FCFLAGS+=\\ \$( ${LIBRSB_CONFIG} --I_opts)
       	LDFLAGS+=\\ \$( ${LIBRSB_CONFIG} --ldflags --extra_libs)

...and LIBRSB_CONFIG is set to:
       	M4_RSB_MOD_ROOTDIR/bin/librsb-config

Official web site:
	$WWWDoc

M4_RSB_MOD_EXTRA_DOC
"
  return 0
}

