dnl tirpc.m4 -- Check for tirpc RPC header and library
dnl

### SunRPC is getting obsolescent/obsolete and has been replaced
### by the tirpc headers and library.    
AC_DEFUN([WO_PROG_TIRPC],
[dnl
PKG_PROG_PKG_CONFIG([0.9.0])
PKG_CHECK_MODULES([TIRPC], [libtirpc],
                  [LIBTIRPC="${TIRPC_LIBS}"
                     RPC_CFLAGS="${TIRPC_CFLAGS} ${TIRPC_LIBS}"
		     AC_MSG_NOTICE([for StdHEP legacy code: using libtirpc headers and library])
		     AC_MSG_NOTICE([                        with $TIRPC_CFLAGS $TIRPC_LIBS])],
		  [  RPC_CFLAGS=
		     AC_MSG_NOTICE([for StdHEP legacy code: using SunRPC headers and library])])
AC_SUBST([RPC_CFLAGS])
])
