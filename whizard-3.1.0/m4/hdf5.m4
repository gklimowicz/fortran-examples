dnl hdf5.m4 -- Check for HDF5 library
dnl

# -------------------------------------------------------------
# Hdf5 library
# -------------------------------------------------------------
##   Provides a --with-hdf5=DIR option and minimum version check for
##   the HDF I/O library. Searches --with-hdf5, $HDF5DIR, and the
##   usual places for HDF5 headers and libraries.
##
##   On success, sets HDF5_CFLAGS, HDF5_LIBS, and #defines HAVE_HDF5.
##   Assumes package is optional unless overridden with $2=yes.
##
AC_DEFUN([WO_PROG_HDF5],
[dnl

AC_ARG_ENABLE([hdf5],
  [AS_HELP_STRING([--enable-hdf5],
    [build WHIZARD with HDF5 support [[yes]]])],
  [], [enable_hdf5="yes"])

if test "x$enable_hdf5" = "xyes"; then
  HAVE_HDF5=0
  AC_ARG_VAR(HDF5_DIR,[root directory of HDF5 installation])
  AC_ARG_WITH(hdf5,
    [AS_HELP_STRING([--with-hdf5=DIR],[root directory of HDF5 installation (default = HDF5_DIR)])],
   dnl action-if-given
  [with_hdf5=$withval
   if test "x${with_hdf5}" != "xyes"; then
      HDF5_PREFIX=$withval
   fi],
   dnl action-if-not-given
   [
   dnl This is "no" if the user did not specify --with-hdf5=foo
   with_hdf5=$withval
   dnl If $HDF5_DIR is set in the user's environment, then treat that
   dnl as though they had said --with-hdf5=$HDF5_DIR.
   if test "x${HDF5_DIR}" != "x"; then
      HDF5_PREFIX=${HDF5_DIR}
      with_hdf5=yes
   fi])

  dnl package requirement; if not specified, the default is to assume that
  dnl the package is optional
  
  dnl GNU-m4 ifelse documentation:
  dnl ifelse (string-1, string-2, equal, [not-equal])
  dnl If string-1 and string-2 are equal (character for character),
  dnl expands to the string in 'equal', otherwise to the string in
  dnl 'not-equal'.
  is_package_required=ifelse([$2], ,no, $2)
  
  if test "x${with_hdf5}" != "xno"; then
     if test -d "${HDF5_PREFIX}/lib"; then
        HDF5_LIBS="-L${HDF5_PREFIX}/lib -lhdf5"
	HDF5_FLIBS="-L${HDF5_PREFIX}/lib -lhdf5_fortran"
        HDF5_CXXLIBS="-L${HDF5_PREFIX}/lib -lhdf5_cpp"
     fi
     dnl If there is an "rpath" flag detected, append it to the various
     dnl LIBS vars.  This avoids hard-coding -Wl,-rpath, in case that is
     dnl not the right approach for some compilers.     
  fi

  if test "x$RPATHFLAG" != "x" && test -d "${HDF5_PREFIX}/lib"; then
        HDF5_LIBS="${HDF5_LIBS} ${RPATHFLAG}${HDF5_PREFIX}/lib"
        HDF5_FLIBS="${HDF5_FLIBS} ${RPATHFLAG}${HDF5_PREFIX}/lib"
        HDF5_CXXLIBS="${HDF5_CXXLIBS} ${RPATHFLAG}${HDF5_PREFIX}/lib"
  fi

  if test -d "${HDF5_PREFIX}/include"; then
        HDF5_CPPFLAGS="-I${HDF5_PREFIX}/include"
  fi

  wo_save_CFLAGS="$CFLAGS"
  wo_save_CPPFLAGS="$CPPFLAGS"
  wo_save_LDFLAGS="$LDFLAGS"
  wo_save_LIBS="$LIBS"

  CFLAGS="${HDF5_CPPFLAGS} ${CFLAGS}"
  CPPFLAGS="${HDF5_CPPFLAGS} ${CPPFLAGS}"
  LDFLAGS="${HDF5_LIBS} ${LDFLAGS}"

  AC_LANG_PUSH([C])
  AC_CHECK_HEADER([hdf5.h],[found_header=yes],[found_header=no])

  dnl ----------------------
  dnl  Minimum version check
  dnl ----------------------
  min_hdf5_version=ifelse([$1], ,1.8.0, $1)
  AC_MSG_CHECKING([for HDF5 version >= $min_hdf5_version])

  dnl Strip the major.minor.micro version numbers out of the min version string
  MAJOR_VER=`echo $min_hdf5_version | sed -e 's/^\([[0-9]]*\).*/\1/'`
  if test "x${MAJOR_VER}" = "x"; then
     MAJOR_VER=0
  fi
  MINOR_VER=`echo $min_hdf5_version | sed -e 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
  if test "x${MINOR_VER}" = "x"; then
     MINOR_VER=0
  fi
  MICRO_VER=`echo $min_hdf5_version | sed -e 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
  if test "x${MICRO_VER}" = "x"; then
     MICRO_VER=0
  fi

  dnl begin additional test(s) if header if available
  succeeded=no
  version_known=no

  if test "x${found_header}" = "xyes"; then
        min_version_succeeded=no
        hdf5_has_cxx=no
  
        dnl Test that HDF5 version is greater than or equal to the required min version.
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[
           @%:@include <hdf5.h>
               ]], [[
               @%:@if H5_VERS_MAJOR > $MAJOR_VER
               /* Sweet nibblets */
               @%:@elif (H5_VERS_MAJOR >= $MAJOR_VER) && (H5_VERS_MINOR >= $MINOR_VER) && (H5_VERS_RELEASE >= $MICRO_VER)
               /* Winner winner, chicken dinner */
               @%:@else
               @%:@  error HDF5 version is too old
               @%:@endif
	       printf("%i.%i.%i",H5_VERS_MAJOR, H5_VERS_MINOR, H5_VERS_RELEASE);
            ]])],[
                min_version_succeeded=yes
		version_known=yes
		HDF5_VERSION=`./conftest`
            ],[
                min_version_succeeded=no
		version_known=no
            ])
            if test "x$min_version_succeeded" = "xno"; then
                  AC_MSG_RESULT([no])
                  if test "x$is_package_required" = "xyes"; then
                      AC_MSG_ERROR([Your HDF5 library version does not meet the minimum version requirement (HDF5 >= $min_hdf5_version). Please use --with-hdf5 to specify the location of a valid installation.])
                  fi
            else
                  AC_MSG_RESULT([yes])
            fi
            dnl Check for -lhdf5
            AC_CHECK_LIB([hdf5],[H5Fopen],[found_library=yes],[found_library=no])
 
            dnl Test for the HDF5 C++ interface by trying to link a test code.
            AC_LANG_PUSH([C++])
 
            AC_MSG_CHECKING([if HDF5 C++ interface is present])
 
            dnl Using the C++ interface requires linking against both the C
            dnl and C++ libs.
            LIBS="${HDF5_LIBS} ${HDF5_CXXLIBS}"
 
            AC_LINK_IFELSE([AC_LANG_PROGRAM([[
              @%:@include <H5Cpp.h>
              @%:@ifndef H5_NO_NAMESPACE
              using namespace H5;
              @%:@endif
                  ]], [[
              H5std_string  fname("test.h5");
              H5File file (fname, H5F_ACC_TRUNC);
              ]])],[
                  hdf5_has_cxx=yes
              ],[
                  hdf5_has_cxx=no
              ])
 
             AC_LANG_POP([C++])
             dnl Not having the C++ interface doesn't disqualify us from using
             dnl the C interface.  We'll set a define if C++ is available, so
             dnl code can conditionally make use of it.
             if test "x$hdf5_has_cxx" = "xyes"; then
                   AC_MSG_RESULT([yes])
                   AC_DEFINE(HAVE_HDF5_CXX, 1, [Define if the HDF5 C++ interface is available])
	     else
                   AC_MSG_RESULT([no])
             fi
             succeeded=no
             if test "x$found_header" = "xyes" && test "x$min_version_succeeded" = "xyes" && test "x$found_library" = "xyes"; then
                   succeeded=yes
             fi 	    
  else
     AC_MSG_RESULT([no])
  fi
  if test "$version_known" = "yes"; then
    AC_MSG_CHECKING([the HDF5 version])
    AC_MSG_RESULT([$HDF5_VERSION])
    AC_SUBST([HDF5_VERSION])
  fi

  dnl Reset variables used by configure tests.
  CFLAGS="$wo_save_CFLAGS"
  CPPFLAGS="$wo_save_CPPFLAGS"
  LDFLAGS="$wo_save_LDFLAGS"
  LIBS="$wo_save_LIBS"

  if test "x$succeeded" = "xno"; then
     if test "x$is_package_required" = "xyes"; then
           AC_MSG_ERROR([HDF5 not found.  Try either --with-hdf5 or setting HDF5_DIR.])
     else	   
           AC_MSG_NOTICE([optional HDF5 library not found, or does not meet version requirements])
     fi
     HDF5_CFLAGS=""
     HDF5_CPPFLAGS=""
     HDF5_LIBS=""
     HDF5_FLIBS=""
     HDF5_CXXLIBS=""
     HDF5_PREFIX=""
  else
     HAVE_HDF5=1
     AC_DEFINE(HAVE_HDF5,1,[define if HDF5 is available])
     AC_SUBST(HDF5_CFLAGS)
     AC_SUBST(HDF5_CPPFLAGS)
     AC_SUBST(HDF5_LIBS)
     AC_SUBST(HDF5_FLIBS)
     AC_SUBST(HDF5_CXXLIBS)
     AC_SUBST(HDF5_PREFIX)
  fi

  AM_CONDITIONAL(HDF5_AVAILABLE,test x$HAVE_HDF5 = x1)

  HDF5_AVAILABLE_FLAG=$succeeded
  AC_SUBST([HDF5_AVAILABLE_FLAG])

  AC_MSG_CHECKING([for HDF5 support])
  if test "x$HAVE_HDF5" = "x0"; then
    enable_hdf5=no
    AC_MSG_RESULT([not found or disabled])
  else
    AC_MSG_RESULT([yes])
  fi
fi
])
