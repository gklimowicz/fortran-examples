dnl lhapdf.m4 -- checks for LHAPDF library
dnl

### Determine paths to LHAPDF components
### Sets LDFLAGS_LHAPDF and the conditional LHAPDF_AVAILABLE if successful
### Also: LHAPDF_ROOT LHAPDF_VERSION LHAPDF_PDFSETS_PATH
AC_DEFUN([WO_PROG_LHAPDF],
[dnl
AC_REQUIRE([AC_PROG_FC])

AC_ARG_ENABLE([lhapdf],
  [AS_HELP_STRING([--enable-lhapdf],
    [enable LHAPDF for structure functions [[yes]]])],
  [], [enable_lhapdf="yes"])

if test "$enable_lhapdf" = "yes"; then
  if test -n "$LHAPDF_DIR"; then
    wo_lhapdf_config_path=$LHAPDF_DIR/bin:$PATH
  else
    wo_lhapdf_config_path=$PATH
  fi
  AC_PATH_PROG([LHAPDF], [lhapdf], [no], 
    [$wo_lhapdf_config_path])

  if test "$LHAPDF" != "no"; then  
    AC_PATH_PROG([LHAPDF_CONFIG], [lhapdf-config], [no], 
      [$wo_lhapdf_config_path])

    if test "$LHAPDF_CONFIG" != "no"; then

       AC_CACHE_CHECK([the LHAPDF version],
          [wo_cv_lhapdf_version],
          [dnl,
             wo_cv_lhapdf_version=[`$LHAPDF_CONFIG --version | $SED -e 's/.*\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/'`]
       ])
       LHAPDF_FULL_VERSION="$wo_cv_lhapdf_version"
       AC_SUBST([LHAPDF_FULL_VERSION])

       AC_CACHE_CHECK([the major version],
       [wo_cv_lhapdf_major_version],
       [wo_cv_lhapdf_major_version=[`echo $wo_cv_lhapdf_version | $SED -e 's/\([0-9][0-9]*\)\..*/\1/'`]
       ])
       LHAPDF_MAJOR_VERSION="$wo_cv_lhapdf_major_version"
       AC_SUBST([LHAPDF_MAJOR_VERSION])

       AC_MSG_CHECKING([the LHAPDF pdfsets path])
       LHAPDF_PDFSETS_PATH=`$LHAPDF_CONFIG --datadir`
       AC_MSG_RESULT([$LHAPDF_PDFSETS_PATH])

       AC_MSG_CHECKING([the standard PDF sets])
       if test -f "$LHAPDF_PDFSETS_PATH/CT10/CT10.info" -a -f "$LHAPDF_PDFSETS_PATH/CT10/CT10_0000.dat" -a -f "$LHAPDF_PDFSETS_PATH/cteq6l1/cteq6l1.info" -a -f "$LHAPDF_PDFSETS_PATH/cteq6l1/cteq6l1_0000.dat"; then
          AC_MSG_RESULT([ all standard PDF sets installed])
       else
          AC_MSG_RESULT([ not all standard PDF sets installed])
          AC_MSG_NOTICE([error: *************************************************************])
          AC_MSG_NOTICE([error: LHAPDF standard PDF sets not installed, please install these ])
          AC_MSG_NOTICE([error:    PDF sets: cteq6l1, CT10.                                  ])
          AC_MSG_NOTICE([error: *************************************************************])
          enable_lhapdf="no"
          AC_MSG_CHECKING([for LHAPDF])
          AC_MSG_RESULT([(disabled)])
       fi
    else
          AC_MSG_NOTICE([error: *****************************************************])
          AC_MSG_NOTICE([error: LHAPDF configure scripts not found or not executable ])
          AC_MSG_NOTICE([error: *****************************************************])
          enable_lhapdf="no"
          AC_MSG_CHECKING([for LHAPDF])
          AC_MSG_RESULT([(disabled)])
    fi


  else

    AC_PATH_PROG([LHAPDF_CONFIG], [lhapdf-config], [no], 
      [$wo_lhapdf_config_path])

    if test "$LHAPDF_CONFIG" != "no"; then
       LHAPDF_ROOT=`$LHAPDF_CONFIG --prefix`

       AC_CACHE_CHECK([the LHAPDF version],
       [wo_cv_lhapdf_version],
       [dnl,
	wo_cv_lhapdf_version=[`$LHAPDF_CONFIG --version | $SED -e 's/.*\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/'`]
    ])
       LHAPDF_FULL_VERSION="$wo_cv_lhapdf_version"
       AC_SUBST([LHAPDF_FULL_VERSION])

       AC_CACHE_CHECK([the major version],
       [wo_cv_lhapdf_major_version],
       [wo_cv_lhapdf_major_version=[`echo $wo_cv_lhapdf_version | $SED -e 's/\([0-9][0-9]*\)\..*/\1/'`]
       ])
       LHAPDF_MAJOR_VERSION="$wo_cv_lhapdf_major_version"
       AC_SUBST([LHAPDF_MAJOR_VERSION])

       AC_MSG_CHECKING([the LHAPDF pdfsets path])
       LHAPDF_PDFSETS_PATH=`$LHAPDF_CONFIG --pdfsets-path`
       if test "$LHAPDF_FULL_VERSION" = "5.5.0"; then
         LHAPDF_PDFSETS_PATH=`$LHAPDF_CONFIG --datarootdir`$LHAPDF_PDFSETS_PATH
       fi
       AC_MSG_RESULT([$LHAPDF_PDFSETS_PATH])

       AC_MSG_CHECKING([the standard PDF sets])
       if test -f "$LHAPDF_PDFSETS_PATH/cteq61.LHpdf" -a -f "$LHAPDF_PDFSETS_PATH/cteq5l.LHgrid" -a -f "$LHAPDF_PDFSETS_PATH/GSG961.LHgrid" -a -f "$LHAPDF_PDFSETS_PATH/cteq6ll.LHpdf"; then
          AC_MSG_RESULT([ all standard PDF sets installed])
       else
          AC_MSG_RESULT([ not all standard PDF sets installed])
          AC_MSG_NOTICE([error: *************************************************************])
          AC_MSG_NOTICE([error: LHAPDF standard PDF sets not installed, please install these ])
          AC_MSG_NOTICE([error:    PDF sets: cteq61.LHpdf, cteq6ll.LHpdf, cteq5l.LHgrid,     ])
          AC_MSG_NOTICE([error:	GSG961.LHgrid.     ])
          AC_MSG_NOTICE([error: *************************************************************])
          enable_lhapdf="no"
          AC_MSG_CHECKING([for LHAPDF])
          AC_MSG_RESULT([(disabled)])
       fi
     else
       enable_lhapdf="no"
     fi
  fi
   
else
  AC_MSG_CHECKING([for LHAPDF])
  AC_MSG_RESULT([(disabled)])
fi

AC_SUBST(LHAPDF_ROOT)
AC_SUBST(LHAPDF_FULL_VERSION)
AC_SUBST(LHAPDF_PDFSETS_PATH)

dnl LHAPDF requires the STD C++ library, when linking statically
if test "$enable_lhapdf" = "yes"; then
   wo_require_stdcpp="yes"
fi

if test "$enable_lhapdf" = "yes"; then
  if test "$LHAPDF_MAJOR_VERSION" = "5"; then
    wo_lhapdf_libdir="-L$LHAPDF_ROOT/lib"
    AC_LANG_PUSH([Fortran])
    AC_CHECK_LIB([LHAPDF],[getxminm],[LDFLAGS_LHAPDF="$wo_lhapdf_libdir -lLHAPDF"],
      [dnl
        AC_MSG_NOTICE([warning:  ********************************************************])
        AC_MSG_NOTICE([warning:  Either your LHAPDF version is too old (you need 5.3.0 or])
        AC_MSG_NOTICE([warning:  higher), or LHAPDF was compiled with a different FORTRAN])
        AC_MSG_NOTICE([warning:  compiler you and forgot to add the proper runtime to    ])
        AC_MSG_NOTICE([warning:  LIBS / LD_LIBRARY_PATH. Disabling LHAPDF support...     ])
        AC_MSG_NOTICE([warning:  ********************************************************])
        enable_lhapdf=no
      ],[$wo_lhapdf_libdir])
    AC_LANG_POP()
  elif test "$LHAPDF_MAJOR_VERSION" = "6"; then
    dnl now see if LHAPDF 6 is functional
    ACX_CHECK_LHAPDF()
  fi
fi
AC_SUBST(LDFLAGS_LHAPDF)


dnl Determine whether we need to stub photon-as-parton related bits of LHAPDF
dnl Only necessary for LHAPDF 5 as all versions of LHAPDF 6 support this
if test "$enable_lhapdf" = "yes"; then
  if test "$LHAPDF_MAJOR_VERSION" = "5"; then
    AC_LANG_PUSH([Fortran])
    AC_CHECK_LIB([LHAPDF],[has_photon],[test],[dnl
       AC_MSG_NOTICE([warning:  ********************************************************])
       AC_MSG_NOTICE([warning:  Your LHAPDF version is not supported for PDF sets like  ])
       AC_MSG_NOTICE([warning:  MRTS2004QED which include the photon as a parton ---    ])
       AC_MSG_NOTICE([warning:  don't try to use those!                                 ])
       AC_MSG_NOTICE([warning:  ********************************************************])
       LHAPDF5_HAS_PHOTON_DUMMY=true
     ],[$wo_lhapdf_libdir])
  fi
fi

AM_CONDITIONAL([LHAPDF6_AVAILABLE], [
   test "$enable_lhapdf" = "yes" -a "$LHAPDF_MAJOR_VERSION" = "6"])
AM_CONDITIONAL([LHAPDF5_AVAILABLE], [
   test "$enable_lhapdf" = "yes" -a "$LHAPDF_MAJOR_VERSION" = "5"])
if test "$enable_lhapdf" = "yes" -a "$LHAPDF_MAJOR_VERSION" = "6"; then
   LHAPDF6_AVAILABLE_FLAG=".true."
else
   LHAPDF6_AVAILABLE_FLAG=".false."
fi
if test "$enable_lhapdf" = "yes" -a "$LHAPDF_MAJOR_VERSION" = "5"; then
   LHAPDF5_AVAILABLE_FLAG=".true."
else
   LHAPDF5_AVAILABLE_FLAG=".false."
fi
AC_SUBST([LHAPDF5_AVAILABLE_FLAG])
AC_SUBST([LHAPDF6_AVAILABLE_FLAG])
AM_CONDITIONAL([LHAPDF5_HAS_PHOTON_DUMMY], [test -n "$LHAPDF5_HAS_PHOTON_DUMMY"])

])


AC_DEFUN([ACX_CHECK_LHAPDF],
[
save_CXXFLAGS="$CXXFLAGS"
save_LIBS="$LIBS"
CXXFLAGS="${CXXFLAGS} `${LHAPDF_CONFIG} --cxxflags`"
LIBS="${LIBS} `${LHAPDF_CONFIG} --ldflags`"
AC_MSG_CHECKING([if LHAPDF is functional])
AC_LANG_PUSH([C++])
AC_LINK_IFELSE(dnl
  [AC_LANG_PROGRAM([[#include "LHAPDF/LHAPDF.h"]], 
    [[using namespace LHAPDF; LHAPDF::PDF* pdf = LHAPDF::mkPDF("CT10nlo", 0); delete pdf;]])], 
  [lhapdfok="yes"], [lhapdfok="no"])
AC_MSG_RESULT([$lhapdfok])
CXXFLAGS="$save_CXXFLAGS"
LIBS="$save_LIBS"
AC_LANG_POP()
AC_MSG_CHECKING(LHAPDF)
if test "${lhapdfok}" = "yes"; then
   LHAPDF_CXXFLAGS="`${LHAPDF_CONFIG} --cxxflags`"
   LDFLAGS_LHAPDF="`${LHAPDF_CONFIG} --ldflags`"
   LHAPDF_LIBS="$LDFLAGS_LHAPDF"
   AC_SUBST(LHAPDF_LIBS)
   AC_MSG_RESULT(yes)
   $1          
else
   AC_MSG_RESULT(no)
   enable_lhapdf=no
   $2
fi

AC_SUBST([LHAPDF_CXXFLAGS])
AC_SUBST([LDFLAGS_LHAPDF])

])
