dnl aux.m4 -- Auxiliary macros for WHIZARD's configure.ac
dnl

dnl Horizontal line for readability:
AC_DEFUN([WO_HLINE],
[AC_MSG_NOTICE([--------------------------------------------------------------])])

dnl Message at the beginning of a configure section
AC_DEFUN([WO_CONFIGURE_SECTION],
[WO_HLINE()
AC_MSG_NOTICE([--- ]$1[ ---])
AC_MSG_NOTICE([])
])


dnl Define a variable and export it
dnl WO_SET(variable, value)
AC_DEFUN([WO_SET],
[$1=$2
AC_SUBST($1)])

dnl Add a list of names to the list of output or executable files
dnl WO_OUTFILES(subdir, files)
AC_DEFUN([WO_OUTFILES],
[for file in $2; do
  OUTFILES="$OUTFILES $1/$file";
done])

dnl WO_EXECUTABLES(subdir, files)
AC_DEFUN([WO_EXECUTABLES],
[for file in $2; do
  OUTFILES="$OUTFILES $1/$file";
  EXECUTABLES="$EXECUTABLES $1/$file"
done])

dnl Add a list of names to the list of automake-controlled files
AC_DEFUN([WO_AUTOMAKE],
[for file in $2; do
  AUTOMAKEFILES="$AUTOMAKEFILES $1/$file";
done])

dnl This adds a package which resides in a subdirectory of 'src'
dnl The `variable' is inserted into AC_SUBST; it will refer to
enl the package subdirectory in Makefiles.  The package
dnl resides in src/XXX and  and optionally has a Makefile (or similar).
dnl WO_PACKAGE(variable, identifier [,Makefiles [,Executables]])
AC_DEFUN([WO_PACKAGE],
[WO_SET($1,src/$2)
ifelse($#, 3,
[WO_OUTFILES($$1, $3)
WO_AUTOMAKE($$1, $3)])
ifelse($#, 4,
[WO_OUTFILES($$1, $3)
WO_EXECUTABLES($$1, $4)])
])

dnl This is like WO_PACKAGE, but it calls AC_ARG_ENABLE in addition.
dnl If the package `id' is disabled or the 'file' is not found,
dnl the variable `var' is set to "no".
dnl WO_PACKAGE_ENABLE(var, id, help, file [,Makefiles [,Executables]])
AC_DEFUN([WO_PACKAGE_ENABLE],
[AC_ARG_ENABLE($2,[$3])
if test "$enable_$2" = "no"; then
AC_MSG_CHECKING([for src/$2/$4])
WO_SET($1, no)
AC_MSG_RESULT([(disabled)])
else
AC_CHECK_FILE(src/$2/$4, enable_$2=yes, enable_$2=no)
if test "$enable_$2" = "no"; then
WO_SET($1, no)
else
ifelse($#, 4, [WO_PACKAGE($1, $2)])
ifelse($#, 5, [WO_PACKAGE($1, $2, $5)])
ifelse($#, 6, [WO_PACKAGE($1, $2, $5, $6)])
fi
fi
])

dnl The same, but disabled by default.
dnl If the package `id' is disabled or the 'file' is not found,
dnl the variable `var' is set to "no".
dnl WO_PACKAGE_DISABLE(var, id, help, file [,Makefiles [,Executables]])
AC_DEFUN([WO_PACKAGE_DISABLE],
[AC_ARG_ENABLE($2,[$3])
if test "$enable_$2" = "yes"; then
AC_CHECK_FILE(src/$2/$4, enable_$2=yes, enable_$2=no)
if test "$enable_$2" = "no"; then
WO_SET($1, no)
else
ifelse($#, 4, [WO_PACKAGE($1, $2)])
ifelse($#, 5, [WO_PACKAGE($1, $2, $5)])
ifelse($#, 6, [WO_PACKAGE($1, $2, $5, $6)])
fi
else
enable_$2="no"
AC_MSG_CHECKING([for src/$2/$4])
WO_SET($1, no)
AC_MSG_RESULT([(disabled)])
fi
])


dnl Extension of AC_PATH_PROG: Search for libraries for which the name
dnl is not exactly known (because it may have the version number in it)
dnl Set output variables $var, $var_DIR, $var_LIB accordingly
dnl WO_PATH_LIB(var, id, name-list, path)
AC_DEFUN([WO_PATH_LIB],
[AC_CACHE_CHECK(for $3, wo_cv_path_$1,
[case "$$1" in
  /*)
  wo_cv_path_$1="$$1" # User-supplied path
  ;;
  *)
  IFS="${IFS=   }"; ac_save_ifs="$IFS"; IFS=":"
  ac_dummy=$4
  for ac_dir in $ac_dummy; do
    test -z "$ac_dir" && ac_dir=.
    unset wo_cv_path_$1
    ac_pwd=`pwd`
    if test -d "$ac_dir"; then
      cd $ac_dir
      for ac_word in $3; do
        test -f "$ac_word" && wo_cv_path_$1="$ac_dir/$ac_word"
      done
      cd $ac_pwd
    fi
    if test -n "$wo_cv_path_$1"; then
      break
    fi
  done
  IFS="$ac_save_ifs"
  if test -z "$wo_cv_path_$1"; then
    wo_cv_path_$1="no"
  fi
  ;;
esac
])
$1=$wo_cv_path_$1
if test "$$1" != "no"; then
$1_DIR=`echo $$1 | sed -e 's|/lib$2.*\.a$||'`
$1_LIB=`echo $$1 | sed -e 's|^.*/lib\($2.*\)\.a$|\1|'`
fi
AC_SUBST($1)
AC_SUBST($1_DIR)
AC_SUBST($1_LIB)
])

dnl WHIZARD summary
AC_DEFUN([WO_SUMMARY],[dnl
WO_VERSION=AC_PACKAGE_VERSION()
function echo_summary () {
echo "|=============================================================================|"
echo "|                                                                             |"
echo "|    WW             WW  WW   WW  WW  WWWWWW      WW      WWWWW    WWWW        |"
echo "|     WW    WW     WW   WW   WW  WW     WW      WWWW     WW  WW   WW  WW      |"
echo "|      WW  WW WW  WW    WWWWWWW  WW    WW      WW  WW    WWWWW    WW   WW     |"
echo "|       WWWW   WWWW     WW   WW  WW   WW      WWWWWWWW   WW  WW   WW  WW      |"
echo "|        WW     WW      WW   WW  WW  WWWWWW  WW      WW  WW   WW  WWWW        |"
echo "|                                                                             |"
echo "|                                                                             |"
echo "|                                        W                                    |"
echo "|                                       sW                                    |"
echo "|                                       WW                                    |"
echo "|                                      sWW                                    |"
echo "|                                      WWW                                    |"
echo "|                                     wWWW                                    |"
echo "|                                    wWWWW                                    |"
echo "|                                    WW WW                                    |"
echo "|                                    WW WW                                    |"
echo "|                                   wWW WW                                    |"
echo "|                                  wWW  WW                                    |"
echo "|                                  WW   WW                                    |"
echo "|                                  WW   WW                                    |"
echo "|                                 WW    WW                                    |"
echo "|                                 WW    WW                                    |"
echo "|                                WW     WW                                    |"
echo "|                                WW     WW                                    |"
echo "|           wwwwww              WW      WW                                    |"
echo "|              WWWWWww          WW      WW                                    |"
echo "|                 WWWWWwwwww   WW       WW                                    |"
echo "|                     wWWWwwwwwWW       WW                                    |"
echo "|                 wWWWWWWWWWWwWWW       WW                                    |"
echo "|                wWWWWW       wW        WWWWWWW                               |"
echo "|                  WWWW       wW        WW  wWWWWWWWwww                       |"
echo "|                   WWWW                      wWWWWWWWwwww                    |"
echo "|                     WWWW                      WWWW     WWw                  |"
echo "|                       WWWWww                   WWWW                         |"
echo "|                           WWWwwww              WWWW                         |"
echo "|                               wWWWWwww       wWWWWW                         |"
echo "|                                     WwwwwwwwwWWW                            |"
echo "|                                                                             |"
echo "|                                                                             |"
echo "|                                                                             |"
echo "|  by:   Wolfgang Kilian, Thorsten Ohl, Juergen Reuter                        |"
echo "|        with contributions from Christian Speckner                           |"
echo "|        Contact: <whizard@desy.de>                                           |"
echo "|                                                                             |"
echo "|  if you use WHIZARD please cite:                                            |"
echo "|        W. Kilian, T. Ohl, J. Reuter,  Eur.Phys.J.C71 (2011) 1742            |"
echo "|                                          @<:@arXiv: 0708.4233 @<:@hep-ph@:>@@:>@        |"
echo "|        M. Moretti, T. Ohl, J. Reuter, arXiv: hep-ph/0102195                 |"
echo "|                                                                             |"
echo "|=============================================================================|"
echo "**************************************************************"
echo "--------------------------------------------------------------"
echo "---      AC_PACKAGE_NAME() CONFIGURATION SUMMARY      ---"
echo "**************************************************************"
echo "Package name: AC_PACKAGE_NAME()"
echo "Version:      AC_PACKAGE_VERSION()"
echo "Date:         $PACKAGE_DATE"
echo "Status:       $PACKAGE_STATUS"
echo "**************************************************************"
echo "---      Compilers      ---"
echo "--------------------------------------------------------------"
echo "Fortran compiler: --- $FC_VENDOR ---"
echo "         Version: --- $FC_VERSION ---"
echo "           Flags: --- $FCFLAGS ---"
echo " float precision: --- $FC_PRECISION ---"
if test "$FC_DEBUG_ON" = ".true."; then
   echo "  debug features: --- on ---"
else
   echo "  debug features: --- off ---"
fi
if test "$FC_OPENMP_OFF" = "!" ; then
   echo "          OpenMP: --- on with max. $FC_OPENMP_DEFAULT_MAX_THREADS threads"
elif test "$FC_OPENMP_ON" = "!" ; then
     echo "          OpenMP: --- off ---"
fi
if test "$MPI_AVAILABLE" = "yes" ; then
     echo "             MPI: --- on ---"
     echo "     MPI Library: --- $MPI_LIBRARY, v$MPI_VERSION ---"
else
     echo "             MPI: --- off ---"
fi
echo "--------------------------------------------------------------"
echo "  OCaml compiler: --- $OCAMLOPT ---"
echo "         Version: --- $OCAMLVERSION ---"
echo "           Flags: --- $OCAMLFLAGS ---"
echo "--------------------------------------------------------------"
echo "    C++ compiler: --- $CXX ---      @<:@interfaces only@:>@"
echo "           Flags: --- $CXXFLAGS ---"
echo "--------------------------------------------------------------"
echo " Python compiler: --- $PYTHON ---   @<:@interfaces only@:>@"
echo "         Version: --- $PYTHON_FULL_VERSION ---"
if test "$PYTHON_API" != "no" ; then
   echo " WhiPy interface: --- yes ---"
else
   echo " WhiPy interface: --- no ---"
fi
echo "**************************************************************"
echo "---      Internal and shipped packages      ---"
echo "--------------------------------------------------------------"
echo "VAMP  (multi-channel adapative integrator) :   yes, v$WO_VERSION"
if test "$OCAMLOPT" != "no" ; then
   echo "O'Mega (matrix element generator)          :   yes, v$WO_VERSION"
else
   echo "O'Mega (matrix element generator)          :   no"
fi
echo "CIRCE1 (lepton beam spectra, parameterized):   yes, v$WO_VERSION"
echo "CIRCE2 (lepton beam spectra, sampled)      :   yes, v$WO_VERSION"
if test "$OCAMLOPT" != "no" ; then
   echo "    incl. tools for generating new spectra :   yes"
else
   echo "    incl. tools for generating new spectra :   no"
fi
echo "--------------------------------------------------------------"
if test "$PYTHIA6_AVAILABLE_FLAG" = ".true." ; then
   echo "PYTHIA6 (parton showering & hadronization) :   yes, v6.427"
   if test "$PYTHIA6_EH_AVAILABLE_FLAG" = "yes" ; then
      echo "        (settings for eh collisions)       :   yes"
   else
      echo "        (settings for eh collisions)       :   no"
   fi
   echo "TAUOLA (tau decays)                        :   yes"
else
   echo "PYTHIA6 (parton showering & hadronization) :   no"
   echo "TAUOLA (tau decays)                        :   no"
fi
echo "StdHEP (event format)                      :   yes, v5.06.01"
echo "--------------------------------------------------------------"
echo "---      External packages      ---"
echo "--------------------------------------------------------------"
if test "$HEPMC_AVAILABLE_FLAG" = "yes" ; then
   echo "HepMC (event format):   yes, v$HEPMC_VERSION"
else
   echo "HepMC (event format):   no"
fi
if test "$HDF5_AVAILABLE_FLAG" = "yes" ; then
   echo "HDF5 (binary format):   yes, v$HDF5_VERSION"
else
   echo "HDF5 (binary format):   no"
fi
if test "$LCIO_AVAILABLE_FLAG" = "yes" ; then
   echo "LCIO (event format) :   yes, v$LCIO_VERSION"
else
   echo "LCIO (event format) :   no"
fi
if test "$LHAPDF5_AVAILABLE_FLAG" = ".true." ; then
   echo "LHAPDF (PDF sets)   :   yes, v$LHAPDF_FULL_VERSION"
elif test "$LHAPDF6_AVAILABLE_FLAG" = ".true." ; then
   echo "LHAPDF (PDF sets)   :   yes, v$LHAPDF_FULL_VERSION"
   echo "        PDF set path:   $LHAPDF_PDFSETS_PATH"
else
   echo "LHAPDF (PDF sets)   :   no"
fi
if test "$HOPPET_AVAILABLE_FLAG" = ".true." ; then
   echo "HOPPET (PDF match.) :   yes, v$HOPPET_VERSION"
else
   echo "HOPPET (PDF match.) :   no"
fi
if test "$FASTJET_AVAILABLE_FLAG" = "yes" ; then
   echo "FastJet (clustering):   yes, v$FASTJET_VERSION"
else
   echo "FastJet (clustering):   no"
fi
if test "$PYTHIA8_AVAILABLE_FLAG" = ".true." ; then
   echo "PYTHIA8 (QCD)       :   yes, v$PYTHIA8_VERSION"
else
   echo "PYTHIA8 (QCD)       :   no"
fi
if test "$GOSAM_AVAILABLE_FLAG" = ".true." ; then
   echo "GoSam (OLP)         :   yes, v$GOSAM_VERSION"
else
   echo "GoSam (OLP)         :   no"
fi
if test "$OPENLOOPS_AVAILABLE_FLAG" = ".true." ; then
   echo "OpenLoops (OLP)     :   yes, v$OPENLOOPS_VERSION"
   echo "           path     :   $OPENLOOPS_DIR"
else
   echo "OpenLoops (OLP)     :   no"
fi
if test "$RECOLA_AVAILABLE_FLAG" = ".true." ; then
   echo "RECOLA (OLP)        :   yes, v$RECOLA_VERSION"
   echo "           path     :   $RECOLA_DIR"   
else
   echo "RECOLA (OLP)        :   no"
fi
if test "$LOOPTOOLS_AVAILABLE_FLAG" = ".true." ; then
   echo "LoopTools           :   yes"
else
   echo "LoopTools           :   no"
fi
echo "--------------------------------------------------------------"
if test "$SIP_ACTIVE" = "yes" ; then
echo "**************************************************************"
echo "***              MAC OS X Darwin system with               ***"
echo "***    Security Integrity Protection (SIP) enabled.        ***"
echo "***  'make check' will not work, and most likely also      ***"
echo "***  'make installcheck' will not work. The installed      ***"
echo "***            WHIZARD will work as intended.              ***"
echo "**************************************************************"
fi
}
echo_summary | tee config-summary.log
])
