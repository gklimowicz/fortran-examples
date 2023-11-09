dnl latex.m4 -- checks for LaTeX programs
dnl

AC_DEFUN([AC_PROG_TEX], [dnl
AC_CHECK_PROGS(PLAINTEX,[tex],no)
AM_CONDITIONAL([TEX_AVAILABLE], [test "$PLAINTEX" != "no"])
AC_SUBST(PLAINTEX)
])

AC_DEFUN([AC_PROG_LATEX], [dnl
AC_CHECK_PROGS(LATEX,[latex elatex lambda],no)
AM_CONDITIONAL([LATEX_AVAILABLE], [test "$LATEX" != "no"])
AC_SUBST(LATEX)
])


AC_DEFUN([_AC_LATEX_TEST], [dnl
AC_REQUIRE([AC_PROG_LATEX])
rm -rf .tmps_latex
mkdir .tmps_latex
cd .tmps_latex
ifelse($#,2,[
$2="no"; export $2;
cat &gt; testconf.tex &lt;&lt; \EOF
$1
EOF
],$#,3,[
echo "\\documentclass{$3}" &gt; testconf.tex
cat &gt;&gt; testconf.tex &lt;&lt; \EOF
$1
EOF
],$#,4,[
echo "\\documentclass{$3}" &gt; testconf.tex
echo "\\usepackage{$4}" &gt; testconf.tex
cat &gt;&gt; testconf.tex &lt;&lt; \EOF
$1
])
cat testconf.tex | $latex 2&gt;&amp;1 1&gt;/dev/null &amp;&amp; $2=yes; export $2;
cd ..
rm -rf .tmps_latex
])

dnl AC_LATEX_CLASSES([book],book)
dnl should set $book="yes"
dnl 
dnl AC_LATEX_CLASSES(allo,book)
dnl should set $book="no"


AC_DEFUN([AC_LATEX_CLASS], [dnl
AC_CACHE_CHECK([for class $1],[ac_cv_latex_class_]translit($1,[-],[_]),[
_AC_LATEX_TEST([
\begin{document}
\end{document}
],[ac_cv_latex_class_]translit($1,[-],[_]),$1)
])
$2=$[ac_cv_latex_class_]translit($1,[-],[_]) ; export $2;
AC_SUBST($2)
])

dnl Checking for dvips

AC_DEFUN([AC_PROG_DVIPS], [dnl
AC_CHECK_PROGS(DVIPS,dvips,no)
AM_CONDITIONAL([DVIPS_AVAILABLE], [test "$DVIPS" != "no"])
AC_SUBST(DVIPS)
])


dnl Checking for ps2pdf and friends

AC_DEFUN([AC_PROG_PS2PDF], [dnl
AC_CHECK_PROGS(PS2PDF,[ps2pdf14 ps2pdf13 ps2pdf12 ps2pdf],no)
AM_CONDITIONAL([PS2PDF_AVAILABLE], [test "$PS2PDF" != "no"])
AC_SUBST(PS2PDF)
if test "$enable_distribution" = "yes"; then
if test "$PDFLATEX" = "no"; then
AC_MSG_NOTICE([error: **********************************])
AC_MSG_NOTICE([error: No way to make documentation PDFs.])
AC_MSG_ERROR([**********************************])
fi
fi
])

dnl Checking for epspdf and epstopdf
dnl epspdf older than 0.4.3 do have problems on the MAC

AC_DEFUN([AC_PROG_EPSPDF], [dnl
AC_CHECK_PROGS(EPSPDF,[epspdf],no)
AM_CONDITIONAL([EPSPDF_AVAILABLE], [test "$EPSPDF" != "no"])
AC_SUBST(EPSPDF)
if test "$EPSPDF" != "no"; then
  EPSPDFVERSION=`$EPSPDF --help 2>&1 | $GREP -i "epspdf.*[[0-9]]\+\.[[0-9]]\+\.[[0-9\]]\+" | tail -n 1 | $SED "s/[[^0123456789\.]]//g"`
  AC_MSG_RESULT([epspdf version is $EPSPDFVERSION])	
  AC_CACHE_VAL([wo_epspdf_cv_integer_version],
    [wo_epspdf_cv_integer_version="`echo $EPSPDFVERSION | \
      $AWK '[$]1 {
        changequote(<<,>>)dnl
          split (<<$>>1, version, "[.+]+");
          printf ("%d%02d%03d", version[1], version[2], version[3])}'`"
        changequote([,])])
  EPSPDFINTEGERVERSION=$wo_epspdf_cv_integer_version
  AC_SUBST([EPSPDFVERSION])
  if test $EPSPDFINTEGERVERSION -ge 004003; then
    AC_MSG_RESULT([epspdf present and newer than 0.4.2; using epspdf for conversion.])
  else
    AC_MSG_RESULT([epspdf version older than 0.4.3; searching for epstopdf.])
  fi
  AM_CONDITIONAL([EPSPDF_043],
	[test $EPSPDFINTEGERVERSION -ge 004003])
else  
  AM_CONDITIONAL([EPSPDF_043],false)  
  EPSPDFINTEGERVERSION=000000
fi
  AC_SUBST([EPSPDFINTEGERVERSION])

])

dnl Catching the buggy version 2.9.8

AC_DEFUN([AC_PROG_EPSTOPDF], [dnl
AC_CHECK_PROGS(EPSTOPDF,[epstopdf],no)
if test "$EPSTOPDF" != "no"; then
  wo_epstopdf_cv_version="`$EPSTOPDF -v 2>&1 | \
     $AWK '{
       if (NR == 1 && toupper([$]1) ~ /EPSTOPDF/) { printf [$]9 } else 
      {if (NR == 2 && toupper([$]1) ~ /EPSTOPDF/) { printf substr([$]2,1,length([$]2)-1) }}}'`" 
  EPSTOPDFVERSION=$wo_epstopdf_cv_version 
  AC_MSG_RESULT([epstopdf version is $EPSTOPDFVERSION])	
  if test "$EPSTOPDFVERSION" = "2.9.8" -o "$EPSTOPDFVERSION" = "2.9.8gw" \
            -o "$EPSTOPDFVERSION" = "2.9.8GW"; then
    AC_MSG_RESULT([********************************************************])	
    AC_MSG_RESULT([epstopdf version 2.9.8 is known to be buggy, be careful.])	
    AC_MSG_RESULT([********************************************************])	
    AM_CONDITIONAL([EPSTOPDF_AVAILABLE], true)
    EPSTOPDF_BUGGY="yes"
  else
    AM_CONDITIONAL([EPSTOPDF_AVAILABLE], true)
    EPSTOPDF_BUGGY="no"
  fi
else
    AM_CONDITIONAL([EPSTOPDF_AVAILABLE], false)
    EPSTOPDF_BUGGY="no"
fi
AC_SUBST(EPSTOPDF)
AC_SUBST(EPSTOPDF_BUGGY)
AC_SUBST(EPSTOPDFVERSION)
if test "$enable_distribution" = "yes"; then
if test "$EPSPDF" = "no" -a "$EPSTOPDF" = "no"; then
AC_MSG_NOTICE([error: ***********************************************])
AC_MSG_NOTICE([error: Neither a viable epspdf nor epstopdf available.])
AC_MSG_ERROR([***********************************************])
fi
fi
])

dnl Checking for supp-pdf.tex (auxiliary for PDF output)

AC_DEFUN([AC_PROG_SUPP_PDF], [dnl
AC_REQUIRE([AC_PROG_TEX])
AC_CACHE_CHECK([for supp-pdf.tex],
[wo_cv_supp_pdf_exists],
[dnl
wo_cv_supp_pdf_exists="no"
if test "$PLAINTEX" != "no"; then
  wo_cmd='echo \\input supp-pdf.tex \\end > conftest.tex'
  eval "$wo_cmd"
  wo_cmd='$PLAINTEX conftest.tex >&5'
  (eval "$wo_cmd") 2>&5 && wo_cv_supp_pdf_exists="yes" 
fi])
AM_CONDITIONAL([SUPP_PDF_AVAILABLE], [test "$wo_cv_supp_pdf_exists" != "no"])
if test "$enable_distribution" = "yes"; then
if test "$wo_cv_supp_pdf_exists" = "no"; then
AC_MSG_NOTICE([error: ********************************************************])
AC_MSG_NOTICE([error: No LaTeX supp-pdf.tex available, please install conTeXt.])
AC_MSG_ERROR([********************************************************])
fi
fi
])

dnl Checking for pdflatex

AC_DEFUN([AC_PROG_PDFLATEX], [dnl
AC_CHECK_PROGS(PDFLATEX,[pdflatex],no)
AM_CONDITIONAL([PDFLATEX_AVAILABLE], [test "$PDFLATEX" != "no"])
AC_SUBST(PDFLATEX)
])

dnl Checking for makeindex

AC_DEFUN([AC_PROG_MAKEINDEX], [dnl
AC_CHECK_PROGS(MAKEINDEX,[makeindex],no)
AM_CONDITIONAL([MAKEINDEX_AVAILABLE], [test "$MAKEINDEX" != "no"])
AC_SUBST(MAKEINDEX)
])

dnl Checking for Metapost

AC_DEFUN([AC_PROG_MPOST], 
[dnl
  AC_CHECK_PROGS([MPOST],[mpost metapost],[no])
  if test "$MPOST" != "no"; then        
        MPOSTVERSIONSTR=`$MPOST -version | $GREP -i 'MetaPost' | head -1` 
	MPOSTVERSION=[`echo $MPOSTVERSIONSTR | $SED -e 's/.*MetaPost \([0-9]*\.[0-9][0-9][0-9]*\).*/\1/'`]
        AC_CACHE_VAL([wo_mpost_cv_integer_version],
          [wo_mpost_cv_integer_version="`echo "$MPOSTVERSION" | \
            $AWK 'NR==1 {
              changequote(<<,>>)dnl
                split (<<$>>1, version, "[.+]+");
                printf ("%d%03d", version[1], version[2])}'`"
              changequote([,])])
        MPOSTINTEGERVERSION=$wo_mpost_cv_integer_version
        AC_MSG_RESULT([MetaPost version is $MPOSTVERSION])
	if test $MPOSTINTEGERVERSION -ge 1750; then
           MPOSTFLAG="--math=scaled"
	else
	   MPOSTFLAG=""
        fi
	MPOST_AVAILABLE_FLAG=".true."
  else
        MPOST_AVAILABLE_FLAG=".false."
  fi
AC_SUBST([MPOSTVERSION])
AC_SUBST([MPOSTINTEGERVERSION])
AM_CONDITIONAL([MPOST_AVAILABLE], [test "$MPOST" != "no"])
AC_SUBST(MPOST_AVAILABLE_FLAG)
AC_SUBST(MPOST)
AC_SUBST(MPOSTFLAG)
])

dnl Putting the above together, check possibilities for online event analysis
AC_DEFUN([WO_CHECK_EVENT_ANALYSIS_METHODS], [dnl
AC_REQUIRE([AC_PROG_LATEX])
AC_REQUIRE([AC_PROG_MPOST])
AC_REQUIRE([AC_PROG_DVIPS])
AC_REQUIRE([AC_PROG_PS2PDF])

AC_CACHE_CHECK([whether we can display event analysis],
[wo_cv_event_analysis],
[dnl
if test "$LATEX" != "no" -a "$MPOST" != "no"; then
  wo_cv_event_analysis="yes"
else
  wo_cv_event_analysis="no"
fi
])
EVENT_ANALYSIS="$wo_cv_event_analysis"
EVENT_ANALYSIS_PS="$EVENT_ANALYSIS"
EVENT_ANALYSIS_PDF="$EVENT_ANALYSIS"

if test "$EVENT_ANALYSIS" != "no"; then
AC_CACHE_CHECK([whether we can display event analysis in PostScript format],
[wo_cv_event_analysis_ps],
[dnl
if test "$DVIPS" != "no"; then
  wo_cv_event_analysis_ps="yes"
else
  wo_cv_event_analysis_ps="no"
fi])
EVENT_ANALYSIS_PS="$wo_cv_event_analysis_ps"
fi

if test "$EVENT_ANALYSIS_PS" != "no"; then
AC_CACHE_CHECK([whether we can display event analysis in PDF format],
[wo_cv_event_analysis_pdf],
[dnl
if test "$PS2PDF" != "no"; then
  wo_cv_event_analysis_pdf="yes"
else
  wo_cv_event_analysis_pdf="no"
fi])
EVENT_ANALYSIS_PDF="$wo_cv_event_analysis_pdf"
fi

AM_CONDITIONAL([EVENT_ANALYSIS_AVAILABLE], [test "$EVENT_ANALYSIS" = "yes" -a "$EVENT_ANALYSIS_PS" = "yes" -a "$EVENT_ANALYSIS_PDF" = "yes"])

AC_SUBST([EVENT_ANALYSIS])
AC_SUBST([EVENT_ANALYSIS_PS])
AC_SUBST([EVENT_ANALYSIS_PDF])
])


dnl Checking for gzip (putting this together with the LaTeX part)

AC_DEFUN([AC_PROG_GZIP],[
AC_CHECK_PROGS(GZIP,[gzip],no)
AM_CONDITIONAL([GZIP_AVAILABLE], [test "$GZIP" != "no"])
if test "$enable_distribution" = "yes"; then
if test "$GZIP" = "no"; then
AC_MSG_NOTICE([error: *********************************************])
AC_MSG_NOTICE([error: Gzip not installed, no distribution possible.])
AC_MSG_ERROR([*********************************************])
fi
fi
AC_SUBST(GZIP)
])
