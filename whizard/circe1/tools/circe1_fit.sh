#! /bin/sh
# mode=${2-slow}
mode=${2-fast}
root=`pwd`
indir=${root}/${3-input}
tmpdir=${root}/tmp
outdir=${root}/output
acc="${1-sband350 sband500 sband800 sband1000 sband1600
         tesla350 tesla500 tesla800 tesla1000 tesla1600
         tesla350-low tesla500-low tesla800-low tesla1000-low tesla1600-low
         xband350 xband500 xband800 xband1000 xband1600}"
xmkdir () {
  for d in "$@"; do
    mkdir $d 2>/dev/null || true
  done
}
rm -fr ${tmpdir}
xmkdir ${outdir} ${tmpdir}
cd ${tmpdir}
cat /dev/null >${outdir}/Params.f90
for a in $acc; do
  case "$a" in
    *1600*) energy=TEV16;;
    *1000*) energy=TEV1;;
     *800*) energy=GEV800;;
     *500*) energy=GEV500;;
  *3[56]0*) energy=GEV350;;
     *170*) energy=GEV170;;
      *90*) energy=GEV090;;
         *) energy=GEV500;;
  esac
  cp ${indir}/${a}_${mode}/lumidiff-??.dat .
  ${root}/circe1_fit.bin
  rm -fr ${outdir}/${a}_${mode}
  mkdir ${outdir}/${a}_${mode}
  cp Slices.mp4 ${outdir}
  cp Errors.tex lumidiff-??x[0-9][0-9].??? ${outdir}/${a}_${mode}
  sed -e "s/@ENERGY@/$energy/g" \
      -e "s/@ACC@/`echo $a | tr a-z A-Z | tr -cd A-Z`/g" Parameters \
    >>${outdir}/Params.f90
done
cd ${root}
rm -fr ${tmpdir}
cat >${outdir}/Params.tex <<'END'
\begin{table}
  \begin{center}
    \renewcommand{\arraystretch}{1.3}
    \begin{tabular}{|c||c|c|c|c|}\hline
      & \texttt{SBAND} & \texttt{TESLA} & \texttt{TESLA'} & \texttt{XBAND}
      \\\hline\hline
END
line () {
  for a in $acc; do
    case $a in
      *350* | *800* | *1000* | *1600*)
          ;;
      *)  echo -n ' & '
          sed -n $1p ${outdir}/${a}_${mode}/Errors.tex
          ;;
    esac
  done
  echo '\\\hline'
}
(echo '$\mathcal{L}/\text{fb}^{-1}\upsilon^{-1}$'; line 1
 echo '$\int d_{e^\pm}$';                          line 2
 echo '$x_{e^\pm}^\alpha$';                        line 3
 echo '$(1-x_{e^\pm})^\alpha$';                    line 4
 echo '$\int d_\gamma$';                           line 5
 echo '$x_\gamma^\alpha$';                         line 6
 echo '$(1-x_\gamma)^\alpha$';                     line 7
) >>${outdir}/Params.tex
cat >>${outdir}/Params.tex <<'END'
    \end{tabular}
  \end{center}
  \caption{\label{tab:param}%
    Version 1, revision 1997 04 16 of the beam spectra at 500 GeV.
    The rows correspond to the luminosity per effective year, the
    integral over the continuum and the powers in the factorized Beta
    distributions~(\ref{eq:beta}).}
\end{table}
END
cat >>${outdir}/Params.tex <<'END'
\begin{table}
  \begin{center}
    \renewcommand{\arraystretch}{1.3}
    \begin{tabular}{|c||c|c|c|c|}\hline
      & \texttt{SBAND} & \texttt{TESLA} & \texttt{TESLA'} & \texttt{XBAND}
      \\\hline\hline
END
line () {
  for a in $acc; do
    case $a in
      *1000*)
        echo -n ' & '
        sed -n $1p ${outdir}/${a}_${mode}/Errors.tex
        ;;
    esac
  done
  echo '\\\hline'
}
(echo '$\mathcal{L}/\text{fb}^{-1}\upsilon^{-1}$'; line 1
 echo '$\int d_{e^\pm}$';                          line 2
 echo '$x_{e^\pm}^\alpha$';                        line 3
 echo '$(1-x_{e^\pm})^\alpha$';                    line 4
 echo '$\int d_\gamma$';                           line 5
 echo '$x_\gamma^\alpha$';                         line 6
 echo '$(1-x_\gamma)^\alpha$';                     line 7
) >>${outdir}/Params.tex
cat >>${outdir}/Params.tex <<'END'
    \end{tabular}
  \end{center}
  \caption{\label{tab:param/TeV}%
    Version 1, revision 1997 04 17 of the beam spectra at 1 TeV.}
\end{table}
END
cat >>${outdir}/Params.tex <<'END'
\begin{table}
  \begin{center}
    \renewcommand{\arraystretch}{1.3}
    \begin{tabular}{|c||c|c|c|c|}\hline
      & 350 GeV & 500 GeV & 800 GeV & 1600 GeV
      \\\hline\hline
END
line () {
  for a in $acc; do
    case $a in
      tesla*-low)
        ;;
      tesla1000)
        ;;
      tesla*)
        echo -n ' & '
        sed -n $1p ${outdir}/${a}_${mode}/Errors.tex
        ;;
    esac
  done
  echo '\\\hline'
}
(echo '$\mathcal{L}/\text{fb}^{-1}\upsilon^{-1}$'; line 1
 echo '$\int d_{e^\pm}$';                          line 2
 echo '$x_{e^\pm}^\alpha$';                        line 3
 echo '$(1-x_{e^\pm})^\alpha$';                    line 4
 echo '$\int d_\gamma$';                           line 5
 echo '$x_\gamma^\alpha$';                         line 6
 echo '$(1-x_\gamma)^\alpha$';                     line 7
) >>${outdir}/Params.tex
cat >>${outdir}/Params.tex <<'END'
    \end{tabular}
  \end{center}
  \caption{\label{tab:param/Tesla}%
    Version 1, revision 1997 04 17 of the beam spectra for TESLA.}
\end{table}
END
exit 0
