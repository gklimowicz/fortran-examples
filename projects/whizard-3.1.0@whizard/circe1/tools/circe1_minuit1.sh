#! /bin/sh
minuit_bin=`pwd`/circe1_minuit1.bin
tmp="$IFS"
  IFS=:
  args=":$*:"
IFS="$tmp"
filter="| \
  awk '/STATUS=(CONVERGED|CALL LIMIT|FAILED)/ { p=1; print }; \
       /@.* \.00000 *fixed/ { next }; \
       /EDM=|CHI2|@/ && p { print }' "
case "$args" in
  *:v:*) filter=;;
esac
(
  cat <<END
  set title
  CIRCE
  parameters
  1  '@ 1        ' 0.00 0.01
  2  '@ lx       ' 0.20 0.01
  3  '@ l(1-x)   ' 0.20 0.01
  4  '@ llx      ' 0.00 0.01
  5  '@ ll(1-x)  ' 0.00 0.01
  6  '@ x        ' 0.00 0.01
  7  '@ lx^2     ' 0.00 0.01
  8  '@ l(1-x)^2 ' 0.00 0.01
  9  '@ llx^2    ' 0.00 0.01
  10 '@ ll(1-x)^2' 0.00 0.01
  11 '@ x^2      ' 0.00 0.01
  12 '@ 1/lx     ' 0.00 0.01
  13 '@ 1/l(1-x) ' 0.00 0.01
  14 '@ 1/llx    ' 0.00 0.01
  15 '@ 1/ll(1-x)' 0.00 0.01
  16 '@ 1/x      ' 0.00 0.01
  17 '@ 1/(1-x)  ' 0.00 0.01

  END
  for p in 1  2  3  4  5  6  7  8  9 10 \
          11 12 13 14 15 16 17; do
    case "$args" in
      *:$p=*:*) val=`echo "$args" | sed 's/.*:'"$p"'=\\([0-9.-]*\\):.*/\\1/'`;
                echo set parameter $p $val;
                echo fix $p;;
        *:$p:*) ;;
             *) echo fix $p;;
    esac
  done
  case "$args" in
    *:S0:*) echo set strategy 0;;
    *:S1:*) echo set strategy 1;;
    *:S2:*) echo set strategy 2;;
  esac
  cat <<END
  migrat 10000 0.01
  stop
  END
) | eval "$minuit_bin $filter"
case "$args" in
  *:p:*) awk '$5 != "" { print $1, $2, $5 }' minuit.fit > chi2
         awk '$5 != "" { print $1, $5 }' minuit.fit > chix
         awk '$5 != "" { print $2, $5 }' minuit.fit > chiy
         gnuplot -geometry -0+0 plot2 >/dev/null 2>&1
esac
exit 0
