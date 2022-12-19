#! /bin/sh -x
########################################################################

# Edited by tests/Makefile using $(SED)
cascade_tests="%%cascade_tests%%"
top_builddir="%%top_builddir%%"/omega
srcdir="%%srcdir%%"
OCAML_NATIVE_EXT="%%OCAML_NATIVE_EXT%%"
SED="%%SED%%"

########################################################################

# Set up variables
omega_SM="$top_builddir/bin/omega_SM$OCAML_NATIVE_EXT"


########################################################################

# Run the tests:
for name in $cascade_tests; do
  file="$srcdir/$name"
  process="`$SED -n 1p $file`"
  cascade="`$SED -n 2p $file`"
  $SED -n '3,$p' $file >$name.expected
  $omega_SM -scatter "$process" -cascade "$cascade" \
	    -summary -forest -quiet 2>$name.result
  diff $name.expected $name.result
  rc=$?
  if test "$rc" -ne 0; then
    exit $rc
  else
    rm -f $name.expected $name.result
  fi
done

exit 0
