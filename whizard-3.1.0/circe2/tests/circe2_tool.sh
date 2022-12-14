#! /bin/sh
########################################################################

src_builddir=../../circe2/src

if test -f OCAML_FLAG; then
    if test -x $src_builddir/circe2_tool.opt; then
	circe2_tool=$src_builddir/circe2_tool.opt
    elif test -x $src_builddir/circe2_tool.bin; then
	circe2_tool=$src_builddir/circe2_tool.bin
    elif test -x $src_builddir/circe2_tool; then
	# ???
	circe2_tool=$src_builddir/circe2_tool
    else
	exit 2
    fi
else
    echo "|=============================================================================|"
    echo "No OCaml for testing tools available, test skipped"
    exit 77
fi


exec $circe2_tool -test
