if test x"${srcdir}" = x ; then srcdir=. ; fi
	./rsbench -oa -Ob --bench -f ${srcdir}/A.mtx -qH -R -n1,4 -T z --verbose --nrhs 1,2 --by-rows || exit 255
	./rsbench -oa -Ob --help || exit 255
	./rsbench --help || exit 255
	./rsbench --version || exit 255
	./rsbench -I || exit 255
	./rsbench -C || exit 255
	./rsbench -oa -Ob --bench -f ${srcdir}/A.mtx --verbose --nrhs 1,4 --by-rows || exit 255
