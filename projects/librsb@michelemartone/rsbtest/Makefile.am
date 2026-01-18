ACLOCAL_AMFLAGS=-I m4
AM_CXXFLAGS=-pipe ${OPENMP_CXXFLAGS}
if WANT_PROGRAMS
bin_PROGRAMS=rsbtest$(EXEEXT)
else
noinst_PROGRAMS=rsbtest$(EXEEXT)
endif
rsbtest_SOURCES=rsbtest.cpp
rsbtest_LDADD=@LIBRSB_LOCAL_LIBS@ $(LIBS)
if WITH_RSBTEST
rsbtest_DEPENDENCIES=@LIBRSB_LOCAL_LIBS@
endif
.PHONY: tests
tests: rsbtest$(EXEEXT)
if HAVE_RUNNABLE_TESTS
	./$< --max-test-time 1800 --quiet
endif

if HAVE_DEV_NULL
RSBT_DEV_NULL=/dev/null
else
RSBT_DEV_NULL=_dev_null.tmp
endif

if HAVE_RSBPP_HPP
RSBT_PP_HPP_PC=
else
RSBT_PP_HPP_PC=!
endif

if HAVE_RSB_HPP
RSBT_HPP_PC=
else
RSBT_HPP_PC=!
endif

RTP=--no-tune --max_t 0.01
STP=--no-tune --max_t 0.01 --max 1 --nrhs 1 --beta 1 --incy 1 --incx 1 --no-trans --alpha 1 --type d --rand # short exec time
BTP=--skip-loading-hermitian-matrices --skip-loading-unsymmetric-matrices --tune-maxt 10 --tune-maxr 10 --verbose-tuning --extra-verbose-interface --min_t 0.01 --max_t 0.01 --mintimes 1 --maxtimes 1 --verbose # bloat, but allow symmetric loading

if HAVE_DD
if HAVE_FULL_DD
IBC= dd if=/dev/stdin of=/dev/stdout bs=16M status=none iflag=fullblock # input buffer command
else
IBC= dd if=/dev/stdin of=/dev/stdout  # input buffer command
endif
endif

.PHONY: qtests
qtests: rsbtest$(EXEEXT)
	if test ! -f A.mtx ; then cp $(abs_top_srcdir)/A.mtx . ; fi
if HAVE_RUNNABLE_TESTS
if HAVE_DD
	./$< --version | $(IBC) | grep -q -i using # diagnostic
	./$< --version | $(IBC) | grep -q librsb
	./$< --help | $(IBC) | grep -q help
	./$<  -V       | $(IBC) | grep -q librsb
	./$< --types all --types abcd --types '?' --verbose | $(IBC) | grep -q configured.to.support.types
	./$< --quiet --types all --only-test-case-n  4        | $(IBC) | grep -q '\<1.*success'
	./$< $(RTP) --serial | $(IBC) | grep -q Building
if WITH_CXX17_FILESYSTEM
	./$< $(STP) --serial . | $(IBC) | grep -q adding
endif
if USE_MKL_SPBLAS
	./$< --mkl A.mkl
else
	! ./$< --mkl A.mkl
endif
	( ! ./$< --unrecognized-option-triggers-abort )
	( ./$< $(RTP) $(BTP) --skip-loading-symmetric-matrices A.mtx | $(IBC) | grep -q skip )
	./$< $(RTP) A_non_existent.mtx | $(IBC) | grep -q problems.opening
	( $(RSBT_HPP_PC) ./$< --render-only $(BTP) A.mtx | $(IBC) | grep -q Rendering )
	( $(RSBT_HPP_PC) ./$< $(RTP) --quiet --types all --nthreads 1,2 --maxtimes 1 -+ A.mtx | $(IBC) | grep -q 2.threads )
	( $(RSBT_HPP_PC) ./$< $(RTP) --quiet --render-only A.mtx > $(RSBT_DEV_NULL) )
	$(RSBT_HPP_PC) ./$< $(RTP) --quiet --max 1 --nrhs 1 --beta 1 --incy 1 --incx 1 --render --no-trans --alpha 1 --type all A.mtx
	( ./$< $(STP) --quiet --skip-loading-if-less-nnz-matrices 4 A.mtx | $(IBC) | grep -q adding )
	( ./$< $(STP) --quiet --skip-loading-if-less-nnz-matrices 7 A.mtx | $(IBC) | grep -q no.matrix )
	./$< $(RTP) --quiet --skip-loading-not-unsymmetric-matrices A.mtx | $(IBC) | grep -q no.matrix
	(  ./$<  --max 1 --nrhs 1,2 --beta 1 --incy 1,2 --incx 1 --alpha 1 --type all              | $(IBC) | grep -q flop )
	(! ./$<  --max 1 --nrhs 1,2 --beta 1 --incy 1,2 --incx 1 --alpha 1 --type all --no-timings | $(IBC) | grep -q flop )
if !WANT_CUSTOM_NUMERICAL_TYPES
	./$< --max 30 --quiet --symmetric   | $(IBC) | grep -q '\<77760.*success'
	./$< --max 30 --quiet --unsymmetric | $(IBC) | grep -q '\<64800.*success'
	./$< --max 30 --quiet --hermitian   | $(IBC) | grep -q '\<63936.*success'
endif
if HAVE_WHICH
	if which rsbench && rsbench -v | $(IBC) | grep -q version ; then rsbench --generate-matrix -r 4 -c 4 -b 4 > test.mtx ; fi
endif
if !WANT_CUSTOM_NUMERICAL_TYPES
if HAVE_WHICH
	if test -f test.mtx && which rsbench -v; then ./$< --no-tune test.mtx --maxtimes 2 | $(IBC) | grep -q '^..\<720.*success' ; fi
endif
	./$< --max 1 --nrhs 1 --beta 1 --incy 1 --incx 1 --no-rectangular | $(IBC) | grep -q '\<144.*success'
endif
	./$< --max 1 --nrhs 1 --beta 1 --incy 1 --incx 1 --no-rectangular --no-trans --alpha 1 --type d  --transA Q | $(IBC) | grep -q no.valid
if !WANT_CUSTOM_NUMERICAL_TYPES
	./$< --max 1 --nrhs 1 --beta 1 --incy 1 --incx 1 --no-rectangular --no-trans --alpha 1 --type d  --transA C | $(IBC) | grep -q '\<2.*success'
	./$< --quiet --types d  --csr-tester --max 30 --no-rand | $(IBC) | grep -q '\<10800.*success'
	./$< --quiet --types d --hcsr-tester --max 30 --no-rand | $(IBC) | grep -q '\<10800.*success'
	./$< --quiet --types d  --coo-tester --max 30 --no-rand | $(IBC) | grep -q '\<17712.*success'
	./$< --quiet --types d --hcoo-tester --max 30 --no-rand | $(IBC) | grep -q '\<17712.*success'
	./$< --quiet --types d --half-indices-tester --max 30 --no-rand | $(IBC) | grep -q '\<17712.*success'
	./$< --quiet --types d --rsb-tester --max 30 --no-rand | $(IBC) | grep -q '\<17712.*success'
	./$< --quiet --types d --recursive-tester --max 30 --no-rand | $(IBC) | grep -q '\<17712.*success'
	./$< --quiet --types d --no-leaf-multivec --nrhs 1 --max 30 --no-rand | $(IBC) | grep -q '\<11808.*success' # note: for nrhs=1 this is bogus
	./$< --quiet --types d --no-leaf-multivec --nrhs 2 --max 30 --tolerate-mismatch | $(IBC) | grep -q 5904 # not bogus, but unfinished feature
	./$< --quiet --types d --max-random-val 100000 --nrhs 2 --max 30 | $(IBC) | grep -q 5904.*successfu
endif
	./$< --quiet --skip-except-every-random-n-test-cases 1000000000 --max 1 | $(IBC) | grep -q 1.compa
	./$< --quiet --skip-except-every-n-test-cases 1000000000 --max 1 | $(IBC) | grep -q 1.compa
	./$< --quiet $(RTP) --skip-loading-if-more-nnz-matrices 5 A.mtx  | $(IBC) | grep -q no.matri
if HAVE_TIMEOUT
	timeout 2 ./$< --max-test-time 1 -q                   | $(IBC) | grep -q '.*success'
endif
if !WANT_CUSTOM_NUMERICAL_TYPES
	./$< --max 1 --coo-stretch-factor 1,2 --quiet         | $(IBC) | grep -q '\<36288.*success'
	./$< --max 1 --coo-stretch-factor 1,2 --quiet --density 1| $(IBC) | grep -q '\<15552.*success'
	./$< --quiet --incx 1,2,5,20 --incy 1,2,5,20 --max 10 -t all  --band-stretch-factor -1,0,1 | $(IBC) | grep -q 'success'
	./$< --quiet --incx 1,2,5,20 --incy 1,2,5,20 --max 10 -t blas --band-stretch-factor -1,0,1 | $(IBC) | grep -q '\<69984.*success'
	./$< --quiet --incx 1,2,5,20 --incy 1,2,5,20 --rectangular --max 10 -t all  | $(IBC) | grep -q 'successfully.$$'
	./$< --quiet --incx 1,2,5,20 --incy 1,2,5,20 --rectangular --max 10 -t blas | $(IBC) | grep -q '\<46656.*success'
	./$< --quiet --incx 1,2,5,20 --incy 1,2,5,20 --rectangular --max 10 -t blas --fortran-indices   | $(IBC) | grep -q '\<46656.*success'
	./$< --quiet --incx 1        --incy 1        --no-rectangular --max 20 -t all --any-hermitian-ok  | $(IBC) | grep -q 'successfully.$$' # note: this is a bit slow
	./$< --quiet --incx 1,2,5,20 --incy 1,2,5,20 --no-rectangular --max 10 -t blas --diagonal-implicit | $(IBC) | grep -q '\<7776.*success'
	./$< --quiet --incx 1,2,5,20 --incy 1,2,5,20 --rectangular --max 10 -t all     --diagonal-implicit | $(IBC) | grep -q 'successfully.$$'
	./$< --quiet --incx 1,2,5,20 --incy 1,2,5,20 --no-rectangular --max 10 -t blas --fortran-indices   | $(IBC) | grep -q '\<7776.*success'
	./$< --quiet --incx 1,2,5,20 --incy 1,2,5,20 --beta -2 --alpha -2 --rectangular --max 10 -t blas| $(IBC) | grep -q '\<3888.*success'
	./$< --quiet --types d --incx 1 --incy 1 --mul 1000000 --inc 1        | $(IBC) | grep -q '\<1512.*success'
	./$< --quiet --types d --incx 1 --incy 1 --mul 100000  --inc 1000000  | $(IBC) | grep -q '\<1512.*success'
	./$< --quiet --types d --only-n-test-cases 17       | $(IBC) | grep -q '\<17.*success'
	./$< --quiet --types d --skip-n-test-cases 4 --only-n-test-cases 17 | $(IBC) | grep -q '\<13.*success'
endif
	./$< --quiet --types d --only-test-case-n  4        | $(IBC) | grep -q '\<1.*success'
	./$< --self-test > $(RSBT_DEV_NULL)
	./$< --quiet --types blas --only-test-case-n  4         | $(IBC) | grep -q success
if !WANT_CUSTOM_NUMERICAL_TYPES
	./$< --quiet --max 3 --inc 1 --mul 2 --nrhs 1,2 -t blas | $(IBC) | grep -q success
	./$< --quiet --max 10 | $(IBC) | grep -q '\<18144.*success'
endif
	./$< --quiet A.mtx $(RTP) --nrhs 1,2 --incx 1,2 --incy 1,2 --report $(RSBT_DEV_NULL) | $(IBC) | grep -q Report
	cp A.mtx A_underscore.mtx # tests LaTeX escaping
	./$< --quiet A_underscore.mtx --report test.tex --no-basename-render $(RTP) --types d --nrhs 1 --incx 1 --incy 1 --no-trans --alpha 1 --beta 1 | $(IBC) | grep LaTeX > $(RSBT_DEV_NULL) && if test $(LATEX) != 'false' ; then latex -interaction=batchmode test.tex ; fi
endif
	@echo " [*] rsbtest tests terminated successfully !"
else
	@echo " [!] Skipping several tests due to unusual configuration."
endif

.PHONY: covrun
covrun: rsbtest$(EXEEXT)
	make check
	make qtests
.PHONY: e
e:
	$(EDITOR) rsbtest.cpp
.PHONY: wc
wc:
	wc *.cpp
.PHONY: t x
t x: tests

.PHONY: check
check: rsbtest$(EXEEXT)
	./$<       $(RTP) --min 1 --max 2
	./$<       $(RTP) --serial --types blas
	./$< A.mtx $(RTP)
if WITH_CXX17_FILESYSTEM
	./$< A.mtx $(RTP) . # iterate matrices in this directory
endif

EXTRA_DIST=A.mtx

.PHONY: coverage
coverage: rsbtest$(EXEEXT)
	make clean
	./configure \
		CXXFLAGS='--coverage -O0 -g -Wcoverage-mismatch -Wall -pedantic ' \
		LIBS='-lgcov' \
		`grep ^LIBRSB_CONFIG Makefile | sed 's/ //g' || true`
	make
	lcov  --directory `pwd` --zerocounters --no-external # deletes *.gcda
	make covrun
#	for f in *.o ; do gcov -f $${f/.o/}  ; done # STL
	gcov -f rsbtest.o
	rm -f *.info
	lcov --capture --directory `pwd` --output-file coverage.info --no-external
	lcov -a coverage.info  -o coverage-total.info
	genhtml coverage-total.info --highlight --legend --no-branch-coverage --function-coverage --branch-coverage  --output-directory coverage-info-dir
if HAVE_WHICH
	if which lynx  ; then lynx  -dump coverage-info-dir/index.html; fi
	if which links ; then links -dump coverage-info-dir/index.html | head -n 30 ; fi
endif

.PHONY: realclean-local
realclean-local: clean
	rm -f *.aux *.dvi *.log *.nav *.out *.snm *.tex *.toc *.vrb \
		*.pdf *.eps
	for f in *.mtx ; do if test $$f != A.mtx; then rm -f $$f; fi; done

.PHONY: distclean-local
distclean-local: realclean-local
	true
