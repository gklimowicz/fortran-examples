#!/bin/bash
#
# Copyright (C) 2008-2022 Michele Martone
# 
# This file is part of librsb.
# 
# librsb is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# librsb is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with librsb; see the file COPYING.
# If not, see <http://www.gnu.org/licenses/>.

# This script is intended for the librsb developer usage.

set -e
set -x
#rm -f *.gcda        *.gcov
lcov           --directory `pwd` --zerocounters --no-external
NOFLUSH='--no-flush-cache-in-iterations --no-flush-cache-around-loop --want-no-memory-benchmark'

if test "$RSB_GCOV_RUN_EXPRESS" != "1" ; then
# misc fast rsbench runs
./rsbench -oa -Ob $NOFLUSH -f pd.mtx --only-lower-triangle --want-no-autotune --times 1 # want_only_lowtri
./rsbench -oa -Ob $NOFLUSH -f pd.mtx --only-upper-triangle --want-no-autotune --times 1 # want_only_upptri
./rsbench -oa -Ob $NOFLUSH -f pd.mtx --times 1 --pre-transpose --matrix-sample-pcnt 50 --no-reuse-io-arrays # want_transpose mtx_sample_rate
./rsbench -oa -Ob $NOFLUSH -f pd.mtx --want-io-only #
./rsbench -oa -Ob $NOFLUSH --dense 10x4 --times 1 # should_generate_dense_nc
./rsbench -oa -Ob $NOFLUSH --dense 4x4 --times 1 --as-hermitian   #
./rsbench -oa -Ob $NOFLUSH --dense 4x4 --times 1 --as-symmetric   #
./rsbench -oa -Ob $NOFLUSH --dense 4x4 --times 1 --as-unsymmetric #
./rsbench -oa -Ob $NOFLUSH --dense 4x4 --times 1 --want-no-ones-fill
./rsbench -oa -Ob $NOFLUSH --dense 4x4 --times 1 --want-autotune=0.1s2xVV # --want-autotune wants '='
./rsbench -oa -Ob $NOFLUSH --dense 4x4 --times 1 --all-flags --dumpout --allow-any-transposition-combination #
./rsbench -oa -Ob $NOFLUSH --dense 4x4 --times 1 --nrhs 1,2 --by-rows --echo-arguments --impatient --matrix-time --repeath-constructor 2 # misc
./rsbench -oa -Ob $NOFLUSH --dense 4x4 --times 1 --nrhs 1,2 --by-columns --echo-arguments --impatient --matrix-time --repeath-constructor 2 # misc

export RSB_MERGE_SF=2.0 RSB_SPLIT_SF=2.0 # TODO: better write unit test for rsb__scale_subm_idx_on_env_var
./rsbench -oa -Ob --dense 10 \
	--update -t 2 \
	--want-print-per-subm-stats --want-getdiag-bench --guess-blocking \
	--in-place-assembly-experimental --want-ancillary-execs --dump-n-lhs-elements 1 \
	--diagonal-dominance-check  \
	--verbose --verbose \
	--all-formats \
	--all-blas-types --one-nonunit-incx-incy-nrhs-per-type $NOFLUSH
rm -f dense-10x10-100nz--?-N-1--base.eps
OMP_NUM_THREADS=1 ./rsbench --plot-matrix --no-submatrix-format-labels -aRd  -f pd.mtx | wc

if ./rsbench -C  | grep -q OpenMP.*on; then
	RSB_USE_HOSTNAME=0 \
	OMP_NUM_THREADS=1 RSB_WANT_SPMV_TRACE=1 ./rsbench -oa -Ob -R cs.mtx --times 1 --type d --notranspose
	! grep -q -F "$HOSTNAME" spmv-trace.eps
	rm spmv-times.eps
	rm spmv-trace.eps
	rm spmv-frame-????.eps
fi
./rsbench --plot-matrix --latex -f cs.mtx | wc
./rsbench --plot-matrix -aRdN -f cs.mtx | wc

./rsbench -oa -Ob --lower 10 -t 1 --less-verbose --less-verbose --less-verbose --all-formats --want-no-autotune $NOFLUSH --all-blas-opts # for coverage of rsb__util_sort_row_major_parallel # perche non copre rsb__BCOR_spmv_sxsa_double_complex_C__tT_r1_c1_uu_sH_dE_uG ?
( set +e ; ./rsbench  -Q10Q & sleep 1 ; kill -INT  `jobs -p`; wait; ) # rsb__sigh
./rsbench -oa -Ob $NOFLUSH --skip-loading-if-less-nnz-matrices 1 --skip-loading-if-more-nnz-matrices 1 --skip-loading-if-matching-substr A.mtx . # rsb__adddir (--skip-loading-if-less-nnz-matrices=1 is to skip the pd.mtx.bin)
fi # RSB_GCOV_RUN_EXPRESS
make qqtests
if test "$RSB_GCOV_RUN_EXPRESS" != "1" ; then
#scripts/devtests.sh
if ./rsbench -oa -Ob $NOFLUSH --setenv RSB_DUMMY=1 || ./rsbench -oa -Ob $NOFLUSH --setenv 'RSB_DUMMY=1' 2>&1 | grep setenv; then true ; else false; fi
if RSB_DUMMY=3 ./rsbench -oa -Ob $NOFLUSH --unsetenv RSB_DUMMY --dense 1 | grep 'RSB_DUMMY=' ; then false ; fi
./rsbench -M || true # it returns non-zero for specific reasons
./rsbench --limits-testing
./rsbench -oa -Ob $NOFLUSH --dense 1  --generate-spacing 2 # rsb__util_coo_array_mul, ...
./rsbench -oa -Ob $NOFLUSH -f us.mtx # rsb__util_coo_upper_to_lower_symmetric

RSB_PR_MBW=1 ./rsbench -oa -Ob --want-memory-benchmark --dense 1 --write-performance-record=test.rpr # rsb__mbw_es_print
./rsbench -oa -Ob -nmb --dense 1 --column-expand -2 --write-no-performance-record
./rsbench -oa -Ob -nmb --dense 1 --column-expand  2 --write-no-performance-record
./rsbench -oa -Ob -nmb --dense 1 --column-expand  2 --write-no-performance-record
./rsbench -oa -Ob -nmb --dense 1 --column-expand  2 --write-no-performance-record -ADs # rsb__sample_program_options_get_flags
rm test.rpr

./rsbench -C | grep -q ZLIB
if ./rsbench -C  | grep -q ZLIB.*on; then
	test -f pg.mtx.gz || gzip < pg.mtx > pg.mtx.gz
	./rsbench -oa -Ob $NOFLUSH -f pg.mtx # test pattern I/O
fi

./rsbench -oa -Ob  --dense 1 --bounded-box 0 # all rsb_do_compute_bounded_boxes
./rsbench --generate-matrix -r 100001 -c 100001 -n 1024 >  /dev/shm/rsb_matrix.mtx
./rsbench -oa -Ob --bench $NOFLUSH -f  /dev/shm/rsb_matrix.mtx --all-blas-opts # for coverage of rsb__util_sort_row_major_parallel
./rsbench --plot-matrix                  -f /dev/shm/rsb_matrix.mtx
./rsbench --plot-matrix -d --ussv-dump   -f /dev/shm/rsb_matrix.mtx
# TODO: need something like : but on --lower 40 as-symmetric/--implicit-diagonal/--all-transposes/formats/-qH/-Fo
#./rsbench -oa -Ob         -f  /dev/shm/rsb_matrix.mtx --want-no-memory-benchmark --types : --inc : --alpha : --beta :
#./rsbench -oa -Ob         -f  /dev/shm/rsb_matrix.mtx --want-no-memory-benchmark --types : --inc : --alpha : --beta :
#./rsbench -oa -Ob --as-symmetric -Fo --lower 40 --want-no-memory-benchmark --types : --inc : --alpha : --beta :
#
#
#./rsbench -oa -Ob --lower 40 --as-symmetric -Fo --types : --beta 1 --alpha 1  --want-no-memory-benchmark  -t 1 
#./rsbench -oa -Ob --lower 40 --as-symmetric -Fo --types : --beta 1 --alpha 1  --want-no-memory-benchmark  -t 1 -qH
#./rsbench -oa -Ob --lower 40 --as-symmetric -Fo --types : --beta 1 --alpha 1  --want-no-memory-benchmark  -t 1 --implicit-diagonal     --all-transposes
#./rsbench -oa -Ob --lower 40 --as-symmetric -Fo --types : --beta 1 --alpha 1  --want-no-memory-benchmark  -t 1 -qH --implicit-diagonal --all-transposes
#
#
#./rsbench -oa -Ob --lower 40                -Fo --types : --want-no-memory-benchmark  -t 1 -qH                     --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40                -Fo --types : --want-no-memory-benchmark  -t 1                         --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40                    --types : --want-no-memory-benchmark  -t 1 -qH                     --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40                    --types : --want-no-memory-benchmark  -t 1                         --all-transposes --alpha : --beta : --inc :
##
#./rsbench -oa -Ob --lower 40 --as-symmetric -Fo --types : --want-no-memory-benchmark  -t 1 -qH                     --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40 --as-symmetric -Fo --types : --want-no-memory-benchmark  -t 1                         --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40 --as-symmetric     --types : --want-no-memory-benchmark  -t 1 -qH                     --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40 --as-symmetric     --types : --want-no-memory-benchmark  -t 1                         --all-transposes --alpha : --beta : --inc :
##
##
#./rsbench -oa -Ob --lower 40                -Fo --types : --want-no-memory-benchmark  -t 1 -qH --implicit-diagonal --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40                -Fo --types : --want-no-memory-benchmark  -t 1     --implicit-diagonal --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40                    --types : --want-no-memory-benchmark  -t 1 -qH --implicit-diagonal --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40                    --types : --want-no-memory-benchmark  -t 1     --implicit-diagonal --all-transposes --alpha : --beta : --inc :
##
#./rsbench -oa -Ob --lower 40 --as-symmetric -Fo --types : --want-no-memory-benchmark  -t 1 -qH --implicit-diagonal --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40 --as-symmetric -Fo --types : --want-no-memory-benchmark  -t 1     --implicit-diagonal --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40 --as-symmetric     --types : --want-no-memory-benchmark  -t 1 -qH --implicit-diagonal --all-transposes --alpha : --beta : --inc :
#./rsbench -oa -Ob --lower 40 --as-symmetric     --types : --want-no-memory-benchmark  -t 1     --implicit-diagonal --all-transposes --alpha : --beta : --inc :

./rsbench -oa -Ob $NOFLUSH -R  --dense 2 --want-getrow-bench --want-print-per-subm-stats # --want-getrow-bench
./rsbench -oa -Ob $NOFLUSH -R  --dense 2 --want-nonzeroes-distplot # obsolescent
./rsbench -oa -Ob $NOFLUSH -R  --dense 2 --want-unordered-coo-test # obsolescent
./rsbench -oa -Ob $NOFLUSH --all-formats --dense 4 --zig-zag  --less-verbose --less-verbose  --less-verbose # coverage of rsb__do_reverse_odd_rows
./rsbench -oa -Ob $NOFLUSH --lower 4 -K # coverage of rsb__util_is_sorted_coo_as_row_major 
./rsbench -oa -Ob $NOFLUSH --lower 4 --z-sorted-coo     # coverage of rsb__do_zsort_coo_submatrices
./rsbench -oa -Ob $NOFLUSH --lower 4 --z-sorted-coo -qH # coverage of rsb__do_zsort_coo_submatrices
./rsbench -I # coverage of rsb_print_mop_maxmins incomplete
./rsbench -Q0.1
./rsbench -oa -Ob $NOFLUSH --lower 4 --ilu0 --types : # coverage of rsb__prec_ilu0
./rsbench -ot -Ob $NOFLUSH --lower 4 --want-ancillary-execs --only-lower-triangle --times 1 --want-no-autotune # rsb_do_spsv_recursive_serial
./rsbench -oa -Ob $NOFLUSH -n1 -R -qH --dense 2 --want-ancillary-execs --times 1 --want-no-autotune # rsb__do_spmv_recursive_serial

 RSB_VERBOSE_TUNING=2 ./rsbench -ot -Ob --lower 10 --want-autotune=0.1s2xVV # coverage of rsb__tattr_dump and verbose tuning
rm lower-10x10-55nz--D-N-1--base.eps
rm lower-10x10-55nz--D-N-1--sv-tuning_trace.gnu
rm lower-10x10-55nz--D-N-1--sv-tuning_trace.dat

./rsbench --plot-matrix -aRzd -f pd.mtx.bin # coverage of rsb_mio.c
./rsbench --matrix-ls       pg.mtx # coverage of rsb_mio.c
./rsbench --matrix-ls       pd.mtx # coverage of rsb_mmls.c
./rsbench --matrix-ls-latex pd.mtx # coverage of rsb_mmls.c
./rsbench -E 0.1s # coverage of rsb_failure_tests.c
./rsbench -oa -Ob $NOFLUSH pd.mtx --alternate-sort 2 -qZ || true # coverage of rsb_msort_up.c
./rsbench -g -r 10 -c 10 -D # coverage of rsb__generate_banded
./rsbench -oa -Ob $NOFLUSH --skip-loading-if-matching-regex cs cs.mtx pd.mtx # coverage of rsb_regexp_match
./rsbench -oS -Ob $NOFLUSH pd.mtx # coverage of rsb__main_block_partitioned_mat_stats
RSB_FPBENCH_MULTITYPE_TIME=0.1 ./rsbench -F || true # coverage of rsb_fpb.c
RSB_BENCHMARK_MIN_SECONDS=0.001  ./rsbench -Oc -f pd.mtx # coverage via rsb__do_completebenchmark
#
./rsbench --generate-matrix -r 100 -c 100 -l 3 | sed s/general/symmetric/g > /dev/shm/rsb_matrix.mtx
RSB_BENCHMARK_MIN_SECONDS=0.001  ./rsbench -Oc -f /dev/shm/rsb_matrix.mtx
#
RSB_BENCHMARK_MIN_SECONDS=0.001  ./rsbench -Oc -f cs.mtx # symmetric coverage via rsb__do_completebenchmark
RSB_BENCHMARK_MIN_SECONDS=0.001  ./rsbench -Oc -f ch.mtx # hermitian coverage via rsb__do_completebenchmark
RSB_BENCHMARK_MIN_SECONDS=0.001  ./rsbench -Oc -f lt.mtx # spsv coverage via rsb__do_completebenchmark
RSB_BENCHMARK_MIN_SECONDS=0.001  ./rsbench -Oc -f ut.mtx # spsv coverage via rsb__do_completebenchmark
PATH="`pwd`:${PATH}" scripts/benchmark.sh --nrhs 1,2,4,8 pd.mtx cs.mtx --want-no-memory-benchmark --all-transposes
TC=`./rsbench -C  | grep types.count`
TC=${TC/types count:/}
test -n ${TC}
TTC=$((2*TC))
./rsbench -oa -Ob --dense 10 --bench --notranspose    --want-no-memory-benchmark --write-performance-record=tmp.rpr | grep record.of.${TC}.samples # more of a test than coverage
./rsbench -oa -Ob --dense 10 --bench --also-transpose --want-no-memory-benchmark --write-performance-record=tmp.rpr | grep record.of.${TTC}.samples # more of a test than coverage
./rsbench -oa -Ob --dense 10 --bench --transpose      --want-no-memory-benchmark --write-performance-record=tmp.rpr | grep record.of.${TC}.samples # more of a test than coverage
./rsbench -oa -Ob --dense 10 --bench --transpose-as N --want-no-memory-benchmark --write-performance-record=tmp.rpr | grep record.of.${TC}.samples # more of a test than coverage
fi # RSB_GCOV_RUN_EXPRESS

RSB_SHORT_TEST_SH=1 scripts/test.sh
make sbtc && ./sbtc
make sbtf && ./sbtf
make ot   && ./ot
if test -d librsbpp; then
	make covrun -C librsbpp
fi
if test -d rsblib; then
	make covrun -C rsblib
fi
if test -d rsbtest; then
	make covrun -C rsbtest
fi
for f in *.o ; do gcov -f ${f/.o/}  ; done
cd examples
#rm -f *.gcda        *.gcov
make tests
make install
for f in *.o ; do gcov -f ${f/.o/}  ; done
cd -

rm -f *.info
lcov --capture --directory `pwd`         --output-file coverage.info --no-external
lcov --capture --directory `pwd`/examples/ --output-file coverage-examples.info --no-external
lcov  -a coverage.info -a coverage-examples.info  -o coverage-total.info
genhtml coverage-total.info --highlight --legend --no-branch-coverage --function-coverage --branch-coverage  --output-directory coverage-info-dir
ls -l coverage-info-dir
which lynx && lynx -dump coverage-info-dir/index.html
which links && links -dump coverage-info-dir/librsb_coverage/index-sort-f.html | head -n 30
echo "[*] Coverage test performed in ${SECONDS}s." 
echo "[*] At next 'make clean', remember to rm -f *.gcov *.gcno" 
