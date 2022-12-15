dnl	Code generator configuration
dnl	Michele Martone
dnl
dnl
dnl	---------------------------------------------------------------------------
dnl	Whether we want OpenMP thread level parallelism (EXPERIMENTAL)
dnl
dnl define(`RSB_M4_WANT_OMP',`ifelse(`@libmmvbr_cv_openmp@',`yes',`1',`')')dnl
define(`RSB_M4_WANT_OMP',`ifelse(`yes',`yes',`1',`')')dnl
define(`RSB_M4_WANT_OMP_IN_RECURSION',`ifelse(`yes',`yes',`1',`')')dnl
define(`RSB_M4_WANT_OMP_IN_KERNELS',`ifelse(`yes',`yes',`',`')')dnl
dnl define(`RSB_M4_WANT_OMP_IN_KERNELS',`ifelse(`yes',`yes',`1',`')')dnl
define(`RSB_M4_MAX_OMP_THREADS',`4')dnl	FIXME : TEMPORARY 
dnl
dnl	---------------------------------------------------------------------------
dnl	Generate code with m4 debug info in it.
dnl
define(`RSB_M4_DEBUG',`1')dnl
dnl
dnl	---------------------------------------------------------------------------
dnl	If 1, enables register blocking, in kernels where this is supported (experimental).
dnl
define(`RSB_M4_WANT_BLOCKING',`1')dnl
dnl
dnl	---------------------------------------------------------------------------
dnl	The number of registers, in case of register blocking (EXPERIMENTAL).
dnl
define(`RSB_M4_REGISTERS',`8')dnl
dnl
dnl	---------------------------------------------------------------------------
define(`RSB_M4_FITTING_SAMPLES',/*12 8*/4)dnl
dnl
dnl	---------------------------------------------------------------------------
define(`RSB_M4_FORTRAN_CONVENTION',`gfortran')dnl
dnl
dnl	---------------------------------------------------------------------------
define(`RSB_M4_BENCHMARK_MIN_SECONDS',/*0.5*/1.0)dnl
dnl
define(`RSB_M4_BENCHMARK_MIN_RUNS',/*5*/10)dnl
dnl
dnl	---------------------------------------------------------------------------
define(`RSB_M4_BUFLEN',128)dnl
dnl
dnl	---------------------------------------------------------------------------
dnl	Implemented types (from 1.3) (same values as RSB_M4_TYPES)
define(`RSB_M4_IMPLEMENTED_TYPES',`double,float,float complex,double complex')dnl
dnl
dnl	---------------------------------------------------------------------------
define(`RSB_M4_USE_RESTRICT',`ifelse(`yes',`yes',`1',`')')dnl
dnl
dnl	---------------------------------------------------------------------------
dnl	Version strings.
dnl
define(`RSB_M4_WANT_LIBRSB_VER_DATE',`202210181826')dnl
define(`RSB_M4_WANT_LIBRSB_VER_MAJOR',`1')dnl
define(`RSB_M4_WANT_LIBRSB_VER_MINOR',`3')dnl
define(`RSB_M4_WANT_LIBRSB_LIBRSB_VER',`10300')dnl
define(`RSB_M4_WANT_LIBRSB_VER_PATCH',`0')dnl
define(`RSB_M4_WANT_LIBRSB_VER_PRERS',`.1')dnl
dnl
dnl	---------------------------------------------------------------------------
define(`RSB_M4_WANT_20110206_BOUNDED_BOX_PATCH',1)dnl
dnl	---------------------------------------------------------------------------
