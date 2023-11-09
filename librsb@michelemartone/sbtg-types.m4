dnl
include(`libspblas_macros.m4')dnl
include(`rsb_fortran_macros.m4')dnl
dnl
`#	'Supported types  :RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES
`#	'Unsupported types:RSB_M4_SPBLAS_MATRIX_UNSUPPORTED_TYPES
blas_type_codes_array_=[dnl
dnl
foreach(`type',RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES,`dnl
singlequote(RSB_M4_SPBLAS_TYPE_CHARCODE(type))`,'dnl
')dnl
'?'];
global blas_type_codes_array=blas_type_codes_array_(1:length(blas_type_codes_array_)-1);dnl
dnl

