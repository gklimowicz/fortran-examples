dnl
include(`rsbpp.m4')dnl
dnl
dnl
#ifndef RSBCK_HPP_INCLUDED
#define RSBCK_HPP_INCLUDED

#include "rsbpp.hpp"

dnl
foreach(`IT',`(rsb_coo_idx_t)',`dnl
foreach(`CT',`(rsb_half_idx_t,rsb_coo_idx_t)',`dnl
foreach(`NT',(`RSBPP_M4_IMPLEMENTED_TYPES'),`dnl
pushdef(`CXX_NT',RSB_M4_MATRIX_TYPE_AS_CXX(NT))`'dnl
// <CXX_NT,IT,CT>
extern template rsb_err_t rsbpp_csr_spmx<CXX_NT,IT,CT>(rsb_flags_t flags, const IT nnz, const IT nr, const IT nc, const CXX_NT* VA, const IT* IP, const CT* JA, const IT nrhs, const IT ldX, const CXX_NT* rhs, const IT ldY, CXX_NT* out, const CXX_NT* alphap, IT incx, IT incy, const rsb_trans_t transA, IT roff, IT coff, bool by_rows);
popdef(`CXX_NT')`'dnl
')dnl
')dnl
')dnl
dnl

dnl
foreach(`IT',`(rsb_coo_idx_t)',`dnl
foreach(`CT',`(rsb_half_idx_t,rsb_coo_idx_t)',`dnl
foreach(`NT',(`RSBPP_M4_IMPLEMENTED_TYPES'),`dnl
pushdef(`CXX_NT',RSB_M4_MATRIX_TYPE_AS_CXX(NT))`'dnl
// <CXX_NT,IT,CT>
extern template rsb_err_t rsbpp_coo_spmx<CXX_NT,IT,CT>(rsb_flags_t flags, const IT nnz, const IT nr, const IT nc, const CXX_NT* VA, const CT* IA, const CT* JA, const IT nrhs, const IT ldX, const CXX_NT* rhs, const IT ldY, CXX_NT* out, const CXX_NT* alphap, IT incx, IT incy, const rsb_trans_t transA, IT roff, IT coff, bool by_rows);
popdef(`CXX_NT')`'dnl
')dnl
')dnl
')dnl
dnl
dnl

#endif /* RSBCK_HPP_INCLUDED */
dnl
