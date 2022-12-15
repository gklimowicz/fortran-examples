/*
Copyright (C) 2020-2022 Michele Martone

This file is part of librsb.

librsb is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

librsb is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with librsb; see the file COPYING.
If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef LIBRSBPP_H_INCLUDED
#define LIBRSBPP_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

/* These functions shall be internal functions of librsb, and therefore be prefixed by 'rsb__'. */
#define rsbpp_coo_spmv rsb__pp_coo_spmv
#define rsbpp_coo_spmm rsb__pp_coo_spmm
#define rsbpp_csr_spmv rsb__pp_csr_spmv
#define rsbpp_csr_spmm rsb__pp_csr_spmm

rsb_err_t rsbpp_coo_spmv(rsb_type_t typecode, rsb_flags_t flags, const rsb_nnz_idx_t nnz, const rsb_nnz_idx_t nr, const rsb_nnz_idx_t nc, const void* VA, const void* IA, const void* JA, const void * rhs, void * out, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff);
rsb_err_t rsbpp_coo_spmm(rsb_type_t typecode, rsb_flags_t flags, const rsb_nnz_idx_t nnz, const rsb_nnz_idx_t nr, const rsb_nnz_idx_t nc, const void* VA, const void* IA, const void* JA, rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const void * rhs, const rsb_coo_idx_t ldY, void * out, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, const rsb_flags_t order);
rsb_err_t rsbpp_csr_spmv(rsb_type_t typecode, rsb_flags_t flags, const rsb_nnz_idx_t nnz, const rsb_nnz_idx_t nr, const rsb_nnz_idx_t nc, const void* VA, const void* IP, const void* JA, const void * rhs, void * out, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff);
rsb_err_t rsbpp_csr_spmm(rsb_type_t typecode, rsb_flags_t flags, const rsb_nnz_idx_t nnz, const rsb_nnz_idx_t nr, const rsb_nnz_idx_t nc, const void* VA, const void* IP, const void* JA, rsb_coo_idx_t nrhs, const rsb_coo_idx_t ldX, const void * rhs, const rsb_coo_idx_t ldY, void * out, const void * alphap, rsb_coo_idx_t incx, rsb_coo_idx_t incy, const rsb_trans_t transA, rsb_coo_idx_t roff, rsb_coo_idx_t coff, const rsb_flags_t order);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* LIBRSBPP_H_INCLUDED */
