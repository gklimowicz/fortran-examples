/*

Copyright (C) 2008-2021 Michele Martone

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
/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief
 * This source file contains matrix I/O functions.
 * */

#ifndef RSB_IO_H_INCLUDED
#define RSB_IO_H_INCLUDED

/*#include "rsb_internals.h"*/
/*#include "rsb.h"*/
#include "rsb_common.h"

#ifdef RSB_WITH_MM
rsb_err_t rsb__util_mm_load_matrix_f(const char *fn, rsb_coo_idx_t ** IA, rsb_coo_idx_t ** JA, void **VA, rsb_coo_idx_t *m, rsb_coo_idx_t *k , rsb_nnz_idx_t *nnz, rsb_type_t typecode, rsb_flags_t flags, rsb_bool_t *is_lowerp, rsb_bool_t *is_upperp);
rsb_err_t rsb__util_mm_load_vector_f(const char *fn, void **VA, rsb_nnz_idx_t *nnz, rsb_type_t typecode);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__util_mm_load_matrix_f_as_csr(const char *fn, rsb_nnz_idx_t ** INDX, rsb_coo_idx_t ** JA, void **VA, rsb_coo_idx_t *m, rsb_coo_idx_t *k , rsb_nnz_idx_t *nnz, rsb_type_t typecode, rsb_flags_t flags);
rsb_err_t rsb__util_mm_load_matrix_f_as_csc(const char *fn, rsb_nnz_idx_t ** INDX, rsb_coo_idx_t ** IA, void **VA, rsb_coo_idx_t *m, rsb_coo_idx_t *k , rsb_nnz_idx_t *nnz, rsb_type_t typecode, rsb_flags_t flags);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
rsb_err_t rsb__do_util_get_matrix_dimensions(const char * filename, size_t * cols, size_t * rows, size_t * nnzp, rsb_flags_t*flagsp);
rsb_err_t rsb__util_mm_info_matrix_f(const char *fn,  rsb_coo_idx_t *m, rsb_coo_idx_t *k , rsb_nnz_idx_t *nnz, rsb_type_t *typecode, rsb_bool_t * is_symmetric, rsb_bool_t * is_hermitian, rsb_bool_t * is_pattern, rsb_bool_t * is_lower, rsb_bool_t * is_upper , rsb_bool_t * is_vector );
rsb_err_t rsb__util_mm_load_coo_matrix(const char *filename, struct rsb_coo_mtx_t * cmp);
#ifdef RSB_OBSOLETE_QUARANTINE_UNUSED
rsb_err_t rsb__util_mm_fill_arrays_for_csc(const char *filename, rsb_nnz_idx_t * INDX, rsb_coo_idx_t * IA, void *VA, rsb_type_t typecode, rsb_flags_t flags);
#endif /* RSB_OBSOLETE_QUARANTINE_UNUSED */
size_t rsb__sys_filesize(const char *filename);
int rsb__fscanf(FILE * fd,const char * fs,rsb_coo_idx_t *IV, rsb_coo_idx_t *JV, void * VAR, void * VAI);
char * rsb__fgets(char* RSB_RESTRICT buf, int len, FILE * RSB_RESTRICT fd);
rsb_err_t rsb__do_file_mtx_get_dims(const char * filename, rsb_coo_idx_t* nrp, rsb_coo_idx_t *ncp, rsb_coo_idx_t *nzp, rsb_flags_t*flagsp);
#if RSB_WANT_EXPERIMENTAL_BINARY_COO
int rsb_getc(FILE *fd, void *ngzfd);
int rsb_ungetc(int cc, FILE *fd, void *ngzfd);
size_t rsb_fread(void* buf, size_t size, size_t nitems, void * fd, void *ngzfd);
rsb_err_t rsb__read_coo_bin_fd(FILE *fd, void *ngzfd, rsb_coo_idx_t * IA, rsb_coo_idx_t * JA, void *VA, const rsb_coo_idx_t nr, const rsb_coo_idx_t nc, const rsb_nnz_idx_t nnz, const rsb_type_t typecode);
#endif /* RSB_WANT_EXPERIMENTAL_BINARY_COO */
#endif /* RSB_WITH_MM */
#endif /* RSB_IO_H_INCLUDED */
/* @endcond */
