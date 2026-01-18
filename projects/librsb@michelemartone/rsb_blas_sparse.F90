! 
! Copyright (C) 2008-2021 Michele Martone
! 
! This file is part of librsb.
! 
! librsb is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published
! by the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
! 
! librsb is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
! License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with librsb; see the file COPYING.
! If not, see <http://www.gnu.org/licenses/>.
! 

!
!> @file
!! @brief Implementation of the Fortran Sparse BLAS interface to \librsb (see \ref rsb_doc_sparse_blas).
!!
!!

#define RSB_HAVE_RSB_KERNELS 1

      MODULE blas_sparse
        !> A Sparse BLAS interface for RSB
        IMPLICIT NONE
PUBLIC

        
        !> inserts a single entry
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE uscr_insert_entry
        MODULE PROCEDURE suscr_insert_entry &
        &, duscr_insert_entry &
        &, cuscr_insert_entry &
        &, zuscr_insert_entry &
        & ;
        END INTERFACE
        
        !> inserts multiple entries
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE uscr_insert_entries
        MODULE PROCEDURE suscr_insert_entries &
        &, duscr_insert_entries &
        &, cuscr_insert_entries &
        &, zuscr_insert_entries &
        & ;
        END INTERFACE
        
        !> inserts a sparse column
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE uscr_insert_col
        MODULE PROCEDURE suscr_insert_col &
        &, duscr_insert_col &
        &, cuscr_insert_col &
        &, zuscr_insert_col &
        & ;
        END INTERFACE
        
        !> inserts a sparse row
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE uscr_insert_row
        MODULE PROCEDURE suscr_insert_row &
        &, duscr_insert_row &
        &, cuscr_insert_row &
        &, zuscr_insert_row &
        & ;
        END INTERFACE
        
        !> inserts a clique
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE uscr_insert_clique
        MODULE PROCEDURE suscr_insert_clique &
        &, duscr_insert_clique &
        &, cuscr_insert_clique &
        &, zuscr_insert_clique &
        & ;
        END INTERFACE
        
        !> inserts a dense block
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE uscr_insert_block
        MODULE PROCEDURE suscr_insert_block &
        &, duscr_insert_block &
        &, cuscr_insert_block &
        &, zuscr_insert_block &
        & ;
        END INTERFACE
        
        !> multiplication  : c <- beta c + alpha A b
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE usmv
        MODULE PROCEDURE susmv &
        &, dusmv &
        &, cusmv &
        &, zusmv &
        & ;
        END INTERFACE
        
        !> triangular solve: b <- alpha A^-1 b
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE ussv
        MODULE PROCEDURE sussv &
        &, dussv &
        &, cussv &
        &, zussv &
        & ;
        END INTERFACE
        
        !> multiplication  : c <- beta c + alpha A b
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE usmm
        MODULE PROCEDURE susmm &
        &, dusmm &
        &, cusmm &
        &, zusmm &
        & ;
        END INTERFACE
        
        !> triangular solve: b <- alpha A^-1 b
        !> \rsb_spblas_f_istat_msg
        !! 
        INTERFACE ussm
        MODULE PROCEDURE sussm &
        &, dussm &
        &, cussm &
        &, zussm &
        & ;
        END INTERFACE
        
        INTEGER, PARAMETER :: blas_sparse_const_success=0
        INTEGER, PARAMETER :: blas_sparse_const_failure=-1 ! value returned by this interface on failure
        INTEGER, PARAMETER :: blas_sparse_const_not_available=-9999 ! value returned by this interface when deactivated
! This file has been auto-generated from blas_enum.h.
        INTEGER,PARAMETER :: blas_rowmajor=101
        INTEGER,PARAMETER :: blas_colmajor=102
        INTEGER,PARAMETER :: blas_no_trans=111
        INTEGER,PARAMETER :: blas_trans=112
        INTEGER,PARAMETER :: blas_conj_trans=113
        INTEGER,PARAMETER :: blas_upper=121
        INTEGER,PARAMETER :: blas_lower=122
        INTEGER,PARAMETER :: blas_non_unit_diag=131
        INTEGER,PARAMETER :: blas_unit_diag=132
        INTEGER,PARAMETER :: blas_left_side=141
        INTEGER,PARAMETER :: blas_right_side=142
        INTEGER,PARAMETER :: blas_base=151
        INTEGER,PARAMETER :: blas_t=152
        INTEGER,PARAMETER :: blas_rnd=153
        INTEGER,PARAMETER :: blas_ieee=154
        INTEGER,PARAMETER :: blas_emin=155
        INTEGER,PARAMETER :: blas_emax=156
        INTEGER,PARAMETER :: blas_eps=157
        INTEGER,PARAMETER :: blas_prec=158
        INTEGER,PARAMETER :: blas_underflow=159
        INTEGER,PARAMETER :: blas_overflow=160
        INTEGER,PARAMETER :: blas_sfmin=161
        INTEGER,PARAMETER :: blas_one_norm=171
        INTEGER,PARAMETER :: blas_real_one_norm=172
        INTEGER,PARAMETER :: blas_two_norm=173
        INTEGER,PARAMETER :: blas_frobenius_norm=174
        INTEGER,PARAMETER :: blas_inf_norm=175
        INTEGER,PARAMETER :: blas_real_inf_norm=176
        INTEGER,PARAMETER :: blas_max_norm=177
        INTEGER,PARAMETER :: blas_real_max_norm=178
        INTEGER,PARAMETER :: blas_increasing_order=181
        INTEGER,PARAMETER :: blas_decreasing_order=182
        INTEGER,PARAMETER :: blas_conj=191
        INTEGER,PARAMETER :: blas_no_conj=192
        INTEGER,PARAMETER :: blas_jrot_inner=201
        INTEGER,PARAMETER :: blas_jrot_outer=202
        INTEGER,PARAMETER :: blas_jrot_sorted=203
        INTEGER,PARAMETER :: blas_prec_single=211
        INTEGER,PARAMETER :: blas_prec_double=212
        INTEGER,PARAMETER :: blas_prec_indigenous=213
        INTEGER,PARAMETER :: blas_prec_extra=214
        INTEGER,PARAMETER :: blas_zero_base=221
        INTEGER,PARAMETER :: blas_one_base=222
        INTEGER,PARAMETER :: blas_general=231
        INTEGER,PARAMETER :: blas_symmetric=232
        INTEGER,PARAMETER :: blas_hermitian=233
        INTEGER,PARAMETER :: blas_triangular=234
        INTEGER,PARAMETER :: blas_lower_triangular=235
        INTEGER,PARAMETER :: blas_upper_triangular=236
        INTEGER,PARAMETER :: blas_lower_symmetric=237
        INTEGER,PARAMETER :: blas_upper_symmetric=238
        INTEGER,PARAMETER :: blas_lower_hermitian=239
        INTEGER,PARAMETER :: blas_upper_hermitian=240
        INTEGER,PARAMETER :: blas_complex=241
        INTEGER,PARAMETER :: blas_real=242
        INTEGER,PARAMETER :: blas_double_precision=243
        INTEGER,PARAMETER :: blas_single_precision=244
        INTEGER,PARAMETER :: blas_num_rows=251
        INTEGER,PARAMETER :: blas_num_cols=252
        INTEGER,PARAMETER :: blas_num_nonzeros=253
        INTEGER,PARAMETER :: blas_invalid_handle=261
        INTEGER,PARAMETER :: blas_new_handle=262
        INTEGER,PARAMETER :: blas_open_handle=263
        INTEGER,PARAMETER :: blas_valid_handle=264
        INTEGER,PARAMETER :: blas_regular=271
        INTEGER,PARAMETER :: blas_irregular=272
        INTEGER,PARAMETER :: blas_block=273
        INTEGER,PARAMETER :: blas_unassembled=274
        INTEGER,PARAMETER :: blas_rsb_spmv_autotuning_on=6660
        INTEGER,PARAMETER :: blas_rsb_spmv_autotuning_off=6661
        INTEGER,PARAMETER :: blas_rsb_spmv_n_autotuning_on=6662
        INTEGER,PARAMETER :: blas_rsb_spmv_n_autotuning_off=6663
        INTEGER,PARAMETER :: blas_rsb_spmv_t_autotuning_on=6664
        INTEGER,PARAMETER :: blas_rsb_spmv_t_autotuning_off=6665
        INTEGER,PARAMETER :: blas_rsb_autotune_next_operation=6666
        INTEGER,PARAMETER :: blas_rsb_rep_rec=9993
        INTEGER,PARAMETER :: blas_rsb_rep_hwi=9994
        INTEGER,PARAMETER :: blas_rsb_rep_rsb=9995
        INTEGER,PARAMETER :: blas_rsb_rep_csr=9996
        INTEGER,PARAMETER :: blas_rsb_rep_coo=9997
        INTEGER,PARAMETER :: blas_rsb_duplicates_ovw=9998
        INTEGER,PARAMETER :: blas_rsb_duplicates_sum=9999

#ifdef RSB_WANT_LONG_IDX_TYPE
        INTEGER,PARAMETER :: RSB_BLAS_IDX_KIND=8
        INTEGER,PARAMETER :: RSB_BLAS_IST_KIND=8 ! for istat
#define C_RSB_INT_KND_ C_INT64_T
#else
        INTEGER,PARAMETER :: RSB_BLAS_IDX_KIND=4
        INTEGER,PARAMETER :: RSB_BLAS_IST_KIND=4 ! for istat
#define C_RSB_INT_KND_ C_INT
#endif

        INTERFACE
          TYPE(C_PTR) FUNCTION &
          &rsb_blas_get_mtx&
          &(A)&
          &BIND(c,NAME = "rsb_blas_get_mtx")
          USE ISO_C_BINDING
          INTEGER(C_INT), VALUE  :: A
          END FUNCTION rsb_blas_get_mtx
        END INTERFACE

CONTAINS

         !> \rsb_spblasl2_ds_msg\rsb_spblas_return_msg
         !> \rsb_spblas_f_istat_msg
         !! 
         
         SUBROUTINE usds(A,istat)
           IMPLICIT NONE
           INTEGER,INTENT(IN)::A
           INTEGER(KIND=RSB_BLAS_IST_KIND)::istat

           istat=blas_sparse_const_success
#if defined(RSB_HAVE_RSB_KERNELS)
           CALL blas_usds(A,istat)
           IF(istat.NE.blas_sparse_const_success)&
            &istat=blas_sparse_const_failure
#else /* RSB_HAVE_RSB_KERNELS */
           istat=blas_sparse_const_not_available
#endif /* RSB_HAVE_RSB_KERNELS */
         END SUBROUTINE
         

         !> \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
         !> \rsb_spblas_f_istat_msg
         !! 
         
         SUBROUTINE uscr_end(A,istat)

           IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
           INTEGER,INTENT(IN)::A

           istat=blas_sparse_const_success
#if defined(RSB_HAVE_RSB_KERNELS)
           CALL blas_uscr_end(A,istat)

           IF(istat.NE.blas_sparse_const_success)&
            &istat=blas_sparse_const_failure
#else /* RSB_HAVE_RSB_KERNELS */
           istat=blas_sparse_const_not_available
#endif /* RSB_HAVE_RSB_KERNELS */
         END SUBROUTINE
         

         !> \rsb_spblasl2_gp_msg\rsb_spblas_return_msg
         !> \rsb_spblas_f_istat_msg
         !! 
         
         SUBROUTINE usgp(A,pname,istat)

           IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
           INTEGER,INTENT(IN)::A
           INTEGER,INTENT(IN)::pname

           istat=blas_sparse_const_success
#if defined(RSB_HAVE_RSB_KERNELS)
           CALL blas_usgp(A,pname,istat)

           !istat does not have the meaning of an error value, here
           !IF(istat.NE.blas_sparse_const_success)istat=blas_sparse_const_failure
#else /* RSB_HAVE_RSB_KERNELS */
           istat=blas_sparse_const_not_available
#endif /* RSB_HAVE_RSB_KERNELS */
         END SUBROUTINE
         

         !> \rsb_spblasl2_sp_msg\rsb_spblas_return_msg
         !> \rsb_spblas_f_istat_msg
         !! 
         
         SUBROUTINE ussp(A,pname,istat)

           IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
           INTEGER,INTENT(IN)::A
           INTEGER,INTENT(IN)::pname

           istat=blas_sparse_const_success
#if defined(RSB_HAVE_RSB_KERNELS)
           CALL blas_ussp(A,pname,istat)

           IF(istat.NE.blas_sparse_const_success)&
            &istat=blas_sparse_const_failure
#else /* RSB_HAVE_RSB_KERNELS */
           istat=blas_sparse_const_not_available
#endif /* RSB_HAVE_RSB_KERNELS */

         END SUBROUTINE
         

        
        
        !> \rsb_spblasl2_cr_begin_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE suscr_begin&
         &(m,n,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: m 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: n 
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_suscr_begin&
           &(m,n,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_begin_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE duscr_begin&
         &(m,n,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: m 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: n 
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_duscr_begin&
           &(m,n,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_begin_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE cuscr_begin&
         &(m,n,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: m 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: n 
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_cuscr_begin&
           &(m,n,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_begin_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE zuscr_begin&
         &(m,n,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: m 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: n 
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_zuscr_begin&
           &(m,n,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_cr_block_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE suscr_block_begin&
         &(Mb,Nb,k,l,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Mb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Nb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: k 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: l 
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_suscr_block_begin&
           &(Mb,Nb,k,l,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_block_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE duscr_block_begin&
         &(Mb,Nb,k,l,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Mb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Nb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: k 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: l 
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_duscr_block_begin&
           &(Mb,Nb,k,l,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_block_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE cuscr_block_begin&
         &(Mb,Nb,k,l,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Mb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Nb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: k 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: l 
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_cuscr_block_begin&
           &(Mb,Nb,k,l,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_block_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE zuscr_block_begin&
         &(Mb,Nb,k,l,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Mb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Nb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: k 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: l 
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_zuscr_block_begin&
           &(Mb,Nb,k,l,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_cr_vbr_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE suscr_variable_block_begin&
         &(Mb,Nb,K,L,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Mb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Nb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: K (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: L (:)
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_suscr_variable_block_begin&
           &(Mb,Nb,K,L,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_vbr_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE duscr_variable_block_begin&
         &(Mb,Nb,K,L,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Mb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Nb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: K (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: L (:)
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_duscr_variable_block_begin&
           &(Mb,Nb,K,L,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_vbr_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE cuscr_variable_block_begin&
         &(Mb,Nb,K,L,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Mb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Nb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: K (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: L (:)
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_cuscr_variable_block_begin&
           &(Mb,Nb,K,L,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_vbr_msg\rsb_spblas_return_mtx_msg
        !> \rsb_spblas_f_istat_msg\rsb_spblasl2_A_msg_ftn

        !! 
        
        SUBROUTINE zuscr_variable_block_begin&
         &(Mb,Nb,K,L,A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Mb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: Nb 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: K (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: L (:)
          INTEGER,INTENT(OUT) :: A


          istat = blas_sparse_const_success
          CALL blas_zuscr_variable_block_begin&
           &(Mb,Nb,K,L,A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE suscr_end&
         &(A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 


          istat = blas_sparse_const_success
          CALL blas_suscr_end&
           &(A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE duscr_end&
         &(A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 


          istat = blas_sparse_const_success
          CALL blas_duscr_end&
           &(A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cuscr_end&
         &(A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 


          istat = blas_sparse_const_success
          CALL blas_cuscr_end&
           &(A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_end_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zuscr_end&
         &(A,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 


          istat = blas_sparse_const_success
          CALL blas_zuscr_end&
           &(A,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE suscr_insert_entry&
         &(A,val,i,j,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          REAL(KIND(1.e0)) :: val 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 


          istat = blas_sparse_const_success
          CALL blas_suscr_insert_entry&
           &(A,val,i,j,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE duscr_insert_entry&
         &(A,val,i,j,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          REAL(KIND(1.d0)) :: val 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 


          istat = blas_sparse_const_success
          CALL blas_duscr_insert_entry&
           &(A,val,i,j,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cuscr_insert_entry&
         &(A,val,i,j,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          COMPLEX(KIND(1.e0)) :: val 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 


          istat = blas_sparse_const_success
          CALL blas_cuscr_insert_entry&
           &(A,val,i,j,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_entry_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zuscr_insert_entry&
         &(A,val,i,j,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          COMPLEX(KIND(1.d0)) :: val 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 


          istat = blas_sparse_const_success
          CALL blas_zuscr_insert_entry&
           &(A,val,i,j,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE suscr_insert_entries&
         &(A,nnz,val,indx,jndx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          REAL(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: jndx (:)


          istat = blas_sparse_const_success
          CALL blas_suscr_insert_entries&
           &(A,nnz,val,indx,jndx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE duscr_insert_entries&
         &(A,nnz,val,indx,jndx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          REAL(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: jndx (:)


          istat = blas_sparse_const_success
          CALL blas_duscr_insert_entries&
           &(A,nnz,val,indx,jndx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cuscr_insert_entries&
         &(A,nnz,val,indx,jndx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          COMPLEX(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: jndx (:)


          istat = blas_sparse_const_success
          CALL blas_cuscr_insert_entries&
           &(A,nnz,val,indx,jndx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_entries_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zuscr_insert_entries&
         &(A,nnz,val,indx,jndx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          COMPLEX(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: jndx (:)


          istat = blas_sparse_const_success
          CALL blas_zuscr_insert_entries&
           &(A,nnz,val,indx,jndx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE suscr_insert_col&
         &(A,j,nnz,val,indx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          REAL(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)


          istat = blas_sparse_const_success
          CALL blas_suscr_insert_col&
           &(A,j,nnz,val,indx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE duscr_insert_col&
         &(A,j,nnz,val,indx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          REAL(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)


          istat = blas_sparse_const_success
          CALL blas_duscr_insert_col&
           &(A,j,nnz,val,indx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cuscr_insert_col&
         &(A,j,nnz,val,indx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          COMPLEX(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)


          istat = blas_sparse_const_success
          CALL blas_cuscr_insert_col&
           &(A,j,nnz,val,indx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_col_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zuscr_insert_col&
         &(A,j,nnz,val,indx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          COMPLEX(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)


          istat = blas_sparse_const_success
          CALL blas_zuscr_insert_col&
           &(A,j,nnz,val,indx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE suscr_insert_row&
         &(A,i,nnz,val,indx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          REAL(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)


          istat = blas_sparse_const_success
          CALL blas_suscr_insert_row&
           &(A,i,nnz,val,indx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE duscr_insert_row&
         &(A,i,nnz,val,indx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          REAL(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)


          istat = blas_sparse_const_success
          CALL blas_duscr_insert_row&
           &(A,i,nnz,val,indx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cuscr_insert_row&
         &(A,i,nnz,val,indx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          COMPLEX(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)


          istat = blas_sparse_const_success
          CALL blas_cuscr_insert_row&
           &(A,i,nnz,val,indx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_row_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zuscr_insert_row&
         &(A,i,nnz,val,indx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nnz 
          COMPLEX(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)


          istat = blas_sparse_const_success
          CALL blas_zuscr_insert_row&
           &(A,i,nnz,val,indx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE suscr_insert_clique&
         &(A,k,l,val,row_stride,col_stride,indx,jndx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: k 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: l 
          REAL(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: row_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: col_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: jndx (:)


          istat = blas_sparse_const_success
          CALL blas_suscr_insert_clique&
           &(A,k,l,val,row_stride,col_stride,indx,jndx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE duscr_insert_clique&
         &(A,k,l,val,row_stride,col_stride,indx,jndx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: k 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: l 
          REAL(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: row_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: col_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: jndx (:)


          istat = blas_sparse_const_success
          CALL blas_duscr_insert_clique&
           &(A,k,l,val,row_stride,col_stride,indx,jndx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cuscr_insert_clique&
         &(A,k,l,val,row_stride,col_stride,indx,jndx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: k 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: l 
          COMPLEX(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: row_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: col_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: jndx (:)


          istat = blas_sparse_const_success
          CALL blas_cuscr_insert_clique&
           &(A,k,l,val,row_stride,col_stride,indx,jndx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_clique_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zuscr_insert_clique&
         &(A,k,l,val,row_stride,col_stride,indx,jndx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: k 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: l 
          COMPLEX(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: row_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: col_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: indx (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: jndx (:)


          istat = blas_sparse_const_success
          CALL blas_zuscr_insert_clique&
           &(A,k,l,val,row_stride,col_stride,indx,jndx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE suscr_insert_block&
         &(A,val,row_stride,col_stride,i,j,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          REAL(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: row_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: col_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 


          istat = blas_sparse_const_success
          CALL blas_suscr_insert_block&
           &(A,val,row_stride,col_stride,i,j,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE duscr_insert_block&
         &(A,val,row_stride,col_stride,i,j,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          REAL(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: row_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: col_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 


          istat = blas_sparse_const_success
          CALL blas_duscr_insert_block&
           &(A,val,row_stride,col_stride,i,j,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cuscr_insert_block&
         &(A,val,row_stride,col_stride,i,j,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          COMPLEX(KIND(1.e0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: row_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: col_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 


          istat = blas_sparse_const_success
          CALL blas_cuscr_insert_block&
           &(A,val,row_stride,col_stride,i,j,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_cr_insert_block_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zuscr_insert_block&
         &(A,val,row_stride,col_stride,i,j,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: A 
          COMPLEX(KIND(1.d0)) :: val (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: row_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: col_stride 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: i 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: j 


          istat = blas_sparse_const_success
          CALL blas_zuscr_insert_block&
           &(A,val,row_stride,col_stride,i,j,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_mv_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE susmv&
         &(transA,alpha,A,x,incx,y,incy,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: transA 
          REAL(KIND(1.e0)) :: alpha 
          INTEGER :: A 
          REAL(KIND(1.e0)) :: x (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incx 
          REAL(KIND(1.e0)) :: y (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incy 


          istat = blas_sparse_const_success
          CALL blas_susmv&
           &(transA,alpha,A,x,incx,y,incy,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_mv_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE dusmv&
         &(transA,alpha,A,x,incx,y,incy,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: transA 
          REAL(KIND(1.d0)) :: alpha 
          INTEGER :: A 
          REAL(KIND(1.d0)) :: x (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incx 
          REAL(KIND(1.d0)) :: y (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incy 


          istat = blas_sparse_const_success
          CALL blas_dusmv&
           &(transA,alpha,A,x,incx,y,incy,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_mv_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cusmv&
         &(transA,alpha,A,x,incx,y,incy,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: transA 
          COMPLEX(KIND(1.e0)) :: alpha 
          INTEGER :: A 
          COMPLEX(KIND(1.e0)) :: x (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incx 
          COMPLEX(KIND(1.e0)) :: y (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incy 


          istat = blas_sparse_const_success
          CALL blas_cusmv&
           &(transA,alpha,A,x,incx,y,incy,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_mv_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zusmv&
         &(transA,alpha,A,x,incx,y,incy,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: transA 
          COMPLEX(KIND(1.d0)) :: alpha 
          INTEGER :: A 
          COMPLEX(KIND(1.d0)) :: x (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incx 
          COMPLEX(KIND(1.d0)) :: y (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incy 


          istat = blas_sparse_const_success
          CALL blas_zusmv&
           &(transA,alpha,A,x,incx,y,incy,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_sv_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE sussv&
         &(transT,alpha,T,x,incx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: transT 
          REAL(KIND(1.e0)) :: alpha 
          INTEGER :: T 
          REAL(KIND(1.e0)) :: x (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incx 


          istat = blas_sparse_const_success
          CALL blas_sussv&
           &(transT,alpha,T,x,incx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_sv_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE dussv&
         &(transT,alpha,T,x,incx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: transT 
          REAL(KIND(1.d0)) :: alpha 
          INTEGER :: T 
          REAL(KIND(1.d0)) :: x (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incx 


          istat = blas_sparse_const_success
          CALL blas_dussv&
           &(transT,alpha,T,x,incx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_sv_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cussv&
         &(transT,alpha,T,x,incx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: transT 
          COMPLEX(KIND(1.e0)) :: alpha 
          INTEGER :: T 
          COMPLEX(KIND(1.e0)) :: x (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incx 


          istat = blas_sparse_const_success
          CALL blas_cussv&
           &(transT,alpha,T,x,incx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_sv_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zussv&
         &(transT,alpha,T,x,incx,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: transT 
          COMPLEX(KIND(1.d0)) :: alpha 
          INTEGER :: T 
          COMPLEX(KIND(1.d0)) :: x (:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: incx 


          istat = blas_sparse_const_success
          CALL blas_zussv&
           &(transT,alpha,T,x,incx,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_mm_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE susmm&
         &(order,transA,nrhs,alpha,A,b,ldb,c,ldc,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: order 
          INTEGER :: transA 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nrhs 
          REAL(KIND(1.e0)) :: alpha 
          INTEGER :: A 
          REAL(KIND(1.e0)) :: b (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldb 
          REAL(KIND(1.e0)) :: c (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldc 


          istat = blas_sparse_const_success
          CALL blas_susmm&
           &(order,transA,nrhs,alpha,A,b,ldb,c,ldc,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_mm_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE dusmm&
         &(order,transA,nrhs,alpha,A,b,ldb,c,ldc,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: order 
          INTEGER :: transA 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nrhs 
          REAL(KIND(1.d0)) :: alpha 
          INTEGER :: A 
          REAL(KIND(1.d0)) :: b (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldb 
          REAL(KIND(1.d0)) :: c (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldc 


          istat = blas_sparse_const_success
          CALL blas_dusmm&
           &(order,transA,nrhs,alpha,A,b,ldb,c,ldc,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_mm_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cusmm&
         &(order,transA,nrhs,alpha,A,b,ldb,c,ldc,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: order 
          INTEGER :: transA 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nrhs 
          COMPLEX(KIND(1.e0)) :: alpha 
          INTEGER :: A 
          COMPLEX(KIND(1.e0)) :: b (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldb 
          COMPLEX(KIND(1.e0)) :: c (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldc 


          istat = blas_sparse_const_success
          CALL blas_cusmm&
           &(order,transA,nrhs,alpha,A,b,ldb,c,ldc,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_mm_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zusmm&
         &(order,transA,nrhs,alpha,A,b,ldb,c,ldc,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: order 
          INTEGER :: transA 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nrhs 
          COMPLEX(KIND(1.d0)) :: alpha 
          INTEGER :: A 
          COMPLEX(KIND(1.d0)) :: b (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldb 
          COMPLEX(KIND(1.d0)) :: c (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldc 


          istat = blas_sparse_const_success
          CALL blas_zusmm&
           &(order,transA,nrhs,alpha,A,b,ldb,c,ldc,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
        
        !> \rsb_spblasl2_sm_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE sussm&
         &(order,transT,nrhs,alpha,T,b,ldb,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: order 
          INTEGER :: transT 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nrhs 
          REAL(KIND(1.e0)) :: alpha 
          INTEGER :: T 
          REAL(KIND(1.e0)) :: b (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldb 


          istat = blas_sparse_const_success
          CALL blas_sussm&
           &(order,transT,nrhs,alpha,T,b,ldb,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_sm_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE dussm&
         &(order,transT,nrhs,alpha,T,b,ldb,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: order 
          INTEGER :: transT 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nrhs 
          REAL(KIND(1.d0)) :: alpha 
          INTEGER :: T 
          REAL(KIND(1.d0)) :: b (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldb 


          istat = blas_sparse_const_success
          CALL blas_dussm&
           &(order,transT,nrhs,alpha,T,b,ldb,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_sm_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE cussm&
         &(order,transT,nrhs,alpha,T,b,ldb,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: order 
          INTEGER :: transT 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nrhs 
          COMPLEX(KIND(1.e0)) :: alpha 
          INTEGER :: T 
          COMPLEX(KIND(1.e0)) :: b (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldb 


          istat = blas_sparse_const_success
          CALL blas_cussm&
           &(order,transT,nrhs,alpha,T,b,ldb,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        !> \rsb_spblasl2_sm_msg\rsb_spblas_return_msg
        !> \rsb_spblas_f_istat_msg
        !! 
        
        SUBROUTINE zussm&
         &(order,transT,nrhs,alpha,T,b,ldb,istat)
          IMPLICIT NONE
          INTEGER(KIND=RSB_BLAS_IST_KIND), INTENT(OUT) ::istat
          INTEGER :: order 
          INTEGER :: transT 
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: nrhs 
          COMPLEX(KIND(1.d0)) :: alpha 
          INTEGER :: T 
          COMPLEX(KIND(1.d0)) :: b (:,:)
          INTEGER(KIND=RSB_BLAS_IDX_KIND) :: ldb 


          istat = blas_sparse_const_success
          CALL blas_zussm&
           &(order,transT,nrhs,alpha,T,b,ldb,istat)

          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
        END SUBROUTINE
        
        
        
      END MODULE blas_sparse
