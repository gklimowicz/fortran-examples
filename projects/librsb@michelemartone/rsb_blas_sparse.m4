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
dnl
include(`rsb_fortran_macros.m4')dnl
dnl
define(`RSB_M4_BLAS_SPARSE_INTERFACE_BEGIN',`ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`INTERFACE',`')')dnl
define(`RSB_M4_BLAS_SPARSE_INTERFACE_END',`ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`END INTERFACE',`')')dnl
dnl
dnl        @author: Michele Martone
dnl
dnl        This macro generates a Sparse BLAS fortran module for librsb.
dnl !! @cond INNERDOC 
dnl ! author: Michele Martone
!
!> @file
!! @brief Implementation of the Fortran Sparse BLAS interface to \librsb (see \ref rsb_doc_sparse_blas).
!!
ifelse(`0',`1',`
!! Supported types are: foreach(`mtype',RSB_M4_SBLAS_MATRIX_SUPPORTED_TYPES,` RSB_M4_C2F_TYPE(mtype)').
!! Supported operations are: foreach(`pmop',RSB_M4_SBLAS_INTERFACE_OPS,` RSB_M4_SBLAS_INTERFACE_IDENTIFIER(pmop)').
')dnl
!!
ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`
#define RSB_HAVE_RSB_KERNELS 1
dnl foreach(`type',RSB_M4_SBLAS_MATRIX_SUPPORTED_TYPES,`dnl
dnl dnl
dnl `#define' RSB_M4_HAVE_TYPE_PREPROCESSOR_SYMBOL(type) 1 /*!< Type type is supported.*/
dnl dnl
dnl ')dnl
')dnl

ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`dnl
      MODULE blas_sparse
        !> A Sparse BLAS interface for RSB
        IMPLICIT NONE
PUBLIC
')dnl

dnl ifelse(RSB_M4_SPBLAS_MATRIX_SUPPORTED_TYPES_LIST_LENGTH,0,`',`dnl
        foreach(`pmop',RSB_M4_SBLAS_GENERIC_OPS,`
        !> RSB_M4_SBLAS_SUBROUTINE_HELP_COMMENT(pmop,`*')
        !> RSB_M4_SBLAS_SUBROUTINE_EXTRA_FORTRAN_HELP_COMMENT(pmop)
dnl        !> RSB_M4_SPBLAS_HELP_INFO(pmop)
        !! 
dnl         MODULE PROCEDURE RSB_M4_INTERFACE_LIST(RSB_M4_COMMA_LIST((RSB_M4_CHOPTRAILINGSPACES(foreach(`mtype',RSB_M4_SBLAS_MATRIX_SUPPORTED_TYPES,`RSB_M4_SBLAS_SUBROUTINE_IDENTIFIER(pmop,mtype) ')))))dnl
        INTERFACE RSB_M4_SBLAS_INTERFACE_IDENTIFIER(pmop)
        ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`PROCEDURE',`MODULE PROCEDURE') RSB_M4_INTERFACE_LIST(RSB_M4_COMMA_LIST((RSB_M4_CHOPTRAILINGSPACES(foreach(`mtype',RSB_M4_SPBLAS_MATRIX_ALL_TYPES,`RSB_M4_SBLAS_SUBROUTINE_IDENTIFIER(pmop,mtype) ')))))dnl
        END INTERFACE
        ')
dnl ')dnl
dnl
ifelse(RSB_M4_want_old_fortran_float_types,`1',`dnl
dnl        ....
')dnl
dnl
        INTEGER, PARAMETER :: blas_sparse_const_success=0
        INTEGER, PARAMETER :: blas_sparse_const_failure=-1 ! value returned by this interface on failure
        INTEGER, PARAMETER :: blas_sparse_const_not_available=-9999 ! value returned by this interface when deactivated
ifelse(`0',`1',`
        INTEGER, PARAMETER :: blas_lower_hermitian=239 ! # FIXME
        INTEGER, PARAMETER :: blas_lower_symmetric=237 ! # FIXME
        !
        INTEGER, PARAMETER :: blas_unit_diag=132 ! # FIXME
        !
        INTEGER, PARAMETER :: blas_lower_triangular=235 ! # FIXME
        INTEGER, PARAMETER :: blas_upper_triangular=236 ! # FIXME
        INTEGER, PARAMETER :: blas_no_trans=111 ! # FIXME
        INTEGER, PARAMETER :: blas_trans=112 ! # FIXME
        INTEGER, PARAMETER :: blas_conj_trans=113 ! # FIXME
        INTEGER, PARAMETER :: blas_rsb_autotuning_on = 666
        INTEGER, PARAMETER :: blas_rsb_autotuning_off = 999
',`dnl
include(`blas_sparse/blas_enum.F90')dnl
')dnl

dnl ifelse(RSB_M4_LONG_IDX,`0',`dnl
dnl         INTEGER,PARAMETER :: RSB_BLAS_IDX_KIND=4
dnl         INTEGER,PARAMETER :: RSB_BLAS_IST_KIND=8 ! for istat
dnl ',`dnl
dnl         INTEGER,PARAMETER :: RSB_BLAS_IDX_KIND=8
dnl         INTEGER,PARAMETER :: RSB_BLAS_IST_KIND=4 ! for istat
dnl ')dnl
dnl dnl
#ifdef RSB_WANT_LONG_IDX_TYPE
        INTEGER,PARAMETER :: RSB_BLAS_IDX_KIND=8
        INTEGER,PARAMETER :: RSB_BLAS_IST_KIND=8 ! for istat
#define C_RSB_INT_KND_ C_INT64_T
#else
        INTEGER,PARAMETER :: RSB_BLAS_IDX_KIND=4
        INTEGER,PARAMETER :: RSB_BLAS_IST_KIND=4 ! for istat
#define C_RSB_INT_KND_ C_INT
#endif
dnl
ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`
        INTERFACE
          TYPE(C_PTR) FUNCTION &
          &rsb_blas_get_mtx&
          &(A)&
          &BIND(c,NAME = "rsb_blas_get_mtx")
          USE ISO_C_BINDING
          INTEGER(C_INT), VALUE  :: A
          END FUNCTION rsb_blas_get_mtx
        END INTERFACE
')dnl
dnl

ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`CONTAINS')

dnl         SUBROUTINE RSB_M4_SBLAS_INTERFACE_RADIX`_'init(istat)
dnl           IMPLICIT NONE
dnl           INTEGER::istat
dnl           istat=blas_sparse_const_success
dnl #ifdef RSB_HAVE_RSB_KERNELS
dnl           CALL RSB_M4_SBLAS2VBR_SUBROUTINE_RADIX`'init(istat)
dnl           IF(istat.NE.blas_sparse_const_success)istat=blas_sparse_const_failure
dnl #else
dnl           istat=blas_sparse_const_not_available
dnl #endif
dnl         END SUBROUTINE
dnl 
         !> RSB_M4_SPBLAS_HELP_INFO(`ds')
         !> RSB_M4_SBLAS_SUBROUTINE_EXTRA_FORTRAN_HELP_COMMENT()
         !! 
         RSB_M4_BLAS_SPARSE_INTERFACE_BEGIN
         SUBROUTINE RSB_M4_SBLAS_INTERFACE_RADIX`'ds(A,istat)
           IMPLICIT NONE
           INTEGER,INTENT(IN)::A
           INTEGER(KIND=RSB_BLAS_IST_KIND)::istat
ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`
           istat=blas_sparse_const_success
`#if defined(RSB_HAVE_RSB_KERNELS)'
           CALL RSB_M4_SBLAS2VBR_SUBROUTINE_IDENTIFIER(`ds',`',`f90')`'(A,istat)
           IF(istat.NE.blas_sparse_const_success)&
            &istat=blas_sparse_const_failure
#else /* RSB_HAVE_RSB_KERNELS */
           istat=blas_sparse_const_not_available
#endif /* RSB_HAVE_RSB_KERNELS */
')dnl
         END SUBROUTINE
         RSB_M4_BLAS_SPARSE_INTERFACE_END

dnl           !> RSB_M4_SBLAS_SUBROUTINE_HELP_COMMENT(`cr_end',`*')
         !> RSB_M4_SPBLAS_HELP_INFO(`cr_end')
         !> RSB_M4_SBLAS_SUBROUTINE_EXTRA_FORTRAN_HELP_COMMENT()
         !! 
         RSB_M4_BLAS_SPARSE_INTERFACE_BEGIN
         SUBROUTINE RSB_M4_SBLAS_SUBROUTINE_IDENTIFIER(`cr_end',`')`'RSB_M4_SBLAS_SUBROUTINE_ARGS(`cr_end',`',`f90')
           IMPLICIT NONE
RSB_M4_SBLAS_SUBROUTINE_INFO_DECLARATION(istat)dnl
           INTEGER,INTENT(IN)::A
ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`
           istat=blas_sparse_const_success
`#if defined(RSB_HAVE_RSB_KERNELS)'
           CALL RSB_M4_SBLAS2VBR_SUBROUTINE_IDENTIFIER(`cr_end',`',`f90')`'dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_SB_INTERFACE(`(RSB_M4_SBLAS_SUBROUTINE_ARGS(`cr_end',`',`f90'))')
           IF(istat.NE.blas_sparse_const_success)&
            &istat=blas_sparse_const_failure
#else /* RSB_HAVE_RSB_KERNELS */
           istat=blas_sparse_const_not_available
#endif /* RSB_HAVE_RSB_KERNELS */
')dnl
         END SUBROUTINE
         RSB_M4_BLAS_SPARSE_INTERFACE_END

dnl           !> RSB_M4_SBLAS_SUBROUTINE_HELP_COMMENT(`gp',`*')
         !> RSB_M4_SPBLAS_HELP_INFO(`gp')
         !> RSB_M4_SBLAS_SUBROUTINE_EXTRA_FORTRAN_HELP_COMMENT()
         !! 
         RSB_M4_BLAS_SPARSE_INTERFACE_BEGIN
         SUBROUTINE RSB_M4_SBLAS_SUBROUTINE_IDENTIFIER(`gp',`')`'RSB_M4_SBLAS_SUBROUTINE_ARGS(`gp',`',`f90')
           IMPLICIT NONE
RSB_M4_SBLAS_SUBROUTINE_INFO_DECLARATION(istat)dnl
           INTEGER,INTENT(IN)::A
           INTEGER,INTENT(IN)::pname
ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`
           istat=blas_sparse_const_success
`#if defined(RSB_HAVE_RSB_KERNELS)'
           CALL RSB_M4_SBLAS2VBR_SUBROUTINE_IDENTIFIER(`gp',`',`f90')`'dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_SB_INTERFACE(`(RSB_M4_SBLAS_SUBROUTINE_ARGS(`gp',`',`f90'))')
           !istat does not have the meaning of an error value, here
           !IF(istat.NE.blas_sparse_const_success)istat=blas_sparse_const_failure
#else /* RSB_HAVE_RSB_KERNELS */
           istat=blas_sparse_const_not_available
#endif /* RSB_HAVE_RSB_KERNELS */
')dnl
         END SUBROUTINE
         RSB_M4_BLAS_SPARSE_INTERFACE_END

         !> RSB_M4_SPBLAS_HELP_INFO(`sp')
dnl           !> RSB_M4_SBLAS_SUBROUTINE_HELP_COMMENT(`sp',`*')
         !> RSB_M4_SBLAS_SUBROUTINE_EXTRA_FORTRAN_HELP_COMMENT()
         !! 
         RSB_M4_BLAS_SPARSE_INTERFACE_BEGIN
         SUBROUTINE RSB_M4_SBLAS_SUBROUTINE_IDENTIFIER(`sp',`')`'RSB_M4_SBLAS_SUBROUTINE_ARGS(`sp',`',`f90')
           IMPLICIT NONE
RSB_M4_SBLAS_SUBROUTINE_INFO_DECLARATION(istat)dnl
           INTEGER,INTENT(IN)::A
           INTEGER,INTENT(IN)::pname
ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`
           istat=blas_sparse_const_success
`#if defined(RSB_HAVE_RSB_KERNELS)'
           CALL RSB_M4_SBLAS2VBR_SUBROUTINE_IDENTIFIER(`sp',`',`f90')`'dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_SB_INTERFACE(`(RSB_M4_SBLAS_SUBROUTINE_ARGS(`sp',`',`f90'))')
           IF(istat.NE.blas_sparse_const_success)&
            &istat=blas_sparse_const_failure
#else /* RSB_HAVE_RSB_KERNELS */
           istat=blas_sparse_const_not_available
#endif /* RSB_HAVE_RSB_KERNELS */
')
         END SUBROUTINE
         RSB_M4_BLAS_SPARSE_INTERFACE_END

        foreach(`pmop',RSB_M4_SBLAS_INTERFACE_OPS,`
dnl        RSB_M4_SBLAS_MATRIX_SUPPORTED_TYPES
        foreach(`mtype',RSB_M4_SPBLAS_MATRIX_ALL_TYPES,`
dnl          !> RSB_M4_SBLAS_SUBROUTINE_HELP_COMMENT(pmop,mtype)
        !> RSB_M4_SPBLAS_HELP_INFO(pmop)
        !> RSB_M4_SBLAS_SUBROUTINE_EXTRA_FORTRAN_HELP_COMMENT(pmop)
        !! 
        RSB_M4_BLAS_SPARSE_INTERFACE_BEGIN
        SUBROUTINE RSB_M4_SBLAS_SUBROUTINE_IDENTIFIER(pmop,mtype)`&
         &'RSB_M4_SBLAS_SUBROUTINE_ARGS(pmop,mtype,`f90')dnl
          IMPLICIT NONE
RSB_M4_SBLAS_SUBROUTINE_INFO_DECLARATION(istat)dnl
RSB_M4_SBLAS_SUBROUTINE_ARGS_DECLARATION(pmop,mtype)dnl
dnl
dnl ifelse(RSB_M4_IS_L3_MOP(pmop),`0',`
dnl RSB_M4_SBLAS_SUBROUTINE_INFO_DECLARATION(istat)dnl
dnl RSB_M4_SBLAS_SUBROUTINE_ARGS_DECLARATION(pmop,mtype)dnl
dnl ')
dnl ifelse(RSB_M4_IS_L3_MOP(pmop),`1',`
dnl ! FIXME: need 2D arrays decls here (pmop).
dnl RSB_M4_SBLAS_SUBROUTINE_INFO_DECLARATION(istat)dnl
dnl RSB_M4_SBLAS_SUBROUTINE_ARGS_DECLARATION(pmop`_2d',mtype)dnl
dnl ')
dnl
ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`
ifelse(`0',`1',`
dnl `#if ( defined(RSB_HAVE_RSB_KERNELS)' && RSB_M4_HAVE_TYPE(mtype))
`#if ( defined(RSB_HAVE_RSB_KERNELS)' && defined(RSB_M4_HAVE_TYPE_PREPROCESSOR_SYMBOL(mtype)))
          istat = blas_sparse_const_success
          CALL RSB_M4_SBLAS2VBR_SUBROUTINE_IDENTIFIER(pmop,mtype,`f90')`'dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_SB_INTERFACE(`(RSB_M4_SBLAS_SUBROUTINE_ARGS(pmop,mtype,`f90'))')
          IF(istat.NE.blas_sparse_const_success)&
           &istat=blas_sparse_const_failure
#else  /* defined(RSB_HAVE_RSB_KERNELS) && RSB_M4_HAVE_TYPE('mtype`) */
          istat=blas_sparse_const_not_available
#endif /* defined(RSB_HAVE_RSB_KERNELS) && RSB_M4_HAVE_TYPE('mtype`) */
          istat=blas_sparse_const_success
',`
          istat = blas_sparse_const_success
          CALL RSB_M4_SBLAS2VBR_SUBROUTINE_IDENTIFIER(pmop,mtype,`f90')`&
           &'dnl
RSB_M4_ARGS_TO_ACTUAL_ARGS_FOR_SB_INTERFACE(`(RSB_M4_SBLAS_SUBROUTINE_ARGS(pmop,mtype,`f90'))')
          IF(istat.NE.blas_sparse_const_success)&
           &istat = blas_sparse_const_failure
')dnl
')dnl
        END SUBROUTINE
        RSB_M4_BLAS_SPARSE_INTERFACE_END
        ')
        ')
dnl
ifelse(RSB_M4_WANT_BLAS_SPARSE_INTERFACE,`1',`',`dnl
      END MODULE blas_sparse
')dnl
dnl
dnl !! @endcond
dnl
