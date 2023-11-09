! 
! Copyright (C) 2008-2022 Michele Martone
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

!> @file.
!! @brief RSB.F90-based usage:
!!        rsb_mtx_alloc_from_coo_const(),
!!        rsb_tune_spmm(),
!!        rsb_file_mtx_save(),
!!        rsb_spmv(),
!!        ...
!! \include fortran_rsb_fi.F90

      SUBROUTINE rsb_mod_example1(res)
! Example using an unsymmetric matrix specified as COO, and autotuning, built at once.
      USE rsb
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER ::res
      INTEGER,TARGET :: istat = 0, i
      INTEGER :: transt = RSB_TRANSPOSITION_N !
      INTEGER(KIND=RSB_IDX_KIND) :: incx = 1, incy = 1
      REAL(KIND=8),TARGET :: alpha = 3, beta = 1
! 1 1
! 1 1
      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz = 4
      INTEGER(KIND=RSB_IDX_KIND) :: nr = 2
      INTEGER(KIND=RSB_IDX_KIND) :: nc = 2
      INTEGER(KIND=RSB_IDX_KIND) :: nrhs = 1
      INTEGER :: order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER ! rhs layout
      INTEGER :: flags = RSB_FLAG_NOFLAGS 
      INTEGER(KIND=RSB_IDX_KIND),TARGET :: IA(4) = (/0, 1, 1,0/)
      INTEGER(KIND=RSB_IDX_KIND),TARGET :: JA(4) = (/0, 0, 1,1/)
      REAL(KIND=8),TARGET :: VA(4) = (/1,1,1,1/)
      REAL(KIND=8),TARGET :: x(2) = (/1, 1/)! reference x 
      REAL(KIND=8),TARGET :: cy(2) = (/9, 9/)! reference cy after 
      REAL(KIND=8),TARGET :: y(2) = (/3, 3/)! y will be overwritten
      TYPE(C_PTR),TARGET :: mtxAp = C_NULL_PTR ! matrix pointer
      REAL(KIND=8) :: tmax = 2.0 ! tuning max time
      INTEGER :: titmax = 2 ! tuning max iterations
      INTEGER,TARGET :: ont = 0     ! optimal number of threads
      !TYPE(C_PTR),PARAMETER :: EO = RSB_NULL_EXIT_OPTIONS
      !TYPE(C_PTR),PARAMETER :: IO = RSB_NULL_INIT_OPTIONS
      ! Note: using C_NULL_PTR instead of the previous lines because of http://gcc.gnu.org/bugzilla/show_bug.cgi?id=59411
      TYPE(C_PTR),PARAMETER :: EO = C_NULL_PTR
      TYPE(C_PTR),PARAMETER :: IO = C_NULL_PTR
      INTEGER,TARGET :: errval

      res = 0

      errval = rsb_lib_init(IO)
      IF (errval.NE.RSB_ERR_NO_ERROR) GOTO 9997

      mtxAp = rsb_mtx_alloc_from_coo_const(C_LOC(VA),C_LOC(IA),C_LOC(JA)&
       &,nnz,&
       & RSB_NUMERICAL_TYPE_DOUBLE,nr,nc,1,1,flags,C_LOC(istat))

      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997

      istat = rsb_file_mtx_save(mtxAp,C_NULL_CHAR)

      ! Structure autotuning:
      istat = rsb_tune_spmm(C_LOC(mtxAp),C_NULL_PTR,C_NULL_PTR,titmax,&
       & tmax,&
       & transt,C_LOC(alpha),C_NULL_PTR,nrhs,order,C_LOC(x),nr,&
       & C_LOC(beta),C_LOC(y),nc)

      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997

      ! Thread count autotuning:
      istat = rsb_tune_spmm(C_NULL_PTR,C_NULL_PTR,C_LOC(ont),titmax,&
       & tmax,&
       & transt,C_LOC(alpha),mtxAp,nrhs,order,C_LOC(x),nr,C_LOC(beta),&
       & C_LOC(y),nc)
      PRINT *, "Optimal number of threads:", ont

      y(:) = (/3, 3/)! restoring reference y (rsb_tune_spmm has overwritten it)
      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997
      
      istat = rsb_file_mtx_save(mtxAp,C_NULL_CHAR)
      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997

      istat = rsb_spmv(transt,C_LOC(alpha),mtxAp,C_LOC(x),incx,&
       & C_LOC(beta),C_LOC(y),incy)
      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997
      DO i = 1, 2
            IF (y(i).NE.cy(i)) PRINT *, "type=d dims=2x2 sym=g diag=g &
      &blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF (y(i).NE.cy(i)) GOTO 9997
      END DO
      PRINT*,"type=d dims=2x2 sym=g diag=g blocks=1x1 usmv alpha= 3&
       & beta= 1 incx=1 incy=1 trans=n is ok"

      IF ( res .NE.RSB_ERR_NO_ERROR) GOTO 9997
      GOTO 9998
9997      res = -1
9998      CONTINUE
      mtxAp = rsb_mtx_free(mtxAp)
      IF (istat.NE.RSB_ERR_NO_ERROR) res = -1 
! 9999      CONTINUE
      errval=rsb_lib_exit(EO)                 ! librsb finalization
      IF (errval.NE.RSB_ERR_NO_ERROR)&
       & STOP "error calling rsb_lib_exit"
      PRINT *, "rsb module fortran test is ok"

      istat = rsb_perror(C_NULL_PTR,istat)
      end SUBROUTINE rsb_mod_example1

      SUBROUTINE rsb_mod_example2(res)
! Example using an unsymmetric matrix specified as COO, and autotuning, built piecewise.
      USE rsb
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER,TARGET :: errval
      INTEGER :: res
      INTEGER :: transt = RSB_TRANSPOSITION_N  ! no transposition
      INTEGER(KIND=RSB_IDX_KIND) :: incX = 1, incB = 1        ! X, B vectors increment
      REAL(KIND=8),TARGET :: alpha = 3,beta = 1
      INTEGER(KIND=RSB_IDX_KIND) :: nnzA = 4, nrA = 3, ncA = 3     ! nonzeroes, rows, columns of matrix A
      INTEGER(KIND=RSB_IDX_KIND),TARGET :: IA(4) = (/1, 2, 3, 3/)  ! row    indices
      INTEGER(KIND=RSB_IDX_KIND),TARGET :: JA(4) = (/1, 2, 1, 3/)  ! column indices
      INTEGER(C_SIGNED_CHAR) :: typecode = RSB_NUMERICAL_TYPE_DOUBLE
      INTEGER :: flags =RSB_FLAG_DEFAULT_MATRIX_FLAGS+RSB_FLAG_SYMMETRIC
      REAL(KIND=8),TARGET :: VA(4) = (/11.0, 22.0, 13.0, 33.0/) ! coefficients
      REAL(KIND=8),TARGET :: X(3) = (/   0,    0,    0/)
      REAL(KIND=8),TARGET :: B(3) = (/-1.0, -2.0, -2.0/)
      TYPE(C_PTR),TARGET  :: mtxAp = C_NULL_PTR
      TYPE(C_PTR)  :: mtxApp = C_NULL_PTR
      REAL(KIND=8),TARGET :: ETIME = 0.0
      !TYPE(C_PTR),PARAMETER :: EO = RSB_NULL_EXIT_OPTIONS
      !TYPE(C_PTR),PARAMETER :: IO = RSB_NULL_INIT_OPTIONS
      ! Note: using C_NULL_PTR instead of the previous lines because of http://gcc.gnu.org/bugzilla/show_bug.cgi?id=59411
      TYPE(C_PTR),PARAMETER :: EO = C_NULL_PTR
      TYPE(C_PTR),PARAMETER :: IO = C_NULL_PTR

      errval = rsb_lib_init(IO)                ! librsb initialization
      IF (errval.NE.RSB_ERR_NO_ERROR) &
       & STOP "error calling rsb_lib_init"
#if defined(__GNUC__) && (__GNUC__ == 4) && (__GNUC_MINOR__ < 5)
#define RSB_SKIP_BECAUSE_OLD_COMPILER 1
#endif
#ifndef RSB_SKIP_BECAUSE_OLD_COMPILER
      mtxAp = rsb_mtx_alloc_from_coo_begin(nnzA,typecode,nrA,ncA,flags,&
       & C_LOC(errval)) ! begin matrix creation
      errval = rsb_mtx_set_vals(mtxAp,&
       & C_LOC(VA),C_LOC(IA),C_LOC(JA),nnzA,flags) ! insert some nonzeroes
      mtxApp = C_LOC(mtxAp) ! Old compilers like e.g.: Gfortran 4.4.7 will NOT compile this.
      IF (errval.NE.RSB_ERR_NO_ERROR) &
       & STOP "error calling rsb_mtx_set_vals"
      errval = rsb_mtx_alloc_from_coo_end(mtxApp)                   ! end matrix creation
      IF (errval.NE.RSB_ERR_NO_ERROR) &
       & STOP "error calling rsb_mtx_alloc_from_coo_end"
      errval = rsb_spmv(transt,C_LOC(alpha),mtxAp,C_LOC(X),&
       & incX,C_LOC(beta),C_LOC(B),incB) ! X := X + (3) * A * B 
      IF (errval.NE.RSB_ERR_NO_ERROR)&
       & STOP "error calling rsb_spmv"
      mtxAp = rsb_mtx_free(mtxAp)                                 ! destroy matrix

      ! The following is optional and depends on configure options, so it is allowed to fail
      errval = rsb_lib_get_opt(RSB_IO_WANT_LIBRSB_ETIME,C_LOC(ETIME))
      IF (errval.EQ.RSB_ERR_NO_ERROR)&
       & PRINT*,"Time spent in librsb is:",ETIME
      ! IF (errval.NE.0)STOP "error calling rsb_lib_get_opt" 
      errval = RSB_ERR_NO_ERROR

      IF (errval.NE.RSB_ERR_NO_ERROR) &
       & STOP "error calling rsb_mtx_free"
#else
      PRINT*,"You have an old Fortran compiler not supporting C_LOC."
      PRINT*,"Skipping a part of the test"
#endif
      errval=rsb_lib_exit(EO)                 ! librsb finalization
      IF (errval.NE.RSB_ERR_NO_ERROR)&
       & STOP "error calling rsb_lib_exit"
      PRINT *, "rsb module fortran test is ok"
      res = errval
      end SUBROUTINE rsb_mod_example2

      SUBROUTINE rsb_mod_example3(res)
! Example using a symmetric matrix specified as CSR, and autotuning, built at once.
      USE rsb
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER ::res
      INTEGER,TARGET :: istat = 0, i
      INTEGER :: transt = RSB_TRANSPOSITION_N !
      INTEGER(KIND=RSB_IDX_KIND) :: incx = 1, incy = 1
      REAL(KIND=8),TARGET :: alpha = 4, beta = 1
! 11 21 
! 21 22 
      ! declaration of VA,IP,JA 
      INTEGER(KIND=RSB_IDX_KIND), PARAMETER :: nnz = 3
      INTEGER(KIND=RSB_IDX_KIND), PARAMETER :: nr = 2
      INTEGER(KIND=RSB_IDX_KIND), PARAMETER :: nc = 2
      INTEGER(KIND=RSB_IDX_KIND), PARAMETER :: nrhs = 1
      INTEGER :: order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER ! rhs layout
      INTEGER :: flags = RSB_FLAG_NOFLAGS + RSB_FLAG_SYMMETRIC + &
       &                 RSB_FLAG_FORTRAN_INDICES_INTERFACE        !  optional flags
      INTEGER(KIND=RSB_IDX_KIND),TARGET :: IP(3) = (/1, 2, 4/)
      INTEGER(KIND=RSB_IDX_KIND),TARGET :: JA(3) = (/1, 1, 2/)
      REAL(KIND=8),TARGET :: VA(3) = (/11,21,22/) ! lower triangle coefficients
      REAL(KIND=8),TARGET :: x(2) = (/1, 2/)! reference x 
      REAL(KIND=8),TARGET :: cy(2) = (/215.0, 264.0/)! reference cy after 
      REAL(KIND=8),TARGET :: y(2) = (/3, 4/)! y will be overwritten
      TYPE(C_PTR),TARGET :: mtxAp = C_NULL_PTR ! matrix pointer
      REAL(KIND=8) :: tmax = 2.0 ! tuning max time
      INTEGER :: titmax = 2 ! tuning max iterations
      INTEGER,TARGET :: ont = 0     ! optimal number of threads
      TYPE(C_PTR),PARAMETER :: EO = C_NULL_PTR
      TYPE(C_PTR),PARAMETER :: IO = C_NULL_PTR
      INTEGER,TARGET :: errval

      errval = rsb_lib_init(IO)                ! librsb initialization
      IF (errval.NE.RSB_ERR_NO_ERROR) &
       & STOP "error calling rsb_lib_init"

      res = 0
      mtxAp = rsb_mtx_alloc_from_csr_const(C_LOC(VA),C_LOC(IP),C_LOC(JA)&
       &,nnz,RSB_NUMERICAL_TYPE_DOUBLE,nr,nc,1,1,flags,C_LOC(istat))

      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997

      istat = rsb_file_mtx_save(mtxAp,C_NULL_CHAR)

      ! Structure autotuning:
      istat = rsb_tune_spmm(C_LOC(mtxAp),C_NULL_PTR,C_NULL_PTR,titmax,&
       & tmax,&
       & transt,C_LOC(alpha),C_NULL_PTR,nrhs,order,C_LOC(x),nr,&
       & C_LOC(beta),C_LOC(y),nc)

      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997

      ! Thread count autotuning:
      istat = rsb_tune_spmm(C_NULL_PTR,C_NULL_PTR,C_LOC(ont),titmax,&
       & tmax,&
       & transt,C_LOC(alpha),mtxAp,nrhs,order,C_LOC(x),nr,C_LOC(beta),&
       & C_LOC(y),nc)
      PRINT *, "Optimal number of threads:", ont

      y(:) = (/3, 4/)! restoring reference y (rsb_tune_spmm has overwritten it)
      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997
      
      istat = rsb_file_mtx_save(mtxAp,C_NULL_CHAR)
      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997

      istat = rsb_spmv(transt,C_LOC(alpha),mtxAp,C_LOC(x),incx,&
       & C_LOC(beta),C_LOC(y),incy)
      PRINT *, y
      IF (istat.NE.RSB_ERR_NO_ERROR) GOTO 9997
      DO i = 1, 2
            IF (y(i).NE.cy(i)) PRINT *, "type=d dims=2x2 sym=s diag=g &
      &blocks=1x1 usmv alpha= 4 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF (y(i).NE.cy(i)) GOTO 9997
      END DO
      PRINT*,"type=d dims=2x2 sym=s diag=g blocks=1x1 usmv alpha= 4&
       & beta= 1 incx=1 incy=1 trans=n is ok"
      GOTO 9998

      errval=rsb_lib_exit(EO)                 ! librsb finalization
      IF (errval.NE.RSB_ERR_NO_ERROR)&
       & STOP "error calling rsb_lib_exit"
      PRINT *, "rsb module fortran test is ok"
9997      res = -1
9998      CONTINUE
      mtxAp = rsb_mtx_free(mtxAp)
      IF (istat.NE.RSB_ERR_NO_ERROR) res = -1 
! 9999      CONTINUE
      istat = rsb_perror(C_NULL_PTR,istat)
      end SUBROUTINE rsb_mod_example3

      PROGRAM main
      USE rsb
      IMPLICIT NONE
      INTEGER :: res = RSB_ERR_NO_ERROR, passed = 0, failed = 0

      CALL rsb_mod_example1(res)
      IF (res.LT.0) failed = failed + 1
      IF (res.EQ.0) passed = passed + 1

      CALL rsb_mod_example2(res)
      IF (res.LT.0) failed = failed + 1
      IF (res.EQ.0) passed = passed + 1

      CALL rsb_mod_example3(res)
      IF (res.LT.0) failed = failed + 1
      IF (res.EQ.0) passed = passed + 1
      
      PRINT *, "FAILED:", failed
      PRINT *, "PASSED:", passed
      IF (failed.GT.0) THEN
       STOP 1
      END IF
      END PROGRAM

