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

!> @file
!! @brief Sparse BLAS-based usage:
!!  uscr_begin(), usgp(), ussp(), usmv(), ...
!! \include fortran.F90

      SUBROUTINE blas_sparse_mod_example(res)
      USE blas_sparse
      USE rsb ! For the second part of the example and RSB_IDX_KIND
      USE iso_c_binding
      IMPLICIT NONE
      INTEGER :: i, j
      INTEGER(KIND=RSB_BLAS_IDX_KIND) :: istat = 0, res
      TYPE(C_PTR),TARGET :: mtxAp = C_NULL_PTR ! matrix pointer
      INTEGER :: A
      INTEGER,PARAMETER :: transn = blas_no_trans
      INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: incx = 1
      INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: incy = 1
      REAL(KIND=8),PARAMETER :: alpha = 3
! Symmetric (declared via lower triangle) matrix based example, e.g.:
! 1 0
! 1 1
      ! declaration of VA,IA,JA 
      !INTEGER,PARAMETER :: nr = 100
      INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: nr = 20
      INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: nc = nr
      INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: nnz = (nr*(nr+1))/2 ! half the square
      INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: nrhs = 1
      INTEGER(KIND=RSB_IDX_KIND) :: nt = 0
      INTEGER :: ic, ir
      ! For compatibility with (or compassion for) flang-13, we turn these three lines off:
      ! REAL(KIND=8),PARAMETER :: VA(nnz) = (/ ((1, ic=1,ir), ir=1,nr ) /) ! (/1, 1, 1/)
      ! INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: IA(nnz) = (/ (((ir), ic=1,ir), ir=1,nr ) /) ! (/1, 2, 2/)
      ! INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: JA(nnz) = (/ (((ic), ic=1,ir), ir=1,nr ) /) ! (/1, 1, 2/)
      ! and use these other ones instead, initializing VA,IA,JA later:
      REAL(KIND=8) :: VA(nnz)
      INTEGER(KIND=RSB_IDX_KIND) :: IA(nnz)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(nnz)
      REAL(KIND=8) :: x(nc,nrhs) = RESHAPE((/((1), ic=1,nc*nrhs)/),[nc,nrhs]) ! reference x ! (/1, 1/)
      REAL(KIND=8),PARAMETER :: cy(nr,nrhs) = RESHAPE((/((alpha+alpha*nr), ir=1,nr*nrhs)/),[nr,nrhs]) ! reference cy after ! (/9, 9/)
      REAL(KIND=8) :: y(nr,nrhs) = RESHAPE((/((alpha), ir=1,nr*nrhs)/),[nr,nrhs]) ! y will be overwritten ! (/3, 3/)
      VA(:) = (/ ((1, ic=1,ir), ir=1,nr ) /) ! (/1, 1, 1/)
      IA(:) = (/ (((ir), ic=1,ir), ir=1,nr ) /) ! (/1, 2, 2/)
      JA(:) = (/ (((ic), ic=1,ir), ir=1,nr ) /) ! (/1, 1, 2/)
      ! First example part: pure blas_sparse code.
      res = 0
      CALL duscr_begin(nr,nc,A,res)
      IF (res.NE.0) GOTO 9999
      CALL ussp(A,blas_lower_symmetric,istat)
      IF (istat.NE.0) GOTO 9997
      CALL ussp(A,blas_rsb_spmv_autotuning_on,istat) ! (experimental) turns auto-tuning + thread setting on
      IF (istat.NE.0) PRINT *,"autotuning returned nonzero:", istat &
       &," ...did you enable autotuning ?"
      !
      ! First style example 
      CALL uscr_insert_entries(A,nnz,VA,IA,JA,istat)
      IF (istat.NE.0) GOTO 9997
      CALL uscr_end(A,istat)
      IF (istat.NE.0) GOTO 9997
      ! CALL ussp(A,blas_rsb_duplicates_sum,istat)
      ! CALL uscr_insert_entries(A,nnz,VA,IA,JA,istat) ! uncomment this to activate add of coefficients to pattern
      CALL usgp(A,blas_rsb_spmv_autotuning_on,nt)  ! (experimental)
      IF (nt.NE.0) PRINT*,"autotuner chose ",nt," threads"
      CALL ussp(A,blas_rsb_spmv_autotuning_off,istat) ! (experimental) turns auto-tuning + thread setting off
      IF (istat.NE.0) GOTO 9997

      DO j = 1, nrhs
        CALL usmv(transn,alpha,A,x(:,j),incx,y(:,j),incy,istat)
      END DO
      IF (istat.NE.0) GOTO 9997
      !
      DO j = 1, nrhs
      DO i = 1, nr
        IF (y(i,j).NE.cy(i,j)) PRINT *, "first check results are not ok"
        IF (y(i,j).NE.cy(i,j)) GOTO 9997
      END DO
      END DO
      !
      y(:,:) = alpha ! reset
      !
      ! Second style example 
      CALL ussp(A,blas_rsb_autotune_next_operation,istat) ! (experimental) turns auto-tuning + thread setting on
      IF (istat.NE.0) GOTO 9997
      DO j = 1, nrhs
        CALL usmv(transn,alpha,A,x(:,j),incx,y(:,j),incy,istat)
      END DO
      CALL usmm(blas_colmajor,transn,nrhs, alpha,A,x,nr,y,nc,istat) ! Equivalent to the above (as long as incx=incy=1).
      CALL usmm(blas_colmajor,transn,nrhs,-alpha,A,x,nr,y,nc,istat) ! Subtract the last usmm call contribution.
      IF (istat.NE.0) GOTO 9997
      !
      DO j = 1, nrhs
      DO i = 1, nr
        IF (y(i,j).NE.cy(i,j)) PRINT *,"second check results are not ok"
        IF (y(i,j).NE.cy(i,j)) GOTO 9997
      END DO
      END DO
      !
      PRINT *, "check results are ok"
      
      ! Second part of the example: access to the rsb.h interface via
      ! the ISO C Binding interface.
      mtxAp = rsb_BLAS_get_mtx(A) ! get pointer to rsb structure (as in the rsb.h API)
      IF(nr.LT.5) istat = rsb_file_mtx_save(mtxAp,C_NULL_CHAR) ! write to stdout (only if matrix small enough)

      GOTO 9998
9997      res = -1
9998      CONTINUE
      CALL usds(A,istat)
      IF (istat.NE.0) res = -1
9999      CONTINUE
      end SUBROUTINE blas_sparse_mod_example

      ! Example loading a matrix from file and measuring SPMM.
      SUBROUTINE blas_sparse_io_example(res)
      USE blas_sparse
      USE rsb ! For rsb_blas_file_mtx_load
      USE iso_c_binding
      IMPLICIT NONE
      INTEGER(KIND=RSB_IDX_KIND) :: res ! note: same as RSB_BLAS_IDX_KIND
      INTEGER :: j
      INTEGER(KIND=RSB_BLAS_IDX_KIND) :: istat = 0
      INTEGER :: A
      INTEGER,PARAMETER :: transn = blas_no_trans
      INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: incx = 1
      INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: incy = 1
      COMPLEX(KIND=8),PARAMETER :: alpha = 3
      INTEGER(KIND=RSB_IDX_KIND) :: nr
      INTEGER(KIND=RSB_IDX_KIND) :: nc
      INTEGER(KIND=RSB_IDX_KIND) :: nz
      INTEGER(KIND=RSB_IDX_KIND) :: st
      INTEGER(KIND=RSB_IDX_KIND),PARAMETER :: nrhs = 4
      COMPLEX(KIND=8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: x
      COMPLEX(KIND=8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: y
      CHARACTER(KIND=C_CHAR,LEN=7),TARGET :: filename = 'pd.mtx'//C_NULL_CHAR
      INTEGER(C_SIGNED_CHAR) :: typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX
      REAL(KIND=C_DOUBLE) :: mvt,mmt,omt
      INTEGER(KIND=C_INT),TARGET::IZERO=0

      res = 0

      A = rsb_blas_file_mtx_load(filename, typecode);
      IF (A.EQ.blas_invalid_handle) GOTO 9997

      CALL usgp(A,blas_num_rows,nr)
      CALL usgp(A,blas_num_cols,nc)
      CALL usgp(A,blas_num_nonzeros,nz)

      PRINT*,"Read matrix ",filename(1:6)," ",nr,"x",nc,":",nz
      CALL usgp(A,blas_general,st)
      IF (st .EQ. 1) PRINT*,"Matrix has no symmetry"
      CALL usgp(A,blas_upper_symmetric,st)
      IF (st .EQ. 1) PRINT*,"Matrix is upper symmetric"
      CALL usgp(A,blas_upper_hermitian,st)
      IF (st .EQ. 1) PRINT*,"Matrix is upper hermitian"
      ! ...

      IF (istat.NE.0) GOTO 9997
      WRITE(*,'(a,i0)') "Using NRHS=",nrhs
      ALLOCATE( x(nc,nrhs))
      ALLOCATE( y(nr,nrhs))
      x = 1.0

      y = 0.0
      mvt = -rsb_time()
      DO j = 1, nrhs
        CALL usmv(transn,alpha,A,x(:,j),incx,y(:,j),incy,istat)
      END DO
      IF (istat.NE.0) GOTO 9997
      mvt = mvt + rsb_time()
      WRITE(*,'(a,e12.4,a)') "Repeated USMV took ",mvt," s"
      
      y = 0.0
      mmt = -rsb_time()
      CALL usmm(blas_colmajor,transn,nrhs, alpha,A,x,nr,y,nc,istat)
      IF (istat.NE.0) GOTO 9997
      mmt = mmt + rsb_time()
      WRITE(*,'(a,e12.4,a)') "A single USMM took ",mmt," s"

      WRITE(*,'(a,g11.4,a)')"USMM-to-USMV speed ratio is is ", mvt/mmt, "x"

      PRINT*,"Call auto-tuning routine.."
      ! Change IZERO value to 1 to get verbose tuning again.
      res = rsb_lib_set_opt(RSB_IO_WANT_VERBOSE_TUNING,C_LOC(IZERO))
      IF (res.NE.0) GOTO 9997
      CALL ussp(A,blas_rsb_autotune_next_operation,istat) ! (experimental) turns auto-tuning + thread setting on
      IF (istat.NE.0) GOTO 9997
      CALL usmm(blas_colmajor,transn,nrhs, alpha,A,x,nr,y,nc,istat)
      IF (istat.NE.0) GOTO 9997

      PRINT*,"Repeat measurement."
      y = 0.0
      omt = -rsb_time()
      CALL usmm(blas_colmajor,transn,nrhs, alpha,A,x,nr,y,nc,istat)
      IF (istat.NE.0) GOTO 9997
      omt = omt + rsb_time()
      WRITE(*,'(a,e12.4,a)') "Tuned USMM took ",omt," s"
      WRITE(*,'(a,g11.4,a)')"Tuned-to-untuned speed ratio is is ",mmt/omt,"x"
      
      GOTO 9998
9997      res = -1
9998      CONTINUE
      CALL usds(A,istat)
      IF (istat.NE.0) res = -1
      end SUBROUTINE blas_sparse_io_example

      PROGRAM main
      USE rsb, ONLY: rsb_lib_init, rsb_lib_exit, C_PTR, C_NULL_PTR,&
       & RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE,RSB_IO_WANT_VERBOSE_TUNING,&
       & rsb_lib_set_opt,RSB_IDX_KIND
      ! USE blas_sparse, only: RSB_BLAS_IDX_KIND ! only if using long indices
      USE iso_c_binding
      IMPLICIT NONE
      INTEGER :: passed = 0, failed = 0
      INTEGER(KIND=RSB_IDX_KIND) :: res ! note: same as RSB_BLAS_IDX_KIND
      !TYPE(C_PTR),PARAMETER :: EO = RSB_NULL_EXIT_OPTIONS
      !TYPE(C_PTR),PARAMETER :: IO = RSB_NULL_INIT_OPTIONS
      ! Note: using C_NULL_PTR instead of the previous lines because of http://gcc.gnu.org/bugzilla/show_bug.cgi?id=59411
      TYPE(C_PTR),PARAMETER :: EO = C_NULL_PTR
      TYPE(C_PTR),PARAMETER :: IO = C_NULL_PTR
      INTEGER(KIND=C_INT),TARGET::IONE=1
      res = rsb_lib_init(IO)
      res = rsb_lib_set_opt(RSB_IO_WANT_VERBOSE_TUNING,C_LOC(IONE))
      
      CALL blas_sparse_mod_example(res)
      IF (res.LT.0) failed = failed + 1
      IF (res.EQ.0) passed = passed + 1

      CALL blas_sparse_io_example(res)
      IF (res.LT.0) failed = failed + 1
      IF (res.EQ.0) passed = passed + 1

      res = rsb_lib_exit(EO)
      
      PRINT *, "FAILED:", failed
      PRINT *, "PASSED:", passed
      IF (failed .GT. 0) THEN
       STOP 1
      END IF
      END PROGRAM
