! /*
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
! */
! /*
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
! */
! Parallel Sparse BLAS fortran interface testing code
!
      module psb_mvsv_tester
      contains
!
! 
      SUBROUTINE ts_sg_de_usmv_2_n_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=3
      REAL*4 :: beta=1
! A =
! 1 1
! 2 1

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      REAL*4 :: VA(4)=&
          &(/1, 1, 2, 1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/9, 12/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE ts_sg_de_usmv_2_n_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_t_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=3
      REAL*4 :: beta=1
! A =
! 1 1
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*4 :: VA(2)=(/1, 1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/6, 6/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE ts_sg_de_usmv_2_t_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_c_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=3
      REAL*4 :: beta=1
! A =
! 1 1
! 2 2

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      REAL*4 :: VA(4)=&
          &(/1, 1, 2, 2/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/12, 12/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE ts_sg_de_usmv_2_c_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=3
      REAL*4 :: beta=0
! A =
! 1 2
! 2 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      REAL*4 :: VA(3)=&
          &(/1, 2, 2/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/9, 6/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE ts_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=3
      REAL*4 :: beta=0
! A =
! 1 2
! 1 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      REAL*4 :: VA(3)=&
          &(/1, 2, 1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/6, 6/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE ts_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=3
      REAL*4 :: beta=0
! A =
! 1 1
! 0 1

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      REAL*4 :: VA(3)=&
          &(/1, 1, 1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/3, 6/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE ts_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_n_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=1
      REAL*4 :: beta=1
! A =
! 1 1
! 0 2

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      REAL*4 :: VA(3)=&
          &(/1, 1, 2/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/5, 5/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE ts_sg_de_usmv_2_n_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_t_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=1
      REAL*4 :: beta=1
! A =
! 1 2
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*4 :: VA(2)=(/1, 2/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/4, 5/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE ts_sg_de_usmv_2_t_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_c_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=1
      REAL*4 :: beta=1
! A =
! 1 0
! 1 3

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 1, 2/)
      REAL*4 :: VA(3)=&
          &(/1, 1, 3/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/5, 6/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE ts_sg_de_usmv_2_c_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=1
      REAL*4 :: beta=0
! A =
! 1 0
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=1
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(1)=(/1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(1)=(/1/)
      REAL*4 :: VA(1)=(/1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/1, 0/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE ts_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=1
      REAL*4 :: beta=0
! A =
! 1 0
! 0 1

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*4 :: VA(2)=(/1, 1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/1, 1/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE ts_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=1
      REAL*4 :: beta=0
! A =
! 1 0
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=1
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(1)=(/1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(1)=(/1/)
      REAL*4 :: VA(1)=(/1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/1, 0/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE ts_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_n_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-1
      REAL*4 :: beta=1
! A =
! 1 1
! 0 2

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      REAL*4 :: VA(3)=&
          &(/1, 1, 2/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/1, 1/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE ts_sg_de_usmv_2_n_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_t_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-1
      REAL*4 :: beta=1
! A =
! 1 0
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=1
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(1)=(/1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(1)=(/1/)
      REAL*4 :: VA(1)=(/1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/2, 3/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE ts_sg_de_usmv_2_t_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_c_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-1
      REAL*4 :: beta=1
! A =
! 1 1
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*4 :: VA(2)=(/1, 1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/2, 2/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE ts_sg_de_usmv_2_c_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-1
      REAL*4 :: beta=0
! A =
! 1 1
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*4 :: VA(2)=(/1, 1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/-2, 0/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE ts_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-1
      REAL*4 :: beta=0
! A =
! 1 1
! 0 3

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      REAL*4 :: VA(3)=&
          &(/1, 1, 3/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/-1, -4/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE ts_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-1
      REAL*4 :: beta=0
! A =
! 1 0
! 3 2

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 1, 2/)
      REAL*4 :: VA(3)=&
          &(/1, 3, 2/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/-4, -2/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE ts_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_n_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-3
      REAL*4 :: beta=1
! A =
! 1 2
! 2 1

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      REAL*4 :: VA(4)=&
          &(/1, 2, 2, 1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/-6, -6/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE ts_sg_de_usmv_2_n_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_t_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-3
      REAL*4 :: beta=1
! A =
! 1 1
! 3 1

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      REAL*4 :: VA(4)=&
          &(/1, 1, 3, 1/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/-9, -3/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE ts_sg_de_usmv_2_t_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_c_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-3
      REAL*4 :: beta=1
! A =
! 1 0
! 1 2

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 1, 2/)
      REAL*4 :: VA(3)=&
          &(/1, 1, 2/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/-3, -3/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE ts_sg_de_usmv_2_c_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-3
      REAL*4 :: beta=0
! A =
! 1 3
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*4 :: VA(2)=(/1, 3/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/-12, 0/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE ts_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-3
      REAL*4 :: beta=0
! A =
! 1 1
! 3 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      REAL*4 :: VA(3)=&
          &(/1, 1, 3/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/-12, -3/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE ts_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE ts_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_sspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*4 :: alpha=-3
      REAL*4 :: beta=0
! A =
! 1 0
! 4 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 1/)
      REAL*4 :: VA(2)=(/1, 4/)
      REAL*4 :: x(2)=(/1, 1/)! reference x 
      REAL*4 :: cy(2)=(/-15, 0/)! reference cy after 
      REAL*4 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=s dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE ts_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_n_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=3
      REAL*8 :: beta=1
! A =
! 1 1
! 0 4

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      REAL*8 :: VA(3)=&
          &(/1, 1, 4/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/9, 15/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE td_sg_de_usmv_2_n_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_t_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=3
      REAL*8 :: beta=1
! A =
! 1 1
! 1 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      REAL*8 :: VA(3)=&
          &(/1, 1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/9, 6/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE td_sg_de_usmv_2_t_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_c_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=3
      REAL*8 :: beta=1
! A =
! 1 1
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*8 :: VA(2)=(/1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/6, 6/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE td_sg_de_usmv_2_c_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=3
      REAL*8 :: beta=0
! A =
! 1 0
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=1
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(1)=(/1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(1)=(/1/)
      REAL*8 :: VA(1)=(/1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/3, 0/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE td_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=3
      REAL*8 :: beta=0
! A =
! 1 0
! 1 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 1/)
      REAL*8 :: VA(2)=(/1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/6, 0/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE td_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=3
      REAL*8 :: beta=0
! A =
! 1 4
! 1 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      REAL*8 :: VA(3)=&
          &(/1, 4, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/6, 12/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE td_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_n_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=1
      REAL*8 :: beta=1
! A =
! 1 3
! 0 2

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      REAL*8 :: VA(3)=&
          &(/1, 3, 2/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/7, 5/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE td_sg_de_usmv_2_n_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_t_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=1
      REAL*8 :: beta=1
! A =
! 1 0
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=1
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(1)=(/1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(1)=(/1/)
      REAL*8 :: VA(1)=(/1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/4, 3/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE td_sg_de_usmv_2_t_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_c_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=1
      REAL*8 :: beta=1
! A =
! 1 0
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=1
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(1)=(/1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(1)=(/1/)
      REAL*8 :: VA(1)=(/1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/4, 3/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE td_sg_de_usmv_2_c_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=1
      REAL*8 :: beta=0
! A =
! 1 3
! 0 2

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      REAL*8 :: VA(3)=&
          &(/1, 3, 2/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/4, 2/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE td_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=1
      REAL*8 :: beta=0
! A =
! 1 0
! 0 1

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*8 :: VA(2)=(/1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/1, 1/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE td_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=1
      REAL*8 :: beta=0
! A =
! 1 3
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*8 :: VA(2)=(/1, 3/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/1, 3/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE td_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_n_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-1
      REAL*8 :: beta=1
! A =
! 1 2
! 1 3

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      REAL*8 :: VA(4)=&
          &(/1, 2, 1, 3/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/0, -1/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE td_sg_de_usmv_2_n_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_t_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-1
      REAL*8 :: beta=1
! A =
! 1 0
! 0 1

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*8 :: VA(2)=(/1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/2, 2/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE td_sg_de_usmv_2_t_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_c_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-1
      REAL*8 :: beta=1
! A =
! 1 0
! 0 1

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*8 :: VA(2)=(/1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/2, 2/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE td_sg_de_usmv_2_c_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-1
      REAL*8 :: beta=0
! A =
! 1 1
! 3 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      REAL*8 :: VA(3)=&
          &(/1, 1, 3/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/-2, -3/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE td_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-1
      REAL*8 :: beta=0
! A =
! 1 0
! 1 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 1/)
      REAL*8 :: VA(2)=(/1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/-2, 0/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE td_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-1
      REAL*8 :: beta=0
! A =
! 1 2
! 0 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*8 :: VA(2)=(/1, 2/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/-1, -2/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE td_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_n_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-3
      REAL*8 :: beta=1
! A =
! 1 0
! 1 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 1/)
      REAL*8 :: VA(2)=(/1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/0, 0/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE td_sg_de_usmv_2_n_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_t_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-3
      REAL*8 :: beta=1
! A =
! 1 0
! 2 3

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 1, 2/)
      REAL*8 :: VA(3)=&
          &(/1, 2, 3/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/-6, -6/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE td_sg_de_usmv_2_t_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_c_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-3
      REAL*8 :: beta=1
! A =
! 1 0
! 1 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 1/)
      REAL*8 :: VA(2)=(/1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/-3, 3/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE td_sg_de_usmv_2_c_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-3
      REAL*8 :: beta=0
! A =
! 1 0
! 2 0

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 1/)
      REAL*8 :: VA(2)=(/1, 2/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/-3, -6/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE td_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-3
      REAL*8 :: beta=0
! A =
! 1 0
! 0 1

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*8 :: VA(2)=(/1, 1/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/-3, -3/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE td_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE td_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_dspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      REAL*8 :: alpha=-3
      REAL*8 :: beta=0
! A =
! 1 0
! 0 4

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      REAL*8 :: VA(2)=(/1, 4/)
      REAL*8 :: x(2)=(/1, 1/)! reference x 
      REAL*8 :: cy(2)=(/-3, -12/)! reference cy after 
      REAL*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=d dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE td_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_n_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=3
      COMPLEX*8 :: beta=1
! A =
! 1+2i 2+1i
! 0+1i 0+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (2.e0,1.e0), (0.e0,1.e0), (0,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(12.e0,9.e0), (3,9)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tc_sg_de_usmv_2_n_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_t_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=3
      COMPLEX*8 :: beta=1
! A =
! 1+2i 0+0i
! 1+0i 2+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 1, 2/)
      COMPLEX*8 :: VA(3)=&
          &(/(1.e0,2.e0), (1.e0,0.e0), (2,0)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(9.e0,6.e0), (9,0)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tc_sg_de_usmv_2_t_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_c_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=3
      COMPLEX*8 :: beta=1
! A =
! 1+2i 0+0i
! 0+0i 1+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      COMPLEX*8 :: VA(2)=(/(1.e0,2.e0), (1,0)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(6.e0,-6.e0), (6,0)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tc_sg_de_usmv_2_c_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=3
      COMPLEX*8 :: beta=0
! A =
! 1+2i 1+5i
! 2+5i 1+6i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (1.e0,5.e0), (2.e0,5.e0), (1,6)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(6.e0,21.e0), (9,33)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tc_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=3
      COMPLEX*8 :: beta=0
! A =
! 1+2i 2+0i
! 0+0i 3+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      COMPLEX*8 :: VA(3)=&
          &(/(1.e0,2.e0), (2.e0,0.e0), (3,0)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(3.e0,6.e0), (15,0)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tc_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=3
      COMPLEX*8 :: beta=0
! A =
! 1+2i 0+5i
! 1+5i 3+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,5.e0), (1.e0,5.e0), (3,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(6.e0,-21.e0), (9,-21)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tc_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_n_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=1
      COMPLEX*8 :: beta=1
! A =
! 1+2i 3+2i
! 3+2i 3+6i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (3.e0,2.e0), (3.e0,2.e0), (3,6)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(7.e0,4.e0), (9,8)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tc_sg_de_usmv_2_n_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_t_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=1
      COMPLEX*8 :: beta=1
! A =
! 1+2i 2+5i
! 1+5i 1+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (2.e0,5.e0), (1.e0,5.e0), (1,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(5.e0,7.e0), (6,7)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tc_sg_de_usmv_2_t_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_c_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=1
      COMPLEX*8 :: beta=1
! A =
! 1+2i 1+0i
! 0+0i 0+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      COMPLEX*8 :: VA(3)=&
          &(/(1.e0,2.e0), (1.e0,0.e0), (0,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(4.e0,-2.e0), (4,-2)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tc_sg_de_usmv_2_c_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=1
      COMPLEX*8 :: beta=0
! A =
! 1+2i 0+0i
! 2+0i 2+6i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 1, 2/)
      COMPLEX*8 :: VA(3)=&
          &(/(1.e0,2.e0), (2.e0,0.e0), (2,6)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(1.e0,2.e0), (4,6)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tc_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=1
      COMPLEX*8 :: beta=0
! A =
! 1+2i 0+1i
! 4+1i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      COMPLEX*8 :: VA(3)=&
          &(/(1.e0,2.e0), (0.e0,1.e0), (4,1)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(5.e0,3.e0), (0,1)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tc_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=1
      COMPLEX*8 :: beta=0
! A =
! 1+2i 2+0i
! 0+0i 0+4i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      COMPLEX*8 :: VA(3)=&
          &(/(1.e0,2.e0), (2.e0,0.e0), (0,4)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(1.e0,-2.e0), (2,-4)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tc_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_n_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-1
      COMPLEX*8 :: beta=1
! A =
! 1+2i 0+3i
! 0+3i 4+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,3.e0), (0.e0,3.e0), (4,0)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(2.e0,-5.e0), (-1,-3)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tc_sg_de_usmv_2_n_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_t_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-1
      COMPLEX*8 :: beta=1
! A =
! 1+2i 1+3i
! 3+3i 0+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (1.e0,3.e0), (3.e0,3.e0), (0,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(-1.e0,-5.e0), (2,-5)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tc_sg_de_usmv_2_t_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_c_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-1
      COMPLEX*8 :: beta=1
! A =
! 1+2i 2+0i
! 0+0i 0+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      COMPLEX*8 :: VA(3)=&
          &(/(1.e0,2.e0), (2.e0,0.e0), (0,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(2.e0,2.e0), (1,2)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tc_sg_de_usmv_2_c_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-1
      COMPLEX*8 :: beta=0
! A =
! 1+2i 1+3i
! 4+3i 1+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (1.e0,3.e0), (4.e0,3.e0), (1,0)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(-2.e0,-5.e0), (-5,-3)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tc_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-1
      COMPLEX*8 :: beta=0
! A =
! 1+2i 0+3i
! 2+3i 4+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,3.e0), (2.e0,3.e0), (4,0)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(-3.e0,-5.e0), (-4,-3)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tc_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-1
      COMPLEX*8 :: beta=0
! A =
! 1+2i 1+0i
! 0+0i 1+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      COMPLEX*8 :: VA(3)=&
          &(/(1.e0,2.e0), (1.e0,0.e0), (1,0)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(-1.e0,2.e0), (-2,0)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tc_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_n_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-3
      COMPLEX*8 :: beta=1
! A =
! 1+2i 0+0i
! 0+0i 2+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      COMPLEX*8 :: VA(2)=(/(1.e0,2.e0), (2,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(0.e0,-6.e0), (-3,-6)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tc_sg_de_usmv_2_n_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_t_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-3
      COMPLEX*8 :: beta=1
! A =
! 1+2i 1+2i
! 0+2i 1+6i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (1.e0,2.e0), (0.e0,2.e0), (1,6)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(0.e0,-12.e0), (-3,-24)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tc_sg_de_usmv_2_t_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_c_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-3
      COMPLEX*8 :: beta=1
! A =
! 1+2i 2+2i
! 1+2i 4+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (2.e0,2.e0), (1.e0,2.e0), (4,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(-3.e0,12.e0), (-15,12)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tc_sg_de_usmv_2_c_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-3
      COMPLEX*8 :: beta=0
! A =
! 1+2i 0+2i
! 2+2i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      COMPLEX*8 :: VA(3)=&
          &(/(1.e0,2.e0), (0.e0,2.e0), (2,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(-3.e0,-12.e0), (-6,-6)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tc_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-3
      COMPLEX*8 :: beta=0
! A =
! 1+2i 0+6i
! 0+6i 1+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,6.e0), (0.e0,6.e0), (1,0)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(-3.e0,-24.e0), (-3,-18)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tc_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE tc_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_cspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*8 :: alpha=-3
      COMPLEX*8 :: beta=0
! A =
! 1+2i 2+5i
! 0+5i 0+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*8 :: VA(4)=&
          &(/(1.e0,2.e0), (2.e0,5.e0), (0.e0,5.e0), (0,2)/)
      COMPLEX*8 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*8 :: cy(2)=(/(-3.e0,21.e0), (-6,21)/)! reference cy after 
      COMPLEX*8 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=c dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tc_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_n_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=3
      COMPLEX*16 :: beta=1
! A =
! 1+2i 0+0i
! 5+0i 1+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 1, 2/)
      COMPLEX*16 :: VA(3)=&
          &(/(1.e0,2.e0), (5.e0,0.e0), (1,0)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(6.e0,6.e0), (21,0)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tz_sg_de_usmv_2_n_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_t_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=3
      COMPLEX*16 :: beta=1
! A =
! 1+2i 0+0i
! 3+0i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 1/)
      COMPLEX*16 :: VA(2)=(/(1.e0,2.e0), (3,0)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(15.e0,6.e0), (3,0)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tz_sg_de_usmv_2_t_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_c_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=3
      COMPLEX*16 :: beta=1
! A =
! 1+2i 0+4i
! 1+4i 1+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,4.e0), (1.e0,4.e0), (1,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(9.e0,-18.e0), (6,-18)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tz_sg_de_usmv_2_c_ap3_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=3
      COMPLEX*16 :: beta=0
! A =
! 1+2i 0+1i
! 0+1i 4+4i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,1.e0), (0.e0,1.e0), (4,4)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(3.e0,9.e0), (12,15)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tz_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=3
      COMPLEX*16 :: beta=0
! A =
! 1+2i 0+3i
! 0+3i 0+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,3.e0), (0.e0,3.e0), (0,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(3.e0,15.e0), (0,15)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tz_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=3
      COMPLEX*16 :: beta=0
! A =
! 1+2i 3+0i
! 0+0i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      COMPLEX*16 :: VA(2)=(/(1.e0,2.e0), (3,0)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(3.e0,-6.e0), (9,0)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tz_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_n_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=1
      COMPLEX*16 :: beta=1
! A =
! 1+2i 3+0i
! 0+0i 2+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 2/)
      COMPLEX*16 :: VA(3)=&
          &(/(1.e0,2.e0), (3.e0,0.e0), (2,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(7.e0,2.e0), (5,2)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tz_sg_de_usmv_2_n_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_t_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=1
      COMPLEX*16 :: beta=1
! A =
! 1+2i 2+2i
! 0+2i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      COMPLEX*16 :: VA(3)=&
          &(/(1.e0,2.e0), (2.e0,2.e0), (0,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(4.e0,4.e0), (5,2)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tz_sg_de_usmv_2_t_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_c_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=1
      COMPLEX*16 :: beta=1
! A =
! 1+2i 0+3i
! 0+3i 1+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,3.e0), (0.e0,3.e0), (1,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(4.e0,-5.e0), (4,-5)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tz_sg_de_usmv_2_c_ap1_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=1
      COMPLEX*16 :: beta=0
! A =
! 1+2i 4+3i
! 2+3i 2+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (4.e0,3.e0), (2.e0,3.e0), (2,0)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(5.e0,5.e0), (4,3)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tz_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=1
      COMPLEX*16 :: beta=0
! A =
! 1+2i 1+2i
! 2+2i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      COMPLEX*16 :: VA(3)=&
          &(/(1.e0,2.e0), (1.e0,2.e0), (2,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(3.e0,4.e0), (1,2)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tz_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=1
      COMPLEX*16 :: beta=0
! A =
! 1+2i 0+1i
! 0+1i 0+6i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,1.e0), (0.e0,1.e0), (0,6)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(1.e0,-3.e0), (0,-7)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tz_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_n_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-1
      COMPLEX*16 :: beta=1
! A =
! 1+2i 0+1i
! 1+1i 0+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,1.e0), (1.e0,1.e0), (0,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(2.e0,-3.e0), (2,-3)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tz_sg_de_usmv_2_n_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_t_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-1
      COMPLEX*16 :: beta=1
! A =
! 1+2i 5+1i
! 0+1i 3+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (5.e0,1.e0), (0.e0,1.e0), (3,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(2.e0,-3.e0), (-5,-3)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tz_sg_de_usmv_2_t_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_c_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-1
      COMPLEX*16 :: beta=1
! A =
! 1+2i 0+0i
! 0+0i 3+4i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      COMPLEX*16 :: VA(2)=(/(1.e0,2.e0), (3,4)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(2.e0,2.e0), (0,4)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tz_sg_de_usmv_2_c_anr1_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-1
      COMPLEX*16 :: beta=0
! A =
! 1+2i 0+2i
! 0+2i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      COMPLEX*16 :: VA(3)=&
          &(/(1.e0,2.e0), (0.e0,2.e0), (0,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(-1.e0,-4.e0), (0,-2)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tz_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-1
      COMPLEX*16 :: beta=0
! A =
! 1+2i 1+3i
! 0+3i 4+6i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (1.e0,3.e0), (0.e0,3.e0), (4,6)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(-1.e0,-5.e0), (-5,-9)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tz_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-1
      COMPLEX*16 :: beta=0
! A =
! 1+2i 1+0i
! 2+0i 2+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (1.e0,0.e0), (2.e0,0.e0), (2,0)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(-3.e0,2.e0), (-3,0)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tz_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_n_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-3
      COMPLEX*16 :: beta=1
! A =
! 1+2i 0+7i
! 2+7i 0+4i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (0.e0,7.e0), (2.e0,7.e0), (0,4)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(0.e0,-27.e0), (-3,-33)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tz_sg_de_usmv_2_n_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_t_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-3
      COMPLEX*16 :: beta=1
! A =
! 1+2i 0+1i
! 0+1i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      COMPLEX*16 :: VA(3)=&
          &(/(1.e0,2.e0), (0.e0,1.e0), (0,1)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(0.e0,-9.e0), (3,-3)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tz_sg_de_usmv_2_t_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_c_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-3
      COMPLEX*16 :: beta=1
! A =
! 1+2i 0+3i
! 0+3i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=3
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(3)=&
          &(/1, 1, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(3)=&
          &(/1, 2, 1/)
      COMPLEX*16 :: VA(3)=&
          &(/(1.e0,2.e0), (0.e0,3.e0), (0,3)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(0.e0,15.e0), (3,9)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tz_sg_de_usmv_2_c_anr3_bp1_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='n'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-3
      COMPLEX*16 :: beta=0
! A =
! 1+2i 1+2i
! 0+2i 0+2i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=4
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(4)=&
          &(/1, 1, 2, 2/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(4)=&
          &(/1, 2, 1, 2/)
      COMPLEX*16 :: VA(4)=&
          &(/(1.e0,2.e0), (1.e0,2.e0), (0.e0,2.e0), (0,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(-6.e0,-12.e0), (0,-12)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is ok"
      END SUBROUTINE tz_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='t'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-3
      COMPLEX*16 :: beta=0
! A =
! 1+2i 0+0i
! 0+0i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=1
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(1)=(/1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(1)=(/1/)
      COMPLEX*16 :: VA(1)=(/(1,2)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(-3.e0,-6.e0), (0,0)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is ok"
      END SUBROUTINE tz_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1 
! 
      SUBROUTINE tz_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      USE psb_base_mod
      IMPLICIT NONE
      CHARACTER(LEN=*) :: afmt
      TYPE(psb_zspmat_type) :: a
      TYPE(psb_desc_type)   :: desc_a
      INTEGER            :: ictxt, iam=-1, np=-1
      INTEGER            :: info=-1
      
      INTEGER::errval,istat=0,i
      CHARACTER::transA='c'
      INTEGER :: incx=1
      INTEGER :: incy=1
      COMPLEX*16 :: alpha=-3
      COMPLEX*16 :: beta=0
! A =
! 1+2i 4+0i
! 0+0i 0+0i

      ! declaration of VA,IA,JA 
      INTEGER(KIND=RSB_IDX_KIND) :: nnz=2
      INTEGER(KIND=RSB_IDX_KIND) :: nr=2
      INTEGER(KIND=RSB_IDX_KIND) :: nc=2
      INTEGER(KIND=RSB_IDX_KIND) :: IA(2)=(/1, 1/)
      INTEGER(KIND=RSB_IDX_KIND) :: JA(2)=(/1, 2/)
      COMPLEX*16 :: VA(2)=(/(1.e0,2.e0), (4,0)/)
      COMPLEX*16 :: x(2)=(/1, 1/)! reference x 
      COMPLEX*16 :: cy(2)=(/(-3.e0,6.e0), (-12,0)/)! reference cy after 
      COMPLEX*16 :: y(2)=(/3, 3/)! y will be overwritten

      errval=0
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            info=-1
            GOTO 9999
      ENDIF
      CALL psb_barrier(ictxt)
      CALL psb_cdall(ictxt,desc_a,info,nl=nr)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spall(a,desc_a,info,nnz=nnz)
      IF (info .NE. 0)GOTO 9996
      CALL psb_barrier(ictxt)
      CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_cdasb(desc_a,info)
      IF (info .NE. 0)GOTO 9996
      CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
      IF(info.NE.0)PRINT *,"matrix assembly failed"
      IF(info.NE.0)GOTO 9996
      
      CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)
      IF(info.NE.0)PRINT *,"psb_spmm failed"
      IF(info.NE.0)GOTO 9996
      DO i=1,2
            IF(y(i).NE.cy(i))PRINT*,"results mismatch:",y,"instead of",cy
            IF(y(i).NE.cy(i))info=-1
            IF(y(i).NE.cy(i))GOTO 9996
      ENDDO
9996      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_spfree(a,desc_a,info)
      IF (info .NE. 0)GOTO 9997
9997      CONTINUE
      IF(info .NE. 0)errval=errval+1
      CALL psb_cdfree(desc_a,info)
      IF (info .NE. 0)GOTO 9998
9998      CONTINUE
      IF(info .NE. 0)errval=errval+1
9999      CONTINUE
      IF(info .NE. 0)errval=errval+1
            IF(errval.NE.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is not ok"
            IF(errval.EQ.0)PRINT*,"type=z dims=2x2 sym=g diag=e blocks=1x1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is ok"
      END SUBROUTINE tz_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1 
!
      end module psb_mvsv_tester
!
