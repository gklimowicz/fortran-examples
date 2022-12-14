      MODULE TRIDIAG_MOD
!@sum TRIDIAG_MOD contains subroutine TRIDIAG
      PRIVATE
      PUBLIC TRIDIAG
      PUBLIC TRIDIAG_CYCLIC
      PUBLIC TRIDIAG_NEW

      Interface Tridiag_new
      Module Procedure tridiag
      Module procedure tridiag_2d_glob
#ifndef OFFLINE_RUN
      Module procedure tridiag_2d_dist_new
      Module procedure tridiag_3d_dist_new
#endif
      End Interface

      contains

#ifdef OPTIMIZED_TRIDIAG
      SUBROUTINE TRIDIAG(A,B,C,R,U,N)
!@sum  TRIDIAG  solves a tridiagonal matrix equation (A,B,C)U=R
!@auth Numerical Recipes
      IMPLICIT NONE
      INTEGER, PARAMETER :: NMAX = 8000  !@var NMAX workspace
      INTEGER, INTENT(IN):: N         !@var N    dimension of arrays
      REAL*8, INTENT(IN) :: A(N)   !@var A    coefficients of u_i-1
      REAL*8, INTENT(IN) :: B(N)   !@var B    coefficients of u_i
      REAL*8, INTENT(IN) :: C(N)   !@var C    coefficients of u_i+1
      REAL*8, INTENT(IN) :: R(N)   !@var R    RHS vector
      REAL*8, INTENT(OUT):: U(N)   !@var U    solution vector
      REAL*8 :: BET, byBET                !@var BET  work variable
      !REAL*8 :: GAM(Size(A))       !@var GAM  work array
      REAL*8 :: GAM(NMAX)       !@var GAM  work array
      INTEGER :: J                 !@var J    loop variable

      IF ( N > NMAX )
     &     call stop_model("TRIDIAG: N > NMAX, increase NMAX",255)

      BET=B(1)
      IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
      byBET = 1d0/BET
      U(1)=R(1) *byBET
      DO J=2,N
        GAM(J)=C(J-1) *byBET
        BET=B(J)-A(J)*GAM(J)
        IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
        byBET = 1d0/BET
        U(J)=(R(J)-A(J)*U(J-1))*byBET
      END DO
      DO J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
      END DO
      RETURN
      END SUBROUTINE TRIDIAG
#else
      SUBROUTINE TRIDIAG(A,B,C,R,U,N)
!@sum  TRIDIAG  solves a tridiagonal matrix equation (A,B,C)U=R
!@auth Numerical Recipes
      IMPLICIT NONE
c      INTEGER, PARAMETER :: NMAX = 8000  !@var NMAX workspace
      INTEGER, INTENT(IN):: N         !@var N    dimension of arrays
      REAL*8, INTENT(IN) :: A(N)   !@var A    coefficients of u_i-1
      REAL*8, INTENT(IN) :: B(N)   !@var B    coefficients of u_i
      REAL*8, INTENT(IN) :: C(N)   !@var C    coefficients of u_i+1
      REAL*8, INTENT(IN) :: R(N)   !@var R    RHS vector
      REAL*8, INTENT(OUT):: U(N)   !@var U    solution vector
      REAL*8 :: BET                !@var BET  work variable
      REAL*8 :: GAM(Size(A))       !@var GAM  work array
      INTEGER :: J                 !@var J    loop variable

c      IF ( N > NMAX )
c     &     call stop_model("TRIDIAG: N > NMAX, increase NMAX",255)
      BET=B(1)
      IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
      U(1)=R(1)/BET
      DO J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
        U(J)=(R(J)-A(J)*U(J-1))/BET
      END DO
      DO J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
      END DO
      RETURN
      END SUBROUTINE TRIDIAG
#endif

      SUBROUTINE TRIDIAG_cyclic(A,B,C,R,U,N)
!@sum TRIDIAG_cyclic solves a cyclic tridiagonal matrix equation (A,B,C)U=R
!@+   having nonzero A(1) and C(N), using the Thomas algorithm and
!@+   Sherman-Morrison formula.
      IMPLICIT NONE
      INTEGER :: N         !@var N    dimension of arrays
      REAL*8  :: A(N)   !@var A    coefficients of u_i-1
      REAL*8  :: B(N)   !@var B    coefficients of u_i
      REAL*8  :: C(N)   !@var C    coefficients of u_i+1
      REAL*8  :: R(N)   !@var R    RHS vector
      INTENT(IN) :: N,A,B,C,R
      REAL*8, INTENT(OUT):: U(N)   !@var U    solution vector
      REAL*8 :: BET                !@var BET  work variable
      REAL*8 :: GAM(Size(A))       !@var GAM  work array
      REAL*8 :: Q(Size(A))         !@var Q    work array
      INTEGER :: J                 !@var J    loop variable
      REAL*8 :: QCOEFF,A1,B1,C1,R1

      A1 = A(1)
      B1 = B(1)
      C1 = C(1)
      R1 = R(1)
      IF(B1.EQ.1D0) THEN
        A1 = A1*2D0
        B1 = B1*2D0
        C1 = C1*2D0
        R1 = R1*2D0
      ENDIF
      BET=B1-1d0
      U(1)=R1/BET
      Q(1)=1d0/BET
      GAM(2)=C1/BET
      DO J=2,N-1
        BET=B(J)-A(J)*GAM(J)
        IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
c        IF (BET.eq.0) stop 'BET==0'
        U(J)=(R(J)-A(J)*U(J-1))/BET
        Q(J)=     -A(J)*Q(J-1) /BET
        GAM(J+1) = C(J)/BET
      END DO
      J=N
        BET=B(J)-A(J)*GAM(J)-A1*C(N)
        IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
c        IF (BET.eq.0) stop 'BET==0'
        U(J)=(R(J)-A(J)*U(J-1))/BET
        Q(J)=(C(N)-A(J)*Q(J-1))/BET
      DO J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
        Q(J)=Q(J)-GAM(J+1)*Q(J+1)
      END DO
      BET=1d0+Q(1)+A1*Q(N)
        IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
c        IF (BET.eq.0) stop 'BET==0'
      QCOEFF = (U(1)+A1*U(N))/BET
      DO J=1,N
        U(J) = U(J) - QCOEFF*Q(J)
      ENDDO
      RETURN
      END SUBROUTINE TRIDIAG_cyclic

      SUBROUTINE TRIDIAG_2D_GLOB(A, B, C, R, U)
!@sum  TRIDIAG  solves an array of tridiagonal matrix equations (A,B,C)U=R
!@auth Numerical Recipes
      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: A(:,:),B(:,:),C(:,:),R(:,:)
      REAL*8, INTENT(OUT)   :: U(:,:)

      REAL*8 :: BET                   !@var BET  work variable
      REAL*8 :: GAM(Size(A,2))        !@var GAM  work array

      Integer :: i, j
      Integer :: N, N_arr

      N = SIZE(A,2)
      n_arr = Size(A,1)
      do i=1,n_arr
        BET=B(i,1)
        IF (BET.eq.0) then
         print*, "TRIDIAG_2D_GLOB: DENOMINATOR = ZERO"
         stop
        end if
        U(i,1)=R(i,1)/BET
        DO J=2,N
          GAM(J)=C(i,J-1)/BET
          BET=B(i,J)-A(i,J)*GAM(J)
          IF (BET.eq.0) then
           print*, "TRIDIAG: DENOMINATOR = ZERO"
           stop
          end if
          U(i,J)=( R(i,J) - A(i,J)*U(i,J-1) )/BET
        END DO
        DO J=N-1,1,-1
          U(i,J)=U(i,J)-GAM(J+1)*U(i,J+1)
        END DO
      end do

      RETURN
      END SUBROUTINE TRIDIAG_2D_GLOB

#ifndef OFFLINE_RUN
      SUBROUTINE TRIDIAG_2D_DIST(A_dist, B_dist, C_dist, R_dist,
     &                           U_dist,grid, j_lower, j_upper )
!@sum  TRIDIAG  solves an array of tridiagonal matrix equations (A,B,C)U=R
!@auth Numerical Recipes
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      USE DOMAIN_DECOMP_1D, ONLY : TRANSP
      IMPLICIT NONE

      Type (DIST_GRID), Intent(IN) :: grid
      REAL*8, INTENT(INOUT) :: A_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: B_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: C_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: R_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(OUT)   :: U_dist(:,grid%j_strt_halo:)
      INTEGER, INTENT(IN)   :: J_LOWER, J_UPPER

      REAL*8, ALLOCATABLE :: A_tr(:,:), B_tr(:,:),
     &                       C_tr(:,:), R_tr(:,:),
     &                       U_tr(:,:)

      REAL*8 :: BET                   !@var BET  work variable
      REAL*8 :: GAM(grid%jm_world)        !@var GAM  work array

      Integer :: i, j
      Integer :: N, N_i


! Determine the size of the global arrays
      N = grid%jm_world
      n_i = grid%ni_loc

! Matrix size consistent with array size?
      if ( J_upper > N ) then
        print*, 'TRIDIAG: upper bound of matrix arrays is too large'
        print*, 'j_upper = ', j_upper, 'jm =', n, ' ( need j_upper<=jm)'
        call stop_model('TRIDIAG: j_upper argument too large', 255)
      end if

! Allocate the transposed arrays
      allocate( a_tr(n_i,n), b_tr(n_i,n), c_tr(n_i,n), r_tr(n_i,n) )
      allocate( u_tr(n_i,n) )

! First get the transpose of A,B,C,R
      call transp(grid, a_dist, a_tr)
      call transp(grid, b_dist, b_tr)
      call transp(grid, c_dist, c_tr)
      call transp(grid, r_dist, r_tr)

! Solve
        do i=1,n_i
          BET=B_tr(i,j_lower)
          IF (BET.eq.0) then
           print*, "TRIDIAG_2D_DIST: DENOM. = ZERO  i,j= ", i,' 1'
           stop
          end if
          U_tr(i,j_lower)=R_tr(i,j_lower)/BET
          DO J=j_lower+1, j_upper
            GAM(J)=C_tr(i,J-1)/BET
            BET=B_tr(i,J)-A_tr(i,J)*GAM(J)
            IF (BET.eq.0) then
             print*, "TRIDIAG_2D_DIST: DENOM. = ZERO i,j= ", i, j
             stop
            end if
          U_tr(i,J)=( R_tr(i,J) - A_tr(i,J)*U_tr(i,J-1) )/BET
          END DO
          DO J=j_upper-1,j_lower,-1
            U_tr(i,J)=U_tr(i,J)-GAM(J+1)*U_tr(i,J+1)
          END DO
        end do


! Transfer the solution to the j-distributed array
        call transp( grid, u_dist, u_tr, reverse=.true.)

      deallocate( a_tr, b_tr, c_tr, r_tr, u_tr )

      RETURN
      END SUBROUTINE TRIDIAG_2D_DIST

      SUBROUTINE TRIDIAG_3D_DIST(A_dist, B_dist, C_dist, R_dist,
     &                           U_dist,grid, j_lower, j_upper )
!@sum  TRIDIAG  solves an array of tridiagonal matrix equations (A,B,C)U=R
!@auth Numerical Recipes
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      USE DOMAIN_DECOMP_1D, ONLY : TRANSP
      IMPLICIT NONE

      Type (DIST_GRID), Intent(IN) :: grid
      REAL*8, INTENT(INOUT) :: A_dist(:,grid%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: B_dist(:,grid%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: C_dist(:,grid%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: R_dist(:,grid%j_strt_halo:,:)
      REAL*8, INTENT(OUT)   :: U_dist(:,grid%j_strt_halo:,:)
      INTEGER, INTENT(IN)   :: J_LOWER, J_UPPER

      REAL*8, ALLOCATABLE :: A_tr(:,:,:), B_tr(:,:,:),
     &                       C_tr(:,:,:), R_tr(:,:,:),
     &                       U_tr(:,:,:)

      REAL*8 :: BET                   !@var BET  work variable
      REAL*8 :: GAM(grid%jm_world)        !@var GAM  work array

      Integer :: i, j, l
      Integer :: N, N_i, N_l


! Determine the size of the global arrays
      N = grid%jm_world
      n_i = grid%ni_loc
      n_l = size(A_dist,3)

! Matrix size consistent with array size?
      if ( J_upper > N ) then
        print*, 'TRIDIAG: upper bound of matrix arrays is too large'
        print*, 'j_upper = ', j_upper, 'jm =', n, ' ( need j_upper<=jm)'
        call stop_model('TRIDIAG: j_upper argument too large', 255)
      end if

! Allocate the transposed arrays
      allocate( a_tr(n_i,n,n_l), b_tr(n_i,n,n_l), c_tr(n_i,n,n_l) )
      allocate( r_tr(n_i,n,n_l), u_tr(n_i,n,n_l) )

! First get the transpose of A,B,C,R
      call transp(grid, a_dist, a_tr)
      call transp(grid, b_dist, b_tr)
      call transp(grid, c_dist, c_tr)
      call transp(grid, r_dist, r_tr)

! Solve
      do l=1,n_l
        do i=1,n_i
          BET=B_tr(i,j_lower,l)
          IF (BET.eq.0) then
           print*, "TRIDIAG_3D_DIST: DENOM. = ZERO  i,j= ", i,' 1'
           stop
          end if
          U_tr(i,j_lower,l)=R_tr(i,j_lower,l)/BET
          DO J=j_lower+1, j_upper
            GAM(J)=C_tr(i,J-1,l)/BET
            BET=B_tr(i,J,l)-A_tr(i,J,l)*GAM(J)
            IF (BET.eq.0) then
             print*, "TRIDIAG_3D_DIST: DENOM. = ZERO i,j= ", i, j
             stop
            end if
          U_tr(i,J,l)=( R_tr(i,J,l) - A_tr(i,J,l)*U_tr(i,J-1,l) )/BET
          END DO
          DO J=j_upper-1,j_lower,-1
            U_tr(i,J,l)=U_tr(i,J,l)-GAM(J+1)*U_tr(i,J+1,l)
          END DO
        end do
      end do

! Transfer the solution to the j-distributed array
        call transp( grid, u_dist, u_tr, reverse=.true.)

      deallocate( a_tr, b_tr, c_tr, r_tr, u_tr )

      RETURN
      END SUBROUTINE TRIDIAG_3D_DIST

      SUBROUTINE TRIDIAG_2D_DIST_new(A_dist, B_dist, C_dist, R_dist,
     &                           U_dist,grid, j_lower, j_upper )
!@sum  TRIDIAG  solves an array of tridiagonal matrix equations (A,B,C)U=R
!@auth Numerical Recipes
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      USE DOMAIN_DECOMP_1D, ONLY : TRANSP,TRANSPOSE_COLUMN
      IMPLICIT NONE

      Type (DIST_GRID), Intent(IN) :: grid
      REAL*8, INTENT(INOUT) :: A_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: B_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: C_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: R_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(OUT)   :: U_dist(:,grid%j_strt_halo:)
      INTEGER, INTENT(IN)   :: J_LOWER, J_UPPER

      REAL*8, ALLOCATABLE :: ABCR(:,:,:,:),ABCR_tr(:,:,:,:),
     &                       U_tr(:,:)

      REAL*8 :: BET
      REAL*8, allocatable :: BYBET(:),GAM(:,:) !@var BET,GAM  work arrays

      Integer :: i, j
      Integer :: N, N_i, IM


! Determine the size of the global arrays
      N = grid%jm_world
      n_i = grid%ni_loc
      allocate( bybet(n_i), gam(n_i,n) )

! Matrix size consistent with array size?
      if ( J_upper > N ) then
        print*, 'TRIDIAG: upper bound of matrix arrays is too large'
        print*, 'j_upper = ', j_upper, 'jm =', n, ' ( need j_upper<=jm)'
        call stop_model('TRIDIAG: j_upper argument too large', 255)
      end if

! Copy j-distributed A,B,C,R into single array to do all tranposes together
      IM = grid%im_world
      allocate(abcr(4,IM,grid%j_strt_halo:grid%j_stop_halo,1))
      do j=grid%j_strt,grid%j_stop
        do i=1,IM
          abcr(1,i,j,1) = a_dist(i,j)
          abcr(2,i,j,1) = b_dist(i,j)
          abcr(3,i,j,1) = c_dist(i,j)
          abcr(4,i,j,1) = r_dist(i,j)
        enddo
      enddo

! Allocate the transposed arrays
      allocate( abcr_tr(4,n_i,n,1) )
      allocate( u_tr(n_i,n) )

! Do the transpose of A,B,C,R
      call transpose_column(grid, abcr, abcr_tr)

! Solve
      do i=1,n_i
        BET=ABCR_tr(2,i,j_lower,1)
        IF (BET.eq.0) then
          print*, "TRIDIAG_2D_DIST: DENOM. = ZERO  i,j= ", i,' 1'
          stop
        end if
        BYBET(I) = 1D0/BET
        U_tr(i,j_lower)=ABCR_tr(4,i,j_lower,1)*BYBET(I)
      enddo
      do J=j_lower+1, j_upper
        do i=1,n_i
          GAM(I,J)=ABCR_tr(3,i,J-1,1)*BYBET(I)
          BET=ABCR_tr(2,i,J,1)-ABCR_tr(1,i,J,1)*GAM(I,J)
          IF (BET.eq.0) then
            print*, "TRIDIAG_2D_DIST: DENOM. = ZERO i,j= ", i, j
            stop
          end if
          BYBET(I) = 1D0/BET
          U_tr(i,J)=( ABCR_tr(4,i,J,1)-ABCR_tr(1,i,J,1)*U_tr(i,J-1) )
     &         *BYBET(I)
        enddo
      enddo
      do J=j_upper-1,j_lower,-1
        do i=1,n_i
          U_tr(i,J)=U_tr(i,J)-GAM(I,J+1)*U_tr(i,J+1)
        enddo
      enddo

! Transfer the solution to the j-distributed array
      call transp( grid, u_dist, u_tr, reverse=.true.)

      deallocate( abcr, abcr_tr, u_tr, bybet, gam )

      RETURN
      END SUBROUTINE TRIDIAG_2D_DIST_new

      SUBROUTINE TRIDIAG_3D_DIST_new(A_dist, B_dist, C_dist, R_dist,
     &                           U_dist,grid, j_lower, j_upper )
!@sum  TRIDIAG  solves an array of tridiagonal matrix equations (A,B,C)U=R
!@auth Numerical Recipes
      USE DOMAIN_DECOMP_1D, ONLY : DIST_GRID
      USE DOMAIN_DECOMP_1D, ONLY : TRANSP,TRANSPOSE_COLUMN
      IMPLICIT NONE

      Type (DIST_GRID), Intent(IN) :: grid
      REAL*8, INTENT(INOUT) :: A_dist(:,grid%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: B_dist(:,grid%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: C_dist(:,grid%j_strt_halo:,:)
      REAL*8, INTENT(INOUT) :: R_dist(:,grid%j_strt_halo:,:)
      REAL*8, INTENT(OUT)   :: U_dist(:,grid%j_strt_halo:,:)
      INTEGER, INTENT(IN)   :: J_LOWER, J_UPPER

      REAL*8, ALLOCATABLE :: ABCR(:,:,:,:),ABCR_tr(:,:,:,:),
     &                       U_tr(:,:,:)

      REAL*8 :: BET
      REAL*8, allocatable :: BYBET(:),GAM(:,:) !@var BET,GAM  work arrays

      Integer :: i, j, l
      Integer :: N, N_i, N_l, IM


! Determine the size of the global arrays
      N = grid%jm_world
      n_i = grid%ni_loc
      n_l = size(A_dist,3)
      allocate( bybet(n_i), gam(n_i,n) )

! Matrix size consistent with array size?
      if ( J_upper > N ) then
        print*, 'TRIDIAG: upper bound of matrix arrays is too large'
        print*, 'j_upper = ', j_upper, 'jm =', n, ' ( need j_upper<=jm)'
        call stop_model('TRIDIAG: j_upper argument too large', 255)
      end if

! Copy j-distributed A,B,C,R into single array to do all tranposes together
      IM = grid%im_world
      allocate(abcr(4,IM,grid%j_strt_halo:grid%j_stop_halo,n_l))
      do l=1,n_l
        do j=grid%j_strt,grid%j_stop
          do i=1,IM
            abcr(1,i,j,l) = a_dist(i,j,l)
            abcr(2,i,j,l) = b_dist(i,j,l)
            abcr(3,i,j,l) = c_dist(i,j,l)
            abcr(4,i,j,l) = r_dist(i,j,l)
          enddo
        enddo
      enddo

! Allocate the transposed arrays
      allocate( abcr_tr(4,n_i,n,n_l) )
      allocate( u_tr(n_i,n,n_l) )

! Do the transpose of A,B,C,R
      call transpose_column(grid, abcr, abcr_tr)

! Solve
      do l=1,n_l
      do i=1,n_i
        BET=ABCR_tr(2,i,j_lower,l)
        IF (BET.eq.0) then
          print*, "TRIDIAG_3D_DIST: DENOM. = ZERO  i,j= ", i,' 1'
          stop
        end if
        BYBET(I) = 1D0/BET
        U_tr(i,j_lower,l)=ABCR_tr(4,i,j_lower,l)*BYBET(I)
      enddo
      do J=j_lower+1, j_upper
        do i=1,n_i
          GAM(I,J)=ABCR_tr(3,i,J-1,l)*BYBET(I)
          BET=ABCR_tr(2,i,J,l)-ABCR_tr(1,i,J,l)*GAM(I,J)
          IF (BET.eq.0) then
            print*, "TRIDIAG_3D_DIST: DENOM. = ZERO i,j= ", i, j
            stop
          end if
          BYBET(I) = 1D0/BET
          U_tr(i,J,l)=(ABCR_tr(4,i,J,l)-ABCR_tr(1,i,J,l)*U_tr(i,J-1,l))
     &         *BYBET(I)
        enddo
      enddo
      do J=j_upper-1,j_lower,-1
        do i=1,n_i
          U_tr(i,J,l)=U_tr(i,J,l)-GAM(I,J+1)*U_tr(i,J+1,l)
        enddo
      enddo
      enddo

! Transfer the solution to the j-distributed array
      call transp( grid, u_dist, u_tr, reverse=.true.)

      deallocate( abcr, abcr_tr, u_tr, bybet, gam )

      RETURN
      END SUBROUTINE TRIDIAG_3D_DIST_new
#endif



      END MODULE TRIDIAG_MOD
