      SUBROUTINE GAUSS(N,X,R,W,ZF)
!----------------------------------------------------------------------------------------------------------------------
!@sum    Perform an n-point quadrature from 2n moments to yield  n abscissas and n weights.
!@auth   Susanne Bauer/Doug Wright
!----------------------------------------------------------------------------------------------------------------------
      use GaussianQuadrature_mod, only: OrthoPolyCoeffs, GaussianCoeffs
      IMPLICIT NONE

      ! Arguments.

      INTEGER, INTENT(IN   ) :: N               ! number of quadrature points
      REAL(8), INTENT(INOUT) :: X(2*N)          ! moments
      REAL(8), INTENT(  OUT) :: R(N)            ! abscissas
      REAL(8), INTENT(  OUT) :: W(N)            ! weights
      REAL(8), INTENT(  OUT) :: ZF              ! =0.0 successful quadrature, =1.0 failed quadrature

      ! Local variables.

      INTEGER :: IFAILTQL
      REAL(8) :: A(N),B(N),ANU(2*N),AMU0

      ZF = 0.0D+00                              ! successful quadrature
      AMU0 = X(1)                               ! normalizing moment
      ANU(:) = X(:)/AMU0                        ! normalize the moments
      CALL OrthoPolyCoeffs(N,ANU,A,B)
      CALL GaussianCoeffs(N,A,B,AMU0,R,W,IFAILTQL)
      IF(     MINVAL(R(:)) .LT. 0.0D+00 ) THEN  ! failed quadrature
        ZF = 1.0D+00
      ELSEIF( MINVAL(W(:)) .LT. 0.0D+00 ) THEN  ! failed quadrature
        ZF = 1.0D+00
      ELSEIF( IFAILTQL .GT. 0 ) THEN            ! failed quadrature
        ZF = 1.0D+00
      ENDIF
      RETURN
      END

      SUBROUTINE GAUSSINV(N,X,W,U)
!----------------------------------------------------------------------------------------------------------------------
!     DLW: 091206: Computes the moments from the abscissas and weights.
!                  This routine is independent of the units used.
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments.

      INTEGER :: N          ! number of quadrature points
      REAL(8) :: U(2*N)     ! moments
      REAL(8) :: X(N)       ! abscissas
      REAL(8) :: W(N)       ! weights
 
      ! Local variables.

      INTEGER :: I 

      DO I=1, 2*N
        U(I) = SUM( W(:) * ( X(:)**(I-1) ) )
      ENDDO

      RETURN
      END SUBROUTINE GAUSSINV


      SUBROUTINE TEST_QUAD
!----------------------------------------------------------------------------------------------------------------------
!     DLW, 091306: Check quadrature routines. 
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: N = 6               ! number of quadrature points
      REAL(8)            :: X(2*N)              ! moments
      REAL(8)            :: R(N)                ! abscissas
      REAL(8)            :: W(N)                ! weights
      REAL(8)            :: ZF                  ! =0.0 successful quadrature, =1.0 failed quadrature
      INTEGER            :: I,K
      REAL(8)            :: N0,DG,SIGMAG,SG 

      ZF = 0.0D+00
      N0 = 1.0D+03
      DG = 0.1D+00
      SIGMAG = 1.6D+00
      SG = EXP( 0.5D+00 * ( LOG(SIGMAG) )**2 )
      DO I=1, 2*N
        K = I-1
        X(I) = N0 * DG**K * SG**(K*K)
      ENDDO
      WRITE(*,90) X(:)
      CALL GAUSS(N,X,R,W,ZF)
      WRITE(*,90) R(:),W(:)
      CALL GAUSSINV(N,R,W,X)
      WRITE(*,90) X(:)

90    FORMAT(50D18.10)
      RETURN
      END SUBROUTINE TEST_QUAD


