!     To calculate weights and abscissas of a Gauss-Legendre quadrature formula
!
!     N : (input) Number of points in the required quadrature formula
!     W : (output) Array of length N, which will contain the weights
!     A : (output) Array of length N containing the abscissas
!     IER : (output) Error parameter, IER=0 implies successful execution
!		IER=302 implies N.LE.0
!		IER=321 implies that some coefficient becomes imaginary
!			during calculations.
!		In both these cases calculations are abandoned.
!		Other values may be set by TQL2
!     WK : Real array of length N*(N+2)+3*(N+1) used as scratch space
!
!     Required routines : GAUSRC, TQL2
 
      SUBROUTINE GAULEG(N,W,A,IER,WK)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION A(N),W(N),WK(3,*)
 
      RI0=2.0
      DO 1500 I=1,N+1
!     Coefficients of recurrence relation
        WK(1,I)=(2*I-1.D0)/I
        WK(2,I)=0.0
        WK(3,I)=(I-1.D0)/I
1500  CONTINUE
 
      CALL GAUSRC(N,W,A,WK,RI0,IER,WK(1,N+2))
      END
