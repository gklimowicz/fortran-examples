!     To calculate weights and abscissas of a Gauss-Laguerre quadrature formula
!	This subroutine requires GAMMA(X), which is calculated to a
!	precision of 1.E-15 only and hence the computed weights and
!	abscissas will also have only that accuracy. For special values
!	of ALP if the integral of weight function RI0 can
!	be calculated accurately it will be possible to get higher
!	accuracy if RI0 is replaced by accurate integral.
!	In particular, for Gauss-Laguerre quarature formula (ALP=0)
!	set RI0=1 to get accurate weights and abscissas.
!
!     N : (input) Number of points in the required quadrature formula
!     ALP : (input) Exponent of x in weight function, w(x)=x**ALP*EXP(-x)
!     W : (output) Array of length N, which will contain the weights
!     A : (output) Array of length N containing the abscissas
!     IER : (output) Error parameter, IER=0 implies successful execution
!		IER=302 implies N.LE.0
!		IER=313 implies ALP.LE.-1
!		IER=321 implies that some coefficient becomes imaginary
!			during calculations.
!		In both these cases calculations are abandoned.
!		Other values may be set by TQL2
!     WK : Real array of length N*(N+2)+3*(N+1) used as scratch space
!
!     Required routines : GAUSRC, TQL2, GAMMA
 
      SUBROUTINE LAGURW(N,ALP,W,A,IER,WK)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      DIMENSION A(N),W(N),WK(3,*)
 
      IF(ALP.LE.-1.0) THEN
        IER=313
        RETURN
      ENDIF

      RI0=GAMMA(ALP+1.0)
      DO 1500 I=1,N+1
!     Coefficients of recurrence relation
        WK(1,I)=-(ALP+I)/I
        WK(2,I)=(ALP+2*I-1)*(ALP+I)/I
        WK(3,I)=(I-1.Q0+ALP)**2*(ALP+I)/I
1500  CONTINUE
 
      CALL GAUSRC(N,W,A,WK,RI0,IER,WK(1,N+2))
      END
