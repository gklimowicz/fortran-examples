!     To calculate weights and abscissas of a Gauss-Jacobi quadrature formula
!
!     N : (input) Number of points in the required quadrature formula
!     ALP : (input) Exponent of 1-x in weight function
!     BETA : (input) Exponent of 1+x in weight function
!		w(x)=(1-x)**ALP*(1+x)**BETA
!     W : (output) Array of length N, which will contain the weights
!     A : (output) Array of length N containing the abscissas
!     IER : (output) Error parameter, IER=0 implies successful execution
!		IER=302 implies N.LE.0
!		IER=313 implies ALP.LE.-1 or BETA.LE.-1
!		IER=321 implies that some coefficient becomes imaginary
!			during calculations.
!		In both these cases calculations are abandoned.
!		Other values may be set by TQL2
!     WK : Real array of length N*(N+2)+3*(N+1) used as scratch space
!
!     Required routines : GAUSRC, TQL2, GAMMA
 
      SUBROUTINE GAUJAC(N,ALP,BETA,W,A,IER,WK)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION A(N),W(N),WK(3,*)
 
      IF(ALP.LE.-1.0.OR.BETA.LE.-1.0) THEN
        IER=313
        RETURN
      ENDIF

      RI0=2.D0**(ALP+BETA+1)*GAMMA(ALP+1.0)*GAMMA(BETA+1.0)/
     1       GAMMA(ALP+BETA+2.0)

!     Coefficients of recurrence relation
      WK(1,1)=1.D0+(ALP+BETA)/2.0
      WK(2,1)=(ALP-BETA)/2.0
!	WK(3,1) is not required and would give 0/0 form in some cases
!	Hence do not calculate it.

      DO 2500 I=2,N+1
        A1=2.0*I*(I+ALP+BETA)*(2*I-2+ALP+BETA)
        WK(1,I)=(2*I-2+ALP+BETA)*(ALP+BETA+2*I-1)*(ALP+BETA+2*I)/A1
        WK(2,I)=(2*I-1.0+ALP+BETA)*(ALP*ALP-BETA*BETA)/A1
        WK(3,I)=2.0*(I-1.D0+ALP)*(I-1.0+BETA)*(ALP+BETA+2*I)/A1
2500  CONTINUE
 
      CALL GAUSRC(N,W,A,WK,RI0,IER,WK(1,N+2))
      END
