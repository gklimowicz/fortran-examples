!	To calculate the integral of a function defined by B-spline expansion
!
!	N : (input) Number of knots for B-spline expansion
!	X : (input) Real array of length N containing the knots.
!		The knots must be distinct and in ascending order.
!       K : (input) Order of B-splines, K=4 for cubic B-splines
!	WT : (input) Real array of length N+K-2 containing the
!		coefficients of B-spline expansion
!	XL : (input) Lower limit of integration
!	XU : (input) Upper limit of integration
!	WK : Scratch array of length 5N+6K+9
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=31 implies that XL is outside the range of table
!		IER=32 implies that XU is outside the range of table
!
!	BSPQD = Integral of \sum WT(I)\phi_I(x) over [XL,XU]
!
!	Required routines : BSPLIN

      FUNCTION BSPQD(N,X,K,WT,XL,XU,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),WT(N+K-2),WK(5*N+6*K+9)
 
      IER=0
      BSPQD=0.0
      K1=K+1
      IF(XL.EQ.XU) RETURN
 
!	Calculate the B-spline basis of order K+1 at XL and XU
      NDERIV=0
      NK=N+K1
      CALL BSPLIN(X,N,K1,XL,NDERIV,WK,WK(NK),WK(2*NK),LF1,IER,WK(3*NK))
      IF(IER.GT.100) RETURN
      CALL BSPLIN(X,N,K1,XU,NDERIV,WK(NK),WK(2*NK),WK(3*NK),LF2,IER,
     1            WK(4*NK))
      IF(IER.GT.100) RETURN
      IF(XL.LT.X(1).OR.XL.GT.X(N)) IER=31
      IF(XU.LT.X(1).OR.XU.GT.X(N)) IER=32

!	Define the extra knots
      K2=2*NK+K
      DO 1200 I=1,N
        WK(K2+I)=X(I)
1200  CONTINUE
      DO 1400 I=1,K
        WK(K2+1-I)=X(1)
        WK(K2+N+I)=X(N)
1400  CONTINUE
 
!	The sum for x=XL
      F1=0.0
      N1=N+K
      DO 2000 I=LF1,LF1+K
        F=0
        DO 1500 J=2,I
          F=F+WT(J-1)*(WK(J+K2)-WK(J-K+K2))
1500    CONTINUE
        F1=F1+F*WK(I)/K
2000  CONTINUE
 
!	The sum for x=XU
      F2=0.0
      DO 2500 I=LF2,LF2+K
        F=0
        DO 2200 J=2,I
          F=F+WT(J-1)*(WK(J+K2)-WK(J-K+K2))
2200    CONTINUE
        F2=F2+F*WK(N1+I)/K
2500  CONTINUE

      BSPQD=F2-F1
      END
