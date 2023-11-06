!	To calculate function value using B-spline expansion
!
!	N : (input) Number of knots to define B-splines
!	X : (input) Real array of length N+2K+1 containing the knots.
!		The knots must be distinct and in ascending order.
!	K : (input) Order of B-splines, K=4 for cubic B-splines
!	NDERIV : (input) Number of derivatives required
!		For NDERIV.LE.0 only function value is calculated
!		For NDERIV=1 first derivative is also calculated
!		For NDERIV>1 both first and second derivatives are calculated
!	WT : (input) Coefficients of B-spline expansion
!	X0 : (input) The point at which expansion has to be evaluated
!	DF : (output) First derivative of function at X0
!	DDF : (output) Second derivative of function at X0
!	WK : Scratch array of length 4N+5K+2
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values of IER may be set by BSPLIN which is called
!
!	BSPEVL = SUM_{i=1}^{N+K-2} WT(I) \phi_i(X0)
!	where \phi_i(x) are B-spline basis functions on knots X
!
!	Required routines : BSPLIN

      FUNCTION BSPEVL(N,X,K,NDERIV,WT,X0,DF,DDF,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N+K),WT(N+K-2),WK(4*N+5*K+2)
 
      BSPEVL=0.0
      NK=(N+K)
      CALL BSPLIN(X,N,K,X0,NDERIV,WK,WK(NK),WK(2*NK),LEFT,IER,WK(3*NK))
      IF(IER.GT.100) RETURN
 
      F=0.0
      DF=0.0
      DDF=0.0
      N1=N+K-1
      N2=2*(N+K)-1
      DO 2000 I=LEFT,LEFT+K-1
        F=F+WT(I)*WK(I)
        DF=DF+WT(I)*WK(N1+I)
        DDF=DDF+WT(I)*WK(N2+I)
2000  CONTINUE
      BSPEVL=F
      END
