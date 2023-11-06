!	To calculate the integral of a function defined by B-spline expansion
!	over 2 dimensions
!
!	NX : (input) Number of knots for B-spline expansion along x direction
!	NY : (input) Number of knots for B-spline expansion along y direction
!	X : (input) Real array of length NX containing the knots along x direction
!		The knots must be distinct and in ascending order.
!	Y : (input) Real array of length NX containing the knots along y direction
!		The knots must be distinct and in ascending order.
!       K : (input) Order of B-splines, K=4 for cubic B-splines
!	WT : (input) Array of length IW*(NY+K-2) containing coefficients
!		of B-spline expansion. 
!	IW : (input) The first dimension of WT as declared in calling program
!	XL : (input) Lower limit of integration along x
!	XU : (input) Upper limit of integration along x
!	YL : (input) Lower limit of integration along y
!	YU : (input) Upper limit of integration along y
!	WK : Scratch array of length 6*MAX(NX,NY)+6*K+9
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values may be set by BSPQD which is called
!
!	BSPQD2 = Integral of \sum WT(I,J)\phi_I(x)\psi_J(y) over [XL,XU] x [YL,YU]
!
!	Required routines : BSPLIN, BSPQD
!

      FUNCTION BSPQD2(NX,NY,X,Y,K,WT,IW,XL,XU,YL,YU,WK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(NX),Y(NY),WT(IW,NY+K-2),WK(*)
 
 
      N1=NY+K
      BSPQD2=0.0
      DO 2000 I=1,NY+K-2
!	Calculate integral along x
        WK(I)= BSPQD(NX,X,K,WT(1,I),XL,XU,WK(N1),IER)
        IF(IER.GT.100) RETURN
2000  CONTINUE
 
!	Calculate integral along y
      BSPQD2=BSPQD(NY,Y,K,WK,YL,YU,WK(N1),IER)
 
      END
 
