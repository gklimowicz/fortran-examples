!	To evaluate a B-spline expansion in 2 dimensions
!
!	NX : (input) Number of knots to define B-splines along 1st dimension
!	NY : (input) Number of knots to define B-splines along 2nd dimension
!	X : (input) Real array of length NX containing the knots.
!		The knots must be distinct and in ascending order.
!	Y : (input) Real array of length NY containing the knots.
!		The knots must be distinct and in ascending order.
!	K : (input) Order of B-splines, K=4 for cubic B-splines
!	NDERIV : (input) Number of derivatives required
!		For NDERIV.LE.0 only function value is calculated
!		For NDERIV=1 first derivative is also calculated
!		For NDERIV>1 both first and second derivatives are calculated
!	WT : (input) Array of length IW*(NY+K-2) containing coefficients
!		of B-spline expansion,
!	IW : (input) First dimension of array WT as defined in the calling program
!	X0,Y0 : (input) The point at which expansion has to be evaluated
!	DFX : (output) First derivative of function w.r.t. X at X0, Y0
!	DFY : (output) First derivative of function w.r.t. Y at X0, Y0
!	DFXX : (output) Second derivative of function w.r.t X,X at X0, Y0
!	DFXY : (output) Second derivative of function w.r.t X,Y at X0, Y0
!	DFYY : (output) Second derivative of function w.r.t Y,Y at X0, Y0
!	WK : Scratch array of length 7*MAX(NX,NY)+8K+2
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values of IER may be set by BSPLIN which is called
!
!	BSPEV2 =SUM_{i=1}^{NX+K-2} SUM_{j=1}^{NY+K-2} WT(i,j)\phi_i(X0)\psi_j(Y0)
!	where \phi_i(x) are B-spline basis functions on knots X
!	and \psi_j(y) are B-spline basis functions on knots Y
!
!	Required routines : BSPLIN

      FUNCTION BSPEV2(NX,NY,X,Y,K,NDERIV,WT,IW,X0,Y0,DFX,DFY,
     1          DFXX,DFXY,DFYY,WK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(NX),Y(NY),WT(IW,NY+K-2),WK(*)
 
      BSPEV2=0.0
      NK=MAX(NX,NY)+K
      CALL BSPLIN(X,NX,K,X0,NDERIV,WK,WK(NK),WK(2*NK),LX,IER,WK(3*NK))
      IF(IER.GT.100) RETURN
      CALL BSPLIN(Y,NY,K,Y0,NDERIV,WK(3*NK),WK(4*NK),WK(5*NK),LY,IER,
     1            WK(6*NK))
      IF(IER.GT.100) RETURN
 
      F=0.0
      DFX=0.0
      DFY=0.0
      DFXX=0.0
      DFYY=0.0
      DFXY=0.0
      N1=NK-1
      N2=2*NK-1
      N3=3*NK-1
      N4=4*NK-1
      N5=5*NK-1
      DO 2000 I=LX,LX+K-1
        DO 2000 J=LY,LY+K-1
          F=F+WT(I,J)*WK(I)*WK(N3+J)
          DFX=DFX+WT(I,J)*WK(N1+I)*WK(N3+J)
          DFY=DFY+WT(I,J)*WK(I)*WK(N4+J)
          DFXX=DFXX+WT(I,J)*WK(N2+I)*WK(N3+J)
          DFXY=DFXY+WT(I,J)*WK(N1+I)*WK(N4+J)
          DFYY=DFYY+WT(I,J)*WK(I)*WK(N5+J)
2000  CONTINUE
      BSPEV2=F
      END
