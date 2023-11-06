!	To calculate coefficients for B-spline interpolation in 2 dimensions
!
!	NX, NY : (input) Number of entries in the table along X, Y directions
!	X, Y : (input) Array of length NX, NY containing the abscissas
!	F : (input) Array of length LA*NY containing the function values
!	K : (input) Order of B-spline required. K=4 gives cubic B-splines
!	AX : (input/output) Real array of length LA*3K containing the
!		triangular decomposition of equation matrix in band form
!		for interpolation along X.
!		For IFLG.GE.2, this array must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	AY : (input/output) Real array of length LA*3K containing the
!		triangular decomposition of equation matrix in band form
!		for interpolation along Y.
!		For IFLG.GE.2, this array must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	LA : (input) the first dimension of arrays AX, AY, F, C as specified
!			 in calling program, LA .GE. MAX(NX,NY)
!	C : (output) Real array of length LA*NY containing coefficients
!		of expansion, which will be calculated by the subroutine
!	XF : (input/output) Real array of size MX, containing
!		the knots used for B-spline calculations along x.
!		The knots must be distinct and in ascending order.
!		For IFLG>1, this array must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	YF : (input/output) Real array of size MY, containing
!		the knots used for B-spline calculations along y.
!		The knots must be distinct and in ascending order.
!		For IFLG>1, this array must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	MX : (input/output) Number of knots for B-splines along X
!		For IFLG>1, this value must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	MY : (input/output) Number of knots for B-splines along Y
!		For IFLG>1, this value must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	IFLG : (input/output) Integer specifying the type of calculation required
!		IFLG<2 The matrix will be calculated and solved for coefficients
!		IFLG>1 The triangular decomposition of matrix is assumed
!			to be available in AX, AY and coefficients C are calculated
!	INTX, INTY : (input/output) Integer arrays of length NX, NY
!		containing information about
!		pivoting during solution of system of linear equations
!		For IFLG>1, these arrays must be supplied, for other values
!		of IFLG they are calculated by the subroutine
!	WK : Scratch array of length NX*LA+NX+NY
!	IER : (output) Error parameter, IER=0 implies successful execution
!		nonzero values may be set by BSPLIN, BSPINT or GAUBND
!
!	Required routines : BSPINT, BSPLIN, GAUBND

      SUBROUTINE BSPINT2(NX,NY,X,Y,F,K,AX,AY,LA,C,XF,YF,MX,MY,IFLG,
     1  INTX,INTY,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NX),Y(NY),AX(LA,3*K),AY(LA,3*K),XF(MX),YF(MY),C(LA,NY)
     1  ,F(LA,NY),WK(NX*LA+NX+NY),INTX(NX),INTY(NY)
 
      IF(IFLG.LE.1) THEN
!	Setup the matrices for interpolation along each direction
        IFLG1=1
        CALL BSPINT(NX,X,F,K,AX,LA,C,XF,MX,IFLG1,INTX,WK,IER)
        IF(IER.GT.100) RETURN
        IFLG1=1
        CALL BSPINT(NY,Y,F,K,AY,LA,C,YF,MY,IFLG1,INTY,WK,IER)
        IF(IER.GT.100) RETURN
        IFLG=2
      ENDIF
 
!	Calculate the interpolation along Y
      DO 2500 I=1,NX
        DO 2500 J=1,NY
          WK(J+(I-1)*LA)=F(I,J)
2500  CONTINUE
      NUM=NX
      IFLG1=2
      KB=K-1
      LN1=LA*NX+1
      CALL GAUBND(NY,KB,NUM,AY,WK,DET,IDET,INTY,LA,IER,IFLG1,WK(LN1))
      IF(IER.GT.100) RETURN
 
!	Calculate the interpolation along X
      DO 3000 J=1,NY
        DO 3000 I=1,NX
          C(I,J)=WK(J+(I-1)*LA)
3000  CONTINUE
      NUM=NY
      IFLG1=2
      CALL GAUBND(NX,KB,NUM,AX,C,DET,IDET,INTX,LA,IER,IFLG1,WK(LN1))
 
      END
