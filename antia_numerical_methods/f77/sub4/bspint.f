!	To calculate coefficients for B-spline interpolation
!
!	N : (input) Number of entries in the table
!	X : (input) Array of length N containing the abscissas
!	F : (input) Array of length N containing the function values
!		F(I) is the tabulated function value at X(I).
!	K : (input) Order of B-spline required. K=4 gives cubic B-splines
!	A : (input/output) Real array of length LA*3K containing the
!		triangular decomposition of equation matrix in band form
!		For IFLG=2, this array must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	LA : (input) The first dimension of A as specified in calling program
!		LA.GE.N
!	C : (output) Coefficients of expansion, which will be calculated
!		provided IFLG.NE.1
!	XF : (input/output) Real array of size NO, containing
!		the knots used for B-spline calculations.
!		The knots must be distinct and in ascending order.
!		For IFLG=2, this array must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	NO : (input/output) Number of knots for B-splines
!		For IFLG=2, this number must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	IFLG : (input/output) Integer specifying the type of calculation required
!		IFLG=0 The matrix will be calculated and solved for coefficients
!		IFLG=1 The matrix will be calculated and triangular decomposition
!			is obtained, but coefficients are not calculated
!		IFLG=2 The triangular decomposition of matrix is assumed
!			to be available in A and coefficients C are calculated
!		IFLG=-1 same as 0, except that no pivoting will be used
!	INC : (input/output) Integer array containing information about
!		pivoting during solution of system of linear equations
!		For IFLG=2, this array must be supplied, for other values
!		of IFLG it is calculated by the subroutine
!	WK : Scratch array of length 3*N+K+7
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=204 implies N<K or K<2
!		other values may be set by BSPLIN or GAUBND
!
!	Required routines : BSPLIN, GAUBND

      SUBROUTINE BSPINT(N,X,F,K,A,LA,C,XF,NO,IFLG,INC,WK,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),A(LA,3*K),WK(3*N+K+7),C(N),F(N),XF(NO),INC(N)
 
      IF(N.LT.K.OR.K.LT.2) THEN
        IER=204
        RETURN
      ENDIF
 
      IF(IFLG.LE.1) THEN
!	set up the knots for B-splines by dropping points near the ends
        XF(1)=X(1)
        KL=(K-2)/2
        KU=(K-1)/2
        DO 2000 I=2+KL,N-1-KU
          XF(I-KL)=X(I)
2000    CONTINUE
        XF(N-KL-KU)=X(N)
        NO=N-KL-KU
        NDB=N+2
        NDERIV=0
 
!	Set up the equation matrix for calculating coefficients of expansion
!	The matrix is in band form A_{i,j} is stored in A(I,J-I+K)
        DO 2500 I=1,N
          XB=X(I)
          CALL BSPLIN(XF,NO,K,XB,NDERIV,WK,WK(NDB),WK(NDB+2),LEFT,IER,
     1                WK(2*NDB+2))
          IF(IER.GT.100) RETURN
          DO 2200 J=MAX(1,I-K+1),MIN(N,I+K-1)
            A(I,J-I+K)=WK(J)
2200      CONTINUE
2500    CONTINUE
      ENDIF

!	Solve the system of equations for a band matrix
      NUM=1
      KB=K-1
      DO 3000 I=1,N
        C(I)=F(I)
3000  CONTINUE
      CALL GAUBND(N,KB,NUM,A,C,DET,IDET,INC,LA,IER,IFLG,WK)
      END
