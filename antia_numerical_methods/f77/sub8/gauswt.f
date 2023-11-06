!	To calculate weights and abscissas of a quadrature formula with
!	specified weight function.
!
!	N : (input) Number of points in the required quadrature formula
!	W : (output) Array of length N, which will contain the weights
!	AB : (input/output) Array of length N containing the abscissas
!		For Gaussian formulas (QGAUS=.TRUE.) AB will be calculated, 
!		while for QGAUS=.FALSE. abscissas must be supplied
!	FMOM : (input) Name of the function routine to calculate the moments
!		FUNCTION FMOM(I) should calculate integral of w(x)x**I
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=303 implies N.LE.0 or N.GE.NPMAX
!		IER=322 implies GAUELM failed to find coefficients of polynomial
!		IER=323 implies POLYR failed to find roots of polynomial
!		IER=324 implies GAUELM failed to find weights 
!	QGAUS : (input) Logical parameter to decide type of formula to be obtained
!		If QGAUS=.TRUE. a Gaussian formula is calculated. In this
!		case both abscissas and weights are calculated.
!		If QGAUS=.FALSE. an interpolatory formula is calculated.
!		In this case only weights are calculated, while abscissas
!		must be supplied.
!
!	Required routines : GAUELM, POLYR, LAGITR, FMOM

      SUBROUTINE GAUSWT(N,W,AB,FMOM,IER,QGAUS)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER (NPMAX=65)
      EXTERNAL FMOM
      DIMENSION A(NPMAX,NPMAX),INC(NPMAX),COF(NPMAX),W(N),AB(N)
      COMPLEX*16  ZERO(NPMAX)

      IER=303
      IF(N.LE.0.OR.N.GE.NPMAX) RETURN
      LJ=NPMAX
      IER=0
      IF(QGAUS) THEN

!	Calculating the coefficients of the orthogonal polynomial
        DO 2000 I=1,N
          DO 1500 J=1,N
1500      A(I,J)=FMOM(I+J-2)
          COF(I)=-FMOM(N+I-1)
2000    CONTINUE
        IFLG=0
        CALL GAUELM(N,1,A,COF,DET,INC,LJ,IER1,IFLG)
        IF(IER1.GT.100) THEN
          IER=322
          RETURN
        ENDIF

!	Find the roots of polynomial, which will be the abscissas
        COF(N+1)=1.
        QREFIN=.TRUE.
!	Array A is used as scratch space
        CALL POLYR(N,COF,ZERO,IER1,QREFIN,A)
        IF(IER1.GT.100) THEN
          IER=323
          RETURN
        ENDIF
        DO 2800 I=1,N
2800    AB(I)=ZERO(I)

      ENDIF

!	Calculate the weights
      DO 4000 I=1,N
        DO 3000 J=1,N
          IF(I.EQ.1) A(I,J)=1.0
          IF(I.GT.1) A(I,J)=AB(J)**(I-1)
3000    CONTINUE
4000  W(I)=FMOM(I-1)
      IFLG=0
      CALL GAUELM(N,1,A,W,DET,INC,LJ,IER1,IFLG)
      IF(IER1.GT.100) IER=324

      END
