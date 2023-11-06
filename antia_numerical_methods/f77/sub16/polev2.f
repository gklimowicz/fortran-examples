!	Evaluating the fitted polynomial and its derivatives at any value
!	of x using known coefficients of orthogonal polynomials in 2 dimensions
!	Should be used to evaluate the polynomial using coefficients calculated
!	by POLFIT2.
!
!	NX : (input) Degree of polynomial in X
!	NY : (input) Degree of polynomial in Y
!	AX : (input) Real array of length LA*2 containing the coefficients
!		alpha and beta for orthogonal polynomials in X
!		AX(I,1) contains alpha and AX(I,2) contains beta
!	AY : (input) Real array of length LA*2 containing the coefficients
!		alpha and beta for orthogonal polynomials in Y
!		AY(I,1) contains alpha and AY(I,2) contains beta
!		The arrays AX and AY can be calculated using POLFIT2
!	LA : (input) First dimension of arrays AX, AY and WT in the calling
!		program. LA > MAX(NX,NY)
!	WT : (input) Real array of length LA*(MY+1) containing the coefficients
!		of the fit. WT(I,J) is the coefficient of
!		PHI_I(X)PSI_J(Y), where PHI_I and PSI_J are orthogonal
!		polynomials in X and Y
!	X0,Y0 : (input) Coordinates of the point at which polynomial
!		needs to be evaluated
!	F : (output) Calculated value of the fitted polynomial at (X0,Y0)
!	DFX : (output) First derivative  dF/dX at X0,Y0
!	DFY : (output) First derivative dF/dY at X0,Y0
!	DFXX : (output) Second derivative d^2F/dXdX at X0,Y0
!	DFXY : (output) Second derivative d^2F/dXdY at X0,Y0
!	DFYY : (output) Second derivative d^2F/dYdY at X0,Y0
!	WK : Real array of length 6(MAX(NX,NY)+2) used as scratch space
!	IER : (output) error parameter, IER=0 implies successful execution
!		IER=604 implies that LA.LE.MAX(NX,NY), in which case no
!			calculations are done.
!	
!	Required routines : POLORT
!
      SUBROUTINE POLEV2(NX,NY,AX,AY,LA,WT,X0,Y0,F,DFX,DFY,
     1      DFXX,DFXY,DFYY,WK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION AX(LA,2),AY(LA,2),WT(LA,NY+1),WK(*)
 
      NK=MAX(NX,NY)+2
      IF(NK-2.GE.LA) THEN
        IER=604
        RETURN
      ENDIF
!	Calculate the orthogonal polynomials along each dimension
      CALL POLORT(NX,AX(1,1),AX(1,2),X0,WK,WK(NK),WK(2*NK))
      CALL POLORT(NY,AY(1,1),AY(1,2),Y0,WK(3*NK),WK(4*NK),WK(5*NK))
 
!	Calculate the fitted polynomial and its derivatives
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
      DO 2000 I=1,NX+1
        DO 2000 J=1,NY+1
          F=F+WT(I,J)*WK(I)*WK(N3+J)
          DFX=DFX+WT(I,J)*WK(N1+I)*WK(N3+J)
          DFY=DFY+WT(I,J)*WK(I)*WK(N4+J)
          DFXX=DFXX+WT(I,J)*WK(N2+I)*WK(N3+J)
          DFXY=DFXY+WT(I,J)*WK(N1+I)*WK(N4+J)
          DFYY=DFYY+WT(I,J)*WK(I)*WK(N5+J)
2000  CONTINUE
      END
