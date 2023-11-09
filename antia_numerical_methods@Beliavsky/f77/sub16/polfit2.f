!	Least squares polynomial fit using orthogonal polynomials in 2 dimensions
!	Weights are assumed to be equal for all points
!
!	NX : (input) Number of data points along X direction
!	NY : (input) Number of data points along Y direction
!	X,Y : (input) Real arrays of length NX,NY containing the coordinates
!		of tabular points
!	F : (input) Real array of length LA*NY containing the function values
!		F(I,J) is the value at X(I),Y(J)
!	AX : (output) Real array of length IC*3 containing information about
!		fit along X direction. AX(I,1), AX(I,2), AX(I,3) will
!		respectively contain the coefficients, alpha, beta, gamma
!	AY : (output) Real array of length IC*3 containing information about
!		fit along Y direction. AY(I,1), AY(I,2), AY(I,3) will
!		respectively contain the coefficients, alpha, beta, gamma
!	LA : (input) First dimension of arrays  F, FY as declared
!		in the calling program. LA .GE. MAX(NX,NY)
!	C : (output) Real array of length IC*(MY+1) containing the fitted
!		coefficients of product of orthogonal polynomials in X & Y
!	IC : (input) First dimension of arrays C, AX, AY as declared in the calling
!		program. IC > MAX(MX,MY)+1
!	MX : (input) Required degree of polynomial in X
!	MY : (input) Required degree of polynomial in Y
!	FY : (output) Real array of length LA*NY containing the values of
!		fitted function at each tabular point	
!	WK : Real array of length MAX[ NX*(NY+MY+1), 6(MAX(MX,MY)+2)]
!		 used as scratch space
!	AW : Real array of length LA*3 used as scratch space
!	CHISQ : (output) the Chi square value for the fit
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=602 implies that IC<MX+1 or IC<MY+1 or LA<NX or LA<NY
!		IER=603 implies that NX<MX+1 or MX<0 or NY < MY+1 or MY<0
!		In both these cases no calculations are done
!		Other values may be set by POLFIT1.
!	
!	The fitted polynomial can be calculated at any value of x using POLEV2
!
!	Required routines : POLFIT1, POLEV2, POLORT
 
      SUBROUTINE POLFIT2(NX,NY,X,Y,F,AX,AY,LA,C,IC,MX,MY,FY,WK,AW,
     1         CHISQ,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(NX),Y(NY),AX(IC,3),AY(IC,3),C(IC,MY+1),F(LA,NY),
     1         FY(LA,NY),WK(*),AW(LA,3)
 
      IF(MAX(NX,NY).GT.LA.OR.IC.LT.MX+1.OR.IC.LT.MY+1) THEN
        IER=602
        RETURN
      ENDIF
 
      IF(MX.LT.0.OR.MY.LT.0.OR.NX.LE.MX.OR.NY.LE.MY) THEN
        IER=603
        RETURN
      ENDIF
 
!     Set the weights to 1
      DO 1000 I=1,MAX(NX,NY)
        AW(I,1)=1.0
1000  CONTINUE
 
!     Set up the RHS for calculating the fits along y-axis
      LJ=NY
      DO 2500 I=1,NX
        DO 2500 J=1,NY
          WK(J+(I-1)*LJ)=F(I,J)
2500  CONTINUE
 
      LN1=LJ*NX+1
      M=MY+1
      NUM=NX
      CALL POLFIT1(NY,MY,NUM,Y,WK,AW(1,1),WK(LN1),AY,AY(1,2),
     1      AY(1,3),AW(1,2),IER)
      IF(IER.GT.100) RETURN
 
!     Set up the RHS for calculating the fits along x-axis
      DO 3000 J=1,M
        DO 3000 I=1,NX
          WK(I+(J-1)*NX)=WK(LN1+J-1+(I-1)*M)
3000  CONTINUE
      NUM=M
      LN1=M*NX+1
      CALL POLFIT1(NX,MX,NUM,X,WK,AW(1,1),WK(LN1),AX,AX(1,2),
     1      AX(1,3),AW(1,2),IER)
      IF(IER.GT.100) RETURN
 
!	Store the calculated coefficients in array C
      M=MX+1
      DO 3500 I=1,MY+1
        DO 3500 J=1,MX+1
          C(J,I)=WK(LN1+J-1+(I-1)*M)
3500  CONTINUE
 
!     Calculate the CHI square
      CHISQ=0.0
      DO 4000 I=1,NY
        DO 4000 J=1,NX
          CALL POLEV2(MX,MY,AX,AY,IC,C,X(J),Y(I),FY(J,I),DFX,DFY,
     1                DFXX,DFXY,DFYY,WK,IER1)
          CHISQ=CHISQ+(F(J,I)-FY(J,I))**2
4000  CONTINUE
 
      END
