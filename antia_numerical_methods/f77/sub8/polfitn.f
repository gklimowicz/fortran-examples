!	Least squares polynomial fit using orthogonal polynomials in n dimensions
!	Weights are assumed to be equal for all points and points are
!	assumed to be on a hyper-rectangular mesh
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the number of
!		data points along each direction
!	X : (input) Real array of length LA*N containing the coordinates
!		of tabular points, X(I,J) contains the Ith point along
!		Jth dimension
!	F : (input) Real array of length NK(1)*NK(2)*...*NK(N) containing
!		the function values. The dimension of F in the calling program
!		must match the size along each dimension, F(NK(1),...,NK(N))
!	AX : (output) Real array of length LA*(3*N+3) containing information about
!		fit along each direction. AX(I,3*J-2), AX(I,3*J-1), AX(I,3*J) will
!		respectively contain the coefficients, alpha, beta, gamma
!		for fit along Jth dimension.
!		The rest of the array is used as scratch space
!	LA : (input) First dimension of arrays X and AX as declared
!		in the calling program. LA .GE. MAX(NK(I))
!	C : (output) Real array of length (MK(1)+1)(MK(2)+1)...(MK(N)+1) containing
!		the fitted coefficients of product of orthogonal polynomials 
!		The dimension of F in the calling program must match the size
!		along each dimension, C(MK(1)+1,MK(2)+1,...,MK(N)+1)
!	MK : (input) Integer array of length N containing the required
!		degree of polynomial in each dimension
!	FY : (output) Real array of same size and shape as F containing the
!		fitted function values at each tabular point	
!	WK : Real array of length 2*NK(1)*NK(2)*...*NK(N)  used as scratch space
!	IWK : Integer array of length 2*N used as scratch space
!	CHISQ : (output) the Chi square value for the fit
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=605 implies that LA < MAX(NK(I))
!		In  this case no calculations are done
!		other values may be set by POLFIT1
!	
!	The fitted polynomial can be calculated at any value of x using POLEVN
!
!	Required routines : POLFIT1, POLEVN, POLORT
 
      SUBROUTINE POLFITN(N,NK,X,F,AX,LA,C,MK,FY,WK,IWK,CHISQ,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(LA,N),AX(LA,3*N+3),C(*),F(*),FY(*),WK(*),NK(N),
     1    MK(N),IWK(2*N)
 
      N1=3*N+1
!     Set the weights to one for fits along each dimension
      DO 1000 I=1,LA
        AX(I,N1)=1.0
1000  CONTINUE
 
      NX=1
      NMAX=NK(N)
      DO 2200 I=1,N-1
        NX=NX*NK(I)
        NMAX=MAX(NMAX,NK(I))
2200  CONTINUE
      IF(NMAX.GT.LA) THEN
        IER=605
        RETURN
      ENDIF
 
      N0=NX*NK(N)+1
      NY=1
 
!     Set up the RHS for fit along Nth dimension
      LJ=NK(N)
      DO 2500 I=1,NX
        DO 2500 J=1,NK(N)
          WK(J+(I-1)*LJ)=F(I+(J-1)*NX)
2500  CONTINUE
 
!     Loop for fit along each dimension
      DO 5000 J1=N,1,-1
        NUM=NX*NY
        LJ=NK(J1)
        M=MK(J1)+1
        NI=1+(J1-1)*3
        CALL POLFIT1(LJ,MK(J1),NUM,X(1,J1),WK,AX(1,N1),WK(N0),
     1          AX(1,NI),AX(1,NI+1),AX(1,NI+2),AX(1,N1+1),IER)
        IF(IER.GT.100) RETURN
 
        IF(J1.GT.1) THEN
!     Set up the RHS for fit along next dimension
          NX1=NX/NK(J1-1)
          NY1=NY*M
          LJ1=NK(J1-1)
          DO 3500 I1=1,NY
            DO 3500 I2=1,M
              DO 3500 I=1,NK(J1-1)
                DO 3500 J=1,NX1
                  WK(I+(J-1)*LJ1+(I2-1)*LJ1*NX1+(I1-1)*NX1*LJ1*M)=
     1            WK(N0+I2-1+(J-1)*M+(I-1)*NX1*M+(I1-1)*NX*M)
3500      CONTINUE
          NX=NX1
          NY=NY1
 
        ELSE
!     Store the fitted coefficients in array C
          DO 3800 I=1,M
            DO 3800 J=1,NY
              C(I+(J-1)*M)=WK(N0+I-1+(J-1)*M)
3800      CONTINUE
        ENDIF
5000  CONTINUE
 
!     Calculate the Chi Square
      CHISQ=0.0
      DO 6000 I=1,N
        IWK(I)=1
6000  CONTINUE
      N0=N+1
 
!     Loop over all points
6200  IJ=IWK(1)
      WK(1)=X(IJ,1)
      ID=NK(1)
      DO 6500 I=2,N
        IJ=IJ+ID*(IWK(I)-1)
        ID=ID*NK(I)
        WK(I)=X(IWK(I),I)
6500  CONTINUE
 
      CALL POLEVN(N,MK,AX,LA,C,WK,FY(IJ),WK(N0),IWK(N0))
      CHISQ=CHISQ+(F(IJ)-FY(IJ))**2
 
!     Choose the next point
      J=1
6800  IF(IWK(J).GE.NK(J)) GO TO 7000
      IWK(J)=IWK(J)+1
      GO TO 6200
 
!     If Jth dimension is exhausted go to next one
7000  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 6800
 
      END
