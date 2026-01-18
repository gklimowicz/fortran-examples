!	To calculate linear least squares fit to B-spline basis functions in 2 dimension
!	version for general weights, but much slower than BSPFIT2
!
!	NX : (input) Number of data points along x-axis
!	NY : (input) Number of data points along y-axis
!	X,Y : (input) Real array of length NX,NY containing the coordinates
!		of points at which function values is available
!	F : (input) Real array of length IC*NY containing the function values
!		F(I,J) should be function value at X(I),Y(J)
!	EF : (input) Real array of length IC*NY containing the error estimates in F(I,J)
!	K : (input) Order of B-splines required, K=4 gives cubic B-splines
!	A : (output) Real array of length LA*(MX+K-2)*(MY+K-2) containing
!		the matrix U of SVD of the design matrix
!	LA : (input) First dimension of array A as declared
!		in the calling program (LA .GE. 3*NX*NY)
!	V : (output) Real array of length IV*(MX+K-2)*(MY+K-2) containing
!		the matrix V of SVD of the design matrix
!	IV : (input) First dimension of V in the calling program
!		IV .GE. (MX+K-2)*(MY+K-2)
!	SIGMA : (output) Real array of length (MX+K-2)*(MY+K-2)
!		containing the singular values of the design matrix
!	C : (output) Real array of length IC*(MY+K-2) containing the
!		fitted coefficients
!	IC : (input) First dimension of arrays C, F, EF, FY as declared
!		in the calling program (IC .GE. NX)
!	XF : (input) Real array of size MX, containing the knots
!		along x-axis used for defining B-spline basis functions.
!	YF : (input) Real array of size MY, containing the knots
!		along y-axis used for defining B-spline basis functions.
!	MX : (input) Number of knots for B-splines along x-axis,
!		the number of basis functions would be MX+K-2
!	MY : (input) Number of knots for B-splines along y-axis,
!		the number of basis functions would be MY+K-2
!	FY : (output) Real array of length IC*NY containing the values of fitted
!		function at each of the tabular points
!	WK : Real array of length 3*NX*NY+(MX+K)*(MY+K) used as scratch space
!	REPS : (input) Required accuracy for solution of equations using SVD
!		singular values less than REPS times maximum will be set to zero
!	RLM : (input) Parameter lambda for smoothing. If RLM.LE.0 no smoothing
!		is applied
!	IDE : (input) Order of derivative to be used for smoothing
!		This is used only when RLM>0. IDE=1 for first derivative
!		and IDE=2 for second derivative smoothing
!	CHISQ : (output) The value of Chi square at minimum
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=608 implies that MX+K-2>NX, MY+K-2>NY or K<2
!		IER=609 implies that RLM>0 and IDE is not acceptable
!		IER=610 implies that EF(I,J).LE.0 for some I,J
!		No calculations are done in all these cases
!		Other values of IER may be set by SVD or BSPLIN
!
!	Required routines : BSPLIN, BSPEV2, SVD, SVDEVL
!
      SUBROUTINE BSPFITW2(NX,NY,X,Y,F,EF,K,A,LA,V,IV,SIGMA,
     1      C,IC,XF,YF,MX,MY,FY,WK,REPS,RLM,IDE,CHISQ,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NX),Y(NY),A(LA,*),XF(MX),YF(MY),C(IC,*),F(IC,NY),
     1      EF(IC,NY),V(IV,*),SIGMA(*),FY(IC,NY),WK(*)
 
      IF(NX.LT.K.OR.NY.LT.K.OR.K.LT.2) THEN
        IER=608
        RETURN
      ENDIF
 
      IF(RLM.GT.0.0.AND.(IDE.LT.1.OR.IDE.GT.2)) THEN
        IER=609
        RETURN
      ENDIF
!     N1 is the number of equations to be solved
      N1=NX*NY
      IF(RLM.GT.0.0) N1=3*N1
 
!     Set up the matrix equation and obtain SVD
!     M is the number of coefficients to be determined
      M1=MX+K-2
      M2=MY+K-2
      M=M1*M2
      NDB=MAX(MX,MY)+5
      NDERIV=0
      IF(RLM.GT.0.0) THEN
        NDERIV=IDE
        NB=NDB-1
        NB1=NDB*4-1
        IF(IDE.EQ.2) NB=2*NDB-1
        IF(IDE.EQ.2) NB1=5*NDB-1
      ENDIF
 
!     Set up the matrix for equations
      DO 2500 I=1,NX
        XB=X(I)
        CALL BSPLIN(XF,MX,K,XB,NDERIV,WK,WK(NDB),WK(NDB*2),LEFT,IER,
     1              WK(NDB*3))
        IF(IER.GT.100) RETURN
        DO 2500 J=1,NY
          YB=Y(J)
          CALL BSPLIN(YF,MY,K,YB,NDERIV,WK(NDB*3),WK(NDB*4),WK(NDB*5),
     1                LEFT,IER,WK(NDB*6))
          IF(IER.GT.100) RETURN
          IF(EF(I,J).LE.0) THEN
            IER=610
            RETURN
          ENDIF
          IJ=I+(J-1)*NX
          DO 2200 K1=1,M1
            DO 2200 K2=1,M2
              IJ1=K1+(K2-1)*M1
              A(IJ,IJ1)=WK(K1)*WK(NDB*3+K2-1)/EF(I,J)
              IF(RLM.GT.0) A(IJ+NX*NY,IJ1)=RLM*WK(NB+K1)*WK(NDB*3+K2-1)
              IF(RLM.GT.0.0) A(IJ+2*NX*NY,IJ1)=RLM*WK(K1)*WK(NB1+K2)
2200      CONTINUE
2500  CONTINUE
      CALL SVD(M,N1,A,V,SIGMA,LA,IV,WK,IER)
      IF(IER.GT.100) RETURN
 
 
 
!     Setup the RHS and solve the equations
      DO 3000 I=1,NX
        DO 3000 J=1,NY
          IJ=I+(J-1)*NX
          WK(IJ)=F(I,J)/EF(I,J)
          IF(RLM.GT.0.0) WK(IJ+NX*NY)=0.0
          IF(RLM.GT.0.0) WK(IJ+2*NX*NY)=0.0
3000  CONTINUE
 
      CALL SVDEVL(M,N1,A,V,SIGMA,LA,IV,WK,WK(N1+2),REPS)
 
      DO 3400 I=1,M1
        DO 3400 J=1,M2
          C(I,J)=WK(I+M1*(J-1))
3400  CONTINUE
 
 
!     Calculate the \chi^2
      CHISQ=0.0
      NDERIV=0
      DO 4000 I=1,NY
        DO 4000 J=1,NX
          FY(J,I)=BSPEV2(MX,MY,XF,YF,K,NDERIV,C,IC,X(J),Y(I),DFX,DFY,
     1                   DFXX,DFXY,DFYY,WK,IER1)
          CHISQ=CHISQ+((F(J,I)-FY(J,I))/EF(J,I))**2
4000  CONTINUE
 
      END
