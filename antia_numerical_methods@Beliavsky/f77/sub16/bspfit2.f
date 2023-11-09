!	To calculate linear least squares fit to B-spline basis functions in 2 dimension
!	Weights are assumed to be unity
!
!	NX : (input) Number of data points along x-axis
!	NY : (input) Number of data points along y-axis
!	X,Y : (input) Real arrays of length NX,NY containing the coordinates
!		of points at which function values are available
!	F : (input) Real array of length LA*NY containing the function values
!		F(I,J) should be function value at X(I),Y(J)
!	K : (input) Order of B-splines required, K=4 gives cubic B-splines
!	AX : (output) Real array of length LA*(MX+K-2) containing the matrix
!		U of SVD of the design matrix for fit along x-axis
!	AY : (output) Real array of length LA*(MY+K-2) containing the matrix
!		U of SVD of the design matrix for fit along y-axis
!	LA : (input) First dimension of arrays F, AX, AY, C, FY as declared
!		in the calling program (LA .GE. 2*MAX(NX,NY))
!	VX : (output) Real array of length IV*(MX+K-2) containing the matrix
!		V of SVD of the design matrix for fit along x-axis
!	VY : (output) Real array of length IV*(MY+K-2) containing the matrix
!		V of SVD of the design matrix for fit along y-axis
!	IV : (input) First dimension of VX, VY in the calling program
!		IV .GE. MAX(MX,MY)+K-2
!	SIGMAX : (output) Real array of length MX+K-2 containing the singular
!		values of the design matrix for fit along x-axis
!	SIGMAY : (output) Real array of length MY+K-2 containing the singular
!		values of the design matrix for fit along y-axis
!	C : (output) Real array of length LA*NY containing the fitted coefficients
!		Note that although the number of coefficients is (MX+K-2)*(MY+K-2)
!		the rest of array is used as scratch space
!	XF : (input) Real array of size MX, containing the knots
!		along x-axis used for defining B-spline basis functions.
!		The knots must be distinct and in ascending order.
!	YF : (input) Real array of size MY, containing the knots
!		along y-axis used for defining B-spline basis functions.
!		The knots must be distinct and in ascending order.
!	MX : (input) Number of knots for B-splines along x-axis,
!		the number of basis functions would be MX+K-2
!	MY : (input) Number of knots for B-splines along y-axis,
!		the number of basis functions would be MY+K-2
!	FY : (output) Real array of length LA*NY containing the values of fitted
!		function at each of the tabular points
!	WK : Real array of length LA*NX+NX+NY used as scratch space
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
!		No calculations are done in all these cases
!		Other values of IER may be set by SVD or BSPLIN
!
!	Required routines : BSPFIT, BSPLIN, BSPEVL, BSPEV2, SVD, SVDEVL
!
      SUBROUTINE BSPFIT2(NX,NY,X,Y,F,K,AX,AY,LA,VX,VY,IV,SIGMAX,
     1        SIGMAY,C,XF,YF,MX,MY,FY,WK,REPS,RLM,IDE,CHISQ,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(NX),Y(NY),AX(LA,MX+K-2),AY(LA,MY+K-2),XF(MX),YF(MY),
     1        C(LA,NY),F(LA,NY),VX(IV,MX+K-2),VY(IV,MY+K-2),
     2        SIGMAX(MX+K-2),SIGMAY(MY+K-2),FY(LA,NY),WK(NX*LA+NX+NY)
 
      N1=4*MAX(NX,NY)+K+12
      N2=N1+MAX(NX,NY)+2
!	Set the weights to 1
      DO 1000 I=0,MAX(NX,NY)
        WK(N1+I)=1.0
1000  CONTINUE

!	Calculate the SVD of matrix for fit along x-axis
      IFLG1=1
      CALL BSPFIT(NX,X,F,WK(N1),K,AX,LA,VX,IV,SIGMAX,C,XF,MX,FY,
     1        IFLG1,WK,REPS,RLM,IDE,CHISQ,WK(N2),IER)
      IF(IER.GT.100) RETURN

!	Calculate the SVD of matrix for fit along y-axis
      IFLG1=1
      CALL BSPFIT(NY,Y,F,WK(N1),K,AY,LA,VY,IV,SIGMAY,C,YF,MY,FY,
     1        IFLG1,WK,REPS,RLM,IDE,CHISQ,WK(N2),IER)
      IF(IER.GT.100) RETURN
 
!	Set up the RHS for calculating the fits along y-axis
      DO 2500 I=1,NX
        DO 2500 J=1,NY
          WK(J+(I-1)*LA)=F(I,J)
          IF(RLM.GT.0.0) WK(J+NY+(I-1)*LA)=0.0
2500  CONTINUE

!	N1 is the number of equations for fit along y-axis
      N1=NY
      IF(RLM.GT.0.0) N1=2*NY
      LN1=LA*NX+1
      M=MY+K-2
      DO 2800 I=1,NX
        NI=1+(I-1)*LA
        CALL SVDEVL(M,N1,AY,VY,SIGMAY,LA,IV,WK(NI),WK(LN1),REPS)
2800  CONTINUE
 
!	Set up the RHS for calculating the fits along x-axis
      DO 3000 J=1,M
        DO 3000 I=1,NX
          C(I,J)=WK(J+(I-1)*LA)
          IF(RLM.GT.0.0) C(I+NX,J)=0.0
3000  CONTINUE
      M1=MX+K-2
      N1=NX
      IF(RLM.GT.0.0) N1=2*NX
      DO 3200 I=1,M
        CALL SVDEVL(M1,N1,AX,VX,SIGMAX,LA,IV,C(1,I),WK,REPS)
3200  CONTINUE
 
!	Calculate the CHI square
      CHISQ=0.0
      NDERIV=0
      DO 4000 I=1,NY
        DO 4000 J=1,NX
          FY(J,I)=BSPEV2(MX,MY,XF,YF,K,NDERIV,C,LA,X(J),Y(I),DFX,DFY,
     1            DFXX,DFXY,DFYY,WK,IER1)
          CHISQ=CHISQ+(F(J,I)-FY(J,I))**2
4000  CONTINUE
 
      END
