!     LEAST SQUARES FIT USING B-SPLINES IN 2 DIMENSIONS

      PROGRAM BSP2
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION YO(90),XO(90),T(900),B(100,90),B1(90),DB(100,52)
      DIMENSION AX(100,90),WK(90000),F(100,90),XF(90),YF(90)
      DIMENSION AY(100,90),INTY(90),ZO(90),ZF(90),VX(90,90),VY(90,90)
      DIMENSION EF(100,90),Y0(100,90),A2(9000,1000),V2(1000,1000)
 
      FUN(X,Y)=SIN(X)*SIN(Y)+1.D-5*RANGAU(SEED)
      F0(X,Y)=SIN(X)*SIN(Y)
      DF0X(X,Y)=COS(X)*SIN(Y)
      DF0Y(X,Y)=SIN(X)*COS(Y)
      DF0XX(X,Y)=-SIN(X)*SIN(Y)
      DF0XY(X,Y)=COS(X)*COS(Y)
      DF0YY(X,Y)=-SIN(X)*SIN(Y)
 
51    FORMAT('   IER =',I4,3X,'NO. OF POINTS =',I5,3X,'NO. OF KNOTS =',
     1       I3,4X,'ORDER OF B-SPLINE =',I3/5X,
     2       'REGULARISATION PARAMETER =',1PE12.4,5X,
     3       'ORDER OF DERIVATIVE =',I2/5X,'CHI SQUARE =',E12.4,5X,
     4       'IT =',I3)
52    FORMAT('   IER =',I4,3X,'X,Y =',1P2E14.6,5X,'F(X) =',E14.6/
     1       5X,4HF' =,2E14.6/5X,5HF'' =,3E14.6)
53    FORMAT('   EXACT VALUES :   F(X) =',1PE14.6/
     1       5X,4HF' =,2E14.6/5X,5HF'' =,3E14.6)
54    FORMAT('   COEFFICIENTS : '/(1P5E14.6))
 
      NX=6
!	For simplicity assume that no. of points and knots are the same
!	along both axes

100   PRINT *,'TYPE NX=NO. OF POINTS, K=ORDER OF B-SPLINE,'
      PRINT *,'     NO=NO. OF KNOTS, RLM=REGULARISATION PARAMETER,'
      PRINT *,'     IDE=ORDER OF DERIVATIVE FOR REGULARISATION'
      READ *,NX,K,NO,RLM,IDE

      PRINT *,'TYPE IT =1/2 FOR BSPFIT2/BSPFITW2'
      PRINT *,'              (QUITS WHEN IT<1)'
      READ *,IT
      IF(IT.LT.1.OR.IT.GT.2) STOP
 
!   Set the seed so that same random numbers are used every time
      SEED=2
!	Set up the table of values with random error added
      H=1.0D0/(NX-1.)
      NY=NX
      DO 1000 I=1,NX
        XO(I)=(I-1)*H
        YO(I)=(I-1)*H
1000  CONTINUE
      H=1.D0/(NO-1.)
      DO 1200 I=1,NO
        XF(I)=(I-1)*H
        YF(I)=(I-1)*H
1200  CONTINUE
        DO 1500 I=1,NY
          DO 1500 J=1,NX
            F(J,I)=FUN(XO(J),YO(I))
            EF(J,I)=1.D-5
1500  CONTINUE
      IC=4
      LA=100
      L2=9000
      IV2=1000
      IV=90
      IFLG=0
      REPS=1.D-7
      IF(IT.EQ.1) THEN
        CALL BSPFIT2(NX,NY,XO,YO,F,K,AX,AY,LA,VX,VY,IV,T,
     1     B1,B,XF,YF,NO,NO,Y0,WK,REPS,RLM,IDE,CHISQ,IER)
      ELSE
        CALL BSPFITW2(NX,NY,XO,YO,F,EF,K,A2,L2,V2,IV2,T,
     1     B,LA,XF,YF,NO,NO,Y0,WK,REPS,RLM,IDE,CHISQ,IER)
      ENDIF
      WRITE(6,51) IER,NX,NO,K,RLM,IDE,CHISQ,IT
      WRITE(6,54) ((B(I,J),I=1,NO+K-2),J=1,NO+K-2)
      IDE=2
 
2000  PRINT *,'TYPE X, Y VALUE AT WHICH FUNCTION IS REQUIRED'
      PRINT *,'                 (QUITS WHEN XX<-100)'
      READ *,XX,YY
      IF(XX.LT.-100) GO TO 100
      F1= BSPEV2(NO,NO,XF,YF,K,IDE,B,LA,XX,YY,DFX,DFY,DFXX,DFXY,
     1           DFYY,WK,IER)
      WRITE(6,52) IER,XX,YY,F1,DFX,DFY,DFXX,DFXY,DFYY
      WRITE(6,53) F0(XX,YY),DF0X(XX,YY),DF0Y(XX,YY),DF0XX(XX,YY),
     1            DF0XY(XX,YY),DF0YY(XX,YY)
      GO TO 2000
 
 
      END
 
 
!     -------------------------------------------------
 
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
!      IMPLICIT REAL*8(A-H,O-Z)
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
 
!     -------------------------------------------------
 
!	To calculate function value using B-spline expansion
!
!	N : (input) Number of knots to define B-splines
!	X : (input) Real array of length N+2K+1 containing the knots.
!		The knots must be distinct and in ascending order.
!	K : (input) Order of B-splines, K=4 for cubic B-splines
!	NDERIV : (input) Number of derivatives required
!		For NDERIV.LE.0 only function value is calculated
!		For NDERIV=1 first derivative is also calculated
!		For NDERIV>1 both first and second derivatives are calculated
!	WT : (input) Coefficients of B-spline expansion
!	X0 : (input) The point at which expansion has to be evaluated
!	DF : (output) First derivative of function at X0
!	DDF : (output) Second derivative of function at X0
!	WK : Scratch array of length 4N+5K+2
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values of IER may be set by BSPLIN which is called
!
!	BSPEVL = SUM_{i=1}^{N+K-2} WT(I) \phi_i(X0)
!	where \phi_i(x) are B-spline basis functions on knots X
!
!	Required routines : BSPLIN

      FUNCTION BSPEVL(N,X,K,NDERIV,WT,X0,DF,DDF,WK,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N+K),WT(N+K-2),WK(4*N+5*K+2)
 
      BSPEVL=0.0
      NK=(N+K)
      CALL BSPLIN(X,N,K,X0,NDERIV,WK,WK(NK),WK(2*NK),LEFT,IER,WK(3*NK))
      IF(IER.GT.100) RETURN
 
      F=0.0
      DF=0.0
      DDF=0.0
      N1=N+K-1
      N2=2*(N+K)-1
      DO 2000 I=LEFT,LEFT+K-1
        F=F+WT(I)*WK(I)
        DF=DF+WT(I)*WK(N1+I)
        DDF=DDF+WT(I)*WK(N2+I)
2000  CONTINUE
      BSPEVL=F
      END
 
!     --------------------------------------------------------
 
!	To calculate linear least squares fit to B-spline basis functions in
!	1 dimension
!
!	N : (input) Number of data points to be fitted
!	X : (input) Real array of length N containing the coordinates
!		of points at which function values is available
!	F : (input) Real array of length N containing the function values
!		F(I) should be function value at X(I)
!	EF : (input) Real array of length N containing the estimated error
!		in F(I). 
!	K : (input) Order of B-splines required, K=4 gives cubic B-splines
!	A : (output) Real array of length LA*(NO+K-2) containing the matrix
!		U of SVD of the design matrix
!	LA : (input) First dimension of A in the calling program (LA.GE.2N)
!	V : (output) Real array of length IV*(NO+K-2) containing the matrix
!		V of SVD of the design matrix
!	IV : (input) First dimension of V, COV in the calling program (IV.GE.NO+K-2)
!	SIGMA : (output) Real array of length NO+K-2 containing the singular
!		values of the design matrix
!	C : (output) Real array of length 2N containing the fitted coefficients
!		Note that although the number of coefficients is NO+K-2, the
!		rest of array is used as scratch space
!	XF : (input) Real array of size NO, containing
!		the knots used for defining B-spline basis functions.
!		The knots must be distinct and in ascending order.
!	NO : (input) Number of knots for B-splines, the number of basis
!		functions would be NO+K-2
!	Y : (output) Real array of length N containing the values of fitted
!		function at each of the tabular points
!	IFLG : (input/output) Integer specifying the type of calculation required
!		IFLG=0 The matrix will be calculated and solved for coefficients
!			the fitted values Y and CHISQ are also calculated
!		IFLG=1 The matrix will be calculated and SVD
!			is obtained, but coefficients are not calculated
!		IFLG=2 The SVD of matrix is assumed to be available in
!			arrays A, V, SIGMA and coefficients C are calculated
!		IFLG=3 The SVD of matrix is assumed to be available in arrays
!			A, V, SIGMA and coefficients C are calculated and in
!			addition fitted values Y and CHISQ are also calculated
!	WK : Real array of length 4*NO+5*K+2 used as scratch space
!	REPS : (input) Required accuracy for solution of equations using SVD
!		singular values less than REPS times maximum will be set to zero
!	RLM : (input) Parameter lambda for smoothing. If RLM.LE.0 no smoothing
!		is applied
!	IDE : (input) Order of derivative to be used for smoothing
!		This is used only when RLM>0. IDE=1 for first derivative
!		and IDE=2 for second derivative smoothing
!	CHISQ : (output) The value of Chi square at minimum
!       COV : (output) Array of length IV*M containing the covariance
!               matrix of the fitted parameters. COV(I,I) will be the
!               variance in C(I).
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=608 implies that NO+K-2>N or K<2
!		IER=609 implies that RLM>0 and IDE is not acceptable
!		IER=610 implies that EF(I).LE.0 for some I
!		No calculations are done in all these cases
!		Other values of IER may be set by SVD or BSPLIN
!
!       THE ARGUMENTS OF THIS FUNCTION HAVE CHANGED FROM THE EARLIER VERSION.
!	NOW THERE IS AN ADDITIONAL ARGUMENT COV TO CALCULATE THE COVARIANCE
!	MATRIX.
!
!	Required routines : BSPLIN, BSPEVL, SVD, SVDEVL
!
      SUBROUTINE BSPFIT(N,X,F,EF,K,A,LA,V,IV,SIGMA,C,XF,NO,Y,IFLG,WK,
     1       REPS,RLM,IDE,CHISQ,COV,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),A(LA,NO+K-2),WK(4*NO+5*K+2),C(N),F(N),XF(NO),EF(N),
     1   Y(N),SIGMA(NO+K-2),V(IV,NO+K-2),COV(IV,NO+K-2)
 
      IF(N.LT.NO+K-2.OR.K.LT.2) THEN
        IER=608
        RETURN
      ENDIF
 
      IF(RLM.GT.0.0.AND.(IDE.LT.1.OR.IDE.GT.2)) THEN
        IER=609
        RETURN
      ENDIF
!	N1 is the number of equations to be solved
      N1=N
      IF(RLM.GT.0.0) N1=2*N
 
      IF(IFLG.LE.1) THEN
!	Set up the matrix equation and obtain SVD
!	M is the number of coefficients to be determined
        M=NO+K-2
        NDB=M+2
        NDERIV=0
        IF(RLM.GT.0.0) THEN
          NDERIV=IDE
          NB=NDB-1
          IF(IDE.EQ.2) NB=2*NDB-1
        ENDIF
 
!	Set up the matrix for equations
        DO 2500 I=1,N
          XB=X(I)
          IF(EF(I).LE.0) THEN
            IER=610
            RETURN
          ENDIF
          CALL BSPLIN(XF,NO,K,XB,NDERIV,WK,WK(NDB),WK(NDB*2),LEFT,IER,
     1                WK(3*NDB))
          IF(IER.GT.100) RETURN
          DO 2200 J=1,M
            A(I,J)=WK(J)/EF(I)
            IF(RLM.GT.0.0) A(I+N,J)=RLM*WK(NB+J)
2200      CONTINUE
2500    CONTINUE
        CALL SVD(M,N1,A,V,SIGMA,LA,IV,WK,IER)
        IF(IER.GT.100) RETURN
 
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
 
      ENDIF
 
!	Setup the RHS and solve the equations
      DO 3000 I=1,N
        C(I)=F(I)/EF(I)
        IF(RLM.GT.0.0) C(I+N)=0.0
3000  CONTINUE
 
      CALL SVDEVL(M,N1,A,V,SIGMA,LA,IV,C,WK,REPS)
 
      IF(IFLG.EQ.2) RETURN
      IFLG=2
 
!	Calculate the \chi^2
      CHISQ=0.0
      NDERIV=0
      DO 4000 I=1,N
        Y(I)=BSPEVL(NO,XF,K,NDERIV,C,X(I),DF,DDF,WK,IER)
        CHISQ=CHISQ+((F(I)-Y(I))/EF(I))**2
4000  CONTINUE

!       Computing the covariance matrix
      SIGMAX=0
      DO I=1,M
        IF(SIGMA(I).GT.SIGMAX) SIGMAX=SIGMA(I)
      ENDDO
      DO I=1,M
      DO J=1,I
        COV(J,I)=0.0
        DO IK=1,M
          IF(SIGMA(IK).GT.REPS*SIGMAX) COV(J,I)=COV(J,I)
     1            +V(J,IK)*V(I,IK)/SIGMA(IK)**2
        ENDDO
        COV(I,J)=COV(J,I)
      ENDDO
      ENDDO

      END
 
!     --------------------------------------------------------
 
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
!      IMPLICIT REAL*8(A-H,O-Z)
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
 
!     ------------------------------------------------------------
 
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
!      IMPLICIT REAL*8(A-H,O-Z)
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
 
!     -------------------------------------------------------------
 
!	To calculate the B-spline basis functions at a specified point
!
!	X : (input) Real array of length NX containing the knots.
!		The knots must be distinct and in ascending order.
!	NX : (input) Number of knots
!	K : (input) Order of B-spline, 0< K <KMAX+1
!		K=4 gives cubic B-splines
!	XB : (input) The point at which B-spline basis functions are to be evaluated
!	NDERIV : (input) Number of derivatives required
!		NDERIV.LE.0 only B-splines are calculated
!		NDERIV=1 first derivative is also calculated
!		NDERIV>1 first and second derivatives are also calculated
!	B : (output) Array of length NX+K-2 containing the value of
!		B-spline basis functions
!	DB : (output) Array of length NX+K-2 containing the value of
!		the first derivative of B-spline basis functions (if NDERIV>0)
!	DDB : (output) Array of length NX+K-2 containing the value of
!		the second derivative of B-spline basis functions (if NDERIV>1)
!	LEFT : (output) XB is located between X(LEFT) and X(LEFT+1)
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=26 implies XB > X(NX)
!		IER=27 implies XB < X(1)
!		IER=203 implies NX<2, K<1 or K>KMAX
!	WK : Real array of length NX+2K+1 used as scratch space
!
!	Required routines : None

      SUBROUTINE BSPLIN(X,NX,K,XB,NDERIV,B,DB,DDB,LEFT,IER,WK)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (KMAX=20)
      DIMENSION X(NX),B(NX+K-2),DR(KMAX),DL(KMAX),DB(NX+K-2),
     1       DDB(NX+K-2),WK(-K:NX+K)
 
      SAVE
      DATA LOW/0/
 
      IF(NX.LE.1.OR.K.LT.1.OR.K.GT.KMAX) THEN
        IER=203
        RETURN
      ENDIF
 
      IER=0
 
      IF(LOW.LT.1.OR.LOW.GE.NX) THEN
!	If the previous value of LOW is inadmissible, set the range to (1,N)
        LOW=1
        IGH=NX
      ELSE
        IGH=LOW+1
      ENDIF
 
1000  IF((XB.LT.X(LOW).AND.XB.LT.X(IGH)).OR.
     1  (XB.GT.X(LOW).AND.XB.GT.X(IGH))) THEN
!	Extend the range
        IF(XB.GT.X(LOW)) THEN
!	Extend the range on higher side
          IF(IGH.GE.NX) THEN
            IER=26
            LOW=NX-1
          ELSE
            NIGH=MIN(NX,IGH+2*(IGH-LOW))
            LOW=IGH
            IGH=NIGH
            GO TO 1000
          ENDIF
        ELSE
!	Extend the range on lower side
          IF(LOW.LE.1) THEN
            IER=27
          ELSE
            NIGH=LOW
            LOW=MAX(1,LOW-2*(IGH-LOW))
            IGH=NIGH
            GO TO 1000
          ENDIF
        ENDIF
      ELSE
 
!	Once the point is bracketed between two tabular points locate it by bisection
1500    IF(IGH-LOW.GT.1.AND.XB.NE.X(LOW)) THEN
          MID=(LOW+IGH)/2
          IF(XB.LE.X(MID).EQV.XB.LE.X(LOW)) THEN
            LOW=MID
          ELSE
            IGH=MID
          ENDIF
          GO TO 1500
        ENDIF
      ENDIF
 
!	Evaluate the B-spline basis functions

!	Define the extra knots on either side of table
!	Note that the program assumes knots from -K+2 to NX+K-1
!	and the B-splines B_{i,k}, i ranges from 1 to NX+K-2
!	The knots are stored in scratch array WK.
      DO 1700 I=1,NX
        WK(I)=X(I)
1700  CONTINUE
      DO 1800 I=1,K
        WK(1-I)=X(1)
        WK(NX+I)=X(NX)
1800  CONTINUE

      DO 1900 I=1,NX+K-2
        B(I)=0.0
        DB(I)=0.0
        DDB(I)=0.0
1900  CONTINUE
      LEFT=LOW
      LX=LOW-1
      J=1
      B(LX+1)=1.
 
!	The recurrence relation for B-splines
      DO 3000 J=1,K-1
        DR(J)=WK(LOW+J)-XB
        DL(J)=XB-WK(LOW+1-J)
        T1=0.0
        DO 2000 I=1,J
          T2=B(LX+I)/(DR(I)+DL(J+1-I))
          B(LX+I)=T1+T2*DR(I)
          T1=T2*DL(J+1-I)
2000    CONTINUE
        B(LX+J+1)=T1

!	Calculate the first derivative using recurrence relations
        IF(J.EQ.K-2.AND.NDERIV.GT.0) THEN
          T1=0.0
          DO 2200 I=1,J+1
            T2=B(LX+I)/(WK(LOW+I)-WK(LOW+I+1-K))
            DB(LX+I)=(K-1)*(T1-T2)
            T1=T2
2200      CONTINUE
          DB(LX+J+2)=(K-1)*T1
        ENDIF

!	Calculate the second derivative using recurrence relations
        IF(J.EQ.K-3.AND.NDERIV.GT.1) THEN
          T2=0.0
          P1=0.0
          DO 2400 I=1,J+1
            T3=B(LX+I)/(WK(LOW+I)-WK(LOW+I+2-K))
            P2=(T2-T3)/(WK(LOW+I)-WK(LOW+I-K+1))
            DDB(LX+I)=(K-2)*(K-1)*(P1-P2)
            T2=T3
            P1=P2
2400      CONTINUE
          P2=T2/(WK(LOW+J+2)-WK(LOW+J+3-K))
          DDB(LX+J+2)=(K-2)*(K-1)*(P1-P2)
          DDB(LX+J+3)=(K-2)*(K-1)*P2
        ENDIF
3000  CONTINUE
 
!	For K=2 the first derivative has to be calculated outside the loop
      IF(K.EQ.2.AND.NDERIV.GT.0) THEN
        T2=1/(WK(LOW+1)-WK(LOW+2-K))
        DB(LX+1)=-T2
        DB(LX+2)=T2
      ENDIF
 
!	For K=3 the second derivative has to be calculated outside the loop
      IF(K.EQ.3.AND.NDERIV.GT.1) THEN
        T3=1./(WK(LOW+1)-WK(LOW+3-K))
        P2= -T3/(WK(LOW+1)-WK(LOW-K+2))
        DDB(LX+1)=-2.*P2
        P1=P2
        P2=T3/(WK(LOW+2)-WK(LOW+3-K))
        DDB(LX+2)=2.*(P1-P2)
        DDB(LX+3)=2.*P2
      ENDIF
 
      END
 
!     ---------------------------------------------------------------
 
!	To calculate the Singular Value Decomposition of a matrix A=U D V-transpose
!
!	N : (input) Number of variables
!	M : (input) Number of equations
!	A : (input/output) Matrix of coefficients of size LA*N
!		After execution it will contain the matrix U
!	V : (output) The matrix V of size LV*N
!	SIGMA : (output) Array of length N, containing the singular values
!	LA : (input) Actual value of first dimension of A in the calling program
!	LV : (input) Actual value of first dimension of V in the calling program
!	E : Scratch array of length N
!	IER : (output) Error parameter, IER=0 if execution is successful
!		IER=12 QR iteration failed to converge to required accuracy
!		IER=105 implies N.LE.0, N.GT.LV, M.LE.0, M.GT.LA, N.GT.M
!
!	Required routines : None

      SUBROUTINE SVD(N,M,A,V,SIGMA,LA,LV,E,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
!      PARAMETER(ITMAX=30,REPS=1.D-16)
!	For REAL*4 use REPS=6.E-8
      PARAMETER(ITMAX=30,REPS=6.E-8)
      DIMENSION A(LA,N),V(LV,N),SIGMA(N),E(N)

      IF(N.GT.M.OR.N.LE.0.OR.M.LE.0.OR.M.GT.LA.OR.N.GT.LV) THEN
        IER=105
        RETURN
      ENDIF

      IER=0
!	Reduction to Bidiagonal form using Householder transformations
      G=0
      RMAX=0

      DO 3000 I=1,N
!	Off-diagonal elements of bidiagonal form
        E(I)=G
        S=0
        DO 1200 J=I,M
1200    S=S+A(J,I)**2

        IF(S.LE.0.0) THEN
!	transformation not required
          G=0
        ELSE
          F=A(I,I)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I)=F-G

          DO 1800 J=I+1,N
            S=0
            DO 1400 K=I,M
1400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 1600 K=I,M
1600        A(K,J)=A(K,J)+F*A(K,I)
1800      CONTINUE
        ENDIF

!	Diagonal elements of bidiagonal form
        SIGMA(I)=G
        S=0
        DO 2000 J=I+1,N
2000    S=S+A(I,J)**2

        IF(S.LE.0.0) THEN
!	Transformation not required
          G=0
        ELSE
          F=A(I,I+1)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I+1)=F-G
          DO 2200 J=I+1,N
!	Temporary storage of intermediate results
2200      E(J)=A(I,J)/H

          DO 2800 J=I+1,M
            S=0
            DO 2400 K=I+1,N
2400        S=S+A(J,K)*A(I,K)
            DO 2600 K=I+1,N
2600        A(J,K)=A(J,K)+S*E(K)
2800      CONTINUE
        ENDIF
        R1=ABS(SIGMA(I))+ABS(E(I))
        IF(R1.GT.RMAX) RMAX=R1
3000  CONTINUE

!	Accumulation of right hand transformation in array V
      DO 4000 I=N,1,-1
        IF(G.NE.0.0) THEN
          H=A(I,I+1)*G
          DO 3200 J=I+1,N
3200      V(J,I)=A(I,J)/H

          DO 3800 J=I+1,N
            S=0
            DO 3400 K=I+1,N
3400        S=S+A(I,K)*V(K,J)
            DO 3600 K=I+1,N
3600        V(K,J)=V(K,J)+S*V(K,I)
3800      CONTINUE
        ENDIF

        DO 3900 J=I+1,N
          V(I,J)=0.0
          V(J,I)=0.0
3900    CONTINUE
        V(I,I)=1
        G=E(I)
4000  CONTINUE

!	Accumulation of left hand transformation overwritten on matrix A
      DO 5000 I=N,1,-1
        G=SIGMA(I)
        DO 4200 J=I+1,N
4200    A(I,J)=0
        IF(G.NE.0.0) THEN
          H=A(I,I)*G

          DO 4700 J=I+1,N
            S=0
            DO 4400 K=I+1,M
4400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 4600 K=I,M
4600        A(K,J)=A(K,J)+F*A(K,I)
4700      CONTINUE

          DO 4800 J=I,M
4800      A(J,I)=A(J,I)/G
        ELSE
          DO 4900 J=I,M
4900      A(J,I)=0.0
        ENDIF
        A(I,I)=A(I,I)+1
5000  CONTINUE

!	Diagonalisation of the bidiagonal form
      AEPS=REPS*RMAX
!	Loop over the singular values
      DO 8000 K=N,1,-1
!	The QR transformation
        DO 7500 ITR=1,ITMAX

!	Test for splitting
          DO 5200 L=K,1,-1
            IF(ABS(E(L)).LT.AEPS) GO TO 6000
            IF(ABS(SIGMA(L-1)).LT.AEPS) GO TO 5400
5200      CONTINUE

!	cancellation of E(L) if L>1
5400      C=0.0
          S=1.0
          DO 5800 I=L,K
            F=S*E(I)
            E(I)=C*E(I)
            IF(ABS(F).LT.AEPS) GO TO 6000
            G=SIGMA(I)
            SIGMA(I)=SQRT(F*F+G*G)
            C=G/SIGMA(I)
            S=-F/SIGMA(I)

            DO 5600 J=1,M
              R1=A(J,L-1)
              R2=A(J,I)
              A(J,L-1)=R1*C+R2*S
              A(J,I)=C*R2-S*R1
5600        CONTINUE
5800      CONTINUE

6000      Z=SIGMA(K)
          IF(L.EQ.K) THEN
!	QR iteration has converged
            IF(Z.LT.0.0) THEN
              SIGMA(K)=-Z
              DO 6200 J=1,N
6200          V(J,K)=-V(J,K)
            ENDIF
            GO TO 8000
          ENDIF

          IF(ITR.EQ.ITMAX) THEN
            IER=12
            GO TO 7500
          ENDIF

!	calculating shift from bottom 2x2 minor
          X=SIGMA(L)
          Y=SIGMA(K-1)
          G=E(K-1)
          H=E(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.*H*Y)
          G=SQRT(1.+F*F)
          IF(F.LT.0.0) G=-G
          F=((X-Z)*(X+Z)+H*(Y/(F+G)-H))/X

!	next QR transformation
          C=1.0
          S=1.0
!	Given's rotation
          DO 7000 I=L+1,K
            G=E(I)
            Y=SIGMA(I)
            H=S*G
            G=C*G
            E(I-1)=SQRT(F*F+H*H)
            C=F/E(I-1)
            S=H/E(I-1)
            F=C*X+S*G
            G=C*G-S*X
            H=S*Y
            Y=C*Y

            DO 6400 J=1,N
              X=V(J,I-1)
              Z=V(J,I)
              V(J,I-1)=C*X+S*Z
              V(J,I)=C*Z-S*X
6400        CONTINUE

            SIGMA(I-1)=SQRT(F*F+H*H)
            IF(SIGMA(I-1).NE.0.0) THEN
              C=F/SIGMA(I-1)
              S=H/SIGMA(I-1)
            ENDIF
            F=C*G+S*Y
            X=C*Y-S*G
            DO 6600 J=1,M
              Y=A(J,I-1)
              Z=A(J,I)
              A(J,I-1)=C*Y+S*Z
              A(J,I)=C*Z-S*Y
6600        CONTINUE
7000      CONTINUE

          E(L)=0
          E(K)=F
          SIGMA(K)=X
7500    CONTINUE
8000  CONTINUE
      END
 
!     -------------------------------------------------------------
 
!	To evaluate the solution of a system of linear equations using SVD
!
!	N : (input) Number of variables
!	M : (input) Number of equations
!	U : (input) array of size LU*N containing the left-hand transformation
!	V : (input) array of size LV*N containing the right-hand transformation
!	SIGMA : (input) array of size N containing the singular values
!	LU : (input) First dimension of array U in the calling program
!	LV : (input) First dimension of array V in the calling program
!	B : (input/output) Array of length M containing the RHS
!		after execution it will contain the solution
!	WK : Scratch array of length N
!	REPS : (input) Relative accuracy.
!               All singular values < REPS*(Max of singular values)
!		will be reduced to zero
!
!	Required routines : None

      SUBROUTINE SVDEVL(N,M,U,V,SIGMA,LU,LV,B,WK,REPS)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(LU,N),V(LV,N),SIGMA(N),B(*),WK(N)

!	Finding the largest singular value
      SMAX=0.0
      DO 2000 I=1,N
        IF(SIGMA(I).GT.SMAX) SMAX=SIGMA(I)
2000  CONTINUE

      AEPS=SMAX*REPS
      DO 3000 I=1,N
        S=0.0
!	Only SIGMA(I) > AEPS contribute to the solution
        IF(SIGMA(I).GT.AEPS) THEN
          DO 2400 J=1,M
2400      S=S+U(J,I)*B(J)
          S=S/SIGMA(I)
        ENDIF
        WK(I)=S
3000  CONTINUE

      DO 4000 I=1,N
        S=0.0
        DO 3400 J=1,N
3400    S=S+V(I,J)*WK(J)
        B(I)=S
4000  CONTINUE
      END
 
!     ------------------------------------------------------
 
!	To generate random numbers with Gaussian probability distribution
!	It generates random numbers with zero mean and variance of 1.
!	
!	SEED : (input/output) real seed, it should be positive and
!		less than AM. It is updated by the routine and should
!		not be modified between two calls, unless a fresh
!		sequence is required
!
!       THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
!       VERSION AS THE SEED IS NOW REAL INSTEAD OF INTEGER.
!
!	Required routines : None

      FUNCTION RANGAU(SEED)
!      IMPLICIT REAL*8(A-H,O-Z)
!	Retain the following declaration even for REAL*4 version
!	otherwise AN and A will be rounded
      REAL*8 AM,A,AC,AN
      PARAMETER(AM=2147483648D0,A=45875D0,AC=453816693D0)
      PARAMETER(PI=3.14159265358979324D0,AN=2147483647D0)

      R1=MOD(A*SEED+AC,AM)
      IF(SEED.EQ.0.0) SEED=0.1D0
      RANGAU=SQRT(2.D0*LOG(AN/SEED))*COS(2.0*PI*R1/AN)
      SEED=R1
      END
