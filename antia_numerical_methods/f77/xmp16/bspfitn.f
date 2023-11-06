!     LEAST SQUARES FIT USING B-SPLINES IN 3 DIMENSIONS

      PROGRAM BSP3
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION T(9000),B(10,10,90),VX(90*90,10)
      DIMENSION AX(90*90,10),WK(90000),F(10,10,20),XF(90,10)
      DIMENSION EF(10,10,20),Y0(10,10,20),A2(9000,2000),V2(2000,2000)
      DIMENSION NK(10),X(90,10),MK(10),IWK(500),X0(10)
 
!	EXERCISE  10.70

      FUN(XX,Y,Z)=SIN(XX)*SIN(Y)*SIN(Z)+1.Q-3*RANGAU(SEED)
      F0(XX,Y,Z)=SIN(XX)*SIN(Y)*SIN(Z)
 
51    FORMAT('   IER =',I4,3X,'NO. OF POINTS =',I5,3X,'NO. OF KNOTS =',
     1       I3,4X,'ORDER OF B-SPLINE =',I3/5X,
     2       'REGULARISATION PARAMETER =',1PD12.4,5X,
     3       'ORDER OF DERIVATIVE =',I2/5X,'CHI SQUARE =',D12.4,5X,
     4       'IT =',I3)
52    FORMAT('   IER =',I4,3X,'X =',1P3D14.6,5X,'F(X) =',D14.6)
53    FORMAT('   EXACT VALUE :   F(X) =',1PD14.6)
 
!	NX is fixed by dimensions
      NX=10
100   PRINT *,'TYPE K=ORDER OF B-SPLINE,'
      PRINT *,'     NO=NO. OF KNOTS, RLM=REGULARISATION PARAMETER,'
      PRINT *,'     IDE=ORDER OF DERIVATIVE FOR REGULARISATION'
      PRINT *,'                 (QUITS WHEN K<1)'
      READ *,K,NO,RLM,IDE
      IF(K.LT.1) STOP

      PRINT *,'TYPE IT =1/2 FOR BSPFITN/BSPFITWN'
      READ *,IT

      ND=3
!     Set the seed here so that same set of random numbers is used every time
      SEED=2

!	Set up the input data table
!	For simplicity assume the no. of knots to be the same along each axis
!	The input values must be such that No. of knots is sufficient for BSPFITWN
      DO 900 I=1,ND
        NK(I)=NX
        MK(I)=NO
900   CONTINUE
      H=1.0Q0/(NX-1.)
      DO 1000 I=1,NX
        DO 1000 J=1,ND
          X(I,J)=(I-1)*H
1000  CONTINUE
      H=1.Q0/(NO-1.)
      DO 1200 I=1,NO
        DO 1200 J=1,ND
          XF(I,J)=(I-1)*H
1200  CONTINUE
      DO 1500 I1=1,NX
        DO 1500 I=1,NX
          DO 1500 J=1,NX
            F(J,I,I1)=FUN(X(J,1),X(I,2),X(I1,3))
            EF(J,I,I1)=1.Q-3
1500  CONTINUE
      LA=90
      L2=9000
      IV2=2000
      IV=90
      IFLG=0
      REPS=1.Q-12
      IF(IT.EQ.1) THEN
        CALL BSPFITN(ND,NK,X,F,K,AX,LA,VX,IV,T,B,XF,MK,Y0,
     1             WK,IWK,REPS,RLM,IDE,CHISQ,IER)
      ELSE
        CALL BSPFITWN(ND,NK,X,F,EF,K,A2,LA,V2,IV,T,B,XF,MK,Y0,
     1              WK,IWK,REPS,RLM,IDE,CHISQ,IER)
      ENDIF
      WRITE(6,51) IER,NX,NO,K,RLS,IDE,CHISQ,IT
      IDE=2

2000  PRINT *,'TYPE X VALUE AT WHICH FUNCTION IS REQUIRED'
      PRINT *,'                 (QUITS WHEN X0(1)<-100)'
      READ *,(X0(I),I=1,ND)
      IF(X0(1).LT.-100) GO TO 100
      F1= BSPEVN(ND,MK,XF,LA,K,B,X0,WK,IWK,IER)
      WRITE(6,52) IER,(X0(I),I=1,ND),F1
      WRITE(6,53) F0(X0(1),X0(2),X0(3))
      GO TO 2000
 
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
      IMPLICIT REAL*16(A-H,O-Z)
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

!     ------------------------------------------------------------
 
!	To calculate function value using B-spline expansion in n-dimensions
!
!	N : (input) Number of dimensions
!	NK : (input) Array of length N containing the number of knots
!		along each dimension
!	X : (input) Real array of length NXD*N containing the knots
!		X(I,J) is the Ith knot along Jth dimension
!	NXD : (input) First dimension of array X in calling program.
!		NXD must be greater than or equal to maximum of NK(I).
!	K : (input) Order of B-splines, K=4 for cubic B-splines
!	WT : (input) Coefficients of B-spline expansion, the coefficients
!		are assumed to be stored in natural FORTRAN order with
!		no gaps in data. The calling program should have
!	 	dimension  WT(NK(1)+K-2,NK(2)+K-2,...,NK(N)+K-2)
!	X0 : (input) Real array of length N, containing the coordinates
!		of the point at which function needs to be evaluated
!	WK : Real scratch array of length (NXD+K)*(3N+1)+K+2
!	IWK : Integer scratch array of length 2N
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values of IER may be set by BSPLIN which is called
!
!	BSPEVN will give the value of the B-spline expansion at X0.
!	This routine does not calculate the derivatives, for that
!	use BSPEVN1 or BSPEVN2.
!
!	Required routines : BSPLIN
 
      FUNCTION BSPEVN(N,NK,X,NXD,K,WT,X0,WK,IWK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION NK(N),X(NXD,N),WT(*),X0(N),WK(*),IWK(2*N)
 
      K1=(NXD+K)*3*N+1
      NDERIV=0
      BSPEVN=0.0
!	Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*(NXD+K)*(I-1)+1
        N2=(NXD+K)*(3*I-2)+1
        N3=(NXD+K)*(3*I-1)+1
        CALL BSPLIN(X(1,I),NK(I),K,X0(I),NDERIV,WK(N1),WK(N2),WK(N3),
     1              IWK(I),IER,WK(K1))
        IF(IER.GT.100) RETURN
        IWK(N+I)=0
1000  CONTINUE
 
!	Calculate the summation over n dimensions
      F=0.0
2000  INDEX=IWK(1)+IWK(N+1)
      TERM=WK(INDEX)
      NDP=NK(1)+K-2
      DO 2200 I=2,N
        N1=(I-1)*3*(NXD+K)
        TERM=TERM*WK(N1+IWK(I)+IWK(N+I))
        INDEX=INDEX+(IWK(I)+IWK(I+N)-1)*NDP
        NDP=NDP*(NK(I)+K-2)
2200  CONTINUE
      F=F+TERM*WT(INDEX)
 
      J=1
!	Choose the next point
2400  IF(IWK(J+N).GE.K-1) GO TO 2600
      IWK(J+N)=IWK(J+N)+1
      GO TO 2000
 
!	If Jth dimension is exhausted go to next one
2600  IWK(J+N)=0
      J=J+1
      IF(J.LE.N) GO TO 2400
 
      BSPEVN=F
 
      END
 
!     -------------------------------------------------------------
 
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
      IMPLICIT REAL*16(A-H,O-Z)
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
 
!     ---------------------------------------------------------------
 
!	To calculate linear least squares fit to B-spline basis functions in N dimension
!	Weights are assumed to be unity
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the number of
!		data points along each dimension
!	X : (input) Real array of length LA*N containing the coordinates
!		of points at which function values is available
!		X(I,J) contains the Ith point along Jth axis
!	F : (input) Real array of length NK(1)NK(2)...NK(N) containing the
!		function values. The dimension of F in the calling program
!		must match the size along each dimension, F(NK(1),...,NK(N))
!	K : (input) Order of B-splines required, K=4 gives cubic B-splines
!	A : (output) Real array of length LA*IV*N containing the matrix
!		U of SVD of the design matrix for fit along each axis
!	LA : (input) First dimension of arrays X and A as declared
!		in the calling program (LA .GE. 2*MAX(NK(I)))
!	V : (output) Real array of length IV*IV*N containing the matrix
!		V of SVD of the design matrix for fit along each axis
!	IV : (input) First dimension of V, XF, SIGMA in the calling program
!		IV .GE. MAX(MK(I))+K-2
!	SIGMA : (output) Real array of length IV*N containing the singular
!		values of the design matrix for fit along each axis
!	C : (output) Real array of length 2NK(1)NK(2)...NK(N) containing
!		the fitted coefficients. Note that although the number of
!		coefficients is (MK(1)+K-2)(MK(2)+K-2)...(MK(N)+K-2)
!		the rest of array is used as scratch space
!		Dimensions of array could be declared as
!		C(MK(1)+K-2, MK(2)+K-2, ..., MK(N-1)+K-2, NX)
!		where the last dimension is increased suitably to accommodate
!		the scratch space required. The last dimension NX should
!		be chosen such that the total size of array exceeds the
!		required value.
!	XF : (input) Real array of size IV*N, containing the knots
!		along each axis used for defining B-spline basis functions.
!		XF(I,J) should contain the Ith knot along Jth dimension
!	MK : (input) Integer array of length N containing the number of
!		knots for B-splines along each axis.
!	FY : (output) Real array of length NK(1)NK(2)...NK(N) containing
!		the values of fitted function at each of the tabular points
!		Dimensions of this array must match those of F.
!	WK : Real array of length 2NK(1)NK(2)...NK(N)+2LA used as scratch space
!	IWK : Integer array of length 3N used as scratch space
!	REPS : (input) Required accuracy for solution of equations using SVD
!		singular values less than REPS times maximum will be set to zero
!	RLM : (input) Parameter lambda for smoothing. If RLM.LE.0 no smoothing
!		is applied
!	IDE : (input) Order of derivative to be used for smoothing
!		This is used only when RLM>0. IDE=1 for first derivative
!		and IDE=2 for second derivative smoothing
!	CHISQ : (output) The value of Chi square at minimum
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=609 implies that RLM>0 and IDE is not acceptable
!		No calculations are done in this case
!		Other values of IER may be set by SVD or BSPFIT
!
!	Required routines : BSPFIT, BSPLIN, BSPEVL, BSPEVN, SVD, SVDEVL
!
 
      SUBROUTINE BSPFITN(N,NK,X,F,K,A,LA,V,IV,SIGMA,C,XF,MK,FY,WK,IWK,
     1         REPS,RLM,IDE,CHISQ,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(LA,N),A(LA*IV,*),XF(IV,N),C(*),F(*),V(IV*IV,*),
     1         SIGMA(IV,N),FY(*),WK(*),NK(N),MK(N),IWK(3*N)
 
      IF(RLM.GT.0.0.AND.(IDE.LT.1.OR.IDE.GT.2)) THEN
        IER=609
        RETURN
      ENDIF
 
      N1=4*LA+K+12
      N2=N1+LA+2
!	Set the weights to one for fits along each dimension
      DO 1000 I=0,LA
        WK(N1+I)=1.0
1000  CONTINUE
!	Calculate the SVD of matrices for fit along each dimension
      DO 2000 I=1,N
        IFLG1=1
        CALL BSPFIT(NK(I),X(1,I),F,WK(N1),K,A(1,I),LA,V(1,I),IV,
     1   SIGMA(1,I),C,XF(1,I),MK(I),FY,IFLG1,WK,REPS,RLM,IDE,CHISQ,
     2             WK(N2),IER)
        IF(IER.GT.100) RETURN
2000  CONTINUE
 
      NX=1
      DO 2200 I=1,N-1
        NX=NX*NK(I)
2200  CONTINUE
      N0=NX*NK(N)+1
      IF(RLM.GT.0.0) N0=2*NX*NK(N)+1
      NY=1
      JU=N
 
      IF(MOD(N,2).EQ.1) THEN
!	If N is odd then fit along the last dimension outside the loop
        N1=NK(N)
        IF(RLM.GT.0.0) N1=2*NK(N)
!	Set up the RHS for fit
        LJ=N1
        DO 2500 I=1,NX
          DO 2500 J=1,NK(N)
            WK(J+(I-1)*LJ)=F(I+(J-1)*NX)
            IF(RLM.GT.0.0) WK(J+NK(N)+(I-1)*LJ)=0.0
2500    CONTINUE
        M=MK(N)+K-2
        DO 2800 I=1,NX
          NI=1+(I-1)*LJ
          CALL SVDEVL(M,N1,A(1,N),V(1,N),SIGMA(1,N),LA,IV,WK(NI),WK(N0)
     1               ,REPS)
2800    CONTINUE
 
!	setup the RHS for fit along the N-1 th dimension
        NX1=NX/NK(N-1)
        NY=M
        LJ1=NK(N-1)
        IF(RLM.GT.0.0) LJ1=2*NK(N-1)
        DO 3000 I1=1,NK(N-1)
          DO 3000 I=1,NX1
            DO 3000 J=1,NY
              C(I1+(I-1)*LJ1+(J-1)*NX1*LJ1)=WK(J+(I-1)*LJ+(I1-1)*NX1*LJ)
              IF(RLM.GT.0.0) C(I1+NK(N-1)+(I-1)*LJ1+(J-1)*NX1*LJ1)=0.0
3000    CONTINUE
        NX=NX1
        JU=N-1
      ELSE
!	setup the RHS for fit along the N th dimension
        LJ=NK(N)
        IF(RLM.GT.0.0) LJ=2*NK(N)
        DO 3200 I=1,NX
          DO 3200 J=1,NK(N)
            C(J+(I-1)*LJ)=F(I+(J-1)*NX)
            IF(RLM.GT.0.0) C(J+NK(N)+(I-1)*LJ)=0.0
3200    CONTINUE
 
      ENDIF
 
!	Loop for fit along each dimension, each pass fits along 2 dimensions
      DO 5000 J1=JU,1,-2
        NUM=NX*NY
        LJ=NK(J1)
        IF(RLM.GT.0.0) LJ=2*NK(J1)
        N1=LJ
        M=MK(J1)+K-2
        DO 3400 I=1,NUM
          NI=1+(I-1)*LJ
          CALL SVDEVL(M,N1,A(1,J1),V(1,J1),SIGMA(1,J1),LA,IV,C(NI),
     1                WK(N0),REPS)
3400    CONTINUE
 
!	Set up the RHS for fit along next dimension
        NX1=NX/NK(J1-1)
        NY1=NY*M
        LJ1=NK(J1-1)
        IF(RLM.GT.0.0) LJ1=2*NK(J1-1)
        DO 3500 I1=1,NY
          DO 3500 I2=1,M
            DO 3500 I=1,NK(J1-1)
              DO 3500 J=1,NX1
                WK(I+(J-1)*LJ1+(I2-1)*LJ1*NX1+(I1-1)*NX1*LJ1*M)=
     1                      C(I2+(J-1)*LJ+(I-1)*NX1*LJ+(I1-1)*NX*LJ)
                IF(RLM.GT.0.0) WK(I+NK(J1-1)+(J-1)*LJ1+(I2-1)*LJ1*NX1+
     1                             (I1-1)*NX1*LJ1*M)=0.0
3500    CONTINUE
        NX=NX1
        NY=NY1
        NUM=NY*NX
        LJ=NK(J1-1)
        IF(RLM.GT.0.0) LJ=2*NK(J1-1)
        M=MK(J1-1)+K-2
        N1=LJ
        DO 3600 I=1,NUM
          NI=1+(I-1)*LJ
          CALL SVDEVL(M,N1,A(1,J1-1),V(1,J1-1),SIGMA(1,J1-1),LA,IV,
     1             WK(NI),WK(N0),REPS)
3600    CONTINUE
 
        IF(J1.EQ.2) THEN
!	Store the fitted coefficients in array C
          N1=NK(1)
          IF(RLM.GT.0.0) N1=2*NK(1)
          DO 3800 I=1,M
            DO 3800 J=1,NY
              C(I+(J-1)*M)=WK(I+(J-1)*LJ)
3800      CONTINUE
        ELSE
 
!	Set up the RHS for fit along next dimension
          LJ1=NK(J1-2)
          IF(RLM.GT.0.0) LJ1=2*NK(J1-2)
          M=MK(J1-1)+K-2
          NX1=NX/NK(J1-2)
          NY1=NY*M
          DO 4000 I1=1,M
            DO 4000 I2=1,NK(J1-2)
              DO 4000 I=1,NX1
                DO 4000 J=1,NY
                  C(I2+(I-1)*LJ1+(I1-1)*LJ1*NX1+(J-1)*LJ1*NX1*M)=
     1                     WK(I1+(I-1)*LJ+(I2-1)*LJ*NX1+(J-1)*LJ*NX)
                  IF(RLM.GT.0.0) C(I2+NK(J1-2)+(I-1)*LJ1+
     1                         (I1-1)*LJ1*NX1+(J-1)*LJ1*NX1*M)=0.0
4000      CONTINUE
          NX=NX1
          NY=NY1
        ENDIF
5000  CONTINUE
 
!	Calculate the Chi Square
      CHISQ=0.0
      DO 6000 I=1,N
        IWK(I)=1
6000  CONTINUE
      N0=N+1
 
!	Loop over all points
6200  IJ=IWK(1)
      WK(1)=X(IJ,1)
      ID=NK(1)
      DO 6500 I=2,N
        IJ=IJ+ID*(IWK(I)-1)
        ID=ID*NK(I)
        WK(I)=X(IWK(I),I)
6500  CONTINUE

      FY(IJ)=BSPEVN(N,MK,XF,IV,K,C,WK,WK(N0),IWK(N0),IER1)
      CHISQ=CHISQ+(F(IJ)-FY(IJ))**2
 
!	Choose the next point
      J=1
6800  IF(IWK(J).GE.NK(J)) GO TO 7000
      IWK(J)=IWK(J)+1
      GO TO 6200
 
!	If Jth dimension is exhausted go to next one
7000  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 6800
 
      END
 
!     -----------------------------------------------------------
 
!	To calculate linear least squares fit to B-spline basis functions in N dimension
!	version for general weights, would require long time to calculate
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the number of
!		data points along each dimension
!	X : (input) Real array of length LA*N containing the coordinates
!		of points at which function values is available
!		X(I,J) contains the Ith point along Jth axis
!	F : (input) Real array of length NK(1)NK(2)...NK(N) containing the
!		function values. The dimension of F in the calling program
!		must match the size along each dimension, F(NK(1),...,NK(N))
!	EF : (input) Real array of length NK(1)NK(2)...NK(N) containing the
!		estimated errors in F
!	K : (input) Order of B-splines required, K=4 gives cubic B-splines
!	A : (output) Real array of length 
!		NK(1)NK(2)...NK(N) * (MK(1)+K-2)(MK(2)+K-2)...(MK(N)+K-2)
!		containing the matrix U of SVD of the design matrix
!		If RLM>0 the size should be (N+1) times this value
!	LA : (input) First dimension of arrays X as declared
!		in the calling program (LA .GE. MAX(NK(I)))
!	V : (output) Real array of length 
!		[(MK(1)+K-2)(MK(2)+K-2)...(MK(N)+K--2)]**2
!		containing the matrix V of SVD of the design matrix
!	IV : (input) First dimension of XF in the calling program
!		IV .GE. MAX(MK(I))
!	SIGMA : (output) Real array of length 
!		(MK(1)+K-2)(MK(2)+K-2)...(MK(N)+K-2)
!		containing the singular values of the design matrix
!	C : (output) Real array of length NK(1)NK(2)...NK(N) containing
!		the fitted coefficients. If RLM>0 then the required size
!		will be N+1 times this value. Note that although the number of
!		coefficients is (MK(1)+K-2)(MK(2)+K-2)...(MK(N)+K-2)
!		the rest of array is used as scratch space
!		Dimensions of array could be declared as
!		C(MK(1)+K-2, MK(2)+K-2, ..., MK(N-1)+K-2, NX)
!		where the last dimension is increased suitably to accommodate
!		the scratch space required. The last dimension NX should
!		be chosen such that the total size of array exceeds the
!		required value.
!	XF : (input) Real array of size IV*N, containing the knots
!		along each axis used for defining B-spline basis functions.
!		XF(I,J) should contain the Ith knot along Jth dimension
!	MK : (input) Integer array of length N containing the number of
!		knots for B-splines along each axis.
!	FY : (output) Real array of length NK(1)NK(2)...NK(N) containing
!		the values of fitted function at each of the tabular points
!		Dimensions of this array must match those of F.
!	WK : Real array of length 3N*(LA+6)+LA+2K+2 used as scratch space
!	IWK : Integer array of length 3N used as scratch space
!	REPS : (input) Required accuracy for solution of equations using SVD
!		singular values less than REPS times maximum will be set to zero
!	RLM : (input) Parameter lambda for smoothing. If RLM.LE.0 no smoothing
!		is applied
!	IDE : (input) Order of derivative to be used for smoothing
!		This is used only when RLM>0. IDE=1 for first derivative
!		and IDE=2 for second derivative smoothing
!	CHISQ : (output) The value of Chi square at minimum
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=608 implies that LA or IV are not large enough
!		IER=609 implies that RLM>0 and IDE is not acceptable
!		No calculations are done in these cases
!		Other values of IER may be set by SVD or BSPLIN
!
!	Required routines : BSPLIN, BSPEVN, SVD, SVDEVL
!
 
      SUBROUTINE BSPFITWN(N,NK,X,F,EF,K,A,LA,V,IV,SIGMA,C,XF,MK,FY,
     1      WK,IWK,REPS,RLM,IDE,CHISQ,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(LA,N),A(*),XF(IV,N),C(*),F(*),V(*),SIGMA(*),FY(*),
     1      WK(*),NK(N),MK(N),IWK(3*N),EF(*)
 
 
      IF(RLM.GT.0.0.AND.(IDE.LT.1.OR.IDE.GT.2)) THEN
        IER=609
        RETURN
      ENDIF

!     N1 is the number of equations to be solved
!     M is the number of coefficients to be determined
      N1=NK(1)
      M=MK(1)+K-2
      NMAX=N1
      MMAX=MK(1)+2*K+2
      DO 2000 I=2,N
        N1=N1*NK(I)
        M=M*(MK(I)+K-2)
        NMAX=MAX(NMAX,NK(I))
        MMAX=MAX(MMAX,MK(I)+2*K+2)
2000  CONTINUE
      IF(NMAX.GT.LA.OR.MMAX.GT.IV) THEN
        IER=608
        RETURN
      ENDIF
 
      N2=N1
      IF(RLM.GT.0.0) N1=(N+1)*N1
 
      NDB=LA+5
      NDERIV=0
      IF(RLM.GT.0.0) THEN
        NDERIV=IDE
        DO 2200 I=1,N
          IWK(I)=NDB+3*(I-1)*NDB
          IF(IDE.EQ.2) IWK(I)=IWK(I)+NDB
2200    CONTINUE
      ENDIF
 
!     Set up the matrix for equations
      DO 2500 I=1,N
        XB=X(1,I)
        N3=3*(I-1)*NDB+1
        CALL BSPLIN(XF(1,I),MK(I),K,XB,NDERIV,WK(N3),WK(N3+NDB),
     1              WK(N3+NDB*2),LEFT,IER,WK(N3+3*NDB))
        IF(IER.GT.100) RETURN
2500  CONTINUE
 
      DO 3000 I=1,N
        IWK(I+N)=1
3000  CONTINUE
      N0=3*N*NDB
 
!	Loop over each tabular point
3200  IJ=IWK(N+1)
      ID=NK(1)
      DO 3300 I=2,N
        IJ=IJ+(IWK(I+N)-1)*ID
        ID=ID*NK(I)
3300  CONTINUE
      C(IJ)=F(IJ)/EF(IJ)
 
!     set each coefficient in the equation
          DO 3510 I=1,N
            IWK(2*N+I)=1
3510      CONTINUE
 
!	Loop over the coefficients
3520      JI=IWK(2*N+1)
          DO 3530 I=1,N+1
            WK(N0+I)=WK(IWK(2*N+1))
3530      CONTINUE
          IF(RLM.GT.0.0) WK(N0+2)=WK(IWK(1)+IWK(2*N+1))

          ID=MK(1)+K-2
          DO 3550 I=2,N
            JI=JI+ID*(IWK(I+2*N)-1)
            ID=ID*(MK(I)+K-2)
            WK1=WK(3*(I-1)*NDB+IWK(I+2*N))
            WK(N0+1)=WK(N0+1)*WK1
            IF(RLM.GT.0.0) THEN
              DO 3540 J=1,N
                IF(I.NE.J) THEN
                  WK(N0+J+1)=WK(N0+J+1)*WK1
                ELSE
                  WK(N0+J+1)=WK(N0+J+1)*WK(IWK(I)+IWK(I+2*N))
                ENDIF
3540          CONTINUE
            ENDIF
3550      CONTINUE

!	setup the coefficient
          A(IJ+(JI-1)*N1)=WK(N0+1)/EF(IJ)
          IF(RLM.GT.0.0) THEN
            DO 3600 I=1,N
              A(IJ+I*N2+(JI-1)*N1)=RLM*WK(N0+1+I)
3600        CONTINUE
          ENDIF
 
!	Go to the next coefficient
          J1=1
3680      IF(IWK(2*N+J1).GE.MK(J1)+K-2) GO TO 3700
          IWK(2*N+J1)=IWK(2*N+J1)+1
          GO TO 3520
 
3700      IWK(2*N+J1)=1
          J1=J1+1
          IF(J1.LE.N) GO TO 3680

!	If all coefficients are done go to the next tabular point
      J=1
3800  IF(IWK(J+N).GE.NK(J)) GO TO 4000
      IWK(J+N)=IWK(J+N)+1
      N3=(3*J-3)*NDB+1
      XB=X(IWK(J+N),J)
      CALL BSPLIN(XF(1,J),MK(J),K,XB,NDERIV,WK(N3),WK(N3+NDB),
     1            WK(N3+NDB*2),LEFT,IER,WK(N0+3*NDB))
      IF(IER.GT.100) RETURN
      GO TO 3200
 
!	If Jth dimension is exhausted go to the next
4000  IWK(J+N)=1
      N3=(3*J-3)*NDB+1
      XB=X(IWK(J+N),J)
      CALL BSPLIN(XF(1,J),MK(J),K,XB,NDERIV,WK(N3),WK(N3+NDB),
     1            WK(N3+NDB*2),LEFT,IER,WK(N0+3*NDB))
      IF(IER.GT.100) RETURN
      J=J+1
      IF(J.LE.N) GO TO 3800
 
!	If all points are done find the SVD of the matrix
      JA=N1
      JV=M
      CALL SVD(M,N1,A,V,SIGMA,JA,JV,WK,IER)
      IF(IER.GT.100) RETURN
 
 
!     Setup the remaining part of RHS and solve the equations
      DO 4100 I=N2+1,N1
        C(I)=0.0
4100  CONTINUE
 
      CALL SVDEVL(M,N1,A,V,SIGMA,JA,JV,C,WK,REPS)
 
 
!     Calculate the \chi^2
      CHISQ=0.0
      DO 4200 I=1,N
        IWK(I)=1
4200  CONTINUE
      N0=N+1
 
!	loop over all points
4300  IJ=IWK(1)
      ID=NK(1)
      WK(1)=X(IJ,1)
      DO 4500 I=2,N
        IJ=IJ+ID*(IWK(I)-1)
        ID=ID*NK(I)
        WK(I)=X(IWK(I),I)
4500  CONTINUE
      FY(IJ)=BSPEVN(N,MK,XF,IV,K,C,WK,WK(N0),IWK(N0),IER1)
      CHISQ=CHISQ+((F(IJ)-FY(IJ))/EF(IJ))**2
 
!	Go to the next point
      J=1
4800  IF(IWK(J).GE.NK(J)) GO TO 5000
      IWK(J)=IWK(J)+1
      GO TO 4300
 
5000  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 4800
 
      END
 
!     --------------------------------------------------------------
 
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
      IMPLICIT REAL*16(A-H,O-Z)
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
 
!     ------------------------------------------------------------
 
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
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(ITMAX=30,REPS=1.Q-33)
!	For REAL*4 use REPS=6.E-8
!      PARAMETER(ITMAX=30,REPS=6.E-8)
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
 
!     --------------------------------------------------------------
 
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
      IMPLICIT REAL*16(A-H,O-Z)
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
 
!     -----------------------------------------------------------
 
!	To generate random numbers with Gaussian probability distribution
!	It generates random numbers with zero mean and variance of 1.
!	
!	SEED : (input/output) real seed, it should be positive and
!		less than AM. It is updated by the routine and should
!		not be modified between two calls, unless a fresh
!		sequence is required
!
!       THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
!       VERSION AS THE SEED IS NOW REAL*16 INSTEAD OF INTEGER.
!
!	Required routines : None

      FUNCTION RANGAU(SEED)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(AM=2147483648Q0,A=45875Q0,AC=453816693Q0,
     1                  AN=2147483647Q0)
      PARAMETER(PI=3.14159265358979323846264338327950288Q0)

      R1=1+MOD(A*SEED+AC,AM)
      IF(SEED.EQ.0.0) SEED=0.1Q0
      RANGAU=SQRT(2.Q0*LOG(AN/SEED))*COS(2.0*PI*R1/AN)
      SEED=R1
      END
