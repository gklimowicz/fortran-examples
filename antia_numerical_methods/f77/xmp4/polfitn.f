!     PROGRAM TO OBTAIN LEAST SQUARES POLYNOMIAL FIT TO A GIVEN DATA
!     USING ORTHOGONAL POLYNOMIALS IN 3 DIMENSIONS
 
      PROGRAM POLYLS
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WK(99000),X(100,9),FB(20,20,20),NK(9),MK(9),XB(9),
     1          IWK(90),DFF(9),D2FF(3,3)
      DIMENSION Y(100),AX(100,16),FXY(20,20,20),C(20020)
 
!     EXAMPLE 10.3
 
!     FUNCTION WITH RANDOM ERROR TO GENERATE DATA
 
      F(X1,Y1,Y2)=((231*X1*Y1-315)*Y2*X1+105)*Y2*Y1-5+RANGAU(SEED)*1.D-3
 
!     THE TRUE FUNCTION AND ITS DERIVATIVES
 
      F1(X1,Y1,Y2)=((231*X1*Y1-315)*Y2*X1+105)*Y2*Y1-5
      DF1(X1,Y1,Y2)=(462*X1*Y1-315)*Y2**2*Y1
      DF2(X1,Y1,Y2)=((462*X1*Y1-315)*X1*Y2+105)*Y2
      DF3(X1,Y1,Y2)=((462*X1*Y1-630)*X1*Y2+105)*Y1
      DF11(X1,Y1,Y2)=462*Y2**2*Y1**2
      DF12(X1,Y1,Y2)=(924*X1*Y1-315)*Y2*Y2
      DF22(X1,Y1,Y2)=462*X1**2*Y2**2
      DF13(X1,Y1,Y2)=(924*X1*Y1-630)*Y1*Y2
      DF23(X1,Y1,Y2)=(924*X1*Y1-630)*X1*Y2+105
      DF33(X1,Y1,Y2)=(462*X1*Y1-630)*Y1*X1
 
51    FORMAT('   NO. OF PTS =',3I4,5X,'DEGREE OF POLYNOMIAL =',3I3/
     1       5X,'IER =',I4,5X,'CHI SQUARE =',1PE14.7)
52    FORMAT('  X =',1P3E14.6,3X,'F =',E14.6)
53    FORMAT(3X,'F =',1PE14.6,5X,4HF' =,3E14.6)
54    FORMAT('   COEFFICIENTS :'/(1P5E14.6))
56    FORMAT(3X,'F =',1PE14.6,5X,4HF' =,3E14.6/3X,
     1       4HF''=,6E14.6)
57    FORMAT('   TRUE VALUES :',6X,'F =',1PE14.6/3X,4HF' =,3E14.6/3X,
     1       4HF''=,6E14.6)
 
      SEED=2
      ND=3
!	No. of data points is fixed by dimension
      N=20
 
!     GENERATING THE DATA WITH RANDOM ERROR
 
      HH=1.D0/(N-1.)
      DO 1000 I=1,N
        X(I,1)=(I-1)*HH
        X(I,2)=(I-1)*HH
        X(I,3)=(I-1)*HH
1000  CONTINUE
      DO 1200 I1=1,N
        DO 1200 I=1,N
          DO 1200 J=1,N
            FXY(I,J,I1)=F(X(I,1),X(J,2),X(I1,3))
1200  CONTINUE
      NK(1)=N
      NK(2)=N
      NK(3)=N
      LA=100
      IC=100
 
      PRINT *,'TYPE M=DEGREE OF POLYNOMIAL TO BE FITTED'
      READ *,(MK(I),I=1,ND)
 
      CALL POLFITN(ND,NK,X,FXY,AX,LA,C,MK,FB,WK,IWK,CHISQ,IER)
 
      WRITE(6,51) (NK(I),I=1,ND),(MK(I),I=1,ND),IER,CHISQ
      WRITE(6,54) (C(I),I=1,(MK(1)+1)*(MK(2)+1)*(MK(3)+1))
 
2000  PRINT *,'TYPE THE X VALUE WHERE FITTED VALUE IS REQUIRED'
      PRINT *,'                 (QUITS WHEN XB(1)<-100)'
      READ *,(XB(I),I=1,ND)
      IF(XB(1).LT.-100) STOP

!	Evaluate the fitted polynomial and derivatives using all three versions of POLEVN
      CALL POLEVN(ND,MK,AX,LA,C,XB,FF,WK,IWK)
      WRITE(6,52)(XB(I),I=1,ND),FF

      CALL POLEVN1(ND,MK,AX,LA,C,XB,FF,DFF,WK,IWK)
      WRITE(6,53) FF,(DFF(I),I=1,ND)

      CALL POLEVN2(ND,MK,AX,LA,C,XB,FF,DFF,D2FF,WK,IWK)
      WRITE(6,56) FF,(DFF(I),I=1,ND),((D2FF(I,J),I=1,J),J=1,ND)

!	The exact value of function and its derivatives
      WRITE(6,57) F1(XB(1),XB(2),XB(3)),DF1(XB(1),XB(2),XB(3)),
     1           DF2(XB(1),XB(2),XB(3)),DF3(XB(1),XB(2),XB(3)),
     1           DF11(XB(1),XB(2),XB(3)),DF12(XB(1),XB(2),XB(3)),
     1           DF22(XB(1),XB(2),XB(3)),DF13(XB(1),XB(2),XB(3)),
     1           DF23(XB(1),XB(2),XB(3)),DF33(XB(1),XB(2),XB(3))
      GO TO 2000
      END
 
!     --------------------------------------------------
 
!	Evaluating the fitted polynomial at any point using known
!	coefficients of orthogonal polynomials in N dimensions.
!	Should be used to evaluate the polynomial using coefficients calculated
!	by POLFITN. It does not calculate the derivatives, for derivatives
!	use POLEVN1 or POLEVN2
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the degree of
!		polynomial in each direction
!	AX : (input) Real array of length LA*(3*N+3) containing the coefficients
!		alpha and beta for orthogonal polynomials in each direction
!		AX(I,3*J-2) contains alpha and AX(I,3*J-1) contains beta
!		for polynomials in Jth dimension
!	LA : (input) First dimension of array AX in the calling program.
!		It must be same as what was used in call to POLFITN while
!		calculating the coefficients.
!	WT : (input) Real array of length (MK(1)+1)(MK(2)+1)...(MK(N)+1)
!		containing the coefficients of the fit. The dimensions of WT in
!		the calling program must match the size along each dimension,
!		WT(MK(1)+1,MK(2)+1,...,MK(N)+1)
!	X0 : (input) Real array of length N containing the coordinates of
!		the point at which polynomial needs to be evaluated
!	F : (output) Calculated value of the fitted polynomial at X0
!	WK : Real array of length 3*N*LA used as scratch space
!	IWK : Integer array of length N used as scratch space
!	
!	Required routines : POLORT
!
      SUBROUTINE POLEVN(N,NK,AX,LA,WT,X0,F,WK,IWK)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),AX(LA,3*N),WT(*),X0(N),WK(LA*3*N),IWK(N)
 
!     Calculate the orthogonal polynomials along each dimension
      DO 1000 I=1,N
        N1=3*LA*(I-1)+1
        N2=3*LA*(I-1)+LA+1
        N3=3*LA*(I-1)+2*LA+1
        NJ=3*(I-1)+1
        XB=X0(I)
        CALL POLORT(NK(I),AX(1,NJ),AX(1,NJ+1),XB,WK(N1),WK(N2),WK(N3))
        IWK(I)=1
1000  CONTINUE
 
!     Calculate the summation over n dimensions
      F=0.0
2000  INDEX=IWK(1)
      TERM=WK(INDEX)
      NDP=NK(1)+1
      DO 2200 I=2,N
        N1=(I-1)*3*LA
        TERM=TERM*WK(N1+IWK(I))
        INDEX=INDEX+(IWK(I)-1)*NDP
        NDP=NDP*(NK(I)+1)
2200  CONTINUE
      F=F+TERM*WT(INDEX)
 
      J=1
!     Choose the next point
2400  IF(IWK(J).GE.NK(J)+1) GO TO 2600
      IWK(J)=IWK(J)+1
      GO TO 2000
 
!     If Jth dimension is exhausted go to next one
2600  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 2400
 
      END
 
!     -------------------------------------------------------
 
!	Evaluating the fitted polynomial and its first derivative at any point
!	using known coefficients of orthogonal polynomials in N dimensions.
!	Should be used to evaluate the polynomial using coefficients calculated
!	by POLFITN. It does not calculate the second derivatives, for that
!	use POLEVN2, for no derivatives use POLEVN
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the degree of
!		polynomial in each direction
!	AX : (input) Real array of length LA*(3*N+3) containing the coefficients
!		alpha and beta for orthogonal polynomials in each direction
!		AX(I,3*J-2) contains alpha and AX(I,3*J-1) contains beta
!		for polynomials in Jth dimension
!	LA : (input) First dimension of array AX in the calling program.
!		It must be same as what was used in call to POLFITN while
!		calculating the coefficients.
!	WT : (input) Real array of length (MK(1)+1)(MK(2)+1)...(MK(N)+1)
!		containing the coefficients of the fit. The dimensions of WT in
!		the calling program must match the size along each dimension,
!		WT(MK(1)+1,MK(2)+1,...,MK(N)+1)
!	X0 : (input) Real array of length N containing the coordinates of
!		the point at which polynomial needs to be evaluated
!	F : (output) Calculated value of the fitted polynomial at X0
!	DF : (output) Real array of length N containing the first derivatives
!		of F at X0
!	WK : Real array of length 3*N*LA+N used as scratch space
!	IWK : Integer array of length N used as scratch space
!	
!	Required routines : POLORT
!
      SUBROUTINE POLEVN1(N,NK,AX,LA,WT,X0,F,DF,WK,IWK)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),AX(LA,3*N),WT(*),X0(N),WK(LA*3*N+N),IWK(N),DF(N)
 
!     Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*LA*(I-1)+1
        N2=3*LA*(I-1)+LA+1
        N3=3*LA*(I-1)+2*LA+1
        NJ=3*(I-1)+1
        XB=X0(I)
        CALL POLORT(NK(I),AX(1,NJ),AX(1,NJ+1),XB,WK(N1),WK(N2),WK(N3))
        IWK(I)=1
1000  CONTINUE
 
!     calculate the sum over all points
      F=0.0
      DO 1800 I=1,N
        DF(I)=0.0
1800  CONTINUE
      N0=LA*3*N
 
2000  INDEX=IWK(1)
      TERM=WK(INDEX)
      WK(N0+1)=WK(LA+INDEX)
      DO 2100 I=2,N
        WK(N0+I)=WK(INDEX)
2100  CONTINUE
      NDP=NK(1)+1
      DO 2300 I=2,N
        N1=(I-1)*3*LA
        N2=IWK(I)
        TERM=TERM*WK(N1+N2)
 
!     terms for first derivatives
        DO 2200 J=1,N
          IF(J.EQ.I) THEN
            WK(N0+J)=WK(N0+J)*WK(N1+N2+LA)
          ELSE
            WK(N0+J)=WK(N0+J)*WK(N1+N2)
          ENDIF
2200    CONTINUE
        INDEX=INDEX+(N2-1)*NDP
        NDP=NDP*(NK(I)+1)
2300  CONTINUE
      F=F+TERM*WT(INDEX)
      DO 2400 I=1,N
        DF(I)=DF(I)+WK(N0+I)*WT(INDEX)
2400  CONTINUE
 
      J=1
!     Go to the next point
2500  IF(IWK(J).GE.NK(J)+1) GO TO 2600
      IWK(J)=IWK(J)+1
      GO TO 2000
 
!     If Jth dimension is exhausted go to the next one
2600  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 2500
 
 
      END
 
!     ---------------------------------------------------
 
!	Evaluating the fitted polynomial and its derivatives at any point
!	using known coefficients of orthogonal polynomials in N dimensions.
!	Should be used to evaluate the polynomial using coefficients calculated
!	by POLFITN. If second derivative is not required use POLEVN1, if
!	no derivatives are required then use POLEVN
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the degree of
!		polynomial in each direction
!	AX : (input) Real array of length LA*(3*N+3) containing the coefficients
!		alpha and beta for orthogonal polynomials in each direction
!		AX(I,3*J-2) contains alpha and AX(I,3*J-1) contains beta
!		for polynomials in Jth dimension
!	LA : (input) First dimension of array AX in the calling program.
!		It must be same as what was used in call to POLFITN while
!		calculating the coefficients.
!	WT : (input) Real array of length (MK(1)+1)(MK(2)+1)...(MK(N)+1)
!		containing the coefficients of the fit. The dimensions of WT in
!		the calling program must match the size along each dimension,
!		WT(MK(1)+1,MK(2)+1,...,MK(N)+1)
!	X0 : (input) Real array of length N containing the coordinates of
!		the point at which polynomial needs to be evaluated
!	F : (output) Calculated value of the fitted polynomial at X0
!	DF : (output) Real array of length N containing the first derivatives
!		of F at X0
!	DDF : (output) Real array of length N*N containing the second derivatives
!		of F at X0, DDF(I,J)=d^2F/(dX(I) dX(J))
!	WK : Real array of length 3*N*LA+N+N*N used as scratch space
!	IWK : Integer array of length N used as scratch space
!	
!	Required routines : POLORT
!
      SUBROUTINE POLEVN2(N,NK,AX,LA,WT,X0,F,DF,DDF,WK,IWK)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),AX(LA,3*N),WT(*),X0(N),WK(*),IWK(N),DF(N),DDF(N,N)
 
!     Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*LA*(I-1)+1
        N2=3*LA*(I-1)+LA+1
        N3=3*LA*(I-1)+2*LA+1
        NJ=3*(I-1)+1
        XB=X0(I)
        CALL POLORT(NK(I),AX(1,NJ),AX(1,NJ+1),XB,WK(N1),WK(N2),WK(N3))
        IWK(I)=1
1000  CONTINUE
 
!     calculate the sum over all points
 
      F=0.0
      DO 1800 I=1,N
        DF(I)=0.0
        DO 1800 J=1,N
          DDF(J,I)=0.0
1800  CONTINUE
      N0=LA*3*N
2000  INDEX=IWK(1)
      TERM=WK(INDEX)

!	Terms for the derivatives
      WK(N0+1)=WK(LA+INDEX)
      DO 2100 I=2,N+N*N
        WK(N0+I)=WK(INDEX)
2100  CONTINUE
      WK(N0+N+1)=WK(2*LA+INDEX)
      DO 2200 I=2,N
        WK(N0+N+I)=WK(INDEX+LA)
2200  CONTINUE
      NDP=NK(1)+1
      DO 2500 I=2,N
        N1=(I-1)*3*LA
        N2=IWK(I)
        TERM=TERM*WK(N1+N2)
        DO 2300 J=1,N
          IF(J.EQ.I) THEN
            WK(N0+J)=WK(N0+J)*WK(N1+N2+LA)
          ELSE
            WK(N0+J)=WK(N0+J)*WK(N1+N2)
          ENDIF
          DO 2300 J1=1,J
            IF(J.EQ.I.AND.J1.EQ.I) THEN
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2+2*LA)
            ELSE IF(J.EQ.I.OR.J1.EQ.I) THEN
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2+LA)
            ELSE
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2)
            ENDIF
2300    CONTINUE
        INDEX=INDEX+(N2-1)*NDP
        NDP=NDP*(NK(I)+1)
2500  CONTINUE
      F=F+TERM*WT(INDEX)
      DO 2600 I=1,N
        DF(I)=DF(I)+WK(N0+I)*WT(INDEX)
        DO 2600 J=1,I
          DDF(I,J)=DDF(I,J)+WK(N0+N+(J-1)*N+I)*WT(INDEX)
2600  CONTINUE
 
      J=1
!	Go to the next point
3000  IF(IWK(J).GE.NK(J)+1) GO TO 3200
      IWK(J)=IWK(J)+1
      GO TO 2000
 
!	If Jth dimension is exhausted, go to the next one
3200  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 3000
 
      DO 3500 I=1,N
        DO 3500 J=I+1,N
          DDF(I,J)=DDF(J,I)
3500  CONTINUE
 
      END
 
!     ---------------------------------------------------
 
!	Least squares polynomial fit using orthogonal polynomials in 1 dimension
!	Modified version of POLFIT to fit multiple sets of function values
!	This routine is meant to be used for fit in multiple dimensions
!
!	N : (input) Number of data points
!	M : (input) Required degree of polynomial
!	NUM : (input) Number of different RHS (function values) to be fitted
!		Each set must be defined over the same abscissas and
!		with same weights.
!	X : (input) Real array of length N containing the abscissas
!	F : (input) Real array of length N*NUM containing the function values
!		for each RHS, F(I,J) is the value at X(I) for Jth set
!		The first dimension of F is assumed to be exactly equal
!		to N to minimise storage requirement.
!	SIG : (input) Real array of length N containing the estimated
!	        error in the function values, SIG(I) is the error in F(I,J).
!		The errors are assumed to be the same for all NUM data sets.
!	A : (output) Real array of length (M+1)*NUM containing the coefficients
!		for the fit for each RHS. The first dimension of A
!		is assumed to be M+1 to minimise storage requirements
!		A(I,J) is the Ith coefficient for Jth RHS.
!	ALP, BETA : (output) Real arrays of length M+1, containing the coefficients
!		required for defining the orthogonal polynomials
!	GAM : (output) Real array of length M+1, containing the quantities
!		\gamma_i for the orthogonal polynomials
!	WK : Real array of length 2N used as scratch space
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=601 implies that N<M+1 or M<0 or N < 1
!		IER=621 implies that GAM(I) vanishes at some I and calculations
!			are abandoned
!	
!	The fitted polynomial can be calculated at any value of x using POLEVL
!
!	Required routines : None
 
      SUBROUTINE POLFIT1(N,M,NUM,X,F,SIG,A,ALP,BETA,GAM,WK,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N,NUM),SIG(N),A(M+1,NUM),ALP(M+1),BETA(M+1),
     1      WK(N,2),GAM(M+1)
 
      IF(M.GE.N.OR.M.LT.0.OR.N.LT.1) THEN
        IER=601
        RETURN
      ENDIF
 
!	initialisation
      I1=1
      I2=2
      GAM0=0.0
      DO 2000 I=1,N
        WK(I,I1)=0.0
        WK(I,I2)=1.0
        GAM0=GAM0+1/SIG(i)**2
2000  CONTINUE
      GAM(1)=GAM0
      BETA(1)=0
      IER=0
 
!	Loop over the degree
      DO 5000 J=0,M
        IF(GAM(J+1).LE.0) THEN
          IER=621
          RETURN
        ENDIF
 
        DO 3200 J1=1,NUM
          S=0.0
          DO 3000 I=1,N
3000      S=S+F(I,J1)*WK(I,I2)/SIG(I)**2
          A(J+1,J1)=S/GAM(J+1)
3200    CONTINUE
        IF(J.EQ.M) RETURN

        S=0
        DO 3400 I=1,N
3400    S=S+X(I)*(WK(I,I2)/SIG(I))**2
        ALP(J+1)=S/GAM(J+1)
        GAM0=0.0
        DO 3600 I=1,N
          WK(I,I1)=(X(I)-ALP(J+1))*WK(I,I2)-BETA(J+1)*WK(I,I1)
          GAM0=GAM0+(WK(I,I1)/SIG(I))**2
3600    CONTINUE
        BETA(J+2)=GAM0/GAM(J+1)
        GAM(J+2)=GAM0
        IT=I1
        I1=I2
        I2=IT
5000  CONTINUE
 
      END
 
!     ---------------------------------------------------
 
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
!      IMPLICIT REAL*8(A-H,O-Z)
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
 
!     ---------------------------------------------------
 
!	Evaluating the orthogonal polynomial basis functions at any value
!	of x using known coefficients
!	Should be used to evaluate the basis using coefficients calculated
!	by POLFIT or POLFIT1.
!
!	M : (input) Degree of polynomial
!	ALP, BETA : (input) Real arrays of length M+1, containing the coefficients
!		required for defining the orthogonal polynomials
!		ALP, BETA could be calculated using POLFIT
!	X : (input) Value of x at which polynomials needs to be evaluated
!	F : (output) Real array of length M+1 containing the value of
!		orthogonal polynomials at X
!	DF : (output) Real array of length M+1 containing first
!		derivative of F at X
!	DDF : (output) Real array of length M+1 containing second
!		derivative of F at X
!	
!	Required routines : None

 
      SUBROUTINE POLORT(M,ALP,BETA,X,F,DF,DDF)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ALP(M+1),BETA(M+1),F(M+1),DF(M+1),DDF(M+1)

      F(1)=1.0
      F(2)=(X-ALP(1))
      DF(1)=0.0
      DF(2)=1.0
      DDF(1)=0.0
      DDF(2)=0.0
 
!	The recurrence relations
      DO 1000 J=3,M+1
        DDF(J)=2.*DF(J-1)+(X-ALP(J-1))*DDF(J-1)-BETA(J-1)*DDF(J-2)
        DF(J)=F(J-1)+(X-ALP(J-1))*DF(J-1)-BETA(J-1)*DF(J-2)
        F(J)=(X-ALP(J-1))*F(J-1)-BETA(J-1)*F(J-2)
1000  CONTINUE
      END
 
!     ---------------------------------------------------
 
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
