!     INTEGRATION IN 3 DIMENSIONS USING B-SPLINE EXPANSION
!     IT PERFORMS B-SPLINE INTERPOLATION TO CALCULATE THE COEFFICIENTS
!     OF EXPANSION AND THEN THESE COEFFICIENTS ARE USED TO EVALUATE
!     THE INTEGRAL OVER SPECIFIED LIMITS

      PROGRAM BSPN
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION B(10,10,10),XL(90),XU(90),AX(90*90,10),INTX(90,10),
     1    WK(90000),XF(90,10),NK(20),F(10,10,10),X(90,10),MK(20)
 
51    FORMAT('  IER =',I4,5X,'K =',I4,5X,' NO. OF KNOTS =',3I4)
52    FORMAT('  IER =',I4,'  INTEGRATION OVER [',1PD14.6,',',D14.6,
     1   '] IN X,'/'  [',D14.6,',',D14.6,'] IN Y AND [',D14.6,',',D14.6,
     2   '] IN Z IS'/D14.6,'  EXACT VALUE =',D14.6)
 
      FUN(XX,Y,Z)=SIN(XX)*SIN(Y)*SIN(Z)
      FX1(XX)=-COS(XX)
 
      NX=10
      N=3
100   PRINT *,'TYPE K = ORDER OF B-SPLINE,  XL,XU, YL,YU, ZL,ZU ='
     1        ,'LIMITS OF INTEGRATION'
      PRINT *,'                        QUITS WHEN (K.LT.2)'
      READ *,K,(XL(I),XU(I),I=1,N)
      IF(K.LT.2) STOP

!	Set up the mesh points in each direction 
      H=1.Q0/(NX-1.)
      DO 1000 I=1,NX
        X(I,1)=(I-1)*H
        X(I,2)=(I-1)*H
        X(I,3)=(I-1)*H
1000  CONTINUE

!	Set up the table of values using a known function
      DO 1500 I1=1,NX
        DO 1500 I=1,NX
          DO 1500 J=1,NX
            F(J,I,I1)=FUN(X(J,1),X(I,2),X(I1,3))
1500  CONTINUE
      DO 1600 I=1,N
        NK(I)=NX
1600  CONTINUE
      LA=90
      IFLG=0

!     Calculate the coefficients of B-splines using interpolation
      CALL BSPINTN(N,NK,X,LA,F,K,AX,B,XF,MK,INTX,WK,IER)
      WRITE(6,51)IER,K,(MK(I),I=1,N)
 
      RINT=BSPQDN(N,MK,XF,LA,K,B,XL,XU,WK,IER)
      REX=(FX1(XU(1))-FX1(XL(1)))*(FX1(XU(2))-FX1(XL(2)))
     1      *(FX1(XU(3))-FX1(XL(3)))
      WRITE(6,52) IER,(XL(I),XU(I),I=1,N),RINT,REX
 
      GO TO 100
 
      END
 
!     -------------------------------------------------
 
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
      IMPLICIT REAL*16(A-H,O-Z)
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
 
!     ---------------------------------------------------------------
 
!	To calculate coefficients for B-spline interpolation in n-dimensions
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N giving the number of
!		tabular points in each dimension
!	X : (input) Array of length NXD*N containing the abscissas
!		X(I,J) is the Ith abscissa in Jth dimension
!	NXD : (input) First dimension of arrays X, XF, INTX as specified
!		in the calling program. NXD must be greater than or equal
!		to the maximum of NK(I)
!	F : (input) Array of length NK(1)*NK(2)* ... *NK(N) containing
!		the table of values. The dimensions of F in the calling
!		program must exactly match the size of table so that
!		there are no gaps in memory allocation.
!	K : (input) Order of B-spline required. K=4 gives cubic B-splines
!	AX : (output) Real array of length NXD*3*K*N containing the
!		triangular decomposition of equation matrix in band form
!		for each dimension.
!	C : (output) Real array of length NK(1)*NK(2)*... *NK(N)
!		containing the coefficients of expansion.
!	XF : (output) Real array of size NXD*N, containing
!		the knots used for B-spline calculations.
!		XF(I,J) is the Ith knot along Jth dimension.
!	MK : (output) Integer array of length N containing number of
!		knots for B-splines in each dimension.
!	INTX : (output) Integer array of length NXD*N containing information
!		about pivoting during solution of system of linear equations
!		for each dimension.
!	WK : Scratch array of length NK(1)*NK(2)*...*NK(N)+3*K
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values may be set by BSPINT, BSPLIN or GAUBND
!
!	Required routines : BSPINT, BSPLIN, GAUBND

 
      SUBROUTINE BSPINTN(N,NK,X,NXD,F,K,AX,C,XF,MK,INTX,WK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(NXD,N),NK(N),F(*),AX(NXD*3*K,N),XF(NXD,N),C(*),
     1        MK(N),WK(*),INTX(NXD,N)
 
!	Calculate the triangular decomposition of matrices for interpolation
!	along each dimension
      DO 2000 I=1,N
        IFLG1=1
        LJ=NK(I)
        CALL BSPINT(NK(I),X(1,I),F,K,AX(1,I),LJ,C,XF(1,I),MK(I),
     1         IFLG1,INTX(1,I),WK,IER)
        IF(IER.GT.100) RETURN
2000  CONTINUE
 
      KB=K-1
      NX=1
      DO 2200 I=1,N-1
        NX=NX*NK(I)
2200  CONTINUE
      N0=NX*NK(N)+1
      NY=1
      JU=N
 
!	If N is odd interpolate along last dimension outside the loop
      IF(MOD(N,2).EQ.1) THEN
        LJ=NK(N)
!	Set up the RHS
        DO 2500 I=1,NX
          DO 2500 J=1,NK(N)
            WK(J+(I-1)*LJ)=F(I+(J-1)*NX)
2500    CONTINUE
        NUM=NX
        IFLG1=2
        CALL GAUBND(NK(N),KB,NUM,AX(1,N),WK,DET,IDET,INTX(1,N),LJ,
     1         IER,IFLG1,WK(N0))
        IF(IER.GT.100) RETURN
 
!	Set up the RHS for next interpolation along N-1 th dimension
        NX1=NX/NK(N-1)
        NY=NK(N)
        LJ1=NK(N-1)
        DO 3000 I1=1,NK(N-1)
          DO 3000 I=1,NX1
            DO 3000 J=1,NY
              C(I1+(I-1)*LJ1+(J-1)*NX1*LJ1)=WK(J+(I-1)*LJ+(I1-1)*NX1*LJ)
3000    CONTINUE
        NX=NX1
        JU=N-1
      ELSE
!	Set up the RHS for interpolation along N th dimension
        LJ=NK(N)
        DO 3200 I=1,NX
          DO 3200 J=1,NK(N)
            C(J+(I-1)*LJ)=F(I+(J-1)*NX)
3200    CONTINUE
 
      ENDIF
 
!	Loop for interpolation in each dimension, each pass
!	interpolates along 2 dimensions
      DO 5000 J1=JU,1,-2
        NUM=NX*NY
        IFLG1=2
        LJ=NK(J1)
        CALL GAUBND(NK(J1),KB,NUM,AX(1,J1),C,DET,IDET,INTX(1,J1),LJ,
     1         IER,IFLG1,WK(N0))
        IF(IER.GT.100) RETURN
 
!	Set up the RHS for interpolation along the next dimension
        NX1=NX/NK(J1-1)
        NY1=NY*NK(J1)
        LJ1=NK(J1-1)
        DO 3500 I1=1,NY
          DO 3500 I2=1,NK(J1)
            DO 3500 I=1,NK(J1-1)
              DO 3500 J=1,NX1
                WK(I+(J-1)*LJ1+(I2-1)*NX+(I1-1)*NX*NK(J1))=
     1          C(I2+(J-1)*LJ+(I-1)*NX1*LJ+(I1-1)*NX*LJ)
3500    CONTINUE
        NX=NX1
        NY=NY1
        NUM=NY*NX
        IFLG1=2
        LJ=NK(J1-1)
        CALL GAUBND(NK(J1-1),KB,NUM,AX(1,J1-1),WK,DET,IDET,
     1         INTX(1,J1-1),LJ,IER,IFLG1,WK(N0))
        IF(IER.GT.100) RETURN
        IF(J1.EQ.2) THEN
!	Store the coefficients in array C
          DO 3800 I=1,NK(1)
            DO 3800 J=1,NY
              C(I+(J-1)*NK(1))=WK(I+(J-1)*LJ)
3800      CONTINUE
        ELSE
 
!	Set up the RHS for interpolation along the next dimension
          LJ1=NK(J1-2)
          NX1=NX/NK(J1-2)
          NY1=NY*NK(J1-1)
          DO 4000 I1=1,NK(J1-1)
            DO 4000 I2=1,NK(J1-2)
              DO 4000 I=1,NX1
                DO 4000 J=1,NY
                  C(I2+(I-1)*LJ1+(I1-1)*NX+(J-1)*NX*NK(J1-1))=
     1            WK(I1+(I-1)*LJ+(I2-1)*LJ*NX1+(J-1)*LJ*NX)
4000      CONTINUE
          NX=NX1
          NY=NY1
        ENDIF
5000  CONTINUE
      END
 
!     -----------------------------------------------------------
 
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
 
!     --------------------------------------------------------
 
!	To calculate the integral of a function defined by B-spline expansion
!
!	N : (input) Number of knots for B-spline expansion
!	X : (input) Real array of length N containing the knots.
!		The knots must be distinct and in ascending order.
!       K : (input) Order of B-splines, K=4 for cubic B-splines
!	WT : (input) Real array of length N+K-2 containing the
!		coefficients of B-spline expansion
!	XL : (input) Lower limit of integration
!	XU : (input) Upper limit of integration
!	WK : Scratch array of length 5N+6K+9
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=31 implies that XL is outside the range of table
!		IER=32 implies that XU is outside the range of table
!
!	BSPQD = Integral of \sum WT(I)\phi_I(x) over [XL,XU]
!
!	Required routines : BSPLIN

      FUNCTION BSPQD(N,X,K,WT,XL,XU,WK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(N),WT(N+K-2),WK(5*N+6*K+9)
 
      IER=0
      BSPQD=0.0
      K1=K+1
      IF(XL.EQ.XU) RETURN
 
!	Calculate the B-spline basis of order K+1 at XL and XU
      NDERIV=0
      NK=N+K1
      CALL BSPLIN(X,N,K1,XL,NDERIV,WK,WK(NK),WK(2*NK),LF1,IER,WK(3*NK))
      IF(IER.GT.100) RETURN
      CALL BSPLIN(X,N,K1,XU,NDERIV,WK(NK),WK(2*NK),WK(3*NK),LF2,IER,
     1            WK(4*NK))
      IF(IER.GT.100) RETURN
      IF(XL.LT.X(1).OR.XL.GT.X(N)) IER=31
      IF(XU.LT.X(1).OR.XU.GT.X(N)) IER=32

!	Define the extra knots
      K2=2*NK+K
      DO 1200 I=1,N
        WK(K2+I)=X(I)
1200  CONTINUE
      DO 1400 I=1,K
        WK(K2+1-I)=X(1)
        WK(K2+N+I)=X(N)
1400  CONTINUE
 
!	The sum for x=XL
      F1=0.0
      N1=N+K
      DO 2000 I=LF1,LF1+K
        F=0
        DO 1500 J=2,I
          F=F+WT(J-1)*(WK(J+K2)-WK(J-K+K2))
1500    CONTINUE
        F1=F1+F*WK(I)/K
2000  CONTINUE
 
!	The sum for x=XU
      F2=0.0
      DO 2500 I=LF2,LF2+K
        F=0
        DO 2200 J=2,I
          F=F+WT(J-1)*(WK(J+K2)-WK(J-K+K2))
2200    CONTINUE
        F2=F2+F*WK(N1+I)/K
2500  CONTINUE

      BSPQD=F2-F1
      END
 
!     --------------------------------------------------------
 
!	To calculate the integral of a function defined by B-spline expansion
!	over hyper-rectangular region in N dimensions
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the number of knots
!		for B-spline representation along each direction
!	X : (input) Real array of length NXD*N containing the knots along each direction
!		X(I,J) is the Ith knot along Jth dimension
!	NXD : (input) First dimension of array X as declared in the calling program
!		NXD must be greater than or equal to maximum of NK(I)
!       K : (input) Order of B-splines, K=4 for cubic B-splines
!	WT : (input) Array of length
!		(NK(1)+K-2)(NK(2)+K-2) ... (NK(N)+K-2)
!		containing coefficients of B-spline expansion. 
!		The coefficients are assumed to be stored in natural
!		Fortran order with no gaps in data. In the calling
!		program the array should have dimension
!		WT(NK(1)+K-2,NK(2)+K-2,...,NK(N)+K-2)
!	XL : (input) Real array of length N containing the lower limits of
!		integration along each direction
!	XU : (input) Real array of length N containing the upper limits of
!		integration along each direction
!	WK : Scratch array of length about twice the size of WT
!		The actual length required is complicated, but for N>2
!		the above limit should be good enough for all practical problems.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values may be set by BSPQD which is called
!
!	BSPQDN = Integral over[XL(1),XU(1)] x ... x [XL(N),XU(N)] of
!	\sum WT(I_1,...,I_N)\phi_I1(x_1) ... \phi_IN(x_N) 
!
!	Required routines : BSPLIN, BSPQD
!
 
      FUNCTION BSPQDN(N,NK,X,NXD,K,WT,XL,XU,WK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(NXD,N),NK(N),WT(*),WK(*),XL(N),XU(N)
 
      N1=1
      DO 1000 I=2,N
        N1=N1*(NK(I)+K-2)
1000  CONTINUE
      N0=N1+4
      N2=N0+N1/(NK(2)+K-2)+4
      BSPQDN=0.0
 
!	Integration along the first dimension
      NK1=NK(1)+K-2
      DO 2000 J=1,N1
        NW=(J-1)*NK1+1
        WK(J)=BSPQD(NK(1),X(1,1),K,WT(NW),XL(1),XU(1),WK(N2),IER)
        IF(IER.GT.100) RETURN
2000  CONTINUE
 
      NA=0
      NB=N0
!	Integrate along the remaining dimensions
      DO 3000 I=2,N
        NK1=NK(I)+K-2
        N1=N1/NK1
        DO 2400 J=1,N1
        NW=(J-1)*NK1+1
        WK(NB+J)=BSPQD(NK(I),X(1,I),K,WK(NW+NA),XL(I),XU(I),WK(N2),IER)
        IF(IER.GT.100) RETURN
2400    CONTINUE
        NT=NA
        NA=NB
        NB=NT
3000  CONTINUE
 
      BSPQDN=WK(NA+1)
      END
 
!     ----------------------------------------------------------
 
!     Solution of a system of linear equations using Gaussian elimination
!     	for a band matrix
!
!     N : (input) Number of equations to be solved
!     KB : (input) Bandwidth of matrix A(I,J)=0 if ABS(I-J)>KB
!     NUM : (input) Number of different sets (each with N equations) of
!		equations to be solved
!     A : (input/output) The matrix of coefficient of size LJ*(3*KB+1)
!		A(I,J-I+KB+1) is the coefficient of x_j in Ith equation
!		at output it will contain the triangular decomposition
!     X : (input/output) The matrix containing right hand sides (size LJ*NUM)
!		X(I,J) is the Ith element of Jth right hand side
!		at output it will contain the solutions
!     DET, IDET : (output) The determinant of the matrix = DET*2**IDET
!     INC : (output) Integer array of length N containing information about
!		interchanges performed during elimination
!     LJ : (input) First dimension of arrays A and X in calling program
!     IER : (output) Error flag, IER=0 signifies successful execution
!		IER=104 implies (N.LE.0 or N.GT.LJ or KB.GT.N)
!		IER=124 implies some pivot turned out to be zero and hence
!	     		matrix must be nearly singular
!     IFLG : (input) Integer variable used as a flag to specify the type
!		of computation required
!     		If IFLG=-1, both elimination and solution are calculated
!     			without pivoting and IFLG is set to 2
!		If IFLG=0, both elimination and solution are computed
!     			with partial pivoting and IFLG is set to 2
!     		If IFLG=1, only elimination is done with pivoting and IFLG is set to 2
!     		If IFLG.GE.2 only solution is calculated, the triangular
!     			decomposition should have been calculated earlier
!     WK : Real array of length 3*KB+1 used as scratch space
!
!	Required routines : None
 
      SUBROUTINE GAUBND(N,KB,NUM,A,X,DET,IDET,INC,LJ,IER,IFLG,WK)
      IMPLICIT REAL*16(A-H,O-Z)
!     For complex matrices use the following statements instead
!      IMPLICIT REAL*8(R)
!      IMPLICIT COMPLEX*16(A-H,S-Z)
 
      DIMENSION A(LJ,3*KB+1),INC(N),X(LJ,NUM),WK(3*KB+1)
 
      IF(N.LE.0.OR.N.GT.LJ.OR.KB.GT.N) THEN
        IER=104
        RETURN
      ENDIF
 
      KB1=KB+1
      IER=124
      IF(IFLG.LE.1) THEN
!     Perform elimination
        DO 2000 I=1,N
          DO 2000 J=2*KB+2,3*KB+1
            A(I,J)=0.0
2000    CONTINUE
 
        DET=1.0
        IDET=0
        DO 2600 K=1,N-1
!     Find the maximum element in the Kth column
          R1=0.0
          KM=K
          IF(IFLG.GE.0) THEN
            DO 2200 L=K,MIN(N,K+KB)
              IF(ABS(A(L,K-L+KB1)).GT.R1) THEN
                R1=ABS(A(L,K-L+KB1))
                KM=L
              ENDIF
2200        CONTINUE
          ENDIF
 
          INC(K)=KM
          IF(KM.NE.K) THEN
!     Interchange the rows if needed
            DO 2300 L=K,MIN(N,2*KB+K)
              WK(L-K+1)=A(K,L-K+KB1)
2300        CONTINUE
            DO 2400 L=K,MIN(N,2*KB+K)
              A(K,L-K+KB1)=A(KM,L-KM+KB1)
2400        A(KM,L-KM+KB1)=WK(L-K+1)
            DET=-DET
          ENDIF
 
          DET=DET*A(K,KB1)
          IF(A(K,KB1).EQ.0.0) RETURN
!     To check for singular or nearly singular matrices replace this
!     statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(K,KB1)).LT.REPS) RETURN
          IF(DET.NE.0.0) THEN
 
!     Scale the value of the determinant DET
2350        IF(ABS(DET).GT.32.) THEN
              DET=DET*0.03125Q0
              IDET=IDET+5
              GO TO 2350
            ENDIF
 
2370        IF(ABS(DET).LT.0.03125Q0) THEN
              DET=DET*32.
              IDET=IDET-5
              GO TO 2370
            ENDIF
          ENDIF
 
          DO 2500 L=K+1,MIN(N,K+KB)
            A(L,K-L+KB1)=A(L,K-L+KB1)/A(K,KB1)
            DO 2500 L1=K+1,MIN(N,2*KB+K)
2500      A(L,L1-L+KB1)=A(L,L1-L+KB1)-A(L,K-L+KB1)*A(K,L1-K+KB1)
2600    CONTINUE
        DET=DET*A(N,KB1)
        INC(N)=N
!     If pivot is zero then return, IER has been set to 124
        IF(A(N,KB1).EQ.0.0) RETURN
!     To check for singular or nearly singular matrices replace this
!     statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(N,KB1)).LT.REPS) RETURN
 
        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF
 
      IER=0
!     Solution for the NUM different right-hand sides
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
!     Forward substitution
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF
          DO 3000 L=K+1,MIN(N,K+KB)
3000    X(L,J)=X(L,J)-A(L,K-L+KB1)*X(K,J)
 
!     back-substitution
        X(N,J)=X(N,J)/A(N,KB1)
        DO 3300 K=N-1,1,-1
          DO 3200 L=MIN(N,K+2*KB),K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L-K+KB1)
3300    X(K,J)=X(K,J)/A(K,KB1)
5000  CONTINUE
      END
