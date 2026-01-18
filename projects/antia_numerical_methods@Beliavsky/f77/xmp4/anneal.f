!	Minimisation using simulated annealing

      PROGRAM EX7
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IP(20),IP1(20),IP0(20)
      COMMON ERR(600),X0(600),F0(600),X(49)
      EXTERNAL RANF

!	EXAMPLE 8.7
 
51    FORMAT(F11.2,F14.4,3F11.6)
52    FORMAT(F11.2)
53    FORMAT('  RANDOM NO. SEED =',I5,5X,'FUNCTION VALUE AT MINIMUM =',
     1         1PE12.4)
54    FORMAT('  SERIAL NO. OF KNOTS FOR MINIMUM ERROR :',15I3)
55    FORMAT('#   THE SELECTED KNOTS :')
56    FORMAT('#    X     INTERPOLATED VALUE  ERROR')
 
!     The number of knots to be fixed
      N=12
      RESMIN=1
!	The output file:
!	The first 50 lines are the input data
!	The next 13 lines contain the selected knots
!	The last 501 lines contain the interpolated values using selected knots
      OPEN(UNIT=21,FILE='titan.out',STATUS='UNKNOWN')

!     Try 2 different sequence of random numbers, -ID is the seed
!     Only the last solution is written out in the file
      DO 6000 ID=1,2
!     Initial distribution of knots
        IP(1)=1
        IP(N)=49
        DO 1000 I=2,N-1
          IP(I)=4*I
          IP1(I)=IP(I)
1000    CONTINUE
        RES0=FCN(IP)
        RESMIN=RES0
        TEMP=MAX(RES0,1.0)
        TEMP0=TEMP
        NCH=500
        IDUM=-ID
 
!       The loop for annealing 
        DO 5000 IT=1,400
          ICH=0
 
!	      Perform 10000 trials at each temperature
          DO 4000 IC=1,10000
            DO 2000 I=1,N
              IP0(I)=IP(I)
2000        CONTINUE
            CALL CHANGE(TEMP,IDUM,IP)
            RES=FCN(IP)
            IF(RES.LT.RESMIN) THEN
!		Store the minimum value
              DO 2100 I=1,N
                IP1(I)=IP(I)
2100          CONTINUE
              RESMIN=RES
            ENDIF
            EXPI=2.0
            IF(ABS((RES0-RES)/TEMP).LT.500.0) THEN
              EXPI=EXP((RES0-RES)/TEMP)
            ELSE IF(RES.GT.RES0) THEN
              EXPI=0.0
            ENDIF
 
            IF(RES.LT.RES0.OR.RANF(IDUM).LT.EXPI) THEN
!			Accept the new point
              RES0=RES
              ICH=ICH+1
            ELSE
              DO 2300 I=1,N
                IP(I)=IP0(I)
2300          CONTINUE
            ENDIF
            IF(ICH.GE.NCH) GO TO 4200
4000      CONTINUE
!         Reduce the temperature slowly
4200      TEMP=TEMP*0.98
 
          IF(RES0-RESMIN.GT.5*TEMP.OR.ICH.EQ.0) THEN
            DO 4400 I=1,N
              IP(I)=IP1(I)
4400        CONTINUE
            RES0=RESMIN
          ENDIF
 
!         The convergence criterion
          IF(ICH.EQ.0.AND.TEMP.LT.4.D-2*RES0) GO TO 5200
          IF(TEMP.LT.1.D-5*RES0) GO TO 5200
5000    CONTINUE

!       Write out the solution
5200    RES=FCN(IP1)
        WRITE(6,53) -ID,RES
        WRITE(6,54) (IP1(I),I=2,N-1)
6000  CONTINUE

!     write the solution in unit 21
      WRITE(21,55)
      WRITE(21,52) (X(IP1(I)),I=1,N)
      WRITE(21,56)
      DO 6400 I=1,501
        WRITE(21,51)X0(I),F0(I),ERR(I)
6400  CONTINUE
      END
 
!     -----------------------------------------------------
 
!	To evaluate the cubic spline interpolant at a specified point
!
!	XB : (input) point at which interpolation is required
!	N : (input) Number of points in the table
!	X : (input) Array of length N, containing the abscissas
!	F : (input) Array of length N, containing the function values at X(I)
!	C : (input) Array of length 3*N containing the spline coefficients
!		which should have been calculated using SPLINE
!	DFB : (output) First derivative of spline at x=XB
!	DDFB : (output) Second derivative of spline at x=XB
!	IER : (output) error parameter, IER=0 if execution is successful
!		IER=24 implies XB is outside the range of table on higher side
!		IER=25 implies XB is outside the range of table on lower side
!		IER=201 implies N<2
!	SPLEVL will be the interpolated value at x=XB
!
!	Required routines : None

      FUNCTION SPLEVL(XB,N,X,F,C,DFB,DDFB,IER)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION X(N),F(N),C(3,N)
      SAVE
      DATA LOW/0/

      SPLEVL=0.0
      IF(N.LE.1) THEN
        IER=201
        RETURN
      ENDIF

!	QASCND is true if table is in ascending order
      QASCND=X(N).GT.X(1)
      IER=0

      IF(LOW.LT.1.OR.LOW.GE.N) THEN
!	If the previous value of LOW is inadmissible, set the range to (1,N)
        LOW=1
        IGH=N
      ELSE
        IGH=LOW+1
      ENDIF

1000  IF((XB.LT.X(LOW).AND.XB.LT.X(IGH)).OR.
     1    (XB.GT.X(LOW).AND.XB.GT.X(IGH))) THEN
!	Extend the range
        IF(XB.GT.X(LOW).EQV.QASCND) THEN
!	Extend the range on higher side
          IF(IGH.GE.N) THEN
            IER=24
            LOW=N-1
          ELSE
            NIGH=MIN(N,IGH+2*(IGH-LOW))
            LOW=IGH
            IGH=NIGH
            GO TO 1000
          ENDIF
        ELSE
!	Extend the range on lower side
          IF(LOW.LE.1) THEN
            IER=25
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

!	Calculate the interpolant and its derivatives
      DX=XB-X(LOW)
      SPLEVL=((C(3,LOW)*DX+C(2,LOW))*DX+C(1,LOW))*DX+F(LOW)
      DFB=(3.*C(3,LOW)*DX+2.*C(2,LOW))*DX+C(1,LOW)
      DDFB=6.*C(3,LOW)*DX+2.*C(2,LOW)
      END
 
!     --------------------------------------------------------
 
!	To calculate coefficients of cubic spline interpolation with
!		not-a-knot boundary conditions
!
!	X : (input) Real array of length N containing x values
!	F : (input) Real array of length N containing values of function at X(I)
!		F(I) is the tabulated function value at X(I).
!	N : (input) Length of table X,F
!	C : (output) Real array of length 3*N containing the spline coefficients
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=201 implies that N<2
!
!	Required routines : None

      SUBROUTINE SPLINE(X,F,N,C,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N),C(3,N)

      IER=0
      IF(N.LE.1) THEN
        IER=201
        RETURN

      ELSE IF(N.EQ.2) THEN
!	Use linear interpolation
        C(1,1)=(F(2)-F(1))/(X(2)-X(1))
        C(2,1)=0.0
        C(3,1)=0.0
        RETURN

      ELSE IF(N.EQ.3) THEN
!	Use quadratic interpolation
        DIV12=(F(2)-F(1))/(X(2)-X(1))
        DIV23=(F(3)-F(2))/(X(3)-X(2))
        C(3,1)=0.0
        C(3,2)=0.0
        C(2,1)=(DIV23-DIV12)/(X(3)-X(1))
        C(2,2)=C(2,1)
        C(1,1)=DIV12+C(2,1)*(X(1)-X(2))
        C(1,2)=DIV23+C(2,1)*(X(2)-X(3))
        RETURN

      ELSE
!	Use cubic splines 

!	Setting up the coefficients of tridiagonal matrix
        C(3,N)=(F(N)-F(N-1))/(X(N)-X(N-1))
        DO 1000 I=N-1,2,-1
          C(3,I)=(F(I)-F(I-1))/(X(I)-X(I-1))
          C(2,I)=2.*(X(I+1)-X(I-1))
!	The right hand sides
1000    C(1,I)=3.*(C(3,I)*(X(I+1)-X(I))+C(3,I+1)*(X(I)-X(I-1)))

!	The not-a-knot boundary conditions
        C1=X(3)-X(1)
        C(2,1)=X(3)-X(2)
        C(1,1)=C(3,2)*C(2,1)*(2.*C1+X(2)-X(1))+C(3,3)*(X(2)-X(1))**2
        C(1,1)=C(1,1)/C1
        CN=X(N)-X(N-2)
        C(2,N)=X(N-1)-X(N-2)
        C(1,N)=C(3,N)*C(2,N)*(2.*CN+X(N)-X(N-1))
        C(1,N)=(C(1,N)+C(3,N-1)*(X(N)-X(N-1))**2)/CN

!	Solving the equation by Gaussian elimination
        G=(X(3)-X(2))/C(2,1)
        C(2,2)=C(2,2)-G*C1
        C(1,2)=C(1,2)-G*C(1,1)
        DO 2000 J=2,N-2
          G=(X(J+2)-X(J+1))/C(2,J)
          C(2,J+1)=C(2,J+1)-G*(X(J)-X(J-1))
2000    C(1,J+1)=C(1,J+1)-G*C(1,J)
        G=CN/C(2,N-1)
        C(2,N)=C(2,N)-G*(X(N-1)-X(N-2))
        C(1,N)=C(1,N)-G*C(1,N-1)

!	The back-substitution
        C(1,N)=C(1,N)/C(2,N)
        DO 3000 I=N-1,2,-1
3000    C(1,I)=(C(1,I)-C(1,I+1)*(X(I)-X(I-1)))/C(2,I)
        C(1,1)=(C(1,1)-C(1,2)*C1)/C(2,1)

!	Calculating the coefficients of cubic spline
        DO 4000 I=1,N-1
          C(2,I)=(3.*C(3,I+1)-2.*C(1,I)-C(1,I+1))/(X(I+1)-X(I))
          C(3,I)=(C(1,I)+C(1,I+1)-2.*C(3,I+1))/(X(I+1)-X(I))**2
4000    CONTINUE
!	Set the coefficients for interval beyond X(N) using continuity
!	of second derivative, although they may not be used.
        C(2,N)=C(2,N-1)+3*(X(N)-X(N-1))*C(3,N-1)
        C(3,N)=0.0
      ENDIF
      END
 
!     --------------------------------------------------------
 
!	To generate uniformly distributed random numbers in interval (0,1)
!
!	ISEED : (input/output) is an integer value used as the seed
!		It should be initialised to negative value before first call
!		and should not be modified between successive calls.
!
!	Required routines : None

      FUNCTION RANF(ISEED)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(M1=714025,IA1=1366,IC1=150889)
      PARAMETER(M2=214326,IA2=3613,IC2=45289)
      PARAMETER(M3=139968,IA3=3877,IC3=29573,ISH=43)
      DIMENSION RAN(ISH)
      SAVE
      DATA IFLG/0/

!	Initialise on first call or when ISEED<0
      IF(ISEED.LT.0.OR.IFLG.EQ.0) THEN
        IFLG=1
        RM1=1.D0/M1
        RM2=1.D0/M2

!	Seeds for the three random number generators
        IS1=MOD(-ISEED,M1)
        IS2=MOD(IA1*IS1+IC1,M1)
        IS3=MOD(IA2*IS2+IC2,M2)
        ISEED=1

!	Store ISH random numbers in the array RAN
        DO 1000 J=1,ISH
          IS1=MOD(IA1*IS1+IC1,M1)
          IS2=MOD(IA2*IS2+IC2,M2)
          RAN(J)=(FLOAT(IS1)+FLOAT(IS2)*RM2)*RM1
1000    CONTINUE
      ENDIF

      IS1=MOD(IA1*IS1+IC1,M1)
      IS2=MOD(IA2*IS2+IC2,M2)
      IS3=MOD(IA3*IS3+IC3,M3)
!	Select a random entry from RAN and store a new number in its place
      I=1+(ISH*IS3)/M3
      RANF=RAN(I)
      RAN(I)=(FLOAT(IS1)+FLOAT(IS2)*RM2)*RM1
      END
 
!     ----------------------------------------------------------
 
!     To calculate the function to be minimised.
!     The function value is the maximum error in spline approximation
!     using the selected knots.
!     IP is an array specifying the serial no. of knots to be used

      FUNCTION FCN(IP)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(49),IP(20),C(3,50),XP(20),FP(20)
      COMMON ERR(600),X0(600),F0(600),X(49)
      DATA IPAS/-999/
!	The input data for the Titanium problem
      DATA F/0.644D0,0.622D0,0.638D0,0.649D0,0.652D0,0.639D0,0.646D0,
     1       0.657D0,0.652D0,0.655D0,0.644D0,0.663D0,0.663D0,0.668D0,
     2       0.676D0,0.676D0,0.686D0,0.679D0,0.678D0,0.683D0,0.694D0,
     3       0.699D0,0.710D0,0.730D0,0.763D0,0.812D0,0.907D0,1.044D0,
     4       1.336D0,1.881D0,2.169D0,2.075D0,1.598D0,1.211D0,0.916D0,
     5       0.746D0,0.672D0,0.627D0,0.615D0,0.607D0,0.606D0,0.609D0,
     6       0.603D0,0.601D0,0.603D0,0.601D0,0.611D0,0.601D0,0.608D0/
      SAVE
 
51    FORMAT(F9.1,F9.3)
52    FORMAT('#    X',8X,'F(X)')
 
      IF(IPAS.EQ.-999) THEN
!	During the first call write out the input data in file titan.out
        NP=49
        WRITE(21,52)
        DO 1000 I=1,NP
          X(I)=575+10*I
          WRITE(21,51) X(I),F(I)
1000    CONTINUE

!	Calculate the cubic spline interpolation at 501 points using full table
        CALL SPLINE(X,F,NP,C,IER)
        N0=501
        H=(X(NP)-X(1))/(N0-1)
        DO 1200 I=1,N0
          XB=X(1)+H*(I-1)
          X0(I)=XB
          F0(I)= SPLEVL(XB,NP,X,F,C,DFB,DDFB,IER)
1200    CONTINUE
        IPAS=1
      ENDIF
 
      N=12
!	Setup the knots using array IP
      DO 1500 I=1,N
        XP(I)=X(IP(I))
        FP(I)=F(IP(I))
1500  CONTINUE
 
!	Calculate the cubic spline interpolation using only N chosen knots
!	and compare with full result to calculate the maximum error
      CALL SPLINE(XP,FP,N,C,IER)
      ERRMAX=0.0
      DO 2000 I=1,N0
        XB=X0(I)
        FB= SPLEVL(XB,N,XP,FP,C,DFB,DDFB,IER)
        DF=ABS(FB-F0(I))
        ERR(I)=FB-F0(I)
        IF(DF.GT.ERRMAX) ERRMAX=DF
2000  CONTINUE
      FCN=ERRMAX
      END
 
!     ----------------------------------------------------------

!	To change the set of knots selected for interpolation
 
      SUBROUTINE CHANGE(TEMP,IDUM,IP)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IP(20)
      EXTERNAL RANF
 
      N=12
      IF(TEMP.GT.0.05) THEN
        DO 2000 I=1,N
!	Shift the index of all knots by at most 1 randomly and
!	maintain the ordering of knots
          IP1=IP(I)+3*RANF(IDUM)-1
          IF(I.GT.1) THEN
            IF(IP1.LE.IP(I-1)+1) IP1=IP(I-1)+2
          ELSE
            IF(IP1.LT.1) IP1=1
          ENDIF
          IF(IP1.GT.49+2*I-2*N) IP1=49+2*I-2*N
          IP(I)=IP1
2000    CONTINUE
      ELSE

!	If the temperature is less than 0.05 then shift only one knot at a time
        I=(N)*RANF(IDUM)+1
        IF(RANF(IDUM).LT.0.5) THEN
          IP1=IP(I)-1
        ELSE
          IP1=IP(I)+1
        ENDIF
        IF(I.GT.1) THEN
          IF(IP1.LE.IP(I-1)) IP1=IP(I-1)+1
        ELSE
          IF(IP1.LT.1) IP1=1
        ENDIF
        IF(I.LT.N) THEN
          IF(IP1.GE.IP(I+1)) IP1=IP(I+1)-1
        ELSE
          IF(IP1.GT.49) IP1=49
        ENDIF
        IP(I)=IP1
      ENDIF
      END
