!       Calculate the incomplete beta function and incomplete gamma function

      PROGRAM IBETA
      IMPLICIT REAL*8(A-H,O-Z)

51    FORMAT('   X =',1PD14.6,'  A =',D14.6,'   B =',D14.6,/,
     1       '  Pgamma(A,X) =',D14.6,'    I_X(A,B) =',D14.6)

100   PRINT *,'TYPE X, A, B   (Quits when X < 0)'
      READ *,X,A,B
      IF(x.LT.0.0) STOP
      F1=GAMMAP(A,X)
      F2=BETAP(A,B,X)
      WRITE(6,51) X,A,B,F1,F2
      GO TO 100
      END

!       ----------------------------------------------------------------------

!	To integrate a function over finite interval using adaptive control
!	of step size
!
!	RINT : (output) Calculated value of the integral
!	XL : (input) The lower limit
!	XU : (input) The upper limit
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	F : (input) Name of the function routine to calculate the integrand
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=31 implies specified accuracy was not achieved on
!			at least one subinterval
!		IER=32 implies that this failure occurred more than IFMAX (=5) times
!		IER=325 implies that subroutine failed to attain required
!			accuracy using NMAX function evaluations
!		In all cases DIF will contain the estimated accuracy
!	NPT : (output) Number of function evaluations used by subroutine
!	NMAX : (input/output) Maximum number of function evaluations to be tried
!		If NMAX.LE.0 it is set to MAXPT (=100000)
!
!		FUNCTION F(X) must be supplied by the user.
!
!	Required routines : KRONRD (or GAUS16), F

      SUBROUTINE ADPINT(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(IPMAX=100,IFMAX=5,MAXPT=100000)
      EXTERNAL F
      DIMENSION XU1(IPMAX)

      IER=0
      IFAIL=0
      RINT=0.0
      DIF=0.0
      IF(XL.EQ.XU) RETURN
      IF(NMAX.LE.0) NMAX=MAXPT
      AEPSL=AEPS
      NPT=0
      RL=XL
      RU=XU
      IU=0

!	To evaluate the integral over [RL,RU]
1000  CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
!1000  CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RM=0.5*(RL+RU)
!	Q=.TRUE. if the interval cannot be divided further
      Q=IU.GE.IPMAX.OR.RM.EQ.RL.OR.RM.EQ.RU

      IF(DIF0.LT.MAX(ABS(FINT)*REPS,AEPSL).OR.Q) THEN
!	Accept the value of FINT if adequate convergence or if the interval
!	cannot be subdivided further
        RINT=RINT+FINT
        DIF=DIF+DIF0
        IF(Q.AND.DIF0.GT.MAX(ABS(RINT)*REPS,AEPSL)) THEN
!	Integration fails to converge on this subinterval. Go to the next subinterval
          IER=31
          IFAIL=IFAIL+1
          IF(IFAIL.GT.IFMAX) THEN
!	If failure is frequent then adjust the convergence criterion.
            IER=32
            AEPSL=DIF*0.5
          ENDIF
        ENDIF

!	If all subintervals are exhausted then return
        IF(IU.LE.0) RETURN

!	otherwise try next subinterval
        RL=RU
        RU=XU1(IU)
        IU=IU-1
      ELSE

!	Subdivide the current interval and try again
        IU=IU+1
        XU1(IU)=RU
        RU=RM
      ENDIF

      IF(NPT.LT.NMAX) GO TO 1000
!	If the number of function evaluations has exceeded the limit then
!	try a last call to estimate the integral over the remaining interval
      IER=325
      RU=XU
      CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
!      CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END

!       ----------------------------------------------------------------------

!	To calculate the incomplete Beta function I_x(a,b) using
!	the integral, called by BETAP
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!               Beta function
!
!	Required routines : GAMMAL, ADPINT, KRONRD, FBETA

      FUNCTION BETAI(A,B,X)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FBETA
      COMMON/BETA/AA,BB

      AA=A
      BB=B
      XL=0
      XU=X
      REPS=1.D-14
      AEPS=1.D-280
      CALL ADPINT(RINT,XL,XU,REPS,AEPS,DIF,FBETA,IER,NPT,NMAX)
      B1=LOG(RINT)+GAMMAL(A+B)-GAMMAL(A)-GAMMAL(B)
      BETAI=EXP(B1)
      IF(BETAI.GT.1.0) BETAI=1.0
      END

!       ----------------------------------------------------------------------

!	To calculate the incomplete Beta function I_x(a,b) for
!	positive real arguments
!	It returns a value of -1 for a.le.0 or b.le.0 or x<0 or x>1
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!                       Beta function
!
!	Required routines : GAMMAL, BETAI, BETSER, BETCON, BETCON1,
!			    ADPINT, KRONRD, FBETA

      FUNCTION BETAP(A,B,X)
      IMPLICIT REAL*8(A-H,O-Z)

51    FORMAT(' ** Error in evaluating Incomplete Beta Function at',/,
     1  ' A=',1PD12.4,' B=',D12.4,' X=',D12.4,'  BETAP=',D12.4)

      IF(A.LE.0.OR.B.LE.0.OR.X.LT.0.OR.X.GT.1) THEN
        BETAP=-1
        RETURN
      ENDIF
      AMAX=MAX(A,B)
      AMIN=MIN(A,B)
      BETAL=(A+B)*LOG(A+B)-A*LOG(A)-B*LOG(B)
      IF(AMAX.LE.30) THEN
        IF(X.LE.0.5) THEN
          BETAP=BETSER(A,B,X)
        ELSE
          BETAP=BETSER(B,A,1-X)
          BETAP=1-BETAP
        ENDIF
      ELSE IF(B.LE.20.AND.X.LE.0.71D0) THEN
        BETAP=BETSER(A,B,X)
      ELSE IF(A.LE.20.AND.X.GE.0.3D0) THEN
        BETAP=BETSER(B,A,1-X)
        BETAP=1-BETAP
      ELSE IF(B.LE.50.AND.X.LE.0.35D0) THEN
        BETAP=BETSER(A,B,X)
      ELSE IF(A.LE.50.AND.X.GE.0.65D0) THEN
        BETAP=BETSER(B,A,1-X)
        BETAP=1-BETAP
      ELSE IF(B.LE.100.AND.X.LE.0.18D0) THEN
        BETAP=BETSER(A,B,X)
      ELSE IF(A.LE.100.AND.X.GE.0.82D0) THEN
        BETAP=BETSER(B,A,1-X)
        BETAP=1-BETAP
      ELSE IF(B.LE.180.AND.X.LE.0.1D0) THEN
        BETAP=BETSER(A,B,X)
      ELSE IF(A.LE.180.AND.X.GE.0.9D0) THEN
        BETAP=BETSER(B,A,1-X)
        BETAP=1-BETAP
      ELSE IF(X.LT.0.5) THEN
        IF(A.LT.2.D0) THEN
          BETAP=BETCON(A,B,X)
        ELSE IF(BETAL.GT.700) THEN
          BETAP=BETCON1(A,B,X)
        ELSE
          BETAP=BETAI(A,B,X)
        ENDIF
      ELSE
        IF(B.LT.2.D0) THEN
          BETAP=BETCON(B,A,1-X)
        ELSE IF(BETAL.GT.700) THEN
          BETAP=BETCON1(B,A,1-X)
        ELSE
          BETAP=BETAI(B,A,1-X)
        ENDIF
        BETAP=1-BETAP
      ENDIF
      
      IF(BETAP.LT.0.0.OR.BETAP.GT.1) PRINT 51,A,B,X,BETAP
      END

!       ----------------------------------------------------------------------

!	To calculate the incomplete Beta function I_x(a,b) using
!	the continued fraction (modified form), called by BETAP
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!                       Beta function
!
!	Required routines : GAMMAL

      FUNCTION BETCON(A,B,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(500),D(500)
           
51    FORMAT(' ** Roundoff error while evaluating the modified ',
     1  'continued fraction for incomplete Beta function at',/,'  A=',
     2  1PD12.4,'  B=',D12.4,'  X=',D12.4,'  CON. FRAC.=',D12.4)

      C(1)=A
      D(1)=A*(1-X*(A+B)/(A+1))
      B1=1+X*(B-1)/(A+1)-X*(A+1)*(A+B+1)/(A+3)
      C(2)=A*B1
      A1=A*X*X*(B-1)*(A+B)/(A+1)**2
      D(2)=B1*D(1)+A1
      C1=C(2)/D(2)
      DO I=2,499
        D1=-X*(A+I-1)*(A+B+I-1)/((A+2*I-2)*(A+2*I-1))
        D3=-X*(A+I)*(A+B+I)/((A+2*I)*(A+2*I+1))
        D2=X*I*(B-I)/((A+2*I-1)*(A+2*I))
        C(I+1)=C(I)*(A+2*I)*(1+D2+D3)-C(I-1)*(A+2*I)*(A+2*I-2)*D2*D1
        D(I+1)=D(I)*(A+2*I)*(1+D2+D3)-D(I-1)*(A+2*I)*(A+2*I-2)*D2*D1
!	Scale the numerator and denominator to prevent underflow/overflow
        IF(ABS(C(I+1)).GT.1.D200) THEN
          C(I+1)=C(I+1)/1.D200
          D(I+1)=D(I+1)/1.D200
          C(I)=C(I)/1.D200
          D(I)=D(I)/1.D200
        ENDIF
        IF(ABS(C(I+1)).LT.1.D-200) THEN
          C(I+1)=C(I+1)*1.D200
          D(I+1)=D(I+1)*1.D200
          C(I)=C(I)*1.D200
          D(I)=D(I)*1.D200
        ENDIF
        C2=C(I+1)/D(I+1)
        IF(ABS(C2-C1).LT.1.D-12) EXIT
        C1=C2
      ENDDO
      IF(C2.LT.0.0) PRINT 51,A,B,X,C2
      B1=A*LOG(X)+B*LOG(1-X)+LOG(C2)-LOG(A)
      B1=B1+GAMMAL(A+B)-GAMMAL(A)-GAMMAL(B)
      BETCON=EXP(B1)
      END

!       ----------------------------------------------------------------------

!	To calculate the incomplete Beta function I_x(a,b) using
!	the continued fraction, called by BETAP
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!                       Beta function
!
!	Required routines : GAMMAL

      FUNCTION BETCON1(A,B,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION c(500),d(500)
      
51    FORMAT(' ** Roundoff error while evaluating the ',
     1  'continued fraction for incomplete Beta function at',/,'  A=',
     2  1PD12.4,'  B=',D12.4,'  X=',D12.4,'  CON. FRAC.=',D12.4)

      C(1)=1
      D(1)=1
      C(2)=1
      D(2)=1-X*(A+B)/(A+1)
      C1=C(2)/D(2)
      DO I=2,498,2
        M=I/2
        C(I+1)=C(I)+M*(B-M)*X*C(I-1)/((A+I-1)*(A+I))
        D(I+1)=D(I)+M*(B-M)*X*D(I-1)/((A+I-1)*(A+I))
        C(I+2)=C(I+1)-(A+M)*(A+B+M)*X*C(I)/((A+I+1)*(A+I))
        D(I+2)=D(I+1)-(A+M)*(A+B+M)*X*D(I)/((A+I+1)*(A+I))
!	Scale the numerator and denominator to prevent underflow/overflow
        IF(ABS(C(I+2)).GT.1.D200) THEN
          C(I+1)=C(I+1)/1.D200
          D(I+1)=D(I+1)/1.D200
          C(I+2)=C(I+2)/1.D200
          D(I+2)=D(I+2)/1.D200
        ENDIF
        IF(ABS(C(I+2)).LT.1.D-200) THEN
          C(I+1)=C(I+1)*1.D200
          D(I+1)=D(I+1)*1.D200
          C(I+2)=C(I+2)*1.D200
          D(I+2)=D(I+2)*1.D200
        ENDIF
        C2=C(I+2)/D(I+2)
        IF(ABS(C2-C1).LT.1.D-12) EXIT
        C1=C2
      ENDDO
      IF(C2.LT.0.0) PRINT 51,A,B,X,C2
      B1=A*LOG(X)+B*LOG(1-X)+LOG(C2)-LOG(A)
      B1=B1+GAMMAL(A+B)-GAMMAL(A)-GAMMAL(B)
      BETCON1=EXP(B1)
      END

!       ----------------------------------------------------------------------

!	To calculate the incomplete Beta function I_x(a,b) using
!	the infinite series, called by BETAP
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!                       Beta function
!
!	Required routines : GAMMAL

      FUNCTION BETSER(A,B,X)
      IMPLICIT REAL*8(A-H,O-Z)

51    FORMAT(' ** Roundoff error while evaluating the infinite series',
     1  ' for incomplete Beta function at',/,'  A=',1PD12.4,'  B=',
     2  D12.4,'  X=',D12.4,'  SUM=',D12.4,'  MAX TERM=',D12.4)

      S3=1
      T=A
      TMAX=A
      DO I=1,500,2
        T1=T*X*(I-B)/I
        T=T1*X*(I+1-B)/(I+1)
        S3=S3+(T1/(A+I)+T/(A+I+1))
        IF(ABS(T/(A+I+1)).GT.TMAX) TMAX=ABS(T/(A+I+1))
        IF(ABS(T/(A+I+1)).LT.1.D-15) EXIT
      ENDDO
      IF((S3).LT.1.D-16*TMAX) PRINT 51,A,B,X,S3,TMAX
      B1=A*LOG(X)+LOG(S3)-LOG(A)
      B1=B1+GAMMAL(A+B)-GAMMAL(A)-GAMMAL(B)
      BETSER=EXP(B1)
      END

!       ----------------------------------------------------------------------

!	To calculate the integrand for incomplete beta function
!       It is used by function BETAI
!
!	A,B : Arguments for the complete Beta function passed through
!               a common block
!	X : (input) Upper limit of integration defining the incomplete
!               Beta function
!
!	Required routines : none

      FUNCTION FBETA(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/BETA/A,B
      FBETA=X**(A-1)*(1-X)**(B-1)
      END

!       ----------------------------------------------------------------------

!	To calculate Gamma function for any real value of XG
!	Use GAMMAL for calculating the logarithm of Gamma function
!	which may be useful for large arguments or when argument is
!	close to a negative integer.
!
!	Required routines : None
 
      FUNCTION GAMMA(XG)
      IMPLICIT REAL*8(A-H,O-Z)
!	PIS is SQRT(2*PI)
      PARAMETER(PI=3.14159265358979323846D0,PIS=2.5066282746310005024D0)
      DIMENSION A(2),B(3),A1(6),B1(7)

!	The coefficients for rational function approximations
      DATA A/1.767971449569122937D+00,  2.909421117928672645D-01/
      DATA B/8.333333333333231537D-02,  1.445531763554246280D-01,
     1       2.012779361583001035D-02/
      DATA A1/3.905731686764559737D+03,  2.204952264401381785D+03,
     1       -1.932467485468849660D+03,  4.643360871045442213D+02,
     1       -4.818088806916028754D+01,  1.896853765546068169D+00/
      DATA B1/3.918055655523400310D+03, -1.088116266563809683D+02,
     1        8.203258626193993149D+02, -9.289402000761705906D+01,
     1        6.521113026294866877D+01, -6.090618615608719044D+00,
     1        1.475909104740280784D+00/
 
      X=ABS(XG)
      IF(X.GT.1000.0) THEN
!	Use asymptotic approximation (Stirling formula)
        GX=(1+1.D0/(12*X)+1./(288*X*X)
     1         -139/(51840*X**3)-571./(2488320D0*X**4))
        GAMMA=X**(X-0.5)*EXP(-X)*PIS*GX

      ELSE IF(X.GT.8.0) THEN
!	Use rational function approximation for Log(Gamma) 
        Y=1./X**2
        RMK=((B(3)*Y+B(2))*Y+B(1))/((A(2)*Y+A(1))*Y+1)
        GAMMA=X**(X-0.5)*EXP(-X)*PIS*EXP(RMK/X)

      ELSE IF(X.GE.2.0) THEN
!	Use rational function approximation for (Gamma) over [2,3]
!	after translating the range if necessary
        F1=1.0
        X1=X
2500    IF(X1.LE.3) GO TO 3000
        F1=F1*(X1-1)
        X1=X1-1
        GO TO 2500
3000    IF(X1.EQ.3) THEN
          GAMMA=F1*2
        ELSE IF(X1.EQ.2) THEN
          GAMMA=F1
        ENDIF

        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1          B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1          A1(1))*X1+1
        GAMMA=F1*FN/FD
      ELSE IF(X.GT.0.0) THEN
!	Use rational function approximation for (Gamma) over [2,3]
!	after translating the range if necessary
        F1=1./X
        X1=X+1
        IF(X.LT.1) THEN
          F1=F1/X1
          X1=X1+1
        ENDIF
        IF(X1.EQ.2) GAMMA=F1

        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1          B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1          A1(1))*X1+1
        GAMMA=F1*FN/FD
 
      ENDIF

      IF(XG.GT.0.0) RETURN
      IX=X
      IF(X.GT.IX) THEN
        GAMMA=PI/(XG*SIN(PI*X)*GAMMA)
      ELSE
        GAMMA=(-1)**IX/0.0
      ENDIF
 
      END

!       ----------------------------------------------------------------------

!	To calculate the Logarithm of Gamma function for a real argument
!	For negative values it give ln(abs(Gamma(x)))
!
!	Required routines : None
 
      FUNCTION GAMMAL(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(4),B(4),A1(6),B1(6)
!	PI2L=LOG(2*PI)/2
      PARAMETER (PI2L=0.918938533204672741D0,PI=3.141592653589793238D0)

!	The coefficients of rational function approximations
      DATA A/ 1.013142782275024216D-2,  7.645657825398191944D-1,
     1        3.381172379819461227D-4,  1.595363637547538209D-2/
      DATA B/0.0, 8.333333333338911768D-2, 8.442856404442060242D-4,
     1            6.093603832366013704D-2/

      DATA A1/ 4.681163846241230144D0,  3.208225429683256526D0,
     1         5.145525793448859216D-1, 1.581117883959157936D-2,
     2        -6.398416804905407512D-5, 5.264566254181773919D-7/
      DATA B1/ 2.938038561191284576D0,  1.489364948862436743D0,
     1        -5.466291543917642961D0,  1.972497734170110410D-1,
     2         7.830146473241555157D-1, 5.756753067834747499D-2/
 
      T=ABS(X)
      IT=T
      IF(T.GE.10.0) THEN
!	Use asymptotic approximation for T>10
        Y=1/T
        FN=((B(4)*Y+B(3))*Y+B(2))*Y+B(1)
        FD=(((A(4)*Y+A(3))*Y+A(2))*Y+A(1))*Y+1
        GAMMAL=FN/FD-T+PI2L+(T-0.5)*LOG(T)
      ELSE
 
!	Use approximation for [4,5]
        FAC=0.0
        T1=T
        IF(T.LT.4) THEN
          DO I=IT,3
            FAC=FAC-LOG(T1)
            T1=T1+1
          ENDDO
        ELSE IF(T.GT.5) THEN
          DO I=IT,5,-1
            T1=T1-1
            FAC=FAC+LOG(T1)
          ENDDO
        ENDIF
        FN=((((B1(6)*T1+B1(5))*T1+B1(4))*T1+B1(3))*T1+B1(2))*T1+B1(1)
        FD=(((((A1(6)*T1+A1(5))*T1+A1(4))*T1+A1(3))*T1+A1(2))*T1+
     1          A1(1))*T1+1
        GAMMAL=FN/FD+FAC
      ENDIF

      IF(X.LE.0) THEN
        IF(T.GT.IT) THEN
          GAMMAL=LOG(PI)-LOG(ABS(SIN(PI*T)*T))-GAMMAL
        ELSE
          GAMMAL=(-1)**IT/0.0
        ENDIF
      ENDIF
 
      END

!       ----------------------------------------------------------------------

!	To calculate the incomplete Gamma function P(a,x) for
!	positive real arguments
!	It returns a value of -1 for a.le.0 or x<0
!
!	A : (input) Argument for the complete Gamma function
!	X : (input) Upper limit of integration defining the incomplete
!                       Gamma function
!
!	Required routines : GAMMA, GAMMAL

      FUNCTION GAMMAP(A,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(500),B(500)

      IF(A.LE.0.OR.X.LT.0) THEN
        GAMMAP=-1
        RETURN
      ENDIF

      IF(X.LT.3) THEN
!	Use power series
        S1=1
        T=A
        DO I=1,500
          T=-T*X/I
          S1=S1+T/(A+I)
          IF(ABS(T/(A+I)).LT.1.D-14) EXIT
        ENDDO
        IF(A.LT.140) THEN
          GAMMAP=X**A*S1/GAMMA(A+1)
        ELSE
          GAMMAP=0.0
        ENDIF
      ELSE IF(A.LT.1.2D0*X) THEN
!	Use continued fraction
        C(1)=1
        B(1)=X
        C(2)=1
        B(2)=X+1-A
        C1=C(2)/B(2)
        DO I=2,498,2
          C(I+1)=X*C(I)+(I/2)*C(I-1)
          B(I+1)=X*B(I)+(I/2)*B(I-1)
          C(I+2)=C(I+1)+(I/2+1-A)*C(I)
          B(I+2)=B(I+1)+(I/2+1-A)*B(I)
          IF(ABS(B(I+2)).GT.1.D200) THEN
            C(I+1)=C(I+1)/1.D200
            B(I+1)=B(I+1)/1.D200
            C(I+2)=C(I+2)/1.D200
            B(I+2)=B(I+2)/1.D200
          ENDIF
          IF(ABS(B(I+2)).LT.1.D-200) THEN
            C(I+1)=C(I+1)*1.D200
            B(I+1)=B(I+1)*1.D200
            C(I+2)=C(I+2)*1.D200
            B(I+2)=B(I+2)*1.D200
          ENDIF
          C2=C(I+2)/B(I+2)
          IF(ABS(C2-C1).LT.1.D-12) EXIT
          C1=C2
        ENDDO
        G1=-X+A*LOG(X)+LOG(C2)-GAMMAL(A)
        GAMMAP=1-EXP(G1)
      ELSE
!	Use the power series for a>x
        S2=1
        T=1
        DO I=1,500
          T=T*X/(A+I)
          S2=S2+T
          IF(ABS(T).LT.1.D-14) EXIT
        ENDDO
        G1=-X+A*LOG(X)+LOG(S2)-GAMMAL(A+1)
        GAMMAP=EXP(G1)
      ENDIF
      END

!       ----------------------------------------------------------------------

!	To integrate a function over a finite interval using Gauss-Kronrod formula
!	For use with ADPINT
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	N : (output) Number of function evaluations used by subroutine
!	F : (input) Name of the function routine to calculate the integrand
!
!	FUNCTION F(X) must be supplied by the user
!
!	Required routines : F

      SUBROUTINE KRONRD(RI,A,B,DIF,N,F)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  W7(4),A7(4),WK7(4),WK15(4),AK15(4)

!	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
!	WK7 are the weights for these points in Kronrod formula
!	WK15 and AK15 are the weights and abscissas for the remaining points
!	in Kronrod formula.
!	Because of symmetry only half the points are given.

      DATA W7  /0.12948496616886969327D0, 0.27970539148927666790D0,
     *          0.38183005050511894495D0, 0.41795918367346938775D0/
      DATA A7  /0.94910791234275852452D0, 0.74153118559939443986D0,
     *          0.40584515137739716690D0, 0.0/
      DATA WK7 /0.06309209262997855329D0, 0.14065325971552591874D0,
     *          0.19035057806478540991D0, 0.20948214108472782801D0/
      DATA WK15/0.02293532201052922496D0, 0.10479001032225018383D0,
     *          0.16900472663926790282D0, 0.20443294007529889241D0/
      DATA AK15/0.99145537112081263920D0, 0.86486442335976907278D0,
     *          0.58608723546769113029D0, 0.20778495500789846760D0/

      AT=(B-A)/2.
      BT=(B+A)/2.
      FBT=F(BT)
      R1=W7(4)*FBT
      RI=WK7(4)*FBT
      DO 2000 K=1,3
        F1=F(AT*A7(K)+BT)
        F2=F(BT-AT*A7(K))
!	7-point Gauss-Legendre formula
        R1=R1+W7(K)*(F1+F2)
!	15-point Kronrod formula
        RI=RI+WK7(K)*(F1+F2)
2000  CONTINUE

      DO 2500 K=1,4
2500  RI=RI+WK15(K)*(F(AT*AK15(K)+BT)+F(BT-AT*AK15(K)))

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=15
      END
