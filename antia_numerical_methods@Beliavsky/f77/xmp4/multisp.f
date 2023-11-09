!     PROGRAM TO CALCULATE MULTIPLE INTEGRALS OVER A HYPERSPHERE IN
!     N DIMENSIONS, THE INTEGRAND IS SPECIFIED IN TERMS OF THE CARTESIAN
!     COORDINATES.
!     THIS PROGRAM USES PRODUCT GAUSS RULES OR COMPOSITE MONOMIAL RULES
!     OR MONTE CARLO METHOD OR EQUIDISTRIBUTED SEQUENCES TO CALCULATE THE
!     INTEGRAL

      PROGRAM MULTI
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=20)
      EXTERNAL SPHND
      DIMENSION A(NMAX),B(NMAX),INT(NMAX),M(NMAX)

!     EXAMPLE 6.15 : INTEGRAL I5
 
51    FORMAT('   IER =',I4,5X,'IT =',I2,5X,'NO. OF FUNCTION',
     1       ' EVALUATIONS =',I8)
52    FORMAT('  INTEGRAL OVER',I3,' DIMENSIONS =',1PE14.6,5X,
     1       'ERROR ESTIMATE =',E14.6)
53    FORMAT('   IER =',I4,5X,'IT =',I2,5X,'M =',I2,5X,
     1       'NO. OF FUNCTION EVALUATIONS =',I8)
54    FORMAT('  INTEGRAL OVER',I3,' DIMENSIONS  S1 =',1PE14.6,5X,
     1       'S2 =',E14.6/9X,'ERROR ESTIMATE =',E14.6)
 
      REPS=1.E-6
      AEPS=1.E-8

100   PRINT *,'TYPE N=NO. OF DIMENSIONS,  IT=1/2/3/4'
      PRINT *,'IT=1  FOR PRODUCT RULES     (MULINT)'  
      PRINT *,'IT=2  FOR MONOMIAL RULES    (STRINT)'
      PRINT *,'IT=3  FOR MONTE CARLO METHOD  (MCARLO)'
      PRINT *,'IT=4  FOR EQUIDISTRIBUTED SEQUENCES (EQUIDS)'
      PRINT *,'        (QUITS WHEN IT.LT.1.OR.IT.GT.4)'

      READ *,N,IT
      IF(IT.LT.1.OR.IT.GT.4) STOP
      PI=2.*ACOS(0.0D0)

!	Set up the limits in each dimension in hyper-spherical coordinates
      DO 500 I=1,NMAX
        A(I)=0.0
        B(I)=PI
        INT(I)=0
        M(I)=0
500   CONTINUE
      B(1)=1.0
      B(N)=2.*PI

      MAXPT=1000000 
      IF(IT.EQ.1) THEN 
        CALL MULINT(A,B,N,M,INT,SPHND,RI,REPS,AEPS,DIF,IER,NUM,MAXPT)
        WRITE(6,51) IER,IT,NUM
        WRITE(6,52) N,RI,DIF
      ELSE IF(IT.EQ.2) THEN
        PRINT *,'TYPE M=1,3 OR 5=DEGREE OF MONOMIAL RULE'
        READ *,MM
        CALL STRINT(A,B,N,MM,INT,SPHND,RI,REPS,AEPS,DIF,IER,NUM,MAXPT)
        WRITE(6,53) IER,IT,MM,NUM
        WRITE(6,52) N,RI,DIF
      ELSE IF(IT.EQ.3) THEN
        CALL MCARLO(A,B,N,MAXPT,SPHND,RI,REPS,AEPS,DIF,NUM,IER)
        WRITE(6,51) IER,IT,NUM
        WRITE(6,52) N,RI,DIF
      ELSE
        CALL EQUIDS(A,B,N,MAXPT,SPHND,S1,S2,REPS,AEPS,DIF,NUM,IER)
        WRITE(6,51) IER,IT,NUM
        WRITE(6,54) N,S1,S2,DIF
      ENDIF
 
      GO TO 100
      END
 
!     -----------------------------------------------------------
 
!	Multiple integration over a hyper-rectangle in n-dimensions
!	using equidistributed sequences
!	Because of large roundoff error in summation it may not be
!	advisable to use this routine in single precision.
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions
!	NPT : (input) Maximum number of function evaluations to be used
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION F(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	S1 : (output) The calculated value of the integral
!	S2 : (output) Another approximation to the value of the integral
!		For smooth functions S2 is expected to be better approximation
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(S2))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	NP : (output) Number of function evaluations used by subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=39 implies specified accuracy was not achieved in
!			which case DIF will contain the estimated accuracy
!		IER=312 implies N<1 or N>21 and no calculations are done
!
!	FUNCTION F(N,X) must be supplied by the user
!	
!	Required routines :  F
!
      SUBROUTINE EQUIDS(A,B,N,NPT,F,S1,S2,REPS,AEPS,DIF,NP,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NCHK=100,NMAX=21)
      DIMENSION A(N),B(N),H(NMAX),XA(NMAX),WT(NMAX),AT(NMAX)
!	The first 21 prime numbers
      DATA AT/2.,3.,5.,7.,11.,13.,17.,19.,23.,29.,31.,37.,41.,
     1        43.,47.,53.,59.,61.,67.,71.,73./

      IER=312
      S1=F(N,A)
      S2=S1
      NP=1
      IF(N.GT.NMAX.OR.N.LT.1) RETURN

      IER=0
      HH=1.0
      DO 1000 I=1,N
        H(I)=B(I)-A(I)
        HH=HH*H(I)
!	The irrational numbers for generating equidistributed sequences
        WT(I)=SQRT(AT(I))
1000  CONTINUE

      RI=0.0
      RI1=0.0
      DIF=0.0
      NPT1=NCHK
      DO 2000 I=1,NPT
!	Generate the abscissas using equidistributed sequences
        DO 1500 J=1,N
          A1=I*WT(J)
          A1=2.*ABS(A1-INT(A1+0.5))*H(J)
          XA(J)=A(J)+A1
1500    CONTINUE
!	Accumulate the sum
        S1=S1+2.*(F(N,XA)-RI)
        S2=S2+S1

        IF(MOD(I,NCHK).EQ.0) THEN
!	To control the roundoff error form partial sums
          SS1=RI+S1/(2*I+1)
          DIFF=S2/((I+1.)**2)
          S2=0.0
          S1=S1-(2*I+1)*DIFF
!	The new approximation to the average value of function
          RI=RI+DIFF

          IF(I.EQ.NPT1) THEN
!	Check for convergence
            DIF1=DIF
            DIF=ABS(RI-RI1)
            RI1=RI
            IF(DIF+DIF1.LT.MAX(AEPS/HH,ABS(RI)*REPS).AND.I.GT.5*NCHK)
     1         THEN
              S1=SS1*HH
              S2=RI*HH
              DIF=(DIF+DIF1)*HH
              NP=I+1
              RETURN
            ENDIF

            NPT1=NPT1*2
          ENDIF
        ENDIF
2000  CONTINUE

!	Integral fails to converge
      IER=39
      S1=(RI+S1/(2*NPT+1))*HH
      S2=(RI+S2/((NPT+1.)**2))*HH
      DIF=(DIF+DIF1)*HH
      NP=NPT+1
      END
 
!     -------------------------------------------------
 
!	Multiple integration over a hyper-rectangle in n-dimensions
!	using Monte-Carlo technique
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions
!	NPT : (input) Maximum number of function evaluations to be used
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	RI : (output) The calculated value of the integral
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	ERR : (output) estimated (absolute) error achieved by the subroutine
!		It is 2.576 times the estimated standard deviation
!		divided by SQRT(NP).
!	NP : (output) Number of function evaluations used by subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=39 implies specified accuracy was not achieved in
!			which case ERR will contain the estimated accuracy
!		IER=311 implies N<1 or N>50 and no calculations are done
!
!	FUNCTION F(N,X) should be supplied by the user
!	
!	Required routines : RANF, F
!
      SUBROUTINE MCARLO(A,B,N,NPT,F,RI,REPS,AEPS,ERR,NP,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=50,NCHK=100)
!	To override the internal RANF if it exists
      EXTERNAL RANF
      DIMENSION A(N),B(N),H(NMAX),XA(NMAX)

      IER=311
      RI=0.0
      NP=0
      IF(N.GT.NMAX.OR.N.LT.1) RETURN

      IER=0
      HH=1.0
      DO 1000 I=1,N
        H(I)=B(I)-A(I)
1000  HH=HH*H(I)

      RI1=0.0
      VAR1=0.0
!	Seed for random number generator, should be changed if another routine is used
      ISEED=-12345
      NPT1=NCHK

      DO 2000 I=1,NPT
!	Generating the abscissas
        DO 1500 J=1,N
1500    XA(J)=A(J)+H(J)*RANF(ISEED)
        F1=F(N,XA)
        RI1=RI1+F1
        VAR1=VAR1+F1*F1

        IF(I.EQ.NPT1) THEN
!	Compute intermediate sums to check for convergence
          RI=RI1/I
          VAR=VAR1/I-RI*RI
          IF(VAR.LT.0.0) VAR=0.0
          ERR=2.576D0*HH*SQRT(VAR/NPT1)
          RI=RI*HH
          NP=I
          IF(ERR.LT.MAX(AEPS/HH,ABS(RI)*REPS)) RETURN
          NPT1=2*NPT1
        ENDIF
2000  CONTINUE

!	Integral fails to converge
      RI=RI1/NPT
      VAR=VAR1/NPT-RI*RI
      ERR=2.576D0*HH*SQRT(VAR/NPT)
      RI=RI*HH
      IF(ERR.GT.MAX(AEPS,ABS(RI)*REPS)) IER=39
      NP=NPT
      END
 
!     --------------------------------------------------

!	Multiple integration over a hyper-rectangle in n-dimensions
!	using product Gauss-Legendre formulas
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions
!	M : (input/output) Integer array of length N specifying the formula
!		to be used along each dimension. M(J)-point formula will
!		be used along Jth dimension, M(J) should be 2,4,8,16 or 32
!		otherwise it will be set to a default value of 2
!	IND : (input/output) Integer array of length N specifying the number
!		of subintervals to be used along each dimension. IND(J)>0
!		otherwise it will be set to a default value of 1
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	RINT : (output) The calculated value of the integral
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=39 implies specified accuracy was not achieved in
!			which case DIF will contain the estimated accuracy
!		IER=305 implies N<1 or N>20 and no calculations are done
!		IER=307 implies that number of points exceeded MAXPT in
!			first attempt and no approximation of RINT is calculated
!	NUM : (output) Number of function evaluations used by subroutine
!	MAXPT : (input/output) Maximum number of function evaluations permitted
!		If MAXPT <1 it is set to a default value of MAXPTS (=1100000)
!
!	FUNCTION F(N,X) must be supplied by the user
!	
!	Required routines : NGAUSS, F
!
      SUBROUTINE MULINT(A,B,N,M,IND,F,RINT,REPS,AEPS,DIF,IER,NUM,MAXPT)
      IMPLICIT LOGICAL(Q)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL F
      PARAMETER(MAXPTS=1100000)
      DIMENSION A(N),B(N),M(N),IND(N)

!	Set M(I) and IND(I) to default values if they are unacceptable
      DO 1000 I=1,N
        IF(M(I).NE.2.AND.M(I).NE.4.AND.M(I).NE.8.AND.M(I).NE.16.
     1     AND.M(I).NE.32) M(I)=2
        IF(IND(I).LT.1) IND(I)=1
1000  CONTINUE
      IF(MAXPT.LE.0) MAXPT=MAXPTS
      NUM=0

!	Evaluate the integral
      CALL NGAUSS(A,B,N,M,IND,F,RINT,IER,NO,MAXPT)
      NUM=NUM+NO
      IF(IER.GT.100) RETURN
      IER=39

!	Iteration to check and improve the accuracy of integral
      DO 3000 I=1,10
        QC=.TRUE.
        DIF=0.0

!	Check for convergence along each dimension by doubling the
!	number of points at which function is evaluated
        DO 2000 J=1,N
          IF(NUM.GT.MAXPT) RETURN
          M1=M(J)
          I1=IND(J)
          IF(M(J).LT.32) THEN
!	If M(J)<32 then double the order of formula
            M(J)=2*M(J)
          ELSE
!	otherwise double the number of subintervals
            IND(J)=2*IND(J)
          ENDIF

          CALL NGAUSS(A,B,N,M,IND,F,RINT1,IER1,NO,MAXPT)
          NUM=NUM+NO
          DIF=DIF+ABS(RINT-RINT1)
          IF(IER1.GT.100) RETURN

          IF(ABS(RINT1-RINT).LT.MAX(AEPS,ABS(RINT1)*REPS)) THEN
!	If satisfactory accuracy is achieved then revert back to old
!	values of M(J) and IND(J)
            M(J)=M1
            IND(J)=I1
          ELSE
!	otherwise use new values
            RINT=RINT1
            QC=.FALSE.
          ENDIF
2000    CONTINUE

        IF(QC) THEN
!	If satisfactory accuracy is achieved for all dimensions then return
          IER=0
          RETURN
        ENDIF
!	If the number of function evaluations exceeds MAXPT then quit
        IF(NUM.GT.MAXPT) RETURN
3000  CONTINUE
      END

!     --------------------------------------------------

!	Multiple integration over a hyper-rectangle in n-dimensions
!	using product Gauss-Legendre formulas with given number of points
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions, N>0 and N<NMAX (=21)
!	M : (input) Integer array of length N specifying the formula
!		to be used along each dimension. M(J)-point formula will
!		be used along Jth dimension, M(J) should be 2,4,8,16 or 32
!		otherwise IER is set to 306 and no calculations are done
!	IND : (input) Integer array of length N specifying the number
!		of subintervals to be used along each dimension. IND(J)>0
!		otherwise IER is set to 306 and no calculations are done
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	RI : (output) The calculated value of the integral
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=305 implies N<1 or N.GE.NMAX, in which case no calculations are done
!		IER=306 implies M(J) is not 2,4,8,16 or 32 or IND(J)<1 for some J
!			in which case no calculations are done
!		IER=307 implies that number of points exceeded MAXPT and
!			no calculations are done
!	NUM : (output) Number of function evaluations used by subroutine
!	MAXPT : (input/output) Maximum number of function evaluations permitted
!		If MAXPT <1 it is set to a default value of MAXPTS (=1100000)
!	
!	Required routines : F
!
      SUBROUTINE NGAUSS(A,B,N,M,IND,F,RI,IER,NUM,MAXPT)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=21,MAXPTS=1100000)
      DIMENSION IP(NMAX),NPT(NMAX),H(NMAX),XA(NMAX),WT(NMAX)
      DIMENSION A(N),B(N),IND(N),M(N),W(31),X(31)

!	Weights and abscissas for Gauss-Legendre quadrature.
!	For N-point formula W(K)=W(N-K+1) and X(K)=-X(N-K+1)
!		For K=1,2,...,N/2. Hence only half points are tabulated.
!	For 2-point W(1); 4-point W(2), W(3); 8-point W(4),...,W(7);
!	16-point W(8),...,W(15); 32-point W(16),...,W(31) are the
!	weights corresponding to abscissas X(I).

      DATA W/1.0D0,
     1       0.34785484513745385737D0, 0.65214515486254614263D0,
     2       0.10122853629037625915D0, 0.22238103445337447054D0,
     3       0.31370664587788728734D0, 0.36268378337836198297D0,
     4       0.02715245941175409485D0, 0.06225352393864789286D0,
     5       0.09515851168249278481D0, 0.12462897125553387205D0,
     6       0.14959598881657673208D0, 0.16915651939500253819D0,
     7       0.18260341504492358887D0, 0.18945061045506849629D0,
     8       0.00701861000947009660D0, 0.01627439473090567061D0,
     9       0.02539206530926205945D0, 0.03427386291302143310D0,
     1       0.04283589802222668066D0, 0.05099805926237617620D0,
     2       0.05868409347853554714D0, 0.06582222277636184684D0,
     3       0.07234579410884850623D0, 0.07819389578707030647D0,
     4       0.08331192422694675522D0, 0.08765209300440381114D0,
     4       0.09117387869576388471D0, 0.09384439908080456564D0,
     5       0.09563872007927485942D0, 0.09654008851472780057D0/

      DATA X/0.57735026918962576451D0,
     1       0.86113631159405257522D0, 0.33998104358485626480D0,
     2       0.96028985649753623168D0, 0.79666647741362673959D0,
     3       0.52553240991632898582D0, 0.18343464249564980494D0,
     4       0.98940093499164993260D0, 0.94457502307323257608D0,
     5       0.86563120238783174388D0, 0.75540440835500303390D0,
     6       0.61787624440264374845D0, 0.45801677765722738634D0,
     7       0.28160355077925891323D0, 0.09501250983763744019D0,
     8       0.99726386184948156354D0, 0.98561151154526833540D0,
     9       0.96476225558750643077D0, 0.93490607593773968917D0,
     1       0.89632115576605212397D0, 0.84936761373256997013D0,
     2       0.79448379596794240696D0, 0.73218211874028968039D0,
     3       0.66304426693021520098D0, 0.58771575724076232904D0,
     4       0.50689990893222939002D0, 0.42135127613063534536D0,
     5       0.33186860228212764978D0, 0.23928736225213707454D0,
     6       0.14447196158279649349D0, 0.04830766568773831623D0/

      IER=305
      RI=0.0
      NUM=0
      IF(N.GE.NMAX.OR.N.LT.1) RETURN
      IF(MAXPT.LE.0) MAXPT=MAXPTS
      IER=307

!	calculate the number of function evaluations required
      NUM=M(1)*IND(1)
      DO 200 I=2,N
200   NUM=NUM*M(I)*IND(I)
      IF(NUM.GT.MAXPT) RETURN

!	Initialisation
      IER=0
      DO 1000 I=1,N
        IP(I)=0
        NPT(I)=M(I)*IND(I)-1
        IF(M(I).NE.2.AND.M(I).NE.4.AND.M(I).NE.8.AND.M(I).NE.16.
     1      AND.M(I).NE.32) IER=306
        IF(IND(I).LT.1) IER=306
1000  H(I)=(B(I)-A(I))/(2*IND(I))
      IF(IER.NE.0) RETURN
      DO 1200 I=N+1,NMAX
        H(I)=1.0
        WT(I)=1.0
1200  CONTINUE

!	Loop for sum over N dimensions
      K=N

3000  DO 3100 I=K,1,-1
        M2=M(I)/2
!	The abscissas are X(NO),...,X(NO+M2-1)
        NO=M2
        H1=H(I)
        J1=IP(I)/M(I)
        J2=IP(I)-J1*M(I)
!	Use the (J2+1)th point in (J1+1)th subinterval
        X1=A(I)+(2*J1+1)*H1
        IF(J2.LT.M2) THEN
!	For the first M2 abscissas
          XA(I)=X1+H1*X(NO+J2)
          WT(I)=W(NO+J2)*WT(I+1)
        ELSE IF(J2-M2.LT.M2) THEN
!	For the next M2 abscissas
          XA(I)=X1-H1*X(NO+J2-M2)
          WT(I)=W(NO+J2-M2)*WT(I+1)
        ELSE
!	For Gaussian formula with odd number of points use the abscissa at x=0
          XA(I)=X1
          WT(I)=W(NO+M2)*WT(I+1)
        ENDIF
3100  CONTINUE

!	Add the new point to running sum
      RI=RI+WT(1)*F(N,XA)
      K=1
3200  IF(IP(K).GE.NPT(K)) GO TO 3400
!	try next point along Kth dimension
      IP(K)=IP(K)+1
      GO TO 3000

!	If Kth dimension is exhausted go to next one
3400  IP(K)=0
      K=K+1
      IF(K.LE.N) GO TO 3200

!	If all points are exhausted compute the value of integral
      DO 4000 I=1,N
4000  RI=RI*H(I)
      END

!     --------------------------------------------------

!	Multiple integration over a hyper-rectangle in n-dimensions
!	using compound monomial rules
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions
!	M : (input/output) Integer specifying the formula to be used
!		M can be 1, 3 or 5, otherwise it will be set to a default value of 3
!		M=1 selects 1-point formula of degree 1
!		M=3 selects 2N-point formula of degree 3 due to Stroud
!		M=5 selects (2N*N+1)-point formula of degree 5
!	IND : (input/output) Integer array of length N specifying the number
!		of subintervals to be used along each dimension. IND(J)>0
!		otherwise it will be set to a default value of 1
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	RINT : (output) The calculated value of the integral
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=39 implies specified accuracy was not achieved in
!			which case DIF will contain the estimated accuracy
!		IER=308 implies that number of points exceeded MAXPT in
!			first attempt and no approximation of RINT is calculated
!		IER=309 implies N<1 or N>50 and no calculations are done
!	NUM : (output) Number of function evaluations used by subroutine
!	MAXPT : (input/output) Maximum number of function evaluations permitted
!		If MAXPT <1 it is set to a default value of MAXPTS (=1000000)
!
!	FUNCTION F(N,X) must be supplied by the user
!	
!	Required routines : STROUD, F
!
      SUBROUTINE STRINT(A,B,N,M,IND,F,RINT,REPS,AEPS,DIF,IER,NUM,MAXPT)
      IMPLICIT LOGICAL(Q)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL F
      PARAMETER(MAXPTS=1100000)
      DIMENSION A(N),B(N),IND(N)

!	set M and IND(I) to default values if they are unacceptable
      IF(M.NE.1.AND.M.NE.3.AND.M.NE.5) M=3
      DO 1000 I=1,N
        IF(IND(I).LT.1) IND(I)=1
1000  CONTINUE
      IF(MAXPT.LE.0) MAXPT=MAXPTS
      NUM=0

!	Evaluate the first approximation to the integral
      CALL STROUD(A,B,N,M,IND,F,RINT,IER,NO,MAXPT)
      NUM=NUM+NO
      IF(IER.GT.100) RETURN
      IER=39

!	Iteration to check and improve the accuracy of integral
      DO 3000 I=1,10
        QC=.TRUE.
        DIF=0.0

        DO 2000 J=1,N
          IF(NUM.GT.MAXPT) RETURN
          I1=IND(J)
!	Double the number of subintervals along the Jth axis
          IND(J)=2*IND(J)

          CALL STROUD(A,B,N,M,IND,F,RINT1,IER1,NO,MAXPT)
          NUM=NUM+NO
          DIF=DIF+ABS(RINT-RINT1)
          IF(IER1.GT.100) RETURN

          IF(ABS(RINT1-RINT).LT.MAX(AEPS,ABS(RINT1)*REPS)) THEN
!	If satisfactory accuracy is achieved then restore the old value of IND(J)
            IND(J)=I1
          ELSE
!	else retain the new value
            RINT=RINT1
            QC=.FALSE.
          ENDIF
2000    CONTINUE

        IF(QC) THEN
!	If satisfactory accuracy is achieved, then return
          IER=0
          RETURN
        ENDIF
        IF(NUM.GT.MAXPT) RETURN
3000  CONTINUE
      END

!     --------------------------------------------------

!	Multiple integration over a hyper-rectangle in n-dimensions
!	using compound monomial rules with given number of points
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions, N>0 and N.LE.NMAX (=50)
!	M : (input) Integer specifying the formula to be used
!		M can be 1,3 or 5, otherwise IER is set to 310 and no
!		calculations are done
!		M=1 selects 1-point formula of degree 1
!		M=3 selects 2N-point formula of degree 3 due to Stroud
!		M=5 selects (2N*N+1)-point formula of degree 5
!	IND : (input) Integer array of length N specifying the number
!		of subintervals to be used along each dimension. IND(J)>0
!		otherwise IER is set to 310 and no calculations are done
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	RI : (output) The calculated value of the integral
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=308 implies that number of points exceeded MAXPT and
!			no calculations are done
!		IER=309 implies N<1 or N>NMAX, in which case no calculations are done
!		IER=310 implies M is not 1,3 or 5 or IND(J)<1 for some J
!			in which case no calculations are done
!	NUM : (output) Number of function evaluations used by subroutine
!	MAXPT : (input/output) Maximum number of function evaluations permitted
!		If MAXPT <1 it is set to a default value of MAXPTS (=1000000)
!
!	FUNCTION F(N,X) must be provided by the user
!	
!	Required routines :  F
!
      SUBROUTINE STROUD(A,B,N,M,IND,F,RI,IER,NUM,MAXPT)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=50,MAXPTS=1100000,PI=3.14159265358979324D0)
      DIMENSION X(NMAX),IP(NMAX),X3(NMAX,NMAX),H(NMAX),XA(NMAX),WT(NMAX)
      DIMENSION A(N),B(N),IND(N)

      IER=309
      RI=0.0
      IF(N.GT.NMAX.OR.N.LT.1) RETURN

!	Calculate the number of function evaluations required
      MPT=0
      IF(M.EQ.1) MPT=1
      IF(M.EQ.3) MPT=2*N
      IF(M.EQ.5) MPT=2*N*N+1
      NUM=MPT
      DO 200 I=1,N
        IF(IND(I).LT.1) IER=310
200   NUM=IND(I)*NUM
      IF(IER.EQ.310) RETURN
      IER=310
      IF(MPT.LE.0) RETURN
      IER=308
      IF(MAXPT.LT.1) MAXPT=MAXPTS
      IF(NUM.GT.MAXPT) RETURN

      IER=0
!	Constants for the (2N*N+1)-point formula of degree 5
      XI=SQRT(0.6D0)
      A0=(25*N*N-115*N+162D0)/(162.D0)
      A1=(70-25*N)/162.D0
      A2=25./324.D0

!	Abscissas for the 2N-point formula of degree 3
      IF(M.EQ.3) THEN
        XI3=SQRT(2./3.D0)
        XI2=1./SQRT(3.D0)
        DO 800 I=1,N
          DO 600 J=1,N-1,2
            AN=J*I*PI/N
            X3(J,I)=XI3*COS(AN)
            X3(J+1,I)=XI3*SIN(AN)
600       CONTINUE
!	When N is odd
          IF((N/2)*2.NE.N) X3(N,I)=XI2*(-1)**I
800     CONTINUE
      ENDIF

      DO 1000 I=1,N
        IP(I)=1
        H(I)=(B(I)-A(I))/(2*IND(I))
!	For abscissas of (2N*N+1)-point formula of degree 5
        WT(I)=H(I)*XI
1000  CONTINUE

!	loop for the sum over all subintervals
      K=N
1300  DO 1400 IN=K,1,-1
        XA(IN)=A(IN)+(2*IP(IN)-1)*H(IN)
1400  CONTINUE

      IF(M.EQ.1) THEN
!	Generalised midpoint rule
        R=F(N,XA)

      ELSE IF(M.EQ.3) THEN
!	Stroud's 2N-point rule of degree 3
        R=0.0
        DO 1600 I=1,N
          DO 1500 J=1,N
1500      X(J)=XA(J)+X3(J,I)*H(J)
          R=R+F(N,X)
          DO 1550 J=1,N
1550      X(J)=XA(J)-X3(J,I)*H(J)
          R=R+F(N,X)
1600    CONTINUE
        R=R/(2*N)

      ELSE IF(M.EQ.5) THEN
!	(2N*N+1)-point rule of degree 5
        R=F(N,XA)*A0
        S1=0.0
        S2=0.0
        DO 1800 I=1,N
1800    X(I)=XA(I)
        DO 2200 I=1,N
          X(I)=X(I)+WT(I)
          S1=S1+F(N,X)
          X(I)=XA(I)-WT(I)
          S1=S1+F(N,X)
          X(I)=XA(I)
          DO 2000 J=I+1,N
            X(I)=XA(I)+WT(I)
            X(J)=XA(J)+WT(J)
            S2=S2+F(N,X)
            X(J)=XA(J)-WT(J)
            S2=S2+F(N,X)
            X(I)=XA(I)-WT(I)
            S2=S2+F(N,X)
            X(J)=XA(J)+WT(J)
            S2=S2+F(N,X)
            X(J)=XA(J)
            X(I)=XA(I)
2000      CONTINUE
2200    CONTINUE
        R=R+A1*S1+A2*S2
      ENDIF

      RI=RI+R
      K=1
3200  IF(IP(K).GE.IND(K)) GO TO 3400
!	Go to the next subinterval along Kth dimension
      IP(K)=IP(K)+1
      GO TO 1300

!	If Kth dimension is exhausted, go to the next one
3400  IP(K)=1
      K=K+1
      IF(K.LE.N) GO TO 3200

!	If all directions are exhausted, compute the value of integral
      DO 4000 I=1,N
4000  RI=2.*RI*H(I)
      END

!     ------------------------------------------

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

!     ------------------------------------------
 
!	Function routine to transform from hyper-spherical coordinates to
!	Cartesian coordinates.
!	It can be used with any of the subroutines for multiple integration
!	for integration over hyper-spherical shells
!
!	N : (input) Number of dimensions
!	X : (input) Real array of length N, giving hyper-spherical coordinates 
!		First coordinate is radial distance r
!	Y : (output) Real array of length N giving the transformed Cartesian
!		coordinates, which is passed on to the FUNCTION FUNSPH(N,Y)
!		for calculating the required function.
!	Since there is no provision to pass error parameter from this
!	routine, if N exceeds NMAX then the execution will terminate.
!
!	The FUNCTION FUNSPH(N,Y) must be supplied by the user. The name of
!	this routine is not passed as argument and hence it has to have
!	the same name as occurring here.
!		
!	Required routines : FUNSPH

      FUNCTION SPHND(N,X)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=50)
      DIMENSION X(N),Y(NMAX)
 
      IF(N.GT.NMAX.OR.N.LT.1) STOP 301

!	The Cartesian coordinates are given by
!	y(1)=x(1)*COS(x(2))
!	y(2)=x(1)*SIN(x(2))*COS(x(3))
!	y(3)=x(1)*SIN(x(2))*SIN(x(3))*COS(x(4))
!	.........
!	y(n-1)=x(1)*SIN(x(2))*SIN(x(3))*...*SIN(x(n-1))*COS(x(n))
!	y(n)=x(1)*SIN(x(2))*SIN(x(3))*...*SIN(x(n))
!
!	P1 is the volume element in hyper-spherical coordinates

      T1=X(1)
      P1=1.0
      DO 1000 I=1,N-1
        Y(I)=T1*COS(X(I+1))
        P1=P1*T1
        T1=T1*SIN(X(I+1))
1000  CONTINUE
      Y(N)=T1

!	FUNSPH should calculate the required integrand using Cartesian coordinates
      SPHND=P1*FUNSPH(N,Y)
 
      END

!     ------------------------------------------
 
      FUNCTION FUNSPH(N,X)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(20)

!	THE INTEGRAND IN CARTESIAN COORDINATES

      XA=0.0
      DO 2500 I=1,N
2500  XA=XA+X(I)**2
      FUNSPH=XA
      END
