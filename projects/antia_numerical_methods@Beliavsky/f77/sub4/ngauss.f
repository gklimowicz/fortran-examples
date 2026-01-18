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
