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
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMAX=50,MAXPTS=1100000)
      PARAMETER(PI=3.14159265358979323846264338327950288Q0)
!      PARAMETER(NMAX=50,MAXPTS=1100000,PI=3.14159265358979324D0)
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
      XI=SQRT(0.6Q0)
      A0=(25*N*N-115*N+162Q0)/(162.Q0)
      A1=(70-25*N)/162.Q0
      A2=25./324.Q0

!	Abscissas for the 2N-point formula of degree 3
      IF(M.EQ.3) THEN
        XI3=SQRT(2./3.Q0)
        XI2=1./SQRT(3.Q0)
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
