!	To find specified eigenvalues of a real symmetric tridiagonal 
!	matrix using bisection on Sturm sequence
!
!	E : (input) Real array of length N containing the off-diagonal elements
!		of the tridiagonal matrix, E(i+1)=A(i,i+1)=A(i+1,i)
!	D : (input) Real array of length N containing the diagonal elements
!		of the tridiagonal matrix, D(i)=A(i,i)
!	N : (input) Order of the matrix
!	M1 : (input) Serial number of lowest eigenvalue to be determined.
!		The eigenvalues are sorted in increasing order
!	M2 : (input) Serial number of highest eigenvalue to be determined.
!		All eigenvalues from M1 to M2 are determined
!	EL : (output) Real array of length M2 containing the calculated
!		lower limit on eigenvalues
!	EU : (output) Real array of length M2 containing the calculated
!		upper limit on eigenvalues
!		The ith eigenvalue is located in interval (EL(i),EU(i))
!	NUM : (output) Number of times Sturm sequence was evaluated to locate
!		the eigenvalues.
!	REPS : (input) Relative accuracy to which eigenvalues are located by bisection
!	WK : Real array of length N used as scratch space
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=110, implies that M1<1 or M2>N, in which case no
!			calculations are done
!
!	Required routines : None

      SUBROUTINE STURM(E,D,N,M1,M2,EL,EU,NUM,REPS,WK,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(XUF=-1.D29,EPS=1.D-30,NBIS=70)
      DIMENSION E(N),D(N),EL(M2),EU(M2),WK(N)

      IF(M1.GT.M2) RETURN
      IF(M1.LT.1.OR.M2.GT.N) THEN
        IER=110
        RETURN
      ENDIF
      E(1)=0.0
!	Finding bounds on eigenvalues using Gerschgorin's theorem
      EMIN=D(N)-ABS(E(N))
      EMAX=D(N)+ABS(E(N))
      WK(1)=0.0
      DO 1000 I=N-1,1,-1
        R=ABS(E(I))+ABS(E(I+1))
        IF(D(I)+R.GT.EMAX) EMAX=D(I)+R
        IF(D(I)-R.LT.EMIN) EMIN=D(I)-R
        WK(I+1)=E(I+1)**2
1000  CONTINUE

!	Initialise the limits to undefined values
      DO 1500 I=M1,M2
        EL(I)=XUF
        EU(I)=XUF
1500  CONTINUE
      NL=0
      NU=N
      NUM=0

!	Loop for each eigenvalue
      DO 5000 M=M1,M2
        IF(EL(M).EQ.XUF) THEN
!	If the lower bound is undefined, then use EMIN
          EL(M)=EMIN
          N1=NL
        ENDIF

        IF(EU(M).EQ.XUF) THEN
!	If upper bound is undefined, use the bound for some higher eigenvalue
!	and if none exists, then use EMAX
          DO 2000 I=M+1,M2
            IF(EU(I).NE.XUF) THEN
              EU(M)=EU(I)
              N2=I
              GO TO 2200
            ENDIF
2000      CONTINUE
          EU(M)=EMAX
          N2=NU
        ENDIF

2200    IFL=0
!	Loop for bisection
        DO 4000 I=1,NBIS
          E1=(EL(M)+EU(M))/2.
          IF(E1.EQ.EL(M).OR.E1.EQ.EU(M)) GO TO 4200
          NUM=NUM+1

!	Count the number of sign changes in the Sturm sequence
          NS=0
          Q1=D(1)-E1
          IF(Q1.LT.0.0) NS=1
          DO 2400 K=2,N
            Q2=D(K)-E1
            IF(Q1.NE.0.0) Q2=Q2-WK(K)/Q1
            IF(Q1.EQ.0.0) Q2=Q2-ABS(E(K))/EPS
            IF(Q2.LT.0.0) NS=NS+1
2400      Q1=Q2

!	If the bounds are two consecutive real number on the machine, then quit
          IF(E1.EQ.EL(M).OR.E1.EQ.EU(M)) GO TO 4200
!	Update the bounds
          IF(NS.GE.M1.AND.NS.LE.M2) EU(NS)=MIN(E1,EU(NS))
          IF(NS.GE.M1-1.AND.NS.LT.M2) EL(NS+1)=MAX(E1,EL(NS+1))
          IF(NS.LT.M) THEN
            N1=MAX(N1,NS)
            IF(E1.GT.EL(M)) EL(M)=E1
          ELSE
            N2=MIN(N2,NS)
            IF(E1.LT.EU(M)) EU(M)=E1
          ENDIF

          IF(N1.EQ.M-1.AND.N2.EQ.M) THEN
!	The eigenvalue is isolated
            IF(IFL.GT.3.AND.ABS(EU(M)-EL(M)).LT.REPS) GO TO 4500
            IF(M.EQ.M1) THEN
              IFL=IFL+1
            ELSE
              IF(EL(M).NE.EU(M-1)) IFL=IFL+1
            ENDIF
          ENDIF
4000    CONTINUE

!	If the eigenvalue cannot be isolated, then set the same bounds
!	for all of them
4200    DO 4400 K=M+1,N2
          EL(K)=EL(M)
4400    EU(K)=EU(M)
4500    N1=M
        N2=MAX(N2,M+1)
5000  CONTINUE
      END
