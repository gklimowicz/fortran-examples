!     PROGRAM TO FIND EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC  MATRIX
!     USING QL ALGORITHM OR STURM'S SEQUENCE
!     NOTE THAT THE NORMALISATION FOR EIGENVECTOR USED BY TQL2 AND TINVIT
!     ARE DIFFERENT

      PROGRAM REALSY
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(30,30),D(30),E(30),EV(30,30),EI(30),WK(210)
      DATA ((A(I,J),J=1,4),I=1,4)/6,4,4,1, 4,6,1,4, 4,1,6,4, 1,4,4,6/

51    FORMAT(10X,13HTHE MATRIX IS)
52    FORMAT(10F8.3)
53    FORMAT('   IER =',I4,5X,'DIAGONAL ELEMENTS =',1P3D14.6/
     1       (2X,5D14.6))
54    FORMAT('   OFF DIAGONAL ELEMENTS =',1P3D14.6/(2X,5D14.6))
55    FORMAT('   TRANSFORMATION MATRIX FOR REDUCTION TO TRIDIAGONAL'
     1       ,' FORM')
56    FORMAT(1P6D13.5)
57    FORMAT('    IER =',I4/'   EIGENVALUE',9X,'EIGENVECTOR')
58    FORMAT(1PD14.6,2X,4D14.6/(5X,5D14.6))
59    FORMAT('   IER =',I4,5X,'M1 =',I3,5X,'M2 =',I3/
     1       '  EIGENVALUES =',1P4D14.6/(2X,5D14.6))
60    FORMAT('   EIGENVALUE',9X,'EIGENVECTOR')

100   PRINT *,'TYPE N=ORDER, IFLG=0/1  FOR QL ALGORITHM/STURM SEQ'
      PRINT *,'           (QUITS WHEN N.LE.0)'
      READ *,N,IFLG
      IF(N.LE.0) STOP

      WRITE(6,51)
      PRINT *,'TYPE IN THE MATRIX ROW-WISE'
        DO 500 I=1,N
        PRINT *,I,'TH ROW'
        READ *,(A(I,J),J=1,N)
        WRITE(6,52) (A(I,J),J=1,N)
500   CONTINUE

!     TOLERANCE FOR TRED2
      IA=30

!     REDUCE THE MATRIX TO TRIDIAGONAL FORM USING HOUSEHOLDER'S METHOD

      CALL TRED2(A,N,IA,D,E,IER)
      WRITE(6,53) IER,(D(I),I=1,N)
      WRITE(6,54) (E(I),I=1,N)
      WRITE(6,55)
      DO 1000 I=1,N
1000  WRITE(6,56) (A(I,J),J=1,N)

      REPS=1.D-14
      IF(IFLG.EQ.0) THEN
!     QL-ALGORITHM FOR TRIDIAGONAL MATRIX

        CALL TQL2(A,N,IA,D,E,REPS,IER)
        WRITE(6,57) IER
        DO 2000 I=1,N
          WRITE(6,58) D(I),(A(J,I),J=1,N)
2000    CONTINUE
      ELSE

!     USE STURM SEQUENCE TO FIND THE REQUIRED EIGENVALUES

        PRINT *,'TYPE  M1,M2= S. NO. OF EIGENVALUES TO DETERMINED'
        READ *,M1,M2
        EPS1=0.01D0
        IV=30
        CALL TRIDIA(E,D,N,M1,M2,EI,EPS1,REPS,EV,IV,WK,IER)
        MZ=M2-M1+1
        WRITE(6,59)IER,M1,M2,(EI(I),I=1,MZ)

!       USE THE FOLLOWING STATEMENTS FOR PRINTING OUT EIGENVECTORS
!       OF TRIDIAGONAL MATRIX
!        DO 2500 I=1,MZ
!          WRITE(6,58) EI(I),(EV(J,I),J=1,N)
!2500    CONTINUE

!       BACK-TRANSFORMING THE MZ=M2-M1+1 EIGENVECTORS

        CALL TRBAK(A,IA,N,EV,IV,MZ)
        WRITE(6,60) 
        DO 3000 I=1,MZ
          WRITE(6,58) EI(I),(EV(J,I),J=1,N)
3000    CONTINUE
      ENDIF
      GO TO 100
      END

!     -----------------------------------------------------

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
      IMPLICIT REAL*8(A-H,O-Z)
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

!     --------------------------------------------------------------

!	To find a specified eigenvalue and eigenvector of a real symmetric
!	tridiagonal matrix using inverse iteration
!
!	E : (input) Real array of length N containing the off-diagonal elements
!		of the tridiagonal matrix, E(i+1)=A(i,i+1)=A(i+1,i)
!	D : (input) Real array of length N containing the diagonal elements
!		of the tridiagonal matrix, D(i)=A(i,i)
!	N : (input) Order of the matrix
!	EL : (input) lower limit on eigenvalue
!	EU : (input) upper limit on eigenvalue
!	EI : (output) The calculated eigenvalue
!	EV : (input/output) Real array of length N containing the calculated
!		eigenvector. If IFLG.NE.0 then EV should contain the previous
!		eigenvector, so that initial guess to EV is orthogonal to it
!	REPS : (input) Relative accuracy to which eigenvalues or eigenvectors
!		are calculated
!	IFLG : (input) Integer variable used as a flag to denote how initial
!		approximation to eigenvector is obtained. If IFLG=0 then
!		initial guess to eigenvector is generated using random numbers.
!		Otherwise, initial guess is chosen to be orthogonal to EV
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=144 implies that inverse iteration failed to converge
!	U : Real array of length N used as scratch space
!	B : Real array of length 4N used as scratch space
!	NUM : (output) Number of iterations required by the subroutine
!
!	FUNCTION RAN1(SEED) is required to generate random numbers.
!	If a different routine is used, the seed should be changed
!	appropriately.
!
!	Required routines : RAN1
!
      SUBROUTINE TINVIT(E,D,N,EL,EU,EI,EV,REPS,IFLG,IER,U,B,NUM)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=100)
      EXTERNAL RAN1
      DIMENSION E(N),D(N),EV(N),B(4,N),U(N)
!	The seed for random number generator
      DATA SEED/12345/

      E(1)=0.0
      IF(IFLG.EQ.0) THEN
!	Choose a random initial guess for eigenvector
        DO 1000 I=1,N
          U(I)=RAN1(SEED)
          EV(I)=U(I)
1000    CONTINUE
      ELSE

!	Choose the initial vector to be orthogonal to EV
        S=0.0
        EVS=0.0
        DO 1200 I=1,N
          U(I)=RAN1(SEED)
          EVS=EVS+EV(I)**2
1200    S=S+U(I)*EV(I)
        S=S/EVS
        DO 1400 I=1,N
          EV(I)=U(I)-EV(I)*S
1400    U(I)=EV(I)
      ENDIF

      EI=(EL+EU)/2.
!	Stretch the bounds to allow for roundoff errors
      EL1=EL-REPS*ABS(EL)
      EU1=EU+REPS*ABS(EU)
      IER=0

!	Loop for inverse iteration
      DO 5000 J=1,NIT
        NUM=J

!	Set up the matrix A-pI
        DO 2000 I=1,N
          B(1,I)=E(I)
          B(2,I)=D(I)-EI
          IF(I.LT.N) B(3,I)=E(I+1)
          B(4,I)=0.0
2000    CONTINUE
        B(3,N)=0.0

!	Gaussian elimination with partial pivoting
        DO 3000 K=1,N-1
          IF(ABS(B(2,K)).LT.ABS(B(1,K+1))) THEN
            DO 2200 I=1,3
              T=B(I,K+1)
              B(I,K+1)=B(I+1,K)
              B(I+1,K)=T
2200        CONTINUE
            T=EV(K+1)
            EV(K+1)=EV(K)
            EV(K)=T
          ENDIF

!	If a pivot is zero make it nonzero
          IF(B(2,K).EQ.0.0) B(2,K)=REPS
          R=-B(1,K+1)/B(2,K)
          B(2,K+1)=B(2,K+1)+R*B(3,K)
          B(3,K+1)=B(3,K+1)+R*B(4,K)
          EV(K+1)=EV(K+1)+R*EV(K)
3000    CONTINUE
        IF(B(2,N).EQ.0.0) B(2,N)=REPS

!	Back-substitution
        EV(N)=EV(N)/B(2,N)
        EV(N-1)=(EV(N-1)-EV(N)*B(3,N-1))/B(2,N-1)
        DO 3400 K=N-2,1,-1
          EV(K)=(EV(K)-EV(K+1)*B(3,K)-EV(K+2)*B(4,K))/B(2,K)
3400    CONTINUE

!	Calculating the Rayleigh quotient (DEI)
        R1=0.0
        KM=0
        S1=0.0
        S2=0.0
        DO 3500 K=1,N
          S1=S1+EV(K)*EV(K)
          S2=S2+EV(K)*U(K)
          IF(R1.LT.ABS(EV(K))) THEN
            R1=ABS(EV(K))
            KM=K
          ENDIF
3500    CONTINUE

        DEI=S2/S1
        UKM=EV(KM)
        DVS=0.0
!	Normalising the eigenvector
        DO 4000 K=1,N
          EV(K)=EV(K)/UKM
          DVS=DVS+ABS(EV(K)-U(K))
4000    CONTINUE

!	The convergence check
        REP=REPS*ABS(EI)
        IF(ABS(DEI).LT.REP.OR.R1*REP.GT.1.0.OR.DVS.LT.REPS) RETURN
        DO 4500 I=1,N
          U(I)=EV(I)
4500    CONTINUE
        IF(EI+DEI.GT.EL1.AND.EI+DEI.LT.EU1) THEN
!	The new shift
          EI=EI+DEI
        ELSE
!	If the new shift is outside the specified bounds, then modify it
          HI=MIN(ABS(EI-EL),ABS(EI-EU))
          EI=EI+SIGN(HI*0.5,DEI)
        ENDIF
5000  CONTINUE

!	Iteration fails to converge
      IER=144
      END

!     ------------------------------------------------------------

!	To find eigenvalues and eigenvectors of Z T Z^T using QL algorithm
!	where T is a symmetric tridiagonal matrix and Z is an orthogonal matrix.
!	If Z is the transformation matrix to reduce original matrix to
!	tridiagonal matrix, it will calculate the eigenvectors of original matrix
!
!	Z : (input/output) Real array of length IZ*N which should contain
!		the transformation matrix required to reduce original real
!		symmetric matrix to tridiagonal form. To find eigenvectors
!		of a symmetric tridiagonal matrix, set Z to unit matrix
!		After execution Z will contain the eigenvector of the original
!		matrix Z T Z^T. Z(i,j) should contain the ith component of
!		jth eigenvector
!	N : (input) Order of matrix
!	IZ : (input) The first dimension of array Z as declared in the
!		calling program. (IZ.GE.N)
!	D : (input/output) Real array of length N, containing the diagonal
!		elements of the tridiagonal matrix, D(i)=T(i,i).
!		After execution it will contain the eigenvalues of the matrix
!	E : (input/output) Real array of length N containing the off-diagonal
!		elements of the tridiagonal matrix, E(i)=T(i,i+1)=T(i+1,i)
!		It is used as scratch space and its contents will be destroyed
!		during execution.
!	REPS : (input) Required tolerance, it should be of order of machine
!		accuracy
!	IER : (output) Error parameter; IER=0 implies successful execution
!		IER=108 implies that N.LE.1 or N>IZ, in which case no
!			calculations are performed
!		IER=143 implies that the QL algorithm failed to converge
!			for some eigenvalue, the calculations are abandoned
!
!	Required routines : None
!	
      SUBROUTINE TQL2(Z,N,IZ,D,E,REPS,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=30)
      DIMENSION Z(IZ,N),D(N),E(N)

      IF(N.LE.0.OR.N.GT.IZ) THEN
        IER=108
        RETURN
      ENDIF
      IER=0
      DO 2000 I=2,N
2000  E(I-1)=E(I)
      E(N)=0
      B=0
      F=0

      DO 4000 L=1,N
        H=REPS*(ABS(D(L))+ABS(E(L)))
        IF(B.LT.H) B=H
!	Look for small off-diagonal elements
        DO 2200 M=L,N
          IF(ABS(E(M)).LE.B) GO TO 2300
2200    CONTINUE
        M=N
2300    IF(M.EQ.L) THEN
!	one eigenvalue is isolated
          D(L)=D(L)+F
          GO TO 4000
        ENDIF

!	Loop for QL transformation 
        DO 3600 IT=1,NIT
!	Find shift
          G=D(L)
          P=(D(L+1)-G)/(2*E(L))
          R=SQRT(P*P+1)
          IF(P.LT.0.0) R=-R
          D(L)=E(L)/(P+R)
          H=G-D(L)
          DO 2600 I=L+1,N
2600      D(I)=D(I)-H
          F=F+H

!	The QL transformation
          P=D(M)
          C=1
          S=0
!	Given's rotations
          DO 3400 I=M-1,L,-1
            G=C*E(I)
            H=C*P
            IF(ABS(P).GE.ABS(E(I))) THEN
              C=E(I)/P
              R=SQRT(C*C+1)
              E(I+1)=S*P*R
              S=C/R
              C=1/R
            ELSE
              C=P/E(I)
              R=SQRT(C*C+1)
              E(I+1)=S*E(I)*R
              S=1./R
              C=C/R
            ENDIF
            P=C*D(I)-S*G
            D(I+1)=H+S*(C*G+S*D(I))

!	Transforming the eigenvectors
            DO 3000 K=1,N
              H=Z(K,I+1)
              Z(K,I+1)=S*Z(K,I)+C*H
              Z(K,I)=C*Z(K,I)-S*H
3000        CONTINUE
3400      CONTINUE

          E(L)=S*P
          D(L)=C*P
          IF(ABS(E(L)).LE.B) THEN
!	One eigenvalue is isolated
            D(L)=D(L)+F
            GO TO 4000
          ENDIF
3600    CONTINUE
!	QL iteration fails to converge
        IER=143
        RETURN
4000  CONTINUE

!	Sort eigenvalues in ascending order by straight selection
      DO 5000 I=1,N-1
        K=I
        P=D(I)
        DO 4200 J=I+1,N
          IF(D(J).LT.P) THEN
            K=J
            P=D(J)
          ENDIF
4200    CONTINUE
        IF(K.NE.I) THEN
!	exchange the eigenvalues and eigenvectors
          D(K)=D(I)
          D(I)=P
          DO 4400 J=1,N
            P=Z(J,I)
            Z(J,I)=Z(J,K)
            Z(J,K)=P
4400      CONTINUE
        ENDIF
5000  CONTINUE
      END

!     --------------------------------------------------------

!	To perform back-transformation on eigenvectors of reduced tridiagonal matrix
!	to obtain eigenvectors of original real symmetric matrix reduced by TRED2
!	This routine is not required if eigenvectors are computed using TQL2
!
!	A : (input/output) Real array of length IA*N containing the transformation
!		matrix. The last column of this array is used as scratch
!		space and will be destroyed during execution
!	IA : (input) The first dimension of array A as declared in the
!		calling program (IA.GE.N)
!	N : (input) Order of the matrix
!	Z : (input/output) Real array of length IZ*NZ containing the
!		eigenvectors of tridiagonal matrix. After execution it will
!		be overwritten by eigenvectors of original matrix
!	IZ : (input) First dimension of array Z as declared in the calling
!		program (IZ.GE.N)
!	NZ : (input) Number of eigenvectors to be transformed
!
!	Required routines : None
!
      SUBROUTINE TRBAK(A,IA,N,Z,IZ,NZ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(IA,N),Z(IZ,NZ)

!	Loop on eigenvectors
      DO 4000 I=1,NZ
!	The matrix multiplication
        DO 3000 J=1,N
          S=0.0
          DO 2200 K=1,N-1
2200      S=S+A(J,K)*Z(K,I)
!	To take care of the last column of A, which is overwritten
          IF(J.EQ.N) S=S+Z(N,I)
          A(J,N)=S
3000    CONTINUE

        DO 3200 J=1,N
3200    Z(J,I)=A(J,N)
4000  CONTINUE
      END

!   ---------------------------------------------------------

!	Reduction of real symmetric matrix to tridiagonal form using
!	Householder's method
!
!	A : (input/output) Real array of length IA*N containing the matrix
!		elements. After execution it will be overwritten by the
!		transformation matrix
!	N : (input) Order of matrix
!	IA : (input) The first dimension of A as specified in the calling program
!	D : (output) Diagonal elements of the transformed tridiagonal matrix
!		D(I) would contain A(I,I)
!	E : (output) Off-diagonal elements of the transformed tridiagonal matrix
!		E(I+1) would contain A(I,I+1)=A(I+1,I)
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=107 implies that N.LE.1 or N>IA, in which case no
!			calculations are done
!
!	Required routines : None

      SUBROUTINE TRED2(A,N,IA,D,E,IER)
      IMPLICIT REAL*8(A-H,O-Z)
!	For REAL*4 use REPS=1.E-30
!      PARAMETER(REPS=1.E-30)
      PARAMETER(REPS=1.D-300)
      DIMENSION A(IA,N),D(N),E(N)

      IF(N.LE.1.OR.N.GT.IA) THEN
        IER=107
        RETURN
      ENDIF
      IER=0

      DO 4000 I=N,2,-1
        F=A(I,I-1)
        G=0
        DO 2000 K=1,I-2
2000    G=G+A(I,K)*A(I,K)
        H=G+F*F
        IF(G.LE.REPS) THEN
!	Skip the transformation
          E(I)=F
          H=0
        ELSE

          G=SQRT(H)
          IF(F.GT.0.0) G=-G
          E(I)=G
          H=H-F*G
          A(I,I-1)=F-G
          F=0
          DO 2600 J=1,I-1
!	Elements of u_i/H_i
            A(J,I)=A(I,J)/H
            G=0
!	Form elements of A_iu_i
            DO 2200 K=1,J
2200        G=G+A(J,K)*A(I,K)
            DO 2400 K=J+1,I-1
2400        G=G+A(K,J)*A(I,K)
!	Components of p_i
            E(J)=G/H
            F=F+G*A(J,I)
2600      CONTINUE

!	calculate u_i^Tp_i/2H_i
          HH=0.5*F/H
          DO 3400 J=1,I-1
            F=A(I,J)
            G=E(J)-HH*F
!	Elements of q_i
            E(J)=G
            DO 3000 K=1,J
3000        A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
3400      CONTINUE
        ENDIF
        D(I)=H
4000  CONTINUE

      D(1)=0
      E(1)=0
!	accumulation of transformation matrix Q
      DO 5000 I=1,N
        IF(D(I).NE.0.0) THEN
          DO 4600 J=1,I-1
            G=0
            DO 4200 K=1,I-1
4200        G=G+A(I,K)*A(K,J)
            DO 4400 K=1,I-1
4400        A(K,J)=A(K,J)-G*A(K,I)
4600      CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1
        DO 4800 J=1,I-1
          A(I,J)=0.0
4800    A(J,I)=0.0
5000  CONTINUE
      END

!  -----------------------------------------------------------------

!	To find specified eigenvalues and eigenvectors of a real symmetric
!	tridiagonal matrix using Sturm sequence coupled with inverse iteration
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
!	EI : (output) Real array of length M2-M1+1 containing the calculated
!		eigenvalues
!	EPS1 : (input) Relative accuracy to which eigenvalues are located by bisection
!		before using inverse iteration. If the inverse iteration
!		does not converge to nearest eigenvalues, EPS1 can be reduced.
!		A value of 0.1-0.01 times typical eigenvalue is generally sufficient.
!	REPS : (input) Desired relative accuracy in eigenvalues and eigenvectors
!	EV : (output) Real array of length IV*(M2-M1+1) containing the
!		eigenvectors. EV(i,j) should contain the ith component of
!		the jth eigenvector
!	IV : (input) The first dimension of array EV as declared in the
!		calling program
!	WK : Real array of length 7N used as scratch space
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=109 implies that N.LE.1 or N>IV or M1<1 or M2>N
!		Other values may be set by TINVIT, only the last value is
!		returned.
!
!	Required routines : STURM, TINVIT, RAN1
!
!
      SUBROUTINE TRIDIA(E,D,N,M1,M2,EI,EPS1,REPS,EV,IV,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION E(N),D(N),WK(N,7),EI(M2-M1+1),EV(IV,M2-M1+1)

      IF(N.LE.1.OR.N.GT.IV.OR.M1.LT.1.OR.M2.GT.N) THEN
        IER=109
        RETURN
      ENDIF
      IER=0
      IF(M1.GT.M2) RETURN
!	Locate the eigenvalues
      CALL STURM(E,D,N,M1,M2,WK(1,2),WK(1,3),NUM,EPS1,WK,IER1)

!	Loop for finding individual eigenvalues and eigenvectors
      DO 2000 I=M1,M2
        IFLG=0
        I1=I-M1+1

        IF(I.GT.M1) THEN
          IF(ABS(WK(I,2)-WK(I-1,2)).LT.3.*REPS*ABS(WK(I,2))) THEN
!	Set the flag for close eigenvalues
            IFLG=1
            DO 1500 J=1,N
1500        EV(J,I1)=EV(J,I1-1)
          ENDIF
        ENDIF
        CALL TINVIT(E,D,N,WK(I,2),WK(I,3),EI(I1),EV(1,I1),REPS,
     1              IFLG,IER1,WK(1,1),WK(1,4),NUM)
        IF(IER1.GT.0) IER=IER1
2000  CONTINUE
      END

!     -------------------------------------------------------------

!	To generate uniformly distributed random numbers in interval (0,1)
!
!	SEED : (input/output) is a real value used as the seed
!		It should be positive during initial call and
!		should not be modified between different calls.
!
!	Required routines : None

      FUNCTION RAN1(SEED)
      IMPLICIT REAL*8(A-H,O-Z)
!	Retain the following declaration even for REAL*4 version
!	otherwise AM and AC will be rounded
      REAL*8 AM,A,AC,AN
      PARAMETER(AM=2147483648D0,A=45875,AC=453816693D0,AN=2147483647D0)

      SEED=MOD(SEED*A+AC,AM)
      RAN1=SEED/AN
      END
