!	To solve a generalised eigenvalue problem in ordinary differential
!	equations using finite difference method
!	Version for complex eigenvalue
!
!	N : (input) Number of mesh points to be used.
!	M : (input) Number of first order differential equations in the system
!	ML : (input) Number of boundary conditions at the first boundary t=T(1)
!	PAR : (input/output) Complex array to be passed on to EQN, EQND, BCS and BCSD
!		for specifying the equations and boundary conditions. PAR(1)
!		is used for passing the eigenvalue and hence should not be
!		used for any other variable. After execution the eigenvalue
!		will be available in PAR(1). Other elements of this array are
!		not used by GEVP_C, but are merely used to pass on any extra parameters
!		that may be required to specify the equations
!	X : (output) Complex array of length M*N containing the eigenvector.
!		X(i,j) is the ith component of solution at jth mesh point.
!		First dimension of X in calling program must be M
!	XC : (output) Complex array of length M*N containing the left eigenvector.
!		It is stored in the same format as X.
!	T : (input) Real array of length N containing the mesh points.
!		These points must be in ascending or descending order.
!		For calculating the deferred correction the mesh spacing
!		must be uniform.
!	E0 : (output) The calculated eigenvalue including deferred correction
!		PAR(1) contains the uncorrected eigenvalue. If deferred
!		correction is not applied E0 is not calculated.
!	EQN : (input) Name of subroutine to specify the differential equation
!		to be solved.
!	BCS : (input) Name of subroutine to calculate the boundary conditions
!		at t=T(1) and T(N)
!	EQND : (input) Name of subroutine to calculate derivative of equation
!		matrix with respect to the eigenvalue.
!	BCS : (input) Name of subroutine to calculate derivative of the boundary
!		conditions with respect to the eigenvalue
!	IWK : Integer array of length M*N used as scratch space
!	WK : Complex array of length (M+ML)*2M*(N+1) used as scratch space
!	IFLAG : (input) Integer variable used as a flag to decide the type of
!		computation required.
!		IFLAG=0 implies that only the eigenvalue is calculated.
!		IFLAG=1 implies that both eigenvalue and eigenvector are
!			calculated.
!		IFLAG=2 implies that first order correction to eigenvalue
!			is also calculated in addition to eigenvector
!			In this case the mesh spacing must be uniform
!	REPS : (input) Required accuracy. This is only used to check convergence
!		of Muller's method for finding zeros of the determinant
!		The truncation error depends on mesh spacing and is not
!		controlled in this routine.
!	RMX : (input) Expected maximum limit on the eigenvalue. This parameter
!		is passed on to MULER2 and iteration is terminated if the
!		magnitude of estimated eigenvalue exceeds RMX.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=704 implies that N<3, M.LE.ML, or ML.LE.0, in which case
!			no calculations are done
!		IER=734 implies that N<5 and deferred correction is requested.
!			In this case deferred correction is not calculated
!		IER=735 implies that the finite difference matrix is singular
!			while calculating the eigenvector. In this case
!			eigenvector is not calculated, but eigenvalue is
!			available in PAR(1).
!		IER=736 implies that mesh spacing is not uniform and
!			deferred correction is not calculated
!		IER=738 implies that eigenvector vanishes.
!		IER=739 implies that inverse iteration for calculating the
!			eigenvector failed to converge.
!		IER=740 implies that inverse iteration for calculating the
!			left eigenvector failed to converge.
!
!	SUBROUTINE EQN, EQND, BCS and BCSD must be supplied by the user.
!	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) calculates the right hand
!		sides for differential equations B y'=A y
!		J is the serial number of mesh point at which calculation
!		is required. M is the number of first order differential
!		equations in the system, ML the number of boundary conditions
!		at the first point, PAR is a Complex array which can be used
!		to pass on any required parameters to define the equations.
!		PAR(1) is the eigenvalue.
!		A and B are Complex arrays of length (M+ML)*M defining the
!		differential equation B Y'= A Y.
!		Y is Complex array of length M specifying the solution at t=T
!		F is a Complex array of length M which should be set to zero.
!		F, A and B must be calculated by the subroutine.
!
!	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) calculates the boundary
!		conditions at both boundaries. M is the number of differential
!		equations, ML is the number of boundary condition at the
!		first mesh point t=T1. PAR is a Complex array which can be used
!		to pass any required parameters. PAR(1) is the eigenvalue.
!		BC is a Complex array of length (M+ML)*M which should contain
!		the coefficients of boundary conditions, BC Y =0.
!		First ML rows will specify the boundary conditions
!		at t=T1, while remaining rows will specify those at t=TN.
!		G is a Complex array of length M which should be set to zero.
!		Y1 and YN are Complex arrays of length M specifying the
!		solution at t=T1 and TN.
!
!	SUBROUTINE EQND(J,M,ML,PAR,A,B,T) calculates the derivatives of
!		matrices A and B (calculated by EQN) with respect to the
!		eigenvalue (PAR(1)). 
!		J is the serial number of mesh point at which calculation
!		is required. M is the number of first order differential
!		equations in the system, ML the number of boundary conditions
!		at the first point, PAR is a Complex array which can be used
!		to pass on any required parameters to define the equations.
!		PAR(1) is the eigenvalue.
!		A and B are Complex arrays of length (M+ML)*M which should
!		give the derivatives of arrays A and B calculated by EQN.
!		T is the value of t at which the derivatives are required.
!
!	SUBROUTINE BCSD(M,ML,PAR,BC,T1,TN) calculates the derivative of
!		matrix BC (calculated by BCS) with respect to the eigenvalue PAR(1).
!		M is the number of differential equations,
!		ML is the number of boundary condition at the
!		first mesh point t=T1. PAR is a Complex array which can be used
!		to pass any required parameters. PAR(1) is the eigenvalue.
!		BC is a Complex array of length (M+ML)*M which should contain
!		the derivative of the matrix BC calculated by BCS.
!		First ML rows will specify the boundary conditions
!		at t=T1, while remaining rows will specify those at t=TN.
!
!	Required routines : SETMAT_C, GAUBLK_C, MULER2, EQN, EQND, BCS, BCSD
!
      SUBROUTINE GEVP_C(N,M,ML,PAR,X,XC,T,E0,EQN,BCS,EQND,BCSD,
     1                IWK,WK,IFLAG,REPS,RMX,IER)
!      IMPLICIT REAL*8(A,B,D-H,O-Z)
!      IMPLICIT COMPLEX*16(C)
!	For complex eigenvalues use the following statements
      IMPLICIT REAL*4(H,O,R,S,T)
      IMPLICIT COMPLEX*8(A-G,P,U-Z)
      EXTERNAL EQN,BCS
      PARAMETER(NIT=20,EPS=1.D-30)
      DIMENSION PAR(*),T(N),X(M,N),XC(M,N),IWK(M,N),WK(M+ML,2*M,N+1),
     1          CZERO(1)

      IF(N.LT.3.OR.M.LE.ML.OR.ML.LE.0) THEN
        IER=704
        RETURN
      ENDIF

      RAPS=0.1D0*REPS
      DE=MAX(1.E-3*ABS(E0),100.*REPS)
      CX1=E0+DE
      CX2=E0-DE
      CX3=E0
      NZ=0
      IER=0

!	First call to MULER2
1000  CALL MULER2(CX1,CX2,CX3,REPS,RAPS,IER,CF,CX,IDET,NZ,CZERO,RMX)
      IF(IER.LT.0) THEN
!	calculate the determinant
        PAR(1)=CX
        CALL SETMAT_C(N,M,ML,WK,WK(1,1,N+1),X,XC,T,PAR,EQN,BCS)
        IFLG=1
        CALL GAUBLK_C(N,M,ML,WK,IFLG,DET,IDET,IWK,XC,IER1)
        CF=DET
!	Call MULER2 again
        GO TO 1000
      ENDIF

!	the eigenvalue
      PAR(1)=CX3
      IF(IER.GT.100.OR.IFLAG.EQ.0) RETURN
      IF(IER1.GT.0) THEN
        IER=735
        RETURN
      ENDIF

!	Initialise the vector for inverse iteration
      DO 1200 I=1,N
        DO 1200 J=1,M
          X(J,I)=1.
1200  XC(J,I)=1.
      IFLG=2

!	Loop for inverse iteration to calculate the eigenvector
      DO 2500 I=1,NIT
        CALL GAUBLK_C(N,M,ML,WK,IFLG,DET,IDET,IWK,X,IER)

!	Normalising the eigenvector
        RMAX=0.0
        DO 1500 J=1,N
          DO 1500 K=1,M
            IF(ABS(X(K,J)).GT.RMAX) THEN
              RMAX=ABS(X(K,J))
              XMAX=X(K,J)
            ENDIF
1500    CONTINUE
        IF(RMAX.EQ.0.0) THEN
!	eigenvector vanishes, hence quit
          IER=738
          RETURN
        ENDIF

        DO 1800 J=1,N
          DO 1800 K=1,M
1800    X(K,J)=X(K,J)/XMAX

!	convergence check
        IF(RMAX*REPS.GT.1.0) GO TO 3000
2500  CONTINUE
      IER=739
      RETURN

3000  IF(IFLAG.LE.1) RETURN
      IF(N.LE.4) THEN
        IER=734
        RETURN
      ENDIF

!	Loop for inverse iteration for the left eigenvector
      DO 5000 IT=1,NIT

!	Forward substitution for U^T
        DO 3500 I=1,N
          DO 3500 J=1,M
            IF(I.GT.1) THEN
              DO 3200 K=1,M
3200          XC(J,I)=XC(J,I)-XC(K,I-1)*WK(K,J+M,I-1)
            ENDIF
              DO 3400 K=1,J-1
3400          XC(J,I)=XC(J,I)-XC(K,I)*WK(K,J,I)
3500    XC(J,I)=XC(J,I)/WK(J,J,I)

!	Back substitution for L^T
        MR=M
        DO 4000 I=N,1,-1
          IF(I.EQ.N-1) MR=M+ML
          DO 4000 J=MIN(M,MR-1),1,-1
            DO 3900 K=J+1,MR
              K0=K
              I1=I
!	Checking for backward interchange
              DO 3600 K1=J+1,K
                IF(K1.GT.M) K0=K-M
                IF(IWK(K1,I).EQ.K0) GO TO 3900
3600          CONTINUE

              K1=K
3700          IF(K1.GT.M) THEN
                K1=K1-M
                I1=I1+1
              ENDIF
              IF(IWK(K1,I1).GT.K1) THEN
                K0=IWK(K1,I1)
                K2=K0
                DO 3800 KK=K1+1,K2
                  IF(KK.GT.M) K0=K2-M
                  IF(IWK(KK,I1).EQ.K0) THEN
                    K1=KK
                    GO TO 3900
                  ENDIF
3800            CONTINUE
                K1=K2
                IF(I1.LE.N) GO TO 3700
!	Hopefully, this statement will not be executed
                STOP
              ENDIF
3900        XC(J,I)=XC(J,I)-XC(K1,I1)*WK(K,J,I)
4000    CONTINUE

!	Exchanges due to interchange matrix N
        DO 4200 I=N,1,-1
          DO 4200 J=M,1,-1
            JJ=IWK(J,I)
            IF(JJ.NE.J) THEN
              XT=XC(J,I)
              XC(J,I)=XC(JJ,I)
              XC(JJ,I)=XT
            ENDIF
4200    CONTINUE

        RMAX=0.0
        DO 4400 I=1,N
          DO 4400 J=1,M
            IF(RMAX.LT.ABS(XC(J,I))) RMAX=ABS(XC(J,I))
4400    CONTINUE
        DO 4600 I=1,N
          DO 4600 J=1,M
4600    XC(J,I)=XC(J,I)/RMAX

        IF(RMAX*REPS.GT.1.0) GO TO 6000
5000  CONTINUE
      IER=740
      RETURN

!	calculating the deferred correction
6000  HH=T(2)-T(1)
      ENUM=0.0
      EDEN=0.0
      DO 7000 J=1,N-1
        TJ=0.5*(T(J)+T(J+1))
        H=T(J+1)-T(J)
        IF(ABS(H-HH).GT.1.D-4*ABS(HH)) THEN
          IER=736
          RETURN
        ENDIF

        JJ=J
        CALL EQN(JJ,M,ML,PAR,WK(1,1,2),WK(1,M+1,2),X(1,J),WK(1,1,3),TJ)
        CALL EQND(JJ,M,ML,PAR,WK(1,1,3),WK(1,M+1,3),TJ)

        DO 6200 I=1,M
          IF(J.EQ.1) THEN
            WK(I,1,1)=(-2*X(I,1)+7*X(I,2)-9*X(I,3)+5*X(I,4)-X(I,5))/24.
            WK(I,2,1)=(43*X(I,1)-112*X(I,2)+102*X(I,3)-40*X(I,4)+
     1                 7*X(I,5))*H/192.
          ELSE IF(J.EQ.N-1) THEN
            WK(I,1,1)=(2*X(I,N)-7*X(I,N-1)+9*X(I,N-2)-5*X(I,N-3)+
     1                  X(I,N-4))/24.
            WK(I,2,1)=(43*X(I,N)-112*X(I,N-1)+102*X(I,N-2)-
     1                 40*X(I,N-3)+7*X(I,N-4))*H/192.
          ELSE
            WK(I,1,1)=(-X(I,J-1)+3.*X(I,J)-3.*X(I,J+1)+X(I,J+2))/24.
            WK(I,2,1)=(X(I,J-1)-X(I,J)-X(I,J+1)+X(I,J+2))*H/16.
          ENDIF
6200    CONTINUE

        DO 6500 I=1,M
          XN=0.0
          XD=0.0
          DO 6400 K=1,M
            XD=XD+WK(I,K+M,3)*(X(K,J+1)-X(K,J))-0.5*H*WK(I,K,3)*
     1         (X(K,J)+X(K,J+1))
            XN=XN+WK(I,K+M,2)*WK(K,1,1)-WK(I,K,2)*WK(K,2,1)
6400      CONTINUE
          ENUM=ENUM+XC(ML+I,J)*XN
          EDEN=EDEN+XC(ML+I,J)*XD
6500    CONTINUE
7000  CONTINUE

      CALL BCSD(M,ML,PAR,WK(1,1,1),T(1),T(N))
      DO 7200 I=1,ML
        XD=0.0
        DO 7100 K=1,M
7100    XD=XD+WK(I,K,1)*X(K,1)
7200  EDEN=EDEN+XD*XC(I,1)
      DO 7500 I=ML+1,M
        XD=0.0
        DO 7400 K=1,M
7400    XD=XD+WK(I,K,1)*X(K,N)
7500  EDEN=EDEN+XD*XC(I,N)

      DEF=ENUM/EDEN
      E0=PAR(1)+DEF
      END
!	--------------------------------
!	
!	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T)
!       IMPLICIT REAL*4(H,O,R,S,T)
!       IMPLICIT COMPLEX*8(A-G,P,U-Z)
!	DIMENSION A(M+ML,M),B(M+ML,M),Y(M),PAR(*),F(M)
!
!	DO 1000 I=1,M
!	F(I)=0.0
!	  DO 1000 K=1,M
!	    A(K,I)=a_{ki}(T,PAR)
!	    B(K,I)=b_{ki}(T,PAR)
! 1000  CONTINUE
!       END
!
!	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN)
!       IMPLICIT REAL*4(H,O,R,S,T)
!       IMPLICIT COMPLEX*8(A-G,P,U-Z)
!	DIMENSION PAR(*),BC(M+ML,M),G(M),Y1(M),YN(M)
!
!	DO 1000 I=1,M
!	    G(I)=0.0
!	  DO 1000 K=1,M
!	  IF(I.LE.ML) THEN
!	    BC(I,K)=bc_{ik}(T1,PAR)
!	  ELSE
!	    BC(I,K)=bc_{ik}(TN,PAR)
!         ENDIF
! 1000  CONTINUE
!       END
!	
!	SUBROUTINE EQND(J,M,ML,PAR,A,B,T)
!       IMPLICIT REAL*4(H,O,R,S,T)
!       IMPLICIT COMPLEX*8(A-G,P,U-Z)
!	DIMENSION A(M+ML,M),B(M+ML,M),PAR(*)
!
!	DO 1000 I=1,M
!	  DO 1000 K=1,M
!	    A(K,I)=da_{ki}/dPAR(1) (T,PAR)
!	    B(K,I)=db_{ki}/dPAR(1) (T,PAR)
! 1000  CONTINUE
!       END
!
!	SUBROUTINE BCSD(M,ML,PAR,BC,T1,TN)
!       IMPLICIT REAL*4(H,O,R,S,T)
!       IMPLICIT COMPLEX*8(A-G,P,U-Z)
!	DIMENSION PAR(*),BC(M+ML,M)
!
!	DO 1000 I=1,M
!	  DO 1000 K=1,M
!	  IF(I.LE.ML) THEN
!	    BC(I,K)=d bc_{ik}/dPAR(1) (T1,PAR)
!	  ELSE
!	    BC(I,K)=d bc_{ik}/dPAR(1) (TN,PAR)
!         ENDIF
! 1000  CONTINUE
!       END
