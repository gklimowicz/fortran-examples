!	To minimise a function of several variables using quasi-Newton method
!
!	N : (input) Number of variables
!	X : (input/output) Real array of length N containing the initial
!		guess for the minimum.
!		After execution it should contain the coordinates of minimiser
!	F : (output) The function value at X
!	G : (output) Real array of length N containing the gradient vector at X
!	H : (output) Real array of length N*N containing the estimated
!		Hessian matrix at X. The first dimension of H should be N.
!	NUM : (output) Number of function evaluations used by the subroutine
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy, all components of the
!		Minimiser will be calculated with accuracy MAX(AEPS, REPS*ABS(X(I)))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=53 implies that Hessian matrix at final point is probably singular
!			The iteration may have converged to a saddle point
!		IER=503 implies that N < 1, in which case no calculations are done
!		IER=526 implies that iteration failed to converge to specified accuracy
!		Other values may be set by LINMIN
!	FCN : (input) Name of the subroutine to calculate the function value
!		and its derivatives
!	WK : Real array of length 3N used as scratch space
!
!	SUBROUTINE FCN(N,X,F,G) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, F is the
!		function value at X and G is the gradient vector. X and G
!		are real arrays of length N. F and G must be calculated by FCN.
!
!	Required routines : LINMIN, FLNM, FCN
!
      SUBROUTINE BFGS(N,X,F,G,H,NUM,REPS,AEPS,IER,FCN,WK)
      IMPLICIT LOGICAL(Q)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL FCN
      PARAMETER(NIT=200)
      DIMENSION X(N),WK(3*N),H(N,N),G(N)

      IER=0
      IF(N.LT.1) THEN
        IER=503
        RETURN
      ENDIF

      DO 2000 I=1,N
!	Initialise the Hessian matrix to unit matrix
        DO 1800 J=1,N
1800    H(J,I)=0.0
        H(I,I)=1.0
2000  CONTINUE
!	If some variable needs to be kept fixed the corresponding
!	diagonal elements, H(I,I) should be set to zero.

      CALL FCN(N,X,F,G)
      DF=ABS(F)
      IF(DF.EQ.0.0) DF=1
      N2=2*N
      NUM=1
      H2=1.

!	The iteration loop
      DO 4000 IT=1,NIT
        DF1=0.0
        H1=H2
        H2=0.0
!	Calculating the search direction WK =S^(k)
        DO 2400 I=1,N
          WK(I)=0.0
          DO 2200 J=1,N
            H2=H2+ABS(H(I,J))
2200      WK(I)=WK(I)-H(I,J)*G(J)
          DF1=DF1+WK(I)*G(I)
2400    CONTINUE

        IF(DF1.EQ.0.0) THEN
!	If gradient vanishes, then quit
!	If Hessian matrix appears to be singular, set the error flag
          IF(ABS(H2/H1).GT.1.3D0) IER=53
          RETURN
        ENDIF

!	Initial guess for line search
        X1=0
        X2=MIN(1.0,-2.*DF/DF1)
        F1=F
        IF(X2.LE.0.0) X2=1
        CALL LINMIN(X1,X2,F1,DF1,REPS,AEPS,IER1,FCN,WK,X,N,NUM)
        IF(IER1.GT.0) IER=IER1
!	If line search fails, then quit
        IF(IER1.GT.100) RETURN

!	The convergence test
        QC=.TRUE.
        DO 2800 I=1,N
          X(I)=WK(N+I)
          WK(N+I)=X1*WK(I)
          IF(ABS(WK(N+I)).GT.MAX(REPS*ABS(X(I)),AEPS)) QC=.FALSE.
          GI=WK(N2+I)
          WK(N2+I)=WK(N2+I)-G(I)
          G(I)=GI
2800    CONTINUE
!	It is possible to apply convergence check on Function value using DF
!	instead of X(I)
        DF=F-F1
        F=F1
        IF(QC) THEN
          IF(ABS(H2/H1).GT.1.3D0) IER=53
          RETURN
        ENDIF

!	Update the matrix using BFGS formula
        DO 3200 I=1,N
          WK(I)=0.0
          DO 3200 J=1,N
3200    WK(I)=WK(I)+H(I,J)*WK(N2+J)
        GHG=0.0
        DG=0.0
        DO 3400 I=1,N
          DG=DG+WK(N2+I)*WK(N+I)
3400    GHG=GHG+WK(I)*WK(N2+I)
        R1=(1.+GHG/DG)/DG
        DO 3600 J=1,N
          DO 3600 I=1,N
3600   H(I,J)=H(I,J)+R1*WK(N+I)*WK(N+J)-(WK(N+I)*WK(J)+WK(I)*WK(N+J))/DG

4000  CONTINUE

!	Iteration fails to converge
      IER=526
      END
