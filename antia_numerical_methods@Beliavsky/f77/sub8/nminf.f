!	To minimise a function of several variables using direction set method
!
!	N : (input) Number of variables
!	X : (input/output) Real array of length N containing the initial
!		guess for the minimum.
!		After execution it should contain the coordinates of minimiser
!	F : (output) The function value at X
!	NUM : (output) Number of function evaluations used by the subroutine
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy, iteration will stop when
!		change in function value is less than MAX(AEPS, REPS*ABS(F))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=504 implies that N < 2, in which case no calculations are done
!		IER=529 implies that iteration failed to converge to specified accuracy
!		Other values may be set by LINMNF
!	FCN : (input) Name of the subroutine to calculate the function value
!	WK : Real array of length N*(2N+2) used as scratch space
!
!	SUBROUTINE FCN(N,X,F) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, F is the
!		function value at X. X is a real array of length N.
!
!	Required routines : LINMNF, FLN, SVD, FCN
!
      SUBROUTINE NMINF(N,X,F,NUM,REPS,AEPS,IER,FCN,WK)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL FCN
      PARAMETER(NIT=200)
      DIMENSION X(N),WK(N,2*(N+1))

      IER=0
      IF(N.LE.1) THEN
        IER=504
        RETURN
      ENDIF

!	Initialise the direction set matrix to identity matrix
      DO 2000 I=1,N
        WK(I,N+3)=0.0
        WK(I,N+4)=0.0
        DO 1800 J=1,N
1800    WK(J,I)=0.0
        WK(I,I)=1.0
2000  CONTINUE

      CALL FCN(N,X,F)
      NUM=1

!	The main iteration loop
      DO 6000 IT=1,NIT
        FI=F
        DO 4000 K=1,N

!	The starting point for line search
          DO 2400 I=1,N
2400      WK(I,N+1)=X(I)
          F0=F
          DFMAX=0.0
          KMAX=1

          DO 2800 I=1,N
            X1=0.0
!	Use previous value as initial approximation to minimum
            X2=WK(I,N+3)
            FUN=F0
            IF(X2.EQ.0.0) X2=1
            CALL LINMNF(X1,X2,F0,REPS,AEPS,IER1,FCN,WK(1,I),WK(1,N+1),
     1                  N,NUM)
            IF(IER1.GT.0) IER=IER1

            WK(I,N+3)=X1
!	Estimate of second derivative along the line
            IF(X1.NE.0.0) WK(I,N+4)=ABS((FUN-F0)/X1**2)
!	The new starting point
            DO 2600 J=1,N
2600        WK(J,N+1)=WK(J,N+1)+X1*WK(J,I)
            IF(FUN-F0.GE.DFMAX.AND.I.LE.N-K+1) THEN
              DFMAX=FUN-F0
              KMAX=I
            ENDIF
2800      CONTINUE

!	Remove the KMAX th direction
          DO 3000 I=KMAX,N-1
            WK(I,N+3)=WK(I+1,N+3)
            WK(I,N+4)=WK(I+1,N+4)
            DO 3000 J=1,N
3000      WK(J,I)=WK(J,I+1)
!	Add a new direction
          DO 3200 I=1,N
3200      WK(I,N)=WK(I,N+1)-X(I)

          X1=0.0
          X2=1
          FUN=F
          WK(N,N+4)=0.0
!	Starting point for the final line search in the loop
          DO 3400 I=1,N
3400      WK(I,N+1)=X(I)
          CALL LINMNF(X1,X2,F,REPS,AEPS,IER1,FCN,WK(1,N),WK(1,N+1),
     1                N,NUM)
          IF(IER1.GT.0) IER=IER1

          WK(N,N+3)=X1
          IF(X1.NE.0.0) WK(N,N+4)=ABS((FUN-F)/X1**2)
          DO 3600 J=1,N
3600      X(J)=X(J)+X1*WK(J,N)
4000    CONTINUE

        IF(ABS(F-FI).LT.MAX(REPS*ABS(F),AEPS)) RETURN

!	The matrix V for SVD
        DO 4500 J=1,N
          IF(WK(J,N+4).GT.0.0) THEN
            DO 4400 I=1,N
4400        WK(I,J)=WK(I,J)/SQRT(WK(J,N+4))
          ENDIF
4500    CONTINUE
        M=N
        LA=N
        LV=N
        CALL SVD(N,M,WK,WK(1,N+3),WK(1,N+1),LA,LV,WK(1,N+2),IER1)
        IF(IER1.GT.0) THEN
          IER=IER1
          RETURN
        ENDIF

        DO 4600 I=1,N
          WK(I,N+3)=SQRT(FUN-F)*WK(I,N+1)
          WK(I,N+4)=0.0
          IF(WK(I,N+1).NE.0.0) WK(I,N+4)=1./WK(I,N+1)**2
4600    CONTINUE

6000  CONTINUE
      IER=529
      END
