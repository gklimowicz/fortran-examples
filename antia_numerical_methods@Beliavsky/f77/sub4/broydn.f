!	Real root of a system of nonlinear equations using Broyden's method
!
!	FCN : (input) Name of the subroutine to calculate the vector function
!	NP : (input) Number of variables and equations to be solved
!	X : (input/output) Real array of length NP containing the starting
!		values. After execution it should contain the computed roots
!	F : (output) Real array of length NP containing the function values
!		at the point specified by array X
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=441 implies that GAUELM failed to solve the required
!			system of linear equations
!		IER=442 implies that iteration did not converge to required
!			accuracy
!	WK : Real array of length MAX(2*NP*NP,NP*(NP+4)) used as scratch space
!	IWK : integer array of length NP used as scratch space
!
!	SUBROUTINE FCN(NP,X,F) must be supplied by the user.
!	Here X and F are real arrays of length NP containing the value of
!	vector function at point X.
!	The system of equations to be solved are F(I)=0
!
!	Required routines : GAUELM, FCN
!
      SUBROUTINE BROYDN(FCN,NP,X,F,REPS,AEPS,IER,WK,IWK)
      IMPLICIT LOGICAL(Q)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION X(NP),F(NP),WK(NP,*),IWK(NP)

      LJ=NP
      IER=0
      CALL FCN(NP,X,F)

!	Generating initial approximation to the Jacobian
      DO 2000 IP=1,NP
        X1=X(IP)
        DK=0.01D0*ABS(X1)
        IF(DK.EQ.0.0) DK=MAX(1.E-6,AEPS*100.)
        X(IP)=X1+DK
        CALL FCN(NP,X,WK)
        DO 1800 J=1,NP
1800    WK(J,NP+IP)=(WK(J,1)-F(J))/DK
        X(IP)=X1
2000  CONTINUE

      NUM=NP
      IFLG=0
      DO 2500 I=1,NP
        DO 2400 J=1,NP
2400    WK(J,I)=0.0
2500  WK(I,I)=1.0

!	To calculate H = B**(-1)
      CALL GAUELM(NP,NUM,WK(1,NP+1),WK,DET,IWK,LJ,IER1,IFLG)
      IF(IER1.GT.0) THEN
!	If GAUELM fails, then quit
        IER=441
        RETURN
      ENDIF

      DO 2800 I=1,NP
        WK(I,NP+1)=F(I)
        WK(I,NP+2)=X(I)
2800  CONTINUE

!	Broyden's iteration
      DO 5000 IPAS=1,200

!	Convergence check
        QCHK=.TRUE.
        DO 3500 J=1,NP
          S1=0.0
          DO 3200 K=1,NP
3200      S1=S1+WK(J,K)*WK(K,NP+1)
          IF(ABS(S1).GT.MAX(REPS*ABS(WK(J,NP+2)),AEPS)) QCHK=.FALSE.
          X(J)=WK(J,NP+2)-S1
          WK(J,NP+3)=-S1
3500    CONTINUE
        IF(QCHK) RETURN

        CALL FCN(NP,X,F)

!	Updating the inverse matrix
        DO 4000 J=1,NP
          WK(J,NP+4)=F(J)-WK(J,NP+1)
          WK(J,NP+1)=F(J)
          WK(J,NP+2)=X(J)
4000    CONTINUE

        SS=0.0
        DO 4300 J=1,NP
          S1=0.0
          S2=0.0
          DO 4100 K=1,NP
            S1=S1+WK(K,NP+3)*WK(K,J)
            S2=S2+WK(J,K)*WK(K,NP+4)
4100      CONTINUE
          F(J)=S1
          X(J)=S2-WK(J,NP+3)
          SS=SS+S1*WK(J,NP+4)
4300    CONTINUE

        DO 4400 J=1,NP
          DO 4400 K=1,NP
            WK(K,J)=WK(K,J)-X(K)*F(J)/SS
4400    CONTINUE
5000  CONTINUE

!	Iteration fails to converge
      IER=442
      END
