!	Real root of a system of nonlinear equations using Newton's method
!
!	FCN : (input) Name of the subroutine to calculate the vector function
!	NP : (input) Number of variables and equations to be solved
!	X : (input/output) Real array of length NP containing the starting
!		values. After execution it should contain the computed roots
!	F : (output) Real array of length NP containing the function values
!		at the point specified by array X
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X(I)))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=441 implies that GAUELM failed to solve the required
!			system of linear equations
!		IER=442 implies that iteration did not converge to required
!			accuracy
!	WK : Real array of length NP*(NP+1) used as scratch space
!	IWK : integer array of length NP used as scratch space
!
!	SUBROUTINE FCN(NP,X,F,DF) must be supplied by the user.
!	Here X and F are real arrays of length NP containing the value of
!	vector function at point X. DF is a real array of length NP*NP 
!	containing the Jacobian. The first dimension of DF in FCN should be NP.
!	The system of equations to be solved are F(I)=0
!
!	Required routines : GAUELM, FCN
!
      SUBROUTINE NEWTON(FCN,NP,X,F,REPS,AEPS,IER,WK,IWK)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      DIMENSION X(NP),F(NP),WK(NP*(NP+1)),IWK(NP)

      LJ=NP
      NUM=1
      IER=0
      NP2=NP*NP

!	Loop for Newton's iteration
      DO 5000 IPAS=1,200
        CALL FCN(NP,X,F,WK)
!	The right hand side of system of equations for Newton's iteration
        DO 2000 I=1,NP
2000    WK(NP2+I)=F(I)

        IFLG=0
        CALL GAUELM(NP,NUM,WK,WK(NP2+1),DET,IWK,LJ,IER1,IFLG)
        IF(IER1.GT.0) THEN
!	If GAUELM fails, then quit
          IER=441
          RETURN
        ENDIF

!	Convergence check
        QCHK=.TRUE.
        DO 3500 J=1,NP
          IF(ABS(WK(NP2+J)).GT.MAX(REPS*ABS(X(J)),AEPS)) QCHK=.FALSE.
          X(J)=X(J)-WK(NP2+J)
3500    CONTINUE
        IF(QCHK) RETURN

5000  CONTINUE

!	Iteration fails to converge
      IER=442
      END
