!	Real root of a system of nonlinear equations using Davidenko's method
!	For use with subroutine BROYDN
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
!		IER=407 implies NP > NMAX or NP < 1, in which case no
!			calculations are done
!		IER=440 implies that THETA is nonzero even after 100 steps
!		Other values may be set by subroutines NEWTON or BROYDN
!	WK : Scratch real array of length NP*(NP+1) (for use with NEWTON) or
!		length MAX(2*NP*NP,NP*(NP+4)) (for use with BROYDN) 
!	IWK : integer array of length NP used as scratch space
!	THETA : which is read by the subroutine is the parameter introduced
!		into the equations. Normally THETA=1 to start with and
!		initial guess should satisfy the equations with this parameter
!		THETA=0 corresponds to the equations that need to be solved.
!		The value of THETA should be gradually reduced from 1 to 0
!
!	SUBROUTINE FCN(NP,X,F) or SUBROUTINE FCN(NP,X,F,DF)
!	must be supplied by the user. Here X and F are real arrays of length
!	NP with F containing the value of vector function at point X. DF is
!	real array of length NP*NP containing the Jacobian which is required
!	by subroutine NEWTON. If subroutine BROYDN is used then DF is not
!	required. The subroutine FCN should parameterise the equations in
!	terms of THETA as explained in section 7.16. A sample routine is
!	appended at the end of this file.
!
!	Required routines : BROYDN (or NEWTON), GAUELM, FCN
!
      SUBROUTINE DAVIDN_B(FCN,NP,X,F,REPS,AEPS,IER,WK,IWK)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMAX=200)
      EXTERNAL FCN
      DIMENSION X(NP),F(NP),WK(*),IWK(NP)
!	To pass on parameters to FCN this common block must be included
!	in subroutine FCN
      COMMON/DNFUN/THETA,X0(NMAX),F0(NMAX)

51    FORMAT(/5X,'THETA =',F6.2/)
52    FORMAT('   IER =',I4,'    X =',1P4D14.6/(6X,5D14.6))

      IF(NP.GT.NMAX.OR.NP.LT.1) THEN
        IER=407
        RETURN
      ENDIF

      PRINT *,' TYPE IN THE INITIAL GUESS'
      READ *,(X(I),I=1,NP)

!	To initialise the parameterised function
      THETA=0.0
!	For NEWTON use the following form
!      CALL FCN(NP,X,F0,WK)
      CALL FCN(NP,X,F0)
      DO 500 I=1,NP
500   X0(I)=X(I)
      REPS0=MAX(1.Q-3,REPS)
      AEPS0=MAX(1.Q-4,AEPS)

      DO 1000 I=1,100
        PRINT *,' TYPE THETA'
        READ *,THETA
        PRINT 51, THETA
        IF(THETA.EQ.0.0) THEN
          REPS0=REPS
          AEPS0=AEPS
        ENDIF

!	Use either NEWTON or BROYDN
!        CALL NEWTON(FCN,NP,X,F,REPS0,AEPS0,IER,WK,IWK)
        CALL BROYDN(FCN,NP,X,F,REPS0,AEPS0,IER,WK,IWK)
        PRINT 52,IER,(X(J),J=1,NP)
        IF(IER.GT.0) RETURN
        IF(THETA.EQ.0.0) RETURN
1000  CONTINUE

      IER=440
      END

!	-------------------
!
!	For use with NEWTON include DF
!
!	SUBROUTINE FCN(NP,X,F)
!	IMPLICIT REAL*16(A-H,O-Z)
!	PARAMETER (NMAX=200)
!	DIMENSION X(NP),F(NP)
!	COMMON/DNFUN/THETA,X0(NMAX),F0(NMAX)
!
!	call the required function to calculate the vector function
!	There is no provision to pass the name FUN, it has to be put explicitly
!	CALL FUN(NP,X,F)
!
!	APPLY PARAMETRISATION IN TERMS OF THETA
!
!	DO 2222 I=1,NP
!	F(I)=F(I)-THETA*F0(I)
!2222	CONTINUE
!	END
