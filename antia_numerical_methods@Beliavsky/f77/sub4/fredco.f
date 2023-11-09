!	To solve linear Fredholm equation of first or second kind using
!	collocation method
!	It can be used to solve Volterra equations by defining the kernel
!	to be zero for t>x
!
!	N : (input) Number of points to be used in the collocation method
!	A : (input) Lower limit of the integral
!	B : (input) Upper limit of the integral
!	F : (output) Real array of length N containing the calculated coefficients
!		of expansion
!		solution = SUM_I F(I)*PHI(I,X)
!	X : (input) Real array of length N containing the points to
!		be used for collocation.
!	REPS : (input) Required relative accuracy to which integrals are to
!		be calculated
!	AEPS : (input) Required absolute accuracy to which integrals are to
!		be calculated
!		REPS and AEPS are passed on to subroutine ADPINT for calculating
!		the integrals when IQ=0. Otherwise these variables are not used.
!	WK : Real array of length N*N used as scratch space.
!	IWK : Integer array of length N used as a scratch space
!	IQ : (input) Integer variable to specify the treatment for integrals
!		PSI(I,X)=Integral[FKER(X,T)*PHI(I,T) dT] over [A,B]
!		If IQ=0, then the integrals are evaluated using subroutine
!			ADPINT, using function routine FUNK to calculate the
!			integrand, which in turn requires, PHI(I,T) and FKER(X,T). 
!		Otherwise the integrals are calculated using a user supplied
!			routine PSI(I,X).
!	IT : (input) Integer variable to specify the type of integral equation
!		If IT=1 Fredholm equation of the first kind is solved
!		If IT=2 Fredholm equation of the second kind is solved
!	IER : (output) The error parameter, IER=0 implies successful execution
!		IER=708 implies N<1, IT>2 or IT.LE.0, No calculations are done.
!		Other values of IER may be set by GAUELM and ADPINT
!
!	FUNCTION FG(X), FUNCTION PHI(I,X) and FUNCTION FKER(X,T) (for IQ=0)
!	or FUNCTION PSI(I,X) (for IQ.NE.0) must be supplied by the user.
!	Names of these function routines are fixed. FG(X) is the right
!	hand side function g(x), FKER(X,T) is the kernel, PHI(I,X) calculates
!	the basis functions phi_i(x), while PSI(I,X) calculates the integrals
!	defined above. The common block ZZFRED is used to pass on the variables
!	to FUNK. XI and II are the values of X and I.
!
!	Required routines : GAUELM, ADPINT, KRONRD, FUNK, FG, FKER, PHI, PSI
!	
      SUBROUTINE FREDCO(N,A,B,F,X,REPS,AEPS,WK,IWK,IQ,IT,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNK
!	To pass arguments to FUNK
      COMMON/ZZFRED/XI,II
      DIMENSION F(N),X(N),WK(N,N),IWK(N)

      IER=0
      IF(N.LT.1.OR.IT.GT.2.OR.IT.LE.0) THEN
        IER=708
        RETURN
      ENDIF
      NMAX=10000

!	Setting up the system of linear equations
      DO 2000 J=1,N
        F(J)=FG(X(J))
        DO 2000 I=1,N
          IF(IQ.EQ.0) THEN
!	Evaluate the integrals numerically, split the range into two 
!	to tackle possible discontinuity at t=X
            XI=X(I)
            II=J
            CALL ADPINT(RI,A,XI,REPS,AEPS,DIF,FUNK,IER,NPT,NMAX)
            IF(IER.GT.100) RETURN
            CALL ADPINT(RI1,XI,B,REPS,AEPS,DIF1,FUNK,IER,NPT1,NMAX)
            WK(I,J)=RI+RI1
            IF(IER.GT.100) RETURN
          ELSE
!	Calculate the integrals PSI(I,X) using user supplied routine
            WK(I,J)=PSI(J,X(I))
          ENDIF
          IF(IT.EQ.2) WK(I,J)=WK(I,J)-PHI(J,X(I))
2000  CONTINUE

!	Solve the resulting system of linear equations
      NUM=1
      LJ=N
      IFLG=0
      CALL GAUELM(N,NUM,WK,F,DET,IWK,LJ,IER,IFLG)
      END
