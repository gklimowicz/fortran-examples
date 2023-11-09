!     PROGRAM TO SOLVE NONLINEAR VOLTERRA EQUATION OF THE SECOND KIND
!     USING THE SIMPSON'S RULE

      PROGRAM INTEQ
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(665),F(665)
      EXTERNAL FG,FKER

!     EXAMPLE 13.6

51    FORMAT('   IER =',I4,5X,'STEP SIZE =',1PD14.6)
52    FORMAT('   X =',1PD14.6,5X,'F(X) =',3D14.6)

      REPS=1.D-8
      A=0.0

100   PRINT *,'TYPE H=STEP SIZE,  N=NO. OF STEPS  (QUITS WHEN N.LE.0)'
      READ *,H,N
      IF(N.LE.0) STOP

      CALL VOLT2(N,A,H,F,X,FG,FKER,REPS,IER)
      WRITE(6,51) IER,H
      DO 1000 I=1,N,5
        WRITE(6,52) X(I),F(I)
1000  CONTINUE
      GO TO 100
      END

!     -------------------------------------------------

!	To solve nonlinear Volterra equation of the second kind using
!	quadrature method with Simpson's rule
!
!	Integral[K(x,t,f(x)) dt] over [0,x] = f(x)+g(x)
!
!
!	N : (input) Number of points at which the solution is required
!	A : (input) Lower limit of the integral
!	H : (input) Uniform spacing to be used between abscissas
!	F : (output) Real array of length N containing the calculated solution
!		F(I) is the value of solution at X(I).
!	X : (output) Real array of length N containing the abscissas used
!		in the trapezoidal rule. The spacing is assumed to be uniform
!	FG : (input) Name of the function routine used to calculate the right
!		hand side g(X). 
!	FKER : (input) Name of the function routine used to calculate the
!		kernel K(x,t,f), note that in this case the "kernel" includes
!		the unknown solution also, because the equation is nonlinear.
!	REPS : (input) The required accuracy to which nonlinear equation is
!		solved using fixed point iteration. It will not control the
!		truncation error which is determined by the mesh spacing.
!	IER : (output) The error parameter, IER=0 implies successful execution
!		IER=712 implies N<3, in which case no calculations are done.
!		IER=752 implies that fixed point iteration failed at some stage.
!			No further calculations are done
!
!	FUNCTION FG(X) and FUNCTION FKER(X,T,F) must be supplied by the user.
!		FG is the right hand side function g(x) and FKER(X,T,F) is 
!		the kernel (the integrand in integral equation). 
!
!	Required routines : FG, FKER
!	
      SUBROUTINE VOLT2(N,A,H,F,X,FG,FKER,REPS,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=20,ETA=1.D-30)
      DIMENSION F(N),X(N)

      IER=0
      IF(N.LT.3) THEN
        IER=712
        RETURN
      ENDIF

!	Generating the starting values, F(1) and F(2)
      X(1)=A
      F(1)=-FG(A)
      G2=FG(A+H)
      FK1=FKER(A+H,A,F(1))
      F21=-G2+H*FK1
      F22=-G2+0.5*H*(FK1+FKER(A+H,A+H,F21))
      F23=0.5*(F(1)+F22)
      F24=-FG(A+0.5*H)+0.25*H*(FKER(A+0.5*H,A,F(1))+
     1     FKER(A+0.5*H,A+0.5*H,F23))
      X(2)=A+H
      F(2)=-G2+H*(FK1+4.*FKER(A+H,A+0.5*H,F24)+
     1      FKER(A+H,A+H,F22))/6.

!	Continuing the solution
      DO 2000 I=3,N,2
        X(I)=A+(I-1)*H
!	Use Simpson's 1/3 rule for odd I
        FI=FKER(X(I),A,F(1))
        DO 1200 J=2,I-2,2
          FI=FI+4.*FKER(X(I),X(J),F(J))
          FI=FI+2.*FKER(X(I),X(J+1),F(J+1))
1200    CONTINUE
        FI=(FI+4.*FKER(X(I),X(I-1),F(I-1)))*H/3.-FG(X(I))

!	The predictor
        F0=2.*F(I-1)-F(I-2)
!	Iteration on the corrector
        DO 1400 J=1,NIT
          F1=FI+H*FKER(X(I),X(I),F0)/3.
          R2=ABS(F1-F0)/(ABS(F1)+ABS(F1-F(I-1))+ETA)
          IF(R2.LT.REPS) GO TO 1500
          F0=F1
1400    CONTINUE
        IER=752
        RETURN

1500    F(I)=F1
        IF(I.GE.N) RETURN

!	For even I+1 use Simpson's 1/3 and 3/8 rule
        X(I+1)=A+I*H
        FI=FKER(X(I+1),A,F(1))
!	Use 1/3 rule on the first I-3 points
        DO 1600 J=2,I-4,2
          FI=FI+4.*FKER(X(I+1),X(J),F(J))
          FI=FI+2.*FKER(X(I+1),X(J+1),F(J+1))
1600    CONTINUE
        IF(I.GT.3) THEN
          FI=(FI+4.*FKER(X(I+1),X(I-3),F(I-3)))/3.
          FI=FI+17.*FKER(X(I+1),X(I-2),F(I-2))/24.
        ELSE
!	If I+1=4 use only 3/8 rule
          FI=3.*FI/8.
        ENDIF
        FI=H*(FI+9.*(FKER(X(I+1),X(I-1),F(I-1))+
     1      FKER(X(I+1),X(I),F(I)))/8.)-FG(X(I+1))

!	The predictor
        F0=2.*F(I)-F(I-1)
!	Iteration on the corrector
        DO 1700 J=1,NIT
          F1=FI+3.*H*FKER(X(I+1),X(I+1),F0)/8.
          R2=ABS(F1-F0)/(ABS(F1)+ABS(F1-F(I))+ETA)
          IF(R2.LT.REPS) GO TO 1800
          F0=F1
1700    CONTINUE
        IER=752
        RETURN

1800    F(I+1)=F1
2000  CONTINUE
      END

!     ------------------------------------------

      FUNCTION FG(X)
      IMPLICIT REAL*8(A-H,O-Z)

!     THE RHS FUNCTION

      FG=12.-13.*EXP(-X)
      END

!     ---------------------------------------

      FUNCTION FKER(X,T,F)
      IMPLICIT REAL*8(A-H,O-Z)

!     THE KERNEL

      FKER=12.*F
      END
