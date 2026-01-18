!     PROGRAM TO SOLVE LINEAR VOLTERRA EQUATIONS USING TRAPEZOIDAL RULE

      PROGRAM INTEQ
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(665),F(665)
      EXTERNAL FG,FKER

!     EXAMPLE 13.7 : VOLTERRA EQUATION OF THE FIRST KIND

51    FORMAT('    IER =',I4,5X,'STEP SIZE =',1PD14.6,5X,'IT =',I3)
52    FORMAT('   X =',1PD14.6,5X,'F(X) =',3D14.6)

      A=0

100   PRINT *,'TYPE H=STEP SIZE,  N=NO. OF STEPS,  IT=+1/-1'
      PRINT *,'IT=-1  FOR SMOOTHED SOLUTION   (QUITS WHEN N.LE.0)'
      READ *,H,N,IT
      IF(N.LE.0) STOP

      CALL VOLT(N,A,H,F,X,FG,FKER,IT,IER)
      WRITE(6,51) IER,H,IT

!     PRINT SOLUTION AT EVERY FIFTH POINT

      DO 1000 I=1,N,5
        WRITE(6,52) X(I),F(I)
1000  CONTINUE
      GO TO 100
      END

!     -------------------------------------------------

!	To solve linear Volterra equations using quadrature method with the
!	trapezoidal rule
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
!		kernel K(x,t)
!	IT : (input) Integer variable to specify the type of integral equation
!		If IT=1 Volterra equation of the first kind is solved
!		If IT=-1 Volterra equation of the first kind is solved and
!			computed values are smoothed as explained in Section 12.8
!		If IT=2 Volterra equation of the second kind is solved
!	IER : (output) The error parameter, IER=0 implies successful execution
!		IER=712 implies N<3, in which case no calculations are done.
!		IER=751 implies that denominator is zero at some stage
!			No further calculations are done
!
!	FUNCTION FG(X) and FUNCTION FKER(X,T) must be supplied by the user.
!		FG is the right hand side function g(x) and FKER(X,T) is 
!		the kernel. 
!
!	Required routines : FG, FKER
!	
      SUBROUTINE VOLT(N,A,H,F,X,FG,FKER,IT,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION F(N),X(N)

      IER=0
      IF(N.LT.3) THEN
        IER=712
        RETURN
      ENDIF

      X(1)=A
      IF(IT.EQ.2) THEN
!	Starting value for equations of second kind
        F(1)=-FG(A)
      ELSE
        DI=H*FKER(A,A)
        IF(DI.EQ.0.0) THEN
          IER=751
          RETURN
        ENDIF
!	Starting value for equations of first kind
        F(1)=(FG(A+H)-FG(A-H))/(2.*DI)
      ENDIF

!	Continue the solution using the trapezoidal rule
      DO 2000 I=2,N
        X(I)=A+(I-1)*H
        FI=0.5*FKER(X(I),A)*F(1)
        DO 1200 J=2,I-1
1200    FI=FI+FKER(X(I),X(J))*F(J)
        FI=FI*H-FG(X(I))

        DI=-H*FKER(X(I),X(I))/2.
        IF(IT.EQ.2) DI=DI+1
        IF(DI.EQ.0.0) THEN
          IER=751
          RETURN
        ENDIF
        F(I)=FI/DI
2000  CONTINUE

      IF(IT.EQ.-1) THEN
!	Apply smoothing for F(2),...,F(N-1)
        FI=0.25*(F(1)+2.*F(2)+F(3))
        DO 3000 I=2,N-2
          FI1=0.25*(F(I)+2.*F(I+1)+F(I+2))
          F(I)=FI
          FI=FI1
3000    CONTINUE
        F(N-1)=FI
      ENDIF
      END

!     ------------------------------------------

      FUNCTION FG(X)
      IMPLICIT REAL*16(A-H,O-Z)

!     RHS FUNCTION

      FG=SINH(X)
      END

!     ---------------------------------------

      FUNCTION FKER(X,T)
      IMPLICIT REAL*16(A-H,O-Z)

!     THE KERNEL

      FKER=EXP(X-T)
      END
