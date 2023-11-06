!	To calculate coefficients of cubic spline interpolation with
!		not-a-knot boundary conditions
!
!	X : (input) Real array of length N containing x values
!	F : (input) Real array of length N containing values of function at X(I)
!		F(I) is the tabulated function value at X(I).
!	N : (input) Length of table X,F
!	C : (output) Real array of length 3*N containing the spline coefficients
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=201 implies that N<2
!
!	Required routines : None

      SUBROUTINE SPLINE(X,F,N,C,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N),C(3,N)

      IER=0
      IF(N.LE.1) THEN
        IER=201
        RETURN

      ELSE IF(N.EQ.2) THEN
!	Use linear interpolation
        C(1,1)=(F(2)-F(1))/(X(2)-X(1))
        C(2,1)=0.0
        C(3,1)=0.0
        RETURN

      ELSE IF(N.EQ.3) THEN
!	Use quadratic interpolation
        DIV12=(F(2)-F(1))/(X(2)-X(1))
        DIV23=(F(3)-F(2))/(X(3)-X(2))
        C(3,1)=0.0
        C(3,2)=0.0
        C(2,1)=(DIV23-DIV12)/(X(3)-X(1))
        C(2,2)=C(2,1)
        C(1,1)=DIV12+C(2,1)*(X(1)-X(2))
        C(1,2)=DIV23+C(2,1)*(X(2)-X(3))
        RETURN

      ELSE
!	Use cubic splines 

!	Setting up the coefficients of tridiagonal matrix
        C(3,N)=(F(N)-F(N-1))/(X(N)-X(N-1))
        DO 1000 I=N-1,2,-1
          C(3,I)=(F(I)-F(I-1))/(X(I)-X(I-1))
          C(2,I)=2.*(X(I+1)-X(I-1))
!	The right hand sides
1000    C(1,I)=3.*(C(3,I)*(X(I+1)-X(I))+C(3,I+1)*(X(I)-X(I-1)))

!	The not-a-knot boundary conditions
        C1=X(3)-X(1)
        C(2,1)=X(3)-X(2)
        C(1,1)=C(3,2)*C(2,1)*(2.*C1+X(2)-X(1))+C(3,3)*(X(2)-X(1))**2
        C(1,1)=C(1,1)/C1
        CN=X(N)-X(N-2)
        C(2,N)=X(N-1)-X(N-2)
        C(1,N)=C(3,N)*C(2,N)*(2.*CN+X(N)-X(N-1))
        C(1,N)=(C(1,N)+C(3,N-1)*(X(N)-X(N-1))**2)/CN

!	Solving the equation by Gaussian elimination
        G=(X(3)-X(2))/C(2,1)
        C(2,2)=C(2,2)-G*C1
        C(1,2)=C(1,2)-G*C(1,1)
        DO 2000 J=2,N-2
          G=(X(J+2)-X(J+1))/C(2,J)
          C(2,J+1)=C(2,J+1)-G*(X(J)-X(J-1))
2000    C(1,J+1)=C(1,J+1)-G*C(1,J)
        G=CN/C(2,N-1)
        C(2,N)=C(2,N)-G*(X(N-1)-X(N-2))
        C(1,N)=C(1,N)-G*C(1,N-1)

!	The back-substitution
        C(1,N)=C(1,N)/C(2,N)
        DO 3000 I=N-1,2,-1
3000    C(1,I)=(C(1,I)-C(1,I+1)*(X(I)-X(I-1)))/C(2,I)
        C(1,1)=(C(1,1)-C(1,2)*C1)/C(2,1)

!	Calculating the coefficients of cubic spline
        DO 4000 I=1,N-1
          C(2,I)=(3.*C(3,I+1)-2.*C(1,I)-C(1,I+1))/(X(I+1)-X(I))
          C(3,I)=(C(1,I)+C(1,I+1)-2.*C(3,I+1))/(X(I+1)-X(I))**2
4000    CONTINUE
!	Set the coefficients for interval beyond X(N) using continuity
!	of second derivative, although they may not be used.
        C(2,N)=C(2,N-1)+3*(X(N)-X(N-1))*C(3,N-1)
        C(3,N)=0.0
      ENDIF
      END
