!     To calculate weights and abscissas of a quadrature formula with
!     specified weight function when recurrence relation for orthogonal
!     polynomial is known.
!
!     N : (input) Number of points in the required quadrature formula
!     W : (output) Array of length N, which will contain the weights
!     AB : (output) Array of length N containing the abscissas
!     COF : (input) Array of length 3*N containing the coefficients of
!		the recurrence relation for orthogonal polynomials
!		P_i(x)=(COF(1,i)*x+COF(2,i))P_{i-1}(x) - COF(3,i)*P_{i-2}(x)
!     RI0 : (input) The integral of weight function over the required
!		interval.
!     IER : (output) Error parameter, IER=0 implies successful execution
!		IER=302 implies N.LE.0
!		IER=321 implies that some coefficient becomes imaginary
!			during calculations.
!		In both these cases calculations are abandoned.
!		Other values may be set by TQL2
!     WK : Real array of length N*(N+2) used as scratch space
!
!     Required routines : TQL2
 
      SUBROUTINE GAUSRC(N,W,AB,COF,RI0,IER,WK)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER (REPS=1.Q-33)
!      PARAMETER (REPS=1.D-15)
!	For REAL*4 use REPS=6.E-8
!      PARAMETER (REPS=6.E-8)
      DIMENSION WK(N,N+2),COF(3,N),W(N),AB(N)
 
      IER=302
      IF(N.LE.0) RETURN
      LJ=N
 
!     Calculate the coefficients of symmetric tridiagonal matrix
      DO 2000 I=1,N
        WK(I,N+1)=-COF(2,I)/COF(1,I)
        IF(I.LT.N) THEN
          R1=COF(3,I+1)/(COF(1,I)*COF(1,I+1))
          IF(R1.GE.0.0) THEN
            WK(I+1,N+2)=SQRT(R1)
          ELSE
            IER=321
            RETURN
          ENDIF
        ENDIF
        DO 1800 J=1,N
          WK(J,I)=0.0
1800    CONTINUE
        WK(I,I)=1.0
2000  CONTINUE
 
!     Find eigenvalues and eigenvectors of the tridiagonal matrix
      CALL TQL2(WK,N,LJ,WK(1,N+1),WK(1,N+2),REPS,IER)
      IF(IER.GT.0) RETURN
 
!     Calculate the abscissas and weights
      DO 3000 I=1,N
        AB(I)=WK(I,N+1)
        W(I)=WK(1,I)**2*RI0
3000  CONTINUE
 
      END
