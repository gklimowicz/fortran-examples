!     To find coefficients of Chebyshev expansion
 
      PROGRAM CHEBEXP
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUN
      DIMENSION C(5000)
 
!     Chebyshev coefficients for ARC TAN(X)
 
51    FORMAT('   IER =',I4,5X,'N =',I5,5X,'COEFFICIENTS :'/(1P5D14.6))
52    FORMAT('   EXACT COEFFICIENTS (ALTERNATE VALUES) :'/(1P5D14.6))
 
!     The number must be somewhat larger than the required coefficients
100   PRINT *,'TYPE NO. OF COEFFICIENTS TO BE CALCULATED'
      PRINT *,'                       (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0) STOP
      CALL CHEBEX(N,C,FUN,IER)
!     Only the first 20 coefficients are printed
      WRITE(6,51) IER,N,(C(I),I=1,20)
!     Calculate the exact coefficients for expansion
      RM=SQRT(2.D0)-1.0
      WRITE(6,52) ((-1)**(I/2)*2*RM**I/I,I=1,20,2)
      GO TO 100
      END
 
!     -------------------------------------------------------
 
!	To calculate the coefficients of Chebyshev expansion of a given
!	function, which can be evaluated at any required point
!	This subroutine uses orthogonality relation over discrete points
!	to compute the coefficients.
!
!	N : (input) Number of coefficients required, this number should
!		be much larger than the actual number required. The
!		accuracy of computed coefficients increases with N.
!	C : (output) Real array of length N containing the required
!		coefficients. The constant term is C(1)/2, while C(I+1)
!		is the coefficient of T_I(x).
!	FUN : (input) Name of function routine to compute the given function
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=613 implies N<10, in which case no calculations are done
!	
!	There is no test for accuracy and this has to be ascertained by
!	doing another calculation with larger N (say 2N).
!
!	FUNCTION FUN(X) must be supplied by the user
!
!	Required routines : FUN
 
      SUBROUTINE CHEBEX(N,C,FUN,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.14159265358979324D0)
      DIMENSION C(N)
 
      IF(N.LT.10) THEN
        IER=613
        RETURN
      ENDIF
      IER=0
 
      DO 2000 I=1,N
        C(I)=0.0
2000  CONTINUE
 
!     evaluate the sums over n points
      DO 3000 I=1,N
        X=COS((2*I-1)*PI/(2*N))
        FX=FUN(X)
        C(1)=C(1)+FX
        C(2)=C(2)+FX*X
        T0=1.0
        T1=X
        DO 2500 J=3,N
          T2=2*X*T1-T0
          C(J)=C(J)+FX*T2
          T0=T1
          T1=T2
2500    CONTINUE
3000  CONTINUE
 
      DO 3500 I=1,N
        C(I)=2.0*C(I)/N
3500  CONTINUE
      END
 
!     ---------------------------------------------------
 
      FUNCTION FUN(X)
      IMPLICIT REAL*8(A-H,O-Z)

!	The required function whose Chebyshev expansion needs to be calculated
      FUN=DATAN(X)
      END
