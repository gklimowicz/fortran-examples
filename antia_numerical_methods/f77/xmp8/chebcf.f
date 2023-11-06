!     PROGRAM TO CONVERT POWER SERIES TO CHEBYSHEV EXPANSION AND VICE VERSA

      PROGRAM CHEBEX
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(50),C(50)

51    FORMAT('  POWER SERIES COEF. :',1P3D14.6/(2X,5D14.6))
52    FORMAT('  COEF. OF CHEBYSHEV EXPANSION :',1P3D14.6/(2X,5D14.6))
53    FORMAT('  N =',I4,5X,'IFLG =',I4)

100   PRINT *,'TYPE N=DEGREE, IFLG=0/1'
      PRINT *,'IFLG=0 TO CONVERT POWER SERIES TO CHEB. EXPANSION'
      PRINT *,'             (QUITS WHEN N.LE.0)'
      READ *,N,IFLG
      IF(N.LE.0) STOP
      IF(IFLG.EQ.0) THEN
        PRINT *,'TYPE POWER SERIES COEFFICIENTS'
        READ *,(P(I),I=1,N+1)
        WRITE(6,51) (P(I),I=1,N+1)
      ELSE
        PRINT *,'TYPE COEFFICIENTS FOR CHEBYSHEV EXPANSION'
        READ *,(C(I),I=1,N+1)
        WRITE(6,52) (C(I),I=1,N+1)
      ENDIF

      CALL CHEBCF(N,C,P,IFLG)
      PRINT 53,N,IFLG
      IF(IFLG.EQ.0) THEN
        WRITE(6,52) (C(I),I=1,N+1)
      ELSE
        WRITE(6,51) (P(I),I=1,N+1)
      ENDIF
      GO TO 100
      END

!     -------------------------------------------------------

!	To calculate the coefficients of Chebyshev expansion from those of
!	power series or vice versa.
!
!	N : (input) The degree of polynomial
!	C : (input/output) Real array of length N+2 containing the coefficients
!		of Chebyshev series, C(I+1) is the coefficient of T_I(x) in
!		expansion. These coefficients would be calculated if IFLG=0,
!		otherwise they must be supplied
!		In the latter case the contents of C are destroyed
!	P : (input/output) Real array of length N+2 containing the coefficients
!		of power series. P(I+1) is the coefficient of x**I.
!		These coefficients must be supplied if IFLG=0,
!		otherwise they would be calculated.
!	IFLG : (input) Flag to decide which coefficients need to be calculated
!		If IFLG=0 then coefficients in Chebyshev expansion are
!		calculated. In this case P(I) must be supplied.
!		otherwise coefficients in power series expansion are
!		calculated. In this case C(I) must be supplied and contents
!		of array C are destroyed during calculations
!
!	Required routines : None
!
      SUBROUTINE CHEBCF(N,C,P,IFLG)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(N+2),P(N+2)

      IF(IFLG.EQ.0) THEN
!	Calculate coefficients of Chebyshev expansion
        DO 2000 I=2,N+2
2000    C(I)=0.0
        C(1)=P(N+1)
!	Nested multiplication
        DO 2400 I=1,N
          CJ=C(2)
          C(2)=C(1)+0.5*C(3)
          C(1)=P(N+1-I)+0.5*CJ
          IF(I.GT.1) THEN
            DO 2200 J=3,I+1
              CJ1=0.5*(CJ+C(J+1))
              CJ=C(J)
              C(J)=CJ1
2200        CONTINUE
          ENDIF
2400    CONTINUE

      ELSE
!	Calculate coefficients in power series
        IF(N.LE.1) THEN
          P(1)=C(1)
          IF(N.EQ.1) P(2)=C(2)
          RETURN
        ENDIF
        DO 4000 I=1,N-2
          CJ=C(N-I)
          C(N-I)=2.*C(N+1-I)
          C(N+1-I)=2.*C(N+2-I)
          IF(I.LE.N-3) THEN
            DO 3400 J=N-I-1,2,-1
              CJ1=2.*CJ-C(J+2)
              CJ=C(J)
              C(J)=CJ1
3400        CONTINUE
          ENDIF
          CJ1=CJ-0.5*C(3)
          CJ=C(1)
          C(1)=CJ1
          P(I)=CJ-0.5*C(2)
4000    CONTINUE
        CJ=C(1)
        C(1)=C(2)
        C(2)=2.*C(3)
        P(N-1)=CJ-0.5*C(2)
        P(N)=C(1)
        P(N+1)=C(2)
      ENDIF
      END
