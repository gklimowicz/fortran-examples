!	Complex zeros of analytic function inside a given circle
!
!	CZ : (input) Complex variable specifying the centre of the required circle
!	RAD : (input) Real variable specifying the radius of the circle.
!		All zeros inside the circle with centre CZ and radius RAD
!		need to be found.
!	CF : (input) Name of the function routine to calculate the function value
!	NZ : (output) Number of zeros located by the subroutine
!	CZERO : (output) Complex array of length NMAX (=5) containing the
!		zeros determined by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=405 implies RAD.LE.0, in which case no calculations are done
!		IER=434 implies that NZ>NMAX and zeros are not determined
!		IER=435 implies that iteration failed to converge to specified
!			accuracy for at least one zero
!		IER=436 implies that POLYC failed to find all roots of the
!			required polynomial
!		Other values of IER may be set by CONTUR
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than max(AEPS, REPS*ABS(CZERO))
!	
!	FUNCTION CF(CX,CDF) must be supplied by the user. Here CF, CX and CDF are
!		all complex variables, and CDF is the first derivative of CF.
!	Required routines : CONTUR, NEWRAC, POLYC, LAGITC, CF
!
      SUBROUTINE DELVES(CZ,RAD,CF,NZ,CZERO,IER,REPS,AEPS)
      IMPLICIT COMPLEX*32(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*16(A,B,D-H,O,P,R-Z)
      PARAMETER(NMAX=5)
      EXTERNAL CF
      DIMENSION CZERO(NMAX),CS(NMAX),COF(NMAX+1),CWK(NMAX+1)

      NZ=0
      IER=405
      IF(RAD.LE.0.0) RETURN

!	To calculate the contour integrals
      CALL CONTUR(CF,CZ,RAD,NZ,CS,IER,NMAX)
      IF(IER.GT.100) RETURN
      IF(NZ.LE.0) RETURN
      IF(NZ.GT.NMAX) THEN
!	If the number of zeros is too large, then quit
        IER=434
        RETURN
      ENDIF

!	Coefficients of the equivalent polynomial
      COF(NZ+1)=1.0
      DO 2000 N=NZ,1,-1
        J=NZ+1-N
        C1=CS(J)
        DO 1500 I=1,J-1
1500    C1=C1+CS(J-I)*COF(NZ+1-I)
        COF(N)=-C1/J
2000  CONTINUE

!	Find zeros of the polynomial
      QREFIN=.FALSE.
      CALL POLYC(NZ,COF,CZERO,IER1,QREFIN,CWK)
      IF(IER1.GT.100) IER=436
      DO 3000 I=1,NZ
3000  CZERO(I)=CZERO(I)+CZ

!	Refine the zeros by iterating on the original function
      DO 4000 I=1,NZ
        CALL NEWRAC(CZERO(I),REPS,AEPS,IER1,CF)
        IF(IER1.GT.100) IER=435
4000  CONTINUE
      END
