!	To find contour integrals over a circle as required by subroutine DELVES
!
!	CF : (input) Name of the function routine to calculate the function value
!	CZ : (input) Complex variable specifying the centre of the required circle
!	RAD : (input) Real variable specifying the radius of the circle.
!	NZ : (output) Number of zeros located by the subroutine inside the circle
!	CS : (output) Complex array of length NMAX containing the values
!		of contour integrals
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=44 implies that the first integration to determine NZ
!			did not converge to satisfactory accuracy, but its
!			magnitude is less than 0.5. NZ is set to zero in this case.
!		IER=45 implies that the first integration to determine NZ
!			did not converge to satisfactory accuracy, but a
!			weaker criterion is satisfied.
!		IER=46 implies that one of the higher integrals to determine
!			CS(I) did not converge to satisfactory accuracy, but a
!			weaker criterion is satisfied.
!		IER=434 implies that NZ>NMAX and no further calculations are done
!		IER=437 implies that the first integration to determine NZ
!			did not converge satisfactorily, no further calculations 
!			are done.
!	NMAX : (input) Maximum number of zeros to be located. If NZ>NMAX
!		no further calculations are done.
!	
!	FUNCTION CF(CX,CDF) must be supplied by the user. Here CF, CX and CDF are
!		complex variables and CDF is the first derivative of CF.
!
!	Required routines : CF

      SUBROUTINE CONTUR(CF,CZ,RAD,NZ,CS,IER,NMAX)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A,B,D-H,O,P,R-Z)
      PARAMETER(REPS0=1D-2,REPS=1D-4,NPT=1024,PI=3.14159265358979324D0)
      PARAMETER(CI=(0.D0,1.D0))
      DIMENSION CFUN(NPT),CZP(NPT),CS(NMAX),CSM(NPT)

      IER=0
      NP=16
      CXP=EXP(CI*2.*PI/NP)
      CZ0=1.
      CSUM=0.0
      CNT=0.0
      II=0

!	Trapezoidal rule for the first integral
!	I=1,7 will try NP from 16 to 1024=NPT
!	To increase number of points increase NPT and the loop count
      DO 2000 I=1,7
        DO 1000 J=1,NP
          II=II+1
          CZP(II)=CZ0
          CFN=CF(RAD*CZ0+CZ,CDF)
          CFUN(II)=CZ0*CDF/CFN
          CZ0=CZ0*CXP
          CSUM=CSUM+CFUN(II)
1000    CONTINUE

        IF(I.EQ.1) THEN
          CZ0=EXP(CI*PI/NP)
        ELSE
          NP=NP*2
          CZ0=EXP(CI*PI/NP)
          CXP=CZ0*CZ0
        ENDIF

!	Latest approximation to the integral
        CINT=CSUM*RAD/NP
        NZ=CINT+0.4D0
        EPS=ABS(CINT-CNT)
        EPS1=ABS(CINT-NZ)
        IF(EPS.LT.REPS0.AND.EPS1.LT.0.1D0) GO TO 2500
        CNT=CINT
2000  CONTINUE

      IF(ABS(CINT).LT.0.5) THEN
!	The integral is approximately zero
        IER=44
        NZ=0
        RETURN
      ENDIF

      IF(EPS1.GT.0.1D0) THEN
!	Integral fails to converge
        IER=437
        RETURN
      ENDIF
!	Accept the value of NZ, although integral hasn't converged properly
      IER=45

2500  NZ=CINT+0.4D0
      IF(NZ.LE.0) RETURN
      IF(NZ.GT.NMAX) THEN
!	Too many zeros inside the circle
        IER=434
        RETURN
      ENDIF

      RMAX=MAX(ABS(CZ),RAD)
      DO 3000 I=1,NZ
!	Contour integrals using trapezoidal rule
        CSUM=0.0
        DO 2800 J=1,NP
          CSUM=CSUM+CFUN(J)*CZP(J)**I
2800    CONTINUE
        CSM(I)=CSUM
        CS(I)=CSUM*RAD**(I+1)/NP
3000  CONTINUE
      II=NP

!	Double the number of points until convergence
3100  CONTINUE
      DO 3400 J=1,NP
        CFN=CF(RAD*CZ0+CZ,CDF)
        CFUNI=CZ0*CDF/CFN
        DO 3200 K=1,NZ
          CSM(K)=CSM(K)+CFUNI*CZ0**K
3200    CONTINUE
        CZ0=CZ0*CXP
3400  CONTINUE

      NP=NP*2
      CZ0=EXP(CI*PI/NP)
      CXP=CZ0*CZ0

!	The convergence test
      QC=.TRUE.
      DO 3600 K=1,NZ
        CINT=CSM(K)*RAD**(K+1)/NP
        EPS=ABS(CINT-CS(K))
        IF(EPS.GT.REPS*MAX(ABS(CINT),RMAX**K)) QC=.FALSE.
        CS(K)=CINT
3600  CONTINUE

      IF(QC) RETURN
      IF(NP.LE.16*NPT) GO TO 3100

!	Integrals fail to converge
      IER=46
      END
