!	To solve initial value problems in ordinary differential equations
!	using extrapolation method
!
!	N : (input) Number of first order differential equations to be solved
!	Y : (input/output) Real array of length N containing the initial
!		values of variables. After execution it will contain the
!		values of variable at the last point where the integration
!		has been successful.
!	DY : (output) Real array of length N containing the derivatives
!		of Y at the last point
!	DIF : (input) Name of subroutine to calculate the right hand side
!		of differential equation y'=f(t,y)
!	H : (input/output) Initial guess for the step size. After execution
!		it will contain the step size used by the program
!	T0 : (input/output) Initial value of independent variable t, at
!		which initial values are specified. After execution it will
!		be set to the point up to which integration has been successful
!	TN : (input) The final value of t at which the solution is required.
!		If integration is successful T0 will be set equal to TN.
!		Intermediate values will not be preserved so if solution
!		is required at intermediate points, TN must be set to first
!		such value and multiple calls will be needed to calculate
!		all required values. For each subsequent call only TN needs
!		to be changed.
!	REPS : (input) Required accuracy in each component of the solution.
!		The subroutine only controls local truncation error and hence
!		actual error could be larger
!	NSTEP : (output) Number of calls to DIF made by the subroutine to
!		complete the integration during each call to EXTP.
!	NMAX : (input/output) Maximum number of function evaluations to be used.
!		If NMAX.LE.0 it will be set to a default value of NMX=100000.
!	IER : (output) Error parameter; IER=0 implies successful execution
!		IER=703 implies N.LE.0, no calculations are done
!		IER=730 implies that step-size has become smaller than
!			REPS*|TN-T0|
!		IER=731 implies that step size is too small for arithmetic used
!		IER=732 implies that integration could not be completed in
!			the specified number of steps.
!		IER=733 implies that denominator for evaluating rational
!			function extrapolation vanished
!	WK : Real array of length 39N used as scratch space
!	IFLG : (input) Integer variable used as a flag to decide the type
!		of extrapolation to be used.
!		If IFLG=0 polynomial extrapolation is used
!		otherwise rational function extrapolation is used
!
!	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
!	the differential equation. T is the value of independent variable,
!	N is the number of variables, Y is a real array of length N containing
!	the values of variables. DY is a real array of length N which should
!	contain the calculated values of derivatives at (T,Y).
!
!	Required routines : DIF
!	
!	
      SUBROUTINE EXTP(N,Y,DY,DIF,H,T0,TN,REPS,NSTEP,NMAX,IER,WK,IFLG)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMX=100000,SFAC=0.8Q0,EPS=1.Q-30)
      DIMENSION Y(N),DY(N),WK(N,39),NSEQ(12),HI(12)
      DATA NSEQ/2,4,6,8,12,16,24,32,48,64,96,128/

      IF(N.LE.0) THEN
        IER=703
        RETURN
      ENDIF
      IF(NMAX.LE.0) NMAX=NMX
      NSTEP=0
      TSTEP=TN-T0
      IER=0
      IF(TSTEP.EQ.0.0) RETURN
      IF(H.EQ.0.0.OR.ABS(H).GT.ABS(TSTEP)) H=TSTEP
      IF(H.LT.0.0.EQV.TN.GT.T0) H=-H

!	Loop for integration
1000  NSTEP=NSTEP+1
      CALL DIF(T0,N,Y,DY)

!	Loop for extrapolation
      DO 2500 IE=1,12
!	Step size for the midpoint method
        H2=H/NSEQ(IE)
        HI(IE)=H2**2
        DO 1200 I=1,N
          WK(I,1)=Y(I)
1200    WK(I,2)=Y(I)+H2*DY(I)

        T1=T0
        DO 1500 J=1,NSEQ(IE)-1
          T1=T1+H2
          CALL DIF(T1,N,WK(1,2),WK(1,3))
          DO 1400 I=1,N
!	The midpoint rule
            C1=WK(I,1)+2.*H2*WK(I,3)
            WK(I,1)=WK(I,2)
1400      WK(I,2)=C1
1500    CONTINUE
        T1=T1+H2
        CALL DIF(T1,N,WK(1,2),WK(1,3))
        NSTEP=NSTEP+NSEQ(IE)
        DO 1800 I=1,N
!	Modified midpoint method
          WK(I,3+IE)=0.5*(WK(I,1)+WK(I,2)+H2*WK(I,3))
          WK(I,15+IE)=WK(I,3+IE)
          WK(I,27+IE)=WK(I,3+IE)
1800    CONTINUE

        IF(IE.GT.1) THEN
!	Perform extrapolation
          DO 2000 J=IE-1,1,-1
            DO 2000 I=1,N
              IF(IFLG.EQ.0) THEN
!	Use polynomial extrapolation
                CN=(WK(I,4+J)-WK(I,15+J))/(HI(J)-HI(IE))
!	D(J,IE-J)
                WK(I,15+J)=HI(IE)*CN
!	C(J,IE-J)
                WK(I,3+J)=HI(J)*CN
              ELSE

!	Use rational function extrapolation
                DEN=HI(J)*WK(I,15+J)/HI(IE)-WK(I,4+J)
                IF(DEN.EQ.0.0) THEN
                  IER=733
                  GO TO 2800
                ENDIF
                CN=(WK(I,4+J)-WK(I,15+J))/DEN
!	C(J,IE-J)
                WK(I,3+J)=HI(J)*WK(I,15+J)*CN/HI(IE)
!	D(J,IE-J)
                WK(I,15+J)=WK(I,4+J)*CN
              ENDIF
              WK(I,27+J)=WK(I,3+J)+WK(I,27+J)
2000      CONTINUE

!	Estimating the truncation error
          ERR=0
          DO 2200 I=1,N
            R2=ABS(Y(I))+ABS(WK(I,28)-Y(I))+EPS
            R1=ABS(WK(I,28)-WK(I,29))/R2
            IF(R1.GT.ERR) ERR=R1
2200      CONTINUE
          ERR=ABS(ERR*TSTEP/H)
          IF(ERR.LT.REPS) GO TO 2800
        ENDIF
2500  CONTINUE

2800  IF(T1.EQ.T0) THEN
!	The step size is too small for the arithmetic used
        IF(IER.EQ.0) IER=731
        RETURN
      ENDIF

      IF(ERR.LT.REPS) THEN
!	Integration is successful
        T0=T0+H
        DO 3000 I=1,N
3000    Y(I)=WK(I,28)
        IF(ABS((TN-T0)/TSTEP).LT.REPS) RETURN

        IF(IE.LT.7) THEN
!	Increase the step size
          H=SFAC*NSEQ(7)*H/NSEQ(IE)
        ELSE IF(IE.GT.7) THEN
!	Decrease the step size
          H=SFAC*H
        ENDIF
        IF(T0+H.GT.TN.EQV.H.GT.0) H=TN-T0
      ELSE
!	If the integration has failed, then decrease the step size
        H=H/32.
      ENDIF

      IF(ABS(H/TSTEP).LT.REPS) THEN
!	The step size is too small
        IF(IER.EQ.0) IER=730
        RETURN
      ENDIF
      IF(NSTEP.LT.NMAX) GO TO 1000
      IER=732
      END
