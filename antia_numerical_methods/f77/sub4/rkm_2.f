!	To solve initial value problems in ordinary differential equations
!	using a second or fourth-order Runge-Kutta method with adaptive
!	step size control
!
!	N : (input) Number of first order differential equations to be solved
!	Y : (input/output) Real array of length N containing the initial
!		values of variables. After execution it will contain the
!		values at the last point where the integration
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
!		to be updated.
!	REPS : (input) Required accuracy in each component of the solution.
!		The subroutine only controls local truncation error and hence
!		actual error could be larger
!	NSTEP : (output) Number of steps required to complete the integration
!		Each step requires 4 or 5 calls to DIF (using RK2)
!	NMAX : (input/output) Maximum number of steps to be used. If NMAX.LE.0
!		it will be set to a default value of NMX=10000.
!	IER : (output) Error parameter; IER=0 implies successful execution
!		IER=701 implies N.LE.0, no calculations are done
!		IER=721 implies that step-size has become smaller than
!			REPS*|TN-T0|
!		IER=722 implies that step size is too small for arithmetic used
!		IER=723 implies that integration could not be completed in
!			the specified number of steps.
!	WK : Real array of length 5N used as scratch space
!	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
!	the differential equation. T is the value of independent variable,
!	N is the number of variables, Y is a real array of length N containing
!	the values of variables. DY is a real array of length N which should
!	contain the calculated values of derivatives at (T,Y).
!
!	Required routines : RK2, DIF
!	
!	
      SUBROUTINE RKM_2(N,Y,DY,DIF,H,T0,TN,REPS,NSTEP,NMAX,IER,WK)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=10000,SFAC=0.9D0,EPS=1.D-30)
!      PARAMETER(E1=0.2D0,E2=0.25D0,CE=15.D0)
!	For use with RK4 use the preceding statement instead of the following one
      PARAMETER(E1=0.33D0,E2=0.5D0,CE=3.D0)
      EXTERNAL DIF
      DIMENSION Y(N),DY(N),WK(N,5)

      IF(N.LE.0) THEN
        IER=701
        RETURN
      ENDIF
      IF(NMAX.LE.0) NMAX=NMX
      NSTEP=0
      TSTEP=TN-T0
      IER=0
      IF(TSTEP.EQ.0.0) RETURN
!	Adjust the initial value of H if needed
      IF(H.EQ.0.0.OR.ABS(H).GT.ABS(TSTEP)) H=TSTEP
      IF(H.LT.0.0.EQV.TN.GT.T0) H=-H
!	Calculate the derivatives at the initial point
      CALL DIF(T0,N,Y,DY)

!	Loop for integration
1000  NSTEP=NSTEP+1
!	Use two steps of h/2
      H2=H/2.
      CALL RK2(N,T0,Y,DY,H2,WK(1,4),DIF,WK)
      T1=T0+H2
      CALL DIF(T1,N,WK(1,4),WK(1,5))
      CALL RK2(N,T1,WK(1,4),WK(1,5),H2,WK(1,3),DIF,WK)

!	Use single step of h
      CALL RK2(N,T0,Y,DY,H,WK(1,4),DIF,WK)

!	Estimate the truncation error
      ERR=0
      DO 2000 I=1,N
        R2=ABS(Y(I))+ABS(WK(I,3)-Y(I))+EPS
        R1=ABS(WK(I,4)-WK(I,3))/R2
        IF(R1.GT.ERR) ERR=R1
2000  CONTINUE
      ERR=ABS(ERR*TSTEP/(CE*H))

      IF(T1.EQ.T0) THEN
!	Step size is too small for arithmetic used
        IER=722
        RETURN
      ENDIF

      IF(ERR.LT.REPS) THEN
!	Integration at this step is successful, update T0 and Y
        T0=T0+H
        DO 3000 I=1,N
3000    Y(I)=WK(I,3)-(WK(I,4)-WK(I,3))/CE
        CALL DIF(T0,N,Y,DY)
        IF(ABS((TN-T0)/TSTEP).LT.REPS) RETURN

!	Adjust the step size
        IF(ERR.EQ.0.0) THEN
          H=2.*H
        ELSE
          H=SFAC*H*(REPS/ERR)**E1
        ENDIF
        IF(T0+H.GT.TN.EQV.H.GT.0) H=TN-T0
      ELSE

!	If the integration is not successful, try again with smaller step size
        H=SFAC*H*(REPS/ERR)**E2
      ENDIF

      IF(ABS(H/TSTEP).LT.REPS) THEN
!	Step size is too small
        IER=721
        RETURN
      ENDIF
      IF(NSTEP.LT.NMAX) GO TO 1000
      IER=723
      END
