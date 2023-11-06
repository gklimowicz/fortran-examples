!	To generate the starting values for fourth order multistep methods
!	in solution of initial value problem in ordinary differential equations.
!	It is used by MSTEP.
!
!	N : (input) Number of first order differential equations to be solved
!	Y : (input/output) Real array of length 4N, the first N elements
!		should contain the initial values of variables. During
!		execution, the solution at the next three points will
!		be stored in this array. 
!	DY : (output) Real array of length 4N, containing the derivatives
!		of Y the points stored in Y.
!	DIF : (input) Name of subroutine to calculate the right hand side
!		of differential equation y'=f(t,y)
!	H : (input/output) Initial guess for the step size. After execution
!		it will contain the step size used by the program
!	T : (input/output) Initial value of independent variable t, at
!		which initial values are specified. After execution it will
!		be set to T+3H if execution is successful.
!	REPS : (input) Required accuracy in each component of the solution.
!		The subroutine only controls local truncation error and hence
!		actual error could be larger
!	IFLG : (input) Integer variable used as flag to decide the
!		type of computation.
!		If IFLG=0 the step size is adjusted to meet the specified accuracy.
!		If IFLG=1 the step size is kept fixed.
!	TSTEP : (input) The size of interval over which integration is requested
!			It is used only for convergence check.
!	NSTP : (input/output) Number of calls to DIF required since the starting
!		of integration. The count is updated by the subroutine.
!	IER : (output) Error parameter; IER=0 implies successful execution
!		IER=724 implies that step size becomes too small
!		IER=725 implies that iteration to generate starting values
!			failed to converge
!	WK : Real array of length 2N used to pass on the value for modifier
!
!	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
!	the differential equation. T is the value of independent variable,
!	N is the number of variables, Y is a real array of length N containing
!	the values of variables. DY is a real array of length N which should
!	contain the calculated values of derivatives at (T,Y).
!
!	Required routines : RK4, DIF
!	
!	
      SUBROUTINE STRT4(N,Y,DY,DIF,H,T,REPS,IFLG,TSTEP,NSTP,IER,WK)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(EPS=1.Q-30,SFAC=0.9Q0,NIT=10)
      EXTERNAL DIF
      DIMENSION Y(N,4),DY(N,4),WK(N,2)

      T0=T
      IER=0
      DO 4000 IT=1,NIT
        CALL DIF(T0,N,Y(1,1),DY(1,1))
!	Generate y_1
        CALL RK4(N,T0,Y(1,1),DY(1,1),H,Y(1,2),DIF,WK)
        T=T0+H
        CALL DIF(T,N,Y(1,2),DY(1,2))
!	Generate y_2
        CALL RK4(N,T,Y(1,2),DY(1,2),H,Y(1,3),DIF,WK)
        H2=2.*H
!	Calculate y_2 using double step
        CALL RK4(N,T0,Y(1,1),DY(1,1),H2,Y(1,4),DIF,WK)
        NSTP=NSTP+11

!	Estimate the truncation error in y_2
        ERR=0
        DO 2000 I=1,N
          R2=ABS(Y(I,1))+ABS(Y(I,2))+ABS(Y(I,3))+EPS
          R1=ABS(Y(I,3)-Y(I,4))/R2
          IF(R1.GT.ERR) ERR=R1
2000    CONTINUE
        ERR=ERR*TSTEP/H

        IF(ERR.LE.REPS.OR.IFLG.GT.0) THEN
!	Accept the computed values
          T=T+H
          CALL DIF(T,N,Y(1,3),DY(1,3))
!	Generate y_3
          CALL RK4(N,T,Y(1,3),DY(1,3),H,Y(1,4),DIF,WK)
          T=T+H
          CALL DIF(T,N,Y(1,4),DY(1,4))
          NSTP=NSTP+5
!	Store Y_3 for use by modifier
          DO 3000 I=1,N
3000      WK(I,2)=Y(I,4)
          RETURN
        ELSE

!	Reduce the step size
          H=SFAC*H*(REPS/ERR)**0.25Q0
          IF(ABS(H/TSTEP).LT.REPS.OR.T.EQ.T0) THEN
!	Step size is too small
            IER=724
            RETURN
          ENDIF
        ENDIF
4000  CONTINUE

!	Iteration to generate starting values fails to converge
      IER=725
      END
