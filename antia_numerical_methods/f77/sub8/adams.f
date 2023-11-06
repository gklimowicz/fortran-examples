!	To perform one step of solution of initial value problems in ordinary
!	differential equations using fourth order Adams-Bashforth-Moulton
!	predictor corrector method. It is called by subroutine MSTEP.
!
!	N : (input) Number of first order differential equations to be solved
!	Y : (input/output) Real array of length 7N containing the solution
!		at last four points. After execution new point will be added.
!	DY : (input/output) Real array of length 7N containing the derivatives
!		of Y at the points in Y
!	DIF : (input) Name of subroutine to calculate the right hand side
!		of differential equation y'=f(t,y)
!	H : (input) Step length to be used.
!	T : (input) Value of independent variable t, at which the solution
!		is required.
!	REPS : (input) Required accuracy in each component of the solution.
!		The subroutine only controls error in iteration on corrector.
!	NSTP : (input/output) Number of calls to DIF required so far.
!		The count is updated by the subroutine.
!	IJ : (input) The index j+1 in corrector formula
!	IJM1 : (input) The index j in corrector formula
!	IJM2 : (input) The index j-1 in corrector formula
!	IJM3 : (input) The index j-2 in corrector formula
!	IJM4 : (input) The index j-3 in corrector formula
!	IER : (output) Error parameter; IER=0 implies successful execution
!		IER=729 implies that iteration on corrector failed to converge
!	WK : Real array of length 2N used to pass the predicted values
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
      SUBROUTINE ADAMS(N,Y,DY,DIF,H,T,REPS,NSTP,IJ,IJM1,IJM2,
     1                 IJM3,IJM4,IER,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(EPS=1.D-30,NIT=10,CFAC=0.2D0)
      DIMENSION Y(N,7),DY(N,7),WK(N,2)

      IER=0
      DO 2000 J=1,N
!	The predictor
        T1=Y(J,IJM1)+H*(55.*DY(J,IJM1)-59.*DY(J,IJM2)+37.*DY(J,IJM3)
     1     -9.*DY(J,IJM4))/24.
!	Modifier
        Y(J,IJ)=T1+251.*(Y(J,IJM1)-WK(J,2))/270.
2000  WK(J,1)=T1

!	Iteration on corrector
      DO 3000 J=1,NIT
        NSTP=NSTP+1
        CALL DIF(T,N,Y(1,IJ),DY(1,IJ))
        ERR=0.0
        DO 2600 K=1,N
!	The corrector
          T2=Y(K,IJM1)+H*(9.*DY(K,IJ)+19.*DY(K,IJM1)-5.*DY(K,IJM2)
     1       +DY(K,IJM3))/24.
          R1=ABS(T2)+ABS(T2-Y(K,IJM1))+EPS
          T1=ABS((Y(K,IJ)-T2)/R1)
          IF(T1.GT.ERR) ERR=T1
2600    Y(K,IJ)=T2
!	The convergence test
        IF(ERR.LT.CFAC*REPS) RETURN
3000  CONTINUE

!	iteration on corrector fails to converge
      IER=729

      END
