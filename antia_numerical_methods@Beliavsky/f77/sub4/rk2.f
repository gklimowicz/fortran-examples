!	To perform one step of integration of ordinary differential equations
!	using a second-order Runge-Kutta method 
!
!	N : (input) Number of first order differential equations to be solved
!	T : (input) Initial value of independent variable t, at
!		which initial values are specified. This value is not updated.
!	Y0 : (input) Real array of length N containing the initial
!		values of variables.
!	DY0 : (input) Real array of length N containing the derivatives
!		of Y at the initial point Y0
!	H : (input) The step size to be used for integration
!	Y1 : (output) Real array of length N containing the solution at t=T+H
!	DIF : (input) Name of subroutine to calculate the right hand side
!		of differential equation y'=f(t,y)
!	WK : Real array of length 2N used as scratch space
!	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
!	the differential equation. T is the value of independent variable,
!	N is the number of variables, Y is a real array of length N containing
!	the values of variables. DY is a real array of length N which should
!	contain the calculated values of derivatives at (T,Y).
!
!	Required routines : DIF
!	
!	
      SUBROUTINE RK2(N,T,Y0,DY0,H,Y1,DIF,WK)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y0(N),DY0(N),Y1(N),WK(N,2)

      DO 1000 I=1,N
        Y1(I)=H*DY0(I)
1000  WK(I,1)=Y0(I)+Y1(I)

      T1=T+H
      CALL DIF(T1,N,WK,WK(1,2))
      DO 1200 I=1,N
1200  Y1(I)=Y0(I)+(Y1(I)+H*WK(I,2))/2.
      END