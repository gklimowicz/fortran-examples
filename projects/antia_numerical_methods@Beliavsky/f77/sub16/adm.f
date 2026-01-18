!	To solve linear parabolic differential equation in 2 dimensions
!	using the alternating direction method
!	The differential equation is assumed to be of the form
!
!	du/dt = Axx(x,y,t)d^2u/dx^2 + Ayy(x,y,t)d^2u/dy^2 + Ax(x,y,t)du/dx
!		+Ay(x,y,t)du/dy + Au(x,y,t)u +A0(x,y,t)
!
!	with Dirichlet boundary conditions
!
!	u(X0,y,t)=BC(1,X0,y,t);		u(XN,y,t)=BC(2,XN,y,t)
!	u(x,Y0,t)=BC(3,x,Y0,t);		u(x,YN,t)=BC(4,x,YN,t)
!
!	T : (input/output) Initial value of "time" where the initial conditions
!		are specified. After execution, it will be replaced by the
!		value of T at the last point where execution is successful
!	DT : (input) The time step to be used for computations. It is kept fixed.
!	X0 : (input) Lower limit on X where the solution is required
!	XN : (input) Upper limit on X where the solution is required.
!		Solution is computed in the interval (X0,XN)
!	Y0 : (input) Lower limit on Y where the solution is required
!	YN : (input) Upper limit on Y where the solution is required.
!		Solution is computed in the interval (Y0,YN)
!	NT : (input) Number of time steps each of length DT  to be executed.
!	NX : (input) Number of mesh points in the X direction.
!	NY : (input) Number of mesh points in the Y direction.
!	X : (output) Real array of length NX containing the mesh points used
!		in X direction. These are calculated by the routine by
!		assuming uniform spacing.
!	Y : (output) Real array of length NY containing the mesh points used
!		in Y direction. These are calculated by the routine by
!		assuming uniform spacing.
!	U : (input/output) Real array of length IU*NY containing the solution
!		at the current time step. It should contain the initial
!		values at the time of calling, if IFLG.NE.0. Otherwise it
!		is computed using Function FIC. After execution it will contain
!		the computed solution at t=T. U(I,J) is the solution at (x_i,y_j)
!	IU : (input) First dimension of array U as declared in the calling program
!	COF : (input) Name of the subroutine to calculate the coefficients
!		in the equation
!	BC : (input) Name of the function routine to calculate the boundary values
!	FIC : (input) Name of the function routine to calculate the initial
!		values when IFLG=0. For other values of IFLG this routine
!		is not used, but a dummy routine may be required by the compiler.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=714 implies that DT=0, XN=X0, YN=Y0, NX<3, NY<3 or IU<NX
!			in which case no calculations are done
!		IER=762 implies that the difference equations are singular
!			and solution cannot be continued further.
!	IFLG : (input/output) Integer variable used as a flag to denote how the
!			initial values are calculated.
!		If IFLG=0 then the initial values are calculated using the
!			function FIC to be supplied by the user. IFLG is set to 1
!		Otherwise the initial values must be supplied in array U.
!	WK : Real array of length MAX(NX,NY)*(4+NY) used as scratch space
!
!	SUBROUTINE COF(X,Y,T,AXX,AYY,AX,AY,AU,A0), FUNCTION BC(IB,X,Y,T)
!	and FUNCTION FIC(X,Y,T) must be supplied by the user 
!	Subroutine COF should calculate the coefficients AXX, AYY, AX, AY
!	AU, A0 defined above for given values of X, Y, T.
!	Function BC should calculate the value of solution at the boundaries
!	as described earlier. Integer IB specifies which boundary is required.
!	Function FIC is required only if IFLG=0, otherwise a dummy routine
!	with this name will suffice. If IFLG=0, function FIC must calculate
!	the initial values at required X,Y,T.
!
!	Required routines : COF, BC, FIC

      SUBROUTINE ADM(T,DT,X0,XN,Y0,YN,NT,NX,NY,X,Y,U,IU,COF,BC,FIC,IER,
     1               IFLG,WK)
      IMPLICIT REAL*16(A-H,O-Z)
!	Some Fortran compilers may not accept MAX(NX,NY) in dimension
!	In that case make it NX+NY, which will occupy more space.
      DIMENSION X(NX),Y(NY),U(IU,NY),WK(MAX(NX,NY),NY+4)

      IF(DT.EQ.0.OR.XN.EQ.X0.OR.YN.EQ.Y0.OR.NX.LE.2.OR.NY.LE.2.
     1   OR.IU.LT.NX) THEN
        IER=714
        RETURN
      ENDIF
      IER=762

!	Setting up the grid points along X and Y
      DX=(XN-X0)/(NX-1)
      DO 1000 I=1,NX
1000  X(I)=X0+(I-1)*DX
      DY=(YN-Y0)/(NY-1)
      DO 1200 I=1,NY
1200  Y(I)=Y0+(I-1)*DY

      IF(IFLG.EQ.0) THEN
!	Calculate the initial values
        DO 1400 J=1,NY
          DO 1400 I=1,NX
            U(I,J)=FIC(X(I),Y(J),T)
1400    CONTINUE
      ENDIF
      IFLG=1
      R=0.5*DT/DX**2
      R1=0.25*DT/DX
      S=0.5*DT/DY**2
      S1=0.25*DT/DY

!	Loop over time steps
      DO 6000 IT=1,NT
        T1=T+DT/2.
        T2=T+DT

!	Setting up the equations for the first half-step
        DO 3500 K=2,NY-1
          DO 2000 J=2,NX-1
            CALL COF(X(J),Y(K),T,AXX,AYY,AX,AY,AU,A0)
            CALL COF(X(J),Y(K),T1,BXX,BYY,BX,BY,BU,B0)
            WK(J,1)=-R*BXX+R1*BX
            WK(J,2)=1+2.*R*BXX-0.25*DT*BU
            WK(J,3)=-R*BXX-R1*BX
            WK(J,4+K)=U(J,K)+S*AYY*(U(J,K+1)-2*U(J,K)+U(J,K-1))+
     1               S1*AY*(U(J,K+1)-U(J,K-1))+0.25*DT*(AU*U(J,K)+A0+B0)
2000      CONTINUE

!	The boundary conditions
          U1=BC(1,X0,Y(K),T1)
          WK(1,4+K)=U1
          WK(2,4+K)=WK(2,4+K)-U1*WK(2,1)
          U1=BC(2,XN,Y(K),T1)
          WK(NX,4+K)=U1
          WK(NX-1,4+K)=WK(NX-1,4+K)-U1*WK(NX-1,3)

!	Gaussian elimination for the tridiagonal system
          DO 2500 J=3,NX-1
            IF(WK(J-1,2).EQ.0) RETURN
            RP=-WK(J,1)/WK(J-1,2)
            WK(J,2)=WK(J,2)+RP*WK(J-1,3)
            WK(J,4+K)=WK(J,4+K)+RP*WK(J-1,4+K)
2500      CONTINUE
          IF(WK(NX-1,2).EQ.0) RETURN

!	Back-substitution
          WK(NX-1,4+K)=WK(NX-1,4+K)/WK(NX-1,2)
          DO 3000 J=NX-2,2,-1
            WK(J,4+K)=(WK(J,4+K)-WK(J,3)*WK(J+1,4+K))/WK(J,2)
3000      CONTINUE
3500    CONTINUE

!	Setting up the equations for the second half step
        DO 4800 K=2,NX-1
          DO 4000 J=2,NY-1
            CALL COF(X(K),Y(J),T1,AXX,AYY,AX,AY,AU,A0)
            CALL COF(X(K),Y(J),T2,BXX,BYY,BX,BY,BU,B0)
            WK(J,1)=-S*BYY+S1*BY
            WK(J,2)=1+2.*S*BYY-0.25*DT*BU
            WK(J,3)=-S*BYY-S1*BY
            WK4=WK(K,J+4)+R*AXX*(WK(K+1,J+4)-2*WK(K,J+4)+WK(K-1,J+4))+
     1      R1*AX*(WK(K+1,J+4)-WK(K-1,J+4))+0.25*DT*(AU*WK(K,J+4)+A0+B0)
            IF(K.GT.2) WK(K-1,J+4)=WK(J,4)
            WK(J,4)=WK4
4000      CONTINUE

!	The boundary conditions
          WK(K,5)=BC(3,X(K),Y0,T2)
          WK(2,4)=WK(2,4)-WK(K,5)*WK(2,1)
          WK(K,NY+4)=BC(4,X(K),YN,T2)
          WK(NY-1,4)=WK(NY-1,4)-WK(K,NY+4)*WK(NY-1,3)

!	Gaussian elimination for the tridiagonal system
          DO 4200 J=3,NY-1
            IF(WK(J-1,2).EQ.0) RETURN
            RP=-WK(J,1)/WK(J-1,2)
            WK(J,2)=WK(J,2)+RP*WK(J-1,3)
            WK(J,4)=WK(J,4)+RP*WK(J-1,4)
4200      CONTINUE
          IF(WK(NY-1,2).EQ.0) RETURN

!	Back-substitution
          WK(NY-1,4)=WK(NY-1,4)/WK(NY-1,2)
          DO 4400 J=NY-2,2,-1
            WK(J,4)=(WK(J,4)-WK(J,3)*WK(J+1,4))/WK(J,2)
4400      CONTINUE
4800    CONTINUE

        DO 5000 J=2,NY-1
5000    WK(NX-1,J+4)=WK(J,4)
        DO 5500 J=1,NY
          U(1,J)=BC(1,X0,Y(J),T2)
          U(NX,J)=BC(2,XN,Y(J),T2)
          DO 5500 K=2,NX-1
5500    U(K,J)=WK(K,J+4)

        T=T2
6000  CONTINUE
      IER=0
      END
