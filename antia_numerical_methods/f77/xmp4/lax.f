!     PROGRAM TO SOLVE NONLINEAR HYPERBOLIC EQUATIONS WRITTEN IN
!     THE CONSERVATION LAW FORM USING LAX-WENDROFF METHOD

      PROGRAM WAVE
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(21),U(5,21),WK(150),IWK(21)
      EXTERNAL FLUX,BC,FIC

!     EXAMPLE 14.3 : THE WAVE EQUATION

51    FORMAT('   IER =',I4,5X,'NX =',I4,5X,'TIME =',1PE14.6,5X,
     1       'TIME STEP =',E14.6/5X,'SOLUTION =',4E14.6/(2X,5E14.6))
52    FORMAT(5X,1P5E14.6)

      N=2
      IU=5
      X0=0
      XN=0.5
      IFLG=0
      PRINT *,'TYPE T=INITIAL TIME,  NX=NO. OF PTS'
      READ *,T,NX

100   PRINT *,'TYPE DT=TIME STEP,  NT=NO. OF TIME STEPS'
     1       ,'       (QUITS WHEN NT.LE.0)'
      READ *,DT,NT
      IF(NT.LE.0) STOP
      CALL LAX(N,T,DT,X0,XN,NT,NX,X,U,IU,FLUX,BC,FIC,IER,IFLG,WK,IWK)
      WRITE(6,51) IER,NX,T,DT,(U(1,I),I=1,NX)
!     WRITE THE SECOND COMPONENT OF THE SOLUTION IN NEXT LINE
      WRITE(6,52) (U(2,I),I=1,NX)
      GO TO 100
      END

!     ----------------------------------------------

!	To solve a system of hyperbolic differential equation using the
!	Lax-Wendroff difference scheme
!	The differential equations are assumed to be of the form
!
!	du_i/dt + df_i/dx =0,   where f_i=f_i(x,t,u)
!
!	while boundary conditions are either Dirichlet form or
!
!	A_i1(t) f_i+A_i2(t) df_i/dx = A_i3(t)
!
!	N : (input) Number of equations in the system
!	T : (input/output) Initial value of "time" where the initial conditions
!		are specified. After execution, it will be replaced by the
!		value of T at the last point where execution is successful
!	DT : (input) The time step to be used for computations. It is kept fixed.
!	X0 : (input) Lower limit on X where the solution is required
!	XN : (input) Upper limit on X where the solution is required.
!		Solution is computed in the interval (X0,XN)
!	NT : (input) Number of time steps each of length DT  to be executed.
!	NX : (input) Number of mesh points in the X direction.
!	X : (output) Real array of length NX containing the mesh points used
!		in X direction. These are calculated by the routine by
!		assuming uniform spacing.
!	U : (input/output) Real array of length IU*NX containing the solution
!		at the current time step. It should contain the initial
!		values at the time of calling, if IFLG.NE.0. Otherwise it
!		is computed using subroutine FIC. After execution it will contain
!		the computed solution at t=T. U(I,J) is the Ith component
!		at Jth mesh point in X.
!	IU : (input) The first dimension of U as declared in the calling
!		program, IU.GE.N
!	FLUX : (input) Name of the subroutine to calculate the "fluxes"
!		f_i(x,t,u)
!	BC : (input) Name of the subroutine to calculate the coefficients
!		in the boundary conditions
!	FIC : (input) Name of the subroutine to calculate the initial
!		values when IFLG=0. For other values of IFLG this routine
!		is not used, but a dummy routine may be required by the compiler.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=715 implies that DT=0, XN=X0, NX<3, or IU<N
!			in which case no calculations are done
!		IER=763 implies that the difference equations are singular
!			and solution cannot be continued further.
!	IFLG : (input/output) Integer variable used as a flag to denote how the
!			initial values are calculated.
!		If IFLG=0 then the initial values are calculated using the
!			subroutine FIC to be supplied by the user. IFLG is set to 1
!		Otherwise the initial values must be supplied in array U.
!	WK : Real array of length N*(NX+5) used as scratch space
!	IWK : Integer array of length N used as scratch space
!
!	SUBROUTINE FLUX(N,X,T,U,F), SUBROUTINE BC(IB,N,X,T,IW,A)
!	and SUBROUTINE FIC(N,X,T,U) must be supplied by the user 
!	Subroutine FLUX should calculate f_i(x,t,u) as
!	defined above for given values of X,T,U. Here X,T are the values
!	of x and t where the calculations are required, N is the number of
!	equations, while U and F are real arrays 
!	of length N containing the solution U and corresponding fluxes F.
!	The subroutine should calculate F(I) using the values of X,T,U.
!
!	Subroutine BC should calculate the coefficients for the boundary
!	conditions for given values of X, T. Here IB is an integer variable
!	used as a flag to denote the boundary at which the boundary conditions
!	are required. IB=1 for boundary at x=X0 and IB=2 for x=XN.
!	N is the number of equations in the system, X and T
!	are the values of x and t where the boundary conditions are required.
!	IW is an integer array which is used as a flag to transmit the
!	type of boundary conditions. This flag also must be set by subroutine
!	BC. IW(I) should be set to zero if the boundary conditions are of
!	Dirichlet form, any other value will imply boundary condition on
!	flux as described above. A is a real array of length N*3 which
!	contains the coefficients required to define the boundary conditions.
!	This array must be dimensioned A(N,3) in the subroutine. For
!	Dirichlet conditions the boundary values are given by
!	A(I,1)U(I)=A(I,3), otherwise the boundary condition is assumed as
!	A(I,1)F(I)+A(I,2)dF(I)/dx=A(I,3)
!	The coefficient A(I,J) and IW(I) for each I, must be set by
!	subroutine BC.
!
!	Subroutine FIC is required only if IFLG=0, otherwise a dummy routine
!	with this name will suffice. If IFLG=0, subroutine FIC must calculate
!	the initial values at required X,T. Here N is the number of equations
!	and U is a real array of length N, which should return the calculated
!	initial values of all components.
!
!	Required routines : FLUX, BC, FIC

      SUBROUTINE LAX(N,T,DT,X0,XN,NT,NX,X,U,IU,FLUX,BC,FIC,IER,IFLG,
     1               WK,IWK)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(IU,NX),X(NX),WK(N,NX+5),IWK(N)

      IF(DT.EQ.0.OR.XN.EQ.X0.OR.NX.LE.2.OR.IU.LT.N) THEN
        IER=715
        RETURN
      ENDIF
      IER=763

      DX=(XN-X0)/(NX-1)
      S=DT/DX
!	Setting the initial values
      DO 1000 J=1,NX
        X(J)=X0+(J-1)*DX
        IF(IFLG.EQ.0) CALL FIC(N,X(J),T,U(1,J))
1000  CONTINUE
      IFLG=1

!	Loop over time steps
      DO 5000 I=1,NT
        T1=T+DT/2.
        TT=T+DT

!	First half-step of Lax-Wendroff scheme
        J1=NX+1
        J2=NX+2
        CALL FLUX(N,X(1),T,U(1,1),WK(1,J1))
        DO 3000 J=2,NX
          CALL FLUX(N,X(J),T,U(1,J),WK(1,J2))
          DO 2500 K=1,N
2500      WK(K,J)=0.5*(U(K,J)+U(K,J-1))-0.5*S*(WK(K,J2)-WK(K,J1))
          JT=J1
          J1=J2
          J2=JT
3000    CONTINUE

!	The boundary conditions at X0
        J1=NX+1
        J2=NX+2
        CALL FLUX(N,X(2),T1,WK(1,2),WK(1,J1))
        CALL BC(1,N,X0,TT,IWK,WK(1,NX+3))
        DO 3200 K=1,N
          IF(IWK(K).EQ.0) THEN
            IF(WK(K,NX+3).EQ.0.0) RETURN
            U(K,1)=WK(K,NX+5)/WK(K,NX+3)
          ELSE
            D1=0.5*WK(K,NX+3)-WK(K,NX+4)/DX
            IF(D1.EQ.0.0) RETURN
            F1=(WK(K,NX+5)-WK(K,J1)*(0.5*WK(K,NX+3)+WK(K,NX+4)/DX))/D1
            U(K,1)=U(K,1)-S*(WK(K,J1)-F1)
          ENDIF
3200    CONTINUE

!	Second half-step of Lax-Wendroff scheme
        DO 4000 J=2,NX-1
          CALL FLUX(N,X(J+1),T1,WK(1,J+1),WK(1,J2))
          DO 3500 K=1,N
3500      U(K,J)=U(K,J)-S*(WK(K,J2)-WK(K,J1))
          JT=J1
          J1=J2
          J2=JT
4000    CONTINUE

!	Boundary conditions at XN
        CALL BC(2,N,XN,TT,IWK,WK(1,NX+3))
        DO 4200 K=1,N
          IF(IWK(K).EQ.0) THEN
            IF(WK(K,NX+3).EQ.0.0) RETURN
            U(K,NX)=WK(K,NX+5)/WK(K,NX+3)
          ELSE
            D1=0.5*WK(K,NX+3)+WK(K,NX+4)/DX
            IF(D1.EQ.0.0) RETURN
            F1=(WK(K,NX+5)-WK(K,J1)*(0.5*WK(K,NX+3)-WK(K,NX+4)/DX))/D1
            U(K,NX)=U(K,NX)-S*(F1-WK(K,J1))
          ENDIF
4200    CONTINUE

        T=TT
5000  CONTINUE
      IER=0
      END

!     --------------------------------------------------

      SUBROUTINE FLUX(N,X,T,U,F)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(N),F(N)

!     THE WAVE EQ. : F1=-U2,  F2=-U1

      F(2)=-U(1)
      F(1)=-U(2)
      END

!     ---------------------------------------------

      SUBROUTINE BC(IB,N,X,T,IC,A)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.1415926535D0)
      DIMENSION IC(N),A(N,3)

!     BOUNDARY CONDITIONS F1'(0,T)=0, F2(0,T)=0, F1(1/2,T)=0, F2'(1/2,T)=0

      DO 1000 I=1,3
      DO 1000 J=1,N
1000  A(J,I)=0.0
      IF(IB.EQ.1) THEN
        A(1,2)=1.0
        A(2,1)=1.0
      ELSE
        A(1,1)=1.0
        A(2,2)=1.0
      ENDIF
      IC(1)=1
      IC(2)=1
      END

!     ---------------------------------------

      SUBROUTINE FIC(N,X,T,U)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.1415926535D0)
      DIMENSION U(N)

!     THE EXACT SOLUTION, USED TO GENERATE STARTING VALUES

      U(1)=-SIN(PI*X)*SIN(PI*T)
      U(2)=COS(PI*X)*COS(PI*T)
      END
