!	To solve linear second order elliptic differential equation using the
!	Alternating direction implicit iterative (ADI) method
!	The differential equation is assumed to be of the form
!
!	Axx(x,y)d^2u/dx^2 + Ayy(x,y)d^2u/dy^2 + Ax(x,y)du/dx
!		+ Ay(x,y)du/dy + A0(x,y)u + F(x,y)=0
!
!	with following boundary conditions on a rectangular region
!
!		A0*u+An*dun=F; where dun is the normal derivative of u
!
!	X0 : (input) Lower limit on X where the solution is required
!	XN : (input) Upper limit on X where the solution is required.
!		Solution is computed in the interval (X0,XN)
!	Y0 : (input) Lower limit on Y where the solution is required
!	YN : (input) Upper limit on Y where the solution is required.
!		Solution is computed in the interval (Y0,YN)
!	KN : (input/output) Parameter k in ADI iteration. The subroutine repeats
!		a cycle of 2^k iteration, but convergence is checked after
!		each iteration. This parameter may be adjusted if it is
!		outside acceptable limits.
!	NX : (input) Number of mesh points in the X direction.
!	NY : (input) Number of mesh points in the Y direction.
!	X : (output) Real array of length NX containing the mesh points used
!		in X direction. These are calculated by the routine by
!		assuming uniform spacing.
!	Y : (output) Real array of length NY containing the mesh points used
!		in Y direction. These are calculated by the routine by
!		assuming uniform spacing.
!	U : (input/output) Real array of length IU*NY containing the solution
!		It should contain the initial values at the time of calling.
!		After execution it will contain	the computed solution.
!		U(I,J) is the solution at (x_i,y_j)
!	IU : (input) The first dimension of U as declared in the calling
!		program, IU.GE.NX
!	COF : (input) Name of the subroutine to calculate the coefficients
!		in the equation
!	BC : (input) Name of the subroutine to calculate the boundary conditions
!	EL : (input/output) Lower limit on the eigenvalues of the partitions
!		S and Y of the finite difference matrix. If EL.LE.0, then
!		the value for Poisson's equation will be used.
!	EU : (input/output) Upper limit on the eigenvalues of the partitions
!		S and Y of the finite difference matrix. If EU.LE.EL, then
!		the value for Poisson's equation will be used.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=717 implies that YN=Y0, XN=X0, NX<3, NY<3 or IU<NX
!			in which case no calculations are done
!		IER=766 implies that the matrix for ADI iteration is singular
!			and calculations have to be abandoned.
!		IER=767 implies that ADI iteration failed to converge to
!			specified accuracy
!	AEPS : (input) Required absolute accuracy. The ADI iteration is
!			continued until the change in all elements is less
!			than AEPS
!	NIT : (output) Number of ADI iterations required by the subroutine.
!	WK : Real array of length MAX(NX,NY)*(2+NY) used as scratch space
!	WKA : Real array of length 7*NX*NY used as scratch space
!
!	SUBROUTINE COF(X,Y,AXX,AYY,AX,AY,A0,F) and SUBROUTINE BC(IB,X,Y,A0,AN,F)
!	must be supplied by the user 
!	Subroutine COF should calculate the coefficients AXX, AYY, 
!	AX, AY, A0, F as defined above for given values of X,Y.
!	Subroutine BC should calculate the Boundary values at each boundary.
!	Here IB is an integer denoting which boundary is being considered.
!	IB=1 implies x=X0, IB=2 implies x=XN, IB=3 implies y=Y0, IB=4 implies y=YN
!
!	Required routines : COF, BC

      SUBROUTINE ADI(X0,XN,Y0,YN,KN,NX,NY,X,Y,U,IU,COF,BC,EL,EU,IER,
     1               AEPS,NIT,WK,WKA)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI2=1.570796326794896619D0,MAXIT=1000,KM=6,NK=2**(KM-1))
!	Some Fortran compilers may not accept MAX(NX,NY) in dimension statement.
!	In that case replace explicitly by maximum of NX, NY or by NX+NY
!	which will waste some memory
      DIMENSION X(NX),Y(NY),U(IU,NY),WK(MAX(NX,NY),*),ALP(KM),BET(KM)
     1          ,SA(NK,KM),WKA(7,NX,NY)

      IF(XN.EQ.X0.OR.YN.EQ.Y0.OR.NX.LE.2.OR.NY.LE.2.OR.IU.LT.NX) THEN
        IER=717
        RETURN
      ENDIF

!	setting up the mesh points
      DX=(XN-X0)/(NX-1)
      DO 1000 I=1,NX
1000  X(I)=X0+(I-1)*DX
      DY=(YN-Y0)/(NY-1)
      DO 1200 I=1,NY
1200  Y(I)=Y0+(I-1)*DY

      DX2=DX*DX
      DY2=DY*DY
      IF(EL.LE.0.0) THEN
!	Use the value for Poisson's equation
        EL1=4.*SIN(PI2/(NX-1.))**2/DX2
        EL2=4.*SIN(PI2/(NY-1.))**2/DY2
        EL=MIN(EL1,EL2)
      ENDIF

      IF(EU.LE.EL) THEN
!	Use the value for Poisson's equation
        EU1=4.*COS(PI2/(NX-1.))**2/DX2
        EU2=4.*COS(PI2/(NY-1.))**2/DY2
        EU=MAX(EU1,EU2)
      ENDIF
      IF(KN.GT.KM-1) KN=KM-1
      IF(KN.LE.2) KN=2

!	Calculating the parameters r_i for ADI iteration
      ALP(1)=EL
      BET(1)=EU
      DO 1400 J=1,KN
        ALP(J+1)=SQRT(ALP(J)*BET(J))
1400  BET(J+1)=0.5*(ALP(J)+BET(J))
      SA(1,1)=SQRT(ALP(KN+1)*BET(KN+1))
      J1=1
      DO 1600 J=1,KN
        PR=ALP(KN-J+1)*BET(KN-J+1)
        DO 1500 L=1,J1
          DISC=SA(L,J)**2-PR
          IF(DISC.LT.0.0) DISC=0.0
          SA(2*L-1,J+1)=SA(L,J)+SQRT(DISC)
          SA(2*L,J+1)=PR/SA(2*L-1,J+1)
1500    CONTINUE
1600  J1=J1*2
      IER=766

!	Setting up the finite difference equations
      DO 2000 I=1,NX
        DO 2000 J=1,NY
          CALL COF(X(I),Y(J),AXX,AYY,AX,AY,A0,F)
          WKA(1,I,J)=-AXX/DX2+0.5*AX/DX
          WKA(2,I,J)=2.*AXX/DX2-0.5*A0
          WKA(3,I,J)=-AXX/DX2-0.5*AX/DX
          WKA(4,I,J)=-AYY/DY2+0.5*AY/DY
          WKA(5,I,J)=2.*AYY/DY2-0.5*A0
          WKA(6,I,J)=-AYY/DY2-0.5*AY/DY
          WKA(7,I,J)=F
2000  CONTINUE

!	The boundary conditions at x=X0
      JL=2
      CALL BC(1,X0,Y0,A0,AX,F)
      IF(AX.NE.0.0) JL=1
      DO 2200 K=1,NY
        CALL BC(1,X0,Y(K),A0,AX,F)
        IF(JL.EQ.1) THEN
          IF(AX.EQ.0.0) RETURN
          WKA(2,1,K)=WKA(2,1,K)+WKA(1,1,K)*2.*DX*A0/AX
          WKA(3,1,K)=WKA(3,1,K)+WKA(1,1,K)
          WKA(7,1,K)=WKA(7,1,K)+WKA(1,1,K)*2.*DX*F/AX
        ELSE
          IF(A0.EQ.0.0) RETURN
          U(1,K)=F/A0
          WKA(7,2,K)=WKA(7,2,K)-WKA(1,2,K)*U(1,K)
        ENDIF
2200  CONTINUE

!	The boundary conditions at x=XN
      JU=NX-1
      CALL BC(2,XN,Y0,A0,AX,F)
      IF(AX.NE.0.0) JU=NX
      DO 2400 K=1,NY
        CALL BC(2,XN,Y(K),A0,AX,F)
        IF(JU.EQ.NX) THEN
          IF(AX.EQ.0.0) RETURN
          WKA(2,NX,K)=WKA(2,NX,K)-WKA(3,NX,K)*2.*DX*A0/AX
          WKA(1,NX,K)=WKA(1,NX,K)+WKA(3,NX,K)
          WKA(7,NX,K)=WKA(7,NX,K)-WKA(3,NX,K)*2.*DX*F/AX
        ELSE
          IF(A0.EQ.0.0) RETURN
          U(NX,K)=F/A0
          WKA(7,NX-1,K)=WKA(7,NX-1,K)-WKA(3,NX-1,K)*U(NX,K)
        ENDIF
2400  CONTINUE

!	The boundary conditions at y=Y0
      KL=2
      CALL BC(3,X0,Y0,A0,AX,F)
      IF(AX.NE.0.0) KL=1
      DO 2600 J=JL,JU
        CALL BC(3,X(J),Y0,A0,AX,F)
        IF(KL.EQ.1) THEN
          IF(AX.EQ.0.0) RETURN
          WKA(5,J,1)=WKA(5,J,1)+WKA(4,J,1)*2.*DY*A0/AX
          WKA(6,J,1)=WKA(6,J,1)+WKA(4,J,1)
          WKA(7,J,1)=WKA(7,J,1)+WKA(4,J,1)*2.*DY*F/AX
        ELSE
          IF(A0.EQ.0.0) RETURN
          U(J,1)=F/A0
          WKA(7,J,2)=WKA(7,J,2)-WKA(4,J,2)*U(J,1)
        ENDIF
2600  CONTINUE

!	The boundary conditions at y=YN
      KU=NY-1
      CALL BC(4,X0,YN,A0,AX,F)
      IF(AX.NE.0.0) KU=NY
      DO 2800 J=JL,JU
        CALL BC(4,X(J),YN,A0,AX,F)
        IF(KU.EQ.NY) THEN
          IF(AX.EQ.0.0) RETURN
          WKA(5,J,NY)=WKA(5,J,NY)-WKA(6,J,NY)*2.*DY*A0/AX
          WKA(4,J,NY)=WKA(4,J,NY)+WKA(6,J,NY)
          WKA(7,J,NY)=WKA(7,J,NY)-WKA(6,J,NY)*2.*DY*F/AX
        ELSE
          IF(A0.EQ.0.0) RETURN
          U(J,NY)=F/A0
          WKA(7,J,NY-1)=WKA(7,J,NY-1)-WKA(6,J,NY-1)*U(J,NY)
        ENDIF
2800  CONTINUE

!	Loop for ADI iteration
      N0=2**KN
      DO 6000 IT=1,MAXIT
        R=SA(MOD(IT-1,N0)+1,KN+1)

!	Gaussian elimination for the first half step
        DO 3500 K=KL,KU
          WK(JL,1)=WKA(2,JL,K)+R
          WK(JL,2+K)=WKA(7,JL,K)-(WKA(5,JL,K)-R)*U(JL,K)
          IF(K.GT.KL) WK(JL,2+K)=WK(JL,2+K)-WKA(4,JL,K)*U(JL,K-1)
          IF(K.LT.KU) WK(JL,2+K)=WK(JL,2+K)-WKA(6,JL,K)*U(JL,K+1)
          DO 3000 J=JL+1,JU
            IF(WK(J-1,1).EQ.0.0) RETURN
            RP=-WKA(1,J,K)/WK(J-1,1)
            WK(J,1)=WKA(2,J,K)+R+RP*WKA(3,J-1,K)
            WK(J,2+K)=WKA(7,J,K)-(WKA(5,J,K)-R)*U(J,K)+RP*WK(J-1,2+K)
            IF(K.GT.KL) WK(J,2+K)=WK(J,2+K)-WKA(4,J,K)*U(J,K-1)
            IF(K.LT.KU) WK(J,2+K)=WK(J,2+K)-WKA(6,J,K)*U(J,K+1)
3000      CONTINUE
          IF(WK(JU,1).EQ.0.0) RETURN

!	Back-substitution
          WK(JU,2+K)=WK(JU,2+K)/WK(JU,1)
          DO 3200 J=JU-1,JL,-1
3200      WK(J,2+K)=(WK(J,2+K)-WKA(3,J,K)*WK(J+1,2+K))/WK(J,1)
3500    CONTINUE

!	Gaussian elimination for the second half-step
        ERR=0.0
        DO 4800 K=JL,JU
          WK(KL,1)=WKA(5,K,KL)+R
          WK(KL,2)=WKA(7,K,KL)-(WKA(2,K,KL)-R)*WK(K,KL+2)
          IF(K.GT.JL) WK(KL,2)=WK(KL,2)-WKA(1,K,KL)*WK(K-1,KL+2)
          IF(K.LT.JU) WK(KL,2)=WK(KL,2)-WKA(3,K,KL)*WK(K+1,KL+2)
          DO 4200 J=KL+1,KU
            IF(WK(J-1,1).EQ.0.0) RETURN
            RP=-WKA(4,K,J)/WK(J-1,1)
            WK(J,1)=WKA(5,K,J)+R+RP*WKA(6,K,J-1)
            WK(J,2)=WKA(7,K,J)-(WKA(2,K,J)-R)*WK(K,J+2)+RP*WK(J-1,2)
            IF(K.GT.JL) WK(J,2)=WK(J,2)-WKA(1,K,J)*WK(K-1,J+2)
            IF(K.LT.JU) WK(J,2)=WK(J,2)-WKA(3,K,J)*WK(K+1,J+2)
4200      CONTINUE
          IF(WK(KU,1).EQ.0.0) RETURN

!	Back-substitution
          U1=WK(KU,2)/WK(KU,1)
          ERR=MAX(ERR,ABS(U1-U(K,KU)))
          U(K,KU)=U1
          DO 4400 J=KU-1,KL,-1
            U1=(WK(J,2)-WKA(6,K,J)*U(K,J+1))/WK(J,1)
            ERR=MAX(ERR,ABS(U1-U(K,J)))
            U(K,J)=U1
4400      CONTINUE

4800    CONTINUE
        IF(ERR.LT.AEPS) THEN
          IER=0
          NIT=IT
          RETURN
        ENDIF

6000  CONTINUE

!	Iteration fails to converge
      NIT=MAXIT
      IER=767
      END
