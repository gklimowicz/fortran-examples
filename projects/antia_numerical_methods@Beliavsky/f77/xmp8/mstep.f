!    PROGRAM TO SOLVE INITIAL VALUE PROBLEMS IN ORDINARY DIFFERENTIAL EQUATIONS
!    USING MULTISTEP METHOD WITH ADAPTIVE STEP SIZE CONTROL

      PROGRAM MULTIS
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL DIF
      DIMENSION X(140),DX(140),WK(40),YF(20),IWK(20)

!     EXERCISE 12.34:  BESSEL'S EQUATION WITH N=0

51    FORMAT('   IER =',I4,5X,'T =',1PD14.6,5X,
     1   'NO. OF FUNCTION EVALUATIONS =',I8/5X,'STEP SIZE =',D14.6,5X
     2   ,'SOLUTION =',2D14.6/(2X,5D14.6))
52    FORMAT('   T0 =',1PD14.6,5X,'INITIAL STEP SIZE =',D14.6,5X,
     1   'IFLG =',I2,5X,'IST =',I2)

      N=2
      REPS=1.D-7
      IFLG=0
      NMAX=100000

      PRINT *,'TYPE T0=INITIAL TIME,  H=INITIAL STEP-SIZE,  IFLG=0/1,',
     1        '  IST=0/1'
      PRINT *,'IFLG=0 FOR ADAPTIVE STEP SIZE CONTROL'
      PRINT *,'IFLG=1 FOR CONSTANT STEP SIZE'
      PRINT *,'IST=0/1 FOR ADAMS/ GEAR'
      READ *,T0,H,IFLG,IST
      WRITE(6,52) T0,H,IFLG,IST

!     TO TACKLE THE SINGULARITY AT T=0 START THE INTEGRATION FROM
!     T0>0, AND USE THE MACLAURIN SERIES TO FIND THE INITIAL VALUES

      X(1)=1
      X(2)=0
      T1=1
      XT=(T0/2)**2
      DO 50 K=1,20
        T1=-T1*XT/(K*K)
        X(1)=X(1)+T1
        X(2)=X(2)+2*K*T1/T0
        IF(ABS(T1).LT.REPS) GO TO 100
50    CONTINUE

!     GO ON TYPING SUCCESSIVE VALUES OF T AT WHICH THE SOLUTION IS REQUIRED

100   PRINT *,'TYPE TN=FINAL VALUE OF TIME     (QUITS WHEN TN.LE.-10)'
      READ *,TN
      IF(TN.LT.-10.0) STOP
      CALL MSTEP(N,X,DX,DIF,H,T0,TN,YF,REPS,NSTEP,NMAX,IER,IFLG,IST,
     1           WK,IWK)
      WRITE(6,51) IER,TN,NSTEP,H,(YF(I),I=1,N)
      IF(IER.EQ.0) GO TO 100
      END

!     ------------------------------------------------------------

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

!     -----------------------------------------------------------

!	Solution of a system of linear equations using Gaussian elimination
!	with partial pivoting
!
!	N : (input) Number of equations to be solved
!	NUM : (input) Number of different sets (each with N equations) of
!	        equations to be solved
!	A : (input/output) The matrix of coefficient of size LJ*N
!	        A(I,J) is the coefficient of x_J in Ith equation
!	     	at output it will contain the triangular decomposition
!	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
!	        X(I,J) is the Ith element of Jth right hand side
!	     	at output it will contain the solutions
!	DET : (output) The determinant of the matrix
!	INC : (output) Integer array of length N containing information about
!		interchanges performed during elimination
!	LJ : (input) First dimension of arrays A and X in calling program
!	IER : (output) Error flag, IER=0 signifies successful execution
!		IER=101 implies (N.LE.0 or N.GT.LJ) 
!		IER=121 implies some pivot turned out to be zero and hence
!			matrix must be nearly singular
!	IFLG : (input) Integer parameter to specify the type of computation required
!		If IFLG.LE.0, both elimination and solution are
!			done and IFLG is set to 2
!		If IFLG=1, only elimination is done and IFLG is set to 2
!		If IFLG.GE.2 only solution is calculated, the triangular
!			decomposition should have been calculated earlier
!
!	Required routines : None

      SUBROUTINE GAUELM(N,NUM,A,X,DET,INC,LJ,IER,IFLG)
      IMPLICIT REAL*8(A-H,O-Z)
!	For complex matrices use the following statements instead
!      IMPLICIT REAL*8(R)
!      IMPLICIT COMPLEX*16(A-H,S-Z)

      DIMENSION A(LJ,N),INC(N),X(LJ,NUM)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=101
        RETURN
      ENDIF

      IER=121
      IF(IFLG.LE.1) THEN
!	Perform elimination

        DET=1.0
        DO 2600 K=1,N-1
!	Find the maximum element in the Kth column
          R1=0.0
          KM=K
          DO 2200 L=K,N
            IF(ABS(A(L,K)).GT.R1) THEN
              R1=ABS(A(L,K))
              KM=L
            ENDIF
2200      CONTINUE

          INC(K)=KM
          IF(KM.NE.K) THEN
!	Interchange the rows if needed
            DO 2300 L=K,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
2300        A(KM,L)=T1
            DET=-DET
          ENDIF

          DET=DET*A(K,K)
          IF(A(K,K).EQ.0.0) RETURN
!	To check for singular or nearly singular matrices replace this
!	statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(K,K)).LT.REPS) RETURN
          DO 2500 L=K+1,N
            A(L,K)=A(L,K)/A(K,K)
            DO 2500 L1=K+1,N
2500      A(L,L1)=A(L,L1)-A(L,K)*A(K,L1)
2600    CONTINUE
        DET=DET*A(N,N)
        INC(N)=N
!	If pivot is zero then return, IER has been set to 121
        IF(A(N,N).EQ.0.0) RETURN
!	To check for singular or nearly singular matrices replace this
!	statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(N,N)).LT.REPS) RETURN

        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
!	Solution for the NUM different right-hand sides
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
!	Forward substitution
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF
          DO 3000 L=K+1,N
3000    X(L,J)=X(L,J)-A(L,K)*X(K,J)

!	back-substitution
        X(N,J)=X(N,J)/A(N,N)
        DO 3300 K=N-1,1,-1
          DO 3200 L=N,K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L)
3300    X(K,J)=X(K,J)/A(K,K)
5000  CONTINUE
      END

!     -----------------------------------------------------------

!	To perform one step of solution of initial value problems in ordinary
!	differential equations using fourth order stiffly stable method.
!	It is called by subroutine MSTEP.
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
!	IFLAG : (input/output) Integer variable used as a flag
!		If IFLAG=0 initial approximation to Jacobian is generated
!			and IFLAG is set to 1
!		Otherwise old approximation to J^{-1} is used.
!	IER : (output) Error parameter; IER=0 implies successful execution
!		IER=729 implies that iteration on corrector failed to converge
!			or cannot be continued because the estimated Jacobian is singular.
!	WK1 : (input/output) Real array of length 3N used to pass the
!		predicted values
!	WK : Real array of length N*MAX(2N, N+4) used as scratch space
!	IWK : Integer array of length N used as scratch space
!
!	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
!	the differential equation. T is the value of independent variable,
!	N is the number of variables, Y is a real array of length N containing
!	the values of variables. DY is a real array of length N which should
!	contain the calculated values of derivatives at (T,Y).
!
!	Required routines : GAUELM, DIF
!	
!	
      SUBROUTINE GEAR(N,Y,DY,DIF,H,T,REPS,NSTP,IJ,IJM1,IJM2,
     1                IJM3,IJM4,IFLAG,IER,WK1,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(EPS=1.D-30,NIT=20,CFAC=1.D-2)
      DIMENSION Y(N,7),DY(N,7),WK1(N,3),WK(N,*),IWK(N)

      IER=0
      DO 1000 J=1,N
!	predictor
        T1=Y(J,IJM1)+H*(55.*DY(J,IJM1)-59.*DY(J,IJM2)+37.*DY(J,IJM3)
     1     -9.*DY(J,IJM4))/24.
!	Modifier
        Y(J,IJ)=T1+6275.*(Y(J,IJM1)-WK1(J,2))/8003.
1000  WK1(J,1)=T1

      LJ=N
      IPAS=0
      IER=0
      NSTP=NSTP+1
      CALL DIF(T,N,Y(1,IJ),DY(1,IJ))
      DO 1200 J=1,N
!	Residual in the corrector
        WK1(J,3)=Y(J,IJ)-(48.*Y(J,IJM1)-36.*Y(J,IJM2)+16.*Y(J,IJM3)-
     1           3.*Y(J,IJM4)+12.*H*DY(J,IJ))/25.
1200  CONTINUE

      IF(IFLAG.EQ.0) THEN
!	Generate initial approximation to Jacobian J 
        DO 2000 IP=1,N
          X1=Y(IP,IJ)
          DK=100.*REPS*ABS(X1)
          IF(DK.EQ.0.0) DK=100.*REPS
          Y(IP,IJ)=X1+DK
          NSTP=NSTP+1
          CALL DIF(T,N,Y(1,IJ),DY(1,IJ))
          DO 1800 J=1,N
            WK(J,1)=Y(J,IJ)-(48.*Y(J,IJM1)-36.*Y(J,IJM2)+16.*Y(J,IJM3)-
     1              3.*Y(J,IJM4)+12.*H*DY(J,IJ))/25.
1800      WK(J,N+IP)=(WK(J,1)-WK1(J,3))/DK
          Y(IP,IJ)=X1
2000    CONTINUE

        NUM=N
        IFLG=0
        DO 2500 I=1,N
          DO 2400 J=1,N
2400      WK(J,I)=0.0
2500    WK(I,I)=1.0

!	Calculate inverse of Jacobian
        CALL GAUELM(N,NUM,WK(1,N+1),WK,DET,IWK,LJ,IER1,IFLG)
        IF(IER1.GT.0) THEN
          IER=729
          RETURN
        ENDIF
!	Reset the flag so that Jacobian is not calculated again
        IFLAG=1
      ENDIF

      DO 2800 I=1,N
        WK(I,N+1)=WK1(I,3)
        WK(I,N+2)=Y(I,IJ)
2800  CONTINUE

!	Iteration on parameters using Broyden's method
      DO 5000 IPAS=1,NIT

!	The convergence check
        ERR=0.0
        DO 3500 J=1,N
          S1=0.0
!	Correction to the Jth component
          DO 3200 K=1,N
3200      S1=S1+WK(J,K)*WK(K,N+1)
          R2=ABS(Y(J,IJ))+ABS(Y(J,IJ)-Y(J,IJM1))+EPS
          ERR=MAX(ERR,ABS(S1/R2))
          Y(J,IJ)=WK(J,N+2)-S1
          WK(J,N+3)=-S1
3500    CONTINUE

        NSTP=NSTP+1
        CALL DIF(T,N,Y(1,IJ),DY(1,IJ))
        IF(ERR.LT.CFAC*REPS) RETURN

!	calculate the residuals
        DO 3800 J=1,N
          WK1(J,3)=Y(J,IJ)-(48.*Y(J,IJM1)-36.*Y(J,IJM2)+16.*Y(J,IJM3)-
     1             3.*Y(J,IJM4)+12.*H*DY(J,IJ))/25.
3800    CONTINUE

!	Update J inverse using Broyden's formula
        DO 4000 J=1,N
          WK(J,N+4)=WK1(J,3)-WK(J,N+1)
          WK(J,N+1)=WK1(J,3)
          WK(J,N+2)=Y(J,IJ)
4000    CONTINUE
        SS=0.0
        DO 4300 J=1,N
          S1=0.0
          S2=0.0
          DO 4100 K=1,N
            S1=S1+WK(K,N+3)*WK(K,J)
            S2=S2+WK(J,K)*WK(K,N+4)
4100      CONTINUE
          WK1(J,3)=S1
          Y(J,IJ)=S2-WK(J,N+3)
          SS=SS+S1*WK(J,N+4)
4300    CONTINUE

        IF(SS.EQ.0.0) THEN
          IER=729
          RETURN
        ENDIF
        DO 4400 J=1,N
          DO 4400 K=1,N
            WK(K,J)=WK(K,J)-Y(K,IJ)*WK1(J,3)/SS
4400    CONTINUE
5000  CONTINUE

!	Iteration fails to converge
      IER=729
      END

!   --------------------------------------------------

!	To solve initial value problems in ordinary differential equations
!	using a fourth-order multistep method with adaptive step size control
!	It can be used with ADAMS (Adams-Bashforth-Moulton predictor corrector method)
!	or with GEAR (Gear's stiffly stable method)
!
!	N : (input) Number of first order differential equations to be solved
!	Y : (input/output) Real array of length 7N, the first N elements
!		should contain the initial values of variables. During
!		execution, the solution at various intermediate points will
!		be stored in this array. The contents of this array must
!		be preserved between two calls to MSTEP if the integration
!		needs to be continued further.
!	DY : (output) Real array of length 7N, containing the derivatives
!		of Y the points stored in Y. This array should also be
!		preserved between two calls to MSTEP.
!	DIF : (input) Name of subroutine to calculate the right hand side
!		of differential equation y'=f(t,y)
!	H : (input/output) Initial guess for the step size. After execution
!		it will contain the step size used by the program
!	T0 : (input/output) Initial value of independent variable t, at
!		which initial values are specified. After execution it will
!		be set to the point up to which integration has been successful
!	TN : (input) The final value of t at which the solution is required.
!		Intermediate values will not be preserved so if solution
!		is required at intermediate points, TN must be set to first
!		such value and multiple calls will be needed to calculate
!		all required values. For each subsequent call only TN needs
!		to be updated. Other variables including the scratch arrays
!		in the call statement must be preserved to their old values.
!	YF : (output) Real array of length N, containing the solution at
!		the last point up to which integration is successful.
!		If execution is successful it will contain the required
!		solution at TN.
!	REPS : (input) Required accuracy in each component of the solution.
!		The subroutine only controls local truncation error and hence
!		actual error could be larger
!	NSTP : (input/output) Number of calls to DIF required since the starting
!		of integration. The count is updated by the subroutine.
!	NMAX : (input/output) Maximum number of calls to DIF to be tried.
!		If NMAX.LE.0 it will be set to a default value of NMX=100000.
!	IER : (output) Error parameter; IER=0 implies successful execution
!		IER=702 implies N.LE.0, no calculations are done
!		IER=724 or 725 implies STRT4 failed to generate starting values
!		IER=726 implies that step-size has become smaller than
!			REPS*|TN-T0|
!		IER=727 implies that step size is too small for arithmetic used
!		IER=728 implies that integration could not be completed in
!			the specified number of steps.
!	IFLG : (input/output) Integer variable used as flag to decide the
!		type of computation.
!		If IFLG=0 the computation is started afresh with starting
!			values being generated and the step size is adjusted
!			to meet the specified accuracy. NSTP is also initialised
!			to zero in the beginning. After execution IFLG is set to 2
!		If IFLG=1 the computation is started afresh with starting
!			values being generated but the step size is kept fixed.
!			NSTP is also initialised to zero in the beginning.
!			After execution IFLG is set to 3
!		If IFLG=2 the integration is continued using already available
!			solution at previous steps. The step size is adjusted
!			to meet the specified accuracy.
!		If IFLG=3 the integration is continued using already available
!			solution at previous steps. The step size is kept fixed.
!		IFLG should not be changed after the first call, unless
!		fresh starting values need to be generated.
!	IST : (input/output) Integer variable used as flag to decide the
!		multistep method to be used:
!		If IST=0 Adams method is used, which may be good for
!			general equations
!		If IST=1 Gear's method is used, which may be used for
!			stiff equations
!	WK : Real array used as scratch space. For use with ADAMS the size
!		of WK should be 2N, while for use with GEAR it should be
!		3N+N*MAX(N+4, 2N). This scratch space must be preserved
!		between two calls to the routine
!	IWK : Integer array of length N used as scratch space.
!
!	Subroutine DIF(T,N,Y,DY) must be supplied by the user to specify
!	the differential equation. T is the value of independent variable,
!	N is the number of variables, Y is a real array of length N containing
!	the values of variables. DY is a real array of length N which should
!	contain the calculated values of derivatives at (T,Y).
!
!	Required routines : ADAMS, GEAR, STRT4, RK4, GAUELM, DIF
!	
!	
      SUBROUTINE MSTEP(N,Y,DY,DIF,H,T0,TN,YF,REPS,NSTP,NMAX,IER,IFLG,
     1                 IST,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL DIF
      PARAMETER(NMX=100000,EPS=1.D-30)
      DIMENSION Y(N,7),DY(N,7),WK(N,*),YF(N),DEL(4),IWK(*)
      SAVE

      IF(N.LE.0) THEN
        IER=702
        RETURN
      ENDIF
      IF(NMAX.LE.0) NMAX=NMX
      IER=0
      TSTEP=TN-T0

      IF(IFLG.LE.1) THEN
        IF(TSTEP.EQ.0.0) RETURN
!	Generate the starting values
        NSTP=0
        IF(H.EQ.0.0) H=TSTEP/8
!	Change the sign of H if needed
        IF(H.LT.0.0.EQV.TN.GT.T0) H=-H

        T=T0
!	Generate the starting values
        CALL STRT4(N,Y,DY,DIF,H,T,REPS,IFLG,TSTEP,NSTP,IER,WK)
        IF(IER.GT.0) RETURN

!	Initialising the array indices
        IJ=4
        IJM1=3
        IJM2=2
        IJM3=1
        T0=T
!	Set the flag for subroutine GEAR
        IFLAG=0
!	Update IFLG
        IFLG=IFLG+2
        NS=4
      ENDIF

      IF(TN.LE.T0.EQV.H.GT.0) GO TO 6000
      IF(TSTEP.EQ.0.0) GO TO 6000
!	Updating the array indices
3400  IJM4=IJM3
      IJM3=IJM2
      IJM2=IJM1
      IJM1=IJ
      IJ=IJ+1
      IF(IJ.GT.7) IJ=1
      T=T0+H
      REPS1=ABS(REPS*H/TSTEP)

!	To perform one step of integration using Adams method
      IF(IST.EQ.0) THEN
        CALL ADAMS(N,Y,DY,DIF,H,T,REPS1,NSTP,IJ,IJM1,IJM2,
     1             IJM3,IJM4,IER,WK)
      ELSE

!	To perform one step of integration using stiffly stable method
        CALL GEAR(N,Y,DY,DIF,H,T,REPS1,NSTP,IJ,IJM1,IJM2,
     1            IJM3,IJM4,IFLAG,IER,WK,WK(1,4),IWK)
      ENDIF

!	Estimating the truncation error
      ERR=0
      DO 4500 I=1,N
        R2=ABS(Y(I,IJ))+ABS(Y(I,IJ)-Y(I,IJM1))+EPS
        R1=ABS(WK(I,1)-Y(I,IJ))/R2
        IF(R1.GT.ERR) ERR=R1
4500  CONTINUE
      ERR=ABS(0.25*ERR*TSTEP/H)

      IF(T.EQ.T0) THEN
!	Step size is too small for the arithmetic used
        IER=727
        DO 4600 I=1,N
4600    YF(I)=Y(I,IJ)
        RETURN
      ENDIF

      IF(IFLG.LE.2) THEN
        IF(ERR.LT.REPS.AND.IER.EQ.0) THEN
!	Integration is successful
          T0=T0+H
          DO 4700 I=1,N
4700      WK(I,2)=WK(I,1)
          NS=NS+1
          IF(TN.LE.T0.EQV.H.GT.0) GO TO 6000
          IF(ABS((TN-T0)/TSTEP).LE.REPS) GO TO 6000

          IF(ERR.LT.REPS/32.AND.NS.GT.7) THEN
!	Double the step size
            H=2.*H
            I2=IJM4-2
            IF(I2.LE.0) I2=I2+7
            DO 4800 I=1,N
              Y(I,IJM1)=Y(I,IJM2)
              Y(I,IJM2)=Y(I,IJM4)
              Y(I,IJM3)=Y(I,I2)
              DY(I,IJM1)=DY(I,IJM2)
              DY(I,IJM2)=DY(I,IJM4)
              DY(I,IJM3)=DY(I,I2)
4800        CONTINUE
            NS=4
          ENDIF

        ELSE
!	If integration has failed, then reduce the step size
          IF(NS.GT.-4.AND.ERR.LT.16.*REPS) THEN
!	H is halved
            DO 5000 I=1,N
              Y1=(45.*Y(I,IJM1)+72.*Y(I,IJM2)+11.*Y(I,IJM3))/128.+
     1            H*(-9.*DY(I,IJM1)+36.*DY(I,IJM2)+3.*DY(I,IJM3))/128.
              Y2=(11.*Y(I,IJM1)+72.*Y(I,IJM2)+45.*Y(I,IJM3))/128.+
     1            H*(-3.*DY(I,IJM1)-36.*DY(I,IJM2)+9.*DY(I,IJM3))/128.
              Y(I,IJ)=Y(I,IJM1)
              Y(I,IJM1)=Y1
              Y(I,IJM3)=Y2
              DY(I,IJ)=DY(I,IJM1)
5000        CONTINUE
            T=T0-H/2
            CALL DIF(T,N,Y(1,IJM1),DY(1,IJM1))
            T=T-H
            CALL DIF(T,N,Y(1,IJM3),DY(1,IJM3))
            NSTP=NSTP+2
            H=H/2.
          ELSE

!	If error is too large or the halving has failed once, then
!	generate fresh starting values with smaller h
            IFLG=IFLG-2
            T=T0
            H=H/8.
            DO 5200 I=1,N
5200        Y(I,1)=Y(I,IJM1)
            CALL STRT4(N,Y,DY,DIF,H,T,REPS,IFLG,TSTEP,NSTP,IER,WK)
            IF(IER.GT.0) THEN
!	If STRT4 fails, then quit
              DO 5300 I=1,N
5300          YF(I)=Y(I,IJM1)
              RETURN
            ENDIF
            IJ=4
            IJM1=3
            IJM2=2
            IJM3=1
            T0=T
            IFLG=IFLG+2
            IFLAG=0
          ENDIF

          NS=-4
          IF(ABS(H/TSTEP).LT.REPS) THEN
!	The step size is too small
            IER=726
            DO 5400 I=1,N
5400        YF(I)=Y(I,IJM1)
            RETURN
          ENDIF

        ENDIF
      ELSE
        IF(IER.GT.0) THEN
!	For integration with fixed H, quit if corrector fails
          DO 5600 I=1,N
5600      YF(I)=Y(I,IJM1)
          RETURN
        ENDIF
        T0=T0+H
        DO 5700 I=1,N
5700    WK(I,2)=WK(I,1)
        IF(TN.LE.T0.EQV.H.GT.0) GO TO 6000
        IF(ABS((TN-T0)/TSTEP).LE.REPS) GO TO 6000
      ENDIF

      IF(NSTP.LT.NMAX) GO TO 3400
!	Quit if the specified number of function evaluations are exceeded
      IER=728
      DO 5800 I=1,N
5800  YF(I)=Y(I,IJ)
      RETURN

6000  TX=(TN-T0)/H
!	Use interpolation to calculate solution at TN
      DO 7000 I=1,N
        DEL(1)=DY(I,IJ)*H
        DEL(2)=Y(I,IJ)-Y(I,IJM1)
        DEL(3)=DY(I,IJM1)*H
        DEL(4)=Y(I,IJM1)-Y(I,IJM2)
        DO 6500 J1=2,3
          DO 6500 J2=4,J1,-1
6500    DEL(J2)=DEL(J2-1)-DEL(J2)
        DEL(4)=DEL(3)-0.5*DEL(4)
        YF(I)=Y(I,IJ)+TX*(DEL(1)+TX*(DEL(2)+(TX+1)*(DEL(3)+
     1        (TX+1)*DEL(4)/2.)))
7000  CONTINUE
      END

!     ----------------------------------------------------

!	To perform one step of integration of ordinary differential equations
!	using a fourth-order Runge-Kutta method 
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
      SUBROUTINE RK4(N,T,Y0,DY0,H,Y1,DIF,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y0(N),DY0(N),Y1(N),WK(N,2)

      H2=H/2.
      DO 1000 I=1,N
        Y1(I)=H2*DY0(I)
!	y_0+0.5k_1
1000  WK(I,1)=Y0(I)+Y1(I)

      T1=T+H2
      CALL DIF(T1,N,WK,WK(1,2))
      DO 1200 I=1,N
        Y1(I)=Y1(I)+H*WK(I,2)
!	y_0+0.5k_2
1200  WK(I,1)=Y0(I)+H2*WK(I,2)

      CALL DIF(T1,N,WK,WK(1,2))
      DO 1400 I=1,N
        Y1(I)=Y1(I)+H*WK(I,2)
!	y_0+k_3
1400  WK(I,1)=Y0(I)+H*WK(I,2)

      T1=T+H
      CALL DIF(T1,N,WK,WK(1,2))
      DO 1600 I=1,N
1600  Y1(I)=Y0(I)+(Y1(I)+H2*WK(I,2))/3.
      END

!     ---------------------------------------------------

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
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(EPS=1.D-30,SFAC=0.9D0,NIT=10)
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
          H=SFAC*H*(REPS/ERR)**0.25D0
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

!     ---------------------------------------------------

      SUBROUTINE DIF(T,N,X,DX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(*),DX(*)

!     BESSEL'S EQUATION FOR N=0
!     Y1'=Y2,   Y2'=-Y2/T-Y1

      DX(1)=X(2)
      DX(2)=-X(2)/T-X(1)
      END
