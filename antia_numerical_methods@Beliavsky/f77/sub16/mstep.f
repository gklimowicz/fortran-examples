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
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL DIF
      PARAMETER(NMX=100000,EPS=1.Q-30)
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
