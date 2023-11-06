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
!      IMPLICIT REAL*8(A-H,O-Z)
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
