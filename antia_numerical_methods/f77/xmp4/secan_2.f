!     PROGRAM TO FIND REAL ROOTS OF A NONLINEAR EQUATION USING SECANT METHOD
!     THE FUNCTION IS CALCULATED IN A SCALED FORM

      PROGRAM ZERO
!      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL CDET
      DIMENSION CZERO(20),CX(20)

!     EXERCISE 7.44:  EIGENVALUES OF HILBERT MATRIX

51    FORMAT('   IER =',I4,3X,' ROOT USING SECAN_2',1PE14.6)
52    FORMAT('   STARTING VALUE =',1PE14.6,4X,'LOWER LIMIT =',E14.6,
     1        3X,'UPPER LIMIT =',E14.6)
53    FORMAT(/'   IER =',I4,4X,' ROOT USING SECANI',1PE14.6,3X,
     1       'FUNCTION VALUE=',E14.6,I7)

      REPS=1.E-6
      AEPS=1.E-17
      NZ=0


100   PRINT *,'TYPE X0=STARTING VALUE,  XL=LOWER LIMIT,  XU=UPPER LIMIT'
      PRINT *,'                            (QUITS WHEN XL.EQ.XU)'
      READ *,X0,XL,XU
      IF(XL.EQ.XU) STOP
      WRITE(6,52) X0,XL,XU
      CALL SECAN_2(X0,XL,XU,X,REPS,AEPS,IER,CDET)
      WRITE(6,51) IER,X

!     USING SECANI WITH REVERSE COMMUNICATION TECHNIQUE

      IER=0
500   CALL SECANI(X0,XL,XU,X,F0,JF,REPS,AEPS,IER)
      IF(IER.LT.0) THEN
        F0=CDET(X,JF)
        GO TO 500
      ENDIF
      WRITE(6,53) IER,X,F0,JF
      GO TO 100
      END

!     ---------------------------------------------

!	Solution of a system of linear equations using Crout's algorithm 
!	with partial pivoting
!
!	N : (input) Number of equations to be solved
!	NUM : (input) Number of different sets (each with N equations) of
!	         equations to be solved
!	A : (input/output) The matrix of coefficient of size LJ*N
!	         A(I,J) is the coefficient of x_J in Ith equation
!	     at output it will contain the triangular decomposition
!	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
!	        X(I,J) is the Ith element of Jth right hand side
!	     	at output it will contain the solutions
!	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
!	INC : (output) Integer array of length N containing information about
!		interchanges performed during elimination
!	LJ : (input) First dimension of arrays A and X in calling program
!	IER : (output) Error flag, IER=0 signifies successful execution
!		IER=102 implies (N.LE.0 or N.GT.LJ) 
!		IER=122 implies some pivot turned out to be zero and hence
!			matrix must be nearly singular
!	IFLG : (input) Integer parameter which determines the type of computation
!		required.
!		If IFLG.LE.0, both elimination and solution are calculated
!			and IFLG is set to 2
!		If IFLG=1, only elimination is done and IFLG is set to 2
!		If IFLG.GE.2 only solution is calculated, the triangular
!		    	decomposition should have been calculated earlier
!	WK : scratch array of length N
!
!	Required routines : None

      SUBROUTINE CROUT(N,NUM,A,X,DET,IDET,INC,LJ,IER,IFLG,WK)
!      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION A(LJ,N),INC(N),X(LJ,NUM),WK(N)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=102
        RETURN
      ENDIF
      IER=122

      IF(IFLG.LE.1) THEN
!	Perform LU decomposition
        DO 1500 I=1,N
          R1=0.0
          DO 1200 J=1,N
            IF(ABS(A(I,J)).GT.R1) R1=ABS(A(I,J))
1200      CONTINUE
          WK(I)=R1
!	If any row is zero, then quit
          IF(R1.EQ.0.0) RETURN
1500    CONTINUE

        DET=1.0
        IDET=0
        DO 2600 K=1,N
          R1=0.0
          KM=K
!	Generate the Kth column of L
          DO 2200 L=K,N
            D1=A(L,K)
            DO 2100 L1=1,K-1
2100        D1=D1-A(L,L1)*A(L1,K)
            A(L,K)=D1

!	Finding the pivot
            R2=ABS(D1/WK(L))
            IF(R2.GT.R1) THEN
              R1=R2
              KM=L
            ENDIF
2200      CONTINUE

          INC(K)=KM
!	Interchange the rows if needed
          IF(KM.NE.K) THEN
            DET=-DET
            DO 2300 L=1,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
2300        A(KM,L)=T1
            T1=WK(K)
            WK(K)=WK(KM)
            WK(KM)=T1
          ENDIF

          DET=DET*A(K,K)
!	If the pivot is zero, then quit
          IF(A(K,K).EQ.0.0) RETURN
!	To check for singular or nearly singular matrices replace this statement by
!         IF(ABS(A(K,K)).LT.REPS) RETURN

          IF(DET.NE.0.0) THEN
!	Scale the value of the determinant DET
2350        IF(ABS(DET).GT.32.) THEN
              DET=DET*0.03125D0
              IDET=IDET+5
              GO TO 2350
            ENDIF

2370        IF(ABS(DET).LT.0.03125D0) THEN
              DET=DET*32.
              IDET=IDET-5
              GO TO 2370
            ENDIF
          ENDIF

!	Generate the Kth row of U
          DO 2500 L=K+1,N
            D1=A(K,L)
            DO 2400 L1=1,K-1
2400        D1=D1-A(K,L1)*A(L1,L)
            A(K,L)=D1/A(K,K)
2500      CONTINUE
2600    CONTINUE
        IER=0

        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
!	Solution for NUM different right-hand sides
      DO 5000 J=1,NUM
!	Forward substitution
        DO 3000 K=1,N
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF

          D1=X(K,J)
          DO 2800 L=1,K-1
2800      D1=D1-A(K,L)*X(L,J)
          X(K,J)=D1/A(K,K)
3000    CONTINUE

!	Back-substitution
        DO 3300 K=N-1,1,-1
          D1=X(K,J)
          DO 3200 L=N,K+1,-1
3200      D1=D1-X(L,J)*A(K,L)
3300    X(K,J)=D1
5000  CONTINUE
      END

!      ------------------------------------------

!	Real zero of a given function using secant iteration
!	Function is calculated as FX*2**JF
!
!	X0 : (input) Initial guess for the zero
!	XL : (input) Lower limit of interval where zero is expected
!	XU : (input) Upper limit of interval where zero is expected
!	X : (output) Computed value of the zero
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=40 implies that function value is equal at two points
!			and it is not possible to continue the iteration
!		IER=402 implies XL>X0 or XU<X0, in which case no calculations are done
!		IER=422 implies that iteration goes outside the specified limits
!		IER=423 implies that iteration failed to converge to specified accuracy
!	FUN : (input) Name of the function routine to calculate the function
!		FUNCTION FUN(X,JF) must be supplied by the user.
!		The function value should be FUN*2**JF
!
!	Required routines : FUN

      SUBROUTINE SECAN_2(X0,XL,XU,X,REPS,AEPS,IER,FUN)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIS=75)

      IER=0
      IF(XL.GT.X0.OR.XU.LT.X0) THEN
!	If X0 is outside the specified interval (XL,XU) then quit
        IER=402
        RETURN
      ENDIF

      X=X0
!	Select the increment for the next point X+DX
      DX=(XU-X0)/100.
      IF(ABS(DX).LT.5.*MAX(REPS*ABS(X),AEPS)) DX=(XL-X0)/5.
      IF(ABS(DX).GT.0.1D0*MAX(ABS(X),100.*AEPS))
     1    DX=SIGN(0.10*MAX(ABS(X),100.*AEPS),DX)
      F1=0.0
      JF1=0

      DO 1000 L=1,NIS
        F=FUN(X,JF)
        DX1=DX
        F1=F1*2.D0**(JF1-JF)

        IF(F1-F.EQ.0.0) THEN
          IF(F.EQ.0.0) RETURN
!	If F1=F and F.NE.0, then quit
          IER=40
          RETURN
        ENDIF

!	The secant iteration
        IF(L.GT.1) DX=DX1*F/(F1-F)
        X=X+DX
        F1=F
        JF1=JF

        IF(ABS(DX).LT.MAX(REPS*ABS(X),AEPS).AND.L.GT.2) RETURN
        IF(X.LT.XL.OR.X.GT.XU) THEN
!	If iteration goes outside the specified limits (XL,XU), then quit
          IER=422
          RETURN
        ENDIF
1000  CONTINUE

!	The iteration fails to converge
      IER=423
      END

!      --------------------------------------------------------

!     Real zero of a given function using secant iteration
!     Function is calculated as FX*2**JF
!     This subroutine uses reverse communication to calculate function
!     values. If IER<0 the function should be evaluated and SECANI should
!     be called back with new function value. Calculation of function
!     value should not change any other variables in the call statement.
!
!     X0 : (input) Initial guess for the zero
!     XL : (input) Lower limit of interval where zero is expected
!     XU : (input) Upper limit of interval where zero is expected
!     X : (output) Value of x at which the function evaluation is required.
!		If IER=0 then it will contain the final value of zero computed
!		by the routine.
!     F : (input) Calculated value of the function at X.
!		If subroutine exits with IER<0, then the calling routine should
!		calculate the function value at X and call SECANI with this value
!		stored in F and JF. Other variables should not be changed.
!     JF : (input) The exponent of function value, the function value
!		should be F*2**JF
!     REPS : (input) Required relative accuracy
!     AEPS : (input) Required absolute accuracy
!     		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!     IER : (input/output) Error parameter, IER=0 implies successful execution
!		Before the first call IER should be set to zero
!		IER<0 implies that execution is not over and the subroutine needs
!			a new function evaluation at X. After calculating the
!			function value SECANI should be called back.
!     		IER=40 implies that function value is equal at two points
!     			and it is not possible to continue the iteration
!     		IER=402 implies XL>X0 or XU<X0, in which case no calculations are done
!     		IER=422 implies that iteration goes outside the specified limits
!     		IER=423 implies that iteration failed to converge to specified accuracy
!
!	Required routines : None (Function is calculated by the calling program)

      SUBROUTINE SECANI(X0,XL,XU,X,F,JF,REPS,AEPS,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIS=75)
      SAVE
 
      IFL=-IER
!     Jump to proper location and continue execution
      IF(IFL.EQ.1) GO TO 1500

!     Initial call to subroutine, start from beginning
      IF(XL.GT.X0.OR.XU.LT.X0) THEN
!     If X0 is outside the specified interval (XL,XU) then quit
        IER=402
        RETURN
      ENDIF
 
      X=X0
!     Select the increment for the next point X+DX
      DX=(XU-X0)/100.
      IF(ABS(DX).LT.5.*MAX(REPS*ABS(X),AEPS)) DX=(XL-X0)/5.
      IF(ABS(DX).GT.0.1D0*MAX(ABS(X),100.*AEPS))
     1       DX=SIGN(0.10*MAX(ABS(X),100.*AEPS),DX)
      F1=0.0
      JF1=0
      L=0
 
1000  L=L+1
      IER=-1
!     To evaluate the function at X
      RETURN

1500  DX1=DX
      F1=F1*2.D0**(JF1-JF)
      IER=0
 
      IF(F1-F.EQ.0.0) THEN
        IF(F.EQ.0.0) RETURN
!     If F1=F and F.NE.0, then quit
        IER=40
        RETURN
      ENDIF
 
!     The secant iteration
      IF(L.GT.1) DX=DX1*F/(F1-F)
      X=X+DX
      F1=F
      JF1=JF
 
      IF(ABS(DX).LT.MAX(REPS*ABS(X),AEPS).AND.L.GT.2) RETURN
      IF(X.LT.XL.OR.X.GT.XU) THEN
!     If iteration goes outside the specified limits (XL,XU), then quit
        IER=422
        RETURN
      ENDIF
      IF(L.LT.NIS) GO TO 1000
 
!     The iteration fails to converge
      IER=423
      END

!      --------------------------------------------------------

      FUNCTION CDET(CX,IDET)
      IMPLICIT LOGICAL(Q)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION CR(20),CW(20),CA(20,20),INT(20)

!     FUNCTION ROUTINE TO CALCULATE THE DETERMINANT OF HILBERT MATRIX

!     THE ORDER OF MATRIX
      N=10
      NUM=1
      LJ=20
      IFLG=1

!     SETTING UP THE HILBERT MATRIX  (H-X I)
      DO 1000 I=1,N
        DO 800 J=1,N
          CA(I,J)=1.D0/(I+J-1.)
800     CONTINUE
        CA(I,I)=CA(I,I)-CX
1000  CONTINUE

      CALL CROUT(N,NUM,CA,CR,CDET,IDET,INT,LJ,IER,IFLG,CW)
      END


