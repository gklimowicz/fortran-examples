!	To solve a linear inversion problem in one dimension using RLS
!	technique with B-spline basis functions
!
!	NK : (input) Number of knots for defining B-splines, the number
!		of basis functions would be NK+K-2
!	XO : (input) Real array of length NK containing the knots
!		used for defining B-spline basis functions.
!		The knots must be distinct and in ascending order.
!	K : (input) Order of B-splines required, K=4 gives cubic B-splines
!	NR : (input) Number of points used in defining the kernels
!	R : (input) Real array of length NR containing the coordinates
!		of points at which kernels are available.
!	RKER : (input) Real array of length IK*NR containing the kernels
!		for the inverse problem. RKER(I,J) should contain the
!		value at R(J) for the Ith kernel. This array
!		must be supplied if IFLG<2, otherwise it is not required
!	IK : (input) First dimension of arrays RKER, AC and A, as specified
!		in the calling program. IK.GE.NM+NS
!	AC : (input/output) Real array of length IK*(NK+K-2) containing
!		the coefficients of matrix defining the inversion problem
!		If IFLG<2, these coefficients are calculating by integrating
!		the kernels with appropriate weights. For IFLG=2,3 these
!		coefficients must be supplied.
!	NM : (input) Number of data points in the inversion problem
!	NS : (input) Number of points to be used for applying regularisation
!		The routine chooses a uniform mesh covering the full interval
!		for applying smoothing.
!	ALP : (input) Regularisation parameter, ALP>0.
!	IDE : (input) Order of derivative to be used for regularisation,
!		IDE should be 1 or 2 for first or second derivative smoothing
!	DI : (input) Real array of length NM, containing the data points
!		for inversion
!	DE : (input) Real array of length NM, containing the estimated
!		error in DI. 
!	DF : (output) Real array of length NM, containing the normalised
!		residuals (DI-DI(fit))/DE for each data point.
!	F : (output) Real array of length NR which will contain the
!		calculated solution at each point in array R.
!	B : (output) Real array of length NM+NS containing the coefficients
!		of basis functions in fitted solution. Although the
!		the number of coefficients is only NK+K-2, the rest
!		of array is used as scratch space
!	IFLG : (input/output) Integer specifying the type of calculation required.
!		IFLG=0 : The matrix coefficients will be calculated using
!			the kernels and then the equations are solved to
!			find the solution for given data points DI.
!			IFLG is set to 4 after calculations.
!		IFLG=1 : The matrix coefficients will be calculated using
!			the kernels and then the SVD of the full matrix is
!			computed, but the solution is not computed.
!			IFLG is set to 4 after calculations.
!		IFLG=2 : The matrix coefficients are assumed to be available
!			in array AC and the matrix is setup and solved
!			to find the solution for given data points DI.
!			IFLG is set to 4 after calculations.
!		IFLG=3 : The matrix coefficients are assumed to be available
!			in array AC and the matrix is setup and the SVD
!			is computed, but the solution is not computed.
!			IFLG is set to 4 after calculations.
!		IFLG=4 : The SVD of matrix is assumed to be available
!			from previous calculations and only the solution
!			for given DI is computed.
!		Since IFLG is set to 4 every-time, it should be reset to
!		0 or 2 before next call when the data or error estimates
!		or smoothing are changed.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=709 implies that NM.LE.NK+K-2 or IK<NM+NS or IV<NK+K-2
!		IER=710 implies that ALP<0 or IDE<1 or IDE>2
!		other values of IER may be set by BSPLIN or SVD or BSPEVL
!	REPS : (input) Required accuracy for solution of equations using
!		SVD. Singular values less than REPS times maximum will be
!		set to zero.
!	CHISQ : (output) The computed value of Chi square for the solution
!	SUMD : (output) The computed value of the smoothing term
!	A : (input/output) Real array of length IK*(NK+K-2) containing
!		the SVD of the matrix of equations. If IFLG<4 this matrix
!		will be calculated, otherwise it must be supplied.
!	AV : (input/output) Real array of length IV*(NK+K-2) containing
!		the matrix V or SVD of the matrix of equations. if IFLG<4
!		this matrix will be calculated, otherwise it must be supplied.
!	IV : (input) The first dimension of AV as declared in the calling
!		program. IV.GE.NK+K-2
!	SIGMA : (input/output) Real array of length NK+K-2 containing the
!		singular values of the matrix A. If IFLG<4 this array will
!		be calculated, otherwise it must be supplied.
!	NSIM : (input) Number of sets to be tried for simulations to
!		calculate the error estimates. If NSIM.LE.1 error estimates
!		are not calculated.
!	FE : (output) Real array of length NR containing the estimated
!		error in F(I). This is calculated only if NSIM>1.
!	WK : Real array of length NR*NSIM+NM+NS+4NK+5K+2, used as scratch space.
!	
!
!	Required routines : BSPLIN, BSPEVL, SVD, SVDEVL, RANGAU
!	
!	
      SUBROUTINE RLS(NK,XO,K,NR,R,RKER,IK,AC,NM,NS,ALP,IDE,DI,DE,DF,F,
     1           B,IFLG,IER,REPS,CHISQ,SUMD,A,AV,IV,SIGMA,NSIM,FE,WK)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION XO(NK),R(NR),RKER(IK,NR),AC(IK,NK+K),DI(NM),DE(NM),
     1   DF(NM),F(NR),B(NK+K),A(IK,NK+K),AV(IV,NK+K),SIGMA(NK+K),D0(2),
     2   FE(NR),WK(NR,*)
 
      NV=NK+K-2
      NE=NM+NS
      IF(NM.LE.NV.OR.IK.LT.NE.OR.IV.LT.NV) THEN
        IER=709
        RETURN
      ENDIF
      IF(ALP.LT.0.0.OR.IDE.LT.1.OR.IDE.GT.2) THEN
        IER=710
        RETURN
      ENDIF
      IER=0
 
      IF(IFLG.LT.2) THEN
!     Setting up the system of linear equations
        NDERIV=0
        DO 3000 I=1,NM
 
          DO 2000 J=1,NV
            SIGMA(J)=0.0
2000      CONTINUE
          H=(R(2)-R(1))/2.
          DO 2500 IR=1,NR
            CALL BSPLIN(XO,NK,K,R(IR),NDERIV,DF,AV,AV(1,2),LEFT,IER,WK)
            IF(IER.GT.100) RETURN
 
            DO 2200 J=1,NV
              SIGMA(J)=SIGMA(J)+H*DF(J)*RKER(I,IR)
2200        CONTINUE
 
            IF(IR.LT.NR-1) THEN
              H=(R(IR+2)-R(IR))/2.0
            ELSE IF(IR.EQ.NR-1) THEN
              H=(R(IR+1)-R(IR))/2.0
            ENDIF
2500      CONTINUE
 
          DO 2800 J=1,NV
            A(I,J)=SIGMA(J)/DE(I)
            AC(I,J)=SIGMA(J)
2800      CONTINUE
3000    CONTINUE
 
      ELSE IF(IFLG.LT.4) THEN
!	The coefficients of matrix are available
        DO 3200 I=1,NM
          DO 3200 J=1,NV
            A(I,J)=AC(I,J)/DE(I)
3200    CONTINUE
 
      ENDIF
 
      IF(IFLG.LT.4) THEN
!     The equations arising from regularisation term
        H=(R(NR)-R(1))/(NS-1)
        NDERIV=IDE
        FA=ALP*SQRT(H)
 
        DO 4000 I=1,NS
          XI=R(1)+H*(I-1)
          CALL BSPLIN(XO,NK,K,XI,NDERIV,DF,AV,AV(1,2),LEFT,IER,WK)
          IF(IER.GT.100) RETURN
          DO 3500 J=1,NV
            A(I+NM,J)=AV(J,IDE)*FA
3500      CONTINUE
4000    CONTINUE
 
        CALL SVD(NV,NE,A,AV,SIGMA,IK,IV,DF,IER)
        IF(IER.GT.0) RETURN
      ENDIF
      IF(IFLG.EQ.1.OR.IFLG.EQ.3) THEN
        IFLG=4
        RETURN
      ENDIF
      IFLG=4
 
!     Set up the RHS of equations
      DO 4200 I=1,NM
        B(I)=DI(I)/DE(I)
4200  CONTINUE
      DO 4300 I=NM+1,NM+NS
        B(I)=0.0
4300  CONTINUE

!	Solve the system of equations using SVD
      CALL SVDEVL(NV,NE,A,AV,SIGMA,IK,IV,B,DF,REPS)
 
      NDERIV=0
      DO 4400 I=1,NR
        F(I)=BSPEVL(NK,XO,K,NDERIV,B,R(I),DF0,DDF,WK,IER)
4400  CONTINUE
 
!	Calculate the smoothing term
      SUMD=0.0
      NDERIV=IDE
      H=(R(NR)-R(1))/(NS-1)
      DO 4600 I=1,NS
        XI=R(1)+(I-1)*H
        FI=BSPEVL(NK,XO,K,NDERIV,B,XI,D0(1),D0(2),WK,IER)
        SUMD=SUMD+D0(IDE)**2
4600  CONTINUE
      SUMD=SUMD*H
 
!	Calculate the Chi square
      CHISQ=0.0
      DO 5000 I=1,NM
        S1=0.0
        DO 4800 J=1,NV
          S1=S1+AC(I,J)*B(J)
4800    CONTINUE
        DF(I)=(DI(I)-S1)/DE(I)
        CHISQ=CHISQ+DF(I)**2
5000  CONTINUE
 
!	Calculate the error estimates in the solution */
      IF(NSIM.LE.1) RETURN
!     Seed for random number generator
      SSD=123
      NDERIV=0
      DO 5100 I=1,NR
        FE(I)=0.0
5100  CONTINUE
 
      DO 6000 IS=1,NSIM
        DO 5200 I=1,NM
          WK(I,NSIM+1)=DI(I)/DE(I)+RANGAU(SSD)
5200    CONTINUE
        DO 5300 I=NM+1,NM+NS
          WK(I,NSIM+1)=0.0
5300    CONTINUE
 
        CALL SVDEVL(NV,NE,A,AV,SIGMA,IK,IV,WK(1,NSIM+1),DF,REPS)
 
        DO 5400 I=1,NR
          WK(I,IS)=BSPEVL(NK,XO,K,NDERIV,WK(1,NSIM+1),R(I),DF0,DDF,
     1            WK(NM+NS+2,NSIM+1),IER)
          FE(I)=FE(I)+WK(I,IS)
5400    CONTINUE
6000  CONTINUE
 
      DO 6500 I=1,NR
        A1=FE(I)/NSIM
        S1=0.0
        DO 6200 IS=1,NSIM
          S1=S1+(WK(I,IS)-A1)**2
6200    CONTINUE
        FE(I)=SQRT(S1/NSIM)
6500  CONTINUE
 
      END
