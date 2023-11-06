!	To calculate linear least squares fit to B-spline basis functions in N dimension
!	version for general weights, would require long time to calculate
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the number of
!		data points along each dimension
!	X : (input) Real array of length LA*N containing the coordinates
!		of points at which function values is available
!		X(I,J) contains the Ith point along Jth axis
!	F : (input) Real array of length NK(1)NK(2)...NK(N) containing the
!		function values. The dimension of F in the calling program
!		must match the size along each dimension, F(NK(1),...,NK(N))
!	EF : (input) Real array of length NK(1)NK(2)...NK(N) containing the
!		estimated errors in F
!	K : (input) Order of B-splines required, K=4 gives cubic B-splines
!	A : (output) Real array of length 
!		NK(1)NK(2)...NK(N) * (MK(1)+K-2)(MK(2)+K-2)...(MK(N)+K-2)
!		containing the matrix U of SVD of the design matrix
!		If RLM>0 the size should be (N+1) times this value
!	LA : (input) First dimension of arrays X as declared
!		in the calling program (LA .GE. MAX(NK(I)))
!	V : (output) Real array of length 
!		[(MK(1)+K-2)(MK(2)+K-2)...(MK(N)+K--2)]**2
!		containing the matrix V of SVD of the design matrix
!	IV : (input) First dimension of XF in the calling program
!		IV .GE. MAX(MK(I))
!	SIGMA : (output) Real array of length 
!		(MK(1)+K-2)(MK(2)+K-2)...(MK(N)+K-2)
!		containing the singular values of the design matrix
!	C : (output) Real array of length NK(1)NK(2)...NK(N) containing
!		the fitted coefficients. If RLM>0 then the required size
!		will be N+1 times this value. Note that although the number of
!		coefficients is (MK(1)+K-2)(MK(2)+K-2)...(MK(N)+K-2)
!		the rest of array is used as scratch space
!		Dimensions of array could be declared as
!		C(MK(1)+K-2, MK(2)+K-2, ..., MK(N-1)+K-2, NX)
!		where the last dimension is increased suitably to accommodate
!		the scratch space required. The last dimension NX should
!		be chosen such that the total size of array exceeds the
!		required value.
!	XF : (input) Real array of size IV*N, containing the knots
!		along each axis used for defining B-spline basis functions.
!		XF(I,J) should contain the Ith knot along Jth dimension
!	MK : (input) Integer array of length N containing the number of
!		knots for B-splines along each axis.
!	FY : (output) Real array of length NK(1)NK(2)...NK(N) containing
!		the values of fitted function at each of the tabular points
!		Dimensions of this array must match those of F.
!	WK : Real array of length 3N*(LA+6)+LA+2K+2 used as scratch space
!	IWK : Integer array of length 3N used as scratch space
!	REPS : (input) Required accuracy for solution of equations using SVD
!		singular values less than REPS times maximum will be set to zero
!	RLM : (input) Parameter lambda for smoothing. If RLM.LE.0 no smoothing
!		is applied
!	IDE : (input) Order of derivative to be used for smoothing
!		This is used only when RLM>0. IDE=1 for first derivative
!		and IDE=2 for second derivative smoothing
!	CHISQ : (output) The value of Chi square at minimum
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=608 implies that LA or IV are not large enough
!		IER=609 implies that RLM>0 and IDE is not acceptable
!		No calculations are done in these cases
!		Other values of IER may be set by SVD or BSPLIN
!
!	Required routines : BSPLIN, BSPEVN, SVD, SVDEVL
!
 
      SUBROUTINE BSPFITWN(N,NK,X,F,EF,K,A,LA,V,IV,SIGMA,C,XF,MK,FY,
     1      WK,IWK,REPS,RLM,IDE,CHISQ,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(LA,N),A(*),XF(IV,N),C(*),F(*),V(*),SIGMA(*),FY(*),
     1      WK(*),NK(N),MK(N),IWK(3*N),EF(*)
 
 
      IF(RLM.GT.0.0.AND.(IDE.LT.1.OR.IDE.GT.2)) THEN
        IER=609
        RETURN
      ENDIF

!     N1 is the number of equations to be solved
!     M is the number of coefficients to be determined
      N1=NK(1)
      M=MK(1)+K-2
      NMAX=N1
      MMAX=MK(1)+2*K+2
      DO 2000 I=2,N
        N1=N1*NK(I)
        M=M*(MK(I)+K-2)
        NMAX=MAX(NMAX,NK(I))
        MMAX=MAX(MMAX,MK(I)+2*K+2)
2000  CONTINUE
      IF(NMAX.GT.LA.OR.MMAX.GT.IV) THEN
        IER=608
        RETURN
      ENDIF
 
      N2=N1
      IF(RLM.GT.0.0) N1=(N+1)*N1
 
      NDB=LA+5
      NDERIV=0
      IF(RLM.GT.0.0) THEN
        NDERIV=IDE
        DO 2200 I=1,N
          IWK(I)=NDB+3*(I-1)*NDB
          IF(IDE.EQ.2) IWK(I)=IWK(I)+NDB
2200    CONTINUE
      ENDIF
 
!     Set up the matrix for equations
      DO 2500 I=1,N
        XB=X(1,I)
        N3=3*(I-1)*NDB+1
        CALL BSPLIN(XF(1,I),MK(I),K,XB,NDERIV,WK(N3),WK(N3+NDB),
     1              WK(N3+NDB*2),LEFT,IER,WK(N3+3*NDB))
        IF(IER.GT.100) RETURN
2500  CONTINUE
 
      DO 3000 I=1,N
        IWK(I+N)=1
3000  CONTINUE
      N0=3*N*NDB
 
!	Loop over each tabular point
3200  IJ=IWK(N+1)
      ID=NK(1)
      DO 3300 I=2,N
        IJ=IJ+(IWK(I+N)-1)*ID
        ID=ID*NK(I)
3300  CONTINUE
      C(IJ)=F(IJ)/EF(IJ)
 
!     set each coefficient in the equation
          DO 3510 I=1,N
            IWK(2*N+I)=1
3510      CONTINUE
 
!	Loop over the coefficients
3520      JI=IWK(2*N+1)
          DO 3530 I=1,N+1
            WK(N0+I)=WK(IWK(2*N+1))
3530      CONTINUE
          IF(RLM.GT.0.0) WK(N0+2)=WK(IWK(1)+IWK(2*N+1))

          ID=MK(1)+K-2
          DO 3550 I=2,N
            JI=JI+ID*(IWK(I+2*N)-1)
            ID=ID*(MK(I)+K-2)
            WK1=WK(3*(I-1)*NDB+IWK(I+2*N))
            WK(N0+1)=WK(N0+1)*WK1
            IF(RLM.GT.0.0) THEN
              DO 3540 J=1,N
                IF(I.NE.J) THEN
                  WK(N0+J+1)=WK(N0+J+1)*WK1
                ELSE
                  WK(N0+J+1)=WK(N0+J+1)*WK(IWK(I)+IWK(I+2*N))
                ENDIF
3540          CONTINUE
            ENDIF
3550      CONTINUE

!	setup the coefficient
          A(IJ+(JI-1)*N1)=WK(N0+1)/EF(IJ)
          IF(RLM.GT.0.0) THEN
            DO 3600 I=1,N
              A(IJ+I*N2+(JI-1)*N1)=RLM*WK(N0+1+I)
3600        CONTINUE
          ENDIF
 
!	Go to the next coefficient
          J1=1
3680      IF(IWK(2*N+J1).GE.MK(J1)+K-2) GO TO 3700
          IWK(2*N+J1)=IWK(2*N+J1)+1
          GO TO 3520
 
3700      IWK(2*N+J1)=1
          J1=J1+1
          IF(J1.LE.N) GO TO 3680

!	If all coefficients are done go to the next tabular point
      J=1
3800  IF(IWK(J+N).GE.NK(J)) GO TO 4000
      IWK(J+N)=IWK(J+N)+1
      N3=(3*J-3)*NDB+1
      XB=X(IWK(J+N),J)
      CALL BSPLIN(XF(1,J),MK(J),K,XB,NDERIV,WK(N3),WK(N3+NDB),
     1            WK(N3+NDB*2),LEFT,IER,WK(N0+3*NDB))
      IF(IER.GT.100) RETURN
      GO TO 3200
 
!	If Jth dimension is exhausted go to the next
4000  IWK(J+N)=1
      N3=(3*J-3)*NDB+1
      XB=X(IWK(J+N),J)
      CALL BSPLIN(XF(1,J),MK(J),K,XB,NDERIV,WK(N3),WK(N3+NDB),
     1            WK(N3+NDB*2),LEFT,IER,WK(N0+3*NDB))
      IF(IER.GT.100) RETURN
      J=J+1
      IF(J.LE.N) GO TO 3800
 
!	If all points are done find the SVD of the matrix
      JA=N1
      JV=M
      CALL SVD(M,N1,A,V,SIGMA,JA,JV,WK,IER)
      IF(IER.GT.100) RETURN
 
 
!     Setup the remaining part of RHS and solve the equations
      DO 4100 I=N2+1,N1
        C(I)=0.0
4100  CONTINUE
 
      CALL SVDEVL(M,N1,A,V,SIGMA,JA,JV,C,WK,REPS)
 
 
!     Calculate the \chi^2
      CHISQ=0.0
      DO 4200 I=1,N
        IWK(I)=1
4200  CONTINUE
      N0=N+1
 
!	loop over all points
4300  IJ=IWK(1)
      ID=NK(1)
      WK(1)=X(IJ,1)
      DO 4500 I=2,N
        IJ=IJ+ID*(IWK(I)-1)
        ID=ID*NK(I)
        WK(I)=X(IWK(I),I)
4500  CONTINUE
      FY(IJ)=BSPEVN(N,MK,XF,IV,K,C,WK,WK(N0),IWK(N0),IER1)
      CHISQ=CHISQ+((F(IJ)-FY(IJ))/EF(IJ))**2
 
!	Go to the next point
      J=1
4800  IF(IWK(J).GE.NK(J)) GO TO 5000
      IWK(J)=IWK(J)+1
      GO TO 4300
 
5000  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 4800
 
      END
