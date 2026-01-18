!	To calculate linear least squares fit to B-spline basis functions in N dimension
!	Weights are assumed to be unity
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
!	K : (input) Order of B-splines required, K=4 gives cubic B-splines
!	A : (output) Real array of length LA*IV*N containing the matrix
!		U of SVD of the design matrix for fit along each axis
!	LA : (input) First dimension of arrays X and A as declared
!		in the calling program (LA .GE. 2*MAX(NK(I)))
!	V : (output) Real array of length IV*IV*N containing the matrix
!		V of SVD of the design matrix for fit along each axis
!	IV : (input) First dimension of V, XF, SIGMA in the calling program
!		IV .GE. MAX(MK(I))+K-2
!	SIGMA : (output) Real array of length IV*N containing the singular
!		values of the design matrix for fit along each axis
!	C : (output) Real array of length 2NK(1)NK(2)...NK(N) containing
!		the fitted coefficients. Note that although the number of
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
!	WK : Real array of length 2NK(1)NK(2)...NK(N)+2LA used as scratch space
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
!		IER=609 implies that RLM>0 and IDE is not acceptable
!		No calculations are done in this case
!		Other values of IER may be set by SVD or BSPFIT
!
!	Required routines : BSPFIT, BSPLIN, BSPEVL, BSPEVN, SVD, SVDEVL
!
 
      SUBROUTINE BSPFITN(N,NK,X,F,K,A,LA,V,IV,SIGMA,C,XF,MK,FY,WK,IWK,
     1         REPS,RLM,IDE,CHISQ,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(LA,N),A(LA*IV,*),XF(IV,N),C(*),F(*),V(IV*IV,*),
     1         SIGMA(IV,N),FY(*),WK(*),NK(N),MK(N),IWK(3*N)
 
      IF(RLM.GT.0.0.AND.(IDE.LT.1.OR.IDE.GT.2)) THEN
        IER=609
        RETURN
      ENDIF
 
      N1=4*LA+K+12
      N2=N1+LA+2
!	Set the weights to one for fits along each dimension
      DO 1000 I=0,LA
        WK(N1+I)=1.0
1000  CONTINUE
!	Calculate the SVD of matrices for fit along each dimension
      DO 2000 I=1,N
        IFLG1=1
        CALL BSPFIT(NK(I),X(1,I),F,WK(N1),K,A(1,I),LA,V(1,I),IV,
     1  SIGMA(1,I),C,XF(1,I),MK(I),FY,IFLG1,WK,REPS,RLM,IDE,CHISQ,
     2        WK(N2),IER)
        IF(IER.GT.100) RETURN
2000  CONTINUE
 
      NX=1
      DO 2200 I=1,N-1
        NX=NX*NK(I)
2200  CONTINUE
      N0=NX*NK(N)+1
      IF(RLM.GT.0.0) N0=2*NX*NK(N)+1
      NY=1
      JU=N
 
      IF(MOD(N,2).EQ.1) THEN
!	If N is odd then fit along the last dimension outside the loop
        N1=NK(N)
        IF(RLM.GT.0.0) N1=2*NK(N)
!	Set up the RHS for fit
        LJ=N1
        DO 2500 I=1,NX
          DO 2500 J=1,NK(N)
            WK(J+(I-1)*LJ)=F(I+(J-1)*NX)
            IF(RLM.GT.0.0) WK(J+NK(N)+(I-1)*LJ)=0.0
2500    CONTINUE
        M=MK(N)+K-2
        DO 2800 I=1,NX
          NI=1+(I-1)*LJ
          CALL SVDEVL(M,N1,A(1,N),V(1,N),SIGMA(1,N),LA,IV,WK(NI),WK(N0)
     1               ,REPS)
2800    CONTINUE
 
!	setup the RHS for fit along the N-1 th dimension
        NX1=NX/NK(N-1)
        NY=M
        LJ1=NK(N-1)
        IF(RLM.GT.0.0) LJ1=2*NK(N-1)
        DO 3000 I1=1,NK(N-1)
          DO 3000 I=1,NX1
            DO 3000 J=1,NY
              C(I1+(I-1)*LJ1+(J-1)*NX1*LJ1)=WK(J+(I-1)*LJ+(I1-1)*NX1*LJ)
              IF(RLM.GT.0.0) C(I1+NK(N-1)+(I-1)*LJ1+(J-1)*NX1*LJ1)=0.0
3000    CONTINUE
        NX=NX1
        JU=N-1
      ELSE
!	setup the RHS for fit along the N th dimension
        LJ=NK(N)
        IF(RLM.GT.0.0) LJ=2*NK(N)
        DO 3200 I=1,NX
          DO 3200 J=1,NK(N)
            C(J+(I-1)*LJ)=F(I+(J-1)*NX)
            IF(RLM.GT.0.0) C(J+NK(N)+(I-1)*LJ)=0.0
3200    CONTINUE
 
      ENDIF
 
!	Loop for fit along each dimension, each pass fits along 2 dimensions
      DO 5000 J1=JU,1,-2
        NUM=NX*NY
        LJ=NK(J1)
        IF(RLM.GT.0.0) LJ=2*NK(J1)
        N1=LJ
        M=MK(J1)+K-2
        DO 3400 I=1,NUM
          NI=1+(I-1)*LJ
          CALL SVDEVL(M,N1,A(1,J1),V(1,J1),SIGMA(1,J1),LA,IV,C(NI),
     1                WK(N0),REPS)
3400    CONTINUE
 
!	Set up the RHS for fit along next dimension
        NX1=NX/NK(J1-1)
        NY1=NY*M
        LJ1=NK(J1-1)
        IF(RLM.GT.0.0) LJ1=2*NK(J1-1)
        DO 3500 I1=1,NY
          DO 3500 I2=1,M
            DO 3500 I=1,NK(J1-1)
              DO 3500 J=1,NX1
                WK(I+(J-1)*LJ1+(I2-1)*LJ1*NX1+(I1-1)*NX1*LJ1*M)=
     1                      C(I2+(J-1)*LJ+(I-1)*NX1*LJ+(I1-1)*NX*LJ)
                IF(RLM.GT.0.0) WK(I+NK(J1-1)+(J-1)*LJ1+(I2-1)*LJ1*NX1+
     1                             (I1-1)*NX1*LJ1*M)=0.0
3500    CONTINUE
        NX=NX1
        NY=NY1
        NUM=NY*NX
        LJ=NK(J1-1)
        IF(RLM.GT.0.0) LJ=2*NK(J1-1)
        M=MK(J1-1)+K-2
        N1=LJ
        DO 3600 I=1,NUM
          NI=1+(I-1)*LJ
          CALL SVDEVL(M,N1,A(1,J1-1),V(1,J1-1),SIGMA(1,J1-1),LA,IV,
     1             WK(NI),WK(N0),REPS)
3600    CONTINUE
 
        IF(J1.EQ.2) THEN
!	Store the fitted coefficients in array C
          N1=NK(1)
          IF(RLM.GT.0.0) N1=2*NK(1)
          DO 3800 I=1,M
            DO 3800 J=1,NY
              C(I+(J-1)*M)=WK(I+(J-1)*LJ)
3800      CONTINUE
        ELSE
 
!	Set up the RHS for fit along next dimension
          LJ1=NK(J1-2)
          IF(RLM.GT.0.0) LJ1=2*NK(J1-2)
          M=MK(J1-1)+K-2
          NX1=NX/NK(J1-2)
          NY1=NY*M
          DO 4000 I1=1,M
            DO 4000 I2=1,NK(J1-2)
              DO 4000 I=1,NX1
                DO 4000 J=1,NY
                  C(I2+(I-1)*LJ1+(I1-1)*LJ1*NX1+(J-1)*LJ1*NX1*M)=
     1                     WK(I1+(I-1)*LJ+(I2-1)*LJ*NX1+(J-1)*LJ*NX)
                  IF(RLM.GT.0.0) C(I2+NK(J1-2)+(I-1)*LJ1+
     1                         (I1-1)*LJ1*NX1+(J-1)*LJ1*NX1*M)=0.0
4000      CONTINUE
          NX=NX1
          NY=NY1
        ENDIF
5000  CONTINUE
 
!	Calculate the Chi Square
      CHISQ=0.0
      DO 6000 I=1,N
        IWK(I)=1
6000  CONTINUE
      N0=N+1
 
!	Loop over all points
6200  IJ=IWK(1)
      WK(1)=X(IJ,1)
      ID=NK(1)
      DO 6500 I=2,N
        IJ=IJ+ID*(IWK(I)-1)
        ID=ID*NK(I)
        WK(I)=X(IWK(I),I)
6500  CONTINUE

      FY(IJ)=BSPEVN(N,MK,XF,IV,K,C,WK,WK(N0),IWK(N0),IER1)
      CHISQ=CHISQ+(F(IJ)-FY(IJ))**2
 
!	Choose the next point
      J=1
6800  IF(IWK(J).GE.NK(J)) GO TO 7000
      IWK(J)=IWK(J)+1
      GO TO 6200
 
!	If Jth dimension is exhausted go to next one
7000  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 6800
 
      END
