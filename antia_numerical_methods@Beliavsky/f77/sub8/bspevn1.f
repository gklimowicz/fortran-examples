!	To calculate function value and first derivative using B-spline
!		expansion in n-dimensions
!
!	N : (input) Number of dimensions
!	NK : (input) Array of length N containing the number of knots
!		along each dimension
!	X : (input) Real array of length NXD*N containing the knots
!		X(I,J) is the Ith knot along Jth dimension
!	NXD : (input) First dimension of array X in calling program.
!		NXD must be greater than or equal to maximum of NK(I).
!	K : (input) Order of B-splines, K=4 for cubic B-splines
!	WT : (input) Coefficients of B-spline expansion, the coefficients
!		are assumed to be stored in natural FORTRAN order with
!		no gaps in data. The calling program should have
!	 	dimension  WT(NK(1)+K-2,NK(2)+K-2,...,NK(N)+K-2)
!	X0 : (input) Real array of length N, containing the coordinates
!		of the point at which function needs to be evaluated
!	DF : (output) Real array of length N which will contain the
!		calculated values of first derivatives.
!		DF(I) will be derivative w.r.t. X_I
!	WK : Real scratch array of length (NXD+K)*(3N+1)+K+N+2
!	IWK : Integer scratch array of length 2N
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values of IER may be set by BSPLIN which is called
!
!	BSPEVN1 will give the value of the B-spline expansion at X0.
!	This routine does not calculate the second derivative, for that
!	use BSPEVN2. To improve efficiency if derivatives are not
!	required then use BSPEVN.
!
!	Required routines : BSPLIN
 
 
      FUNCTION BSPEVN1(N,NK,X,NXD,K,WT,X0,DF,WK,IWK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),X(NXD,N),WT(*),X0(N),DF(N),WK(*),IWK(2*N)
 
      K1=(NXD+K)*3*N+1
      NDERIV=1
      BSPEVN1=0.0
!	Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*(NXD+K)*(I-1)+1
        N2=(NXD+K)*(3*I-2)+1
        N3=(NXD+K)*(3*I-1)+1
        CALL BSPLIN(X(1,I),NK(I),K,X0(I),NDERIV,WK(N1),WK(N2),WK(N3),
     1              IWK(I),IER,WK(K1))
        IF(IER.GT.100) RETURN
        IWK(N+I)=0
1000  CONTINUE
 
!	calculate the sum over all points
      F=0.0
      DO 1800 I=1,N
        DF(I)=0.0
1800  CONTINUE
      N0=(NXD+K)*3*N

2000  INDEX=IWK(1)+IWK(N+1)
      TERM=WK(INDEX)
      WK(N0+1)=WK(NXD+K+INDEX)
      DO 2100 I=2,N
        WK(N0+I)=WK(INDEX)
2100  CONTINUE
      NDP=NK(1)+K-2
      DO 2300 I=2,N
        N1=(I-1)*3*(NXD+K)
        N2=IWK(I)+IWK(N+I)
        TERM=TERM*WK(N1+N2)

!	terms for first derivatives
        DO 2200 J=1,N
          IF(J.EQ.I) THEN
            WK(N0+J)=WK(N0+J)*WK(N1+N2+NXD+K)
          ELSE
            WK(N0+J)=WK(N0+J)*WK(N1+N2)
          ENDIF
2200    CONTINUE
        INDEX=INDEX+(N2-1)*NDP
        NDP=NDP*(NK(I)+K-2)
2300  CONTINUE
      F=F+TERM*WT(INDEX)
      DO 2400 I=1,N
        DF(I)=DF(I)+WK(N0+I)*WT(INDEX)
2400  CONTINUE
 
      J=1
!	Go to the next point
2500  IF(IWK(J+N).GE.K-1) GO TO 2600
      IWK(J+N)=IWK(J+N)+1
      GO TO 2000
 
!	If Jth dimension is exhausted go to the next one
2600  IWK(J+N)=0
      J=J+1
      IF(J.LE.N) GO TO 2500
 
      BSPEVN1=F
 
      END
