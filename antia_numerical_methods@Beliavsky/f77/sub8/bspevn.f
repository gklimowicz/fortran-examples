!	To calculate function value using B-spline expansion in n-dimensions
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
!	WK : Real scratch array of length (NXD+K)*(3N+1)+K+2
!	IWK : Integer scratch array of length 2N
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values of IER may be set by BSPLIN which is called
!
!	BSPEVN will give the value of the B-spline expansion at X0.
!	This routine does not calculate the derivatives, for that
!	use BSPEVN1 or BSPEVN2.
!
!	Required routines : BSPLIN
 
      FUNCTION BSPEVN(N,NK,X,NXD,K,WT,X0,WK,IWK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),X(NXD,N),WT(*),X0(N),WK(*),IWK(2*N)
 
      K1=(NXD+K)*3*N+1
      NDERIV=0
      BSPEVN=0.0
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
 
!	Calculate the summation over n dimensions
      F=0.0
2000  INDEX=IWK(1)+IWK(N+1)
      TERM=WK(INDEX)
      NDP=NK(1)+K-2
      DO 2200 I=2,N
        N1=(I-1)*3*(NXD+K)
        TERM=TERM*WK(N1+IWK(I)+IWK(N+I))
        INDEX=INDEX+(IWK(I)+IWK(I+N)-1)*NDP
        NDP=NDP*(NK(I)+K-2)
2200  CONTINUE
      F=F+TERM*WT(INDEX)
 
      J=1
!	Choose the next point
2400  IF(IWK(J+N).GE.K-1) GO TO 2600
      IWK(J+N)=IWK(J+N)+1
      GO TO 2000
 
!	If Jth dimension is exhausted go to next one
2600  IWK(J+N)=0
      J=J+1
      IF(J.LE.N) GO TO 2400
 
      BSPEVN=F
 
      END
