!	To calculate the integral of a function defined by B-spline expansion
!	over hyper-rectangular region in N dimensions
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the number of knots
!		for B-spline representation along each direction
!	X : (input) Real array of length NXD*N containing the knots along each direction
!		X(I,J) is the Ith knot along Jth dimension
!	NXD : (input) First dimension of array X as declared in the calling program
!		NXD must be greater than or equal to maximum of NK(I)
!       K : (input) Order of B-splines, K=4 for cubic B-splines
!	WT : (input) Array of length
!		(NK(1)+K-2)(NK(2)+K-2) ... (NK(N)+K-2)
!		containing coefficients of B-spline expansion. 
!		The coefficients are assumed to be stored in natural
!		Fortran order with no gaps in data. In the calling
!		program the array should have dimension
!		WT(NK(1)+K-2,NK(2)+K-2,...,NK(N)+K-2)
!	XL : (input) Real array of length N containing the lower limits of
!		integration along each direction
!	XU : (input) Real array of length N containing the upper limits of
!		integration along each direction
!	WK : Scratch array of length about twice the size of WT
!		The actual length required is complicated, but for N>2
!		the above limit should be good enough for all practical problems.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		Nonzero values may be set by BSPQD which is called
!
!	BSPQDN = Integral over[XL(1),XU(1)] x ... x [XL(N),XU(N)] of
!	\sum WT(I_1,...,I_N)\phi_I1(x_1) ... \phi_IN(x_N) 
!
!	Required routines : BSPLIN, BSPQD
!
 
      FUNCTION BSPQDN(N,NK,X,NXD,K,WT,XL,XU,WK,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(NXD,N),NK(N),WT(*),WK(*),XL(N),XU(N)
 
      N1=1
      DO 1000 I=2,N
        N1=N1*(NK(I)+K-2)
1000  CONTINUE
      N0=N1+4
      N2=N0+N1/(NK(2)+K-2)+4
      BSPQDN=0.0
 
!	Integration along the first dimension
      NK1=NK(1)+K-2
      DO 2000 J=1,N1
        NW=(J-1)*NK1+1
        WK(J)=BSPQD(NK(1),X(1,1),K,WT(NW),XL(1),XU(1),WK(N2),IER)
        IF(IER.GT.100) RETURN
2000  CONTINUE
 
      NA=0
      NB=N0
!	Integrate along the remaining dimensions
      DO 3000 I=2,N
        NK1=NK(I)+K-2
        N1=N1/NK1
        DO 2400 J=1,N1
        NW=(J-1)*NK1+1
        WK(NB+J)=BSPQD(NK(I),X(1,I),K,WK(NW+NA),XL(I),XU(I),WK(N2),IER)
        IF(IER.GT.100) RETURN
2400    CONTINUE
        NT=NA
        NA=NB
        NB=NT
3000  CONTINUE
 
      BSPQDN=WK(NA+1)
      END
