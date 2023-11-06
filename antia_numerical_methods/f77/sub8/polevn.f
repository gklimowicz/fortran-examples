!	Evaluating the fitted polynomial at any point using known
!	coefficients of orthogonal polynomials in N dimensions.
!	Should be used to evaluate the polynomial using coefficients calculated
!	by POLFITN. It does not calculate the derivatives, for derivatives
!	use POLEVN1 or POLEVN2
!
!	N : (input) Number of dimensions
!	NK : (input) Integer array of length N containing the degree of
!		polynomial in each direction
!	AX : (input) Real array of length LA*(3*N+3) containing the coefficients
!		alpha and beta for orthogonal polynomials in each direction
!		AX(I,3*J-2) contains alpha and AX(I,3*J-1) contains beta
!		for polynomials in Jth dimension
!	LA : (input) First dimension of array AX in the calling program.
!		It must be same as what was used in call to POLFITN while
!		calculating the coefficients.
!	WT : (input) Real array of length (MK(1)+1)(MK(2)+1)...(MK(N)+1)
!		containing the coefficients of the fit. The dimensions of WT in
!		the calling program must match the size along each dimension,
!		WT(MK(1)+1,MK(2)+1,...,MK(N)+1)
!	X0 : (input) Real array of length N containing the coordinates of
!		the point at which polynomial needs to be evaluated
!	F : (output) Calculated value of the fitted polynomial at X0
!	WK : Real array of length 3*N*LA used as scratch space
!	IWK : Integer array of length N used as scratch space
!	
!	Required routines : POLORT
!
      SUBROUTINE POLEVN(N,NK,AX,LA,WT,X0,F,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),AX(LA,3*N),WT(*),X0(N),WK(LA*3*N),IWK(N)
 
!     Calculate the orthogonal polynomials along each dimension
      DO 1000 I=1,N
        N1=3*LA*(I-1)+1
        N2=3*LA*(I-1)+LA+1
        N3=3*LA*(I-1)+2*LA+1
        NJ=3*(I-1)+1
        XB=X0(I)
        CALL POLORT(NK(I),AX(1,NJ),AX(1,NJ+1),XB,WK(N1),WK(N2),WK(N3))
        IWK(I)=1
1000  CONTINUE
 
!     Calculate the summation over n dimensions
      F=0.0
2000  INDEX=IWK(1)
      TERM=WK(INDEX)
      NDP=NK(1)+1
      DO 2200 I=2,N
        N1=(I-1)*3*LA
        TERM=TERM*WK(N1+IWK(I))
        INDEX=INDEX+(IWK(I)-1)*NDP
        NDP=NDP*(NK(I)+1)
2200  CONTINUE
      F=F+TERM*WT(INDEX)
 
      J=1
!     Choose the next point
2400  IF(IWK(J).GE.NK(J)+1) GO TO 2600
      IWK(J)=IWK(J)+1
      GO TO 2000
 
!     If Jth dimension is exhausted go to next one
2600  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 2400
 
      END
