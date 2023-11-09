!	Evaluating the fitted polynomial and its derivatives at any point
!	using known coefficients of orthogonal polynomials in N dimensions.
!	Should be used to evaluate the polynomial using coefficients calculated
!	by POLFITN. If second derivative is not required use POLEVN1, if
!	no derivatives are required then use POLEVN
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
!	DF : (output) Real array of length N containing the first derivatives
!		of F at X0
!	DDF : (output) Real array of length N*N containing the second derivatives
!		of F at X0, DDF(I,J)=d^2F/(dX(I) dX(J))
!	WK : Real array of length 3*N*LA+N+N*N used as scratch space
!	IWK : Integer array of length N used as scratch space
!	
!	Required routines : POLORT
!
      SUBROUTINE POLEVN2(N,NK,AX,LA,WT,X0,F,DF,DDF,WK,IWK)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NK(N),AX(LA,3*N),WT(*),X0(N),WK(*),IWK(N),DF(N),DDF(N,N)
 
!     Calculate the B-spline basis functions along each dimension
      DO 1000 I=1,N
        N1=3*LA*(I-1)+1
        N2=3*LA*(I-1)+LA+1
        N3=3*LA*(I-1)+2*LA+1
        NJ=3*(I-1)+1
        XB=X0(I)
        CALL POLORT(NK(I),AX(1,NJ),AX(1,NJ+1),XB,WK(N1),WK(N2),WK(N3))
        IWK(I)=1
1000  CONTINUE
 
!     calculate the sum over all points
 
      F=0.0
      DO 1800 I=1,N
        DF(I)=0.0
        DO 1800 J=1,N
          DDF(J,I)=0.0
1800  CONTINUE
      N0=LA*3*N
2000  INDEX=IWK(1)
      TERM=WK(INDEX)

!	Terms for the derivatives
      WK(N0+1)=WK(LA+INDEX)
      DO 2100 I=2,N+N*N
        WK(N0+I)=WK(INDEX)
2100  CONTINUE
      WK(N0+N+1)=WK(2*LA+INDEX)
      DO 2200 I=2,N
        WK(N0+N+I)=WK(INDEX+LA)
2200  CONTINUE
      NDP=NK(1)+1
      DO 2500 I=2,N
        N1=(I-1)*3*LA
        N2=IWK(I)
        TERM=TERM*WK(N1+N2)
        DO 2300 J=1,N
          IF(J.EQ.I) THEN
            WK(N0+J)=WK(N0+J)*WK(N1+N2+LA)
          ELSE
            WK(N0+J)=WK(N0+J)*WK(N1+N2)
          ENDIF
          DO 2300 J1=1,J
            IF(J.EQ.I.AND.J1.EQ.I) THEN
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2+2*LA)
            ELSE IF(J.EQ.I.OR.J1.EQ.I) THEN
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2+LA)
            ELSE
              WK(N0+N+(J1-1)*N+J)=WK(N0+N+(J1-1)*N+J)*WK(N1+N2)
            ENDIF
2300    CONTINUE
        INDEX=INDEX+(N2-1)*NDP
        NDP=NDP*(NK(I)+1)
2500  CONTINUE
      F=F+TERM*WT(INDEX)
      DO 2600 I=1,N
        DF(I)=DF(I)+WK(N0+I)*WT(INDEX)
        DO 2600 J=1,I
          DDF(I,J)=DDF(I,J)+WK(N0+N+(J-1)*N+I)*WT(INDEX)
2600  CONTINUE
 
      J=1
!	Go to the next point
3000  IF(IWK(J).GE.NK(J)+1) GO TO 3200
      IWK(J)=IWK(J)+1
      GO TO 2000
 
!	If Jth dimension is exhausted, go to the next one
3200  IWK(J)=1
      J=J+1
      IF(J.LE.N) GO TO 3000
 
      DO 3500 I=1,N
        DO 3500 J=I+1,N
          DDF(I,J)=DDF(J,I)
3500  CONTINUE
 
      END
