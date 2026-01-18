!	Linear interpolation in N dimensions
!
!	N : (input) Number of variables
!	XB : (input) Real array of length N containing the coordinates
!		of the point at which interpolation is required
!	X : (input) Real array of length NXD*N containing the abscissas
!	F : (input) Real array of dimension exactly F(NDIM(1),...,NDIM(N))
!		containing the function values
!		F(I1,I2,...,IN)=f(X(I1,1),X(I2,2),...,X(IN,N))
!	NP : (input) Integer array of length N containing the number of tabular points
!		NP(I) is the number of tabular point along Ith dimension
!	FB : (output) Real variable giving the interpolated value
!	NDIM : (input) Integer array of length N, containing the
!		dimension of F as declared in calling program
!	NXD : (input) First dimension of array X, NXD .GE. Max(NP(i))
!	IN : Integer array of length 2N used as scratch space
!	HN : Real array of length N used as scratch space
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=207 implies NP(I)>NDIM(I) or NP(I)<2 for some I
!
!	Required routines : LOCATE

      SUBROUTINE LINRN(N,XB,X,F,NP,FB,NDIM,NXD,IN,HN,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(*),XB(N),F(*),NP(N),IN(2*N),HN(N),NDIM(N)

      IER=0
      ND=0
!	Locate the point within a hypercube in the table
      DO 2000 I=1,N
        IF(NP(I).GT.NDIM(I).OR.NP(I).GT.NXD.OR.NP(I).LT.2) THEN
          IER=207
          RETURN
        ENDIF

        LOW=LOCATE(XB(I),X(ND+1),NP(I))
        IN(I+N)=LOW
        HN(I)=(XB(I)-X(LOW+ND))/(X(LOW+ND+1)-X(LOW+ND))
        ND=ND+NXD
        IN(I)=0
2000  CONTINUE

      FB=0.0
3000  INDEX=IN(1)+IN(1+N)
      NDP=NDIM(1)
      DO 3200 I=2,N
        INDEX=INDEX+(IN(I)+IN(I+N)-1)*NDP
        NDP=NDP*NDIM(I)
3200  CONTINUE
!	F(IN(1)+IN(1+N), IN(2)+IN(2+N), ... ,IN(N)+IN(2N))
      TERM=F(INDEX)

      DO 3400 I=1,N
        IF(IN(I).EQ.0) THEN
          TERM=TERM*(1.-HN(I))
        ELSE
          TERM=TERM*HN(I)
        ENDIF
3400  CONTINUE

      FB=FB+TERM

!	select the next point
      K=1
4200  IF(IN(K).GE.1) GO TO 4400
      IN(K)=IN(K)+1
      GO TO 3000

!	If the Kth dimension is exhausted, go to the next one
4400  IN(K)=0
      K=K+1
      IF(K.LE.N) GO TO 4200

      END
