!     PROGRAM TO PERFORM LINEAR INTERPOLATION IN N DIMENSIONS

      PROGRAM LININT
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(8,8,8,8),X1(20,9),XB(10),N1(10),IWK(20),WK(20)
     1          ,NDIM(20)

!     SAMPLE LINEAR FUNCTION FOR GENERATING TABLE OF VALUES.
!     INTERPOLATION SHOULD GIVE EXACT VALUE APART FROM ROUNDOFF

      FN(X,Y,Z,V)=(1.+X)*(2.+Y)*(3.+Z)*(4.+V)
 
51    FORMAT('  F(',3(1PD14.6,','),D14.6,') =',D14.6/
     1       8X,'  IER =',I4,5X,'EXACT VALUE =',D14.6)

!     GENERATING THE TABLE OF VALUES IN FOUR DIMENSIONS
!     WITH UNIFORM SPACING OF H=0.1 AND CONTAINING N=5 POINTS
!     ALONG EACH AXES SPANNING THE INTERVAL [0.1,0.5]

      H=0.1D0
      N=5
      M=4
      DO 1000 J=1,M
        N1(J)=N
        NDIM(J)=8
        DO 1000 I=1,N
1000  X1(I,J)=I*H
      DO 1400 I=1,N
        DO 1400 J=1,N
          DO 1400 K=1,N
            DO 1400 L=1,N
              F(L,K,J,I)=FN(X1(L,1),X1(K,2),X1(J,3),X1(I,4))
1400  CONTINUE

      NXD=20 

2000  PRINT *,'TYPE XB(I), I=1,4  THE COORDINATES OF THE REQUIRED POINT'
      PRINT *,'               (QUITS WHEN XB(1)<-100)'
      READ *,(XB(I),I=1,M)
      IF(XB(1).LT.-100) STOP
      CALL LINRN(M,XB,X1,F,N1,FB,NDIM,NXD,IWK,WK,IER)
      WRITE(6,51) (XB(I),I=1,M),FB,IER,FN(XB(1),XB(2),XB(3),XB(4))
      GO TO 2000
      END
 
!     --------------------------------------------------
 
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
      IMPLICIT REAL*8(A-H,O-Z)
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
 
!     ------------------------------------------------------
 
!	To locate a given point between two points of an ordered table
!
!	XB : (input) The point which needs to be located
!	X : (input) Real array of length NP containing the ordered table
!	NP : (input) length of table
!	LOCATE should give the value such that XB is between X(LOCATE) and 
!		X(LOCATE+1)
!
!	Required routines : None

      FUNCTION LOCATE(XB,X,NP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NP)

      LOW=1
      IGH=NP
      IF(XB.LT.X(LOW).EQV.XB.LT.X(IGH)) THEN
        IF(ABS(XB-X(LOW)).GT.ABS(XB-X(IGH))) LOW=IGH-1
      ELSE

!	If the point is within the range of table locate it by bisection
1500    IF(IGH-LOW.GT.1) THEN
          MID=(LOW+IGH)/2
          IF(XB.LT.X(MID).EQV.XB.LT.X(LOW)) THEN
            LOW=MID
          ELSE
            IGH=MID
          ENDIF
          GO TO 1500
        ENDIF
      ENDIF

      LOCATE=LOW
      END
