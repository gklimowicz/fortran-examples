!	To calculate polynomial interpolation in 2 dimensions on a
!	rectangular mesh of points
!
!	(XB1,XB2) : (input) is the point at which interpolation is required
!	X1 : (input) Real array of length N1 containing the abscissas
!	X2 : (input) Real array of length N2 containing the abscissas
!	F : (input) Real array of length NDIM*N2 containing the function values
!		F(I,J)=f(X1(I),X2(J))
!	NDIM : (input) First dimension of array F as specified in calling program
!	N1 : (input) Length of array X1, i.e Number of points along first dimension
!	N2 : (input) Length of array X2, i.e Number of points along second dimension
!	NP1 : (input) Number of points to be used for interpolation along X1
!	NP2 : (input) Number of points to be used for interpolation along X2
!	FB : (output) Interpolated value of function at (XB1,XB2)
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=206 implies N1 > NDIM, in which case no calculations
!			are done
!
!	Required routines : DIVDIF0, NEARST

      SUBROUTINE POLY2(XB1,XB2,X1,X2,F,NDIM,N1,N2,NP1,NP2,FB,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=10)
      DIMENSION X1(N1),X2(N2),F(NDIM,N2),XN(NMAX),XD(NMAX),FB1(NMAX)

      IF(N1.GT.NDIM) THEN
        IER=206
        RETURN
      ENDIF

      IER=0
      REPS=0.0
!	Find the point nearest to XB2 in X2
      NEXT=NEARST(XB2,X2,N2)

      NUSE=NP1
      IFLG=0
      CALL DIVDIF0(XB1,X1,F(1,NEXT),NUSE,N1,FB1,REPS,IER1,IFLG,IF1)
!	Set IFLG=1 so that next time DIVDIF0 does not try to locate
!	the point again in the table.
      IFLG=1
      FB=FB1(NUSE)
      XD(1)=FB
      XN(1)=X2(NEXT)
      PX=1.0

      IP=NEXT
      IN=NEXT

!	The number of points to be used along X2
      NIT=MIN(NMAX,NP2,N2)
      IF(NP2.LT.1) NIT=MIN(4,N2,NMAX)

!	Calculate the successive interpolation polynomials
      DO 5000 J=2,NIT

!	Find the next nearest point in X2
        IF(IN.LE.1) GO TO 2200
        IF(IP.GE.N2) GO TO 2000
        IF(ABS(XB2-X2(IP+1)).LT.ABS(XB2-X2(IN-1))) GO TO 2200
2000    IN=IN-1
        NEXT=IN
        GO TO 2800
2200    IP=IP+1
        NEXT=IP

2800    NUSE=NP1
!	interpolate along X1 to calculate function value at (XB1,X2(NEXT))
        CALL DIVDIF0(XB1,X1,F(1,NEXT),NUSE,N1,FB1,REPS,IER1,IFLG,IF1)
        XD(J)=FB1(NUSE)
        XN(J)=X2(NEXT)
!	Calculate the divided difference for interpolation in X2
        DO 3000 K=J-1,1,-1
3000    XD(K)=(XD(K+1)-XD(K))/(XN(J)-XN(K))

        PX=PX*(XB2-XN(J-1))
        ERR=XD(1)*PX
        FB=FB+ERR

5000  CONTINUE
      END
