!	Interpolation using Newton's divided difference formula
!	simplified version of DIVDIF without derivative calculation
!
!	XB : (input) value of x at which interpolation is required
!	X : (input) real array of length NTAB containing x values
!	F : (input) real array of length NTAB containing function values
!		F(I) is the tabulated function value at X(I).
!	NUSE : (input/output) Number of points to be used for interpolation
!		After execution it will contain the number actually used
!	NTAB : (input) Number of points in the table
!	FB : (output) Real array containing interpolated values
!	       FB(I) should contain interpolation using I points
!	       FB(NUSE) should be the final value
!	AEPS : (input) Required accuracy
!	IER : Error parameter, IER=0 if the execution is successful
!		IER=21 implies NUSE<1, in which case it is set to MIN(6,NTAB)
!		IER=22 implies NUSE>NTAB or NMAX, in which case it is reduced
!		IER=23 implies interpolation has not converged to specified accuracy
!	IFLG : (input) Flag to decide whether nearest point has to be found
!		If IFLG=0 find the nearest point to XB to start interpolation
!		otherwise if IF1 is admissible use X(IF1) as the first point
!	IF1 : (input/output) The first point to be used for interpolation
!		when IFLG.NE.0
!		If IFLG=0 then IF1 is set to the index of nearest point in X
!
!	Required routines : NEARST

      SUBROUTINE DIVDIF0(XB,X,F,NUSE,NTAB,FB,AEPS,IER,IFLG,IF1)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=10)
      DIMENSION X(NTAB),F(NTAB),FB(*),XN(NMAX),XD(NMAX)

!	Choose the first point
      IF(IFLG.EQ.0.OR.IF1.LT.1.OR.IF1.GT.NTAB) THEN
        NEXT=NEARST(XB,X,NTAB)
        IF1=NEXT
      ELSE
        NEXT=IF1
      ENDIF
      FB(1)=F(NEXT)
      XD(1)=F(NEXT)
      XN(1)=X(NEXT)
      IER=0
      PX=1.0

!	Points between IN and IP are used for interpolation
      IP=NEXT
      IN=NEXT

!	Maximum number of points to be used for interpolation
      NIT=MIN(NMAX,NUSE,NTAB)
      IF(NUSE.GT.NMAX.OR.NUSE.GT.NTAB) IER=22
      IF(NUSE.LT.1) THEN
        IER=21
        NIT=MIN(6,NTAB,NMAX)
      ENDIF
      NUSE=1

!	Calculate successive interpolation polynomial
      DO 5000 J=2,NIT

!	Choose the next nearest point to XB
        IF(IN.LE.1) GO TO 2200
        IF(IP.GE.NTAB) GO TO 2000
        IF(ABS(XB-X(IP+1)).LT.ABS(XB-X(IN-1))) GO TO 2200
2000    IN=IN-1
        NEXT=IN
        GO TO 2800
2200    IP=IP+1
        NEXT=IP

!	Calculating the divided differences
2800    XD(J)=F(NEXT)
        XN(J)=X(NEXT)
        DO 3000 K=J-1,1,-1
3000    XD(K)=(XD(K+1)-XD(K))/(XN(J)-XN(K))

        PX=PX*(XB-XN(J-1))
        ERR=XD(1)*PX
        FB(J)=FB(J-1)+ERR
        NUSE=J

        IF(ABS(ERR).LT.AEPS) RETURN
5000  CONTINUE

      IER=23
      END
