!	To calculate rational function interpolation
!
!	XB : (input) x value at which interpolation is required
!	X : (input) Array of length NTAB containing the abscissas
!	F : (input) Array of length NTAB containing the function values at X(I)
!	NUSE : (input/output) Number of points to be used for interpolation
!		after execution NUSE will contain the number of points actually used
!	NTAB : (input) Number of points in the table
!	FB : (output) The interpolated value at x=XB
!	AEPS : (input) The required accuracy in interpolation
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=21 implies NUSE <2 in which case it is increased
!		IER=22 implies NUSE>NTAB or NMAX, in which case it is reduced
!		IER=23 implies interpolation did not converge to specified accuracy
!		IER=205 implies denominator turns out to be zero and
!			interpolation cannot proceed
!
!	Required routines : NEARST

      SUBROUTINE RATNAL(XB,X,F,NUSE,NTAB,FB,AEPS,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMAX=10)
      DIMENSION X(NTAB),F(NTAB),XN(NMAX),D(NMAX),C(NMAX)

!	To find the entry nearest to XB
      NEXT=NEARST(XB,X,NTAB)
      FB=F(NEXT)
      D(1)=F(NEXT)
      C(1)=D(1)
      XN(1)=X(NEXT)
      IER=0

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
!	If XB coincides with a tabular point, then RETURN
      IF(XB.EQ.XN(1)) RETURN

!	Calculate the successive rational function interpolation
      DO 5000 J=2,NIT

!	Choosing the next nearest point
        IF(IN.LE.1) GO TO 2200
        IF(IP.GE.NTAB) GO TO 2000
        IF(ABS(XB-X(IP+1)).LT.ABS(XB-X(IN-1))) GO TO 2200
2000    IN=IN-1
        NEXT=IN
        GO TO 2800
2200    IP=IP+1
        NEXT=IP

2800    XN(J)=X(NEXT)
        C(J)=F(NEXT)
        D(J)=F(NEXT)

!	Using the recurrences to calculate the differences C(I) and D(I)
        DO 3000 II=J-1,1,-1
          W=C(II+1)-D(II)
          RX=(XB-XN(II))/(XB-XN(J))
          FAC=(RX*D(II)-C(II+1))

          IF(FAC.EQ.0.0) THEN
            IER=205
            RETURN
          ENDIF

          FAC=W/FAC
          C(II)=RX*D(II)*FAC
          D(II)=C(II+1)*FAC
3000    CONTINUE

        FB=FB+C(1)
        NUSE=J
        IF(ABS(C(1)).LT.AEPS) RETURN
5000  CONTINUE

      IER=23
      RETURN
      END
