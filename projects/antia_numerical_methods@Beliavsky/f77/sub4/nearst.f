!	To locate the nearest point in an ordered table using bisection
!
!	XB : (input) given value of x for which nearest point is needed
!	X : (input) array of length NTAB containing table of values
!	NTAB : (input) length of table
!	After execution X(NEARST) is the tabular point closest to XB 
!
!	Required routines : None

      FUNCTION NEARST(XB,X,NTAB)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NTAB)

      LOW=1
      IGH=NTAB
      IF(.NOT.(XB.LT.X(LOW).EQV.XB.LT.X(IGH))) THEN

!	If the point is within the range of table, then locate it by bisection

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

      IF(ABS(XB-X(LOW)).LT.ABS(XB-X(IGH))) THEN
        NEARST=LOW
      ELSE
        NEARST=IGH
      ENDIF
      END
