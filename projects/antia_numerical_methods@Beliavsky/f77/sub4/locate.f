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
!      IMPLICIT REAL*8(A-H,O-Z)
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
