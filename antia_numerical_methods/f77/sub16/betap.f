!	To calculate the incomplete Beta function I_x(a,b) for
!	positive real arguments
!	It returns a value of -1 for a.le.0 or b.le.0 or x<0 or x>1
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!                       Beta function
!
!	Required routines : GAMMAL, BETAI, BETSER, BETCON, BETCON1,
!			    ADPINT, KRONRD, FBETA

      FUNCTION BETAP(A,B,X)
      IMPLICIT REAL*16(A-H,O-Z)

51    FORMAT(' ** Error in evaluating Incomplete Beta Function at',/,
     1  ' A=',1PD12.4,' B=',D12.4,' X=',D12.4,'  BETAP=',D12.4)

      IF(A.LE.0.OR.B.LE.0.OR.X.LT.0.OR.X.GT.1) THEN
        BETAP=-1
        RETURN
      ENDIF
      AMAX=MAX(A,B)
      AMIN=MIN(A,B)
      BETAL=(A+B)*LOG(A+B)-A*LOG(A)-B*LOG(B)
      IF(AMAX.LE.30) THEN
        IF(X.LE.0.5) THEN
          BETAP=BETSER(A,B,X)
        ELSE
          BETAP=BETSER(B,A,1-X)
          BETAP=1-BETAP
        ENDIF
      ELSE IF(B.LE.20.AND.X.LE.0.71Q0) THEN
        BETAP=BETSER(A,B,X)
      ELSE IF(A.LE.20.AND.X.GE.0.3Q0) THEN
        BETAP=BETSER(B,A,1-X)
        BETAP=1-BETAP
      ELSE IF(B.LE.50.AND.X.LE.0.35Q0) THEN
        BETAP=BETSER(A,B,X)
      ELSE IF(A.LE.50.AND.X.GE.0.65Q0) THEN
        BETAP=BETSER(B,A,1-X)
        BETAP=1-BETAP
      ELSE IF(B.LE.100.AND.X.LE.0.18Q0) THEN
        BETAP=BETSER(A,B,X)
      ELSE IF(A.LE.100.AND.X.GE.0.82Q0) THEN
        BETAP=BETSER(B,A,1-X)
        BETAP=1-BETAP
      ELSE IF(B.LE.180.AND.X.LE.0.1Q0) THEN
        BETAP=BETSER(A,B,X)
      ELSE IF(A.LE.180.AND.X.GE.0.9Q0) THEN
        BETAP=BETSER(B,A,1-X)
        BETAP=1-BETAP
      ELSE IF(X.LT.0.5) THEN
        IF(A.LT.2.Q0) THEN
          BETAP=BETCON(A,B,X)
        ELSE IF(BETAL.GT.700) THEN
          BETAP=BETCON1(A,B,X)
        ELSE
          BETAP=BETAI(A,B,X)
        ENDIF
      ELSE
        IF(B.LT.2.Q0) THEN
          BETAP=BETCON(B,A,1-X)
        ELSE IF(BETAL.GT.700) THEN
          BETAP=BETCON1(B,A,1-X)
        ELSE
          BETAP=BETAI(B,A,1-X)
        ENDIF
        BETAP=1-BETAP
      ENDIF
      
      IF(BETAP.LT.0.0.OR.BETAP.GT.1) PRINT 51,A,B,X,BETAP
      END
