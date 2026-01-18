!	To calculate the incomplete Beta function I_x(a,b) using
!	the infinite series, called by BETAP
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!                       Beta function
!
!	Required routines : GAMMAL

      FUNCTION BETSER(A,B,X)
      IMPLICIT REAL*16(A-H,O-Z)

51    FORMAT(' ** Roundoff error while evaluating the infinite series',
     1  ' for incomplete Beta function at',/,'  A=',1PD12.4,'  B=',
     2  D12.4,'  X=',D12.4,'  SUM=',D12.4,'  MAX TERM=',D12.4)

      S3=1
      T=A
      TMAX=A
      DO I=1,500,2
        T1=T*X*(I-B)/I
        T=T1*X*(I+1-B)/(I+1)
        S3=S3+(T1/(A+I)+T/(A+I+1))
        IF(ABS(T/(A+I+1)).GT.TMAX) TMAX=ABS(T/(A+I+1))
        IF(ABS(T/(A+I+1)).LT.1.Q-15) EXIT
      ENDDO
      IF((S3).LT.1.Q-16*TMAX) PRINT 51,A,B,X,S3,TMAX
      B1=A*LOG(X)+LOG(S3)-LOG(A)
      B1=B1+GAMMAL(A+B)-GAMMAL(A)-GAMMAL(B)
      BETSER=EXP(B1)
      END
