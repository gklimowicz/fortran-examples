!	To calculate the incomplete Beta function I_x(a,b) using
!	the continued fraction, called by BETAP
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!                       Beta function
!
!	Required routines : GAMMAL

      FUNCTION BETCON1(A,B,X)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION c(500),d(500)
      
51    FORMAT(' ** Roundoff error while evaluating the ',
     1  'continued fraction for incomplete Beta function at',/,'  A=',
     2  1PE12.4,'  B=',E12.4,'  X=',E12.4,'  CON. FRAC.=',E12.4)

      C(1)=1
      D(1)=1
      C(2)=1
      D(2)=1-X*(A+B)/(A+1)
      C1=C(2)/D(2)
      DO I=2,498,2
        M=I/2
        C(I+1)=C(I)+M*(B-M)*X*C(I-1)/((A+I-1)*(A+I))
        D(I+1)=D(I)+M*(B-M)*X*D(I-1)/((A+I-1)*(A+I))
        C(I+2)=C(I+1)-(A+M)*(A+B+M)*X*C(I)/((A+I+1)*(A+I))
        D(I+2)=D(I+1)-(A+M)*(A+B+M)*X*D(I)/((A+I+1)*(A+I))
!	Scale the numerator and denominator to prevent underflow/overflow
        IF(ABS(C(I+2)).GT.1.D30) THEN
          C(I+1)=C(I+1)/1.D30
          D(I+1)=D(I+1)/1.D30
          C(I+2)=C(I+2)/1.D30
          D(I+2)=D(I+2)/1.D30
        ENDIF
        IF(ABS(C(I+2)).LT.1.D-30) THEN
          C(I+1)=C(I+1)*1.D30
          D(I+1)=D(I+1)*1.D30
          C(I+2)=C(I+2)*1.D30
          D(I+2)=D(I+2)*1.D30
        ENDIF
        C2=C(I+2)/D(I+2)
        IF(ABS(C2-C1).LT.1.E-9) EXIT
        C1=C2
      ENDDO
      IF(C2.LT.0.0) PRINT 51,A,B,X,C2
      B1=A*LOG(X)+B*LOG(1-X)+LOG(C2)-LOG(A)
      B1=B1+GAMMAL(A+B)-GAMMAL(A)-GAMMAL(B)
      BETCON1=EXP(B1)
      END
