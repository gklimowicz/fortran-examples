!	To calculate the incomplete Beta function I_x(a,b) using
!	the continued fraction (modified form), called by BETAP
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!                       Beta function
!
!	Required routines : GAMMAL

      FUNCTION BETCON(A,B,X)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(500),D(500)
           
51    FORMAT(' ** Roundoff error while evaluating the modified ',
     1  'continued fraction for incomplete Beta function at',/,'  A=',
     2  1PE12.4,'  B=',E12.4,'  X=',E12.4,'  CON. FRAC.=',E12.4)

      C(1)=A
      D(1)=A*(1-X*(A+B)/(A+1))
      B1=1+X*(B-1)/(A+1)-X*(A+1)*(A+B+1)/(A+3)
      C(2)=A*B1
      A1=A*X*X*(B-1)*(A+B)/(A+1)**2
      D(2)=B1*D(1)+A1
      C1=C(2)/D(2)
      DO I=2,499
        D1=-X*(A+I-1)*(A+B+I-1)/((A+2*I-2)*(A+2*I-1))
        D3=-X*(A+I)*(A+B+I)/((A+2*I)*(A+2*I+1))
        D2=X*I*(B-I)/((A+2*I-1)*(A+2*I))
        C(I+1)=C(I)*(A+2*I)*(1+D2+D3)-C(I-1)*(A+2*I)*(A+2*I-2)*D2*D1
        D(I+1)=D(I)*(A+2*I)*(1+D2+D3)-D(I-1)*(A+2*I)*(A+2*I-2)*D2*D1
!	Scale the numerator and denominator to prevent underflow/overflow
        IF(ABS(C(I+1)).GT.1.D30) THEN
          C(I+1)=C(I+1)/1.D30
          D(I+1)=D(I+1)/1.D30
          C(I)=C(I)/1.D30
          D(I)=D(I)/1.D30
        ENDIF
        IF(ABS(C(I+1)).LT.1.D-30) THEN
          C(I+1)=C(I+1)*1.D30
          D(I+1)=D(I+1)*1.D30
          C(I)=C(I)*1.D30
          D(I)=D(I)*1.D30
        ENDIF
        C2=C(I+1)/D(I+1)
        IF(ABS(C2-C1).LT.1.E-9) EXIT
        C1=C2
      ENDDO
      IF(C2.LT.0.0) PRINT 51,A,B,X,C2
      B1=A*LOG(X)+B*LOG(1-X)+LOG(C2)-LOG(A)
      B1=B1+GAMMAL(A+B)-GAMMAL(A)-GAMMAL(B)
      BETCON=EXP(B1)
      END
