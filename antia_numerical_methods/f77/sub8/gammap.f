!	To calculate the incomplete Gamma function P(a,x) for
!	positive real arguments
!	It returns a value of -1 for a.le.0 or x<0
!
!	A : (input) Argument for the complete Gamma function
!	X : (input) Upper limit of integration defining the incomplete
!                       Gamma function
!
!	Required routines : GAMMA, GAMMAL

      FUNCTION GAMMAP(A,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(500),B(500)

      IF(A.LE.0.OR.X.LT.0) THEN
        GAMMAP=-1
        RETURN
      ENDIF

      IF(X.LT.3) THEN
!	Use power series
        S1=1
        T=A
        DO I=1,500
          T=-T*X/I
          S1=S1+T/(A+I)
          IF(ABS(T/(A+I)).LT.1.D-14) EXIT
        ENDDO
        IF(A.LT.140) THEN
          GAMMAP=X**A*S1/GAMMA(A+1)
        ELSE
          GAMMAP=0.0
        ENDIF
      ELSE IF(A.LT.1.2D0*X) THEN
!	Use continued fraction
        C(1)=1
        B(1)=X
        C(2)=1
        B(2)=X+1-A
        C1=C(2)/B(2)
        DO I=2,498,2
          C(I+1)=X*C(I)+(I/2)*C(I-1)
          B(I+1)=X*B(I)+(I/2)*B(I-1)
          C(I+2)=C(I+1)+(I/2+1-A)*C(I)
          B(I+2)=B(I+1)+(I/2+1-A)*B(I)
          IF(ABS(B(I+2)).GT.1.D200) THEN
            C(I+1)=C(I+1)/1.D200
            B(I+1)=B(I+1)/1.D200
            C(I+2)=C(I+2)/1.D200
            B(I+2)=B(I+2)/1.D200
          ENDIF
          IF(ABS(B(I+2)).LT.1.D-200) THEN
            C(I+1)=C(I+1)*1.D200
            B(I+1)=B(I+1)*1.D200
            C(I+2)=C(I+2)*1.D200
            B(I+2)=B(I+2)*1.D200
          ENDIF
          C2=C(I+2)/B(I+2)
          IF(ABS(C2-C1).LT.1.D-12) EXIT
          C1=C2
        ENDDO
        G1=-X+A*LOG(X)+LOG(C2)-GAMMAL(A)
        GAMMAP=1-EXP(G1)
      ELSE
!	Use the power series for a>x
        S2=1
        T=1
        DO I=1,500
          T=T*X/(A+I)
          S2=S2+T
          IF(ABS(T).LT.1.D-14) EXIT
        ENDDO
        G1=-X+A*LOG(X)+LOG(S2)-GAMMAL(A+1)
        GAMMAP=EXP(G1)
      ENDIF
      END
