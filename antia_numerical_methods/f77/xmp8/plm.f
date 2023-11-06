!     TO EVALUATE ASSOCIATED LEGENDRE FUNCTION

      PROGRAM LEGEND
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(5000)
 
51    FORMAT('   L =',I6,5X,'M =',I6,5X,'X =',1PD14.6/5X,
     1       'PL(X) =',D14.6,5X,'PLM(X) =',D14.6)

 
100   PRINT *,'TYPE X, L, M      (QUITS WHEN X<-1000)'
      READ *,X,L,M
      IF(X.LT.-1000) STOP

      CALL PLEG(L,X,P)
      F1=P(L+1)
      CALL PLM(L,M,X,P)
      WRITE(6,51) L,M,X,F1,P(L+1)
      GO TO 100
      END
 
!     ----------------------------------------------------
 
!	To calculate Legendre polynomial P_l(X)
!
!	L : (input) Order of polynomial (L.GE.0)
!	X : (input) Argument at which the value of polynomial is required
!	P : (output) Real array of length L+1, which will contain the
!		calculated values of polynomials. P(j+1) will contain P_j(X)
!
!	Required routines : None
 
      SUBROUTINE PLEG(L,X,P)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(L+1)
 
      IF(L.LT.0) RETURN
      P(1)=1
      P(2)=X
      DO 2000 N=2,L
        P(N+1)=((2*N-1)*X*P(N)-(N-1)*P(N-1))/N
2000  CONTINUE
      END
 
!     ----------------------------------------------------
 
!	To calculate the associated Legendre functions P_lm(X)
!
!	L,M : (input) Order of function (L.GE.0), ABS(M).LE.L
!	X : (input) Argument at which the value of polynomial is required
!	P : (output) Real array of length L+1, which will contain the
!		calculated values of polynomials. P(j+1) will contain
!		P_jM(X) for j.GE.M
!
!	Required routines : None
 
      SUBROUTINE PLM(L,M,X,P)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(L+1)
 
      IF(L.LT.0.OR.ABS(M).GT.L) RETURN
!	First compute P_MM(x)
      MM=ABS(M)
      PM=1.
      DO 2000 I=1,2*MM-1,2
        PM=PM*I
2000  CONTINUE

      IF(M.LT.0) THEN
!	Modify the normalisation factor
        DO 2100 I=1,2*MM
          PM=PM/I
2100    CONTINUE
        IF(MOD(MM,2).EQ.1) PM=-PM
      ENDIF
 
      RM=MM/2.D0
      IF(MM.EQ.0) THEN
        P(MM+1)=1.0
      ELSE IF(ABS(X).LT.1) THEN
        P(MM+1)=PM*(1-X*X)**RM
      ELSE
        P(MM+1)=0.0
      ENDIF

!	Use the recurrence relation to compute P_nM(x)
      P(MM+2)=X*P(MM+1)*(2*MM+1)/(MM+1-M)
      DO 2200 N=MM+2,L
        P(N+1)=((2*N-1)*X*P(N)-(N-1+M)*P(N-1))/(N-M)
2200  CONTINUE
      END
