!	To calculate the probability that two uncorrelated sequences
!	of length (n+2) will give a correlation coefficient exceeding |x|
!	For large even values of n the series can give large roundoff
!	in which case the distribution is approximated by Normal distribution
!
!       N  : (input) Length of data sets should be N+2
!       XX : (input) The function will calculate the probability for
!               correlation exceeding |x| for uncorrelated sequences
!
!	Required routines : GAMMAL, ERF

      FUNCTION PCOR(N,XX)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(PI=3.14159265358979323846Q0,PIS=1.7724538509055144Q0)

      X=ABS(XX)
      PCOR=0.0
      IF(X.GT.1.0) RETURN
      IF(2*(N/2).EQ.N) THEN
!	if n is even
        L=(N-2)/2
        SS=X
        P=X
        PMAX=P
        DO K=1,l
          P=-(L-K+1)*P*X*X/K
          SS=SS+P/(2*K+1.Q0)
          IF(ABS(P).GT.PMAX) PMAX=ABS(P)
        ENDDO
        PL=GAMMAL((N+1.Q0)/2)-GAMMAL(N/2.Q0)
        PCOR=2.*EXP(PL)*SS/PIS
        IF(PMAX.GT.1.Q5.OR.PCOR.GT.1.0) PCOR=ERF(X*SQRT(N/2.Q0))
      ELSE
!	if n is odd
        L=(N-3)/2
        SS=SQRT(1.-X*X)
        P=SS
        IF(N.EQ.1) SS=0.0
        DO K=1,L
          P=P*(2*K)*(1-X*X)/(2*K+1.Q0)
          SS=SS+P
        ENDDO
        PCOR=(ASIN(X)+X*SS)*2/PI
      ENDIF
      PCOR=1-PCOR
      END
