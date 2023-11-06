!	To calculate the Logarithm of Gamma function for a real argument
!	For negative values it give ln(abs(Gamma(x)))
!
!	Required routines : None
 
      FUNCTION GAMMAL(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(4),B(4),A1(6),B1(6)
!	PI2L=LOG(2*PI)/2
      PARAMETER (PI2L=0.918938533204672741D0,PI=3.141592653589793238D0)

!	The coefficients of rational function approximations
      DATA A/ 1.013142782275024216D-2,  7.645657825398191944D-1,
     1        3.381172379819461227D-4,  1.595363637547538209D-2/
      DATA B/0.0, 8.333333333338911768D-2, 8.442856404442060242D-4,
     1            6.093603832366013704D-2/

      DATA A1/ 4.681163846241230144D0,  3.208225429683256526D0,
     1         5.145525793448859216D-1, 1.581117883959157936D-2,
     2        -6.398416804905407512D-5, 5.264566254181773919D-7/
      DATA B1/ 2.938038561191284576D0,  1.489364948862436743D0,
     1        -5.466291543917642961D0,  1.972497734170110410D-1,
     2         7.830146473241555157D-1, 5.756753067834747499D-2/
 
      T=ABS(X)
      IT=T
      IF(T.GE.10.0) THEN
!	Use asymptotic approximation for T>10
        Y=1/T
        FN=((B(4)*Y+B(3))*Y+B(2))*Y+B(1)
        FD=(((A(4)*Y+A(3))*Y+A(2))*Y+A(1))*Y+1
        GAMMAL=FN/FD-T+PI2L+(T-0.5)*LOG(T)
      ELSE
 
!	Use approximation for [4,5]
        FAC=0.0
        T1=T
        IF(T.LT.4) THEN
          DO I=IT,3
            FAC=FAC-LOG(T1)
            T1=T1+1
          ENDDO
        ELSE IF(T.GT.5) THEN
          DO I=IT,5,-1
            T1=T1-1
            FAC=FAC+LOG(T1)
          ENDDO
        ENDIF
        FN=((((B1(6)*T1+B1(5))*T1+B1(4))*T1+B1(3))*T1+B1(2))*T1+B1(1)
        FD=(((((A1(6)*T1+A1(5))*T1+A1(4))*T1+A1(3))*T1+A1(2))*T1+
     1          A1(1))*T1+1
        GAMMAL=FN/FD+FAC
      ENDIF

      IF(X.LE.0) THEN
        IF(T.GT.IT) THEN
          GAMMAL=LOG(PI)-LOG(ABS(SIN(PI*T)*T))-GAMMAL
        ELSE
          GAMMAL=(-1)**IT/0.0
        ENDIF
      ENDIF
 
      END
