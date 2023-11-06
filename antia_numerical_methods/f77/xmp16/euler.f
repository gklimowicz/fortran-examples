!	To sum an alternating series using Euler's transform

      PROGRAM EULERSUM
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL TERM

!     EXERCISE 6.37
 
51    FORMAT(5X,'N =',I6,5X,'M1 =',I5,5X,'M2 =',I5)
52    FORMAT('IER =',I4,5X,'N1 =',I5,5X,'N2 =',I5/5X,'SUM =',
     1       1PD14.6,5X,'ESTIMATED ERROR =',D14.6)
 
      REPS=1.Q-32
      AEPS=1.Q-39
      A0=1.0
100   PRINT *,'TYPE  N=NO. OF TERMS,  M1, M2'
      PRINT *,'    M1,M2=NO. OF TERMS TO BE SUMMED DIRECTLY AT TWO ENDS'
      PRINT *,'                                 (QUITS WHEN N < 0)'
      READ *,N,M1,M2
      IF(N.LT.0) STOP
      WRITE(6,51) N,M1,M2
      CALL EULER(N,M1,M2,A0,REPS,AEPS,DIF,N1,N2,SUM,IER,TERM)
      WRITE(6,52) IER,N1,N2,SUM,DIF
      GO TO 100
      END
 
!     -------------------------------------------------------------
 
!	Summation of alternating series using Euler transform
!
!	N : (input) Number of terms to be summed.
!		If N<1 it is assumed to be infinite
!	M1 : (input/output) Number of terms in the beginning to be summed
!		separately. If M1+M2+2*NMAX .GE. N then M1 is set to N and
!		the sum is evaluated by direct summation
!	M2 : (input) Number of terms at the end to be summed separately 
!		The Euler transform is applied to terms from M1+1 to N-M2
!	A0 : (input) The sign of first term. The sum assuming the first
!		term to be positive is multiplied by A0
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(SUM))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	N1 : (output) The number of terms actually used at the beginning
!	N2 : (output) The number of terms actually used at the end
!		N2 is relevant only for finite series
!	SUM : (output) The calculated value of the sum
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=31 implies M1+M2+2*NMAX .GE. N, in which case direct
!			sum is calculated
!		IER=32 implies that the transformed series did not converge to
!			specified accuracy at lower end
!		IER=34 implies that the transformed series did not converge to
!			specified accuracy at upper end
!		IER=36 implies that the transformed series did not converge to
!			specified accuracy at both ends
!	TERM : (input) Name of the function routine to calculate the terms of
!		series, FUNCTION TERM(I) should calculate the Ith term
!		without sign. The Ith term of the series will be
!		A0*TERM(I)*(-1)**(I-1). Here I ranges from 1 to N
!
!	FUNCTION TERM(I) must be supplied by the user
!
!	Required routines : TERM

      SUBROUTINE EULER(N,M1,M2,A0,REPS,AEPS,DIF,N1,N2,SUM,IER,TERM)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER (NMAX=20)
      DIMENSION D(NMAX,NMAX)
 
      IER=0
      DIF=0.0
      N2=0
      IF(M1+M2+2*NMAX.GE.N.AND.N.GT.0) THEN
!	Sum the series directly
        M1=N
        N1=N
        IER=31
      ENDIF

!	Sum the first M1 terms separately
      S1=0.0
      A=A0
      DO 1000 I=1,M1
        S1=S1+A*TERM(I)
        A=-A
1000  CONTINUE
      SUM=S1
      IF(M1.GE.N.AND.N.GT.0) RETURN
 
!	Sum the remaining terms using Euler transform
      S2=0.0
      T1=-1
      DO 2000 I=M1+1,M1+NMAX
        D(I-M1,1)=TERM(I)
        DO 1800 J=2,I-M1
!	The differences
          D(I-M1,J)=D(I-M1,J-1)-D(I-M1-1,J-1)
1800    CONTINUE
        T1=-T1/2.
        TS=T1*D(I-M1,I-M1)
        S2=S2+TS
        IF(ABS(TS).LT.MAX(AEPS,REPS*ABS(S1+S2))) GO TO 2500
2000  CONTINUE
!	Sum of transformed series does not converge
      IER=32
      I=M1+NMAX

2500  N1=I
!	Sum of the infinite series
      SUM=S1+A*S2
      DIF=ABS(TS)
      IF(N.LE.0) RETURN
 
!	For finite series calculate the contribution from the other end also
      S1=0
      A=A0*(-1)**(N-1)
!	Sum the last M2 terms separately
      DO 3000 I=1,M2
        S1=S1+A*TERM(N-I+1)
        A=-A
3000  CONTINUE
 
!	Sum the remaining terms using Euler transform
      S2=0.0
      T1=-1
      DO 4000 I=M2+1,M2+NMAX
        D(I-M2,1)=TERM(N-I+1)
        DO 3800 J=2,I-M2
          D(I-M2,J)=D(I-M2,J-1)-D(I-M2-1,J-1)
3800    CONTINUE
        T1=-T1/2.
        TS=T1*D(I-M2,I-M2)
        S2=S2+TS
        IF(ABS(TS).LT.MAX(AEPS,REPS*ABS(SUM+S1+S2))) GO TO 4500
4000  CONTINUE
!	Sum of transformed series does not converge
      IER=IER+4
      I=M2+NMAX

4500  N2=I
!	Add the two parts to get the sum of finite series
      DIF=DIF+ABS(TS)
      SUM=SUM+S1+A*S2
      IF(IER.GT.0.AND.IER.LT.30) IER=IER+30
      END
 
!     ----------------------------------------------------------
 
!	TO CALCULATE TERMS OF SERIES WITHOUT SIGN 

      FUNCTION TERM(I)
      IMPLICIT REAL*16(A-H,O-Z)
 
      TERM=1.Q0/(2.*I-1.Q0)
      END
