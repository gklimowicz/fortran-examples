!     TO EVALUATE SPHERICAL HARMONIC

      PROGRAM SPHAR
      IMPLICIT REAL*8(A-H,O-X)
      IMPLICIT COMPLEX*16(Y,Z)
 
51    FORMAT('   L =',I6,5X,'M =',I6,5X,'THETA =',1PD14.6,5X,
     1       'PHI =',D14.6/5X,'COS(THETA) =',D14.6/5X,
     2       'YLM(X) =',2D14.6,5X,2D14.6)

 
100   PRINT *,'TYPE L, M, THETA, PHI      (QUITS WHEN L<0)'
      READ *,L,M,THETA,PHI
      IF(L.LT.0) STOP

      Y1=YLM(L,M,THETA,PHI)

      X=COS(THETA)
!	Calculate YLM using COS(THETA) also and print both values
!	These should be identical
      Y2=YLM_X(L,M,X,PHI)
      WRITE(6,51) L,M,THETA,PHI,X,Y1,Y2
      GO TO 100
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
 
!     ----------------------------------------------------
 
!	To compute spherical harmonic Y_lm(THETA,PHI)
!
!	L : (input) Degree of spherical harmonic
!	M : (input) Azimuthal order of spherical harmonic
!	THETA, PHI : (input) Real variables specifying the angular coordinates 
!		at which the spherical harmonic needs to be evaluated
!	YLM is complex value of spherical harmonic and must be declared
!		to be complex in the calling routine
!
!	Required routines : PLM
 
      FUNCTION YLM(L,M,THETA,PHI)
      IMPLICIT REAL*8(A-H,O-X)
      IMPLICIT COMPLEX*16(Y,Z)
      PARAMETER(LMAX=5001,PI=3.14159265358979324D0)
      DIMENSION P(LMAX)
 
      YLM=0.0
      IF(L.LT.0.OR.ABS(M).GT.L.OR.L.GE.LMAX) RETURN
      MM=ABS(M)
!	To use X instead of THETA in argument comment out this line
      X=COS(THETA)
      CALL PLM(L,MM,X,P)

      CLM=(2*L+1.)/(4.*PI)
      DO 200 I=L-MM+1,L+MM
        CLM=CLM/I
200   CONTINUE
      CLM=SQRT(CLM)

      IF(MOD(MM,2).EQ.1.AND.M.GE.0) CLM=-CLM
      YI=(0.D0,1.D0)
      YLM=CLM*P(L+1)*EXP(YI*M*PHI)
      END
 
!     ----------------------------------------------------
 
!	To compute spherical harmonic Y_lm(COS(THETA),PHI)
!	Version of YLM with argument as COS(THETA) instead of THETA
!
!	L : (input) Degree of spherical harmonic
!	M : (input) Azimuthal order of spherical harmonic
!	X, PHI : (input) Real variables specifying the values of COS(THETA)
!		and PHI at which the spherical harmonic needs to be evaluated
!	YLM_X is complex value of spherical harmonic and must be declared
!		to be complex in the calling routine
!
!	Required routines : PLM
 
      FUNCTION YLM_X(L,M,X,PHI)
      IMPLICIT REAL*8(A-H,O-X)
      IMPLICIT COMPLEX*16(Y,Z)
      PARAMETER(LMAX=5001,PI=3.14159265358979324D0)
      DIMENSION P(LMAX)
 
      YLM_X=0.0
      IF(L.LT.0.OR.ABS(M).GT.L.OR.L.GE.LMAX) RETURN
      MM=ABS(M)
!	To use THETA instead of X in argument uncomment  this line
!      X=COS(THETA)
      CALL PLM(L,MM,X,P)

      CLM=(2*L+1.)/(4.*PI)
      DO 200 I=L-MM+1,L+MM
        CLM=CLM/I
200   CONTINUE
      CLM=SQRT(CLM)

      IF(MOD(MM,2).EQ.1.AND.M.GE.0) CLM=-CLM
      YI=(0.D0,1.D0)
      YLM_X=CLM*P(L+1)*EXP(YI*M*PHI)
      END
