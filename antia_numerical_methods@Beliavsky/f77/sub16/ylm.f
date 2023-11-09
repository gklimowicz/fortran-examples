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
      IMPLICIT REAL*16(A-H,O-X)
      IMPLICIT COMPLEX*32(Y,Z)
!      PARAMETER(LMAX=5001,PI=3.14159265358979324D0)
      PARAMETER(LMAX=5001,PI=3.14159265358979323846264338327950288Q0)
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
      YI=(0.Q0,1.Q0)
      YLM=CLM*P(L+1)*EXP(YI*M*PHI)
      END
