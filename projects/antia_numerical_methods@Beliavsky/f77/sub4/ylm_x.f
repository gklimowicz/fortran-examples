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
      IMPLICIT REAL*4(A-H,O-X)
      IMPLICIT COMPLEX*8(Y,Z)
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
