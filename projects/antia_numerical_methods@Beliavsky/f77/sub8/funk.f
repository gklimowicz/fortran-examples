!	Function routine to calculate the integrand for calculating
!	PSI(I,X) as required by subroutine FREDCO.
!
!	Function FKER(X,T) is the kernel K(x,t) and PHI(I,T) is the Ith
!	basis function, phi_i(t). The argument X and I are passed through
!	the common block.
!
!	FUNCTION FKER(X,T) and FUNCTION PHI(I,T) must be supplied by the user
!
!	Required routines : FKER, PHI

      FUNCTION FUNK(T)
      IMPLICIT REAL*8(A-H,O-Z)
!	To pass parameters from subroutine FREDCO
      COMMON/ZZFRED/X,I

      FUNK=FKER(X,T)*PHI(I,T)
      END
