!	To calculate the integrand for incomplete beta function
!       It is used by function BETAI
!
!	A,B : Arguments for the complete Beta function passed through
!               a common block
!	X : (input) Upper limit of integration defining the incomplete
!               Beta function
!
!	Required routines : none

      FUNCTION FBETA(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/BETA/A,B
      FBETA=X**(A-1)*(1-X)**(B-1)
      END
