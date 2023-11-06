!	Function routine to transform from hyper-spherical coordinates to
!	Cartesian coordinates.
!	It can be used with any of the subroutines for multiple integration
!	for integration over hyper-spherical shells
!
!	N : (input) Number of dimensions
!	X : (input) Real array of length N, giving hyper-spherical coordinates 
!		First coordinate is radial distance r
!	Y : (output) Real array of length N giving the transformed Cartesian
!		coordinates, which is passed on to the FUNCTION FUNSPH(N,Y)
!		for calculating the required function.
!	Since there is no provision to pass error parameter from this
!	routine, if N exceeds NMAX then the execution will terminate.
!
!	The FUNCTION FUNSPH(N,Y) must be supplied by the user. The name of
!	this routine is not passed as argument and hence it has to have
!	the same name as occurring here.
!		
!	Required routines : FUNSPH

      FUNCTION SPHND(N,X)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=50)
      DIMENSION X(N),Y(NMAX)
 
      IF(N.GT.NMAX.OR.N.LT.1) STOP 301

!	The Cartesian coordinates are given by
!	y(1)=x(1)*COS(x(2))
!	y(2)=x(1)*SIN(x(2))*COS(x(3))
!	y(3)=x(1)*SIN(x(2))*SIN(x(3))*COS(x(4))
!	.........
!	y(n-1)=x(1)*SIN(x(2))*SIN(x(3))*...*SIN(x(n-1))*COS(x(n))
!	y(n)=x(1)*SIN(x(2))*SIN(x(3))*...*SIN(x(n))
!
!	P1 is the volume element in hyper-spherical coordinates

      T1=X(1)
      P1=1.0
      DO 1000 I=1,N-1
        Y(I)=T1*COS(X(I+1))
        P1=P1*T1
        T1=T1*SIN(X(I+1))
1000  CONTINUE
      Y(N)=T1

!	FUNSPH should calculate the required integrand using Cartesian coordinates
      SPHND=P1*FUNSPH(N,Y)
 
      END
