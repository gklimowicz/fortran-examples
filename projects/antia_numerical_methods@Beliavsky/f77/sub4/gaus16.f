!	To integrate a function over a finite interval using 16 point
!	Gauss-Legendre formula, for use with ADPINT
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	N : (output) Number of function evaluations used by subroutine
!	F : (input) Name of the function routine to calculate the integrand
!
!	FUNCTION F(X) must be supplied by the user
!
!	Required routines : F

      SUBROUTINE GAUS16(RI,A,B,DIF,N,F)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  W8(4),A8(4),W16(8),A16(8)

!	W8 and A8 are the weights and abscissas for the 8-point Gauss formula
!	W16 and A16 are the weights and abscissas for the 16-point Gauss formula
!	Because of symmetry only half the points are given.

      DATA W8 /0.10122853629037625915D0, 0.22238103445337447054D0,
     *         0.31370664587788728734D0, 0.36268378337836198297D0/
      DATA A8 /0.96028985649753623168D0, 0.79666647741362673959D0,
     *         0.52553240991632898582D0, 0.18343464249564980494D0/

      DATA W16/0.02715245941175409485D0, 0.06225352393864789286D0,
     *         0.09515851168249278481D0, 0.12462897125553387205D0,
     *         0.14959598881657673208D0, 0.16915651939500253819D0,
     *         0.18260341504492358887D0, 0.18945061045506849629D0/
      DATA A16/0.98940093499164993260D0, 0.94457502307323257608D0,
     *         0.86563120238783174388D0, 0.75540440835500303390D0,
     *         0.61787624440264374845D0, 0.45801677765722738634D0,
     *         0.28160355077925891323D0, 0.09501250983763744019D0/

      AT=(B-A)/2.
      BT=(B+A)/2.
      R1=0.0
!	8-point Gauss-Legendre formula
      DO 2000 K=1,4
        R1=R1+W8(K)*(F(AT*A8(K)+BT)+F(BT-AT*A8(K)))
2000  CONTINUE

      RI=0.0
!	16-point Gauss-Legendre formula
      DO 2500 K=1,8
        RI=RI+W16(K)*(F(AT*A16(K)+BT)+F(BT-AT*A16(K)))
2500  CONTINUE

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=24
      END
