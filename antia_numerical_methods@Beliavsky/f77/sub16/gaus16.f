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
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION  W8(4),A8(4),W16(8),A16(8)

!	W8 and A8 are the weights and abscissas for the 8-point Gauss formula
!	W16 and A16 are the weights and abscissas for the 16-point Gauss formula
!	Because of symmetry only half the points are given.

      DATA W8 /  1.0122853629037625915253135430996220Q-1,
     *           2.2238103445337447054435599442624126Q-1,
     *           3.1370664587788728733796220198660110Q-1,
     *           3.6268378337836198296515044927719567Q-1/
      DATA A8 /  9.6028985649753623168356086856947230Q-1,
     *           7.9666647741362673959155393647583006Q-1,
     *           5.2553240991632898581773904918924604Q-1,
     *           1.8343464249564980493947614236018303Q-1/

      DATA W16/  2.7152459411754094851780572456019197Q-2,
     *           6.2253523938647892862843836994377013Q-2,
     *           9.5158511682492784809925107602247148Q-2,
     *           1.2462897125553387205247628219201544Q-1,
     *           1.4959598881657673208150173054747880Q-1,
     *           1.6915651939500253818931207903036093Q-1,
     *           1.8260341504492358886676366796921865Q-1,
     *           1.8945061045506849628539672320828440Q-1/
      DATA A16/  9.8940093499164993259615417345033306Q-1,
     *           9.4457502307323257607798841553460820Q-1,
     *           8.6563120238783174388046789771239397Q-1,
     *           7.5540440835500303389510119484744268Q-1,
     *           6.1787624440264374844667176404879140Q-1,
     *           4.5801677765722738634241944298357810Q-1,
     *           2.8160355077925891323046050146049710Q-1,
     *           9.5012509837637440185319335424958114Q-2/

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
