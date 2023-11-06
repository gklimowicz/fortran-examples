!	To integrate a function over a finite interval using Gauss-Kronrod formula
!	For use with ADPINT
!	Since the weights and abscissas are not accurate to REAL*16
!	accuracy, this routine will not achive the maximum accuracy
!	permissible by arithmetic. Use GAUS16 instead of KRONRD
!	if high accuracy is required.
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

      SUBROUTINE KRONRD(RI,A,B,DIF,N,F)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION  W7(4),A7(4),WK7(4),WK15(4),AK15(4)

!	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
!	WK7 are the weights for these points in Kronrod formula
!	WK15 and AK15 are the weights and abscissas for the remaining points
!	in Kronrod formula.
!	Because of symmetry only half the points are given.

      DATA W7  /0.1294849661688696932706114326790820183Q0,
     *          0.2797053914892766679014677714237795825Q0,
     *          0.3818300505051189449503697754889751339Q0,
     *          0.4179591836734693877551020408163265306Q0/
      DATA A7  /0.9491079123427585245261896840478512624Q0,
     *          0.7415311855993944398638647732807884070Q0,
     *          0.4058451513773971669066064120769614633Q0, 0.0/
      DATA WK7 /0.0630920926299785532907006631892042866Q0,
     *          0.1406532597155259187451895905102379204Q0,
     *          0.1903505780647854099132564024210136828Q0,
     *          0.2094821410847278280129991748917142637Q0/
      DATA WK15/0.0229353220105292249637320080589695920Q0,
     *          0.1047900103222501838398763225415180174Q0,
     *          0.1690047266392679028265834265985502841Q0,
     *          0.2044329400752988924141619992346490847Q0/
      DATA AK15/0.9914553711208126392068546975263285166Q0,
     *          0.8648644233597690727897127886409262012Q0,
     *          0.5860872354676911302941448382587295984Q0,
     *          0.2077849550078984676006894037732449134Q0/

      AT=(B-A)/2.
      BT=(B+A)/2.
      FBT=F(BT)
      R1=W7(4)*FBT
      RI=WK7(4)*FBT
      DO 2000 K=1,3
        F1=F(AT*A7(K)+BT)
        F2=F(BT-AT*A7(K))
!	7-point Gauss-Legendre formula
        R1=R1+W7(K)*(F1+F2)
!	15-point Kronrod formula
        RI=RI+WK7(K)*(F1+F2)
2000  CONTINUE

      DO 2500 K=1,4
2500  RI=RI+WK15(K)*(F(AT*AK15(K)+BT)+F(BT-AT*AK15(K)))

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=15
      END
