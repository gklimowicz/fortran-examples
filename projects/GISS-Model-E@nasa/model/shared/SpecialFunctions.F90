!------------------------------------------------------------------------------
Module SpecialFunctions_mod
!------------------------------------------------------------------------------
!@sum Module that contains methods to calculate special functions
!@auth SSSO ASTG
  implicit none
  private
  public besselI0

  public erfGam ! used by AMP tracers
  public erf
  interface erf
    module procedure erf_single
    module procedure erf_double
  end interface
  
  integer, parameter :: DP = selected_real_kind(14)
  INTEGER, PARAMETER :: maxNumIterations=10000
  REAL(kind=DP), PARAMETER :: epsCF=3.0D-07
  REAL(kind=DP), PARAMETER :: epsSR=3.0D-09

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  ! Based upon Abromowitz and Stegun 7.1.26
  ! http://en.wikipedia.org/wiki/Abramowitz_and_Stegun
  ! Original source:
  ! C. Hastings, Jr.  Approximations for Digital Computers. Princeton University Press, 1955.

  real(kind=DP) function erf_double(x)
    real(kind=DP), intent(in) :: x

    real(kind=DP), parameter ::  p = +0.3275911
    real(kind=DP), parameter :: a1 = +0.254829592
    real(kind=DP), parameter :: a2 = -0.284496736
    real(kind=DP), parameter :: a3 = +1.421413741
    real(kind=DP), parameter :: a4 = -1.453152027
    real(kind=DP), parameter :: a5 = +1.061405429

    real(kind=DP) :: t

    t = 1/(1 + p *abs(x))
    erf_double = 1 - (a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5)*exp(-x**2)

    if (x < 0) erf_double = -erf_double

  end function erf_double

  real function erf_single(x)
    real, intent(in) :: x

    erf_single = erf(real(x,kind=DP))

  end function erf_single

  FUNCTION GammaContinuedFraction(A,X)
!@sum Computes the incomplete gamma function P(a,x) using the continued
!@+ fraction representations 
    REAL(kind=DP) :: a, x
    REAL(kind=DP) :: GammaContinuedFraction
!
    ! smallest representable floating-point number.
    REAL(kind=DP), PARAMETER :: smallestFP=1.0D-30 
    REAL(kind=DP) :: an,b,c,d,del,h
    INTEGER :: i

    B = X + 1.0d0 - A
    C = 1.0D0/smallestFP
    D = 1.0D0/B
    H = D
    DO I = 1,maxNumIterations
      AN = -I*(I - A)
      B = B + 2.0D0
      D = AN*D + B
      IF (ABS(D) < smallestFP) D = smallestFP
      C=B+AN/C
      IF (ABS(C) < smallestFP) C = smallestFP
      D = 1.0D0/D
      DEL = D*C
      H = H*DEL
      IF(ABS(DEL - 1.0D0) < epsCF) exit
    END DO
    if (i > maxNumIterations) &
      WRITE(*,*)'GammaContinuedFraction: A TOO LARGE, maxNumIterations TOO SMALL', &
      GammaContinuedFraction,A,X,logGamma(A)
    GammaContinuedFraction = EXP(-X + A*LOG(X) - logGamma(A))*H

  END FUNCTION GammaContinuedFraction


  FUNCTION GammaSeriesFraction(A,X)
!@sum Computes the incomplete gamma function P(a,x) using the series
!@+ fraction representations 
    REAL(kind=DP) :: A,X
    REAL(kind=DP) :: GammaSeriesFraction
    INTEGER :: N
    REAL(kind=DP) :: AP,DEL,SUM

    IF (X >= 0.D+00) THEN
      IF (X < 0.) STOP 'GammaSeriesFraction: X < 0'
      GammaSeriesFraction = 0.D0
      RETURN
    END IF
    AP = A
    SUM = 1.D0/A
    DEL = SUM
    DO N = 1,maxNumIterations
      AP = AP + 1.D0
      DEL = DEL*X/AP
      SUM = SUM + DEL
      IF (ABS(DEL) < ABS(SUM)*epsSR) exit
    END DO
    if (n > maxNumIterations) &
      WRITE(*,*)'GammaSeriesFraction: A TOO LARGE, maxNumIterations TOO SMALL'
    GammaSeriesFraction = SUM*EXP(-X + A*LOG(X) - logGamma(A))

  END FUNCTION GammaSeriesFraction


  REAL(kind=DP) FUNCTION LogGamma(XX)
!@sum The natural log of the gamma function.
    REAL(kind=DP) , intent(in) :: XX
    INTEGER J
    REAL(kind=DP) SER,TMP,X,Y
    REAL(kind=DP) :: stp = 2.5066282746310005D0
    REAL(kind=DP), DIMENSION(6) :: cof = (/76.18009172947146, &
      -86.50532032941677, 24.01409824083091, &
      -1.231739572450155, 0.1208650973866179e-2, &
      -0.5395239384953e-5/)

    X = XX
    Y = X
    TMP = X + 5.5D0
    TMP = (X + 0.5D0)*LOG(TMP) - TMP
    SER = 1.000000000190015D0
    DO J = 1,6
      Y = Y + 1.D0
      SER = SER + COF(J)/Y
    END DO
    logGamma = TMP + LOG(STP*SER/X)

  END FUNCTION LogGamma

  REAL(kind=DP) FUNCTION erfGam(X)
!@sum  The error function (2/sqrt(pi))*Integrate[Exp[-t^2],{t,0,x}] 
    REAL(kind=DP) :: X

    erfGam = 0.d0
    IF(X < 0.0D0)THEN
      erfGam=-incompleteGamma(0.5D0,X**2)
    ELSE
      erfGam= incompleteGamma(0.5D0,X**2)
    END IF

  END FUNCTION erfGam


  REAL(kind=DP) FUNCTION incompleteGamma(A,X)
!@sum  Computes the incomplete gamma function P(a,x) by selecting between 
!@+ series and continued fraction representations based on the size of 
!@+ the arguments.
    REAL(kind=DP) :: A,X

    IF(X < 0.0D0 .OR. A <= 0.0D+00)THEN
      WRITE(*,*)'incompleteGamma: BAD ARGUMENTS'
    END IF
    IF(X < A + 1.0D0)THEN
      incompleteGamma = GammaSeriesFraction(A,X)
    ELSE
      incompleteGamma = 1.0D0 - GammaContinuedFraction(A,X)
    END IF

  END FUNCTION incompleteGamma

  REAL(kind=DP) FUNCTION BesselI0(x)
!@sum Modified Bessel's function I0
    REAL(kind=DP), intent(in) :: x
    REAL(kind=DP) :: ax,y
    REAL(kind=DP), parameter :: p1=1.0d0, p2=3.5156229d0, p3=3.0899424d0, &
      p4=1.2067492d0, p5=0.2659732d0, p6=0.360768d-1, &
      p7=0.45813d-2
    REAL(kind=DP), parameter :: q1=0.39894228d0,q2=0.1328592d-1, &
      q3=0.225319d-2, q4=-0.157565d-2, q5=0.916281d-2, &
      q6=-0.2057706d-1, q7=0.2635537d-1, q8=-0.1647633d-1, &
      q9=0.392377d-2

    if (abs(x) < 3.75d0) then
      y = (x/3.75d0)**2.d0
      BesselI0 = (p1 + y*(p2 + y*(p3 + y*(p4 + y*(p5 + y*(p6 + y*p7))))))
    else
      ax = abs(x)
      y = 3.75d0/ax
      BesselI0 = (exp(ax)/sqrt(ax))*(q1 + y*(q2 + y*(q3 + y*(q4 &
        + y*(q5 + y*(q6 + y*(q7 + y*(q8 + y*q9))))))))
    end if

  END FUNCTION BesselI0

end Module SpecialFunctions_mod
