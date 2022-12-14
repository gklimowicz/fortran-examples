! Largely basd on discussions at http://en.wikipedia.org/wiki/Kepler's_laws_of_planetary_motion
module OrbitUtilities_mod
  use KindParameters_mod, only: DP
  implicit none

  private

  public :: computeTrueAnomaly
  public :: computeMeanAnomaly
  public :: computeDistance

  interface computeTrueAnomaly
     module procedure computeTrueAnomaly_main
     module procedure computeTrueAnomaly_parameterized
  end interface computeTrueAnomaly

contains

  ! Uses Newton-Raphson method to invert
  ! M = E - e sin E
  function computeEccentricAnomaly(meanAnomaly, eccentricity) result(eccentricAnomaly)
    real(kind=DP) :: eccentricAnomaly
    real(kind=DP), intent(in) :: meanAnomaly
    real(kind=DP), intent(in) :: eccentricity

    real(kind=DP) :: f, df, x, dx

    associate (MA => meanAnomaly, eps => eccentricity, EA => eccentricAnomaly)
      ! first guess
      x = MA + eps*(sin(MA) + eps*sin(2*MA)/2)
      do
         f = x - eps*sin(x) - MA
         df = 1 - eps*cos(x)
         dx = - f/df
         x = x + dx
         if (abs(dx) < 1.d-10) exit
      end do
      EA = x
    end associate
              
  end function computeEccentricAnomaly


  function computeTrueAnomaly_main(meanAnomaly, eccentricity) result(trueAnomaly)
    real(kind=DP) :: trueAnomaly
    real(kind=DP), intent(in) :: meanAnomaly
    real(kind=DP), intent(in) :: eccentricity

    real(kind=DP) :: eccentricAnomaly
    real(kind=DP), parameter :: pi = 2*asin(1.d0)

    associate(EA=>eccentricAnomaly, e=>eccentricity, nu=>trueAnomaly, M=>meanAnomaly)
      EA = computeEccentricAnomaly(M, e)
      nu = computeTrueAnomaly(M, e, EA)
    end associate

  end function computeTrueAnomaly_main

  function computeTrueAnomaly_parameterized( &
       & meanAnomaly, eccentricity, eccentricAnomaly) result(trueAnomaly)
    real(kind=DP) :: trueAnomaly
    real(kind=DP), intent(in) :: meanAnomaly ! unused - except for disambiguation
    real(kind=DP), intent(in) :: eccentricity
    real(kind=DP), intent(in) :: eccentricAnomaly

    real(kind=DP), parameter :: PI = 2*asin(1.d0)

    associate(EA=>eccentricAnomaly, e=>eccentricity, nu=>trueAnomaly)
      nu = modulo(2*atan(sqrt((1+e)/(1-e)) * tan(EA/2)), 2*PI)
    end associate

  end function computeTrueAnomaly_parameterized


  function computeMeanAnomaly(trueAnomaly, eccentricity) result(meanAnomaly)
    real(kind=DP) :: meanAnomaly
    real(kind=DP), intent(in) :: trueAnomaly
    real(kind=DP), intent(in) :: eccentricity

    real(kind=DP) :: eccentricAnomaly

    associate(EA => eccentricAnomaly, nu => trueAnomaly, M => meanAnomaly, e => eccentricity)
      EA = 2*atan(sqrt((1-e)/(1+e)) * tan(nu/2))
      M = EA - e * sin(EA)
    end associate
  end function computeMeanAnomaly


  function computeDistance(meanDistance, eccentricity, meanAnomaly) result(distance)
    real(kind=DP) :: distance
    real(kind=DP), intent(in) :: meanDistance
    real(kind=DP), intent(in) :: eccentricity
    real(kind=DP), intent(in) :: meanAnomaly

    real(kind=DP) :: eccentricAnomaly
    
    associate (r=>distance, a=>meanDistance, e=>eccentricity, &
         & ma=>meanAnomaly, ea=>eccentricAnomaly)
      ea = computeEccentricAnomaly(ma, e)
      r = a*(1 - e*cos(ea))
    end associate

  end function computeDistance


end module OrbitUtilities_mod
