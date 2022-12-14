#include "rundeck_opts.h"

module constant
   use MathematicalConstants_mod, only: PI, TWOPI
   use MathematicalConstants_mod, only: RADIAN
   use MathematicalConstants_mod, only: zero, one
   use MathematicalConstants_mod, only: rt2, byrt2
   use MathematicalConstants_mod, only: rt3, byrt3
   use MathematicalConstants_mod, only: rt12, byrt12
   use MathematicalConstants_mod, only: by3, by6, by9, by12
   use PlanetaryParams_mod, only: PlanetaryParams
#ifdef PLANET_PARAMS
  use PlanetParams_mod, only : exoPlanetParams=>PlanetParams
#endif
!@sum  CONSTANT definitions for physical constants and useful numbers
!@auth G. Schmidt
  implicit none
  save
  !**** Conventions: 'by' implies reciprocal, 'rt' implies square root

  !**** Numerical constants

  real*8,parameter :: UNDEF_VAL=huge(1.d0)
  integer,parameter :: IUNDEF_VAL=huge(1)

!@param undef Missing value
  real*8,parameter :: undef=-1.d30
!@param teeny  small positive value used in num/(den+teeny) to avoid 0/0
  real*8,parameter :: teeny=1.d-30
  integer*8,parameter :: intNaN=-1  ! i.e. = Z'FFFFFFFFFFFFFFFF'
!@param NaN NaN
#if (defined COMPILER_PGI || defined COMPILER_NAG)
  real*8,parameter :: NaN=1d30
#else
  real*8,parameter :: NaN=transfer(intNaN,1.d0)
#endif

  !**** Physical constants

!@param stbo Stefan-Boltzmann constant (W/m^2 K^4)
! real*8,parameter :: stbo =5.6692d-8 ! used in radiation til 2015
! real*8,parameter :: stbo =5.67051d-8 ! year 2000 best estimate 
  real*8,parameter :: stbo =5.67037321d-8 ! CODATA recommended in 2010
! real*8,parameter :: stbo =5.670367d-8 ! CODATA recommended in 2015
!                           +-0.000013d-8 uncertainty 

  !**** Latent heats:
  !**** Note that for energy conservation the efective latent heat at any
  !**** temperature must follow these formulae (assuming a reference
  !**** temperature of 0 Celcius, and constant specific heats).
  !**** If specific heats vary as a function of temperature, the extra
  !**** term becomes an integral
  !**** lhe(T) = lhe(0) + (shv-shw) T (in C)
  !**** lhm(T) = lhm(0) + (shw-shi) T (in C)
  !**** lhs(T) = lhs(0) + (shv-shi) T (in C)
!@param lhe   latent heat of evap at 0 C 2.5008d6 [J kg-1]
  real*8,parameter :: lhe = 2.5d6
!@param lhm   latent heat of melt at 0 C 334590 [J kg-1]
  real*8,parameter :: lhm = 3.34d5
!@param bylhm  1/lhm [kg J-1]
  real*8,parameter :: bylhm = 1./lhm
!@param lhs  latent heat of sublimation at 0 C [J kg-1]
  real*8,parameter :: lhs = lhe+lhm

!@param rhow density of pure water [kg m-3]
  real*8,parameter :: rhow = 1d3
!@param rhows density of average sea water [kg m-3]
  real*8,parameter :: rhows = 1030d0
!@param byrhows recip. density of average sea water [m^3 kg-1]
  real*8,parameter :: byrhows = 1d0/rhows
!@param rhoi density of pure ice [kg m-3]
  real*8,parameter :: rhoi = 916.6d0
!@param byrhoi 1/rhoi [m^3 kg-1]
  real*8,parameter :: byrhoi = 1d0/rhoi

!@param tf freezing point of water at 1 atm [K]
  real*8,parameter :: tf = 273.15d0
!@param bytf 1/tf [K-1]
  real*8,parameter :: bytf = 1d0/tf

!@param shw heat capacity of water (at 20 C) [J kg-1 K-1]
  real*8,parameter :: shw  = 4185.
!@param byshw 1/shw [kg K J-1]
  real*8,parameter :: byshw = 1d0/shw

!@param shi heat capacity of pure ice (at 0 C) [J kg-1 K-1]
  real*8,parameter :: shi  = 2060.
!@param byshi 1/shi [kg K J-1]
  real*8,parameter :: byshi = 1d0/shi

!@param fraction of N2 in dry air (0-1)
      real*8,parameter :: pN2 = 0.780840d0
!@param fraction of O2 in dry air (0-1)
      real*8,parameter :: pO2 = 0.209476d0

  !**** RGAS = R/M_A = 1000* 8.314510 J/mol K /28.9655 g/mol
  !**** For values of CO2 much larger than present day (> 4x conc)
  !**** the molar mass of dry air M_A could change.
  !**** Assume that M_O2 = 31.9988 and M_CO2 = 44.00995
  !**** and current percentages 20.946% and 0.0350% (US Stand. Atm.)
  !**** Assuming CO2 displaces other gases equally M_A=28.9602 + n*0.00527
  !**** where n is multiple of present day CO2 conc (350 ppm)
  !**** For 4xCO2  M_A = 28.9813  => rgas = 286.89
  !**** For 10xCO2 M_A = 29.0129  => rgas = 286.58
!@param gasc  gas constant [J K-1 mol-1]
  real*8,parameter :: gasc = 8.314510d0
!@param bygasc  1/gasc [K mol J-1]
  real*8,parameter :: bygasc = 1./gasc
!@param mair molar mass of dry air (28.9655 g/mol)
#ifdef PLANET_PARAMS
  real*8,parameter :: mair = exoPlanetParams%mair
#else
  real*8,parameter :: mair = 28.9655d0
#endif
!@param rgas gas constant (287.05 J/K kg)
  real*8,parameter :: rgas = 1d3 * gasc / mair ! = 287.05...

!@param mwat molar mass of water vapour [g mol-1]
  real*8,parameter :: mwat = 18.015d0
!@param rvap  gas constant for water vapour (461.5) [J K-1 kg-1]
  !**** defined as R/M_W = 1000* 8.314510 J/mol K /18.015 g/mol
  real*8,parameter :: rvap = 1d3 * gasc / mwat ! = 461.5...

!@param mrat  mass ratio of air to water vapour (0.62197) [1]
  real*8,parameter :: mrat = mwat/mair    ! = 0.62197....
!@param bymrat 1/mrat (1.6078) [1]
  real*8,parameter :: bymrat = 1./mrat    ! = 1.6078....
!@param deltx coeff. of humidity in virtual temperature defn. (0.6078) [1]
  real*8,parameter :: deltx = bymrat-1.   ! = 0.6078....

!@param srat ratio of specific heats at const. press. and vol. (=1.401)
#ifdef PLANET_PARAMS
  real*8,parameter :: srat = exoPlanetParams%srat
#else
  real*8,parameter :: srat = 1.401d0
#endif
!@param kapa ideal gas law exponent for dry air (.2862)
  !**** kapa = (g-1)/g where g=1.401 = c_p/c_v
  real*8,parameter :: kapa = (srat - 1.)/srat  ! =.2862....
!@param bykapa,bykapap1,bykapap2 various useful reciprocals of kapa [1]
  real*8,parameter :: bykapa = 1./kapa
  real*8,parameter :: bykapap1 = 1./(kapa+1.)
  real*8,parameter :: bykapap2 = 1./(kapa+2.)

!@param sha specific heat of dry air (const. pres.) (rgas/kapa) [J kg-1 K-1]
  real*8,parameter :: sha = rgas/kapa
!@param bysha 1/sha [kg K J-1]
  real*8,parameter :: bysha = 1./sha

!@param shv specific heat of water vapour (const. pres.) [J kg-1 K-1]
  !**** shv is currently assumed to be zero to aid energy conservation in
  !**** the atmosphere. Once the heat content associated with water
  !**** vapour is included, this can be set to the standard value
  !**** Literature values are 1911 (Arakawa), 1952 (Wallace and Hobbs)
  !**** Smithsonian Met Tables = 4*rvap + delta = 1858--1869 ????
  !     real*8,parameter :: shv = 4.*rvap  ????
  real*8,parameter :: shv = 0.

  !**** air viscosity - temperature independent
!@var visc_air0 dynamic viscosity of air [kg m-1 s-1]
  real*8,parameter :: visc_air0 = 1.7d-5

!@var visc_air_kin0 kinematic viscosity of air (1 bar 15 deg C) [m^2 s-1]
  real*8,parameter :: visc_air_kin0 = 1.46d-5

!@var visc_wtr_kin kinematic viscosity of water (35 psu, 20 deg C) [m^2 s-1]
  real*8,parameter :: visc_wtr_kin = 1.05d-6

!@var avog Avogadro's constant (molecules/mole)
!@var byavog 1 over Avogadro's constant (molecules/mole)^-1
  real*8,parameter :: avog=6.02214129d23
  real*8,parameter :: byavog=1.d0/avog

!@var loschmidt_constant Loschmidt constant (cm-3 at STP)
  real*8,parameter :: loschmidt_constant = 2.6867805D+19

!**** Astronomical constants
!@param daysPerYear number of solar days per orbital period
  real*8, protected :: daysPerYear
!@param omega earth's rotation rate (7.29 s^-1)
  real*8, protected :: omega
!@param omega2 2*omega
  real*8, protected :: omega2

!@param radius radius of the earth (6371000 m, IUGG)
#ifdef PLANET_PARAMS
  real*8,parameter :: radius = exoPlanetParams%radius
#else
  real*8,parameter :: radius = 6371000.
#endif
!@param areag surface area of the earth (m^2)
  real*8,parameter :: areag = 4.*pi*radius*radius

!@param grav gravitaional accelaration (9.80665) [m s-2]
  !**** SI reference gravity (at 45 deg) = 9.80665
#ifdef PLANET_PARAMS
  real*8,parameter :: grav = exoPlanetParams%grav
#else
  real*8,parameter :: grav = 9.80665d0
#endif
!@param bygrav 1/grav
  real*8,parameter :: bygrav = 1d0/grav

  !**** lapse rate related variables
!@param GAMD dry adiabatic lapse rate (0.0098) [K m-1]
  real*8, parameter :: gamd = grav*kapa/rgas
!@param BMOIST moist adiabatic lapse rate [K m-1]
  real*8, parameter :: bmoist = 0.0065d0
!@param BBYG moist adiabatic lapse rate divided by grav [K s^2 m-2]
  real*8, parameter :: bbyg = bmoist*bygrav
!@param GBYRB grav divided by rgas and bmoist [kg m^2 s-1 J-1]
  real*8, parameter :: gbyrb = grav/(rgas*bmoist)

  !**** Useful conversion factors

!@param kg2mb conversion from kg/m^2 to milli-bars [mbar m^2 kg-1]
!@param mb2kg conversion from milli-bars to kg/m^2 [kg m-2 mbar-2]
  real*8,parameter :: kg2mb = 1d-2*grav, mb2kg = 1d2*bygrav
!@param kgpa2mm conversion from kg/m^2 water to mm [mm m^2 kg-1]
!@param mm2kgpa conversion from mm water to kg/m^2 [kg mm-1 m-2]
  real*8,parameter :: kgpa2mm = 1d0, mm2kgpa = 1d0

#ifdef PLANET_PARAMS
  character(len=16), parameter :: planet_name=exoPlanetParams%name
#else
  character(len=16), parameter :: planet_name='Earth'
#endif


  type (PlanetaryParams) :: planetParams
  logical :: init = .false.

contains

   subroutine initializeConstants()
      ! Some qty's can no longer be Fortran PARAMETERs due to the
      ! need to support exoplanet run-time configuration

      real*8 :: rotationPeriod
      real*8 :: orbitalPeriod
      real*8 :: rotationsPerYear

      if (init) return

      init = .true.

      planetParams = PlanetaryParams() ! read from rundeck
      rotationPeriod = planetParams%getSiderealRotationPeriod()
      orbitalPeriod = planetParams%getSiderealOrbitalPeriod()
      omega = 2*pi/rotationPeriod
      omega2 = 2*omega

      rotationsPerYear = orbitalPeriod / rotationPeriod

      daysPerYear = rotationsPerYear - 1 ! minus the revolution contribution

!!$      syr = orbitalPeriod


   end subroutine initializeConstants

  real*8 function visc_air(T)
!@sum visc_air dynamic viscosity of air (function of T) (kg/m s)
!@auth Sutherland formula
    real*8, intent(in) :: T  ! temperature (K)
    real*8, parameter :: n0=1.827d-5, T0=291.15d0, C=120d0

    visc_air = n0*sqrt((T/T0)**3)*(T0+C)/(T+C)

    return
  end function visc_air

  real*8 function visc_air_kin(T)
!@sum visc_air_kin kinematic viscosity of air (function of T) (m2/s)
!@auth COARE formula - Andreas (1989) CRREL Rep. 89-11
    real*8, intent(in) :: T  ! temperature (K)
    real*8, parameter :: nu0=1.326d-5, a0=6.542d-3, b0=8.301d-6, &
         &     c0=4.84d-9
    real*8 :: Tc  ! temperature in deg C

    Tc=T-tf
    visc_air_kin = nu0*(1.+Tc*(a0+Tc*(b0-c0*Tc)))   !m2/s

    return
  end function visc_air_kin


end module constant
