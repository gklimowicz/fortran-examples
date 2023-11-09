#include "rundeck_opts.h"
module VerticalRes
!@sum Vertical Resolution file, 20 layers, top at .1 mb
!@auth Original Development Team
  use constant, only : grav,mb2kg
#ifdef PLANET_PARAMS
  use PlanetParams_mod, only : PlanetParams
#endif
  Implicit None
!@var LM    = number of dynamical layers
!@var LS1   = lowest layer of strtosphere
  Integer*4,Parameter :: LM=20

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
!@var PLbot pressure levels at bottom of layers (mb)

#ifdef PLANET_PARAMS
  real*8, parameter :: PSF = PlanetParams%psf
#else
  real*8, parameter :: PSF = 984d0
#endif
  real*8, parameter, private :: pratio = PSF/984d0

  real*8, parameter :: &
    PLBOT(1:LM+1) = pratio* &
    (/ 984d0, 964d0, 934d0, 884d0, 810d0, &  ! Pbot L=1,..
    710d0, 550d0, 390d0, 285d0, 210d0, &  !      L=...
    150d0,                             &  !      L=LS1
    110d0,  80d0,  55d0,  35d0,  20d0, &  !      L=...
    10d0,   3d0,   1d0,  .3d0, .1d0 /)    !      L=..,LM+1

  integer, private :: iii ! iterator

!@var delp nominal pressure thicknesses of layers (mb)
  real*8, dimension(lm), parameter, private :: delp = plbot(1:lm)-plbot(2:lm+1)

#ifndef PLANET_PARAMS
  integer, parameter :: ls1=11
#else
! Calculate ls1 from PlanetParams%ptop using allowed compile-time operations.
! This clunky coding will be enthusiastically removed once vertical grid
! information is set at runtime.
  real*8, parameter, private :: ptop0 = PlanetParams%ptop
  integer, dimension(lm), parameter, private :: lev=(/ (iii,iii=1,lm) /)
! factor 1d6 is simply to make all pdiff magnitudes greater than 1
  real*8, parameter, dimension(lm+1), private :: pdiff = 1d6*(plbot-ptop0)
  integer, parameter, dimension(lm), private :: &
                                ! lprod is nonzero if ptop0 lies between two elements of plbot
    lprod = lev*int(min(1d0,max(0d0,-pdiff(1:lm)*pdiff(2:lm+1)))), &
                                ! ldiff is nonzero if ptop0 is exactly equal to an element of plbot
    ldiff = lev*(1-int(min(1d0,abs(pdiff(1:lm)))))
! summation to get the smaller candiate for ls1
  integer, parameter, private :: ls1_lower = &
    lprod( 1)+lprod( 2)+lprod( 3)+lprod( 4)+lprod( 5) &
    +lprod( 6)+lprod( 7)+lprod( 8)+lprod( 9)+lprod(10) &
    +lprod(11)+lprod(12)+lprod(13)+lprod(14)+lprod(15) &
    +lprod(16)+lprod(17)+lprod(18)+lprod(19)+lprod(20) &
                                !
    +ldiff( 1)+ldiff( 2)+ldiff( 3)+ldiff( 4)+ldiff( 5) &
    +ldiff( 6)+ldiff( 7)+ldiff( 8)+ldiff( 9)+ldiff(10) &
    +ldiff(11)+ldiff(12)+ldiff(13)+ldiff(14)+ldiff(15) &
    +ldiff(16)+ldiff(17)+ldiff(18)+ldiff(19)+ldiff(20)
! choose ls1 by rounding ptop to the closest element of plbot
  integer, parameter, private :: &
    wt = 1-nint((plbot(ls1_lower)-ptop0)/delp(ls1_lower))
  integer, parameter :: ls1 = wt*ls1_lower + (1-wt)*(ls1_lower+1)
#endif

  integer, parameter :: ls1_nominal=ls1

  Real*8,Parameter :: &
    PTOP = plbot(ls1), &
    PMTOP = plbot(lm+1), &
    PSFMPT = PSF-PTOP, &
    PSTRAT = PTOP-PMTOP

!@var tropomask unity from layers 1-ls1-1, zero above
!@var stratmask zero from layers 1-ls1-1, unity above
  real*8, dimension(lm), parameter, private :: &
    tropomask = (/ (1d0, iii=1,ls1-1), (0d0, iii=ls1,lm) /), &
    stratmask = 1d0-tropomask

!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV
!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var AM(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)  = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L) = fraction of variable mass in each layer = DSIG(L)
  Real*8,Parameter :: MDRYA = psf*mb2kg, MTOP = pmtop*mb2kg, MFIXs = pstrat*mb2kg, &
    MFIX(LM)  = delp*stratmask*mb2kg, &
    MFRAC(LM) = delp*tropomask/psfmpt

!**** Vertival resolution
!****                         ---MSURF=10034.0---    ---MSURF=5781.8----
!**** Layer   MFIX   MFRAC    MVAR     AM     SUM    MVAR     AM     SUM
!**** =====   ====   =====    ====     ==     ===    ====     ==     ===
!****   Top                                   1.0                    1.0
!****    20    2.0    0        0.0    2.0             0.0    2.0
!****    19    7.1    0        0.0    7.1             0.0    7.1
!****    18   20.4    0        0.0   20.4             0.0   20.4
!****
!****    12  305.9    0        0.0  203.9             0.0  305.9
!****    11  407.9    0        0.0  224.3             0.0  407.9
!****    10    0.0  60/834   611.8  611.8  1529.6   305.9  305.9  1529.6
!****     9    0.0  70/834
!****
!****     3    0.0  50/834
!****     2    0.0  30/834
!****     1    0.0  20/834   203.9  203.9 10034.0   102.0  102.0  5781.8
!****         ----  ------   -----  -----           -----  -----
!****       1528.6 834/834  8504.4 10033.0         4252.2 5780.8  


End Module VerticalRes
