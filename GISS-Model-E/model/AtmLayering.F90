#include "rundeck_opts.h"

module VerticalRes
!@sum Vertical Resolution and Layering file
  use constant, only : grav,mb2kg
  implicit none

!@var PSF global mean surface pressure  (mb)
  real*8, parameter :: PSF = 984d0

  integer, private :: iii ! iterator
  integer, parameter, private :: lmbig=200 ! needed until PDTs are supported

  type layering_t
!@var LM    = number of dynamical layers
    integer :: lm
    !integer :: ls1
    integer :: ls1_nominal
!@var PLbot nominal pressure levels at bottom of layers (mb)
    real*8, dimension(lmbig) :: plbot
  end type layering_t

!
! Define a set of layering options
!

! 12 layers, top at 10 mb
  type(layering_t), parameter, private :: L12 = layering_t( &
    lm = 12, &
    !ls1 = 9, &
    ls1_nominal = 9, &
    plbot = (/  &
      PSF, 934d0, 854d0, 720d0, 550.d0,   &  ! Pbot L=1,5
      390d0, 285d0, 210d0,                &  !      L=...
      150d0,                              &  !      L=LS1
      100d0,  60d0,  30d0, 10d0           &  !      L=..,LM+1
      , (0d0, iii=12+2,lmbig) /) &
      )

! 20 layers, top at .1 mb
  type(layering_t), parameter, private :: L20 = layering_t( &
    lm = 20, &
    !ls1 = 11, &
    ls1_nominal = 11, &
    plbot = (/  &
      PSF, 964d0, 934d0, 884d0, 810d0,   &  ! Pbot L=1,..
      710d0, 550d0, 390d0, 285d0, 210d0, &  !      L=...
      150d0,                             &  !      L=LS1
      110d0,  80d0,  55d0,  35d0,  20d0, &  !      L=...
      10d0,   3d0,   1d0,  .3d0, .1d0    &  !      L=..,LM+1
      , (0d0, iii=20+2,lmbig) /) &
      )

! 40 layers, top at .1 mb
  type(layering_t), parameter, private :: L40 = layering_t( &
    lm = 40, &
    !ls1 = 24, &
    ls1_nominal = 24, &
    plbot = (/  &
      PSF,   964d0, 942d0, 917d0, 890d0, 860d0, 825d0, &  !  L=1,..
      785d0, 740d0, 692d0, 642d0, 591d0, 539d0, 489d0, &  !  L=...
      441d0, 396d0, 354d0, 316d0, 282d0, 251d0, 223d0, &
      197d0, 173d0,                                    &
      150d0,                                           &  !  L=LS1
      128d0, 108d0,  90d0,  73d0,  57d0,  43d0,  31d0, &  !  L=...
      20d0,  10d0,5.62d0,3.16d0,1.78d0,  1.d0,         &
      .562d0,.316d0,.178d0, .1d0                       &  !  L=..,LM+1
      , (0d0, iii=40+2,lmbig) /) &
      )

! 96 layers, top at .1 mb
  type(layering_t), parameter, private :: L96 = layering_t( &
       lm = 96, &
       !ls1 = 53, &
       ls1_nominal = 53, &
       plbot = (/  &
       PSF,   969d0, 954d0, 939d0, 924d0, 909d0, 894d0, &
       879d0, 864d0, 849d0, 834d0, 819d0, 804d0, 789d0, &
       774d0, 759d0, 744d0, 729d0, 714d0, 697d0, 680d0, &
       658d0, 636d0, 614d0, 591d0, 568d0, 545d0, 522d0, &
       500d0, 478d0, 457d0, 437d0, 417d0, 399d0, 381d0, &
       364d0, 348d0, 333d0, 318d0, 304d0, 290d0, 277d0, &
       264d0, 252d0, 240d0, 228d0, 217d0, 206d0, 195d0, &
       184d0, 174d0, 164d0, 154d0, 144d0, 135d0, 126d0, &
       117d0, 108d0, 100d0,  92d0,  85d0,  78d0,  71d0, &
       65d0,  59d0,  53d0,  48d0,  43d0,  38d0,  34d0, &
       30d0,  26d0,  23d0,  20d0,  17d0,  14d0,  12d0, &
       10d0,  8.552d0,  7.104d0,  6.075d0,  5.047d0, &
       4.316d0,  3.586d0,  3.066d0,  2.547d0,  2.178d0, &
       1.810d0,  1.527d0,  1.262d0,  1.032d0,  0.832d0, &
       0.650d0,  0.476d0,  0.316d0,  0.178d0,  .1d0 &
      , (0d0, iii=96+2,lmbig) /) &
       )

  ! alternate L96:  4x 1-2-1 smoothing of layer thickness (weaker @ upper boundary)
  type(layering_t), parameter, private :: L96smoothed = layering_t( &
       lm = 96, &
       !ls1 = 53, &
       ls1_nominal = 53, &
       plbot = (/  &
       PSF,969d0,954d0,939d0,924d0,909d0, &
       894d0,879d0,864d0,849d0,834d0,819d0, &
       804d0,789d0,774d0,758.9921875d0,743.921875d0,728.61328125d0, &
       712.7109375d0,695.71484375d0,677.18359375d0,657.03515625d0,635.62109375d0,613.43359375d0, &
       590.8203125d0,568d0,545.18359375d0,522.58984375d0,500.40625d0,478.77734375d0, &
       457.8125d0,437.59375d0,418.1796875d0,399.59765625d0,381.84765625d0,364.92578125d0, &
       348.80859375d0,333.4140625d0,318.62890625d0,304.375d0,290.625d0,277.37109375d0, &
       264.58984375d0,252.2265625d0,240.22265625d0,228.55078125d0,217.1875d0,206.078125d0, &
       195.1875d0,184.546875d0,174.1875d0,164.078125d0,154.1875d0,144.546875d0, &
       135.1875d0,126.078125d0,117.19140625d0,108.5859375d0,100.3671875d0,92.58984375d0, &
       85.2265625d0,78.22265625d0,71.5546875d0,65.22265625d0,59.22265625d0,53.5546875d0, &
       48.22265625d0,43.22265625d0,38.5546875d0,34.22265625d0,30.22265625d0,26.55078125d0, &
       23.1875d0,20.078125d0,17.18965625d0,14.5684375d0,12.28657421875d0,10.35730859375d0, &
       8.735375d0,7.3664921875d0,6.210015625d0,5.23436328125d0,4.4119296875d0,3.7186796875d0, &
       3.1341478515625d0,2.64043282890625d0,2.2207877225d0, &
       1.8587369305253123d0,1.54771247844d0,1.2782114910409375d0, &
       1.042731898456875d0,0.83677161146875d0,0.6514690681375d0, &
       0.47725106352312496d0,0.31688139972875d0,0.17859880642375d0, &
       .1d0 &
      , (0d0, iii=96+2,lmbig) /) &
       )

! 102 layers = L96 plus 6 layers extending the top to .002 mb
  type(layering_t), parameter, private :: L102 = layering_t( &
       lm = 102, &
       !ls1 = L96%ls1, &
       ls1_nominal = L96%ls1_nominal, &
       plbot = (/ (L96%plbot(iii), iii=1,97),  &
          0.056d0,  0.032d0,  0.018d0,  0.010d0,  0.005d0, .002d0 &
          , (0d0, iii=102+2,lmbig) /) &
       )

  type(layering_t), parameter, private :: L102smoothed = layering_t( &
       lm = 102, &
       !ls1 = L96smoothed%ls1, &
       ls1_nominal = L96smoothed%ls1_nominal, &
       plbot = (/ (L96smoothed%plbot(iii), iii=1,97),  &
                  (L102%plbot(iii), iii=98,103)  &
          , (0d0, iii=102+2,lmbig) /) &
       )

! Select layering using rundeck CPP symbol ATM_LAYERING
  type(layering_t), parameter, private :: layering = ATM_LAYERING

! Expose the components of the layering selection
  integer, parameter :: lm = layering%lm

  integer, parameter :: ls1_nominal = layering%ls1_nominal

  real*8, parameter :: &
    plbot(1:lm+1) = layering%plbot(1:lm+1)

!@var delp nominal pressure thicknesses of layers (mb)
  real*8, dimension(lm), parameter, private :: delp = plbot(1:lm)-plbot(2:lm+1)

!@var PMTOP model top pressure  (mb)
!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV

  real*8, parameter :: PMTOP = plbot(lm+1)
  real*8, parameter :: MDRYA = psf*mb2kg, MTOP = pmtop*mb2kg

#ifndef STDHYB

! discontinuous-transition sigma-and-CP version of mfix and mfrac:

!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var AM(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)  = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L) = fraction of variable mass in each layer = DSIG(L)

!@var tropomask unity from layers 1-ls1-1, zero above
!@var stratmask zero from layers 1-ls1-1, unity above
  real*8, dimension(lm), parameter, private :: &
    tropomask = (/ (1d0, iii=1,ls1_nominal-1), (0d0, iii=ls1_nominal,lm) /), &
    stratmask = 1d0-tropomask

  real*8, parameter :: MFIXs = (plbot(ls1_nominal)-pmtop)*mb2kg, &
    MFIX(LM)  = delp*stratmask*mb2kg, &
    MFRAC(LM) = delp*tropomask/(psf-plbot(ls1_nominal))

#else

! A standard hybrid (STDHYB) coordinate expresses the k-th layer edge pressure
! pe(k) as the sum of
!   (1) a layer-dependent fixed pressure
!   (2) a layer-dependent coefficient times surface pressure psurf
! Here, this formulation is expressed as
!   pe(k) = wtfix(k)*plbot(k) + (1-wtfix(k))*psige(k)
! where the "sigma-layer" pressure is
!   psige(k) = sige(k)*(psurf-pmtop) + pmtop
! Instead of layer-edge pressure parameters, ModelE currently requires per-layer
! mass parameters mfix and mfrac:
! mfix(k)  = mass in each layer independent of mcol (kg/m^2)
! mfrac(k) = proportionality coeff. for mass in each layer that varies with mcol
! mcol  = spatially and temporally varying column-integrated mass (kg/m^2)
!         including dynamically inactive (radiation-only) mass above model top
! am(k) = air mass of layer k = mfix(k) + mcol*mfrac(k) (kg/m^2)
!
! so the calculations below separate pe(k)-pe(k+1) into psurf-dependent and
! psurf-independent parts.

! STDHYB differs from the non-STDHYB case in that the latter always
! subtracts a fixed mass from MCOL; STDHYB absorbs the subtracted
! fixed mass into mfix.
! non-STDHYB:
!   AM(K) = MFIX(K) + MFRAC(K)*(MCOL - MFIXs - MTOP)  ! subtract every calculation
! STDHYB (mfix, mfrac not the same as those on previous line):
!   AM(K) = MFIX(K)  + MFRAC(K)*MCOL
!   MFIX(K) = MFIX0(K)-coeff*MTOP  ! subtract during initialization only

  ! edge sigmas as if entire atmosphere were sigma (for psige defined above)
  real*8, dimension(lm+1), parameter, private :: sige_=(plbot-plbot(lm+1))/(plbot(1)-plbot(lm+1))

  ! pf0 (hPa) nominal midpoint pressure of transition from troposphere to
  !    stratosphere in wtfix formula below.
  ! pfw (hPa) nominal width of transition in wtfix formula below.  Note that the
  !    effective pfw is smaller on the stratospheric side of the transition.

  real*8, parameter, private :: &
       pf0=100d0,pfw=70d0
       ! pf0=450d0,pfw=225d0

  ! wtfix is the weight for plbot in the pe(k) expression above; it asymptotes to
  ! unity in the stratosphere and zero at the surface.
  ! The shape of wtfix is an engineering choice and is not fundamental to STDHYB.
  real*8, dimension(lm+1), parameter, private :: wtfix0 = &
       .5d0*(1d0+tanh( &
       (pf0-plbot)/(pfw*(1d0-max(0d0,(pf0-plbot)/(1.3d0*pf0)))) &
       ) )

  ! rescale to span 0 to 1
  real*8, dimension(lm+1), parameter, private :: &
       wtfix = (wtfix0-wtfix0(1))/(wtfix0(lm+1)-wtfix0(1))

  real*8, dimension(lm+1), parameter, private :: pwtfix = plbot*wtfix

  real*8, dimension(lm+1), parameter, private :: cvar = (1d0-wtfix)*sige_

  real*8, dimension(lm), parameter, private :: dcvar = cvar(1:lm)-cvar(2:lm+1)

  real*8, dimension(lm), parameter, private :: &
       dpfix = (pwtfix(1:lm) - pwtfix(2:lm+1)) &
       - plbot(lm+1)*(dcvar(1:lm) + (wtfix(1:lm)-wtfix(2:lm+1)))

  real*8, dimension(lm), parameter :: mfix = dpfix*mb2kg, mfrac=dcvar

! The discontinuous-transition sigma-and-CP version in the non-STDHYB block can be
! expressed in STDHYB form; here is some commented-out code to do so if STDHYB ever
! becomes the default. If-test logic can be used instead of wt1 when a verticalres
! init routine is introduced.

!!@var varmask unity from layers 1-ls1-1, zero above
!!@var fixmask zero in layers 1-ls1-1, unity above
!  real*8, dimension(lm), parameter, private :: &
!    varmask = (/ (1d0, iii=1,layering%ls1-1), (0d0, iii=layering%ls1,lm) /), &
!    fixmask = 1d0-varmask
!
!  real*8, dimension(lm), parameter, private :: &
!      mfrac1 = varmask*delp/(plbot(1)-plbot(layering%ls1)), &
!      mfix1 = mb2kg*(-plbot(layering%ls1)*mfrac1 + delp*fixmask)
!
!  real*8, dimension(lm), parameter :: mfix2 = dpfix*mb2kg, mfrac2=dcvar
!
! If ls1 is set to lm+1, mfix,mfrac will be mfix1,mfrac1
!  real*8, parameter :: wt1 = 1d0 - layering%ls1/(lm+1)
!  real*8, parameter :: &
!       mfix(lm)  = wt1*mfix1 + (1d0-wt1)*mfix2
!  real*8, parameter :: &
!       mfrac(lm) = wt1*mfrac1 + (1d0-wt1)*mfrac2

#endif /* STDHYB or not */

End Module VerticalRes
