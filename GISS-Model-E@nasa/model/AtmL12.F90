module VerticalRes
!@sum Vertical Resolution file 12 layers, top at 10 mb
!@auth Original Development Team
  Implicit None
!@var LM    = number of dynamical layers
!@var LS1   = lowest layer of strtosphere
  Integer*4,Parameter :: LM=12, LS1=9, LS1_NOMINAL=LS1

!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV
!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var AM(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)  = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L) = fraction of variable mass in each layer = DSIG(L)
  Real*8,Parameter :: MDRYA = 98400/9.80665d0, MTOP = 1000/9.80665d0, MFIXs = 14000/9.80665d0, &
    MFIX(LM) = (/ 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0, &
    5000/9.80665d0, 4000/9.80665d0, 3000/9.80665d0, 2000/9.80665d0 /), &
    MFRAC(LM) = (/ 50/834d0, 80/834d0, 134/834d0, 170/834d0, 160/834d0, &
    105/834d0, 75/834d0,  60/834d0, 0d0,0d0,0d0,0d0 /)

!**** Vertival resolution
!****                         ---MSURF=10034.0---    ---MSURF=5781.8----
!**** Layer   MFIX   MFRAC    MVAR     AM     SUM    MVAR     AM     SUM
!**** =====   ====   =====    ====     ==     ===    ====     ==     ===
!****   Top                                 102.0                  102.0
!****    12  203.9    0        0.0  203.9             0.0  203.9
!****    11  305.9    0        0.0  305.9             0.0  305.9
!****    10  407.9    0        0.0  407.9             0.0  407.9
!****     9  509.9    0        0.0  509.9  1529.6     0.0  509.9  1529.6
!****     8    0.0  60/834   611.8  611.8           305.9  305.9
!****     7    0.0  75/834
!****
!****     3    0.0 134/834
!****     2    0.0  80/834
!****     1    0.0  50/834   509.9  509.9 10034.0   102.0  102.0  5781.8
!****         ----  ------   -----  -----           -----  -----
!****       1427.6 834/834  8504.4 9932.0          4252.2 5679.8  

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
!@var PLbot pressure levels at bottom of layers (mb)
  Real*8,Parameter :: PSF=984.d0, PTOP = 150.d0, PMTOP = 10d0, PSFMPT = PSF-PTOP, &
    PSTRAT = PTOP-PMTOP, &
    PLBOT(1:LM+1) = (/ PSF, 934d0, 854d0, 720d0, 550.d0, &  ! Pbot L=1,5
    390d0, 285d0, 210d0,                &  !      L=...
    PTOP,                              &  !      L=LS1
    100d0,  60d0,  30d0, PMTOP /)          !      L=..,LM+1

End Module VerticalRes


