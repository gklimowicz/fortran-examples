module VerticalRes
!@sum Vertical Resolution file, 23 layers, top at .20576514 Pa
!@ver 2013/03/21
!@auth Original Development Team
  Implicit None
!@var LM    = number of dynamical layers
!@var LS1   = lowest layer of strtosphere
  Integer*4,Parameter :: LM=23, LS1=12, LS1_NOMINAL=LS1

!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV
!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var AM(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)  = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L) = fraction of variable mass in each layer = DSIG(L)
  Real*8,Parameter :: MDRYA = 98400/9.80665d0, MTOP = .20576514d0/9.80665d0, MFIXs = 15000/9.80665d0-MTOP, &
    MFIX(LM) = (/ 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0, &
    3300/9.80665d0, 3080/9.80665d0, 3000/9.80665d0, 2460/9.80665d0, 1380/9.80665d0, &
    780/9.80665d0,  537/9.80665d0,  317/9.80665d0, 99.9004742d0/9.80665d0, &
    31.59957612/9.80665d0, 11.37992166d0/9.80665d0, 2.91426288d0/9.80665d0 /), &
    MFRAC(LM) = (/ 24/834d0, 31/834d0, 45/834d0, 65/834d0,109/834d0, 140/834d0, &
    145/834d0,111/834d0, 69/834d0, 53/834d0, 42/834d0, &
    0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0 /)

!**** Vertival resolution
!****                         ---MSURF=10034.0---    ---MSURF=5781.8----
!**** Layer   MFIX   MFRAC    MVAR     AM     SUM    MVAR     AM     SUM
!**** =====   ====   =====    ====     ==     ===    ====     ==     ===
!****   Top                                   .02                    .02
!****    23    .30    0        0.0    .30             0.0    .30
!****    22   1.16    0        0.0   1.16             0.0   1.16
!****
!****    12  336.5    0        0.0  336.5  1529.6     0.0  336.5  1529.6
!****    11    0.0  42/834   428.3  428.3           214.1  214.1
!****
!****     1    0.0  24/834   244.7  244.7 10034.0   122.4  122.4  5781.8
!****         ----  ------   -----  -----           -----  -----
!****       1528.6 834/834  8504.4 10033.0         4252.2 5780.8  

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
!@var PLbot pressure levels at bottom of layers (mb)
  Real*8,Parameter :: PSF=984.d0, PTOP = 150.d0, PMTOP = .0020576514d0, PSFMPT = PSF-PTOP, &
    PSTRAT = PTOP-PMTOP, &
    PLBOT(1:LM+1) = (/ PSF,    960d0,  929d0,  884d0,  819d0, &  ! Pbot L=1,5
    710d0,  570d0,  425d0,  314d0,  245d0, &  !      L=6,10
    192d0,                                 &  !      L=11
    PTOP,                                  &  !      L=LS1
    117d0, 86.2d0, 56.2d0, 31.6d0, 17.8d0, &  !      L=13,17
    10d0, 4.63d0, 1.46d0,  0.460995258d0, &  !      L=18,21
    .1449994968d0, .0312002802d0,  PMTOP /)   !      L=22,LM+1

End Module VerticalRes

