module VerticalRes
!@sum Vertical Resolution file, 53 layers, top at .20576514 Pa
!@ver 2013/03/20
!@auth Original Development Team
  Implicit None
!@var LM    = number of dynamical layers
!@var LS1   = lowest layer of strtosphere
  Integer*4,Parameter :: LM=53, LS1=27, LS1_NOMINAL=LS1

!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV
!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var AM(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)  = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L) = fraction of variable mass in each layer = DSIG(L)
  Real*8,Parameter :: MDRYA = 98400/9.80665d0, MTOP = .20576514d0/9.80665d0, &
    MFIXs = 15400/9.80665d0-MTOP, &
    MFIX(LM) = (/ 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, &
    0d0,0d0,0d0,0d0,0d0, 0d0, &
    1900/9.80665d0, 1800/9.80665d0, 1700/9.80665d0, 1600/9.80665d0, 1500/9.80665d0, &
    1400/9.80665d0, 1300/9.80665d0, 1200/9.80665d0, 1100/9.80665d0,  900/9.80665d0, &
    289.57d0/9.80665d0, 205.719d0/9.80665d0, 146.148d0/9.80665d0, 103.829d0/9.80665d0, &
    73.763d0/9.80665d0,  52.404d0/9.80665d0, 37.2288d0/9.80665d0, 35.1388d0/9.80665d0, &
    24.5991d0/9.80665d0, 13.8008d0/9.80665d0,  7.7995d0/9.80665d0, 4.37955d0/9.80665d0, &
    2.46066d0/9.80665d0, 1.37946d0/9.80665d0,  .78033d0/9.80665d0, .481068d0/9.80665d0, &
    .31316686d0/9.80665d0 /), &
    MFRAC(LM) = (/ 30/830d0, 30/830d0, 30/830d0, 30/830d0, 30/830d0, 30/830d0, &
    30/830d0, 30/830d0, 30/830d0, 30/830d0, 34/830d0, 44/830d0, &
    45/830d0, 46/830d0, 45/830d0, 43/830d0, 40/830d0, 36/830d0, &
    33/830d0, 28/830d0, 26/830d0, 24/830d0, 23/830d0, 22/830d0, &
    21/830d0, 20/830d0, &
    0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, &
    0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0 /)

!**** Vertival resolution
!****                         ---MSURF=10034.0---    ---MSURF=5802.2----
!**** Layer   MFIX   MFRAC    MVAR     AM     SUM    MVAR     AM     SUM
!**** =====   ====   =====    ====     ==     ===    ====     ==     ===
!****   Top                                   .02                    .02
!****    53    .03    0        0.0    .03             0.0    .03
!****
!****    27  193.7    0        0.0  193.7  1570.4     0.0  193.7  1570.4
!****    26    0.0  20/830   203.9  203.9           102.0  102.0
!****
!****     1    0.0  30/830   305.9  305.9 10034.0   153.0  153.0  5802.2
!****         ----  ------   -----  -----           -----  -----
!****       1570.3 830/830  8463.6 10034.0         4231.8 5802.2  

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
!@var PLbot pressure levels at bottom of layers (mb)
  Real*8,Parameter :: PSF=984.d0, PTOP = 154.d0, PMTOP = .0020576514d0, PSFMPT = PSF-PTOP, &
    PSTRAT = PTOP-PMTOP, &
    PLBOT(1:LM+1) = (/ &
    PSF,   954d0,  924d0,  894d0,  864d0,  834d0,  804d0,  774d0, &
    744d0, 714d0,  680d0,  636d0,  591d0,  545d0,  500d0,  457d0, &
    417d0, 381d0,  348d0,  318d0,  290d0,  264d0,  240d0,  217d0, &
    195d0, 174d0, &                              ! Pbot L=1,26
    PTOP, &                                      !      L=LS1
    135d0, 117d0,  100d0,   84d0,   69d0,   55d0,   42d0,   30d0, &
    19d0,  10d0, 7.1043d0, 5.04711d0, 3.58563d0, 2.54734d0, &
    1.80971d0, 1.28567d0, 0.913382d0, 0.561994d0, 0.316003d0, &
    0.177995d0, 0.1d0, 0.0562045d0, 0.0315979d0, 0.0178033d0, &
    0.01d0, 0.00518932d0, PMTOP /)              !      L=..,LM+1

End Module VerticalRes

