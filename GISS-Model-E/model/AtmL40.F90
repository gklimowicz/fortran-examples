#include "rundeck_opts.h"
module VerticalRes
!@sum Vertical Resolution file, 40 layers, top at .1 mb
!@auth Original Development Team
  Implicit None
!@var LM    = number of dynamical layers
!@var LS1   = lowest layer of strtosphere
  Integer*4,Parameter :: LM=40, LS1=24, LS1_NOMINAL=LS1

!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV
!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var AM(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)  = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L) = fraction of variable mass in each layer = DSIG(L)
  Real*8,Parameter :: MDRYA = 98400/9.80665d0, MTOP = 10/9.80665d0, MFIXs = 14990/9.80665d0, &
    MFIX(LM) = (/ 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0,  &
    0d0,0d0,0d0, &
    2200/9.80665d0, 2000/9.80665d0, 1800/9.80665d0, 1700/9.80665d0, 1600/9.80665d0, &
    1400/9.80665d0, 1200/9.80665d0, 1100/9.80665d0, 1000/9.80665d0,  438/9.80665d0, &
    246/9.80665d0,  138/9.80665d0,   78/9.80665d0, 43.8d0/9.80665d0, 24.6d0/9.80665d0, &
    13.8d0/9.80665d0,7.8d0/9.80665d0 /), &
    MFRAC(LM) = (/ 20/834d0, 22/834d0, 25/834d0, 27/834d0, 30/834d0, &
    35/834d0, 40/834d0, 45/834d0, 48/834d0, 50/834d0, &
    51/834d0, 52/834d0, 50/834d0, 48/834d0, 45/834d0, &
    42/834d0, 38/834d0, 34/834d0, 31/834d0, 28/834d0, &
    26/834d0, 24/834d0, 23/834d0, &   
    0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0 /)

!**** Vertival resolution
!****                         ---MSURF=10034.0---    ---MSURF=5781.8----
!**** Layer   MFIX   MFRAC    MVAR     AM     SUM    MVAR     AM     SUM
!**** =====   ====   =====    ====     ==     ===    ====     ==     ===
!****   Top                                   1.0                    1.0
!****    40     .8    0        0.0     .8     1.8     0.0     .8     1.8
!****    39    1.4    0        0.0    1.4     3.2     0.0    1.4     3.2
!****    38    2.5    0        0.0    2.5     5.7     0.0    2.5     5.7
!****    37    4.5    0        0.0    4.5    10.2     0.0    4.5    10.2
!****    36    8.0    0        0.0    8.0    18.2     0.0    8.0    18.2
!****    35   14.1    0        0.0   14.1    32.2     0.0   14.1    32.2
!****    34   25.1    0        0.0   25.1    57.3     0.0   25.1    57.3
!****    33   44.7    0        0.0   44.7   102.0     0.0   44.7   102.0
!****    32  102.0    0        0.0  102.0   203.9     0.0  102.0   203.9
!****    31  112.2    0        0.0  112.2   316.1     0.0  112.2   316.1
!****
!****    30  122.4    0        0.0  122.4   438.5     0.0  122.4   438.5
!****    29  142.8    0        0.0  142.8   581.2     0.0  142.8   581.2
!****    28  163.2    0        0.0  163.2   744.4     0.0  163.2   744.4
!****    27  173.4    0        0.0  173.4   917.7     0.0  173.4   917.7
!****    26  183.5    0        0.0  183.5  1101.3     0.0  183.5  1101.3
!****    25  203.9    0        0.0  203.9  1305.2     0.0  203.9  1305.2
!****    24  224.3    0        0.0  224.3  1529.6     0.0  224.3  1529.6
!****    23    0.0  23/834   234.5  234.5  1764.1   117.3  117.3  1646.8
!****    22    0.0  24/834   244.7  244.7  2008.8   122.4  122.4  1769.2
!****
!****    10    0.0  50/834   509.9  509.9  7545.9   254.9  254.9  4293.0
!****     9    0.0  48/834   489.5  489.5  8004.8   244.7  244.7  4537.7
!****     8    0.0  45/834   458.9  458.9  8412.7   229.4  229.4  4767.2
!****     7    0.0  40/834   407.9  407.9  8769.6   203.9  203.9  4971.1
!****     6    0.0  35/834   356.9  356.9  9075.5   178.5  178.5  5149.6
!****     5    0.0  30/834   305.9  305.9  9605.7   153.0  153.0  5302.5
!****     4    0.0  27/834   275.3  275.3  9350.8   137.7  137.7  5440.2
!****     3    0.0  25/834   254.9  254.9  9605.7   127.5  127.5  5567.7
!****     2    0.0  22/834   224.3  224.3  9830.1   112.2  112.2  5679.8
!****     1    0.0  20/834   203.9  203.9 10034.0   102.0  102.0  5781.8
!****         ----  ------   -----  -----           -----  -----
!****       1528.6 834/834  8504.4 10033.0         4252.2 5780.8  

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
!@var PLbot pressure levels at bottom of layers (mb)
  Real*8,Parameter :: PSF=984.d0, PTOP = 150.d0, PMTOP = .1d0, PSFMPT = PSF-PTOP, &
    PSTRAT = PTOP-PMTOP, &
    PLBOT(1:LM+1) = (/ PSF,   964d0, 942d0, 917d0, 890d0, 860d0, 825d0, &  !  L=1,..   
    785d0, 740d0, 692d0, 642d0, 591d0, 539d0, 489d0, &  !  L=...    
    441d0, 396d0, 354d0, 316d0, 282d0, 251d0, 223d0, &
    197d0, 173d0,                                    &
    PTOP,                                            &  !  L=LS1    
    128d0, 108d0,  90d0,  73d0,  57d0,  43d0,  31d0, &  !  L=...    
     20d0,  10d0,5.62d0,3.16d0,1.78d0,  1.d0,        &
    .562d0,.316d0,.178d0, PMTOP /)                      !  L=..,LM+1

End Module VerticalRes
