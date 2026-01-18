#include "rundeck_opts.h"

      Module VerticalRes
!**** AtmL40A.F90     Vertical Resolution file, 40 layers, top at .1 mb
      Implicit None

!****                                  ----MSURF=10026----     ----MSURF=5674-----
!**** Layer   MFIX   MFIXS   MFRAC     MVAR     MA     SUM     MVAR     AM     SUM
!**** =====   ====   =====   =====     ====     ==     ===     ====     ==     ===
!****   Top                                              1                       1

!****    40      1      1      0          0      1       2        0      1       2
!****    39      4      5      0          0      4       6        0      4       6
!****    38      9     14      0          0      9      15        0      9      15
!****    37     16     30      0          0     16      31        0     16      31
!****    36     25     55      0          0     25      56        0     25      56
!****    35     36     91      0          0     36      92        0     36      92
!****    34     49    140      0          0     49     141        0     49     141
!****    33     64    204      0          0     64     205        0     64     205
!****    32     81    285      0          0     81     286        0     81     286
!****    31    100    385      0          0    100     386        0    100     386

!****    30    120    505      0          0    120     506        0    120     506
!****    29    140    645      0          0    140     646        0    140     646
!****    28    160    805      0          0    160     806        0    160     806
!****    27    180    985      0          0    180     986        0    180     986
!****    26    200   1185    0/256        0    200    1186        0    200    1186
!****    25    200   1385    1/256       15    215    1401       -2    198    1384
!****    24    200   1585    3/256       45    245    1646       -6    194    1578
!****    23    200   1785    6/256       90    290    1936      -12    188    1766
!****    22    290   1985   10/256      150    350    2286      -20    180    1946
!****    21    200   2185   13/256      195    395    2681      -26    174    2120

!****    20    200   2385   15/256      225    425    3106      -30    170    2290
!****    19    200   2585   16/256      240    440    3546      -32    168    2458
!****    18    200   2785   16/256      240    440    3986      -32    168    2626
!****    17    200   2985   16/256      240    440    4426      -32    168    2794
!****    16    200   3185   16/256      240    440    4866      -32    168    2962
!****    15    200   3385   16/256      240    440    5306      -32    168    3130
!****    14    200   3585   16/256      240    440    5746      -32    168    3298
!****    13    200   3785   16/256      240    440    6186      -32    168    3466
!****    12    200   3985   16/256      240    440    6626      -32    168    3634
!****    11    200   4185   16/256      240    440    7066      -32    168    3802

!****    10    200   4385   16/256      240    440    7506      -32    168    3970
!****     9    200   4585   15/256      225    425    7931      -30    170    4140
!****     8    200   4785   13/256      195    395    8326      -26    174    4314
!****     7    200   4985   10/256      150    350    8676      -20    180    4494
!****     6    200   5185    6/256       90    290    8966      -12    188    4682
!****     5    200   5385    3/256       45    245    9211       -6    194    4876
!****     4    200   5585    1/256       15    215    9426       -2    198    5074
!****     3    200   5785    0/256        0    200    9626        0    200    5274
!****     2    200   5985      0          9    200    9826        0    200    5474
!****     1    200   6185      0          0    200   10026        0    200    5674
!****         ----         -------     ----  -----             ----   ----
!****         6185         256/256     3840  10025             -512   5673

!@var LM    = number of dynamical layers
!@var LS1   = lowest layer from top with MFRAC(LS1) = 0
!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV
!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var MA(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)      = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L)     = fraction of variable mass in each layer = DSIG(L)
!@var MFIXSOLD(L)  = Sum [MFIX(L+1:LM)]    defined from 0:LM
!@var MFRACSOLD(L) = Sum [MFRAC(L+1:LM)]   defined from 0:LM

      Integer :: LLL
      Integer,Parameter :: LM=40, LS1=26, LS1_NOMINAL=LS1
      Real*8,Parameter  :: MDRYA = 10034, MTOP = 1, &
         MFIXSOLD(0:LM) = (/ 6185,  5985,5785,5585,5385,5185,  4985,4785,4585,4385,4185, &
                                    3985,3785,3585,3385,3185,  2985,2785,2585,2385,2185, &
                                    1985,1785,1585,1385,1195,   985, 805, 645, 505, 385, &
                                     285, 204, 140,  91,  55,    30,  14,   5,   1,   0 /), &
         MFRACSOLD(0:LM) = (/ 1d0,  1d0,       1d0,       256/256d0, 255/256d0, 252/256d0, &
                                    246/256d0, 236/256d0, 223/256d0, 208/256d0, 192/256d0, &
                                    176/256d0, 160/256d0, 144/256d0, 128/256d0, 112/256d0, &
                                     96/256d0,  80/256d0,  64/256d0,  48/256d0,  33/256d0, &
                                     20/256d0,  10/256d0,   4/256d0,   1/256d0,   0/256d0, &
                                     0d0,0d0,0d0,0d0,0d0,  0d0,0d0,0d0,0d0,0d0,  0d0,0d0,0d0,0d0,0d0 /),&
         MFIXs    = MFIXSOLD(0), &
         MFIX(LM) = (/ ( MFIXSOLD(LLL-1) -  MFIXSOLD(LLL), LLL=1,LM) /), &
        MFRAC(LM) = (/ (MFRACSOLD(LLL-1) - MFRACSOLD(LLL), LLL=1,LM) /)

!@var MVARSMEAN = MDRYA - MFIXs - MTOP
!@var PSF       = global mean surface pressure (mb) = MDRYA * GRAV / 100
!@var PMTOP     = pressure of dynamical top (mb) = MTOP * GRAV / 100
!@var PTOP      = pressure at bottom of layer LS1 (mb) =
!@var PSTRAT    = pressure of stratosphere (mb) = PTOP - PMTOP
!@var PSFMPT    = mean pressure of troposhere (mb) = PSF - PTOP
!@var PLBOT(L)  = mean pressure at bottom of layer L (mb) = [MTOP + MFIXS(L-1) + MFRACS(L-1)*MVARSMEAN]*GRAV / 100

      Real*8,Parameter :: PSF = MDRYA*9.80665d0 / 100, &
                        PMTOP =  MTOP*9.80665d0 / 100, &
                         PTOP = (MTOP + MFIXSOLD(LS1-1))*9.80665d0 / 100, &
                       PSTRAT = PTOP - PMTOP, &
                       PSFMPT = PSF - PTOP, &
                    MVARSMEAN = MDRYA - MFIXs - MTOP, &
                PLBOT(1:LM+1) = (/ ((MTOP + MFIXSOLD(LLL) + MFRACSOLD(LLL)*MVARSMEAN)*9.80665d0 / 100, LLL=0,LM) /)
      EndModule VerticalRes
