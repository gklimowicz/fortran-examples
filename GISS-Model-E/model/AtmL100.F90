#include "rundeck_opts.h"                                                            

      Module VerticalRes
!**** AtmL100.F90   Atmospheric Vertical Resolution file, 100 layers, dynamical top at .02 (kg/m^2)                        
      Implicit None                                                                      

!****  L     L                                      -----MSURF=9952.00------      -----MSURF=5830.00------
!**** New   Old      MFIX       MFIXS   MFRAC       MVAR       MA        SUM      MVAR       MA        SUM
!**** ===   ===      ====       =====   =====       ====       ==        ===      ====       ==        ===
!****   0   101                                                      41/2048                         .0200

!****   1   100     35/2048     35/2048   0           0      .0171     .0371        0      .0171     .0371
!****   2    99     66/2048    101/2048   0           0      .0322     .0693        0      .0322     .0693
!****   3    98    122/2048    223/2048   0           0      .0596     .1289        0      .0596     .1289
!****   4    97    228/2048    451/2048   0           0      .1113     .2402        0      .1113     .2402
!****   5    96    424/2048    875/2048   0           0      .2070     .4473        0      .2070     .4473
!****   6    95    790/2048   1665/2048   0           0      .3857     .8330        0      .3857     .8330
!****   7    94   1469/2048   3134/2048   0           0      .7173    1.5503        0      .7173    1.5503
!****   8    93   2735/2048   5869/2048   0           0     1.3354    2.8857        0     1.3354    2.8857
!****   9    92   5092/2048  10961/2048   0           0     2.4863    5.3721        0     2.4863    5.3721
!****  10    91   9478/2048  20439/2048   0           0     4.6279   10.0000        0     4.6279   10.0000

!****  11    90      6          15.98     0           0       6.00     16.00        0       6.00     16.00
!****  12    89      7.00       22.98     0           0       7.00     23.00        0       7.00     23.00
!****  13    88      8.00       30.98     0           0       8.00     31.00        0       8.00     31.00
!****  14    87      9.00       39.98     0           0       9.00     40.00        0       9.00     40.00
!****  15    86     10.00       49.98     0           0      10.00     50.00        0      10.00     50.00
!****  16    85     12.00       61.98     0           0      12.00     62.00        0      12.00     62.00
!****  17    84     14.00       75.98     0           0      14.00     76.00        0      14.00     76.00
!****  18    83     16.00       91.98     0           0      16.00     92.00        0      16.00     92.00
!****  19    82     18.00      109.98     0           0      18.00    110.00        0      18.00    110.00
!****  20    81     20.00      129.98     0           0      20.00    130.00        0      20.00    130.00

!****  21    80     22.00      151.98     0           0      22.00    152.00        0      22.00    152.00
!****  22    79     24.00      175.98     0           0      24.00    176.00        0      24.00    176.00
!****  23    78     26.00      201.98     0           0      26.00    202.00        0      26.00    202.00
!****  24    77     28.00      229.98     0           0      28.00    230.00        0      28.00    230.00
!****  25    76     30.00      259.98     0           0      30.00    260.00        0      30.00    260.00
!****  26    75     32.00      291.98     0           0      32.00    292.00        0      32.00    292.00
!****  27    74     34.00      325.98     0           0      34.00    326.00        0      34.00    326.00
!****  28    73     36.00      361.98     0           0      36.00    362.00        0      36.00    362.00
!****  29    72     38.00      399.98     0           0      38.00    400.00        0      38.00    400.00
!****  30    71     40.00      439.98     0           0      40.00    440.00        0      40.00    440.00

!****  31    70     42.00      481.98     0           0      42.00    482.00        0      42.00    482.00
!****  32    69     44.00      525.98     0           0      44.00    526.00        0      44.00    526.00
!****  33    68     46.00      571.98     0           0      46.00    572.00        0      46.00    572.00
!****  34    67     48.00      619.98     0           0      48.00    620.00        0      48.00    620.00
!****  35    66     50.00      669.98     0           0      50.00    670.00        0      50.00    670.00
!****  36    65     52.00      721.98     0           0      52.00    722.00        0      52.00    722.00
!****  37    64     54.00      775.98     0           0      54.00    776.00        0      54.00    776.00
!****  38    63     56.00      831.98     0           0      56.00    832.00        0      56.00    832.00
!****  39    62     58.00      889.98     0           0      58.00    890.00        0      58.00    890.00
!****  40    61     60.00      949.98     0           0      60.00    950.00        0      60.00    950.00

!****  41    60     62.00     1011.98     0           0      62.00   1012.00        0      62.00   1012.00
!****  42    59     64.00     1075.98    0/2048       0      64.00   1076.00        0      64.00   1076.00
!****  43    58     66.00     1141.98    1/2048      1.50    67.50   1143.50       -.50    65.50   1141.50
!****  44    57     68.00     1209.98    3/2048      4.50    72.50   1216.00      -1.50    66.50   1208.00
!****  45    56     70.00     1279.98    6/2048      9.00    79.00   1295.00      -3.00    67.00   1275.00
!****  46    55     72.00     1351.98   10/2048     15.00    87.00   1382.00      -5.00    67.00   1342.00
!****  47    54     74.00     1425.98   15/2048     22.50    96.50   1478.50      -7.50    66.50   1408.50
!****  48    53     76.00     1501.98   21/2048     31.50   107.50   1586.00     -10.50    65.50   1474.00
!****  49    52     78.00     1579.98   28/2048     42.00   120.00   1706.00     -14.00    64.00   1538.00
!****  50    51     80.00     1659.98   36/2048     54.00   134.00   1840.00     -18.00    62.00   1600.00

!****  51    50     81.00     1740.98   43/2048     64.50   145.50   1985.50     -21.50    59.50   1659.50
!****  52    49     82.00     1822.98   49/2048     73.50   155.50   2141.00     -24.50    57.50   1717.00
!****  53    48     83.00     1905.98   54/2048     81.00   164.00   2305.00     -27.00    56.00   1773.00
!****  54    47     84.00     1989.98   58/2048     87.00   171.00   2476.00     -29.00    55.00   1828.00
!****  55    46     85.00     2074.98   61/2048     91.50   176.50   2652.50     -30.50    54.50   1882.50
!****  56    45     86.00     2160.98   63/2048     94.50   180.50   2833.00     -31.50    54.50   1937.00
!****  57    44     87.00     2247.98   64/2048     96.00   183.00   3016.00     -32.00    55.00  1992.00
!****  58    43     88.00     2335.98   64/2048     96.00   184.00   3200.00     -32.00    56.00   2048.00
!****  59    42     89.00     2424.98   64/2048     96.00   185.00   3385.00     -32.00    57.00   2105.00
!****  60    41     90.00     2514.98   64/2048     96.00   186.00   3571.00     -32.00    58.00   2163.00

!****  61    40     91.00     2605.98   64/2048     96.00   187.00   3758.00     -32.00    59.00   2222.00
!****  62    39     92.00     2697.98   64/2048     96.00   188.00   3946.00     -32.00    60.00   2282.00
!****  63    38     93.00     2790.98   64/2048     96.00   189.00   4135.00     -32.00    61.00   2343.00
!****  64    37     94.00     2884.98   64/2048     96.00   190.00   4325.00     -32.00    62.00   2405.00
!****  65    36     95.00     2979.98   64/2048     96.00   191.00   4516.00     -32.00    63.00   2468.00
!****  66    35     96.00     3075.98   64/2048     96.00   192.00   4708.00     -32.00    64.00   2532.00
!****  67    34     97.00     3172.98   64/2048     96.00   193.00   4901.00     -32.00    65.00   2597.00
!****  68    33     98.00     3270.98   64/2048     96.00   194.00   5095.00     -32.00    66.00   2663.00
!****  69    32     99.00     3369.98   64/2048     96.00   195.00   5290.00     -32.00    67.00   2730.00
!****  70    31    100.00     3469.98   64/2048     96.00   196.00   5486.00     -32.00    68.00   2798.00

!****  71    30    101.00     3570.98   64/2048     96.00   197.00   5683.00     -32.00    69.00   2867.00
!****  72    29    102.00     3672.98   64/2048     96.00   198.00   5881.00     -32.00    70.00   2937.00
!****  73    28    103.00     3775.98   64/2048     96.00   199.00   6080.00     -32.00    71.00   3008.00
!****  74    27    104.00     3879.98   64/2048     96.00   200.00   6280.00     -32.00    72.00   3080.00
!****  75    26    105.00     3984.98   63/2048     94.50   199.50   6479.50     -31.50    73.50   3153.50
!****  76    25    106.00     4090.98   61/2048     91.50   197.50   6677.00     -30.50    75.50   3229.00
!****  77    24    107.00     4197.98   58/2048     87.00   194.00   6871.00     -29.00    78.00   3307.00
!****  78    23    108.00     4305.98   54/2048     81.00   189.00   7060.00     -27.00    81.00   3388.00
!****  79    22    109.00     4414.98   49/2048     73.50   182.50   7342.50     -24.50    84.50   3472.50
!****  80    21    110.00     4524.98   43/2048     64.50   174.50   7417.00     -21.50    88.50   3561.00

!****  81    20    111.00     4635.98   36/2048     54.00   165.00   7582.00     -18.00    93.00   3654.00
!****  82    19    112.00     4747.98   28/2048     42.00   154.00   7736.00     -14.00    98.00   3752.00
!****  83    18    113.00     4860.98   21/2048     31.50   144.50   7880.50     -10.50   102.50   3854.50
!****  84    17    114.00     4974.98   15/2048     22.50   136.50   8017.00      -7.50   106.50   3961.00
!****  85    16    115.00     5089.98   10/2048     15.00   130.00   8147.00      -5.00   110.00   4071.00
!****  86    15    116.00     5205.98    6/2048      9.00   125.00   8272.00      -3.00   113.00   4184.00
!****  87    14    117.00     5322.98    3/2048      4.50   121.50   8393.50      -1.50   115.50   4299.50
!****  88    13    118.00     5440.98    1/2048      1.50   119.50   8513.00       -.50   117.50   4417.00
!****  89    12    119.00     5559.98    0/2048       0     119.00   8632.00        0     119.00   4536.00
!****  90    11    120.00     5679.98     0           0     120.00   8752.00        0     120.00   4656.00

!****  91    10    120.00     5799.98     0           0     120.00   8872.00        0     120.00   4776.00
!****  92     9    120.00     5919.98     0           0     120.00   8992.00        0     120.00   4896.00
!****  93     8    120.00     6039.98     0           0     120.00   9112.00        0     120.00   5016.00
!****  94     7    120.00     6159.98     0           0     120.00   9232.00        0     120.00   5136.00
!****  95     6    120.00     6279.98     0           0     120.00   9352.00        0     120.00   5256.00
!****  96     5    120.00     6399.98     0           0     120.00   9472.00        0     120.00   5376.00
!****  97     4    120.00     6519.98     0           0     120.00   9592.00        0     120.00   5496.00
!****  98     3    120.00     6639.98     0           0     120.00   9712.00        0     120.00   5616.00
!****  99     2    120.00     6759.98     0           0     120.00   9832.00        0     120.00   5736.00
!**** 100     1    120.00     6879.98     0           0     120.00   9952.00        0     120.00   5856.00
!****              ------             --------    -------   ------              -------   ------
!****             6879.98             2048/2048   3072.00  9951.98             -1024.00  5855.98  

!**** In this section new layering is used where Lnew is counted from the top down
!@var LM        = number of dynamical layers
!@var MTOP      = mass above dynamical top = .02 (kg/m^2)
!@var MFIX(L)   = fixed mass of each layer (kg/m^2)
!@var MFIXS(L)  = Sum [MFIX(1:L)]
!@var MFRAC(L)  = fixed fraction of variable mass in each layer (dimensionless) = MVAR(L) / MVARS
!@var MFRACS(L) = Sum [MFRAC(1:L)] ;   Sum [MFRAC(1:LM)] = MFRACS(LM) = 1  
!@var MVARS     = Sum [MVAR(1:LM)] = spatially and temporally varying column mass (kg/m^2), may be negative
!@var MVAR(L)   = MVARS * MFRAC(L) = spatially and temporally varying mass of each layer (kg/m^2), may be negative
!@var MSURF     = MTOP + MFIXS(LM) + MVARS = spatially and temporally varying total column mass (kg/m^2)
!@var MA(L)     = MFIX(L) + MVARS*MFRAC(L) = mass of each layer (kg/m^2)
!@var MDRYA     = global mean dry atmospheric mass = 10034.00 (kg/m^2), water vapor has no mass at present

!**** Using exponential layering, top 10 mass thicknesses plus MTOP sum to 10 (kg/m^2).
!**** Counting from top down, mass of each layer is A * X^L where A, X are constants with A > 0, X > 1.
!**** (1):  MTOP = 41/2048 = A * Sum[X^(-INFINITY:0)] = A * [X^1 - X^(-INFINITY)] / (X-1) = A * X / (X-1)
!**** (2):  MTOP + MFIXS(10) = 10 = A * Sum[X^(-INFINITY:10)] = A * X^11 / (X-1)
!**** (2) / (1):  10 * 2048/41 = 20480/41 = X^10  =>  X = (20480/41)^(1/10)  ;  (1):  A = (41/2048) * (X-1) / X
!**** MFIX(L)  = A * X^L = (41/2048) * (X-1) * X^(L-1)
!**** MFIXS(L) = A * X^(L+1) / (X-1) - MTOP = (41/2048) * X^L - (41/2048) 

      Integer :: LLL
      Integer,Parameter :: LM = 100
      Real*8, Parameter :: MDRYA = 10034, MTOP = 41/2048d0, XXX = (10/MTOP)**(1/10d0), & 
         MFIXSNEW(0:LM) = (/ 0d0, 35/2048d0, 101/2048d0, 223/2048d0,  451/2048d0,  875/2048d0, &
                                1665/2048d0,3134/2048d0,5869/2048d0,10961/2048d0,20439/2048d0, &
                             ((LLL-5 + (LLL-5)**2)/2 - 5 - MTOP, LLL=11,15), &
                             (LLL-10 + (LLL-10)**2 + 20 - MTOP, LLL=16,50), &
                             ((LLL+30 + (LLL+30)**2)/2 - 1580 - MTOP, LLL=51,90), &
                             (LLL*120 - 5120 - MTOP, LLL=91,100) /), &
         MFRACSNEW(0:LM) = (/ (0d0, LLL=0,41), &
                         0/2048d0,   1/2048d0,   4/2048d0,  10/2048d0,  20/2048d0,  35/2048d0,  56/2048d0,  84/2048d0, 120/2048d0, &
           163/2048d0, 212/2048d0, 266/2048d0, 324/2048d0, 385/2048d0, 448/2048d0, 512/2048d0, 576/2048d0, 640/2048d0, 704/2048d0, &
           768/2048d0, 832/2048d0, 896/2048d0, 960/2048d0,1024/2048d0,1088/2048d0,1152/2048d0,1216/2048d0,1280/2048d0,1344/2048d0, &
          1408/2048d0,1472/2048d0,1536/2048d0,1600/2048d0,1663/2048d0,1724/2048d0,1782/2048d0,1836/2048d0,1885/2048d0,1928/2048d0, &
          1964/2048d0,1992/2048d0,2013/2048d0,2028/2048d0,2038/2048d0,2044/2048d0,2047/2048d0,2048/2048d0,2048/2048d0, &
                              (1d0, LLL=90,100) /), &
!        MFIXNEW(LM) = (/ 35/2048d0,  66/2048d0, 122/2048d0, 228/2048d0, 424/2048d0, &
!                        790/2048d0,1469/2048d0,2735/2048d0,5092/2048d0,9478/2048d0, &
!                        (LLL-5d0, LLL=11,15), (2*(LLL-10d0), LLL=16,50), (LLL+30d0, LLL=51,90), (120d0, LLL=91,100) /), &
!         MFRACNEW(LM) = (/ (0d0, LLL=1,41), &
!                            0/2048d0, 1/2048d0, 3/2048d0, 6/2048d0,10/2048d0,15/2048d0,21/2048d0,28/2048d0,36/2048d0, &
!                 43/2048d0,49/2048d0,54/2048d0,58/2048d0,61/2048d0,63/2048d0,64/2048d0,64/2048d0,64/2048d0,64/2048d0, &
!                 64/2048d0,64/2048d0,64/2048d0,64/2048d0,64/2048d0,64/2048d0,64/2048d0,64/2048d0,64/2048d0,64/2048d0, &
!                 64/2048d0,64/2048d0,64/2048d0,64/2048d0,63/2048d0,61/2048d0,58/2048d0,54/2048d0,49/2048d0,43/2048d0, &
!                 36/2048d0,28/2048d0,21/2048d0,15/2048d0,10/2048d0, 6/2048d0, 3/2048d0, 1/2048d0, 0/2048d0, &
!                           (0d0, LLL=90,100) /), &
           MFIXNEW(LM) = (/ ( MFIXSNEW(LLL) -  MFIXSNEW(LLL-1), LLL=1,LM) /), &
          MFRACNEW(LM) = (/ (MFRACSNEW(LLL) - MFRACSNEW(LLL-1), LLL=1,LM) /)

!**** In this section old layering is used where Lold is counted from the bottom up 
      Real*8,Parameter :: &
          MFIXs       =      MFIXSNEW(LM), &
!         MFIXS(0:LM) = (/ ( MFIXSNEW(LM  -LLL), LLL=0,LM) /), &  !  note that MFIXS  is defined 0:LM, not 1:LM+1
         MFRACS(0:LM) = (/ (MFRACSNEW(LM  -LLL), LLL=0,LM) /), &  !  note that MFRACS is defined 0:LM, not 1:LM+1
             MFIX(LM) = (/ (  MFIXNEW(LM+1-LLL), LLL=1,LM) /), &
            MFRAC(LM) = (/ ( MFRACNEW(LM+1-LLL), LLL=1,LM) /) 

!@var LS1 = lowest layer of strtosphere = layer of lowest elevation with MFRAC(LS1:LM) = 0
!@var MFIX(L)   = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L)  = fraction of variable mass in each layer = DSIG(L)
!@var MVARSMEAN = MDRYA - MFIXS(LM) - MTOP
!@var PSF    = global mean surface pressure (mb) = MDRYA * GRAV / 100 
!@var PMTOP  = pressure of dynamical top (mb) = MTOP * GRAV / 100
!@var PTOP   = pressure at bottom of layer LS1 (mb) = 
!@var PSTRAT = pressure of stratosphere (mb) = PTOP - PMTOP
!@var PSFMPT = mean pressure of troposhere (mb) = PSF - PTOP
!@var PLBOT(L) = mean pressure at bottom of layer L (mb) = [MTOP + MFIXS(L-1) + MFRACS(L-1)*MVARSMEAN]*GRAV / 100
      Integer,Parameter :: LS1 = 59, LS1_NOMINAL=LS1
      Real*8, Parameter :: PSF = MDRYA*9.80665d0 / 100, &
                         PMTOP =  MTOP*9.80665d0 / 100, &
                          PTOP = (MTOP + MFIXSNEW(LM-LS1+1))*9.80665d0 / 100, &
!                         PTOP = (MTOP + MFIXS(LS1-1))*9.80665d0 / 100, &
                        PSTRAT = PTOP - PMTOP, &
                        PSFMPT = PSF - PTOP, &
                     MVARSMEAN = MDRYA - MFIXs - MTOP, &
!                    MVARSMEAN = MDRYA - MFIXS(0) - MTOP, &
                 PLBOT(1:LM+1) = (/ ((MTOP + MFIXSNEW(LM-LLL) + MFRACS(LLL)*MVARSMEAN)*9.80665d0 / 100, LLL=0,LM) /)   
!                PLBOT(1:LM+1) = (/ ((MTOP + MFIXS(LLL) + MFRACS(LLL)*MVARSMEAN)*9.80665d0 / 100, LLL=0,LM) /)   
     EndModule VerticalRes

