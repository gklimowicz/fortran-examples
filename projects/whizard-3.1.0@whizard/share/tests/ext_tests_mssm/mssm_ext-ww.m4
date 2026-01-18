dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-ww.m4 > mssm_ext-ww.sin
dnl   whizard -r mssm_ext-ww.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-ww.m4 -1   $, mssm_ext_ww_)
! --------------------------------------------------------------------------
model = MSSM
read_slha ("sps1a.slha")
?vis_history = false
! --------------------------------------------------------------------------
me = 0
mW = 80.419 
wW = 2.048
mZ = 91.1880
wZ = 2.446
mtop = 178.0
mb = 4.6
GF = 1.16639e-5
alphas = 0.118
show (al_h,mu_h,tanb_h)
! --------------------------------------------------------------------------
iterations = 3:200000
! --------------------------------------------------------------------------
PROCESS(`Wm,Wp', `se1,SE1'      ,  500 GeV,   192.145      ,  0.001     , 4)
PROCESS(`Wm,Wp', `se2,SE2'      ,  500 GeV,    14.2151     ,  0.0004    , 3)
PROCESS(`Wm,Wp', `smu1,SMU1'    ,  500 GeV,   192.146      ,  0.001     , 3)
PROCESS(`Wm,Wp', `smu2,SMU2'    ,  500 GeV,    14.2145     ,  0.0004    , 4)
PROCESS(`Wm,Wp', `stau1,STAU1'  ,  500 GeV,     7.9266     ,  0.0002    , 4)
dnl Rounding error in error estimate of WHIZARD 1 CATPISS number
PROCESS(`Wm,Wp', `stau2,STAU2'  ,  500 GeV,   168.046      ,  0.001     , 5)
PROCESS(`Wm,Wp', `stau1,STAU2'  ,  500 GeV,    17.8521     ,  0.0001    , 4)
PROCESS(`Wm,Wp', `sn11,SN11'    ,  500 GeV,   157.809      ,  0.003     , 3)
PROCESS(`Wm,Wp', `sn21,SN21'    ,  500 GeV,   157.806      ,  0.003     , 3)
PROCESS(`Wm,Wp', `sn31,SN31'    ,  500 GeV,   152.509      ,  0.003     , 3)
PROCESS(`Wm,Wp', `neu1,neu1'    ,  500 GeV,     3.8824     ,  0.0001    , 3)
PROCESS(`Wm,Wp', `neu1,neu2'    ,  500 GeV,   121.2925     ,  0.0007    , 4)
PROCESS(`Wm,Wp', `neu1,neu3'    ,  500 GeV,     6.8934     ,  0.0002    , 3)
PROCESS(`Wm,Wp', `neu1,neu4'    ,  500 GeV,     1.4973     ,  0.0006    , 3)
dnl Underestimated error in the WHIZARD 1 CATPISS value
PROCESS(`Wm,Wp', `neu2,neu2'    ,  500 GeV,  5996.57       ,  0.02      , 5)
PROCESS(`Wm,Wp', `"ch1+","ch1-"',  500 GeV,  3771.6        ,  0.4       , 3)
PROCESS(`Wm,Wp', `h,h'          ,  500 GeV,  6024.7        ,  0.4       , 3)
PROCESS(`Wm,Wp', `Z,h'          ,  500 GeV, 75539.0        ,  7.0       , 3)
PROCESS(`Wm,Wp', `Z,HH'         ,  500 GeV,     1.70944    ,  0.00008   , 3)
PROCESS(`Wm,Wp', `se1,SE1'      , 2000 GeV,    26.5380     ,  0.0006    , 4)
PROCESS(`Wm,Wp', `se2,SE2'      , 2000 GeV,     1.02966    ,  0.00004   , 3)
PROCESS(`Wm,Wp', `smu1,SMU1'    , 2000 GeV,    26.5380     ,  0.0006    , 4)
PROCESS(`Wm,Wp', `smu2,SMU2'    , 2000 GeV,     1.02972    ,  0.00004   , 3)
PROCESS(`Wm,Wp', `stau1,STAU1'  , 2000 GeV,     0.83284    ,  0.00003   , 3)
PROCESS(`Wm,Wp', `stau2,STAU2'  , 2000 GeV,    22.4195     ,  0.0005    , 3)
PROCESS(`Wm,Wp', `stau1,STAU2'  , 2000 GeV,     2.32935    ,  0.00005   , 3)
PROCESS(`Wm,Wp', `sn11,SN11'    , 2000 GeV,    23.486      ,  0.001     , 3)
PROCESS(`Wm,Wp', `sn21,SN21'    , 2000 GeV,    23.487      ,  0.001     , 3)
PROCESS(`Wm,Wp', `sn31,SN31'    , 2000 GeV,    23.429      ,  0.001     , 3)
PROCESS(`Wm,Wp', `su1,su1c'     , 2000 GeV,    41.590      ,  0.001     , 3)
PROCESS(`Wm,Wp', `su2,su2c'     , 2000 GeV,     1.07608    ,  0.00003   , 3)
PROCESS(`Wm,Wp', `sc1,sc1c'     , 2000 GeV,    41.588      ,  0.001     , 4)
PROCESS(`Wm,Wp', `sc2,sc2c'     , 2000 GeV,     1.07603    ,  0.00003   , 3)
PROCESS(`Wm,Wp', `st1,st1c'     , 2000 GeV,   180.637      ,  0.004     , 3)
PROCESS(`Wm,Wp', `st2,st2c'     , 2000 GeV,   204.461      ,  0.003     , 3)
PROCESS(`Wm,Wp', `st1,st2c'     , 2000 GeV,    85.178      ,  0.002     , 3)
PROCESS(`Wm,Wp', `sd1,sd1c'     , 2000 GeV,    39.0067     ,  0.0004    , 4)
PROCESS(`Wm,Wp', `sd2,sd2c'     , 2000 GeV,     0.269305   ,  0.000008  , 3)
PROCESS(`Wm,Wp', `ss1,ss1c'     , 2000 GeV,    39.0062     ,  0.0004    , 4)
PROCESS(`Wm,Wp', `ss2,ss2c'     , 2000 GeV,     0.269291   ,  0.000008  , 3)
PROCESS(`Wm,Wp', `sb1,sb1c'     , 2000 GeV,   141.457      ,  0.002     , 3)
PROCESS(`Wm,Wp', `sb2,sb2c'     , 2000 GeV,    19.7133     ,  0.0004    , 3)
PROCESS(`Wm,Wp', `sb1,sb2c'     , 2000 GeV,    61.090      ,  0.001     , 3)
PROCESS(`Wm,Wp', `neu1,neu1'    , 2000 GeV,     1.27423    ,  0.00008   , 3)
PROCESS(`Wm,Wp', `neu1,neu2'    , 2000 GeV,    24.472      ,  0.003     , 3)
PROCESS(`Wm,Wp', `neu1,neu3'    , 2000 GeV,    12.8790     ,  0.0008    , 4)
PROCESS(`Wm,Wp', `neu1,neu4'    , 2000 GeV,     9.7064     ,  0.0007    , 3)
PROCESS(`Wm,Wp', `neu2,neu2'    , 2000 GeV,  1041.50       ,  0.05      , 3)
PROCESS(`Wm,Wp', `neu2,neu3'    , 2000 GeV,   365.615      ,  0.006     , 4)
PROCESS(`Wm,Wp', `neu2,neu4'    , 2000 GeV,   467.775      ,  0.008     , 3)
PROCESS(`Wm,Wp', `neu3,neu3'    , 2000 GeV,    82.347      ,  0.008     , 3)
PROCESS(`Wm,Wp', `neu3,neu4'    , 2000 GeV,   138.18       ,  0.01      , 3)
PROCESS(`Wm,Wp', `neu4,neu4'    , 2000 GeV,   117.80       ,  0.01      , 3)
PROCESS(`Wm,Wp', `"ch1+","ch1-"', 2000 GeV,   944.2        ,  0.1       , 3)
PROCESS(`Wm,Wp', `"ch2+","ch2-"', 2000 GeV,   258.37       ,  0.04      , 3)
PROCESS(`Wm,Wp', `"ch1+","ch2-"', 2000 GeV,   130.98       ,  0.02      , 3)
PROCESS(`Wm,Wp', `h,h'          , 2000 GeV,  6061.0        ,  1.0       , 3)
PROCESS(`Wm,Wp', `h,HH'         , 2000 GeV,     2.1752     ,  0.0006    , 3)
PROCESS(`Wm,Wp', `HH,HH'        , 2000 GeV,     6.7509     ,  0.0011    , 3)
PROCESS(`Wm,Wp', `HA,HA'        , 2000 GeV,     6.7273     ,  0.0004    , 3)
PROCESS(`Wm,Wp', `Z,h'          , 2000 GeV, 86200.0        , 20.0       , 3)
PROCESS(`Wm,Wp', `Z,HH'         , 2000 GeV,    16.3939     ,  0.0037    , 3)
PROCESS(`Wm,Wp', `HA,h'         , 2000 GeV,     0.0060123  ,  0.0000007 , 3)
PROCESS(`Wm,Wp', `HA,HH'        , 2000 GeV,     3.4708     ,  0.0007    , 3)
PROCESS(`Wm,Wp', `Hp,Hm'        , 2000 GeV,    19.6060     ,  0.0023    , 3)
! --------------------------------------------------------------------------
END_TESTSUITE 
