dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-wz.m4 > mssm_ext-wz.sin
dnl   whizard -r mssm_ext-wz.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-wz.m4 -1   $, mssm_ext_wz_)
! --------------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                    -----
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
PROCESS(`Wm,Z', `se1,SN11'      ,  500 GeV,    96.639      ,  0.002     , 3)
PROCESS(`Wm,Z', `smu1,SN21'     ,  500 GeV,    96.638      ,  0.002     , 3)
PROCESS(`Wm,Z', `stau1,SN31'    ,  500 GeV,    14.952      ,  0.001     , 3)
PROCESS(`Wm,Z', `stau2,SN31'    ,  500 GeV,    85.875      ,  0.002     , 3)
dnl Bad WHIZARD 1 integration for the CATPISS paper
PROCESS(`Wm,Z', `"ch1-",neu1'   ,  500 GeV,    61.626      ,  0.003     , 4)
PROCESS(`Wm,Z', `"ch1-",neu2'   ,  500 GeV,  2835.0        ,  0.3       , 3)
PROCESS(`Wm,Z', `"ch2-",neu1'   ,  500 GeV,    11.7619     ,  0.0007    , 3)
PROCESS(`Wm,Z', `Wm,h'          ,  500 GeV,76213.0         ,  6.0       , 3)
PROCESS(`Wm,Z', `Wm,HH'         ,  500 GeV,    4.2446      ,  0.0002    , 3)
PROCESS(`Wm,Z', `Wm,HA'         ,  500 GeV,    1.07037     ,  0.00001   , 3)
PROCESS(`Wm,Z', `Z,Hm'          ,  500 GeV,    0.17723     ,  0.00002   , 3)
PROCESS(`Wm,Z', `se1,SN11'      , 2000 GeV,    15.728      ,  0.005     , 3)
PROCESS(`Wm,Z', `smu1,SN21'     , 2000 GeV,    15.727      ,  0.005     , 3)
PROCESS(`Wm,Z', `stau1,SN31'    , 2000 GeV,     1.4268     ,  0.0002    , 3)
PROCESS(`Wm,Z', `stau2,SN31'    , 2000 GeV,    14.478      ,  0.002     , 3)
PROCESS(`Wm,Z', `sd1,su1c'      , 2000 GeV,    24.220      ,  0.001     , 3)
PROCESS(`Wm,Z', `ss1,sc1c'      , 2000 GeV,    24.221      ,  0.001     , 3)
PROCESS(`Wm,Z', `sb1,st1c'      , 2000 GeV,    40.676      ,  0.004     , 3)
PROCESS(`Wm,Z', `sb2,st2c'      , 2000 GeV,     8.3706     ,  0.0007    , 3)
PROCESS(`Wm,Z', `sb1,st2c'      , 2000 GeV,    63.592      ,  0.006     , 3)
PROCESS(`Wm,Z', `sb2,st1c'      , 2000 GeV,     3.9236     ,  0.005     , 3)
PROCESS(`Wm,Z', `"ch1-",neu1'   , 2000 GeV,    16.389      ,  0.001     , 3)
PROCESS(`Wm,Z', `"ch1-",neu2'   , 2000 GeV,   668.1        ,  0.1       , 3)
PROCESS(`Wm,Z', `"ch1-",neu3'   , 2000 GeV,   278.53       ,  0.01      , 3)
PROCESS(`Wm,Z', `"ch1-",neu4'   , 2000 GeV,   270.97       ,  0.02      , 3)
PROCESS(`Wm,Z', `"ch2-",neu1'   , 2000 GeV,    12.380      ,  0.001     , 3)
PROCESS(`Wm,Z', `"ch2-",neu2'   , 2000 GeV,   218.38       ,  0.02      , 3)
PROCESS(`Wm,Z', `"ch2-",neu3'   , 2000 GeV,    76.494      ,  0.005     , 3)
PROCESS(`Wm,Z', `"ch2-",neu4'   , 2000 GeV,    97.693      ,  0.007     , 3)
PROCESS(`Wm,Z', `h,Hm'          , 2000 GeV,    0.0044399   ,  0.0000005 , 3)
PROCESS(`Wm,Z', `HH,Hm'         , 2000 GeV,    6.1592      ,  0.0002    , 3)
PROCESS(`Wm,Z', `HA,Hm'         , 2000 GeV,    5.9726      ,  0.0005    , 3)
PROCESS(`Wm,Z', `Wm,h'          , 2000 GeV,82890.0         , 20.0       , 3)
PROCESS(`Wm,Z', `Wm,HH'         , 2000 GeV,   15.783       ,  0.003     , 3)
PROCESS(`Wm,Z', `Wm,HA'         , 2000 GeV,    0.24815     ,  0.00007   , 3)
PROCESS(`Wm,Z', `Z,Hm'          , 2000 GeV,    0.25403     ,  0.00007   , 3)
! --------------------------------------------------------------------------
END_TESTSUITE 
