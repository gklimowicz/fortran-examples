dnl Process this with
dnl
dnl   m4 testuite.m4 mssmtest-wa.m4 > mssmtest-wa.sin
dnl   whizard -r mssmtest-wa.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssmtest-wa.m4 -1   $, mssm_ext_wa_)
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
PROCESS(`Wm,A', `se1,SN11'      ,  500 GeV,    92.927      ,  0.007     , 3)
PROCESS(`Wm,A', `smu1,SN21'     ,  500 GeV,    92.942      ,  0.007     , 3)
PROCESS(`Wm,A', `stau1,SN31'    ,  500 GeV,    12.100      ,  0.001     , 3)
PROCESS(`Wm,A', `stau2,SN31'    ,  500 GeV,    85.167      ,  0.007     , 3)
PROCESS(`Wm,A', `"ch1-",neu1'   ,  500 GeV,    15.821      ,  0.002     , 3)
PROCESS(`Wm,A', `"ch1-",neu2'   ,  500 GeV,  1223.5        ,  0.1       , 3)
PROCESS(`Wm,A', `"ch2-",neu1'   ,  500 GeV,     3.03742    ,  0.00007   , 3)
PROCESS(`Wm,A', `Wm,h'          ,  500 GeV,12855.0         ,  3.0       , 3)
PROCESS(`Wm,A', `Wm,HH'         ,  500 GeV,    0.54011     ,  0.00006   , 3)
PROCESS(`Wm,A', `se1,SN11'      , 2000 GeV,    14.477      ,  0.003     , 3)
PROCESS(`Wm,A', `smu1,SN21'     , 2000 GeV,    14.479      ,  0.003     , 3)
PROCESS(`Wm,A', `stau1,SN31'    , 2000 GeV,     1.2566     ,  0.0003    , 3)
PROCESS(`Wm,A', `stau2,SN31'    , 2000 GeV,    13.372      ,  0.003     , 3)
PROCESS(`Wm,A', `sd1,su1c'      , 2000 GeV,     6.260      ,  0.001     , 3)
PROCESS(`Wm,A', `ss1,sc1c'      , 2000 GeV,     6.262      ,  0.001     , 3)
PROCESS(`Wm,A', `sb1,st1c'      , 2000 GeV,     5.528      ,  0.001     , 3)
PROCESS(`Wm,A', `sb2,st2c'      , 2000 GeV,     0.5417     ,  0.0001    , 3)
PROCESS(`Wm,A', `sb1,st2c'      , 2000 GeV,     6.267      ,  0.001     , 3)
PROCESS(`Wm,A', `sb2,st1c'      , 2000 GeV,     0.8595     ,  0.002     , 3)
PROCESS(`Wm,A', `"ch1-",neu1'   , 2000 GeV,     3.8332     ,  0.0006    , 3)
PROCESS(`Wm,A', `"ch1-",neu2'   , 2000 GeV,   303.04       ,  0.05      , 3)
PROCESS(`Wm,A', `"ch1-",neu3'   , 2000 GeV,    50.902      ,  0.008     , 3)
PROCESS(`Wm,A', `"ch1-",neu4'   , 2000 GeV,    52.648      ,  0.008     , 3)
PROCESS(`Wm,A', `"ch2-",neu1'   , 2000 GeV,     6.5764     ,  0.0009    , 3)
PROCESS(`Wm,A', `"ch2-",neu2'   , 2000 GeV,    34.003      ,  0.005     , 3)
PROCESS(`Wm,A', `"ch2-",neu3'   , 2000 GeV,    47.719      ,  0.007     , 3)
PROCESS(`Wm,A', `"ch2-",neu4'   , 2000 GeV,    59.636      ,  0.008     , 3)
PROCESS(`Wm,A', `h,Hm'          , 2000 GeV,    0.0045192   ,  0.0000008 , 3)
PROCESS(`Wm,A', `HH,Hm'         , 2000 GeV,    4.9610      ,  0.0009    , 3)
PROCESS(`Wm,A', `HA,Hm'         , 2000 GeV,    4.9671      ,  0.0009    , 3)
PROCESS(`Wm,A', `Wm,h'          , 2000 GeV,15811.0         ,  4.0       , 3)
PROCESS(`Wm,A', `Wm,HH'         , 2000 GeV,    3.0172      ,  0.0007    , 3)
! --------------------------------------------------------------------------
END_TESTSUITE 
