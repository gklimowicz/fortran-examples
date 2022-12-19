dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-bt.m4 > mssm_ext-bt.sin
dnl   whizard -r mssm_ext-bt.sin
dnl
dnl ---------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-bt.m4 -1   $, mssm_ext_bt_)
! -----------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                 -----
! -----------------------------------------------------------------------
model = MSSM
read_slha ("sps1a.slha")
?vis_history = false
! -----------------------------------------------------------------------
me = 0
mW = 80.419 
wW = 2.048
mZ = 91.1880
wZ = 2.446
mtop = 178.0
mb = 4.6
GF = 1.16639e-5
alphas = 0.118
?fatal_beam_decay=false
show (al_h,mu_h,tanb_h)
! -----------------------------------------------------------------------
iterations = 3:200000
! -----------------------------------------------------------------------
PROCESS(`b,T', `"ch1-",neu1'    ,  500 GeV,  37.75442   , 0.00007    , 3)
PROCESS(`b,T', `"ch1-",neu2'    ,  500 GeV, 171.6667    , 0.0006     , 4)
PROCESS(`b,T', `"ch2-",neu1'    ,  500 GeV,  17.91595   , 0.00004    , 3)
PROCESS(`b,T', `Z,Hm'           ,  500 GeV,  20.3544    , 0.0001     , 3)
PROCESS(`b,T', `sb1,st1c'       , 2000 GeV, 667.45      , 0.02       , 4)
PROCESS(`b,T', `sb1,st2c'       , 2000 GeV, 609.52      , 0.02       , 3)
PROCESS(`b,T', `sb2,st1c'       , 2000 GeV, 692.66      , 0.02       , 3)
PROCESS(`b,T', `sb2,st2c'       , 2000 GeV, 775.71      , 0.02       , 3)
PROCESS(`b,T', `"ch1-",neu1'    , 2000 GeV,   0.584741  , 0.000006   , 3)
PROCESS(`b,T', `"ch1-",neu2'    , 2000 GeV,   6.1435    , 0.0002     , 3)
PROCESS(`b,T', `"ch1-",neu3'    , 2000 GeV,   7.20626   , 0.00009    , 3)
PROCESS(`b,T', `"ch1-",neu4'    , 2000 GeV,   9.7429    , 0.0001     , 3)
PROCESS(`b,T', `"ch2-",neu1'    , 2000 GeV,   2.89723   , 0.00003    , 3)
PROCESS(`b,T', `"ch2-",neu2'    , 2000 GeV,   8.10775   , 0.00008    , 4)
PROCESS(`b,T', `"ch2-",neu3'    , 2000 GeV,  54.046     , 0.001      , 4)
PROCESS(`b,T', `"ch2-",neu4'    , 2000 GeV,  48.0844    , 0.0009     , 3)
PROCESS(`b,T', `Hm,h'           , 2000 GeV,  26.660     , 0.001      , 4)
PROCESS(`b,T', `Hm,HH'          , 2000 GeV,   2.00611   , 0.00008    , 3)
PROCESS(`b,T', `Hm,HA'          , 2000 GeV,   1.90817   , 0.00008    , 3)
PROCESS(`b,T', `Z,Hm'           , 2000 GeV,  34.766     , 0.001      , 3)
! -----------------------------------------------------------------------
END_TESTSUITE 
