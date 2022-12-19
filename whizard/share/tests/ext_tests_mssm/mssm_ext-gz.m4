dnl Process this with
dnl
dnl   m4 testuite.m4 mssmtest-gz.m4 > mssmtest-gz.sin
dnl   whizard -r mssmtest-gz.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssmtest-gz.m4 -1   $, mssm_ext_gz_)
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
PROCESS(`g,Z', `su1,su1c'       , 2000 GeV,    59.7774     ,  0.0009    , 3)
PROCESS(`g,Z', `su2,su2c'       , 2000 GeV,    10.9085     ,  0.0002    , 3)
PROCESS(`g,Z', `sc1,sc1c'       , 2000 GeV,    59.7772     ,  0.0009    , 3)
PROCESS(`g,Z', `sc2,sc2c'       , 2000 GeV,    10.9084     ,  0.0002    , 3)
PROCESS(`g,Z', `st1,st1c'       , 2000 GeV,     0.0074179  ,  0.0000002 , 3)
PROCESS(`g,Z', `st2,st2c'       , 2000 GeV,    20.2091     ,  0.0003    , 3)
PROCESS(`g,Z', `st1,st2c'       , 2000 GeV,   111.986      ,  0.004     , 3)
dnl Bad WHIZARD 1 integration for the CATPISS paper
PROCESS(`g,Z', `sd1,sd1c'       , 2000 GeV,    86.9615     ,  0.0001    , 5)
PROCESS(`g,Z', `sd2,sd2c'       , 2000 GeV,     2.73090    ,  0.00004   , 4)
PROCESS(`g,Z', `ss1,ss1c'       , 2000 GeV,    86.959      ,  0.0001    , 4)
PROCESS(`g,Z', `ss2,ss2c'       , 2000 GeV,     2.73090    ,  0.00004   , 4)
dnl
PROCESS(`g,Z', `sb1,sb1c'       , 2000 GeV,    73.360      ,  0.001     , 3)
PROCESS(`g,Z', `sb2,sb2c'       , 2000 GeV,     0.268554   ,  0.000004  , 3)
PROCESS(`g,Z', `sb1,sb2c'       , 2000 GeV,    12.5214     ,  0.0002    , 3)
! --------------------------------------------------------------------------
END_TESTSUITE 
