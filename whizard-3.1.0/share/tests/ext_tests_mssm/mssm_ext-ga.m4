dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-ga.m4 > mssm_ext-ga.sin
dnl   whizard -r mssm_ext-ga.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-ga.m4 -1   $, mssm_ext_ga_)
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
PROCESS(`g,A', `su1,su1c'       , 2000 GeV,    55.4290     ,  0.0008    , 3)
PROCESS(`g,A', `su2,su2c'       , 2000 GeV,    57.0184     ,  0.0009    , 3)
PROCESS(`g,A', `sc1,sc1c'       , 2000 GeV,    55.4288     ,  0.0008    , 3)
PROCESS(`g,A', `sc2,sc2c'       , 2000 GeV,    57.0175     ,  0.0009    , 3)
PROCESS(`g,A', `st1,st1c'       , 2000 GeV,    73.382      ,  0.002     , 3)
PROCESS(`g,A', `st2,st2c'       , 2000 GeV,    54.1136     ,  0.0008    , 3)
PROCESS(`g,A', `sd1,sd1c'       , 2000 GeV,    13.7569     ,  0.0002    , 3)
PROCESS(`g,A', `sd2,sd2c'       , 2000 GeV,    14.2737     ,  0.0002    , 4)
PROCESS(`g,A', `ss1,ss1c'       , 2000 GeV,    13.7568     ,  0.0002    , 3)
PROCESS(`g,A', `ss2,ss2c'       , 2000 GeV,    14.2735     ,  0.0002    , 3)
PROCESS(`g,A', `sb1,sb1c'       , 2000 GeV,    14.9546     ,  0.0003    , 3)
PROCESS(`g,A', `sb2,sb2c'       , 2000 GeV,    14.3171     ,  0.0002    , 3)
! --------------------------------------------------------------------------
END_TESTSUITE 
