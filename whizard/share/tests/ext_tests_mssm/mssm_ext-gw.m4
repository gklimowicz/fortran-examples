dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-gw.m4 > mssm_ext-gw.sin
dnl   whizard -r mssm_ext-gw.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-gw.m4 -1   $, mssm_ext_gw_)
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
PROCESS(`g,Wm', `sd1,su1c'      , 2000 GeV,   187.616      ,  0.003     , 3)
PROCESS(`g,Wm', `ss1,sc1c'      , 2000 GeV,   187.617      ,  0.003     , 3)
PROCESS(`g,Wm', `sb1,st1c'      , 2000 GeV,   138.625      ,  0.004     , 3)
PROCESS(`g,Wm', `sb2,st2c'      , 2000 GeV,    16.5095     ,  0.0003    , 3)
PROCESS(`g,Wm', `sb1,st2c'      , 2000 GeV,   195.692      ,  0.004     , 3)
PROCESS(`g,Wm', `sb2,st1c'      , 2000 GeV,    20.7532     ,  0.0007    , 3)
! --------------------------------------------------------------------------
END_TESTSUITE 
