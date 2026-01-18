dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-uuckm.m4 > mssm_ext-uuckm.sin
dnl   whizard -r mssm_ext-uuckm.sin
dnl
dnl ---------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-uuckm.m4 -1   $, mssm_ext_uuckm_)
! -----------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                 -----
! -----------------------------------------------------------------------
! ----- Note that these values deviate from those in the CATPISS paper 
! ----- as we are using here a different set of inputs for the CKM matrix
! -----------------------------------------------------------------------
model = MSSM_CKM
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
show (al_h,mu_h,tanb_h)
! -----------------------------------------------------------------------
iterations = 3:200000
! -----------------------------------------------------------------------
PROCESS(`u,U', `sd1,sd1c'     , 2000 GeV,  160.231     , 0.012      , 3)
PROCESS(`u,U', `ss1,ss1c'     , 2000 GeV,  168.352     , 0.012      , 3)
! -----------------------------------------------------------------------
END_TESTSUITE 
