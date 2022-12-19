dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-ee2.m4 > mssm_ext-ee2.sin
dnl   whizard -r mssm_ext-ee2.sin
dnl
dnl ---------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-ee2.m4 -1   $, mssm_ext_ee2_)
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
show (al_h,mu_h,tanb_h)
! -----------------------------------------------------------------------
iterations = 3:200000
! -----------------------------------------------------------------------
PROCESS(`e1,e1', `se1,se1'      ,  500 GeV, 520.31      , 0.03       , 3)
PROCESS(`e1,e1', `se2,se2'      ,  500 GeV, 459.59      , 0.01       , 3)
dnl Bad WHIZARD 1 integration for the CATPISS paper
PROCESS(`e1,e1', `se1,se2'      ,  500 GeV, 159.96      , 0.02       , 5)
PROCESS(`e1,e1', `se1,se1'      , 2000 GeV,  36.836     , 0.002      , 3)
PROCESS(`e1,e1', `se2,se2'      , 2000 GeV,  28.650     , 0.003      , 3)
dnl Bad WHIZARD 1 integration for the CATPISS paper
PROCESS(`e1,e1', `se1,se2'      , 2000 GeV,  56.522     , 0.008      , 4)
! -----------------------------------------------------------------------
END_TESTSUITE 
