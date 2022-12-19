dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-gg.m4 > mssm_ext-gg.sin
dnl   whizard -r mssm_ext-gg.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-gg.m4 -1   $, mssm_ext_gg_)
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
dnl There are sometimes numerical instabilities in this integration
PROCESS(`g,g', `sgl,sgl'        , 2000 GeV, 13575.6        ,  0.1       , 5)
PROCESS(`g,g', `su1,su1c'       , 2000 GeV,   185.615      ,  0.003     , 4)
PROCESS(`g,g', `su2,su2c'       , 2000 GeV,   191.590      ,  0.003     , 4)
PROCESS(`g,g', `sc1,sc1c'       , 2000 GeV,   185.612      ,  0.003     , 3)
PROCESS(`g,g', `sc2,sc2c'       , 2000 GeV,   191.588      ,  0.003     , 3)
PROCESS(`g,g', `st1,st1c'       , 2000 GeV,   250.71       ,  0.01      , 4)
PROCESS(`g,g', `st2,st2c'       , 2000 GeV,   180.541      ,  0.003     , 3)
PROCESS(`g,g', `sd1,sd1c'       , 2000 GeV,   184.081      ,  0.003     , 4)
PROCESS(`g,g', `sd2,sd2c'       , 2000 GeV,   191.875      ,  0.003     , 4)
PROCESS(`g,g', `ss1,ss1c'       , 2000 GeV,   184.079      ,  0.003     , 4)
PROCESS(`g,g', `ss2,ss2c'       , 2000 GeV,   191.873      ,  0.003     , 3)
PROCESS(`g,g', `sb1,sb1c'       , 2000 GeV,   201.884      ,  0.004     , 3)
dnl Bad WHIZARD 1 integration in the CATPISS paper
PROCESS(`g,g', `sb2,sb2c'       , 2000 GeV,   192.516      ,  0.003     , 5)
! --------------------------------------------------------------------------
END_TESTSUITE 
