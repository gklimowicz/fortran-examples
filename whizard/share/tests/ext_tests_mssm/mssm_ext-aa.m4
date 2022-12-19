dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-aa.m4 > mssmtest_ext-aa.sin
dnl   whizard -r mssm_ext-aa.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-aa.m4 -1   $, mssm_ext_aa_)
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
!
seed = 0
!
show (al_h,mu_h,tanb_h)
! --------------------------------------------------------------------------
iterations = 3:200000
! --------------------------------------------------------------------------
PROCESS(`A,A', `se1,SE1'        ,  500 GeV,   210.005      ,  0.007     , 3)
PROCESS(`A,A', `se2,SE2'        ,  500 GeV,   250.321      ,  0.011     , 3)
PROCESS(`A,A', `smu1,SMU1'      ,  500 GeV,   210.005      ,  0.007     , 3)
PROCESS(`A,A', `smu2,SMU2'      ,  500 GeV,   250.322      ,  0.011     , 3)
PROCESS(`A,A', `stau1,STAU1'    ,  500 GeV,   263.362      ,  0.013     , 3)
PROCESS(`A,A', `stau2,STAU2'    ,  500 GeV,   207.618      ,  0.007     , 3)
PROCESS(`A,A', `"ch1+","ch1-"'  ,  500 GeV,  1459.04       ,  0.06      , 3)
dnl Misprint in the CATPISS paper: 29.056 instead of 20.056 !!!
PROCESS(`A,A', `se1,SE1'        , 2000 GeV,    29.056      ,  0.005     , 3)
PROCESS(`A,A', `se2,SE2'        , 2000 GeV,    31.381      ,  0.006     , 3)
PROCESS(`A,A', `smu1,SMU1'      , 2000 GeV,    29.056      ,  0.005     , 3)
PROCESS(`A,A', `smu2,SMU2'      , 2000 GeV,    31.379      ,  0.006     , 3)
PROCESS(`A,A', `stau1,STAU1'    , 2000 GeV,    31.714      ,  0.006     , 3)
PROCESS(`A,A', `stau2,STAU2'    , 2000 GeV,    28.897      ,  0.005     , 3)
PROCESS(`A,A', `su1,su1c'       , 2000 GeV,     9.4536     ,  0.0004    , 3)
PROCESS(`A,A', `su2,su2c'       , 2000 GeV,     9.7244     ,  0.0005    , 3)
PROCESS(`A,A', `sc1,sc1c'       , 2000 GeV,     9.4534     ,  0.0004    , 3)
dnl Very bad WHIZARD 1 integration for the CATPISS paper
PROCESS(`A,A', `sc2,sc2c'       , 2000 GeV,     9.7230     ,  0.0005    , 4)
PROCESS(`A,A', `st1,st1c'       , 2000 GeV,    12.5159     ,  0.0009    , 3)
PROCESS(`A,A', `st2,st2c'       , 2000 GeV,     9.2298     ,  0.0004    , 3)
PROCESS(`A,A', `sd1,sd1c'       , 2000 GeV,     0.58655    ,  0.00003   , 3)
PROCESS(`A,A', `sd2,sd2c'       , 2000 GeV,     0.60853    ,  0.00003   , 3)
PROCESS(`A,A', `ss1,ss1c'       , 2000 GeV,     0.58656    ,  0.00003   , 3)
PROCESS(`A,A', `ss2,ss2c'       , 2000 GeV,     0.60863    ,  0.00003   , 3)
PROCESS(`A,A', `sb1,sb1c'       , 2000 GeV,     0.63761    ,  0.00003   , 3)
PROCESS(`A,A', `sb2,sb2c'       , 2000 GeV,     0.61045    ,  0.00003   , 3)
dnl Bad WHIZARD 1 integration for the CATPISS paper
PROCESS(`A,A', `"ch1+","ch1-"'  , 2000 GeV,   274.020      ,  0.009     , 4)
PROCESS(`A,A', `"ch2+","ch2-"'  , 2000 GeV,   181.542      ,  0.006     , 3)
dnl Bad WHIZARD 1 integration for the CATPISS paper
PROCESS(`A,A', `Hp,Hm'          , 2000 GeV,    20.644      ,  0.002     , 4)
! --------------------------------------------------------------------------
END_TESTSUITE 
