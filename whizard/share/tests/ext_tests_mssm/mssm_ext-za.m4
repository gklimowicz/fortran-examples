dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-za.m4 > mssm_ext-za.sin
dnl   whizard -r mssm_ext-za.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-za.m4 -1   $, mssm_ext_za_)
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
PROCESS(`Z,A', `se1,SE1'        ,  500 GeV,    64.0633     ,  0.0007    , 3)
PROCESS(`Z,A', `se2,SE2'        ,  500 GeV,    50.7284     ,  0.0007    , 3)
PROCESS(`Z,A', `smu1,SMU1'      ,  500 GeV,    64.0628     ,  0.0007    , 3)
PROCESS(`Z,A', `smu2,SMU2'      ,  500 GeV,    50.7284     ,  0.0007    , 3)
PROCESS(`Z,A', `stau1,STAU1'    ,  500 GeV,    36.4567     ,  0.0006    , 3)
PROCESS(`Z,A', `stau2,STAU2'    ,  500 GeV,    46.6053     ,  0.0005    , 3)
PROCESS(`Z,A', `stau1,STAU2'    ,  500 GeV,    24.0446     ,  0.0003    , 3)
PROCESS(`Z,A', `se1,SE1'        , 2000 GeV,     8.7333     ,  0.0005    , 3)
PROCESS(`Z,A', `se2,SE2'        , 2000 GeV,     6.0451     ,  0.0004    , 3)
PROCESS(`Z,A', `smu1,SMU1'      , 2000 GeV,     8.7327     ,  0.0005    , 3)
PROCESS(`Z,A', `smu2,SMU2'      , 2000 GeV,     6.0455     ,  0.0004    , 3)
PROCESS(`Z,A', `stau1,STAU1'    , 2000 GeV,     4.1336     ,  0.0003    , 3)
PROCESS(`Z,A', `stau2,STAU2'    , 2000 GeV,     6.3907     ,  0.0004    , 3)
PROCESS(`Z,A', `stau1,STAU2'    , 2000 GeV,     2.3102     ,  0.0001    , 3)
PROCESS(`Z,A', `su1,su1c'       , 2000 GeV,    10.1949     ,  0.0002    , 3)
PROCESS(`Z,A', `su2,su2c'       , 2000 GeV,     1.86042    ,  0.00003   , 3)
PROCESS(`Z,A', `sc1,sc1c'       , 2000 GeV,    10.1949     ,  0.0002    , 3)
PROCESS(`Z,A', `sc2,sc2c'       , 2000 GeV,     1.86039    ,  0.00003   , 3)
PROCESS(`Z,A', `st1,st1c'       , 2000 GeV,     0.00126510 ,  0.00000003, 3)
PROCESS(`Z,A', `st2,st2c'       , 2000 GeV,     3.44660    ,  0.00005   , 3)
PROCESS(`Z,A', `st1,st2c'       , 2000 GeV,    19.0982     ,  0.0007    , 3)
PROCESS(`Z,A', `sd1,sd1c'       , 2000 GeV,     3.70757    ,  0.00005   , 3)
PROCESS(`Z,A', `sd2,sd2c'       , 2000 GeV,     0.116438   ,  0.000002  , 3)
PROCESS(`Z,A', `ss1,ss1c'       , 2000 GeV,     3.70774    ,  0.00006   , 3)
PROCESS(`Z,A', `ss2,ss2c'       , 2000 GeV,     0.116435   ,  0.000002  , 3)
PROCESS(`Z,A', `sb1,sb1c'       , 2000 GeV,     3.12782    ,  0.00006   , 3)
PROCESS(`Z,A', `sb2,sb2c'       , 2000 GeV,     0.0114501  ,  0.0000002 , 3)
PROCESS(`Z,A', `sb1,sb2c'       , 2000 GeV,     0.533885   ,  0.000009  , 3)
PROCESS(`Z,A', `Hp,Hm'          , 2000 GeV,     6.1849     ,  0.00005   , 4)
! --------------------------------------------------------------------------
END_TESTSUITE 
