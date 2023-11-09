dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-zz.m4 > mssm_ext-zz.sin
dnl   whizard -r mssm_ext-zz.sin
dnl
dnl ------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-zz.m4 -1   $, mssm_ext_zz_)
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
PROCESS(`Z,Z', `se1,SE1'        ,  500 GeV,    35.7923     ,  0.0004    , 3)
PROCESS(`Z,Z', `se2,SE2'        ,  500 GeV,    22.9508     ,  0.0004    , 3)
PROCESS(`Z,Z', `smu1,SMU1'      ,  500 GeV,    35.7920     ,  0.0004    , 3)
PROCESS(`Z,Z', `smu2,SMU2'      ,  500 GeV,    22.9509     ,  0.0004    , 3)
PROCESS(`Z,Z', `stau1,STAU1'    ,  500 GeV,    19.7282     ,  0.0003    , 3)
PROCESS(`Z,Z', `stau2,STAU2'    ,  500 GeV,    30.0574     ,  0.0002    , 3)
PROCESS(`Z,Z', `stau1,STAU2'    ,  500 GeV,     0.51455    ,  0.00002   , 3)
PROCESS(`Z,Z', `sn11,SN11'      ,  500 GeV,   232.517      ,  0.003     , 3)
PROCESS(`Z,Z', `sn21,SN21'      ,  500 GeV,   232.515      ,  0.003     , 3)
PROCESS(`Z,Z', `sn31,SN31'      ,  500 GeV,   233.341      ,  0.003     , 3)
PROCESS(`Z,Z', `h,h'            ,  500 GeV,  7887.5        ,  0.1       , 3)
PROCESS(`Z,Z', `se1,SE1'        , 2000 GeV,     3.8011     ,  0.0002    , 3)
dnl Bad WHIZARD 1 integration in the CATPISS paper
PROCESS(`Z,Z', `se2,SE2'        , 2000 GeV,     1.9234     ,  0.0001    , 5)
PROCESS(`Z,Z', `smu1,SMU1'      , 2000 GeV,     3.8008     ,  0.0002    , 3)
PROCESS(`Z,Z', `smu2,SMU2'      , 2000 GeV,     1.9239     ,  0.0001    , 3)
PROCESS(`Z,Z', `stau1,STAU1'    , 2000 GeV,     1.99985    ,  0.00008   , 3)
PROCESS(`Z,Z', `stau2,STAU2'    , 2000 GeV,     3.6161     ,  0.0001    , 3)
PROCESS(`Z,Z', `stau1,STAU2'    , 2000 GeV,     0.057456   ,  0.000002  , 3)
PROCESS(`Z,Z', `sn11,SN11'      , 2000 GeV,    32.037      ,  0.002     , 3)
PROCESS(`Z,Z', `sn21,SN21'      , 2000 GeV,    32.037      ,  0.002     , 3)
PROCESS(`Z,Z', `sn31,SN31'      , 2000 GeV,    32.072      ,  0.002     , 3)
PROCESS(`Z,Z', `su1,su1c'       , 2000 GeV,    15.6792     ,  0.0003    , 3)
PROCESS(`Z,Z', `su2,su2c'       , 2000 GeV,     1.20948    ,  0.00002   , 3)
PROCESS(`Z,Z', `sc1,sc1c'       , 2000 GeV,    15.6791     ,  0.0003    , 3)
PROCESS(`Z,Z', `sc2,sc2c'       , 2000 GeV,     1.20949    ,  0.00002   , 3)
PROCESS(`Z,Z', `st1,st1c'       , 2000 GeV,   262.155      ,  0.008     , 3)
PROCESS(`Z,Z', `st2,st2c'       , 2000 GeV,   475.11       ,  0.01      , 3)
PROCESS(`Z,Z', `st1,st2c'       , 2000 GeV,    10.7125     ,  0.0002    , 3)
PROCESS(`Z,Z', `sd1,sd1c'       , 2000 GeV,    30.5474     ,  0.0005    , 3)
PROCESS(`Z,Z', `sd2,sd2c'       , 2000 GeV,     0.238127   ,  0.000006  , 3)
PROCESS(`Z,Z', `ss1,ss1c'       , 2000 GeV,    30.5475     ,  0.0005    , 3)
PROCESS(`Z,Z', `ss2,ss2c'       , 2000 GeV,     0.238115   ,  0.000006  , 3)
PROCESS(`Z,Z', `sb1,sb1c'       , 2000 GeV,    20.7329     ,  0.0002    , 3)
PROCESS(`Z,Z', `sb2,sb2c'       , 2000 GeV,    10.6865     ,  0.0002    , 3)
PROCESS(`Z,Z', `sb1,sb2c'       , 2000 GeV,    18.6455     ,  0.0002    , 3)
PROCESS(`Z,Z', `h,h'            , 2000 GeV,  7802.5        ,  0.3       , 3)
PROCESS(`Z,Z', `h,HH'           , 2000 GeV,     2.7726     ,  0.0002    , 3)
PROCESS(`Z,Z', `HH,HH'          , 2000 GeV,    11.5209     ,  0.0004    , 3)
PROCESS(`Z,Z', `HA,HA'          , 2000 GeV,    11.3528     ,  0.0004    , 3)
PROCESS(`Z,Z', `Hp,Hm'          , 2000 GeV,     3.17136    ,  0.00005   , 3)
! --------------------------------------------------------------------------
END_TESTSUITE 
