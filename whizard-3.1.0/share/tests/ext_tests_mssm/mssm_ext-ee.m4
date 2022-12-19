dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-ee.m4 > mssm_ext-ee.sin
dnl   whizard -r mssm_ext-ee.sin
dnl
dnl ----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-ee.m4 1775 2010-02-12 17:55:32Z jr_reuter $, mssm_ext_ee_)
! ------------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                  -----
! ------------------------------------------------------------------------
model = MSSM
read_slha ("sps1a.slha")
?vis_history = false
! ------------------------------------------------------------------------
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
! ------------------------------------------------------------------------
iterations = 3:200000
! ------------------------------------------------------------------------
PROCESS(`e1,E1', `se1,SE1'      ,  500 GeV,  54.687      , 0.003      , 3)
PROCESS(`e1,E1', `se2,SE2'      ,  500 GeV, 274.682      , 0.001      , 4)
PROCESS(`e1,E1', `se1,SE2'      ,  500 GeV,  75.167      , 0.003      , 3)
PROCESS(`e1,E1', `smu1,SMU1'    ,  500 GeV,  22.5478     , 0.0009     , 3)
PROCESS(`e1,E1', `smu2,SMU2'    ,  500 GeV,  51.837      , 0.002      , 3)
PROCESS(`e1,E1', `stau1,STAU1'  ,  500 GeV,  55.580      , 0.002      , 3)
PROCESS(`e1,E1', `stau2,STAU2'  ,  500 GeV,  19.0174     , 0.0007     , 3)
PROCESS(`e1,E1', `stau1,STAU2'  ,  500 GeV,   1.41191    , 0.00005    , 3)
PROCESS(`e1,E1', `sn11,SN11'    ,  500 GeV, 493.38       , 0.02       , 3)
PROCESS(`e1,E1', `sn21,SN21'    ,  500 GeV,  14.8638     , 0.0006     , 3)
PROCESS(`e1,E1', `sn31,SN31'    ,  500 GeV,  15.1394     , 0.0008     , 3)
PROCESS(`e1,E1', `neu1,neu1'    ,  500 GeV, 240.636      , 0.007      , 3)
PROCESS(`e1,E1', `neu1,neu2'    ,  500 GeV,  62.374      , 0.002      , 3)
PROCESS(`e1,E1', `neu1,neu3'    ,  500 GeV,   7.78131    , 0.00004    , 3)
PROCESS(`e1,E1', `neu1,neu4'    ,  500 GeV,   1.03460    , 0.00003    , 3)
PROCESS(`e1,E1', `neu2,neu2'    ,  500 GeV,  70.730      , 0.003      , 3)
PROCESS(`e1,E1', `"ch1+","ch1-"',  500 GeV, 162.788      , 0.007      , 3)
PROCESS(`e1,E1', `Z,h'          ,  500 GeV,  59.376      , 0.001      , 3)
PROCESS(`e1,E1', `Z,HH'         ,  500 GeV,   6.179180E-4, 0.000005E-4, 3)
PROCESS(`e1,E1', `se1,SE1'      , 2000 GeV,  78.866      , 0.004      , 3)
PROCESS(`e1,E1', `se2,SE2'      , 2000 GeV,  91.776      , 0.005      , 3)
PROCESS(`e1,E1', `se1,SE2'      , 2000 GeV,   7.2372     , 0.0004     , 3)
PROCESS(`e1,E1', `smu1,SMU1'    , 2000 GeV,   6.8265     , 0.0003     , 3)
PROCESS(`e1,E1', `smu2,SMU2'    , 2000 GeV,   5.8105     , 0.0002     , 3)
PROCESS(`e1,E1', `stau1,STAU1'  , 2000 GeV,   5.7141     , 0.0002     , 3)
PROCESS(`e1,E1', `stau2,STAU2'  , 2000 GeV,   6.5045     , 0.0003     , 3)
PROCESS(`e1,E1', `stau1,STAU2'  , 2000 GeV,   0.214058   , 0.000008   , 3)
PROCESS(`e1,E1', `sn11,SN11'    , 2000 GeV, 272.15       , 0.01       , 3)
PROCESS(`e1,E1', `sn21,SN21'    , 2000 GeV,   2.9232     , 0.0001     , 3)
PROCESS(`e1,E1', `sn31,SN31'    , 2000 GeV,   2.9245     , 0.0001     , 3)
PROCESS(`e1,E1', `su1,su1c'     , 2000 GeV,   7.6188     , 0.0003     , 3)
PROCESS(`e1,E1', `su2,su2c'     , 2000 GeV,   4.6935     , 0.0002     , 3)
PROCESS(`e1,E1', `sc1,sc1c'     , 2000 GeV,   7.6182     , 0.0003     , 3)
PROCESS(`e1,E1', `sc2,sc2c'     , 2000 GeV,   4.6933     , 0.0002     , 3)
PROCESS(`e1,E1', `st1,st1c'     , 2000 GeV,   5.9847     , 0.0002     , 3)
PROCESS(`e1,E1', `st2,st2c'     , 2000 GeV,   5.3792     , 0.0002     , 3)
PROCESS(`e1,E1', `st1,st2c'     , 2000 GeV,   1.24264    , 0.00005    , 3)
PROCESS(`e1,E1', `sd1,sd1c'     , 2000 GeV,   5.2059     , 0.0002     , 3)
PROCESS(`e1,E1', `sd2,sd2c'     , 2000 GeV,   1.17595    , 0.00005    , 3)
PROCESS(`e1,E1', `ss1,ss1c'     , 2000 GeV,   5.2058     , 0.0002     , 3)
PROCESS(`e1,E1', `ss2,ss2c'     , 2000 GeV,   1.17585    , 0.00005    , 3)
PROCESS(`e1,E1', `sb1,sb1c'     , 2000 GeV,   4.9387     , 0.0002     , 3)
PROCESS(`e1,E1', `sb2,sb2c'     , 2000 GeV,   1.12946    , 0.00004    , 3)
PROCESS(`e1,E1', `sb1,sb2c'     , 2000 GeV,   0.516432   , 0.000009   , 3)
PROCESS(`e1,E1', `neu1,neu1'    , 2000 GeV,  26.3087     , 0.0009     , 3)
PROCESS(`e1,E1', `neu1,neu2'    , 2000 GeV,   9.9475     , 0.0004     , 3)
PROCESS(`e1,E1', `neu1,neu3'    , 2000 GeV,   0.64796    , 0.00001    , 3)
PROCESS(`e1,E1', `neu1,neu4'    , 2000 GeV,   1.36564    , 0.00005    , 3)
PROCESS(`e1,E1', `neu2,neu2'    , 2000 GeV,  18.6845     , 0.0008     , 3)
PROCESS(`e1,E1', `neu2,neu3'    , 2000 GeV,   1.85590    , 0.00004    , 3)
PROCESS(`e1,E1', `neu2,neu4'    , 2000 GeV,   3.03951    , 0.00009    , 3)
PROCESS(`e1,E1', `neu3,neu3'    , 2000 GeV,   0.0042214  , 0.0000002  , 3)
PROCESS(`e1,E1', `neu3,neu4'    , 2000 GeV,   9.9362     , 0.0003     , 3)
PROCESS(`e1,E1', `neu4,neu4'    , 2000 GeV,   0.135482   , 0.000005   , 3)
PROCESS(`e1,E1', `"ch1+","ch1-"', 2000 GeV,  45.080      , 0.002      , 3)
PROCESS(`e1,E1', `"ch2+","ch2-"', 2000 GeV,  26.9864     , 0.0006     , 3)
PROCESS(`e1,E1', `"ch1+","ch2-"', 2000 GeV,   4.01053    , 0.00009    , 3)
PROCESS(`e1,E1', `Z,h'          , 2000 GeV,   3.11492    , 0.00009    , 3)
PROCESS(`e1,E1', `Z,HH'         , 2000 GeV,   5.5058E-4  , 0.0002E-4  , 3)
PROCESS(`e1,E1', `HA,h'         , 2000 GeV,   5.3433E-4  , 0.0002E-4  , 3)
PROCESS(`e1,E1', `HA,HH'        , 2000 GeV,   2.37434    , 0.00009    , 3)
PROCESS(`e1,E1', `Hp,Hm'        , 2000 GeV,   5.5339     , 0.0002     , 3)
! ------------------------------------------------------------------------
END_TESTSUITE
