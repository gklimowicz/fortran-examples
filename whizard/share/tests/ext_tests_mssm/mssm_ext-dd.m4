dnl Process this with
dnl
dnl   m4 testuite.m4 mssmtest-dd.m4 > mssmtest-dd.sin
dnl   whizard -r mssmtest-dd.sin
dnl
dnl ----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssmtest-dd.m4 1764 2010-02-11 23:33:52Z jr_reuter $, mssm_ext_dd_)
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
PROCESS(`d,D', `se1,SE1'        ,  500 GeV,   3.3472     , 0.0006     , 3)
PROCESS(`d,D', `se2,SE2'        ,  500 GeV,   2.0047     , 0.0003     , 3)
PROCESS(`d,D', `smu1,SMU1'      ,  500 GeV,   3.3465     , 0.0006     , 3)
PROCESS(`d,D', `smu2,SMU2'      ,  500 GeV,   2.0041     , 0.0003     , 3)
PROCESS(`d,D', `stau1,STAU1'    ,  500 GeV,   1.7271     , 0.0003     , 3)
PROCESS(`d,D', `stau2,STAU2'    ,  500 GeV,   2.4570     , 0.0004     , 3)
PROCESS(`d,D', `stau1,STAU2'    ,  500 GeV,   0.6950     , 0.0001     , 3)
PROCESS(`d,D', `sn11,SN11'      ,  500 GeV,   7.318      , 0.001      , 3)
dnl Bad WHIZARD 1 integration for the CATPISS paper
PROCESS(`d,D', `sn21,SN21'      ,  500 GeV,   7.314      , 0.001      , 4)
PROCESS(`d,D', `sn31,SN31'      ,  500 GeV,   7.454      , 0.001      , 3)
PROCESS(`d,D', `neu1,neu1'      ,  500 GeV,   0.1189331  , 0.0000007  , 3)
PROCESS(`d,D', `neu1,neu2'      ,  500 GeV,   0.249935   , 0.000001   , 3)
PROCESS(`d,D', `neu1,neu3'      ,  500 GeV,   0.817225   , 0.000004   , 3)
PROCESS(`d,D', `neu1,neu4'      ,  500 GeV,   0.0212673  , 0.0000002  , 3)
PROCESS(`d,D', `neu2,neu2'      ,  500 GeV,   1.939907   , 0.000009   , 3)
PROCESS(`d,D', `"ch1+","ch1-"'  ,  500 GeV, 137.161      , 0.003      , 3)
PROCESS(`d,D', `Z,h'            ,  500 GeV,  29.235      , 0.003      , 3)
PROCESS(`d,D', `Z,HH'           ,  500 GeV,   0.00030421 , 0.00000002 , 3)
PROCESS(`d,D', `sgl,sgl'        , 2000 GeV,1133.4        , 0.2        , 3)
PROCESS(`d,D', `se1,SE1'        , 2000 GeV,   0.9845     , 0.0002     , 3)
PROCESS(`d,D', `se2,SE2'        , 2000 GeV,   0.21578    , 0.00004    , 3)
PROCESS(`d,D', `smu1,SMU1'      , 2000 GeV,   0.9843     , 0.0002     , 3)
PROCESS(`d,D', `smu2,SMU2'      , 2000 GeV,   0.21572    , 0.00004    , 3)
PROCESS(`d,D', `stau1,STAU1'    , 2000 GeV,   0.17264    , 0.00003    , 3)
PROCESS(`d,D', `stau2,STAU2'    , 2000 GeV,   0.8171     , 0.0001     , 4)
PROCESS(`d,D', `stau1,STAU2'    , 2000 GeV,   0.10538    , 0.00002    , 3)
PROCESS(`d,D', `sn11,SN11'      , 2000 GeV,   1.4391     , 0.0002     , 3)
PROCESS(`d,D', `sn21,SN21'      , 2000 GeV,   1.4385     , 0.0003     , 3)
PROCESS(`d,D', `sn31,SN31'      , 2000 GeV,   1.4400     , 0.0002     , 3)
PROCESS(`d,D', `su1,su1c'       , 2000 GeV, 174.67       , 0.02       , 3)
PROCESS(`d,D', `su2,su2c'       , 2000 GeV, 185.19       , 0.03       , 3)
PROCESS(`d,D', `sc1,sc1c'       , 2000 GeV, 178.10       , 0.02       , 3)
PROCESS(`d,D', `sc2,sc2c'       , 2000 GeV, 185.21       , 0.02       , 3)
PROCESS(`d,D', `st1,st1c'       , 2000 GeV, 244.45       , 0.03       , 3)
PROCESS(`d,D', `st2,st2c'       , 2000 GeV, 168.80       , 0.02       , 3)
PROCESS(`d,D', `st1,st2c'       , 2000 GeV,   0.61183    , 0.00008    , 3)
PROCESS(`d,D', `sd1,sd1c'       , 2000 GeV, 790.3        , 0.1        , 3)
PROCESS(`d,D', `sd2,sd2c'       , 2000 GeV, 926.9        , 0.1        , 3)
PROCESS(`d,D', `ss1,ss1c'       , 2000 GeV, 175.92       , 0.02       , 3)
PROCESS(`d,D', `ss2,ss2c'       , 2000 GeV, 185.47       , 0.02       , 3)
PROCESS(`d,D', `sb1,sb1c'       , 2000 GeV, 200.54       , 0.03       , 3)
PROCESS(`d,D', `sb2,sb2c'       , 2000 GeV, 186.37       , 0.02       , 3)
PROCESS(`d,D', `sb1,sb2c'       , 2000 GeV,   0.25429    , 0.00005    , 3)
PROCESS(`d,D', `neu1,neu1'      , 2000 GeV,   0.079125   , 0.000004   , 3)
PROCESS(`d,D', `neu1,neu2'      , 2000 GeV,   0.34310    , 0.00002    , 3)
PROCESS(`d,D', `neu1,neu3'      , 2000 GeV,   0.173875   , 0.000003   , 3)
PROCESS(`d,D', `neu1,neu4'      , 2000 GeV,   0.140020   , 0.000003   , 3)
PROCESS(`d,D', `neu2,neu2'      , 2000 GeV,   3.1011     , 0.0002     , 3)
PROCESS(`d,D', `neu2,neu3'      , 2000 GeV,   1.07909    , 0.00002    , 3)
PROCESS(`d,D', `neu2,neu4'      , 2000 GeV,   1.16852    , 0.00006    , 3)
PROCESS(`d,D', `neu3,neu3'      , 2000 GeV,   0.00266298 , 0.00000004 , 3)
PROCESS(`d,D', `neu3,neu4'      , 2000 GeV,   4.76810    , 0.00009    , 3)
PROCESS(`d,D', `neu4,neu4'      , 2000 GeV,   0.087994   , 0.000006   , 3)
PROCESS(`d,D', `"ch1+","ch1-"'  , 2000 GeV,  10.504      , 0.002      , 3)
dnl Very bad fluctuation in WHIZARD 1 CATPISS value, corrected value
dnl for the cross section, 4.4954 -> 4.4960, cross-checked with WHIZARD 1
PROCESS(`d,D', `"ch2+","ch2-"'  , 2000 GeV,   4.4960     , 0.0001     , 3)
PROCESS(`d,D', `"ch1+","ch2-"'  , 2000 GeV,   0.77407    , 0.00005    , 4)
PROCESS(`d,D', `Z,h'            , 2000 GeV,   1.5337     , 0.0003     , 3)
PROCESS(`d,D', `Z,HH'           , 2000 GeV,   2.7112E-4  , 0.0005E-4  , 3)
! ------------------------------------------------------------------------
END_TESTSUITE
