dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-uu.m4 > mssm_ext-uu.sin
dnl   whizard -r mssm_ext-uu.sin
dnl
dnl ----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-uu.m4 1764 2010-02-11 23:33:52Z jr_reuter $, mssm_ext_uu_)
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
PROCESS(`u,U', `se1,SE1'        ,  500 GeV,   5.1698     , 0.0009     , 3)
PROCESS(`u,U', `se2,SE2'        ,  500 GeV,   6.538      , 0.001      , 3)
PROCESS(`u,U', `smu1,SMU1'      ,  500 GeV,   5.1687     , 0.0009     , 3)
PROCESS(`u,U', `smu2,SMU2'      ,  500 GeV,   6.536      , 0.001      , 4)
PROCESS(`u,U', `stau1,STAU1'    ,  500 GeV,   6.992      , 0.001      , 3)
dnl Misprint in the CATPISS paper: 4.1264 instead of 4.1246 !!!
PROCESS(`u,U', `stau2,STAU2'    ,  500 GeV,   4.1264     , 0.0007     , 3)
PROCESS(`u,U', `stau1,STAU2'    ,  500 GeV,   0.54193    , 0.00009    , 3)
PROCESS(`u,U', `sn11,SN11'      ,  500 GeV,   5.706      , 0.001      , 3)
PROCESS(`u,U', `sn21,SN21'      ,  500 GeV,   5.704      , 0.001      , 5)
PROCESS(`u,U', `sn31,SN31'      ,  500 GeV,   5.813      , 0.001      , 3)
PROCESS(`u,U', `neu1,neu1'      ,  500 GeV,   2.24829    , 0.00002    , 3)
PROCESS(`u,U', `neu1,neu2'      ,  500 GeV,   0.0538560  , 0.0000009  , 3)
PROCESS(`u,U', `neu1,neu3'      ,  500 GeV,   0.524526   , 0.000003   , 3)
PROCESS(`u,U', `neu1,neu4'      ,  500 GeV,   0.00982339 , 0.00000008 , 3)
PROCESS(`u,U', `neu2,neu2'      ,  500 GeV,   3.66472    , 0.00003    , 3)
PROCESS(`u,U', `"ch1+","ch1-"'  ,  500 GeV, 153.97       , 0.02       , 3)
PROCESS(`u,U', `Z,h'            ,  500 GeV,  22.797      , 0.001      , 3)
PROCESS(`u,U', `Z,HH'           ,  500 GeV,   0.000237224, 0.000000001, 4)
PROCESS(`u,U', `sgl,sgl'        , 2000 GeV,1137.8        , 0.2        , 3)
PROCESS(`u,U', `se1,SE1'        , 2000 GeV,   1.5469     , 0.0002     , 3)
PROCESS(`u,U', `se2,SE2'        , 2000 GeV,   0.7318     , 0.0001     , 3)
PROCESS(`u,U', `smu1,SMU1'      , 2000 GeV,   1.5466     , 0.0003     , 3)
PROCESS(`u,U', `smu2,SMU2'      , 2000 GeV,   0.7316     , 0.0001     , 4)
PROCESS(`u,U', `stau1,STAU1'    , 2000 GeV,   0.7194     , 0.0001     , 3)
PROCESS(`u,U', `stau2,STAU2'    , 2000 GeV,   1.3957     , 0.0002     , 6)
PROCESS(`u,U', `stau1,STAU2'    , 2000 GeV,   0.08217    , 0.00001    , 3)
PROCESS(`u,U', `sn11,SN11'      , 2000 GeV,   1.1222     , 0.0002     , 3)
PROCESS(`u,U', `sn21,SN21'      , 2000 GeV,   1.1217     , 0.0002     , 6)
PROCESS(`u,U', `sn31,SN31'      , 2000 GeV,   1.1229     , 0.0002     , 3)
PROCESS(`u,U', `su1,su1c'       , 2000 GeV, 799.6        , 0.1        , 3)
PROCESS(`u,U', `su2,su2c'       , 2000 GeV, 879.7        , 0.1        , 3)
PROCESS(`u,U', `su1,su2c'       , 2000 GeV, 784.16       , 0.03       , 3)
PROCESS(`u,U', `sc1,sc1c'       , 2000 GeV, 178.39       , 0.02       , 3)
PROCESS(`u,U', `sc2,sc2c'       , 2000 GeV, 185.62       , 0.02       , 3)
PROCESS(`u,U', `st1,st1c'       , 2000 GeV, 245.11       , 0.03       , 3)
PROCESS(`u,U', `st2,st2c'       , 2000 GeV, 169.22       , 0.02       , 3)
PROCESS(`u,U', `st1,st2c'       , 2000 GeV,   0.47714    , 0.00008    , 3)
PROCESS(`u,U', `sd1,sd1c'       , 2000 GeV, 166.60       , 0.02       , 3)
PROCESS(`u,U', `sd2,sd2c'       , 2000 GeV, 185.56       , 0.03       , 3)
PROCESS(`u,U', `ss1,ss1c'       , 2000 GeV, 175.68       , 0.02       , 3)
PROCESS(`u,U', `ss2,ss2c'       , 2000 GeV, 185.58       , 0.02       , 3)
PROCESS(`u,U', `sb1,sb1c'       , 2000 GeV, 200.364      , 0.008      , 3)
PROCESS(`u,U', `sb2,sb2c'       , 2000 GeV, 186.500      , 0.007      , 3)
PROCESS(`u,U', `sb1,sb2c'       , 2000 GeV,   0.198272   , 0.000008   , 3)
PROCESS(`u,U', `neu1,neu1'      , 2000 GeV,   1.2165     , 0.0001     , 3)
PROCESS(`u,U', `neu1,neu2'      , 2000 GeV,   0.10850    , 0.00001    , 3)
PROCESS(`u,U', `neu1,neu3'      , 2000 GeV,   0.096752   , 0.000005   , 4)
PROCESS(`u,U', `neu1,neu4'      , 2000 GeV,   0.067293   , 0.000006   , 4)
PROCESS(`u,U', `neu2,neu2'      , 2000 GeV,   4.2296     , 0.0004     , 3)
PROCESS(`u,U', `neu2,neu3'      , 2000 GeV,   0.211458   , 0.000008   , 6)
PROCESS(`u,U', `neu2,neu4'      , 2000 GeV,   0.55025    , 0.00008    , 3)
PROCESS(`u,U', `neu3,neu3'      , 2000 GeV,   0.00033843 , 0.00000001 , 3)
PROCESS(`u,U', `neu3,neu4'      , 2000 GeV,   4.4433     , 0.0002     , 4)
PROCESS(`u,U', `neu4,neu4'      , 2000 GeV,   0.016389   , 0.000003   , 3)
PROCESS(`u,U', `"ch1+","ch1-"'  , 2000 GeV,  10.734      , 0.002      , 3)
PROCESS(`u,U', `"ch2+","ch2-"'  , 2000 GeV,   5.0401     , 0.0002     , 3)
PROCESS(`u,U', `"ch1+","ch2-"'  , 2000 GeV,   1.5362     , 0.0002     , 3)
PROCESS(`u,U', `Z,h'            , 2000 GeV,   1.1960     , 0.0002     , 3)
PROCESS(`u,U', `Z,HH'           , 2000 GeV,   2.1142E-4  , 0.0004E-4  , 3)
! ------------------------------------------------------------------------
END_TESTSUITE
