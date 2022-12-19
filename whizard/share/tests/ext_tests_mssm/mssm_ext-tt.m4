dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-tt.m4 > mssm_ext-tt.sin
dnl   whizard -r mssm_ext-tt.sin
dnl
dnl ----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-tt.m4 -1   $, mssm_ext_tt_)
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
PROCESS(`e3,E3', `stau1,STAU1'  ,  500 GeV,  257.32      , 0.01       , 3)
PROCESS(`e3,E3', `stau2,STAU2'  ,  500 GeV,   46.368     , 0.002      , 3)
PROCESS(`e3,E3', `stau1,STAU2'  ,  500 GeV,   81.720     , 0.003      , 3)
PROCESS(`e3,E3', `sn31,SN31'    ,  500 GeV,  502.27      , 0.02       , 3)
PROCESS(`e3,E3', `neu1,neu1'    ,  500 GeV,  249.954     , 0.009      , 3)
PROCESS(`e3,E3', `neu1,neu2'    ,  500 GeV,   69.969     , 0.002      , 3)
PROCESS(`e3,E3', `neu1,neu3'    ,  500 GeV,   17.0394    , 0.0001     , 3)
PROCESS(`e3,E3', `neu1,neu4'    ,  500 GeV,    7.01414   , 0.00006    , 3)
PROCESS(`e3,E3', `neu2,neu2'    ,  500 GeV,   82.353     , 0.003      , 3)
PROCESS(`e3,E3', `"ch1+","ch1-"',  500 GeV,  185.093     , 0.006      , 3)
PROCESS(`e3,E3', `h,h'          ,  500 GeV,    0.35339   , 0.00002    , 3)
PROCESS(`e3,E3', `Z,h'          ,  500 GeV,   59.589     , 0.003      , 3)
PROCESS(`e3,E3', `Z,HH'         ,  500 GeV,    2.83169   , 0.00009    , 4)
PROCESS(`e3,E3', `Z,HA'         ,  500 GeV,    2.99162   , 0.00009    , 3)
PROCESS(`e3,E3', `stau1,STAU1'  , 2000 GeV,   79.636     , 0.004      , 3)
PROCESS(`e3,E3', `stau2,STAU2'  , 2000 GeV,   66.862     , 0.003      , 3)
PROCESS(`e3,E3', `stau1,STAU2'  , 2000 GeV,   18.9588    , 0.0008     , 3)
PROCESS(`e3,E3', `sn31,SN31'    , 2000 GeV,  272.01      , 0.01       , 3)
PROCESS(`e3,E3', `neu1,neu1'    , 2000 GeV,   26.431     , 0.001      , 3)
PROCESS(`e3,E3', `neu1,neu2'    , 2000 GeV,    9.8940    , 0.0004     , 3)
PROCESS(`e3,E3', `neu1,neu3'    , 2000 GeV,    0.79136   , 0.00002    , 3)
PROCESS(`e3,E3', `neu1,neu4'    , 2000 GeV,    1.5075    , 0.0005     , 3)
PROCESS(`e3,E3', `neu2,neu2'    , 2000 GeV,   18.8879    , 0.0009     , 3)
PROCESS(`e3,E3', `neu2,neu3'    , 2000 GeV,    1.75884   , 0.00005    , 3)
PROCESS(`e3,E3', `neu2,neu4'    , 2000 GeV,    2.9640    , 0.0001     , 3)
PROCESS(`e3,E3', `neu3,neu3'    , 2000 GeV,    0.0469966 , 0.0000009  , 3)
PROCESS(`e3,E3', `neu3,neu4'    , 2000 GeV,    8.5857    , 0.0003     , 3)
PROCESS(`e3,E3', `neu4,neu4'    , 2000 GeV,    0.264389  , 0.000005   , 3)
PROCESS(`e3,E3', `"ch1+","ch1-"', 2000 GeV,   45.147     , 0.002      , 3)
PROCESS(`e3,E3', `"ch2+","ch2-"', 2000 GeV,   26.5162    , 0.0006     , 3)
PROCESS(`e3,E3', `"ch1+","ch2-"', 2000 GeV,    4.21267   , 0.00009    , 3)
PROCESS(`e3,E3', `h,h'          , 2000 GeV,    0.00012422, 0.00000003 , 3)
PROCESS(`e3,E3', `h,HH'         , 2000 GeV,    0.0051669 , 0.0000003  , 3)
PROCESS(`e3,E3', `HH,HH'        , 2000 GeV,    0.079301  , 0.000006   , 3)
PROCESS(`e3,E3', `HA,HA'        , 2000 GeV,    0.079758  , 0.000006   , 3)
PROCESS(`e3,E3', `Z,h'          , 2000 GeV,    3.1802    , 0.0001     , 3)
PROCESS(`e3,E3', `Z,HH'         , 2000 GeV,    4.6706    , 0.0003     , 3)
PROCESS(`e3,E3', `Z,HA'         , 2000 GeV,    4.6821    , 0.0003     , 3)
PROCESS(`e3,E3', `HA,h'         , 2000 GeV,    0.0051434 , 0.0000003  , 3)
PROCESS(`e3,E3', `HA,HH'        , 2000 GeV,    1.48793   , 0.00009    , 3)
PROCESS(`e3,E3', `Hp,Hm'        , 2000 GeV,    5.2344    , 0.0002     , 3)
! ------------------------------------------------------------------------
END_TESTSUITE 
