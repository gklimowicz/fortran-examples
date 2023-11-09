dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-bb.m4 > mssm_ext-bb.sin
dnl   whizard -r mssm_ext-bb.sin
dnl
dnl ----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-bb.m4 1764 2010-02-11 23:33:52Z jr_reuter $, mssm_ext_bb_)
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
!
seed = 0
!
show (al_h,mu_h,tanb_h)
! ------------------------------------------------------------------------
iterations = 3:200000
! ------------------------------------------------------------------------
PROCESS(`b,B', `neu1,neu1'      ,  500 GeV,   6.078898   , 0.000004   , 3)
PROCESS(`b,B', `neu1,neu2'      ,  500 GeV,  27.52342    , 0.00002    , 3)
PROCESS(`b,B', `neu1,neu3'      ,  500 GeV,  11.191450   , 0.000004   , 3)
PROCESS(`b,B', `neu1,neu4'      ,  500 GeV,   4.487316   , 0.000001   , 3)
PROCESS(`b,B', `neu2,neu2'      ,  500 GeV,  31.52604    , 0.00001    , 3)
PROCESS(`b,B', `"ch1+","ch1-"'  ,  500 GeV, 162.816      , 0.005      , 3)
PROCESS(`b,B', `h,h'            ,  500 GeV,   0.79711    , 0.00004    , 3)
dnl Bad integration in WHIZARD 1/truncated error in the CATPISS paper
PROCESS(`b,B', `Z,h'            ,  500 GeV,  30.487      , 0.001      , 5)
PROCESS(`b,B', `Z,HH'           ,  500 GeV,  50.838      , 0.001      , 3)
PROCESS(`b,B', `Z,HA'           ,  500 GeV,  52.025      , 0.001      , 3)
PROCESS(`b,B', `sb1,sb1c'       , 2000 GeV, 896.92       , 0.03       , 3)
PROCESS(`b,B', `sb2,sb2c'       , 2000 GeV, 933.08       , 0.03       , 3)
PROCESS(`b,B', `sb1,sb2c'       , 2000 GeV, 742.46       , 0.02       , 4)
PROCESS(`b,B', `st1,st1c'       , 2000 GeV, 475.02       , 0.02       , 3)
PROCESS(`b,B', `st2,st2c'       , 2000 GeV, 178.057      , 0.007      , 3)
PROCESS(`b,B', `st1,st2c'       , 2000 GeV,  50.581      , 0.002      , 3)
PROCESS(`b,B', `neu1,neu1'      , 2000 GeV,   0.096786   , 0.000002   , 3)
PROCESS(`b,B', `neu1,neu2'      , 2000 GeV,   0.445637   , 0.000009   , 3)
PROCESS(`b,B', `neu1,neu3'      , 2000 GeV,   0.1367346  , 0.0000007  , 3)
PROCESS(`b,B', `neu1,neu4'      , 2000 GeV,   0.1064429  , 0.0000007  , 3)
PROCESS(`b,B', `neu2,neu2'      , 2000 GeV,   3.54561    , 0.00009    , 3)
PROCESS(`b,B', `neu2,neu3'      , 2000 GeV,   0.928660   , 0.000006   , 3)
PROCESS(`b,B', `neu2,neu4'      , 2000 GeV,   1.08823    , 0.00001    , 3)
PROCESS(`b,B', `neu3,neu3'      , 2000 GeV,   0.264224   , 0.000009   , 3)
PROCESS(`b,B', `neu3,neu4'      , 2000 GeV,   2.78541    , 0.00003    , 3)
PROCESS(`b,B', `neu4,neu4'      , 2000 GeV,   0.46432    , 0.00001    , 3)
PROCESS(`b,B', `"ch1+","ch1-"'  , 2000 GeV,  13.9123     , 0.0006     , 3)
PROCESS(`b,B', `"ch2+","ch2-"'  , 2000 GeV, 104.774      , 0.002      , 3)
PROCESS(`b,B', `"ch1+","ch2-"'  , 2000 GeV,   6.78942    , 0.00009    , 3)
PROCESS(`b,B', `h,h'            , 2000 GeV,   0.00076246 , 0.00000001 , 3)
PROCESS(`b,B', `h,HH'           , 2000 GeV,   0.061079   , 0.000003   , 3)
PROCESS(`b,B', `HH,HH'          , 2000 GeV,   1.18500    , 0.00009    , 3)
PROCESS(`b,B', `HA,HA'          , 2000 GeV,   1.19373    , 0.00009    , 3)
PROCESS(`b,B', `HA,h'           , 2000 GeV,   0.076825   , 0.000004   , 3)
PROCESS(`b,B', `HA,HH'          , 2000 GeV,   2.4064     , 0.0001     , 3)
PROCESS(`b,B', `Z,h'            , 2000 GeV,   1.78212    , 0.00005    , 3)
PROCESS(`b,B', `Z,HH'           , 2000 GeV,  16.9839     , 0.0008     , 3)
dnl Upward fluctuation in WHIZARD 1, hence modified error
PROCESS(`b,B', `Z,HA'           , 2000 GeV,  17.0182     , 2*0.0008   , 4)
PROCESS(`b,B', `Hp,Hm'          , 2000 GeV,   2.31882    , 0.00009    , 3)
! ------------------------------------------------------------------------
END_TESTSUITE
