dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-tn.m4 > mssm_ext-tn.sin
dnl   whizard -r mssm_ext-tn.sin
dnl
dnl -------------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-tn.m4 -1   $, mssm_ext_tn_)
! ---------------------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                           -----
! ---------------------------------------------------------------------------------
! ----- Note that ALL these cross sections from all three different programs
! ----- in the CATPISS paper are too small by a factor of two because of the
! ----- wrong chirality assignment for the neutrinos in all three programs.
! ---------------------------------------------------------------------------------
model = MSSM
read_slha ("sps1a.slha")
?vis_history = false
! ---------------------------------------------------------------------------------
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
! ---------------------------------------------------------------------------------
iterations = 3:200000
! ---------------------------------------------------------------------------------
PROCESS(`e3,N3', `stau1,SN31'  ,  500 GeV,  2*84.129     ,        2*0.003      , 3)
PROCESS(`e3,N3', `stau2,SN31'  ,  500 GeV, 2*139.852     ,        2*0.006      , 3)
PROCESS(`e3,N3', `"ch1-",neu1' ,  500 GeV, 2*146.265     ,        2*0.004      , 3)
PROCESS(`e3,N3', `"ch1-",neu2' ,  500 GeV,  2*56.217     ,        2*0.001      , 3)
PROCESS(`e3,N3', `"ch2-",neu1' ,  500 GeV,   2*7.52316   ,        2*0.00006    , 3)
PROCESS(`e3,N3', `Wm,h'        ,  500 GeV, 2*133.478     ,        2*0.005      , 3)
PROCESS(`e3,N3', `Wm,HH'       ,  500 GeV,  2*57.989     ,        2*0.002      , 3)
PROCESS(`e3,N3', `Wm,HA'       ,  500 GeV,  2*58.583     ,        2*0.002      , 3)
PROCESS(`e3,N3', `Z,Hm'        ,  500 GeV,  2*17.9860    ,        2*0.0005     , 3)
PROCESS(`e3,N3', `stau1,SN31'  , 2000 GeV,  2*12.2724    ,        2*0.0004     , 3)
PROCESS(`e3,N3', `stau2,SN31'  , 2000 GeV,  2*61.463     ,        2*0.003      , 3)
PROCESS(`e3,N3', `"ch1-",neu1' , 2000 GeV,  2*21.3863    ,        2*0.0009     , 3)
PROCESS(`e3,N3', `"ch1-",neu2' , 2000 GeV,  2*21.336     ,        2*0.001      , 4)
PROCESS(`e3,N3', `"ch1-",neu3' , 2000 GeV,   2*2.2046    ,        2*0.0002     , 3)
PROCESS(`e3,N3', `"ch1-",neu4' , 2000 GeV,   2*3.44365   ,        2*0.00007    , 3)
PROCESS(`e3,N3', `"ch2-",neu1' , 2000 GeV,   2*1.95691   ,        2*0.00006    , 4)
PROCESS(`e3,N3', `"ch2-",neu2' , 2000 GeV,   2*3.49538   ,        2*0.00006    , 3)
PROCESS(`e3,N3', `"ch2-",neu3' , 2000 GeV,  2*25.8690    ,        2*0.0007     , 3)
PROCESS(`e3,N3', `"ch2-",neu4' , 2000 GeV,  2*23.1989    ,        2*0.0008     , 3)
PROCESS(`e3,N3', `Hm,h'        , 2000 GeV,   2*0.0024223 ,        2*0.0000001  , 3)
PROCESS(`e3,N3', `Hm,HH'       , 2000 GeV,   2*4.8560    ,        2*0.0003     , 3)
PROCESS(`e3,N3', `Hm,HA'       , 2000 GeV,   2*4.8578    ,        2*0.0003     , 3)
PROCESS(`e3,N3', `Wm,h'        , 2000 GeV,   2*7.5376    ,        2*0.0003     , 3)
PROCESS(`e3,N3', `Wm,HH'       , 2000 GeV,   2*8.5432    ,        2*0.0004     , 3)
PROCESS(`e3,N3', `Wm,HA'       , 2000 GeV,   2*8.5672    ,        2*0.0004     , 4)
PROCESS(`e3,N3', `Z,Hm'        , 2000 GeV,  2*13.9881    ,        2*0.0006     , 3)
! ---------------------------------------------------------------------------------
END_TESTSUITE 
