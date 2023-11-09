dnl Process this with
dnl
dnl   m4 testuite.m4 mssm_ext-en.m4 > mssm_ext-en.sin
dnl   whizard -r mssm_ext-en.sin
dnl
dnl --------------------------------------------------------------------------
BEGIN_TESTSUITE($Id: mssm_ext-en.m4 -1   $, mssm_ext_en_)
! ----------------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                      -----
! ----------------------------------------------------------------------------
! ----- Note that ALL these cross sections from all three different programs
! ----- in the CATPISS paper are too small by a factor of two because of the
! ----- wrong chirality assignment for the neutrinos in all three programs.
! ----------------------------------------------------------------------------
model = MSSM
read_slha ("sps1a.slha")
?vis_history = false
! ----------------------------------------------------------------------------
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
! ----------------------------------------------------------------------------
iterations = 3:200000
! ----------------------------------------------------------------------------
PROCESS(`e1,N1', `se1,SN11'    ,  500 GeV, 2*158.694   ,        2*0.007   , 3)
PROCESS(`e1,N1', `se2,SN11'    ,  500 GeV,  2*68.513   ,        2*0.003   , 3)
PROCESS(`e1,N1', `smu1,SN21'   ,  500 GeV,  2*58.491   ,        2*0.002   , 3)
PROCESS(`e1,N1', `stau1,SN31'  ,  500 GeV,   2*8.5021  ,        2*0.0003  , 3)
PROCESS(`e1,N1', `stau2,SN31'  ,  500 GeV,  2*51.790   ,        2*0.002   , 3)
PROCESS(`e1,N1', `"ch1-",neu1' ,  500 GeV, 2*137.416   ,        2*0.004   , 3)
PROCESS(`e1,N1', `"ch1-",neu2' ,  500 GeV,  2*58.795   ,        2*0.001   , 4)
PROCESS(`e1,N1', `"ch2-",neu1' ,  500 GeV,   2*1.16072 ,        2*0.00003 , 3)
PROCESS(`e1,N1', `se1,SN11'    , 2000 GeV,  2*67.095   ,        2*0.003   , 3)
PROCESS(`e1,N1', `se2,SN11'    , 2000 GeV,   2*6.5470  ,        2*0.0004  , 3)
PROCESS(`e1,N1', `smu1,SN21'   , 2000 GeV,  2*13.8935  ,        2*0.0005  , 3)
PROCESS(`e1,N1', `stau1,SN31'  , 2000 GeV,   2*1.11690 ,        2*0.00004 , 3)
PROCESS(`e1,N1', `stau2,SN31'  , 2000 GeV,  2*12.7836  ,        2*0.0005  , 3)
PROCESS(`e1,N1', `"ch1-",neu1' , 2000 GeV,  2*21.4203  ,        2*0.0009  , 3)
PROCESS(`e1,N1', `"ch1-",neu2' , 2000 GeV,  2*21.283   ,        2*0.001   , 3)
PROCESS(`e1,N1', `"ch1-",neu3' , 2000 GeV,   2*2.26760 ,        2*0.00007 , 3)
PROCESS(`e1,N1', `"ch1-",neu4' , 2000 GeV,   2*3.51046 ,        2*0.00006 , 3)
PROCESS(`e1,N1', `"ch2-",neu1' , 2000 GeV,   2*1.73607 ,        2*0.00006 , 3)
PROCESS(`e1,N1', `"ch2-",neu2' , 2000 GeV,   2*3.61122 ,        2*0.00006 , 3)
PROCESS(`e1,N1', `"ch2-",neu3' , 2000 GeV,  2*26.9511  ,        2*0.0007  , 3)
PROCESS(`e1,N1', `"ch2-",neu4' , 2000 GeV,  2*24.0223  ,        2*0.0008  , 3)
! ----------------------------------------------------------------------------
END_TESTSUITE 
