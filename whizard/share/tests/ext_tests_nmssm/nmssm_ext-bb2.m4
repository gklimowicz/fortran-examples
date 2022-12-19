dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-bb2.m4 > nmssm_ext-bb2.sin
dnl   whizard -r nmssm_ext-bb2.sin
dnl
dnl ---------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-bb2.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_bb2_)
! -----------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                 -----
! -----------------------------------------------------------------------
model = NMSSM
read_slha ("nmssm.slha")
?vis_history = false
! -----------------------------------------------------------------------
helicity_selection_threshold = 1.E7
helicity_selection_cutoff = 20000
show(real)
! -----------------------------------------------------------------------
iterations = 3:200000
! -----------------------------------------------------------------------
PROCESS(`b,B', `n1,N1'	      ,  5000 GeV,  4.7765207E-01,  1.85E-04,  3)
PROCESS(`b,B', `n2,N2'	      ,  5000 GeV,  4.7749355E-01,  1.86E-04,  3)
PROCESS(`b,B', `n3,N3'	      ,  5000 GeV,  4.7759490E-01,  1.85E-04,  3)
PROCESS(`b,B', `t,T'	      ,  5000 GeV,  5.3007275E+04,  8.66E+01,  3)
PROCESS(`b,B', `Wp,Wm'	      ,  5000 GeV,  1.5071085E+02,  1.19E-01,  3)
PROCESS(`b,B', `Z,Z'	      ,  5000 GeV,  1.1978725E+01,  7.52E-03,  3)
PROCESS(`b,B', `se1,SE1'      ,  5000 GeV,  1.5979867E-01,  7.86E-05,  3)
PROCESS(`b,B', `se2,SE2'      ,  5000 GeV,  3.8953888E-02,  1.89E-05,  3)
PROCESS(`b,B', `smu1,SMU1'    ,  5000 GeV,  1.5981656E-01,  7.89E-05,  3)
PROCESS(`b,B', `smu2,SMU2'    ,  5000 GeV,  3.8910119E-02,  1.92E-05,  3)
PROCESS(`b,B', `stau1,STAU1'  ,  5000 GeV,  3.9434404E-02,  1.94E-05,  3)
PROCESS(`b,B', `stau2,STAU2'  ,  5000 GeV,  4.1068162E-02,  2.03E-05,  3)
PROCESS(`b,B', `stau1,STAU2'  ,  5000 GeV,  5.9068465E-02,  2.92E-05,  3)
PROCESS(`b,B', `sn11,SN11'    ,  5000 GeV,  2.3699732E-01,  1.15E-04,  3)
PROCESS(`b,B', `sn21,SN21'    ,  5000 GeV,  2.3660115E-01,  1.16E-04,  3)
PROCESS(`b,B', `sn31,SN31'    ,  5000 GeV,  2.3662859E-01,  1.17E-04,  3)
PROCESS(`b,B', `su1,su1c'     ,  5000 GeV,  3.8403744E+01,  5.25E-04,  3)
PROCESS(`b,B', `su2,su2c'     ,  5000 GeV,  3.8024882E+01,  5.47E-04,  3)
PROCESS(`b,B', `sc1,sc1c'     ,  5000 GeV,  3.8403366E+01,  5.53E-02,  3)
PROCESS(`b,B', `sc2,sc2c'     ,  5000 GeV,  3.8025307E+01,  5.47E-04,  3)
PROCESS(`b,B', `st1,st1c'     ,  5000 GeV,  1.2576918E+02,  1.29E-02,  3)
PROCESS(`b,B', `st2,st2c'     ,  5000 GeV,  8.1485236E+01,  9.51E-03,  3)
PROCESS(`b,B', `st1,st2c'     ,  5000 GeV,  2.0233286E+01,  2.17E-02,  3)
PROCESS(`b,B', `sd1,sd1c'     ,  5000 GeV,  3.8421662E+01,  5.53E-04,  3)
PROCESS(`b,B', `sd2,sd2c'     ,  5000 GeV,  3.7982339E+01,  5.47E-04,  3)
PROCESS(`b,B', `ss1,ss1c'     ,  5000 GeV,  3.8421769E+01,  5.52E-04,  3)
PROCESS(`b,B', `ss2,ss2c'     ,  5000 GeV,  3.7982339E+01,  5.47E-04,  3)
PROCESS(`b,B', `sb1,sb1c'     ,  5000 GeV,  2.8019022E+02,  1.65E-02,  3)
PROCESS(`b,B', `sb2,sb2c'     ,  5000 GeV,  2.7761095E+02,  1.46E-02,  3)
PROCESS(`b,B', `sb1,sb2c'     ,  5000 GeV,  2.1167281E+02,  1.51E-02,  3)
PROCESS(`b,B', `neu1,neu1'    ,  5000 GeV,  2.0906104E-02,  9.81E-06,  3)
PROCESS(`b,B', `neu1,neu2'    ,  5000 GeV,  1.1688598E-02,  4.93E-06,  3)
PROCESS(`b,B', `neu1,neu3'    ,  5000 GeV,  1.2569256E-01,  6.02E-05,  3)
PROCESS(`b,B', `neu1,neu4'    ,  5000 GeV,  6.2096698E-03,  2.57E-06,  3)
PROCESS(`b,B', `neu1,neu5'    ,  5000 GeV,  1.7448150E-02,  4.96E-06,  3)
PROCESS(`b,B', `neu2,neu2'    ,  5000 GeV,  8.4919673E-03,  3.14E-06,  3)
PROCESS(`b,B', `neu2,neu3'    ,  5000 GeV,  6.6996609E-02,  3.09E-05,  3)
PROCESS(`b,B', `neu2,neu4'    ,  5000 GeV,  4.3768934E-02,  1.39E-05,  3)
PROCESS(`b,B', `neu2,neu5'    ,  5000 GeV,  4.0873055E-02,  1.40E-05,  3)
PROCESS(`b,B', `neu3,neu3'    ,  5000 GeV,  6.4260693E-01,  3.06E-04,  3)
PROCESS(`b,B', `neu3,neu4'    ,  5000 GeV,  9.9687798E-02,  4.12E-05,  3)
PROCESS(`b,B', `neu3,neu5'    ,  5000 GeV,  3.7698011E-02,  1.11E-05,  3)
PROCESS(`b,B', `neu4,neu4'    ,  5000 GeV,  1.9844804E-03,  9.56E-07,  3)
PROCESS(`b,B', `neu4,neu5'    ,  5000 GeV,  8.4493230E-01,  2.39E-04,  3)
PROCESS(`b,B', `neu5,neu5'    ,  5000 GeV,  2.4990446E-03,  9.29E-07,  3)
PROCESS(`b,B', `"ch1+","ch1-"',  5000 GeV,  2.4757316E+00,  1.59E-03,  3)
PROCESS(`b,B', `"ch2+","ch2-"',  5000 GeV,  2.4105629E+01,  1.22E-02,  3)
PROCESS(`b,B', `"ch1+","ch2-"',  5000 GeV,  5.8246733E-01,  2.34E-04,  3)
PROCESS(`b,B', `sgl,sgl'      ,  5000 GeV,  3.4542771E+02,  1.68E-02,  3)
PROCESS(`b,B', `Z,h01'        ,  5000 GeV,  2.6781510E-01,  2.46E-04,  3)
PROCESS(`b,B', `Z,h02'        ,  5000 GeV,  5.4636957E-03,  5.91E-06,  3)
PROCESS(`b,B', `Z,h03'        ,  5000 GeV,  1.0978407E-01,  9.33E-05,  3)
PROCESS(`b,B', `A01,h01'      ,  5000 GeV,  9.0211307E-04,  8.08E-07,  3)
PROCESS(`b,B', `A01,h02'      ,  5000 GeV,  9.7361647E-05,  8.94E-08,  3)
PROCESS(`b,B', `A01,h03'      ,  5000 GeV,  1.4450248E-03,  1.58E-06,  3)
PROCESS(`b,B', `A02,h01'      ,  5000 GeV,  3.9552475E-04,  2.38E-07,  3)
PROCESS(`b,B', `A02,h02'      ,  5000 GeV,  2.4803214E-03,  2.34E-06,  3)
PROCESS(`b,B', `A02,h03'      ,  5000 GeV,  1.3573632E-01,  1.53E-04,  3)
PROCESS(`b,B', `Hp,Hm'        ,  5000 GeV,  3.1027993E+00,  2.56E-03,  4)
PROCESS(`b,B', `Z,A01'        ,  5000 GeV,  2.0097950E-03,  4.93E-07,  3)
PROCESS(`b,B', `Z,A02'        ,  5000 GeV,  1.0928482E-01,  9.29E-05,  3)
PROCESS(`b,B', `h01,h01'      ,  5000 GeV,  5.0737808E-05,  4.02E-08,  3)
PROCESS(`b,B', `h01,h02'      ,  5000 GeV,  4.0581740E-04,  3.65E-07,  3)
PROCESS(`b,B', `h01,h03'      ,  5000 GeV,  3.2958899E-04,  2.14E-07,  3)
PROCESS(`b,B', `h02,h02'      ,  5000 GeV,  9.5993746E-05,  8.68E-08,  3)
PROCESS(`b,B', `h02,h03'      ,  5000 GeV,  1.8444098E-03,  1.68E-06,  3)
PROCESS(`b,B', `h03,h03'      ,  5000 GeV,  1.5150707E-04,  1.46E-07,  3)
PROCESS(`b,B', `A01,A01'      ,  5000 GeV,  4.8460825E-07,  3.52E-10,  3)
PROCESS(`b,B', `A01,A02'      ,  5000 GeV,  2.0853402E-04,  1.72E-07,  3)
PROCESS(`b,B', `A02,A02'      ,  5000 GeV,  1.3506794E-04,  1.45E-07,  3)
! ----------------------------------------------------------------------
END_TESTSUITE
