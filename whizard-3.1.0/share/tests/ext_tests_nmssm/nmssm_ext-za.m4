dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-za.m4 > nmssm_ext-za.sin
dnl   whizard -r nmssm_ext-za.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-za.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_za_)
! -------------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                   -----
! -------------------------------------------------------------------------
model = NMSSM
read_slha ("nmssm.slha")
?vis_history = false
! -------------------------------------------------------------------------
helicity_selection_threshold = 1.E7
helicity_selection_cutoff = 20000
!show(real)
! -------------------------------------------------------------------------
iterations = 3:50000
! -------------------------------------------------------------------------
PROCESS(`A,Z', `t,T'		,  3000 GeV,  4.7148089E+02,  1.62E-01,  3)
PROCESS(`A,Z', `Wp,Wm'		,  3000 GeV,  2.3815974E+05,  3.30E+02,  3)
PROCESS(`A,Z', `se1,SE1'	,  3000 GeV,  4.0288953E+00,  5.99E-03,  3)
PROCESS(`A,Z', `se2,SE2'	,  3000 GeV,  3.0902050E+00,  4.63E-03,  3)
PROCESS(`A,Z', `smu1,SMU1'	,  3000 GeV,  4.0221098E+00,  6.04E-03,  3)
PROCESS(`A,Z', `smu2,SMU2'	,  3000 GeV,  3.0939320E+00,  4.57E-03,  3)
PROCESS(`A,Z', `stau1,STAU1'	,  3000 GeV,  9.2711947E-03,  1.42E-05,  3)
PROCESS(`A,Z', `stau2,STAU2'	,  3000 GeV,  2.2349162E-02,  3.30E-05,  3)
PROCESS(`A,Z', `stau1,STAU2'	,  3000 GeV,  3.5742883E+00,  5.31E-03,  3)
PROCESS(`A,Z', `su1,su1c'	,  3000 GeV,  3.9894974E+00,  1.37E-03,  3)
PROCESS(`A,Z', `su2,su2c'	,  3000 GeV,  8.1671773E-01,  2.81E-04,  3)
PROCESS(`A,Z', `sc1,sc1c'	,  3000 GeV,  3.9905985E+00,  1.35E-03,  3)
PROCESS(`A,Z', `sc2,sc2c'	,  3000 GeV,  8.1688763E-01,  2.91E-04,  3)
PROCESS(`A,Z', `st1,st1c'	,  3000 GeV,  3.5578476E-01,  1.86E-04,  3)
PROCESS(`A,Z', `st2,st2c'	,  3000 GeV,  2.7587370E-01,  8.06E-05,  3)
PROCESS(`A,Z', `st1,st2c'	,  3000 GeV,  3.8472182E+01,  3.86E-02,  3)
PROCESS(`A,Z', `sd1,sd1c'	,  3000 GeV,  1.4995206E+00,  5.08E-04,  3)
PROCESS(`A,Z', `sd2,sd2c'	,  3000 GeV,  5.1052372E-02,  1.77E-05,  3)
PROCESS(`A,Z', `ss1,ss1c'	,  3000 GeV,  1.4988856E+00,  4.89E-04,  3)
PROCESS(`A,Z', `ss2,ss2c'	,  3000 GeV,  5.1062779E-02,  1.71E-05,  3)
PROCESS(`A,Z', `sb1,sb1c'	,  3000 GeV,  1.8633006E-01,  6.37E-05,  3)
PROCESS(`A,Z', `sb2,sb2c'	,  3000 GeV,  3.2120074E-01,  1.09E-04,  3)
PROCESS(`A,Z', `sb1,sb2c'	,  3000 GeV,  5.2511465E-01,  1.76E-04,  3)
PROCESS(`A,Z', `"ch1+","ch1-"'	,  3000 GeV,  3.7103808E+02,  9.87E-02,  3)
PROCESS(`A,Z', `"ch1+","ch2-"'	,  3000 GeV,  5.6722565E+01,  3.19E-02,  3)
PROCESS(`A,Z', `"ch2+","ch2-"'	,  3000 GeV,  2.8241363E+01,  7.79E-03,  3)
PROCESS(`A,Z', `Hp,Hm'		,  3000 GeV,  1.5728086E+00,  4.75E-04,  3)
PROCESS(`A,Z', `t,T'		,  5000 GeV,  2.0106882E+02,  4.20E-02,  3)
PROCESS(`A,Z', `Wp,Wm'		,  5000 GeV,  2.3891931E+05,  2.17E+02,  3)
PROCESS(`A,Z', `se1,SE1'	,  5000 GeV,  1.5236563E+00,  1.44E-03,  3)
PROCESS(`A,Z', `se2,SE2'	,  5000 GeV,  1.1724302E+00,  1.11E-03,  3)
PROCESS(`A,Z', `smu1,SMU1'	,  5000 GeV,  1.5207938E+00,  1.44E-03,  3)
PROCESS(`A,Z', `smu2,SMU2'	,  5000 GeV,  1.1714825E+00,  1.11E-03,  3)
PROCESS(`A,Z', `stau1,STAU1'	,  5000 GeV,  3.4901967E-03,  3.35E-06,  3)
PROCESS(`A,Z', `stau2,STAU2'	,  5000 GeV,  8.5176294E-03,  8.06E-06,  3)
PROCESS(`A,Z', `stau1,STAU2'	,  5000 GeV,  1.3465546E+00,  1.25E-03,  3)
PROCESS(`A,Z', `su1,su1c'	,  5000 GeV,  2.0959712E+00,  9.74E-04,  3)
PROCESS(`A,Z', `su2,su2c'	,  5000 GeV,  4.2917180E-01,  1.97E-04,  3)
PROCESS(`A,Z', `sc1,sc1c'	,  5000 GeV,  2.0982522E+00,  9.56E-04,  3)
PROCESS(`A,Z', `sc2,sc2c'	,  5000 GeV,  4.2898610E-01,  1.97E-04,  3)
PROCESS(`A,Z', `st1,st1c'	,  5000 GeV,  1.8461567E-01,  1.09E-04,  3)
PROCESS(`A,Z', `st2,st2c'	,  5000 GeV,  1.3603745E-01,  4.86E-05,  3)
PROCESS(`A,Z', `st1,st2c'	,  5000 GeV,  1.3518376E+01,  7.64E-03,  3)
PROCESS(`A,Z', `sd1,sd1c'	,  5000 GeV,  7.8737267E-01,  3.62E-04,  3)
PROCESS(`A,Z', `sd2,sd2c'	,  5000 GeV,  2.6811458E-02,  1.24E-05,  3)
PROCESS(`A,Z', `ss1,ss1c'	,  5000 GeV,  7.8722473E-01,  3.63E-04,  3)
PROCESS(`A,Z', `ss2,ss2c'	,  5000 GeV,  2.6806936E-02,  1.24E-05,  3)
PROCESS(`A,Z', `sb1,sb1c'	,  5000 GeV,  9.7963012E-02,  4.47E-05,  3)
PROCESS(`A,Z', `sb2,sb2c'	,  5000 GeV,  1.6855604E-01,  7.61E-05,  3)
PROCESS(`A,Z', `sb1,sb2c'	,  5000 GeV,  2.7516986E-01,  1.26E-04,  3)
PROCESS(`A,Z', `"ch1+","ch1-"'	,  5000 GeV,  1.5833221E+02,  2.48E-02,  3)
PROCESS(`A,Z', `"ch1+","ch2-"'	,  5000 GeV,  2.6033927E+01,  8.81E-03,  3)
PROCESS(`A,Z', `"ch2+","ch2-"'	,  5000 GeV,  1.3523317E+01,  2.00E-03,  3)
PROCESS(`A,Z', `Hp,Hm'		,  5000 GeV,  7.7371343E-01,  2.48E-04,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
