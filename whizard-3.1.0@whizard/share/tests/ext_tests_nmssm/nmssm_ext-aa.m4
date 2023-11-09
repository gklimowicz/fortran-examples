dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-aa.m4 > nmssm_ext-aa.sin
dnl   whizard -r nmssm_ext-aa.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-aa.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_aa_)
! -------------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                   -----
! -------------------------------------------------------------------------
model = NMSSM
read_slha ("nmssm.slha")
?vis_history = false
! -------------------------------------------------------------------------
helicity_selection_threshold = 1.E7
helicity_selection_cutoff = 2000
! show(real)
! -------------------------------------------------------------------------
iterations = 3:50000
! -------------------------------------------------------------------------
PROCESS(`A,A', `t,T'		,  3000 GeV,  9.4255447e+01,  2.45e-02,  3)
PROCESS(`A,A', `Wp,Wm'		,  3000 GeV,  9.3864054E+04,  1.28E+02,  3)
PROCESS(`A,A', `se1, SE1'	,  3000 GeV,  1.5139082E+01,  2.28E-02,  3)
PROCESS(`A,A', `se2, SE2'	,  3000 GeV,  1.5125771E+01,  2.29E-02,  3)
PROCESS(`A,A', `smu1, SMU1'	,  3000 GeV,  1.5106579E+01,  2.27E-02,  3)
PROCESS(`A,A', `smu2, SMU2'	,  3000 GeV,  1.5173115E+01,  2.27E-02,  3)
PROCESS(`A,A', `stau1, STAU1'	,  3000 GeV,  1.5320611E+01,  2.36E-02,  3)
PROCESS(`A,A', `stau2, STAU2'	,  3000 GeV,  1.5017961E+01,  2.18E-02,  3)
PROCESS(`A,A', `su1,su1c'	,  3000 GeV,  4.0170998E+00,  1.39E-03,  3)
PROCESS(`A,A', `su2,su2c'	,  3000 GeV,  4.0164032E+00,  1.34E-03,  3)
PROCESS(`A,A', `sc1,sc1c'	,  3000 GeV,  4.0171434E+00,  1.38E-03,  3)
PROCESS(`A,A', `sc2,sc2c'	,  3000 GeV,  4.0175098E+00,  1.40E-03,  3)
PROCESS(`A,A', `st1,st1c'	,  3000 GeV,  4.7473881E+00,  2.54E-03,  3)
PROCESS(`A,A', `st2,st2c'	,  3000 GeV,  3.7282666E+00,  1.09E-03,  3)
PROCESS(`A,A', `sd1,sd1c'	,  3000 GeV,  2.5079842E-01,  8.54E-05,  3)
PROCESS(`A,A', `sd2,sd2c'	,  3000 GeV,  2.5103269E-01,  8.38E-05,  3)
PROCESS(`A,A', `ss1,ss1c'	,  3000 GeV,  2.5086139E-01,  8.64E-05,  3)
PROCESS(`A,A', `ss2,ss2c'	,  3000 GeV,  2.5069643E-01,  8.71E-05,  3)
PROCESS(`A,A', `sb1,sb1c'	,  3000 GeV,  2.5148809E-01,  8.49E-05,  3)
PROCESS(`A,A', `sb2,sb2c'	,  3000 GeV,  2.5036105E-01,  8.58E-05,  3)
PROCESS(`A,A', `"ch1+","ch1-"'	,  3000 GeV,  1.7727630E+02,  4.62E-02,  3)
PROCESS(`A,A', `"ch2+","ch2-"'	,  3000 GeV,  8.9185999E+01,  2.37E-02,  3)
PROCESS(`A,A', `Hp,Hm'		,  3000 GeV,  5.9536447E+00,  1.89E-03,  3)
PROCESS(`A,A', `t,T'		,  5000 GeV,  4.0928233E+01,  6.28E-03,  3)
PROCESS(`A,A', `Wp,Wm'		,  5000 GeV,  9.3923228E+04,  8.28E+01,  3)
PROCESS(`A,A', `se1,SE1'	,  5000 GeV,  5.7483575E+00,  5.46E-03,  3)
PROCESS(`A,A', `se2,SE2'	,  5000 GeV,  5.7456901E+00,  5.45E-03,  3)
PROCESS(`A,A', `smu1,SMU1'	,  5000 GeV,  5.7473263E+00,  5.46E-03,  3)
PROCESS(`A,A', `smu2,SMU2'	,  5000 GeV,  5.7569416E+00,  5.47E-03,  3)
PROCESS(`A,A', `stau1,STAU1'	,  5000 GeV,  5.7793061E+00,  5.54E-03,  3)
PROCESS(`A,A', `stau2,STAU2'	,  5000 GeV,  5.7190065E+00,  5.40E-03,  3)
PROCESS(`A,A', `su1,su1c'	,  5000 GeV,  2.1092243E+00,  9.77E-04,  3)
PROCESS(`A,A', `su2,su2c'	,  5000 GeV,  2.1103144E+00,  9.58E-04,  3)
PROCESS(`A,A', `sc1,sc1c'	,  5000 GeV,  2.1100129E+00,  9.56E-04,  3)
PROCESS(`A,A', `sc2,sc2c'	,  5000 GeV,  2.1119405E+00,  9.77E-04,  3)
PROCESS(`A,A', `st1,st1c'	,  5000 GeV,  2.4641905E+00,  1.45E-03,  3)
PROCESS(`A,A', `st2,st2c'	,  5000 GeV,  1.8382527E+00,  6.52E-04,  3)
PROCESS(`A,A', `sd1,sd1c'	,  5000 GeV,  1.3168148E-01,  6.04E-05,  3)
PROCESS(`A,A', `sd2,sd2c'	,  5000 GeV,  1.3177093E-01,  6.09E-05,  3)
PROCESS(`A,A', `ss1,ss1c'	,  5000 GeV,  1.3169584E-01,  6.03E-05,  3)
PROCESS(`A,A', `ss2,ss2c'	,  5000 GeV,  1.3172183E-01,  6.05E-05,  3)
PROCESS(`A,A', `sb1,sb1c'	,  5000 GeV,  1.3224558E-01,  5.99E-05,  3)
PROCESS(`A,A', `sb2,sb2c'	,  5000 GeV,  1.3135316E-01,  6.01E-05,  3)
PROCESS(`A,A', `"ch1+","ch1-"'	,  5000 GeV,  7.5745777E+01,  1.18E-02,  3)
PROCESS(`A,A', `"ch2+","ch2-"'	,  5000 GeV,  4.2769180E+01,  6.33E-03,  3)
PROCESS(`A,A', `Hp,Hm'		,  5000 GeV,  2.9244027E+00,  9.48E-04,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
