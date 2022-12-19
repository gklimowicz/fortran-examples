dnl Process this with
dnl
dnl   m4 testuite.m4 testsuite_sample.m4 > testsuite_sample.sin
dnl   whizard -r testsuite_sample
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: testsuite_sample.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_zz2_)
! -------------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                   -----
! -------------------------------------------------------------------------
model = NMSSM
read_slha ("nmssm.slha")
?vis_history = false
! -------------------------------------------------------------------------
!helicity_selection_threshold = 1.E7
!helicity_selection_cutoff = 20000
!show(real)
! -------------------------------------------------------------------------
iterations = 3:50000
! -------------------------------------------------------------------------
PROCESS(`Z,Z', `t,T'		,  5000 GeV,  6.3663452E+02,  4.10E-01,  3)
PROCESS(`Z,Z', `Wp,Wm'	   	,  5000 GeV,  6.0641622E+05,  8.23E+02,  3)
PROCESS(`Z,Z', `Z,Z'	   	,  5000 GeV,  3.0759878E+01,  2.28E-02,  3)
PROCESS(`Z,Z', `se1,SE1'	,  5000 GeV,  4.6827560E-01,  3.35E-04,  3)
PROCESS(`Z,Z', `se2,SE2'	,  5000 GeV,  2.8928802E-01,  2.06E-04,  3)
PROCESS(`Z,Z', `smu1,SMU1'	,  5000 GeV,  4.6861168E-01,  3.34E-04,  3)
PROCESS(`Z,Z', `smu2,SMU2'	,  5000 GeV,  2.8899838E-01,  2.07E-04,  3)
PROCESS(`Z,Z', `stau1,STAU1'	,  5000 GeV,  3.7406350E-01,  2.67E-04,  3)
PROCESS(`Z,Z', `stau2,STAU2'	,  5000 GeV,  3.7479827E-01,  2.65E-04,  3)
PROCESS(`Z,Z', `stau1,STAU2'	,  5000 GeV,  5.6162694E-03,  3.96E-06,  3)
PROCESS(`Z,Z', `sn11,SN11'	,  5000 GeV,  5.2441341E+00,  5.02E-03,  3)
PROCESS(`Z,Z', `sn21,SN21'	,  5000 GeV,  5.2533816E+00,  5.01E-03,  3)
PROCESS(`Z,Z', `sn31,SN31'	,  5000 GeV,  5.2487743E+00,  5.00E-03,  3)
PROCESS(`Z,Z', `su1,su1c'	,  5000 GeV,  2.3849179E+00,  8.17E-04,  3)
PROCESS(`Z,Z', `su2,su2c'	,  5000 GeV,  1.4935063E-01,  5.12E-05,  3)
PROCESS(`Z,Z', `sc1,sc1c'	,  5000 GeV,  2.3858442E+00,  8.08E-04,  3)
PROCESS(`Z,Z', `sc2,sc2c'	,  5000 GeV,  1.4919940E-01,  5.00E-05,  3)
PROCESS(`Z,Z', `st1,st1c'	,  5000 GeV,  6.1429364E+01,  1.14E-02,  3)
PROCESS(`Z,Z', `st2,st2c'	,  5000 GeV,  1.2138685E+02,  9.05E-02,  3)
PROCESS(`Z,Z', `st1,st2c'	,  5000 GeV,  2.1936014E+00,  9.64E-04,  3)
PROCESS(`Z,Z', `sd1,sd1c'	,  5000 GeV,  5.1618647E+00,  1.77E-03,  3)
PROCESS(`Z,Z', `sd2,sd2c'	,  5000 GeV,  2.1004911E-02,  7.36E-06,  3)
PROCESS(`Z,Z', `ss1,ss1c'	,  5000 GeV,  5.1639219E+00,  1.75E-03,  3)
PROCESS(`Z,Z', `ss2,ss2c'	,  5000 GeV,  2.0978904E-02,  7.40E-06,  3)
PROCESS(`Z,Z', `sb1,sb1c'	,  5000 GeV,  1.4029722E+00,  4.75E-04,  3)
PROCESS(`Z,Z', `sb2,sb2c'	,  5000 GeV,  3.4385312E+00,  1.19E-03,  3)
PROCESS(`Z,Z', `sb1,sb2c'	,  5000 GeV,  1.1244023E+00,  3.83E-04,  3)
PROCESS(`Z,Z', `neu1,neu1'	,  5000 GeV,  7.1985290E-01,  5.09E-04,  3)
PROCESS(`Z,Z', `neu1,neu2'	,  5000 GeV,  9.8048153E-01,  7.99E-04,  3)
PROCESS(`Z,Z', `neu1,neu3'	,  5000 GeV,  1.4873207E+00,  1.11E-03,  3)
PROCESS(`Z,Z', `neu1,neu4'	,  5000 GeV,  7.2680241E+00,  4.11E-03,  3)
PROCESS(`Z,Z', `neu1,neu5'	,  5000 GeV,  1.5521212E-01,  8.45E-05,  3)
PROCESS(`Z,Z', `neu2,neu2'	,  5000 GeV,  7.2208348E+00,  6.13E-03,  3)
PROCESS(`Z,Z', `neu2,neu3'	,  5000 GeV,  3.6724409E+00,  2.97E-03,  3)
PROCESS(`Z,Z', `neu2,neu4'	,  5000 GeV,  4.3908722E+00,  2.47E-03,  3)
PROCESS(`Z,Z', `neu2,neu5'	,  5000 GeV,  2.3428498E+01,  1.34E-02,  3)
PROCESS(`Z,Z', `neu3,neu3'	,  5000 GeV,  1.5521031E+00,  1.48E-03,  3)
PROCESS(`Z,Z', `neu3,neu4'	,  5000 GeV,  8.2447296E+00,  4.58E-03,  3)
PROCESS(`Z,Z', `neu3,neu5'	,  5000 GeV,  6.3749503E+00,  3.67E-03,  3)
PROCESS(`Z,Z', `neu4,neu4'	,  5000 GeV,  3.7928006E+01,  1.88E-02,  3)
PROCESS(`Z,Z', `neu4,neu5'	,  5000 GeV,  8.9942889E-01,  6.64E-04,  3)
PROCESS(`Z,Z', `neu5,neu5'	,  5000 GeV,  2.7159161E+01,  1.38E-02,  3)
PROCESS(`Z,Z', `"ch1+","ch1-"'	,  5000 GeV,  3.4206675E+02,  1.35E-01,  3)
PROCESS(`Z,Z', `"ch1+","ch2-"'	,  5000 GeV,  7.3337377E+01,  1.97E-02,  3)
PROCESS(`Z,Z', `"ch2+","ch2-"'	,  5000 GeV,  1.7238883E+01,  8.17E-03,  3)
dnl Minimal fluctuations in the Higgs sector
PROCESS(`Z,Z', `h01,h01'	,  5000 GeV,  7.4490352E+03,  9.67E+00,  3)
PROCESS(`Z,Z', `h01,h02'	,  5000 GeV,  3.2080731E+02,  4.17E-01,  4)
PROCESS(`Z,Z', `h01,h03'	,  5000 GeV,  1.1647312E+00,  8.22E-04,  3)
PROCESS(`Z,Z', `h02,h02'	,  5000 GeV,  5.3238961E+00,  4.51E-03,  3)
PROCESS(`Z,Z', `h02,h03'	,  5000 GeV,  5.6636899E-01,  4.43E-04,  3)
PROCESS(`Z,Z', `h03,h03'	,  5000 GeV,  2.0857873E+00,  1.20E-03,  3)
PROCESS(`Z,Z', `A01,A01'	,  5000 GeV,  3.2147391E+00,  2.56E-03,  3)
PROCESS(`Z,Z', `A01,A02'	,  5000 GeV,  7.9844798E-01,  6.09E-04,  3)
PROCESS(`Z,Z', `A02,A02'	,  5000 GeV,  1.9136559E+00,  1.26E-03,  3)
PROCESS(`Z,Z', `Z,A01'		,  5000 GeV,  7.5040141E-04,  5.91E-07,  3)
PROCESS(`Z,Z', `Z,A02'	   	,  5000 GeV,  1.2193386E-01,  8.75E-05,  3)
PROCESS(`Z,Z', `Hp,Hm'	   	,  5000 GeV,  2.0506880E-01,  4.34E-05,  4)
! -------------------------------------------------------------------------
END_TESTSUITE
