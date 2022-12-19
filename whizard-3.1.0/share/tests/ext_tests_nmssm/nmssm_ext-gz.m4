dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-gz.m4 > nmssm_ext-gz.sin
dnl   whizard -r nmssm_ext-gz.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-gz.m4 $, nmssm_ext_gz_)
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
PROCESS(`gl,Z', `t,T'		,  3000 GeV,  2.6515166E+03,  9.53E-01,  3)
PROCESS(`gl,Z', `su1,su1c'	,  3000 GeV,  2.2437656E+01,  7.82E-03,  3)
PROCESS(`gl,Z', `su2,su2c'	,  3000 GeV,  4.5932203E+00,  1.59E-03,  3)
PROCESS(`gl,Z', `sc1,sc1c'	,  3000 GeV,  2.2427425E+01,  7.93E-03,  3)
PROCESS(`gl,Z', `sc2,sc2c'	,  3000 GeV,  4.5934159E+00,  1.55E-03,  3)
PROCESS(`gl,Z', `st1,st1c'	,  3000 GeV,  2.0002664E+00,  1.08E-03,  3)
PROCESS(`gl,Z', `st2,st2c'	,  3000 GeV,  1.5509284E+00,  4.67E-04,  3)
PROCESS(`gl,Z', `st1,st2c'	,  3000 GeV,  2.1652096E+02,  2.16E-01,  3)
PROCESS(`gl,Z', `sd1,sd1c'	,  3000 GeV,  3.3719697E+01,  1.11E-02,  3)
PROCESS(`gl,Z', `sd2,sd2c'	,  3000 GeV,  1.1478395E+00,  3.69E-04,  3)
PROCESS(`gl,Z', `ss1,ss1c'	,  3000 GeV,  3.3732569E+01,  1.08E-02,  3)
PROCESS(`gl,Z', `ss2,ss2c'	,  3000 GeV,  1.1474302E+00,  4.06E-04,  3)
PROCESS(`gl,Z', `sb1,sb1c'	,  3000 GeV,  4.1904564E+00,  1.42E-03,  3)
PROCESS(`gl,Z', `sb2,sb2c'	,  3000 GeV,  7.2280282E+00,  2.42E-03,  3)
PROCESS(`gl,Z', `sb1,sb2c'	,  3000 GeV,  1.1808205E+01,  4.29E-03,  3)
PROCESS(`gl,Z', `t,T'		,  5000 GeV,  1.1306609E+03,  4.03E-01,  3)
PROCESS(`gl,Z', `su1,su1c'	,  5000 GeV,  1.1783452E+01,  9.31E-03,  3)
PROCESS(`gl,Z', `su2,su2c'	,  5000 GeV,  2.4117499E+00,  1.91E-03,  3)
PROCESS(`gl,Z', `sc1,sc1c'	,  5000 GeV,  1.1801173E+01,  9.22E-03,  3)
PROCESS(`gl,Z', `sc2,sc2c'	,  5000 GeV,  2.4149850E+00,  1.89E-03,  3)
PROCESS(`gl,Z', `st1,st1c'	,  5000 GeV,  1.0383736E+00,  1.02E-03,  3)
PROCESS(`gl,Z', `st2,st2c'	,  5000 GeV,  7.6521954E-01,  4.53E-04,  3)
PROCESS(`gl,Z', `st1,st2c'	,  5000 GeV,  7.5870194E+01,  7.56E-02,  3)
PROCESS(`gl,Z', `sd1,sd1c'	,  5000 GeV,  1.7696149E+01,  1.43E-02,  3)
PROCESS(`gl,Z', `sd2,sd2c'	,  5000 GeV,  6.0260176E-01,  4.79E-04,  3)
PROCESS(`gl,Z', `ss1,ss1c'	,  5000 GeV,  1.7708484E+01,  1.41E-02,  3)
PROCESS(`gl,Z', `ss2,ss2c'	,  5000 GeV,  6.0299570E-01,  4.78E-04,  3)
PROCESS(`gl,Z', `sb1,sb1c'	,  5000 GeV,  2.2030422E+00,  1.73E-03,  3)
PROCESS(`gl,Z', `sb2,sb2c'	,  5000 GeV,  3.7950412E+00,  2.94E-03,  3)
PROCESS(`gl,Z', `sb1,sb2c'	,  5000 GeV,  6.1859449E+00,  4.92E-03,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
