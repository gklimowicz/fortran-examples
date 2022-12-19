dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-bt.m4 > nmssm_ext-bt.sin
dnl   whizard -r nmssm_ext-bt.sin
dnl
dnl ----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-bt.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_bt_)
! ------------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                  -----
! ------------------------------------------------------------------------
model = NMSSM
read_slha ("nmssm.slha")
?vis_history = false
! ------------------------------------------------------------------------
helicity_selection_threshold = 1.E7
helicity_selection_cutoff = 20000
!show(real)
?fatal_beam_decay = false
! ------------------------------------------------------------------------
iterations = 3:50000
! ------------------------------------------------------------------------
PROCESS(`b,T', `se1,SN11'	,  3000 GeV,  2.0727263E+00,  1.74E-03,  3)
PROCESS(`b,T', `smu1,SN21'	,  3000 GeV,  2.0730556E+00,  1.74E-03,  3)
PROCESS(`b,T', `stau1,SN31'	,  3000 GeV,  1.0240231E+00,  8.48E-04,  3)
PROCESS(`b,T', `stau2,SN31'	,  3000 GeV,  1.0536333E+00,  8.81E-04,  3)
PROCESS(`b,T', `sd1,su1c'	,  3000 GeV,  2.5559031E+00,  2.14E-03,  3)
PROCESS(`b,T', `ss1,sc1c'	,  3000 GeV,  2.5594234E+00,  2.11E-03,  3)
PROCESS(`b,T', `sb1,st1c'	,  3000 GeV,  3.0220812E+02,  2.54E-01,  3)
PROCESS(`b,T', `sb2,st2c'	,  3000 GeV,  2.0255609E+02,  1.59E-01,  3)
PROCESS(`b,T', `sb1,st2c'	,  3000 GeV,  2.1366914E+02,  1.68E-01,  3)
PROCESS(`b,T', `sb2,st1c'	,  3000 GeV,  2.9504922E+02,  2.40E-01,  3)
PROCESS(`b,T', `"ch1-",neu1'	,  3000 GeV,  2.8882343E-01,  1.02E-04,  3)
PROCESS(`b,T', `"ch1-",neu2'	,  3000 GeV,  2.4993128E+00,  9.48E-04,  3)
PROCESS(`b,T', `"ch1-",neu3'	,  3000 GeV,  2.9919995E+00,  1.70E-03,  3)
PROCESS(`b,T', `"ch1-",neu4'	,  3000 GeV,  3.8485118E+00,  1.62E-03,  3)
PROCESS(`b,T', `"ch1-",neu5'	,  3000 GeV,  2.9211304E+00,  1.81E-03,  3)
PROCESS(`b,T', `"ch2-",neu1'	,  3000 GeV,  2.8173847E+00,  9.54E-04,  3)
PROCESS(`b,T', `"ch2-",neu2'	,  3000 GeV,  3.8842796E+01,  1.69E-02,  3)
PROCESS(`b,T', `"ch2-",neu3'	,  3000 GeV,  5.0155712E+00,  2.31E-03,  3)
PROCESS(`b,T', `"ch2-",neu4'	,  3000 GeV,  1.8425725E+01,  1.45E-02,  3)
PROCESS(`b,T', `"ch2-",neu5'	,  3000 GeV,  2.0593468E+01,  1.56E-02,  3)
PROCESS(`b,T', `Hm,h01'		,  3000 GeV,  3.5677510E+01,  4.71E-02,  3)
PROCESS(`b,T', `Hm,h02'		,  3000 GeV,  8.9384651E+00,  9.83E-03,  3)
PROCESS(`b,T', `Hm,h03'		,  3000 GeV,  4.6962074E-01,  6.82E-04,  3)
PROCESS(`b,T', `Hm,A01'		,  3000 GeV,  2.0539684E+00,  2.70E-03,  4)
PROCESS(`b,T', `Hm,A02'		,  3000 GeV,  3.8865176E-01,  5.77E-04,  3)
PROCESS(`b,T', `Wm,h01'		,  3000 GeV,  2.1084693E+02,  4.38E-01,  3)
PROCESS(`b,T', `Wm,h02'		,  3000 GeV,  7.4807976E+00,  9.39E-03,  3)
PROCESS(`b,T', `Wm,h03'		,  3000 GeV,  1.6053746E+02,  5.28E+01,  3)
PROCESS(`b,T', `Wm,A01'		,  3000 GeV,  1.5279928E+00,  1.05E-03,  3)
PROCESS(`b,T', `Wm,A02'		,  3000 GeV,  1.7114660E+02,  4.93E+01,  3)
PROCESS(`b,T', `Hm,Z'		,  3000 GeV,  6.9312808E+01,  6.46E-02,  4)
! ------------------------------------------------------------------------
END_TESTSUITE
