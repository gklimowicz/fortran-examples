dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-gw.m4 > nmssm_ext-gw.sin
dnl   whizard -r nmssm_ext-gw.sin
dnl
dnl ------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-gw.m4 $, nmssm_ext_gw_)
! --------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds              -----
! --------------------------------------------------------------------
model = NMSSM
read_slha ("nmssm.slha")
?vis_history = false
! --------------------------------------------------------------------
!helicity_selection_threshold = 1.E7
!helicity_selection_cutoff = 20000
!show(real)
! --------------------------------------------------------------------
iterations = 3:50000
! --------------------------------------------------------------------
PROCESS(`gl,Wm', `sd1,su1c',  3000 GeV,  7.2549248E+01,  2.50E-02,  3)
PROCESS(`gl,Wm', `ss1,sc1c',  3000 GeV,  7.2573674E+01,  2.38E-02,  3)
PROCESS(`gl,Wm', `sb1,st1c',  3000 GeV,  1.2367371E+02,  1.14E-01,  3)
PROCESS(`gl,Wm', `sb2,st2c',  3000 GeV,  9.6998656E+01,  8.22E-02,  3)
PROCESS(`gl,Wm', `sb1,st2c',  3000 GeV,  8.7173874E+01,  7.45E-02,  3)
PROCESS(`gl,Wm', `sb2,st1c',  3000 GeV,  1.5973470E+02,  1.49E-01,  3)
PROCESS(`gl,Wm', `sd1,su1c',  5000 GeV,  3.8113282E+01,  2.30E-02,  3)
PROCESS(`gl,Wm', `ss1,sc1c',  5000 GeV,  3.8073622E+01,  3.00E-02,  3)
PROCESS(`gl,Wm', `sb1,st1c',  5000 GeV,  3.9736256E+01,  3.87E-02,  3)
PROCESS(`gl,Wm', `sb2,st2c',  5000 GeV,  4.4305077E+01,  4.15E-02,  3)
PROCESS(`gl,Wm', `sb1,st2c',  5000 GeV,  3.9313200E+01,  3.75E-02,  3)
PROCESS(`gl,Wm', `sb2,st1c',  5000 GeV,  5.1364531E+01,  4.89E-02,  3)
! --------------------------------------------------------------------
END_TESTSUITE
