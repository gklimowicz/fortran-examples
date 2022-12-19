dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-gg.m4 > nmssm_ext-gg.sin
dnl   whizard -r nmssm_ext-gg.sin
dnl
dnl --------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-gg.m4 $, nmssm_ext_gg_)
! ----------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                -----
! ----------------------------------------------------------------------
model = NMSSM
read_slha ("nmssm.slha")
?vis_history = false
! ----------------------------------------------------------------------
helicity_selection_threshold = 1.E7
helicity_selection_cutoff = 20000
show(real)
! ----------------------------------------------------------------------
iterations = 3:50000
! ----------------------------------------------------------------------
PROCESS(`gl,gl', `t,T'	   ,  3000 GeV,	 2.5081146E+03,  1.95E+00,    3)
PROCESS(`gl,gl', `sgl,sgl' ,  3000 GeV,	 9.1938444E+03,  1.91E+00,    3)
PROCESS(`gl,gl', `su1,su1c',  3000 GeV,	 7.0232438E+01,  6.04E-03,    3)
PROCESS(`gl,gl', `su2,su2c',  3000 GeV,	 7.0202127E+01,  7.07E-03,    4)
PROCESS(`gl,gl', `sc1,sc1c',  3000 GeV,	 7.0238595E+01,  5.97E-03,    3)
PROCESS(`gl,gl', `sc2,sc2c',  3000 GeV,	 7.0211786E+01,  6.42E-03,    3)
PROCESS(`gl,gl', `st1,st1c',  3000 GeV,	 8.6198986E+01,  1.66E-02,    3)
PROCESS(`gl,gl', `st2,st2c',  3000 GeV,	 6.0369126E+01,  5.81E-03,    4)
PROCESS(`gl,gl', `sd1,sd1c',  3000 GeV,	 7.0127359E+01,  5.97E-03,    3)
PROCESS(`gl,gl', `sd2,sd2c',  3000 GeV,	 7.0163549E+01,  7.01E-03,    4)
PROCESS(`gl,gl', `ss1,ss1c',  3000 GeV,	 7.0133800E+01,  6.95E-03,    3)
PROCESS(`gl,gl', `ss2,ss2c',  3000 GeV,	 7.0182724E+01,  5.20E-03,    4)
PROCESS(`gl,gl', `sb1,sb1c',  3000 GeV,	 7.0400739E+01,  6.53E-03,    3)
PROCESS(`gl,gl', `sb2,sb2c',  3000 GeV,	 6.9920822E+01,  5.72E-03,    3)
PROCESS(`gl,gl', `t,T'	   ,  5000 GeV,	 1.1247848E+03,  8.53E-02,    3)
dnl Fluctuations in the photon polarization vectors for high energies
PROCESS(`gl,gl', `sgl,sgl' ,  5000 GeV,	 4.6574401E+03,  2.21E-01,    4)
PROCESS(`gl,gl', `su1,su1c',  5000 GeV,	 3.8821548E+01,  1.39E-03,    3)
PROCESS(`gl,gl', `su2,su2c',  5000 GeV,	 3.8810184E+01,  1.40E-03,    4)
PROCESS(`gl,gl', `sc1,sc1c',  5000 GeV,	 3.8819174E+01,  1.38E-03,    4)
PROCESS(`gl,gl', `sc2,sc2c',  5000 GeV,	 3.8807131E+01,  1.38E-03,    3)
PROCESS(`gl,gl', `st1,st1c',  5000 GeV,	 4.5800799E+01,  2.10E-03,    3)
PROCESS(`gl,gl', `st2,st2c',  5000 GeV,	 3.3557935E+01,  9.27E-04,    4)
PROCESS(`gl,gl', `sd1,sd1c',  5000 GeV,	 3.8766493E+01,  1.40E-03,    3)
PROCESS(`gl,gl', `sd2,sd2c',  5000 GeV,	 3.8793026E+01,  1.38E-03,    4)
PROCESS(`gl,gl', `ss1,ss1c',  5000 GeV,	 3.8769950E+01,  1.40E-03,    4)
PROCESS(`gl,gl', `ss2,ss2c',  5000 GeV,	 3.8790550E+01,  1.40E-03,    4)
PROCESS(`gl,gl', `sb1,sb1c',  5000 GeV,	 3.8911090E+01,  1.41E-03,    4)
PROCESS(`gl,gl', `sb2,sb2c',  5000 GeV,	 3.8650275E+01,  1.39E-03,    3)
! ----------------------------------------------------------------------
END_TESTSUITE
