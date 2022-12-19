dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-ga.m4 > nmssm_ext-ga.sin
dnl   whizard -r nmssm_ext-ga
dnl
dnl ---------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-ga.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_ga_)
! -----------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                 -----
! -----------------------------------------------------------------------
model = NMSSM
read_slha ("nmssm.slha")
?vis_history = false
! -----------------------------------------------------------------------
helicity_selection_threshold = 1.E7
helicity_selection_cutoff = 20000
! show(real)
! -----------------------------------------------------------------------
iterations = 3:50000
! -----------------------------------------------------------------------
PROCESS(`gl,A', `t,T'	     ,  3000 GeV,   5.2971809E+02, 1.37E-01,   3)
PROCESS(`gl,A', `su1,su1c'   ,  3000 GeV,   2.2573580E+01, 7.92E-03,   3)
PROCESS(`gl,A', `su2,su2c'   ,  3000 GeV,   2.2572249E+01, 7.84E-03,   3)
PROCESS(`gl,A', `sc1,sc1c'   ,  3000 GeV,   2.2579369E+01, 7.71E-03,   3)
PROCESS(`gl,A', `sc2,sc2c'   ,  3000 GeV,   2.2578124E+01, 8.09E-03,   3)
PROCESS(`gl,A', `st1,st1c'   ,  3000 GeV,   2.6687660E+01, 1.35E-02,   3)
PROCESS(`gl,A', `st2,st2c'   ,  3000 GeV,   2.0965880E+01, 6.24E-03,   3)
PROCESS(`gl,A', `sd1,sd1c'   ,  3000 GeV,   5.6383495E+00, 1.89E-03,   3)
PROCESS(`gl,A', `sd2,sd2c'   ,  3000 GeV,   5.6459148E+00, 1.89E-03,   3)
PROCESS(`gl,A', `ss1,ss1c'   ,  3000 GeV,   5.6394649E+00, 1.93E-03,   3)
PROCESS(`gl,A', `ss2,ss2c'   ,  3000 GeV,   5.6431251E+00, 1.93E-03,   3)
PROCESS(`gl,A', `sb1,sb1c'   ,  3000 GeV,   5.6568356E+00, 1.93E-03,   3)
PROCESS(`gl,A', `sb2,sb2c'   ,  3000 GeV,   5.6280896E+00, 1.89E-03,   3)
PROCESS(`gl,A', `t,T'	     ,  5000 GeV,   2.3006541E+02, 6.72E-03,   3)
dnl Fluctuations in the photon polarization vectors for high energies
PROCESS(`gl,A', `su1,su1c'   ,  5000 GeV,   1.1863544E+01, 5.21E-04,   3)
PROCESS(`gl,A', `su2,su2c'   ,  5000 GeV,   1.1858754E+01, 5.32E-04,   4)
PROCESS(`gl,A', `sc1,sc1c'   ,  5000 GeV,   1.1864777E+01, 5.10E-04,   5)
PROCESS(`gl,A', `sc2,sc2c'   ,  5000 GeV,   1.1860867E+01, 5.18E-04,   3)
PROCESS(`gl,A', `st1,st1c'   ,  5000 GeV,   1.3848759E+01, 6.41E-04,   3)
PROCESS(`gl,A', `st2,st2c'   ,  5000 GeV,   1.0334309E+01, 4.70E-04,   4)
PROCESS(`gl,A', `sd1,sd1c'   ,  5000 GeV,   2.9620062E+00, 1.27E-04,   3)
PROCESS(`gl,A', `sd2,sd2c'   ,  5000 GeV,   2.9639171E+00, 1.29E-04,   4)
PROCESS(`gl,A', `ss1,ss1c'   ,  5000 GeV,   2.9622609E+00, 1.28E-04,   4)
PROCESS(`gl,A', `ss2,ss2c'   ,  5000 GeV,   2.9639384E+00, 1.30E-04,   4)
PROCESS(`gl,A', `sb1,sb1c'   ,  5000 GeV,   2.9723032E+00, 1.30E-04,   3)
PROCESS(`gl,A', `sb2,sb2c'   ,  5000 GeV,   2.9535715E+00, 1.27E-04,   4)
! -----------------------------------------------------------------------
END_TESTSUITE
