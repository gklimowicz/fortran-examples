dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-en.m4 > nmssm_ext-en.sin
dnl   whizard -r nmssm_ext-en.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-en.m4 $, nmssm_ext_en_)
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
PROCESS(`e1,N1', `Z,Wm'		,  3000 GeV,  7.1180098E+02,  1.28E+00,  3)
PROCESS(`e1,N1', `se1,SN11'	,  3000 GeV,  7.1516406E+01,  9.76E-02,  3)
PROCESS(`e1,N1', `se2,SN11'	,  3000 GeV,  6.6334194E+00,  9.64E-03,  3)
PROCESS(`e1,N1', `smu1,SN21'	,  3000 GeV,  1.2385832E+01,  1.38E-02,  3)
PROCESS(`e1,N1', `stau1,SN31'	,  3000 GeV,  6.1062640E+00,  5.22E-03,  3)
PROCESS(`e1,N1', `stau2,SN31'	,  3000 GeV,  6.2641030E+00,  5.28E-03,  3)
PROCESS(`e1,N1', `neu1,"ch1-"'	,  3000 GeV,  2.0380748E+01,  2.40E-02,  3)
PROCESS(`e1,N1', `neu2,"ch1-"'	,  3000 GeV,  2.7937816E+00,  4.08E-03,  3)
PROCESS(`e1,N1', `neu3,"ch1-"'	,  3000 GeV,  2.0499984E+01,  3.02E-02,  3)
PROCESS(`e1,N1', `neu4,"ch1-"'	,  3000 GeV,  1.1242076E+00,  9.14E-04,  3)
PROCESS(`e1,N1', `neu5,"ch1-"'	,  3000 GeV,  6.9004088E-01,  6.54E-04,  3)
PROCESS(`e1,N1', `neu1,"ch2-"'	,  3000 GeV,  9.9103088E-01,  1.64E-03,  3)
PROCESS(`e1,N1', `neu2,"ch2-"'	,  3000 GeV,  6.9465414E-01,  6.98E-04,  3)
PROCESS(`e1,N1', `neu3,"ch2-"'	,  3000 GeV,  1.0220287E+00,  1.64E-03,  3)
PROCESS(`e1,N1', `neu4,"ch2-"'	,  3000 GeV,  2.3385982E+01,  2.32E-02,  3)
PROCESS(`e1,N1', `neu5,"ch2-"'	,  3000 GeV,  2.3691116E+01,  2.16E-02,  3)
PROCESS(`e1,N1', `Z,Wm'		,  5000 GeV,  3.0162326E+02,  6.36E-01,  3)
PROCESS(`e1,N1', `se1,SN11'	,  5000 GeV,  2.9475778E+01,  3.56E-02,  3)
PROCESS(`e1,N1', `se2,SN11'	,  5000 GeV,  2.4708784E+00,  3.74E-03,  3)
PROCESS(`e1,N1', `smu1,SN21'	,  5000 GeV,  4.5180620E+00,  3.92E-03,  3)
PROCESS(`e1,N1', `stau1,SN31'	,  5000 GeV,  2.2310762E+00,  1.90E-03,  3)
PROCESS(`e1,N1', `stau2,SN31'	,  5000 GeV,  2.2900994E+00,  1.98E-03,  3)
PROCESS(`e1,N1', `neu1,"ch1-"'	,  5000 GeV,  7.5554683E+00,  9.36E-03,  3)
PROCESS(`e1,N1', `neu2,"ch1-"'	,  5000 GeV,  1.0443521E+00,  1.60E-03,  3)
PROCESS(`e1,N1', `neu3,"ch1-"'	,  5000 GeV,  7.8134296E+00,  1.20E-02,  3)
PROCESS(`e1,N1', `neu4,"ch1-"'	,  5000 GeV,  3.5623320E-01,  3.18E-04,  3)
PROCESS(`e1,N1', `neu5,"ch1-"'	,  5000 GeV,  2.5540244E-01,  2.74E-04,  3)
PROCESS(`e1,N1', `neu1,"ch2-"'	,  5000 GeV,  3.5589178E-01,  4.24E-04,  3)
PROCESS(`e1,N1', `neu2,"ch2-"'	,  5000 GeV,  2.5129162E-01,  2.78E-04,  3)
PROCESS(`e1,N1', `neu3,"ch2-"'	,  5000 GeV,  3.3349600E-01,  3.00E-04,  3)
PROCESS(`e1,N1', `neu4,"ch2-"'	,  5000 GeV,  8.6121750E+00,  9.70E-03,  3)
PROCESS(`e1,N1', `neu5,"ch2-"'	,  5000 GeV,  8.5984160E+00,  9.08E-03,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
