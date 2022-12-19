dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-zz1.m4 > nmssm_ext-zz1.sin
dnl   whizard -r nmssm_ext-zz1.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-zz1.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_zz1_)
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
PROCESS(`Z,Z', `t,T'		,  3000 GeV,  1.3967134E+03,  1.46E+00,  3)
PROCESS(`Z,Z', `Wp,Wm'		,  3000 GeV,  6.0765347E+05,  1.35E+03,  3)
PROCESS(`Z,Z', `Z,Z'		,  3000 GeV,  8.4984692E+01,  9.92E-02,  3)
PROCESS(`Z,Z', `se1,SE1'	,  3000 GeV,  1.2485340E+00,  1.41E-03,  3)
PROCESS(`Z,Z', `se2,SE2'	,  3000 GeV,  7.7167582E-01,  8.68E-04,  3)
PROCESS(`Z,Z', `smu1,SMU1'	,  3000 GeV,  1.2465830E+00,  1.42E-03,  3)
PROCESS(`Z,Z', `smu2,SMU2'	,  3000 GeV,  7.7408572E-01,  8.51E-04,  4)
PROCESS(`Z,Z', `stau1,STAU1'	,  3000 GeV,  1.0003363E+00,  1.14E-03,  3)
PROCESS(`Z,Z', `stau2,STAU2'	,  3000 GeV,  9.9911458E-01,  1.10E-03,  3)
PROCESS(`Z,Z', `stau1,STAU2'	,  3000 GeV,  1.4934907E-02,  1.67E-05,  3)
PROCESS(`Z,Z', `sn11,SN11'	,  3000 GeV,  1.3975055E+01,  2.07E-02,  3)
PROCESS(`Z,Z', `sn21,SN21'	,  3000 GeV,  1.4030331E+01,  2.08E-02,  3)
PROCESS(`Z,Z', `sn31,SN31'	,  3000 GeV,  1.3961302E+01,  2.10E-02,  3)
PROCESS(`Z,Z', `su1,su1c'	,  3000 GeV,  4.6361968E+00,  1.02E-03,  3)
PROCESS(`Z,Z', `su2,su2c'	,  3000 GeV,  3.0433837E-01,  8.05E-05,  3)
PROCESS(`Z,Z', `sc1,sc1c'	,  3000 GeV,  4.6341941E+00,  9.69E-04,  3)
PROCESS(`Z,Z', `sc2,sc2c'	,  3000 GeV,  3.0390402E-01,  8.31E-05,  3)
PROCESS(`Z,Z', `st1,st1c'	,  3000 GeV,  4.7421516E+02,  8.10E-02,  3)
PROCESS(`Z,Z', `st2,st2c'	,  3000 GeV,  4.3095314E+01,  4.58E-02,  3)
PROCESS(`Z,Z', `st1,st2c'	,  3000 GeV,  6.1038412E+00,  4.78E-03,  3)
PROCESS(`Z,Z', `sd1,sd1c'	,  3000 GeV,  1.0001657E+01,  1.95E-03,  3)
PROCESS(`Z,Z', `sd2,sd2c'	,  3000 GeV,  4.5302150E-02,  1.57E-05,  3)
PROCESS(`Z,Z', `ss1,ss1c'	,  3000 GeV,  1.0002384E+01,  2.05E-03,  3)
PROCESS(`Z,Z', `ss2,ss2c'	,  3000 GeV,  4.5282252E-02,  1.61E-05,  3)
PROCESS(`Z,Z', `sb1,sb1c'	,  3000 GeV,  2.7847591E+00,  6.54E-04,  3)
PROCESS(`Z,Z', `sb2,sb2c'	,  3000 GeV,  7.1965945E+00,  2.16E-03,  3)
PROCESS(`Z,Z', `sb1,sb2c'	,  3000 GeV,  2.1561519E+00,  4.51E-04,  3)
PROCESS(`Z,Z', `neu1,neu1'	,  3000 GeV,  1.4992185E+00,  1.61E-03,  3)
PROCESS(`Z,Z', `neu1,neu2'	,  3000 GeV,  1.8048862E+00,  2.24E-03,  3)
PROCESS(`Z,Z', `neu1,neu3'	,  3000 GeV,  3.3510535E+00,  3.59E-03,  3)
PROCESS(`Z,Z', `neu1,neu4'	,  3000 GeV,  1.3589698E+01,  1.30E-02,  3)
PROCESS(`Z,Z', `neu1,neu5'	,  3000 GeV,  2.5461034E-01,  2.32E-04,  3)
PROCESS(`Z,Z', `neu2,neu2'	,  3000 GeV,  1.2924366E+01,  1.74E-02,  3)
PROCESS(`Z,Z', `neu2,neu3'	,  3000 GeV,  9.5629250E+00,  1.02E-02,  3)
PROCESS(`Z,Z', `neu2,neu4'	,  3000 GeV,  8.3508117E+00,  7.67E-03,  3)
PROCESS(`Z,Z', `neu2,neu5'	,  3000 GeV,  4.5508369E+01,  4.33E-02,  3)
PROCESS(`Z,Z', `neu3,neu3'	,  3000 GeV,  2.4555420E+00,  3.88E-03,  3)
PROCESS(`Z,Z', `neu3,neu4'	,  3000 GeV,  1.5611183E+01,  1.44E-02,  3)
PROCESS(`Z,Z', `neu3,neu5'	,  3000 GeV,  1.2425948E+01,  1.18E-02,  3)
PROCESS(`Z,Z', `neu4,neu4'	,  3000 GeV,  7.8120553E+01,  6.12E-02,  3)
PROCESS(`Z,Z', `neu4,neu5'	,  3000 GeV,  1.5341093E+00,  1.92E-03,  3)
PROCESS(`Z,Z', `neu5,neu5'	,  3000 GeV,  5.5808973E+01,  5.28E-02,  3)
PROCESS(`Z,Z', `"ch1+","ch1-"'	,  3000 GeV,  7.9878873E+02,  4.52E-01,  3)
PROCESS(`Z,Z', `"ch1+","ch2-"'	,  3000 GeV,  1.6371703E+02,  5.88E-02,  3)
PROCESS(`Z,Z', `"ch2+","ch2-"'	,  3000 GeV,  3.5201258E+01,  2.70E-02,  3)
PROCESS(`Z,Z', `h01,h01'	,  3000 GeV,  7.4261286E+03,  1.60E+01,  3)
PROCESS(`Z,Z', `h01,h02'	,  3000 GeV,  3.2236214E+02,  6.81E-01,  3)
PROCESS(`Z,Z', `h01,h03'	,  3000 GeV,  2.7714982E+00,  3.39E-03,  3)
PROCESS(`Z,Z', `h02,h02'	,  3000 GeV,  8.1417077E+00,  1.17E-02,  3)
PROCESS(`Z,Z', `h02,h03'	,  3000 GeV,  1.5115585E+00,  2.00E-03,  3)
PROCESS(`Z,Z', `h03,h03'	,  3000 GeV,  3.2205750E+00,  3.70E-03,  3)
PROCESS(`Z,Z', `A01,A01'	,  3000 GeV,  6.5038707E+00,  8.38E-03,  3)
PROCESS(`Z,Z', `A01,A02'	,  3000 GeV,  2.0536918E+00,  2.68E-03,  3)
PROCESS(`Z,Z', `A02,A02'	,  3000 GeV,  2.9014956E+00,  3.19E-03,  3)
PROCESS(`Z,Z', `Z,A01'		,  3000 GeV,  1.8464988E-03,  2.28E-06,  3)
PROCESS(`Z,Z', `Z,A02'	   	,  3000 GeV,  2.9668334E-01,  3.34E-04,  3)
PROCESS(`Z,Z', `Hp,Hm'	   	,  3000 GeV,  4.1675255E-01,  7.64E-05,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
