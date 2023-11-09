dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-tt2.m4 > nmssm_ext-tt2.sin
dnl   whizard -r nmssm_ext-tt2.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-tt1.m4 $, nmssm_ext_tt2_)
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
iterations = 5:100000
! -------------------------------------------------------------------------
PROCESS(`e3,E3', `n1,N1'	,  5000 GeV,  9.7637058E-01,  9.17E-05,  3)
PROCESS(`e3,E3', `n2,N2'	,  5000 GeV,  9.7631483E-01,  9.33E-05,  3)
PROCESS(`e3,E3', `n3,N3'	,  5000 GeV,  5.3558298E+04,  8.51E+01,  3)
PROCESS(`e3,E3', `t,T'		,  5000 GeV,  7.1134414E+00,  1.59E-03,  3)
PROCESS(`e3,E3', `Wp,Wm'	,  5000 GeV,  1.9178002E+02,  2.69E-01,  3)
PROCESS(`e3,E3', `Z,Z'		,  5000 GeV,  1.0334306E+01,  9.01E-03,  3)
PROCESS(`e3,E3', `se1,SE1'	,  5000 GeV,  1.1718506E+00,  5.76E-04,  3)
PROCESS(`e3,E3', `se2,SE2'	,  5000 GeV,  1.0505128E+00,  5.15E-04,  3)
PROCESS(`e3,E3', `smu1,SMU1'	,  5000 GeV,  1.1714728E+00,  5.79E-04,  3)
PROCESS(`e3,E3', `smu2,SMU2'	,  5000 GeV,  1.0498219E+00,  5.19E-04,  3)
PROCESS(`e3,E3', `stau1,STAU1'	,  5000 GeV,  1.4191398E+01,  1.41E-02,  3)
PROCESS(`e3,E3', `stau2,STAU2'	,  5000 GeV,  1.3934237E+01,  1.37E-02,  3)
PROCESS(`e3,E3', `stau1,STAU2'	,  5000 GeV,  1.5722796E+01,  1.38E-02,  3)
PROCESS(`e3,E3', `sn11,SN11'	,  5000 GeV,  4.8398879E-01,  2.38E-04,  3)
PROCESS(`e3,E3', `sn21,SN21'	,  5000 GeV,  4.8366323E-01,  2.39E-04,  3)
PROCESS(`e3,E3', `sn31,SN31'	,  5000 GeV,  7.6642306E+01,  5.85E-02,  3)
PROCESS(`e3,E3', `su1,su1c'	,  5000 GeV,  1.6398867E+00,  8.07E-04,  3)
PROCESS(`e3,E3', `su2,su2c'	,  5000 GeV,  1.0786987E+00,  5.33E-04,  3)
PROCESS(`e3,E3', `sc1,sc1c'	,  5000 GeV,  1.6400928E+00,  8.05E-04,  3)
PROCESS(`e3,E3', `sc2,sc2c'	,  5000 GeV,  1.0795631E+00,  5.30E-04,  3)
PROCESS(`e3,E3', `st1,st1c'	,  5000 GeV,  1.2061169E+00,  5.92E-04,  3)
PROCESS(`e3,E3', `st2,st2c'	,  5000 GeV,  9.4671576E-01,  4.65E-04,  3)
PROCESS(`e3,E3', `st1,st2c'	,  5000 GeV,  2.7852701E-01,  1.37E-04,  3)
PROCESS(`e3,E3', `sd1,sd1c'	,  5000 GeV,  1.1073849E+00,  5.46E-04,  3)
PROCESS(`e3,E3', `sd2,sd2c'	,  5000 GeV,  2.6978007E-01,  1.32E-04,  3)
PROCESS(`e3,E3', `ss1,ss1c'	,  5000 GeV,  1.1081091E+00,  5.42E-04,  3)
PROCESS(`e3,E3', `ss2,ss2c'	,  5000 GeV,  2.6988195E-01,  1.32E-04,  3)
PROCESS(`e3,E3', `sb1,sb1c'	,  5000 GeV,  3.7400016E-01,  1.82E-04,  3)
PROCESS(`e3,E3', `sb2,sb2c'	,  5000 GeV,  4.5060324E-01,  2.21E-04,  3)
PROCESS(`e3,E3', `sb1,sb2c'	,  5000 GeV,  2.7686787E-01,  1.35E-04,  3)
PROCESS(`e3,E3', `neu1,neu1'	,  5000 GeV,  4.4616682E+00,  3.21E-03,  3)
PROCESS(`e3,E3', `neu1,neu2'	,  5000 GeV,  6.2984215E-01,  4.62E-04,  3)
PROCESS(`e3,E3', `neu1,neu3'	,  5000 GeV,  1.2191838E+00,  8.81E-04,  3)
PROCESS(`e3,E3', `neu1,neu4'	,  5000 GeV,  8.6703517E-02,  5.90E-05,  3)
PROCESS(`e3,E3', `neu1,neu5'	,  5000 GeV,  3.8920326E-02,  2.43E-05,  3)
PROCESS(`e3,E3', `neu2,neu2'	,  5000 GeV,  9.6075630E-02,  6.95E-05,  3)
PROCESS(`e3,E3', `neu2,neu3'	,  5000 GeV,  9.2296444E-01,  6.67E-04,  3)
PROCESS(`e3,E3', `neu2,neu4'	,  5000 GeV,  6.4886500E-02,  4.11E-05,  3)
PROCESS(`e3,E3', `neu2,neu5'	,  5000 GeV,  3.7132733E-02,  2.43E-05,  3)
PROCESS(`e3,E3', `neu3,neu3'	,  5000 GeV,  3.1159502E+00,  2.28E-03,  3)
PROCESS(`e3,E3', `neu3,neu4'	,  5000 GeV,  1.8261317E-01,  1.28E-04,  3)
PROCESS(`e3,E3', `neu3,neu5'	,  5000 GeV,  6.1401855E-02,  3.78E-05,  3)
PROCESS(`e3,E3', `neu4,neu4'	,  5000 GeV,  1.5868859E-03,  7.48E-07,  3)
PROCESS(`e3,E3', `neu4,neu5'	,  5000 GeV,  1.7621701E+00,  1.12E-03,  3)
PROCESS(`e3,E3', `neu5,neu5'	,  5000 GeV,  1.3020049E-03,  9.32E-07,  3)
PROCESS(`e3,E3', `"ch1+","ch1-"',  5000 GeV,  8.4877855E+00,  6.57E-03,  3)
PROCESS(`e3,E3', `"ch2+","ch2-"',  5000 GeV,  4.5935212E+00,  2.46E-03,  3)
PROCESS(`e3,E3', `"ch1+","ch2-"',  5000 GeV,  2.6416303E-01,  1.45E-04,  3)
PROCESS(`e3,E3', `"ch2+","ch1-"',  5000 GeV,  2.6443569E-01,  1.43E-04,  3)
# This is a highly fluctuating process, therefore we take a rather high
# error. 
PROCESS(`e3,E3', `h01,h01'	,  5000 GeV,  4.9290295E-06,  2.00E-07,  3)
PROCESS(`e3,E3', `h01,h02'	,  5000 GeV,  2.1625199E-04,  1.96E-07,  3)
PROCESS(`e3,E3', `h01,h03'	,  5000 GeV,  7.8024065E-05,  5.43E-08,  3)
PROCESS(`e3,E3', `h02,h02'	,  5000 GeV,  5.1112918E-05,  4.66E-08,  3)
PROCESS(`e3,E3', `h02,h03'	,  5000 GeV,  9.8440561E-04,  8.90E-07,  3)
PROCESS(`e3,E3', `h03,h03'	,  5000 GeV,  2.2536319E-05,  1.90E-08,  3)
PROCESS(`e3,E3', `A01,A01'	,  5000 GeV,  2.3956259E-07,  2.03E-10,  3)
PROCESS(`e3,E3', `A01,A02'	,  5000 GeV,  1.0783899E-04,  9.59E-08,  3)
PROCESS(`e3,E3', `A02,A02'	,  5000 GeV,  1.4882435E-05,  1.40E-08,  3)
PROCESS(`e3,E3', `Z,h01'	,  5000 GeV,  4.8921997E-01,  5.28E-04,  3)
PROCESS(`e3,E3', `Z,h02'	,  5000 GeV,  1.0467581E-02,  1.17E-05,  3)
PROCESS(`e3,E3', `Z,h03'	,  5000 GeV,  3.6806884E-02,  3.69E-05,  3)
PROCESS(`e3,E3', `Z,A01'	,  5000 GeV,  9.0419491E-04,  5.38E-07,  3)
PROCESS(`e3,E3', `Z,A02'	,  5000 GeV,  3.6660921E-02,  3.68E-05,  3)
PROCESS(`e3,E3', `A01,h01'	,  5000 GeV,  4.8084173E-04,  4.30E-07,  3)
PROCESS(`e3,E3', `A01,h02'	,  5000 GeV,  6.0110499E-05,  5.61E-08,  3)
PROCESS(`e3,E3', `A01,h03'	,  5000 GeV,  2.8205463E-03,  3.13E-06,  3)
PROCESS(`e3,E3', `A02,h01'	,  5000 GeV,  1.0814750E-04,  8.61E-08,  3)
PROCESS(`e3,E3', `A02,h02'	,  5000 GeV,  2.2514974E-03,  2.24E-06,  3)
PROCESS(`e3,E3', `A02,h03'	,  5000 GeV,  2.9419677E-01,  3.26E-04,  3)
PROCESS(`e3,E3', `Hp,Hm'	,  5000 GeV,  7.3254962E-01,  6.95E-04,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
