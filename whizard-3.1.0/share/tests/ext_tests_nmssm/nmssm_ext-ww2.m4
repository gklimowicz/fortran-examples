dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-ww2.m4 > nmssm_ext-ww2.sin
dnl   whizard -r nmssm_ext-ww2.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-ww2.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_ww2_)
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
iterations = 3:50000
! -------------------------------------------------------------------------
dnl too strong fluctuations for that process for 5 TeV for this number of 
dnl calls and iterations, leaving this out
dnl PROCESS(`Wm,Wp', `t,T'  	,  5000 GeV,  2.2234859E+03,  3.11E+01,  3)
PROCESS(`Wm,Wp', `Z,Z'		,  5000 GeV,  3.0256862E+05,  4.09E+02,  3)
PROCESS(`Wm,Wp', `se1,SE1'	,  5000 GeV,  4.2961815E+00,  2.69E-03,  3)
PROCESS(`Wm,Wp', `se2,SE2'	,  5000 GeV,  9.7748764E-02,  6.04E-05,  3)
dnl Wrong ordering of smuons in the original WHIZARD/FeynRules comparison
PROCESS(`Wm,Wp', `smu1,SMU1'	,  5000 GeV,  4.2245897E+00,  2.69E-03,  3)
PROCESS(`Wm,Wp', `smu2,SMU2'	,  5000 GeV,  9.7880662E-02,  5.98E-05,  3)
PROCESS(`Wm,Wp', `stau1,STAU1'	,  5000 GeV,  1.0729103E+00,  6.48E-04,  3)
PROCESS(`Wm,Wp', `stau2,STAU2'	,  5000 GeV,  1.1320350E+00,  6.86E-04,  3)
PROCESS(`Wm,Wp', `stau1,STAU2'	,  5000 GeV,  1.1026023E+00,  7.03E-04,  3)
PROCESS(`Wm,Wp', `sn11,SN11'	,  5000 GeV,  4.0040052E+00,  3.01E-03,  3)
dnl Wrong ordering of smuons in the original WHIZARD/FeynRules comparison
PROCESS(`Wm,Wp', `sn21,SN21'	,  5000 GeV,  3.9698134E+00,  3.08E-03,  3)
PROCESS(`Wm,Wp', `sn31,SN31'	,  5000 GeV,  3.9973031E+00,  3.07E-03,  3)
PROCESS(`Wm,Wp', `su1,su1c'	,  5000 GeV,  7.5989113E+00,  4.34E-03,  3)
PROCESS(`Wm,Wp', `su2,su2c'	,  5000 GeV,  1.1032293E-01,  6.65E-05,  3)
PROCESS(`Wm,Wp', `sc1,sc1c'	,  5000 GeV,  7.6028783E+00,  4.34E-03,  3)
PROCESS(`Wm,Wp', `sc2,sc2c'	,  5000 GeV,  1.0994796E-01,  6.66E-05,  3)
PROCESS(`Wm,Wp', `st1,st1c'	,  5000 GeV,  6.8687413E+01,  2.31E-02,  3)
PROCESS(`Wm,Wp', `st2,st2c'	,  5000 GeV,  1.3780698E+02,  5.47E-02,  3)
PROCESS(`Wm,Wp', `st1,st2c'	,  5000 GeV,  3.2163755E+01,  1.43E-02,  3)
PROCESS(`Wm,Wp', `sd1,sd1c'	,  5000 GeV,  7.3351390E+00,  2.76E-03,  3)
PROCESS(`Wm,Wp', `sd2,sd2c'	,  5000 GeV,  2.7574850E-02,  1.66E-05,  3)
PROCESS(`Wm,Wp', `ss1,ss1c'	,  5000 GeV,  7.3280538E+00,  2.80E-03,  3)
PROCESS(`Wm,Wp', `ss2,ss2c'	,  5000 GeV,  2.7589163E-02,  1.65E-05,  3)
PROCESS(`Wm,Wp', `sb1,sb1c'	,  5000 GeV,  1.5746729E+01,  1.21E-02,  3)
PROCESS(`Wm,Wp', `sb2,sb2c'	,  5000 GeV,  2.1558563E+01,  1.38E-02,  3)
PROCESS(`Wm,Wp', `sb1,sb2c'	,  5000 GeV,  1.8252757E+01,  1.27E-02,  3)
PROCESS(`Wm,Wp', `neu1,neu1'	,  5000 GeV,  4.0046375E-01,  1.34E-04,  3)
PROCESS(`Wm,Wp', `neu1,neu2'	,  5000 GeV,  4.4419852E+00,  1.34E-03,  3)
PROCESS(`Wm,Wp', `neu1,neu3'	,  5000 GeV,  1.1718684E+01,  4.49E-03,  3)
PROCESS(`Wm,Wp', `neu1,neu4'	,  5000 GeV,  3.6958497E+00,  1.07E-03,  3)
PROCESS(`Wm,Wp', `neu1,neu5'	,  5000 GeV,  4.1553157E+00,  1.59E-03,  3)
PROCESS(`Wm,Wp', `neu2,neu2'	,  5000 GeV,  1.1407746E+01,  6.90E-03,  3)
PROCESS(`Wm,Wp', `neu2,neu3'	,  5000 GeV,  5.2925007E+01,  1.32E-02,  3)
PROCESS(`Wm,Wp', `neu2,neu4'	,  5000 GeV,  2.9786186E+01,  9.96E-03,  3)
PROCESS(`Wm,Wp', `neu2,neu5'	,  5000 GeV,  2.9080333E+01,  9.27E-03,  3)
PROCESS(`Wm,Wp', `neu3,neu3'	,  5000 GeV,  2.5259643E+02,  1.01E-01,  3)
PROCESS(`Wm,Wp', `neu3,neu4'	,  5000 GeV,  9.6820168E+01,  2.20E-02,  3)
PROCESS(`Wm,Wp', `neu3,neu5'	,  5000 GeV,  9.2469130E+01,  2.20E-02,  3)
PROCESS(`Wm,Wp', `neu4,neu4'	,  5000 GeV,  1.8883780E+01,  7.82E-03,  3)
PROCESS(`Wm,Wp', `neu4,neu5'	,  5000 GeV,  2.8424277E+01,  1.27E-02,  3)
PROCESS(`Wm,Wp', `neu5,neu5'	,  5000 GeV,  1.6508763E+01,  6.74E-03,  3)
PROCESS(`Wm,Wp', `"ch1+","ch1-"',  5000 GeV,  2.7365872E+02,  2.38E-01,  3)
PROCESS(`Wm,Wp', `"ch2+","ch2-"',  5000 GeV,  9.0179187E+01,  7.29E-02,  3)
PROCESS(`Wm,Wp', `"ch1+","ch2-"',  5000 GeV,  3.9432532E+01,  3.21E-02,  3)
PROCESS(`Wm,Wp', `h01,h01'	,  5000 GeV,  5.7542757E+03,  5.41E+01,  3)
PROCESS(`Wm,Wp', `h01,h02'	,  5000 GeV,  2.4598879E+02,  2.31E+00,  3)
PROCESS(`Wm,Wp', `h01,h03'	,  5000 GeV,  1.1767950E+00,  5.60E-04,  3)
PROCESS(`Wm,Wp', `h02,h02'	,  5000 GeV,  4.7975263E+00,  2.97E-02,  3)
PROCESS(`Wm,Wp', `h02,h03'	,  5000 GeV,  4.2294600E-01,  1.35E-04,  3)
PROCESS(`Wm,Wp', `h03,h03'	,  5000 GeV,  3.1540532E+00,  8.58E-04,  3)
PROCESS(`Wm,Wp', `A01,A01'	,  5000 GeV,  3.2330117E+00,  1.44E-03,  3)
PROCESS(`Wm,Wp', `A01,A02'	,  5000 GeV,  5.6717388E-01,  2.31E-04,  3)
PROCESS(`Wm,Wp', `A02,A02'	,  5000 GeV,  4.5513645E+00,  1.26E-03,  3)
PROCESS(`Wm,Wp', `Z,h01'	,  5000 GeV,  8.2984240E+04,  1.12E+02,  3)
PROCESS(`Wm,Wp', `Z,h02'	,  5000 GeV,  1.7813470E+03,  2.42E+00,  3)
PROCESS(`Wm,Wp', `Z,h03'	,  5000 GeV,  4.1611376E-02,  5.51E-05,  3)
PROCESS(`Wm,Wp', `A01,h01'	,  5000 GeV,  1.3628187E-03,  6.92E-07,  3)
PROCESS(`Wm,Wp', `A01,h02'	,  5000 GeV,  5.3713419E-02,  2.73E-05,  3)
PROCESS(`Wm,Wp', `A01,h03'	,  5000 GeV,  3.9356374E-01,  1.80E-04,  3)
PROCESS(`Wm,Wp', `A02,h01'	,  5000 GeV,  4.4113338E-03,  2.04E-06,  3)
PROCESS(`Wm,Wp', `A02,h02'	,  5000 GeV,  1.7295811E-01,  8.07E-05,  3)
PROCESS(`Wm,Wp', `A02,h03'	,  5000 GeV,  3.9868079E-01,  2.30E-04,  3)
PROCESS(`Wm,Wp', `Hp,Hm'	,  5000 GeV,  6.8585415E+00,  3.08E-03,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
