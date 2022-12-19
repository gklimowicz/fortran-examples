dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-dd1.m4 > nmssm_ext-dd1.sin
dnl   whizard -r nmssm_ext-dd1
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-dd1.m4 $, nmssm_ext_dd1_)
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
PROCESS(`d,D', `n1,N1'		,  3000 GeV,  1.3300240E+00,  8.71E-04,  3)
PROCESS(`d,D', `n2,N2'		,  3000 GeV,  1.3287159E+00,  8.95E-04,  3)
PROCESS(`d,D', `n3,N3'		,  3000 GeV,  1.3285066E+00,  8.96E-04,  3)
PROCESS(`d,D', `t,T'		,  3000 GeV,  5.5660195E+02,  9.05E-02,  3)
PROCESS(`d,D', `Wp,Wm'		,  3000 GeV,  1.5065884E+02,  3.09E-01,  3)
PROCESS(`d,D', `Z,Z'		,  3000 GeV,  3.2931023E+01,  1.40E-02,  3)
PROCESS(`d,D', `se1,SE1'	,  3000 GeV,  4.3665895E-01,  3.70E-04,  3)
PROCESS(`d,D', `se2,SE2'	,  3000 GeV,  1.0630514E-01,  9.02E-05,  3)
PROCESS(`d,D', `smu1,SMU1'	,  3000 GeV,  4.3644812E-01,  3.70E-04,  3)
PROCESS(`d,D', `smu2,SMU2'	,  3000 GeV,  1.0630732E-01,  9.10E-05,  3)
PROCESS(`d,D', `stau1,STAU1'	,  3000 GeV,  1.0800505E-01,  9.18E-05,  3)
PROCESS(`d,D', `stau2,STAU2'	,  3000 GeV,  1.1168450E-01,  9.67E-05,  3)
PROCESS(`d,D', `stau1,STAU2'	,  3000 GeV,  1.6136740E-01,  1.38E-04,  3)
PROCESS(`d,D', `sn11,SN11'	,  3000 GeV,  6.4651612E-01,  5.64E-04,  3)
PROCESS(`d,D', `sn21,SN21'	,  3000 GeV,  6.4745457E-01,  5.53E-04,  3)
PROCESS(`d,D', `sn31,SN31'	,  3000 GeV,  6.4829456E-01,  5.45E-04,  3)
PROCESS(`d,D', `su1,su1c'	,  3000 GeV,  5.5503948E+01,  6.68E-02,  3)
PROCESS(`d,D', `su2,su2c'	,  3000 GeV,  5.5268553E+01,  7.66E-02,  3)
PROCESS(`d,D', `sc1,sc1c'	,  3000 GeV,  5.5676113E+01,  4.85E-02,  3)
PROCESS(`d,D', `sc2,sc2c'	,  3000 GeV,  5.5134978E+01,  4.73E-02,  3)
PROCESS(`d,D', `st1,st1c'	,  3000 GeV,  8.4218398E+01,  7.19E-02,  3)
PROCESS(`d,D', `st2,st2c'	,  3000 GeV,  2.8178709E+01,  2.40E-02,  3)
PROCESS(`d,D', `st1,st2c'	,  3000 GeV,  1.9813135E-01,  1.66E-04,  3)
PROCESS(`d,D', `sd1,sd1c'	,  3000 GeV,  3.8551545E+02,  4.17E-01,  3)
PROCESS(`d,D', `sd2,sd2c'	,  3000 GeV,  4.1841590E+02,  4.53E-01,  3)
PROCESS(`d,D', `ss1,ss1c'	,  3000 GeV,  5.5616447E+01,  4.75E-02,  3)
PROCESS(`d,D', `ss2,ss2c'	,  3000 GeV,  5.5102350E+01,  4.66E-02,  3)
PROCESS(`d,D', `sb1,sb1c'	,  3000 GeV,  5.5758175E+01,  4.69E-02,  3)
PROCESS(`d,D', `sb2,sb2c'	,  3000 GeV,  5.4491858E+01,  4.73E-02,  3)
PROCESS(`d,D', `sb1,sb2c'	,  3000 GeV,  1.9624038E-01,  1.68E-04,  3)
PROCESS(`d,D', `neu1,neu1'	,  3000 GeV,  4.3557183E-02,  2.37E-05,  3)
PROCESS(`d,D', `neu1,neu2'	,  3000 GeV,  2.2946499E-02,  1.11E-05,  3)
PROCESS(`d,D', `neu1,neu3'	,  3000 GeV,  2.5937369E-01,  1.48E-04,  3)
PROCESS(`d,D', `neu1,neu4'	,  3000 GeV,  9.0690319E-03,  6.05E-06,  3)
PROCESS(`d,D', `neu1,neu5'	,  3000 GeV,  4.2697498E-02,  1.15E-05,  3)
PROCESS(`d,D', `neu2,neu2'	,  3000 GeV,  8.2451345E-03,  3.33E-06,  3)
PROCESS(`d,D', `neu2,neu3'	,  3000 GeV,  1.4244790E-01,  7.47E-05,  3)
PROCESS(`d,D', `neu2,neu4'	,  3000 GeV,  3.5662005E-02,  7.32E-06,  3)
PROCESS(`d,D', `neu2,neu5'	,  3000 GeV,  2.9431217E-02,  8.14E-06,  3)
PROCESS(`d,D', `neu3,neu3'	,  3000 GeV,  1.3246643E+00,  7.55E-04,  3)
PROCESS(`d,D', `neu3,neu4'	,  3000 GeV,  2.1462193E-01,  1.06E-04,  3)
PROCESS(`d,D', `neu3,neu5'	,  3000 GeV,  9.0713931E-02,  2.73E-05,  3)
PROCESS(`d,D', `neu4,neu4'	,  3000 GeV,  1.7783210E-03,  1.50E-06,  3)
PROCESS(`d,D', `neu4,neu5'	,  3000 GeV,  2.3826792E+00,  6.13E-04,  3)
PROCESS(`d,D', `neu5,neu5'	,  3000 GeV,  2.1030030E-05,  3.72E-09,  3)
PROCESS(`d,D', `"ch1+","ch1-"'  ,  3000 GeV,  4.9522093E+00,  5.02E-03,  3)
PROCESS(`d,D', `"ch2+","ch2-"'	,  3000 GeV,  1.8083640E+00,  3.59E-04,  3)
PROCESS(`d,D', `"ch1+","ch2-"'	,  3000 GeV,  1.2000862E-01,  6.63E-05,  3)
PROCESS(`d,D', `sgl,sgl'	,  3000 GeV,  6.4966039E+02,  3.52E-01,  3)
PROCESS(`d,D', `Z,h01'		,  3000 GeV,  6.5162015E-01,  5.44E-04,  3)
PROCESS(`d,D', `Z,h02'		,  3000 GeV,  1.3952756E-02,  1.17E-05,  3)
PROCESS(`d,D', `Z,h03'		,  3000 GeV,  1.8883861E-07,  1.55E-10,  3)
PROCESS(`d,D', `A01,h01'	,  3000 GeV,  3.7778612E-07,  3.20E-10,  3)
PROCESS(`d,D', `A01,h02'	,  3000 GeV,  1.4881599E-05,  1.25E-08,  3)
PROCESS(`d,D', `A01,h03'	,  3000 GeV,  2.4865178E-03,  2.08E-06,  3)
PROCESS(`d,D', `A02,h01'	,  3000 GeV,  2.8250123E-05,  2.40E-08,  3)
PROCESS(`d,D', `A02,h02'	,  3000 GeV,  1.1097687E-03,  9.50E-07,  3)
PROCESS(`d,D', `A02,h03'	,  3000 GeV,  7.9049726E-02,  6.62E-05,  3)
PROCESS(`d,D', `Hp,Hm'		,  3000 GeV,  5.4147842E-02,  4.60E-05,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
