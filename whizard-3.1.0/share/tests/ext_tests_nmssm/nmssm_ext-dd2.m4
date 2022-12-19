dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-dd2.m4 > nmssm_ext-dd2.sin
dnl   whizard -r nmssm_ext-dd2.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-dd2.m4 $, nmssm_ext_dd2_)
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
PROCESS(`d,D', `n1,N1'		,  5000 GeV,  4.7784112E-01,  1.85E-04,  3)
PROCESS(`d,D', `n2,N2'		,  5000 GeV,  4.7758343E-01,  1.86E-04,  3)
PROCESS(`d,D', `n3,N3'		,  5000 GeV,  4.7769266E-01,  1.85E-04,  3)
PROCESS(`d,D', `t,T'		,  5000 GeV,  2.0032293E+02,  1.99E-02,  3)
PROCESS(`d,D', `Wp,Wm'		,  5000 GeV,  6.3693894E+01,  8.97E-02,  3)
PROCESS(`d,D', `Z,Z'		,  5000 GeV,  1.3862681E+01,  8.38E-03,  3)
PROCESS(`d,D', `se1,SE1'	,  5000 GeV,  1.5981999E-01,  7.89E-05,  3)
PROCESS(`d,D', `se2,SE2'	,  5000 GeV,  3.8912530E-02,  1.92E-05,  3)
PROCESS(`d,D', `smu1,SMU1'	,  5000 GeV,  1.5984163E-01,  7.87E-05,  3)
PROCESS(`d,D', `smu2,SMU2'	,  5000 GeV,  3.8917434E-02,  1.90E-05,  3)
PROCESS(`d,D', `stau1,STAU1'	,  5000 GeV,  3.9433110E-02,  1.95E-05,  3)
PROCESS(`d,D', `stau2,STAU2'	,  5000 GeV,  4.1099326E-02,  2.02E-05,  3)
PROCESS(`d,D', `stau1,STAU2'	,  5000 GeV,  5.9175389E-02,  2.88E-05,  3)
PROCESS(`d,D', `sn11,SN11'	,  5000 GeV,  2.3666839E-01,  1.16E-04,  3)
PROCESS(`d,D', `sn21,SN21'	,  5000 GeV,  2.3650576E-01,  1.17E-04,  3)
PROCESS(`d,D', `sn31,SN31'	,  5000 GeV,  2.3667068E-01,  1.16E-04,  3)
PROCESS(`d,D', `su1,su1c'	,  5000 GeV,  5.4188804E+01,  4.00E-02,  3)
PROCESS(`d,D', `su2,su2c'	,  5000 GeV,  3.7995133E+01,  3.50E-02,  3)
PROCESS(`d,D', `sc1,sc1c'	,  5000 GeV,  3.8389675E+01,  2.07E-02,  3)
PROCESS(`d,D', `sc2,sc2c'	,  5000 GeV,  3.8041949E+01,  2.05E-02,  3)
PROCESS(`d,D', `st1,st1c'	,  5000 GeV,  4.2508975E+01,  2.28E-02,  3)
PROCESS(`d,D', `st2,st2c'	,  5000 GeV,  3.3363312E+01,  1.78E-02,  3)
PROCESS(`d,D', `st1,st2c'	,  5000 GeV,  1.3602127E-01,  7.38E-05,  3)
PROCESS(`d,D', `sd1,sd1c'	,  5000 GeV,  3.5251577E+02,  2.42E-01,  3)
PROCESS(`d,D', `sd2,sd2c'	,  5000 GeV,  3.6971079E+02,  2.60E-01,  3)
PROCESS(`d,D', `ss1,ss1c'	,  5000 GeV,  3.8420980E+01,  2.07E-02,  3)
PROCESS(`d,D', `ss2,ss2c'	,  5000 GeV,  3.7988366E+01,  2.04E-02,  3)
PROCESS(`d,D', `sb1,sb1c'	,  5000 GeV,  3.8112499E+01,  2.07E-02,  3)
PROCESS(`d,D', `sb2,sb2c'	,  5000 GeV,  3.8003883E+01,  2.05E-02,  3)
PROCESS(`d,D', `sb1,sb2c'	,  5000 GeV,  1.3541762E-01,  7.38E-05,  3)
PROCESS(`d,D', `neu1,neu1'	,  5000 GeV,  2.0813272E-02,  9.72E-06,  3)
PROCESS(`d,D', `neu1,neu2'	,  5000 GeV,  1.0399708E-02,  4.47E-06,  3)
PROCESS(`d,D', `neu1,neu3'	,  5000 GeV,  1.2552355E-01,  6.00E-05,  3)
PROCESS(`d,D', `neu1,neu4'	,  5000 GeV,  4.5749537E-03,  2.34E-06,  3)
PROCESS(`d,D', `neu1,neu5'	,  5000 GeV,  1.6089088E-02,  4.72E-06,  3)
PROCESS(`d,D', `neu2,neu2'	,  5000 GeV,  3.4847042E-03,  1.34E-06,  3)
PROCESS(`d,D', `neu2,neu3'	,  5000 GeV,  6.7061824E-02,  3.08E-05,  3)
PROCESS(`d,D', `neu2,neu4'	,  5000 GeV,  1.2740050E-02,  2.96E-06,  3)
PROCESS(`d,D', `neu2,neu5'	,  5000 GeV,  1.1002325E-02,  3.25E-06,  3)
PROCESS(`d,D', `neu3,neu3'	,  5000 GeV,  6.4192058E-01,  3.08E-04,  3)
PROCESS(`d,D', `neu3,neu4'	,  5000 GeV,  9.7546660E-02,  4.11E-05,  3)
PROCESS(`d,D', `neu3,neu5'	,  5000 GeV,  3.4778709E-02,  1.07E-05,  3)
PROCESS(`d,D', `neu4,neu4'	,  5000 GeV,  9.8420993E-04,  5.65E-07,  3)
PROCESS(`d,D', `neu4,neu5'	,  5000 GeV,  8.6151187E-01,  2.42E-04,  3)
PROCESS(`d,D', `neu5,neu5'	,  5000 GeV,  8.1787827E-06,  1.75E-09,  3)
PROCESS(`d,D', `"ch1+","ch1-"'	,  5000 GeV,  2.1980122E+00,  1.40E-03,  3)
PROCESS(`d,D', `"ch2+","ch2-"'	,  5000 GeV,  6.5329750E-01,  1.47E-04,  3)
PROCESS(`d,D', `"ch1+","ch2-"'	,  5000 GeV,  5.7847960E-02,  2.51E-05,  3)
PROCESS(`d,D', `sgl,sgl'	,  5000 GeV,  3.4525322E+02,  1.76E-01,  3)
PROCESS(`d,D', `Z,h01'		,  5000 GeV,  2.3396073E-01,  1.14E-04,  3)
PROCESS(`d,D', `Z,h02'		,  5000 GeV,  5.0127366E-03,  2.46E-06,  3)
PROCESS(`d,D', `Z,h03'		,  5000 GeV,  1.0204811E-07,  4.99E-11,  3)
PROCESS(`d,D', `A01,h01'	,  5000 GeV,  1.3656842E-07,  6.71E-11,  3)
PROCESS(`d,D', `A01,h02'	,  5000 GeV,  5.3767519E-06,  2.64E-09,  3)
PROCESS(`d,D', `A01,h03'	,  5000 GeV,  1.3573170E-03,  6.67E-07,  3)
PROCESS(`d,D', `A02,h01'	,  5000 GeV,  1.5569966E-05,  7.67E-09,  3)
PROCESS(`d,D', `A02,h02'	,  5000 GeV,  6.1291346E-04,  3.01E-07,  3)
PROCESS(`d,D', `A02,h03'	,  5000 GeV,  1.4677142E-01,  7.18E-05,  3)
PROCESS(`d,D', `Hp,Hm'		,  5000 GeV,  1.0027830E-01,  4.91E-05,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
