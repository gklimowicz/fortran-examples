dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-wz.m4 > nmssm_ext-wz.sin
dnl   whizard -r nmssm_ext-wz.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-wz.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_wz_)
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
PROCESS(`Wm,Z', `se1,SN11'	,  3000 GeV,  6.1823728E+00,  6.36E-03,  3)
dnl Wrong ordering of smuons in the original WHIZARD/FeynRules comparison
PROCESS(`Wm,Z', `smu1,SN21'	,  3000 GeV,  6.0174555E+00,  6.36E-03,  3)
PROCESS(`Wm,Z', `stau1,SN31'	,  3000 GeV,  3.0583758E+00,  3.16E-03,  3)
PROCESS(`Wm,Z', `stau2,SN31'	,  3000 GeV,  3.1677196E+00,  3.14E-03,  3)
PROCESS(`Wm,Z', `sd1,su1c'	,  3000 GeV,  5.8539943E+00,  3.56E-03,  3)
PROCESS(`Wm,Z', `ss1,sc1c'	,  3000 GeV,  5.8484252E+00,  3.57E-03,  3)
PROCESS(`Wm,Z', `sb1,st1c'	,  3000 GeV,  5.0997897E+01,  3.27E-02,  3)
PROCESS(`Wm,Z', `sb2,st2c'	,  3000 GeV,  1.0032500E+02,  7.33E-02,  3)
PROCESS(`Wm,Z', `sb1,st2c'	,  3000 GeV,  7.8599193E+01,  5.22E-02,  3)
PROCESS(`Wm,Z', `sb2,st1c'	,  3000 GeV,  5.7834983E+01,  3.44E-02,  3)
PROCESS(`Wm,Z', `"ch1-",neu1'	,  3000 GeV,  1.2799495E+01,  1.21E-02,  3)
PROCESS(`Wm,Z', `"ch1-",neu2'	,  3000 GeV,  5.7350699E+01,  4.23E-02,  3)
PROCESS(`Wm,Z', `"ch1-",neu3'	,  3000 GeV,  4.2624862E+02,  5.06E-01,  3)
PROCESS(`Wm,Z', `"ch1-",neu4'	,  3000 GeV,  1.1574215E+02,  7.16E-02,  3)
PROCESS(`Wm,Z', `"ch1-",neu5'	,  3000 GeV,  1.9356826E+02,  8.62E-02,  3)
PROCESS(`Wm,Z', `"ch2-",neu1'	,  3000 GeV,  8.6745758E+00,  5.77E-03,  3)
PROCESS(`Wm,Z', `"ch2-",neu2'	,  3000 GeV,  3.6788337E+01,  2.52E-02,  3)
PROCESS(`Wm,Z', `"ch2-",neu3'	,  3000 GeV,  1.1242986E+02,  8.93E-02,  3)
PROCESS(`Wm,Z', `"ch2-",neu4'	,  3000 GeV,  7.5919004E+01,  6.71E-02,  3)
PROCESS(`Wm,Z', `"ch2-",neu5'	,  3000 GeV,  4.4208533E+01,  3.18E-02,  3)
PROCESS(`Wm,Z', `Hm,h01'	,  3000 GeV,  8.3326664E-03,  9.93E-06,  3)
PROCESS(`Wm,Z', `Hm,h02'	,  3000 GeV,  3.2535350E-01,  3.96E-04,  3)
PROCESS(`Wm,Z', `Hm,h03'	,  3000 GeV,  4.2075824E+00,  3.56E-03,  3)
PROCESS(`Wm,Z', `Hm,A01'	,  3000 GeV,  7.3823137E-01,  8.43E-04,  3)
PROCESS(`Wm,Z', `Hm,A02'	,  3000 GeV,  2.3423849E+00,  1.68E-03,  3)
PROCESS(`Wm,Z', `Wm,h01'	,  3000 GeV,  7.9377656E+04,  1.76E+02,  3)
PROCESS(`Wm,Z', `Wm,h02'	,  3000 GeV,  1.6924572E+03,  3.78E+00,  3)
PROCESS(`Wm,Z', `Wm,h03'	,  3000 GeV,  3.5044992E-02,  7.04E-05,  3)
PROCESS(`Wm,Z', `Wm,A01'	,  3000 GeV,  1.3415206E-04,  2.55E-07,  3)
PROCESS(`Wm,Z', `Wm,A02'	,  3000 GeV,  3.4721284E-02,  4.63E-05,  3)
PROCESS(`Wm,Z', `Hm,Z'		,  3000 GeV,  3.5638173E-02,  4.73E-05,  3)
PROCESS(`Wm,Z', `se1,SN11'	,  5000 GeV,  2.3723458E+00,  2.79E-03,  3)
dnl Wrong ordering of smuons in the original WHIZARD/FeynRules comparison
PROCESS(`Wm,Z', `smu1,SN21'	,  5000 GeV,  2.3079440E+00,  1.77E-03,  3)
PROCESS(`Wm,Z', `stau1,SN31'	,  5000 GeV,  1.1719223E+00,  1.39E-03,  3)
PROCESS(`Wm,Z', `stau2,SN31'	,  5000 GeV,  1.2065173E+00,  1.40E-03,  3)
PROCESS(`Wm,Z', `sd1,su1c'	,  5000 GeV,  3.5116862E+00,  2.53E-03,  3)
PROCESS(`Wm,Z', `ss1,sc1c'	,  5000 GeV,  3.5130672E+00,  2.51E-03,  3)
PROCESS(`Wm,Z', `sb1,st1c'	,  5000 GeV,  1.1442296E+01,  6.51E-03,  3)
PROCESS(`Wm,Z', `sb2,st2c'	,  5000 GeV,  2.3471568E+01,  1.78E-02,  3)
PROCESS(`Wm,Z', `sb1,st2c'	,  5000 GeV,  1.7936197E+01,  1.37E-02,  3)
PROCESS(`Wm,Z', `sb2,st1c'	,  5000 GeV,  1.2779349E+01,  7.51E-03,  3)
PROCESS(`Wm,Z', `"ch1-",neu1'	,  5000 GeV,  5.5856617E+00,  3.51E-03,  3)
PROCESS(`Wm,Z', `"ch1-",neu2'	,  5000 GeV,  2.5887192E+01,  1.20E-02,  3)
PROCESS(`Wm,Z', `"ch1-",neu3'	,  5000 GeV,  1.8741212E+02,  1.32E-01,  3)
PROCESS(`Wm,Z', `"ch1-",neu4'	,  5000 GeV,  5.2929252E+01,  2.04E-02,  3)
PROCESS(`Wm,Z', `"ch1-",neu5'	,  5000 GeV,  8.6822394E+01,  2.60E-02,  3)
PROCESS(`Wm,Z', `"ch2-",neu1'	,  5000 GeV,  4.2531655E+00,  1.72E-03,  3)
PROCESS(`Wm,Z', `"ch2-",neu2'	,  5000 GeV,  1.8621651E+01,  8.07E-03,  3)
PROCESS(`Wm,Z', `"ch2-",neu3'	,  5000 GeV,  5.0990585E+01,  2.50E-02,  3)
PROCESS(`Wm,Z', `"ch2-",neu4'	,  5000 GeV,  3.8780175E+01,  2.23E-02,  3)
PROCESS(`Wm,Z', `"ch2-",neu5'	,  5000 GeV,  2.3030112E+01,  9.19E-03,  3)
PROCESS(`Wm,Z', `Hm,h01'	,  5000 GeV,  3.0693647E-03,  2.05E-06,  3)
PROCESS(`Wm,Z', `Hm,h02'	,  5000 GeV,  1.2069146E-01,  8.00E-05,  3)
PROCESS(`Wm,Z', `Hm,h03'	,  5000 GeV,  2.1415142E+00,  1.16E-03,  3)
PROCESS(`Wm,Z', `Hm,A01'	,  5000 GeV,  2.7179804E-01,  1.73E-04,  3)
PROCESS(`Wm,Z', `Hm,A02'	,  5000 GeV,  1.2822785E+00,  4.44E-04,  3)
PROCESS(`Wm,Z', `Wm,h01'	,  5000 GeV,  7.9407356E+04,  1.08E+02,  3)
PROCESS(`Wm,Z', `Wm,h02'	,  5000 GeV,  1.7051486E+03,  2.28E+00,  3)
PROCESS(`Wm,Z', `Wm,h03'	,  5000 GeV,  3.9902925E-02,  5.23E-05,  3)
PROCESS(`Wm,Z', `Wm,A01'	,  5000 GeV,  6.9977224E-05,  7.93E-08,  3)
PROCESS(`Wm,Z', `Wm,A02'	,  5000 GeV,  1.3608286E-02,  1.25E-05,  3)
PROCESS(`Wm,Z', `Hm,Z'		,  5000 GeV,  1.3989255E-02,  1.25E-05,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
