dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-ee1.m4 > nmssm_ext-ee1.sin
dnl   whizard -r nmssm_ext1-ee.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-ee1.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_ee1_)
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
PROCESS(`e1,E1', `n1,N1'	,  3000 GeV,  5.3521893E+04,  1.45E+02,  3)
PROCESS(`e1,E1', `n2,N2'	,  3000 GeV,  2.7154546E+00,  4.15E-04,  3)
PROCESS(`e1,E1', `n3,N3'	,  3000 GeV,  2.7160902E+00,  4.32E-04,  3)
PROCESS(`e1,E1', `t,T'		,  3000 GeV,  1.9711937E+01,  8.03E-03,  3)
PROCESS(`e1,E1', `Wp,Wm'	,  3000 GeV,  4.5554211E+02,  9.22E-01,  3)
PROCESS(`e1,E1', `Z,Z'		,  3000 GeV,  2.4906532E+01,  1.06E-02,  3)
PROCESS(`e1,E1', `se1,SE1'	,  3000 GeV,  5.8568844E+01,  8.84E-02,  3)
PROCESS(`e1,E1', `se2,SE2'	,  3000 GeV,  6.8229901E+01,  1.13E-01,  3)
PROCESS(`e1,E1', `se1,SE2'	,  3000 GeV,  3.3698124E+00,  4.78E-03,  3)
PROCESS(`e1,E1', `smu1,SMU1'	,  3000 GeV,  3.1962762E+00,  2.77E-03,  3)
PROCESS(`e1,E1', `smu2,SMU2'	,  3000 GeV,  2.8669086E+00,  2.44E-03,  3)
PROCESS(`e1,E1', `stau1,STAU1'	,  3000 GeV,  2.7144097E+00,  2.28E-03,  3)
PROCESS(`e1,E1', `stau2,STAU2'	,  3000 GeV,  2.6915841E+00,  2.31E-03,  3)
PROCESS(`e1,E1', `stau1,STAU2'	,  3000 GeV,  3.2991385E-01,  2.81E-04,  3)
PROCESS(`e1,E1', `sn11,SN11'	,  3000 GeV,  1.7332945E+02,  2.12E-01,  3)
PROCESS(`e1,E1', `sn21,SN21'	,  3000 GeV,  1.3243721E+00,  1.12E-03,  3)
PROCESS(`e1,E1', `sn31,SN31'	,  3000 GeV,  1.3263646E+00,  1.11E-03,  3)
PROCESS(`e1,E1', `su1,su1c'	,  3000 GeV,  2.3819439E+00,  1.99E-03,  3)
PROCESS(`e1,E1', `su2,su2c'	,  3000 GeV,  1.5652821E+00,  1.33E-03,  3)
PROCESS(`e1,E1', `sc1,sc1c'	,  3000 GeV,  2.3819439E+00,  2.03E-03,  3)
PROCESS(`e1,E1', `sc2,sc2c'	,  3000 GeV,  1.5665966E+00,  1.33E-03,  3)
PROCESS(`e1,E1', `st1,st1c'	,  3000 GeV,  2.3914396E+00,  2.04E-03,  3)
PROCESS(`e1,E1', `st2,st2c'	,  3000 GeV,  7.9873153E-01,  6.75E-04,  3)
PROCESS(`e1,E1', `st1,st2c'	,  3000 GeV,  4.0418989E-01,  3.47E-04,  3)
PROCESS(`e1,E1', `sd1,sd1c'	,  3000 GeV,  1.6048294E+00,  1.36E-03,  3)
PROCESS(`e1,E1', `sd2,sd2c'	,  3000 GeV,  3.9136515E-01,  3.30E-04,  3)
PROCESS(`e1,E1', `ss1,ss1c'	,  3000 GeV,  1.6044757E+00,  1.37E-03,  3)
PROCESS(`e1,E1', `ss2,ss2c'	,  3000 GeV,  3.9098748E-01,  3.32E-04,  3)
PROCESS(`e1,E1', `sb1,sb1c'	,  3000 GeV,  5.4571906E-01,  4.70E-04,  3)
PROCESS(`e1,E1', `sb2,sb2c'	,  3000 GeV,  6.4790435E-01,  5.46E-04,  3)
PROCESS(`e1,E1', `sb1,sb2c'	,  3000 GeV,  4.0084687E-01,  3.45E-04,  3)
PROCESS(`e1,E1', `neu1,neu1'	,  3000 GeV,  1.2048218E+01,  1.39E-02,  3)
PROCESS(`e1,E1', `neu1,neu2'	,  3000 GeV,  1.7009856E+00,  1.96E-03,  3)
PROCESS(`e1,E1', `neu1,neu3'	,  3000 GeV,  3.2924906E+00,  3.83E-03,  3)
PROCESS(`e1,E1', `neu1,neu4'	,  3000 GeV,  2.4130666E-01,  2.52E-04,  3)
PROCESS(`e1,E1', `neu1,neu5'	,  3000 GeV,  1.0399155E-01,  9.60E-05,  3)
PROCESS(`e1,E1', `neu2,neu2'	,  3000 GeV,  2.5286887E-01,  2.91E-04,  3)
PROCESS(`e1,E1', `neu2,neu3'	,  3000 GeV,  2.4889892E+00,  2.90E-03,  3)
PROCESS(`e1,E1', `neu2,neu4'	,  3000 GeV,  1.3103275E-01,  1.17E-04,  3)
PROCESS(`e1,E1', `neu2,neu5'	,  3000 GeV,  5.7754667E-02,  5.36E-05,  3)
PROCESS(`e1,E1', `neu3,neu3'	,  3000 GeV,  8.4042542E+00,  9.74E-03,  3)
PROCESS(`e1,E1', `neu3,neu4'	,  3000 GeV,  4.9416895E-01,  5.44E-04,  3)
PROCESS(`e1,E1', `neu3,neu5'	,  3000 GeV,  1.6679747E-01,  1.49E-04,  3)
PROCESS(`e1,E1', `neu4,neu4'	,  3000 GeV,  3.2678974E-03,  1.64E-06,  3)
PROCESS(`e1,E1', `neu4,neu5'	,  3000 GeV,  4.8893282E+00,  4.57E-03,  3)
PROCESS(`e1,E1', `neu5,neu5'	,  3000 GeV,  5.8502305E-05,  3.19E-08,  3)
PROCESS(`e1,E1', `"ch1+","ch1-"',  3000 GeV,  2.2850785E+01,  2.89E-02,  3)
PROCESS(`e1,E1', `"ch2+","ch2-"',  3000 GeV,  1.2598006E+01,  1.03E-02,  3)
PROCESS(`e1,E1', `"ch1+","ch2-"',  3000 GeV,  7.5164320E-01,  6.16E-04,  3)
PROCESS(`e1,E1', `Z,h01'	,  3000 GeV,  1.3302479E+00,  1.13E-03,  3)
PROCESS(`e1,E1', `Z,h02'	,  3000 GeV,  2.8504985E-02,  2.40E-05,  3)
PROCESS(`e1,E1', `Z,h03'	,  3000 GeV,  3.8558359E-07,  3.21E-10,  3)
PROCESS(`e1,E1', `A01,h01'	,  3000 GeV,  7.7307473E-07,  6.59E-10,  3)
PROCESS(`e1,E1', `A01,h02'	,  3000 GeV,  3.0387620E-05,  2.61E-08,  3)
PROCESS(`e1,E1', `A01,h03'	,  3000 GeV,  5.0789448E-03,  4.34E-06,  3)
PROCESS(`e1,E1', `A02,h01'	,  3000 GeV,  5.7728434E-05,  4.95E-08,  3)
PROCESS(`e1,E1', `A02,h02'	,  3000 GeV,  2.2687853E-03,  1.92E-06,  3)
PROCESS(`e1,E1', `A02,h03'	,  3000 GeV,  1.6113881E-01,  1.38E-04,  3)
PROCESS(`e1,E1', `Hp,Hm'	,  3000 GeV,  3.9652879E-01,  3.38E-04,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
