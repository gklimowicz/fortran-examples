dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-tn.m4 > nmssm_ext-tn.sin
dnl   whizard -r nmssm_ext-tn.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-tn.m4 $, nmssm_ext_tn_)
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
PROCESS(`e3,N3', `Z,Wm'		,  3000 GeV,  7.0808692E+02,  1.22E+00,  3)
PROCESS(`e3,N3', `se1,SN11'	,  3000 GeV,  1.2377155E+01,  1.05E-02,  3)
PROCESS(`e3,N3', `smu1,SN21'	,  3000 GeV,  1.2351754E+01,  1.07E-02,  3)
PROCESS(`e3,N3', `stau1,SN31'	,  3000 GeV,  3.9002258E+01,  3.82E-02,  3)
PROCESS(`e3,N3', `stau2,SN31'	,  3000 GeV,  3.9080046E+01,  3.96E-02,  3)
PROCESS(`e3,N3', `"ch1-",neu1'	,  3000 GeV,  2.0379344E+01,  2.42E-02,  3)
PROCESS(`e3,N3', `"ch1-",neu2'	,  3000 GeV,  2.8031970E+00,  2.80E-03,  3)
PROCESS(`e3,N3', `"ch1-",neu3'	,  3000 GeV,  2.0592216E+01,  3.00E-02,  3)
PROCESS(`e3,N3', `"ch1-",neu4'	,  3000 GeV,  1.1275763E+00,  8.26E-04,  3)
PROCESS(`e3,N3', `"ch1-",neu5'	,  3000 GeV,  6.9514150E-01,  6.62E-04,  3)
PROCESS(`e3,N3', `"ch2-",neu1'	,  3000 GeV,  1.0088133E+00,  1.08E-03,  3)
dnl Heavy downward fluctuation of WHIZARD for the FeynRules comparison
PROCESS(`e3,N3', `"ch2-",neu2'	,  3000 GeV,  8.1168194E-01,  5.10E-04,  3)
PROCESS(`e3,N3', `"ch2-",neu3'	,  3000 GeV,  1.0359709E+00,  7.66E-04,  3)
PROCESS(`e3,N3', `"ch2-",neu4'	,  3000 GeV,  2.3329018E+01,  2.34E-02,  3)
PROCESS(`e3,N3', `"ch2-",neu5'	,  3000 GeV,  2.3673460E+01,  2.24E-02,  3)
PROCESS(`e3,N3', `Wm,h01'	,  3000 GeV,  6.3382376E+00,  9.94E-03,  3)
PROCESS(`e3,N3', `Wm,h02'	,  3000 GeV,  1.3776369E-01,  2.22E-04,  3)
PROCESS(`e3,N3', `Wm,h03'	,  3000 GeV,  3.2809450E-01,  5.36E-04,  3)
PROCESS(`e3,N3', `Wm,A01'	,  3000 GeV,  1.3482623E-02,  1.25E-05,  3)
PROCESS(`e3,N3', `Wm,A02'	,  3000 GeV,  3.2837822E-01,  5.24E-04,  3)
PROCESS(`e3,N3', `Hm,h01'	,  3000 GeV,  6.1301638E-04,  8.00E-07,  3)
PROCESS(`e3,N3', `Hm,h02'	,  3000 GeV,  3.0119642E-02,  4.08E-05,  3)
PROCESS(`e3,N3', `Hm,h03'	,  3000 GeV,  7.5923322E-01,  1.22E-03,  3)
PROCESS(`e3,N3', `Hm,A01'	,  3000 GeV,  2.4992130E-02,  3.96E-05,  3)
PROCESS(`e3,N3', `Hm,A02'	,  3000 GeV,  7.4864918E-01,  1.20E-03,  3)
PROCESS(`e3,N3', `Z,Hm'		,  3000 GeV,  5.6412838E-01,  6.56E-04,  3)
PROCESS(`e3,N3', `Z,Wm'		,  5000 GeV,  3.0075190E+02,  6.00E-01,  3)
PROCESS(`e3,N3', `se1,SN11'	,  5000 GeV,  4.5270622E+00,  3.84E-03,  3)
PROCESS(`e3,N3', `smu1,SN21'	,  5000 GeV,  4.5264314E+00,  3.82E-03,  3)
PROCESS(`e3,N3', `stau1,SN31'	,  5000 GeV,  1.5884301E+01,  1.58E-02,  3)
PROCESS(`e3,N3', `stau2,SN31'	,  5000 GeV,  1.6052889E+01,  1.58E-02,  3)
PROCESS(`e3,N3', `"ch1-",neu1'	,  5000 GeV,  7.5353728E+00,  9.50E-03,  3)
PROCESS(`e3,N3', `"ch1-",neu2'	,  5000 GeV,  1.0480540E+00,  1.58E-03,  3)
PROCESS(`e3,N3', `"ch1-",neu3'	,  5000 GeV,  7.8071594E+00,  1.22E-02,  3)
PROCESS(`e3,N3', `"ch1-",neu4'	,  5000 GeV,  3.5652074E-01,  3.20E-04,  3)
PROCESS(`e3,N3', `"ch1-",neu5'	,  5000 GeV,  2.5527442E-01,  2.78E-04,  3)
PROCESS(`e3,N3', `"ch2-",neu1'	,  5000 GeV,  3.5872782E-01,  4.34E-04,  3)
PROCESS(`e3,N3', `"ch2-",neu2'	,  5000 GeV,  2.8608974E-01,  3.24E-04,  3)
PROCESS(`e3,N3', `"ch2-",neu3'	,  5000 GeV,  3.3569090E-01,  3.12E-04,  3)
PROCESS(`e3,N3', `"ch2-",neu4'	,  5000 GeV,  8.6069194E+00,  9.64E-03,  3)
PROCESS(`e3,N3', `"ch2-",neu5'	,  5000 GeV,  8.5530408E+00,  9.42E-03,  3)
PROCESS(`e3,N3', `Wm,h01'	,  5000 GeV,  2.2808574E+00,  3.70E-03,  3)
PROCESS(`e3,N3', `Wm,h02'	,  5000 GeV,  4.8330392E-02,  8.00E-05,  3)
PROCESS(`e3,N3', `Wm,h03'	,  5000 GeV,  1.3055033E-01,  2.32E-04,  3)
PROCESS(`e3,N3', `Wm,A01'	,  5000 GeV,  2.3628748E-03,  3.24E-06,  3)
PROCESS(`e3,N3', `Wm,A02'	,  5000 GeV,  1.3035782E-01,  2.36E-04,  3)
PROCESS(`e3,N3', `Hm,h01'	,  5000 GeV,  1.6391130E-04,  2.52E-07,  4)
PROCESS(`e3,N3', `Hm,h02'	,  5000 GeV,  7.9993172E-03,  1.22E-05,  3)
PROCESS(`e3,N3', `Hm,h03'	,  5000 GeV,  1.3996335E+00,  2.34E-03,  3)
PROCESS(`e3,N3', `Hm,A01'	,  5000 GeV,  1.3173624E-02,  2.14E-05,  3)
PROCESS(`e3,N3', `Hm,A02'	,  5000 GeV,  1.3976411E+00,  2.30E-03,  3)
PROCESS(`e3,N3', `Z,Hm'		,  5000 GeV,  2.1952906E-01,  2.80E-04,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
