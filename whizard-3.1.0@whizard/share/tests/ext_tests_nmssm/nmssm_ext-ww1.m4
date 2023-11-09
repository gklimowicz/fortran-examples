dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-ww1.m4 > nmssm_ext-ww1.sin
dnl   whizard -r nmssm_ext-ww1.sin
dnl
dnl ----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-ww1.m4 $, nmssm_ext_ww1_)
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
PROCESS(`Wm,Wp', `t,T'		,  3000 GeV,  5.2959795E+03,  6.12E+00,  3)
PROCESS(`Wm,Wp', `Z,Z'		,  3000 GeV,  3.0247597E+05,  6.89E+02,  3)
PROCESS(`Wm,Wp', `se1,SE1'	,  3000 GeV,  1.1329904E+01,  1.10E-02,  3)
PROCESS(`Wm,Wp', `se2,SE2'	,  3000 GeV,  2.7012580E-01,  2.80E-04,  3)
dnl Wrong ordering of smuons in the original WHIZARD/FeynRules comparison
PROCESS(`Wm,Wp', `smu1,SMU1'	,  3000 GeV,  1.1135872E+01,  7.39E-04,  3)
PROCESS(`Wm,Wp', `smu2,SMU2'	,  3000 GeV,  2.7000110E-01,  2.83E-04,  3)
PROCESS(`Wm,Wp', `stau1,STAU1'	,  3000 GeV,  2.8487449E+00,  2.57E-03,  3)
PROCESS(`Wm,Wp', `stau2,STAU2'	,  3000 GeV,  2.9971964E+00,  2.77E-03,  3)
PROCESS(`Wm,Wp', `stau1,STAU2'	,  3000 GeV,  2.9221293E+00,  2.87E-03,  3)
PROCESS(`Wm,Wp', `sn11,SN11'	,  3000 GeV,  1.0571424E+01,  1.29E-02,  3)
dnl Wrong ordering of smuons in the original WHIZARD/FeynRules comparison
PROCESS(`Wm,Wp', `sn21,SN21'	,  3000 GeV,  1.0459269E+01,  9.72E-04,  3)
PROCESS(`Wm,Wp', `sn31,SN31'	,  3000 GeV,  1.0537059E+01,  1.32E-02,  3)
PROCESS(`Wm,Wp', `su1,su1c'	,  3000 GeV,  1.3906157E+01,  1.03E-02,  3)
PROCESS(`Wm,Wp', `su2,su2c'	,  3000 GeV,  2.0989466E-01,  2.18E-04,  3)
PROCESS(`Wm,Wp', `sc1,sc1c'	,  3000 GeV,  1.3893551E+01,  1.05E-02,  3)
PROCESS(`Wm,Wp', `sc2,sc2c'	,  3000 GeV,  2.0939722E-01,  2.21E-04,  3)
PROCESS(`Wm,Wp', `st1,st1c'	,  3000 GeV,  3.4799086E+02,  1.62E-01,  3)
PROCESS(`Wm,Wp', `st2,st2c'	,  3000 GeV,  1.1205768E+02,  5.20E-02,  3)
PROCESS(`Wm,Wp', `st1,st2c'	,  3000 GeV,  1.3950058E+02,  1.39E-01,  3)
PROCESS(`Wm,Wp', `sd1,sd1c'	,  3000 GeV,  1.3417266E+01,  5.98E-03,  3)
PROCESS(`Wm,Wp', `sd2,sd2c'	,  3000 GeV,  5.2428831E-02,  5.48E-05,  3)
PROCESS(`Wm,Wp', `ss1,ss1c'	,  3000 GeV,  1.3419887E+01,  5.94E-03,  3)
PROCESS(`Wm,Wp', `ss2,ss2c'	,  3000 GeV,  5.2525094E-02,  5.35E-05,  3)
PROCESS(`Wm,Wp', `sb1,sb1c'	,  3000 GeV,  1.0074353E+02,  9.23E-02,  3)
PROCESS(`Wm,Wp', `sb2,sb2c'	,  3000 GeV,  1.0224847E+02,  1.20E-01,  3)
PROCESS(`Wm,Wp', `sb1,sb2c'	,  3000 GeV,  1.0507622E+02,  1.14E-01,  3)
PROCESS(`Wm,Wp', `neu1,neu1'	,  3000 GeV,  8.2937872E-01,  3.66E-04,  3)
PROCESS(`Wm,Wp', `neu1,neu2'	,  3000 GeV,  9.3181396E+00,  4.43E-03,  3)
PROCESS(`Wm,Wp', `neu1,neu3'	,  3000 GeV,  2.7325608E+01,  1.47E-02,  3)
PROCESS(`Wm,Wp', `neu1,neu4'	,  3000 GeV,  6.9325380E+00,  3.14E-03,  3)
PROCESS(`Wm,Wp', `neu1,neu5'	,  3000 GeV,  8.1595608E+00,  4.14E-03,  3)
PROCESS(`Wm,Wp', `neu2,neu2'	,  3000 GeV,  2.3267567E+01,  1.93E-02,  3)
PROCESS(`Wm,Wp', `neu2,neu3'	,  3000 GeV,  1.1910560E+02,  3.95E-02,  3)
PROCESS(`Wm,Wp', `neu2,neu4'	,  3000 GeV,  5.9961958E+01,  3.10E-02,  3)
PROCESS(`Wm,Wp', `neu2,neu5'	,  3000 GeV,  5.8876960E+01,  2.82E-02,  3)
PROCESS(`Wm,Wp', `neu3,neu3'	,  3000 GeV,  5.9103675E+02,  3.23E-01,  3)
PROCESS(`Wm,Wp', `neu3,neu4'	,  3000 GeV,  2.1384803E+02,  6.39E-02,  3)
PROCESS(`Wm,Wp', `neu3,neu5'	,  3000 GeV,  2.0551958E+02,  6.53E-02,  3)
PROCESS(`Wm,Wp', `neu4,neu4'	,  3000 GeV,  3.9232297E+01,  2.82E-02,  3)
PROCESS(`Wm,Wp', `neu4,neu5'	,  3000 GeV,  5.5256662E+01,  4.26E-02,  3)
PROCESS(`Wm,Wp', `neu5,neu5'	,  3000 GeV,  3.4148657E+01,  2.44E-02,  3)
PROCESS(`Wm,Wp', `"ch1+","ch1-"',  3000 GeV,  6.2084243E+02,  8.54E-01,  3)
PROCESS(`Wm,Wp', `"ch2+","ch2-"',  3000 GeV,  1.8256030E+02,  2.32E-01,  3)
PROCESS(`Wm,Wp', `"ch1+","ch2-"',  3000 GeV,  8.7939214E+01,  1.14E-01,  3)
PROCESS(`Wm,Wp', `h01,h01'	,  3000 GeV,  5.7425857E+03,  7.39E+01,  3)
PROCESS(`Wm,Wp', `h01,h02'	,  3000 GeV,  2.4757009E+02,  1.99E+00,  3)
PROCESS(`Wm,Wp', `h01,h03'	,  3000 GeV,  2.8255891E+00,  1.06E-03,  3)
PROCESS(`Wm,Wp', `h02,h02'	,  3000 GeV,  7.9026180E+00,  4.88E-02,  3)
PROCESS(`Wm,Wp', `h02,h03'	,  3000 GeV,  1.1093989E+00,  4.87E-04,  3)
PROCESS(`Wm,Wp', `h03,h03'	,  3000 GeV,  5.4109278E+00,  1.52E-03,  3)
PROCESS(`Wm,Wp', `A01,A01'	,  3000 GeV,  6.5571884E+00,  3.38E-03,  3)
PROCESS(`Wm,Wp', `A01,A02'	,  3000 GeV,  1.4810855E+00,  9.89E-04,  3)
PROCESS(`Wm,Wp', `A02,A02'	,  3000 GeV,  7.6283505E+00,  2.04E-03,  3)
PROCESS(`Wm,Wp', `Z,h01'	,  3000 GeV,  8.2814733E+04,  1.86E+02,  3)
PROCESS(`Wm,Wp', `Z,h02'	,  3000 GeV,  1.7735179E+03,  3.99E+00,  3)
PROCESS(`Wm,Wp', `Z,h03'	,  3000 GeV,  3.5660809E-02,  7.34E-05,  3)
PROCESS(`Wm,Wp', `A01,h01'	,  3000 GeV,  5.8535711E-03,  4.74E-06,  3)
PROCESS(`Wm,Wp', `A01,h02'	,  3000 GeV,  2.2991862E-01,  1.83E-04,  3)
PROCESS(`Wm,Wp', `A01,h03'	,  3000 GeV,  1.0786432E+00,  9.25E-04,  3)
PROCESS(`Wm,Wp', `A02,h01'	,  3000 GeV,  1.1968002E-02,  1.04E-05,  3)
PROCESS(`Wm,Wp', `A02,h02'	,  3000 GeV,  4.6867557E-01,  4.02E-04,  3)
PROCESS(`Wm,Wp', `A02,h03'	,  3000 GeV,  2.9370314E-01,  3.08E-04,  3)
PROCESS(`Wm,Wp', `Hp,Hm'	,  3000 GeV,  1.4338390E+01,  1.24E-02,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
