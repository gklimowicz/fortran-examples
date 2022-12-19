dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-uu2.m4 > nmssm_ext-uu2.sin
dnl   whizard -r nmssm_ext-uu2.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-uu2.m4 ohl $, nmssm_ext_uu2_)
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
PROCESS(`u,U', `n1,N1'		,  5000 GeV,  3.7009049E-01,  1.04E-04,  3)
PROCESS(`u,U', `n2,N2'		,  5000 GeV,  3.6989807E-01,  1.05E-04,  3)
PROCESS(`u,U', `n3,N3'		,  5000 GeV,  3.7031276E-01,  9.85E-05,  3)
PROCESS(`u,U', `t,T'		,  5000 GeV,  2.0073151E+02,  3.13E-02,  3)
PROCESS(`u,U', `Wp,Wm'		,  5000 GeV,  6.3697675E+01,  1.55E-01,  3)
PROCESS(`u,U', `Z,Z'		,  5000 GeV,  6.3953137E+00,  6.77E-03,  3)
PROCESS(`u,U', `se1,SE1'	,  5000 GeV,  2.5945571E-01,  2.23E-04,  3)
PROCESS(`u,U', `se2,SE2'	,  5000 GeV,  1.3241413E-01,  1.11E-04,  3)
PROCESS(`u,U', `smu1,SMU1'	,  5000 GeV,  2.5959150E-01,  2.21E-04,  3)
PROCESS(`u,U', `smu2,SMU2'	,  5000 GeV,  1.3209038E-01,  1.14E-04,  3)
PROCESS(`u,U', `stau1,STAU1'	,  5000 GeV,  1.4968850E-01,  1.25E-04,  3)
PROCESS(`u,U', `stau2,STAU2'	,  5000 GeV,  1.5091592E-01,  1.28E-04,  3)
PROCESS(`u,U', `stau1,STAU2'	,  5000 GeV,  4.5791472E-02,  3.89E-05,  3)
PROCESS(`u,U', `sn11,SN11'	,  5000 GeV,  1.8322852E-01,  1.57E-04,  3)
PROCESS(`u,U', `sn21,SN21'	,  5000 GeV,  1.8340824E-01,  1.55E-04,  3)
PROCESS(`u,U', `sn31,SN31'	,  5000 GeV,  1.8326038E-01,  1.57E-04,  3)
PROCESS(`u,U', `su1,su1c'	,  5000 GeV,  3.5220701E+02,  3.86E-01,  3)
PROCESS(`u,U', `su2,su2c'	,  5000 GeV,  3.6018003E+02,  3.96E-01,  3)
PROCESS(`u,U', `su1,su2c'	,  5000 GeV,  1.2986802E+02,  1.58E-01,  3)
PROCESS(`u,U', `sc1,sc1c'	,  5000 GeV,  3.8502204E+01,  3.26E-02,  3)
PROCESS(`u,U', `sc2,sc2c'	,  5000 GeV,  3.8090527E+01,  3.27E-02,  3)
PROCESS(`u,U', `st1,st1c'	,  5000 GeV,  4.2655592E+01,  3.60E-02,  3)
PROCESS(`u,U', `st2,st2c'	,  5000 GeV,  3.3401485E+01,  2.91E-02,  3)
PROCESS(`u,U', `st1,st2c'	,  5000 GeV,  1.0536113E-01,  9.07E-05,  3)
PROCESS(`u,U', `sd1,sd1c'	,  5000 GeV,  5.3493499E+01,  6.24E-02,  3)
PROCESS(`u,U', `sd2,sd2c'	,  5000 GeV,  3.8018002E+01,  5.50E-02,  3)
PROCESS(`u,U', `ss1,ss1c'	,  5000 GeV,  3.8371479E+01,  3.25E-02,  3)
PROCESS(`u,U', `ss2,ss2c'	,  5000 GeV,  3.8009791E+01,  3.21E-02,  3)
PROCESS(`u,U', `sb1,sb1c'	,  5000 GeV,  3.8156960E+01,  3.26E-02,  3)
PROCESS(`u,U', `sb2,sb2c'	,  5000 GeV,  3.8049929E+01,  3.19E-02,  3)
PROCESS(`u,U', `sb1,sb2c'	,  5000 GeV,  1.0476507E-01,  9.08E-05,  3)
PROCESS(`u,U', `neu1,neu1'	,  5000 GeV,  2.2487077E-01,  1.84E-04,  3)
PROCESS(`u,U', `neu1,neu2'	,  5000 GeV,  3.0076160E-02,  2.39E-05,  3)
PROCESS(`u,U', `neu1,neu3'	,  5000 GeV,  8.1339629E-03,  7.08E-06,  3)
PROCESS(`u,U', `neu1,neu4'	,  5000 GeV,  3.7161300E-03,  2.93E-06,  3)
PROCESS(`u,U', `neu1,neu5'	,  5000 GeV,  9.3125034E-03,  4.44E-06,  3)
PROCESS(`u,U', `neu2,neu2'	,  5000 GeV,  3.9093471E-03,  3.71E-06,  3)
PROCESS(`u,U', `neu2,neu3'	,  5000 GeV,  1.2237136E-01,  1.02E-04,  3)
PROCESS(`u,U', `neu2,neu4'	,  5000 GeV,  3.5578949E-02,  2.01E-05,  3)
PROCESS(`u,U', `neu2,neu5'	,  5000 GeV,  5.7623756E-03,  2.60E-06,  3)
PROCESS(`u,U', `neu3,neu3'	,  5000 GeV,  7.6764068E-01,  6.16E-04,  3)
PROCESS(`u,U', `neu3,neu4'	,  5000 GeV,  2.5811051E-02,  2.47E-05,  3)
PROCESS(`u,U', `neu3,neu5'	,  5000 GeV,  8.5993279E-03,  3.78E-06,  3)
PROCESS(`u,U', `neu4,neu4'	,  5000 GeV,  8.1071844E-04,  7.10E-07,  3)
PROCESS(`u,U', `neu4,neu5'	,  5000 GeV,  6.9266636E-01,  3.45E-04,  3)
PROCESS(`u,U', `neu5,neu5'	,  5000 GeV,  1.9205888E-05,  1.06E-08,  3)
PROCESS(`u,U', `"ch1+","ch1-"'	,  5000 GeV,  2.1870334E+00,  2.44E-03,  3)
PROCESS(`u,U', `"ch2+","ch2-"'	,  5000 GeV,  9.7441659E-01,  3.85E-04,  3)
PROCESS(`u,U', `"ch1+","ch2-"'	,  5000 GeV,  1.0317157E-01,  8.36E-05,  3)
PROCESS(`u,U', `sgl,sgl'	,  5000 GeV,  3.4541427E+02,  2.79E-01,  3)
PROCESS(`u,U', `Z,h01'		,  5000 GeV,  1.8150194E-01,  1.51E-04,  3)
PROCESS(`u,U', `Z,h02'		,  5000 GeV,  3.8870524E-03,  3.26E-06,  3)
PROCESS(`u,U', `Z,h03'		,  5000 GeV,  7.9075716E-08,  6.75E-11,  3)
PROCESS(`u,U', `A01,h01'	,  5000 GeV,  1.0577662E-07,  9.08E-11,  3)
PROCESS(`u,U', `A01,h02'	,  5000 GeV,  4.1616917E-06,  3.55E-09,  3)
PROCESS(`u,U', `A01,h03'	,  5000 GeV,  1.0511429E-03,  8.98E-07,  3)
PROCESS(`u,U', `A02,h01'	,  5000 GeV,  1.2072132E-05,  1.02E-08,  3)
PROCESS(`u,U', `A02,h02'	,  5000 GeV,  4.7454137E-04,  4.10E-07,  3)
PROCESS(`u,U', `A02,h03'	,  5000 GeV,  1.1354355E-01,  9.80E-05,  3)
PROCESS(`u,U', `Hp,Hm'		,  5000 GeV,  1.6299905E-01,  1.38E-04,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
