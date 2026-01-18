dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-tt1.m4 > nmssm_ext-tt1.sin
dnl   whizard -r nmssm_ext-tt1.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-tt1.m4 $, nmssm_ext_tt1_)
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
PROCESS(`e3,E3', `n1,N1'	,  3000 GeV,  2.7145945E+00,  4.72E-04,  3)
PROCESS(`e3,E3', `n2,N2'	,  3000 GeV,  2.7156229E+00,  4.46E-04,  3)
PROCESS(`e3,E3', `n3,N3'	,  3000 GeV,  5.3273824E+04,  1.46E+02,  3)
PROCESS(`e3,E3', `t,T'	 	,  3000 GeV,  1.9730048E+01,  7.70E-03,  3)
PROCESS(`e3,E3', `Wp,Wm'	,  3000 GeV,  4.5525676E+02,  9.19E-01,  3)
PROCESS(`e3,E3', `Z,Z'		,  3000 GeV,  2.4815330E+01,  2.74E-02,  3)
PROCESS(`e3,E3', `se1,SE1'	,  3000 GeV,  3.1958583E+00,  2.77E-03,  3)
PROCESS(`e3,E3', `se2,SE2'	,  3000 GeV,  2.8655474E+00,  2.48E-03,  3)
PROCESS(`e3,E3', `smu1,SMU1'	,  3000 GeV,  3.1952694E+00,  2.75E-03,  3)
PROCESS(`e3,E3', `smu2,SMU2'	,  3000 GeV,  2.8687442E+00,  2.43E-03,  3)
PROCESS(`e3,E3', `stau1,STAU1'	,  3000 GeV,  3.1478840E+01,  5.01E-02,  3)
PROCESS(`e3,E3', `stau2,STAU2'	,  3000 GeV,  3.0631791E+01,  4.88E-02,  3)
PROCESS(`e3,E3', `stau1,STAU2'	,  3000 GeV,  3.5761187E+01,  4.89E-02,  3)
PROCESS(`e3,E3', `sn11,SN11'	,  3000 GeV,  1.3231619E+00,  1.14E-03,  3)
PROCESS(`e3,E3', `sn21,SN21'	,  3000 GeV,  1.3238730E+00,  1.14E-03,  3)
PROCESS(`e3,E3', `sn31,SN31'	,  3000 GeV,  1.7333889E+02,  2.11E-01,  3)
PROCESS(`e3,E3', `su1,su1c'	,  3000 GeV,  2.3832865E+00,  2.03E-03,  3)
PROCESS(`e3,E3', `su2,su2c'	,  3000 GeV,  1.5642847E+00,  1.36E-03,  3)
PROCESS(`e3,E3', `sc1,sc1c'	,  3000 GeV,  2.3804980E+00,  2.05E-03,  3)
PROCESS(`e3,E3', `sc2,sc2c'	,  3000 GeV,  1.5651887E+00,  1.34E-03,  3)
PROCESS(`e3,E3', `st1,st1c'	,  3000 GeV,  2.3961737E+00,  2.04E-03,  3)
PROCESS(`e3,E3', `st2,st2c'	,  3000 GeV,  8.0062838E-01,  6.83E-04,  3)
PROCESS(`e3,E3', `st1,st2c'	,  3000 GeV,  4.0595365E-01,  3.42E-04,  3)
PROCESS(`e3,E3', `sd1,sd1c'	,  3000 GeV,  1.6059001E+00,  1.36E-03,  3)
PROCESS(`e3,E3', `sd2,sd2c'	,  3000 GeV,  3.9077742E-01,  3.31E-04,  3)
PROCESS(`e3,E3', `ss1,ss1c'	,  3000 GeV,  1.6040548E+00,  1.37E-03,  3)
PROCESS(`e3,E3', `ss2,ss2c'	,  3000 GeV,  3.9050065E-01,  3.38E-04,  3)
PROCESS(`e3,E3', `sb1,sb1c'	,  3000 GeV,  5.4589993E-01,  4.70E-04,  3)
PROCESS(`e3,E3', `sb2,sb2c'	,  3000 GeV,  6.4672358E-01,  5.56E-04,  3)
PROCESS(`e3,E3', `sb1,sb2c'	,  3000 GeV,  4.0165701E-01,  3.43E-04,  3)
PROCESS(`e3,E3', `neu1,neu1'	,  3000 GeV,  1.2073446E+01,  1.38E-02,  3)
PROCESS(`e3,E3', `neu1,neu2'	,  3000 GeV,  1.7025344E+00,  2.01E-03,  3)
PROCESS(`e3,E3', `neu1,neu3'	,  3000 GeV,  3.2959035E+00,  3.86E-03,  3)
PROCESS(`e3,E3', `neu1,neu4'	,  3000 GeV,  2.4716418E-01,  2.67E-04,  3)
PROCESS(`e3,E3', `neu1,neu5'	,  3000 GeV,  1.1051274E-01,  1.05E-04,  3)
PROCESS(`e3,E3', `neu2,neu2'	,  3000 GeV,  2.6260245E-01,  3.04E-04,  3)
PROCESS(`e3,E3', `neu2,neu3'	,  3000 GeV,  2.4921765E+00,  2.92E-03,  3)
PROCESS(`e3,E3', `neu2,neu4'	,  3000 GeV,  1.8711468E-01,  1.78E-04,  3)
PROCESS(`e3,E3', `neu2,neu5'	,  3000 GeV,  1.1176564E-01,  1.13E-04,  3)
PROCESS(`e3,E3', `neu3,neu3'	,  3000 GeV,  8.3805914E+00,  1.02E-02,  3)
PROCESS(`e3,E3', `neu3,neu4'	,  3000 GeV,  4.9803015E-01,  5.62E-04,  3)
PROCESS(`e3,E3', `neu3,neu5'	,  3000 GeV,  1.7336478E-01,  1.56E-04,  3)
PROCESS(`e3,E3', `neu4,neu4'	,  3000 GeV,  4.8029379E-03,  2.77E-06,  3)
PROCESS(`e3,E3', `neu4,neu5'	,  3000 GeV,  4.8487241E+00,  4.75E-03,  3)
PROCESS(`e3,E3', `neu5,neu5'	,  3000 GeV,  4.2191208E-03,  4.67E-06,  3)
PROCESS(`e3,E3', `"ch1+","ch1-"',  3000 GeV,  2.2812760E+01,  2.90E-02,  3)
PROCESS(`e3,E3', `"ch2+","ch2-"',  3000 GeV,  1.2589218E+01,  1.03E-02,  3)
PROCESS(`e3,E3', `"ch1+","ch2-"',  3000 GeV,  7.5689892E-01,  6.34E-04,  3)
PROCESS(`e3,E3', `h01,h01'	,  3000 GeV,  1.2810537E-05,  1.30E-08,  3)
PROCESS(`e3,E3', `h01,h02'	,  3000 GeV,  2.1913209E-03,  3.29E-06,  3)
PROCESS(`e3,E3', `h01,h03'	,  3000 GeV,  5.5587321E-04,  7.99E-07,  3)
PROCESS(`e3,E3', `h02,h02'	,  3000 GeV,  4.9996639E-04,  7.82E-07,  4)
PROCESS(`e3,E3', `h02,h03'	,  3000 GeV,  9.1787245E-03,  1.40E-05,  3)
PROCESS(`e3,E3', `h03,h03'	,  3000 GeV,  7.3340440E-05,  1.12E-07,  3)
PROCESS(`e3,E3', `A01,A01'	,  3000 GeV,  1.9560453E-06,  2.90E-09,  3)
PROCESS(`e3,E3', `A01,A02'	,  3000 GeV,  7.2787463E-04,  1.10E-06,  3)
PROCESS(`e3,E3', `A02,A02'	,  3000 GeV,  9.6466011E-06,  1.46E-08,  3)
PROCESS(`e3,E3', `Z,h01'	,  3000 GeV,  1.3600289E+00,  2.42E-03,  3)
PROCESS(`e3,E3', `Z,h02'	,  3000 GeV,  3.0832631E-02,  5.67E-05,  3)
PROCESS(`e3,E3', `Z,h03'	,  3000 GeV,  8.7948382E-02,  1.48E-04,  3)
PROCESS(`e3,E3', `Z,A01'	,  3000 GeV,  5.9808786E-03,  6.72E-06,  3)
PROCESS(`e3,E3', `Z,A02'	,  3000 GeV,  8.8054045E-02,  1.46E-04,  3)
PROCESS(`e3,E3', `A01,h01'	,  3000 GeV,  4.8926322E-03,  7.42E-06,  3)
PROCESS(`e3,E3', `A01,h02'	,  3000 GeV,  5.3028073E-04,  8.09E-07,  3)
PROCESS(`e3,E3', `A01,h03'	,  3000 GeV,  5.6160530E-03,  1.03E-05,  3)
PROCESS(`e3,E3', `A02,h01'	,  3000 GeV,  8.3238915E-04,  1.23E-06,  3)
PROCESS(`e3,E3', `A02,h02'	,  3000 GeV,  1.1460292E-02,  1.82E-05,  3)
PROCESS(`e3,E3', `A02,h03'	,  3000 GeV,  1.5714039E-01,  2.99E-04,  3)
PROCESS(`e3,E3', `Hp,Hm'	,  3000 GeV,  3.9487341E-01,  6.36E-04,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
