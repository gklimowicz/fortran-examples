dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-ee.m4 > nmssm_ext-ee.sin
dnl   whizard -r nmssm_ext-ee
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-ee.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_ee2_)
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
PROCESS(`e1,E1', `n1,N1'	,  5000 GeV,  5.3542631E+04,  8.50E+01,  3)
PROCESS(`e1,E1', `n2,N2'	,  5000 GeV,  9.7653103E-01,  9.23E-05,  3)
PROCESS(`e1,E1', `n3,N3'	,  5000 GeV,  9.7635885E-01,  9.36E-05,  3)
PROCESS(`e1,E1', `t,T'		,  5000 GeV,  7.1159051E+00,  1.61E-03,  3)
PROCESS(`e1,E1', `Wp,Wm'	,  5000 GeV,  1.9122473E+02,  2.68E-01,  3)
PROCESS(`e1,E1', `Z,Z'		,  5000 GeV,  1.0478015E+01,  6.32E-03,  3)
PROCESS(`e1,E1', `se1,SE1'	,  5000 GeV,  2.7041505E+01,  2.57E-02,  3)
PROCESS(`e1,E1', `se2,SE2'	,  5000 GeV,  2.9977109E+01,  3.06E-02,  3)
PROCESS(`e1,E1', `se1,SE2'	,  5000 GeV,  1.2502181E+00,  1.06E-03,  3)
PROCESS(`e1,E1', `smu1,SMU1'	,  5000 GeV,  1.1715446E+00,  5.82E-04,  3)
PROCESS(`e1,E1', `smu2,SMU2'	,  5000 GeV,  1.0504992E+00,  5.15E-04,  3)
PROCESS(`e1,E1', `stau1,STAU1'	,  5000 GeV,  9.9168895E-01,  4.83E-04,  3)
PROCESS(`e1,E1', `stau2,STAU2'	,  5000 GeV,  9.8981565E-01,  4.88E-04,  3)
PROCESS(`e1,E1', `stau1,STAU2'	,  5000 GeV,  1.2073167E-01,  6.00E-05,  3)
PROCESS(`e1,E1', `sn11,SN11'	,  5000 GeV,  7.6639759E+01,  5.87E-02,  3)
PROCESS(`e1,E1', `sn21,SN21'	,  5000 GeV,  4.8370603E-01,  2.38E-04,  3)
PROCESS(`e1,E1', `sn31,SN31'	,  5000 GeV,  4.8377582E-01,  2.38E-04,  3)
PROCESS(`e1,E1', `sc2,sc2c'	,  5000 GeV,  1.0799307E+00,  5.27E-04,  3)
PROCESS(`e1,E1', `sc1,sc1c'	,  5000 GeV,  1.6385992E+00,  8.16E-04,  3)
PROCESS(`e1,E1', `su2,su2c'	,  5000 GeV,  1.0794030E+00,  5.28E-04,  3)
PROCESS(`e1,E1', `su1,su1c'	,  5000 GeV,  1.6399474E+00,  8.07E-04,  3)
PROCESS(`e1,E1', `st1,st1c'	,  5000 GeV,  1.2051313E+00,  5.97E-04,  3)
PROCESS(`e1,E1', `st2,st2c'	,  5000 GeV,  9.4604008E-01,  4.63E-04,  3)
PROCESS(`e1,E1', `st1,st2c'	,  5000 GeV,  2.7850741E-01,  1.35E-04,  3)
PROCESS(`e1,E1', `sd1,sd1c'	,  5000 GeV,  1.1071735E+00,  5.48E-04,  3)
PROCESS(`e1,E1', `ss1,ss1c'	,  5000 GeV,  1.1073977E+00,  5.48E-04,  3)
PROCESS(`e1,E1', `sb2,sb2c'	,  5000 GeV,  4.5086848E-01,  2.18E-04,  3)
PROCESS(`e1,E1', `ss2,ss2c'	,  5000 GeV,  2.6945558E-01,  1.34E-04,  3)
PROCESS(`e1,E1', `sb1,sb1c'	,  5000 GeV,  3.7348139E-01,  1.85E-04,  3)
PROCESS(`e1,E1', `sd2,sd2c'	,  5000 GeV,  2.6997913E-01,  1.31E-04,  3)
PROCESS(`e1,E1', `sb1,sb2c'	,  5000 GeV,  2.7681447E-01,  1.35E-04,  3)
PROCESS(`e1,E1', `neu1,neu1'	,  5000 GeV,  4.4534250E+00,  3.20E-03,  3)
PROCESS(`e1,E1', `neu1,neu2'	,  5000 GeV,  6.2935770E-01,  4.51E-04,  3)
PROCESS(`e1,E1', `neu1,neu3'	,  5000 GeV,  1.2201185E+00,  8.70E-04,  3)
PROCESS(`e1,E1', `neu1,neu4'	,  5000 GeV,  8.5488270E-02,  5.76E-05,  3)
PROCESS(`e1,E1', `neu1,neu5'	,  5000 GeV,  3.7669965E-02,  2.36E-05,  3)
PROCESS(`e1,E1', `neu2,neu2'	,  5000 GeV,  9.3520991E-02,  6.66E-05,  3)
PROCESS(`e1,E1', `neu2,neu3'	,  5000 GeV,  9.2255931E-01,  6.60E-04,  3)
PROCESS(`e1,E1', `neu2,neu4'	,  5000 GeV,  4.8420385E-02,  2.87E-05,  3)
PROCESS(`e1,E1', `neu2,neu5'	,  5000 GeV,  2.1101591E-02,  1.31E-05,  3)
PROCESS(`e1,E1', `neu3,neu3'	,  5000 GeV,  3.1178635E+00,  2.26E-03,  3)
PROCESS(`e1,E1', `neu3,neu4'	,  5000 GeV,  1.8161434E-01,  1.26E-04,  3)
PROCESS(`e1,E1', `neu3,neu5'	,  5000 GeV,  5.9780766E-02,  3.70E-05,  3)
PROCESS(`e1,E1', `neu4,neu4'	,  5000 GeV,  1.1128203E-03,  4.94E-07,  3)
PROCESS(`e1,E1', `neu4,neu5'	,  5000 GeV,  1.7730710E+00,  1.10E-03,  3)
PROCESS(`e1,E1', `neu5,neu5'	,  5000 GeV,  2.2572280E-05,  1.01E-08,  3)
PROCESS(`e1,E1', `"ch1+","ch1-"',  5000 GeV,  8.4802746E+00,  6.58E-03,  3)
PROCESS(`e1,E1', `"ch2+","ch2-"',  5000 GeV,  4.6006705E+00,  2.45E-03,  3)
PROCESS(`e1,E1', `"ch1+","ch2-"',  5000 GeV,  2.6356797E-01,  1.43E-04,  3)
PROCESS(`e1,E1', `Z,h01'	,  5000 GeV,  4.7787333E-01,  2.36E-04,  3)
PROCESS(`e1,E1', `Z,h02'	,  5000 GeV,  1.0244223E-02,  5.04E-06,  3)
PROCESS(`e1,E1', `Z,h03'	,  5000 GeV,  2.0848498E-07,  1.03E-10,  3)
PROCESS(`e1,E1', `A01,h01'	,  5000 GeV,  2.7917432E-07,  1.38E-10,  3)
PROCESS(`e1,E1', `A01,h02'	,  5000 GeV,  1.0978578E-05,  5.42E-09,  3)
PROCESS(`e1,E1', `A01,h03'	,  5000 GeV,  2.7734579E-03,  1.37E-06,  3)
PROCESS(`e1,E1', `A02,h01'	,  5000 GeV,  3.1827326E-05,  1.57E-08,  3)
PROCESS(`e1,E1', `A02,h02'	,  5000 GeV,  1.2530522E-03,  6.18E-07,  3)
PROCESS(`e1,E1', `A02,h03'	,  5000 GeV,  2.9985928E-01,  1.47E-04,  3)
PROCESS(`e1,E1', `Hp,Hm'	,  5000 GeV,  7.3524903E-01,  3.60E-04,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
