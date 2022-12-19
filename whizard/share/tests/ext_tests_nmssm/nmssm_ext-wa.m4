dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-wa.m4 > nmssm_ext-wa.sin
dnl   whizard -r nmssm_ext-wa.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-wa.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_wa_)
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
PROCESS(`Wm,A', `Z,Wm'		,  3000 GeV,  1.1853089E+05,  3.22E+02,  3)
PROCESS(`Wm,A', `se1,SN11'	,  3000 GeV,  7.1292368E+00,  8.76E-03,  3)
PROCESS(`Wm,A', `smu1,SN21'	,  3000 GeV,  7.1170086E+00,  8.88E-03,  3)
PROCESS(`Wm,A', `stau1,SN31'	,  3000 GeV,  3.5287938E+00,  4.47E-03,  3)
PROCESS(`Wm,A', `stau2,SN31'	,  3000 GeV,  3.6216916E+00,  4.43E-03,  3)
PROCESS(`Wm,A', `sd1,su1c'	,  3000 GeV,  2.2987199E+00,  2.57E-03,  3)
PROCESS(`Wm,A', `ss1,sc1c'	,  3000 GeV,  2.2974001E+00,  2.58E-03,  3)
PROCESS(`Wm,A', `sb1,st1c'	,  3000 GeV,  4.5548471E+00,  5.77E-03,  3)
PROCESS(`Wm,A', `sb2,st2c'	,  3000 GeV,  1.9043009E+00,  2.22E-03,  3)
PROCESS(`Wm,A', `sb1,st2c'	,  3000 GeV,  1.7018845E+00,  1.98E-03,  3)
PROCESS(`Wm,A', `sb2,st1c'	,  3000 GeV,  5.9037120E+00,  7.52E-03,  3)
PROCESS(`Wm,A', `"ch1-",neu1'	,  3000 GeV,  4.6487827E+00,  6.50E-03,  3)
PROCESS(`Wm,A', `"ch1-",neu2'	,  3000 GeV,  1.4741834E+01,  2.08E-02,  3)
PROCESS(`Wm,A', `"ch1-",neu3'	,  3000 GeV,  1.9765262E+02,  2.76E-01,  3)
PROCESS(`Wm,A', `"ch1-",neu4'	,  3000 GeV,  2.9964591E+01,  3.79E-02,  3)
PROCESS(`Wm,A', `"ch1-",neu5'	,  3000 GeV,  3.1435849E+01,  4.07E-02,  3)
PROCESS(`Wm,A', `"ch2-",neu1'	,  3000 GeV,  3.2625723E+00,  3.32E-03,  3)
PROCESS(`Wm,A', `"ch2-",neu2'	,  3000 GeV,  3.0312131E+01,  3.03E-02,  3)
PROCESS(`Wm,A', `"ch2-",neu3'	,  3000 GeV,  1.6451678E+01,  1.70E-02,  3)
PROCESS(`Wm,A', `"ch2-",neu4'	,  3000 GeV,  2.5181546E+01,  2.53E-02,  3)
PROCESS(`Wm,A', `"ch2-",neu5'	,  3000 GeV,  2.2475369E+01,  2.24E-02,  3)
PROCESS(`Wm,A', `Hm,h01'	,  3000 GeV,  8.5704569E-03,  8.89E-06,  3)
PROCESS(`Wm,A', `Hm,h02'	,  3000 GeV,  3.3508330E-01,  3.49E-04,  3)
PROCESS(`Wm,A', `Hm,h03'	,  3000 GeV,  1.1566947E+00,  4.67E-04,  3)
PROCESS(`Wm,A', `Hm,A01'	,  3000 GeV,  7.6103816E-01,  7.86E-04,  3)
PROCESS(`Wm,A', `Hm,A02'	,  3000 GeV,  1.1513994E+00,  4.56E-04,  3)
PROCESS(`Wm,A', `Wm,h01'	,  3000 GeV,  1.6181633E+04,  4.42E+01,  3)
PROCESS(`Wm,A', `Wm,h02'	,  3000 GeV,  3.4644355E+02,  9.43E-01,  3)
PROCESS(`Wm,A', `Wm,h03'	,  3000 GeV,  7.6833179E-03,  1.92E-05,  3)
PROCESS(`Wm,A', `Z,Wm'	 	,  5000 GeV,  1.1915254E+05,  1.91E+02,  3)
PROCESS(`Wm,A', `se1,SN11'	,  5000 GeV,  2.7194113E+00,  2.06E-03,  3)
PROCESS(`Wm,A', `smu1,SN21'	,  5000 GeV,  2.7198052E+00,  2.04E-03,  3)
PROCESS(`Wm,A', `stau1,SN31'	,  5000 GeV,  1.3434103E+00,  1.03E-03,  3)
PROCESS(`Wm,A', `stau2,SN31'	,  5000 GeV,  1.3810638E+00,  1.03E-03,  3)
PROCESS(`Wm,A', `sd1,su1c'	,  5000 GeV,  1.4030154E+00,  8.74E-04,  3)
PROCESS(`Wm,A', `ss1,sc1c'	,  5000 GeV,  1.4027770E+00,  8.74E-04,  3)
PROCESS(`Wm,A', `sb1,st1c'	,  5000 GeV,  1.9177728E+00,  1.47E-03,  3)
PROCESS(`Wm,A', `sb2,st2c'	,  5000 GeV,  1.6145271E+00,  1.21E-03,  3)
PROCESS(`Wm,A', `sb1,st2c'	,  5000 GeV,  1.4329251E+00,  1.08E-03,  3)
PROCESS(`Wm,A', `sb2,st1c'	,  5000 GeV,  2.4901742E+00,  1.91E-03,  3)
PROCESS(`Wm,A', `"ch1-",neu1'	,  5000 GeV,  2.0334158E+00,  1.74E-03,  3)
PROCESS(`Wm,A', `"ch1-",neu2'	,  5000 GeV,  6.4478331E+00,  5.63E-03,  3)
PROCESS(`Wm,A', `"ch1-",neu3'	,  5000 GeV,  8.6281156E+01,  7.66E-02,  3)
PROCESS(`Wm,A', `"ch1-",neu4'	,  5000 GeV,  1.3516987E+01,  1.10E-02,  3)
PROCESS(`Wm,A', `"ch1-",neu5'	,  5000 GeV,  1.4266261E+01,  1.15E-02,  3)
PROCESS(`Wm,A', `"ch2-",neu1'	,  5000 GeV,  1.7136337E+00,  1.09E-03,  3)
PROCESS(`Wm,A', `"ch2-",neu2'	,  5000 GeV,  1.5702365E+01,  9.96E-03,  3)
PROCESS(`Wm,A', `"ch2-",neu3'	,  5000 GeV,  8.7299458E+00,  5.59E-03,  3)
PROCESS(`Wm,A', `"ch2-",neu4'	,  5000 GeV,  1.3166412E+01,  8.69E-03,  3)
PROCESS(`Wm,A', `"ch2-",neu5'	,  5000 GeV,  1.1743778E+01,  7.72E-03,  3)
PROCESS(`Wm,A', `Hm,h01'	,  5000 GeV,  3.1698409E-03,  2.03E-06,  3)
PROCESS(`Wm,A', `Hm,h02'	,  5000 GeV,  1.2417811E-01,  7.94E-05,  3)
PROCESS(`Wm,A', `Hm,h03'	,  5000 GeV,  6.4781044E-01,  3.44E-04,  3)
PROCESS(`Wm,A', `Hm,A01'	,  5000 GeV,  2.7977773E-01,  1.78E-04,  3)
PROCESS(`Wm,A', `Hm,A02'	,  5000 GeV,  6.4532326E-01,  3.39E-04,  3)
PROCESS(`Wm,A', `Wm,h01'	,  5000 GeV,  1.6275348E+04,  2.58E+01,  3)
PROCESS(`Wm,A', `Wm,h02'	,  5000 GeV,  3.5003047E+02,  5.52E-01,  3)
PROCESS(`Wm,A', `Wm,h03'	,  5000 GeV,  8.2874985E-03,  1.29E-05,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
