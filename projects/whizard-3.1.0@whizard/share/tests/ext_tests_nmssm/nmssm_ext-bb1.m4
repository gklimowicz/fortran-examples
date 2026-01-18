dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-bb1.m4 > nmssm_ext-bb1.sin
dnl   whizard -r nmssm_ext-bb1.sin
dnl
dnl ---------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-bb1.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_bb1_)
! -----------------------------------------------------------------------
! ----- Passed 5 checks with different random seeds                 -----
! -----------------------------------------------------------------------
model = NMSSM
read_slha ("nmssm.slha")
?vis_history = false
! -----------------------------------------------------------------------
helicity_selection_threshold = 1.E7
helicity_selection_cutoff = 20000
show(real)
! -----------------------------------------------------------------------
iterations = 3:200000
! -----------------------------------------------------------------------
PROCESS(`b,B', `n1,N1'	      ,  3000 GeV,  1.3267831E+00,  9.09E-04,  3)
PROCESS(`b,B', `n2,N2'	      ,  3000 GeV,  1.3293584E+00,  8.88E-04,  3)
PROCESS(`b,B', `n3,N3'	      ,  3000 GeV,  1.3273127E+00,  9.08E-04,  3)
PROCESS(`b,B', `t,T'	      ,  3000 GeV,  5.2144408E+04,  1.48E+02,  3)
PROCESS(`b,B', `Wp,Wm'	      ,  3000 GeV,  3.3903881E+02,  4.20E-01,  3)
PROCESS(`b,B', `Z,Z'	      ,  3000 GeV,  2.9872868E+01,  3.21E-02,  3)
PROCESS(`b,B', `se1,SE1'      ,  3000 GeV,  4.3638043E-01,  3.71E-04,  3)
PROCESS(`b,B', `se2,SE2'      ,  3000 GeV,  1.0634283E-01,  9.06E-05,  3)
PROCESS(`b,B', `smu1,SMU1'    ,  3000 GeV,  4.3661518E-01,  3.70E-04,  3)
PROCESS(`b,B', `smu2,SMU2'    ,  3000 GeV,  1.0632861E-01,  9.00E-05,  3)
PROCESS(`b,B', `stau1,STAU1'  ,  3000 GeV,  1.0801117E-01,  9.13E-05,  3)
PROCESS(`b,B', `stau2,STAU2'  ,  3000 GeV,  1.1187495E-01,  9.47E-05,  3)
PROCESS(`b,B', `stau1,STAU2'  ,  3000 GeV,  1.6168606E-01,  1.36E-04,  3)
PROCESS(`b,B', `sn11,SN11'    ,  3000 GeV,  6.4716889E-01,  5.62E-04,  3)
PROCESS(`b,B', `sn21,SN21'    ,  3000 GeV,  6.4761916E-01,  5.56E-04,  3)
PROCESS(`b,B', `sn31,SN31'    ,  3000 GeV,  6.4714271E-01,  5.58E-04,  3)
PROCESS(`b,B', `su1,su1c'     ,  3000 GeV,  5.5742041E+01,  4.78E-02,  3)
PROCESS(`b,B', `su2,su2c'     ,  3000 GeV,  5.5159550E+01,  4.75E-02,  3)
PROCESS(`b,B', `sc1,sc1c'     ,  3000 GeV,  5.5797214E+01,  4.74E-02,  3)
PROCESS(`b,B', `sc2,sc2c'     ,  3000 GeV,  5.5156147E+01,  4.74E-02,  3)
PROCESS(`b,B', `st1,st1c'     ,  3000 GeV,  1.7348251E+02,  2.06E-01,  3)
PROCESS(`b,B', `st2,st2c'     ,  3000 GeV,  4.3466051E+01,  5.73E-02,  3)
PROCESS(`b,B', `st1,st2c'     ,  3000 GeV,  2.0700610E+01,  2.94E-02,  3)
PROCESS(`b,B', `sd1,sd1c'     ,  3000 GeV,  5.5636486E+01,  4.74E-02,  3)
PROCESS(`b,B', `sd2,sd2c'     ,  3000 GeV,  5.5081261E+01,  4.65E-02,  3)
PROCESS(`b,B', `ss1,ss1c'     ,  3000 GeV,  5.5611371E+01,  4.75E-02,  3)
PROCESS(`b,B', `ss2,ss2c'     ,  3000 GeV,  5.5025180E+01,  4.75E-02,  3)
PROCESS(`b,B', `sb1,sb1c'     ,  3000 GeV,  3.8501685E+02,  3.29E-01,  3)
dnl The WHIZARD 1 cross section seems to be fluctuating up
PROCESS(`b,B', `sb2,sb2c'     ,  3000 GeV,  3.7789934E+02,  3.11E-01,  4)
PROCESS(`b,B', `sb1,sb2c'     ,  3000 GeV,  2.8687609E+02,  2.60E-01,  3)
PROCESS(`b,B', `neu1,neu1'    ,  3000 GeV,  4.3922016E-02,  2.44E-05,  3)
PROCESS(`b,B', `neu1,neu2'    ,  3000 GeV,  2.7778895E-02,  1.30E-05,  3)
PROCESS(`b,B', `neu1,neu3'    ,  3000 GeV,  2.6056866E-01,  1.45E-04,  3)
PROCESS(`b,B', `neu1,neu4'    ,  3000 GeV,  1.7624684E-02,  7.28E-06,  3)
PROCESS(`b,B', `neu1,neu5'    ,  3000 GeV,  5.0275637E-02,  1.40E-05,  3)
PROCESS(`b,B', `neu2,neu2'    ,  3000 GeV,  2.5812694E-02,  1.05E-05,  3)
PROCESS(`b,B', `neu2,neu3'    ,  3000 GeV,  1.4246248E-01,  7.53E-05,  3)
PROCESS(`b,B', `neu2,neu4'    ,  3000 GeV,  1.3964712E-01,  4.78E-05,  3)
PROCESS(`b,B', `neu2,neu5'    ,  3000 GeV,  1.2972948E-01,  5.02E-05,  3)
PROCESS(`b,B', `neu3,neu3'    ,  3000 GeV,  1.3271518E+00,  7.58E-04,  3)
PROCESS(`b,B', `neu3,neu4'    ,  3000 GeV,  2.2735489E-01,  1.02E-04,  3)
PROCESS(`b,B', `neu3,neu5'    ,  3000 GeV,  1.0490707E-01,  2.97E-05,  3)
PROCESS(`b,B', `neu4,neu4'    ,  3000 GeV,  4.8621366E-03,  2.60E-06,  3)
PROCESS(`b,B', `neu4,neu5'    ,  3000 GeV,  2.3416438E+00,  6.18E-04,  3)
PROCESS(`b,B', `neu5,neu5'    ,  3000 GeV,  8.2010808E-03,  3.39E-06,  3)
PROCESS(`b,B', `"ch1+","ch1-"',  3000 GeV,  5.5580555E+00,  5.75E-03,  3)
PROCESS(`b,B', `"ch2+","ch2-"',  3000 GeV,  4.8828820E+01,  3.62E-02,  3)
PROCESS(`b,B', `"ch1+","ch2-"',  3000 GeV,  1.7992620E+00,  1.07E-03,  3)
PROCESS(`b,B', `sgl,sgl'      ,  3000 GeV,  6.5001243E+02,  3.54E-01,  3)
PROCESS(`b,B', `Z,h01'        ,  3000 GeV,  7.3219932E-01,  1.13E-03,  3)
PROCESS(`b,B', `Z,h02'        ,  3000 GeV,  1.8317299E-02,  3.18E-05,  3)
PROCESS(`b,B', `Z,h03'        ,  3000 GeV,  2.6999897E-01,  3.76E-04,  3)
PROCESS(`b,B', `A01,h01'      ,  3000 GeV,  9.1885371E-03,  1.39E-05,  3)
PROCESS(`b,B', `A01,h02'      ,  3000 GeV,  9.5026501E-04,  1.46E-06,  3)
PROCESS(`b,B', `A01,h03'      ,  3000 GeV,  3.4858710E-03,  6.06E-06,  3)
PROCESS(`b,B', `A02,h01'      ,  3000 GeV,  1.8883303E-03,  2.54E-06,  3)
PROCESS(`b,B', `A02,h02'      ,  3000 GeV,  1.8364076E-02,  2.84E-05,  3)
PROCESS(`b,B', `A02,h03'      ,  3000 GeV,  7.2089059E-02,  1.37E-04,  3)
PROCESS(`b,B', `Hp,Hm'        ,  3000 GeV,  1.0213147E+00,  1.37E-03,  3)
PROCESS(`b,B', `Z,A01'        ,  3000 GeV,  1.2044804E-02,  1.13E-05,  3)
PROCESS(`b,B', `Z,A02'        ,  3000 GeV,  2.6926666E-01,  3.72E-04,  3)
PROCESS(`b,B', `h01,h01'      ,  3000 GeV,  1.2010048E-04,  1.57E-07,  3)
PROCESS(`b,B', `h01,h02'      ,  3000 GeV,  4.1045455E-03,  6.18E-06,  3)
PROCESS(`b,B', `h01,h03'      ,  3000 GeV,  1.3237779E-03,  1.67E-06,  3)
PROCESS(`b,B', `h02,h02'      ,  3000 GeV,  9.4370065E-04,  1.42E-06,  3)
PROCESS(`b,B', `h02,h03'      ,  3000 GeV,  1.7181224E-02,  2.64E-05,  3)
PROCESS(`b,B', `h03,h03'      ,  3000 GeV,  1.5900755E-04,  2.40E-07,  3)
PROCESS(`b,B', `A01,A01'      ,  3000 GeV,  3.7700752E-06,  5.02E-09,  3)
PROCESS(`b,B', `A01,A02'      ,  3000 GeV,  1.3715377E-03,  2.07E-06,  3)
PROCESS(`b,B', `A02,A02'      ,  3000 GeV,  3.5838597E-05,  5.76E-08,  3)
! ----------------------------------------------------------------------
END_TESTSUITE
