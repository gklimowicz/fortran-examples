dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-uu1.m4 > nmssm_ext-uu1.sin
dnl   whizard -r nmssm_ext-uu1
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-uu1.m4 1758 2010-02-11 13:27:04Z ohl $, nmssm_ext_uu1_)
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
PROCESS(`u,U', `n1,N1'		,  3000 GeV,  1.0295549E+00,  2.79E-04,  3)
PROCESS(`u,U', `n2,N2'		,  3000 GeV,  1.0296720E+00,  2.81E-04,  3)
PROCESS(`u,U', `n3,N3'		,  3000 GeV,  1.0294710E+00,  2.79E-04,  3)
PROCESS(`u,U', `t,T'		,  3000 GeV,  5.5736643E+02,  9.40E-02,  3)
PROCESS(`u,U', `Wp,Wm'		,  3000 GeV,  1.5109945E+02,  3.08E-01,  3)
PROCESS(`u,U', `Z,Z'		,  3000 GeV,  1.5163758E+01,  6.37E-03,  3)
PROCESS(`u,U', `se1,SE1'	,  3000 GeV,  7.0902694E-01,  6.06E-04,  3)
PROCESS(`u,U', `se2,SE2'	,  3000 GeV,  3.6118292E-01,  3.07E-04,  3)
PROCESS(`u,U', `smu1,SMU1'	,  3000 GeV,  7.0889908E-01,  6.02E-04,  3)
PROCESS(`u,U', `smu2,SMU2'	,  3000 GeV,  3.6049419E-01,  3.12E-04,  3)
PROCESS(`u,U', `stau1,STAU1'	,  3000 GeV,  4.0963673E-01,  3.44E-04,  3)
PROCESS(`u,U', `stau2,STAU2'	,  3000 GeV,  4.0964995E-01,  3.57E-04,  3)
PROCESS(`u,U', `stau1,STAU2'	,  3000 GeV,  1.2520584E-01,  1.05E-04,  3)
PROCESS(`u,U', `sn11,SN11'	,  3000 GeV,  5.0157143E-01,  4.28E-04,  3)
PROCESS(`u,U', `sn21,SN21'	,  3000 GeV,  5.0092094E-01,  4.39E-04,  3)
PROCESS(`u,U', `sn31,SN31'	,  3000 GeV,  5.0127891E-01,  4.31E-04,  3)
PROCESS(`u,U', `su1,su1c'	,  3000 GeV,  3.8744139E+02,  4.19E-01,  3)
PROCESS(`u,U', `su2,su2c'	,  3000 GeV,  4.0183903E+02,  4.39E-01,  3)
PROCESS(`u,U', `su1,su2c'	,  3000 GeV,  2.6821989E+02,  2.79E-01,  3)
PROCESS(`u,U', `sc1,sc1c'	,  3000 GeV,  5.5842665E+01,  4.78E-02,  3)
PROCESS(`u,U', `sc2,sc2c'	,  3000 GeV,  5.5327820E+01,  4.73E-02,  3)
PROCESS(`u,U', `st1,st1c'	,  3000 GeV,  8.4683885E+01,  7.09E-02,  3)
PROCESS(`u,U', `st2,st2c'	,  3000 GeV,  2.8244460E+01,  2.43E-02,  3)
PROCESS(`u,U', `st1,st2c'	,  3000 GeV,  1.5337503E-01,  1.31E-04,  3)
PROCESS(`u,U', `sd1,sd1c'	,  3000 GeV,  5.4688073E+01,  6.60E-02,  3)
PROCESS(`u,U', `sd2,sd2c'	,  3000 GeV,  5.5149023E+01,  7.69E-02,  3)
PROCESS(`u,U', `ss1,ss1c'	,  3000 GeV,  5.5478563E+01,  4.81E-02,  3)
PROCESS(`u,U', `ss2,ss2c'	,  3000 GeV,  5.5051000E+01,  4.72E-02,  3)
PROCESS(`u,U', `sb1,sb1c'	,  3000 GeV,  5.5706205E+01,  4.79E-02,  3)
PROCESS(`u,U', `sb2,sb2c'	,  3000 GeV,  5.4500326E+01,  4.70E-02,  3)
PROCESS(`u,U', `sb1,sb2c'	,  3000 GeV,  1.5216858E-01,  1.29E-04,  3)
PROCESS(`u,U', `neu1,neu1'	,  3000 GeV,  4.6905088E-01,  2.65E-04,  3)
PROCESS(`u,U', `neu1,neu2'	,  3000 GeV,  6.4045330E-02,  3.42E-05,  3)
PROCESS(`u,U', `neu1,neu3'	,  3000 GeV,  1.6539942E-02,  1.01E-05,  3)
PROCESS(`u,U', `neu1,neu4'	,  3000 GeV,  7.8552408E-03,  4.25E-06,  3)
PROCESS(`u,U', `neu1,neu5'	,  3000 GeV,  2.5230575E-02,  6.64E-06,  3)
PROCESS(`u,U', `neu2,neu2'	,  3000 GeV,  7.2574068E-03,  5.51E-06,  3)
PROCESS(`u,U', `neu2,neu3'	,  3000 GeV,  2.5140902E-01,  1.47E-04,  3)
PROCESS(`u,U', `neu2,neu4'	,  3000 GeV,  8.9064573E-02,  2.81E-05,  3)
PROCESS(`u,U', `neu2,neu5'	,  3000 GeV,  1.5866250E-02,  3.96E-06,  3)
PROCESS(`u,U', `neu3,neu3'	,  3000 GeV,  1.6021020E+00,  9.04E-04,  3)
PROCESS(`u,U', `neu3,neu4'	,  3000 GeV,  4.7601740E-02,  3.76E-05,  3)
PROCESS(`u,U', `neu3,neu5'	,  3000 GeV,  2.4267001E-02,  5.60E-06,  3)
PROCESS(`u,U', `neu4,neu4'	,  3000 GeV,  1.5711764E-03,  1.06E-06,  3)
PROCESS(`u,U', `neu4,neu5'	,  3000 GeV,  1.9131520E+00,  5.37E-04,  3)
PROCESS(`u,U', `neu5,neu5'	,  3000 GeV,  4.4289624E-05,  1.38E-08,  3)
PROCESS(`u,U', `"ch1+","ch1-"'	,  3000 GeV,  5.0006803E+00,  4.98E-03,  3)
PROCESS(`u,U', `"ch2+","ch2-"'	,  3000 GeV,  2.7003611E+00,  5.82E-04,  3)
PROCESS(`u,U', `"ch1+","ch2-"'	,  3000 GeV,  1.9853248E-01,  1.28E-04,  3)
PROCESS(`u,U', `sgl,sgl'	,  3000 GeV,  6.5039937E+02,  3.48E-01,  3)
PROCESS(`u,U', `Z,h01'		,  3000 GeV,  5.0433542E-01,  4.27E-04,  3)
PROCESS(`u,U', `Z,h02'		,  3000 GeV,  1.0812240E-02,  9.05E-06,  3)
PROCESS(`u,U', `Z,h03'		,  3000 GeV,  1.4609808E-07,  1.22E-10,  3)
PROCESS(`u,U', `A01,h01'	,  3000 GeV,  2.9261777E-07,  2.51E-10,  3)
PROCESS(`u,U', `A01,h02'	,  3000 GeV,  1.1518655E-05,  9.82E-09,  3)
PROCESS(`u,U', `A01,h03'	,  3000 GeV,  1.9252619E-03,  1.64E-06,  3)
PROCESS(`u,U', `A02,h01'	,  3000 GeV,  2.1873214E-05,  1.86E-08,  3)
PROCESS(`u,U', `A02,h02'	,  3000 GeV,  8.5942209E-04,  7.31E-07,  3)
PROCESS(`u,U', `A02,h03'	,  3000 GeV,  6.1103851E-02,  5.23E-05,  3)
PROCESS(`u,U', `Hp,Hm'		,  3000 GeV,  8.7753658E-02,  7.59E-05,  3)
! -------------------------------------------------------------------------
END_TESTSUITE
