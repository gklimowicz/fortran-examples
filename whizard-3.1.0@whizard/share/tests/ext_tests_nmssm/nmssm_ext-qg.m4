dnl Process this with
dnl
dnl   m4 testuite.m4 nmssm_ext-qg.m4 > nmssm_ext-qg.sin
dnl   whizard -r nmssm_ext-qg.sin
dnl
dnl -----------------------------------------------------------------------
BEGIN_TESTSUITE($Id: nmssm_ext-qg.m4 $, nmssm_ext_qg_)
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
PROCESS(`u,gl', `su1,sgl'	,  3000 GeV,  1.8992459E+03,  1.29E+00,  3)
PROCESS(`u,gl', `su2,sgl'	,  3000 GeV,  1.9003712E+03,  1.27E+00,  3)
PROCESS(`U,gl', `su1c,sgl'	,  3000 GeV,  1.8990060E+03,  1.31E+00,  3)
PROCESS(`U,gl', `su2c,sgl'	,  3000 GeV,  1.8991037E+03,  1.29E+00,  3)
PROCESS(`d,gl', `sd1,sgl'	,  3000 GeV,  1.8976483E+03,  1.30E+00,  3)
PROCESS(`d,gl', `sd2,sgl'	,  3000 GeV,  1.9001410E+03,  1.28E+00,  3)
PROCESS(`D,gl', `sd1c,sgl'	,  3000 GeV,  1.8966650E+03,  1.31E+00,  3)
PROCESS(`D,gl', `sd2c,sgl'	,  3000 GeV,  1.9007399E+03,  1.29E+00,  3)
PROCESS(`u,gl', `su1,sgl'	,  5000 GeV,  1.0844597E+03,  9.65E-01,  3)
PROCESS(`u,gl', `su2,sgl'	,  5000 GeV,  1.0858668E+03,  9.42E-01,  3)
PROCESS(`U,gl', `su1c,sgl'	,  5000 GeV,  1.0848775E+03,  9.55E-01,  3)
PROCESS(`U,gl', `su2c,sgl'	,  5000 GeV,  1.0855153E+03,  9.54E-01,  3)
PROCESS(`d,gl', `sd1,sgl'	,  5000 GeV,  1.0863452E+03,  9.39E-01,  3)
PROCESS(`d,gl', `sd2,sgl'	,  5000 GeV,  1.0860693E+03,  9.48E-01,  3)
PROCESS(`D,gl', `sd1c,sgl'	,  5000 GeV,  1.0858148E+03,  9.49E-01,  3)
PROCESS(`D,gl', `sd2c,sgl'	,  5000 GeV,  1.0876225E+03,  9.31E-01,  3)
! -------------------------------------------------------------------------
	END_TESTSUITE
