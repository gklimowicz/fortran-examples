#include "rundeck_opts.h"
!------------------------------------------------------------------------------
module AmpTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  AMPTracersMetadata_mod encapsulates the TRACERS_AMP metadata
!@auth NCCS ASTG
  use sharedTracersMetadata_mod, only: DMS_setspec, &
    SO2_setspec, H2O2_s_setspec, NH3_setspec
  USE CONSTANT, only: mwat
  USE AERO_PARAM, only: &
    SOLU_DD1, SOLU_DD2, SOLU_AKK, SOLU_ACC, &
    SOLU_DS1, SOLU_DS2, SOLU_SSA, SOLU_SSC, &
    SOLU_SSS, SOLU_OCC, SOLU_BC1, SOLU_BC2, SOLU_BC3, &
    SOLU_DBC, SOLU_BOC, SOLU_BCS, SOLU_OCS, SOLU_MXX
  USE AERO_PARAM, only: &
    DG_DD1, DG_DD2, DG_AKK,                    &
    DG_DS1, DG_DS2, DG_SSA, DG_SSC, DG_ACC,    &
    DG_SSS, DG_OCC, DG_BC1, DG_BC2, DG_BC3,    &
    DG_DBC, DG_BOC, DG_BCS, DG_OCS, DG_MXX
  USE AERO_ACTV, only: &
    DENS_SULF, DENS_DUST, DENS_SEAS, DENS_BCAR, DENS_OCAR
  use Tracer_com, only: &
    n_M_NO3,   n_M_NH4,   n_M_H2O,   n_M_AKK_SU, &
    n_N_AKK_1, n_M_ACC_SU,n_N_ACC_1, n_M_DD1_SU, &
    n_M_DD1_DU,n_N_DD1_1, n_M_DS1_SU,n_M_DS1_DU, &
    n_N_DS1_1 ,n_M_DD2_SU,n_M_DD2_DU,n_N_DD2_1 , &
    n_M_DS2_SU,n_M_DS2_DU,n_N_DS2_1 ,n_M_SSA_SU, &
    n_M_SSA_SS,n_M_SSC_SS,                       &
    n_M_OCC_SU,n_M_OCC_OC,n_N_OCC_1 ,            &
    n_M_BC1_SU,n_M_BC1_BC,n_N_BC1_1 ,n_M_BC2_SU, &
    n_M_BC2_BC,n_N_BC2_1 ,n_M_BC3_SU,n_M_BC3_BC, &
    n_N_BC3_1 ,n_M_DBC_SU,n_M_DBC_BC,n_M_DBC_DU, &
    n_N_DBC_1 ,n_M_BOC_SU,n_M_BOC_BC,n_M_BOC_OC, &
    n_N_BOC_1, n_M_BCS_SU,n_M_BCS_BC,n_N_BCS_1 , &
    n_M_MXX_SU,n_M_MXX_BC,n_M_MXX_OC,n_M_MXX_DU, &
    n_M_MXX_SS,n_N_MXX_1 ,n_M_OCS_SU,n_M_OCS_OC, &
    n_N_OCS_1,n_M_SSS_SS,n_M_SSS_SU,             &
    n_H2SO4, n_N_SSA_1, n_N_SSC_1
  use RunTimeControls_mod, only: &
    tracers_nitrate, tracers_aerosols_koch, tracers_aerosols_seasalt, &
    tracers_amp_m1, tracers_amp_m2,         &
    tracers_amp_m3, tracers_amp_m4,         &
    tracers_amp_m5, tracers_amp_m6,         &
    tracers_amp_m7, tracers_amp_m8,         &
    tracers_special_shindell
  use Tracer_com, only: ntmAMPi, ntmAMPe, ntmAMP, ntm_chem, coupled_chem
  use OldTracer_mod, only: set_needtrs
  use OldTracer_mod, only: nPart
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: set_has_chemistry
  use OldTracer_mod, only: oldAddTracer
  use Tracer_mod, only: Tracer

  implicit none
  private

  public AMP_initMetadata

#ifdef TRACERS_AMP_M1
  integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,1 ,1,  & !AKK
    2 ,2 ,3 ,3 ,3,  & !ACC,DD1
    4 ,4 ,4 ,5 ,5,  & !DS1,DD2
    5 ,6 ,6 ,6 ,7,  & !DD2,DS2,SSA
    7 ,7, 8, 8,     & !SSA,SSC
    9 ,9 ,9 ,10,10, & !OCC,BC1
    10,11,11,11,12, & !BC1,BC2,BC3
    12,12,13,13,13, & !BC3,DBC
    13,14,14,14,14, & !DBC,BOC
    15,15,15,16,16, & !BCS,MXX
    16,16,16,16/)
  integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,0 ,1,  & !AKK
    0 ,2 ,0 ,0 ,3,  & !ACC,DD1
    0 ,0 ,4 ,0 ,0,  & !DS1,DD2
    5 ,0 ,0 ,6 ,0,  & !DD2,DS2,SSA
    0 ,7 ,0 ,8 ,    & !SSA,SSC
    0 ,0 ,9 ,0 ,0 , & !OCC,BC1
    10,0 ,0 ,11,0 , & !BC1,BC2,BC3
    0 ,12,0 ,0 ,0 , & !BC3,DBC
    13,0 ,0 ,0 ,14, & !DBC,BOC
    0 ,0 ,15,0 ,0 , & !BCS,MXX
    0 ,0 , 0,16/)
  integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/ &
    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10, &
    11,12,13,14,15,16,17,18,19,20, &
    21,22,   24,25,26,27,28,29,30, &
    31,32,33,34,35,36,37,38,39,40, &
    41,42,43,44,45,46,47,48,49,50, &
    51,52,53,54   /)
  integer, parameter :: AMP_trm_nm1(ntmAMP)=(/ &
    0 ,         0 ,         0 ,ntm_chem+4 ,ntm_chem+4,  & !AKK
    ntm_chem+6 ,ntm_chem+6 ,ntm_chem+8 ,ntm_chem+8 ,ntm_chem+8,  & !ACC,DD1
    ntm_chem+11,ntm_chem+11,ntm_chem+11,ntm_chem+14,ntm_chem+14, & !DS1,DD2
    ntm_chem+14,ntm_chem+17,ntm_chem+17,ntm_chem+17,ntm_chem+20, & !DD2,DS2,SSA
    ntm_chem+20,ntm_chem+20,ntm_chem+23,ntm_chem+23,             & !SSA,SSC
    ntm_chem+25,ntm_chem+25,ntm_chem+25,ntm_chem+28,ntm_chem+28, & !OCC,BC1
    ntm_chem+28,ntm_chem+31,ntm_chem+31,ntm_chem+31,ntm_chem+34, & !BC1,BC2,BC3
    ntm_chem+34,ntm_chem+34,ntm_chem+37,ntm_chem+37,ntm_chem+37, & !BC3,DBC
    ntm_chem+37,ntm_chem+41,ntm_chem+41,ntm_chem+41,ntm_chem+41, & !DBC,BOC
    ntm_chem+45,ntm_chem+45,ntm_chem+45,ntm_chem+48,ntm_chem+48, & !BCS,MXX
    ntm_chem+48,ntm_chem+48,ntm_chem+48,ntm_chem+48/)
  integer, parameter :: AMP_trm_nm2(ntmAMP)=(/ &
    0 ,         0 ,         0 ,ntm_chem+4 ,ntm_chem+4,  & !AKK
    ntm_chem+6 ,ntm_chem+6 ,ntm_chem+9 ,ntm_chem+9 ,ntm_chem+9,  & !ACC,DD1
    ntm_chem+12,ntm_chem+12,ntm_chem+12,ntm_chem+15,ntm_chem+15, & !DS1,DD2
    ntm_chem+15,ntm_chem+18,ntm_chem+18,ntm_chem+18,ntm_chem+21, & !DD2,DS2,SSA
    ntm_chem+21,ntm_chem+21,ntm_chem+23,ntm_chem+23,             & !SSA,SSC
    ntm_chem+26,ntm_chem+26,ntm_chem+26,ntm_chem+29,ntm_chem+29, & !OCC,BC1
    ntm_chem+29,ntm_chem+32,ntm_chem+32,ntm_chem+32,ntm_chem+35, & !BC1,BC2,BC3
    ntm_chem+35,ntm_chem+35,ntm_chem+39,ntm_chem+39,ntm_chem+39, & !BC3,DBC
    ntm_chem+39,ntm_chem+43,ntm_chem+43,ntm_chem+43,ntm_chem+43, & !DBC,BOC
    ntm_chem+46,ntm_chem+46,ntm_chem+46,ntm_chem+52,ntm_chem+52, & !BCS,MXX
    ntm_chem+52,ntm_chem+52,ntm_chem+52,ntm_chem+52/)
#endif
#ifdef TRACERS_AMP_M2
  integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,1 ,1,  & !AKK
    2 ,2 ,3 ,3 ,3,  & !ACC,DD1
    4 ,4 ,4 ,5 ,5,  & !DS1,DD2
    5 ,6 ,6 ,6 ,7,  & !DD2,DS2,SSA
    7 ,8,           & !SSA,SSC
    9 ,9 ,9 ,10,10, & !OCC,BC1
    10,11,11,11,12, & !BC1,BC2,OSC
    12,12,13,13,13, & !BC3,DBC
    13,14,14,14,14, & !DBC,BOC
    15,15,15,16,16, & !BCS,MXX
    16,16,16,16/)
  integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,0 ,1,  & !AKK
    0 ,2 ,0 ,0 ,3,  & !ACC,DD1
    0 ,0 ,4 ,0 ,0,  & !DS1,DD2
    5 ,0 ,0 ,6 ,0,  & !DD2,DS2,SSA
    0 ,0,           & !SSA,SSC
    0 ,0 ,9 ,0 , 0, & !OCC,BC1
    10,0 ,0 ,11,0 , & !BC1,BC2,OSC
    0 ,12,0 ,0 ,0 , & !DBC
    13,0 ,0 ,0 ,14, & !DBC,BOC
    0 ,0 ,15,0 ,0 , & !BCS,MXX
    0 ,0 ,0 ,16/)
  integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/ &
    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10, &
    11,12,13,14,15,16,17,18,19,20, &
    21,      24,  26,27,28,29,30, &
    31,32,33,34,35,36,37,38,39,40, &
    41,42,43,44,45,46,47,48,49,50, &
    51,52,53,54   /)

  integer, parameter :: AMP_trm_nm1(ntmAMP)=(/ &
    0 ,0 ,0 ,4 ,4,  & !AKK
    6 ,6 ,8 ,8 ,8,  & !ACC,DD1
    11,11,11,14,14, & !DS1,DD2
    14,17,17,17,20, & !DD2,DS2,SSA
    20,22,          & !SSA,SSC
    23,23,23,26,26, & !OCC,BC1
    26,29,29,29,32, & !BC1,BC2
    32,32,35,35,35, & !OSC,DBC
    35,39,39,39,39, & !DBC,BOC
    43,43,43,46,46, & !BCS,MXX
    46,46,46,46/)
  integer, parameter :: AMP_trm_nm2(ntmAMP)=(/ &
    0 ,0 ,0 ,4 ,4,  & !AKK
    6 ,6 ,9 ,9 ,9,  & !ACC,DD1
    12,12,12,15,15, & !DS1,DD2
    15,18,18,18,21, & !DD2,DS2,SSA
    21,22,          & !SSA,SSC
    24,24,24,27,27, & !OCC,BC1
    27,30,30,30,33, & !BC1,BC2
    33,33,37,37,37, & !OCS,DBC
    37,41,41,41,41, & !DBC,BOC
    44,44,44,50,50, & !BCS,MXX
    50,50,50,50/)
#endif
#ifdef TRACERS_AMP_M3
  integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,1 ,1,  & !AKK
    2 ,2 ,3 ,3 ,3,  & !ACC,DD1
    4 ,4 ,4 ,5 ,5,  & !DS1,DD2
    5 ,6 ,6 ,6 ,7,  & !DD2,DS2,SSA
    7 ,8,           & !SSA,SSC
    9 ,9 ,9 ,10,10, & !OCC,BC1
    10,11,11,11,12, & !BC1,BC2,BOc
    12,12,12,13,13, & !BOC,MXX
    13,13,13,13/)
  integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,0 ,1,  & !AKK
    0 ,2 ,0 ,0 ,3,  & !ACC,DD1
    0 ,0 ,4 ,0 ,0,  & !DS1,DD2
    5 ,0 ,0 ,6 ,0,  & !DD2,DS2,SSA
    0 ,0,           & !SSA,SSC
    0 ,0 ,9 ,0 ,0 , & !OCC,BC1
    10,0 ,0 ,11,0 , & !BC1,BC2,BOc
    0 ,0 ,12,0 ,0 , & !BOC,MXX
    0 ,0 ,0 ,13/)
  integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/ &
    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10, &
    11,12,13,14,15      ,18   ,20, &
    21,22,23,24,25,26,27,28,29,30, &
    31,32,33,34,35,36,37,38,39,40, &
    41,42,43,44  /)

  integer, parameter :: AMP_trm_nm1(ntmAMP)=(/ &
    0 ,0 ,0 ,4 ,4,  & !AKK
    6 ,6 ,8 ,8 ,8,  & !ACC,DD1
    11,11,11,14,14, & !DS1,DD2
    14,17,17,17,20, & !DD2,DS2,SSA
    20,22,          & !SSA,SSC
    23,23,23,26,26, & !OCC,BC1
    26,29,29,29,32, & !BC1,BC2,BOC
    32,32,32,36,36, & !BOC,MXX
    36,36,36,36/)

  integer, parameter :: AMP_trm_nm2(ntmAMP)=(/ &
    0 ,0 ,0 ,4 ,4,  & !AKK
    6 ,6 ,9 ,9 ,9,  & !ACC,DD1
    12,12,12,15,15, & !DS1,DD2
    15,18,18,18,21, & !DD2,DS2,SSA
    21,22,          & !SSA,SSC
    24,24,24,27,27, & !OCC,BC1
    27,30,30,30,34, & !BC1,BC2,BOC
    34,34,34,40,40, & !BOC,MXX
    40,40,40,40/)

#endif
#ifdef TRACERS_AMP_M4
  integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,1 ,1,  & !ACC
    2 ,2 ,2 ,3 ,3,  & !DD1,DS1
    3 ,4 ,4 ,4 ,5,  & !DS1,DD2,DS2
    5 ,5 ,6 ,6,     & !DS2,SSS
    7 ,7 ,7 ,8 ,8,  & !OCC,BC1
    8 ,9 ,9 ,9 ,10, & !BC1,BC2,MXX
    10,10,10,10,10  & !MXX
    /)
  integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,0 ,1,  & !ACC
    0 ,0 ,2 ,0 ,0,  & !DD1,DS1
    3 ,0 ,0 ,4 ,0,  & !DS1,DD2,DS2
    0 ,5 ,0 ,0,     & !DS2,SSS
    0 ,0 ,6 ,0 ,0,  & !OCC,BC1
    7 ,0 ,0 ,8 ,0,  & !BC1,BC2,MXX
    0, 0, 0 ,0 ,9  & !MXX
    /)
  integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/ &
    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10, &
    11,12,13,14,15,16,17,18,19, &
    21,22,23,24,25,26,27,28,29,30, &
    31,32,33,34,35   /)
  integer, parameter :: AMP_trm_nm1(ntmAMP)=(/ &
    0 ,0 ,0 ,4 ,4,  & !ACC
    6 ,6 ,6 ,9 ,9,  & !DD1
    9 ,12,12,12,15, & !DS1,DD2
    15,15,18,18,    & !DD2,DS2,SSA
    20,20,20,23,23, & !SSA,SSC
    23,26,26,26,29, & !OCC,BC1
    29,29,29,29,29/)

  integer, parameter :: AMP_trm_nm2(ntmAMP)=(/ &
    0 ,0 ,0 ,4 ,4,  & !ACC
    7 ,7 ,7 ,10,10,  & !DD1
    10,13,13,13,16, & !DS1,DD2
    16,16,19,19,    & !DD2,DS2,SSA
    21,21,21,24,24, & !SSA,SSC
    24,27,27,27,33, & !OCC,BC1
    33,33,33,33,33/)

#endif
#ifdef TRACERS_AMP_M5
  integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,1 ,1,  & !AKK
    2 ,2 ,3 ,3 ,3,  & !ACC,DD1
    4 ,4 ,4 ,       & !DS1
    5 ,             & !SSA
    5 ,6 ,          & !SSA,SSC
    7 ,7 ,7, 8, 8,  & !OCC,BC1
    8, 9, 9, 9,10,  & !BC1,BC2,BC3
    10,10,11,11,11, & !BC3,DBC
    11,12,12,12,12, & !DBC,BOC
    13,13,13,14,14, & !BCS,MXX
    14,14,14,14/)
  integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,0 ,1,  & !AKK
    0 ,2 ,0 ,0 ,3,  & !ACC,DD1
    0 ,0 ,4 ,       & !DS1
    0 ,             & !SSA
    0 ,0 ,          & !SSA,SSC
    0 ,0 ,7, 0, 0,  & !OCC,BC1
    8, 0, 0, 9, 0,  & !BC1,BC2,BC3
    0, 10,0 ,0 ,0, & !BC3,DBC
    11,0 ,0 ,0 ,12, & !DBC,BOC
    0, 0, 13,0 ,0, & !BCS,MXX
    0,0,0,14/)
  integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/ &
    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10, &
    11,12,13,14,15,      18,   20, &
    21,22,23,24,25,26,27,28,29,30, &
    31,32,33,34,35,36,37,38,39,40, &
    41,42,43,44,45,46,47,48   /)
#endif
#ifdef TRACERS_AMP_M6
  integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,1 ,1,  & !AKK
    2 ,2 ,3 ,3 ,3,  & !ACC,DD1
    4 ,4 ,4 ,       & !DS1
    5 ,             & !SSA
    5 ,6 ,          & !SSA,SSC
    7 ,7 ,7, 8, 8,  & !OCC,BC1
    8, 9, 9, 9,10,  & !BC1,BC2,OCS
    10,10,11,11,11, & !OCS,DBC
    11,12,12,12,12, & !DBC,BOC
    13,13,13,14,14, & !BCS,MXX
    14,14,14,14/)
  integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,0 ,1,  & !AKK
    0 ,2 ,0 ,0 ,3,  & !ACC,DD1
    0 ,0 ,4 ,       & !DS1
    0 ,             & !SSA
    0 ,0 ,          & !SSA,SSC
    0 ,0 ,7, 0, 0,  & !OCC,BC1
    8, 0, 0, 9, 0,  & !BC1,BC2,OCS
    0 ,10,0 ,0 ,0 , & !OCS,DBC
    11,0 ,0 ,0 ,12, & !DBC,BOC
    0 ,0 ,13,0 ,0 , & !BCS,MXX
    0 ,0 ,0 ,14/)
  integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/ &
    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10, &
    11,12,13,14,15,      18,   20, &
    21,22,23,24,25,26,27,28,29,30, &
    31,32,33,34,35,36,37,38,39,40, &
    41,42,43,44,45,46,47,48   /)
#endif
#ifdef TRACERS_AMP_M7
  integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,1 ,1,  & !AKK
    2 ,2 ,3 ,3 ,3,  & !ACC,DD1
    4 ,4 ,4 ,       & !DS1
    5 ,             & !SSA
    5 ,6 ,          & !SSA,SSC
    7 ,7 ,7, 8, 8,  & !OCC,BC1
    8, 9, 9, 9,     & !BC1,BC2
    10,10,10,10,    & !OCS,DBC
    11,11,          & !MXX
    11,11,11,11/)   !MXX
  integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,0 ,1,  & !AKK
    0 ,2 ,0 ,0 ,3,  & !ACC,DD1
    0 ,0 ,4 ,       & !DS1
    0 ,             & !SSA
    0 ,0 ,          & !SSA,SSC
    0 ,0 ,7, 0, 0,  & !OCC,BC1
    8, 0, 0, 9,     & !BC1,BC2
    0, 0, 0,10,    & !OCS,DBC
    0, 0,          & !MXX
    0, 0, 0,11/)   !MXX
  integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/ &
    1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10, &
    11,12,13,14,15,      18,   20, &
    21,22,23,24,25,26,27,28,29,30, &
    31,32,33,34,35,36,37,38/)
#endif
#ifdef TRACERS_AMP_M8
  integer, parameter :: AMP_MODES_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,       & !
    1 ,1 ,2 ,2 ,2,  & !ACC,DD1
    3 ,3 ,3 ,       & !DS1
    4 ,             & !SSS
    4 ,5 ,          & !SSS
    5 ,5 ,5, 6, 6,  & !OCC,BC1
    6, 7, 7, 7, 8,  & !BC1,BC2
    8,8,8,8,8/)    !MXX
  integer, parameter :: AMP_NUMB_MAP(ntmAMP)=(/ &
    0 ,0 ,0 ,       & !
    0 ,1 ,0 ,0 ,2,  & !ACC,DD1
    0 ,0 ,3 ,       & !DS1
    0 ,             & !SSS
    0 ,0 ,          & !SSS
    0 ,0 ,5, 0, 0,  & !OCC,BC1
    6, 0, 0, 7, 0,  & !BC1,BC2
    0,0,0,0,8/)    !MXX
  integer, parameter :: AMP_AERO_MAP(ntmAMP)=(/ &
    1 ,2 ,3 ,4 ,5 , &
    6 ,7 ,8 ,9 ,10, &
    11,12,13,15,16, &
    17,18,19,20,21, &
    22,23,24,25,26, &
    27,28,29   /)
#endif

  public AMP_MODES_MAP
  public AMP_NUMB_MAP
  public AMP_AERO_MAP
  public AMP_TRM_NM1
  public AMP_TRM_NM2

  real(8), parameter :: microns2meters = 1.0d-6
  REAL(8), PARAMETER :: RG_AKK = microns2meters*DG_AKK/2.0d0 
  REAL(8), PARAMETER :: RG_ACC = microns2meters*DG_ACC/2.0d0
  REAL(8), PARAMETER :: RG_DD1 = microns2meters*DG_DD1/2.0d0
  REAL(8), PARAMETER :: RG_DD2 = microns2meters*DG_DD2/2.0d0
  REAL(8), PARAMETER :: RG_DS1 = microns2meters*DG_DS1/2.0d0
  REAL(8), PARAMETER :: RG_DS2 = microns2meters*DG_DS2/2.0d0
  REAL(8), PARAMETER :: RG_SSA = microns2meters*DG_SSA/2.0d0
  REAL(8), PARAMETER :: RG_SSC = microns2meters*DG_SSC/2.0d0
  REAL(8), PARAMETER :: RG_SSS = microns2meters*DG_SSS/2.0d0
  REAL(8), PARAMETER :: RG_OCC = microns2meters*DG_OCC/2.0d0
  REAL(8), PARAMETER :: RG_BC1 = microns2meters*DG_BC1/2.0d0
  REAL(8), PARAMETER :: RG_BC2 = microns2meters*DG_BC2/2.0d0
  REAL(8), PARAMETER :: RG_BC3 = microns2meters*DG_BC3/2.0d0
  REAL(8), PARAMETER :: RG_DBC = microns2meters*DG_DBC/2.0d0
  REAL(8), PARAMETER :: RG_BOC = microns2meters*DG_BOC/2.0d0
  REAL(8), PARAMETER :: RG_BCS = microns2meters*DG_BCS/2.0d0
  REAL(8), PARAMETER :: RG_OCS = microns2meters*DG_OCS/2.0d0
  REAL(8), PARAMETER :: RG_MXX = microns2meters*DG_MXX/2.0d0

  real(8), parameter :: SULF_MolecMass = 96.0d0
  real(8), parameter :: DUST_MolecMass = 1.0d0
  real(8), parameter :: SEAS_MolecMass = 75.0d0
  real(8), parameter :: BCAR_MolecMass = 12.0d0
  real(8), parameter :: OCAR_MolecMass = 12.0d0

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine AMP_initMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    !**** Tracers for Scheme AMP: Aerosol Microphysics (Mechanism M1 - M8)
    call  M_NO3_setSpec('M_NO3')
    call  M_NH4_setSpec('M_NH4')
    call  M_H2O_setSpec('M_H2O')

    if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3 &
      .or. tracers_amp_m5 .or. tracers_amp_m6 .or.  &
      tracers_amp_m7) then
      n_M_AKK_SU = AMP_SetSpec('AKK', 'SU')
      n_N_AKK_1  = AMP_SetSpec('AKK', '1' )
    end if

    n_M_ACC_SU = AMP_SetSpec('ACC', 'SU')
    n_N_ACC_1  = AMP_SetSpec('ACC', '1' )
    n_M_DD1_SU = AMP_SetSpec('DD1', 'SU')
    n_M_DD1_DU = AMP_SetSpec('DD1', 'DU')
    n_N_DD1_1  = AMP_SetSpec('DD1', '1' )
    n_M_DS1_SU = AMP_SetSpec('DS1', 'SU')
    n_M_DS1_DU = AMP_SetSpec('DS1', 'DU')
    n_N_DS1_1  = AMP_SetSpec('DS1', '1' )

    if (tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3 &
      .or. tracers_amp_m4) then
      n_M_DD2_SU = AMP_SetSpec('DD2','SU')
      n_M_DD2_DU = AMP_SetSpec('DD2','DU')
      n_N_DD2_1  = AMP_SetSpec('DD2','1' )
      n_M_DS2_SU = AMP_SetSpec('DS2','SU')
      n_M_DS2_DU = AMP_SetSpec('DS2','DU')
      n_N_DS2_1  = AMP_SetSpec('DS2','1' )
    end if

    if (tracers_amp_m1 .or.  tracers_amp_m2 .or. tracers_amp_m3 &
      .or. tracers_amp_m5 .or.  tracers_amp_m6 .or.  &
      tracers_amp_m7) then
      n_M_SSA_SU = AMP_SetSpec('SSA','SU')
      n_M_SSA_SS = AMP_SetSpec('SSA','SS')
      n_N_SSA_1  = AMP_SetSpec('SSA','1' )
      n_M_SSC_SS = AMP_SetSpec('SSC','SS')
      n_N_SSC_1  = AMP_SetSpec('SSC','1' )
    end if

    if ( tracers_amp_m4 .or. tracers_amp_m8) then ! cases not tested
      n_M_SSS_SU = AMP_SetSpec('SSS','SU') ! need rundeck with these
      n_M_SSS_SS = AMP_SetSpec('SSS','SS') ! settings
    end if

    n_M_OCC_SU = AMP_SetSpec('OCC','SU')
    n_M_OCC_OC = AMP_SetSpec('OCC','OC')
    n_N_OCC_1  = AMP_SetSpec('OCC','1' )
    n_M_BC1_SU = AMP_SetSpec('BC1','SU')
    n_M_BC1_BC = AMP_SetSpec('BC1','BC')
    n_N_BC1_1  = AMP_SetSpec('BC1','1' )
    n_M_BC2_SU = AMP_SetSpec('BC2','SU')
    n_M_BC2_BC = AMP_SetSpec('BC2','BC')
    n_N_BC2_1  = AMP_SetSpec('BC2','1' )

    if ( tracers_amp_m1 .or. tracers_amp_m5) then
      n_M_BC3_SU = AMP_SetSpec('BC3','SU')
      n_M_BC3_BC = AMP_SetSpec('BC3','BC')
      n_N_BC3_1  = AMP_SetSpec('BC3','1' )
    end if

    if ( tracers_amp_m2 .or. tracers_amp_m6) then
      n_M_OCS_SU = AMP_SetSpec('OCS','SU')
      n_M_OCS_OC = AMP_SetSpec('OCS','OC')
      n_N_OCS_1  = AMP_SetSpec('OCS','1' )
    end if

    if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m6) then
      n_M_DBC_SU = AMP_SetSpec('DBC','SU') 
      n_M_DBC_BC = AMP_SetSpec('DBC','BC')
      n_M_DBC_DU = AMP_SetSpec('DBC','DU') 
      n_N_DBC_1  = AMP_SetSpec('DBC','1' ) 
    end if

    if (tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3 &
      .or. tracers_amp_m6 .or. tracers_amp_m7) then
      n_M_BOC_SU = AMP_SetSpec('BOC','SU')
      n_M_BOC_BC = AMP_SetSpec('BOC','BC')
      n_M_BOC_OC = AMP_SetSpec('BOC','OC')
      n_N_BOC_1  = AMP_SetSpec('BOC','1' )
    end if

    if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m5 &
      .or. TRACERS_AMP_M6) then
      n_M_BCS_SU = AMP_SetSpec('BCS','SU')
      n_M_BCS_BC = AMP_SetSpec('BCS','BC')
      n_N_BCS_1  = AMP_SetSpec('BCS','1' )
    end if

    n_M_MXX_SU = AMP_SetSpec('MXX','SU')
    n_M_MXX_BC = AMP_SetSpec('MXX','BC')
    n_M_MXX_OC = AMP_SetSpec('MXX','OC')
    n_M_MXX_DU = AMP_SetSpec('MXX','DU')
    n_M_MXX_SS = AMP_SetSpec('MXX','SS')
    n_N_MXX_1  = AMP_SetSpec('MXX','1' )

    call  H2SO4_setSpec('H2SO4')
    call  DMS_setSpec('DMS')  ! duplicate with Koch
    call  SO2_setSpec('SO2')  ! duplicate with Koch
    if (.not. tracers_special_shindell .or. coupled_chem.eq.0) then
      call  H2O2_s_setSpec('H2O2_s') ! duplicate with Koch
    endif
    call  NH3_setSpec('NH3')  ! duplicate with nitrate
    if (tracers_aerosols_koch.or.tracers_aerosols_seasalt) then
      call stop_model('contradictory tracer specs', 255)
    end if
    if (tracers_nitrate) then
      call stop_model('contradictory tracer specs', 255)
    end if

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine H2SO4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_H2SO4 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 98.d0)
      call set_trpdens(n, DENS_SULF)
      call set_trradius(n, DG_ACC * .5d-6)
      call set_fq_aer(n, SOLU_ACC)
      call set_tr_wd_type(n, npart)
      call set_has_chemistry(n, .true.)
    end subroutine H2SO4_setSpec

    subroutine M_NO3_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_NO3 = n
      ntmAMPi=n                 ! always the first tracer in AMP
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 62.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0)  !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_has_chemistry(n, .true.)
    end subroutine M_NO3_setSpec

    subroutine M_NH4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_NH4 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 18.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.d+0)
      call set_tr_wd_type(n, npart)
      call set_has_chemistry(n, .true.)
    end subroutine M_NH4_setSpec

    subroutine M_H2O_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_M_H2O = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, mwat)
      call set_trpdens(n, 1.d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.d+0)
      call set_tr_wd_type(n, npart) !nWater
    end subroutine M_H2O_setSpec

!------------------------------------------------------------------------------
    function AMP_setSpec(mode, component) result (tracerIndex)
!------------------------------------------------------------------------------
      use OldTracer_mod, only: om2oc, set_om2oc
      use Dictionary_mod, only: sync_param
      use RunTimeControls_mod, only: dynamic_biomass_burning
      implicit none
      character(len=*), intent(in) :: mode
      character(len=*), intent(in) :: component
      ! local variables
      character(len=1) :: prefix
      character(len=64) :: tracerName
      integer :: tracerIndex
      real*8 :: tmp

      prefix = getTracerPrefix(component)
      tracerName = prefix // "_" // mode // "_" // component
      tracerIndex = oldAddTracer(trim(tracerName))

      if (trim(tracerName) == 'N_MXX_1') then
        ntmAMPe = tracerIndex     ! always the last tracer in AMP
        if (ntmAMPi + ntmAMP-1 /= ntmAMPe) &
          call stop_model( 'ntmAMPi+ntmAMP-1 /= ntmAMPe', 255 )
      end if


      if (trim(component) == 'OC') then
        call set_om2oc(tracerIndex, 1.4d0)
        tmp = om2oc(tracerIndex)
        call sync_param(trim(tracerName)//"_om2oc",tmp)
        call set_om2oc(tracerIndex, tmp)
      endif
      call set_ntm_power(tracerIndex, -11)
      if (trim(component) == 'OC') then
        tmp = getMolecularMass(component) * om2oc(tracerIndex)
      else
        tmp = getMolecularMass(component)
      endif
      call set_tr_mm(tracerIndex, tmp)
      call set_trpdens(tracerIndex, getDensity(component))
      call set_trradius(tracerIndex, getRadius(mode))
      call set_fq_aer(tracerIndex, getSolubility(mode))
      call set_tr_wd_type(tracerIndex, nPART)

    end function AMP_setSpec

!------------------------------------------------------------------------------
    function getTracerPrefix(component) result (unitPrefix)
!------------------------------------------------------------------------------
      implicit none
      character(len=*), intent(in) :: component
      character(len=1) :: unitPrefix

      select case (trim(component))
      case ('1')
        unitPrefix = 'N'        ! number density
      case default
        unitPrefix = 'M'        ! mixing ratio
      end select

    end function getTracerPrefix

!------------------------------------------------------------------------------
    function getMolecularMass(component) result (molecularMass)
!------------------------------------------------------------------------------
      implicit none
      character(len=*), intent(in) :: component
      ! local variables
      real(8) :: molecularMass

      select case (trim(component))
      case ('SU')
        molecularMass = SULF_MolecMass
      case ('DU')
        molecularMass = DUST_MolecMass
      case ('SS')
        molecularMass = SEAS_MolecMass
      case ('BC')
        molecularMass = BCAR_MolecMass
      case ('OC')
        molecularMass = OCAR_MolecMass
      case ('1')
        molecularMass = 1.0d0
      case default
        call stop_model('Incorrect molecular mass choice', 255)
      end select

    end function getMolecularMass

!------------------------------------------------------------------------------
    function getDensity(component)  result (density)
!------------------------------------------------------------------------------
      implicit none
      character(len=*), intent(in) :: component
      real(8) :: density

      select case (trim(component))
      case ('SU')
        density = DENS_SULF
      case ('DU')
        density = DENS_DUST
      case ('SS')
        density = DENS_SEAS
      case ('BC')
        density = DENS_BCAR
      case ('OC')
        density = DENS_OCAR
      case ('1')
        density = 1.0d0
      case default
        call stop_model('Incorrect density choice', 255)
      end select

    end function getDensity

!------------------------------------------------------------------------------
    function getRadius(mode)  result (radius)
!------------------------------------------------------------------------------
      implicit none
      character(len=*), intent(in) :: mode
      real(8) :: radius

      select case (trim(mode))
      case ('AKK')
        radius = RG_AKK
      case ('ACC')
        radius = RG_ACC
      case ('DD1')
        radius = RG_DD1
      case ('DS1')
        radius = RG_DS1
      case ('DD2')
        radius = RG_DD2
      case ('DS2')
        radius = RG_DS2
      case ('SSA')
        radius = RG_SSA
      case ('SSC')
        radius = RG_SSC
      case ('SSS')
        radius = RG_SSS
      case ('OCC')
        radius = RG_OCC
      case ('BC1')
        radius = RG_BC1
      case ('BC2')
        radius = RG_BC2
      case ('BC3')
        radius = RG_BC3
      case ('DBC')
        radius = RG_DBC
      case ('BOC')
        radius = RG_BOC
      case ('BCS')
        radius = RG_BCS
      case ('MXX')
        radius = RG_MXX
      case ('OCS')
        radius = RG_OCS
      case default
        call stop_model('Incorrect radius choice', 255)
      end select

    end function getRadius

!------------------------------------------------------------------------------
    function getSolubility(mode)  result (solubility)
!------------------------------------------------------------------------------
      implicit none
      character(len=*), intent(in) :: mode
      real(8) :: solubility

      select case (trim(mode))
      case ('AKK')
        solubility = SOLU_AKK
      case ('ACC')
        solubility = SOLU_ACC
      case ('DD1')
        solubility = SOLU_DD1
      case ('DS1')
        solubility = SOLU_DS1
      case ('DD2')
        solubility = SOLU_DD2
      case ('DS2')
        solubility = SOLU_DS2
      case ('SSA')
        solubility = SOLU_SSA
      case ('SSC')
        solubility = SOLU_SSC
      case ('SSS')
        solubility = SOLU_SSS
      case ('OCC')
        solubility = SOLU_OCC
      case ('BC1')
        solubility = SOLU_BC1
      case ('BC2')
        solubility = SOLU_BC2
      case ('BC3')
        solubility = SOLU_BC3
      case ('DBC')
        solubility = SOLU_DBC
      case ('BOC')
        solubility = SOLU_BOC
      case ('BCS')
        solubility = SOLU_BCS
      case ('MXX')
        solubility = SOLU_MXX
      case ('OCS')
        solubility = SOLU_OCS
      case default
        call stop_model('Incorrect solubility choice', 255)
      end select

    end function getSolubility

  end subroutine AMP_initMetadata

end module AmpTracersMetadata_mod
