module @ID@_top_real_decay
  use kinds
  use diagnostics
  use omega95
  use parameters_SM_tt_threshold
  implicit none
  private
  public :: calculate_amplitude

  integer, parameter :: n_prt = 4

  type(momentum) :: p1, p2, p3, p4
  type(momentum) :: p12, p14
  type(spinor) :: owf_u3_1__1_0
  type(conjspinor) :: owf_d3b__1_3_0
  type(vector) :: owf_gl___4_0, owf_wm_2_0
  type(spinor) :: owf_d3_1__12_0, owf_u3_1__14_0_X1
  complex(default) :: oks_u3_1_wpd3_1_gl__

contains

  function calculate_amplitude (k_decay, s, top_width) result (amp)
    complex(default) :: amp
    type(momentum), dimension(:), intent(in) :: k_decay
    integer, dimension(n_prt), intent(in) :: s
    real(default), intent(in) :: top_width
    real(default) :: dynamic_top_mass
    p1 = - k_decay(1) ! incoming
    p2 =   k_decay(2) ! outgoing
    p3 =   k_decay(3) ! outgoing
    p4 =   k_decay(4) ! outgoing
    p12 = p1 + p2
    p14 = p1 + p4
    dynamic_top_mass = sqrt (p1 * p1)
    owf_u3_1__1_0 = u (dynamic_top_mass, - p1, s(1))
    owf_wm_2_0 = conjg (eps (mass(24), p2, s(2)))
    owf_d3b__1_3_0 = ubar (mass(5), p3, s(3))
    owf_gl___4_0 = conjg (eps (mass(21), p4, s(4)))
    owf_d3_1__12_0 = pr_psi(p12,mass(5),wd_tl(p12,width(5)), .false., &
       + f_vlf(gcc,owf_wm_2_0,owf_u3_1__1_0))
    owf_u3_1__14_0_X1 = pr_psi(p14,dynamic_top_mass,wd_tl(p14,top_width), &
       .false., + f_vf((-gs),owf_gl___4_0,owf_u3_1__1_0))
    oks_u3_1_wpd3_1_gl__ = ( &
       + f_fv((-gs),owf_d3b__1_3_0,owf_gl___4_0))*owf_d3_1__12_0
    oks_u3_1_wpd3_1_gl__ = oks_u3_1_wpd3_1_gl__ + ( &
       + f_fvl(gcc,owf_d3b__1_3_0,owf_wm_2_0))*owf_u3_1__14_0_X1
    amp = - oks_u3_1_wpd3_1_gl__ ! 2 vertices, 1 propagators
  end function calculate_amplitude

end module @ID@_top_real_decay

module @ID@_anti_top_real_decay
  use kinds
  use diagnostics
  use omega95
  use parameters_SM_tt_threshold
  implicit none
  private
  public :: calculate_amplitude

  integer, parameter :: n_prt = 4

  type(momentum) :: p1, p2, p3, p4
  type(momentum) :: p12, p14
  type(spinor) :: owf_d3_2__3_0
  type(conjspinor) :: owf_u3b__1_1_0
  type(vector) :: owf_gl_1_2_4_0, owf_wp_2_0
  type(conjspinor) :: owf_d3b__1_12_0, owf_u3b__2_14_0
  complex(default) :: oks_u3b__1wmd3b__2gl_2_1

contains

  function calculate_amplitude (k_decay, s, top_width) result (amp)
    complex(default) :: amp
    type(momentum), dimension(:), intent(in) :: k_decay
    integer, dimension(n_prt), intent(in) :: s
    real(default), intent(in) :: top_width
    real(default) :: dynamic_top_mass
    p1 = - k_decay(1) ! incoming
    p2 =   k_decay(2) ! outgoing
    p3 =   k_decay(3) ! outgoing
    p4 =   k_decay(4) ! outgoing
    p12 = p1 + p2
    p14 = p1 + p4
    dynamic_top_mass = sqrt (p1 * p1)
    owf_u3b__1_1_0 = vbar (dynamic_top_mass, - p1, s(1))
    owf_wp_2_0 = conjg (eps (mass(24), p2, s(2)))
    owf_d3_2__3_0 = v (mass(5), p3, s(3))
    owf_gl_1_2_4_0 = conjg (eps (mass(21), p4, s(4)))
    owf_d3b__1_12_0 = pr_psibar(p12,mass(5),wd_tl(p12,width(5)), .false., &
       + f_fvl(gcc,owf_u3b__1_1_0,owf_wp_2_0))
    owf_u3b__2_14_0 = pr_psibar(p14,dynamic_top_mass,wd_tl(p14,top_width), &
       .false., + f_fv((-gs),owf_u3b__1_1_0,owf_gl_1_2_4_0))
    oks_u3b__1wmd3b__2gl_2_1 = owf_d3b__1_12_0 * &
       f_vf((-gs),owf_gl_1_2_4_0,owf_d3_2__3_0)
    oks_u3b__1wmd3b__2gl_2_1 = oks_u3b__1wmd3b__2gl_2_1 + owf_u3b__2_14_0*( &
       + f_vlf(gcc,owf_wp_2_0,owf_d3_2__3_0))
    amp = - oks_u3b__1wmd3b__2gl_2_1 ! 2 vertices, 1 propagators
  end function calculate_amplitude

end module @ID@_anti_top_real_decay

module @ID@_threshold
  use kinds
  use diagnostics
  use numeric_utils
  use physics_defs, only: THR_POS_WP, THR_POS_WM, THR_POS_B, THR_POS_BBAR, THR_POS_GLUON
  use physics_defs, only: ass_boson, ass_quark
  use physics_defs, only: PROC_MODE_UNDEFINED, PROC_MODE_TT, PROC_MODE_WBWB
  use constants
  use lorentz
  use omega95
  use parameters_SM_tt_threshold
  use ttv_formfactors
  use @ID@_top_real_decay, top_real_decay_calculate_amplitude => calculate_amplitude
  use @ID@_anti_top_real_decay, anti_top_real_decay_calculate_amplitude => calculate_amplitude
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none
  private
  public :: init, calculate_blob, compute_born, &
       compute_production_owfs, &
       compute_decay_owfs, table_spin_states, compute_production_me, &
       top_decay_born, anti_top_decay_born, top_propagators, compute_real, abs2, &
       apply_boost, compute_projected_momenta, convert_to_mom_and_invert_sign, &
       handle_onshell_test_point, set_production_factors, get_selected_beam_helicities

  logical, parameter, public :: test_ward = .false.
  logical, parameter, public :: test_onshell = .false.

  integer, parameter, public :: n_prt = 6
  integer, parameter :: n_prt_OS = 4
  integer, parameter, public :: n_in = 2
  integer, parameter, public :: n_out = 4
  integer, parameter, public :: n_cindex = 2
  integer, parameter, public :: n_flv = 1
  integer, parameter, public :: n_hel = 144
  integer, parameter, public :: n_hel_OS = 16

  integer, public :: process_mode = PROC_MODE_UNDEFINED

  integer, parameter :: OFS = 1
  integer, parameter :: ONS = 2
  integer, parameter :: ONS_BOOST = 3


  ! NB: you MUST NOT change the value of N_ here!!!
  !     It is defined here for convenience only and must be
  !     compatible with hardcoded values in the amplitude!
  real(default), parameter, public :: N_ = 3
  real(default), public :: production_factors

  integer, dimension(n_prt_OS,n_hel_OS), save, protected :: table_spin_states_OS
  data table_spin_states_OS(:,   1) / -1, -1, -1, -1 /
  data table_spin_states_OS(:,   2) / -1, -1, -1,  1 /
  data table_spin_states_OS(:,   3) / -1, -1,  1, -1 /
  data table_spin_states_OS(:,   4) / -1, -1,  1,  1 /
  data table_spin_states_OS(:,   5) / -1,  1, -1, -1 /
  data table_spin_states_OS(:,   6) / -1,  1, -1,  1 /
  data table_spin_states_OS(:,   7) / -1,  1,  1, -1 /
  data table_spin_states_OS(:,   8) / -1,  1,  1,  1 /
  data table_spin_states_OS(:,   9) /  1, -1, -1, -1 /
  data table_spin_states_OS(:,  10) /  1, -1, -1,  1 /
  data table_spin_states_OS(:,  11) /  1, -1,  1, -1 /
  data table_spin_states_OS(:,  12) /  1, -1,  1,  1 /
  data table_spin_states_OS(:,  13) /  1,  1, -1, -1 /
  data table_spin_states_OS(:,  14) /  1,  1, -1,  1 /
  data table_spin_states_OS(:,  15) /  1,  1,  1, -1 /
  data table_spin_states_OS(:,  16) /  1,  1,  1,  1 /

  integer, dimension(n_prt,n_hel), save, protected :: table_spin_states
  data table_spin_states(:,   1) / -1, -1, -1, -1, -1, -1 /
  data table_spin_states(:,   2) / -1, -1, -1, -1, -1,  1 /
  data table_spin_states(:,   3) / -1, -1, -1, -1,  1, -1 /
  data table_spin_states(:,   4) / -1, -1, -1, -1,  1,  1 /
  data table_spin_states(:,   5) / -1, -1, -1,  0, -1, -1 /
  data table_spin_states(:,   6) / -1, -1, -1,  0, -1,  1 /
  data table_spin_states(:,   7) / -1, -1, -1,  0,  1, -1 /
  data table_spin_states(:,   8) / -1, -1, -1,  0,  1,  1 /
  data table_spin_states(:,   9) / -1, -1, -1,  1, -1, -1 /
  data table_spin_states(:,  10) / -1, -1, -1,  1, -1,  1 /
  data table_spin_states(:,  11) / -1, -1, -1,  1,  1, -1 /
  data table_spin_states(:,  12) / -1, -1, -1,  1,  1,  1 /
  data table_spin_states(:,  13) / -1, -1,  0, -1, -1, -1 /
  data table_spin_states(:,  14) / -1, -1,  0, -1, -1,  1 /
  data table_spin_states(:,  15) / -1, -1,  0, -1,  1, -1 /
  data table_spin_states(:,  16) / -1, -1,  0, -1,  1,  1 /
  data table_spin_states(:,  17) / -1, -1,  0,  0, -1, -1 /
  data table_spin_states(:,  18) / -1, -1,  0,  0, -1,  1 /
  data table_spin_states(:,  19) / -1, -1,  0,  0,  1, -1 /
  data table_spin_states(:,  20) / -1, -1,  0,  0,  1,  1 /
  data table_spin_states(:,  21) / -1, -1,  0,  1, -1, -1 /
  data table_spin_states(:,  22) / -1, -1,  0,  1, -1,  1 /
  data table_spin_states(:,  23) / -1, -1,  0,  1,  1, -1 /
  data table_spin_states(:,  24) / -1, -1,  0,  1,  1,  1 /
  data table_spin_states(:,  25) / -1, -1,  1, -1, -1, -1 /
  data table_spin_states(:,  26) / -1, -1,  1, -1, -1,  1 /
  data table_spin_states(:,  27) / -1, -1,  1, -1,  1, -1 /
  data table_spin_states(:,  28) / -1, -1,  1, -1,  1,  1 /
  data table_spin_states(:,  29) / -1, -1,  1,  0, -1, -1 /
  data table_spin_states(:,  30) / -1, -1,  1,  0, -1,  1 /
  data table_spin_states(:,  31) / -1, -1,  1,  0,  1, -1 /
  data table_spin_states(:,  32) / -1, -1,  1,  0,  1,  1 /
  data table_spin_states(:,  33) / -1, -1,  1,  1, -1, -1 /
  data table_spin_states(:,  34) / -1, -1,  1,  1, -1,  1 /
  data table_spin_states(:,  35) / -1, -1,  1,  1,  1, -1 /
  data table_spin_states(:,  36) / -1, -1,  1,  1,  1,  1 /
  data table_spin_states(:,  37) / -1,  1, -1, -1, -1, -1 /
  data table_spin_states(:,  38) / -1,  1, -1, -1, -1,  1 /
  data table_spin_states(:,  39) / -1,  1, -1, -1,  1, -1 /
  data table_spin_states(:,  40) / -1,  1, -1, -1,  1,  1 /
  data table_spin_states(:,  41) / -1,  1, -1,  0, -1, -1 /
  data table_spin_states(:,  42) / -1,  1, -1,  0, -1,  1 /
  data table_spin_states(:,  43) / -1,  1, -1,  0,  1, -1 /
  data table_spin_states(:,  44) / -1,  1, -1,  0,  1,  1 /
  data table_spin_states(:,  45) / -1,  1, -1,  1, -1, -1 /
  data table_spin_states(:,  46) / -1,  1, -1,  1, -1,  1 /
  data table_spin_states(:,  47) / -1,  1, -1,  1,  1, -1 /
  data table_spin_states(:,  48) / -1,  1, -1,  1,  1,  1 /
  data table_spin_states(:,  49) / -1,  1,  0, -1, -1, -1 /
  data table_spin_states(:,  50) / -1,  1,  0, -1, -1,  1 /
  data table_spin_states(:,  51) / -1,  1,  0, -1,  1, -1 /
  data table_spin_states(:,  52) / -1,  1,  0, -1,  1,  1 /
  data table_spin_states(:,  53) / -1,  1,  0,  0, -1, -1 /
  data table_spin_states(:,  54) / -1,  1,  0,  0, -1,  1 /
  data table_spin_states(:,  55) / -1,  1,  0,  0,  1, -1 /
  data table_spin_states(:,  56) / -1,  1,  0,  0,  1,  1 /
  data table_spin_states(:,  57) / -1,  1,  0,  1, -1, -1 /
  data table_spin_states(:,  58) / -1,  1,  0,  1, -1,  1 /
  data table_spin_states(:,  59) / -1,  1,  0,  1,  1, -1 /
  data table_spin_states(:,  60) / -1,  1,  0,  1,  1,  1 /
  data table_spin_states(:,  61) / -1,  1,  1, -1, -1, -1 /
  data table_spin_states(:,  62) / -1,  1,  1, -1, -1,  1 /
  data table_spin_states(:,  63) / -1,  1,  1, -1,  1, -1 /
  data table_spin_states(:,  64) / -1,  1,  1, -1,  1,  1 /
  data table_spin_states(:,  65) / -1,  1,  1,  0, -1, -1 /
  data table_spin_states(:,  66) / -1,  1,  1,  0, -1,  1 /
  data table_spin_states(:,  67) / -1,  1,  1,  0,  1, -1 /
  data table_spin_states(:,  68) / -1,  1,  1,  0,  1,  1 /
  data table_spin_states(:,  69) / -1,  1,  1,  1, -1, -1 /
  data table_spin_states(:,  70) / -1,  1,  1,  1, -1,  1 /
  data table_spin_states(:,  71) / -1,  1,  1,  1,  1, -1 /
  data table_spin_states(:,  72) / -1,  1,  1,  1,  1,  1 /
  data table_spin_states(:,  73) /  1, -1, -1, -1, -1, -1 /
  data table_spin_states(:,  74) /  1, -1, -1, -1, -1,  1 /
  data table_spin_states(:,  75) /  1, -1, -1, -1,  1, -1 /
  data table_spin_states(:,  76) /  1, -1, -1, -1,  1,  1 /
  data table_spin_states(:,  77) /  1, -1, -1,  0, -1, -1 /
  data table_spin_states(:,  78) /  1, -1, -1,  0, -1,  1 /
  data table_spin_states(:,  79) /  1, -1, -1,  0,  1, -1 /
  data table_spin_states(:,  80) /  1, -1, -1,  0,  1,  1 /
  data table_spin_states(:,  81) /  1, -1, -1,  1, -1, -1 /
  data table_spin_states(:,  82) /  1, -1, -1,  1, -1,  1 /
  data table_spin_states(:,  83) /  1, -1, -1,  1,  1, -1 /
  data table_spin_states(:,  84) /  1, -1, -1,  1,  1,  1 /
  data table_spin_states(:,  85) /  1, -1,  0, -1, -1, -1 /
  data table_spin_states(:,  86) /  1, -1,  0, -1, -1,  1 /
  data table_spin_states(:,  87) /  1, -1,  0, -1,  1, -1 /
  data table_spin_states(:,  88) /  1, -1,  0, -1,  1,  1 /
  data table_spin_states(:,  89) /  1, -1,  0,  0, -1, -1 /
  data table_spin_states(:,  90) /  1, -1,  0,  0, -1,  1 /
  data table_spin_states(:,  91) /  1, -1,  0,  0,  1, -1 /
  data table_spin_states(:,  92) /  1, -1,  0,  0,  1,  1 /
  data table_spin_states(:,  93) /  1, -1,  0,  1, -1, -1 /
  data table_spin_states(:,  94) /  1, -1,  0,  1, -1,  1 /
  data table_spin_states(:,  95) /  1, -1,  0,  1,  1, -1 /
  data table_spin_states(:,  96) /  1, -1,  0,  1,  1,  1 /
  data table_spin_states(:,  97) /  1, -1,  1, -1, -1, -1 /
  data table_spin_states(:,  98) /  1, -1,  1, -1, -1,  1 /
  data table_spin_states(:,  99) /  1, -1,  1, -1,  1, -1 /
  data table_spin_states(:, 100) /  1, -1,  1, -1,  1,  1 /
  data table_spin_states(:, 101) /  1, -1,  1,  0, -1, -1 /
  data table_spin_states(:, 102) /  1, -1,  1,  0, -1,  1 /
  data table_spin_states(:, 103) /  1, -1,  1,  0,  1, -1 /
  data table_spin_states(:, 104) /  1, -1,  1,  0,  1,  1 /
  data table_spin_states(:, 105) /  1, -1,  1,  1, -1, -1 /
  data table_spin_states(:, 106) /  1, -1,  1,  1, -1,  1 /
  data table_spin_states(:, 107) /  1, -1,  1,  1,  1, -1 /
  data table_spin_states(:, 108) /  1, -1,  1,  1,  1,  1 /
  data table_spin_states(:, 109) /  1,  1, -1, -1, -1, -1 /
  data table_spin_states(:, 110) /  1,  1, -1, -1, -1,  1 /
  data table_spin_states(:, 111) /  1,  1, -1, -1,  1, -1 /
  data table_spin_states(:, 112) /  1,  1, -1, -1,  1,  1 /
  data table_spin_states(:, 113) /  1,  1, -1,  0, -1, -1 /
  data table_spin_states(:, 114) /  1,  1, -1,  0, -1,  1 /
  data table_spin_states(:, 115) /  1,  1, -1,  0,  1, -1 /
  data table_spin_states(:, 116) /  1,  1, -1,  0,  1,  1 /
  data table_spin_states(:, 117) /  1,  1, -1,  1, -1, -1 /
  data table_spin_states(:, 118) /  1,  1, -1,  1, -1,  1 /
  data table_spin_states(:, 119) /  1,  1, -1,  1,  1, -1 /
  data table_spin_states(:, 120) /  1,  1, -1,  1,  1,  1 /
  data table_spin_states(:, 121) /  1,  1,  0, -1, -1, -1 /
  data table_spin_states(:, 122) /  1,  1,  0, -1, -1,  1 /
  data table_spin_states(:, 123) /  1,  1,  0, -1,  1, -1 /
  data table_spin_states(:, 124) /  1,  1,  0, -1,  1,  1 /
  data table_spin_states(:, 125) /  1,  1,  0,  0, -1, -1 /
  data table_spin_states(:, 126) /  1,  1,  0,  0, -1,  1 /
  data table_spin_states(:, 127) /  1,  1,  0,  0,  1, -1 /
  data table_spin_states(:, 128) /  1,  1,  0,  0,  1,  1 /
  data table_spin_states(:, 129) /  1,  1,  0,  1, -1, -1 /
  data table_spin_states(:, 130) /  1,  1,  0,  1, -1,  1 /
  data table_spin_states(:, 131) /  1,  1,  0,  1,  1, -1 /
  data table_spin_states(:, 132) /  1,  1,  0,  1,  1,  1 /
  data table_spin_states(:, 133) /  1,  1,  1, -1, -1, -1 /
  data table_spin_states(:, 134) /  1,  1,  1, -1, -1,  1 /
  data table_spin_states(:, 135) /  1,  1,  1, -1,  1, -1 /
  data table_spin_states(:, 136) /  1,  1,  1, -1,  1,  1 /
  data table_spin_states(:, 137) /  1,  1,  1,  0, -1, -1 /
  data table_spin_states(:, 138) /  1,  1,  1,  0, -1,  1 /
  data table_spin_states(:, 139) /  1,  1,  1,  0,  1, -1 /
  data table_spin_states(:, 140) /  1,  1,  1,  0,  1,  1 /
  data table_spin_states(:, 141) /  1,  1,  1,  1, -1, -1 /
  data table_spin_states(:, 142) /  1,  1,  1,  1, -1,  1 /
  data table_spin_states(:, 143) /  1,  1,  1,  1,  1, -1 /
  data table_spin_states(:, 144) /  1,  1,  1,  1,  1,  1 /

  integer, public :: nhel_max

  type(spinor) :: owf_t_4, owf_b_6, owf_e_2
  type(conjspinor) :: owf_t_3, owf_b_5, owf_e_1
  type(vector) :: owf_Wp_3, owf_Wm_4
  type(spinor) :: owf_wb_46
  type(conjspinor) :: owf_wb_35
  type(vector) :: owf_A_12, owf_Z_12

contains

  subroutine init (par, scheme)
    real(default), dimension(*), intent(in) :: par
    integer, intent(in) :: scheme
    call import_from_whizard (par, scheme)
  end subroutine init

  subroutine compute_production_owfs (p_ofs, hi, spins)
    type(momentum), intent(in), dimension(:) :: p_ofs
    integer, intent(in), optional :: hi
    integer, dimension(2), intent(in), optional :: spins
    integer, dimension(n_prt_OS) :: s_OS
    integer, dimension(2) :: s
    type(momentum) :: p12
    p12 = p_ofs(1) + p_ofs(2)
    if (process_mode == PROC_MODE_TT) then
       s_OS = table_spin_states_OS(:,hi)
       owf_e_1 = vbar (mass(11), - p_ofs(1), s_OS(1))
       owf_e_2 = u (mass(11), - p_ofs(2), s_OS(2))
       owf_t_3 = ubar (ttv_mtpole (p12 * p12), p_ofs(3), s_OS(3))
       owf_t_4 = v (ttv_mtpole (p12 * p12), p_ofs(4), s_OS(4))
       owf_A_12 = pr_feynman (p12, v_ff (qlep, owf_e_1, owf_e_2))
       if (.not. threshold%settings%Z_disabled) then
          owf_Z_12 = pr_unitarity (p12, mass(23), wd_tl (p12, width(23)), &
               .false., + va_ff (gnclep(1), gnclep(2), owf_e_1, owf_e_2))
       end if
    else
       if (present (hi)) then
          s = table_spin_states(1:2,hi)
       else
          if (present (spins)) then
             s = spins
          else
             call msg_fatal ("compute_production_owfs: " // &
                  "Please give either helicity index or spins")
          end if
       end if
       owf_e_1 = vbar (mass(11), - p_ofs(1), s(1))
       owf_e_2 = u (mass(11), - p_ofs(2), s(2))
       owf_A_12 = pr_feynman (p12, v_ff (qlep, owf_e_1, owf_e_2))
       if (.not. threshold%settings%Z_disabled) then
          owf_Z_12 = pr_unitarity (p12, mass(23), wd_tl (p12, width(23)), &
               .false., + va_ff (gnclep(1), gnclep(2), owf_e_1, owf_e_2))
       end if
    end if
  end subroutine compute_production_owfs

  subroutine compute_decay_owfs (p_ofs, hi, spins)
    type(momentum), dimension(:), intent(in) :: p_ofs
    integer, intent(in), optional :: hi
    integer, dimension(3:6), intent(in), optional :: spins
    integer, dimension(3:6) :: s
    if (present (hi)) then
       s = table_spin_states(3:6,hi)
    else
       if (present (spins)) then
          s = spins
       else
          call msg_fatal ("compute_decay_owfs: " // &
               "Please give either helicity index or spins")
       end if
    end if
    owf_Wp_3 = conjg (eps (mass(24), p_ofs(3), s(THR_POS_WP)))
    owf_Wm_4 = conjg (eps (mass(24), p_ofs(4), s(THR_POS_WM)))
    owf_b_5 = ubar (mass(5), p_ofs(5), s(THR_POS_B))
    owf_b_6 = v (mass(5), p_ofs(6), s(THR_POS_BBAR))
  end subroutine compute_decay_owfs

  function calculate_blob (ffi, p12, ptop_ofs, h_t, h_tbar, ptop_ons, tree_contrib) result (amp)
    complex(default) :: amp
    integer, intent(in) :: ffi
    type(momentum), intent(in) :: p12
    type(momentum), intent(in), dimension(2) :: ptop_ofs
    type(momentum), intent(in), dimension(2), optional :: ptop_ons
    integer, intent(in), optional :: h_t, h_tbar
    logical, intent(in), optional :: tree_contrib
    complex(default) :: blob_Z_vec, blob_Z_ax, ttv_vec, ttv_ax
    type(momentum) :: ptop, ptopbar
    real(default) :: mtop, top_width, add_tree
    integer :: u
    logical :: tc
    u = output_unit
    if (present (tree_contrib)) then
       tc = tree_contrib
    else
       tc = .not. (threshold%settings%interference .or. &
            threshold%settings%force_minus_one)
    end if
    if (tc) then
       add_tree = one
    else
       add_tree = zero
    end if
    if (process_mode == PROC_MODE_TT) then
       blob_Z_vec = gncup(1) * (ttv_formfactor (ptop_ofs(1), ptop_ofs(2), 1) + add_tree)
       blob_Z_ax = gncup(2) * (ttv_formfactor (ptop_ofs(1), ptop_ofs(2), 2) + add_tree)
       if (.not. threshold%settings%Z_disabled) then
          amp = owf_Z_12 * va_ff (blob_Z_vec, blob_Z_ax, owf_t_3, owf_t_4)
       else
          amp = zero
       end if
       amp = amp + owf_A_12 * v_ff (qup, owf_t_3, owf_t_4) * &
            (ttv_formfactor (ptop_ofs(1), ptop_ofs(2), 1) + add_tree)
    else if (process_mode == PROC_MODE_WBWB) then
       ttv_vec = ttv_formfactor (ptop_ofs(1), ptop_ofs(2), 1, ffi) + add_tree
       ttv_ax = ttv_formfactor (ptop_ofs(1), ptop_ofs(2), 2, ffi) + add_tree
       blob_Z_vec = gncup(1) * ttv_vec
       blob_Z_ax = gncup(2) * ttv_ax
       mtop = ttv_mtpole (p12 * p12)
       if (threshold%settings%factorized_computation) then
          if (threshold%settings%onshell_projection%production) then
             if (debug_active (D_THRESHOLD)) then
                call assert_equal (u, sqrt (ptop_ons(1) * ptop_ons(1)), &
                     mtop, "Production: ptop is projected", exit_on_fail=.true., &
                     rel_smallness=tiny_07)
                call assert_equal (u, sqrt (ptop_ons(2) * ptop_ons(2)), &
                     mtop, "Production: ptopbar is projected", exit_on_fail=.true., &
                     rel_smallness=tiny_07)
             end if
             ptop = ptop_ons(1)
             ptopbar = ptop_ons(2)
          else
             ptop = ptop_ofs(1)
             ptopbar = ptop_ofs(2)
          end if
          owf_t_3 = ubar (sqrt (ptop * ptop), ptop, h_t)
          owf_t_4 = v (sqrt (ptopbar * ptopbar), ptopbar, h_tbar)
          if (.not. threshold%settings%Z_disabled) then
             amp = owf_Z_12 * va_ff (blob_Z_vec, blob_Z_ax, owf_t_3, owf_t_4)
          else
             amp = zero
          end if
          amp = amp + owf_A_12 * v_ff (qup, owf_t_3, owf_t_4) * ttv_vec
       else
          top_width = ttv_wtpole (p12*p12, ffi)
          owf_wb_35 = pr_psibar (ptop_ofs(1), mtop, wd_tl (ptop_ofs(1), top_width), .false., &
               + f_fvl (gccq33, owf_b_5, owf_Wp_3))
          owf_wb_46 = pr_psi (ptop_ofs(2), mtop, wd_tl (ptop_ofs(2), top_width), .false., &
               + f_vlf (gccq33, owf_Wm_4, owf_b_6))
          if (.not. threshold%settings%Z_disabled) then
             amp = owf_Z_12 * va_ff (blob_Z_vec, blob_Z_ax, owf_wb_35, owf_wb_46)
          else
             amp = zero
          end if
          amp = amp + owf_A_12 * v_ff (qup, owf_wb_35, owf_wb_46) * ttv_vec
       end if
    else
       call msg_fatal ("Undefined process mode!")
    end if
  end function calculate_blob

  function top_propagators (ffi, p12, p_ofs) result(one_over_p)
    complex(default) :: one_over_p
    integer, intent(in) :: ffi
    type(momentum), intent(in) :: p12
    type(momentum), intent(in), dimension(2) :: p_ofs
    real(default) :: top_mass, top_width
    top_mass = ttv_mtpole (p12 * p12)
    if (threshold%settings%onshell_projection%width) then
      top_width = ttv_wtpole (p12*p12, ffi)
      one_over_p = one / cmplx (p_ofs(1) * p_ofs(1) - top_mass**2, &
           top_mass*top_width, kind=default)
      one_over_p = one_over_p / cmplx (p_ofs(2) * p_ofs(2) - top_mass**2, &
           top_mass * top_width, kind=default)
    else
      call msg_fatal ("on-shell projection width: You should really keep it on-shell. Not supported anymore.")
    end if
  end function top_propagators

  function top_decay_born (p_ons, h_t, h_W, h_b) result (me)
    complex(default) :: me
    type(momentum), intent(in), dimension(:) :: p_ons
    integer, intent(in) :: h_t
    integer, intent(in) :: h_W, h_b
    type(momentum) :: pwp, pb, ptop
    pwp = p_ons(THR_POS_WP)
    pb = p_ons(THR_POS_B)
    ptop = pwp + pb
    owf_Wp_3 = conjg (eps (mass(24), pwp, h_W))
    owf_b_5 = ubar (mass(5), pb, h_b)
    me = f_fvl (gccq33, owf_b_5, owf_Wp_3) * u (sqrt(ptop*ptop), ptop, h_t)
  end function top_decay_born

  function momentum_mode () result (mode)
    integer :: mode
    if (threshold%settings%onshell_projection%decay) then
      if (threshold%settings%onshell_projection%boost_decay) then
         mode = ONS_BOOST
      else
         mode = ONS
      end if
    else
       mode = OFS
    end if
  end function momentum_mode

  function anti_top_decay_born (p_ons, h_tbar, h_W, h_b) result (me)
    complex(default) :: me
    type(momentum), intent(in), dimension(:) :: p_ons
    integer, intent(in) :: h_tbar
    integer, intent(in) :: h_W, h_b
    type(momentum) :: pwm, pbbar, ptopbar
    pwm = p_ons(THR_POS_WM)
    pbbar = p_ons(THR_POS_BBAR)
    ptopbar = pwm + pbbar
    owf_Wm_4 = conjg (eps (mass(24), pwm, h_W))
    owf_b_6 = v (mass(5), pbbar, h_b)
    me = vbar (sqrt(ptopbar*ptopbar), ptopbar, h_tbar) * &
         f_vlf (gccq33, owf_Wm_4, owf_b_6)
  end function anti_top_decay_born

  subroutine compute_born (n_tot, p_ofs, ffi, sel_hel_beam, amp, tree_contrib)
    integer, intent(in) :: n_tot
    real(default), dimension(0:3,*), intent(in) :: p_ofs
    integer, intent(in) :: ffi, sel_hel_beam
    complex(default), dimension(:), intent(inout) :: amp
    logical, intent(in), optional :: tree_contrib
    complex(default), dimension(-1:1,-1:1,-1:1,-1:1) :: production_me
    complex(default), dimension(-1:1,-1:1,-1:1,1:2) :: born_decay_me
    complex(default) :: propagators
    integer :: hi
    type(momentum), dimension(2) :: ptop_ofs, ptop_ons
    type(momentum), dimension(:), allocatable :: mom_ofs, mom_ons, p_decay
    type(momentum) :: p12
    type(lorentz_transformation_t), dimension(2) :: lt
    integer, dimension(2) :: sel_hel
    logical :: this_helicity_selected
    amp = zero
    allocate (mom_ofs (n_tot), mom_ons (n_tot), p_decay(n_tot))
    call convert_to_mom_and_invert_sign (p_ofs, n_tot, mom_ofs)
    p12 = mom_ofs(1) + mom_ofs(2)
    ptop_ofs = get_top_momenta_offshell (p_ofs)
    if (threshold%settings%onshell_projection%active ()) &
         call compute_projections ()
    if (threshold%settings%factorized_computation) then
         propagators = top_propagators (ffi, p12, ptop_ofs)
         call compute_partial_matrix_elements ()
    end if
    sel_hel = get_selected_beam_helicities (sel_hel_beam)
    do hi = 1, nhel_max
       this_helicity_selected = sel_hel_beam < 0 .or. &
            (table_spin_states (1, hi) == sel_hel(1) &
            .and. table_spin_states (2, hi) == sel_hel(2))
       if (.not. this_helicity_selected) cycle
       if (threshold%settings%factorized_computation) then
          amp(hi) = multiply_partial_matrix_elements ( &
               production_me, born_decay_me, propagators, hi)
       else
          call compute_production_owfs (mom_ofs, hi)
          if (process_mode == PROC_MODE_WBWB) call compute_decay_owfs (mom_ofs, hi)
          amp(hi) = - calculate_blob (ffi, p12, ptop_ofs, tree_contrib=tree_contrib)
          if (process_mode == PROC_MODE_TT .and. &
               threshold%settings%helicity_approximation%ultra) then
             if (threshold%settings%sel_hel_top /= table_spin_states_OS(3,hi) .or. &
                threshold%settings%sel_hel_topbar /= table_spin_states_OS(4,hi)) then
                amp(hi) = zero
             end if
          end if
       end if
    end do
  contains
    subroutine compute_projections ()
      call compute_projected_momenta (mom_ofs, mom_ons)
      ptop_ons(1) = mom_ons(THR_POS_WP) + mom_ons(THR_POS_B)
      ptop_ons(2) = mom_ons(THR_POS_WM) + mom_ons(THR_POS_BBAR)
      if (.not. threshold%settings%onshell_projection%boost_decay) then
         lt = inverse (boost_to_cms (ptop_ons, ttv_mtpole (p12 * p12)))
      else
         lt = identity
      end if
    end subroutine compute_projections
    subroutine compute_partial_matrix_elements ()
      integer :: i
      production_me = compute_production_me (ffi, mom_ofs, ptop_ofs, ptop_ons, tree_contrib)
      select case (momentum_mode ())
      case (OFS)
         p_decay = mom_ofs
      case (ONS)
         p_decay(1:2) = mom_ons(1:2)
         do i = 3, n_tot
            p_decay(i) = apply_boost ((lt(ass_leg(i))), mom_ons(i))
         end do
      case (ONS_BOOST)
         p_decay = mom_ons
      end select
      if (.not. threshold%settings%onshell_projection%boost_decay &
             .and. threshold%settings%onshell_projection%decay &
             .and. debug2_active (D_THRESHOLD)) &
           call check_rest_frame (ttv_mtpole (p12 * p12), p_decay)
      born_decay_me = compute_decay_me (p_decay)
    end subroutine compute_partial_matrix_elements
  end subroutine compute_born

  function multiply_partial_matrix_elements ( &
         production_me, born_decay_me, propagators, hi) result (amp)
    complex(default) :: amp
    complex(default), dimension(-1:1,-1:1,-1:1,-1:1), intent(in) :: production_me
    complex(default), dimension(-1:1,-1:1,-1:1,1:2), intent(in) :: born_decay_me
    complex(default), intent(in) :: propagators
    integer, intent(in) :: hi
    complex(default) :: prod, dec1, dec2
    integer :: h_t, h_tbar
    integer, dimension(n_prt) :: s
    s = table_spin_states(:,hi)
    if (threshold%settings%helicity_approximation%ultra) then
       h_t = threshold%settings%sel_hel_top
       h_tbar = threshold%settings%sel_hel_topbar
       prod = production_me(s(1), s(2), h_t, h_tbar)
       dec1 = born_decay_me(s(ass_quark(2)), s(ass_boson(2)), h_tbar, 1)
       dec2 = born_decay_me(s(ass_quark(1)), s(ass_boson(1)), h_t, 2)
       amp = abs2 (prod) * abs2 (propagators) * abs2 (dec1) * abs2 (dec2)
    else if (threshold%settings%helicity_approximation%extra) then
       prod = zero
       do h_t = -1, 1, 2
          do h_tbar = -1, 1, 2
             prod = prod + abs2(production_me(s(1), s(2), h_t, h_tbar))
          end do
       end do
       dec1 = zero
       do h_tbar = -1, 1, 2
          dec1 = dec1 + abs2(born_decay_me(s(ass_quark(2)), s(ass_boson(2)), h_tbar, 1))
       end do
       dec2 = zero
       do h_t = -1, 1, 2
          dec2 = dec2 + abs2(born_decay_me(s(ass_quark(1)), s(ass_boson(1)), h_t, 2))
       end do
       amp = prod * abs2 (propagators) * dec1 * dec2 / 4
    else if (threshold%settings%helicity_approximation%simple) then
       amp = zero
       do h_t = -1, 1, 2
          do h_tbar = -1, 1, 2
             prod = production_me(s(1), s(2), h_t, h_tbar)
             dec1 = born_decay_me(s(ass_quark(2)), s(ass_boson(2)), h_tbar, 1)
             dec2 = born_decay_me(s(ass_quark(1)), s(ass_boson(1)), h_t, 2)
             amp = amp + abs2 (prod) * abs2 (propagators) * &
                  abs2 (dec1) * abs2 (dec2)
          end do
       end do
    else
       amp = zero
       do h_t = -1, 1, 2
          do h_tbar = -1, 1, 2
             prod = production_me(s(1), s(2), h_t, h_tbar)
             dec1 = born_decay_me(s(ass_quark(2)), s(ass_boson(2)), h_tbar, 1)
             dec2 = born_decay_me(s(ass_quark(1)), s(ass_boson(1)), h_t, 2)
             amp = amp + prod * propagators * dec1 * dec2
          end do
       end do
    end if
  end function multiply_partial_matrix_elements

  function ass_leg (i_particle)
    integer :: ass_leg
    integer, intent(in) :: i_particle
    if (i_particle == 3 .or. i_particle == 5) then
       ass_leg = 1
    else if (i_particle == 4 .or. i_particle == 6) then
       ass_leg = 2
    else
       call msg_fatal ("ass_leg called with invalid argument!")
    end if
  end function ass_leg

  subroutine check_rest_frame (mtop, p_decay)
    real(default), intent(in) :: mtop
    type(momentum), dimension(:), intent(in) :: p_decay
    integer :: u
    type(momentum) :: p_test
    u = output_unit
    p_test = p_decay (THR_POS_WP) + p_decay (THR_POS_B)
    call assert_equal (u, mtop, sqrt (p_test * p_test), &
         "Born phase-space: Check if top quark is in rest frame: ", abs_smallness = tiny_07, &
         rel_smallness = tiny_07, exit_on_fail = .true.)
    p_test = p_decay (THR_POS_WM) + p_decay (THR_POS_BBAR)
    call assert_equal (u, mtop, sqrt (p_test * p_test), &
         "Born phase-space: Check if anti-top quark is in rest frame: ", abs_smallness = tiny_07, &
         rel_smallness = tiny_07, exit_on_fail = .true.)
  end subroutine check_rest_frame

  function compute_decay_me (p_decay) result (born_decay_me)
    complex(default), dimension(-1:1,-1:1,-1:1,1:2) :: born_decay_me
    type(momentum), intent(in), dimension(:) :: p_decay
    procedure(top_decay_born), pointer :: top_decay_born_
    integer :: h_t, h_W, h_b, leg
    do leg = 1, 2
       if (leg == 1) then
          top_decay_born_ => anti_top_decay_born
       else
          top_decay_born_ => top_decay_born
       end if
       do h_t = -1, 1, 2
          do h_W = -1, 1, 1
             do h_b = -1, 1, 2
                born_decay_me(h_b, h_W, h_t, leg) = top_decay_born_ (p_decay, h_t, h_W, h_b)
             end do
          end do
       end do
    end do
  end function compute_decay_me

  function get_top_momenta_offshell (k, leg) result (p_top)
     type(momentum), dimension(2) :: p_top
     real(default), dimension(0:3,*), intent(in) :: k
     integer, intent(in), optional :: leg
     type(vector4_t), dimension(2) :: p_tmp
     if (process_mode /= PROC_MODE_WBWB) then
        p_tmp(1)%p = k(:,3)
        p_tmp(2)%p = k(:,4)
     else
        p_tmp(1)%p = k(:,THR_POS_WP) + k(:,THR_POS_B)
        p_tmp(2)%p = k(:,THR_POS_WM) + k(:,THR_POS_BBAR)
        if (present (leg)) then
           if (leg == 1) then
              p_tmp(1)%p = p_tmp(1)%p + k(:,THR_POS_GLUON)
           else
              p_tmp(2)%p = p_tmp(2)%p + k(:,THR_POS_GLUON)
           end if
        end if
     end if
     p_top(1) = p_tmp(1)%p
     p_top(2) = p_tmp(2)%p
  end function get_top_momenta_offshell

  subroutine convert_to_mom_and_invert_sign (p, n, mom)
    real(default), intent(in), dimension(0:3,*) :: p
    integer, intent(in) :: n
    type(momentum), intent(out), dimension(:), allocatable :: mom
    integer :: i
    allocate (mom(n))
    do i = 1, n
       if (i <= 2) then
          mom(i) = - p(:,i)
       else
          mom(i) = p(:,i)
       end if
    end do
  end subroutine convert_to_mom_and_invert_sign

  subroutine compute_projected_momenta (p_ofs, p_ons)
     type(momentum), intent(in), dimension(:) :: p_ofs
     type(momentum), intent(out), dimension(:) :: p_ons
     type(momentum), dimension(2) :: ptop_ons
     real(default), dimension(4) :: tmp, test
     type(momentum) :: p12, p35
     type(lorentz_transformation_t) :: lt
     if (threshold%settings%onshell_projection%active ()) then
        p12 = p_ofs(1) + p_ofs(2); p35 = p_ofs(3) + p_ofs(5)
        call compute_projected_top_momenta (p12, p35, ptop_ons, lt)
        call compute_projected_top_decay_products (p12, lt, p_ofs, p_ons)
        if (debug_active (D_THRESHOLD)) then
           if (- p12%t > 2 * ttv_mtpole (p12*p12)) then
              !!! No sum (...) function for type(momentum), need to do this explicitly
              tmp = p_ofs(THR_POS_WP) + p_ofs(THR_POS_WM) + p_ofs(THR_POS_B) + p_ofs(THR_POS_BBAR)
              test = - p12
              call assert_equal (output_unit, tmp, test, &
                   "overall: momentum conservation", &
                   abs_smallness=10*tiny_07, &
                   rel_smallness=tiny_07, &
                   exit_on_fail=.true.)
           end if
        end if
     end if
  end subroutine compute_projected_momenta

  elemental function boost_to_cms (ptop_ons, mtop) result (lt)
    type(lorentz_transformation_t) :: lt
    type(momentum), intent(in) :: ptop_ons
    real(default), intent(in) :: mtop
    real(default), dimension(4) :: tmp
    type(vector4_t) :: v4_tmp
    !!! There is no direct conversion between type(momentum) and type(vector4_t)
    !!! So, we need a temporary real array to store the momentum values
    tmp = ptop_ons; v4_tmp = tmp
    lt = boost (v4_tmp, mtop)
  end function boost_to_cms

  !!! The boost to the center-of-mass system only has a reasonable meaning
  !!! above the threshold. Below the threshold, we do not apply boost at all, so
  !!! that the top quarks stay in the rest frame. However, with top quarks exactly
  !!! at rest, problems arise in the matrix elements (e.g. in the computation
  !!! of angles). Therefore, we apply a boost which is not exactly 1, but has a
  !!! tiny value differing from that.

  subroutine compute_projected_top_momenta (p12, p35, ptop_ons, lt)
    type(momentum), intent(in) :: p12, p35
    type(momentum), intent(out), dimension(2) :: ptop_ons
    type(momentum), dimension(2) :: ptop_ons_rest
    type(lorentz_transformation_t), intent(out) :: lt
    real(default) :: sqrts, scale_factor, mtop, s_minus_threshold
    real(default), dimension(1:3) :: unit_vec
    real(default), dimension(4) :: tmp, test
    integer :: u
    mtop = ttv_mtpole (p12 * p12)
    sqrts = - p12%t
    s_minus_threshold = sqrts**2 - four * mtop**2
    if (s_minus_threshold > zero) then
       scale_factor = sqrt (s_minus_threshold) / 2
    else
       scale_factor = tiny_10
    end if
    unit_vec = p35%x / sqrt (dot_product(p35%x, p35%x))
    ptop_ons(1) = [sqrts / 2, scale_factor * unit_vec]
    ptop_ons(2) = [sqrts / 2, - scale_factor * unit_vec]
    ptop_ons_rest(1) = [mtop, zero, zero, zero]
    ptop_ons_rest(2) = [mtop, zero, zero, zero]
    lt = boost_to_cms (ptop_ons(1), mtop)
    if (debug_active (D_THRESHOLD)) then
       u = output_unit
       if (s_minus_threshold > zero) then
          tmp = apply_boost (inverse (lt), ptop_ons(1))
          test = ptop_ons_rest(1)
          call assert_equal (u, tmp, test, &
               "verify that we have the right boost", abs_smallness=tiny_07 * 10, &
               rel_smallness=tiny_07, exit_on_fail=.true.)
          tmp = apply_boost (lt, ptop_ons_rest(1))
          test = ptop_ons(1)
          call assert_equal(u, test, tmp, "test the inverse boost", &
               exit_on_fail=.true.)
          call assert (u, p12 == - (ptop_ons(1) + ptop_ons(2)), &
               "momentum conservation with a flip", exit_on_fail=.true.)
          call assert_equal (u, ptop_ons(1) * ptop_ons(1), mtop**2, &
               "mass onshell", abs_smallness=tiny_07, &
                rel_smallness=tiny_07, exit_on_fail=.true.)
          call assert_equal (u, ptop_ons(2) * ptop_ons(2), mtop**2, &
               "mass onshell", abs_smallness=tiny_07, &
                rel_smallness=tiny_07, exit_on_fail=.true.)
          call assert_equal (u, dot_product(unit_vec, unit_vec), one, &
              "unit vector length", exit_on_fail=.true.)
       end if
    end if
  end subroutine compute_projected_top_momenta

  subroutine compute_projected_top_decay_products (p12, lt, p_ofs, p_ons)
    type(momentum), intent(in) :: p12
    type(lorentz_transformation_t), intent(in) :: lt
    type(momentum), intent(in), dimension(:) :: p_ofs
    type(momentum), intent(out), dimension(*) :: p_ons
    type(momentum), dimension(3:6) :: p_ons_rest
    real(default) :: mtop, mw2, mb2, en_w, en_b
    type(vector4_t) :: p_tmp_1, p_tmp_2
    type(vector4_t), dimension(3) :: p_decay
    integer :: u
    logical :: momenta_already_onshell
    u = output_unit
    mtop = ttv_mtpole (p12*p12)
    mw2 = mass(24)**2
    mb2 = mass(5)**2
    en_w = (mtop**2 + mw2 - mb2) / (2 * mtop)
    en_b = (mtop**2 - mw2 + mb2) / (2 * mtop)
    momenta_already_onshell = process_mode == PROC_MODE_WBWB
    if (process_mode == PROC_MODE_TT) then
       momenta_already_onshell = .true.
    else
       momenta_already_onshell = &
            onshell_tops (p_ofs, [THR_POS_WP, THR_POS_B], [THR_POS_WM, THR_POS_BBAR])
    end if
    p_ons(1:2) = p_ofs(1:2)
    if (.not. momenta_already_onshell) then
       p_tmp_1%p = apply_boost (inverse (lt), p_ofs(THR_POS_B))
       p_tmp_2%p = apply_boost (inverse (lt), p_ofs(THR_POS_WP))
       p_decay = create_two_particle_decay (mtop**2, p_tmp_1, p_tmp_2)
       p_ons_rest(THR_POS_B) = p_decay(2)%p
       p_ons_rest(THR_POS_WP) = p_decay(3)%p
       p_ons(THR_POS_WP) = apply_boost (lt, p_ons_rest(THR_POS_WP))
       p_ons(THR_POS_B) = apply_boost (lt, p_ons_rest(THR_POS_B))
       p_tmp_1%p = apply_boost (lt, p_ofs(THR_POS_BBAR))
       p_tmp_2%p = apply_boost (lt, p_ofs(THR_POS_WM))
       p_decay = create_two_particle_decay (mtop**2, p_tmp_1, p_tmp_2)
       p_ons_rest(THR_POS_BBAR) = p_decay(2)%p
       p_ons_rest(THR_POS_WM) = p_decay(3)%p
       p_ons(THR_POS_WM) = apply_boost (inverse(lt), p_ons_rest(THR_POS_WM))
       p_ons(THR_POS_BBAR) = apply_boost (inverse(lt), p_ons_rest(THR_POS_BBAR))
    else
       p_ons(3:6) = p_ofs(3:6)
    end if
    if (debug_active (D_THRESHOLD) .and. .not. momenta_already_onshell) then
       call assert_equal (u, en_w + en_b, mtop, "top energy", &
            exit_on_fail=.true.)
       call assert_equal (u, en_w + en_b, mtop, "top energy", &
          exit_on_fail=.true.)
       call assert_equal (u, p_ons_rest(THR_POS_WM) * p_ons_rest(THR_POS_WM), mw2, &
            "W- mass onshell", rel_smallness=tiny_07, &
            exit_on_fail=.true.)
       call assert_equal (u, p_ons_rest(THR_POS_BBAR) * p_ons_rest(THR_POS_BBAR), mb2, &
            "bbar mass onshell", rel_smallness=tiny_07, &
            exit_on_fail=.true.)
       call assert_equal (u, p_ons_rest(THR_POS_WP) * p_ons_rest(THR_POS_WP), mw2, &
            "W+ mass onshell", rel_smallness=tiny_07, &
            exit_on_fail=.true.)
       call assert_equal (u, p_ons_rest(THR_POS_B) * p_ons_rest(THR_POS_B), mb2, &
            "b mass onshell", rel_smallness=tiny_07, &
            exit_on_fail=.true.)
    end if
  end subroutine compute_projected_top_decay_products

  elemental function apply_boost (boost_in, mom) result (mom_result)
    type(momentum) :: mom_result
    type(momentum), intent(in) :: mom
    type(lorentz_transformation_t), intent(in) :: boost_in
    type(vector4_t) :: tmp_v4
    real(default), dimension(4) :: tmp
    tmp = mom
    tmp_v4 = tmp
    tmp = boost_in * tmp_v4
    mom_result = tmp
  end function apply_boost

  pure function check_if_onshell (p1, p2, m) result (onshell)
    logical :: onshell
    type(momentum), intent(in) :: p1, p2
    real(default), intent(in) :: m
    type(momentum) :: pp
    real(default) :: mm
    pp = p1 + p2
    mm = sqrt (pp * pp)
    onshell = nearly_equal (mm, m)
  end function check_if_onshell

  function compute_production_me (ffi, p_ofs, ptop_ofs, ptop_ons, &
         tree_contrib) result (production_me)
    complex(default), dimension(-1:1,-1:1,-1:1,-1:1) :: production_me
    integer, intent(in) :: ffi
    type(momentum), intent(in), dimension(:) :: p_ofs
    type(momentum), intent(in), dimension(2) :: ptop_ofs
    type(momentum), intent(in), dimension(2), optional :: ptop_ons
    logical, intent(in), optional :: tree_contrib
    type(momentum) :: p12
    integer :: h_t, h_tbar, h_pos, h_el
    p12 = p_ofs(1) + p_ofs(2)
    do h_tbar = -1, 1, 2
    do h_t = -1, 1, 2
    do h_pos = -1, 1, 2
    do h_el = -1, 1, 2
       call compute_production_owfs (p_ofs, spins = [h_el, h_pos])
       production_me(h_el, h_pos, h_t, h_tbar) = &
            calculate_blob (ffi, p12, ptop_ofs, h_t, h_tbar, ptop_ons, tree_contrib=tree_contrib)
    end do
    end do
    end do
    end do
  end function compute_production_me

  subroutine check_for_consistent_flags_for_real (settings)
    type(settings_t), intent(in) :: settings
    if (.not. settings%factorized_computation)  &
         call msg_fatal ('compute_real: OFFSHELL_STRATEGY is not '&
         &'factorized (activate with 2')
    if (.not. (settings%helicity_approximation%simple .or. &
              settings%helicity_approximation%ultra)) &
         call msg_fatal ('compute_real: OFFSHELL_STRATEGY is not '&
         &'helicity-approximated (activate with 32)')
  end subroutine check_for_consistent_flags_for_real

  function skip (settings, h_t, h_tbar)
    logical :: skip
    type(settings_t), intent(in) :: settings
    integer, intent(in) :: h_t, h_tbar
    skip = settings%helicity_approximation%ultra .and. &
         (h_t /= settings%sel_hel_top .or. h_tbar /= settings%sel_hel_topbar)
  end function skip

  function compute_real (n_tot, p_ofs, p_ons, leg, ffi, sel_hel_beam) result (amp2)
    real(default) :: amp2
    integer, intent(in) :: n_tot
    real(default), dimension(0:3,*), intent(in) :: p_ofs
    real(default), dimension(0:3,*), intent(in) :: p_ons
    integer, intent(in) :: leg, ffi, sel_hel_beam
    complex(default), dimension(-1:1,-1:1,-1:1,-1:1) :: production_me
    complex(default), dimension(-1:1,-1:1,-1:1,-1:1) :: real_decay_me
    complex(default), dimension(-1:1,-1:1,-1:1) :: born_decay_me
    complex(default) :: top_propagators_
    complex(default) :: real_dec, born_dec, prod
    real(default) :: total
    integer, dimension(2) :: h_ass_t
    integer, dimension(n_prt) :: s
    integer :: hi, h_t, h_tbar, h_gl
    type(momentum), dimension(:), allocatable :: mom_ofs, mom_ons
    integer :: other_leg
    integer, dimension(2) :: sel_hel
    logical :: this_helicity_selected
    call check_for_consistent_flags_for_real (threshold%settings)
    call convert_to_mom_and_invert_sign (p_ofs, n_tot, mom_ofs)
    call convert_to_mom_and_invert_sign (p_ons, n_tot, mom_ons)
    call compute_real_amplitudes (ffi, mom_ofs, mom_ons, leg, &
         production_me, real_decay_me, born_decay_me, top_propagators_)
    total = zero
    sel_hel = get_selected_beam_helicities (sel_hel_beam)
    do hi = 1, nhel_max
       s = table_spin_states(:,hi)
       this_helicity_selected = sel_hel_beam < 0 .or. &
            (s(1) == sel_hel(1) .and. s(2) == sel_hel(2))
       if (.not. this_helicity_selected) cycle
       do h_t = -1, 1, 2
       do h_tbar = -1, 1, 2
          if (skip (threshold%settings, h_t, h_tbar)) cycle
          h_ass_t = [h_t, h_tbar]
          other_leg = 3 - leg
          prod = production_me(s(1), s(2), h_t, h_tbar)
          born_dec = born_decay_me(s(ass_quark(other_leg)), &
               s(ass_boson(other_leg)), h_ass_t(other_leg))
          do h_gl = -1, 1, 2
             real_dec = real_decay_me(h_gl, s(ass_quark(leg)), &
                  s(ass_boson(leg)), h_ass_t(leg))
             total = total + abs2 (prod) * abs2 (real_dec) * abs2 (born_dec) * &
                  abs2(top_propagators_)
          end do
       end do
       end do
    end do
    !!! Add color factor. Real ~ N^2 - 1; Born(has to be divided out) ~ N
    amp2 = total * (N_**2 - one) / N_

  end function compute_real

  subroutine compute_real_amplitudes (ffi, p_ofs, p_ons, leg, &
       production_me, real_decay_me, born_decay_me, top_propagators_)
    integer, intent(in) :: ffi
    type(momentum), intent(in), dimension(:) :: p_ofs, p_ons
    integer, intent(in) :: leg
    complex(default), intent(out), dimension(-1:1,-1:1,-1:1,-1:1) :: production_me
    complex(default), intent(out), dimension(-1:1,-1:1,-1:1,-1:1) :: real_decay_me
    complex(default), intent(out), dimension(-1:1,-1:1,-1:1) :: born_decay_me
    complex(default), intent(out) :: top_propagators_
    procedure(top_real_decay_calculate_amplitude), pointer :: top_decay_real
    procedure(top_decay_born), pointer :: top_decay_born_
    type(momentum), dimension(2) :: ptop_ofs, ptop_ons
    type(momentum) :: p12
    type(momentum), dimension(7) :: p_ons_work
    type(momentum), dimension(4) :: p_ons_real_decay
    type(lorentz_transformation_t) :: lt
    integer :: h_t, h_b, h_W, h_gl
    real(default) :: mtop
    integer :: other_leg
    p12 = p_ofs(1) + p_ofs(2)
    ptop_ons(1) = p_ons(THR_POS_WP) + p_ons(THR_POS_B)
    ptop_ons(2) = p_ons(THR_POS_WM) + p_ons(THR_POS_BBAR)
    ptop_ons(leg) = ptop_ons(leg) + p_ons(THR_POS_GLUON)
    ptop_ofs(1) = p_ofs(THR_POS_WP) + p_ofs(THR_POS_B)
    ptop_ofs(2) = p_ofs(THR_POS_WM) + p_ofs(THR_POS_BBAR)
    ptop_ofs(leg) = ptop_ofs(leg) + p_ofs(THR_POS_GLUON)
    production_me = compute_production_me (ffi, p_ofs, ptop_ofs, ptop_ons)
    if (leg == 1) then
       top_decay_real => top_real_decay_calculate_amplitude
       top_decay_born_ => anti_top_decay_born
    else
       top_decay_real => anti_top_real_decay_calculate_amplitude
       top_decay_born_ => top_decay_born
    end if
    p_ons_real_decay = [ptop_ons(leg), p_ons(ass_boson(leg)), &
         p_ons(ass_quark(leg)), p_ons(THR_POS_GLUON)]
    mtop = ttv_mtpole (p12 * p12)
    if (.not. threshold%settings%onshell_projection%boost_decay) then
       p_ons_work([1,2]) = p_ons([1,2])
       lt = inverse (boost_to_cms (ptop_ons(leg), mtop))
       p_ons_real_decay = apply_boost (lt, p_ons_real_decay)
       p_ons_work([ass_boson(leg), ass_quark(leg), THR_POS_GLUON]) = &
               apply_boost (lt, p_ons([ass_boson(leg), ass_quark(leg), THR_POS_GLUON]))
       other_leg = 3 - leg
       p_ons_work([ass_boson(other_leg), ass_quark(other_leg)]) = &
               apply_boost (inverse (lt), p_ons([ass_boson(other_leg), ass_quark(other_leg)]))
       if (debug2_active (D_THRESHOLd)) call check_rest_frame ()
    else
       p_ons_work = p_ons
    end if
    if (debug2_active (D_THRESHOLD)) call check_phase_space_point ()
    do h_t = -1, 1, 2
    do h_W = -1, 1, 1
    do h_b = -1, 1, 2
       born_decay_me(h_b, h_W, h_t) = top_decay_born_ (p_ons_work, h_t, h_W, h_b)
       if (.not. test_ward) then
          do h_gl = -1, 1, 2
             real_decay_me(h_gl, h_b, h_W, h_t) = top_decay_real &
                  (p_ons_real_decay, [h_t, h_W, h_b, h_gl], zero)
          end do
       else
          do h_gl = -1, 1, 2
             real_decay_me(h_gl, h_b, h_W, h_t) = top_decay_real &
                  (p_ons_real_decay, [h_t, h_W, h_b, 4], zero)
          end do
       end if
    end do
    end do
    end do
    top_propagators_ = top_propagators (ffi, p12, ptop_ofs)
  contains
    subroutine check_rest_frame ()
      integer :: u, other_leg
      type(momentum) :: p_test
      u = output_unit
      call assert_equal (u, mtop, p_ons_real_decay(1)%t, &
           'Real top is in rest frame?', abs_smallness = tiny_07, &
           rel_smallness = tiny_07, exit_on_fail = .true.)
      other_leg = 3 - leg
      p_test = p_ons_work (ass_boson(other_leg)) + p_ons_work (ass_quark(other_leg))
      call assert_equal (u, mtop, p_test%t, 'Born top is in rest frame?', &
           abs_smallness = tiny_07, rel_smallness = tiny_07, exit_on_fail = .true.)
    end subroutine check_rest_frame

    subroutine check_phase_space_point ()
      real(default) :: sqrts, E
      integer :: u, i, other_leg
      type(momentum) :: p_test
      if (debug_active (D_THRESHOLD)) then
         u = output_unit
         sqrts = sqrt(p12 * p12)
         if (threshold%settings%onshell_projection%decay) then
           call assert_equal (u, mtop, &
                sqrt (p_ons_real_decay(1) * p_ons_real_decay(1)), 'Real-decay Top is on-shell', &
                abs_smallness = tiny_07, rel_smallness = tiny_07, exit_on_fail = .true.)
           if (threshold%settings%onshell_projection%boost_decay) &
                call assert_equal (u, sqrts / two, sum (p_ons_real_decay(2:4)%t), &
                     'Real decay momenta have E = sqrts / 2', abs_smallness = tiny_07, &
                     rel_smallness = tiny_07, exit_on_fail = .true.)
         end if
         !!! Use absoulate value to handle slightly negative values
         call assert_equal (u, zero, sqrt (abs (p_ons_real_decay(4)* p_ons_real_decay(4))), &
              'Gluon is on-shell', abs_smallness = 1E-5_default, &
              rel_smallness = 1E-5_default, exit_on_fail = .true.)
         call assert_equal (u, mass(5), sqrt (p_ons_real_decay(3) * p_ons_real_decay(3)), &
              'Real-decay Bottom is on-shell', abs_smallness = tiny_07, &
              rel_smallness = tiny_07, exit_on_fail = .true.)
         call assert_equal (u, mass(24), sqrt (p_ons_real_decay(2) * p_ons_real_decay(2)), &
              'Real-decay W is on-shell', abs_smallness = tiny_07, &
              rel_smallness = tiny_07, exit_on_fail = .true.)
         other_leg = 3 - leg
         if (threshold%settings%onshell_projection%decay) then
              p_test = p_ons_work(ass_quark(other_leg)) + p_ons_work(ass_boson(other_leg))
              call assert_equal (u, mtop, sqrt (p_test * p_test), &
                   'Born-decay Top is on-shell', abs_smallness = tiny_07, &
                   rel_smallness = tiny_07, exit_on_fail = .true.)
         end if
         p_test = p_ons_work(ass_quark(other_leg))
         call assert_equal (u, mass(5), sqrt (p_test * p_test), &
              'Born-decay Bottom is on-shell', abs_smallness = tiny_07, &
              rel_smallness = tiny_07, exit_on_fail = .true.)
         p_test = p_ons_work(ass_boson(other_leg))
         call assert_equal (u, mass(24), sqrt (p_test * p_test), &
              'Born-decy W is on-shell', abs_smallness = tiny_07, &
              rel_smallness = tiny_07, exit_on_fail = .true.)
         do i = 1, 3
            E = sum (p_ons_real_decay (2:4)%x(i)) + p_ons_work(ass_quark(other_leg))%x(i) &
                 + p_ons_work(ass_boson(other_leg))%x(i)
            call assert_equal (u, zero, E, 'Total momentum vanishes', &
                 abs_smallness = 1E-5_default, rel_smallness = 1E-5_default, &
                 exit_on_fail = .true.)
         end do
      end if
    end subroutine check_phase_space_point

  end subroutine compute_real_amplitudes

  subroutine handle_onshell_test_point (n_tot, p_ofs, p_ofs_out)
    integer, intent(in) :: n_tot
    real(c_default_float), dimension(0:3,*), intent(in) :: p_ofs
    real(c_default_float), dimension(0:3,n_tot), intent(out) :: p_ofs_out
    type(spinor) :: test_psi, test_spinor1, test_spinor2
    type(conjspinor) :: test_psibar, test_conjspinor1, test_conjspinor2
    type(momentum) :: p35, p46, p3, p4, p5, p6
    type(vector) :: vp
    complex(kind=default) :: c_one
    real(default) :: mtop
    integer :: i
    p_ofs_out(:,1:n_tot) = p_ofs(:,1:n_tot)
    if (test_onshell) then
       !!! This is above threshold: sqrts = 350 GeV
       p_ofs_out(:,3) =  [  95.237818224532234, -39.152407665077284, &
                            28.564227360034039,  15.943661703665313]
       p_ofs_out(:,4) =  [ 101.19883586743811,   39.626666946167894, &
                           -10.149946667584741,  45.833335786383174]
       p_ofs_out(:,5) =  [ 79.762181775467766,   66.551897152939105, &
                          -28.767162575206349,  -32.979705643001502]
       p_ofs_out(:,6) =  [ 73.801164132561894,  -67.026156434029716, &
                           10.352881882757050,  -28.797291847046992]
       c_one = one
       mtop = mass(6)     !!! We assume a fixed mtpole = m1S
       p3 = p_ofs_out(:,3)
       p4 = p_ofs_out(:,4)
       p5 = p_ofs_out(:,5)
       p6 = p_ofs_out(:,6)
       p35 = p_ofs_out(:,3) + p_ofs_out(:,5)
       p46 = p_ofs_out(:,4) + p_ofs_out(:,6)
       call check_spinor_sum ()
    end if
  contains
    subroutine check_spinor_sum ()
      integer :: s3, s4, s5, s6
      vp = p35
      test_psi%a = [one, two, three, four]
      test_spinor1 = f_vf (c_one, vp, test_psi) + mtop * test_psi
      test_spinor2 = u (mtop, p35, +1) * (ubar (mtop, p35, +1) * test_psi) + &
                     u (mtop, p35, -1) * (ubar (mtop, p35, -1) * test_psi)
      do i = 1, 4
        call assert_equal (output_unit, test_spinor1%a(i), test_spinor2%a(i), &
             "(p+m)1=(sum u ubar)1", abs_smallness = tiny_07, &
             rel_smallness = tiny_07, exit_on_fail = .true.)
      end do
      !!! top
      do s3 = -1, 1 ; do s5 = -1, 1, 2
         test_psibar = + f_fvl(gccq33,ubar (mass(5), p5, s5),conjg (eps (mass(24), p3, s3)))
         test_conjspinor1 = pr_psibar (p35,mtop,wd_tl(p35,width(6)),.false., test_psibar)
         test_conjspinor2 = (test_psibar * u (mtop, p35, +1)) * ubar (mtop, p35, +1) + &
                            (test_psibar * u (mtop, p35, -1)) * ubar (mtop, p35, -1)
         test_conjspinor2 = test_conjspinor2 * (one / cmplx (p35*p35 - mtop**2, &
              mtop*width(6), kind=default))
         do i = 1, 4
            call assert_equal (output_unit, test_conjspinor1%a(i), test_conjspinor2%a(i), &
                 "(p+m)/(p^2-m^2)=(sum u ubar)/(p^2-m^2)", abs_smallness = tiny_07, &
                 rel_smallness = tiny_07, exit_on_fail = .true.)
         end do
      end do; end do
      !!! topbar
      do s4 = -1, 1 ; do s6 = -1, 1, 2
         test_psi = + f_vlf(gccq33,conjg (eps (mass(24), p4, s4)), v (mass(5), p6, s6))
         test_spinor1 = - pr_psi (p46,mtop,wd_tl(p46,width(6)),.false., test_psi)
         test_spinor2 = v (mtop, p46, +1) * (vbar (mtop, p46, +1) * test_psi) + &
                        v (mtop, p46, -1) * (vbar (mtop, p46, -1) * test_psi)
         test_spinor2 = test_spinor2 * (one / cmplx (p46*p46 - mtop**2, &
              mtop*width(6), kind=default))
         do i = 1, 4
           call assert_equal (output_unit, test_spinor1%a(i), test_spinor2%a(i), &
                "- (-p+m)/(p^2-m^2)=(sum v vbar)/(p^2-m^2)", abs_smallness = tiny_07, &
                rel_smallness = tiny_07, exit_on_fail = .true.)
         end do
      end do; end do
    end subroutine check_spinor_sum
  end subroutine handle_onshell_test_point

  function get_selected_beam_helicities (sel_hel_beam) result (sel_hel)
    integer, dimension(2) :: sel_hel
    integer, intent(in) :: sel_hel_beam
    integer :: i
    if (sel_hel_beam >= 0) then
       do i = 0, 1
          if (btest (sel_hel_beam, i)) then
             sel_hel(i+1) = 1
          else
             sel_hel(i+1) = -1
          end if
       end do
    else
       sel_hel = 0
    end if
  end function get_selected_beam_helicities


  subroutine set_production_factors (sel_hel_beam)
    integer, intent(in) :: sel_hel_beam
    !!! Colour factors: N_ colors can be produced
    if (sel_hel_beam >= 0) then
       production_factors = N_
    else
       !!! Helicity factors: Average over incoming helicities
       production_factors = N_ / four
    end if
  end subroutine

end module @ID@_threshold


subroutine @ID@_threshold_init (par, scheme) bind(C)
  use iso_c_binding
  use kinds
  use @ID@_threshold
  implicit none
  real(c_default_float), dimension(*), intent(in) :: par
  integer(c_int), intent(in) :: scheme
  call init (par, scheme)
end subroutine @ID@_threshold_init

subroutine @ID@_set_process_mode (mode) bind(C)
  use iso_c_binding
  use kinds
  use @ID@_threshold
  use physics_defs, only: PROC_MODE_TT
  implicit none
  integer, intent(in) :: mode
  process_mode = mode
  if (process_mode == PROC_MODE_TT) then
     nhel_max = n_hel_OS
  else
     nhel_max = n_hel
  end if
end subroutine @ID@_set_process_mode

!!! p_ons is supplied correctly for the real computation
subroutine @ID@_get_amp_squared (amp2, p_ofs, p_ons, leg, n_tot, sel_hel_beam) bind(C)
  use iso_c_binding
  use kinds
  use constants
  use numeric_utils
  use diagnostics
  use physics_defs, only: THR_POS_WP, THR_POS_WM, THR_POS_B, THR_POS_BBAR
  use physics_defs, only: PROC_MODE_WBWB
  use, intrinsic :: iso_fortran_env, only: output_unit
  use opr_@ID@, full_proc_new_event => new_event
  use opr_@ID@, full_proc_get_amplitude => get_amplitude
  use opr_@ID@, full_proc_number_spin_states => number_spin_states
  use opr_@ID@, full_proc_number_particles_out => number_particles_out
  use opr_@ID@, full_proc_flavor_states => flavor_states
  use @ID@_threshold
  use parameters_sm_tt_threshold
  use ttv_formfactors
  implicit none
  real(c_default_float), intent(out) :: amp2
  real(c_default_float), dimension(0:3,*), intent(in) :: p_ofs
  real(c_default_float), dimension(0:3,*), intent(in) :: p_ons
  integer, intent(in) :: leg, n_tot, sel_hel_beam
  complex(default), dimension(:), save, allocatable :: amp_with_FF, &
       amp_no_FF, amp_omega_full
  logical :: real_computation
  integer :: hi, n_total_hel
  real_computation = full_proc_number_particles_out () == 5
  if (.not. allocated (amp_omega_full)) then
     if (real_computation) then
        n_total_hel = n_hel * 2 ! times 2 due to the gluon
     else
        n_total_hel = n_hel
     end if
     call allocate_amps ()
  end if
  call set_production_factors (sel_hel_beam)
  if (real_computation) then
     amp2 = compute_real (n_tot, p_ofs, p_ons, leg, FF, sel_hel_beam)
  else
     amp2 = compute_born_special_cases (n_tot, p_ofs, leg, sel_hel_beam)
  end if
  amp2 = amp2 * production_factors

contains

  subroutine allocate_amps ()
    !!! We use the `save` attribute to avoid reallocation for every PS point
    allocate (amp_omega_full(n_total_hel))
    allocate (amp_with_FF(n_total_hel))
    allocate (amp_no_FF(n_total_hel))
    amp_omega_full = zero
    amp_with_FF = zero
    amp_no_FF = zero
  end subroutine allocate_amps

  subroutine handle_test_onshell (n_tot, p_ofs_work)
    integer, intent(in) :: n_tot
    real(c_default_float), dimension(0:3,n_tot) :: p_ofs_work
    if (test_onshell) then
       call compute_born (n_tot, p_ofs_work, TREE, -1, amp_no_FF, tree_contrib=.true.)
       do hi = 1, size(amp_omega_full)
          call assert_equal (output_unit, amp_omega_full(hi), &
               amp_no_FF(hi), "Signal \= Factorized", exit_on_fail=.true., &
               rel_smallness=tiny_07)
       end do
       stop
    end if
  end subroutine handle_test_onshell

  function compute_born_special_cases (n_tot, p_ofs, leg, sel_hel_beam) result (amp2)
    real(c_default_float) :: amp2
    integer, intent(in) :: n_tot, leg
    integer, intent(in) :: sel_hel_beam
    real(c_default_float), dimension(0:3,*), intent(in) :: p_ofs
    real(c_default_float), dimension(0:3,n_tot) :: p_ofs_work
    integer, dimension(6,1) :: table_flavor_states
    integer :: use_ff
    logical :: simple, extra, ultra
    call handle_onshell_test_point (n_tot, p_ofs, p_ofs_work)
    if (threshold%settings%interference) then
       call threshold%formfactor%disable ()
       call full_proc_new_event (p_ofs_work)
       do hi = 1, full_proc_number_spin_states()
          amp_omega_full(hi) = full_proc_get_amplitude (1, hi, 1)
       end do
       if (process_mode == PROC_MODE_WBWB) then
          call full_proc_flavor_states (table_flavor_states)
          if (table_flavor_states(1,1) /= -11 .or. &
              table_flavor_states(2,1) /= +11 .or. &
              table_flavor_states(THR_POS_WP,1) /= +24 .or. &
              table_flavor_states(THR_POS_WM,1) /= -24 .or. &
              table_flavor_states(THR_POS_B,1) /= +5 .or. &
              table_flavor_states(THR_POS_BBAR,1) /= -5) then
             call msg_fatal ("The factorized computation requires " // &
                  "'E1, e1 => Wp, Wm, b, B' as final state (in this order!)")
          end if
       end if
    end if
    call threshold%formfactor%activate ()
    call handle_test_onshell (n_tot, p_ofs_work)
    select case (FF)
    case (EXPANDED_HARD, EXPANDED_SOFT, EXPANDED_SOFT_SWITCHOFF, &
            EXPANDED_NOTSOHARD)
       call compute_born (n_tot, p_ofs_work, FF, sel_hel_beam, &
            amp_with_FF, tree_contrib=.false.)
       if (threshold%settings%interference) then
          amp_no_FF = amp_omega_full
       else
          call compute_born (n_tot, p_ofs_work, TREE, sel_hel_beam, &
               amp_no_FF, tree_contrib=.true.)
       end if
       amp2 = real (sum (abs2 (amp_no_FF))) + &
            2 * sum (real (amp_no_FF * conjg (amp_with_FF)))   !!! No |FFtilde|^2 here
    case (MATCHED, MATCHED_NOTSOHARD)
       call compute_born (n_tot, p_ofs_work, RESUMMED, sel_hel_beam, &
            amp_with_FF, tree_contrib=.false.)
       if (threshold%settings%factorized_computation .and. ( &
            threshold%settings%helicity_approximation%simple .or. &
            threshold%settings%helicity_approximation%extra .or. &
            threshold%settings%helicity_approximation%ultra)) then
          amp2 = real (sum (amp_with_FF))
       else
          amp2 = real (sum (abs2 (amp_with_FF)))                 !!! |FFtilde|^2
       end if
       select case (leg)
       case (-1)    !!! BORN
          simple = threshold%settings%helicity_approximation%simple
          threshold%settings%helicity_approximation%simple = .false.
          extra = threshold%settings%helicity_approximation%extra
          threshold%settings%helicity_approximation%extra = .false.
          ultra = threshold%settings%helicity_approximation%ultra
          threshold%settings%helicity_approximation%ultra = .false.
          if (FF == MATCHED) then
             use_ff = MATCHED_EXPANDED
          else
             use_ff = MATCHED_EXPANDED_NOTSOHARD
          end if
          call compute_born (n_tot, p_ofs_work, use_ff, sel_hel_beam, &
               amp_with_FF, tree_contrib=.false.)
          amp2 = amp2 + &                                      !!! 2 RealOf{resummed - expanded}{full}
               2 * sum (real (amp_with_FF * conjg (amp_omega_full)))
          threshold%settings%helicity_approximation%simple = simple
          threshold%settings%helicity_approximation%extra = extra
          threshold%settings%helicity_approximation%ultra = ultra
       case (0)     !!! SUBTRACTION
       case default
          call msg_bug ("threshold: compute_born_special_cases: this should not happen!")
       end select
    case default
       call compute_born (n_tot, p_ofs_work, FF, sel_hel_beam, amp_with_FF)
       if (threshold%settings%interference) then
          amp_no_FF = amp_omega_full
          if (threshold%settings%flip_relative_sign) &
              amp_no_FF = - amp_no_FF
          amp2 = real (sum (abs2 (amp_omega_full) + abs2 (amp_with_FF) + &
               2 * real (amp_no_FF * conjg (amp_with_FF))))
       else
          if (threshold%settings%factorized_computation .and. ( &
               threshold%settings%helicity_approximation%simple .or. &
               threshold%settings%helicity_approximation%extra .or. &
               threshold%settings%helicity_approximation%ultra)) then
             amp2 = real (sum (amp_with_FF))
          else
             amp2 = real (sum (abs2 (amp_with_FF)))
          end if
       end if
    end select
    if (test_ward)  amp2 = 0
    if (threshold%settings%only_interference_term) then
       if (FF == TREE) then
          call compute_born (n_tot, p_ofs_work, FF, sel_hel_beam,  &
               amp_with_FF, tree_contrib=.true.)
       else
          call compute_born (n_tot, p_ofs_work, FF, sel_hel_beam, &
               amp_with_FF, tree_contrib=.false.)
       end if
       amp2 = sum (2 * real (amp_omega_full * conjg (amp_with_FF)))
    end if
  end function compute_born_special_cases

end subroutine @ID@_get_amp_squared
