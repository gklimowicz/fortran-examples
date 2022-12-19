! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

submodule (ttv_formfactors) ttv_formfactors_s

  use numeric_utils
  use string_utils
  use physics_defs, only: CF, CA, TR
  use io_units
  use system_dependencies

  implicit none

  integer, parameter :: VECTOR = 1
  integer, parameter :: AXIAL = 2
  real(default), parameter :: NF = 5.0_default

  real(default), parameter :: z3 = 1.20205690315959428539973816151_default
  real(default), parameter :: A1 = 31./9.*CA - 20./9.*TR*NF
  real(default), parameter :: A2 = (4343./162. + 4.*pi**2 - pi**4/4. + &
       22./3.*z3)*CA**2 - (1798./81. + 56./3.*z3)*CA*TR*NF - &
       (55./3. - 16.*z3)*CF*TR*NF + (20./9.*TR*NF)**2
  complex(default), parameter :: ieps = imago*tiny_10


  interface char
    module procedure int_to_char, real_to_char, complex_to_char, logical_to_char
  end interface char


contains

  module subroutine onshell_projection_debug_write (onshell_projection)
    class(onshell_projection_t), intent(in) :: onshell_projection
    if (debug_on) call msg_debug (D_THRESHOLD, "onshell_projection%production", &
         onshell_projection%production)
    if (debug_on) call msg_debug (D_THRESHOLD, "onshell_projection%decay", &
         onshell_projection%decay)
    if (debug_on) call msg_debug (D_THRESHOLD, "onshell_projection%width", &
         onshell_projection%width)
    if (debug_on) call msg_debug (D_THRESHOLD, "onshell_projection%boost_decay", &
         onshell_projection%boost_decay)
  end subroutine onshell_projection_debug_write

  pure module subroutine onshell_projection_set_all (onshell_projection, flag)
    class(onshell_projection_t), intent(inout) :: onshell_projection
    logical, intent(in) :: flag
    onshell_projection%production = flag
    onshell_projection%decay = flag
  end subroutine onshell_projection_set_all

  pure module function onshell_projection_active (onshell_projection) result (active)
    logical :: active
    class(onshell_projection_t), intent(in) :: onshell_projection
    active = onshell_projection%production .or. &
         onshell_projection%decay
  end function onshell_projection_active

  ! TODO: (bcn 2016-03-21) break this up into a part regarding the
  ! FF grid and a part regarding the settings
  module subroutine settings_setup_flags (settings, ff_in, offshell_strategy_in, &
           top_helicity_selection)
    class(settings_t), intent(inout) :: settings
    integer, intent(in) :: ff_in, offshell_strategy_in, top_helicity_selection
    logical :: bit_top, bit_topbar
    !!! RESUMMED_SWITCHOFF = - 2
    !!! MATCHED = -1, &
    SWITCHOFF_RESUMMED                 = ff_in < 0
    TOPPIK_RESUMMED                    = ff_in <= 1
    settings%nlo = btest(offshell_strategy_in, 0)
    settings%factorized_computation = btest(offshell_strategy_in, 1)
    settings%interference = btest(offshell_strategy_in, 2)
    call settings%onshell_projection%set_all(btest(offshell_strategy_in, 3))
    settings%no_nlo_width_in_signal_propagators = btest(offshell_strategy_in, 4)
    settings%helicity_approximation%simple = btest(offshell_strategy_in, 5)
    if (.not. settings%onshell_projection%active ()) then
       settings%onshell_projection%production = btest(offshell_strategy_in, 6)
       settings%onshell_projection%decay = btest(offshell_strategy_in, 7)
    end if
    settings%onshell_projection%width = .not. btest(offshell_strategy_in, 8)
    settings%onshell_projection%boost_decay = btest(offshell_strategy_in, 9)
    settings%helicity_approximation%extra = btest(offshell_strategy_in, 10)
    settings%force_minus_one = btest(offshell_strategy_in, 11)
    settings%flip_relative_sign = btest(offshell_strategy_in, 12)
    if (top_helicity_selection > -1) then
       settings%helicity_approximation%ultra = .true.
       bit_top = btest (top_helicity_selection, 0)
       bit_topbar = btest (top_helicity_selection, 1)
       if (bit_top) then
          settings%sel_hel_top = 1
       else
          settings%sel_hel_top = -1
       end if
       if (bit_topbar) then
          settings%sel_hel_topbar = 1
       else
          settings%sel_hel_topbar = -1
       end if
    end if
    settings%only_interference_term = btest(offshell_strategy_in, 14)
    settings%Z_disabled = btest(offshell_strategy_in, 15)
    if (ff_in == MATCHED .or. ff_in == MATCHED_NOTSOHARD) then
       settings%onshell_projection%width = .true.
       settings%onshell_projection%production = .true.
       settings%onshell_projection%decay = .true.
       settings%factorized_computation = .true.
       settings%interference = .true.
       settings%onshell_projection%boost_decay = .true.
    end if
    if (debug_on) call msg_debug (D_THRESHOLD, "SWITCHOFF_RESUMMED", SWITCHOFF_RESUMMED)
    if (debug_on) call msg_debug (D_THRESHOLD, "TOPPIK_RESUMMED", TOPPIK_RESUMMED)
    if (debug_active (D_THRESHOLD)) &
         call settings%write ()
  end subroutine settings_setup_flags

  module subroutine settings_write (settings, unit)
    class(settings_t), intent(in) :: settings
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, '(A,L1)') "settings%helicity_approximation%simple = ", &
         settings%helicity_approximation%simple
    write (u, '(A,L1)') "settings%helicity_approximation%extra = ", &
         settings%helicity_approximation%extra
    write (u, '(A,L1)') "settings%helicity_approximation%ultra = ", &
         settings%helicity_approximation%ultra
    write (u, '(A,L1)') "settings%initialized_parameters = ", &
         settings%initialized_parameters
    write (u, '(A,L1)') "settings%initialized_ps = ", &
         settings%initialized_ps
    write (u, '(A,L1)') "settings%initialized_ff = ", &
         settings%initialized_ff
    write (u, '(A,L1)') "settings%mpole_dynamic = ", &
         settings%mpole_dynamic
    write (u, '(A,I5)') "settings%offshell_strategy = ", &
         settings%offshell_strategy
    write (u, '(A,L1)') "settings%factorized_computation = ", &
         settings%factorized_computation
    write (u, '(A,L1)') "settings%interference = ", settings%interference
    write (u, '(A,L1)') "settings%only_interference_term = ", &
         settings%only_interference_term
    write (u, '(A,L1)') "settings%Z_disabled = ", &
         settings%Z_disabled
    write (u, '(A,L1)') "settings%nlo = ", settings%nlo
    write (u, '(A,L1)') "settings%no_nlo_width_in_signal_propagators = ", &
         settings%no_nlo_width_in_signal_propagators
    write (u, '(A,L1)') "settings%force_minus_one = ", settings%force_minus_one
    write (u, '(A,L1)') "settings%flip_relative_sign = ", settings%flip_relative_sign
    call settings%onshell_projection%debug_write ()
  end subroutine settings_write

  pure module function settings_use_nlo_width (settings, ff) result (nlo)
    logical :: nlo
    class(settings_t), intent(in) :: settings
    integer, intent(in) :: ff
    nlo = settings%nlo
  end function settings_use_nlo_width

  pure module subroutine formfactor_activate (formfactor)
    class(formfactor_t), intent(inout) :: formfactor
    formfactor%active = .true.
  end subroutine formfactor_activate

  pure module subroutine formfactor_disable (formfactor)
    class(formfactor_t), intent(inout) :: formfactor
    formfactor%active = .false.
  end subroutine formfactor_disable

  module function formfactor_compute (formfactor, ps, vec_type, FF_mode) result (FF)
    complex(default) :: FF
    class(formfactor_t), intent(in) :: formfactor
    type(phase_space_point_t), intent(in) :: ps
    integer, intent(in) :: vec_type, FF_mode
    real(default) :: f
    if (threshold%settings%initialized_parameters .and. formfactor%active) then
       select case (FF_mode)
       case (MATCHED, MATCHED_NOTSOHARD, RESUMMED, RESUMMED_SWITCHOFF)
          FF = resummed_formfactor (ps, vec_type) - one
       case (MATCHED_EXPANDED)
          f = f_switch_off (v_matching (ps%sqrts, GAM_M1S))
          FF = - expanded_formfactor (f * AS_HARD, f * AS_HARD, ps, vec_type) &
               + resummed_formfactor (ps, vec_type)
       case (MATCHED_EXPANDED_NOTSOHARD)
          f = f_switch_off (v_matching (ps%sqrts, GAM_M1S))
          FF = - expanded_formfactor (f * alphas_notsohard (ps%sqrts), f * &
               alphas_notsohard (ps%sqrts), ps, vec_type) &
               + resummed_formfactor (ps, vec_type)
       case (EXPANDED_HARD)
          FF = expanded_formfactor (AS_HARD, AS_HARD, ps, vec_type) - one
       case (EXPANDED_NOTSOHARD)
          FF = expanded_formfactor (alphas_notsohard (ps%sqrts), &
               alphas_notsohard (ps%sqrts), ps, vec_type) - one
       case (EXPANDED_SOFT)
          FF = expanded_formfactor (AS_HARD, alphas_soft (ps%sqrts), ps, &
               vec_type) - one
       case (EXPANDED_SOFT_SWITCHOFF)
          f = f_switch_off (v_matching (ps%sqrts, GAM_M1S))
          FF = expanded_formfactor (f * AS_HARD, &
               f * alphas_soft (ps%sqrts), ps, vec_type) - one
       case (RESUMMED_ANALYTIC_LL)
          FF = formfactor_LL_analytic (alphas_soft (ps%sqrts), ps%sqrts, &
               ps%p, vec_type) - one
       case (TREE)
          FF = zero
       case default
          FF = zero
       end select
    else
       FF = zero
    end if
    if (debug2_active (D_THRESHOLD)) then
       call update_global_sqrts_dependent_variables (ps%sqrts)
       call msg_debug2 (D_THRESHOLD, "threshold%settings%initialized_parameters", &
            threshold%settings%initialized_parameters)
       call msg_debug2 (D_THRESHOLD, "formfactor%active", formfactor%active)
       call msg_debug2 (D_THRESHOLD, "FF_mode", FF_mode)
       call msg_debug2 (D_THRESHOLD, "FF", FF)
       call msg_debug2 (D_THRESHOLD, "v", sqrts_to_v (ps%sqrts, GAM))
       call msg_debug2 (D_THRESHOLD, "vec_type", vec_type)
       call ps%write ()
    end if
  end function formfactor_compute

  pure module subroutine width_init (width, aemi, sw, mw, mb, vtb, gam_inv)
    class(width_t), intent(inout) :: width
    real(default), intent(in) :: aemi, sw, mw, mb, vtb, gam_inv
    width%aem = one / aemi
    width%sw = sw
    width%mw = mw
    width%mb = mb
    width%vtb = vtb
    width%gam_inv = gam_inv
  end subroutine width_init

  pure module function width_compute (width, top_mass, sqrts, initial) result (gamma)
    real(default) :: gamma
    class(width_t), intent(in) :: width
    real(default), intent(in) :: top_mass, sqrts
    logical, intent(in), optional :: initial
    real(default) :: alphas
    logical :: ini
    ini = .false.;  if (present (initial))  ini = initial
    if (ini) then
       alphas = AS_HARD
    else
       alphas = alphas_notsohard (sqrts)
    end if
    if (threshold%settings%nlo) then
       gamma = top_width_sm_qcd_nlo_jk (width%aem, width%sw, width%vtb, &
            top_mass, width%mw, width%mb, alphas) + width%gam_inv
    else
       gamma = top_width_sm_lo (width%aem, width%sw, width%vtb, top_mass, &
            width%mw, width%mb) + width%gam_inv
    end if
  end function width_compute

  pure module subroutine phase_space_point_init_rel (ps_point, p2, k2, q2, m)
    class(phase_space_point_t), intent(inout) :: ps_point
    real(default), intent(in) :: p2
    real(default), intent(in) :: k2
    real(default), intent(in) :: q2
    real(default), intent(in), optional :: m
    ps_point%p2 = p2
    ps_point%k2 = k2
    ps_point%q2 = q2
    call rel_to_nonrel (p2, k2, q2, ps_point%sqrts, ps_point%p, ps_point%p0)
    ps_point%mpole = m1s_to_mpole (ps_point%sqrts)
    ps_point%en = sqrts_to_en (ps_point%sqrts)
    ps_point%inside_grid = sqrts_within_range (ps_point%sqrts)
    if ( present(m) ) ps_point%onshell = ps_point%is_onshell (m)
  end subroutine phase_space_point_init_rel

  pure module subroutine phase_space_point_init_nonrel (ps_point, sqrts, p, p0, m)
    class(phase_space_point_t), intent(inout) :: ps_point
    real(default), intent(in) :: sqrts
    real(default), intent(in) :: p
    real(default), intent(in) :: p0
    real(default), intent(in), optional :: m
    ps_point%sqrts = sqrts
    ps_point%p = p
    ps_point%p0 = p0
    call nonrel_to_rel (sqrts, p, p0, ps_point%p2, ps_point%k2, ps_point%q2)
    ps_point%mpole = m1s_to_mpole (sqrts)
    ps_point%en = sqrts_to_en (sqrts, ps_point%mpole)
    ps_point%inside_grid = sqrts_within_range (sqrts)
    if ( present(m) ) ps_point%onshell = ps_point%is_onshell (m)
  end subroutine phase_space_point_init_nonrel

  !!! convert squared 4-momenta into sqrts, p0 = E_top-sqrts/2 and abs. 3-momentum p
  pure subroutine rel_to_nonrel (p2, k2, q2, sqrts, p, p0)
    real(default), intent(in) :: p2
    real(default), intent(in) :: k2
    real(default), intent(in) :: q2
    real(default), intent(out) :: sqrts
    real(default), intent(out) :: p
    real(default), intent(out) :: p0
    sqrts = sqrt(q2)
    p0 = abs(p2 - k2) / (2. * sqrts)
    p = sqrt (0.5_default * (- p2 - k2 + sqrts**2/2. + 2.* p0**2))
  end subroutine rel_to_nonrel

  !!! convert sqrts, p0 = E_top-sqrts/2 and abs. 3-momentum p into squared 4-momenta
  pure subroutine nonrel_to_rel (sqrts, p, p0, p2, k2, q2)
    real(default), intent(in) :: sqrts
    real(default), intent(in) :: p
    real(default), intent(in) :: p0
    real(default), intent(out) :: p2
    real(default), intent(out) :: k2
    real(default), intent(out) :: q2
    p2 = (sqrts/2.+p0)**2 - p**2
    k2 = (sqrts/2.-p0)**2 - p**2
    q2 = sqrts**2
  end subroutine nonrel_to_rel

  pure function complex_m2 (m, w) result (m2c)
    real(default), intent(in) :: m
    real(default), intent(in) :: w
    complex(default) :: m2c
    m2c = m**2 - imago*m*w
  end function complex_m2

  pure module function phase_space_point_is_onshell (ps_point, m) result (flag)
    logical :: flag
    class(phase_space_point_t), intent(in) :: ps_point
    real(default), intent(in) :: m
    flag = nearly_equal (ps_point%p2 , m**2, rel_smallness=1E-5_default) .and. &
         nearly_equal (ps_point%k2 , m**2, rel_smallness=1E-5_default)
  end function phase_space_point_is_onshell

  module subroutine phase_space_point_write (psp, unit)
    class(phase_space_point_t), intent(in) :: psp
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, '(A)') char ("p2 = " // str (psp%p2))
    write (u, '(A)') char ("k2 = " // str (psp%k2))
    write (u, '(A)') char ("q2 = " // str (psp%q2))
    write (u, '(A)') char ("sqrts = " // str (psp%sqrts))
    write (u, '(A)') char ("p = " // str (psp%p))
    write (u, '(A)') char ("p0 = " // str (psp%p0))
    write (u, '(A)') char ("mpole = " // str (psp%mpole))
    write (u, '(A)') char ("en = " // str (psp%en))
    write (u, '(A)') char ("inside_grid = " // str (psp%inside_grid))
    write (u, '(A)') char ("onshell = " // str (psp%onshell))
  end subroutine phase_space_point_write

  function set_nrqcd_order (nrqcd_order_in) result (nrqcdorder)
    integer :: nrqcdorder
    real(default), intent(in) :: nrqcd_order_in
    nrqcdorder = 1
    if ( int(nrqcd_order_in) > nrqcdorder ) then
      call msg_warning ("reset to highest available NRQCD_ORDER = " // char(nrqcdorder))
    else
      nrqcdorder = int(nrqcd_order_in)
    end if
  end function set_nrqcd_order

  module subroutine init_parameters (mpole_out, gam_out, m1s_in, Vtb, gam_inv, &
         aemi, sw, az, mz, mw, mb, h_in, f_in, nrqcd_order_in, ff_in, &
         offshell_strategy_in, v1_in, v2_in, scan_sqrts_min, &
         scan_sqrts_max, scan_sqrts_stepsize, mpole_fixed, top_helicity_selection)
    real(default), intent(out) :: mpole_out
    real(default), intent(out) :: gam_out
    real(default), intent(in) :: m1s_in
    real(default), intent(in) :: Vtb
    real(default), intent(in) :: gam_inv
    real(default), intent(in) :: aemi
    real(default), intent(in) :: sw
    real(default), intent(in) :: az
    real(default), intent(in) :: mz
    real(default), intent(in) :: mw
    real(default), intent(in) :: mb
    real(default), intent(in) :: h_in
    real(default), intent(in) :: f_in
    real(default), intent(in) :: nrqcd_order_in
    real(default), intent(in) :: ff_in
    real(default), intent(in) :: offshell_strategy_in
    real(default), intent(in) :: v1_in
    real(default), intent(in) :: v2_in
    real(default), intent(in) :: scan_sqrts_min
    real(default), intent(in) :: scan_sqrts_max
    real(default), intent(in) :: scan_sqrts_stepsize
    logical, intent(in) :: mpole_fixed
    real(default), intent(in) :: top_helicity_selection
    if (debug_active (D_THRESHOLD))  call show_input()
    threshold%settings%initialized_parameters = .false.
    M1S = m1s_in
    threshold%settings%mpole_dynamic = .not. mpole_fixed
    threshold%settings%offshell_strategy = int (offshell_strategy_in)
    call threshold%settings%setup_flags (int(ff_in), &
         threshold%settings%offshell_strategy, &
         int (top_helicity_selection))
    NRQCD_ORDER = set_nrqcd_order (nrqcd_order_in)
    v1 = v1_in
    v2 = v2_in
    sqrts_min = scan_sqrts_min
    sqrts_max = scan_sqrts_max
    sqrts_it = scan_sqrts_stepsize
    !!! global hard parameters incl. hard alphas used in all form factors
    RESCALE_H = h_in
    MU_HARD   = M1S * RESCALE_H
    AS_MZ     = az
    MASS_Z    = mz
    AS_HARD   = running_as (MU_HARD, az, mz, 2, NF)
    call threshold%width%init (aemi, sw, mw, mb, vtb, gam_inv)
    GAM_M1S = threshold%width%compute (M1S, zero, initial=.true.)
    call compute_global_auxiliary_numbers ()
    !!! soft parameters incl. mtpole
    !!! (depend on sqrts: initialize with sqrts ~ 2*M1S)
    NUSTAR_FIXED = - one
    NUSTAR_DYNAMIC = NUSTAR_FIXED  < zero
    RESCALE_F = f_in
    call update_global_sqrts_dependent_variables (2. * M1S)
    mtpole_init = MTPOLE
    mpole_out = mtpole_init
    gam_out = GAM
    threshold%settings%initialized_parameters = .true.
  contains
      subroutine show_input()
        if (debug_on) call msg_debug (D_THRESHOLD, "init_parameters")
        if (debug_on) call msg_debug (D_THRESHOLD, "m1s_in", m1s_in)
        if (debug_on) call msg_debug (D_THRESHOLD, "Vtb", Vtb)
        if (debug_on) call msg_debug (D_THRESHOLD, "gam_inv", gam_inv)
        if (debug_on) call msg_debug (D_THRESHOLD, "aemi", aemi)
        if (debug_on) call msg_debug (D_THRESHOLD, "sw", sw)
        if (debug_on) call msg_debug (D_THRESHOLD, "az", az)
        if (debug_on) call msg_debug (D_THRESHOLD, "mz", mz)
        if (debug_on) call msg_debug (D_THRESHOLD, "mw", mw)
        if (debug_on) call msg_debug (D_THRESHOLD, "mb", mb)
        if (debug_on) call msg_debug (D_THRESHOLD, "h_in", h_in)
        if (debug_on) call msg_debug (D_THRESHOLD, "f_in", f_in)
        if (debug_on) call msg_debug (D_THRESHOLD, "nrqcd_order_in", nrqcd_order_in)
        if (debug_on) call msg_debug (D_THRESHOLD, "ff_in", ff_in)
        if (debug_on) call msg_debug (D_THRESHOLD, "offshell_strategy_in", offshell_strategy_in)
        if (debug_on) call msg_debug (D_THRESHOLD, "top_helicity_selection", top_helicity_selection)
        if (debug_on) call msg_debug (D_THRESHOLD, "v1_in", v1_in)
        if (debug_on) call msg_debug (D_THRESHOLD, "v2_in", v2_in)
        if (debug_on) call msg_debug (D_THRESHOLD, "scan_sqrts_min", scan_sqrts_min)
        if (debug_on) call msg_debug (D_THRESHOLD, "scan_sqrts_max", scan_sqrts_max)
        if (debug_on) call msg_debug (D_THRESHOLD, "scan_sqrts_stepsize", scan_sqrts_stepsize)
        if (debug_on) call msg_debug (D_THRESHOLD, "AS_HARD", AS_HARD)
      end subroutine show_input

  end subroutine init_parameters

  subroutine compute_global_auxiliary_numbers ()
    !!! auxiliary numbers needed later
    !!! current coefficients Ai(S,L,J), cf. arXiv:hep-ph/0609151, Eqs. (63)-(64)
    !!! 3S1 coefficients (s-wave, vector current)
    B0 = coeff_b0(NF) * (4.*pi)
    B1 = coeff_b1(NF) * (4.*pi)**2
    aa2(1) = (CF*(CA*CF*(9.*CA - 100.*CF) - &
              B0*(26.*CA**2 + 19.*CA*CF - 32.*CF**2)))/(26.*B0**2 *CA)
    aa3(1) = CF**2/( B0**2 *(6.*B0 - 13.*CA)*(B0 - 2.*CA)) * &
              (CA**2 *(9.*CA - 100.*CF) + B0*CA*(74.*CF - CA*16.) - &
              6.*B0**2 *(2.*CF - CA))
    aa4(1) = (24.*CF**2 * (11.*CA - 3.*B0)*(5.*CA + 8.*CF)) / &
              (13.*CA*(6.*B0 - 13.*CA)**2)
    aa5(1) =  (CF**2 * (CA*(15.-28) + B0*5.))/(6.*(B0-2.*CA)**2)
    aa8(1) = zero
    aa0(1) = -((8.*CF*(CA + CF)*(CA + 2.*CF))/(3.*B0**2))
    !!! 3P1 coefficients (p-wave, axial vector current)
    aa2(2) = -1./3. * (CF*(CA+2.*CF)/B0 - CF**2/(4.*B0) )
    aa3(2) =  zero
    aa4(2) =  zero
    aa5(2) =  1./3. * CF**2/(4.*(B0-2.*CA))
    aa8(2) = -1./3. * CF**2/(B0-CA)
    aa0(2) = -1./3. * 8.*CA*CF*(CA+4.*CF)/(3.*B0**2)
  end subroutine compute_global_auxiliary_numbers

  module subroutine init_threshold_grids (test)
    real(default), intent(in) :: test
    if (debug_active (D_THRESHOLD)) then
       call msg_debug (D_THRESHOLD, "init_threshold_grids")
       call msg_debug (D_THRESHOLD, "TOPPIK_RESUMMED", TOPPIK_RESUMMED)
    end if
    if (test > zero) then
      call msg_message ("TESTING ONLY: Skip threshold initialization and use tree-level SM.")
      return
    end if
    if (.not. threshold%settings%initialized_parameters) &
         call msg_fatal ("init_threshold_grid: parameters not initialized!")
    !!! !!! !!! MAC OS X and BSD don't load the global module with parameter values stored
    !!! if (parameters_ref == parameters_string ()) return
    call dealloc_grids ()
    if (TOPPIK_RESUMMED) call init_formfactor_grid ()
    parameters_ref = parameters_string ()
  end subroutine init_threshold_grids

  !!! LL/NLL resummation of nonrelativistic Coulomb potential
  pure function resummed_formfactor (ps, vec_type) result (c)
    type(phase_space_point_t), intent(in) :: ps
    integer, intent(in) :: vec_type
    complex(default) :: c
    c = one
    if (.not. threshold%settings%initialized_ff .or. .not. ps%inside_grid) return
    if (POINTS_SQ > 1) then
       call interpolate_linear (sq_grid, p_grid, ff_grid(:,:,1,vec_type), ps%sqrts, ps%p, c)
    else
       call interpolate_linear (p_grid, ff_grid(1,:,1,vec_type), ps%p, c)
    end if
  end function resummed_formfactor

  !!! leading nonrelativistic O(alphas^1) contribution (-> expansion of resummation)
  function expanded_formfactor (alphas_hard, alphas_soft, ps, vec_type) result (FF)
    complex(default) :: FF
    real(default), intent(in) :: alphas_hard, alphas_soft
    type(phase_space_point_t), intent(in) :: ps
    integer, intent(in) :: vec_type
    real(default) :: shift_from_hard_current
    complex(default) :: v, contrib_from_potential
    FF = one
    if (.not. threshold%settings%initialized_parameters) return
    call update_global_sqrts_dependent_variables (ps%sqrts)
    v = sqrts_to_v (ps%sqrts, GAM)
    if (NRQCD_ORDER == 1) then
       if (vec_type == AXIAL) then
          shift_from_hard_current = - CF / pi
       else
          shift_from_hard_current = - two * CF / pi
       end if
    else
       shift_from_hard_current = zero
    end if
    if (ps%onshell) then
       contrib_from_potential = CF * ps%mpole * Pi / (4 * ps%p)
    else
       if (vec_type == AXIAL) then
          contrib_from_potential = - CF * ps%mpole / (two * ps%p) * &
               (imago * ps%mpole * v / ps%p + &
               (ps%mpole**2 * v**2 + (ps%p)**2) / (4 *Pi * (ps%p)**2) * ( &
               (log (- ps%mpole * v - ps%p))**2 - &
               (log (- ps%mpole * v + ps%p))**2 + &
               (log (ps%mpole * v - ps%p))**2 - &
               (log (ps%mpole * v + ps%p))**2 ))
       else
          contrib_from_potential = imago * CF * ps%mpole * &
               log ((ps%p + ps%mpole * v) / &
               (-ps%p + ps%mpole * v) + ieps) / (two * ps%p)
       end if
    end if
    FF = one + alphas_soft * contrib_from_potential + &
         alphas_hard * shift_from_hard_current
  end function expanded_formfactor

  subroutine init_formfactor_grid ()
    type(string_t) :: ff_file
    if (debug_on) call msg_debug (D_THRESHOLD, "init_formfactor_grid")
    threshold%settings%initialized_ff = .false.
    ff_file = "SM_tt_threshold.grid"
    call msg_message ()
    call msg_message ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call msg_message (" Initialize e+e- => ttbar threshold resummation:")
    call msg_message (" Use analytic (LL) or TOPPIK (NLL) form factors for ttA/ttZ vector")
    call msg_message (" and axial vector couplings (S/P-wave) in the threshold region.")
    call msg_message (" Cf. threshold shapes from A. Hoang et al.: [arXiv:hep-ph/0107144],")
    call msg_message (" [arXiv:1309.6323].")
    if (NRQCD_ORDER > 0) then
      call msg_message (" Numerical NLL solutions calculated with TOPPIK [arXiv:hep-ph/9904468]")
      call msg_message (" by M. Jezabek, T. Teubner.")
    end if
    call msg_message ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call msg_message ()
    call read_formfactor_grid (ff_file)
    if (.not. threshold%settings%initialized_ff) then
      if (.not. threshold%settings%initialized_ps) call init_threshold_phase_space_grid ()
      call scan_formfactor_over_phase_space_grid ()
      call write_formfactor_grid (ff_file)
    end if
  end subroutine init_formfactor_grid

  subroutine read_formfactor_grid (ff_file)
    type(string_t), intent(in) :: ff_file
    complex(single), dimension(:,:,:,:), allocatable :: ff_grid_sp
    character(len(parameters_ref)) :: parameters
    integer :: u, st
    logical :: ex
    integer, dimension(4) :: ff_shape
    if (debug_on) call msg_debug (D_THRESHOLD, "read_formfactor_grid")
    inquire (file=char(ff_file), exist=ex)
    if (.not. ex) return
    u = free_unit ()
    call msg_message ("Opening grid file: " // char(ff_file))
    open (unit=u, status='old', file=char(ff_file), form='unformatted', iostat=st)
    if (st /= 0) call msg_fatal ("iostat = " // char(st))
    read (u) parameters
    read (u) ff_shape
    if (ff_shape(4) /= 2)  call msg_fatal ("read_formfactor_grid: i = " // char(ff_shape(4)))
    if (parameters /= parameters_string ()) then
       call msg_message ("Threshold setup has changed: recalculate threshold grid.")
       close (unit=u, status='delete')
       return
    end if
    call msg_message ("Threshold setup unchanged: reusing existing threshold grid.")
    POINTS_SQ = ff_shape(1)
    POINTS_P = ff_shape(2)
    if (debug_active (D_THRESHOLD)) then
       call msg_debug (D_THRESHOLD, "ff_shape(1) (POINTS_SQ)", ff_shape(1))
       call msg_debug (D_THRESHOLD, "ff_shape(2)", ff_shape(2))
       call msg_debug (D_THRESHOLD, "ff_shape(3) (POINTS_P0)", ff_shape(3))
       call msg_debug (D_THRESHOLD, "ff_shape(4) (==2)", ff_shape(4))
    end if
    allocate (sq_grid(POINTS_SQ))
    read (u) sq_grid
    allocate (p_grid(POINTS_P))
    read (u) p_grid
    POINTS_P0 = ff_shape(3)
    allocate (ff_grid_sp(POINTS_SQ,POINTS_P,POINTS_P0,2))
    read (u) ff_grid_sp
    allocate (ff_grid(POINTS_SQ,POINTS_P,POINTS_P0,2))
    ff_grid = cmplx (ff_grid_sp, kind=default)
    close (u, iostat=st)
    if (st > 0)  call msg_fatal ("close " // char(ff_file) // ": iostat = " // char(st))
    threshold%settings%initialized_ps = .true.
    threshold%settings%initialized_ff = .true.
  end subroutine read_formfactor_grid

  subroutine write_formfactor_grid (ff_file)
    type(string_t), intent(in) :: ff_file
    integer :: u, st
    if (.not. threshold%settings%initialized_ff) then
      call msg_warning ("write_formfactor_grid: no grids initialized!")
      return
    end if
    u = free_unit ()
    open (unit=u, status='replace', file=char(ff_file), form='unformatted', iostat=st)
    if (st /= 0)  call msg_fatal ("open " // char(ff_file) // ": iostat = " // char(st))
    write (u) parameters_string ()
    write (u) shape(ff_grid)
    write (u) sq_grid
    write (u) p_grid
    write (u) cmplx(ff_grid, kind=single)
    close (u, iostat=st)
    if (st > 0)  call msg_fatal ("close " // char(ff_file) // ": iostat = " // char(st))
  end subroutine write_formfactor_grid

  pure function parameters_string () result (str)
    character(len(parameters_ref)) :: str
    str = char(M1S) // " " // char(GAM_M1S) &
         // " " // char(NRQCD_ORDER) &
         // " " // char(RESCALE_H) &
         // " " // char(RESCALE_F) &
         //  " " // char(sqrts_min) &
         // " " // char(sqrts_max) &
         // " " // char(sqrts_it)
  end function parameters_string

  subroutine update_global_sqrts_dependent_variables (sqrts)
    real(default), intent(in) :: sqrts
    real(default) :: nu_soft, f
    logical :: only_once_for_fixed_nu, already_done
    real(default), save :: last_sqrts = - one
    if (debug_on) call msg_debug (D_THRESHOLD, "update_global_sqrts_dependent_variables")
    if (debug_on) call msg_debug (D_THRESHOLD, "sqrts", sqrts)
    if (debug_on) call msg_debug (D_THRESHOLD, "last_sqrts", last_sqrts)
    already_done = threshold%settings%initialized_parameters .and. &
         nearly_equal (sqrts, last_sqrts, rel_smallness=1E-6_default)
    if (debug_on) call msg_debug (D_THRESHOLD, "already_done", already_done)
    only_once_for_fixed_nu = .not. NUSTAR_DYNAMIC .and. MTPOLE > zero
    if (debug_on) call msg_debug (D_THRESHOLD, "only_once_for_fixed_nu", only_once_for_fixed_nu)
    if (only_once_for_fixed_nu .or. already_done) return
    last_sqrts = sqrts
    nu_soft = RESCALE_F * nustar (sqrts)
    MU_SOFT = M1S * RESCALE_H * nu_soft
    MU_USOFT = M1S * RESCALE_H * nu_soft**2
    AS_SOFT = running_as (MU_SOFT, AS_HARD, MU_HARD, NRQCD_ORDER, NF)
    AS_LL_SOFT = running_as (MU_SOFT, AS_HARD, MU_HARD, 0, NF)
    AS_USOFT = running_as (MU_USOFT, AS_HARD, MU_HARD, 0, NF) !!! LL here
    if (SWITCHOFF_RESUMMED) then
       f = f_switch_off (v_matching (sqrts, GAM_M1S))
       AS_SOFT = AS_SOFT * f
       AS_LL_SOFT = AS_LL_SOFT * f
       AS_USOFT = AS_USOFT * f
    end if
    MTPOLE = m1s_to_mpole (sqrts)
    GAM = threshold%width%compute (MTPOLE, sqrts)
    if (debug_on) call msg_debug (D_THRESHOLD, "GAM", GAM)
    if (debug_on) call msg_debug (D_THRESHOLD, "nu_soft", nu_soft)
    if (debug_on) call msg_debug (D_THRESHOLD, "MTPOLE", MTPOLE)
    if (debug_on) call msg_debug (D_THRESHOLD, "AS_SOFT", AS_SOFT)
    if (debug_on) call msg_debug (D_THRESHOLD, "AS_LL_SOFT", AS_LL_SOFT)
    if (debug_on) call msg_debug (D_THRESHOLD, "AS_USOFT", AS_USOFT)
  end subroutine update_global_sqrts_dependent_variables

  !!! Coulomb potential coefficients needed by TOPPIK
  pure function xc (a_soft, i_xc) result (xci)
    real(default), intent(in) :: a_soft
    integer, intent(in) :: i_xc
    real(default) :: xci
    xci = zero
    select case (i_xc)
      case (0)
        xci = one
        if ( NRQCD_ORDER>0 ) xci = xci + a_soft/(4.*pi) * A1
        if ( NRQCD_ORDER>1 ) xci = xci + (a_soft/(4.*pi))**2 * A2
      case (1)
        if ( NRQCD_ORDER>0 ) xci = xci + a_soft/(4.*pi) * B0
        if ( NRQCD_ORDER>1 ) xci = xci + (a_soft/(4.*pi))**2 * (B1 + 2*B0*A1)
      case (2)
        if ( NRQCD_ORDER>1 ) xci = xci + (a_soft/(4.*pi))**2 * B0**2
      case default
        return
    end select
  end function xc

  function current_coeff (a_hard, a_soft, a_usoft, i) result (coeff)
    real(default), intent(in) :: a_hard, a_soft, a_usoft
    integer, intent(in) :: i
    real(default) :: coeff
    real(default) :: matching_c, c1
    real(default) :: z, w
    if (debug_on) call msg_debug (D_THRESHOLD, "current_coeff")
    coeff = one
    if (NRQCD_ORDER == 0) return
    z = a_soft / a_hard
    w = a_usoft / a_soft
    !!! hard s/p-wave 1-loop matching coefficients, cf. arXiv:hep-ph/0604072
    select case (i)
      case (1)
        matching_c = one - 2.*(CF/pi) * a_hard
      case (2)
        matching_c = one -    (CF/pi) * a_hard
     case default
        call msg_fatal ("current_coeff: unknown coeff i = " // char(i))
    end select
    !!! current coefficient c1, cf. arXiv:hep-ph/0609151, Eq. (62)
    c1 = exp( a_hard * pi * ( aa2(i)*(1.-z) + aa3(i)*log(z) + &
         aa4(i)*(1.-z**(1.-13.*CA/(6.*B0))) + aa5(i)*(1.-z**(1.-2.*CA/B0)) + &
         aa8(i)*(1.-z**(1.-CA/B0)) + aa0(i)*(z-1.-log(w)/w) ))
    coeff = matching_c * c1
  end function current_coeff

  pure module function v_matching (sqrts, gamma) result (v)
    real(default) :: v
    real(default), intent(in) :: sqrts, gamma
    v = abs (sqrts_to_v_1S (sqrts, gamma))
  end function v_matching

  pure module function f_switch_off (v) result (fval)
    real(default), intent(in) :: v
    real(default) :: fval
    real(default) :: vm, f1, f2, x
    f1 = one
    f2 = zero + tiny_10
    vm = (v1+v2) / 2.
    if ( v < v1 ) then
      fval = f1
    else if (v < v2) then
       x = (v - v1) / (v2 - v1)
       fval = 1 - x**2 * (3 - 2 * x)
    else
      fval = f2
    end if
  end function f_switch_off

  function formfactor_LL_analytic (a_soft, sqrts, p, vec_type) result (c)
    real(default), intent(in) :: a_soft
    real(default), intent(in) :: sqrts
    real(default), intent(in) :: p
    integer, intent(in) :: vec_type
    complex(default) :: c
    real(default) :: en
    c = one
    if (.not. threshold%settings%initialized_parameters) return
    call update_global_sqrts_dependent_variables (sqrts)
    en = sqrts_to_en (sqrts, MTPOLE)
    select case (vec_type)
      case (1)
        c = G0p (CF*a_soft, en, p, MTPOLE, GAM) / G0p_tree (en, p, MTPOLE, GAM)
      case (2)
        c = G0p_ax (CF*a_soft, en, p, MTPOLE, GAM) / G0p_tree (en, p, MTPOLE, GAM)
      case default
        call msg_fatal ("unknown ttZ/ttA vertex component, vec_type = " // char(vec_type))
    end select
  end function formfactor_LL_analytic

  !!! Max's LL nonrelativistic threshold Green's function
  function G0p (a, en, p, m, w) result (c)
    real(default), intent(in) :: a
    real(default), intent(in) :: en
    real(default), intent(in) :: p
    real(default), intent(in) :: m
    real(default), intent(in) :: w
    complex(default) :: c
    complex(default) :: k, ipk, la, z1, z2
    complex(default) :: one, two, cc, dd
    k   = sqrt( -m*en -imago*m*w )
    ipk = imago * p / k
    la  = a * m / 2. / k
    one = cmplx (1., kind=default)
    two = cmplx (2., kind=default)
    cc  = 2. - la
    dd  = ( 1. + ipk ) / 2.
    z1  = nr_hypgeo (two, one, cc, dd)
    dd  = ( 1. - ipk ) / 2.
    z2  = nr_hypgeo (two, one, cc, dd)
    c   = - imago * m / (4.*p*k) / (1.-la) * ( z1 - z2 )
  end function G0p

  !!! tree level version: a_soft -> 0
  pure function G0p_tree (en, p, m, w) result (c)
    real(default), intent(in) :: en
    real(default), intent(in) :: p
    real(default), intent(in) :: m
    real(default), intent(in) :: w
    complex(default) :: c
    c = m / (p**2 - m*(en+imago*w))
  end function G0p_tree

  !!! Peter Poier's LL nonrelativistic axial threshold Green's function
  function G0p_ax (a, en, p, m, w) result (c)
    real(default), intent(in) :: a
    real(default), intent(in) :: en
    real(default), intent(in) :: p
    real(default), intent(in) :: m
    real(default), intent(in) :: w
    complex(default) :: c
    complex(default) :: k, ipk, la, z1, z2, z3, z4
    complex(default) :: zero, two, three, cc, ddp, ddm
    k   = sqrt( -m*en -imago*m*w )
    ipk = imago * p / k
    la  = a * m / 2. / k
    zero = cmplx (0., kind=default)
    two = cmplx (2., kind=default)
    three = cmplx (3., kind=default)
    cc  = 1. - la
    ddp = ( 1. + ipk ) / 2.
    ddm = ( 1. - ipk ) / 2.
    z1  = nr_hypgeo (zero, two, cc, ddp)
    z2  = nr_hypgeo (zero, two, cc, ddm)
    cc  = 2. - la
    z3  = nr_hypgeo (zero, three, cc, ddm)
    z4  = nr_hypgeo (zero, three, cc, ddp)
    c   = m / 2. / p**3 * ( 2.*p + imago*k*(1.-la)*(z1-z2) + imago*k*(z3-z4) )
  end function G0p_ax

  pure function nustar (sqrts) result (nu)
    real(default), intent(in) :: sqrts
    real(default) :: nu
    real(default), parameter :: nustar_offset = 0.05_default
    complex(default) :: arg
    if (NUSTAR_DYNAMIC) then
      !!! from [arXiv:1309.6323], Eq. (3.2) (other definitions possible)
      arg = ( sqrts - 2.*M1S + imago*GAM_M1S ) / M1S
      nu  = nustar_offset + abs(sqrt(arg))
    else
      nu  = NUSTAR_FIXED
    end if
  end function nustar

  pure function alphas_soft (sqrts) result (a_soft)
    real(default) :: a_soft
    real(default), intent(in) :: sqrts
    real(default) :: mu_soft, nusoft
    nusoft = RESCALE_F * nustar (sqrts)
    mu_soft = RESCALE_H * M1S * nusoft
    a_soft = running_as (mu_soft, AS_HARD, MU_HARD, NRQCD_ORDER, NF)
  end function alphas_soft

  pure module function alphas_notsohard (sqrts) result (a_soft)
    real(default) :: a_soft
    real(default), intent(in) :: sqrts
    real(default) :: mu_notsohard
    ! complex(default) :: v
    ! v = sqrts_to_v_1S (sqrts, GAM_M1S)
    ! mu_notsohard = RESCALE_H * M1S * sqrt(abs(v))
    mu_notsohard = RESCALE_H * M1S * sqrt(nustar (sqrts))
    a_soft = running_as (mu_notsohard, AS_MZ, MASS_Z, 2, NF)
  end function alphas_notsohard

  pure module function m1s_to_mpole (sqrts) result (mpole)
    real(default), intent(in) :: sqrts
    real(default) :: mpole
    mpole = mtpole_init
    if (threshold%settings%mpole_dynamic) then
       mpole = M1S * ( 1. + deltaM(sqrts) )
    else
       mpole = M1S
    end if
  end function m1s_to_mpole

  !pure
  !function mpole_to_M1S (mpole, sqrts, nl) result (m)
    !real(default), intent(in) :: mpole
    !real(default), intent(in) :: sqrts
    !integer, intent(in) :: nl
    !real(default) :: m
    !m = mpole * ( 1. - deltaM(sqrts, nl) )
  !end function mpole_to_M1S

  pure function deltaM (sqrts) result (del)
    real(default), intent(in) :: sqrts
    real(default) :: del
    real(default) :: ac
    ac  = CF * alphas_soft (sqrts)
    del = ac**2 / 8.
    if (NRQCD_ORDER > 0) then
      del = del + ac**3 / (8. * pi * CF) * &
           (B0 * (log (RESCALE_H * RESCALE_F * nustar (sqrts) / ac) + one) + A1 / 2.)
    end if
  end function deltaM

  pure function sqrts_within_range (sqrts) result (flag)
    real(default), intent(in) :: sqrts
    logical :: flag
    flag = ( sqrts >= sqrts_min - tiny_07 .and. sqrts <= sqrts_max + tiny_07 )
  end function

  ! The mapping is such that even for min=max, we get three points:
  ! min - it , min, min + it
  pure function sqrts_iter (i_sq) result (sqrts)
    integer, intent(in) :: i_sq
    real(default) :: sqrts
    if (POINTS_SQ > 1) then
       sqrts = sqrts_min - sqrts_it + &
               (sqrts_max - sqrts_min + two * sqrts_it) * &
               real(i_sq - 1) / real(POINTS_SQ - 1)
    else
       sqrts = sqrts_min
    end if
  end function sqrts_iter

  function scan_formfactor_over_p_LL_analytic (a_soft, sqrts, vec_type) result (ff_analytic)
    real(default), intent(in) :: a_soft
    real(default), intent(in) :: sqrts
    integer, intent(in) :: vec_type
    complex(default), dimension(POINTS_P) :: ff_analytic
    integer :: i_p
    ff_analytic = [(formfactor_LL_analytic (a_soft, sqrts, p_grid(i_p), vec_type), i_p=1, POINTS_P)]
  end function scan_formfactor_over_p_LL_analytic

  !!! tttoppik wrapper
  subroutine scan_formfactor_over_p_TOPPIK (a_soft, sqrts, vec_type, p_grid_out, mpole_in, ff_toppik)
    real(default), intent(in) :: a_soft
    real(default), intent(in) :: sqrts
    integer, intent(in) :: vec_type
    real(default), dimension(POINTS_P), intent(out), optional :: p_grid_out
    real(default), intent(in), optional :: mpole_in
    complex(default), dimension(POINTS_P), optional :: ff_toppik
    integer :: i_p
    real(default) :: mpole, alphas_hard, f
    real(default), dimension(POINTS_P) :: p_toppik
    type(nr_spline_t) :: toppik_spline
    real*8 :: xenergy, xtm, xtg, xalphas, xscale, xc0, xc1, xc2, xim, xdi, &
        xcutn, xcutv, xkincm, xkinca, xkincv, xcdeltc, &
        xcdeltl, xcfullc, xcfulll, xcrm2
    integer, parameter :: nmax=900
    real*8 :: xdsdp(nmax), xpp(nmax), xww(nmax)
    complex*16 :: zff(nmax)
    integer :: np, jknflg, jgcflg, jvflg
    if (debug_on) call msg_debug (D_THRESHOLD, "scan_formfactor_over_p_TOPPIK")
    if (POINTS_P > nmax-40) call msg_fatal ("TOPPIK: POINTS_P must be <=" // char(nmax-40))
    if (debug_on) call msg_debug (D_THRESHOLD, "POINTS_P", POINTS_P)
    if (present (ff_toppik))  ff_toppik = zero
    mpole = MTPOLE;  if (present (mpole_in)) mpole = mpole_in
    xenergy = sqrts_to_en (sqrts, MTPOLE)
    xtm     = mpole
    xtg     = GAM
    xalphas = a_soft
    xscale  = MU_SOFT
    xcutn   = 175.E6
    xcutv   = 175.E6
    xc0     = xc (a_soft, 0)
    xc1     = xc (a_soft, 1)
    xc2     = xc (a_soft, 2)
    xcdeltc = 0.
    xcdeltl = 0.
    xcfullc = 0.
    xcfulll = 0.
    xcrm2   = 0.
    xkincm  = 0.
    xkinca  = 0.
    jknflg  = 0
    jgcflg  = 0
    xkincv  = 0.
    jvflg   = 0
    select case (vec_type)
      case (VECTOR)
         if (debug_on) call msg_debug (D_THRESHOLD, "calling tttoppik")
         call tttoppik &
                (xenergy,xtm,xtg,xalphas,xscale,xcutn,xcutv,xc0,xc1,xc2, &
                 xcdeltc,xcdeltl,xcfullc,xcfulll,xcrm2,xkincm,xkinca,jknflg, &
                 jgcflg, xkincv,jvflg,xim,xdi,np,xpp,xww,xdsdp,zff)
      case (AXIAL)
         if (debug_on) call msg_debug (D_THRESHOLD, "calling tttoppikaxial")
         call tttoppikaxial &
                (xenergy,xtm,xtg,xalphas,xscale,xcutn,xcutv,xc0,xc1,xc2, &
                 xcdeltc,xcdeltl,xcfullc,xcfulll,xcrm2,xkincm,xkinca,jknflg, &
                 jgcflg, xkincv,jvflg,xim,xdi,np,xpp,xww,xdsdp,zff)
         !!! 1st ~10 TOPPIK p-wave entries are ff_unstable: discard them
         zff(1:10) = [(zff(11), i_p=1, 10)]
      case default
         call msg_fatal ("unknown ttZ/ttA vertex component, vec_type = " // char(vec_type))
    end select
    if (present (p_grid_out)) p_grid_out = xpp(1:POINTS_P)
    if (.not. present (ff_toppik)) return
    !!! keep track of TOPPIK instabilities and try to repair later
    if (np < 0) then
      ff_toppik(1) = 2.d30
      if (debug_active (D_THRESHOLD)) then
         call msg_warning ("caught TOPPIK instability at sqrts = " // char(sqrts))
      end if
      return
    end if
    p_toppik = xpp(1:POINTS_P)
    ff_toppik = zff(1:POINTS_P)
    !!! TOPPIK output p-grid scales with en above ~ 4 GeV:
    !!! interpolate for global sqrts/p grid
    if (.not. nearly_equal (p_toppik(42), p_grid(42), rel_smallness=1E-6_default)) then
      call toppik_spline%init (p_toppik, ff_toppik)
      ff_toppik(2:POINTS_P) = [(toppik_spline%interpolate (p_grid(i_p)), i_p=2, POINTS_P)]
      call toppik_spline%dealloc ()
    end if
    !!! TOPPIK output includes tree level ~ 1, a_soft @ LL in current coefficient!
    if (SWITCHOFF_RESUMMED) then
       f = f_switch_off (v_matching (sqrts, GAM_M1S))
       alphas_hard = AS_HARD * f
    else
       alphas_hard = AS_HARD
    end if
    ff_toppik = ff_toppik * current_coeff (alphas_hard, AS_LL_SOFT, AS_USOFT, vec_type)
    if (debug_on) call msg_debug (D_THRESHOLD, &
         "current_coeff (alphas_hard, AS_LL_SOFT, AS_USOFT, vec_type)", &
         current_coeff (alphas_hard, AS_LL_SOFT, AS_USOFT, vec_type))
  end subroutine scan_formfactor_over_p_TOPPIK

  function scan_formfactor_over_p (sqrts, vec_type) result (ff)
    real(default), intent(in) :: sqrts
    integer, intent(in) :: vec_type
    complex(default), dimension(POINTS_P) :: ff
    if (debug_on) call msg_debug (D_THRESHOLD, "scan_formfactor_over_p")
    select case (NRQCD_ORDER)
      case (0)
       ! ff = scan_formfactor_over_p_LL_analytic (AS_SOFT, sqrts, vec_type)
        call scan_formfactor_over_p_TOPPIK (AS_SOFT, sqrts, vec_type, ff_toppik=ff)
      case (1)
        call scan_formfactor_over_p_TOPPIK (AS_SOFT, sqrts, vec_type, ff_toppik=ff)
      case default
        call msg_fatal ("NRQCD_ORDER = " // char(NRQCD_ORDER))
    end select
  end function scan_formfactor_over_p

  subroutine scan_formfactor_over_phase_space_grid ()
    integer :: i_sq, vec_type, unstable_loop
    logical, dimension(:,:), allocatable :: ff_unstable
    real(default) :: t1, t2, t3, t_toppik, t_p0_dep
    if (debug_on) call msg_debug (D_THRESHOLD, "scan_formfactor_over_phase_space_grid")
    allocate (ff_grid(POINTS_SQ,POINTS_P,POINTS_P0,2))
    allocate (ff_unstable(POINTS_SQ,2))
    t_toppik = zero
    t_p0_dep = zero
    write (msg_buffer, "(3(A,F7.3,1X),A)") "Scanning from ", &
         sqrts_min - sqrts_it, "GeV to ", &
         sqrts_max + sqrts_it, "GeV in steps of ", sqrts_it, "GeV"
    call msg_message ()
    ENERGY_SCAN: do i_sq = 1, POINTS_SQ
      if (signal_is_pending ())  return
      call update_global_sqrts_dependent_variables (sq_grid(i_sq))
      !!! vector and axial vector
      do vec_type = VECTOR, AXIAL
        call cpu_time (t1)
        unstable_loop = 0
        UNTIL_STABLE: do
           ff_grid(i_sq,:,1,vec_type) = scan_formfactor_over_p (sq_grid(i_sq), vec_type)
           ff_unstable(i_sq,vec_type) = abs(ff_grid(i_sq,1,1,vec_type)) > 1.d30
           unstable_loop = unstable_loop + 1
           if (ff_unstable(i_sq,vec_type) .and. unstable_loop < 10) then
              cycle
           else
              exit
           end if
        end do UNTIL_STABLE
        call cpu_time (t2)
        !!!  include p0 dependence by an integration over the p0-independent FF
        call cpu_time (t3)
        t_toppik = t_toppik + t2 - t1
        t_p0_dep = t_p0_dep + t3 - t2
      end do
      call msg_show_progress (i_sq, POINTS_SQ)
    end do ENERGY_SCAN
    if (debug_active (D_THRESHOLD)) then
       print *, "time for TOPPIK call:   ", t2 - t1, " seconds."
       print *, "time for p0 dependence: ", t3 - t2, " seconds."
    end if
    if (any (ff_unstable))  call handle_TOPPIK_instabilities (ff_grid, ff_unstable)
    if (allocated(Vmatrix))  deallocate(Vmatrix)
    if (allocated(q_grid))  deallocate(q_grid)
    threshold%settings%initialized_ff = .true.
  end subroutine scan_formfactor_over_phase_space_grid

  subroutine init_threshold_phase_space_grid ()
    integer :: i_sq
    if (debug_on) call msg_debug (D_THRESHOLD, "init_threshold_phase_space_grid")
    if (sqrts_it > tiny_07) then
       POINTS_SQ = int ((sqrts_max - sqrts_min) / sqrts_it + tiny_07) + 3
    else
       POINTS_SQ = 1
    end if
    if (debug_on) call msg_debug (D_THRESHOLD, "Number of sqrts grid points: POINTS_SQ", POINTS_SQ)
    if (debug_on) call msg_debug (D_THRESHOLD, "sqrts_max", sqrts_max)
    if (debug_on) call msg_debug (D_THRESHOLD, "sqrts_min", sqrts_min)
    if (debug_on) call msg_debug (D_THRESHOLD, "sqrts_it", sqrts_it)
    allocate (sq_grid(POINTS_SQ))
    sq_grid = [(sqrts_iter (i_sq), i_sq=1, POINTS_SQ)]
    POINTS_P = 600
    allocate (p_grid(POINTS_P))
    p_grid = p_grid_from_TOPPIK ()
    POINTS_P0 = 1
    threshold%settings%initialized_ps = .true.
  end subroutine init_threshold_phase_space_grid

  subroutine init_p0_grid (p_in, n)
    real(default), dimension(:), allocatable, intent(in) :: p_in
    integer, intent(in) :: n
    if (debug_on) call msg_debug (D_THRESHOLD, "init_p0_grid")
    if (debug_on) call msg_debug (D_THRESHOLD, "n", n)
    if (debug_on) call msg_debug (D_THRESHOLD, "size(p_in)", size(p_in))
    if (.not. allocated (p_in))  call msg_fatal ("init_p0_grid: p_in not allocated!")
    if (allocated (p0_grid))  deallocate (p0_grid)
    allocate (p0_grid(n))
    p0_grid(1) = zero
    p0_grid(2:n) = p_in(1:n-1)
  end subroutine init_p0_grid

  !!! Andre's procedure to refine an existing grid
  pure subroutine finer_grid (gr, fgr, n_in)
    real(default), dimension(:), intent(in) :: gr
    real(default), dimension(:), allocatable, intent(inout) :: fgr
    integer, intent(in), optional :: n_in
    integer :: n, i, j
    real(default), dimension(:), allocatable :: igr
    n = 4
    if ( present(n_in) ) n = n_in
    allocate( igr(n) )
    if ( allocated(fgr) ) deallocate( fgr )
    allocate( fgr(n*(size(gr)-1)+1) )
    do i=1, size(gr)-1
      do j=0, n-1
        igr(j+1) = gr(i) + real(j)*(gr(i+1)-gr(i))/real(n)
      end do
      fgr((i-1)*n+1:i*n) = igr
    end do
    fgr(size(fgr)) = gr(size(gr))
    deallocate( igr )
  end subroutine finer_grid

  subroutine dealloc_grids ()
    if ( allocated(sq_grid) ) deallocate( sq_grid )
    if ( allocated( p_grid) ) deallocate(  p_grid )
    if ( allocated(p0_grid) ) deallocate( p0_grid )
    if ( allocated(ff_grid) ) deallocate( ff_grid )
    threshold%settings%initialized_ps = .false.
    threshold%settings%initialized_ff = .false.
  end subroutine dealloc_grids

  subroutine trim_p_grid (n_p_new)
    integer, intent(in) :: n_p_new
    real(default), dimension(n_p_new) :: p_save
    complex(default), dimension(POINTS_SQ,n_p_new,POINTS_P0,2) :: ff_save
    if (n_p_new > POINTS_P) then
      call msg_fatal ("trim_p_grid: new size larger than old size.")
      return
    end if
    p_save = p_grid(1:n_p_new)
    ff_save = ff_grid(:,1:n_p_new,:,:)
    deallocate( p_grid, ff_grid )
    allocate( p_grid(n_p_new), ff_grid(POINTS_SQ,n_p_new,POINTS_P0,2) )
    p_grid = p_save
    ff_grid = ff_save
  end subroutine trim_p_grid

  !!! try to repair TOPPIK instabilities by interpolation of adjacent sq_grid points
  subroutine handle_TOPPIK_instabilities (ff, nan)
    complex(default), dimension(:,:,:,:), intent(inout) :: ff
    logical, dimension(:,:), intent(in) :: nan
    integer :: i, i_sq, n_nan
    logical :: interrupt
    n_nan = sum (merge ([(1, i=1, 2*POINTS_SQ)], &
         [(0, i=1, 2*POINTS_SQ)], reshape (nan, [2*POINTS_SQ])) )
    interrupt = n_nan > 3
    do i = 1, 2
      if (interrupt ) exit
      if (.not. any (nan(:,i))) cycle
      do i_sq = 2, POINTS_SQ - 1
        if (.not. nan(i_sq,i)) cycle
        if (nan(i_sq+1,i) .or. nan(i_sq-1,i)) then
          interrupt = .true.
          exit
        end if
        ff(i_sq,:,:,i) = (ff(i_sq-1,:,:,i) + ff(i_sq+1,:,:,i)) / two
      end do
    end do
    if (.not. interrupt) return
    call msg_fatal ("Too many TOPPIK instabilities! Check your parameter setup " &
                     // "or slightly vary the scales sh and/or sf.")
  end subroutine handle_TOPPIK_instabilities

  pure function sqrts_to_v (sqrts, gamma) result (v)
    complex(default) :: v
    real(default), intent(in) :: sqrts, gamma
    real(default) :: m
    m = m1s_to_mpole (sqrts)
    v = sqrt ((sqrts - two * m + imago * gamma) / m)
  end function sqrts_to_v

  pure function sqrts_to_v_1S (sqrts, gamma) result (v)
    complex(default) :: v
    real(default), intent(in) :: sqrts, gamma
    v = sqrt ((sqrts - two * M1S + imago * gamma) / M1S)
  end function sqrts_to_v_1S

  pure function v_to_sqrts (v) result (sqrts)
    real(default), intent(in) :: v
    real(default) :: sqrts
    real(default) :: m
    m = mtpole_init
    sqrts = 2.*m + m*v**2
  end function v_to_sqrts

  !!! -q^2 times the Coulomb potential V at LO resp. NLO
  function minus_q2_V (a, q, p, p0r, vec_type) result (v)
    real(default), intent(in) :: a
    real(default), intent(in) :: q
    real(default), intent(in) :: p
    real(default), intent(in) :: p0r
    integer, intent(in) :: vec_type
    complex(default) :: p0, log_mppp, log_mmpm, log_mu_s, v
    p0 = abs(p0r) + ieps
    log_mppp = log( (p-p0+q) * (p+p0+q) )
    log_mmpm = log( (p-p0-q) * (p+p0-q) )
    select case (vec_type)
      case (1)
        select case (NRQCD_ORDER)
          case (0)
            v = CF*a * 2.*pi*(log_mppp-log_mmpm) * q/p
          case (1)
            log_mu_s = 2.*log(MU_SOFT)
            v = CF*a * (2.*(4.*pi+A1*a)*(log_mppp-log_mmpm) &
                      + B0*a*((log_mmpm-log_mu_s)**2-(log_mppp-log_mu_s)**2)) * q/(4.*p)
          case default
            call msg_fatal ("NRQCD_ORDER = " // char(NRQCD_ORDER))
        end select
      case (2)
        !!! not implemented yet
        v = zero
      case default
        call msg_fatal ("unknown ttZ/ttA vertex component, vec_type = " // char(vec_type))
    end select
  end function minus_q2_V

  !!! compute support points (~> q-grid) for numerical integration: trim p-grid and
  !!! merge with singular points of integrand: q = p, |p-p0|, p+p0, sqrt(mpole*E)
  subroutine compute_support_points (en, i_p, i_p0, n_trim)
    real(default), intent(in) :: en
    integer, intent(in) :: i_p
    integer, intent(in) :: i_p0
    integer, intent(in) :: n_trim
    real(default) :: p, p0
    real(default), dimension(4) :: sing_vals
    integer :: n_sing, i_q
    if (mod (POINTS_P, n_trim) /= 0) call msg_fatal ("trim p-grid for q-integration: POINTS_P = " &
                                  // char(POINTS_P) // " and n_trim = " // char(n_trim))
    n_q = POINTS_P / n_trim + merge(0,1,n_trim==1)
    p = p_grid(i_p)
    p0 = p0_grid(i_p0)
    n_sing = 0
    if ( i_p /= 1 .and. mod(i_p,n_trim) /= 0 ) then
      n_sing = n_sing+1
      sing_vals(n_sing) = p
    end if
    if ( i_p0 /= 1 ) then
      n_sing = n_sing+1
      sing_vals(n_sing) = p0 + p
      if ( i_p0 /= i_p+1 ) then
        n_sing = n_sing+1
        sing_vals(n_sing) = abs( p0 - p )
      end if
    end if
    if ( en > 0. ) then
      n_sing = n_sing+1
      sing_vals(n_sing) = sqrt( MTPOLE * en )
    end if
    if ( allocated(q_grid) ) deallocate( q_grid )
    allocate( q_grid(n_q+n_sing) )
    q_grid(1) = p_grid(1)
    q_grid(2:n_q) = [(p_grid(i_q), i_q=max(n_trim,2), POINTS_P, n_trim)]
    if (n_sing > 0 ) q_grid(n_q+1:n_q+n_sing) = sing_vals(1:n_sing)
    call nr_sort (q_grid)
  end subroutine compute_support_points

  !!! cf. arXiv:hep-ph/9503238, validated against arXiv:hep-ph/0008171
  pure function formfactor_ttv_relativistic_nlo (alphas, ps, J0) result (c)
    real(default), intent(in) :: alphas
    type(phase_space_point_t), intent(in) :: ps
    complex(default), intent(in) :: J0
    complex(default) :: c
    real(default) :: p2, k2, q2, kp, pq, kq
    complex(default) :: D2, chi, ln1, ln2, L1, L2, z, S, m2, m
    complex(default) :: JA, JB, JC, JD, JE, IA, IB, IC, ID, IE
    complex(default) :: CCmsbar
    complex(default) :: dF1, dF2, dM1, dM2
    complex(default), dimension(12) :: P1
    complex(default), parameter :: ximo = zero
    p2 = ps%p2
    k2 = ps%k2
    q2 = ps%q2
    m2 =  complex_m2 (ps%mpole, GAM)
    !!! kinematic abbreviations
    kp = 0.5_default * (-q2 + p2 + k2)
    pq = 0.5_default * ( k2 - p2 - q2)
    kq = 0.5_default * (-p2 + k2 + q2)
    D2 = kp**2 - k2*p2
    chi = p2*k2*q2 + 2.*m2*((p2 + k2)*kp - 2.*p2*k2) + m2**2 * q2
    ln1 = log( (1. - p2/m2)*(1,0) + ieps )
    ln2 = log( (1. - k2/m2)*(1,0) + ieps )
    L1 = (1. - m2/p2) * ln1
    L2 = (1. - m2/k2) * ln2
    z = sqrt( (1.-4.*m2/q2)*(1,0) )
    S = 0.5_default * z * log( (z+1.)/(z-1.) + ieps )
    m = sqrt(m2)

    !!! loop integrals in terms of J0
    JA = 1./D2 * (J0/2.*(-m2*pq - p2*kq) + kp*L2 - p2*L1 - 2.*pq*S)
    JB = 1./D2 * (J0/2.*( m2*kq + k2*pq) + kp*L1 - k2*L2 + 2.*kq*S)
    JC = 1/(4.*D2) * (2.*p2 + 2*kp*m2/k2 - 4.*kp*S + 2.*kp*(1. - m2/k2)*L2 + &
            (2.*kp*(p2 - m2) + 3.*p2*(m2 - k2))*JA + p2*(m2 - p2)*JB)
    JD = 1./(4.*D2) * (2.*kp*((k2 - m2)*JA + (p2 - m2)*JB - 1.) - k2*(2.*m2/k2 &
            - 2.*S + (1. - m2/k2)*L2 + (p2 - m2)*JA) - p2*(-2.*S + (1. - &
            m2/p2)*L1 + (k2 - m2)*JB))
    JE = 1./(4.*D2) * (2.*k2 + 2*kp*m2/p2 - 4.*kp*S + 2.*kp*(1. - m2/p2)*L1 + &
            (2.*kp*(k2 - m2) + 3.*k2*(m2 - p2))*JB + k2*(m2 - k2)*JA)
    IA = 1./D2 * (-(kq/2.)*J0 - 2.*q2/chi *((m2 - p2)*k2 - (m2 - k2)*kp)*S + &
            1./(m2 - p2)*(p2 - kp + p2*q2/chi *(k2 - m2)*(m2 + kp))*L1 + &
            k2*q2/chi *(m2 + kp)*L2)
    IB = 1./D2 * ( (pq/2.)*J0 - 2.*q2/chi *((m2 - k2)*p2 - (m2 - p2)*kp)*S + &
            1./(m2 - k2)*(k2 - kp + k2*q2/chi *(p2 - m2)*(m2 + kp))*L2 + &
            p2*q2/chi *(m2 + kp)*L1)
    IC = 1./(4.*D2) * (2.*p2*J0 - 4.*kp/k2*(1. + m2/(k2 - m2)*L2) + (2.*kp - &
            3.*p2)*JA - p2*JB + (-2.*kp*(m2 - p2) + 3.*p2*(m2 - k2))*IA + &
            p2*(m2 - p2)*IB)
    ID = 1./(4.*D2) * (-2.*kp*J0 + 2.*(1. + m2/(k2 - m2)*L2) + 2.*(1. + &
            m2/(p2 - m2)*L1) + (2.*kp - k2)*JA + (2.*kp - p2)*JB + (k2*(m2 - &
            p2) - 2.*kp*(m2 - k2))*IA + (p2*(m2 - k2) - 2.*kp*(m2 - p2))*IB)
    IE = 1./(4.*D2) * (2.*k2*J0 - 4.*kp/p2*(1. + m2/(p2 - m2)*L1) + (2.*kp - &
            3.*k2)*JB - k2*JA + (-2.*kp*(m2 - k2) + 3.*k2*(m2 - p2))*IB + &
            k2*(m2 - k2)*IA)

    !!! divergent part ~ 1/epsilon: depends on subtraction scheme
    CCmsbar = -2.0_default * log(RESCALE_H)

    ! real top mass in the loop numerators
!    m2 = cmplx(real(m2), kind=default)
!    m  = sqrt(m2)

    !!! quark self energies
    dF1 = - (ximo+1.) * (CCmsbar + (1.+m2/p2)*(1.-L1))
    dF2 = - (ximo+1.) * (CCmsbar + (1.+m2/k2)*(1.-L2))
    dM1 = m/p2 * ( (ximo+1.)*(1.+m2/p2*ln1) - 3.*ln1 )
    dM2 = m/k2 * ( (ximo+1.)*(1.+m2/k2*ln2) - 3.*ln2 )

    !!! coefficient list: vertex function Gamma_mu (k,p) = sum_i( Vi_mu * Pi )
    P1(1)  =  2.*JA - 2.*JC + ximo*(m2*IC + p2*ID)
    P1(2)  =  2.*JB - 2.*JE + ximo*(k2*ID + m2*IE)
    P1(3)  = -2.*J0 + 2.*JA + 2.*JB - 2.*JD + ximo*(-J0/2. - k2/2.*IC - &
                 kp*ID + m2*ID + p2/2.*IE + JA)
    P1(4)  = -2.*JD + ximo*(k2*IC + m2*ID - JA)
    P1(5)  = J0 - JA - JB + ximo*(J0/4. + k2/4.*IC + kp/2.*ID + p2/4.*IE - &
                 1./2.*JA - 1./2.*JB)
    P1(6)  = -m2*J0 - k2*JA - p2*JB + k2/2.*JC + kp*JD + p2/2.*JE + &
                 (1./2. + CCmsbar - 2.*S) &
                 + ximo*(-m2*J0/4. - m2/4.*k2*IC - m2/2.*kp*ID - m2/4.*p2*IE &
                 - k2/2.*JA - p2/2.*JB + (CCmsbar + 2.))
    P1(7)  =  2.*m*J0 - 4.*m*JA + ximo*m*(J0/2. - 2.*kp*IC + k2/2.*IC - &
                 p2*ID - kp*ID - p2/2.*IE - JA)
    P1(8)  =  2.*m*J0 - 4.*m*JB + ximo*m*(J0/2. + k2/2.*IC - kp*ID + k2*ID - &
                 p2/2.*IE - JB)
    P1(9)  =  ximo*m*(ID + IE)
    P1(10) =  ximo*m*(ID + IC)
    P1(11) =  ximo*m*( p2*ID + kp*IC + p2/2.*IE - k2/2.*IC) + dM2
                                 !!! self energy contribution: ~ gamma_mu.k_slash = V11
    P1(12) =  ximo*m*(-k2*ID - kp*IE + p2/2.*IE - k2/2.*IC) + dM1
                                 !!! self energy contribution: ~ gamma_mu.p_slash = V12

    !!! leading form factor: V6 = gamma_mu, V5 = gamma_mu.k_slash.p_slash ~> -m^2*gamma_mu
    c = one + alphas * CF / (4.*pi) * ( P1(6) - m2*P1(5) &
                 !!! self energy contributions ~ gamma^mu
                 + dF1 + dF2 + m*( dM1 + dM2 ) )
                 !!! on-shell subtraction: UV divergence cancels
!                 + 0.5_default*( dF1 + dF2 + m*( dM1 + dM2 ) )
  end function formfactor_ttv_relativistic_nlo

  pure function sqrts_to_en (sqrts, mpole_in) result (en)
    real(default), intent(in) :: sqrts
    real(default), intent(in), optional :: mpole_in
    real(default) :: mpole, en
    if (present (mpole_in)) then
      mpole = mpole_in
    else
      mpole = m1s_to_mpole (sqrts)
    end if
    en = sqrts - two * mpole
  end function sqrts_to_en

  function p_grid_from_TOPPIK (mpole_in) result (p_toppik)
    real(default), intent(in), optional :: mpole_in
    real(default), dimension(POINTS_P) :: p_toppik
    real(default) :: mpole
    if (debug_on) call msg_debug (D_THRESHOLD, "p_grid_from_TOPPIK")
    mpole = MTPOLE;  if (present (mpole_in))  mpole = mpole_in
    call scan_formfactor_over_p_TOPPIK &
                 (alphas_soft(2. * M1S), 2. * M1S, 1, p_toppik, mpole)
    if (.not. strictly_monotonous (p_toppik)) &
      call msg_fatal ("p_grid NOT strictly monotonous!")
  end function p_grid_from_TOPPIK

  pure module function int_to_char (i) result (c)
    integer, intent(in) :: i
    character(len=len(trim(int2fixed(i)))) :: c
    c = int2char (i)
  end function int_to_char

  pure module function real_to_char (r) result (c)
    real(default), intent(in) :: r
    character(len=len(trim(real2fixed(r)))) :: c
    c = real2char (r)
  end function real_to_char

  pure module function complex_to_char (z) result (c)
    complex(default), intent(in) :: z
    character(len=len(trim(real2fixed(real(z))))+len(trim(real2fixed(aimag(z))))+5) :: c
    character(len=len(trim(real2fixed(real(z))))) :: re
    character(len=len(trim(real2fixed(aimag(z))))) :: im
    re = real_to_char (real(z))
    im = real_to_char (aimag(z))
    if (nearly_equal (aimag(z), zero)) then
      c = re
    else
      c = re // " + " // im // "*I"
    end if
  end function complex_to_char

  pure module function logical_to_char (l) result (c)
    logical, intent(in) :: l
    character(len=1) :: c
    write (c, '(l1)') l
  end function logical_to_char

  subroutine get_rest_frame (p1_in, p2_in, p1_out, p2_out)
    type(vector4_t), intent(in) :: p1_in, p2_in
    type(vector4_t), intent(out) :: p1_out, p2_out
    type(lorentz_transformation_t) :: L
    L = inverse (boost (p1_in + p2_in, (p1_in + p2_in)**1))
    p1_out = L * p1_in; p2_out = L * p2_in
  end subroutine get_rest_frame

  function shift_momentum (p_in, E, p) result (p_out)
    type(vector4_t) :: p_out
    type(vector4_t), intent(in) :: p_in
    real(default), intent(in) :: E, p
    type(vector3_t) :: vec
    vec = p_in%p(1:3) / space_part_norm (p_in)
    p_out = vector4_moving (E, p * vec)
  end function shift_momentum

  subroutine evaluate_one_to_two_splitting_threshold (p_origin, &
      p1_in, p2_in, p1_out, p2_out, msq_in, jac)
    type(vector4_t), intent(in) :: p_origin
    type(vector4_t), intent(in) :: p1_in, p2_in
    type(vector4_t), intent(inout) :: p1_out, p2_out
    real(default), intent(in), optional :: msq_in
    real(default), intent(inout), optional :: jac
    type(lorentz_transformation_t) :: L
    type(vector4_t) :: p1_rest, p2_rest
    real(default) :: msq, msq1, msq2
    real(default) :: m
    real(default) :: E1, E2, E_max
    real(default) :: p, lda
    real(default), parameter :: E_offset = 0.001_default
    !!! (TODO-cw-2016-10-13) Find a better way to get masses
    real(default), parameter :: mb = 4.2_default
    real(default), parameter :: mw = 80.419_default

    call get_rest_frame (p1_in, p2_in, p1_rest, p2_rest)

    msq = p_origin**2; m = sqrt(msq)
    msq1 = p1_in**2
    msq2 = m * (m - two * p1_rest%p(0))
    E1 = (msq + msq1 - msq2) / (two * m)
    E_max = (msq - (mb + mw)**2) / (two * m)
    E_max = E_max - E_offset
    if (E1 > E_max) then
       E1 = E_max
       msq2 = m * (m - two * E_max)
    end if

    lda = lambda (msq, msq1, msq2)
    if (lda < zero) call msg_fatal &
         ("Threshold Splitting: lambda < 0 encountered! Use a higher offset.")
    p = sqrt(lda) / (two * m)

    E1 = sqrt (msq1 + p**2)
    E2 = sqrt (msq2 + p**2)

    p1_out = shift_momentum (p1_rest, E1, p)
    p2_out = shift_momentum (p2_rest, E2, p)

    L = boost (p_origin, p_origin**1)
    p1_out = L  * p1_out
    p2_out = L  * p2_out
  end subroutine evaluate_one_to_two_splitting_threshold

  module subroutine generate_on_shell_decay_threshold (p_decay, p_top, p_decay_onshell)
    !!! Gluon must be on first position in this array
    type(vector4_t), intent(in), dimension(:) :: p_decay
    type(vector4_t), intent(inout) :: p_top
    type(vector4_t), intent(inout), dimension(:) :: p_decay_onshell
    procedure(evaluate_one_to_two_splitting_special), pointer :: ppointer
    ppointer => evaluate_one_to_two_splitting_threshold
    call generate_on_shell_decay (p_top, p_decay, p_decay_onshell, 1, &
         evaluate_special = ppointer)
  end subroutine generate_on_shell_decay_threshold


end submodule ttv_formfactors_s

