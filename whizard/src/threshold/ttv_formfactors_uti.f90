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

module ttv_formfactors_uti

  use kinds, only: default
  use debug_master, only: debug_on
  use constants
  use ttv_formfactors
  use diagnostics
  use sm_physics, only: running_as
  use numeric_utils

  implicit none
  private

  public :: ttv_formfactors_1
  public :: ttv_formfactors_2

contains

  subroutine ttv_formfactors_1 (u)
    integer, intent(in) :: u
    real(default) :: m1s, Vtb, wt_inv, alphaemi, sw, alphas_mz, mz, &
         mw, mb, sh, sf, NRQCD_ORDER, FF, offshell_strategy, v1, v2, &
         scan_sqrts_max, sqrts, scan_sqrts_min, scan_sqrts_stepsize, &
         test, gam_out, mpole
    type(formfactor_t) :: formfactor
    type(phase_space_point_t) :: ps
    logical :: mpole_fixed
    integer :: top_helicity_selection
    write (u, "(A)")  "* Test output: ttv_formfactors_1"
    write (u, "(A)")  "*   Purpose: Basic setup"
    write (u, "(A)")

    m1s = 172.0_default
    Vtb = one
    wt_inv = zero
    alphaemi = 125.0_default
    alphas_mz = 0.118_default
    mz = 91.1876_default
    mw = 80.399_default
    sw = sqrt(one - mw**2 / mz**2)
    mb = 4.2_default
    sh = one
    sf = one
    NRQCD_ORDER = one
    FF = MATCHED
    offshell_strategy = 0
    top_helicity_selection = -1
    v1 = 0.3_default
    v2 = 0.5_default
    scan_sqrts_stepsize = 0.0_default
    test = - one
    write (u, "(A)") "Check high energy behavior"
    sqrts = 500.0_default
    scan_sqrts_min = sqrts
    scan_sqrts_max = sqrts
    write (u, "(A)") "Check that the mass is not fixed"
    mpole_fixed = .false.

    call init_parameters &
         (mpole, gam_out, m1s, Vtb, wt_inv, &
          alphaemi, sw, alphas_mz, mz, mw, &
          mb, sh, sf, NRQCD_ORDER, FF, offshell_strategy, &
          v1, v2, scan_sqrts_min, scan_sqrts_max, &
          scan_sqrts_stepsize, mpole_fixed, real(top_helicity_selection, default))
    call init_threshold_grids (test)
    call threshold%formfactor%activate ()
    call formfactor%activate ()
    call assert (u, m1s_to_mpole (350.0_default) > m1s + 0.1_default, &
         "m1s_to_mpole (350.0_default) > m1s")
    write (u, "(A)")

    ! For simplicity we test on-shell back-to-back tops
    call ps%init (m1s**2, m1s**2, sqrts**2, mpole)
    call assert_equal (u, f_switch_off (v_matching (ps%sqrts, GAM_M1S)), tiny_10, &
         "f_switch_off (v_matching (ps%sqrts, GAM_M1S))")
    call assert (u, &
         abs (formfactor%compute (ps, 1, EXPANDED_HARD)) > &
         abs (formfactor%compute (ps, 1, RESUMMED)), &
         "expansion with hard alphas should be larger " // &
         "than resummed (with switchoff)")
    call assert_equal (u, &
         abs (formfactor%compute (ps, 1, RESUMMED)), zero, &
         "resummed (with switchoff) should be zero", abs_smallness=tiny_10)
    call assert_equal (u, &
         abs (formfactor%compute (ps, 1, EXPANDED_SOFT_SWITCHOFF)), zero, &
         "expanded (with switchoff) should be zero", abs_smallness=tiny_10)
    write (u, "(A)") ""

    write (u, "(A)") "Check global variables"
    call assert_equal (u, AS_HARD, &
         running_as (m1s, alphas_mz, mz, 2, 5.0_default), "hard alphas")
    call assert_equal (u, AS_SOFT, zero, "soft alphas", abs_smallness=tiny_10)
    call assert_equal (u, AS_USOFT, zero, "ultrasoft alphas", abs_smallness=tiny_10)
    call assert_equal (u, AS_LL_SOFT, zero, "LL soft alphas", abs_smallness=tiny_10)

    !!! care: the formfactor contains the tree level that we usually subtract again
    write (u, "(A)") "Check low energy behavior"
    sqrts = 2 * m1s + 0.01_default
    scan_sqrts_min = sqrts
    scan_sqrts_max = sqrts
    write (u, "(A)") "Check that the mass is fixed"
    mpole_fixed = .true.
    call init_parameters &
         (mpole, gam_out, m1s, Vtb, wt_inv, &
          alphaemi, sw, alphas_mz, mz, mw, &
          mb, sh, sf, NRQCD_ORDER, FF, offshell_strategy, &
          v1, v2, scan_sqrts_min, scan_sqrts_max, &
          scan_sqrts_stepsize, mpole_fixed, real(top_helicity_selection, default))
    call init_threshold_grids (test)

    call ps%init (m1s**2, m1s**2, sqrts**2, mpole)
    call assert_equal (u, m1s_to_mpole (350.0_default), m1s, &
         "m1s_to_mpole (350.0_default) == m1s")
    call assert_equal (u, m1s_to_mpole (550.0_default), m1s, &
         "m1s_to_mpole (550.0_default) == m1s")
    write (u, "(A)") ""
    call assert_equal (u, f_switch_off (v_matching (ps%sqrts, GAM_M1S)), &
         one, "f_switch_off (v_matching (ps%sqrts, GAM_M1S))")
    call formfactor%disable ()
    call assert_equal (u, &
         abs(formfactor%compute (ps, 1, 1)), &
         zero, &
         "disabled formfactor should return zero")
    call formfactor%activate ()
    call assert_equal (u, &
         formfactor%compute (ps, 1, EXPANDED_SOFT_SWITCHOFF), &
         formfactor%compute (ps, 1, EXPANDED_SOFT), &
         "switchoff function should do nothing here")
    write (u, "(A)") ""

    write (u, "(A)")  "* Test output end: ttv_formfactors_1"
  end subroutine ttv_formfactors_1

  subroutine ttv_formfactors_2 (u)
    integer, intent(in) :: u
    write (u, "(A)")  "* Test output: ttv_formfactors_2"
    write (u, "(A)")  "*   Purpose: Test flags"
    write (u, "(A)")

    write (u, "(A)") "RESUMMED_SWITCHOFF + NLO"
    call threshold%settings%setup_flags (-2, 1, -1)
    call assert (u, SWITCHOFF_RESUMMED, "SWITCHOFF_RESUMMED")
    call assert (u, TOPPIK_RESUMMED, "TOPPIK_RESUMMED")
    call assert (u, threshold%settings%nlo, "threshold%settings%nlo")
    call assert (u, .not. threshold%settings%factorized_computation, &
         ".not. threshold%settings%factorized_computation")
    call assert (u, .not. threshold%settings%interference, &
         ".not. threshold%settings%interference")
    call assert (u, .not. threshold%settings%no_nlo_width_in_signal_propagators, &
         ".not. threshold%settings%no_nlo_width_in_signal_propagators")

    write (u, "(A)") "MATCHED + FACTORIZATION"
    call threshold%settings%setup_flags (-1, 0+2, -1)
    call assert (u, .not. threshold%settings%nlo, ".not. threshold%settings%nlo")
    call assert (u, TOPPIK_RESUMMED, "TOPPIK_RESUMMED")
    call assert (u, threshold%settings%factorized_computation, &
         "threshold%settings%factorized_computation")

    write (u, "(A)") "RESUMMED + INTERFERENCE"
    call threshold%settings%setup_flags (1, 0+0+4, -1)
    call assert (u, .not. SWITCHOFF_RESUMMED, ".not. SWITCHOFF_RESUMMED")
    call assert (u, TOPPIK_RESUMMED, "TOPPIK_RESUMMED")
    call assert (u, .not. threshold%settings%nlo, ".not. threshold%settings%nlo")
    call assert (u, .not. threshold%settings%factorized_computation, &
         ".not. threshold%settings%factorized_computation")
    call assert (u, threshold%settings%interference, "threshold%settings%interference")

    write (u, "(A)") "EXPANDED_HARD"
    call threshold%settings%setup_flags (4, 0+2+4, -1)
    call assert (u, .not. SWITCHOFF_RESUMMED, ".not. SWITCHOFF_RESUMMED")
    call assert (u, .not. TOPPIK_RESUMMED, ".not. TOPPIK_RESUMMED")
    call assert (u, .not. threshold%settings%nlo, ".not. threshold%settings%nlo")
    call assert (u, threshold%settings%factorized_computation, &
         "threshold%settings%factorized_computation")
    call assert (u, threshold%settings%interference, "threshold%settings%interference")

    write (u, "(A)") "EXPANDED_SOFT"
    call threshold%settings%setup_flags (5, 1+2+4, -1)
    call assert (u, .not. SWITCHOFF_RESUMMED, ".not. SWITCHOFF_RESUMMED")
    call assert (u, .not. TOPPIK_RESUMMED, ".not. TOPPIK_RESUMMED")
    call assert (u, threshold%settings%nlo, "threshold%settings%nlo")
    call assert (u, threshold%settings%factorized_computation, &
         "threshold%settings%factorized_computation")
    call assert (u, threshold%settings%interference, &
         "threshold%settings%interference")

    write (u, "(A)") "EXPANDED_SOFT_SWITCHOFF"
    call threshold%settings%setup_flags (6, 0+0+0+8, -1)
    call assert (u, .not. SWITCHOFF_RESUMMED, "SWITCHOFF_RESUMMED")
    call assert (u, .not. TOPPIK_RESUMMED, ".not. TOPPIK_RESUMMED")
    call assert (u, .not. threshold%settings%nlo, "threshold%settings%nlo")
    call assert (u, .not. threshold%settings%factorized_computation, &
         "threshold%settings%factorized_computation")
    call assert (u, .not. threshold%settings%interference, &
         "threshold%settings%interference")

    write (u, "(A)") "RESUMMED_ANALYTIC_LL"
    call threshold%settings%setup_flags (7, 0+0+4+8, -1)
    call assert (u, .not. SWITCHOFF_RESUMMED, "SWITCHOFF_RESUMMED")
    call assert (u, .not. TOPPIK_RESUMMED, ".not. TOPPIK_RESUMMED")
    call assert (u, .not. threshold%settings%nlo, "threshold%settings%nlo")
    call assert (u, .not. threshold%settings%factorized_computation, &
         "threshold%settings%factorized_computation")
    call assert (u, threshold%settings%interference, "threshold%settings%interference")
    call assert (u, threshold%settings%onshell_projection%production, &
         "threshold%settings%onshell_projection%production")

    write (u, "(A)") "EXPANDED_SOFT_HARD"
    call threshold%settings%setup_flags (8, 0+2+0+128, -1)
    call assert (u, .not. SWITCHOFF_RESUMMED, "SWITCHOFF_RESUMMED")
    call assert (u, .not. TOPPIK_RESUMMED, ".not. TOPPIK_RESUMMED")
    call assert (u, .not. threshold%settings%nlo, "threshold%settings%nlo")
    call assert (u, threshold%settings%factorized_computation, &
         "threshold%settings%factorized_computation")
    call assert (u, .not. threshold%settings%interference, "threshold%settings%interference")
    call assert (u, .not. threshold%settings%onshell_projection%production, &
         "threshold%settings%onshell_projection%production")
    call assert (u, threshold%settings%onshell_projection%decay, &
         "threshold%settings%onshell_projection%decay")

    write (u, "(A)") "EXTRA_TREE"
    call threshold%settings%setup_flags (9, 1+0+0+16+64, -1)
    call assert (u, .not. SWITCHOFF_RESUMMED, "SWITCHOFF_RESUMMED")
    call assert (u, .not. TOPPIK_RESUMMED, ".not. TOPPIK_RESUMMED")
    call assert (u, threshold%settings%nlo, "threshold%settings%nlo")
    call assert (u, .not. threshold%settings%factorized_computation, &
         "threshold%settings%factorized_computation")
    call assert (u, .not. threshold%settings%interference, "threshold%settings%interference")
    call assert (u, threshold%settings%onshell_projection%production, &
         "threshold%settings%onshell_projection%production")
    call assert (u, .not. threshold%settings%onshell_projection%decay, &
         "threshold%settings%onshell_projection%decay")
    call assert (u, threshold%settings%no_nlo_width_in_signal_propagators, &
         "threshold%settings%no_nlo_width_in_signal_propagators")

    write (u, "(A)") "test projection of width"
    call threshold%settings%setup_flags (9, 0+0+0+0+256, -1)
    call assert (u, .not. threshold%settings%onshell_projection%production, &
         "threshold%settings%onshell_projection%production")
    call assert (u, .not. threshold%settings%onshell_projection%decay, &
         "threshold%settings%onshell_projection%decay")
    call assert (u, .not. threshold%settings%onshell_projection%width, &
         "threshold%settings%onshell_projection%width")

    write (u, "(A)") "test boost of decay momenta"
    call threshold%settings%setup_flags (9, 512, -1)
    if (debug_on) call msg_debug (D_THRESHOLD, &
         "threshold%settings%onshell_projection%boost_decay", &
         threshold%settings%onshell_projection%boost_decay)
    call threshold%settings%setup_flags (9, 0, -1)
    if (debug_on) call msg_debug (D_THRESHOLD, &
         ".not. threshold%settings%onshell_projection%boost_decay", &
         .not. threshold%settings%onshell_projection%boost_decay)

    write (u, "(A)") "test helicity approximations"
    call threshold%settings%setup_flags (9, 32, -1)
    call assert (u, threshold%settings%helicity_approximation%simple, &
         "threshold%settings%helicity_approximation%simple")
    call assert (u, .not. threshold%settings%helicity_approximation%extra, &
         ".not. threshold%settings%helicity_approximation%extra")
    call assert (u, .not. threshold%settings%helicity_approximation%ultra, &
         ".not. threshold%settings%helicity_approximation%ultra")
    call threshold%settings%setup_flags (9, 1024, -1)
    call assert (u, .not. threshold%settings%helicity_approximation%simple, &
         ".not. threshold%settings%helicity_approximation%simple")
    call assert (u, threshold%settings%helicity_approximation%extra, &
         "threshold%settings%helicity_approximation%extra")

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: ttv_formfactors_2"
  end subroutine ttv_formfactors_2


end module ttv_formfactors_uti
