module @ID@_virtual
  use iso_c_binding
  use kinds
  use diagnostics
  use openloops
  use parameters_sm_tt_threshold
  implicit none
  integer, dimension(4) :: id
end module @ID@_virtual

subroutine @ID@_start_openloops () bind(C)
  use @ID@_virtual
  use openloops_blha, only: olp_printparameter
  integer, dimension(12), parameter :: pdgs = [1,2,3,4,5,6,11,13,15,23,24,25]
  character(32) :: buffer
  integer :: i
  !!! Need to set this parameter to switch off OpenLoops' internal CS I-operators,
  !!! which are not implemented for decays.
  call set_parameter ("ir_on", 0)
  !!! coupling order alpha_ew^1, implies QCD correction for loop process
  call set_parameter("order_ew", 1)
  call set_parameter("model", "sm")
  do i = 1, size(pdgs)
     if (pdgs(i) > 4) then
        write (buffer, "(I0)")  pdgs(i)
        call set_parameter("mass(" // trim (adjustl (buffer)) // ")", mass(pdgs(i)))
        call set_parameter("width(" // trim (adjustl (buffer)) // ")", width(pdgs(i)))
     end if
  end do
  call set_parameter("alpha", real(1 / alphaemi, kind=double))
  ! Increase verbosity level to list loaded libraries
  call set_parameter("verbose", 3)
  call set_parameter("psp_tolerance", 10e-7_double)
  call set_parameter("use_cms", 0)
  call set_parameter("me_cache", 0)
  ! 1 for tree-like matrix elements (tree, color and spin correlations),
  ! 11 for loop, 12 for loop^2
  id(1) = register_process("6(-1) -> 24 5", 11)
  id(2) = register_process("6(+1) -> 24 5", 11)
  !!! Anti-particles need opposite sign for the helicity
  !!! Different conventions in OMega and OpenLoops
  id(3) = register_process("-6(+1) -> -24 -5", 11)
  id(4) = register_process("-6(-1) -> -24 -5", 11)
  call start()
  call olp_printparameter ("parameters.ol")
end subroutine @ID@_start_openloops

subroutine @ID@_olp_eval2 (i_flv, alpha_s_c, p_ofs, mu_c, &
       sel_hel_beam, sqme_c, acc_c) bind(C)
  use @ID@_threshold
  use @ID@_virtual
  use debug_master, only: debug_on
  use physics_defs, only: ass_boson, ass_quark
  use physics_defs, only: THR_POS_WP, THR_POS_WM
  use physics_defs, only: THR_POS_B, THR_POS_BBAR
  use ttv_formfactors
  use omega95
  implicit none
  integer(c_int), intent(in) :: i_flv
  real(c_default_float), intent(in) :: alpha_s_c
  real(c_default_float), dimension(0:3,*), intent(in) :: p_ofs
  real(c_default_float), intent(in) :: mu_c
  integer, intent(in) :: sel_hel_beam
  real(c_default_float), dimension(4), intent(out) :: sqme_c
  real(c_default_float), intent(out) :: acc_c
  type(momentum), dimension(:), allocatable :: mom_ofs, mom_ons
  type(momentum), dimension(2) :: ptop_ofs, ptop_ons
  type(momentum) :: p12
  complex(default), dimension(-1:1,-1:1,-1:1,-1:1) :: production_me
  integer :: h_el, h_pos, h_t, h_tbar
  integer :: leg, this_id, h_b, h_W
  real(double) :: virt_decay(0:3), total(0:3), acc
  real(double) :: p_decay(0:3,3,2)
  real(double) :: mu, alpha_s, dynamic_top_mass
  complex(default) :: born_decay_me, bw
  real(default) :: prod2, born_decay_me2
  logical :: eval_this_beam_helicities
  integer, dimension(2) :: sel_hel
  if (debug_on) call msg_debug (D_ME_METHODS, "@ID@_olp_eval2")
  if (i_flv /= 1)  call msg_fatal ("i_flv /= 1, threshold interface was not built for this")
  if (any (id <= 0))  call msg_fatal ("Could not register process in OpenLoops")
  if (.not. threshold%settings%factorized_computation)  &
       call msg_fatal ("@ID@_olp_eval2: OFFSHELL_STRATEGY is not factorized")
  alpha_s = alpha_s_c
  mu = mu_c
  allocate (mom_ofs(6), mom_ons(6))
  call convert_to_mom_and_invert_sign (p_ofs, 6, mom_ofs)
  call set_parameter("alpha_s", alpha_s)
  call set_parameter("mu", mu)
  total = 0

  sel_hel = get_selected_beam_helicities (sel_hel_beam)
  call set_production_factors (sel_hel_beam)

  call compute_projected_momenta (mom_ofs, mom_ons)
  call set_parameter("width(6)", zero)
  call set_parameter("width(24)", zero)
  ptop_ons(1) = mom_ons(THR_POS_WP) + mom_ons(THR_POS_B)
  ptop_ons(2) = mom_ons(THR_POS_WM) + mom_ons(THR_POS_BBAR)
  ptop_ofs(1) = mom_ofs(THR_POS_WP) + mom_ofs(THR_POS_B)
  ptop_ofs(2) = mom_ofs(THR_POS_WM) + mom_ofs(THR_POS_BBAR)
  p_decay(:,1,1) = ptop_ons(1)
  p_decay(:,2,1) = mom_ons(THR_POS_WP)
  p_decay(:,3,1) = mom_ons(THR_POS_B)
  p_decay(:,1,2) = ptop_ons(2)
  p_decay(:,2,2) = mom_ons(THR_POS_WM)
  p_decay(:,3,2) = mom_ons(THR_POS_BBAR)
  if (.not. threshold%settings%onshell_projection%boost_decay &
         .and. debug2_active (D_THRESHOLD)) then
     p12 = mom_ofs(1) + mom_ofs(2)
     call check_rest_frame (ttv_mtpole (p12 * p12))
  end if
  bw = top_propagators (FF, mom_ofs(1) + mom_ofs(2), ptop_ofs)
  production_me = compute_production_me (FF, mom_ofs, ptop_ofs, ptop_ons)
  prod2 = zero
  do h_t = -1, 1, 2
  do h_tbar = -1, 1, 2
     if (skip (h_t, h_tbar)) cycle
     do h_el = -1, 1, 2
     do h_pos = -1, 1, 2
        eval_this_beam_helicities = sel_hel_beam < 0 .or. &
            (sel_hel(1) == h_el .and. sel_hel(2) == h_pos)
        if (.not. eval_this_beam_helicities) cycle
        prod2 = abs2 (production_me(h_el, h_pos, h_t, h_tbar))
     do leg = 1, 2
        dynamic_top_mass = sqrt (ptop_ons(leg) * ptop_ons(leg))
        call set_parameter("mass(6)", dynamic_top_mass)
        born_decay_me2 = zero
        do h_W = -1, 1
        do h_b = -1, 1, 2
           if (leg == 1) then
              born_decay_me = anti_top_decay_born (mom_ons, h_tbar, h_W, h_b)
              this_id = (3 + h_t) / 2
           else
              born_decay_me = top_decay_born (mom_ons, h_t, h_W, h_b)
              this_id = (3 + h_tbar) / 2 + 2
           end if
           born_decay_me2 = born_decay_me2 + abs2 (born_decay_me)
        end do
        end do
        call evaluate_loop (id(this_id), p_decay(:,:,leg), &
             virt_decay(3), virt_decay(0:2), acc)
        acc_c = max (acc_c, acc)
        total = total + prod2 * born_decay_me2 * virt_decay * abs2 (bw)
     end do
  end do
  end do
  end do
  end do
  !!! OpenLoops applies the averaging factor of two regardless whether
  !!! the MEs are polarized or not
  total = total * production_factors
  sqme_c = [total(2), total(1), total(0), total(3)]
  if (debug2_active (D_ME_METHODS)) &
       print *, 'sqme_c =    ', sqme_c !!! Debugging
contains
  function skip (h_t, h_tbar)
    logical :: skip
    integer, intent(in) :: h_t, h_tbar
    associate (s => threshold%settings)
         skip = s%helicity_approximation%ultra &
              .and. (h_t /= s%sel_hel_top .or. h_tbar /= s%sel_hel_topbar)
    end associate
  end function skip

  subroutine check_rest_frame (mtop)
    use, intrinsic :: iso_fortran_env, only: output_unit
    use constants, only: tiny_07
    use numeric_utils, only: assert_equal
    real(default), intent(in) :: mtop
    call assert_equal (output_unit, mtop, sqrt (ptop_ons(1) * ptop_ons(1)), &
         'Top quark in rest frame?', abs_smallness = tiny_07, &
         rel_smallness = tiny_07, exit_on_fail = .true.)
    call assert_equal (output_unit, mtop, sqrt (ptop_ons(2) * ptop_ons(2)), &
         'Anti-top quark in rest frame?', abs_smallness = tiny_07, &
         rel_smallness = tiny_07, exit_on_fail = .true.)
  end subroutine check_rest_frame
end subroutine @ID@_olp_eval2

subroutine @ID@_stop_openloops () bind(C)
  use iso_c_binding
  use openloops
  call finish()
end subroutine @ID@_stop_openloops
