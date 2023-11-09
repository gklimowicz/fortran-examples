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

submodule (kinematics) kinematics_s

  use debug_master, only: debug_on
  use format_utils, only: write_separator
  use diagnostics
  use io_units
  use phs_points, only: assignment(=), size
  use interactions
  use phs_fks
  use ttv_formfactors, only: m1s_to_mpole

  implicit none

contains

  module subroutine kinematics_write (object, unit)
    class(kinematics_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, c
    u = given_output_unit (unit)
    if (object%f_allocated) then
       write (u, "(1x,A)")  "Flux * PHS volume:"
       write (u, "(2x,ES19.12)")  object%phs_factor
       write (u, "(1x,A)")  "Jacobian factors per channel:"
       do c = 1, size (object%f)
          write (u, "(3x,I0,':',1x,ES14.7)", advance="no")  c, object%f(c)
          if (c == object%selected_channel) then
             write (u, "(1x,A)")  "[selected]"
          else
             write (u, *)
          end if
       end do
    end if
    if (object%sf_chain_allocated) then
       call write_separator (u)
       call object%sf_chain%write (u)
    end if
    if (object%phs_allocated) then
       call write_separator (u)
       call object%phs%write (u)
    end if
  end subroutine kinematics_write

  module subroutine kinematics_final (object)
    class(kinematics_t), intent(inout) :: object
    if (object%sf_chain_allocated) then
       call object%sf_chain%final ()
       deallocate (object%sf_chain)
       object%sf_chain_allocated = .false.
    end if
    if (object%phs_allocated) then
       call object%phs%final ()
       deallocate (object%phs)
       object%phs_allocated = .false.
    end if
    if (object%f_allocated) then
       deallocate (object%f)
       object%f_allocated = .false.
    end if
  end subroutine kinematics_final

  module subroutine kinematics_configure (kin, pcm, pcm_work, &
       sf_chain, beam_config, phs_config, nlo_type, is_i_sub)
    class(kinematics_t), intent(out) :: kin
    class(pcm_t), intent(inout) :: pcm
    class(pcm_workspace_t), intent(in) :: pcm_work
    type(sf_chain_t), intent(in), target :: sf_chain
    type(process_beam_config_t), intent(in), target :: beam_config
    class(phs_config_t), intent(in), target :: phs_config
    integer, intent(in) :: nlo_type
    logical, intent(in) :: is_i_sub
    logical :: extended_sf
    extended_sf = nlo_type == NLO_DGLAP .or. &
         (nlo_type == NLO_REAL .and. is_i_sub)
    call kin%init_sf_chain (sf_chain, beam_config, &
         extended_sf = pcm%has_pdfs .and. extended_sf)
    !!! Add one for additional Born matrix element
    call kin%init_phs (phs_config)
    call kin%set_nlo_info (nlo_type)
    select type (phs => kin%phs)
    type is (phs_fks_t)
       call phs%allocate_momenta (phs_config, .not. (nlo_type == NLO_REAL))
       select type (pcm)
       type is (pcm_nlo_t)
          call pcm%region_data%init_phs_identifiers (phs%phs_identifiers)
          !!! The triple select type pyramid of doom
          select type (pcm_work)
          type is (pcm_nlo_workspace_t)
             if (allocated (pcm_work%real_kinematics%alr_to_i_phs)) &
                  call pcm%region_data%set_alr_to_i_phs (phs%phs_identifiers, &
                       pcm_work%real_kinematics%alr_to_i_phs)
          end select
       end select
    end select
  end subroutine kinematics_configure

  module subroutine kinematics_set_nlo_info (k, nlo_type)
    class(kinematics_t), intent(inout) :: k
    integer, intent(in) :: nlo_type
    if (nlo_type == NLO_VIRTUAL)  k%only_cm_frame = .true.
  end subroutine kinematics_set_nlo_info

  module subroutine kinematics_set_threshold (kin, factorization_mode)
    class(kinematics_t), intent(inout) :: kin
    integer, intent(in) :: factorization_mode
    kin%threshold = factorization_mode == FACTORIZATION_THRESHOLD
  end subroutine kinematics_set_threshold

  module subroutine kinematics_init_sf_chain (k, sf_chain, config, extended_sf)
    class(kinematics_t), intent(inout) :: k
    type(sf_chain_t), intent(in), target :: sf_chain
    type(process_beam_config_t), intent(in) :: config
    logical, intent(in), optional :: extended_sf
    integer :: n_strfun, n_channel
    integer :: c
    k%n_in = config%data%get_n_in ()
    n_strfun = config%n_strfun
    n_channel = config%n_channel
    allocate (k%sf_chain)
    k%sf_chain_allocated = .true.
    call k%sf_chain%init (sf_chain, n_channel)
    if (n_strfun /= 0) then
       do c = 1, n_channel
          call k%sf_chain%set_channel (c, config%sf_channel(c))
       end do
    end if
    call k%sf_chain%link_interactions ()
    call k%sf_chain%exchange_mask ()
    call k%sf_chain%init_evaluators (extended_sf = extended_sf)
  end subroutine kinematics_init_sf_chain

  module subroutine kinematics_init_phs (k, config)
    class(kinematics_t), intent(inout) :: k
    class(phs_config_t), intent(in), target :: config
    k%n_channel = config%get_n_channel ()
    call config%allocate_instance (k%phs)
    call k%phs%init (config)
    k%phs_allocated = .true.
    allocate (k%f (k%n_channel))
    k%f = 0
    k%f_allocated = .true.
  end subroutine kinematics_init_phs

  module subroutine kinematics_evaluate_radiation_kinematics (k, r_in)
    class(kinematics_t), intent(inout) :: k
    real(default), intent(in), dimension(:) :: r_in
    select type (phs => k%phs)
    type is (phs_fks_t)
       if (phs%mode == PHS_MODE_ADDITIONAL_PARTICLE) then
          call phs%generate_radiation_variables &
               (r_in(phs%n_r_born + 1 : phs%n_r_born + 3), &
                threshold = k%threshold)
          call phs%compute_cms_energy ()
       end if
    end select
  end subroutine kinematics_evaluate_radiation_kinematics

  module subroutine kinematics_generate_fsr_in (kin)
    class(kinematics_t), intent(inout) :: kin
    select type (phs => kin%phs)
    type is (phs_fks_t)
       call phs%generate_fsr_in ()
    end select
  end subroutine kinematics_generate_fsr_in

  module subroutine kinematics_compute_xi_ref_momenta (k, reg_data, nlo_type)
    class(kinematics_t), intent(inout) :: k
    type(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: nlo_type
    logical :: use_contributors
    use_contributors = allocated (reg_data%alr_contributors)
    select type (phs => k%phs)
    type is (phs_fks_t)
       if (use_contributors) then
          call phs%compute_xi_ref_momenta (contributors = reg_data%alr_contributors)
       else if (k%threshold) then
          if (.not. is_subtraction_component (k%emitter, nlo_type)) &
               call phs%compute_xi_ref_momenta_threshold ()
       else
          call phs%compute_xi_ref_momenta ()
       end if
    end select
  end subroutine kinematics_compute_xi_ref_momenta

  module subroutine kinematics_compute_selected_channel &
       (k, mci_work, phs_channel, p, success)
    class(kinematics_t), intent(inout) :: k
    type(mci_work_t), intent(in) :: mci_work
    integer, intent(in) :: phs_channel
    type(vector4_t), dimension(:), intent(out) :: p
    logical, intent(out) :: success
    integer :: sf_channel
    k%selected_channel = phs_channel
    sf_channel = k%phs%config%get_sf_channel (phs_channel)
    call k%sf_chain%compute_kinematics (sf_channel, mci_work%get_x_strfun ())
    call k%sf_chain%get_out_momenta (p(1:k%n_in))
    call k%phs%set_incoming_momenta (p(1:k%n_in))
    call k%phs%compute_flux ()
    call k%phs%select_channel (phs_channel)
    call k%phs%evaluate_selected_channel (phs_channel, &
         mci_work%get_x_process ())

    select type (phs => k%phs)
    type is (phs_fks_t)
       if (debug_on)  call msg_debug2 (D_REAL, "phase space is phs_FKS")
       if (phs%q_defined) then
          call phs%get_born_momenta (p)
          if (debug_on) then
             call msg_debug2 (D_REAL, "q is defined")
             call msg_debug2 (D_REAL, "get_born_momenta called")
          end if
          k%phs_factor = phs%get_overall_factor ()
          success = .true.
      else
         k%phs_factor = 0
         success = .false.
      end if
    class default
      if (phs%q_defined) then
         call k%phs%get_outgoing_momenta (p(k%n_in + 1 :))
         k%phs_factor = k%phs%get_overall_factor ()
         success = .true.
      else
         k%phs_factor = 0
         success = .false.
      end if
    end select
  end subroutine kinematics_compute_selected_channel

  module subroutine kinematics_redo_sf_chain (kin, mci_work, phs_channel)
    class(kinematics_t), intent(inout) :: kin
    type(mci_work_t), intent(in) :: mci_work
    integer, intent(in) :: phs_channel
    real(default), dimension(:), allocatable :: x
    integer :: sf_channel, n
    real(default) :: xi, y
    n = size (mci_work%get_x_strfun ())
    if (n > 0) then
       allocate (x(n))
       x = mci_work%get_x_strfun ()
       sf_channel = kin%phs%config%get_sf_channel (phs_channel)
       call kin%sf_chain%compute_kinematics (sf_channel, x)
    end if
  end subroutine kinematics_redo_sf_chain

  module subroutine kinematics_compute_other_channels (k, mci_work, phs_channel)
    class(kinematics_t), intent(inout) :: k
    type(mci_work_t), intent(in) :: mci_work
    integer, intent(in) :: phs_channel
    integer :: c, c_sf
    call k%phs%evaluate_other_channels (phs_channel)
    do c = 1, k%n_channel
       c_sf = k%phs%config%get_sf_channel (c)
       k%f(c) = k%sf_chain%get_f (c_sf) * k%phs%get_f (c)
    end do
  end subroutine kinematics_compute_other_channels

  module subroutine kinematics_get_incoming_momenta (k, p)
    class(kinematics_t), intent(in) :: k
    type(vector4_t), dimension(:), intent(out) :: p
    type(interaction_t), pointer :: int
    integer :: i
    int => k%sf_chain%get_out_int_ptr ()
    do i = 1, k%n_in
       p(i) = int%get_momentum (k%sf_chain%get_out_i (i))
    end do
  end subroutine kinematics_get_incoming_momenta

  module function kinematics_get_boost_to_lab (kin) result (lt)
    type(lorentz_transformation_t) :: lt
    class(kinematics_t), intent(in) :: kin
    lt = kin%phs%get_lorentz_transformation ()
  end function kinematics_get_boost_to_lab

  module function kinematics_get_boost_to_cms (kin) result (lt)
    type(lorentz_transformation_t) :: lt
    class(kinematics_t), intent(in) :: kin
    lt = inverse (kin%phs%get_lorentz_transformation ())
  end function kinematics_get_boost_to_cms

  module subroutine kinematics_recover_mcpar (k, mci_work, phs_channel, p)
    class(kinematics_t), intent(inout) :: k
    type(mci_work_t), intent(inout) :: mci_work
    integer, intent(in) :: phs_channel
    type(vector4_t), dimension(:), intent(in) :: p
    integer :: c, c_sf
    real(default), dimension(:), allocatable :: x_sf, x_phs
    c = phs_channel
    c_sf = k%phs%config%get_sf_channel (c)
    k%selected_channel = c
    call k%sf_chain%recover_kinematics (c_sf)
    call k%phs%set_incoming_momenta (p(1:k%n_in))
    call k%phs%compute_flux ()
    call k%phs%set_outgoing_momenta (p(k%n_in+1:))
    call k%phs%inverse ()
    do c = 1, k%n_channel
       c_sf = k%phs%config%get_sf_channel (c)
       k%f(c) = k%sf_chain%get_f (c_sf) * k%phs%get_f (c)
    end do
    k%phs_factor = k%phs%get_overall_factor ()
    c = phs_channel
    c_sf = k%phs%config%get_sf_channel (c)
    allocate (x_sf (k%sf_chain%config%get_n_bound ()))
    allocate (x_phs (k%phs%config%get_n_par ()))
    call k%phs%select_channel (c)
    call k%sf_chain%get_mcpar (c_sf, x_sf)
    call k%phs%get_mcpar (c, x_phs)
    call mci_work%set_x_strfun (x_sf)
    call mci_work%set_x_process (x_phs)
  end subroutine kinematics_recover_mcpar

  module subroutine kinematics_recover_sfchain (k, channel, p)
    class(kinematics_t), intent(inout) :: k
    integer, intent(in) :: channel
    type(vector4_t), dimension(:), intent(in) :: p
    k%selected_channel = channel
    call k%sf_chain%recover_kinematics (channel)
  end subroutine kinematics_recover_sfchain

  module subroutine kinematics_get_mcpar (k, phs_channel, r)
    class(kinematics_t), intent(in) :: k
    integer, intent(in) :: phs_channel
    real(default), dimension(:), intent(out) :: r
    integer :: sf_channel, n_par_sf, n_par_phs
    sf_channel = k%phs%config%get_sf_channel (phs_channel)
    n_par_phs = k%phs%config%get_n_par ()
    n_par_sf = k%sf_chain%config%get_n_bound ()
    if (n_par_sf > 0) then
       call k%sf_chain%get_mcpar (sf_channel, r(1:n_par_sf))
    end if
    if (n_par_phs > 0) then
       call k%phs%get_mcpar (phs_channel, r(n_par_sf+1:))
    end if
  end subroutine kinematics_get_mcpar

  module subroutine kinematics_evaluate_sf_chain &
       (k, fac_scale, negative_sf, sf_rescale)
    class(kinematics_t), intent(inout) :: k
    real(default), intent(in) :: fac_scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(inout), optional :: sf_rescale
    select case (k%sf_chain%get_status ())
    case (SF_DONE_KINEMATICS)
       call k%sf_chain%evaluate (fac_scale, negative_sf = negative_sf, &
            sf_rescale = sf_rescale)
    end select
  end subroutine kinematics_evaluate_sf_chain

  module subroutine kinematics_return_beam_momenta (k)
    class(kinematics_t), intent(in) :: k
    call k%sf_chain%return_beam_momenta ()
  end subroutine kinematics_return_beam_momenta

  module function kinematics_lab_is_cm (k) result (lab_is_cm)
     logical :: lab_is_cm
     class(kinematics_t), intent(in) :: k
     lab_is_cm = k%phs%config%lab_is_cm
  end function kinematics_lab_is_cm

  module subroutine kinematics_modify_momenta_for_subtraction (k, p_in, p_out)
    class(kinematics_t), intent(inout) :: k
    type(vector4_t), intent(in), dimension(:) :: p_in
    type(vector4_t), intent(out), dimension(:), allocatable :: p_out
    allocate (p_out (size (p_in)))
    if (k%threshold) then
       select type (phs => k%phs)
       type is (phs_fks_t)
          p_out = phs%get_onshell_projected_momenta ()
       end select
    else
       p_out = p_in
    end if
  end subroutine kinematics_modify_momenta_for_subtraction

  module subroutine kinematics_threshold_projection (k, pcm_work, nlo_type)
    class(kinematics_t), intent(inout) :: k
    type(pcm_nlo_workspace_t), intent(inout) :: pcm_work
    integer, intent(in) :: nlo_type
    real(default) :: sqrts, mtop
    type(lorentz_transformation_t) :: L_to_cms
    type(vector4_t), dimension(:), allocatable :: p_tot, p_onshell
    integer :: n_tot
    n_tot = k%phs%get_n_tot ()
    allocate (p_tot (size (pcm_work%real_kinematics%p_born_cms%phs_point(1))))
    select type (phs => k%phs)
    type is (phs_fks_t)
       p_tot = pcm_work%real_kinematics%p_born_cms%phs_point(1)
    class default
       p_tot(1 : k%n_in) = phs%p
       p_tot(k%n_in + 1 : n_tot) = phs%q
    end select
    sqrts = sum (p_tot (1:k%n_in))**1
    mtop = m1s_to_mpole (sqrts)
    L_to_cms = get_boost_for_threshold_projection (p_tot, sqrts, mtop)
    call pcm_work%real_kinematics%p_born_cms%set_momenta (1, p_tot)
    p_onshell = pcm_work%real_kinematics%p_born_onshell%phs_point(1)
    call threshold_projection_born (mtop, L_to_cms, p_tot, p_onshell)
    pcm_work%real_kinematics%p_born_onshell%phs_point(1) = p_onshell
    if (debug2_active (D_THRESHOLD)) then
       print *, 'On-shell projected Born: '
       call vector4_write_set (p_onshell)
    end if
  end subroutine kinematics_threshold_projection

  module subroutine kinematics_evaluate_radiation (k, p_in, p_out, success)
    class(kinematics_t), intent(inout) :: k
    type(vector4_t), intent(in), dimension(:) :: p_in
    type(vector4_t), intent(out), dimension(:), allocatable :: p_out
    logical, intent(out) :: success
    type(vector4_t), dimension(:), allocatable :: p_real
    type(vector4_t), dimension(:), allocatable :: p_born
    real(default) :: xi_max_offshell, xi_offshell, y_offshell, jac_rand_dummy, phi
    select type (phs => k%phs)
    type is (phs_fks_t)
       allocate (p_born (size (p_in)))
       if (k%threshold) then
          p_born = phs%get_onshell_projected_momenta ()
       else
          p_born = p_in
       end if
       if (.not. k%phs%lab_is_cm () .and. .not. k%threshold) then
            p_born = inverse (k%phs%lt_cm_to_lab) * p_born
       end if
       call phs%compute_xi_max (p_born, k%threshold)
       if (k%emitter >= 0) then
          allocate (p_real (size (p_born) + 1))
          allocate (p_out (size (p_born) + 1))
          if (k%emitter <= k%n_in) then
             call phs%generate_isr (k%i_phs, p_real)
          else
             if (k%threshold) then
                jac_rand_dummy = 1._default
                call compute_y_from_emitter (phs%generator%real_kinematics%x_rad (I_Y), &
                     phs%generator%real_kinematics%p_born_cms%get_momenta(1), &
                     k%n_in, k%emitter, .false., phs%generator%y_max, jac_rand_dummy, &
                     y_offshell)
                call phs%compute_xi_max (k%emitter, k%i_phs, y_offshell, &
                     phs%generator%real_kinematics%p_born_cms%get_momenta(1), &
                     xi_max_offshell)
                xi_offshell = xi_max_offshell * phs%generator%real_kinematics%xi_tilde
                phi = phs%generator%real_kinematics%phi
                call phs%generate_fsr (k%emitter, k%i_phs, p_real, &
                     xi_y_phi = [xi_offshell, y_offshell, phi], no_jacobians = .true.)
                call phs%generator%real_kinematics%p_real_cms%set_momenta (k%i_phs, p_real)
                call phs%generate_fsr_threshold (k%emitter, k%i_phs, p_real)
                if (debug2_active (D_SUBTRACTION)) &
                     call generate_fsr_threshold_for_other_emitters (k%emitter, k%i_phs)
             else if (k%i_con > 0) then
                call phs%generate_fsr (k%emitter, k%i_phs, p_real, k%i_con)
             else
                call phs%generate_fsr (k%emitter, k%i_phs, p_real)
             end if
          end if
          success = check_scalar_products (p_real)
          if (debug2_active (D_SUBTRACTION)) then
             call msg_debug2 (D_SUBTRACTION, "Real phase-space: ")
             call vector4_write_set (p_real)
          end if
          p_out = p_real
       else
          allocate (p_out (size (p_in))); p_out = p_in
          success = .true.
       end if
    end select
  contains
    subroutine generate_fsr_threshold_for_other_emitters (emitter, i_phs)
      integer, intent(in) :: emitter, i_phs
      integer :: ii_phs, this_emitter
      select type (phs => k%phs)
      type is (phs_fks_t)
         do ii_phs = 1, size (phs%phs_identifiers)
            this_emitter = phs%phs_identifiers(ii_phs)%emitter
            if (ii_phs /= i_phs .and. this_emitter /= emitter) &
                 call phs%generate_fsr_threshold (this_emitter, i_phs)
         end do
      end select
    end subroutine
  end subroutine kinematics_evaluate_radiation


end submodule kinematics_s

