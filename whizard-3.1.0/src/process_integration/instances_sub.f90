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

submodule (instances) instances_s

  use debug_master, only: debug_on
  use io_units
  use format_utils, only: write_separator
  use constants
  use diagnostics
  use numeric_utils
  use helicities
  use flavors
  use pdg_arrays, only: is_quark, is_charged_lepton, flv_eqv_expr_class

  !!! We should depend less on these modules (move it to pcm_nlo_t e.g.)
  use phs_wood, only: phs_wood_t
  use phs_fks
  use blha_olp_interfaces, only: prc_blha_t
  use blha_config, only: BLHA_AMP_COLOR_C
  use prc_omega, only: prc_omega_t, omega_state_t
  use prc_external, only: prc_external_t, prc_external_state_t
  use prc_threshold, only: prc_threshold_t
  use blha_olp_interfaces, only: blha_result_array_size
  use prc_openloops, only: prc_openloops_t, openloops_state_t
  use prc_recola, only: prc_recola_t
  use blha_olp_interfaces, only: blha_color_c_fill_offdiag, blha_color_c_fill_diag
  use ttv_formfactors, only: m1s_to_mpole

  implicit none

contains

  module subroutine term_instance_write &
       (term, unit, kin, show_eff_state, testflag)
    class(term_instance_t), intent(in) :: term
    integer, intent(in), optional :: unit
    type(kinematics_t), intent(in), optional :: kin
    logical, intent(in), optional :: show_eff_state
    logical, intent(in), optional :: testflag
    real(default) :: fac_scale, ren_scale
    integer :: u
    logical :: state
    u = given_output_unit (unit)
    state = .true.;  if (present (show_eff_state))  state = show_eff_state
    if (term%active) then
       if (associated (term%config)) then
          write (u, "(1x,A,I0,A,I0,A)")  "Term #", term%config%i_term, &
               " (component #", term%config%i_component, ")"
       else
          write (u, "(1x,A)")  "Term [undefined]"
       end if
    else
       write (u, "(1x,A,I0,A)")  "Term #", term%config%i_term, &
            " [inactive]"
    end if
    if (term%checked) then
       write (u, "(3x,A,L1)")      "passed cuts           = ", term%passed
    end if
    if (term%passed) then
       write (u, "(3x,A,ES19.12)")  "overall scale         = ", term%scale
       write (u, "(3x,A,ES19.12)")  "factorization scale   = ", term%get_fac_scale ()
       write (u, "(3x,A,ES19.12)")  "renormalization scale = ", term%get_ren_scale ()
       if (allocated (term%alpha_qcd_forced)) then
          write (u, "(3x,A,ES19.12)")  "alpha(QCD) forced     = ", &
               term%alpha_qcd_forced
       end if
       write (u, "(3x,A,ES19.12)")  "reweighting factor    = ", term%weight
    end if
    !!! This used to be a member of term_instance
    if (present (kin)) then
       call kin%write (u)
    end if
    call write_separator (u)
    write (u, "(1x,A)")  "Amplitude (transition matrix of the &
         &hard interaction):"
    call write_separator (u)
    call term%int_hard%basic_write (u, testflag = testflag)
    if (state .and. term%isolated%has_trace) then
       call write_separator (u)
       write (u, "(1x,A)")  "Evaluators for the hard interaction:"
       call term%isolated%write (u, testflag = testflag)
    end if
    if (state .and. term%connected%has_trace) then
       call write_separator (u)
       write (u, "(1x,A)")  "Evaluators for the connected process:"
       call term%connected%write (u, testflag = testflag)
    end if
  end subroutine term_instance_write

  module subroutine term_instance_final (term)
    class(term_instance_t), intent(inout) :: term
    if (allocated (term%amp)) deallocate (term%amp)
    if (allocated (term%core_state)) deallocate (term%core_state)
    if (allocated (term%ren_scale))  deallocate (term%ren_scale)
    if (allocated (term%fac_scale))  deallocate (term%fac_scale)
    if (allocated (term%es_scale))  deallocate (term%es_scale)
    if (allocated (term%alpha_qcd_forced)) &
       deallocate (term%alpha_qcd_forced)
    if (allocated (term%p_seed)) deallocate(term%p_seed)
    if (allocated (term%p_hard)) deallocate (term%p_hard)
    call term%connected%final ()
    call term%isolated%final ()
    call term%int_hard%final ()
    term%pcm => null ()
    term%pcm_work => null ()
  end subroutine term_instance_final

  module subroutine term_instance_configure &
       (term_instance, process, i, pcm_work, sf_chain, kin)
    class(term_instance_t), intent(out), target :: term_instance
    type(process_t), intent(in), target :: process
    integer, intent(in) :: i
    class(pcm_workspace_t), intent(in), target :: pcm_work
    type(sf_chain_t), intent(in), target :: sf_chain
    type(kinematics_t), intent(inout), target :: kin
    type(process_term_t) :: term
    integer :: i_component
    logical :: requires_extended_sf
    term = process%get_term_ptr (i)
    i_component = term%i_component
    if (i_component /= 0) then
       call term_instance%init &
            (process%get_pcm_ptr (), pcm_work, process%get_nlo_type_component (i_component))
       requires_extended_sf = term_instance%nlo_type == NLO_DGLAP .or. &
              (term_instance%nlo_type == NLO_REAL .and. process%get_i_sub (i) == i)
       call term_instance%setup_dynamics (process, i, kin, &
            real_finite = process%component_is_real_finite (i_component))
       select type (phs => kin%phs)
       type is (phs_fks_t)
          call term_instance%set_emitter (kin)
          call term_instance%setup_fks_kinematics (kin, &
               process%get_var_list_ptr (), &
               process%get_beam_config_ptr ())
       end select
       select type (pcm => term_instance%pcm)
       type is (pcm_nlo_t)
          call kin%set_threshold (pcm%settings%factorization_mode)
       end select
       call term_instance%setup_expressions (process%get_meta (), process%get_config ())
    end if
  end subroutine term_instance_configure

  module subroutine term_instance_init (term_instance, pcm, pcm_work, nlo_type)
    class(term_instance_t), intent(out) :: term_instance
    class(pcm_t), intent(in), target :: pcm
    class(pcm_workspace_t), intent(in), target :: pcm_work
    integer, intent(in) :: nlo_type
    term_instance%pcm => pcm
    term_instance%pcm_work => pcm_work
    term_instance%nlo_type = nlo_type
  end subroutine term_instance_init

  module subroutine term_instance_setup_dynamics &
       (term, process, i_term, kin, real_finite)
    class(term_instance_t), intent(inout), target :: term
    type(process_t), intent(in), target:: process
    integer, intent(in) :: i_term
    type(kinematics_t), intent(in) :: kin
    logical, intent(in), optional :: real_finite
    class(prc_core_t), pointer :: core => null ()
    type(process_beam_config_t) :: beam_config
    type(interaction_t), pointer :: sf_chain_int
    type(interaction_t), pointer :: src_int
    type(quantum_numbers_mask_t), dimension(:), allocatable :: mask_in
    type(state_matrix_t), pointer :: state_matrix
    type(flavor_t), dimension(:), allocatable :: flv_int, flv_src, f_in, f_out
    integer, dimension(:,:), allocatable :: flv_born, flv_real
    type(flavor_t), dimension(:,:), allocatable :: flv_pdf
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn_pdf
    integer :: n_in, n_vir, n_out, n_tot, n_sub
    integer :: n_flv_born, n_flv_real, n_flv_total
    integer :: i, j
    logical :: me_already_squared, keep_fs_flavors
    logical :: decrease_n_tot
    logical :: requires_extended_sf
    me_already_squared = .false.
    keep_fs_flavors = .false.
    term%config => process%get_term_ptr (i_term)
    term%int_hard = term%config%int
    core => process%get_core_term (i_term)
    term%negative_sf = process%get_negative_sf ()
    call core%allocate_workspace (term%core_state)
    select type (core)
    class is (prc_external_t)
       call reduce_interaction (term%int_hard, &
            core%includes_polarization (), .true., .false.)
       me_already_squared = .true.
       allocate (term%amp (term%int_hard%get_n_matrix_elements ()))
    class default
       allocate (term%amp (term%config%n_allowed))
    end select
    if (allocated (term%core_state)) then
       select type (core_state => term%core_state)
       type is (openloops_state_t)
          call core_state%init_threshold (process%get_model_ptr ())
       end select
    end if
    term%amp = cmplx (0, 0, default)
    decrease_n_tot = term%nlo_type == NLO_REAL .and. &
         term%config%i_term_global /= term%config%i_sub
    if (present (real_finite)) then
       if (real_finite) decrease_n_tot = .false.
    end if
    if (decrease_n_tot) then
       allocate (term%p_seed (term%int_hard%get_n_tot () - 1))
    else
       allocate (term%p_seed (term%int_hard%get_n_tot ()))
    end if
    allocate (term%p_hard (term%int_hard%get_n_tot ()))
    sf_chain_int => kin%sf_chain%get_out_int_ptr ()
    n_in = term%int_hard%get_n_in ()
    do j = 1, n_in
       i = kin%sf_chain%get_out_i (j)
       call term%int_hard%set_source_link (j, sf_chain_int, i)
    end do
    call term%isolated%init (kin%sf_chain, term%int_hard)
    allocate (mask_in (n_in))
    mask_in = kin%sf_chain%get_out_mask ()
    select type (phs => kin%phs)
      type is (phs_wood_t)
         if (me_already_squared) then
            call term%isolated%setup_identity_trace &
                 (core, mask_in, .true., .false.)
         else
            call term%isolated%setup_square_trace &
                 (core, mask_in, term%config%col, .false.)
         end if
      type is (phs_fks_t)
         select case (phs%mode)
         case (PHS_MODE_ADDITIONAL_PARTICLE)
            if (me_already_squared) then
               call term%isolated%setup_identity_trace &
                    (core, mask_in, .true., .false.)
            else
               keep_fs_flavors = term%config%data%n_flv > 1
               call term%isolated%setup_square_trace &
                    (core, mask_in, term%config%col, &
                    keep_fs_flavors)
            end if
         case (PHS_MODE_COLLINEAR_REMNANT)
            if (me_already_squared) then
               call term%isolated%setup_identity_trace &
                    (core, mask_in, .true., .false.)
            else
               call term%isolated%setup_square_trace &
                    (core, mask_in, term%config%col, .false.)
            end if
         end select
      class default
         call term%isolated%setup_square_trace &
              (core, mask_in, term%config%col, .false.)
    end select
    if (term%nlo_type == NLO_VIRTUAL .or. (term%nlo_type == NLO_REAL .and. &
         term%config%i_term_global == term%config%i_sub) .or. &
         term%nlo_type == NLO_MISMATCH) then
       n_sub = term%get_n_sub ()
    else if (term%nlo_type == NLO_DGLAP) then
       n_sub = n_beams_rescaled + term%get_n_sub ()
    else
       !!! No integration of real subtraction in interactions yet
       n_sub = 0
    end if
    keep_fs_flavors = keep_fs_flavors .or. me_already_squared
    requires_extended_sf = term%nlo_type == NLO_DGLAP .or. &
         (term%is_subtraction () .and. process%pcm_contains_pdfs ())
    call term%connected%setup_connected_trace (term%isolated, &
         undo_helicities = undo_helicities (core, me_already_squared), &
         keep_fs_flavors = keep_fs_flavors, &
         requires_extended_sf = requires_extended_sf)
    associate (int_eff => term%isolated%int_eff)
      state_matrix => int_eff%get_state_matrix_ptr ()
      n_tot = int_eff%get_n_tot  ()
      flv_int = quantum_numbers_get_flavor &
           (state_matrix%get_quantum_number (1))
      allocate (f_in (n_in))
      f_in = flv_int(1:n_in)
      deallocate (flv_int)
    end associate
    n_in = term%connected%trace%get_n_in ()
    n_vir = term%connected%trace%get_n_vir ()
    n_out = term%connected%trace%get_n_out ()
    allocate (f_out (n_out))
    do j = 1, n_out
       call term%connected%trace%find_source &
            (n_in + n_vir + j, src_int, i)
       if (associated (src_int)) then
          state_matrix => src_int%get_state_matrix_ptr ()
          flv_src = quantum_numbers_get_flavor &
               (state_matrix%get_quantum_number (1))
          f_out(j) = flv_src(i)
          deallocate (flv_src)
       end if
    end do

    beam_config = process%get_beam_config ()
    select type (pcm => term%pcm)
    type is (pcm_nlo_t)
       term%flv_dep_cut_eval = pcm%settings%nlo_correction_type == "EW" &
                               .and. pcm%region_data%alphas_power > 0 &
                               .and. any(is_charged_lepton(f_out%get_pdg()))
    end select

    call term%connected%setup_subevt (term%isolated%sf_chain_eff, &
         beam_config%data%flv, f_in, f_out)
    call term%connected%setup_var_list &
         (process%get_var_list_ptr (), beam_config%data)
    ! Does connected%trace never have any helicity qn?
    call term%init_interaction_qn_index (core, term%connected%trace, n_sub, &
         process%get_model_ptr (), is_polarized = .false.)
    call term%init_interaction_qn_index &
         (core, term%int_hard, n_sub, process%get_model_ptr ())
    call term%init_eqv_expr_classes ()
    if (requires_extended_sf) then
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          n_in = pcm%region_data%get_n_in ()
          flv_born = pcm%region_data%get_flv_states_born ()
          flv_real = pcm%region_data%get_flv_states_real ()
          n_flv_born = pcm%region_data%get_n_flv_born ()
          n_flv_real = pcm%region_data%get_n_flv_real ()
          n_flv_total = n_flv_born + n_flv_real
          allocate (flv_pdf(n_in, n_flv_total), &
               qn_pdf(n_in, n_flv_total))
          call flv_pdf(:, :n_flv_born)%init (flv_born(:n_in, :))
          call flv_pdf(:, n_flv_born + 1:n_flv_total)%init (flv_real(:n_in, :))
          call qn_pdf%init (flv_pdf)
          call sf_chain_int%init_qn_index (qn_pdf, n_flv_born, n_flv_real)
       end select
    end if
  contains

   function undo_helicities (core, me_squared) result (val)
     logical :: val
     class(prc_core_t), intent(in) :: core
     logical, intent(in) :: me_squared
     select type (core)
     class is (prc_external_t)
        val = me_squared .and. .not. core%includes_polarization ()
     class default
        val = .false.
     end select
   end function undo_helicities

   subroutine reduce_interaction (int, polarized_beams, keep_fs_flavors, &
      keep_colors)
     type(interaction_t), intent(inout) :: int
     logical, intent(in) :: polarized_beams
     logical, intent(in) :: keep_fs_flavors, keep_colors
     type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask
     logical, dimension(:), allocatable :: mask_f, mask_c, mask_h
     integer :: n_tot, n_in
     n_in = int%get_n_in (); n_tot = int%get_n_tot ()
     allocate (qn_mask (n_tot))
     allocate (mask_f (n_tot), mask_c (n_tot), mask_h (n_tot))
     mask_c = .not. keep_colors
     mask_f (1 : n_in) = .false.
     if (keep_fs_flavors) then
        mask_f (n_in + 1 : ) = .false.
     else
        mask_f (n_in + 1 : ) = .true.
     end if
     if (polarized_beams) then
        mask_h (1 : n_in) = .false.
     else
        mask_h (1 : n_in) = .true.
     end if
     mask_h (n_in + 1 : ) = .true.
     call qn_mask%init (mask_f, mask_c, mask_h)
     call int%reduce_state_matrix (qn_mask, keep_order = .true.)
   end subroutine reduce_interaction
end subroutine term_instance_setup_dynamics

  module subroutine setup_interaction_qn_index &
       (int, data, qn_config, n_sub, is_polarized)
    class(interaction_t), intent(inout) :: int
    class(process_constants_t), intent(in) :: data
    type(quantum_numbers_t), dimension(:, :), intent(in) :: qn_config
    integer, intent(in) :: n_sub
    logical, intent(in) :: is_polarized
    integer :: i
    type(quantum_numbers_t), dimension(:, :), allocatable :: qn_hel
    if (is_polarized) then
       call setup_interaction_qn_hel (int, data, qn_hel)
       call int%init_qn_index (qn_config, n_sub, qn_hel)
       call int%set_qn_index_helicity_flip (.true.)
    else
       call int%init_qn_index (qn_config, n_sub)
    end if
  end subroutine setup_interaction_qn_index

  module subroutine setup_interaction_qn_hel (int, data, qn_hel)
    class(interaction_t), intent(in) :: int
    class(process_constants_t), intent(in) :: data
    type(quantum_numbers_t), dimension(:, :), allocatable, intent(out) :: &
         qn_hel
    type(helicity_t), dimension(:), allocatable :: hel
    integer, dimension(:), allocatable :: index_table
    integer, dimension(:, :), allocatable :: hel_state
    integer :: i, j, n_hel_unique
    associate (n_in => int%get_n_in (), n_tot => int%get_n_tot ())
      allocate (hel_state (n_tot, data%get_n_hel ()), &
           source = data%hel_state)
      allocate (index_table (data%get_n_hel ()), &
           source = 0)
      forall (j=1:data%get_n_hel (), i=n_in+1:n_tot) hel_state(i, j) = 0
      n_hel_unique = 0
      HELICITY: do i = 1, data%get_n_hel ()
         do j = 1, data%get_n_hel ()
            if (index_table (j) == 0) then
               index_table(j) = i; n_hel_unique = n_hel_unique + 1
               cycle HELICITY
            else if (all (hel_state(:, i) == &
                 hel_state(:, index_table(j)))) then
               cycle HELICITY
            end if
         end do
      end do HELICITY
      allocate (qn_hel (n_tot, n_hel_unique))
      allocate (hel (n_tot))
      do j = 1, n_hel_unique
         call hel%init (hel_state(:, index_table(j)))
         call qn_hel(:, j)%init (hel)
      end do
    end associate
  end subroutine setup_interaction_qn_hel

  module subroutine term_instance_init_eqv_expr_classes (term)
    class(term_instance_t), intent(inout), target :: term
    type(interaction_t), pointer :: src_int
    type(state_matrix_t), pointer :: state_matrix
    type(flavor_t), dimension(:), allocatable :: flv_src
    logical, dimension(:,:,:), allocatable :: eqv_expr_class
    logical, dimension (:), allocatable :: evaluated
    integer :: n_in, n_vir, n_out
    integer :: k, j, i
    n_in = term%connected%trace%get_n_in ()
    n_vir = term%connected%trace%get_n_vir ()
    n_out = term%connected%trace%get_n_out ()
    allocate (eqv_expr_class (3, n_out, &
         term%connected%trace%get_qn_index_n_flv ()))
    do k = 1, term%connected%trace%get_qn_index_n_flv ()
       do j = 1, n_out
          call term%connected%trace%find_source &
               (n_in + n_vir + j, src_int, i)
          if (associated (src_int)) then
             state_matrix => src_int%get_state_matrix_ptr ()
             flv_src = quantum_numbers_get_flavor &
                  (state_matrix%get_quantum_number (k))
             eqv_expr_class (:, j, k) = flv_eqv_expr_class (flv_src(i)%get_pdg())
             deallocate (flv_src)
          end if
       end do
    end do
    if (term%flv_dep_cut_eval) then
       allocate (evaluated (term%connected%trace%get_qn_index_n_flv ()))
       evaluated = .false.
       allocate (term%i_flv_to_i_flv_rep (term%connected%trace%get_qn_index_n_flv ()))
       do i = 1, term%connected%trace%get_qn_index_n_flv ()
          if (.not. evaluated (i)) then
             do k = i, term%connected%trace%get_qn_index_n_flv ()
                if (same_eqv_expr_class(eqv_expr_class (:,:,i), eqv_expr_class (:,:,k))) then
                   term%i_flv_to_i_flv_rep (k) = i
                   evaluated (k) = .true.
                end if
             end do
          end if
       end do
    end if

  contains

     function same_eqv_expr_class (flv_mask1, flv_mask2) result (same)
        logical, dimension (:,:), intent(in) :: flv_mask1, flv_mask2
        logical :: same
        integer :: l
        same = .true.
        do l = 1, size (flv_mask1, dim = 2)
           same = same .and. all (flv_mask1(:,l) .eqv. flv_mask2(:,l))
        end do
     end function same_eqv_expr_class
  end subroutine term_instance_init_eqv_expr_classes

  module subroutine term_instance_init_interaction_qn_index (term, core, &
       int, n_sub, model, is_polarized)
    class(term_instance_t), intent(inout), target :: term
    class(prc_core_t), intent(in) :: core
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: n_sub
    class(model_data_t), intent(in) :: model
    logical, intent(in), optional :: is_polarized
    logical :: polarized
    type(quantum_numbers_t), dimension(:, :), allocatable :: qn_config
    integer, dimension(:,:), allocatable :: flv_born
    type(flavor_t), dimension(:), allocatable :: flv
    integer :: i
    select type (core)
    class is (prc_external_t)
       if (present (is_polarized)) then
          polarized = is_polarized
       else
          polarized = core%includes_polarization ()
       end if
       select type (pcm_work => term%pcm_work)
       type is (pcm_nlo_workspace_t)
          associate (is_born => .not. (term%nlo_type == NLO_REAL .and. &
                  .not. term%is_subtraction ()))
             select type (pcm => term%pcm)
             type is (pcm_nlo_t)
                qn_config = pcm%get_qn (is_born)
             end select
             call setup_interaction_qn_index (int, term%config%data, &
                  qn_config, n_sub, polarized)
          end associate
       class default
          call term%config%data%get_flv_state (flv_born)
          allocate (flv (size (flv_born, dim = 1)))
          allocate (qn_config (size (flv_born, dim = 1), size (flv_born, dim = 2)))
          do i = 1, core%data%n_flv
             call flv%init (flv_born(:,i), model)
             call qn_config(:, i)%init (flv)
          end do
          call setup_interaction_qn_index (int, term%config%data, &
               qn_config, n_sub, polarized)
       end select
    class default
       call int%init_qn_index ()
    end select
  end subroutine term_instance_init_interaction_qn_index

  module subroutine term_instance_setup_fks_kinematics &
       (term, kin, var_list, beam_config)
    class(term_instance_t), intent(inout), target :: term
    type(kinematics_t), intent(inout) :: kin
    type(var_list_t), intent(in) :: var_list
    type(process_beam_config_t), intent(in) :: beam_config
    integer :: mode
    logical :: singular_jacobian
    if (.not. (term%nlo_type == NLO_REAL .or. term%nlo_type == NLO_DGLAP .or. &
       term%nlo_type == NLO_MISMATCH)) return
    singular_jacobian = var_list%get_lval &
         (var_str ("?powheg_use_singular_jacobian"))
    if (term%nlo_type == NLO_REAL) then
       mode = check_generator_mode (GEN_REAL_PHASE_SPACE)
    else if (term%nlo_type == NLO_MISMATCH) then
       mode = check_generator_mode (GEN_SOFT_MISMATCH)
    else
       mode = PHS_MODE_UNDEFINED
    end if
    select type (phs => kin%phs)
    type is (phs_fks_t)
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          select type (pcm_work => term%pcm_work)
          type is (pcm_nlo_workspace_t)
             call pcm%setup_phs_generator (pcm_work, &
                  phs%generator, phs%config%sqrts, mode, singular_jacobian)
             if (beam_config%has_structure_function ()) then
                pcm_work%isr_kinematics%isr_mode = SQRTS_VAR
             else
                pcm_work%isr_kinematics%isr_mode = SQRTS_FIXED
             end if
             if (debug_on)  call msg_debug &
                  (D_PHASESPACE, "isr_mode: ", pcm_work%isr_kinematics%isr_mode)
          end select
       end select
    class default
       call msg_fatal ("Phase space should be an FKS phase space!")
    end select
  contains
    function check_generator_mode (gen_mode_default) result (gen_mode)
       integer :: gen_mode
       integer, intent(in) :: gen_mode_default
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          associate (settings => pcm%settings)
             if (settings%test_coll_limit .and. settings%test_anti_coll_limit) &
                call msg_fatal ("You cannot check the collinear and anti-collinear limit "&
                     &"at the same time!")
             if (settings%test_soft_limit .and. .not. settings%test_coll_limit &
                  .and. .not. settings%test_anti_coll_limit) then
                gen_mode = GEN_SOFT_LIMIT_TEST
             else if (.not. settings%test_soft_limit .and. settings%test_coll_limit) then
                gen_mode = GEN_COLL_LIMIT_TEST
             else if (.not. settings%test_soft_limit .and. settings%test_anti_coll_limit) then
                gen_mode = GEN_ANTI_COLL_LIMIT_TEST
             else if (settings%test_soft_limit .and. settings%test_coll_limit) then
                gen_mode = GEN_SOFT_COLL_LIMIT_TEST
             else if (settings%test_soft_limit .and. settings%test_anti_coll_limit) then
                gen_mode = GEN_SOFT_ANTI_COLL_LIMIT_TEST
             else
                gen_mode = gen_mode_default
             end if
          end associate
       end select
    end function check_generator_mode
  end subroutine term_instance_setup_fks_kinematics

  module subroutine term_instance_compute_seed_kinematics &
       (term, kin, mci_work, phs_channel, success)
    class(term_instance_t), intent(inout), target :: term
    type(kinematics_t), intent(inout) :: kin
    type(mci_work_t), intent(in) :: mci_work
    integer, intent(in) :: phs_channel
    logical, intent(out) :: success
    call kin%compute_selected_channel &
         (mci_work, phs_channel, term%p_seed, success)
  end subroutine term_instance_compute_seed_kinematics

  module subroutine term_instance_evaluate_projections (term, kin)
    class(term_instance_t), intent(inout) :: term
    type(kinematics_t), intent(inout) :: kin
    if (kin%threshold .and. term%nlo_type > BORN) then
       if (debug2_active (D_THRESHOLD)) &
            print *, 'Evaluate on-shell projection: ', &
            char (component_status (term%nlo_type))
       select type (pcm_work => term%pcm_work)
       type is (pcm_nlo_workspace_t)
          call kin%threshold_projection (pcm_work, term%nlo_type)
       end select
    end if
  end subroutine term_instance_evaluate_projections

  module subroutine term_instance_compute_hard_kinematics &
       (term, kin, recover, skip_term, success)
    class(term_instance_t), intent(inout) :: term
    type(kinematics_t), intent(inout) :: kin
    integer, intent(in), optional :: skip_term
    logical, intent(in), optional :: recover
    logical, intent(out) :: success
    type(vector4_t), dimension(:), allocatable :: p
    if (allocated (term%core_state)) &
       call term%core_state%reset_new_kinematics ()
    if (present (skip_term)) then
       if (term%config%i_term_global == skip_term) return
    end if

    if (present (recover)) then
       if (recover) return
    end if
    if (term%nlo_type == NLO_REAL .and. kin%emitter >= 0) then
       call kin%evaluate_radiation (term%p_seed, p, success)
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          if (pcm%dalitz_plot%active) then
             if (kin%emitter > kin%n_in) then
                if (p(kin%emitter)**2 > tiny_07) &
                     call pcm%register_dalitz_plot (kin%emitter, p)
             end if
          end if
       end select
    else if (is_subtraction_component (kin%emitter, term%nlo_type)) then
       call kin%modify_momenta_for_subtraction (term%p_seed, p)
       success = .true.
    else
       allocate (p (size (term%p_seed))); p = term%p_seed
       success = .true.
    end if
    call term%int_hard%set_momenta (p)
    if (debug_on) then
       call msg_debug2 (D_REAL, "inside compute_hard_kinematics")
       if (debug2_active (D_REAL))  call vector4_write_set (p)
    end if
  end subroutine term_instance_compute_hard_kinematics

  module subroutine term_instance_recover_seed_kinematics &
       (term, kin, p_seed_ref)
    class(term_instance_t), intent(inout) :: term
    type(kinematics_t), intent(in) :: kin
    integer :: n_in
    type(vector4_t), dimension(:), intent(in), optional :: p_seed_ref
    n_in = kin%n_in
    call kin%get_incoming_momenta (term%p_seed(1:n_in))
    associate (int_eff => term%isolated%int_eff)
       call int_eff%set_momenta (term%p_seed(1:n_in), outgoing = .false.)
       if (present (p_seed_ref)) then
          term%p_seed(n_in + 1 : ) = p_seed_ref
       else
          term%p_seed(n_in + 1 : ) = int_eff%get_momenta (outgoing = .true.)
       end if
    end associate
    call term%isolated%receive_kinematics ()
  end subroutine term_instance_recover_seed_kinematics

  module subroutine term_instance_apply_real_partition (term, kin)
    class(term_instance_t), intent(inout) :: term
    type(kinematics_t), intent(in) :: kin
    real(default) :: f, sqme
    integer :: i_component
    integer :: i_amp, n_amps, qn_index
    logical :: is_subtraction
    i_component = term%config%i_component
    select type (pcm => term%pcm)
    type is (pcm_nlo_t)
       if (pcm%component_selected (i_component) .and. &
           pcm%nlo_type (i_component) == NLO_REAL) then
          is_subtraction = pcm%component_type (i_component) == &
               COMP_REAL_SING .and. kin%emitter < 0
          if (is_subtraction) return
          select case (pcm%component_type (i_component))
          case (COMP_REAL_FIN)
             call term%connected%trace%set_duplicate_flv_zero()
          end select
          f = pcm%real_partition%get_f (term%p_hard)
          n_amps = term%connected%trace%get_n_matrix_elements ()
          do i_amp = 1, n_amps
             qn_index = term%connected%trace%get_qn_index (i_amp, i_sub = 0)
             if (term%passed_array(i_amp) .or. .not. term%passed) then
                sqme = real (term%connected%trace%get_matrix_element (qn_index))
             else
                sqme = zero
             end if
             if (debug_on)  call msg_debug2 &
                  (D_PROCESS_INTEGRATION, "term_instance_apply_real_partition")
             select case (pcm%component_type (i_component))
             case (COMP_REAL_FIN)
                if (debug_on)  call msg_debug2 &
                     (D_PROCESS_INTEGRATION, "Real finite")
                sqme = sqme * (one - f)
             case (COMP_REAL_SING)
                if (debug_on)  call msg_debug2 &
                     (D_PROCESS_INTEGRATION, "Real singular")
                sqme = sqme * f
             end select
             if (debug_on)  call msg_debug2 &
                  (D_PROCESS_INTEGRATION, "apply_damping: sqme", sqme)
             call term%connected%trace%set_matrix_element &
                  (qn_index, cmplx (sqme, zero, default))
          end do
       end if
    end select
  end subroutine term_instance_apply_real_partition

  pure module function term_instance_get_p_hard (term_instance) result (p_hard)
    type(vector4_t), dimension(:), allocatable :: p_hard
    class(term_instance_t), intent(in) :: term_instance
    allocate (p_hard (size (term_instance%p_hard)))
    p_hard = term_instance%p_hard
  end function term_instance_get_p_hard

  module subroutine term_instance_set_emitter (term, kin)
    class(term_instance_t), intent(inout) :: term
    type(kinematics_t), intent(inout) :: kin
    integer :: i_phs
    logical :: set_emitter
    select type (pcm => term%pcm)
    type is (pcm_nlo_t)
       select type (phs => kin%phs)
       type is (phs_fks_t)
          !!! Without resonances, i_alr = i_phs
          i_phs = term%config%i_term
          kin%i_phs = i_phs
          set_emitter = i_phs <= pcm%region_data%n_phs .and. &
               term%nlo_type == NLO_REAL
          if (set_emitter) then
             kin%emitter = phs%phs_identifiers(i_phs)%emitter
             select type (pcm => term%pcm)
             type is (pcm_nlo_t)
                if (allocated (pcm%region_data%i_phs_to_i_con)) &
                   kin%i_con = pcm%region_data%i_phs_to_i_con (i_phs)
             end select
          end if
       end select
    end select
  end subroutine term_instance_set_emitter

  module subroutine term_instance_setup_expressions (term, meta, config)
    class(term_instance_t), intent(inout), target :: term
    type(process_metadata_t), intent(in), target :: meta
    type(process_config_data_t), intent(in) :: config
    if (allocated (config%ef_cuts)) &
         call term%connected%setup_cuts (config%ef_cuts)
    if (allocated (config%ef_scale)) &
         call term%connected%setup_scale (config%ef_scale)
    if (allocated (config%ef_fac_scale)) &
         call term%connected%setup_fac_scale (config%ef_fac_scale)
    if (allocated (config%ef_ren_scale)) &
         call term%connected%setup_ren_scale (config%ef_ren_scale)
    if (allocated (config%ef_weight)) &
         call term%connected%setup_weight (config%ef_weight)
  end subroutine term_instance_setup_expressions

  module subroutine term_instance_setup_event_data (term, kin, core, model)
    class(term_instance_t), intent(inout), target :: term
    type(kinematics_t), intent(in) :: kin
    class(prc_core_t), intent(in) :: core
    class(model_data_t), intent(in), target :: model
    integer :: n_in
    logical :: mask_color
    type(quantum_numbers_mask_t), dimension(:), allocatable :: mask_in
    n_in = term%int_hard%get_n_in ()
    allocate (mask_in (n_in))
    mask_in = kin%sf_chain%get_out_mask ()
    call setup_isolated (term%isolated, core, model, mask_in, term%config%col)
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       mask_color = pcm_work%is_fixed_order_nlo_events ()
    class default
       mask_color = .false.
    end select
    call setup_connected (term%connected, term%isolated, core, &
         term%nlo_type, mask_color)
  contains
    subroutine setup_isolated (isolated, core, model, mask, color)
      type(isolated_state_t), intent(inout), target :: isolated
      class(prc_core_t), intent(in) :: core
      class(model_data_t), intent(in), target :: model
      type(quantum_numbers_mask_t), intent(in), dimension(:) :: mask
      integer, intent(in), dimension(:) :: color
      select type (core)
      class is (prc_blha_t)
         call isolated%matrix%init_identity(isolated%int_eff)
         isolated%has_matrix = .true.
      class default
         call isolated%setup_square_matrix (core, model, mask, color)
      end select
      !!! TODO (PS-09-10-20) We should not square the flows
      !!! if they come from BLHA either
      call isolated%setup_square_flows (core, model, mask)
    end subroutine setup_isolated

    subroutine setup_connected (connected, isolated, core, nlo_type, mask_color)
      type(connected_state_t), intent(inout), target :: connected
      type(isolated_state_t), intent(in), target :: isolated
      class(prc_core_t), intent(in) :: core
      integer, intent(in) :: nlo_type
      logical, intent(in) :: mask_color
      type(quantum_numbers_mask_t), dimension(:), allocatable :: mask
      call connected%setup_connected_matrix (isolated)
      if (term%nlo_type == NLO_VIRTUAL .or. (term%nlo_type == NLO_REAL &
           .and. term%config%i_term_global == term%config%i_sub) &
           .or. term%nlo_type == NLO_DGLAP) then
         !!! We do not care about the subtraction matrix elements in
         !!! connected%matrix, because all entries there are supposed
         !!! to be squared. To be able to match with flavor quantum numbers,
         !!! we remove the subtraction quantum entries from the state matrix.
         allocate (mask (connected%matrix%get_n_tot()))
         call mask%set_sub (1)
         call connected%matrix%reduce_state_matrix (mask, keep_order = .true.)
      end if
      call term%init_interaction_qn_index (core, connected%matrix, 0, model, &
           is_polarized = .false.)
      select type (core)
      class is (prc_blha_t)
         call connected%setup_connected_flows &
              (isolated, mask_color = mask_color)
      class default
         call connected%setup_connected_flows (isolated)
      end select
      call connected%setup_state_flv (isolated%get_n_out ())
    end subroutine setup_connected
  end subroutine term_instance_setup_event_data

  module subroutine term_instance_evaluate_color_correlations (term, core)
    class(term_instance_t), intent(inout) :: term
    class(prc_core_t), intent(inout) :: core
    integer :: i_flv_born
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          if (debug_on) call msg_debug2 (D_SUBTRACTION, &
               "term_instance_evaluate_color_correlations: " // &
               "use_internal_color_correlations:", &
               pcm%settings%use_internal_color_correlations)
          if (debug_on) call msg_debug2 (D_SUBTRACTION, "fac_scale", term%get_fac_scale ())
          do i_flv_born = 1, pcm%region_data%n_flv_born
             select case (term%nlo_type)
             case (NLO_REAL)
                call transfer_me_array_to_bij (pcm, i_flv_born, &
                     pcm_work%real_sub%sqme_born (i_flv_born), &
                     pcm_work%real_sub%sqme_born_color_c (:, :, i_flv_born))
             case (NLO_MISMATCH)
                call transfer_me_array_to_bij (pcm, i_flv_born, &
                     pcm_work%soft_mismatch%sqme_born (i_flv_born), &
                     pcm_work%soft_mismatch%sqme_born_color_c (:, :, i_flv_born))
             case (NLO_VIRTUAL)
                !!! This is just a copy of the above with a different offset and can for sure be unified
                call transfer_me_array_to_bij (pcm, i_flv_born, &
                     -one, pcm_work%virtual%sqme_color_c (:, :, i_flv_born))
              case (NLO_DGLAP)
                call transfer_me_array_to_bij (pcm, i_flv_born, &
                     pcm_work%dglap_remnant%sqme_born (i_flv_born), &
                     pcm_work%dglap_remnant%sqme_color_c_extra (:, :, i_flv_born))
             end select
          end do
       end select
    end select
  contains
    function get_trivial_cf_factors (n_tot, flv, factorization_mode) result (beta_ij)
      integer, intent(in) :: n_tot, factorization_mode
      integer, intent(in), dimension(:) :: flv
      real(default), dimension(n_tot, n_tot) :: beta_ij
      if (factorization_mode == NO_FACTORIZATION) then
         beta_ij = get_trivial_cf_factors_default (n_tot, flv)
      else
         beta_ij = get_trivial_cf_factors_threshold (n_tot, flv)
      end if
    end function get_trivial_cf_factors

    function get_trivial_cf_factors_default (n_tot, flv) result (beta_ij)
      integer, intent(in) :: n_tot
      integer, intent(in), dimension(:) :: flv
      real(default), dimension(n_tot, n_tot) :: beta_ij
      integer :: i, j
      beta_ij = zero
      if (count (is_quark (flv)) == 2) then
         do i = 1, n_tot
            do j = 1, n_tot
               if (is_quark(flv(i)) .and. is_quark(flv(j))) then
                  if (i == j) then
                     beta_ij(i,j)= -cf
                  else
                     beta_ij(i,j) = cf
                  end if
               end if
            end do
         end do
      end if
    end function get_trivial_cf_factors_default

    function get_trivial_cf_factors_threshold (n_tot, flv) result (beta_ij)
      integer, intent(in) :: n_tot
      integer, intent(in), dimension(:) :: flv
      real(default), dimension(n_tot, n_tot) :: beta_ij
      integer :: i
      beta_ij = zero
      do i = 1, 4
         beta_ij(i,i) = -cf
      end do
      beta_ij(1,2) = cf; beta_ij(2,1) = cf
      beta_ij(3,4) = cf; beta_ij(4,3) = cf
    end function get_trivial_cf_factors_threshold

    subroutine transfer_me_array_to_bij (pcm, i_flv, &
         sqme_born, sqme_color_c)
      type(pcm_nlo_t), intent(in) :: pcm
      integer, intent(in) :: i_flv
      real(default), intent(in) :: sqme_born
      real(default), dimension(:,:), intent(inout) :: sqme_color_c
      logical :: special_case_interferences
      integer :: i_color_c, i_sub, n_offset, i_qn
      real(default), dimension(:), allocatable :: sqme
      real(default) :: sqme_born_c
      if (debug_on) call msg_debug2 (D_PROCESS_INTEGRATION, "transfer_me_array_to_bij")
      if (pcm%settings%use_internal_color_correlations) then
         !!! A negative value for sqme_born indicates that the Born matrix
         !!! element is multiplied at a different place, e.g. in the case
         !!! of the virtual component
         sqme_color_c = get_trivial_cf_factors &
              (pcm%region_data%get_n_legs_born (), &
              pcm%region_data%get_flv_states_born (i_flv), &
              pcm%settings%factorization_mode)
         if (sqme_born > zero) then
            sqme_color_c = sqme_born * sqme_color_c
         else if (sqme_born == zero) then
            sqme_color_c = zero
         end if
      else
         special_case_interferences = &
              pcm%region_data%nlo_correction_type == "EW"
         n_offset = 0
         if (term%nlo_type == NLO_VIRTUAL) then
            n_offset = 1
         else if (pcm%has_pdfs .and. (term%is_subtraction () &
                .or. term%nlo_type == NLO_DGLAP)) then
            n_offset = n_beams_rescaled
         end if
         allocate (sqme (term%get_n_sub_color ()), source = zero)
         do i_sub = 1, term%get_n_sub_color ()
            i_qn = term%connected%trace%get_qn_index (i_flv, i_sub = i_sub + n_offset)
            if (term%passed_array(i_flv) .or. .not. term%passed) then
               sqme(i_sub) = real(term%connected%trace%get_matrix_element (i_qn), default)
            else
               sqme(i_sub) = zero
            end if
         end do
         call blha_color_c_fill_offdiag (pcm%region_data%n_legs_born, &
              sqme, sqme_color_c)
         i_qn = term%connected%trace%get_qn_index (i_flv, i_sub = 0)
         if (term%passed_array(i_flv) .or. .not. term%passed) then
            sqme_born_c = real(term%connected%trace%get_matrix_element (i_qn), default)
         else
            sqme_born_c = zero
         end if
         call blha_color_c_fill_diag (sqme_born_c, &
              pcm%region_data%get_flv_states_born (i_flv), &
              sqme_color_c, special_case_interferences)
      end if
    end subroutine transfer_me_array_to_bij
  end subroutine term_instance_evaluate_color_correlations

  module subroutine term_instance_evaluate_charge_correlations (term, core)
    class(term_instance_t), intent(inout) :: term
    class(prc_core_t), intent(inout) :: core
    integer :: i_flv_born
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          do i_flv_born = 1, pcm%region_data%n_flv_born
             select case (term%nlo_type)
             case (NLO_REAL)
                call transfer_me_array_to_bij (pcm, i_flv_born, &
                     pcm_work%real_sub%sqme_born (i_flv_born), &
                     pcm_work%real_sub%sqme_born_charge_c (:, :, i_flv_born))
             case (NLO_MISMATCH)
                call transfer_me_array_to_bij (pcm, i_flv_born, &
                     pcm_work%soft_mismatch%sqme_born (i_flv_born), &
                     pcm_work%soft_mismatch%sqme_born_charge_c (:, :, i_flv_born))
             case (NLO_VIRTUAL)
                call transfer_me_array_to_bij (pcm, i_flv_born, &
                     one, pcm_work%virtual%sqme_charge_c (:, :, i_flv_born))
             end select
          end do
       end select
    end select
  contains
    subroutine transfer_me_array_to_bij (pcm, i_flv, sqme_born, sqme_charge_c)
      type(pcm_nlo_t), intent(in) :: pcm
      integer, intent(in) :: i_flv
      real(default), intent(in) :: sqme_born
      real(default), dimension(:,:), intent(inout) :: sqme_charge_c
      integer :: n_legs_born, i, j
      real(default), dimension(:), allocatable :: sigma
      real(default), dimension(:), allocatable :: Q
      n_legs_born = pcm%region_data%n_legs_born
      associate (flv_born => pcm%region_data%flv_born(i_flv))
         allocate (sigma (n_legs_born), Q (size (flv_born%charge)))
         Q = flv_born%charge
         sigma(1:flv_born%n_in) = -one
         sigma(flv_born%n_in + 1: ) = one
      end associate
      do i = 1, n_legs_born
         do j = 1, n_legs_born
            sqme_charge_c(i, j) = sigma(i) * sigma(j) * Q(i) * Q(j) * (-one)
         end do
      end do
      sqme_charge_c = sqme_charge_c * sqme_born
    end subroutine transfer_me_array_to_bij
  end subroutine term_instance_evaluate_charge_correlations

  module subroutine term_instance_evaluate_spin_correlations (term, core)
    class(term_instance_t), intent(inout) :: term
    class(prc_core_t), intent(inout) :: core
    integer :: i_flv, i_sub, i_emitter, emitter, i_qn
    integer :: n_flv, n_sub_color, n_sub_spin, n_offset,i,j
    real(default), dimension(1:3, 1:3) :: sqme_spin_c
    real(default), dimension(:), allocatable :: sqme_spin_c_all
    real(default), dimension(:), allocatable :: sqme_spin_c_arr
    if (debug_on) call msg_debug2 (D_PROCESS_INTEGRATION, &
         "term_instance_evaluate_spin_correlations")
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       if (pcm_work%real_sub%requires_spin_correlations () &
            .and. term%nlo_type == NLO_REAL) then
          select type (core)
          type is (prc_openloops_t)
             select type (pcm => term%pcm)
             type is (pcm_nlo_t)
                n_flv = term%connected%trace%get_qn_index_n_flv ()
                n_sub_color = term%get_n_sub_color ()
                n_sub_spin = term%get_n_sub_spin ()
                n_offset = 0; if(pcm%has_pdfs) n_offset = n_beams_rescaled
                allocate (sqme_spin_c_arr(6))
                do i_flv = 1, n_flv
                   allocate (sqme_spin_c_all(n_sub_spin))
                   do i_sub = 1, n_sub_spin
                      i_qn = term%connected%trace%get_qn_index (i_flv, &
                             i_sub = i_sub + n_offset + n_sub_color)
                      if (term%passed_array(i_flv) .or. .not. term%passed) then
                         sqme_spin_c_all(i_sub) = &
                            real(term%connected%trace%get_matrix_element (i_qn), default)
                      else
                         sqme_spin_c_all(i_sub) = zero
                      end if
                   end do
                   do i_emitter = 1, pcm%region_data%n_emitters
                      emitter = pcm%region_data%emitters(i_emitter)
                      if (emitter > 0) then
                         call split_array (sqme_spin_c_all, sqme_spin_c_arr)
                         do j = 1, size (sqme_spin_c, dim=2)
                            do i = j, size (sqme_spin_c, dim=1)
                               !!! Restoring the symmetric matrix packed into a 1-dim array
                               !!! c.f. [[prc_openloops_compute_sqme_spin_c]]
                               sqme_spin_c(i,j) = sqme_spin_c_arr(j + i * (i - 1) / 2)
                               if (i /= j) sqme_spin_c(j,i) = sqme_spin_c(i,j)
                            end do
                         end do
                         pcm_work%real_sub%sqme_born_spin_c(:,:,emitter,i_flv) = sqme_spin_c
                      end if
                   end do
                   deallocate (sqme_spin_c_all)
                end do
             end select
          class default
             call msg_fatal &
                  ("Spin correlations so far only supported by OpenLoops.")
          end select
       end if
    end select
  end subroutine term_instance_evaluate_spin_correlations

  module subroutine term_instance_apply_fks &
       (term, kin, alpha_s_sub, alpha_qed_sub)
    class(term_instance_t), intent(inout) :: term
    class(kinematics_t), intent(inout) :: kin
    real(default), intent(in) :: alpha_s_sub, alpha_qed_sub
    real(default), dimension(:), allocatable :: sqme
    integer :: i, i_phs, emitter, i_qn
    logical :: is_subtraction
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          if (term%connected%has_matrix) then
             allocate (sqme (pcm%get_n_alr ()))
          else
             allocate (sqme (1))
          end if
          sqme = zero
          select type (phs => kin%phs)
          type is (phs_fks_t)
             if (pcm%has_pdfs .and. &
                  pcm%settings%use_internal_color_correlations) then
                call msg_fatal ("Color correlations for proton processes " // &
                     "so far only supported by OpenLoops.")
             end if
             call pcm_work%set_real_and_isr_kinematics &
                  (phs%phs_identifiers, kin%phs%get_sqrts ())
             if (kin%emitter < 0) then
                call pcm_work%set_subtraction_event ()
                do i_phs = 1, pcm%region_data%n_phs
                   emitter = phs%phs_identifiers(i_phs)%emitter
                   call pcm_work%real_sub%compute (emitter, &
                        i_phs, alpha_s_sub, alpha_qed_sub, term%connected%has_matrix, sqme)
                end do
             else
                call pcm_work%set_radiation_event ()
                emitter = kin%emitter; i_phs = kin%i_phs
                do i = 1, term%connected%trace%get_qn_index_n_flv ()
                   i_qn = term%connected%trace%get_qn_index (i)
                   if (term%passed_array(i) .or. .not. term%passed) then
                      pcm_work%real_sub%sqme_real_non_sub (i, i_phs) = &
                        real (term%connected%trace%get_matrix_element (i_qn))
                   else
                      pcm_work%real_sub%sqme_real_non_sub (i, i_phs) = zero
                   end if
                end do
                call pcm_work%real_sub%compute (emitter, i_phs, alpha_s_sub, &
                     alpha_qed_sub, term%connected%has_matrix, sqme)
             end if
          end select
       end select
    end select
    if (term%connected%has_trace) &
         call term%connected%trace%set_only_matrix_element &
              (1, cmplx (sum(sqme), 0, default))
    select type (pcm => term%pcm)
    type is (pcm_nlo_t)
       is_subtraction = kin%emitter < 0
       if (term%connected%has_matrix) &
            call refill_evaluator (cmplx (sqme * term%weight, 0, default), &
                 pcm%get_qn (is_subtraction), &
                 pcm%region_data%get_flavor_indices (is_subtraction), &
                 term%connected%matrix)
       if (term%connected%has_flows) &
            call refill_evaluator (cmplx (sqme * term%weight, 0, default), &
                 pcm%get_qn (is_subtraction), &
                 pcm%region_data%get_flavor_indices (is_subtraction), &
                 term%connected%flows)
    end select
  end subroutine term_instance_apply_fks

  module subroutine term_instance_evaluate_sqme_virt (term, alpha_s, alpha_qed)
    class(term_instance_t), intent(inout) :: term
    real(default), intent(in) :: alpha_s, alpha_qed
    real(default), dimension(2) :: alpha_coupling
    type(vector4_t), dimension(:), allocatable :: p_born
    real(default), dimension(:), allocatable :: sqme_virt
    integer :: i_flv, i_qn_born, i_qn_virt
    if (term%nlo_type /= NLO_VIRTUAL)  call msg_fatal ("Trying to " // &
         "evaluate virtual matrix element with unsuited term_instance.")
    if (debug2_active (D_VIRTUAL)) then
       call msg_debug2 &
            (D_VIRTUAL, "Evaluating virtual-subtracted matrix elements")
       print *, 'ren_scale: ', term%get_ren_scale ()
       print *, 'fac_scale: ', term%get_fac_scale ()
       if (allocated (term%es_scale)) then
          print *, 'ES  scale: ', term%es_scale
       else
          print *, 'ES  scale: ', term%get_ren_scale ()
       end if
    end if
    select type (pcm => term%pcm)
    type is (pcm_nlo_t)
       select type (pcm_work => term%pcm_work)
       type is (pcm_nlo_workspace_t)
          alpha_coupling = [alpha_s, alpha_qed]
          if (debug2_active (D_VIRTUAL)) then
             print *, 'alpha_s: ', alpha_coupling (1)
             print *, 'alpha_qed: ', alpha_coupling (2)
          end if
          allocate (p_born (pcm%region_data%n_legs_born))
          if (pcm%settings%factorization_mode == FACTORIZATION_THRESHOLD) then
             p_born = pcm_work%real_kinematics%p_born_onshell%get_momenta(1)
          else
             p_born = term%int_hard%get_momenta ()
          end if
          call pcm_work%set_momenta_and_scales_virtual &
               (p_born, term%ren_scale, term%get_fac_scale (), &
               term%es_scale)
          select type (pcm_work => term%pcm_work)
          type is (pcm_nlo_workspace_t)
             associate (virtual => pcm_work%virtual)
               do i_flv = 1, term%connected%trace%get_qn_index_n_flv ()
                  i_qn_born = term%connected%trace%get_qn_index (i_flv, i_sub = 0)
                  i_qn_virt = term%connected%trace%get_qn_index (i_flv, i_sub = 1)
                  if (term%passed_array(i_flv) .or. .not. term%passed) then
                     virtual%sqme_born(i_flv) = &
                          real (term%connected%trace%get_matrix_element (i_qn_born))
                     virtual%sqme_virt_fin(i_flv) = &
                          real (term%connected%trace%get_matrix_element (i_qn_virt))
                  else
                     virtual%sqme_born(i_flv) = zero
                     virtual%sqme_virt_fin(i_flv) = zero
                  end if
               end do
             end associate
          end select
          call pcm_work%compute_sqme_virt (term%pcm, term%p_hard, &
               alpha_coupling, term%connected%has_matrix, sqme_virt)
          call term%connected%trace%set_only_matrix_element &
               (1, cmplx (sum(sqme_virt), 0, default))
          if (term%connected%has_matrix) &
               call refill_evaluator (cmplx (sqme_virt * term%weight, &
                   0, default), pcm%get_qn (.true.), &
                   remove_duplicates_from_int_array ( &
                   pcm%region_data%get_flavor_indices (.true.)), &
                   term%connected%matrix)
          if (term%connected%has_flows) &
               call refill_evaluator (cmplx (sqme_virt * term%weight, &
                    0, default), pcm%get_qn (.true.), &
                    remove_duplicates_from_int_array ( &
                    pcm%region_data%get_flavor_indices (.true.)), &
                    term%connected%flows)
       end select
    end select
  end subroutine term_instance_evaluate_sqme_virt

  module subroutine term_instance_evaluate_sqme_mismatch (term, alpha_s)
    class(term_instance_t), intent(inout) :: term
    real(default), intent(in) :: alpha_s
    real(default), dimension(:), allocatable :: sqme_mism
    if (term%nlo_type /= NLO_MISMATCH) call msg_fatal &
       ("Trying to evaluate soft mismatch with unsuited term_instance.")
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       call pcm_work%compute_sqme_mismatch &
            (term%pcm, alpha_s, term%connected%has_matrix, sqme_mism)
    end select
    call term%connected%trace%set_only_matrix_element &
         (1, cmplx (sum (sqme_mism) * term%weight, 0, default))
    if (term%connected%has_matrix) then
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          if (term%connected%has_matrix) &
               call refill_evaluator (cmplx (sqme_mism * term%weight, 0, default), &
                    pcm%get_qn (.true.), &
                    remove_duplicates_from_int_array ( &
                    pcm%region_data%get_flavor_indices (.true.)), &
                    term%connected%matrix)
          if (term%connected%has_flows) &
               call refill_evaluator (cmplx (sqme_mism * term%weight, 0, default), &
                    pcm%get_qn (.true.), &
                    remove_duplicates_from_int_array ( &
                    pcm%region_data%get_flavor_indices (.true.)), &
                    term%connected%flows)
       end select
    end if
  end subroutine term_instance_evaluate_sqme_mismatch

  module subroutine term_instance_evaluate_sqme_dglap (term, alpha_s, alpha_qed)
    class(term_instance_t), intent(inout) :: term
    real(default), intent(in) :: alpha_s, alpha_qed
    real(default), dimension(2) :: alpha_coupling
    real(default), dimension(:), allocatable :: sqme_dglap
    integer :: i_flv
    if (term%nlo_type /= NLO_DGLAP) call msg_fatal &
       ("Trying to evaluate DGLAP remnant with unsuited term_instance.")
    if (debug_on)  call msg_debug2 &
         (D_PROCESS_INTEGRATION, "term_instance_evaluate_sqme_dglap")
    select type (pcm => term%pcm)
    type is (pcm_nlo_t)
       select type (pcm_work => term%pcm_work)
       type is (pcm_nlo_workspace_t)
          alpha_coupling = [alpha_s,alpha_qed]
          if (debug2_active (D_PROCESS_INTEGRATION)) then
             associate (n_flv => pcm_work%dglap_remnant%reg_data%n_flv_born)
               print *, "size(sqme_born) = ", &
                    size (pcm_work%dglap_remnant%sqme_born)
               call term%connected%trace%write ()
             end associate
          end if
          call pcm_work%compute_sqme_dglap_remnant (pcm, alpha_coupling, &
               term%connected%has_matrix, sqme_dglap)
       end select
    end select
    call term%connected%trace%set_only_matrix_element &
         (1, cmplx (sum (sqme_dglap) * term%weight, 0, default))
    if (term%connected%has_matrix) then
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          call refill_evaluator (cmplx (sqme_dglap * term%weight, 0, default), &
               pcm%get_qn (.true.), &
               remove_duplicates_from_int_array ( &
               pcm%region_data%get_flavor_indices (.true.)), &
               term%connected%matrix)
          if (term%connected%has_flows) then
             call refill_evaluator &
                  (cmplx (sqme_dglap * term%weight, 0, default), &
                  pcm%get_qn (.true.), &
                  remove_duplicates_from_int_array ( &
                  pcm%region_data%get_flavor_indices (.true.)), &
                  term%connected%flows)
          end if
       end select
    end if
  end subroutine term_instance_evaluate_sqme_dglap

  module subroutine term_instance_reset (term)
    class(term_instance_t), intent(inout) :: term
    call term%connected%reset_expressions ()
    if (allocated (term%alpha_qcd_forced))  deallocate (term%alpha_qcd_forced)
    term%active = .false.
  end subroutine term_instance_reset

  module subroutine term_instance_set_alpha_qcd_forced (term, alpha_qcd)
    class(term_instance_t), intent(inout) :: term
    real(default), intent(in) :: alpha_qcd
    if (allocated (term%alpha_qcd_forced)) then
       term%alpha_qcd_forced = alpha_qcd
    else
       allocate (term%alpha_qcd_forced, source = alpha_qcd)
    end if
  end subroutine term_instance_set_alpha_qcd_forced

  module subroutine term_instance_compute_eff_kinematics (term)
    class(term_instance_t), intent(inout) :: term
    term%checked = .false.
    term%passed = .false.
    call term%isolated%receive_kinematics ()
    call term%connected%receive_kinematics ()
  end subroutine term_instance_compute_eff_kinematics

  module subroutine term_instance_recover_hard_kinematics (term)
    class(term_instance_t), intent(inout) :: term
    term%checked = .false.
    term%passed = .false.
    call term%connected%send_kinematics ()
    call term%isolated%send_kinematics ()
  end subroutine term_instance_recover_hard_kinematics

  module subroutine term_instance_evaluate_expressions &
       (term, config, scale_forced)
    class(term_instance_t), intent(inout) :: term
    type(process_beam_config_t), intent(in) :: config
    real(default), intent(in), allocatable, optional :: scale_forced
    real(default) :: scale = 0
    real(default) :: weight = 1
    real(default), allocatable :: fac_scale, ren_scale
    type(interaction_t), pointer :: src_int
    type(state_matrix_t), pointer :: state_matrix
    type(flavor_t), dimension(:), allocatable :: flv_int, flv_src, f_in, f_out
    logical :: passed
    integer :: n_in, n_vir, n_out, n_tot, n_flv
    integer :: i, j, k
    n_flv = term%connected%trace%get_qn_index_n_flv ()
    if (.not. allocated (term%passed_array)) allocate (term%passed_array(n_flv))
    if (term%flv_dep_cut_eval) then
       do k = 1, n_flv
          if (k == term%i_flv_to_i_flv_rep(k)) then
             n_in = term%int_hard%get_n_in ()
             associate (int_eff => term%isolated%int_eff)
               state_matrix => int_eff%get_state_matrix_ptr ()
               n_tot = int_eff%get_n_tot  ()
               flv_int = quantum_numbers_get_flavor &
                    (state_matrix%get_quantum_number (k))
               allocate (f_in (n_in))
               f_in = flv_int(1:n_in)
               deallocate (flv_int)
             end associate
             n_in = term%connected%trace%get_n_in ()
             n_vir = term%connected%trace%get_n_vir ()
             n_out = term%connected%trace%get_n_out ()
             allocate (f_out (n_out))
             do j = 1, n_out
                call term%connected%trace%find_source &
                     (n_in + n_vir + j, src_int, i)
                if (associated (src_int)) then
                   state_matrix => src_int%get_state_matrix_ptr ()
                   flv_src = quantum_numbers_get_flavor &
                        (state_matrix%get_quantum_number (k))
                   f_out(j) = flv_src(i)
                   deallocate (flv_src)
                end if
             end do

             call term%connected%renew_flv_content_subevt &
                  (term%isolated%sf_chain_eff, &
                   config%data%flv, f_in, f_out)
             call term%connected%evaluate_expressions (passed, &
               scale, fac_scale, ren_scale, weight, &
               scale_forced, force_evaluation = .true.)
             if (k == 1) then
                term%scale = scale
                if (allocated (fac_scale)) then
                   if (.not. allocated (term%fac_scale)) then
                      allocate (term%fac_scale, source = fac_scale)
                   else
                      term%fac_scale = fac_scale
                   end if
                end if
                if (allocated (ren_scale)) then
                   if (.not. allocated (term%ren_scale)) then
                      allocate (term%ren_scale, source = ren_scale)
                   else
                      term%ren_scale = ren_scale
                   end if
                end if
                term%weight = weight
             end if
             term%passed_array(k) = passed
             deallocate (f_in)
             deallocate (f_out)
          else
             term%passed_array(k) = term%passed_array(term%i_flv_to_i_flv_rep(k))
          end if
       end do
       term%passed = any (term%passed_array)
    else
       call term%connected%evaluate_expressions (term%passed, &
            term%scale, term%fac_scale, term%ren_scale, term%weight, &
            scale_forced, force_evaluation = .true.)
       term%passed_array = term%passed
    end if
    term%checked = .true.
  end subroutine term_instance_evaluate_expressions

  module subroutine term_instance_evaluate_interaction (term, core, kin)
    class(term_instance_t), intent(inout) :: term
    class(prc_core_t), intent(in), pointer :: core
    type(kinematics_t), intent(inout) :: kin
    if (debug_on) call msg_debug2 (D_PROCESS_INTEGRATION, &
         "term_instance_evaluate_interaction")
    if (kin%only_cm_frame .and. (.not. kin%lab_is_cm())) then
         term%p_hard = kin%get_boost_to_cms () * term%int_hard%get_momenta ()
    else
         term%p_hard = term%int_hard%get_momenta ()
    end if
    select type (core)
    class is (prc_external_t)
       call term%evaluate_interaction_external (core, kin)
    class default
       call term%evaluate_interaction_default (core)
    end select
    call term%int_hard%set_matrix_element (term%amp)
  end subroutine term_instance_evaluate_interaction

  module subroutine term_instance_evaluate_interaction_default (term, core)
    class(term_instance_t), intent(inout) :: term
    class(prc_core_t), intent(in) :: core
    real(default) :: fac_scale, ren_scale
    integer :: i
    if (allocated (term%fac_scale)) then
       fac_scale = term%fac_scale
    else
       fac_scale = term%scale
    end if
    if (allocated (term%ren_scale)) then
       ren_scale = term%ren_scale
    else
       ren_scale = term%scale
    end if
    do i = 1, term%config%n_allowed
       term%amp(i) = core%compute_amplitude (term%config%i_term, term%p_hard, &
            term%config%flv(i), term%config%hel(i), term%config%col(i), &
            fac_scale, ren_scale, term%alpha_qcd_forced, &
            term%core_state)
    end do
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       call pcm_work%set_fac_scale (fac_scale)
    end select
  end subroutine term_instance_evaluate_interaction_default

  module subroutine term_instance_evaluate_interaction_external &
       (term, core, kin)
    class(term_instance_t), intent(inout) :: term
    class(prc_core_t), intent(inout) :: core
    type(kinematics_t), intent(inout) :: kin
    if (debug_on) call msg_debug2 (D_PROCESS_INTEGRATION, &
         "term_instance_evaluate_interaction_external")
    select type (core_state => term%core_state)
    type is (openloops_state_t)
       select type (core)
       type is (prc_openloops_t)
          call core%compute_alpha_s (core_state, term%get_ren_scale ())
          if (allocated (core_state%threshold_data)) &
               call evaluate_threshold_parameters (core_state, core, kin%phs%get_sqrts ())
       end select
    class is (prc_external_state_t)
       select type (core)
       class is (prc_external_t)
          call core%compute_alpha_s (core_state, term%get_ren_scale ())
       end select
    end select
    call evaluate_threshold_interaction ()
    if (term%nlo_type == NLO_VIRTUAL) then
       call term%evaluate_interaction_external_loop (core)
    else
       call term%evaluate_interaction_external_tree (core)
    end if
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       call pcm_work%set_fac_scale (term%get_fac_scale ())
    end select

  contains
    subroutine evaluate_threshold_parameters (core_state, core, sqrts)
       type(openloops_state_t), intent(inout) :: core_state
       type(prc_openloops_t), intent(inout) :: core
       real(default), intent(in) :: sqrts
       real(default) :: mtop, wtop
       mtop = m1s_to_mpole (sqrts)
       wtop = core_state%threshold_data%compute_top_width &
              (mtop, core_state%alpha_qcd)
       call core%set_mass_and_width (6, mtop, wtop)
    end subroutine

    subroutine evaluate_threshold_interaction ()
       integer :: leg
       select type (core)
       type is (prc_threshold_t)
          if (term%nlo_type > BORN) then
             select type (pcm_work => term%pcm_work)
             type is (pcm_nlo_workspace_t)
                if (kin%emitter >= 0) then
                   call core%set_offshell_momenta &
                        (pcm_work%real_kinematics%p_real_cms%get_momenta(term%config%i_term))
                   leg = thr_leg (kin%emitter)
                   call core%set_leg (leg)
                   call core%set_onshell_momenta &
                        (pcm_work%real_kinematics%p_real_onshell(leg)%get_momenta(term%config%i_term))
                else
                   call core%set_leg (0)
                   call core%set_offshell_momenta &
                        (pcm_work%real_kinematics%p_born_cms%get_momenta(1))
                end if
             end select
          else
             call core%set_leg (-1)
             call core%set_offshell_momenta (term%p_hard)
          end if
       end select
    end subroutine evaluate_threshold_interaction
  end subroutine term_instance_evaluate_interaction_external

  module subroutine term_instance_evaluate_interaction_external_tree &
       (term, core)
    class(term_instance_t), intent(inout) :: term
    class(prc_core_t), intent(inout) :: core
    real(default) :: sqme
    real(default), dimension(:), allocatable :: sqme_color_c
    real(default), dimension(:), allocatable :: sqme_spin_c
    real(default), dimension(6) :: sqme_spin_c_tmp
    integer :: n_flv, n_hel, n_sub_color, n_sub_spin, n_pdf_off
    integer :: i_flv, i_hel, i_sub, i_color_c, i_color_c_eqv, &
         i_spin_c, i_spin_c_eqv
    integer :: i_flv_eqv, i_hel_eqv
    integer :: emitter, i_emitter
    logical :: bad_point, bp
    logical, dimension(:,:), allocatable :: eqv_me_evaluated
    if (debug_on) call msg_debug2 (D_PROCESS_INTEGRATION, &
         "term_instance_evaluate_interaction_external_tree")
    allocate (sqme_color_c (blha_result_array_size &
         (term%int_hard%get_n_tot (), BLHA_AMP_COLOR_C)))
    n_flv = term%int_hard%get_qn_index_n_flv ()
    n_hel = term%int_hard%get_qn_index_n_hel ()
    n_sub_color = term%get_n_sub_color ()
    n_sub_spin = term%get_n_sub_spin ()
    allocate (eqv_me_evaluated(n_flv,n_hel))
    eqv_me_evaluated = .false.
    do i_flv = 1, n_flv
       if (.not. term%passed_array(i_flv) .and. term%passed) cycle
       do i_hel = 1, n_hel
          i_flv_eqv = core%data%eqv_flv_index(i_flv)
          i_hel_eqv = core%data%eqv_hel_index(i_hel)
          if (.not. eqv_me_evaluated(i_flv_eqv, i_hel_eqv)) then
             select type (core)
             class is (prc_external_t)
                call core%update_alpha_s (term%core_state, term%get_ren_scale ())
                call core%compute_sqme (i_flv, i_hel, term%p_hard, &
                     term%get_ren_scale (), sqme, bad_point)
                call term%pcm_work%set_bad_point (bad_point)
                associate (i_int => term%int_hard%get_qn_index &
                     (i_flv = i_flv, i_hel = i_hel, i_sub = 0))
                  term%amp(i_int) = cmplx (sqme, 0, default)
                end associate
             end select
             n_pdf_off = 0
             if (term%pcm%has_pdfs .and. &
                  (term%is_subtraction () .or. term%nlo_type == NLO_DGLAP)) then
                n_pdf_off = n_pdf_off + n_beams_rescaled
                do i_sub = 1, n_pdf_off
                   term%amp(term%int_hard%get_qn_index (i_flv, i_hel, i_sub)) = &
                        term%amp(term%int_hard%get_qn_index (i_flv, i_hel, i_sub = 0))
                end do
             end if
             if (term%pcm%has_pdfs .and. term%nlo_type == NLO_DGLAP) then
                sqme_color_c = zero
                select type (pcm => term%pcm)
                type is (pcm_nlo_t)
                   if (pcm%settings%nlo_correction_type == "EW" .and. &
                      pcm%region_data%alphas_power > 0) then
                      select type (core)
                      class is (prc_blha_t)
                         call core%compute_sqme_color_c_raw (i_flv, i_hel, &
                              term%p_hard, term%get_ren_scale (), sqme_color_c, &
                              bad_point)
                         call term%pcm_work%set_bad_point (bad_point)
                      class is (prc_recola_t)
                         call core%compute_sqme_color_c_raw (i_flv, i_hel, &
                              term%p_hard, term%get_ren_scale (), sqme_color_c, &
                              bad_point)
                         call term%pcm_work%set_bad_point (bad_point)
                      end select
                   end if
                end select
                do i_sub = 1, n_sub_color
                   i_color_c = term%int_hard%get_qn_index &
                        (i_flv, i_hel, i_sub + n_pdf_off)
                   term%amp(i_color_c) = cmplx (sqme_color_c(i_sub), 0, default)
                end do
             end if
             if ((term%nlo_type == NLO_REAL .and. term%is_subtraction ()) .or. &
                  term%nlo_type == NLO_MISMATCH) then
                sqme_color_c = zero
                select type (core)
                class is (prc_blha_t)
                   call core%compute_sqme_color_c_raw (i_flv, i_hel, &
                        term%p_hard, term%get_ren_scale (), sqme_color_c, bad_point)
                   call term%pcm_work%set_bad_point (bad_point)
                class is (prc_recola_t)
                   call core%compute_sqme_color_c_raw (i_flv, i_hel, &
                        term%p_hard, term%get_ren_scale (), sqme_color_c, bad_point)
                   call term%pcm_work%set_bad_point (bad_point)
                end select
                do i_sub = 1, n_sub_color
                   i_color_c = term%int_hard%get_qn_index &
                        (i_flv, i_hel, i_sub + n_pdf_off)
                   term%amp(i_color_c) = cmplx (sqme_color_c(i_sub), 0, default)
                end do
                if (n_sub_spin > 0) then
                   bad_point = .false.
                   allocate (sqme_spin_c(0))
                   select type (core)
                   type is (prc_openloops_t)
                      select type (pcm => term%pcm)
                      type is (pcm_nlo_t)
                         do i_emitter = 1, pcm%region_data%n_emitters
                            emitter = pcm%region_data%emitters(i_emitter)
                            if (emitter > 0) then
                               call core%compute_sqme_spin_c &
                                    (i_flv, &
                                    i_hel, &
                                    emitter, &
                                    term%p_hard, &
                                    term%get_ren_scale (), &
                                    sqme_spin_c_tmp, &
                                    bp)
                               sqme_spin_c = [sqme_spin_c, sqme_spin_c_tmp]
                               bad_point = bad_point .or. bp
                            end if
                         end do
                      end select
                      do i_sub = 1, n_sub_spin
                         i_spin_c = term%int_hard%get_qn_index (i_flv, i_hel, &
                              i_sub + n_pdf_off + n_sub_color)
                         term%amp(i_spin_c) = cmplx &
                              (sqme_spin_c(i_sub), 0, default)
                      end do
                   end select
                   deallocate (sqme_spin_c)
                end if
             end if
             eqv_me_evaluated(i_flv_eqv, i_hel_eqv) = .true.
          else
             associate (i_int => term%int_hard%get_qn_index &
                     (i_flv = i_flv, i_hel = i_hel, i_sub = 0), &
                     i_int_eqv => term%int_hard%get_qn_index &
                     (i_flv = i_flv_eqv, i_hel = i_hel_eqv, i_sub = 0))
                term%amp(i_int) = term%amp(i_int_eqv)
             end associate
             n_pdf_off = 0
             if (term%pcm%has_pdfs .and. &
                  (term%is_subtraction () .or. term%nlo_type == NLO_DGLAP)) then
                n_pdf_off = n_pdf_off + n_beams_rescaled
                do i_sub = 1, n_pdf_off
                   term%amp(term%int_hard%get_qn_index (i_flv, i_hel, i_sub)) = &
                        term%amp(term%int_hard%get_qn_index (i_flv, i_hel, i_sub = 0))
                end do
             end if
             if (term%pcm%has_pdfs .and. term%nlo_type == NLO_DGLAP) then
                do i_sub = 1, n_sub_color
                   i_color_c = term%int_hard%get_qn_index &
                        (i_flv, i_hel, i_sub + n_pdf_off)
                   i_color_c_eqv = term%int_hard%get_qn_index &
                        (i_flv_eqv, i_hel_eqv, i_sub + n_pdf_off)
                   term%amp(i_color_c) = term%amp(i_color_c_eqv)
                end do
             end if
             if ((term%nlo_type == NLO_REAL .and. term%is_subtraction ()) .or. &
                  term%nlo_type == NLO_MISMATCH) then
                do i_sub = 1, n_sub_color
                   i_color_c = term%int_hard%get_qn_index &
                        (i_flv, i_hel, i_sub + n_pdf_off)
                   i_color_c_eqv = term%int_hard%get_qn_index &
                        (i_flv_eqv, i_hel_eqv, i_sub + n_pdf_off)
                   term%amp(i_color_c) = term%amp(i_color_c_eqv)
                end do
                do i_sub = 1, n_sub_spin
                   i_spin_c = term%int_hard%get_qn_index (i_flv, i_hel, &
                        i_sub + n_pdf_off + n_sub_color)
                   i_spin_c_eqv = term%int_hard%get_qn_index (i_flv_eqv, i_hel_eqv, &
                        i_sub + n_pdf_off + n_sub_color)
                   term%amp(i_spin_c) = term%amp(i_spin_c_eqv)
                end do
             end if
          end if
       end do
    end do
  end subroutine term_instance_evaluate_interaction_external_tree

  module subroutine term_instance_evaluate_interaction_external_loop &
       (term, core)
    class(term_instance_t), intent(inout) :: term
    class(prc_core_t), intent(in) :: core
    integer :: n_hel, n_sub, n_flv
    integer :: i, i_flv, i_hel, i_sub, i_virt, i_color_c, i_color_c_eqv
    integer :: i_flv_eqv, i_hel_eqv
    real(default), dimension(4) :: sqme_virt
    real(default), dimension(:), allocatable :: sqme_color_c
    real(default) :: es_scale
    logical :: bad_point
    logical, dimension(:,:), allocatable :: eqv_me_evaluated
    if (debug_on) call msg_debug (D_PROCESS_INTEGRATION, &
         "term_instance_evaluate_interaction_external_loop")
    allocate (sqme_color_c (blha_result_array_size &
         (term%int_hard%get_n_tot (), BLHA_AMP_COLOR_C)))
    n_flv = term%int_hard%get_qn_index_n_flv ()
    n_hel = term%int_hard%get_qn_index_n_hel ()
    n_sub = term%int_hard%get_qn_index_n_sub ()
    allocate (eqv_me_evaluated(n_flv,n_hel))
    eqv_me_evaluated = .false.
    i_virt = 1
    do i_flv = 1, n_flv
       if (.not. term%passed_array(i_flv) .and. term%passed) cycle
       do i_hel = 1, n_hel
          i_flv_eqv = core%data%eqv_flv_index(i_flv)
          i_hel_eqv = core%data%eqv_hel_index(i_hel)
          if (.not. eqv_me_evaluated(i_flv_eqv, i_hel_eqv)) then
             select type (core)
             class is (prc_external_t)
                if (allocated (term%es_scale)) then
                   es_scale = term%es_scale
                else
                   es_scale = term%get_ren_scale ()
                end if
                call core%compute_sqme_virt (i_flv, i_hel, term%p_hard, &
                     term%get_ren_scale (), es_scale, &
                     term%pcm%blha_defaults%loop_method, &
                     sqme_virt, bad_point)
                call term%pcm_work%set_bad_point (bad_point)
             end select
             associate (i_born => term%int_hard%get_qn_index (i_flv, i_hel = i_hel, i_sub = 0), &
                  i_loop => term%int_hard%get_qn_index (i_flv, i_hel = i_hel, i_sub = i_virt))
               term%amp(i_loop) = cmplx (sqme_virt(3), 0, default)
               term%amp(i_born) = cmplx (sqme_virt(4), 0, default)
             end associate
             select type (pcm => term%pcm)
             type is (pcm_nlo_t)
                select type (core)
                class is (prc_blha_t)
                   call core%compute_sqme_color_c_raw (i_flv, i_hel, &
                        term%p_hard, term%get_ren_scale (), &
                        sqme_color_c, bad_point)
                   call term%pcm_work%set_bad_point (bad_point)
                   do i_sub = 1 + i_virt, n_sub
                      i_color_c = term%int_hard%get_qn_index &
                           (i_flv, i_hel = i_hel, i_sub = i_sub)
                      ! Index shift: i_sub - i_virt
                      term%amp(i_color_c) = &
                           cmplx (sqme_color_c(i_sub - i_virt), 0, default)
                   end do
                type is (prc_recola_t)
                   call core%compute_sqme_color_c_raw (i_flv, i_hel, &
                        term%p_hard, term%get_ren_scale (), sqme_color_c, bad_point)
                   call term%pcm_work%set_bad_point (bad_point)
                   do i_sub = 1 + i_virt, n_sub
                      i_color_c = term%int_hard%get_qn_index &
                           (i_flv, i_hel = i_hel, i_sub = i_sub)
                      ! Index shift: i_sub - i_virt
                      term%amp(i_color_c) = &
                           cmplx (sqme_color_c(i_sub - i_virt), 0, default)
                   end do
                end select
             end select
             eqv_me_evaluated(i_flv_eqv, i_hel_eqv) = .true.
          else
             associate (i_born => term%int_hard%get_qn_index (i_flv, i_hel = i_hel, i_sub = 0), &
                  i_loop => term%int_hard%get_qn_index (i_flv, i_hel = i_hel, i_sub = i_virt), &
                  i_born_eqv => term%int_hard%get_qn_index &
                  (i_flv_eqv, i_hel = i_hel_eqv, i_sub = 0), &
                  i_loop_eqv => term%int_hard%get_qn_index &
                  (i_flv_eqv, i_hel = i_hel_eqv, i_sub = 1))
               term%amp(i_loop) = term%amp(i_loop_eqv)
               term%amp(i_born) = term%amp(i_born_eqv)
             end associate
             do i_sub = 1 + i_virt, n_sub
                i_color_c = term%int_hard%get_qn_index &
                     (i_flv, i_hel = i_hel, i_sub = i_sub)
                i_color_c_eqv = term%int_hard%get_qn_index &
                     (i_flv_eqv, i_hel = i_hel_eqv, i_sub = i_sub)
                ! Index shift: i_sub - i_virt
                term%amp(i_color_c) = term%amp(i_color_c_eqv)
             end do
          end if
       end do
    end do
   end subroutine term_instance_evaluate_interaction_external_loop

  module subroutine term_instance_evaluate_trace (term, kin)
    class(term_instance_t), intent(inout) :: term
    type(kinematics_t), intent(inout) :: kin
    real(default) :: fac_scale
    if (allocated (term%fac_scale)) then
       fac_scale = term%fac_scale
    else
       fac_scale = term%scale
    end if
    call kin%evaluate_sf_chain (fac_scale, term%negative_sf)
    call term%evaluate_scaled_sf_chains (kin)
    call term%isolated%evaluate_sf_chain (fac_scale)
    call term%isolated%evaluate_trace ()
    call term%connected%evaluate_trace ()
  end subroutine term_instance_evaluate_trace

  module subroutine term_instance_evaluate_event_data (term)
    class(term_instance_t), intent(inout) :: term
    logical :: only_momenta
    only_momenta = term%nlo_type > BORN
    call term%isolated%evaluate_event_data (only_momenta)
    call term%connected%evaluate_event_data (only_momenta)
  end subroutine term_instance_evaluate_event_data

  module subroutine term_instance_set_fac_scale (term, fac_scale)
    class(term_instance_t), intent(inout) :: term
    real(default), intent(in) :: fac_scale
    term%fac_scale = fac_scale
  end subroutine term_instance_set_fac_scale

  module function term_instance_get_fac_scale (term) result (fac_scale)
    class(term_instance_t), intent(in) :: term
    real(default) :: fac_scale
    if (allocated (term%fac_scale)) then
       fac_scale = term%fac_scale
    else
       fac_scale = term%scale
    end if
  end function term_instance_get_fac_scale

  module function term_instance_get_ren_scale (term) result (ren_scale)
    class(term_instance_t), intent(in) :: term
    real(default) :: ren_scale
    if (allocated (term%ren_scale)) then
       ren_scale = term%ren_scale
    else
       ren_scale = term%scale
    end if
  end function term_instance_get_ren_scale

  module function term_instance_get_alpha_s (term, core) result (alpha_s)
    class(term_instance_t), intent(in) :: term
    class(prc_core_t), intent(in) :: core
    real(default) :: alpha_s
    alpha_s = core%get_alpha_s (term%core_state)
    if (alpha_s < zero)  alpha_s = term%config%alpha_s
  end function term_instance_get_alpha_s

  module subroutine term_instance_get_helicities_for_openloops &
       (term, helicities)
    class(term_instance_t), intent(in) :: term
    integer, dimension(:,:), allocatable, intent(out) :: helicities
    type(helicity_t), dimension(:), allocatable :: hel
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn
    type(quantum_numbers_mask_t) :: qn_mask
    integer :: h, i, j, n_in
    call qn_mask%set_sub (1)
    call term%isolated%trace%get_quantum_numbers_mask (qn_mask, qn)
    n_in = term%int_hard%get_n_in ()
    allocate (helicities (size (qn, dim=1), n_in))
    allocate (hel (n_in))
    do i = 1, size (qn, dim=1)
       do j = 1, n_in
          hel(j) = qn(i, j)%get_helicity ()
          call hel(j)%diagonalize ()
          call hel(j)%get_indices (h, h)
          helicities (i, j) = h
       end do
    end do
  end subroutine term_instance_get_helicities_for_openloops

  elemental module function term_instance_get_i_term_global &
       (term) result (i_term)
    integer :: i_term
    class(term_instance_t), intent(in) :: term
    i_term = term%config%i_term_global
  end function term_instance_get_i_term_global

  elemental module function term_instance_is_subtraction (term) result (sub)
    logical :: sub
    class(term_instance_t), intent(in) :: term
    sub = term%config%i_term_global == term%config%i_sub
  end function term_instance_is_subtraction

  module function term_instance_get_n_sub (term) result (n_sub)
    integer :: n_sub
    class(term_instance_t), intent(in) :: term
    n_sub = term%config%n_sub
  end function term_instance_get_n_sub

  module function term_instance_get_n_sub_color (term) result (n_sub_color)
    integer :: n_sub_color
    class(term_instance_t), intent(in) :: term
    n_sub_color = term%config%n_sub_color
  end function term_instance_get_n_sub_color

  module function term_instance_get_n_sub_spin (term) result (n_sub_spin)
    integer :: n_sub_spin
    class(term_instance_t), intent(in) :: term
    n_sub_spin = term%config%n_sub_spin
  end function term_instance_get_n_sub_spin

  module subroutine process_instance_write_header (object, unit, testflag)
    class(process_instance_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    call write_separator (u, 2)
    if (associated (object%process)) then
       call object%process%write_meta (u, testflag)
    else
       write (u, "(1x,A)") "Process instance [undefined process]"
       return
    end if
    write (u, "(3x,A)", advance = "no")  "status = "
    select case (object%evaluation_status)
    case (STAT_INITIAL);            write (u, "(A)")  "initialized"
    case (STAT_ACTIVATED);          write (u, "(A)")  "activated"
    case (STAT_BEAM_MOMENTA);       write (u, "(A)")  "beam momenta set"
    case (STAT_FAILED_KINEMATICS);  write (u, "(A)")  "failed kinematics"
    case (STAT_SEED_KINEMATICS);    write (u, "(A)")  "seed kinematics"
    case (STAT_HARD_KINEMATICS);    write (u, "(A)")  "hard kinematics"
    case (STAT_EFF_KINEMATICS);     write (u, "(A)")  "effective kinematics"
    case (STAT_FAILED_CUTS);        write (u, "(A)")  "failed cuts"
    case (STAT_PASSED_CUTS);        write (u, "(A)")  "passed cuts"
    case (STAT_EVALUATED_TRACE);    write (u, "(A)")  "evaluated trace"
       call write_separator (u)
       write (u, "(3x,A,ES19.12)")  "sqme   = ", object%sqme
    case (STAT_EVENT_COMPLETE);   write (u, "(A)")  "event complete"
       call write_separator (u)
       write (u, "(3x,A,ES19.12)")  "sqme   = ", object%sqme
       write (u, "(3x,A,ES19.12)")  "weight = ", object%weight
       if (.not. vanishes (object%excess)) &
            write (u, "(3x,A,ES19.12)")  "excess = ", object%excess
    case default;                 write (u, "(A)")  "undefined"
    end select
    if (object%i_mci /= 0) then
       call write_separator (u)
       call object%mci_work(object%i_mci)%write (u, testflag)
    end if
    call write_separator (u, 2)
  end subroutine process_instance_write_header

  module subroutine process_instance_write (object, unit, testflag)
    class(process_instance_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u, i
    u = given_output_unit (unit)
    call object%write_header (u)
    if (object%evaluation_status >= STAT_BEAM_MOMENTA) then
       call object%sf_chain%write (u)
       call write_separator (u, 2)
       if (object%evaluation_status >= STAT_SEED_KINEMATICS) then
          if (object%evaluation_status >= STAT_HARD_KINEMATICS) then
             call write_separator (u, 2)
             write (u, "(1x,A)") "Active terms:"
             if (any (object%term%active)) then
                do i = 1, size (object%term)
                   if (object%term(i)%active) then
                      call write_separator (u)
                      call object%term(i)%write (u, &
                           kin = object%kin(i), &
                           show_eff_state = &
                           object%evaluation_status >= STAT_EFF_KINEMATICS, &
                           testflag = testflag)
                   end if
                end do
             end if
          end if
          call write_separator (u, 2)
       end if
    end if
  end subroutine process_instance_write

  module subroutine process_instance_init (instance, process)
    class(process_instance_t), intent(out), target :: instance
    type(process_t), intent(inout), target :: process
    integer :: i
    class(pcm_t), pointer :: pcm
    type(process_term_t), pointer :: term
    type(var_list_t), pointer :: var_list
    integer :: i_born, i_real, i_real_fin, i_component
    if (debug_on)  call msg_debug &
         (D_PROCESS_INTEGRATION, "process_instance_init")
    instance%process => process
    instance%pcm => process%get_pcm_ptr ()
    call instance%process%check_library_sanity ()
    call instance%setup_sf_chain (process%get_beam_config_ptr ())
    allocate (instance%mci_work (process%get_n_mci ()))
    do i = 1, size (instance%mci_work)
       call instance%process%init_mci_work (instance%mci_work(i), i)
    end do
    call instance%process%reset_selected_cores ()
    pcm => instance%process%get_pcm_ptr ()
    call pcm%allocate_workspace (instance%pcm_work)
    select type (pcm)
    type is (pcm_nlo_t)
       !!! The process is kept when the integration is finalized, but not the
       !!! process_instance. Thus, we check whether pcm has been initialized
       !!! but set up the pcm_work each time.
       i_real_fin = process%get_associated_real_fin (1)
       if (.not. pcm%initialized) then
          i_born = pcm%get_i_core (pcm%i_born)
          i_real = pcm%get_i_core (pcm%i_real)
          call pcm%init_qn (process%get_model_ptr ())
          if (i_real_fin > 0) call pcm%allocate_ps_matching ()
          var_list => process%get_var_list_ptr ()
          if (var_list%get_sval (var_str ("$dalitz_plot")) /= var_str ('')) &
               call pcm%activate_dalitz_plot (var_list%get_sval (var_str ("$dalitz_plot")))
       end if
       pcm%initialized = .true.
       select type (pcm_work => instance%pcm_work)
       type is (pcm_nlo_workspace_t)
          call pcm_work%init_config (pcm, &
               process%component_can_be_integrated (), &
               process%get_nlo_type_component (), process%get_energy (), &
               i_real_fin, process%get_model_ptr ())
       end select
    end select
    ! TODO wk-03-01 n_terms will eventually acquire a different meaning
    allocate (instance%kin (process%get_n_terms ()))
    do i = 1, process%get_n_terms ()
       term => process%get_term_ptr (i)
       i_component = term%i_component
       call instance%kin(i)%configure (pcm, instance%pcm_work, &
            instance%sf_chain, &
            process%get_beam_config_ptr (), &
            process%get_phs_config (i_component), &
            process%get_nlo_type_component (i_component), &
            term%i_sub == i)
    end do
    ! TODO wk-03-01 n_terms will eventually acquire a different meaning
    allocate (instance%term (process%get_n_terms ()))
    do i = 1, process%get_n_terms ()
       call instance%term(i)%configure (process, i, instance%pcm_work, &
            instance%sf_chain, instance%kin(i))
    end do
    call instance%set_i_mci_to_real_component ()
    call instance%find_same_kinematics ()
    instance%evaluation_status = STAT_INITIAL
  end subroutine process_instance_init

  module subroutine process_instance_final (instance)
    class(process_instance_t), intent(inout) :: instance
    class(process_instance_hook_t), pointer :: current
    integer :: i
    instance%process => null ()
    if (allocated (instance%mci_work)) then
       do i = 1, size (instance%mci_work)
          call instance%mci_work(i)%final ()
       end do
       deallocate (instance%mci_work)
    end if
    call instance%sf_chain%final ()
    if (allocated (instance%kin)) then
       do i = 1, size (instance%kin)
          call instance%kin(i)%final ()
       end do
       deallocate (instance%kin)
    end if
    if (allocated (instance%term)) then
       do i = 1, size (instance%term)
          call instance%term(i)%final ()
       end do
       deallocate (instance%term)
    end if
    call instance%pcm_work%final ()
    instance%evaluation_status = STAT_UNDEFINED
    do while (associated (instance%hook))
       current => instance%hook
       call current%final ()
       instance%hook => current%next
       deallocate (current)
    end do
    instance%hook => null ()
  end subroutine process_instance_final

  module subroutine process_instance_reset (instance, reset_mci)
    class(process_instance_t), intent(inout), target :: instance
    logical, intent(in), optional :: reset_mci
    integer :: i
    call instance%process%reset_selected_cores ()
    do i = 1, size (instance%term)
       call instance%term(i)%reset ()
    end do
    instance%term%checked = .false.
    instance%term%passed = .false.
    instance%kin%new_seed = .true.
    if (present (reset_mci)) then
       if (reset_mci)  instance%i_mci = 0
    end if
    instance%selected_channel = 0
    instance%evaluation_status = STAT_INITIAL
  end subroutine process_instance_reset

  module subroutine process_instance_sampler_test (instance, i_mci, n_calls)
    class(process_instance_t), intent(inout), target :: instance
    integer, intent(in) :: i_mci
    integer, intent(in) :: n_calls
    integer :: i_mci_work
    i_mci_work = instance%process%get_i_mci_work (i_mci)
    call instance%choose_mci (i_mci_work)
    call instance%reset_counter ()
    call instance%process%sampler_test (instance, n_calls, i_mci_work)
    call instance%process%set_counter_mci_entry (i_mci_work, instance%get_counter ())
  end subroutine process_instance_sampler_test

  module subroutine process_instance_generate_weighted_event (instance, i_mci)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_mci
    integer :: i_mci_work
    i_mci_work = instance%process%get_i_mci_work (i_mci)
    call instance%choose_mci (i_mci_work)
    associate (mci_work => instance%mci_work(i_mci_work))
       call instance%process%generate_weighted_event &
          (i_mci_work, mci_work, instance, &
           instance%keep_failed_events ())
    end associate
  end subroutine process_instance_generate_weighted_event

  module subroutine process_instance_generate_unweighted_event (instance, i_mci)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_mci
    integer :: i_mci_work
    i_mci_work = instance%process%get_i_mci_work (i_mci)
    call instance%choose_mci (i_mci_work)
    associate (mci_work => instance%mci_work(i_mci_work))
       call instance%process%generate_unweighted_event &
          (i_mci_work, mci_work, instance)
    end associate
  end subroutine process_instance_generate_unweighted_event

  module subroutine process_instance_recover_event (instance)
    class(process_instance_t), intent(inout) :: instance
    integer :: i_mci
    i_mci = instance%i_mci
    call instance%process%set_i_mci_work (i_mci)
    associate (mci_instance => instance%mci_work(i_mci)%mci)
      call mci_instance%fetch (instance, instance%selected_channel)
    end associate
  end subroutine process_instance_recover_event

  module subroutine process_instance_activate (instance)
    class(process_instance_t), intent(inout) :: instance
    integer :: i, j
    integer, dimension(:), allocatable :: i_term
    associate (mci_work => instance%mci_work(instance%i_mci))
       call instance%process%select_components &
            (mci_work%get_active_components ())
    end associate
    associate (process => instance%process)
       do i = 1, instance%process%get_n_components ()
          if (instance%process%component_is_selected (i)) then
             allocate (i_term (size (process%get_component_i_terms (i))))
             i_term = process%get_component_i_terms (i)
             do j = 1, size (i_term)
                instance%term(i_term(j))%active = .true.
             end do
          end if
          if (allocated (i_term)) deallocate (i_term)
       end do
    end associate
    instance%evaluation_status = STAT_ACTIVATED
  end subroutine process_instance_activate

  module subroutine process_instance_find_same_kinematics (instance)
    class(process_instance_t), intent(inout) :: instance
    integer :: i_term1, i_term2, k, n_same
    do i_term1 = 1, size (instance%term)
       if (.not. allocated (instance%term(i_term1)%same_kinematics)) then
          n_same = 1 !!! Index group includes the index of its term_instance
          do i_term2 = 1, size (instance%term)
             if (i_term1 == i_term2) cycle
             if (compare_md5s (i_term1, i_term2)) n_same = n_same + 1
          end do
          allocate (instance%term(i_term1)%same_kinematics (n_same))
          associate (same_kinematics1 => instance%term(i_term1)%same_kinematics)
             same_kinematics1 = 0
             k = 1
             do i_term2 = 1, size (instance%term)
                if (compare_md5s (i_term1, i_term2)) then
                   same_kinematics1(k) = i_term2
                   k = k + 1
                end if
             end do
             do k = 1, size (same_kinematics1)
                if (same_kinematics1(k) == i_term1) cycle
                i_term2 = same_kinematics1(k)
                allocate (instance%term(i_term2)%same_kinematics (n_same))
                instance%term(i_term2)%same_kinematics = same_kinematics1
             end do
          end associate
       end if
    end do
  contains
    function compare_md5s (i, j) result (same)
      logical :: same
      integer, intent(in) :: i, j
      character(32) :: md5sum_1, md5sum_2
      integer :: mode_1, mode_2
      mode_1 = 0; mode_2 = 0
      select type (phs => instance%kin(i)%phs%config)
      type is (phs_fks_config_t)
         md5sum_1 = phs%md5sum_born_config
         mode_1 = phs%mode
      class default
         md5sum_1 = phs%md5sum_phs_config
      end select
      select type (phs => instance%kin(j)%phs%config)
      type is (phs_fks_config_t)
         md5sum_2 = phs%md5sum_born_config
         mode_2 = phs%mode
      class default
         md5sum_2 = phs%md5sum_phs_config
      end select
      same = (md5sum_1 == md5sum_2) .and. (mode_1 == mode_2)
    end function compare_md5s
  end subroutine process_instance_find_same_kinematics

  module subroutine process_instance_transfer_same_kinematics (instance, i_term)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_term
    integer :: i, i_term_same
    associate (same_kinematics => instance%term(i_term)%same_kinematics)
      do i = 1, size (same_kinematics)
         i_term_same = same_kinematics(i)
         instance%term(i_term_same)%p_seed = instance%term(i_term)%p_seed
         associate (phs => instance%kin(i_term_same)%phs)
           call phs%set_lorentz_transformation &
                (instance%kin(i_term)%phs%get_lorentz_transformation ())
           select type (phs)
           type is (phs_fks_t)
              call phs%set_momenta (instance%term(i_term_same)%p_seed)
              if (i_term_same /= i_term) then
                 call phs%set_reference_frames (.false.)
              end if
           end select
         end associate
         instance%kin(i_term_same)%new_seed = .false.
      end do
    end associate
  end subroutine process_instance_transfer_same_kinematics

  module subroutine process_instance_redo_sf_chains &
       (instance, i_term, phs_channel)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in), dimension(:) :: i_term
    integer, intent(in) :: phs_channel
    integer :: i
    do i = 1, size (i_term)
       call instance%kin(i_term(i))%redo_sf_chain &
            (instance%mci_work(instance%i_mci), phs_channel)
    end do
  end subroutine process_instance_redo_sf_chains

  module subroutine process_instance_integrate (instance, i_mci, &
       n_it, n_calls, adapt_grids, adapt_weights, final, pacify)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_mci
    integer, intent(in) :: n_it
    integer, intent(in) :: n_calls
    logical, intent(in), optional :: adapt_grids
    logical, intent(in), optional :: adapt_weights
    logical, intent(in), optional :: final, pacify
    integer :: nlo_type, i_mci_work
    nlo_type = instance%process%get_component_nlo_type (i_mci)
    i_mci_work = instance%process%get_i_mci_work (i_mci)
    call instance%choose_mci (i_mci_work)
    call instance%reset_counter ()
    associate (mci_work => instance%mci_work(i_mci_work), &
               process => instance%process)
       call process%integrate (i_mci_work, mci_work, &
            instance, n_it, n_calls, adapt_grids, adapt_weights, &
            final, pacify, nlo_type = nlo_type)
       call process%set_counter_mci_entry (i_mci_work, instance%get_counter ())
    end associate
  end subroutine process_instance_integrate

  module subroutine process_instance_setup_sf_chain (instance, config)
    class(process_instance_t), intent(inout) :: instance
    type(process_beam_config_t), intent(in), target :: config
    integer :: n_strfun
    n_strfun = config%n_strfun
    if (n_strfun /= 0) then
       call instance%sf_chain%init (config%data, config%sf)
    else
       call instance%sf_chain%init (config%data)
    end if
    if (config%sf_trace) then
       call instance%sf_chain%setup_tracing (config%sf_trace_file)
    end if
  end subroutine process_instance_setup_sf_chain

  module subroutine process_instance_setup_event_data (instance, model, i_core)
    class(process_instance_t), intent(inout), target :: instance
    class(model_data_t), intent(in), optional, target :: model
    integer, intent(in), optional :: i_core
    class(model_data_t), pointer :: current_model
    integer :: i
    class(prc_core_t), pointer :: core => null ()
    if (present (model)) then
       current_model => model
    else
       current_model => instance%process%get_model_ptr ()
    end if
    do i = 1, size (instance%term)
       associate (term => instance%term(i), kin => instance%kin(i))
         if (associated (term%config)) then
            core => instance%process%get_core_term (i)
            call term%setup_event_data (kin, core, current_model)
         end if
       end associate
    end do
    core => null ()
  end subroutine process_instance_setup_event_data

  module subroutine process_instance_choose_mci (instance, i_mci)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_mci
    instance%i_mci = i_mci
    call instance%reset ()
  end subroutine process_instance_choose_mci

  module subroutine process_instance_set_mcpar (instance, x, warmup_flag)
    class(process_instance_t), intent(inout) :: instance
    real(default), dimension(:), intent(in) :: x
    logical, intent(in), optional :: warmup_flag
    logical :: activate
    activate = .true.; if (present (warmup_flag)) activate = .not. warmup_flag
    if (instance%evaluation_status == STAT_INITIAL) then
       associate (mci_work => instance%mci_work(instance%i_mci))
          call mci_work%set (x)
       end associate
       if (activate) call instance%activate ()
    end if
  end subroutine process_instance_set_mcpar

  module subroutine process_instance_receive_beam_momenta (instance)
    class(process_instance_t), intent(inout) :: instance
    if (instance%evaluation_status >= STAT_INITIAL) then
       call instance%sf_chain%receive_beam_momenta ()
       instance%evaluation_status = STAT_BEAM_MOMENTA
    end if
  end subroutine process_instance_receive_beam_momenta

  module subroutine process_instance_set_beam_momenta (instance, p)
    class(process_instance_t), intent(inout) :: instance
    type(vector4_t), dimension(:), intent(in) :: p
    if (instance%evaluation_status >= STAT_INITIAL) then
       call instance%sf_chain%set_beam_momenta (p)
       instance%evaluation_status = STAT_BEAM_MOMENTA
    end if
  end subroutine process_instance_set_beam_momenta

  module subroutine process_instance_recover_beam_momenta (instance, i_term)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_term
    if (.not. instance%process%lab_is_cm ()) then
       if (instance%evaluation_status >= STAT_EFF_KINEMATICS) then
          call instance%kin(i_term)%return_beam_momenta ()
       end if
    end if
  end subroutine process_instance_recover_beam_momenta

  module subroutine process_instance_select_channel (instance, channel)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: channel
    instance%selected_channel = channel
  end subroutine process_instance_select_channel

  module subroutine process_instance_compute_seed_kinematics &
       (instance, recover, skip_term)
    class(process_instance_t), intent(inout) :: instance
    logical, intent(in), optional :: recover
    integer, intent(in), optional :: skip_term
    integer :: channel, skip_component, i, j
    logical :: success
    integer, dimension(:), allocatable :: i_term
    channel = instance%selected_channel
    if (channel == 0) then
       call msg_bug ("Compute seed kinematics: undefined integration channel")
    end if
    if (present (skip_term)) then
       skip_component = instance%term(skip_term)%config%i_component
    else
       skip_component = 0
    end if
    if (present (recover)) then
       if (recover) return
    end if
    if (instance%evaluation_status >= STAT_ACTIVATED) then
       success = .true.
       do i = 1, instance%process%get_n_components ()
          if (i == skip_component)  cycle
          if (instance%process%component_is_selected (i)) then
             allocate (i_term (size (instance%process%get_component_i_terms (i))))
             i_term = instance%process%get_component_i_terms (i)
             do j = 1, size (i_term)
                associate (term => instance%term(i_term(j)), kin => instance%kin(i_term(j)))
                  if (kin%new_seed) then
                     call term%compute_seed_kinematics (kin, &
                          instance%mci_work(instance%i_mci), channel, success)
                     call instance%transfer_same_kinematics (i_term(j))
                  end if
                  if (.not. success)  exit
                  select type (pcm => instance%pcm)
                  class is (pcm_nlo_t)
                     call term%evaluate_projections (kin)
                     call kin%evaluate_radiation_kinematics &
                          (instance%mci_work(instance%i_mci)%get_x_process ())
                     call kin%generate_fsr_in ()
                     call kin%compute_xi_ref_momenta (pcm%region_data, term%nlo_type)
                  end select
                end associate
             end do
          end if
          if (allocated (i_term)) deallocate (i_term)
       end do
       if (success) then
          instance%evaluation_status = STAT_SEED_KINEMATICS
       else
          instance%evaluation_status = STAT_FAILED_KINEMATICS
       end if
    end if
    associate (mci_work => instance%mci_work(instance%i_mci))
       select type (pcm_work => instance%pcm_work)
       class is (pcm_nlo_workspace_t)
          call pcm_work%set_x_rad (mci_work%get_x_process ())
       end select
    end associate
  end subroutine process_instance_compute_seed_kinematics

  pure module function process_instance_get_x_process (instance) result (x)
    real(default), dimension(:), allocatable :: x
    class(process_instance_t), intent(in) :: instance
    allocate (x(size (instance%mci_work(instance%i_mci)%get_x_process ())))
    x = instance%mci_work(instance%i_mci)%get_x_process ()
  end function process_instance_get_x_process

  pure module function process_instance_get_active_component_type &
       (instance) result (nlo_type)
    integer :: nlo_type
    class(process_instance_t), intent(in) :: instance
    nlo_type = instance%process%get_component_nlo_type (instance%i_mci)
  end function process_instance_get_active_component_type

  module subroutine process_instance_recover_mcpar (instance, i_term)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_term
    integer :: channel, i
    if (instance%evaluation_status >= STAT_EFF_KINEMATICS) then
       channel = instance%selected_channel
       if (channel == 0) then
          call msg_bug ("Recover MC parameters: undefined integration channel")
       end if
       call instance%kin(i_term)%recover_mcpar &
            (instance%mci_work(instance%i_mci), channel, instance%term(i_term)%p_seed)
       if (instance%term(i_term)%nlo_type == NLO_REAL) then
          do i = 1, size (instance%term)
             if (i /= i_term .and. instance%term(i)%nlo_type == NLO_REAL) then
                if (instance%term(i)%active) then
                   call instance%kin(i)%recover_mcpar &
                        (instance%mci_work(instance%i_mci), channel, &
                         instance%term(i)%p_seed)
                end if
             end if
          end do
       end if
    end if
  end subroutine process_instance_recover_mcpar

  module subroutine process_instance_recover_sfchain (instance, i_term)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_term
    integer :: channel
    if (instance%evaluation_status >= STAT_EFF_KINEMATICS) then
       channel = instance%selected_channel
       if (channel == 0) then
          call msg_bug ("Recover sfchain: undefined integration channel")
       end if
       call instance%kin(i_term)%recover_sfchain &
            (channel, instance%term(i_term)%p_seed)
    end if
  end subroutine process_instance_recover_sfchain

  module subroutine process_instance_compute_hard_kinematics &
       (instance, recover, skip_term)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in), optional :: skip_term
    logical, intent(in), optional :: recover
    integer :: i
    logical :: success
    success = .true.
    if (instance%evaluation_status >= STAT_SEED_KINEMATICS) then
       do i = 1, size (instance%term)
          associate (term => instance%term(i), kin => instance%kin(i))
            if (term%active) then
               call term%compute_hard_kinematics &
                    (kin, recover, skip_term, success)
               if (.not. success) exit
               !!! Ren scale is zero when this is commented out! Understand!
               if (term%nlo_type == NLO_REAL) &
                  call kin%redo_sf_chain (instance%mci_work(instance%i_mci), &
                      instance%selected_channel)
            end if
          end associate
       end do
       if (success) then
          instance%evaluation_status = STAT_HARD_KINEMATICS
       else
          instance%evaluation_status = STAT_FAILED_KINEMATICS
       end if
    end if
  end subroutine process_instance_compute_hard_kinematics

  module subroutine process_instance_recover_seed_kinematics (instance, i_term)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_term
    type(vector4_t), dimension(:), allocatable :: p_seed_ref
    integer :: i
    if (instance%evaluation_status >= STAT_EFF_KINEMATICS) then
       call instance%term(i_term)%recover_seed_kinematics (instance%kin(i_term))
       if (instance%term(i_term)%nlo_type == NLO_REAL) then
          allocate (p_seed_ref &
               (instance%term(i_term)%isolated%int_eff%get_n_out ()))
          p_seed_ref = instance%term(i_term)%isolated%int_eff%get_momenta &
               (outgoing = .true.)
          do i = 1, size (instance%term)
             if (i /= i_term .and. instance%term(i)%nlo_type == NLO_REAL) then
                if (instance%term(i)%active) then
                   call instance%term(i)%recover_seed_kinematics &
                        (instance%kin(i), p_seed_ref)
                end if
             end if
          end do
       end if
    end if
  end subroutine process_instance_recover_seed_kinematics

  module subroutine process_instance_compute_eff_kinematics &
       (instance, skip_term)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in), optional :: skip_term
    integer :: i
    if (instance%evaluation_status >= STAT_HARD_KINEMATICS) then
       do i = 1, size (instance%term)
          if (present (skip_term)) then
             if (i == skip_term)  cycle
          end if
          if (instance%term(i)%active) then
             call instance%term(i)%compute_eff_kinematics ()
          end if
       end do
       instance%evaluation_status = STAT_EFF_KINEMATICS
    end if
  end subroutine process_instance_compute_eff_kinematics

  module subroutine process_instance_recover_hard_kinematics (instance, i_term)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_term
    integer :: i
    if (instance%evaluation_status >= STAT_EFF_KINEMATICS) then
       call instance%term(i_term)%recover_hard_kinematics ()
       do i = 1, size (instance%term)
          if (i /= i_term) then
             if (instance%term(i)%active) then
                call instance%term(i)%compute_eff_kinematics ()
             end if
          end if
       end do
       instance%evaluation_status = STAT_EFF_KINEMATICS
    end if
  end subroutine process_instance_recover_hard_kinematics

  module subroutine process_instance_evaluate_expressions &
       (instance, scale_forced)
    class(process_instance_t), intent(inout) :: instance
    real(default), intent(in), allocatable, optional :: scale_forced
    integer :: i
    logical :: passed_real
    if (instance%evaluation_status >= STAT_EFF_KINEMATICS) then
       do i = 1, size (instance%term)
          if (instance%term(i)%active) then
             call instance%term(i)%evaluate_expressions &
                  (instance%process%get_beam_config (), scale_forced)
          end if
       end do
       call evaluate_real_scales_and_cuts ()
       call set_ellis_sexton_scale ()
       if (.not. passed_real) then
          instance%evaluation_status = STAT_FAILED_CUTS
       else
          if (any (instance%term%passed)) then
             instance%evaluation_status = STAT_PASSED_CUTS
          else
             instance%evaluation_status = STAT_FAILED_CUTS
          end if
       end if
    end if
  contains
    subroutine evaluate_real_scales_and_cuts ()
      integer :: i
      passed_real = .true.
      select type (pcm => instance%pcm)
      type is (pcm_nlo_t)
         do i = 1, size (instance%term)
            if (instance%term(i)%active .and. instance%term(i)%nlo_type == NLO_REAL) then
               if (pcm%settings%cut_all_real_sqmes) &
                    passed_real = passed_real .and. instance%term(i)%passed
               if (pcm%settings%use_born_scale) &
                    call replace_scales (instance%term(i))
            end if
         end do
      end select
    end subroutine evaluate_real_scales_and_cuts

    subroutine replace_scales (this_term)
      type(term_instance_t), intent(inout) :: this_term
      integer :: i_sub
      i_sub = this_term%config%i_sub
      if (this_term%config%i_term_global /= i_sub .and. i_sub > 0) then
         this_term%ren_scale = instance%term(i_sub)%ren_scale
         this_term%fac_scale = instance%term(i_sub)%fac_scale
      end if
    end subroutine replace_scales

    subroutine set_ellis_sexton_scale ()
      real(default) :: es_scale
      type(var_list_t), pointer :: var_list
      integer :: i
      var_list => instance%process%get_var_list_ptr ()
      es_scale = var_list%get_rval (var_str ("ellis_sexton_scale"))
      do i = 1, size (instance%term)
         if (instance%term(i)%active .and. instance%term(i)%nlo_type == NLO_VIRTUAL) then
            if (es_scale > zero) then
               if (allocated (instance%term(i)%es_scale)) then
                  instance%term(i)%es_scale = es_scale
               else
                  allocate (instance%term(i)%es_scale, source=es_scale)
               end if
            end if
         end if
      end do
    end subroutine set_ellis_sexton_scale
  end subroutine process_instance_evaluate_expressions

  module subroutine process_instance_compute_other_channels &
       (instance, skip_term)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in), optional :: skip_term
    integer :: channel, skip_component, i, j
    integer, dimension(:), allocatable :: i_term
    channel = instance%selected_channel
    if (channel == 0) then
       call msg_bug ("Compute other channels: undefined integration channel")
    end if
    if (present (skip_term)) then
       skip_component = instance%term(skip_term)%config%i_component
    else
       skip_component = 0
    end if
    if (instance%evaluation_status >= STAT_PASSED_CUTS) then
       do i = 1, instance%process%get_n_components ()
          if (i == skip_component)  cycle
          if (instance%process%component_is_selected (i)) then
             allocate (i_term (size (instance%process%get_component_i_terms (i))))
             i_term = instance%process%get_component_i_terms (i)
             do j = 1, size (i_term)
                call instance%kin(i_term(j))%compute_other_channels &
                     (instance%mci_work(instance%i_mci), channel)
             end do
          end if
          if (allocated (i_term)) deallocate (i_term)
       end do
    end if
  end subroutine process_instance_compute_other_channels

  module subroutine process_instance_reset_core_kinematics (instance)
    class(process_instance_t), intent(inout) :: instance
    integer :: i
    if (instance%evaluation_status >= STAT_PASSED_CUTS) then
       do i = 1, size (instance%term)
          associate (term => instance%term(i))
            if (term%active .and. term%passed) then
               if (allocated (term%core_state)) &
                    call term%core_state%reset_new_kinematics ()
            end if
          end associate
       end do
    end if
  end subroutine process_instance_reset_core_kinematics

  module subroutine process_instance_evaluate_trace (instance, recover)
    class(process_instance_t), intent(inout) :: instance
    logical, intent(in), optional :: recover
    class(prc_core_t), pointer :: core => null ()
    integer :: i, i_real_fin, i_core, i_qn, i_flv
    real(default) :: alpha_s, alpha_qed, pt
    class(prc_core_t), pointer :: core_sub => null ()
    class(model_data_t), pointer :: model => null ()
    logical :: has_pdfs
    if (debug_on) call msg_debug2 (D_PROCESS_INTEGRATION, "process_instance_evaluate_trace")
    has_pdfs = instance%process%pcm_contains_pdfs ()
    instance%sqme = zero
    select type (pcm_work => instance%pcm_work)
    type is (pcm_nlo_workspace_t)
       if (allocated(pcm_work%real_sub%sqme_real_arr)) then
          pcm_work%real_sub%sqme_real_arr = zero
       end if
    end select
    call instance%reset_matrix_elements ()
    if (instance%evaluation_status >= STAT_PASSED_CUTS) then
       do i = 1, size (instance%term)
          associate (term => instance%term(i), kin => instance%kin(i))
            if (term%active .and. term%passed) then
               core => instance%process%get_core_term (i)
               select type (pcm => instance%process%get_pcm_ptr ())
               class is (pcm_nlo_t)
                  i_core = pcm%get_i_core (pcm%i_sub)
                  core_sub => instance%process%get_core_ptr (i_core)
               end select
               call term%evaluate_interaction (core, kin)
               call term%evaluate_trace (kin)
               i_real_fin = instance%process%get_associated_real_fin (1)
               if (instance%process%uses_real_partition ()) &
                    call term%apply_real_partition (kin)
               if (term%config%i_component == i_real_fin) then
                  if (term%nlo_type == NLO_REAL .and. .not. term%is_subtraction ()) then
                     !!! Force the scale pT into the events for the real finite
                     associate (p_hard => term%p_hard)
                        !!! This is only the correct pt for ISR
                        pt = transverse_part(p_hard(size(p_hard)))
                        call term%set_fac_scale (pt)
                        select type (core)
                        class is (prc_external_t)
                           select type (core_state => term%core_state)
                           class is (prc_external_state_t)
                              core_state%alpha_qcd = core%qcd%alpha%get (pt)
                           end select
                        type is (prc_omega_t)
                           select type (core_state => term%core_state)
                           type is (omega_state_t)
                              core_state%alpha_qcd = core%qcd%alpha%get (pt)
                           end select
                        end select
                     end associate
                  end if
               else
                  if (term%nlo_type == BORN) then
                     do i_flv = 1, term%connected%trace%get_qn_index_n_flv ()
                        i_qn = term%connected%trace%get_qn_index (i_flv, i_sub = 0)
                        if (.not. term%passed_array(i_flv)) then
                           call term%connected%trace%set_matrix_element &
                                (i_qn, cmplx (zero, zero, default))
                        end if
                     end do
                  end if
                  if ((term%nlo_type == NLO_REAL .and. kin%emitter < 0) &
                       .or. term%nlo_type == NLO_MISMATCH &
                       .or. term%nlo_type == NLO_DGLAP) &
                       call term%set_born_sqmes (core)
                  if (term%is_subtraction () .or. &
                       term%nlo_type == NLO_DGLAP) &
                       call term%set_sf_factors (kin, has_pdfs)
                  if (term%nlo_type > BORN) then
                     if (.not. (term%nlo_type == NLO_REAL .and. &
                          kin%emitter >= 0)) then
                        select type (pcm => term%pcm)
                        type is (pcm_nlo_t)
                           if (char (pcm%settings%nlo_correction_type) == "QCD" .or. &
                                char (pcm%settings%nlo_correction_type) == "Full") &
                                call term%evaluate_color_correlations (core_sub)
                           if (char (pcm%settings%nlo_correction_type) == "EW" .or. &
                                char (pcm%settings%nlo_correction_type) == "Full") then
                              call term%evaluate_charge_correlations (core_sub)
                              select type (pcm => term%pcm)
                              type is (pcm_nlo_t)
                                 associate (reg_data => pcm%region_data)
                                   if (reg_data%alphas_power > 0) &
                                        call term%evaluate_color_correlations (core_sub)
                                 end associate
                              end select
                           end if
                        end select
                     end if
                     if (term%is_subtraction ()) then
                        call term%evaluate_spin_correlations (core_sub)
                     end if
                  end if
                  alpha_s = core%get_alpha_s (term%core_state)
                  alpha_qed = core%get_alpha_qed (term%core_state)
                  if (term%nlo_type > BORN) then
                     select type (pcm => term%pcm)
                     type is (pcm_nlo_t)
                        if (alpha_qed == -1 .and. (&
                             char (pcm%settings%nlo_correction_type) == "EW" .or. &
                             char (pcm%settings%nlo_correction_type) == "Full")) then
                           call msg_bug("Attempting to compute EW corrections with alpha_qed = -1")
                        end if
                     end select
                  end if
                  if (present (recover)) then
                     if (recover)  return
                  end if
                  select case (term%nlo_type)
                  case (NLO_REAL)
                     call term%apply_fks (kin, alpha_s, alpha_qed)
                  case (NLO_VIRTUAL)
                     call term%evaluate_sqme_virt (alpha_s, alpha_qed)
                  case (NLO_MISMATCH)
                     call term%evaluate_sqme_mismatch (alpha_s)
                  case (NLO_DGLAP)
                     call term%evaluate_sqme_dglap (alpha_s, alpha_qed)
                  end select
               end if
            end if
            core_sub => null ()
            instance%sqme = instance%sqme + real (sum (&
                 term%connected%trace%get_matrix_element () * &
                 term%weight))
          end associate
       end do
       core => null ()
       if (instance%pcm_work%is_valid ()) then
          instance%evaluation_status = STAT_EVALUATED_TRACE
       else
          instance%evaluation_status = STAT_FAILED_KINEMATICS
       end if
    else
       !!! Failed kinematics or failed cuts: set sqme to zero
       instance%sqme = zero
    end if
  end subroutine process_instance_evaluate_trace

  module subroutine term_instance_set_born_sqmes (term, core)
    class(term_instance_t), intent(inout) :: term
    class(prc_core_t), intent(in) :: core
    integer :: i_flv, ii_flv
    real(default) :: sqme
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       do i_flv = 1, term%connected%trace%get_qn_index_n_flv ()
          ii_flv = term%connected%trace%get_qn_index (i_flv, i_sub = 0)
          if (term%passed_array (i_flv) .or. .not. term%passed) then
             sqme = real (term%connected%trace%get_matrix_element (ii_flv))
          else
             sqme = zero
          end if
          select case (term%nlo_type)
          case (NLO_REAL)
             pcm_work%real_sub%sqme_born(i_flv) = sqme
          case (NLO_MISMATCH)
             pcm_work%soft_mismatch%sqme_born(i_flv) = sqme
          case (NLO_DGLAP)
             pcm_work%dglap_remnant%sqme_born(i_flv) = sqme
          end select
       end do
    end select
  end subroutine term_instance_set_born_sqmes

  module subroutine term_instance_set_sf_factors (term, kin, has_pdfs)
    class(term_instance_t), intent(inout) :: term
    type(kinematics_t), intent(inout) :: kin
    logical, intent(in) :: has_pdfs
    type(interaction_t), pointer :: sf_chain_int
    real(default) :: factor_born, factor_real
    integer :: n_in, alr, em
    integer :: i_born, i_real
    select type (pcm_work => term%pcm_work)
    type is (pcm_nlo_workspace_t)
       if (.not. has_pdfs) then
          pcm_work%real_sub%sf_factors = one
          return
       end if
       select type (pcm => term%pcm)
       type is (pcm_nlo_t)
          sf_chain_int => kin%sf_chain%get_out_int_ptr ()
          associate (reg_data => pcm%region_data)
             n_in = reg_data%get_n_in ()
             do alr = 1, reg_data%n_regions
                em = reg_data%regions(alr)%emitter
                if (em <= n_in) then
                   i_born = reg_data%regions(alr)%uborn_index
                   i_real = reg_data%regions(alr)%real_index
                   factor_born = sf_chain_int%get_matrix_element &
                        (sf_chain_int%get_sf_qn_index_born (i_born, i_sub = 0))
                   factor_real = sf_chain_int%get_matrix_element &
                        (sf_chain_int%get_sf_qn_index_real (i_real, i_sub = em))
                   call set_factor (pcm_work, alr, em, factor_born, factor_real)
                   if (em == 0) then
                      do em = 1, 2
                         factor_real = sf_chain_int%get_matrix_element &
                              (sf_chain_int%get_sf_qn_index_real (i_real, i_sub = em))
                         call set_factor (pcm_work, alr, em, factor_born, factor_real)
                      end do
                   else
                      factor_real = sf_chain_int%get_matrix_element &
                           (sf_chain_int%get_sf_qn_index_real (i_real, i_sub = 0))
                      call set_factor (pcm_work, alr, 0, factor_born, factor_real)
                   end if
                end if
             end do
          end associate
       end select
    end select
  contains
    subroutine set_factor (pcm_work, alr, em, factor_born, factor_real)
      type(pcm_nlo_workspace_t), intent(inout), target :: pcm_work
      integer, intent(in) :: alr, em
      real(default), intent(in) :: factor_born, factor_real
      real(default) :: factor
      if (any (vanishes ([factor_real, factor_born], tiny(1._default), tiny(1._default)))) then
         factor = zero
      else
         factor = factor_real / factor_born
      end if
      select case (term%nlo_type)
      case (NLO_REAL)
         pcm_work%real_sub%sf_factors(alr, em) = factor
      case (NLO_DGLAP)
         pcm_work%dglap_remnant%sf_factors(alr, em) = factor
      end select
    end subroutine
  end subroutine term_instance_set_sf_factors

  module subroutine process_instance_apply_real_partition (instance)
    class(process_instance_t), intent(inout) :: instance
    integer :: i_component, i_term
    integer, dimension(:), allocatable :: i_terms
    associate (process => instance%process)
       i_component = process%get_first_real_component ()
       if (process%component_is_selected (i_component) .and. &
              process%get_component_nlo_type (i_component) == NLO_REAL) then
          allocate (i_terms, source=process%get_component_i_terms (i_component))
          do i_term = 1, size (i_terms)
             call instance%term(i_terms(i_term))%apply_real_partition ( &
                  instance%kin(i_terms(i_term)))
          end do
       end if
       if (allocated (i_terms)) deallocate (i_terms)
    end associate
  end subroutine process_instance_apply_real_partition

  module subroutine process_instance_set_i_mci_to_real_component (instance)
    class(process_instance_t), intent(inout) :: instance
    integer :: i_mci, i_component
    type(process_component_t), pointer :: component => null ()
    select type (pcm_work => instance%pcm_work)
    type is (pcm_nlo_workspace_t)
       if (allocated (pcm_work%i_mci_to_real_component)) then
          call msg_warning &
               ("i_mci_to_real_component already allocated - replace it")
          deallocate (pcm_work%i_mci_to_real_component)
       end if
       allocate (pcm_work%i_mci_to_real_component (size (instance%mci_work)))
       do i_mci = 1, size (instance%mci_work)
          do i_component = 1, instance%process%get_n_components ()
             component => instance%process%get_component_ptr (i_component)
             if (component%i_mci /= i_mci) cycle
             select case (component%component_type)
             case (COMP_MASTER, COMP_REAL)
                pcm_work%i_mci_to_real_component (i_mci) = &
                     component%config%get_associated_real ()
             case (COMP_REAL_FIN)
                pcm_work%i_mci_to_real_component (i_mci) = &
                     component%config%get_associated_real_fin ()
             case (COMP_REAL_SING)
                pcm_work%i_mci_to_real_component (i_mci) = &
                     component%config%get_associated_real_sing ()
             end select
          end do
       end do
       component => null ()
    end select
  end subroutine process_instance_set_i_mci_to_real_component

  module subroutine process_instance_evaluate_event_data (instance, weight)
    class(process_instance_t), intent(inout) :: instance
    real(default), intent(in), optional :: weight
    integer :: i
    if (instance%evaluation_status >= STAT_EVALUATED_TRACE) then
       do i = 1, size (instance%term)
          associate (term => instance%term(i))
            if (term%active) then
               call term%evaluate_event_data ()
            end if
          end associate
       end do
       if (present (weight)) then
          instance%weight = weight
       else
          instance%weight = &
               instance%mci_work(instance%i_mci)%mci%get_event_weight ()
          instance%excess = &
               instance%mci_work(instance%i_mci)%mci%get_event_excess ()
       end if
       instance%n_dropped = &
            instance%mci_work(instance%i_mci)%mci%get_n_event_dropped ()
       instance%evaluation_status = STAT_EVENT_COMPLETE
    else
       !!! failed kinematics etc.: set weight to zero
       instance%weight = zero
       !!! Maybe we want to process and keep the event nevertheless
       if (instance%keep_failed_events ()) then
          do i = 1, size (instance%term)
             associate (term => instance%term(i))
               if (term%active) then
                  call term%evaluate_event_data ()
               end if
             end associate
          end do
!           do i = 1, size (instance%term)
!              instance%term(i)%fac_scale = zero
!           end do
          instance%evaluation_status = STAT_EVENT_COMPLETE
       end if
    end if
  end subroutine process_instance_evaluate_event_data

  module subroutine process_instance_compute_sqme_rad (instance, &
       i_term, i_phs, is_subtraction, alpha_s_external, scale_forced)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_term, i_phs
    logical, intent(in) :: is_subtraction
    real(default), intent(in), optional :: alpha_s_external
    real(default), intent(in), allocatable, optional :: scale_forced
    class(prc_core_t), pointer :: core
    integer :: i_real_fin
    logical :: has_pdfs
    has_pdfs = instance%process%pcm_contains_pdfs ()
    select type (pcm_work => instance%pcm_work)
    type is (pcm_nlo_workspace_t)
       if (allocated(pcm_work%real_sub%sqme_real_arr)) then
          pcm_work%real_sub%sqme_real_arr = zero
       end if
    end select
    if (debug_on) call msg_debug2 (D_PROCESS_INTEGRATION, "process_instance_compute_sqme_rad")
    select type (pcm_work => instance%pcm_work)
    type is (pcm_nlo_workspace_t)
       associate (term => instance%term(i_term), kin => instance%kin(i_term))
          core => instance%process%get_core_term (i_term)
          if (is_subtraction) then
             call pcm_work%set_subtraction_event ()
          else
             call pcm_work%set_radiation_event ()
          end if
          call term%int_hard%set_momenta (pcm_work%get_momenta &
               (term%pcm, i_phs = i_phs, born_phsp = is_subtraction))
          if (allocated (term%core_state)) &
               call term%core_state%reset_new_kinematics ()
          if (present (alpha_s_external)) then
             call term%set_alpha_qcd_forced (alpha_s_external)
          end if
          call term%compute_eff_kinematics ()
          call term%evaluate_expressions &
               (instance%process%get_beam_config (), scale_forced)
          call term%evaluate_interaction (core, kin)
          call term%evaluate_trace (kin)
          if (term%is_subtraction ()) then
             call term%set_sf_factors (kin, has_pdfs)
             select type (pcm => instance%pcm)
             type is (pcm_nlo_t)
                if (char (pcm%settings%nlo_correction_type) == "QCD" .or. &
                     char (pcm%settings%nlo_correction_type) == "Full") &
                     call term%evaluate_color_correlations (core)
                if (char (pcm%settings%nlo_correction_type) == "EW" .or. &
                     char (pcm%settings%nlo_correction_type) == "Full") &
                     call term%evaluate_charge_correlations (core)
             end select
             call term%evaluate_spin_correlations (core)
          end if
          i_real_fin = instance%process%get_associated_real_fin (1)
          if (term%config%i_component /= i_real_fin) &
               call term%apply_fks (kin, core%get_alpha_s (term%core_state), &
                                    core%get_alpha_qed (term%core_state))
          if (instance%process%uses_real_partition ()) &
               call instance%apply_real_partition ()
       end associate
    end select
    core => null ()
  end subroutine process_instance_compute_sqme_rad

  module subroutine process_instance_normalize_weight (instance)
    class(process_instance_t), intent(inout) :: instance
    if (.not. vanishes (instance%weight)) then
       instance%weight = sign (1._default, instance%weight)
    end if
  end subroutine process_instance_normalize_weight

  module subroutine process_instance_evaluate_sqme (instance, channel, x)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: channel
    real(default), dimension(:), intent(in) :: x
    call instance%reset ()
    call instance%set_mcpar (x)
    call instance%select_channel (channel)
    call instance%compute_seed_kinematics ()
    call instance%compute_hard_kinematics ()
    call instance%compute_eff_kinematics ()
    call instance%evaluate_expressions ()
    call instance%compute_other_channels ()
    call instance%evaluate_trace ()
  end subroutine process_instance_evaluate_sqme

  module subroutine process_instance_recover &
       (instance, channel, i_term, update_sqme, recover_phs, scale_forced)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: channel
    integer, intent(in) :: i_term
    logical, intent(in) :: update_sqme
    logical, intent(in) :: recover_phs
    real(default), intent(in), allocatable, optional :: scale_forced
    logical :: skip_phs, recover
    call instance%activate ()
    instance%evaluation_status = STAT_EFF_KINEMATICS
    call instance%recover_hard_kinematics (i_term)
    call instance%recover_seed_kinematics (i_term)
    call instance%select_channel (channel)
    recover = instance%pcm_work%is_nlo ()
    if (recover_phs) then
       call instance%recover_mcpar (i_term)
       call instance%recover_beam_momenta (i_term)
       call instance%compute_seed_kinematics &
            (recover = recover, skip_term = i_term)
       call instance%compute_hard_kinematics &
            (recover = recover, skip_term = i_term)
       call instance%compute_eff_kinematics (i_term)
       call instance%compute_other_channels (i_term)
    else
       call instance%recover_sfchain (i_term)
    end if
    call instance%evaluate_expressions (scale_forced)
    if (update_sqme) then
       call instance%reset_core_kinematics ()
       call instance%evaluate_trace (recover)
    end if
  end subroutine process_instance_recover

  module subroutine process_instance_evaluate (sampler, c, x_in, val, x, f)
    class(process_instance_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    call sampler%evaluate_sqme (c, x_in)
    if (sampler%is_valid ()) then
       call sampler%fetch (val, x, f)
    end if
    call sampler%record_call ()
    call sampler%evaluate_after_hook ()
  end subroutine process_instance_evaluate

  module function process_instance_is_valid (sampler) result (valid)
    class(process_instance_t), intent(in) :: sampler
    logical :: valid
    valid = sampler%evaluation_status >= STAT_PASSED_CUTS
  end function process_instance_is_valid

  module subroutine process_instance_append_after_hook (sampler, new_hook)
    class(process_instance_t), intent(inout), target :: sampler
    class(process_instance_hook_t), intent(inout), target :: new_hook
    class(process_instance_hook_t), pointer :: last
    if (associated (new_hook%next)) then
       call msg_bug ("process_instance_append_after_hook: " // &
            "reuse of SAME hook object is forbidden.")
    end if
    if (associated (sampler%hook)) then
       last => sampler%hook
       do while (associated (last%next))
          last => last%next
       end do
       last%next => new_hook
    else
       sampler%hook => new_hook
    end if
  end subroutine process_instance_append_after_hook

  module subroutine process_instance_evaluate_after_hook (sampler)
    class(process_instance_t), intent(in) :: sampler
    class(process_instance_hook_t), pointer :: current
    current => sampler%hook
    do while (associated(current))
       call current%evaluate (sampler)
       current => current%next
    end do
  end subroutine process_instance_evaluate_after_hook

  module subroutine process_instance_rebuild (sampler, c, x_in, val, x, f)
    class(process_instance_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x_in
    real(default), intent(in) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    call msg_bug ("process_instance_rebuild not implemented yet")
    x = 0
    f = 0
  end subroutine process_instance_rebuild

  module subroutine process_instance_fetch (sampler, val, x, f)
    class(process_instance_t), intent(in) :: sampler
    real(default), intent(out) :: val
    real(default), dimension(:,:), intent(out) :: x
    real(default), dimension(:), intent(out) :: f
    integer, dimension(:), allocatable :: i_terms
    integer :: i, i_term_base, cc
    integer :: n_channel

    val = 0
    associate (process => sampler%process)
       FIND_COMPONENT: do i = 1, process%get_n_components ()
         if (sampler%process%component_is_selected (i)) then
            allocate (i_terms (size (process%get_component_i_terms (i))))
            i_terms = process%get_component_i_terms (i)
            i_term_base = i_terms(1)
            associate (k => sampler%kin(i_term_base))
              n_channel = k%n_channel
              do cc = 1, n_channel
                 call k%get_mcpar (cc, x(:,cc))
              end do
              f = k%f
              val = sampler%sqme * k%phs_factor
            end associate
            if (allocated (i_terms)) deallocate (i_terms)
            exit FIND_COMPONENT
         end if
       end do FIND_COMPONENT
    end associate
  end subroutine process_instance_fetch

  module subroutine process_instance_init_simulation (instance, i_mci, &
     safety_factor, keep_failed_events)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_mci
    real(default), intent(in), optional :: safety_factor
    logical, intent(in), optional :: keep_failed_events
    call instance%mci_work(i_mci)%init_simulation &
         (safety_factor, keep_failed_events)
  end subroutine process_instance_init_simulation

  module subroutine process_instance_final_simulation (instance, i_mci)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_mci
    call instance%mci_work(i_mci)%final_simulation ()
  end subroutine process_instance_final_simulation

  module subroutine process_instance_get_mcpar (instance, channel, x)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: channel
    real(default), dimension(:), intent(out) :: x
    integer :: i
    if (instance%evaluation_status >= STAT_SEED_KINEMATICS) then
       do i = 1, size (instance%term)
          if (instance%term(i)%active) then
             call instance%kin(i)%get_mcpar (channel, x)
             return
          end if
       end do
       call msg_bug ("Process instance: get_mcpar: no active channels")
    else
       call msg_bug ("Process instance: get_mcpar: no seed kinematics")
    end if
  end subroutine process_instance_get_mcpar

  module function process_instance_has_evaluated_trace (instance) result (flag)
    class(process_instance_t), intent(in) :: instance
    logical :: flag
    flag = instance%evaluation_status >= STAT_EVALUATED_TRACE
  end function process_instance_has_evaluated_trace

  module function process_instance_is_complete_event (instance) result (flag)
    class(process_instance_t), intent(in) :: instance
    logical :: flag
    flag = instance%evaluation_status >= STAT_EVENT_COMPLETE
  end function process_instance_is_complete_event

  module function process_instance_select_i_term (instance) result (i_term)
    integer :: i_term
    class(process_instance_t), intent(in) :: instance
    integer :: i_mci
    i_mci = instance%i_mci
    i_term = instance%process%select_i_term (i_mci)
  end function process_instance_select_i_term

  module function process_instance_get_beam_int_ptr (instance) result (ptr)
    class(process_instance_t), intent(in), target :: instance
    type(interaction_t), pointer :: ptr
    ptr => instance%sf_chain%get_beam_int_ptr ()
  end function process_instance_get_beam_int_ptr

  module function process_instance_get_trace_int_ptr &
       (instance, i_term) result (ptr)
    class(process_instance_t), intent(in), target :: instance
    integer, intent(in) :: i_term
    type(interaction_t), pointer :: ptr
    ptr => instance%term(i_term)%connected%get_trace_int_ptr ()
  end function process_instance_get_trace_int_ptr

  module function process_instance_get_matrix_int_ptr &
       (instance, i_term) result (ptr)
    class(process_instance_t), intent(in), target :: instance
    integer, intent(in) :: i_term
    type(interaction_t), pointer :: ptr
    ptr => instance%term(i_term)%connected%get_matrix_int_ptr ()
  end function process_instance_get_matrix_int_ptr

  module function process_instance_get_flows_int_ptr &
       (instance, i_term) result (ptr)
    class(process_instance_t), intent(in), target :: instance
    integer, intent(in) :: i_term
    type(interaction_t), pointer :: ptr
    ptr => instance%term(i_term)%connected%get_flows_int_ptr ()
  end function process_instance_get_flows_int_ptr

  module function process_instance_get_state_flv &
       (instance, i_term) result (state_flv)
    class(process_instance_t), intent(in) :: instance
    integer, intent(in) :: i_term
    type(state_flv_content_t) :: state_flv
    state_flv = instance%term(i_term)%connected%get_state_flv ()
  end function process_instance_get_state_flv

  module function process_instance_get_isolated_state_ptr &
       (instance, i_term) result (ptr)
    class(process_instance_t), intent(in), target :: instance
    integer, intent(in) :: i_term
    type(isolated_state_t), pointer :: ptr
    ptr => instance%term(i_term)%isolated
  end function process_instance_get_isolated_state_ptr

  module function process_instance_get_connected_state_ptr &
       (instance, i_term) result (ptr)
    class(process_instance_t), intent(in), target :: instance
    integer, intent(in) :: i_term
    type(connected_state_t), pointer :: ptr
    ptr => instance%term(i_term)%connected
  end function process_instance_get_connected_state_ptr

  module subroutine process_instance_get_beam_index (instance, i_term, i_beam)
    class(process_instance_t), intent(in) :: instance
    integer, intent(in) :: i_term
    integer, dimension(:), intent(out) :: i_beam
    call instance%term(i_term)%connected%get_beam_index (i_beam)
  end subroutine process_instance_get_beam_index

  module subroutine process_instance_get_in_index (instance, i_term, i_in)
    class(process_instance_t), intent(in) :: instance
    integer, intent(in) :: i_term
    integer, dimension(:), intent(out) :: i_in
    call instance%term(i_term)%connected%get_in_index (i_in)
  end subroutine process_instance_get_in_index

  module function process_instance_get_sqme (instance, i_term) result (sqme)
    real(default) :: sqme
    class(process_instance_t), intent(in) :: instance
    integer, intent(in), optional :: i_term
    if (instance%evaluation_status >= STAT_EVALUATED_TRACE) then
       if (present (i_term)) then
          sqme = instance%term(i_term)%connected%trace%get_matrix_element (1)
       else
          sqme = instance%sqme
       end if
    else
       sqme = 0
    end if
  end function process_instance_get_sqme

  module function process_instance_get_weight (instance) result (weight)
    real(default) :: weight
    class(process_instance_t), intent(in) :: instance
    if (instance%evaluation_status >= STAT_EVENT_COMPLETE) then
       weight = instance%weight
    else
       weight = 0
    end if
  end function process_instance_get_weight

  module function process_instance_get_excess (instance) result (excess)
    real(default) :: excess
    class(process_instance_t), intent(in) :: instance
    if (instance%evaluation_status >= STAT_EVENT_COMPLETE) then
       excess = instance%excess
    else
       excess = 0
    end if
  end function process_instance_get_excess

  module function process_instance_get_n_dropped (instance) result (n_dropped)
    integer :: n_dropped
    class(process_instance_t), intent(in) :: instance
    if (instance%evaluation_status >= STAT_EVENT_COMPLETE) then
       n_dropped = instance%n_dropped
    else
       n_dropped = 0
    end if
  end function process_instance_get_n_dropped

  module function process_instance_get_channel (instance) result (channel)
    integer :: channel
    class(process_instance_t), intent(in) :: instance
    channel = instance%selected_channel
  end function process_instance_get_channel

  module subroutine process_instance_set_fac_scale (instance, fac_scale)
    class(process_instance_t), intent(inout) :: instance
    real(default), intent(in) :: fac_scale
    integer :: i_term
    i_term = 1
    call instance%term(i_term)%set_fac_scale (fac_scale)
  end subroutine process_instance_set_fac_scale

  module function process_instance_get_fac_scale &
       (instance, i_term) result (fac_scale)
    class(process_instance_t), intent(in) :: instance
    integer, intent(in) :: i_term
    real(default) :: fac_scale
    fac_scale = instance%term(i_term)%get_fac_scale ()
  end function process_instance_get_fac_scale

  module function process_instance_get_alpha_s &
       (instance, i_term) result (alpha_s)
    real(default) :: alpha_s
    class(process_instance_t), intent(in) :: instance
    integer, intent(in) :: i_term
    class(prc_core_t), pointer :: core => null ()
    core => instance%process%get_core_term (i_term)
    alpha_s = instance%term(i_term)%get_alpha_s (core)
    core => null ()
  end function process_instance_get_alpha_s

  module function process_instance_get_qcd (process_instance) result (qcd)
    type(qcd_t) :: qcd
    class(process_instance_t), intent(in) :: process_instance
    qcd = process_instance%process%get_qcd ()
  end function process_instance_get_qcd

  module subroutine process_instance_reset_counter (process_instance)
    class(process_instance_t), intent(inout) :: process_instance
    call process_instance%mci_work(process_instance%i_mci)%reset_counter ()
  end subroutine process_instance_reset_counter

  module subroutine process_instance_record_call (process_instance)
    class(process_instance_t), intent(inout) :: process_instance
    call process_instance%mci_work(process_instance%i_mci)%record_call &
         (process_instance%evaluation_status)
  end subroutine process_instance_record_call

  pure module function process_instance_get_counter &
       (process_instance) result (counter)
    class(process_instance_t), intent(in) :: process_instance
    type(process_counter_t) :: counter
    counter = process_instance%mci_work(process_instance%i_mci)%get_counter ()
  end function process_instance_get_counter

  pure module function process_instance_get_actual_calls_total &
       (process_instance) result (n)
    class(process_instance_t), intent(in) :: process_instance
    integer :: n
    integer :: i
    type(process_counter_t) :: counter
    n = 0
    do i = 1, size (process_instance%mci_work)
       counter = process_instance%mci_work(i)%get_counter ()
       n = n + counter%total
    end do
  end function process_instance_get_actual_calls_total

  module subroutine process_instance_reset_matrix_elements (instance)
    class(process_instance_t), intent(inout) :: instance
    integer :: i_term
    do i_term = 1, size (instance%term)
       call instance%term(i_term)%connected%trace%set_matrix_element &
            (cmplx (0, 0, default))
       call instance%term(i_term)%connected%matrix%set_matrix_element &
            (cmplx (0, 0, default))
    end do
  end subroutine process_instance_reset_matrix_elements

  module subroutine process_instance_get_test_phase_space_point (instance, &
         i_component, i_core, p)
    type(vector4_t), dimension(:), allocatable, intent(out) :: p
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_component, i_core
    real(default), dimension(:), allocatable :: x
    logical :: success
    integer :: i_term
    instance%i_mci = i_component
    i_term = instance%process%get_i_term (i_core)
    associate (term => instance%term(i_term), kin => instance%kin(i_term))
       allocate (x (instance%mci_work(i_component)%config%n_par))
       x = 0.5_default
       call instance%set_mcpar (x, .true.)
       call instance%select_channel (1)
       call term%compute_seed_kinematics &
            (kin, instance%mci_work(i_component), 1, success)
       call kin%evaluate_radiation_kinematics &
            (instance%mci_work(instance%i_mci)%get_x_process ())
       call term%compute_hard_kinematics (kin, success = success)
       allocate (p (size (term%p_hard)))
       p = term%int_hard%get_momenta ()
    end associate
  end subroutine process_instance_get_test_phase_space_point

  pure module function process_instance_get_p_hard &
       (process_instance, i_term) result (p_hard)
    type(vector4_t), dimension(:), allocatable :: p_hard
    class(process_instance_t), intent(in) :: process_instance
    integer, intent(in) :: i_term
    allocate (p_hard (size (process_instance%term(i_term)%get_p_hard ())))
    p_hard = process_instance%term(i_term)%get_p_hard ()
  end function process_instance_get_p_hard

  module function process_instance_get_first_active_i_term &
       (instance) result (i_term)
    integer :: i_term
    class(process_instance_t), intent(in) :: instance
    integer :: i
    i_term = 0
    do i = 1, size (instance%term)
       if (instance%term(i)%active) then
          i_term = i
          exit
       end if
    end do
  end function process_instance_get_first_active_i_term

  module function process_instance_get_real_of_mci (instance) result (i_real)
    integer :: i_real
    class(process_instance_t), intent(in) :: instance
    select type (pcm_work => instance%pcm_work)
    type is (pcm_nlo_workspace_t)
       i_real = pcm_work%i_mci_to_real_component (instance%i_mci)
    end select
  end function process_instance_get_real_of_mci

  module function process_instance_get_connected_states &
       (instance, i_component) result (connected)
    type(connected_state_t), dimension(:), allocatable :: connected
    class(process_instance_t), intent(in) :: instance
    integer, intent(in) :: i_component
    connected = instance%process%get_connected_states (i_component, &
        instance%term(:)%connected)
  end function process_instance_get_connected_states

  module function process_instance_get_sqrts (instance) result (sqrts)
    class(process_instance_t), intent(in) :: instance
    real(default) :: sqrts
    sqrts = instance%process%get_sqrts ()
  end function process_instance_get_sqrts

  module function process_instance_get_polarization (instance) result (pol)
    class(process_instance_t), intent(in) :: instance
    real(default), dimension(:), allocatable :: pol
    pol = instance%process%get_polarization ()
  end function process_instance_get_polarization

  module function process_instance_get_beam_file (instance) result (file)
    class(process_instance_t), intent(in) :: instance
    type(string_t) :: file
    file = instance%process%get_beam_file ()
  end function process_instance_get_beam_file

  module function process_instance_get_process_name (instance) result (name)
    class(process_instance_t), intent(in) :: instance
    type(string_t) :: name
    name = instance%process%get_id ()
  end function process_instance_get_process_name

  module subroutine process_instance_get_trace &
       (instance, pset, i_term, n_incoming)
    class(process_instance_t), intent(in), target :: instance
    type(particle_set_t), intent(out) :: pset
    integer, intent(in) :: i_term
    integer, intent(in), optional :: n_incoming
    type(interaction_t), pointer :: int
    logical :: ok
    int => instance%get_trace_int_ptr (i_term)
    call pset%init (ok, int, int, FM_IGNORE_HELICITY, &
         [0._default, 0._default], .false., .true., n_incoming)
  end subroutine process_instance_get_trace

  module subroutine process_instance_set_trace &
       (instance, pset, i_term, recover_beams, check_match, success)
    class(process_instance_t), intent(inout), target :: instance
    type(particle_set_t), intent(in) :: pset
    integer, intent(in) :: i_term
    logical, intent(in), optional :: recover_beams, check_match
    logical, intent(out), optional :: success
    type(interaction_t), pointer :: int
    integer :: n_in
    int => instance%get_trace_int_ptr (i_term)
    n_in = instance%process%get_n_in ()
    call pset%fill_interaction (int, n_in, &
         recover_beams = recover_beams, &
         check_match = check_match, &
         state_flv = instance%get_state_flv (i_term), &
         success = success)
  end subroutine process_instance_set_trace

  module subroutine process_instance_set_alpha_qcd_forced &
       (instance, i_term, alpha_qcd)
    class(process_instance_t), intent(inout) :: instance
    integer, intent(in) :: i_term
    real(default), intent(in) :: alpha_qcd
    call instance%term(i_term)%set_alpha_qcd_forced (alpha_qcd)
  end subroutine process_instance_set_alpha_qcd_forced

  module function process_instance_has_nlo_component (instance) result (nlo)
    class(process_instance_t), intent(in) :: instance
    logical :: nlo
    nlo = instance%process%is_nlo_calculation ()
  end function process_instance_has_nlo_component

  module function process_instance_keep_failed_events (instance) result (keep)
    logical :: keep
    class(process_instance_t), intent(in) :: instance
    keep = instance%mci_work(instance%i_mci)%keep_failed_events
  end function process_instance_keep_failed_events

  module function process_instance_get_term_indices &
       (instance, nlo_type) result (i_term)
    integer, dimension(:), allocatable :: i_term
    class(process_instance_t), intent(in) :: instance
    integer :: nlo_type
    allocate (i_term (count (instance%term%nlo_type == nlo_type)))
    i_term = pack (instance%term%get_i_term_global (), &
         instance%term%nlo_type == nlo_type)
  end function process_instance_get_term_indices

  module function process_instance_get_boost_to_lab &
       (instance, i_term) result (lt)
    type(lorentz_transformation_t) :: lt
    class(process_instance_t), intent(in) :: instance
    integer, intent(in) :: i_term
    lt = instance%kin(i_term)%get_boost_to_lab ()
  end function process_instance_get_boost_to_lab

  module function process_instance_get_boost_to_cms &
       (instance, i_term) result (lt)
    type(lorentz_transformation_t) :: lt
    class(process_instance_t), intent(in) :: instance
    integer, intent(in) :: i_term
    lt = instance%kin(i_term)%get_boost_to_cms ()
  end function process_instance_get_boost_to_cms

  module function process_instance_lab_is_cm &
       (instance, i_term) result (lab_is_cm)
    logical :: lab_is_cm
    class(process_instance_t), intent(in) :: instance
    integer, intent(in) :: i_term
    lab_is_cm = instance%kin(i_term)%phs%lab_is_cm ()
  end function process_instance_lab_is_cm

  module subroutine pacify_process_instance (instance)
    type(process_instance_t), intent(inout) :: instance
    integer :: i
    do i = 1, size (instance%kin)
       call pacify (instance%kin(i)%phs)
    end do
  end subroutine pacify_process_instance


end submodule instances_s

