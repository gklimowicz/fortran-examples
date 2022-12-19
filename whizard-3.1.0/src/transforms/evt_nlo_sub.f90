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

submodule (evt_nlo) evt_nlo_s

  use debug_master, only: debug_on
  use io_units, only: given_output_unit
  use diagnostics
  use format_utils, only: write_separator
  use numeric_utils, only: nearly_equal
  use physics_defs, only: BORN, NLO_REAL
  use lorentz
  use interactions, only: interaction_t
  use pcm, only: pcm_nlo_t, pcm_nlo_workspace_t
  use prc_core, only: prc_core_t
  use prc_external, only: prc_external_t
  use phs_fks, only: SQRTS_FIXED, SQRTS_VAR

  implicit none

contains

  module subroutine evt_nlo_write_name (evt, unit)
    class(evt_nlo_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)") "Event transform: NLO"
  end subroutine evt_nlo_write_name

  module subroutine evt_nlo_write (evt, unit, verbose, more_verbose, testflag)
    class(evt_nlo_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, more_verbose, testflag
    integer :: u, i
    u = given_output_unit (unit)
    call write_separator (u, 2)
    call evt%write_name (u)
    call write_separator (u)
    call evt%base_write (u, testflag = testflag, show_set = .true.)
    write (u,'(A,ES16.9)')  "sqme_rad = ", evt%sqme_rad
    write (u, "(3x,A,I0)")  "i_evaluation = ", evt%i_evaluation
    call write_separator (u)
    write (u, "(1x,A)")  "Radiated particle sets:"
    do i = 1, size (evt%particle_set_nlo)
       call evt%particle_set_nlo(i)%write (u, testflag = testflag)
       call write_separator (u)
    end do
  end subroutine evt_nlo_write

  module subroutine evt_nlo_connect &
       (evt, process_instance, model, process_stack)
    class(evt_nlo_t), intent(inout), target :: evt
    type(process_instance_t), intent(in), target :: process_instance
    class(model_data_t), intent(in), target :: model
    type(process_stack_t), intent(in), optional :: process_stack
    if (debug_on) call msg_debug (D_TRANSFORMS, "evt_nlo_connect")
    call evt%base_connect (process_instance, model, process_stack)
    select type (pcm_work => process_instance%pcm_work)
    class is (pcm_nlo_workspace_t)
       select type (pcm => process_instance%pcm)
       type is (pcm_nlo_t)
          call pcm%setup_phs_generator (pcm_work, evt%phs_fks_generator, &
               process_instance%get_sqrts ())
          call evt%set_i_evaluation_mappings (pcm%region_data, &
               pcm_work%real_kinematics%alr_to_i_phs)
       end select
    end select
    call evt%set_mode (process_instance)
    call evt%setup_general_event_kinematics (process_instance)
    if (evt%mode > EVT_NLO_SEPARATE_BORNLIKE) &
         call evt%setup_real_event_kinematics (process_instance)
    if (debug_on) call msg_debug2 (D_TRANSFORMS, "evt_nlo_connect: success")
  end subroutine evt_nlo_connect

  module subroutine evt_nlo_set_i_evaluation_mappings &
       (evt, reg_data, alr_to_i_phs)
    class(evt_nlo_t), intent(inout) :: evt
    type(region_data_t), intent(in) :: reg_data
    integer, intent(in), dimension(:) :: alr_to_i_phs
    integer :: n_phs, alr
    integer :: i_evaluation, i_phs, emitter
    logical :: checked
    type(registered_triple_t), allocatable, target :: check_list
    i_evaluation = 1
    n_phs = reg_data%n_phs
    allocate (evt%i_evaluation_to_i_phs (n_phs), source = 0)
    allocate (evt%i_evaluation_to_emitter (n_phs), source = -1)
    allocate (evt%i_evaluation_to_i_term (0 : n_phs), source = 0)
    do alr = 1, reg_data%n_regions
       i_phs = alr_to_i_phs (alr)
       emitter = reg_data%regions(alr)%emitter
       call search_check_list (checked)
       if (.not. checked) then
          evt%i_evaluation_to_i_phs (i_evaluation) = i_phs
          evt%i_evaluation_to_emitter (i_evaluation) = emitter
          i_evaluation = i_evaluation + 1
       end if
    end do
    call fill_i_evaluation_to_i_term ()
    if (.not. (all (evt%i_evaluation_to_i_phs > 0) &
       .and. all (evt%i_evaluation_to_emitter > -1))) then
       call msg_fatal ("evt_nlo: Inconsistent mappings!")
    else
       if (debug2_active (D_TRANSFORMS)) then
          print *, 'evt_nlo Mappings, i_evaluation -> '
          print *, 'i_phs: ', evt%i_evaluation_to_i_phs
          print *, 'emitter: ', evt%i_evaluation_to_emitter
       end if
    end if
  contains
    subroutine fill_i_evaluation_to_i_term ()
      integer :: i_term, i_evaluation, term_emitter
      !!! First find subtraction component
      i_evaluation = 1
      do i_term = 1, evt%process%get_n_terms ()
         if (evt%process_instance%term(i_term)%nlo_type /= NLO_REAL) cycle
         term_emitter = evt%process_instance%kin(i_term)%emitter
         if (term_emitter < 0) then
            evt%i_evaluation_to_i_term (0) = i_term
         else if (evt%i_evaluation_to_emitter(i_evaluation) == term_emitter) then
            evt%i_evaluation_to_i_term (i_evaluation) = i_term
            i_evaluation = i_evaluation + 1
         end if
      end do
    end subroutine fill_i_evaluation_to_i_term

    subroutine search_check_list (found)
      logical, intent(out) :: found
      type(registered_triple_t), pointer :: current_triple => null ()
      if (allocated (check_list)) then
         current_triple => check_list
         do
            if (all (current_triple%phs_em == [i_phs, emitter])) then
               found = .true.
               exit
            end if
            if (.not. associated (current_triple%next)) then
               allocate (current_triple%next)
               current_triple%next%phs_em = [i_phs, emitter]
               found = .false.
               exit
            else
               current_triple => current_triple%next
            end if
         end do
      else
         allocate (check_list)
         check_list%phs_em = [i_phs, emitter]
         found = .false.
      end if
    end subroutine search_check_list
  end subroutine evt_nlo_set_i_evaluation_mappings

  module function evt_nlo_get_i_phs (evt) result (i_phs)
    integer :: i_phs
    class(evt_nlo_t), intent(in) :: evt
    i_phs = evt%i_evaluation_to_i_phs (evt%i_evaluation)
  end function evt_nlo_get_i_phs

  module function evt_nlo_get_emitter (evt) result (emitter)
    integer :: emitter
    class(evt_nlo_t), intent(in) :: evt
    emitter = evt%i_evaluation_to_emitter (evt%i_evaluation)
  end function evt_nlo_get_emitter

  module function evt_nlo_get_i_term (evt) result (i_term)
    integer :: i_term
    class(evt_nlo_t), intent(in) :: evt
    if (evt%mode >= EVT_NLO_SEPARATE_REAL) then
       i_term = evt%i_evaluation_to_i_term (evt%i_evaluation)
    else
       i_term = evt%process_instance%get_first_active_i_term ()
    end if
  end function evt_nlo_get_i_term

  module subroutine evt_nlo_generate_weighted (evt, probability)
    class(evt_nlo_t), intent(inout) :: evt
    real(default), intent(inout) :: probability
    real(default) :: sqme
    call print_debug_info ()
    sqme = probability
    if (evt%mode > EVT_NLO_SEPARATE_BORNLIKE) then
       if (evt%i_evaluation == 0) then
          call evt%reset_phs_identifiers ()
          call evt%evaluate_real_kinematics ()
          if (evt%mode == EVT_NLO_SEPARATE_REAL) then
             sqme = evt%compute_subtraction_sqmes ()
          else
             sqme = sqme - evt%compute_all_sqme_rad ()
          end if
       else
          call evt%compute_real ()
          sqme = evt%sqme_rad
       end if
    end if
    probability = sqme
    if (debug_on)  call msg_debug &
         (D_TRANSFORMS, "probability (after)", probability)
  contains
    function status_code_to_string (mode) result (smode)
      type(string_t) :: smode
      integer, intent(in) :: mode
      select case (mode)
      case (EVT_NLO_UNDEFINED)
         smode = var_str ("Undefined")
      case (EVT_NLO_SEPARATE_BORNLIKE)
         smode = var_str ("Born-like")
      case (EVT_NLO_SEPARATE_REAL)
         smode = var_str ("Real")
      case (EVT_NLO_COMBINED)
         smode = var_str ("Combined")
      end select
    end function status_code_to_string

    subroutine print_debug_info ()
       if (debug_on)  call msg_debug (D_TRANSFORMS, "evt_nlo_generate_weighted")
       if (debug_on)  call msg_debug &
            (D_TRANSFORMS, char ("mode: " // status_code_to_string (evt%mode)))
       if (debug_on)  call msg_debug &
            (D_TRANSFORMS, "probability (before)", probability)
       if (debug_on)  call msg_debug &
            (D_TRANSFORMS, "evt%i_evaluation", evt%i_evaluation)
       if (debug2_active (D_TRANSFORMS)) then
          if (evt%mode > EVT_NLO_SEPARATE_BORNLIKE) then
             if (evt%i_evaluation == 0) then
                print *, 'Evaluate subtraction component'
             else
                print *, 'Evaluate radiation component'
             end if
          end if
       end if
    end subroutine print_debug_info
  end subroutine evt_nlo_generate_weighted

  module subroutine evt_nlo_reset_phs_identifiers (evt)
     class(evt_nlo_t), intent(inout) :: evt
     evt%event_deps%phs_identifiers%evaluated = .false.
  end subroutine evt_nlo_reset_phs_identifiers

  module subroutine evt_nlo_connected_set_real_IS_momenta (evt)
    class(evt_nlo_t), intent(inout) :: evt
    type(vector4_t) :: p_hard, p_beam, p_remn
    type(interaction_t), pointer :: int_matrix
    integer :: i, i_term, n_in, i_in_beam, i_in_hard, i_in_remn
    i_term = evt%get_i_term ()
    int_matrix => evt%process_instance%get_matrix_int_ptr (i_term)
    n_in = evt%particle_set%get_n_in ()
    do i = 1, n_in
       i_in_beam = i
       i_in_hard = n_in + i
       i_in_remn = 2 * n_in + i
       p_hard = evt%process_instance%term(i_term)%int_hard%get_momentum (i)
       p_beam = int_matrix%get_momentum (i_in_beam)
       p_remn = p_beam - p_hard
       call int_matrix%set_momentum (p_hard , i_in_hard)
       call int_matrix%set_momentum (p_remn , i_in_remn)
    end do
  end subroutine evt_nlo_connected_set_real_IS_momenta

  module subroutine evt_nlo_make_particle_set &
       (evt, factorization_mode, keep_correlations, r)
    class(evt_nlo_t), intent(inout) :: evt
    integer, intent(in) :: factorization_mode
    logical, intent(in) :: keep_correlations
    real(default), dimension(:), intent(in), optional :: r
    if (evt%mode >= EVT_NLO_SEPARATE_BORNLIKE) then
       call make_factorized_particle_set (evt, factorization_mode, &
            keep_correlations, r, evt%get_i_term (), &
            evt%get_selected_quantum_numbers (evt%selected_i_flv))
    else
       call make_factorized_particle_set (evt, factorization_mode, &
            keep_correlations, r)
    end if
  end subroutine evt_nlo_make_particle_set

  module subroutine evt_nlo_evaluate_real_kinematics (evt)
    class(evt_nlo_t), intent(inout) :: evt
    integer :: alr, i_phs, i_con, emitter
    real(default), dimension(3) :: x_rad
    logical :: use_contributors
    integer :: n_regions
    integer :: i_term
    type(vector4_t), dimension(:), allocatable :: p_real

    select type (pcm_work => evt%process_instance%pcm_work)
    class is (pcm_nlo_workspace_t)
       x_rad = pcm_work%real_kinematics%x_rad
       associate (event_deps => evt%event_deps)
          i_term = evt%get_i_term ()
          event_deps%p_born_lab%phs_point(1) = &
               evt%process_instance%term(i_term)%p_seed
          event_deps%p_born_cms%phs_point(1) &
               = evt%boost_to_cms (event_deps%p_born_lab%phs_point(1))
          call evt%phs_fks_generator%set_sqrts_hat &
               (event_deps%p_born_cms%get_energy (1, 1))
          use_contributors = allocated (event_deps%contributors)
          select type (pcm => evt%process_instance%pcm)
          type is (pcm_nlo_t)
             n_regions = pcm%region_data%n_regions
          end select
          do alr = 1, n_regions
             i_phs = pcm_work%real_kinematics%alr_to_i_phs(alr)
             if (event_deps%phs_identifiers(i_phs)%evaluated) cycle
             emitter = event_deps%phs_identifiers(i_phs)%emitter
             associate (generator => evt%phs_fks_generator)
                if (emitter <= evt%process%get_n_in ()) then
                   call generator%prepare_generation (x_rad, i_phs, emitter, &
                        event_deps%p_born_cms%phs_point(1)%get (), &
                        event_deps%phs_identifiers)
                   ! TODO wk 19-02-28: intent of p_real (also below)?
                   p_real = event_deps%p_real_lab%phs_point(i_phs)
                   select case (generator%isr_kinematics%isr_mode)
                   case (SQRTS_FIXED)
                      call generator%generate_isr_fixed_beam_energy (i_phs, &
                           event_deps%p_born_cms%phs_point(1)%get (), &
                           p_real)
                   case (SQRTS_VAR)
                      call generator%generate_isr (i_phs, &
                           event_deps%p_born_lab%phs_point(1)%get (), &
                           p_real)
                   end select
                   event_deps%p_real_lab%phs_point(i_phs) = p_real
                   event_deps%p_real_cms%phs_point(i_phs) &
                        = evt%boost_to_cms (event_deps%p_real_lab%phs_point(i_phs))
                else
                   if (use_contributors) then
                      i_con = event_deps%alr_to_i_con(alr)
                      call generator%prepare_generation (x_rad, i_phs, emitter, &
                           event_deps%p_born_cms%phs_point(1)%get (), &
                           event_deps%phs_identifiers, event_deps%contributors, i_con)
                      p_real = event_deps%p_real_cms%phs_point(i_phs)
                      call generator%generate_fsr (emitter, i_phs, i_con, &
                           event_deps%p_born_cms%phs_point(1)%get (), &
                           p_real)
                      event_deps%p_real_cms%phs_point(i_phs) = p_real
                   else
                      call generator%prepare_generation (x_rad, i_phs, emitter, &
                           event_deps%p_born_cms%phs_point(1)%get (), &
                           event_deps%phs_identifiers)
                      p_real = event_deps%p_real_cms%phs_point(i_phs)
                      call generator%generate_fsr (emitter, i_phs, &
                           event_deps%p_born_cms%phs_point(1)%get (), &
                           p_real)
                      event_deps%p_real_cms%phs_point(i_phs) = p_real
                   end if
                   event_deps%p_real_lab%phs_point(i_phs) &
                        = evt%boost_to_lab (event_deps%p_real_cms%phs_point(i_phs))
                end if
             end associate
             call pcm_work%set_momenta &
                  (event_deps%p_born_lab%phs_point(1)%get (), &
                  event_deps%p_real_lab%phs_point(i_phs)%get (), &
                  i_phs)
             call pcm_work%set_momenta &
                  (event_deps%p_born_cms%phs_point(1)%get (), &
                  event_deps%p_real_cms%phs_point(i_phs)%get (), &
                  i_phs, cms = .true.)
             event_deps%phs_identifiers(i_phs)%evaluated = .true.
          end do
       end associate
    end select
  end subroutine evt_nlo_evaluate_real_kinematics

  module function evt_nlo_compute_subtraction_sqmes (evt) result (sqme)
    class(evt_nlo_t), intent(inout) :: evt
    real(default) :: sqme
    integer :: i_phs, i_term
    if (debug_on)  call msg_debug &
         (D_TRANSFORMS, "evt_nlo_compute_subtraction_sqmes")
    sqme = zero
    associate (event_deps => evt%event_deps)
      i_phs = 1; i_term = evt%i_evaluation_to_i_term(0)
      call evt%process_instance%compute_sqme_rad &
           (i_term, i_phs, is_subtraction = .true.)
      sqme = sqme + evt%process_instance%get_sqme (i_term)
    end associate
  end function evt_nlo_compute_subtraction_sqmes

  module subroutine evt_nlo_compute_real (evt)
    class(evt_nlo_t), intent(inout) :: evt
    integer :: i_phs, i_term
    if (debug_on) call msg_debug (D_TRANSFORMS, "evt_nlo_compute_real")
    i_phs = evt%get_i_phs ()
    i_term = evt%i_evaluation_to_i_term (evt%i_evaluation)
    associate (event_deps => evt%event_deps)
      call evt%process_instance%compute_sqme_rad (i_term, i_phs, &
           is_subtraction = .false.)
      evt%sqme_rad = evt%process_instance%get_sqme (i_term)
    end associate
  end subroutine evt_nlo_compute_real

  module function evt_nlo_compute_all_sqme_rad (evt) result (sqme)
    class(evt_nlo_t), intent(inout) :: evt
    real(default) :: sqme
    integer :: i_phs, i_term
    if (debug_on) call msg_debug (D_TRANSFORMS, "evt_nlo_compute_all_sqme_rad")
    sqme = zero
    do i_term = 1, size (evt%process_instance%term)
       if (evt%is_valid_event (i_term)) then
          associate (term => evt%process_instance%term(i_term))
            if (term%nlo_type == NLO_REAL .and. &
                 .not. term%is_subtraction ()) then
               i_phs = evt%process_instance%kin(i_term)%i_phs
               call evt%process_instance%compute_sqme_rad ( &
                    i_term, i_phs, is_subtraction = .false.)
               sqme = sqme + evt%process_instance%get_sqme (i_term)
            end if
          end associate
       end if
    end do
  end function evt_nlo_compute_all_sqme_rad

  module function evt_nlo_boost_to_cms (evt, p_lab) result (p_cms)
    type(phs_point_t), intent(in) :: p_lab
    class(evt_nlo_t), intent(in) :: evt
    type(phs_point_t) :: p_cms
    type(vector4_t) :: p0, p1
    type(lorentz_transformation_t) :: lt_lab_to_cms, lt
    real(default) :: sqrts_hat
    integer :: i_boost, n_legs_born
    if (evt%event_deps%lab_is_cm) then
       lt_lab_to_cms = identity
    else
       n_legs_born = size (evt%event_deps%p_born_lab%phs_point(1))
       if (size (p_lab) == n_legs_born) then
          i_boost = evt%get_i_term ()
          lt_lab_to_cms = evt%process_instance%get_boost_to_cms (i_boost)
       else
          sqrts_hat = (p_lab%select (1) + p_lab%select (2))**1
          p0 = p_lab%select (1) + p_lab%select (2)
          lt = boost (p0, sqrts_hat)
          p1 = inverse(lt) * p_lab%select (1)
          lt_lab_to_cms = inverse (lt * rotation_to_2nd (3, space_part (p1)))
       end if
    end if
    p_cms = lt_lab_to_cms * p_lab
  end function evt_nlo_boost_to_cms

  module function evt_nlo_boost_to_lab (evt, p_cms) result (p_lab)
    type(phs_point_t) :: p_lab
    class(evt_nlo_t), intent(in) :: evt
    type(phs_point_t), intent(in) :: p_cms
    type(lorentz_transformation_t) :: lt_cms_to_lab
    integer :: i_boost
    if (evt%event_deps%lab_is_cm) then
       lt_cms_to_lab = identity
    else
       i_boost = evt%get_i_term ()
       lt_cms_to_lab = evt%process_instance%get_boost_to_lab (i_boost)
    end if
    p_lab = lt_cms_to_lab * p_cms
  end function evt_nlo_boost_to_lab

  module subroutine evt_nlo_setup_general_event_kinematics &
       (evt, process_instance)
    class(evt_nlo_t), intent(inout) :: evt
    type(process_instance_t), intent(in) :: process_instance
    integer :: n_born
    associate (event_deps => evt%event_deps)
       event_deps%lab_is_cm = process_instance%lab_is_cm (1)
       select type (pcm => process_instance%pcm)
       type is (pcm_nlo_t)
          n_born = pcm%region_data%n_legs_born
       end select
       call event_deps%p_born_cms%init (n_born, 1)
       call event_deps%p_born_lab%init (n_born, 1)
    end associate
  end subroutine evt_nlo_setup_general_event_kinematics

  module subroutine evt_nlo_setup_real_event_kinematics (evt, process_instance)
    class(evt_nlo_t), intent(inout) :: evt
    type(process_instance_t), intent(in) :: process_instance
    integer :: n_real, n_phs
    integer :: i_real
    associate (event_deps => evt%event_deps)
       select type (pcm => process_instance%pcm)
       class is (pcm_nlo_t)
          n_real = pcm%region_data%n_legs_real
       end select
       i_real = evt%process%get_first_real_term ()
       select type (phs => process_instance%kin(i_real)%phs)
       type is (phs_fks_t)
          event_deps%phs_identifiers = phs%phs_identifiers
       end select
       n_phs = size (event_deps%phs_identifiers)
       call event_deps%p_real_cms%init (n_real, n_phs)
       call event_deps%p_real_lab%init (n_real, n_phs)
       select type (pcm => process_instance%pcm)
       type is (pcm_nlo_t)
          if (allocated (pcm%region_data%alr_contributors)) then
             allocate (event_deps%contributors &
                  (size (pcm%region_data%alr_contributors)))
             event_deps%contributors = pcm%region_data%alr_contributors
          end if
          if (allocated (pcm%region_data%alr_to_i_contributor)) then
             allocate (event_deps%alr_to_i_con &
                  (size (pcm%region_data%alr_to_i_contributor)))
             event_deps%alr_to_i_con = pcm%region_data%alr_to_i_contributor
          end if
       end select
    end associate
  end subroutine evt_nlo_setup_real_event_kinematics

  module subroutine evt_nlo_set_mode (evt, process_instance)
    class(evt_nlo_t), intent(inout) :: evt
    type(process_instance_t), intent(in) :: process_instance
    integer :: i_real
    select type (pcm => process_instance%pcm)
    type is (pcm_nlo_t)
       if (pcm%settings%combined_integration) then
          evt%mode = EVT_NLO_COMBINED
       else
          i_real = evt%process%get_first_real_component ()
          if (i_real == evt%process%extract_active_component_mci ()) then
             evt%mode = EVT_NLO_SEPARATE_REAL
          else
             evt%mode = EVT_NLO_SEPARATE_BORNLIKE
          end if
       end if
    end select
  end subroutine evt_nlo_set_mode

  module function evt_nlo_is_valid_event (evt, i_term) result (valid)
    logical :: valid
    class(evt_nlo_t), intent(in) :: evt
    integer, intent(in) :: i_term
    valid = evt%process_instance%term(i_term)%passed
  end function evt_nlo_is_valid_event

  module function evt_nlo_get_selected_quantum_numbers &
       (evt, i_flv) result (qn_select)
    class(evt_nlo_t), intent(in) :: evt
    integer, intent(in) :: i_flv
    type(quantum_numbers_t), dimension(:), allocatable :: qn_select
    integer :: i_term, index
    i_term = evt%get_i_term ()
    associate (term => evt%process_instance%term(i_term))
       index = term%connected%matrix%get_qn_index (i_flv, i_sub = 0)
       qn_select = term%connected%matrix%get_quantum_numbers (index)
    end associate
  end function evt_nlo_get_selected_quantum_numbers

  module subroutine evt_nlo_prepare_new_event (evt, i_mci, i_term)
    class(evt_nlo_t), intent(inout) :: evt
    integer, intent(in) :: i_mci, i_term
    real(default) :: s, x
    real(default) :: sqme_total
    real(default), dimension(:), allocatable :: sqme_flv
    integer :: i, i_flv, i_core, emitter, n_in
    logical, save :: warn_once = .true.
    class(prc_core_t), pointer :: core => null ()
    call evt%reset ()
    call evt%rng%generate (x)
    do i = 1, size (evt%process_instance%term)
       associate (term => evt%process_instance%term(i))
          if (evt%i_evaluation == 0) then
             if (term%nlo_type == BORN) then
                allocate (sqme_flv (term%config%data%n_flv))
                exit
             end if
          else
             if (term%nlo_type == NLO_REAL .and. .not. term%is_subtraction()) then
                allocate (sqme_flv (term%config%data%n_flv))
                exit
             end if
          end if
       end associate
    end do
    sqme_total = zero
    sqme_flv = zero
    i_core = evt%process%get_i_core (i_term)
    core => evt%process%get_core_ptr (i_core)
    do i = 1, size (evt%process_instance%term)
       associate (term => evt%process_instance%term(i))
          if (i == evt%i_evaluation + 1 .and. (term%nlo_type == BORN .or. &
               (term%nlo_type == NLO_REAL .and. .not. term%is_subtraction())) ) then
             sqme_total = sqme_total + real (sum ( term%connected%matrix%get_matrix_element ()))
             !!! TODO (VR 2020-02-19) figure out why this select type is needed for prc_omega_t
             !!! For NLO and prc_omega_t the connected trace seems to be set up incorrectly!
             !!! (PS 2020-11-05) This leads to real events of processes with structure functions
             !!! having a wrong flavor distribution if computed with O'Mega.
             !!! The flavor distributions are identical with and also without the special case
             !!! for O'Mega and wrong in both cases.
             !!! However, this case it is not critical as long as O'Mega does not provide matrix elements
             !!! exclusive in coupling orders and is thus only rarely used for NLO applications anyways
             select type (core)
             class is (prc_external_t)
                do i_flv = 1, size (sqme_flv)
                   if (allocated (term%passed_array)) then
                      if (term%passed .and. .not. term%passed_array(i_flv)) cycle
                   end if
                   sqme_flv(i_flv) = sqme_flv(i_flv) &
                        + real (term%connected%matrix%get_matrix_element ( &
                        term%connected%matrix%get_qn_index (i_flv, i_sub = 0)))
                end do
             class default
                sqme_flv = sqme_flv &
                     + real (term%connected%matrix%get_matrix_element ())
                emitter = evt%process_instance%kin(i)%emitter
                n_in = evt%process_instance%kin(i)%n_in
                if (warn_once .and. term%nlo_type == NLO_REAL .and. emitter <= n_in) then
                   warn_once = .false.
                   call msg_warning("evt_nlo_prepare_new_event: fNLO flavor&
                        & distributions with O'Mega are wrong.")
                end if
             end select
          end if
       end associate
    end do
    if (debug2_active (D_TRANSFORMS)) then
       if (.not. nearly_equal(sqme_total, sum (sqme_flv))) then
          call msg_warning ("evt_nlo_prepare_new_event: &
          &sum over flavored sqmes does not match total sqme.")
       end if
    end if
    !!! Need absolute values to take into account negative weights
    x = x * abs (sqme_total)
    s = abs (sqme_flv (1))
    evt%selected_i_flv = 1
    if (s < x) then
       do i_flv = 2, size (sqme_flv)
          s = s + abs (sqme_flv (i_flv))
          if (s > x) then
             evt%selected_i_flv = i_flv
             exit
          end if
       end do
    end if
    if (debug2_active (D_TRANSFORMS)) then
       call msg_print_color ("Selected i_flv: ", COL_GREEN)
       print *, evt%selected_i_flv
    end if
  end subroutine evt_nlo_prepare_new_event


end submodule evt_nlo_s

