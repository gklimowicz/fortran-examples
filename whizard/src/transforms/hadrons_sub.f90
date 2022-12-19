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

submodule (hadrons) hadrons_s

  use debug_master, only: debug_on
  use constants
  use diagnostics
  use format_utils, only: write_separator
  use helicities
  use hep_common
  use io_units
  use numeric_utils, only: vanishes
  use subevents

  implicit none

contains

  elemental module function hadrons_method_of_string (string) result (i)
    integer :: i
    type(string_t), intent(in) :: string
    select case (char(string))
    case ("WHIZARD")
       i = HADRONS_WHIZARD
    case ("PYTHIA6")
       i = HADRONS_PYTHIA6
    case ("PYTHIA8")
       i = HADRONS_PYTHIA8
    case default
       i = HADRONS_UNDEFINED
    end select
  end function hadrons_method_of_string

  elemental module function hadrons_method_to_string (i) result (string)
    type(string_t) :: string
    integer, intent(in) :: i
    select case (i)
    case (HADRONS_WHIZARD)
       string = "WHIZARD"
    case (HADRONS_PYTHIA6)
       string = "PYTHIA6"
    case (HADRONS_PYTHIA8)
       string = "PYTHIA8"
    case default
       string = "UNDEFINED"
    end select
  end function hadrons_method_to_string

  module subroutine hadron_settings_init (hadron_settings, var_list)
    class(hadron_settings_t), intent(out) :: hadron_settings
    type(var_list_t), intent(in) :: var_list
    hadron_settings%active = &
         var_list%get_lval (var_str ("?hadronization_active"))
    hadron_settings%method = hadrons_method_of_string ( &
         var_list%get_sval (var_str ("$hadronization_method")))
    hadron_settings%enhanced_fraction = &
         var_list%get_rval (var_str ("hadron_enhanced_fraction"))
    hadron_settings%enhanced_width = &
         var_list%get_rval (var_str ("hadron_enhanced_width"))
  end subroutine hadron_settings_init

  module subroutine hadron_settings_write (settings, unit)
    class(hadron_settings_t), intent(in) :: settings
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)")  "Hadronization settings:"
    call write_separator (u)
    write (u, "(1x,A)")  "Master switches:"
    write (u, "(3x,A,1x,L1)") &
         "active                = ", settings%active
    write (u, "(1x,A)")  "General settings:"
    if (settings%active) then
       write (u, "(3x,A)") &
            "hadron_method         =  " // &
            char (hadrons_method_to_string (settings%method))
    else
       write (u, "(3x,A)") " [Hadronization off]"
    end if
    write (u, "(1x,A)")  "pT generation parameters"
    write (u, "(3x,A,1x,ES19.12)") &
         "enhanced_fraction     = ", settings%enhanced_fraction
    write (u, "(3x,A,1x,ES19.12)") &
         "enhanced_width        = ", settings%enhanced_width
  end subroutine hadron_settings_write

  pure module subroutine hadrons_import_rng (hadrons, rng)
    class(hadrons_t), intent(inout) :: hadrons
    class(rng_t), intent(inout), allocatable :: rng
    call move_alloc (from = rng, to = hadrons%rng)
  end subroutine hadrons_import_rng

  module subroutine hadrons_hadrons_init &
       (hadrons, shower_settings, hadron_settings, model_hadrons)
    class(hadrons_hadrons_t), intent(out) :: hadrons
    type(shower_settings_t), intent(in) :: shower_settings
    type(hadron_settings_t), intent(in) :: hadron_settings
    type(model_t), intent(in), target :: model_hadrons
    hadrons%model => model_hadrons
    hadrons%shower_settings = shower_settings
    hadrons%hadron_settings = hadron_settings
    call msg_message &
         ("Hadronization: WHIZARD model for hadronization and decays")
  end subroutine hadrons_hadrons_init

  module subroutine hadrons_hadrons_hadronize (hadrons, particle_set, valid)
    class(hadrons_hadrons_t), intent(inout) :: hadrons
    type(particle_set_t), intent(in) :: particle_set
    logical, intent(out) :: valid
    integer, dimension(:), allocatable :: cols, acols, octs
    integer :: n
    if (signal_is_pending ()) return
    if (debug_on) call msg_debug (D_TRANSFORMS, "hadrons_hadrons_hadronize")
    call particle_set%write (6, compressed=.true.)
    n = particle_set%get_n_tot ()
    allocate (cols (n), acols (n), octs (n))
    call extract_color_systems (particle_set, cols, acols, octs)
    print *, "size(cols)  = ", size (cols)
    if (size(cols) > 0) then
       print *, "cols  = ", cols
    end if
    print *, "size(acols) = ", size(acols)
    if (size(acols) > 0) then
       print *, "acols = ", acols
    end if
    print *, "size(octs)  = ", size(octs)
    if (size (octs) > 0) then
       print *, "octs  = ", octs
    end if
    !!! if all arrays are empty, i.e. zero particles found, nothing to do
  end subroutine hadrons_hadrons_hadronize

  module subroutine lund_pt_init (lund_pt, settings)
    class (lund_pt_t), intent(out) :: lund_pt
    type(hadron_settings_t), intent(in) :: settings
  end subroutine lund_pt_init

  module subroutine hadrons_hadrons_make_particle_set &
         (hadrons, particle_set, model, valid)
    class(hadrons_hadrons_t), intent(in) :: hadrons
    type(particle_set_t), intent(inout) :: particle_set
    class(model_data_t), intent(in), target :: model
    logical, intent(out) :: valid
    if (signal_is_pending ()) return
    valid = .false.
    if (valid) then
    else
       call msg_fatal ("WHIZARD hadronization not yet implemented")
    end if
  end subroutine hadrons_hadrons_make_particle_set

  subroutine extract_color_systems (p_set, cols, acols, octs)
    type(particle_set_t), intent(in) :: p_set
    integer, dimension(:), allocatable, intent(out) :: cols, acols, octs
    logical, dimension(:), allocatable :: mask
    integer :: i, n, n_cols, n_acols, n_octs
    n = p_set%get_n_tot ()
    allocate (mask (n))
    do i = 1, n
       mask(i) = p_set%prt(i)%col%get_col () /= 0 .and. &
            p_set%prt(i)%col%get_acl () == 0 .and. &
            p_set%prt(i)%get_status () == PRT_OUTGOING
    end do
    n_cols = count (mask)
    allocate (cols (n_cols))
    cols = p_set%get_indices (mask)
    do i = 1, n
       mask(i) = p_set%prt(i)%col%get_col () == 0 .and. &
            p_set%prt(i)%col%get_acl () /= 0 .and. &
            p_set%prt(i)%get_status () == PRT_OUTGOING
    end do
    n_acols = count (mask)
    allocate (acols (n_acols))
    acols = p_set%get_indices (mask)
    do i = 1, n
       mask(i) = p_set%prt(i)%col%get_col () /= 0 .and. &
            p_set%prt(i)%col%get_acl () /= 0 .and. &
            p_set%prt(i)%get_status () == PRT_OUTGOING
    end do
    n_octs = count (mask)
    allocate (octs (n_octs))
    octs = p_set%get_indices (mask)
  end subroutine extract_color_systems

  module subroutine hadrons_pythia6_init &
       (hadrons, shower_settings, hadron_settings, model_hadrons)
    class(hadrons_pythia6_t), intent(out) :: hadrons
    type(shower_settings_t), intent(in) :: shower_settings
    type(hadron_settings_t), intent(in) :: hadron_settings
    type(model_t), intent(in), target :: model_hadrons
    logical :: pygive_not_set_by_shower
    hadrons%model => model_hadrons
    hadrons%shower_settings = shower_settings
    hadrons%hadron_settings = hadron_settings
    pygive_not_set_by_shower = .not. (shower_settings%method == PS_PYTHIA6 &
         .and. (shower_settings%isr_active .or. shower_settings%fsr_active))
    if (pygive_not_set_by_shower) then
       call pythia6_set_verbose (shower_settings%verbose)
       call pythia6_set_config (shower_settings%pythia6_pygive)
    end if
    call msg_message &
         ("Hadronization: Using PYTHIA6 interface for hadronization and decays")
  end subroutine hadrons_pythia6_init

  module subroutine hadrons_pythia6_hadronize (hadrons, particle_set, valid)
    class(hadrons_pythia6_t), intent(inout) :: hadrons
    type(particle_set_t), intent(in) :: particle_set
    logical, intent(out) :: valid
    integer :: N, NPAD, K
    real(double) :: P, V
    common /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)
    save /PYJETS/
    if (signal_is_pending ()) return
    if (debug_on) call msg_debug (D_TRANSFORMS, "hadrons_pythia6_hadronize")
    call pygive ("MSTP(111)=1")    !!! Switch on hadronization and decays
    call pygive ("MSTJ(1)=1")      !!! String fragmentation
    call pygive ("MSTJ(21)=2")     !!! String fragmentation keeping resonance momentum
    call pygive ("MSTJ(28)=0")     !!! Switch off tau decays
    if (debug_active (D_TRANSFORMS)) then
       call msg_debug (D_TRANSFORMS, "N", N)
       call pylist(2)
       print *, ' line 7 : ', k(7,1:5), p(7,1:5)
    end if
    call pyedit (12)
    call pythia6_set_last_treated_line (N)
    call pyexec ()
    call pyedit (12)
    valid = .true.
  end subroutine hadrons_pythia6_hadronize

  module subroutine hadrons_pythia6_make_particle_set &
       (hadrons, particle_set, model, valid)
    class(hadrons_pythia6_t), intent(in) :: hadrons
    type(particle_set_t), intent(inout) :: particle_set
    class(model_data_t), intent(in), target :: model
    logical, intent(out) :: valid
    if (signal_is_pending ()) return
    valid = pythia6_handle_errors ()
    if (valid) then
       call pythia6_combine_with_particle_set &
            (particle_set, model, hadrons%model, hadrons%shower_settings)
    end if
  end subroutine hadrons_pythia6_make_particle_set

  module subroutine hadrons_pythia8_init &
       (hadrons, shower_settings, hadron_settings, model_hadrons)
    class(hadrons_pythia8_t), intent(out) :: hadrons
    type(shower_settings_t), intent(in) :: shower_settings
    type(hadron_settings_t), intent(in) :: hadron_settings
    type(model_t), intent(in), target :: model_hadrons
    hadrons%model => model_hadrons
    hadrons%shower_settings = shower_settings
    hadrons%hadron_settings = hadron_settings
    call msg_message ("Hadronization: Using PYTHIA8 interface " // &
         "for hadronization and decays.")
    ! TODO sbrass which verbose?
    call hadrons%pythia%init (verbose = shower_settings%verbose)
    call hadrons%lhaup%init ()
  end subroutine hadrons_pythia8_init

  module subroutine hadrons_pythia8_transfer_settings (hadrons)
    class(hadrons_pythia8_t), intent(inout), target :: hadrons
    real(default) :: r
    if (debug_on)  call msg_debug &
         (D_TRANSFORMS, "hadrons_pythia8_transfer_settings")
    if (debug_on)  call msg_debug2 &
         (D_TRANSFORMS, "pythia_initialized", hadrons%pythia_initialized)
    if (hadrons%pythia_initialized) return
    call hadrons%pythia%import_rng (hadrons%rng)
    call hadrons%pythia%parse_and_set_config &
         (hadrons%shower_settings%pythia8_config)
    if (len (hadrons%shower_settings%pythia8_config_file) > 0) &
         call hadrons%pythia%read_file &
         (hadrons%shower_settings%pythia8_config_file)
    call hadrons%pythia%read_string (var_str ("Beams:frameType = 5"))
    call hadrons%pythia%read_string (var_str ("ProcessLevel:all = off"))
    if (.not. hadrons%shower_settings%verbose) then
       call hadrons%pythia%read_string (var_str ("Print:quiet = on"))
    end if
    call hadrons%pythia%set_lhaup_ptr (hadrons%lhaup)
    call hadrons%pythia%init_pythia ()
    hadrons%pythia_initialized = .true.
  end subroutine hadrons_pythia8_transfer_settings

  module subroutine hadrons_pythia8_set_user_process (hadrons, pset)
    class(hadrons_pythia8_t), intent(inout) :: hadrons
    type(particle_set_t), intent(in) :: pset
    integer, dimension(2) :: beam_pdg
    real(default), dimension(2) :: beam_energy
    integer, parameter :: process_id = 0, n_processes = 0
    if (debug_on)  call msg_debug &
         (D_TRANSFORMS, "hadrons_pythia8_set_user_process")
    beam_pdg = [pset%prt(1)%get_pdg (), pset%prt(2)%get_pdg ()]
    beam_energy = [energy(pset%prt(1)%p), energy(pset%prt(2)%p)]
    call hadrons%lhaup%set_init (beam_pdg, beam_energy, &
         n_processes, unweighted = .false., negative_weights = .false.)
    call hadrons%lhaup%set_process_parameters (process_id = process_id, &
         cross_section = one, error = one)
  end subroutine hadrons_pythia8_set_user_process

  module subroutine hadrons_pythia8_import_particle_set (hadrons, particle_set)
    class(hadrons_pythia8_t), target, intent(inout) :: hadrons
    type(particle_set_t), intent(in) :: particle_set
    integer, parameter :: PROCESS_ID = 1
    if (debug_on)  call msg_debug &
         (D_TRANSFORMS, "hadrons_pythia8_import_particle_set")
    if (.not. hadrons%user_process_set) then
       call hadrons%set_user_process (particle_set)
       hadrons%user_process_set = .true.
    end if
    call hadrons%lhaup%set_event_process (process_id = PROCESS_ID, &
         scale = -one, alpha_qcd = -one, alpha_qed = -one, weight = -one)
    call hadrons%lhaup%set_event (process_id = PROCESS_ID, &
         particle_set = particle_set, polarization = .true.)
    if (debug_active (D_TRANSFORMS)) then
       call hadrons%lhaup%list_init ()
    end if
  end subroutine hadrons_pythia8_import_particle_set

  module subroutine hadrons_pythia8_hadronize (hadrons, particle_set, valid)
    class(hadrons_pythia8_t), intent(inout) :: hadrons
    type(particle_set_t), intent(in) :: particle_set
    logical, intent(out) :: valid
    if (signal_is_pending ()) return
    call hadrons%import_particle_set (particle_set)
    if (.not. hadrons%pythia_initialized) &
         call hadrons%transfer_settings ()
    call hadrons%pythia%next (valid)
    if (debug_active (D_TRANSFORMS)) then
       call hadrons%pythia%list_event ()
       call particle_set%write (summary=.true., compressed=.true.)
    end if
  end subroutine hadrons_pythia8_hadronize

  module subroutine hadrons_pythia8_make_particle_set &
         (hadrons, particle_set, model, valid)
    class(hadrons_pythia8_t), intent(in) :: hadrons
    type(particle_set_t), intent(inout) :: particle_set
    class(model_data_t), intent(in), target :: model
    logical, intent(out) :: valid
    type(particle_t), dimension(:), allocatable :: beam
    if (debug_on)  call msg_debug &
         (D_TRANSFORMS, "hadrons_pythia8_make_particle_set")
    if (signal_is_pending ()) return
    associate (settings => hadrons%shower_settings)
      if (debug_active (D_TRANSFORMS)) then
         call msg_debug (D_TRANSFORMS, 'Combine PYTHIA8 with particle set')
         call msg_debug (D_TRANSFORMS, 'Particle set before replacing')
         call particle_set%write (summary=.true., compressed=.true.)
         call hadrons%pythia%list_event ()
         call msg_debug (D_TRANSFORMS, string = "settings%hadron_collision", &
              value = settings%hadron_collision)
      end if
      call hadrons%pythia%get_hadron_particles (&
           model, hadrons%model, particle_set, &
           helicity = PRT_DEFINITE_HELICITY)
    end associate
    if (debug_active (D_TRANSFORMS)) then
       print *, 'Particle set after replacing'
       call particle_set%write (summary=.true., compressed=.true.)
    end if
    valid = .true.
  end subroutine hadrons_pythia8_make_particle_set

  module subroutine evt_hadrons_init (evt, model_hadrons)
    class(evt_hadrons_t), intent(out) :: evt
    type(model_t), intent(in), target :: model_hadrons
    evt%model_hadrons => model_hadrons
    evt%is_first_event = .true.
  end subroutine evt_hadrons_init

  module subroutine evt_hadrons_write_name (evt, unit)
    class(evt_hadrons_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Event transform: hadronization"
  end subroutine evt_hadrons_write_name

  module subroutine evt_hadrons_write &
       (evt, unit, verbose, more_verbose, testflag)
    class(evt_hadrons_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, more_verbose, testflag
    integer :: u
    u = given_output_unit (unit)
    call write_separator (u, 2)
    call evt%write_name (u)
    call write_separator (u)
    call evt%base_write (u, testflag = testflag, show_set = .false.)
    if (evt%particle_set_exists)  &
         call evt%particle_set%write &
         (u, summary = .true., compressed = .true., testflag = testflag)
    call write_separator (u)
    call evt%hadrons%shower_settings%write (u)
    call write_separator (u)
    call evt%hadrons%hadron_settings%write (u)
  end subroutine evt_hadrons_write

  module subroutine evt_hadrons_first_event (evt)
    class(evt_hadrons_t), intent(inout) :: evt
    if (debug_on) call msg_debug (D_TRANSFORMS, "evt_hadrons_first_event")
    associate (settings => evt%hadrons%shower_settings)
       settings%hadron_collision = .false.
       if (all (evt%particle_set%prt(1:2)%flv%get_pdg_abs () <= 39)) then
          settings%hadron_collision = .false.
       else if (all (evt%particle_set%prt(1:2)%flv%get_pdg_abs () >= 100)) then
          settings%hadron_collision = .true.
       else
          call msg_fatal ("evt_hadrons didn't recognize beams setup")
       end if
       if (debug_on)  call msg_debug &
            (D_TRANSFORMS, "hadron_collision", settings%hadron_collision)
       if (.not. (settings%isr_active .or. settings%fsr_active)) then
          call msg_fatal ("Hadronization without shower is not supported")
       end if
    end associate
    evt%is_first_event = .false.
  end subroutine evt_hadrons_first_event

  module subroutine evt_hadrons_generate_weighted (evt, probability)
    class(evt_hadrons_t), intent(inout) :: evt
    real(default), intent(inout) :: probability
    logical :: valid
    if (signal_is_pending ())  return
    evt%particle_set = evt%previous%particle_set
    if (evt%is_first_event) then
       call evt%first_event ()
    end if
    call evt%hadrons%hadronize (evt%particle_set, valid)
    probability = 1
    evt%particle_set_exists = valid
  end subroutine evt_hadrons_generate_weighted

  module subroutine evt_hadrons_make_particle_set &
       (evt, factorization_mode, keep_correlations, r)
    class(evt_hadrons_t), intent(inout) :: evt
    integer, intent(in) :: factorization_mode
    logical, intent(in) :: keep_correlations
    real(default), dimension(:), intent(in), optional :: r
    logical :: valid
    call evt%hadrons%make_particle_set (evt%particle_set, evt%model, valid)
    evt%particle_set_exists = evt%particle_set_exists .and. valid
  end subroutine evt_hadrons_make_particle_set

  module subroutine evt_hadrons_connect &
       (evt, process_instance, model, process_stack)
    class(evt_hadrons_t), intent(inout), target :: evt
    type(process_instance_t), intent(in), target :: process_instance
    class(model_data_t), intent(in), target :: model
    type(process_stack_t), intent(in), optional :: process_stack
    call evt%base_connect (process_instance, model, process_stack)
    call evt%make_rng (evt%process)
  end subroutine evt_hadrons_connect

  module subroutine evt_hadrons_make_rng (evt, process)
    class(evt_hadrons_t), intent(inout) :: evt
    type(process_t), intent(inout) :: process
    class(rng_t), allocatable :: rng
    call process%make_rng (rng)
    call evt%hadrons%import_rng (rng)
  end subroutine evt_hadrons_make_rng

  module subroutine evt_hadrons_prepare_new_event (evt, i_mci, i_term)
    class(evt_hadrons_t), intent(inout) :: evt
    integer, intent(in) :: i_mci, i_term
    call evt%reset ()
  end subroutine evt_hadrons_prepare_new_event


end submodule hadrons_s

