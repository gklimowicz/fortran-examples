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

submodule (shower_pythia8) shower_pythia8_s

  use debug_master, only: debug_on
  use constants
  use numeric_utils, only: vanishes
  use io_units
  use physics_defs
  use diagnostics
  use helicities
  use subevents

  implicit none

contains

  module subroutine shower_pythia8_init &
       (shower, settings, taudec_settings, pdf_data, os_data)
    class(shower_pythia8_t), intent(out) :: shower
    type(shower_settings_t), intent(in) :: settings
    type(taudec_settings_t), intent(in) :: taudec_settings
    type(pdf_data_t), intent(in) :: pdf_data
    type(os_data_t), intent(in) :: os_data
    if (debug_on) call msg_debug (D_SHOWER, "shower_pythia8_init")
    shower%settings = settings
    shower%taudec_settings = taudec_settings
    shower%os_data = os_data
    call shower%pdf_data%init (pdf_data)
    shower%name = "PYTHIA8"
    call shower%write_msg ()
    call shower%pythia%init (verbose = settings%verbose)
    call shower%lhaup%init ()
  end subroutine shower_pythia8_init

  module subroutine shower_pythia8_set_user_process (shower, pset)
    class(shower_pythia8_t), intent(inout) :: shower
    type(particle_set_t), intent(in) :: pset
    integer, dimension(2) :: beam_pdg
    real(default), dimension(2) :: beam_energy
    integer, parameter :: process_id = 1, n_processes = 1
    if (debug_on) call msg_debug (D_SHOWER, "shower_pythia8_set_user_process")
    ! TODO sbrass find correct beam entries, fallback would be first two entries
    beam_pdg = [pset%prt(1)%get_pdg (), pset%prt(2)%get_pdg ()]
    beam_energy = [energy(pset%prt(1)%p), energy(pset%prt(2)%p)]
    call shower%lhaup%set_init (beam_pdg, beam_energy, &
         n_processes, unweighted = .false., negative_weights = .false.)
    call shower%lhaup%set_process_parameters (process_id = process_id, &
         cross_section = one, error = one)
  end subroutine shower_pythia8_set_user_process
  module subroutine shower_pythia8_import_particle_set &
       (shower, particle_set)
    class(shower_pythia8_t), target, intent(inout) :: shower
    type(particle_set_t), intent(in) :: particle_set
    type(particle_set_t) :: pset_reduced
    integer, parameter :: PROCESS_ID = 1
    logical :: keep_beams
    if (debug_on) call msg_debug (D_SHOWER, "shower_pythia8_import_particle_set")
    if (.not. shower%user_process_set) then
       call shower%set_user_process (particle_set)
       shower%user_process_set = .true.
    end if
    if (debug_active (D_SHOWER)) then
       call particle_set%write (summary=.true., compressed=.true.)
    end if
    call shower%lhaup%set_event_process (process_id = PROCESS_ID, scale = shower%fac_scale, &
         alpha_qcd = shower%alpha_s, alpha_qed = -one, weight = -one)
    call shower%lhaup%set_event (process_id = PROCESS_ID, particle_set = particle_set, &
         keep_beams = .false., keep_remnants = .true., polarization = .true.)
    if (debug_active (D_SHOWER)) then
       call shower%lhaup%list_init ()
       call shower%lhaup%list_event ()
    end if
  end subroutine shower_pythia8_import_particle_set

  module subroutine shower_pythia8_generate_emissions &
        (shower, valid, number_of_emissions)
    class(shower_pythia8_t), intent(inout), target :: shower
    logical, intent(out) :: valid
    integer, optional, intent(in) :: number_of_emissions
    if (signal_is_pending ()) return
    call shower%transfer_settings ()
    call shower%pythia%next (valid)
  end subroutine shower_pythia8_generate_emissions

  module subroutine shower_pythia8_make_particle_set &
       (shower, particle_set, model, model_hadrons)
    class(shower_pythia8_t), intent(in) :: shower
    type(particle_set_t), intent(inout) :: particle_set
    class(model_data_t), intent(in), target :: model
    class(model_data_t), intent(in), target :: model_hadrons
    type(particle_t), dimension(:), allocatable :: beam
    integer :: n_whizard, n_tot_pythia
    if (debug_on) call msg_debug (D_SHOWER, "shower_pythia8_make_particle_set")
    if (signal_is_pending ()) return
    associate (settings => shower%settings)
      if (debug_active (D_SHOWER)) then
         call msg_debug (D_SHOWER, 'Combine PYTHIA8 with particle set')
         call msg_debug (D_SHOWER, 'Particle set before replacing')
         call particle_set%write (summary=.true., compressed=.true.)
         call shower%pythia%list_event ()
         call msg_debug (D_SHOWER, string = "settings%hadron_collision", &
              value = settings%hadron_collision)
      end if
    end associate
    call shower%pythia%get_shower_particles (&
         model, model_hadrons, particle_set, &
         helicity = PRT_DEFINITE_HELICITY, &
         recover_beams = shower%settings%hadron_collision)
    if (debug_active (D_SHOWER)) then
       print *, 'Particle set after replacing'
       call particle_set%write (summary=.true., compressed=.true.)
    end if
  end subroutine shower_pythia8_make_particle_set

  module subroutine shower_pythia8_transfer_settings (shower)
    class(shower_pythia8_t), intent(inout), target :: shower
    if (debug_on)  call msg_debug &
         (D_SHOWER, "shower_pythia8_transfer_settings")
    if (debug_on)  call msg_debug2 &
         (D_SHOWER, "pythia_initialized", shower%pythia_initialized)
    if (shower%pythia_initialized) return
    associate (pythia => shower%pythia)
      call pythia%set_lhaup_ptr (shower%lhaup)
      call pythia%import_rng (shower%rng)
      call shower%pythia%parse_and_set_config (shower%settings%pythia8_config)
      if (len (shower%settings%pythia8_config_file) > 0) &
           call pythia%read_file (shower%settings%pythia8_config_file)
      call pythia%read_string (var_str ("Beams:frameType = 5"))
      ! call pythia%read_string (var_str ("ParticleDecays:allowPhotonRadiation = off"))
      call pythia%read_string (var_str ("HadronLevel:all = off"))
      if (.not. shower%settings%verbose) then
         call pythia%read_string (var_str ("Print:quiet = on"))
      end if
      if (.not. shower%settings%isr_active) then
         call pythia%read_string (var_str ("PartonLevel:ISR = off"))
      else
         call pythia%read_string (var_str ("PartonLevel:ISR = on"))
      end if
      if (.not. shower%settings%fsr_active) then
         call pythia%read_string (var_str ("PartonLevel:FSR = off"))
      else
         call pythia%read_string (var_str ("PartonLevel:FSR = on"))
      end if
      call pythia%init_pythia ()
    end associate
    shower%pythia_initialized = .true.
  end subroutine shower_pythia8_transfer_settings

  module subroutine shower_pythia8_get_final_colored_ME_momenta &
         (shower, momenta)
    class(shower_pythia8_t), intent(in) :: shower
    type(vector4_t), dimension(:), allocatable, intent(out) :: momenta
    if (debug_on)  call msg_debug &
         (D_MATCHING, "shower_pythia8_get_final_colored_ME_momenta")
    call shower%pythia%get_final_colored_ME_momenta (momenta)
  end subroutine shower_pythia8_get_final_colored_ME_momenta


end submodule shower_pythia8_s

