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

module dispatch_transforms

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use process
  use variables
  use system_defs, only: LF
  use system_dependencies, only: LHAPDF6_AVAILABLE
  use sf_lhapdf, only: lhapdf_initialize
  use pdf, only: pdf_data_t
  use diagnostics
  use models
  use os_interface
  use beam_structures
  use resonances, only: resonance_history_set_t
  use instances, only: process_instance_t, process_instance_hook_t

  use event_base, only: event_callback_t, event_callback_nop_t
  use hepmc_interface, only: HEPMC3_MODE_HEPMC2, HEPMC3_MODE_HEPMC3
  use hepmc_interface, only: HEPMC3_MODE_ROOT, HEPMC3_MODE_ROOTTREE
  use hepmc_interface, only: HEPMC3_MODE_HEPEVT
  use eio_base
  use eio_raw
  use eio_checkpoints
  use eio_callback
  use eio_lhef
  use eio_hepmc
  use eio_lcio
  use eio_stdhep
  use eio_ascii
  use eio_weights
  use eio_dump

  use event_transforms
  use resonance_insertion
  use isr_epa_handler
  use decays
  use shower_base
  use shower_core
  use shower
  use shower_pythia6
  use shower_pythia8
  use hadrons
  use mlm_matching
  use powheg_matching
  use ckkw_matching
  use tauola_interface !NODEP!
  use evt_nlo

  implicit none
  private

  public :: dispatch_evt_nlo
  public :: dispatch_evt_resonance
  public :: dispatch_evt_isr_epa_handler
  public :: dispatch_evt_decay
  public :: dispatch_evt_shower
  public :: dispatch_evt_shower_hook
  public :: dispatch_matching
  public :: dispatch_evt_hadrons
  public :: dispatch_eio

contains

  subroutine dispatch_evt_nlo (evt, keep_failed_events)
    class(evt_t), intent(out), pointer :: evt
    logical, intent(in) :: keep_failed_events
    call msg_message ("Simulate: activating fixed-order NLO events")
    allocate (evt_nlo_t :: evt)
    evt%only_weighted_events = .true.
    select type (evt)
    type is (evt_nlo_t)
       evt%i_evaluation = 0
       evt%keep_failed_events = keep_failed_events
    end select
  end subroutine dispatch_evt_nlo

  subroutine dispatch_evt_resonance (evt, var_list, res_history_set, libname)
    class(evt_t), intent(out), pointer :: evt
    type(var_list_t), intent(in) :: var_list
    type(resonance_history_set_t), dimension(:), intent(in) :: res_history_set
    type(string_t), intent(in) :: libname
    logical :: resonance_history
    resonance_history = var_list%get_lval (var_str ("?resonance_history"))
    if (resonance_history) then
       allocate (evt_resonance_t :: evt)
       call msg_message ("Simulate: activating resonance insertion")
       select type (evt)
       type is (evt_resonance_t)
          call evt%set_resonance_data (res_history_set)
          call evt%set_library (libname)
       end select
    else
       evt => null ()
    end if
  end subroutine dispatch_evt_resonance

  subroutine dispatch_evt_isr_epa_handler (evt, var_list)
    class(evt_t), intent(out), pointer :: evt
    type(var_list_t), intent(in) :: var_list
    logical :: isr_recoil
    logical :: epa_recoil
    logical :: isr_handler_active
    logical :: epa_handler_active
    type(string_t) :: isr_handler_mode
    type(string_t) :: epa_handler_mode
    logical :: isr_keep_mass
    real(default) :: sqrts
    real(default) :: isr_q_max
    real(default) :: epa_q_max
    real(default) :: isr_mass
    real(default) :: epa_mass
    isr_handler_active = var_list%get_lval (var_str ("?isr_handler"))
    if (isr_handler_active) then
       call msg_message ("Simulate: activating ISR handler")
       isr_recoil = &
            var_list%get_lval (var_str ("?isr_recoil"))
       isr_handler_mode = &
            var_list%get_sval (var_str ("$isr_handler_mode"))
       isr_keep_mass = &
            var_list%get_lval (var_str ("?isr_handler_keep_mass"))
       if (isr_recoil) then
          call msg_fatal ("Simulate: ISR handler is incompatible &
               &with ?isr_recoil=true")
       end if
    end if
    epa_handler_active = var_list%get_lval (var_str ("?epa_handler"))
    if (epa_handler_active) then
       call msg_message ("Simulate: activating EPA handler")
       epa_recoil = var_list%get_lval (var_str ("?epa_recoil"))
       epa_handler_mode = var_list%get_sval (var_str ("$epa_handler_mode"))
       if (epa_recoil) then
          call msg_fatal ("Simulate: EPA handler is incompatible &
               &with ?epa_recoil=true")
       end if
    end if
    if (isr_handler_active .and. epa_handler_active) then
       if (isr_handler_mode /= epa_handler_mode) then
          call msg_fatal ("Simulate: ISR/EPA handler: modes must coincide")
       end if
    end if
    if (isr_handler_active .or. epa_handler_active) then
       allocate (evt_isr_epa_t :: evt)
       select type (evt)
       type is (evt_isr_epa_t)
          if (isr_handler_active) then
             call evt%set_mode_string (isr_handler_mode)
          else
             call evt%set_mode_string (epa_handler_mode)
          end if
          sqrts = var_list%get_rval (var_str ("sqrts"))
          if (isr_handler_active) then
             isr_q_max = var_list%get_rval (var_str ("isr_q_max"))
             isr_mass = var_list%get_rval (var_str ("isr_mass"))
             call evt%set_data_isr (sqrts, isr_q_max, isr_mass, isr_keep_mass)
          end if
          if (epa_handler_active) then
             epa_q_max = var_list%get_rval (var_str ("epa_q_max"))
             epa_mass = var_list%get_rval (var_str ("epa_mass"))
             call evt%set_data_epa (sqrts, epa_q_max, epa_mass)
          end if
          call msg_message ("Simulate: ISR/EPA handler mode: " &
               // char (evt%get_mode_string ()))
       end select
    else
       evt => null ()
    end if
  end subroutine dispatch_evt_isr_epa_handler

  subroutine dispatch_evt_decay (evt, var_list)
    class(evt_t), intent(out), pointer :: evt
    type(var_list_t), intent(in), target :: var_list
    logical :: allow_decays
    allow_decays = var_list%get_lval (var_str ("?allow_decays"))
    if (allow_decays) then
       allocate (evt_decay_t :: evt)
       call msg_message ("Simulate: activating decays")
       select type (evt)
       type is (evt_decay_t)
          call evt%set_var_list (var_list)
       end select
    else
       evt => null ()
    end if
  end subroutine dispatch_evt_decay

  subroutine dispatch_evt_shower (evt, var_list, model, fallback_model, &
       os_data, beam_structure, process)
    class(evt_t), intent(out), pointer :: evt
    type(var_list_t), intent(in) :: var_list
    type(model_t), pointer, intent(in) :: model, fallback_model
    type(os_data_t), intent(in) :: os_data
    type(beam_structure_t), intent(in) :: beam_structure
    type(process_t), intent(in), optional :: process
    type(string_t) :: lhapdf_file, lhapdf_dir, process_name
    integer :: lhapdf_member
    type(shower_settings_t) :: settings
    type(taudec_settings_t) :: taudec_settings
    call msg_message ("Simulate: activating parton shower")
    allocate (evt_shower_t :: evt)
    call settings%init (var_list)
    if (associated (model)) then
       call taudec_settings%init (var_list, model)
    else
       call taudec_settings%init (var_list, fallback_model)
    end if
    if (present (process)) then
       process_name = process%get_id ()
    else
       process_name = 'dispatch_testing'
    end if
    select type (evt)
    type is (evt_shower_t)
       call evt%init (fallback_model, os_data)
       lhapdf_member = &
            var_list%get_ival (var_str ("lhapdf_member"))
       if (LHAPDF6_AVAILABLE) then
          lhapdf_dir = &
               var_list%get_sval (var_str ("$lhapdf_dir"))
          lhapdf_file = &
               var_list%get_sval (var_str ("$lhapdf_file"))
          call lhapdf_initialize &
               (1, lhapdf_dir, lhapdf_file, lhapdf_member, evt%pdf_data%pdf)
       end if
       if (present (process))  call evt%pdf_data%setup ("Shower", &
            beam_structure, lhapdf_member, process%get_pdf_set ())
       select case (settings%method)
       case (PS_WHIZARD)
          allocate (shower_t :: evt%shower)
       case (PS_PYTHIA6)
          allocate (shower_pythia6_t :: evt%shower)
       case (PS_PYTHIA8)
          allocate (shower_pythia8_t :: evt%shower)
       case default
          call msg_fatal ('Shower: Method ' // &
            char (var_list%get_sval (var_str ("$shower_method"))) // &
            'not implemented!')
       end select
       call evt%shower%init (settings, taudec_settings, evt%pdf_data, os_data)
       call dispatch_matching (evt, settings, var_list, process_name, evt%pdf_data)
    class default
       call dispatch_matching (evt, settings, var_list, process_name)
    end select
  end subroutine dispatch_evt_shower

  subroutine dispatch_evt_shower_hook (hook, var_list, process_instance, beam_structure, pdf_set)
    class(process_instance_hook_t), pointer, intent(out) :: hook
    type(var_list_t), intent(in) :: var_list
    class(process_instance_t), intent(in), target :: process_instance
    type(beam_structure_t), intent(in) :: beam_structure
    integer, intent(in) :: pdf_set
    type(pdf_data_t) :: pdf_data
    type(string_t) :: lhapdf_file, lhapdf_dir
    integer :: lhapdf_member
    if (var_list%get_lval (var_str ('?powheg_matching'))) then
       call msg_message ("Integration hook: add POWHEG hook")
       allocate (powheg_matching_hook_t :: hook)
       select type (hook)
       type is (powheg_matching_hook_t)
          lhapdf_member = var_list%get_ival (var_str ("lhapdf_member"))
          if (LHAPDF6_AVAILABLE) then
             lhapdf_dir = var_list%get_sval (var_str ("$lhapdf_dir"))
             lhapdf_file = var_list%get_sval (var_str ("$lhapdf_file"))
             call lhapdf_initialize (1, lhapdf_dir, lhapdf_file, lhapdf_member, pdf_data%pdf)
          end if
          call pdf_data%setup ("Shower", beam_structure, lhapdf_member, pdf_set)
          call hook%init (var_list, process_instance, pdf_data)
       end select
    else
       hook => null ()
    end if
  end subroutine dispatch_evt_shower_hook

  subroutine dispatch_matching (evt, settings, var_list, process_name, pdf_data)
    class(evt_t), intent(inout) :: evt
    type(shower_settings_t), intent(in) :: settings
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_name
    type(pdf_data_t), intent(in), optional :: pdf_data
    select type (evt)
    type is (evt_shower_t)
       if (settings%mlm_matching .and. settings%ckkw_matching) then
          call msg_fatal ("Both MLM and CKKW matching activated," // &
               LF // "     aborting simulation")
       end if
       if (settings%powheg_matching) then
          call msg_message ("Simulate: applying POWHEG matching")
          allocate (powheg_matching_t :: evt%matching)
       end if
       if (settings%mlm_matching) then
          call msg_message ("Simulate: applying MLM matching")
          allocate (mlm_matching_t :: evt%matching)
       end if
       if (settings%ckkw_matching) then
          call msg_warning ("Simulate: CKKW(-L) matching not yet supported")
          allocate (ckkw_matching_t :: evt%matching)
       end if
       if (allocated (evt%matching)) then
          call evt%matching%init (var_list, process_name)
          if (present(pdf_data)) then
             select type (matching => evt%matching)
             type is (powheg_matching_t)
                matching%process_deps%pdf_data = pdf_data
             end select
          end if
       end if
    end select
  end subroutine dispatch_matching

  subroutine dispatch_evt_hadrons (evt, var_list, fallback_model)
    class(evt_t), intent(out), pointer :: evt
    type(var_list_t), intent(in) :: var_list
    type(model_t), pointer, intent(in) :: fallback_model
    type(shower_settings_t) :: shower_settings
    type(hadron_settings_t) :: hadron_settings
    allocate (evt_hadrons_t :: evt)
    call msg_message ("Simulate: activating hadronization")
    call shower_settings%init (var_list)
    call hadron_settings%init (var_list)
    select type (evt)
    type is (evt_hadrons_t)
       call evt%init (fallback_model)
       select case (hadron_settings%method)
       case (HADRONS_WHIZARD)
          allocate (hadrons_hadrons_t :: evt%hadrons)
       case (HADRONS_PYTHIA6)
          allocate (hadrons_pythia6_t :: evt%hadrons)
       case (HADRONS_PYTHIA8)
          allocate (hadrons_pythia8_t :: evt%hadrons)
       case default
          call msg_fatal ('Hadronization: Method ' // &
            char (var_list%get_sval (var_str ("hadronization_method"))) // &
            'not implemented!')
       end select
       call evt%hadrons%init &
            (shower_settings, hadron_settings, fallback_model)
    end select
  end subroutine dispatch_evt_hadrons

  subroutine dispatch_eio (eio, method, var_list, fallback_model, &
         event_callback)
    class(eio_t), allocatable, intent(inout) :: eio
    type(string_t), intent(in) :: method
    type(var_list_t), intent(in) :: var_list
    type(model_t), target, intent(in) :: fallback_model
    class(event_callback_t), allocatable, intent(in) :: event_callback
    logical :: check, keep_beams, keep_remnants, recover_beams
    logical :: use_alphas_from_file, use_scale_from_file
    logical :: fixed_order_nlo_events
    logical :: write_sqme_prc, write_sqme_ref, write_sqme_alt
    logical :: output_cross_section, ensure_order
    type(string_t) :: lhef_version, lhef_extension, raw_version
    type(string_t) :: extension_default, debug_extension, dump_extension, &
         extension_hepmc, &
         extension_lha, extension_hepevt, extension_ascii_short, &
         extension_ascii_long, extension_athena, extension_mokka, &
         extension_stdhep, extension_stdhep_up, extension_stdhep_ev4, &
         extension_raw, extension_hepevt_verb, extension_lha_verb, &
         extension_lcio
    integer :: checkpoint
    integer :: lcio_run_id, hepmc3_mode
    logical :: show_process, show_transforms, show_decay, verbose, pacified
    logical :: dump_weights, dump_compressed, dump_summary, dump_screen
    logical :: proc_as_run_id, hepmc3_write_flows
    keep_beams = &
         var_list%get_lval (var_str ("?keep_beams"))
    keep_remnants = &
         var_list%get_lval (var_str ("?keep_remnants"))
    ensure_order = &
         var_list%get_lval (var_str ("?hepevt_ensure_order"))
    recover_beams = &
         var_list%get_lval (var_str ("?recover_beams"))
    use_alphas_from_file = &
         var_list%get_lval (var_str ("?use_alphas_from_file"))
    use_scale_from_file = &
         var_list%get_lval (var_str ("?use_scale_from_file"))
    fixed_order_nlo_events = &
         var_list%get_lval (var_str ("?fixed_order_nlo_events"))
    select case (char (method))
    case ("raw")
       allocate (eio_raw_t :: eio)
       select type (eio)
       type is (eio_raw_t)
          check = &
               var_list%get_lval (var_str ("?check_event_file"))
          raw_version = &
               var_list%get_sval (var_str ("$event_file_version"))
          extension_raw = &
               var_list%get_sval (var_str ("$extension_raw"))
          call eio%set_parameters (check, use_alphas_from_file, &
               use_scale_from_file, fixed_order_nlo_events, &
               raw_version, extension_raw)
       end select
    case ("checkpoint")
       allocate (eio_checkpoints_t :: eio)
       select type (eio)
       type is (eio_checkpoints_t)
          checkpoint = &
               var_list%get_ival (var_str ("checkpoint"))
          pacified = &
               var_list%get_lval (var_str ("?pacify"))
          call eio%set_parameters (checkpoint, blank = pacified)
       end select
    case ("callback")
       allocate (eio_callback_t :: eio)
       select type (eio)
       type is (eio_callback_t)
          checkpoint = &
               var_list%get_ival (var_str ("event_callback_interval"))
          if (allocated (event_callback)) then
             call eio%set_parameters (event_callback, checkpoint)
          else
             call eio%set_parameters (event_callback_nop_t (), 0)
          end if
       end select
    case ("lhef")
       allocate (eio_lhef_t :: eio)
       select type (eio)
       type is (eio_lhef_t)
          lhef_version = &
               var_list%get_sval (var_str ("$lhef_version"))
          lhef_extension = &
               var_list%get_sval (var_str ("$lhef_extension"))
          write_sqme_prc = &
               var_list%get_lval (var_str ("?lhef_write_sqme_prc"))
          write_sqme_ref = &
               var_list%get_lval (var_str ("?lhef_write_sqme_ref"))
          write_sqme_alt = &
               var_list%get_lval (var_str ("?lhef_write_sqme_alt"))
          call eio%set_parameters ( &
               keep_beams, keep_remnants, recover_beams, &
               use_alphas_from_file, use_scale_from_file, &
               char (lhef_version), lhef_extension, &
               write_sqme_ref, write_sqme_prc, write_sqme_alt)
       end select
    case ("hepmc")
       allocate (eio_hepmc_t :: eio)
       select type (eio)
       type is (eio_hepmc_t)
          output_cross_section = &
               var_list%get_lval (var_str ("?hepmc_output_cross_section"))
          extension_hepmc = &
               var_list%get_sval (var_str ("$extension_hepmc"))
          hepmc3_write_flows = &
               var_list%get_lval (var_str ("?hepmc3_write_flows"))
          select case (char (var_list%get_sval (var_str ("$hepmc3_mode"))))
          case ("HepMC2")
             hepmc3_mode = HEPMC3_MODE_HEPMC2
          case ("HepMC3")
             hepmc3_mode = HEPMC3_MODE_HEPMC3
          case ("Root")
             hepmc3_mode = HEPMC3_MODE_ROOT
             if (extension_hepmc /= "root") then
                call msg_message ("Events: HepMC3 Root mode, using " // &
                     "event sample extension 'root'")
                extension_hepmc = "root"
             end if
          case ("RootTree")
             hepmc3_mode = HEPMC3_MODE_ROOTTREE
             if (extension_hepmc /= "root") then
                call msg_message ("Events: HepMC3 RootTree mode, using " // &
                     "event sample extension 'root'")
                extension_hepmc = "root"
             end if
          case ("HepEVT")
             hepmc3_mode = HEPMC3_MODE_HEPEVT
          case default
             call msg_fatal ("Only supported HepMC3 modes are: 'HepMC2', " // &
                  "'HepMC3', 'HepEVT', 'Root', and 'RootTree'.")
          end select
          call eio%set_parameters (recover_beams, &
               use_alphas_from_file, use_scale_from_file, &
               extension = extension_hepmc, &
               output_cross_section = output_cross_section, &
               hepmc3_mode = hepmc3_mode, &
               hepmc3_write_flows = hepmc3_write_flows)
       end select
    case ("lcio")
       allocate (eio_lcio_t :: eio)
       select type (eio)
       type is (eio_lcio_t)
          extension_lcio = &
               var_list%get_sval (var_str ("$extension_lcio"))
          proc_as_run_id = &
               var_list%get_lval (var_str ("?proc_as_run_id"))
          lcio_run_id = &
               var_list%get_ival (var_str ("lcio_run_id"))
          call eio%set_parameters (recover_beams, &
               use_alphas_from_file, use_scale_from_file, &
               extension_lcio, proc_as_run_id = proc_as_run_id, &
               lcio_run_id = lcio_run_id)
       end select
    case ("stdhep")
       allocate (eio_stdhep_hepevt_t :: eio)
       select type (eio)
       type is (eio_stdhep_hepevt_t)
          extension_stdhep = &
               var_list%get_sval (var_str ("$extension_stdhep"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, recover_beams, &
                use_alphas_from_file, use_scale_from_file, extension_stdhep)
       end select
    case ("stdhep_up")
       allocate (eio_stdhep_hepeup_t :: eio)
       select type (eio)
       type is (eio_stdhep_hepeup_t)
          extension_stdhep_up = &
               var_list%get_sval (var_str ("$extension_stdhep_up"))
          call eio%set_parameters (keep_beams, keep_remnants, ensure_order, &
               recover_beams, use_alphas_from_file, &
               use_scale_from_file, extension_stdhep_up)
       end select
    case ("stdhep_ev4")
       allocate (eio_stdhep_hepev4_t :: eio)
       select type (eio)
       type is (eio_stdhep_hepev4_t)
          extension_stdhep_ev4 = &
               var_list%get_sval (var_str ("$extension_stdhep_ev4"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, recover_beams, &
                use_alphas_from_file, use_scale_from_file, extension_stdhep_ev4)
       end select
    case ("ascii")
       allocate (eio_ascii_ascii_t :: eio)
       select type (eio)
       type is (eio_ascii_ascii_t)
          extension_default = &
               var_list%get_sval (var_str ("$extension_default"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, extension_default)
       end select
    case ("athena")
       allocate (eio_ascii_athena_t :: eio)
       select type (eio)
       type is (eio_ascii_athena_t)
          extension_athena = &
               var_list%get_sval (var_str ("$extension_athena"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, extension_athena)
       end select
    case ("debug")
       allocate (eio_ascii_debug_t :: eio)
       select type (eio)
       type is (eio_ascii_debug_t)
          debug_extension = &
               var_list%get_sval (var_str ("$debug_extension"))
          show_process = &
               var_list%get_lval (var_str ("?debug_process"))
          show_transforms = &
               var_list%get_lval (var_str ("?debug_transforms"))
          show_decay = &
               var_list%get_lval (var_str ("?debug_decay"))
          verbose = &
               var_list%get_lval (var_str ("?debug_verbose"))
          call eio%set_parameters ( &
               extension = debug_extension, &
               show_process = show_process, &
               show_transforms = show_transforms, &
               show_decay = show_decay, &
               verbose = verbose)
       end select
    case ("dump")
       allocate (eio_dump_t :: eio)
       select type (eio)
       type is (eio_dump_t)
          dump_extension = &
               var_list%get_sval (var_str ("$dump_extension"))
          pacified = &
               var_list%get_lval (var_str ("?pacify"))
          dump_weights = &
               var_list%get_lval (var_str ("?dump_weights"))
          dump_compressed = &
               var_list%get_lval (var_str ("?dump_compressed"))
          dump_summary = &
               var_list%get_lval (var_str ("?dump_summary"))
          dump_screen = &
               var_list%get_lval (var_str ("?dump_screen"))
          call eio%set_parameters ( &
               extension = dump_extension, &
               pacify = pacified, &
               weights = dump_weights, &
               compressed = dump_compressed, &
               summary = dump_summary, &
               screen = dump_screen)
       end select
    case ("hepevt")
       allocate (eio_ascii_hepevt_t :: eio)
       select type (eio)
       type is (eio_ascii_hepevt_t)
          extension_hepevt = &
               var_list%get_sval (var_str ("$extension_hepevt"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, extension_hepevt)
       end select
    case ("hepevt_verb")
       allocate (eio_ascii_hepevt_verb_t :: eio)
       select type (eio)
       type is (eio_ascii_hepevt_verb_t)
          extension_hepevt_verb = &
               var_list%get_sval (var_str ("$extension_hepevt_verb"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, extension_hepevt_verb)
       end select
    case ("lha")
       allocate (eio_ascii_lha_t :: eio)
       select type (eio)
       type is (eio_ascii_lha_t)
          extension_lha = &
               var_list%get_sval (var_str ("$extension_lha"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, extension_lha)
       end select
    case ("lha_verb")
       allocate (eio_ascii_lha_verb_t :: eio)
       select type (eio)
       type is (eio_ascii_lha_verb_t)
          extension_lha_verb = var_list%get_sval ( &
               var_str ("$extension_lha_verb"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, extension_lha_verb)
       end select
    case ("long")
       allocate (eio_ascii_long_t :: eio)
       select type (eio)
       type is (eio_ascii_long_t)
          extension_ascii_long = &
               var_list%get_sval (var_str ("$extension_ascii_long"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, extension_ascii_long)
       end select
    case ("mokka")
       allocate (eio_ascii_mokka_t :: eio)
       select type (eio)
       type is (eio_ascii_mokka_t)
          extension_mokka = &
               var_list%get_sval (var_str ("$extension_mokka"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, extension_mokka)
       end select
    case ("short")
       allocate (eio_ascii_short_t :: eio)
       select type (eio)
       type is (eio_ascii_short_t)
          extension_ascii_short = &
               var_list%get_sval (var_str ("$extension_ascii_short"))
          call eio%set_parameters &
               (keep_beams, keep_remnants, ensure_order, extension_ascii_short)
       end select
    case ("weight_stream")
       allocate (eio_weights_t :: eio)
       select type (eio)
       type is (eio_weights_t)
          pacified = &
               var_list%get_lval (var_str ("?pacify"))
          call eio%set_parameters (pacify = pacified)
       end select
    case default
       call msg_fatal ("Event I/O method '" // char (method) &
            // "' not implemented")
    end select
    call eio%set_fallback_model (fallback_model)
  end subroutine dispatch_eio


end module dispatch_transforms
