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

submodule (phs_wood) phs_wood_s

  use io_units
  use constants
  use numeric_utils
  use diagnostics
  use physics_defs
  use md5
  use process_constants
  use sf_mappings
  use sf_base

  implicit none

contains

  module subroutine phs_wood_config_final (object)
    class(phs_wood_config_t), intent(inout) :: object
    logical :: opened
    if (object%io_unit /= 0) then
       inquire (unit = object%io_unit, opened = opened)
       if (opened)  close (object%io_unit)
    end if
    call object%clear_phase_space ()
    call object%forest%final ()
  end subroutine phs_wood_config_final

  module subroutine phs_wood_config_increase_n_par (phs_config)
    class(phs_wood_config_t), intent(inout) :: phs_config
    if (phs_config%is_combined_integration) then
       phs_config%n_par = phs_config%n_par + 3
    end if
  end subroutine phs_wood_config_increase_n_par

  module subroutine phs_wood_config_write (object, unit, include_id)
    class(phs_wood_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: include_id
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)") &
         "Partonic phase-space configuration (phase-space forest):"
    call object%base_write (unit)
    write (u, "(1x,A)")    "Phase-space configuration parameters:"
    call object%par%write (u)
    call object%mapping_defaults%write (u)
    write (u, "(3x,A,A,A)")  "Run ID: '", char (object%run_id), "'"
  end subroutine phs_wood_config_write

  module subroutine phs_wood_config_write_forest (object, unit)
    class(phs_wood_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    call object%forest%write (u)
  end subroutine phs_wood_config_write_forest

  module subroutine phs_wood_config_set_parameters (phs_config, par)
    class(phs_wood_config_t), intent(inout) :: phs_config
    type(phs_parameters_t), intent(in) :: par
    phs_config%par = par
  end subroutine phs_wood_config_set_parameters

  module subroutine phs_wood_config_enable_equivalences (phs_config)
    class(phs_wood_config_t), intent(inout) :: phs_config
    phs_config%use_equivalences = .true.
  end subroutine phs_wood_config_enable_equivalences

  module subroutine phs_wood_config_set_mapping_defaults &
       (phs_config, mapping_defaults)
    class(phs_wood_config_t), intent(inout) :: phs_config
    type(mapping_defaults_t), intent(in) :: mapping_defaults
    phs_config%mapping_defaults = mapping_defaults
  end subroutine phs_wood_config_set_mapping_defaults

  module subroutine phs_wood_config_set_input (phs_config, unit)
    class(phs_wood_config_t), intent(inout) :: phs_config
    integer, intent(in) :: unit
    phs_config%io_unit = unit
    rewind (unit)
  end subroutine phs_wood_config_set_input

  module subroutine phs_wood_config_generate_phase_space (phs_config)
    class(phs_wood_config_t), intent(inout) :: phs_config
    integer :: off_shell, extra_off_shell
    logical :: valid
    integer :: unit_fds
    type(string_t) :: file_name
    logical :: file_exists
    call msg_message ("Phase space: generating configuration ...")
    off_shell = phs_config%par%off_shell
    if (phs_config%use_cascades2) then
       file_name = char (phs_config%id) // ".fds"
       inquire (file=char (file_name), exist=file_exists)
       if (.not. file_exists) call msg_fatal &
            ("The O'Mega input file " // char (file_name) // &
            " does not exist. " // "Please make sure that the " // &
            "variable ?omega_write_phs_output has been set correctly.")
       unit_fds = free_unit ()
       open (unit=unit_fds, file=char(file_name), status='old', action='read')
       do extra_off_shell = 0, max (phs_config%n_tot - 3, 0)
          phs_config%par%off_shell = off_shell + extra_off_shell
          allocate (phs_config%feyngraph_set)
          call feyngraph_set_generate (phs_config%feyngraph_set, &
               phs_config%model, phs_config%n_in, phs_config%n_out, &
               phs_config%flv, &
               phs_config%par, phs_config%fatal_beam_decay, unit_fds, &
               phs_config%vis_channels)
          if (feyngraph_set_is_valid (phs_config%feyngraph_set)) then
             exit
          else
             call msg_message ("Phase space: ... failed.  &
                  &Increasing phs_off_shell ...")
             call phs_config%feyngraph_set%final ()
             deallocate (phs_config%feyngraph_set)
          end if
       end do
       close (unit_fds)
    else
       allocate (phs_config%cascade_set)
       do extra_off_shell = 0, max (phs_config%n_tot - 3, 0)
          phs_config%par%off_shell = off_shell + extra_off_shell
          call cascade_set_generate (phs_config%cascade_set, &
               phs_config%model, phs_config%n_in, phs_config%n_out, &
               phs_config%flv, &
               phs_config%par, phs_config%fatal_beam_decay)
          if (cascade_set_is_valid (phs_config%cascade_set)) then
             exit
          else
             call msg_message ("Phase space: ... failed.  &
                  &Increasing phs_off_shell ...")
          end if
       end do
    end if
    if (phs_config%use_cascades2) then
       valid = feyngraph_set_is_valid (phs_config%feyngraph_set)
    else
       valid = cascade_set_is_valid (phs_config%cascade_set)
    end if
    if (valid) then
       call msg_message ("Phase space: ... success.")
    else
       call msg_fatal ("Phase-space: generation failed")
    end if
  end subroutine phs_wood_config_generate_phase_space

  module subroutine phs_wood_config_write_phase_space (phs_config, &
       filename_vis, unit)
    class(phs_wood_config_t), intent(in) :: phs_config
    integer, intent(in), optional :: unit
    type(string_t), intent(in), optional :: filename_vis
    type(string_t) :: setenv_tex, setenv_mp, pipe, pipe_dvi
    integer :: u, unit_tex, unit_dev, status
    if (allocated (phs_config%cascade_set) .or. &
         allocated (phs_config%feyngraph_set)) then
       if (present (unit)) then
          u = unit
       else
          u = phs_config%io_unit
       end if
       write (u, "(1x,A,A)") "process ", char (phs_config%id)
       write (u, "(A)")
       if (phs_config%use_cascades2) then
          call feyngraph_set_write_process_bincode_format (phs_config%feyngraph_set, u)
       else
          call cascade_set_write_process_bincode_format (phs_config%cascade_set, u)
       end if
       write (u, "(A)")
       write (u, "(3x,A,A,A32,A)") "md5sum_process    = ", &
            '"', phs_config%md5sum_process, '"'
       write (u, "(3x,A,A,A32,A)") "md5sum_model_par  = ", &
            '"', phs_config%md5sum_model_par, '"'
       write (u, "(3x,A,A,A32,A)") "md5sum_phs_config = ", &
            '"', phs_config%md5sum_phs_config, '"'
       call phs_config%par%write (u)
       if (phs_config%use_cascades2) then
          call feyngraph_set_write_file_format (phs_config%feyngraph_set, u)
       else
          call cascade_set_write_file_format (phs_config%cascade_set, u)
       end if
       if (phs_config%vis_channels) then
          unit_tex = free_unit ()
          open (unit=unit_tex, file=char(filename_vis // ".tex"), &
               action="write", status="replace")
          if (phs_config%use_cascades2) then
             call feyngraph_set_write_graph_format (phs_config%feyngraph_set, &
                  filename_vis // "-graphs", phs_config%id, unit_tex)
          else
             call cascade_set_write_graph_format (phs_config%cascade_set, &
                  filename_vis // "-graphs", phs_config%id, unit_tex)
          end if
          close (unit_tex)
          call msg_message ("Phase space: visualizing channels in file " &
               // char(trim(filename_vis)) // "...")
          if (phs_config%os_data%event_analysis_ps) then
             BLOCK: do
                unit_dev = free_unit ()
                open (file = "/dev/null", unit = unit_dev, &
                     action = "write", iostat = status)
                if (status /= 0) then
                   pipe = ""
                   pipe_dvi = ""
                else
                   pipe = " > /dev/null"
                   pipe_dvi = " 2>/dev/null 1>/dev/null"
                end if
                close (unit_dev)
                if (phs_config%os_data%whizard_texpath /= "") then
                   setenv_tex = "TEXINPUTS=" // &
                        phs_config%os_data%whizard_texpath // ":$TEXINPUTS "
                   setenv_mp = "MPINPUTS=" // &
                        phs_config%os_data%whizard_texpath // ":$MPINPUTS "
                else
                   setenv_tex = ""
                   setenv_mp = ""
                end if
                call os_system_call (setenv_tex // &
                     phs_config%os_data%latex // " " // &
                     filename_vis // ".tex " // pipe, status)
                if (status /= 0)  exit BLOCK
                if (phs_config%os_data%mpost /= "") then
                   call os_system_call (setenv_mp // &
                        phs_config%os_data%mpost // " " // &
                        filename_vis // "-graphs.mp" // pipe, status)
                else
                   call msg_fatal ("Could not use MetaPOST.")
                end if
                if (status /= 0)  exit BLOCK
                call os_system_call (setenv_tex // &
                     phs_config%os_data%latex // " " // &
                     filename_vis // ".tex" // pipe, status)
                if (status /= 0)  exit BLOCK
                call os_system_call &
                     (phs_config%os_data%dvips // " -o " // filename_vis &
                     // ".ps " // filename_vis // ".dvi" // pipe_dvi, status)
                if (status /= 0)  exit BLOCK
                if (phs_config%os_data%event_analysis_pdf) then
                   call os_system_call (phs_config%os_data%ps2pdf // " " // &
                        filename_vis // ".ps", status)
                   if (status /= 0)  exit BLOCK
                end if
                exit BLOCK
             end do BLOCK
             if (status /= 0) then
                call msg_error ("Unable to compile analysis output file")
             end if
          end if
       end if
    else
       call msg_fatal ("Phase-space configuration: &
            &no phase space object generated")
    end if
  end subroutine phs_wood_config_write_phase_space

  module subroutine phs_wood_config_clear_phase_space (phs_config)
    class(phs_wood_config_t), intent(inout) :: phs_config
    if (allocated (phs_config%cascade_set)) then
       call cascade_set_final (phs_config%cascade_set)
       deallocate (phs_config%cascade_set)
    end if
    if (allocated (phs_config%feyngraph_set)) then
       call phs_config%feyngraph_set%final ()
       deallocate (phs_config%feyngraph_set)
    end if
  end subroutine phs_wood_config_clear_phase_space

  module subroutine phs_wood_config_extract_resonance_history_set &
       (phs_config, res_set, include_trivial)
    class(phs_wood_config_t), intent(in) :: phs_config
    type(resonance_history_set_t), intent(out) :: res_set
    logical, intent(in), optional :: include_trivial
    call phs_config%forest%extract_resonance_history_set &
         (res_set, include_trivial)
  end subroutine phs_wood_config_extract_resonance_history_set

  module subroutine phs_wood_config_configure (phs_config, sqrts, &
       sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
       ignore_mismatch, nlo_type, subdir)
    class(phs_wood_config_t), intent(inout) :: phs_config
    real(default), intent(in) :: sqrts
    logical, intent(in), optional :: sqrts_fixed
    logical, intent(in), optional :: lab_is_cm
    logical, intent(in), optional :: azimuthal_dependence
    logical, intent(in), optional :: rebuild
    logical, intent(in), optional :: ignore_mismatch
    integer, intent(in), optional :: nlo_type
    type(string_t), intent(in), optional :: subdir
    type(string_t) :: filename, filename_vis
    logical :: variable_limits
    logical :: ok, exist, found, check, match, rebuild_phs
    integer :: g, c0, c1, n
    if (present (nlo_type)) then
      phs_config%nlo_type = nlo_type
    else
      phs_config%nlo_type = BORN
    end if
    phs_config%sqrts = sqrts
    phs_config%par%sqrts = sqrts
    if (present (sqrts_fixed)) &
         phs_config%sqrts_fixed = sqrts_fixed
    if (present (lab_is_cm)) &
         phs_config%lab_is_cm = lab_is_cm
    if (present (azimuthal_dependence)) &
         phs_config%azimuthal_dependence = azimuthal_dependence
    if (present (rebuild)) then
       rebuild_phs = rebuild
    else
       rebuild_phs = .true.
    end if
    if (present (ignore_mismatch)) then
       check = .not. ignore_mismatch
       if (ignore_mismatch) &
            call msg_warning ("Reading phs file: MD5 sum check disabled")
    else
       check = .true.
    end if
    phs_config%md5sum_forest = ""
    call phs_config%compute_md5sum (include_id = .false.)
    if (phs_config%io_unit == 0) then
       filename = phs_config%make_phs_filename (subdir)
       filename_vis = phs_config%make_phs_filename (subdir) // "-vis"
       if (.not. rebuild_phs) then
          if (check) then
             call phs_config%read_phs_file (exist, found, match, subdir=subdir)
             rebuild_phs = .not. (exist .and. found .and. match)
          else
             call phs_config%read_phs_file (exist, found, subdir=subdir)
             rebuild_phs = .not. (exist .and. found)
          end if
       end if
       if (.not. mpi_is_comm_master ()) then
          rebuild_phs = .false.
          call msg_message ("MPI: Workers do not build phase space configuration.")
       end if
       if (rebuild_phs) then
          call phs_config%generate_phase_space ()
          phs_config%io_unit = free_unit ()
          if (phs_config%id /= "") then
             call msg_message ("Phase space: writing configuration file '" &
                  // char (filename) // "'")
             open (phs_config%io_unit, file = char (filename), &
                  status = "replace", action = "readwrite")
          else
             open (phs_config%io_unit, status = "scratch", action = "readwrite")
          end if
          call phs_config%write_phase_space (filename_vis)
          rewind (phs_config%io_unit)
       else
          call msg_message ("Phase space: keeping configuration file '" &
               // char (filename) // "'")
       end if
    end if
    if (phs_config%io_unit == 0) then
       ok = .true.
    else
      call phs_config%forest%read (phs_config%io_unit, phs_config%id, &
           phs_config%n_in, phs_config%n_out, phs_config%model, ok)
       if (.not. phs_config%io_unit_keep_open) then
          close (phs_config%io_unit)
          phs_config%io_unit = 0
       end if
    end if
    if (ok) then
       call phs_config%forest%set_flavors (phs_config%flv(:,1))
       variable_limits = .not. phs_config%lab_is_cm
       call phs_config%forest%set_parameters (phs_config%mapping_defaults, &
            variable_limits)
       call phs_config%forest%setup_prt_combinations ()
       phs_config%n_channel = phs_config%forest%get_n_channels ()
       phs_config%n_par = phs_config%forest%get_n_parameters ()
       allocate (phs_config%channel (phs_config%n_channel))
       if (phs_config%use_equivalences) then
          call phs_config%forest%set_equivalences ()
          call phs_config%forest%get_equivalences (phs_config%channel, &
               phs_config%azimuthal_dependence)
          phs_config%provides_equivalences = .true.
       end if
       call phs_config%forest%set_s_mappings ()
       call phs_config%record_on_shell ()
       if (phs_config%mapping_defaults%enable_s_mapping) then
          call phs_config%record_s_mappings ()
       end if
       allocate (phs_config%chain (phs_config%n_channel), source = 0)
       do g = 1, phs_config%forest%get_n_groves ()
          call phs_config%forest%get_grove_bounds (g, c0, c1, n)
          phs_config%chain (c0:c1) = g
       end do
       phs_config%provides_chains = .true.
       call phs_config%compute_md5sum_forest ()
    else
       write (msg_buffer, "(A,A,A)") &
            "Phase space: process '", &
            char (phs_config%id), "' not found in configuration file"
       call msg_fatal ()
    end if
  end subroutine phs_wood_config_configure

  module subroutine phs_wood_config_compute_md5sum_forest (phs_config)
    class(phs_wood_config_t), intent(inout) :: phs_config
    integer :: u
    u = free_unit ()
    open (u, status = "scratch", action = "readwrite")
    call phs_config%write_forest (u)
    rewind (u)
    phs_config%md5sum_forest = md5sum (u)
    close (u)
  end subroutine phs_wood_config_compute_md5sum_forest

  module function phs_wood_make_phs_filename &
       (phs_config, subdir) result (filename)
    class(phs_wood_config_t), intent(in) :: phs_config
    type(string_t), intent(in), optional :: subdir
    type(string_t) :: filename
    type(string_t) :: basename, suffix, comp_code, comp_index
    basename = phs_config%id
    call split (basename, suffix, "_", back=.true.)
    comp_code = extract (suffix, 1, 1)
    comp_index = extract (suffix, 2)
    if (comp_code == "i" .and. verify (comp_index, "1234567890") == 0) then
       suffix = "." // comp_code // comp_index
    else
       basename = phs_config%id
       suffix = ""
    end if
    if (phs_config%run_id /= "") then
       filename = basename // "." // phs_config%run_id // suffix // ".phs"
    else
       filename = basename // suffix // ".phs"
    end if
    if (present (subdir)) then
       filename = subdir // "/" // filename
    end if
  end function phs_wood_make_phs_filename

  module subroutine phs_wood_config_reshuffle_flavors &
       (phs_config, reshuffle, flv_extra)
    class(phs_wood_config_t), intent(inout) :: phs_config
    integer, intent(in), dimension(:), allocatable :: reshuffle
    type(flavor_t), intent(in) :: flv_extra
    call phs_config%forest%set_flavors (phs_config%flv(:,1), reshuffle, &
         flv_extra)
  end subroutine phs_wood_config_reshuffle_flavors

  module subroutine phs_wood_config_set_momentum_links (phs_config, reshuffle)
    class(phs_wood_config_t), intent(inout) :: phs_config
    integer, intent(in), dimension(:), allocatable :: reshuffle
    call phs_config%forest%set_momentum_links (reshuffle)
  end subroutine phs_wood_config_set_momentum_links

  module subroutine phs_wood_config_record_s_mappings (phs_config)
    class(phs_wood_config_t), intent(inout) :: phs_config
    logical :: flag
    real(default) :: mass, width
    integer :: c
    do c = 1, phs_config%n_channel
       call phs_config%forest%get_s_mapping (c, flag, mass, width)
       if (flag) then
          if (mass == 0) then
             call msg_fatal ("Phase space: s-channel resonance " &
                  // " has zero mass")
          end if
          if (width == 0) then
             call msg_fatal ("Phase space: s-channel resonance " &
                  // " has zero width")
          end if
          call phs_config%channel(c)%set_resonant (mass, width)
       end if
    end do
  end subroutine phs_wood_config_record_s_mappings

  module subroutine phs_wood_config_record_on_shell (phs_config)
    class(phs_wood_config_t), intent(inout) :: phs_config
    logical :: flag
    real(default) :: mass
    integer :: c
    do c = 1, phs_config%n_channel
       call phs_config%forest%get_on_shell (c, flag, mass)
       if (flag) then
          call phs_config%channel(c)%set_on_shell (mass)
       end if
    end do
  end subroutine phs_wood_config_record_on_shell

  module function phs_wood_config_get_md5sum (phs_config) result (md5sum)
    class(phs_wood_config_t), intent(in) :: phs_config
    character(32) :: md5sum
    if (phs_config%md5sum_forest /= "") then
       md5sum = phs_config%md5sum_forest
    else
       md5sum = phs_config%md5sum_phs_config
    end if
  end function phs_wood_config_get_md5sum

  module subroutine phs_wood_read_phs_file &
       (phs_config, exist, found, match, subdir)
    class(phs_wood_config_t), intent(inout) :: phs_config
    logical, intent(out) :: exist
    logical, intent(out) :: found
    logical, intent(out), optional :: match
    type(string_t), intent(in), optional :: subdir
    type(string_t) :: filename
    integer :: u
    filename = phs_config%make_phs_filename (subdir)
    inquire (file = char (filename), exist = exist)
    if (exist) then
       u = free_unit ()
       open (u, file = char (filename), action = "read", status = "old")
       call phs_config%forest%read (u, phs_config%id, phs_config%n_in, &
            phs_config%n_out, phs_config%model, found, &
            phs_config%md5sum_process, phs_config%md5sum_model_par, &
            phs_config%md5sum_phs_config, match = match)
       close (u)
    else
       found = .false.
       if (present (match))  match = .false.
    end if
  end subroutine phs_wood_read_phs_file

  module subroutine phs_wood_config_startup_message (phs_config, unit)
    class(phs_wood_config_t), intent(in) :: phs_config
    integer, intent(in), optional :: unit
    integer :: n_groves, n_eq
    n_groves = phs_config%forest%get_n_groves ()
    n_eq = phs_config%forest%get_n_equivalences ()
    call phs_config%base_startup_message (unit)
    if (phs_config%n_channel == 1) then
       write (msg_buffer, "(A,2(I0,A))") &
            "Phase space: found ", phs_config%n_channel, &
            " channel, collected in ", n_groves, &
            " grove."
    else if (n_groves == 1) then
       write (msg_buffer, "(A,2(I0,A))") &
            "Phase space: found ", phs_config%n_channel, &
            " channels, collected in ", n_groves, &
            " grove."
       else
       write (msg_buffer, "(A,2(I0,A))") &
            "Phase space: found ", phs_config%n_channel, &
            " channels, collected in ", n_groves, &
            " groves."
    end if
    call msg_message (unit = unit)
    if (phs_config%use_equivalences) then
       if (n_eq == 1) then
          write (msg_buffer, "(A,I0,A)") &
               "Phase space: Using ", n_eq, &
               " equivalence between channels."
       else
          write (msg_buffer, "(A,I0,A)") &
               "Phase space: Using ", n_eq, &
               " equivalences between channels."
       end if
    else
       write (msg_buffer, "(A)") &
            "Phase space: no equivalences between channels used."
    end if
    call msg_message (unit = unit)
    write (msg_buffer, "(A,2(1x,I0,1x,A))") &
         "Phase space: wood"
    call msg_message (unit = unit)
  end subroutine phs_wood_config_startup_message

  module subroutine phs_wood_write (object, unit, verbose)
    class(phs_wood_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    u = given_output_unit (unit)
    call object%base_write (u)
  end subroutine phs_wood_write

  module subroutine phs_wood_write_forest (object, unit)
    class(phs_wood_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    call object%forest%write (u)
  end subroutine phs_wood_write_forest

  module subroutine phs_wood_final (object)
    class(phs_wood_t), intent(inout) :: object
    call object%forest%final ()
  end subroutine phs_wood_final

  module subroutine phs_wood_init (phs, phs_config)
    class(phs_wood_t), intent(out) :: phs
    class(phs_config_t), intent(in), target :: phs_config
    call phs%base_init (phs_config)
    select type(phs_config)
    type is (phs_wood_config_t)
       phs%forest = phs_config%forest
       if (phs_config%is_combined_integration) then
          phs%n_r_born = phs_config%n_par - 3
       end if
    end select
  end subroutine phs_wood_init

  module subroutine phs_wood_evaluate_selected_channel (phs, c_in, r_in)
    class(phs_wood_t), intent(inout) :: phs
    integer, intent(in) :: c_in
    real(default), intent(in), dimension(:) :: r_in
    logical :: ok
    phs%q_defined = .false.
    if (phs%p_defined) then
       call phs%forest%set_prt_in (phs%p)
       phs%r(:,c_in) = r_in
       call phs%forest%evaluate_selected_channel (c_in, phs%active_channel, &
            phs%sqrts_hat, phs%r, phs%f, phs%volume, ok)
       select type (config => phs%config)
       type is (phs_wood_config_t)
          if (config%is_combined_integration) then
             if (phs%n_r_born >= 0) then
                phs%r_real = r_in (phs%n_r_born + 1 : phs%n_r_born + 3)
             else
                call msg_fatal ("n_r_born should be larger than 0!")
             end if
          end if
       end select
       if (ok) then
          phs%q = phs%forest%get_momenta_out ()
          phs%q_defined = .true.
       end if
    end if
  end subroutine phs_wood_evaluate_selected_channel

  module subroutine phs_wood_evaluate_other_channels (phs, c_in)
    class(phs_wood_t), intent(inout) :: phs
    integer, intent(in) :: c_in
    integer :: c
    if (phs%q_defined) then
       call phs%forest%evaluate_other_channels (c_in, phs%active_channel, &
            phs%sqrts_hat, phs%r, phs%f, combine=.true.)
       select type (config => phs%config)
       type is (phs_wood_config_t)
          if (config%is_combined_integration) then
             if (phs%n_r_born >= 0) then
                do c = 1, size (phs%r, 2)
                   phs%r(phs%n_r_born + 1 : phs%n_r_born + 3, c) = phs%r_real
                end do
             else
                phs%r_defined = .false.
             end if
          end if
       end select
       phs%r_defined = .true.
    end if
  end subroutine phs_wood_evaluate_other_channels

  module subroutine phs_wood_inverse (phs)
    class(phs_wood_t), intent(inout) :: phs
    if (phs%p_defined .and. phs%q_defined) then
       call phs%forest%set_prt_in (phs%p)
       call phs%forest%set_prt_out (phs%q)
       call phs%forest%recover_channel (1, phs%sqrts_hat, phs%r, &
            phs%f, phs%volume)
       call phs%forest%evaluate_other_channels (1, phs%active_channel, &
            phs%sqrts_hat, phs%r, phs%f, combine=.false.)
       phs%r_defined = .true.
    end if
  end subroutine phs_wood_inverse


end submodule phs_wood_s

