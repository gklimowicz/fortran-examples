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

submodule (process_config) process_config_s

  use format_utils, only: write_separator
  use io_units
  use diagnostics
  use md5
  use physics_defs
  use helicities
  use colors
  use quantum_numbers
  use state_matrices
  use prc_external
  use prc_openloops, only: prc_openloops_t
  use prc_threshold, only: prc_threshold_t
  use blha_olp_interfaces, only: prc_blha_t

  implicit none

contains

  module function flagged (v_list, id, def) result (flag)
    logical :: flag
    integer, dimension(:), intent(in) :: v_list
    integer, intent(in) :: id
    logical, intent(in), optional :: def
    logical :: default_result
    default_result = .false.;  if (present (def))  default_result = def
    if (default_result) then
       flag = all (v_list /= -id)
    else
       flag = all (v_list /= -id) .and. any (v_list == id)
    end if
  end function flagged

  module subroutine set_flag (v_list, value, flag)
    integer, dimension(:), intent(inout), allocatable :: v_list
    integer, intent(in) :: value
    logical, intent(in), optional :: flag
    if (present (flag)) then
       if (flag) then
          v_list = [v_list, value]
       else
          v_list = [v_list, -value]
       end if
    end if
  end subroutine set_flag

  module subroutine process_config_data_write &
       (config, u, counters, model, expressions)
    class(process_config_data_t), intent(in) :: config
    integer, intent(in) :: u
    logical, intent(in) :: counters
    logical, intent(in) :: model
    logical, intent(in) :: expressions
    write (u, "(1x,A)") "Configuration data:"
    if (counters) then
       write (u, "(3x,A,I0)") "Number of incoming particles = ", &
            config%n_in
       write (u, "(3x,A,I0)") "Number of process components = ", &
            config%n_components
       write (u, "(3x,A,I0)") "Number of process terms      = ", &
            config%n_terms
       write (u, "(3x,A,I0)") "Number of MCI configurations = ", &
            config%n_mci
    end if
    if (associated (config%model)) then
       write (u, "(3x,A,A)")  "Model = ", char (config%model_name)
       if (model) then
          call write_separator (u)
          call config%model%write (u)
          call write_separator (u)
       end if
    else
       write (u, "(3x,A,A,A)")  "Model = ", char (config%model_name), &
            " [not associated]"
    end if
    call config%qcd%write (u, show_md5sum = .false.)
    call write_separator (u)
    if (expressions) then
       if (allocated (config%ef_cuts)) then
          call write_separator (u)
          write (u, "(3x,A)") "Cut expression:"
          call config%ef_cuts%write (u)
       end if
       if (allocated (config%ef_scale)) then
          call write_separator (u)
          write (u, "(3x,A)") "Scale expression:"
          call config%ef_scale%write (u)
       end if
       if (allocated (config%ef_fac_scale)) then
          call write_separator (u)
          write (u, "(3x,A)") "Factorization scale expression:"
          call config%ef_fac_scale%write (u)
       end if
       if (allocated (config%ef_ren_scale)) then
          call write_separator (u)
          write (u, "(3x,A)") "Renormalization scale expression:"
          call config%ef_ren_scale%write (u)
       end if
       if (allocated (config%ef_weight)) then
          call write_separator (u)
          write (u, "(3x,A)") "Weight expression:"
          call config%ef_weight%write (u)
       end if
    else
       call write_separator (u)
       write (u, "(3x,A)") "Expressions (cut, scales, weight): [not shown]"
    end if
    if (config%md5sum /= "") then
       call write_separator (u)
       write (u, "(3x,A,A,A)")  "MD5 sum (config)  = '", config%md5sum, "'"
    end if
  end subroutine process_config_data_write

  module function process_config_data_get_qcd (config) result (qcd)
    class(process_config_data_t), intent(in) :: config
    type(qcd_t) :: qcd
    qcd = config%qcd
  end function process_config_data_get_qcd

  module subroutine process_config_data_compute_md5sum (config)
    class(process_config_data_t), intent(inout) :: config
    integer :: u
    if (config%md5sum == "") then
       u = free_unit ()
       open (u, status = "scratch", action = "readwrite")
       call config%write (u, counters = .false., &
            model = .true., expressions = .true.)
       rewind (u)
       config%md5sum = md5sum (u)
       close (u)
    end if
  end subroutine process_config_data_compute_md5sum

  pure module function process_config_data_get_md5sum (config) result (md5)
    character(32) :: md5
    class(process_config_data_t), intent(in) :: config
    md5 = config%md5sum
  end function process_config_data_get_md5sum

  module subroutine process_environment_final (env)
    class(process_environment_t), intent(inout) :: env
    if (associated (env%model)) then
       call env%model%final ()
       deallocate (env%model)
    end if
    if (associated (env%var_list)) then
       call env%var_list%final (follow_link=.true.)
       deallocate (env%var_list)
    end if
  end subroutine process_environment_final

  module subroutine process_environment_write (env, unit, &
       show_var_list, show_model, show_lib, show_beams, show_os_data)
    class(process_environment_t), intent(in) :: env
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_var_list
    logical, intent(in), optional :: show_model
    logical, intent(in), optional :: show_lib
    logical, intent(in), optional :: show_beams
    logical, intent(in), optional :: show_os_data
    integer :: u, iostat
    integer, dimension(:), allocatable :: v_list
    character(0) :: iomsg
    u = given_output_unit (unit)
    allocate (v_list (0))
    call set_flag (v_list, F_SHOW_VAR_LIST, show_var_list)
    call set_flag (v_list, F_SHOW_MODEL, show_model)
    call set_flag (v_list, F_SHOW_LIB, show_lib)
    call set_flag (v_list, F_SHOW_BEAMS, show_beams)
    call set_flag (v_list, F_SHOW_OS_DATA, show_os_data)
    call env%write_formatted (u, "LISTDIRECTED", v_list, iostat, iomsg)
  end subroutine process_environment_write

  module subroutine process_environment_write_formatted &
       (dtv, unit, iotype, v_list, iostat, iomsg)
    class(process_environment_t), intent(in) :: dtv
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, dimension(:), intent(in) :: v_list
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    associate (env => dtv)
      if (flagged (v_list, F_SHOW_VAR_LIST, .true.)) then
         write (unit, "(1x,A)")  "Variable list:"
         if (associated (env%var_list)) then
            call write_separator (unit)
            call env%var_list%write (unit)
         else
            write (unit, "(3x,A)")  "[not allocated]"
         end if
         call write_separator (unit)
      end if
      if (flagged (v_list, F_SHOW_MODEL, .true.)) then
         write (unit, "(1x,A)")  "Model:"
         if (associated (env%model)) then
            call write_separator (unit)
            call env%model%write (unit)
         else
            write (unit, "(3x,A)")  "[not allocated]"
         end if
         call write_separator (unit)
      end if
      if (flagged (v_list, F_SHOW_LIB, .true.)) then
         write (unit, "(1x,A)")  "Process library:"
         if (associated (env%lib)) then
            call write_separator (unit)
            call env%lib%write (unit)
         else
            write (unit, "(3x,A)")  "[not allocated]"
         end if
      end if
      if (flagged (v_list, F_SHOW_BEAMS, .true.)) then
         call write_separator (unit)
         call env%beam_structure%write (unit)
      end if
      if (flagged (v_list, F_SHOW_OS_DATA, .true.)) then
         write (unit, "(1x,A)")  "Operating-system data:"
         call write_separator (unit)
         call env%os_data%write (unit)
      end if
    end associate
    iostat = 0
  end subroutine process_environment_write_formatted

  module subroutine process_environment_init &
       (env, model, lib, os_data, var_list, beam_structure)
    class(process_environment_t), intent(out) :: env
    type(model_t), intent(in), target :: model
    type(process_library_t), intent(in), target :: lib
    type(os_data_t), intent(in) :: os_data
    type(var_list_t), intent(in), target, optional :: var_list
    type(beam_structure_t), intent(in), optional :: beam_structure
    allocate (env%model)
    call env%model%init_instance (model)
    env%lib => lib
    env%os_data = os_data
    allocate (env%var_list)
    if (present (var_list)) then
       call env%var_list%init_snapshot (var_list, follow_link=.true.)
       env%var_list_is_set = .true.
    end if
    if (present (beam_structure)) then
       env%beam_structure = beam_structure
    end if
  end subroutine process_environment_init

  module function process_environment_got_var_list (env) result (flag)
    class(process_environment_t), intent(in) :: env
    logical :: flag
    flag = env%var_list_is_set
  end function process_environment_got_var_list

  module function process_environment_get_var_list_ptr (env) result (var_list)
    class(process_environment_t), intent(in) :: env
    type(var_list_t), pointer :: var_list
    var_list => env%var_list
  end function process_environment_get_var_list_ptr

  module function process_environment_get_model_ptr (env) result (model)
    class(process_environment_t), intent(in) :: env
    type(model_t), pointer :: model
    model => env%model
  end function process_environment_get_model_ptr

  module function process_environment_get_lib_ptr (env) result (lib)
    class(process_environment_t), intent(inout) :: env
    type(process_library_t), pointer :: lib
    lib => env%lib
  end function process_environment_get_lib_ptr

  module subroutine process_environment_reset_lib_ptr (env)
    class(process_environment_t), intent(inout) :: env
    env%lib => null ()
  end subroutine process_environment_reset_lib_ptr

  module subroutine process_environment_check_lib_sanity (env, meta)
    class(process_environment_t), intent(in) :: env
    type(process_metadata_t), intent(in) :: meta
    if (associated (env%lib)) then
       if (env%lib%get_update_counter () /= meta%lib_update_counter) then
          call msg_fatal ("Process '" // char (meta%id) &
               // "': library has been recompiled after integration")
       end if
    end if
  end subroutine process_environment_check_lib_sanity

  module subroutine process_environment_fill_process_constants &
       (env, id, i_component, data)
    class(process_environment_t), intent(in) :: env
    type(string_t), intent(in) :: id
    integer, intent(in) :: i_component
    type(process_constants_t), intent(out) :: data
    call env%lib%fill_constants (id, i_component, data)
  end subroutine process_environment_fill_process_constants

  module function process_environment_get_beam_structure &
       (env) result (beam_structure)
    class(process_environment_t), intent(in) :: env
    type(beam_structure_t) :: beam_structure
    beam_structure = env%beam_structure
  end function process_environment_get_beam_structure

  module function process_environment_has_pdfs (env) result (flag)
    class(process_environment_t), intent(in) :: env
    logical :: flag
    flag = env%beam_structure%has_pdf ()
  end function process_environment_has_pdfs

  module function process_environment_has_polarized_beams (env) result (flag)
    class(process_environment_t), intent(in) :: env
    logical :: flag
    flag = env%beam_structure%has_polarized_beams ()
  end function process_environment_has_polarized_beams

  module function process_environment_get_os_data (env) result (os_data)
    class(process_environment_t), intent(in) :: env
    type(os_data_t) :: os_data
    os_data = env%os_data
  end function process_environment_get_os_data

  module subroutine process_metadata_write (meta, u, screen)
    class(process_metadata_t), intent(in) :: meta
    integer, intent(in) :: u
    logical, intent(in) :: screen
    integer :: i
    select case (meta%type)
    case (PRC_UNKNOWN)
       if (screen) then
          write (msg_buffer, "(A)") "Process [undefined]"
       else
          write (u, "(1x,A)") "Process [undefined]"
       end if
       return
    case (PRC_DECAY)
       if (screen) then
          write (msg_buffer, "(A,1x,A,A,A)") "Process [decay]:", &
               "'", char (meta%id), "'"
       else
          write (u, "(1x,A)", advance="no") "Process [decay]:"
       end if
    case (PRC_SCATTERING)
       if (screen) then
          write (msg_buffer, "(A,1x,A,A,A)") "Process [scattering]:", &
               "'", char (meta%id), "'"
       else
          write (u, "(1x,A)", advance="no") "Process [scattering]:"
       end if
    case default
       call msg_bug ("process_write: undefined process type")
    end select
    if (screen)  then
       call msg_message ()
    else
       write (u, "(1x,A,A,A)") "'", char (meta%id), "'"
    end if
    if (meta%num_id /= 0) then
       if (screen) then
          write (msg_buffer, "(2x,A,I0)") "ID (num)      = ", meta%num_id
          call msg_message ()
       else
          write (u, "(3x,A,I0)") "ID (num)      = ", meta%num_id
       end if
    end if
    if (screen) then
       if (meta%run_id /= "") then
          write (msg_buffer, "(2x,A,A,A)") "Run ID        = '", &
               char (meta%run_id), "'"
          call msg_message ()
       end if
    else
       write (u, "(3x,A,A,A)") "Run ID        = '", char (meta%run_id), "'"
    end if
    if (allocated (meta%lib_name)) then
       if (screen) then
          write (msg_buffer, "(2x,A,A,A)")  "Library name  = '", &
               char (meta%lib_name), "'"
          call msg_message ()
       else
          write (u, "(3x,A,A,A)")  "Library name  = '", &
               char (meta%lib_name), "'"
       end if
    else
       if (screen) then
          write (msg_buffer, "(2x,A)")  "Library name  = [not associated]"
          call msg_message ()
       else
          write (u, "(3x,A)")  "Library name  = [not associated]"
       end if
    end if
    if (screen) then
       write (msg_buffer, "(2x,A,I0)")  "Process index = ", meta%lib_index
       call msg_message ()
    else
       write (u, "(3x,A,I0)")  "Process index = ", meta%lib_index
    end if
    if (allocated (meta%component_id)) then
       if (screen) then
          if (any (meta%active)) then
             write (msg_buffer, "(2x,A)")  "Process components:"
          else
             write (msg_buffer, "(2x,A)")  "Process components: [none]"
          end if
          call msg_message ()
       else
          write (u, "(3x,A)")  "Process components:"
       end if
       do i = 1, size (meta%component_id)
          if (.not. meta%active(i))  cycle
          if (screen) then
             write (msg_buffer, "(4x,I0,9A)")  i, ": '", &
                  char (meta%component_id (i)), "':   ", &
                  char (meta%component_description (i))
             call msg_message ()
          else
             write (u, "(5x,I0,9A)")  i, ": '", &
                  char (meta%component_id (i)), "':   ", &
                  char (meta%component_description (i))
          end if
       end do
    end if
    if (screen) then
       write (msg_buffer, "(A)")  repeat ("-", 72)
       call msg_message ()
    else
       call write_separator (u)
    end if
  end subroutine process_metadata_write

  module subroutine process_metadata_show (meta, u, model_name)
    class(process_metadata_t), intent(in) :: meta
    integer, intent(in) :: u
    type(string_t), intent(in) :: model_name
    integer :: i
    select case (meta%type)
    case (PRC_UNKNOWN)
       write (u, "(A)") "Process: [undefined]"
       return
    case default
       write (u, "(A)", advance="no") "Process:"
    end select
    write (u, "(1x,A)", advance="no") char (meta%id)
    select case (meta%num_id)
    case (0)
    case default
       write (u, "(1x,'(',I0,')')", advance="no") meta%num_id
    end select
    select case (char (model_name))
    case ("")
    case default
       write (u, "(1x,'[',A,']')", advance="no")  char (model_name)
    end select
    write (u, *)
    if (allocated (meta%component_id)) then
       do i = 1, size (meta%component_id)
          if (meta%active(i)) then
             write (u, "(2x,I0,':',1x,A)")  i, &
                  char (meta%component_description (i))
          end if
       end do
    end if
  end subroutine process_metadata_show

  module subroutine process_metadata_init (meta, id, lib, var_list)
    class(process_metadata_t), intent(out) :: meta
    type(string_t), intent(in) :: id
    type(process_library_t), intent(in), target :: lib
    type(var_list_t), intent(in) :: var_list
    select case (lib%get_n_in (id))
    case (1);  meta%type = PRC_DECAY
    case (2);  meta%type = PRC_SCATTERING
    case default
       call msg_bug ("Process '" // char (id) // "': impossible n_in")
    end select
    meta%id = id
    meta%run_id = var_list%get_sval (var_str ("$run_id"))
    allocate (meta%lib_name)
    meta%lib_name = lib%get_name ()
    meta%lib_update_counter = lib%get_update_counter ()
    if (lib%contains (id)) then
       meta%lib_index = lib%get_entry_index (id)
       meta%num_id = lib%get_num_id (id)
       call lib%get_component_list (id, meta%component_id)
       meta%n_components = size (meta%component_id)
       call lib%get_component_description_list &
            (id, meta%component_description)
       allocate (meta%active (meta%n_components), source = .true.)
    else
       call msg_fatal ("Process library does not contain process '" &
            // char (id) // "'")
    end if
    if (.not. lib%is_active ()) then
       call msg_bug ("Process init: inactive library not handled yet")
    end if
  end subroutine process_metadata_init

  module subroutine process_metadata_deactivate_component (meta, i)
    class(process_metadata_t), intent(inout) :: meta
    integer, intent(in) :: i
    call msg_message ("Process component '" &
         // char (meta%component_id(i)) // "': matrix element vanishes")
    meta%active(i) = .false.
  end subroutine process_metadata_deactivate_component

  module subroutine process_phs_config_write (phs_config, unit)
    class(process_phs_config_t), intent(in) :: phs_config
    integer, intent(in), optional :: unit
    integer :: u, iostat
    integer, dimension(:), allocatable :: v_list
    character(0) :: iomsg
    u = given_output_unit (unit)
    allocate (v_list (0))
    call phs_config%write_formatted (u, "LISTDIRECTED", v_list, iostat, iomsg)
  end subroutine process_phs_config_write

  module subroutine process_phs_config_write_formatted &
       (dtv, unit, iotype, v_list, iostat, iomsg)
    class(process_phs_config_t), intent(in) :: dtv
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, dimension(:), intent(in) :: v_list
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    associate (phs_config => dtv)
      write (unit, "(1x, A)")  "Phase-space configuration entry:"
      call phs_config%phs_par%write (unit)
      call phs_config%mapping_defs%write (unit)
    end associate
    iostat = 0
  end subroutine process_phs_config_write_formatted

  module subroutine process_beam_config_write (object, unit, verbose)
    class(process_beam_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u, i, c
    u = given_output_unit (unit)
    call object%data%write (u, verbose = verbose)
    if (object%data%initialized) then
       write (u, "(3x,A,L1)")  "Azimuthal dependence    = ", &
            object%azimuthal_dependence
       write (u, "(3x,A,L1)")  "Lab frame is c.m. frame = ", &
            object%lab_is_cm
       if (object%md5sum /= "") then
          write (u, "(3x,A,A,A)")  "MD5 sum (beams/strf) = '", &
               object%md5sum, "'"
       end if
       if (allocated (object%sf)) then
          do i = 1, size (object%sf)
             call object%sf(i)%write (u)
          end do
          if (any_sf_channel_has_mapping (object%sf_channel)) then
             write (u, "(1x,A,L1)")  "Structure-function mappings per channel:"
             do c = 1, object%n_channel
                write (u, "(3x,I0,':')", advance="no")  c
                call object%sf_channel(c)%write (u)
             end do
          end if
       end if
    end if
  end subroutine process_beam_config_write

  module subroutine process_beam_config_final (object)
    class(process_beam_config_t), intent(inout) :: object
    call object%data%final ()
  end subroutine process_beam_config_final

  module subroutine process_beam_config_init_beam_structure &
       (beam_config, beam_structure, sqrts, model, decay_rest_frame)
    class(process_beam_config_t), intent(out) :: beam_config
    type(beam_structure_t), intent(in) :: beam_structure
    logical, intent(in), optional :: decay_rest_frame
    real(default), intent(in) :: sqrts
    class(model_data_t), intent(in), target :: model
    call beam_config%data%init_structure (beam_structure, &
         sqrts, model, decay_rest_frame)
    beam_config%lab_is_cm = beam_config%data%lab_is_cm
  end subroutine process_beam_config_init_beam_structure

  module subroutine process_beam_config_init_scattering &
       (beam_config, flv_in, sqrts, beam_structure)
    class(process_beam_config_t), intent(out) :: beam_config
    type(flavor_t), dimension(2), intent(in) :: flv_in
    real(default), intent(in) :: sqrts
    type(beam_structure_t), intent(in), optional :: beam_structure
    if (present (beam_structure)) then
       if (beam_structure%polarized ()) then
          call beam_config%data%init_sqrts (sqrts, flv_in, &
               beam_structure%get_smatrix (), beam_structure%get_pol_f ())
       else
          call beam_config%data%init_sqrts (sqrts, flv_in)
       end if
    else
       call beam_config%data%init_sqrts (sqrts, flv_in)
    end if
  end subroutine process_beam_config_init_scattering

  module subroutine process_beam_config_init_decay &
       (beam_config, flv_in, rest_frame, beam_structure)
    class(process_beam_config_t), intent(out) :: beam_config
    type(flavor_t), dimension(1), intent(in) :: flv_in
    logical, intent(in), optional :: rest_frame
    type(beam_structure_t), intent(in), optional :: beam_structure
    if (present (beam_structure)) then
       if (beam_structure%polarized ()) then
          call beam_config%data%init_decay (flv_in, &
               beam_structure%get_smatrix (), beam_structure%get_pol_f (), &
               rest_frame = rest_frame)
       else
          call beam_config%data%init_decay (flv_in, rest_frame = rest_frame)
       end if
    else
       call beam_config%data%init_decay (flv_in, &
            rest_frame = rest_frame)
    end if
    beam_config%lab_is_cm = beam_config%data%lab_is_cm
  end subroutine process_beam_config_init_decay

  module subroutine process_beam_config_startup_message &
       (beam_config, unit, beam_structure)
    class(process_beam_config_t), intent(in) :: beam_config
    integer, intent(in), optional :: unit
    type(beam_structure_t), intent(in), optional :: beam_structure
    integer :: u
    u = free_unit ()
    open (u, status="scratch", action="readwrite")
    if (present (beam_structure)) then
       call beam_structure%write (u)
    end if
    call beam_config%data%write (u)
    rewind (u)
    do
       read (u, "(1x,A)", end=1)  msg_buffer
       call msg_message ()
    end do
1   continue
    close (u)
  end subroutine process_beam_config_startup_message

  module subroutine process_beam_config_init_sf_chain &
       (beam_config, sf_config, sf_trace_file)
    class(process_beam_config_t), intent(inout) :: beam_config
    type(sf_config_t), dimension(:), intent(in) :: sf_config
    type(string_t), intent(in), optional :: sf_trace_file
    integer :: i
    beam_config%n_strfun = size (sf_config)
    allocate (beam_config%sf (beam_config%n_strfun))
    do i = 1, beam_config%n_strfun
       associate (sf => sf_config(i))
         call beam_config%sf(i)%init (sf%i, sf%data)
         if (.not. sf%data%is_generator ()) then
            beam_config%n_sfpar = beam_config%n_sfpar + sf%data%get_n_par ()
         end if
       end associate
    end do
    if (present (sf_trace_file)) then
       beam_config%sf_trace = .true.
       beam_config%sf_trace_file = sf_trace_file
    end if
  end subroutine process_beam_config_init_sf_chain

  module subroutine process_beam_config_allocate_sf_channels &
       (beam_config, n_channel)
    class(process_beam_config_t), intent(inout) :: beam_config
    integer, intent(in) :: n_channel
    beam_config%n_channel = n_channel
    call allocate_sf_channels (beam_config%sf_channel, &
         n_channel = n_channel, &
         n_strfun = beam_config%n_strfun)
  end subroutine process_beam_config_allocate_sf_channels

  module subroutine process_beam_config_set_sf_channel &
       (beam_config, c, sf_channel)
    class(process_beam_config_t), intent(inout) :: beam_config
    integer, intent(in) :: c
    type(sf_channel_t), intent(in) :: sf_channel
    beam_config%sf_channel(c) = sf_channel
  end subroutine process_beam_config_set_sf_channel

  module subroutine process_beam_config_sf_startup_message &
       (beam_config, sf_string, unit)
    class(process_beam_config_t), intent(in) :: beam_config
    type(string_t), intent(in) :: sf_string
    integer, intent(in), optional :: unit
    if (beam_config%n_strfun > 0) then
       call msg_message ("Beam structure: " // char (sf_string), unit = unit)
       write (msg_buffer, "(A,3(1x,I0,1x,A))") &
            "Beam structure:", &
            beam_config%n_channel, "channels,", &
            beam_config%n_sfpar, "dimensions"
       call msg_message (unit = unit)
       if (beam_config%sf_trace) then
          call msg_message ("Beam structure: tracing &
               &values in '" // char (beam_config%sf_trace_file) // "'")
       end if
    end if
  end subroutine process_beam_config_sf_startup_message

  module function process_beam_config_get_pdf_set (beam_config) result (pdf_set)
    class(process_beam_config_t), intent(in) :: beam_config
    integer :: pdf_set
    integer :: i
    pdf_set = 0
    if (allocated (beam_config%sf)) then
       do i = 1, size (beam_config%sf)
          pdf_set = beam_config%sf(i)%get_pdf_set ()
          if (pdf_set /= 0)  return
       end do
    end if
  end function process_beam_config_get_pdf_set

  module function process_beam_config_get_beam_file (beam_config) result (file)
    class(process_beam_config_t), intent(in) :: beam_config
    type(string_t) :: file
    integer :: i
    file = ""
    if (allocated (beam_config%sf)) then
       do i = 1, size (beam_config%sf)
          file = beam_config%sf(i)%get_beam_file ()
          if (file /= "")  return
       end do
    end if
  end function process_beam_config_get_beam_file

  module subroutine process_beam_config_compute_md5sum (beam_config)
    class(process_beam_config_t), intent(inout) :: beam_config
    integer :: u
    if (beam_config%md5sum == "") then
       u = free_unit ()
       open (u, status = "scratch", action = "readwrite")
       call beam_config%write (u, verbose=.true.)
       rewind (u)
       beam_config%md5sum = md5sum (u)
       close (u)
    end if
  end subroutine process_beam_config_compute_md5sum

  pure module function process_beam_config_get_md5sum (beam_config) result (md5)
    character(32) :: md5
    class(process_beam_config_t), intent(in) :: beam_config
    md5 = beam_config%md5sum
  end function process_beam_config_get_md5sum

  pure module function process_beam_config_has_structure_function &
       (beam_config) result (has_sf)
    logical :: has_sf
    class(process_beam_config_t), intent(in) :: beam_config
    has_sf = beam_config%n_strfun > 0
  end function process_beam_config_has_structure_function

  module subroutine process_component_final (object)
    class(process_component_t), intent(inout) :: object
    if (allocated (object%phs_config)) then
       call object%phs_config%final ()
    end if
  end subroutine process_component_final

  module subroutine process_component_write (object, unit)
    class(process_component_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%config)) then
       write (u, "(1x,A,I0)")  "Component #", object%index
       call object%config%write (u)
       if (object%md5sum_phs /= "") then
          write (u, "(3x,A,A,A)")  "MD5 sum (phs)       = '", &
               object%md5sum_phs, "'"
       end if
    else
       write (u, "(1x,A)") "Process component: [not allocated]"
    end if
    if (.not. object%active) then
       write (u, "(1x,A)") "[Inactive]"
       return
    end if
    write (u, "(1x,A)") "Referenced data:"
    if (allocated (object%i_term)) then
       write (u, "(3x,A,999(1x,I0))") "Terms                    =", &
            object%i_term
    else
       write (u, "(3x,A)") "Terms                    = [undefined]"
    end if
    if (object%i_mci /= 0) then
       write (u, "(3x,A,I0)") "MC dataset               = ", object%i_mci
    else
       write (u, "(3x,A)") "MC dataset               = [undefined]"
    end if
    if (allocated (object%phs_config)) then
       call object%phs_config%write (u)
    end if
  end subroutine process_component_write

  module subroutine process_component_init (component, &
       i_component, env, meta, config, &
       active, &
       phs_config_template)
    class(process_component_t), intent(out) :: component
    integer, intent(in) :: i_component
    type(process_environment_t), intent(in) :: env
    type(process_metadata_t), intent(in) :: meta
    type(process_config_data_t), intent(in) :: config
    logical, intent(in) :: active
    class(phs_config_t), intent(in), allocatable :: phs_config_template

    type(process_constants_t) :: data

    component%index = i_component
    component%config => &
         config%process_def%get_component_def_ptr (i_component)

    component%active = active
    if (component%active) then
       allocate (component%phs_config, source = phs_config_template)
       call env%fill_process_constants (meta%id, i_component, data)
       call component%phs_config%init (data, config%model)
    end if
  end subroutine process_component_init

  elemental module function process_component_is_active &
       (component) result (active)
    logical :: active
    class(process_component_t), intent(in) :: component
    active = component%active
  end function process_component_is_active

  module subroutine process_component_configure_phs &
       (component, sqrts, beam_config, rebuild, &
        ignore_mismatch, subdir)
    class(process_component_t), intent(inout) :: component
    real(default), intent(in) :: sqrts
    type(process_beam_config_t), intent(in) :: beam_config
    logical, intent(in), optional :: rebuild
    logical, intent(in), optional :: ignore_mismatch
    type(string_t), intent(in), optional :: subdir
    logical :: no_strfun
    integer :: nlo_type
    no_strfun = beam_config%n_strfun == 0
    nlo_type = component%config%get_nlo_type ()
    call component%phs_config%configure (sqrts, &
         azimuthal_dependence = beam_config%azimuthal_dependence, &
         sqrts_fixed = no_strfun, &
         lab_is_cm = beam_config%lab_is_cm .and. no_strfun, &
         rebuild = rebuild, ignore_mismatch = ignore_mismatch, &
         nlo_type = nlo_type, &
         subdir = subdir)
  end subroutine process_component_configure_phs

  module subroutine process_component_compute_md5sum (component)
    class(process_component_t), intent(inout) :: component
    component%md5sum_phs = component%phs_config%get_md5sum ()
  end subroutine process_component_compute_md5sum

  module subroutine process_component_collect_channels (component, coll)
    class(process_component_t), intent(inout) :: component
    type(phs_channel_collection_t), intent(inout) :: coll
    call component%phs_config%collect_channels (coll)
  end subroutine process_component_collect_channels

  module function process_component_get_config (component) &
       result (config)
    type(process_component_def_t) :: config
    class(process_component_t), intent(in) :: component
    config = component%config
  end function process_component_get_config

  pure module function process_component_get_md5sum (component) result (md5)
    type(string_t) :: md5
    class(process_component_t), intent(in) :: component
    md5 = component%config%get_md5sum () // component%md5sum_phs
  end function process_component_get_md5sum

  module function process_component_get_n_phs_par (component) result (n_par)
    class(process_component_t), intent(in) :: component
    integer :: n_par
    n_par = component%phs_config%get_n_par ()
  end function process_component_get_n_phs_par

  module subroutine process_component_get_phs_config (component, phs_config)
    class(process_component_t), intent(in), target :: component
    class(phs_config_t), intent(out), pointer :: phs_config
    phs_config => component%phs_config
  end subroutine process_component_get_phs_config

  elemental module function process_component_get_nlo_type &
       (component) result (nlo_type)
    integer :: nlo_type
    class(process_component_t), intent(in) :: component
    nlo_type = component%config%get_nlo_type ()
  end function process_component_get_nlo_type

  module function process_component_needs_mci_entry &
       (component, combined_integration) result (value)
    logical :: value
    class(process_component_t), intent(in) :: component
    logical, intent(in), optional :: combined_integration
    value = component%active
    if (present (combined_integration)) then
       if (combined_integration) &
            value = value .and. component%component_type <= COMP_MASTER
    end if
  end function process_component_needs_mci_entry

  elemental module function process_component_can_be_integrated &
       (component) result (active)
    logical :: active
    class(process_component_t), intent(in) :: component
    active = component%config%can_be_integrated ()
  end function process_component_can_be_integrated

  module subroutine process_term_write (term, unit)
    class(process_term_t), intent(in) :: term
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A,I0)")  "Term #", term%i_term_global
    write (u, "(3x,A,I0)")  "Process component index      = ", &
         term%i_component
    write (u, "(3x,A,I0)")  "Term index w.r.t. component  = ", &
         term%i_term
    call write_separator (u)
    write (u, "(1x,A)")  "Hard interaction:"
    call write_separator (u)
    call term%int%basic_write (u)
  end subroutine process_term_write

  module subroutine process_term_write_state_summary (term, core, unit)
    class(process_term_t), intent(in) :: term
    class(prc_core_t), intent(in) :: core
    integer, intent(in), optional :: unit
    integer :: u, i, f, h, c
    type(state_iterator_t) :: it
    character :: sgn
    u = given_output_unit (unit)
    write (u, "(1x,A,I0)")  "Term #", term%i_term_global
    call it%init (term%int%get_state_matrix_ptr ())
    do while (it%is_valid ())
       i = it%get_me_index ()
       f = term%flv(i)
       h = term%hel(i)
       if (allocated (term%col)) then
          c = term%col(i)
       else
          c = 1
       end if
       if (core%is_allowed (term%i_term, f, h, c)) then
          sgn = "+"
       else
          sgn = " "
       end if
       write (u, "(1x,A1,1x,I0,2x)", advance="no")  sgn, i
       call quantum_numbers_write (it%get_quantum_numbers (), u)
       write (u, *)
       call it%advance ()
    end do
  end subroutine process_term_write_state_summary

  module subroutine process_term_final (term)
    class(process_term_t), intent(inout) :: term
    call term%int%final ()
  end subroutine process_term_final

  module subroutine process_term_init &
       (term, i_term_global, i_component, i_term, core, model, &
        nlo_type, use_beam_pol, subtraction_method, &
        has_pdfs, n_emitters)
    class(process_term_t), intent(inout), target :: term
    integer, intent(in) :: i_term_global
    integer, intent(in) :: i_component
    integer, intent(in) :: i_term
    class(prc_core_t), intent(inout) :: core
    class(model_data_t), intent(in), target :: model
    integer, intent(in), optional :: nlo_type
    logical, intent(in), optional :: use_beam_pol
    type(string_t), intent(in), optional :: subtraction_method
    logical, intent(in), optional :: has_pdfs
    integer, intent(in), optional :: n_emitters
    class(modelpar_data_t), pointer :: alpha_s_ptr
    logical :: use_internal_color
    term%i_term_global = i_term_global
    term%i_component = i_component
    term%i_term = i_term
    call core%get_constants (term%data, i_term)
    alpha_s_ptr => model%get_par_data_ptr (var_str ("alphas"))
    if (associated (alpha_s_ptr)) then
       term%alpha_s = alpha_s_ptr%get_real ()
    else
       term%alpha_s = -1
    end if
    use_internal_color = .false.
    if (present (subtraction_method)) &
         use_internal_color = (char (subtraction_method) == 'omega') &
         .or. (char (subtraction_method) == 'threshold')
    call term%setup_interaction (core, model, nlo_type = nlo_type, &
         pol_beams = use_beam_pol, use_internal_color = use_internal_color, &
         has_pdfs = has_pdfs, n_emitters = n_emitters)
  end subroutine process_term_init

  module subroutine process_term_setup_interaction (term, core, model, &
     nlo_type, pol_beams, has_pdfs, use_internal_color, n_emitters)
    class(process_term_t), intent(inout) :: term
    class(prc_core_t), intent(inout) :: core
    class(model_data_t), intent(in), target :: model
    logical, intent(in), optional :: pol_beams
    logical, intent(in), optional :: has_pdfs
    integer, intent(in), optional :: nlo_type
    logical, intent(in), optional :: use_internal_color
    integer, intent(in), optional :: n_emitters
    integer :: n, n_tot
    type(flavor_t), dimension(:), allocatable :: flv
    type(color_t), dimension(:), allocatable :: col
    type(helicity_t), dimension(:), allocatable :: hel
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    logical :: is_pol, use_color
    integer :: nlo_t, n_sub
    is_pol = .false.; if (present (pol_beams)) is_pol = pol_beams
    nlo_t = BORN; if (present (nlo_type)) nlo_t = nlo_type
    n_tot = term%data%n_in + term%data%n_out
    call count_number_of_states ()
    term%n_allowed = n
    call compute_n_sub (n_emitters, has_pdfs)
    call fill_quantum_numbers ()
    call term%int%basic_init &
         (term%data%n_in, 0, term%data%n_out, set_relations = .true.)
    select type (core)
    class is (prc_blha_t)
       call setup_states_blha_olp ()
    type is (prc_threshold_t)
       call setup_states_threshold ()
    class is (prc_external_t)
       call setup_states_other_prc_external ()
    class default
       call setup_states_omega ()
    end select
    call term%int%freeze ()
  contains
    subroutine count_number_of_states ()
      integer :: f, h, c
      n = 0
      select type (core)
      class is (prc_external_t)
         do f = 1, term%data%n_flv
            do h = 1, term%data%n_hel
               do c = 1, term%data%n_col
                  n = n + 1
               end do
            end do
         end do
      class default !!! Omega and all test cores
         do f = 1, term%data%n_flv
            do h = 1, term%data%n_hel
               do c = 1, term%data%n_col
                  if (core%is_allowed (term%i_term, f, h, c))  n = n + 1
               end do
            end do
         end do
      end select
    end subroutine count_number_of_states

    subroutine compute_n_sub (n_emitters, has_pdfs)
      integer, intent(in), optional :: n_emitters
      logical, intent(in), optional :: has_pdfs
      logical :: can_have_sub
      integer :: n_sub_color, n_sub_spin
      use_color = .false.; if (present (use_internal_color)) &
           use_color = use_internal_color
      can_have_sub = nlo_t == NLO_VIRTUAL .or. &
           (nlo_t == NLO_REAL .and. term%i_term_global == term%i_sub) .or. &
           nlo_t == NLO_MISMATCH .or. nlo_t == NLO_DGLAP
      n_sub_color = 0; n_sub_spin = 0
      if (can_have_sub) then
         if (.not. use_color) n_sub_color = n_tot * (n_tot - 1) / 2
         if (nlo_t == NLO_REAL) then
            if (present (n_emitters)) then
               n_sub_spin = 6 * n_emitters
            end if
         end if
      end if
      n_sub = n_sub_color + n_sub_spin
      !!! For the virtual subtraction we also need the finite virtual contribution
      !!! corresponding to the $\epsilon^0$-pole
      if (nlo_t == NLO_VIRTUAL)  n_sub = n_sub + 1
      if (present (has_pdfs)) then
         if (has_pdfs &
              .and. ((nlo_t == NLO_REAL .and. can_have_sub) &
              .or. nlo_t == NLO_DGLAP)) then
            !!! necessary dummy, needs refactoring,
            !!! c.f. [[term_instance_evaluate_interaction_external_tree]]
            n_sub = n_sub + n_beams_rescaled
         end if
      end if
      term%n_sub = n_sub
      term%n_sub_color = n_sub_color
      term%n_sub_spin = n_sub_spin
    end subroutine compute_n_sub

    subroutine fill_quantum_numbers ()
      integer :: nn
      logical :: can_have_sub
      select type (core)
      class is (prc_external_t)
         can_have_sub = nlo_t == NLO_VIRTUAL .or. &
              (nlo_t == NLO_REAL .and. term%i_term_global == term%i_sub) .or. &
              nlo_t == NLO_MISMATCH .or. nlo_t == NLO_DGLAP
         if (can_have_sub) then
            nn = (n_sub + 1) * n
         else
            nn = n
         end if
      class default
         nn = n
      end select
      allocate (term%flv (nn), term%col (nn), term%hel (nn))
      allocate (flv (n_tot), col (n_tot), hel (n_tot))
      allocate (qn (n_tot))
    end subroutine fill_quantum_numbers

    subroutine setup_states_blha_olp ()
      integer :: s, f, c, h, i
      i = 0
      associate (data => term%data)
         do s = 0, n_sub
             do f = 1, data%n_flv
                do h = 1, data%n_hel
                   do c = 1, data%n_col
                      i = i + 1
                      term%flv(i) = f
                      term%hel(i) = h
                      !!! Dummy-initialization of color
                      term%col(i) = c
                      call flv%init (data%flv_state (:,f), model)
                      call color_init_from_array (col, &
                           data%col_state(:,:,c), data%ghost_flag(:,c))
                      call col(1:data%n_in)%invert ()
                      if (is_pol) then
                         select type (core)
                         type is (prc_openloops_t)
                            call hel%init (data%hel_state (:,h))
                            call qn%init (flv, hel, col, s)
                         class default
                            call msg_fatal ("Polarized beams only supported by OpenLoops")
                         end select
                      else
                         call qn%init (flv, col, s)
                      end if
                      call qn%tag_hard_process ()
                      call term%int%add_state (qn)
                  end do
               end do
             end do
         end do
      end associate
    end subroutine setup_states_blha_olp

    subroutine setup_states_threshold ()
      integer :: s, f, c, h, i
      i = 0
      n_sub = 0; if (nlo_t == NLO_VIRTUAL) n_sub = 1
      associate (data => term%data)
         do s = 0, n_sub
            do f = 1, term%data%n_flv
               do h = 1, data%n_hel
                  do c = 1, data%n_col
                     i = i + 1
                     term%flv(i) = f
                     term%hel(i) = h
                     !!! Dummy-initialization of color
                     term%col(i) = 1
                     call flv%init (term%data%flv_state (:,f), model)
                     if (is_pol) then
                        call hel%init (data%hel_state (:,h))
                        call qn%init (flv, hel, s)
                     else
                        call qn%init (flv, s)
                     end if
                     call qn%tag_hard_process ()
                     call term%int%add_state (qn)
                  end do
               end do
            end do
         end do
      end associate
    end subroutine setup_states_threshold

    subroutine setup_states_other_prc_external ()
      integer :: s, f, i, c, h
      if (is_pol) &
         call msg_fatal ("Polarized beams only supported by OpenLoops")
      i = 0
      !!! n_sub = 0; if (nlo_t == NLO_VIRTUAL) n_sub = 1
      associate (data => term%data)
        do s = 0, n_sub
           do f = 1, data%n_flv
              do h = 1, data%n_hel
                 do c = 1, data%n_col
                    i = i + 1
                    term%flv(i) = f
                    term%hel(i) = h
                    !!! Dummy-initialization of color
                    term%col(i) = c
                    call flv%init (data%flv_state (:,f), model)
                    call color_init_from_array (col, &
                         data%col_state(:,:,c), data%ghost_flag(:,c))
                    call col(1:data%n_in)%invert ()
                    call qn%init (flv, col, s)
                    call qn%tag_hard_process ()
                    call term%int%add_state (qn)
                 end do
              end do
           end do
        end do
      end associate
    end subroutine setup_states_other_prc_external

    subroutine setup_states_omega ()
      integer :: f, h, c, i
      i = 0
      associate (data => term%data)
         do f = 1, data%n_flv
            do h = 1, data%n_hel
              do c = 1, data%n_col
                 if (core%is_allowed (term%i_term, f, h, c)) then
                    i = i + 1
                    term%flv(i) = f
                    term%hel(i) = h
                    term%col(i) = c
                    call flv%init (data%flv_state(:,f), model)
                    call color_init_from_array (col, &
                         data%col_state(:,:,c), &
                         data%ghost_flag(:,c))
                    call col(:data%n_in)%invert ()
                    call hel%init (data%hel_state(:,h))
                    call qn%init (flv, col, hel)
                    call qn%tag_hard_process ()
                    call term%int%add_state (qn)
                 end if
              end do
            end do
         end do
      end associate
    end subroutine setup_states_omega

  end subroutine process_term_setup_interaction

  module subroutine process_term_get_process_constants &
       (term, prc_constants)
    class(process_term_t), intent(inout) :: term
    type(process_constants_t), intent(out) :: prc_constants
    prc_constants = term%data
  end subroutine process_term_get_process_constants


end submodule process_config_s

