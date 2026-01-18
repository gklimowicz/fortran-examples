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

submodule (whizard) whizard_s

  use system_defs, only: VERSION_STRING
  use system_defs, only: EOF, BACKSLASH
  use diagnostics
  use eval_trees
  use models
  use phs_forests
  use prclib_stacks
  use slha_interface
  use commands

  implicit none

contains

  module subroutine pt_stack_final (pt_stack)
    class(pt_stack_t), intent(inout) :: pt_stack
    type(pt_entry_t), pointer :: current
    do while (associated (pt_stack%last))
       current => pt_stack%last
       pt_stack%last => current%previous
       call parse_tree_final (current%parse_tree_t)
       deallocate (current)
    end do
  end subroutine pt_stack_final

  module subroutine pt_stack_push (pt_stack, parse_tree)
    class(pt_stack_t), intent(inout) :: pt_stack
    type(parse_tree_t), intent(out), pointer :: parse_tree
    type(pt_entry_t), pointer :: current
    allocate (current)
    parse_tree => current%parse_tree_t
    current%previous => pt_stack%last
    pt_stack%last => current
  end subroutine pt_stack_push

  module subroutine whizard_init (whizard, options, paths, logfile)
    class(whizard_t), intent(out), target :: whizard
    type(whizard_options_t), intent(in) :: options
    type(paths_t), intent(in), optional :: paths
    type(string_t), intent(in), optional :: logfile
    call init_syntax_tables ()
    whizard%options = options
    call whizard%global%global_init (paths, logfile)
    call whizard%init_job_id ()
    call whizard%init_rebuild_flags ()
    call whizard%unpack_files ()
    call whizard%preload_model ()
    call whizard%preload_library ()
    call whizard%global%init_fallback_model &
         (var_str ("SM_hadrons"), var_str ("SM_hadrons.mdl"))
  end subroutine whizard_init

  module subroutine whizard_final (whizard)
    class(whizard_t), intent(inout), target :: whizard
    call whizard%global%final ()
    call whizard%pt_stack%final ()
    call whizard%pack_files ()
    call final_syntax_tables ()
  end subroutine whizard_final

  module subroutine whizard_init_job_id (whizard)
    class(whizard_t), intent(inout), target :: whizard
    associate (var_list => whizard%global%var_list, options => whizard%options)
      if (options%job_id /= "") then
         call var_list%set_string (var_str ("$job_id"), &
              options%job_id, is_known=.true.)
      end if
    end associate
  end subroutine whizard_init_job_id

  module subroutine whizard_init_rebuild_flags (whizard)
    class(whizard_t), intent(inout), target :: whizard
    associate (var_list => whizard%global%var_list, options => whizard%options)
      call var_list%append_log (var_str ("?rebuild_library"), &
           options%rebuild_library, intrinsic=.true.)
      call var_list%append_log (var_str ("?recompile_library"), &
           options%recompile_library, intrinsic=.true.)
      call var_list%append_log (var_str ("?rebuild_phase_space"), &
           options%rebuild_phs, intrinsic=.true.)
      call var_list%append_log (var_str ("?rebuild_grids"), &
           options%rebuild_grids, intrinsic=.true.)
      call var_list%append_log (var_str ("?rebuild_events"), &
           options%rebuild_events, intrinsic=.true.)
    end associate
  end subroutine whizard_init_rebuild_flags

  module subroutine whizard_pack_files (whizard)
    class(whizard_t), intent(in), target :: whizard
    logical :: exist
    integer :: i
    type(string_t) :: file
    if (allocated (whizard%options%pack_args)) then
       do i = 1, size (whizard%options%pack_args)
          file = whizard%options%pack_args(i)
          call msg_message ("Packing file/dir '" // char (file) // "'")
          exist = os_file_exist (file) .or. os_dir_exist (file)
          if (exist) then
             call os_pack_file (whizard%options%pack_args(i), &
                  whizard%global%os_data)
          else
             call msg_error ("File/dir '" // char (file) // "' not found")
          end if
       end do
    end if
  end subroutine whizard_pack_files

  module subroutine whizard_unpack_files (whizard)
    class(whizard_t), intent(in), target :: whizard
    logical :: exist
    integer :: i
    type(string_t) :: file
    if (allocated (whizard%options%unpack_args)) then
       do i = 1, size (whizard%options%unpack_args)
          file = whizard%options%unpack_args(i)
          call msg_message ("Unpacking file '" // char (file) // "'")
          exist = os_file_exist (file)
          if (exist) then
             call os_unpack_file (whizard%options%unpack_args(i), &
                  whizard%global%os_data)
          else
             call msg_error ("File '" // char (file) // "' not found")
          end if
       end do
    end if
  end subroutine whizard_unpack_files

  module subroutine whizard_preload_model (whizard)
    class(whizard_t), intent(inout), target :: whizard
    type(string_t) :: model_name
    model_name = whizard%options%preload_model
    if (model_name /= "") then
       call whizard%global%read_model (model_name, whizard%global%preload_model)
       whizard%global%model => whizard%global%preload_model
       if (associated (whizard%global%model)) then
          call whizard%global%model%link_var_list (whizard%global%var_list)
          call whizard%global%var_list%set_string (var_str ("$model_name"), &
               model_name, is_known = .true.)
          call msg_message ("Preloaded model: " &
               // char (model_name))
       else
          call msg_fatal ("Preloading model " // char (model_name) &
               // " failed")
       end if
    else
       call msg_message ("No model preloaded")
    end if
  end subroutine whizard_preload_model

  module subroutine whizard_preload_library (whizard)
    class(whizard_t), intent(inout), target :: whizard
    type(string_t) :: library_name, libs
    type(string_t), dimension(:), allocatable :: libname_static
    type(prclib_entry_t), pointer :: lib_entry
    integer :: i
    call get_prclib_static (libname_static)
    do i = 1, size (libname_static)
       allocate (lib_entry)
       call lib_entry%init_static (libname_static(i))
       call whizard%global%add_prclib (lib_entry)
    end do
    libs = adjustl (whizard%options%preload_libraries)
    if (libs == "" .and. whizard%options%default_lib /= "") then
       allocate (lib_entry)
       call lib_entry%init (whizard%options%default_lib)
       call whizard%global%add_prclib (lib_entry)
       call msg_message ("Preloaded library: " // &
            char (whizard%options%default_lib))
    end if
    SCAN_LIBS: do while (libs /= "")
       call split (libs, library_name, " ")
       if (library_name /= "") then
          allocate (lib_entry)
          call lib_entry%init (library_name)
          call whizard%global%add_prclib (lib_entry)
          call msg_message ("Preloaded library: " // char (library_name))
       end if
    end do SCAN_LIBS
  end subroutine whizard_preload_library

  module subroutine init_syntax_tables ()
    if (.not. syntax_tables_exist) then
       call syntax_model_file_init ()
       call syntax_phs_forest_init ()
       call syntax_pexpr_init ()
       call syntax_slha_init ()
       call syntax_cmd_list_init ()
       syntax_tables_exist = .true.
    end if
  end subroutine init_syntax_tables

  module subroutine final_syntax_tables ()
    if (syntax_tables_exist) then
       call syntax_model_file_final ()
       call syntax_phs_forest_final ()
       call syntax_pexpr_final ()
       call syntax_slha_final ()
       call syntax_cmd_list_final ()
       syntax_tables_exist = .false.
    end if
  end subroutine final_syntax_tables

  module subroutine write_syntax_tables ()
    integer :: unit
    character(*), parameter :: file_model = "whizard.model_file.syntax"
    character(*), parameter :: file_phs = "whizard.phase_space_file.syntax"
    character(*), parameter :: file_pexpr = "whizard.prt_expressions.syntax"
    character(*), parameter :: file_slha = "whizard.slha.syntax"
    character(*), parameter :: file_sindarin = "whizard.sindarin.syntax"
    if (.not. syntax_tables_exist)  call init_syntax_tables ()
    unit = free_unit ()
    print *, "Writing file '" // file_model // "'"
    open (unit=unit, file=file_model, status="replace", action="write")
    write (unit, "(A)")  VERSION_STRING
    write (unit, "(A)")  "Syntax definition file: " // file_model
    call syntax_model_file_write (unit)
    close (unit)
    print *, "Writing file '" // file_phs // "'"
    open (unit=unit, file=file_phs, status="replace", action="write")
    write (unit, "(A)")  VERSION_STRING
    write (unit, "(A)")  "Syntax definition file: " // file_phs
    call syntax_phs_forest_write (unit)
    close (unit)
    print *, "Writing file '" // file_pexpr // "'"
    open (unit=unit, file=file_pexpr, status="replace", action="write")
    write (unit, "(A)")  VERSION_STRING
    write (unit, "(A)")  "Syntax definition file: " // file_pexpr
    call syntax_pexpr_write (unit)
    close (unit)
    print *, "Writing file '" // file_slha // "'"
    open (unit=unit, file=file_slha, status="replace", action="write")
    write (unit, "(A)")  VERSION_STRING
    write (unit, "(A)")  "Syntax definition file: " // file_slha
    call syntax_slha_write (unit)
    close (unit)
    print *, "Writing file '" // file_sindarin // "'"
    open (unit=unit, file=file_sindarin, status="replace", action="write")
    write (unit, "(A)")  VERSION_STRING
    write (unit, "(A)")  "Syntax definition file: " // file_sindarin
    call syntax_cmd_list_write (unit)
    close (unit)
  end subroutine write_syntax_tables

  module subroutine whizard_process_ifile (whizard, ifile, quit, quit_code)
    class(whizard_t), intent(inout), target :: whizard
    type(ifile_t), intent(in) :: ifile
    logical, intent(out) :: quit
    integer, intent(out) :: quit_code
    type(lexer_t), target :: lexer
    type(stream_t), target :: stream
    call msg_message ("Reading commands given on the command line")
    call lexer_init_cmd_list (lexer)
    call stream_init (stream, ifile)
    call whizard%process_stream (stream, lexer, quit, quit_code)
    call stream_final (stream)
    call lexer_final (lexer)
  end subroutine whizard_process_ifile

  module subroutine whizard_process_stdin (whizard, quit, quit_code)
    class(whizard_t), intent(inout), target :: whizard
    logical, intent(out) :: quit
    integer, intent(out) :: quit_code
    type(lexer_t), target :: lexer
    type(stream_t), target :: stream
    call msg_message ("Reading commands from standard input")
    call lexer_init_cmd_list (lexer)
    call stream_init (stream, 5)
    call whizard%process_stream (stream, lexer, quit, quit_code)
    call stream_final (stream)
    call lexer_final (lexer)
  end subroutine whizard_process_stdin

  module subroutine whizard_process_file (whizard, file, quit, quit_code)
    class(whizard_t), intent(inout), target :: whizard
    type(string_t), intent(in) :: file
    logical, intent(out) :: quit
    integer, intent(out) :: quit_code
    type(lexer_t), target :: lexer
    type(stream_t), target :: stream
    logical :: exist
    call msg_message ("Reading commands from file '" // char (file) // "'")
    inquire (file=char(file), exist=exist)
    if (exist) then
       call lexer_init_cmd_list (lexer)
       call stream_init (stream, char (file))
       call whizard%process_stream (stream, lexer, quit, quit_code)
       call stream_final (stream)
       call lexer_final (lexer)
    else
       call msg_error ("File '" // char (file) // "' not found")
    end if
  end subroutine whizard_process_file

  module subroutine whizard_process_stream &
       (whizard, stream, lexer, quit, quit_code)
    class(whizard_t), intent(inout), target :: whizard
    type(stream_t), intent(inout), target :: stream
    type(lexer_t), intent(inout), target :: lexer
    logical, intent(out) :: quit
    integer, intent(out) :: quit_code
    type(parse_tree_t), pointer :: parse_tree
    type(command_list_t), target :: command_list
    call lexer_assign_stream (lexer, stream)
    call whizard%pt_stack%push (parse_tree)
    call parse_tree_init (parse_tree, syntax_cmd_list, lexer)
    if (associated (parse_tree%get_root_ptr ())) then
       whizard%global%lexer => lexer
       call command_list%compile (parse_tree%get_root_ptr (), &
            whizard%global)
    end if
    call whizard%global%activate ()
    call command_list%execute (whizard%global)
    call command_list%final ()
    quit = whizard%global%quit
    quit_code = whizard%global%quit_code
  end subroutine whizard_process_stream

  module subroutine whizard_shell (whizard, quit_code)
    class(whizard_t), intent(inout), target :: whizard
    integer, intent(out) :: quit_code
    type(lexer_t), target :: lexer
    type(stream_t), target :: stream
    type(string_t) :: prompt1
    type(string_t) :: prompt2
    type(string_t) :: input
    type(string_t) :: extra
    integer :: last
    integer :: iostat
    logical :: mask_tmp
    logical :: quit
    call msg_message ("Launching interactive shell")
    call lexer_init_cmd_list (lexer)
    prompt1 = "whish? "
    prompt2 = "     > "
    COMMAND_LOOP: do
       call put (6, prompt1)
       call get (5, input, iostat=iostat)
       if (iostat > 0 .or. iostat == EOF) exit COMMAND_LOOP
       CONTINUE_INPUT: do
          last = len_trim (input)
          if (extract (input, last, last) /= BACKSLASH)  exit CONTINUE_INPUT
          call put (6, prompt2)
          call get (5, extra, iostat=iostat)
          if (iostat > 0) exit COMMAND_LOOP
          input = replace (input, last, extra)
       end do CONTINUE_INPUT
       call stream_init (stream, input)
       mask_tmp = mask_fatal_errors
       mask_fatal_errors = .true.
       call whizard%process_stream (stream, lexer, quit, quit_code)
       msg_count = 0
       mask_fatal_errors = mask_tmp
       call stream_final (stream)
       if (quit)  exit COMMAND_LOOP
    end do COMMAND_LOOP
    print *
    call lexer_final (lexer)
  end subroutine whizard_shell


end submodule whizard_s

