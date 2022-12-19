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

module whizard

  use io_units
  use iso_varying_string, string_t => varying_string
  use os_interface
  use ifiles
  use lexers
  use parser
  use rt_data

  implicit none
  private

  public :: whizard_options_t
  public :: whizard_t
  public :: init_syntax_tables
  public :: final_syntax_tables
  public :: write_syntax_tables

  type :: whizard_options_t
     type(string_t) :: job_id
     type(string_t), dimension(:), allocatable :: pack_args
     type(string_t), dimension(:), allocatable :: unpack_args
     type(string_t) :: preload_model
     type(string_t) :: default_lib
     type(string_t) :: preload_libraries
     logical :: rebuild_library = .false.
     logical :: recompile_library = .false.
     logical :: rebuild_phs = .false.
     logical :: rebuild_grids = .false.
     logical :: rebuild_events = .false.
  end type whizard_options_t

  type, extends (parse_tree_t) :: pt_entry_t
     type(pt_entry_t), pointer :: previous => null ()
  end type pt_entry_t

  type :: pt_stack_t
     type(pt_entry_t), pointer :: last => null ()
   contains
     procedure :: final => pt_stack_final
     procedure :: push => pt_stack_push
  end type pt_stack_t

  type :: whizard_t
     type(whizard_options_t) :: options
     type(rt_data_t) :: global
     type(pt_stack_t) :: pt_stack
   contains
     procedure :: init => whizard_init
     procedure :: final => whizard_final
     procedure :: init_job_id => whizard_init_job_id
     procedure :: init_rebuild_flags => whizard_init_rebuild_flags
     procedure :: pack_files => whizard_pack_files
     procedure :: unpack_files => whizard_unpack_files
     procedure :: preload_model => whizard_preload_model
     procedure :: preload_library => whizard_preload_library
     procedure :: process_ifile => whizard_process_ifile
     procedure :: process_stdin => whizard_process_stdin
     procedure :: process_file => whizard_process_file
     procedure :: process_stream => whizard_process_stream
     procedure :: shell => whizard_shell
  end type whizard_t


  logical :: syntax_tables_exist = .false.

  save

  interface
    module subroutine pt_stack_final (pt_stack)
      class(pt_stack_t), intent(inout) :: pt_stack
    end subroutine pt_stack_final
    module subroutine pt_stack_push (pt_stack, parse_tree)
      class(pt_stack_t), intent(inout) :: pt_stack
      type(parse_tree_t), intent(out), pointer :: parse_tree
    end subroutine pt_stack_push
    module subroutine whizard_init (whizard, options, paths, logfile)
      class(whizard_t), intent(out), target :: whizard
      type(whizard_options_t), intent(in) :: options
      type(paths_t), intent(in), optional :: paths
      type(string_t), intent(in), optional :: logfile
    end subroutine whizard_init
    module subroutine whizard_final (whizard)
      class(whizard_t), intent(inout), target :: whizard
    end subroutine whizard_final
    module subroutine whizard_init_job_id (whizard)
      class(whizard_t), intent(inout), target :: whizard
    end subroutine whizard_init_job_id
    module subroutine whizard_init_rebuild_flags (whizard)
      class(whizard_t), intent(inout), target :: whizard
    end subroutine whizard_init_rebuild_flags
    module subroutine whizard_pack_files (whizard)
      class(whizard_t), intent(in), target :: whizard
    end subroutine whizard_pack_files
    module subroutine whizard_unpack_files (whizard)
      class(whizard_t), intent(in), target :: whizard
    end subroutine whizard_unpack_files
    module subroutine whizard_preload_model (whizard)
      class(whizard_t), intent(inout), target :: whizard
    end subroutine whizard_preload_model
    module subroutine whizard_preload_library (whizard)
      class(whizard_t), intent(inout), target :: whizard
    end subroutine whizard_preload_library
    module subroutine init_syntax_tables ()
    end subroutine init_syntax_tables
    module subroutine final_syntax_tables ()
    end subroutine final_syntax_tables
    module subroutine write_syntax_tables ()
    end subroutine write_syntax_tables
    module subroutine whizard_process_ifile (whizard, ifile, quit, quit_code)
      class(whizard_t), intent(inout), target :: whizard
      type(ifile_t), intent(in) :: ifile
      logical, intent(out) :: quit
      integer, intent(out) :: quit_code
    end subroutine whizard_process_ifile
    module subroutine whizard_process_stdin (whizard, quit, quit_code)
      class(whizard_t), intent(inout), target :: whizard
      logical, intent(out) :: quit
      integer, intent(out) :: quit_code
    end subroutine whizard_process_stdin
    module subroutine whizard_process_file (whizard, file, quit, quit_code)
      class(whizard_t), intent(inout), target :: whizard
      type(string_t), intent(in) :: file
      logical, intent(out) :: quit
      integer, intent(out) :: quit_code
    end subroutine whizard_process_file
    module subroutine whizard_process_stream &
         (whizard, stream, lexer, quit, quit_code)
      class(whizard_t), intent(inout), target :: whizard
      type(stream_t), intent(inout), target :: stream
      type(lexer_t), intent(inout), target :: lexer
      logical, intent(out) :: quit
      integer, intent(out) :: quit_code
    end subroutine whizard_process_stream
    module subroutine whizard_shell (whizard, quit_code)
      class(whizard_t), intent(inout), target :: whizard
      integer, intent(out) :: quit_code
    end subroutine whizard_shell
  end interface

end module whizard
