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

module restricted_subprocesses

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use particle_specifiers, only: prt_spec_t
  use resonances, only: resonance_history_t, resonance_history_set_t
  use variables, only: var_list_t
  use models, only: model_t
  use event_transforms, only: evt_t
  use resonance_insertion, only: evt_resonance_t
  use rt_data, only: rt_data_t
  use process_configurations, only: process_configuration_t
  use process, only: process_t, process_ptr_t
  use instances, only: process_instance_t, process_instance_ptr_t

  implicit none
  private

  public :: restricted_process_configuration_t
  public :: resonant_subprocess_set_t
  public :: get_libname_res
  public :: spawn_resonant_subprocess_libraries

  type, extends (process_configuration_t) :: restricted_process_configuration_t
     private
   contains
     procedure :: init_resonant_process
  end type restricted_process_configuration_t

  type :: resonant_subprocess_set_t
     private
     integer, dimension(:), allocatable :: n_history
     type(resonance_history_set_t), dimension(:), allocatable :: res_history_set
     logical :: lib_active = .false.
     type(string_t) :: libname
     type(string_t), dimension(:), allocatable :: proc_id
     type(process_ptr_t), dimension(:), allocatable :: subprocess
     type(process_instance_ptr_t), dimension(:), allocatable :: instance
     logical :: filled = .false.
     type(evt_resonance_t), pointer :: evt => null ()
   contains
     procedure :: write => resonant_subprocess_set_write
     procedure :: init => resonant_subprocess_set_init
     procedure :: fill_resonances => resonant_subprocess_set_fill_resonances
     procedure :: get_resonance_history_set &
          => resonant_subprocess_set_get_resonance_history_set
     procedure :: create_library => resonant_subprocess_set_create_library
     procedure :: add_to_library => resonant_subprocess_set_add_to_library
     procedure :: freeze_library => resonant_subprocess_set_freeze_library
     procedure :: compile_library => resonant_subprocess_set_compile_library
     procedure :: is_active => resonant_subprocess_set_is_active
     procedure :: get_n_process => resonant_subprocess_set_get_n_process
     procedure :: get_libname => resonant_subprocess_set_get_libname
     procedure :: get_proc_id => resonant_subprocess_set_get_proc_id
     procedure :: prepare_process_objects &
          => resonant_subprocess_set_prepare_process_objects
     procedure :: prepare_process_instances &
          => resonant_subprocess_set_prepare_process_instances
     procedure :: connect_transform => &
          resonant_subprocess_set_connect_transform
     procedure :: set_on_shell_limit => resonant_subprocess_set_on_shell_limit
     procedure :: set_on_shell_turnoff => resonant_subprocess_set_on_shell_turnoff
     procedure :: set_background_factor &
          => resonant_subprocess_set_background_factor
     procedure :: dump_instances => resonant_subprocess_set_dump_instances
     procedure :: fill_momenta => resonant_subprocess_set_fill_momenta
     procedure :: determine_on_shell_histories &
          => resonant_subprocess_set_determine_on_shell_histories
     procedure :: evaluate_subprocess &
          => resonant_subprocess_set_evaluate_subprocess
     procedure :: get_master_sqme &
          => resonant_subprocess_set_get_master_sqme
     procedure :: get_subprocess_sqme &
          => resonant_subprocess_set_get_subprocess_sqme
     procedure :: compute_probabilities &
          => resonant_subprocess_set_compute_probabilities
  end type resonant_subprocess_set_t




  interface
    module subroutine init_resonant_process &
         (prc_config, prc_name, prt_in, prt_out, res_history, model, var_list)
      class(restricted_process_configuration_t), intent(out) :: prc_config
      type(string_t), intent(in) :: prc_name
      type(prt_spec_t), dimension(:), intent(in) :: prt_in
      type(prt_spec_t), dimension(:), intent(in) :: prt_out
      type(resonance_history_t), intent(in) :: res_history
      type(model_t), intent(in), target :: model
      type(var_list_t), intent(in), target :: var_list
    end subroutine init_resonant_process
    module subroutine resonant_subprocess_set_write (prc_set, unit, testflag)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine resonant_subprocess_set_write
    module subroutine resonant_subprocess_set_init (prc_set, n_component)
      class(resonant_subprocess_set_t), intent(out) :: prc_set
      integer, intent(in) :: n_component
    end subroutine resonant_subprocess_set_init
    module subroutine resonant_subprocess_set_fill_resonances (prc_set, &
         res_history_set, i_component)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      type(resonance_history_set_t), intent(in) :: res_history_set
      integer, intent(in) :: i_component
    end subroutine resonant_subprocess_set_fill_resonances
    module function resonant_subprocess_set_get_resonance_history_set &
         (prc_set) result (res_history_set)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      type(resonance_history_set_t), dimension(:), allocatable :: &
           res_history_set
    end function resonant_subprocess_set_get_resonance_history_set
    elemental module function get_libname_res (proc_id) result (libname)
      type(string_t), intent(in) :: proc_id
      type(string_t) :: libname
    end function get_libname_res
    module subroutine spawn_resonant_subprocess_libraries &
         (libname, local, global, libname_res)
      type(string_t), intent(in) :: libname
      type(rt_data_t), intent(inout), target :: local
      type(rt_data_t), intent(inout), target :: global
      type(string_t), dimension(:), allocatable, intent(inout) :: libname_res
    end subroutine spawn_resonant_subprocess_libraries
    module subroutine resonant_subprocess_set_create_library (prc_set, &
         libname, global, exist)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      type(string_t), intent(in) :: libname
      type(rt_data_t), intent(inout), target :: global
      logical, intent(out) :: exist
    end subroutine resonant_subprocess_set_create_library
    module subroutine resonant_subprocess_set_add_to_library (prc_set, &
         i_component, prt_in, prt_out, global)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      integer, intent(in) :: i_component
      type(prt_spec_t), dimension(:), intent(in) :: prt_in
      type(prt_spec_t), dimension(:), intent(in) :: prt_out
      type(rt_data_t), intent(inout), target :: global
    end subroutine resonant_subprocess_set_add_to_library
    module subroutine resonant_subprocess_set_freeze_library (prc_set, global)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      type(rt_data_t), intent(inout), target :: global
    end subroutine resonant_subprocess_set_freeze_library
    module subroutine resonant_subprocess_set_compile_library (prc_set, global)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      type(rt_data_t), intent(inout), target :: global
    end subroutine resonant_subprocess_set_compile_library
    module function resonant_subprocess_set_is_active (prc_set) result (flag)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      logical :: flag
    end function resonant_subprocess_set_is_active
    module function resonant_subprocess_set_get_n_process (prc_set) result (n)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      integer :: n
    end function resonant_subprocess_set_get_n_process
    module function resonant_subprocess_set_get_libname &
         (prc_set) result (libname)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      type(string_t) :: libname
    end function resonant_subprocess_set_get_libname
    module function resonant_subprocess_set_get_proc_id &
         (prc_set, i) result (proc_id)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      integer, intent(in) :: i
      type(string_t) :: proc_id
    end function resonant_subprocess_set_get_proc_id
    module subroutine resonant_subprocess_set_prepare_process_objects &
         (prc_set, local, global)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      type(rt_data_t), intent(inout), target :: local
      type(rt_data_t), intent(inout), optional, target :: global
    end subroutine resonant_subprocess_set_prepare_process_objects
    module subroutine resonant_subprocess_set_prepare_process_instances &
         (prc_set, global)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      type(rt_data_t), intent(in), target :: global
    end subroutine resonant_subprocess_set_prepare_process_instances
    module subroutine resonant_subprocess_set_connect_transform (prc_set, evt)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      class(evt_t), intent(in), target :: evt
    end subroutine resonant_subprocess_set_connect_transform
    module subroutine resonant_subprocess_set_on_shell_limit &
         (prc_set, on_shell_limit)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      real(default), intent(in) :: on_shell_limit
    end subroutine resonant_subprocess_set_on_shell_limit
    module subroutine resonant_subprocess_set_on_shell_turnoff &
         (prc_set, on_shell_turnoff)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      real(default), intent(in) :: on_shell_turnoff
    end subroutine resonant_subprocess_set_on_shell_turnoff
    module subroutine resonant_subprocess_set_background_factor &
         (prc_set, background_factor)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      real(default), intent(in) :: background_factor
    end subroutine resonant_subprocess_set_background_factor
    module subroutine resonant_subprocess_set_dump_instances &
         (prc_set, unit, testflag)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine resonant_subprocess_set_dump_instances
    module subroutine resonant_subprocess_set_fill_momenta (prc_set)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
    end subroutine resonant_subprocess_set_fill_momenta
    module subroutine resonant_subprocess_set_determine_on_shell_histories &
         (prc_set, i_component, index_array)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      integer, intent(in) :: i_component
      integer, dimension(:), allocatable, intent(out) :: index_array
    end subroutine resonant_subprocess_set_determine_on_shell_histories
    module subroutine resonant_subprocess_set_evaluate_subprocess &
         (prc_set, index_array)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      integer, dimension(:), intent(in) :: index_array
    end subroutine resonant_subprocess_set_evaluate_subprocess
    module function resonant_subprocess_set_get_master_sqme &
         (prc_set) result (sqme)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      real(default) :: sqme
    end function resonant_subprocess_set_get_master_sqme
    module subroutine resonant_subprocess_set_get_subprocess_sqme &
         (prc_set, sqme)
      class(resonant_subprocess_set_t), intent(in) :: prc_set
      real(default), dimension(:), intent(inout) :: sqme
    end subroutine resonant_subprocess_set_get_subprocess_sqme
    module subroutine resonant_subprocess_set_compute_probabilities &
         (prc_set, prob_array)
      class(resonant_subprocess_set_t), intent(inout) :: prc_set
      real(default), dimension(:), allocatable, intent(out) :: prob_array
    end subroutine resonant_subprocess_set_compute_probabilities
  end interface

end module restricted_subprocesses
