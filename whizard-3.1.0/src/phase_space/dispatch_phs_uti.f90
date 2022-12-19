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

module dispatch_phs_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use variables
  use io_units, only: free_unit
  use os_interface, only: os_data_t
  use process_constants
  use model_data
  use models
  use phs_base
  use phs_none
  use phs_forests
  use phs_wood
  use mappings
  use dispatch_phase_space

  implicit none
  private

  public :: dispatch_phs_1
  public :: dispatch_phs_2

contains

  subroutine dispatch_phs_1 (u)
    integer, intent(in) :: u
    type(var_list_t) :: var_list
    class(phs_config_t), allocatable :: phs
    type(phs_parameters_t) :: phs_par
    type(os_data_t) :: os_data
    type(mapping_defaults_t) :: mapping_defs

    write (u, "(A)")  "* Test output: dispatch_phs_1"
    write (u, "(A)")  "*   Purpose: select phase-space configuration method"
    write (u, "(A)")

    call var_list%init_defaults (0)

    write (u, "(A)")  "* Allocate PHS as phs_none_t"
    write (u, "(A)")

    call var_list%set_string (&
         var_str ("$phs_method"), &
         var_str ("none"), is_known = .true.)
    call dispatch_phs (phs, var_list, os_data, var_str ("dispatch_phs_1"))
    call phs%write (u)

    call phs%final ()
    deallocate (phs)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate PHS as phs_single_t"
    write (u, "(A)")

    call var_list%set_string (&
         var_str ("$phs_method"), &
         var_str ("single"), is_known = .true.)
    call dispatch_phs (phs, var_list, os_data, var_str ("dispatch_phs_1"))
    call phs%write (u)

    call phs%final ()
    deallocate (phs)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate PHS as phs_wood_t"
    write (u, "(A)")

    call var_list%set_string (&
         var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call dispatch_phs (phs, var_list, os_data, var_str ("dispatch_phs_1"))
    call phs%write (u)

    call phs%final ()
    deallocate (phs)

    write (u, "(A)")
    write (u, "(A)")  "* Setting parameters for phs_wood_t"
    write (u, "(A)")

    phs_par%m_threshold_s = 123
    phs_par%m_threshold_t = 456
    phs_par%t_channel = 42
    phs_par%off_shell = 17
    phs_par%keep_nonresonant = .false.
    mapping_defs%energy_scale = 987
    mapping_defs%invariant_mass_scale = 654
    mapping_defs%momentum_transfer_scale = 321
    mapping_defs%step_mapping = .false.
    mapping_defs%step_mapping_exp = .false.
    mapping_defs%enable_s_mapping = .true.
    call dispatch_phs (phs, var_list, os_data, var_str ("dispatch_phs_1"), &
         mapping_defs, phs_par)
    call phs%write (u)

    call phs%final ()

    call var_list%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_phs_1"

  end subroutine dispatch_phs_1

  subroutine dispatch_phs_2 (u)
    use phs_base_ut, only: init_test_process_data
    use phs_wood_ut, only: write_test_phs_file
    use phs_forests
    integer, intent(in) :: u
    type(var_list_t) :: var_list
    type(os_data_t) :: os_data
    type(process_constants_t) :: process_data
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    class(phs_config_t), allocatable :: phs
    integer :: u_phs

    write (u, "(A)")  "* Test output: dispatch_phs_2"
    write (u, "(A)")  "*   Purpose: select 'wood' phase-space &
         &for a test process"
    write (u, "(A)")  "*            and read phs configuration from file"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize a process"
    write (u, "(A)")

    call var_list%init_defaults (0)
    call os_data%init ()
    call syntax_model_file_init ()
    call model_list%read_model &
         (var_str ("Test"), var_str ("Test.mdl"), os_data, model)

    call syntax_phs_forest_init ()

    call init_test_process_data (var_str ("dispatch_phs_2"), process_data)

    write (u, "(A)")  "* Write phase-space file"

    u_phs = free_unit ()
    open (u_phs, file = "dispatch_phs_2.phs", action = "write", status = "replace")
    call write_test_phs_file (u_phs, var_str ("dispatch_phs_2"))
    close (u_phs)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate PHS as phs_wood_t"
    write (u, "(A)")

    call var_list%set_string (&
         var_str ("$phs_method"), &
         var_str ("wood"), is_known = .true.)
    call var_list%set_string (&
         var_str ("$phs_file"), &
         var_str ("dispatch_phs_2.phs"), is_known = .true.)
    call dispatch_phs (phs, var_list, os_data, var_str ("dispatch_phs_2"))

    call phs%init (process_data, model)
    call phs%configure (sqrts = 1000._default)

    call phs%write (u)
    write (u, "(A)")
    select type (phs)
    type is (phs_wood_config_t)
       call phs%write_forest (u)
    end select

    call phs%final ()

    call var_list%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_phs_2"

  end subroutine dispatch_phs_2


end module dispatch_phs_uti
