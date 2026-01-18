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

module phs_wood_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use os_interface
  use lorentz
  use flavors
  use model_data
  use process_constants
  use mappings
  use phs_base
  use phs_forests

  use phs_wood

  use phs_base_ut, only: init_test_process_data, init_test_decay_data

  implicit none
  private

  public :: write_test_phs_file

  public :: phs_wood_1
  public :: phs_wood_2
  public :: phs_wood_3
  public :: phs_wood_4
  public :: phs_wood_5
  public :: phs_wood_6
  public :: phs_wood_vis_1

contains

  subroutine phs_wood_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    class(phs_config_t), allocatable :: phs_data
    type(mapping_defaults_t) :: mapping_defaults
    real(default) :: sqrts
    integer :: u_phs, iostat
    character(32) :: buffer

    write (u, "(A)")  "* Test output: phs_wood_1"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &phase-space configuration data"
    write (u, "(A)")

    call model%init_test ()

    call syntax_phs_forest_init ()

    write (u, "(A)")  "* Initialize a process"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_wood_1"), process_data)

    write (u, "(A)")  "* Create a scratch phase-space file"
    write (u, "(A)")

    u_phs = free_unit ()
    open (u_phs, status = "scratch", action = "readwrite")
    call write_test_phs_file (u_phs, var_str ("phs_wood_1"))
    rewind (u_phs)
    do
       read (u_phs, "(A)", iostat = iostat)  buffer
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Setup phase-space configuration object"
    write (u, "(A)")

    mapping_defaults%step_mapping = .false.

    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_input (u_phs)
       call phs_data%set_mapping_defaults (mapping_defaults)
    end select

    sqrts = 1000._default
    call phs_data%configure (sqrts)

    call phs_data%write (u)
    write (u, "(A)")

    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%write_forest (u)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    close (u_phs)
    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_wood_1"

  end subroutine phs_wood_1

  subroutine phs_wood_2 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(process_constants_t) :: process_data
    real(default) :: sqrts, E
    class(phs_config_t), allocatable, target :: phs_data
    class(phs_t), pointer :: phs => null ()
    type(vector4_t), dimension(2) :: p, q
    integer :: u_phs

    write (u, "(A)")  "* Test output: phs_wood_2"
    write (u, "(A)")  "*   Purpose: test simple single-channel phase space"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)

    write (u, "(A)")  "* Initialize a process and a matching &
         &phase-space configuration"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_wood_2"), process_data)
    u_phs = free_unit ()
    open (u_phs, status = "scratch", action = "readwrite")
    call write_test_phs_file (u_phs, var_str ("phs_wood_2"))
    rewind (u_phs)

    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_input (u_phs)
    end select

    sqrts = 1000._default
    call phs_data%configure (sqrts)

    call phs_data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize the phase-space instance"
    write (u, "(A)")

    call phs_data%allocate_instance (phs)
    call phs%init (phs_data)

    call phs%write (u, verbose=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Set incoming momenta"
    write (u, "(A)")

    E = sqrts / 2
    p(1) = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    p(2) = vector4_moving (E,-sqrt (E**2 - flv%get_mass ()**2), 3)

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute phase-space point &
         &for x = 0.125, 0.5"
    write (u, "(A)")

    call phs%evaluate_selected_channel (1, [0.125_default, 0.5_default])
    call phs%evaluate_other_channels (1)
    call phs%write (u)
    write (u, "(A)")
    select type (phs)
    type is (phs_wood_t)
       call phs%write_forest (u)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Inverse kinematics"
    write (u, "(A)")

    call phs%get_outgoing_momenta (q)
    call phs%final ()
    deallocate (phs)

    call phs_data%allocate_instance (phs)
    call phs%init (phs_data)

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()
    call phs%set_outgoing_momenta (q)

    call phs%inverse ()
    call phs%write (u)
    write (u, "(A)")
    select type (phs)
    type is (phs_wood_t)
       call phs%write_forest (u)
    end select

    call phs%final ()
    deallocate (phs)

    close (u_phs)
    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_wood_2"

  end subroutine phs_wood_2

  subroutine phs_wood_3 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    type(phs_parameters_t) :: phs_par
    class(phs_config_t), allocatable :: phs_data
    integer :: iostat
    character(80) :: buffer

    write (u, "(A)")  "* Test output: phs_wood_3"
    write (u, "(A)")  "*   Purpose: generate a phase-space configuration"
    write (u, "(A)")

    call model%init_test ()

    call syntax_phs_forest_init ()

    write (u, "(A)")  "* Initialize a process and phase-space parameters"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_wood_3"), process_data)
    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)

    phs_par%sqrts = 1000
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_parameters (phs_par)
       phs_data%io_unit_keep_open = .true.
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Generate a scratch phase-space file"
    write (u, "(A)")

    call phs_data%configure (phs_par%sqrts)

    select type (phs_data)
    type is (phs_wood_config_t)
       rewind (phs_data%io_unit)
       do
          read (phs_data%io_unit, "(A)", iostat = iostat)  buffer
          if (iostat /= 0)  exit
          write (u, "(A)") trim (buffer)
       end do
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_wood_3"

  end subroutine phs_wood_3

  subroutine phs_wood_4 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    type(phs_parameters_t) :: phs_par
    class(phs_config_t), allocatable, target :: phs_data
    integer :: iostat
    character(80) :: buffer
    class(phs_t), pointer :: phs => null ()
    real(default) :: E, pL
    type(vector4_t), dimension(2) :: p
    type(vector4_t), dimension(3) :: q

    write (u, "(A)")  "* Test output: phs_wood_4"
    write (u, "(A)")  "*   Purpose: generate a phase-space configuration"
    write (u, "(A)")

    call model%init_test ()

    call syntax_phs_forest_init ()

    write (u, "(A)")  "* Initialize a process and phase-space parameters"
    write (u, "(A)")

    process_data%id = "phs_wood_4"
    process_data%model_name = "Test"
    process_data%n_in = 2
    process_data%n_out = 3
    process_data%n_flv = 1
    allocate (process_data%flv_state (process_data%n_in + process_data%n_out, &
         process_data%n_flv))
    process_data%flv_state(:,1) = [25, 25, 25, 6, -6]

    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)

    phs_par%sqrts = 1000
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_parameters (phs_par)
       phs_data%io_unit_keep_open = .true.
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Generate a scratch phase-space file"
    write (u, "(A)")

    call phs_data%configure (phs_par%sqrts)

    select type (phs_data)
    type is (phs_wood_config_t)
       rewind (phs_data%io_unit)
       do
          read (phs_data%io_unit, "(A)", iostat = iostat)  buffer
          if (iostat /= 0)  exit
          write (u, "(A)") trim (buffer)
       end do
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Initialize the phase-space instance"
    write (u, "(A)")

    call phs_data%allocate_instance (phs)
    call phs%init (phs_data)

    write (u, "(A)")  "* Set incoming momenta"
    write (u, "(A)")

    select type (phs_data)
    type is (phs_wood_config_t)
       E = phs_data%sqrts / 2
       pL = sqrt (E**2 - phs_data%flv(1,1)%get_mass ()**2)
    end select
    p(1) = vector4_moving (E, pL, 3)
    p(2) = vector4_moving (E, -pL, 3)

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()

    write (u, "(A)")  "* Compute phase-space point &
         &for x = 0.1, 0.2, 0.3, 0.4, 0.5"
    write (u, "(A)")

    call phs%evaluate_selected_channel (1, &
         [0.1_default, 0.2_default, 0.3_default, 0.4_default, 0.5_default])
    call phs%evaluate_other_channels (1)
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Inverse kinematics"
    write (u, "(A)")

    call phs%get_outgoing_momenta (q)
    call phs%final ()
    deallocate (phs)

    call phs_data%allocate_instance (phs)
    call phs%init (phs_data)

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()
    call phs%set_outgoing_momenta (q)

    call phs%inverse ()
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call phs%final ()
    deallocate (phs)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_wood_4"

  end subroutine phs_wood_4

  subroutine phs_wood_5 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    type(phs_parameters_t) :: phs_par
    class(phs_config_t), allocatable :: phs_data

    write (u, "(A)")  "* Test output: phs_wood_5"
    write (u, "(A)")  "*   Purpose: generate a phase-space configuration"
    write (u, "(A)")

    call model%init_test ()

    call syntax_phs_forest_init ()

    write (u, "(A)")  "* Initialize a process and phase-space parameters"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_wood_5"), process_data)
    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)

    phs_par%sqrts = 1000
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_parameters (phs_par)
       call phs_data%enable_equivalences ()
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Generate a scratch phase-space file"
    write (u, "(A)")

    call phs_data%configure (phs_par%sqrts)
    call phs_data%write (u)
    write (u, "(A)")

    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%write_forest (u)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_wood_5"

  end subroutine phs_wood_5

  subroutine phs_wood_6 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    type(phs_parameters_t) :: phs_par
    class(phs_config_t), allocatable :: phs_data
    logical :: exist, found, match
    integer :: u_phs
    character(*), parameter :: filename = "phs_wood_6_p.phs"

    write (u, "(A)")  "* Test output: phs_wood_6"
    write (u, "(A)")  "*   Purpose: generate and check  phase-space file"
    write (u, "(A)")

    call model%init_test ()

    call syntax_phs_forest_init ()

    write (u, "(A)")  "* Initialize a process and phase-space parameters"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_wood_6"), process_data)
    process_data%id = "phs_wood_6_p"
    process_data%md5sum = "1234567890abcdef1234567890abcdef"
    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)

    phs_par%sqrts = 1000
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_parameters (phs_par)
    end select

    write (u, "(A)")  "* Remove previous phs file, if any"
    write (u, "(A)")

    inquire (file = filename, exist = exist)
    if (exist) then
       u_phs = free_unit ()
       open (u_phs, file = filename, action = "write")
       close (u_phs, status = "delete")
    end if

    write (u, "(A)")  "* Check phase-space file (should fail)"
    write (u, "(A)")

    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%read_phs_file (exist, found, match)
       write (u, "(1x,A,L1)")  "exist = ", exist
       write (u, "(1x,A,L1)")  "found = ", found
       write (u, "(1x,A,L1)")  "match = ", match
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Generate a phase-space file"
    write (u, "(A)")

    call phs_data%configure (phs_par%sqrts)

    write (u, "(1x,A,A,A)")  "MD5 sum (process)    = '", &
         phs_data%md5sum_process, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (model par)  = '", &
         phs_data%md5sum_model_par, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (phs config) = '", &
         phs_data%md5sum_phs_config, "'"

    write (u, "(A)")
    write (u, "(A)")  "* Check MD5 sum"
    write (u, "(A)")

    call phs_data%final ()
    deallocate (phs_data)
    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)
    phs_par%sqrts = 1000
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_parameters (phs_par)
       phs_data%sqrts = phs_par%sqrts
       phs_data%par%sqrts = phs_par%sqrts
    end select
    call phs_data%compute_md5sum ()

    write (u, "(1x,A,A,A)")  "MD5 sum (process)    = '", &
         phs_data%md5sum_process, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (model par)  = '", &
         phs_data%md5sum_model_par, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (phs config) = '", &
         phs_data%md5sum_phs_config, "'"

    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%read_phs_file (exist, found, match)
       write (u, "(1x,A,L1)")  "exist = ", exist
       write (u, "(1x,A,L1)")  "found = ", found
       write (u, "(1x,A,L1)")  "match = ", match
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Modify sqrts and check MD5 sum"
    write (u, "(A)")

    call phs_data%final ()
    deallocate (phs_data)
    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)
    phs_par%sqrts = 500
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_parameters (phs_par)
       phs_data%sqrts = phs_par%sqrts
       phs_data%par%sqrts = phs_par%sqrts
    end select
    call phs_data%compute_md5sum ()

    write (u, "(1x,A,A,A)")  "MD5 sum (process)    = '", &
         phs_data%md5sum_process, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (model par)  = '", &
         phs_data%md5sum_model_par, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (phs config) = '", &
         phs_data%md5sum_phs_config, "'"

    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%read_phs_file (exist, found, match)
       write (u, "(1x,A,L1)")  "exist = ", exist
       write (u, "(1x,A,L1)")  "found = ", found
       write (u, "(1x,A,L1)")  "match = ", match
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Modify process and check MD5 sum"
    write (u, "(A)")

    call phs_data%final ()
    deallocate (phs_data)
    process_data%md5sum = "77777777777777777777777777777777"
    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)
    phs_par%sqrts = 1000
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_parameters (phs_par)
       phs_data%sqrts = phs_par%sqrts
       phs_data%par%sqrts = phs_par%sqrts
    end select
    call phs_data%compute_md5sum ()

    write (u, "(1x,A,A,A)")  "MD5 sum (process)    = '", &
         phs_data%md5sum_process, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (model par)  = '", &
         phs_data%md5sum_model_par, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (phs config) = '", &
         phs_data%md5sum_phs_config, "'"

    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%read_phs_file (exist, found, match)
       write (u, "(1x,A,L1)")  "exist = ", exist
       write (u, "(1x,A,L1)")  "found = ", found
       write (u, "(1x,A,L1)")  "match = ", match
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Modify phs parameter and check MD5 sum"
    write (u, "(A)")

    call phs_data%final ()
    deallocate (phs_data)
    allocate (phs_wood_config_t :: phs_data)
    process_data%md5sum = "1234567890abcdef1234567890abcdef"
    call phs_data%init (process_data, model)
    phs_par%sqrts = 1000
    phs_par%off_shell = 17
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_parameters (phs_par)
       phs_data%sqrts = phs_par%sqrts
       phs_data%par%sqrts = phs_par%sqrts
    end select
    call phs_data%compute_md5sum ()

    write (u, "(1x,A,A,A)")  "MD5 sum (process)    = '", &
         phs_data%md5sum_process, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (model par)  = '", &
         phs_data%md5sum_model_par, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (phs config) = '", &
         phs_data%md5sum_phs_config, "'"

    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%read_phs_file (exist, found, match)
       write (u, "(1x,A,L1)")  "exist = ", exist
       write (u, "(1x,A,L1)")  "found = ", found
       write (u, "(1x,A,L1)")  "match = ", match
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Modify model parameter and check MD5 sum"
    write (u, "(A)")

    call phs_data%final ()
    deallocate (phs_data)
    allocate (phs_wood_config_t :: phs_data)
    call model%set_par (var_str ("ms"), 100._default)
    call phs_data%init (process_data, model)
    phs_par%sqrts = 1000
    phs_par%off_shell = 1
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_parameters (phs_par)
       phs_data%sqrts = phs_par%sqrts
       phs_data%par%sqrts = phs_par%sqrts
    end select
    call phs_data%compute_md5sum ()

    write (u, "(1x,A,A,A)")  "MD5 sum (process)    = '", &
         phs_data%md5sum_process, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (model par)  = '", &
         phs_data%md5sum_model_par, "'"
    write (u, "(1x,A,A,A)")  "MD5 sum (phs config) = '", &
         phs_data%md5sum_phs_config, "'"

    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%read_phs_file (exist, found, match)
       write (u, "(1x,A,L1)")  "exist = ", exist
       write (u, "(1x,A,L1)")  "found = ", found
       write (u, "(1x,A,L1)")  "match = ", match
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_wood_6"

  end subroutine phs_wood_6

  subroutine phs_wood_vis_1 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    class(phs_config_t), allocatable :: phs_data
    type(mapping_defaults_t) :: mapping_defaults
    type(string_t) :: vis_file, pdf_file, ps_file
    real(default) :: sqrts
    logical :: exist, exist_pdf, exist_ps
    integer :: u_phs, iostat, u_vis
    character(95) :: buffer

    write (u, "(A)")  "* Test output: phs_wood_vis_1"
    write (u, "(A)")  "*   Purpose: visualizing the &
         &phase-space configuration"
    write (u, "(A)")

    call os_data%init ()
    call model%init_test ()

    call syntax_phs_forest_init ()

    write (u, "(A)")  "* Initialize a process"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_wood_vis_1"), process_data)

    write (u, "(A)")  "* Create a scratch phase-space file"
    write (u, "(A)")

    u_phs = free_unit ()
    open (u_phs, status = "scratch", action = "readwrite")
    call write_test_phs_file (u_phs, var_str ("phs_wood_vis_1"))
    rewind (u_phs)
    do
       read (u_phs, "(A)", iostat = iostat)  buffer
       if (iostat /= 0)  exit
       write (u, "(A)") trim (buffer)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Setup phase-space configuration object"
    write (u, "(A)")

    mapping_defaults%step_mapping = .false.

    allocate (phs_wood_config_t :: phs_data)
    call phs_data%init (process_data, model)
    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%set_input (u_phs)
       call phs_data%set_mapping_defaults (mapping_defaults)
       phs_data%os_data = os_data
       phs_data%io_unit = 0
       phs_data%io_unit_keep_open = .true.
       phs_data%vis_channels = .true.
    end select

    sqrts = 1000._default
    call phs_data%configure (sqrts)

    call phs_data%write (u)
    write (u, "(A)")

    select type (phs_data)
    type is (phs_wood_config_t)
       call phs_data%write_forest (u)
    end select

    vis_file = "phs_wood_vis_1.phs-vis.tex"
    ps_file  = "phs_wood_vis_1.phs-vis.ps"
    pdf_file = "phs_wood_vis_1.phs-vis.pdf"
    inquire (file = char (vis_file), exist = exist)
    if (exist) then
       u_vis = free_unit ()
       open (u_vis, file = char (vis_file), action = "read", status = "old")
       iostat = 0
       do while (iostat == 0)
          read (u_vis, "(A)", iostat = iostat)  buffer
          if (iostat == 0)  write (u, "(A)")  trim (buffer)
       end do
       close (u_vis)
    else
       write (u, "(A)")  "[Visualize LaTeX file is missing]"
    end if
    inquire (file = char (ps_file), exist = exist_ps)
    if (exist_ps) then
       write (u, "(A)")  "[Visualize Postscript file exists and is nonempty]"
    else
       write (u, "(A)")  "[Visualize Postscript file is missing/non-regular]"
    end if
    inquire (file = char (pdf_file), exist = exist_pdf)
    if (exist_pdf) then
       write (u, "(A)")  "[Visualize PDF file exists and is nonempty]"
    else
       write (u, "(A)")  "[Visualize PDF file is missing/non-regular]"
    end if

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    close (u_phs)
    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_wood_vis_1"

  end subroutine phs_wood_vis_1


  subroutine write_test_phs_file (u_phs, procname)
    integer, intent(in) :: u_phs
    type(string_t), intent(in), optional :: procname
    if (present (procname)) then
       write (u_phs, "(A,A)")  "process ", char (procname)
    else
       write (u_phs, "(A)")  "process testproc"
    end if
    write (u_phs, "(A,A)")  "   md5sum_process    = ", '""'
    write (u_phs, "(A,A)")  "   md5sum_model_par  = ", '""'
    write (u_phs, "(A,A)")  "   md5sum_phs_config = ", '""'
    write (u_phs, "(A)")  "   sqrts         = 1000"
    write (u_phs, "(A)")  "   m_threshold_s =   50"
    write (u_phs, "(A)")  "   m_threshold_t =  100"
    write (u_phs, "(A)")  "   off_shell = 2"
    write (u_phs, "(A)")  "   t_channel = 6"
    write (u_phs, "(A)")  "   keep_nonresonant = T"
    write (u_phs, "(A)")  "  grove #1"
    write (u_phs, "(A)")  "    tree 3"
  end subroutine write_test_phs_file


end module phs_wood_uti
