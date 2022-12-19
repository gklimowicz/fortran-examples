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

module phs_rambo_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use flavors
  use lorentz
  use model_data
  use process_constants
  use phs_base

  use phs_rambo

  use phs_base_ut, only: init_test_process_data, init_test_decay_data

  implicit none
  private

  public :: phs_rambo_1
  public :: phs_rambo_2
  public :: phs_rambo_3
  public :: phs_rambo_4

contains

  subroutine phs_rambo_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    class(phs_config_t), allocatable :: phs_data
    real(default) :: sqrts

    write (u, "(A)")  "* Test output: phs_rambo_1"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &phase-space configuration data"
    write (u, "(A)")

    call model%init_test ()

    write (u, "(A)")  "* Initialize a process and a matching &
         &phase-space configuration"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_rambo_1"), process_data)

    allocate (phs_rambo_config_t :: phs_data)
    call phs_data%init (process_data, model)

    sqrts = 1000._default
    call phs_data%configure (sqrts)

    call phs_data%write (u)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_rambo_1"

  end subroutine phs_rambo_1

  subroutine phs_rambo_2 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(process_constants_t) :: process_data
    real(default) :: sqrts, E
    class(phs_config_t), allocatable, target :: phs_data
    class(phs_t), pointer :: phs => null ()
    type(vector4_t), dimension(2) :: p, q

    write (u, "(A)")  "* Test output: phs_rambo_2"
    write (u, "(A)")  "*   Purpose: test simple two-channel phase space"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)

    write (u, "(A)")  "* Initialize a process and a matching &
         &phase-space configuration"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_rambo_2"), process_data)

    allocate (phs_rambo_config_t :: phs_data)
    call phs_data%init (process_data, model)

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
         &for x = 0.5, 0.125"
    write (u, "(A)")

    call phs%evaluate_selected_channel (1, [0.5_default, 0.125_default])
    call phs%evaluate_other_channels (1)
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Inverse kinematics"
    write (u, "(A)")

    call phs%get_outgoing_momenta (q)
    deallocate (phs)
    call phs_data%allocate_instance (phs)
    call phs%init (phs_data)

    sqrts = 1000._default
    call phs_data%configure (sqrts)

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()
    call phs%set_outgoing_momenta (q)

    call phs%inverse ()
    call phs%write (u)

    call phs%final ()
    deallocate (phs)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_rambo_2"

  end subroutine phs_rambo_2

  subroutine phs_rambo_3 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(process_constants_t) :: process_data
    real(default) :: sqrts, E
    class(phs_config_t), allocatable, target :: phs_data
    class(phs_t), pointer :: phs => null ()
    type(vector4_t), dimension(2) :: p, q
    type(lorentz_transformation_t) :: lt

    write (u, "(A)")  "* Test output: phs_rambo_3"
    write (u, "(A)")  "*   Purpose: phase-space evaluation in lab frame"
    write (u, "(A)")

    call model%init_test ()
    call flv%init (25, model)

    write (u, "(A)")  "* Initialize a process and a matching &
         &phase-space configuration"
    write (u, "(A)")

    call init_test_process_data (var_str ("phs_rambo_3"), process_data)

    allocate (phs_rambo_config_t :: phs_data)
    call phs_data%init (process_data, model)

    sqrts = 1000._default
    call phs_data%configure (sqrts, lab_is_cm=.false., sqrts_fixed=.false.)

    call phs_data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize the phase-space instance"
    write (u, "(A)")

    call phs_data%allocate_instance (phs)
    call phs%init (phs_data)

    call phs%write (u, verbose=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Set incoming momenta in lab system"
    write (u, "(A)")

    lt = boost (0.1_default, 1) * boost (0.3_default, 3)

    E = sqrts / 2
    p(1) = lt * vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    p(2) = lt * vector4_moving (E,-sqrt (E**2 - flv%get_mass ()**2), 3)

    call vector4_write (p(1), u)
    call vector4_write (p(2), u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute phase-space point &
         &for x = 0.5, 0.125"
    write (u, "(A)")

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()

    call phs%evaluate_selected_channel (1, [0.5_default, 0.125_default])
    call phs%evaluate_other_channels (1)
    call pacify (phs)
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Extract outgoing momenta in lab system"
    write (u, "(A)")

    call phs%get_outgoing_momenta (q)
    call vector4_write (q(1), u)
    call vector4_write (q(2), u)

    write (u, "(A)")
    write (u, "(A)")  "* Inverse kinematics"
    write (u, "(A)")

    deallocate (phs)
    call phs_data%allocate_instance (phs)
    call phs%init (phs_data)

    sqrts = 1000._default
    call phs_data%configure (sqrts)

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()
    call phs%set_outgoing_momenta (q)

    call phs%inverse ()
    call pacify (phs)
    call phs%write (u)

    call phs%final ()
    deallocate (phs)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_rambo_3"

  end subroutine phs_rambo_3

  subroutine phs_rambo_4 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(process_constants_t) :: process_data
    class(phs_config_t), allocatable, target :: phs_data
    class(phs_t), pointer :: phs => null ()
    type(vector4_t), dimension(1) :: p
    type(vector4_t), dimension(2) :: q

    write (u, "(A)")  "* Test output: phs_rambo_4"
    write (u, "(A)")  "*   Purpose: test simple two-channel phase space"
    write (u, "(A)")

    call model%init_test ()

    call model%set_par (var_str ("ff"), 0.4_default)
    call model%set_par (var_str ("mf"), &
         model%get_real (var_str ("ff")) * model%get_real (var_str ("ms")))
    call flv%init (25, model)

    write (u, "(A)")  "* Initialize a decay and a matching &
         &phase-space configuration"
    write (u, "(A)")

    call init_test_decay_data (var_str ("phs_rambo_4"), process_data)

    allocate (phs_rambo_config_t :: phs_data)
    call phs_data%init (process_data, model)

    call phs_data%configure (flv%get_mass ())

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

    p(1) = vector4_at_rest (flv%get_mass ())

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute phase-space point &
         &for x = 0.5, 0.125"
    write (u, "(A)")

    call phs%evaluate_selected_channel (1, [0.5_default, 0.125_default])
    call phs%evaluate_other_channels (1)
    call phs%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Inverse kinematics"
    write (u, "(A)")

    call phs%get_outgoing_momenta (q)
    deallocate (phs)
    call phs_data%allocate_instance (phs)
    call phs%init (phs_data)

    call phs_data%configure (flv%get_mass ())

    call phs%set_incoming_momenta (p)
    call phs%compute_flux ()
    call phs%set_outgoing_momenta (q)

    call phs%inverse ()
    call phs%write (u)

    call phs%final ()
    deallocate (phs)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_rambo_4"

  end subroutine phs_rambo_4


end module phs_rambo_uti
