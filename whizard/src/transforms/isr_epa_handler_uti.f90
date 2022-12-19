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

module isr_epa_handler_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use format_utils, only: write_separator
  use os_interface
  use lorentz, only: vector4_t, vector4_moving, operator(*)
  use rng_base, only: rng_t
  use models, only: syntax_model_file_init, syntax_model_file_final
  use models, only: model_list_t, model_t
  use particles, only: particle_set_t, pacify

  use event_transforms
  use isr_epa_handler, only: evt_isr_epa_t

  use rng_base_ut, only: rng_test_t

  implicit none
  private

  public :: isr_handler_1
  public :: isr_handler_2
  public :: isr_handler_3
  public :: epa_handler_1
  public :: epa_handler_2
  public :: epa_handler_3

contains

  subroutine isr_handler_1 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(evt_trivial_t), target :: evt_trivial
    type(evt_isr_epa_t), target :: evt_isr_epa
    type(vector4_t), dimension(8) :: p
    real(default) :: sqrts
    real(default), dimension(2) :: x, xb
    real(default) :: probability

    write (u, "(A)")  "* Test output: isr_handler_1"
    write (u, "(A)")  "*   Purpose: apply photon handler trivially (no-op)"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment"
    write (u, "(A)")

    call syntax_model_file_init ()
    call os_data%init ()
    call model_list%read_model &
         (var_str ("SM"), var_str ("SM.mdl"), &
         os_data, model)

    write (u, "(A)")  "* Initialize particle set"
    write (u, "(A)")

    call pset%init_direct (n_beam = 2, n_in = 2, n_rem = 2, n_vir = 0, n_out = 2, &
         pdg = [11, -11, 11, -11, 22, 22, 13, -13], model = model)

    sqrts = 100._default
    x = [0.6_default, 0.9_default]
    xb= 1 - x

    p(1) = vector4_moving (sqrts/2, sqrts/2, 3)
    p(2) = vector4_moving (sqrts/2,-sqrts/2, 3)
    p(3:4) = x * p(1:2)
    p(5:6) = xb  * p(1:2)
    p(7:8) = p(3:4)

    call pset%set_momentum (p, on_shell = .false.)

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call pacify (evt_trivial%particle_set)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize ISR handler transform"
    write (u, "(A)")

    evt_trivial%next => evt_isr_epa
    evt_isr_epa%previous => evt_trivial

    call evt_isr_epa%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Fill ISR handler transform"
    write (u, "(A)")

    call evt_isr_epa%prepare_new_event (1, 1)
    call evt_isr_epa%generate_weighted (probability)
    call evt_isr_epa%make_particle_set (0, .false.)

    call evt_isr_epa%write (u)

    write (u, "(A)")
    write (u, "(A,1x,F8.5)")  "Event probability =", probability

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_isr_epa%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: isr_handler_1"

  end subroutine isr_handler_1

  subroutine isr_handler_2 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(evt_trivial_t), target :: evt_trivial
    type(evt_isr_epa_t), target :: evt_isr_epa
    type(vector4_t), dimension(8) :: p
    real(default) :: sqrts
    real(default), dimension(2) :: x, xb
    class(rng_t), allocatable :: rng
    real(default) :: probability

    write (u, "(A)")  "* Test output: isr_handler_2"
    write (u, "(A)")  "*   Purpose: apply photon handler with two-photon recoil"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment"
    write (u, "(A)")

    call syntax_model_file_init ()
    call os_data%init ()
    call model_list%read_model &
         (var_str ("SM"), var_str ("SM.mdl"), &
         os_data, model)

    write (u, "(A)")  "* Initialize particle set"
    write (u, "(A)")

    call pset%init_direct (n_beam = 2, n_in = 2, n_rem = 2, n_vir = 0, n_out = 2, &
         pdg = [11, -11, 11, -11, 22, 22, 13, -13], model = model)

    sqrts = 100._default
    x = [0.6_default, 0.9_default]
    xb= 1 - x

    p(1) = vector4_moving (sqrts/2, sqrts/2, 3)
    p(2) = vector4_moving (sqrts/2,-sqrts/2, 3)
    p(3:4) = x * p(1:2)
    p(5:6) = xb  * p(1:2)
    p(7:8) = p(3:4)

    call pset%set_momentum (p, on_shell = .false.)

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call pacify (evt_trivial%particle_set)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize ISR handler transform"
    write (u, "(A)")

    evt_trivial%next => evt_isr_epa
    evt_isr_epa%previous => evt_trivial

    call evt_isr_epa%set_mode_string (var_str ("recoil"))
    call evt_isr_epa%set_data_isr ( &
         sqrts = sqrts, &
         q_max = sqrts, &
         m = 511.e-3_default, &
         keep_mass = .false. &
         )

    allocate (rng_test_t :: rng)
    call rng%init (3)  ! default would produce pi for azimuthal angle
    call evt_isr_epa%import_rng (rng)

    call evt_isr_epa%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Fill ISR handler transform"
    write (u, "(A)")

    call evt_isr_epa%prepare_new_event (1, 1)
    call evt_isr_epa%generate_weighted (probability)
    call evt_isr_epa%make_particle_set (0, .false.)

    call evt_isr_epa%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A,1x,F8.5)")  "Event probability =", probability

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_isr_epa%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: isr_handler_2"

  end subroutine isr_handler_2

  subroutine isr_handler_3 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(evt_trivial_t), target :: evt_trivial
    type(evt_isr_epa_t), target :: evt_isr_epa
    type(vector4_t), dimension(8) :: p
    real(default) :: sqrts
    real(default), dimension(2) :: x0
    real(default), dimension(2) :: x, xb
    class(rng_t), allocatable :: rng
    real(default) :: probability

    write (u, "(A)")  "* Test output: isr_handler_3"
    write (u, "(A)")  "*   Purpose: apply photon handler for boosted beams &
         &and two-photon recoil"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment"
    write (u, "(A)")

    call syntax_model_file_init ()
    call os_data%init ()
    call model_list%read_model &
         (var_str ("SM"), var_str ("SM.mdl"), &
         os_data, model)

    write (u, "(A)")  "* Initialize particle set"
    write (u, "(A)")

    call pset%init_direct (n_beam = 2, n_in = 2, n_rem = 2, n_vir = 0, n_out = 2, &
         pdg = [11, -11, 11, -11, 22, 22, 13, -13], model = model)

    write (u, "(A)")  "* Event data"
    write (u, "(A)")

    sqrts = 100._default
    write (u, "(A,2(1x,F12.7))")  "sqrts   =", sqrts

    x0 = [0.9_default, 0.4_default]
    write (u, "(A,2(1x,F12.7))")  "x0      =", x0

    write (u, "(A)")
    write (u, "(A,2(1x,F12.7))")  "sqs_hat =", sqrts * sqrt (product (x0))

    x = [0.6_default, 0.9_default]
    xb= 1 - x
    write (u, "(A,2(1x,F12.7))")  "x       =", x

    write (u, "(A)")
    write (u, "(A,2(1x,F12.7))")  "x0 * x  =", x0 * x

    p(1) = x0(1) * vector4_moving (sqrts/2, sqrts/2, 3)
    p(2) = x0(2) * vector4_moving (sqrts/2,-sqrts/2, 3)
    p(3:4) = x * p(1:2)
    p(5:6) = xb  * p(1:2)
    p(7:8) = p(3:4)

    call pset%set_momentum (p, on_shell = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call pacify (evt_trivial%particle_set)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize ISR handler transform"
    write (u, "(A)")

    evt_trivial%next => evt_isr_epa
    evt_isr_epa%previous => evt_trivial

    call evt_isr_epa%set_mode_string (var_str ("recoil"))
    call evt_isr_epa%set_data_isr ( &
         sqrts = sqrts, &
         q_max = sqrts, &
         m = 511.e-3_default, &
         keep_mass = .false. &
         )

    allocate (rng_test_t :: rng)
    call rng%init (3)  ! default would produce pi for azimuthal angle
    call evt_isr_epa%import_rng (rng)

    call evt_isr_epa%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Fill ISR handler transform"
    write (u, "(A)")

    call evt_isr_epa%prepare_new_event (1, 1)
    call evt_isr_epa%generate_weighted (probability)
    call evt_isr_epa%make_particle_set (0, .false.)

    call evt_isr_epa%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A,1x,F8.5)")  "Event probability =", probability

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_isr_epa%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: isr_handler_3"

  end subroutine isr_handler_3

  subroutine epa_handler_1 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(evt_trivial_t), target :: evt_trivial
    type(evt_isr_epa_t), target :: evt_isr_epa
    type(vector4_t), dimension(8) :: p
    real(default) :: sqrts
    real(default), dimension(2) :: x, xb
    real(default) :: probability

    write (u, "(A)")  "* Test output: epa_handler_1"
    write (u, "(A)")  "*   Purpose: apply beam handler trivially (no-op)"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment"
    write (u, "(A)")

    call syntax_model_file_init ()
    call os_data%init ()
    call model_list%read_model &
         (var_str ("SM"), var_str ("SM.mdl"), &
         os_data, model)

    write (u, "(A)")  "* Initialize particle set"
    write (u, "(A)")

    call pset%init_direct &
         (n_beam = 2, n_in = 2, n_rem = 2, n_vir = 0, n_out = 2, &
         pdg = [11, -11, 22, 22, 11, -11, 13, -13], &
         model = model)

    sqrts = 100._default
    x = [0.6_default, 0.9_default]
    xb= 1 - x

    p(1) = vector4_moving (sqrts/2, sqrts/2, 3)
    p(2) = vector4_moving (sqrts/2,-sqrts/2, 3)
    p(3:4) = x * p(1:2)
    p(5:6) = xb  * p(1:2)
    p(7:8) = p(3:4)

    call pset%set_momentum (p, on_shell = .false.)

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call pacify (evt_trivial%particle_set)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize EPA handler transform"
    write (u, "(A)")

    evt_trivial%next => evt_isr_epa
    evt_isr_epa%previous => evt_trivial

    call evt_isr_epa%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Fill EPA handler transform"
    write (u, "(A)")

    call evt_isr_epa%prepare_new_event (1, 1)
    call evt_isr_epa%generate_weighted (probability)
    call evt_isr_epa%make_particle_set (0, .false.)

    call evt_isr_epa%write (u)

    write (u, "(A)")
    write (u, "(A,1x,F8.5)")  "Event probability =", probability

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_isr_epa%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: epa_handler_1"

  end subroutine epa_handler_1

  subroutine epa_handler_2 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(evt_trivial_t), target :: evt_trivial
    type(evt_isr_epa_t), target :: evt_isr_epa
    type(vector4_t), dimension(8) :: p
    real(default) :: sqrts
    real(default), dimension(2) :: x, xb
    class(rng_t), allocatable :: rng
    real(default) :: probability

    write (u, "(A)")  "* Test output: epa_handler_2"
    write (u, "(A)")  "*   Purpose: apply beam handler with two-beam recoil"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment"
    write (u, "(A)")

    call syntax_model_file_init ()
    call os_data%init ()
    call model_list%read_model &
         (var_str ("SM"), var_str ("SM.mdl"), &
         os_data, model)

    write (u, "(A)")  "* Initialize particle set"
    write (u, "(A)")

    call pset%init_direct (n_beam = 2, n_in = 2, n_rem = 2, n_vir = 0, n_out = 2, &
         pdg = [11, -11, 22, 22, 11, -11, 13, -13], model = model)

    sqrts = 100._default
    x = [0.6_default, 0.9_default]
    xb= 1 - x

    p(1) = vector4_moving (sqrts/2, sqrts/2, 3)
    p(2) = vector4_moving (sqrts/2,-sqrts/2, 3)
    p(3:4) = x * p(1:2)
    p(5:6) = xb  * p(1:2)
    p(7:8) = p(3:4)

    call pset%set_momentum (p, on_shell = .false.)

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call pacify (evt_trivial%particle_set)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize EPA handler transform"
    write (u, "(A)")

    evt_trivial%next => evt_isr_epa
    evt_isr_epa%previous => evt_trivial

    call evt_isr_epa%set_mode_string (var_str ("recoil"))
    call evt_isr_epa%set_data_epa ( &
         sqrts = sqrts, &
         q_max = sqrts, &
         m = 511.e-3_default &
         )

    allocate (rng_test_t :: rng)
    call rng%init (3)  ! default would produce pi for azimuthal angle
    call evt_isr_epa%import_rng (rng)

    call evt_isr_epa%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Fill EPA handler transform"
    write (u, "(A)")

    call evt_isr_epa%prepare_new_event (1, 1)
    call evt_isr_epa%generate_weighted (probability)
    call evt_isr_epa%make_particle_set (0, .false.)

    call evt_isr_epa%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A,1x,F8.5)")  "Event probability =", probability

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_isr_epa%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: epa_handler_2"

  end subroutine epa_handler_2

  subroutine epa_handler_3 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(evt_trivial_t), target :: evt_trivial
    type(evt_isr_epa_t), target :: evt_isr_epa
    type(vector4_t), dimension(8) :: p
    real(default) :: sqrts
    real(default), dimension(2) :: x0
    real(default), dimension(2) :: x, xb
    class(rng_t), allocatable :: rng
    real(default) :: probability

    write (u, "(A)")  "* Test output: epa_handler_3"
    write (u, "(A)")  "*   Purpose: apply beam handler for boosted beams &
         &and two-beam recoil"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment"
    write (u, "(A)")

    call syntax_model_file_init ()
    call os_data%init ()
    call model_list%read_model &
         (var_str ("SM"), var_str ("SM.mdl"), &
         os_data, model)

    write (u, "(A)")  "* Initialize particle set"
    write (u, "(A)")

    call pset%init_direct (n_beam = 2, n_in = 2, n_rem = 2, n_vir = 0, n_out = 2, &
         pdg = [11, -11, 22, 22, 11, -11, 13, -13], model = model)

    write (u, "(A)")  "* Event data"
    write (u, "(A)")

    sqrts = 100._default
    write (u, "(A,2(1x,F12.7))")  "sqrts   =", sqrts

    x0 = [0.9_default, 0.4_default]
    write (u, "(A,2(1x,F12.7))")  "x0      =", x0

    write (u, "(A)")
    write (u, "(A,2(1x,F12.7))")  "sqs_hat =", sqrts * sqrt (product (x0))

    x = [0.6_default, 0.9_default]
    xb= 1 - x
    write (u, "(A,2(1x,F12.7))")  "x       =", x

    write (u, "(A)")
    write (u, "(A,2(1x,F12.7))")  "x0 * x  =", x0 * x

    p(1) = x0(1) * vector4_moving (sqrts/2, sqrts/2, 3)
    p(2) = x0(2) * vector4_moving (sqrts/2,-sqrts/2, 3)
    p(3:4) = x * p(1:2)
    p(5:6) = xb  * p(1:2)
    p(7:8) = p(3:4)

    call pset%set_momentum (p, on_shell = .false.)

    write (u, "(A)")
    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call pacify (evt_trivial%particle_set)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize EPA handler transform"
    write (u, "(A)")

    evt_trivial%next => evt_isr_epa
    evt_isr_epa%previous => evt_trivial

    call evt_isr_epa%set_mode_string (var_str ("recoil"))
    call evt_isr_epa%set_data_epa ( &
         sqrts = sqrts, &
         q_max = sqrts, &
         m = 511.e-3_default &
         )

    allocate (rng_test_t :: rng)
    call rng%init (3)  ! default would produce pi for azimuthal angle
    call evt_isr_epa%import_rng (rng)

    call evt_isr_epa%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A)")  "* Fill EPA handler transform"
    write (u, "(A)")

    call evt_isr_epa%prepare_new_event (1, 1)
    call evt_isr_epa%generate_weighted (probability)
    call evt_isr_epa%make_particle_set (0, .false.)

    call evt_isr_epa%write (u, testflag=.true.)

    write (u, "(A)")
    write (u, "(A,1x,F8.5)")  "Event probability =", probability

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_isr_epa%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: epa_handler_3"

  end subroutine epa_handler_3


end module isr_epa_handler_uti

