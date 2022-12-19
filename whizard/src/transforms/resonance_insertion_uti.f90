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

module resonance_insertion_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use format_utils, only: write_separator
  use os_interface
  use lorentz
  use rng_base, only: rng_t
  use flavors, only: flavor_t
  use colors, only: color_t
  use models, only: syntax_model_file_init, syntax_model_file_final
  use models, only: model_list_t, model_t
  use particles, only: particle_t, particle_set_t

  use resonances, only: resonance_info_t
  use resonances, only: resonance_history_t
  use resonances, only: resonance_history_set_t

  use event_transforms
  use resonance_insertion

  use rng_base_ut, only: rng_test_t

  implicit none
  private

  public :: resonance_insertion_1
  public :: resonance_insertion_2
  public :: resonance_insertion_3
  public :: resonance_insertion_4
  public :: resonance_insertion_5
  public :: resonance_insertion_6

contains

  subroutine resonance_insertion_1 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(evt_trivial_t), target :: evt_trivial
    type(evt_resonance_t), target :: evt_resonance
    type(flavor_t) :: fw
    type(color_t) :: col
    real(default) :: mw, ew, pw
    type(vector4_t), dimension(5) :: p
    class(rng_t), allocatable :: rng
    real(default) :: probability
    integer, dimension(:), allocatable :: i_invalid
    type(particle_t), dimension(:), allocatable :: prt_invalid
    integer :: i

    write (u, "(A)")  "* Test output: resonance_insertion_1"
    write (u, "(A)")  "*   Purpose: apply simple resonance insertion"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize environment"
    write (u, "(A)")

    call syntax_model_file_init ()
    call os_data%init ()
    call model_list%read_model &
         (var_str ("SM"), var_str ("SM.mdl"), &
         os_data, model)
    ! reset slightly in order to avoid a rounding ambiguity
    call model%set_real (var_str ("mW"), 80.418_default)

    write (u, "(A)")  "* Initialize particle set"
    write (u, "(A)")

    call pset%init_direct (n_beam = 0, n_in = 2, n_rem = 0, n_vir = 0, n_out = 3, &
         pdg = [1, -1, 1, -2, 24], model = model)

    call fw%init (24, model)

    mw = fw%get_mass ()
    ew = 200._default
    pw = sqrt (ew**2 - mw**2)

    p(1) = vector4_moving (ew, ew, 3)
    p(2) = vector4_moving (ew,-ew, 3)
    p(3) = vector4_moving (ew/2, vector3_moving ([pw/2, mw/2, 0._default]))
    p(4) = vector4_moving (ew/2, vector3_moving ([pw/2,-mw/2, 0._default]))
    p(5) = vector4_moving (ew, vector3_moving ([-pw, 0._default, 0._default]))

    call pset%set_momentum (p, on_shell = .true.)

    call col%init_col_acl (1,0)
    call pset%set_color (1, col)

    call col%init_col_acl (0,1)
    call pset%set_color (2, col)

    call col%init_col_acl (2,0)
    call pset%set_color (3, col)

    call col%init_col_acl (0,2)
    call pset%set_color (4, col)

    call col%init_col_acl (0,0)
    call pset%set_color (5, col)

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Prepare resonance history set"
    write (u, "(A)")

    call res_history_set(1)%init ()

    call res_info%init (3, -24, model, 2)
    call res_history%add_resonance (res_info)
    call res_history_set(1)%enter (res_history)
    call res_history%clear ()

    call res_history_set(1)%freeze ()

    write (u, "(A)")  "* Initialize resonance insertion transform"
    write (u, "(A)")

    evt_trivial%next => evt_resonance
    evt_resonance%previous => evt_trivial

    allocate (rng_test_t :: rng)
    call evt_resonance%import_rng (rng)

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)
    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Fill resonance insertion transform"
    write (u, "(A)")

    call evt_resonance%prepare_new_event (1, 1)
    call evt_resonance%init_selector ([1._default])
    call evt_resonance%generate_weighted (probability)
    call evt_resonance%make_particle_set (0, .false.)

    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A,1x,F8.5)")  "Event probability =", probability

    write (u, "(A)")
    call evt_resonance%find_prt_invalid_color (i_invalid, prt_invalid)
    write (u, "(A)")  "Particles with invalid color:"
    select case (size (prt_invalid))
    case (0)
       write (u, "(2x,A)")  "[none]"
    case default
       do i = 1, size (prt_invalid)
          write (u, "(1x,A,1x,I0)", advance="no")  "Particle", i_invalid(i)
          call prt_invalid(i)%write (u)
       end do
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_resonance%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonance_insertion_1"

  end subroutine resonance_insertion_1

  subroutine resonance_insertion_2 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(evt_trivial_t), target :: evt_trivial
    type(evt_resonance_t), target :: evt_resonance
    type(color_t) :: col
    class(rng_t), allocatable :: rng
    real(default) :: probability
    type(particle_t), dimension(:), allocatable :: prt_invalid
    integer, dimension(:), allocatable :: i_invalid
    integer :: i

    write (u, "(A)")  "* Test output: resonance_insertion_2"
    write (u, "(A)")  "*   Purpose: resonance insertion with color mismatch"
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

    call pset%init_direct (n_beam = 0, n_in = 2, n_rem = 0, n_vir = 0, n_out = 3, &
         pdg = [1, -1, 1, -2, 24], model = model)

    call col%init_col_acl (1,0)
    call pset%set_color (1, col)

    call col%init_col_acl (0,2)
    call pset%set_color (2, col)

    call col%init_col_acl (1,0)
    call pset%set_color (3, col)

    call col%init_col_acl (0,2)
    call pset%set_color (4, col)

    call col%init_col_acl (0,0)
    call pset%set_color (5, col)

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Prepare resonance history set"
    write (u, "(A)")

    call res_history_set(1)%init ()

    call res_info%init (3, -24, model, 2)
    call res_history%add_resonance (res_info)
    call res_history_set(1)%enter (res_history)
    call res_history%clear ()

    call res_history_set(1)%freeze ()

    write (u, "(A)")  "* Initialize resonance insertion transform"
    write (u, "(A)")

    evt_trivial%next => evt_resonance
    evt_resonance%previous => evt_trivial

    allocate (rng_test_t :: rng)
    call evt_resonance%import_rng (rng)

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)
    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Fill resonance insertion transform"
    write (u, "(A)")

    call evt_resonance%prepare_new_event (1, 1)
    call evt_resonance%init_selector ([1._default])
    call evt_resonance%generate_weighted (probability)
    call evt_resonance%make_particle_set (0, .false.)

    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A,1x,F8.5)")  "Event probability =", probability

    write (u, "(A)")
    call evt_resonance%find_prt_invalid_color (i_invalid, prt_invalid)
    write (u, "(A)")  "Particles with invalid color:"
    select case (size (prt_invalid))
    case (0)
       write (u, "(2x,A)")  "[none]"
    case default
       do i = 1, size (prt_invalid)
          write (u, "(1x,A,1x,I0)", advance="no")  "Particle", i_invalid(i)
          call prt_invalid(i)%write (u)
       end do
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_resonance%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonance_insertion_2"

  end subroutine resonance_insertion_2

  subroutine resonance_insertion_3 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(evt_trivial_t), target :: evt_trivial
    type(evt_resonance_t), target :: evt_resonance
    type(color_t) :: col
    class(rng_t), allocatable :: rng
    real(default) :: probability
    type(particle_t), dimension(:), allocatable :: prt_invalid
    integer, dimension(:), allocatable :: i_invalid
    integer :: i

    write (u, "(A)")  "* Test output: resonance_insertion_3"
    write (u, "(A)")  "*   Purpose: resonance insertion with color mismatch"
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

    call pset%init_direct (n_beam = 0, n_in = 2, n_rem = 0, n_vir = 0, n_out = 6, &
         pdg = [2, -2, 24, 5, 5, -5, -24, -5], model = model)

    call col%init_col_acl (1,0)
    call pset%set_color (1, col)

    call col%init_col_acl (0,2)
    call pset%set_color (2, col)

    call col%init_col_acl (0,0)
    call pset%set_color (3, col)

    call col%init_col_acl (1,0)
    call pset%set_color (4, col)

    call col%init_col_acl (3,0)
    call pset%set_color (5, col)

    call col%init_col_acl (0,3)
    call pset%set_color (6, col)

    call col%init_col_acl (0,0)
    call pset%set_color (7, col)

    call col%init_col_acl (0,2)
    call pset%set_color (8, col)

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Prepare resonance history set"
    write (u, "(A)")

    call res_history_set(1)%init ()

    call res_info%init (3, 6, model, 6)
    call res_history%add_resonance (res_info)
    call res_info%init (12, 25, model, 6)
    call res_history%add_resonance (res_info)
    call res_info%init (60, -6, model, 6)
    call res_history%add_resonance (res_info)
    call res_history_set(1)%enter (res_history)
    call res_history%clear ()

    call res_history_set(1)%freeze ()

    write (u, "(A)")  "* Initialize resonance insertion transform"
    write (u, "(A)")

    evt_trivial%next => evt_resonance
    evt_resonance%previous => evt_trivial

    allocate (rng_test_t :: rng)
    call evt_resonance%import_rng (rng)

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)
    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Fill resonance insertion transform"
    write (u, "(A)")

    call evt_resonance%prepare_new_event (1, 1)
    call evt_resonance%init_selector ([1._default])
    call evt_resonance%generate_weighted (probability)
    call evt_resonance%make_particle_set (0, .false.)

    call evt_resonance%write (u)

    write (u, "(A)")
    call evt_resonance%find_prt_invalid_color (i_invalid, prt_invalid)
    write (u, "(A)")  "Particles with invalid color:"
    select case (size (prt_invalid))
    case (0)
       write (u, "(2x,A)")  "[none]"
    case default
       do i = 1, size (prt_invalid)
          write (u, "(1x,A,1x,I0)", advance="no")  "Particle", i_invalid(i)
          call prt_invalid(i)%write (u)
       end do
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_resonance%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonance_insertion_3"

  end subroutine resonance_insertion_3

  subroutine resonance_insertion_4 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(evt_trivial_t), target :: evt_trivial
    type(evt_resonance_t), target :: evt_resonance
    type(color_t) :: col
    class(rng_t), allocatable :: rng
    real(default) :: probability
    integer :: i

    write (u, "(A)")  "* Test output: resonance_insertion_4"
    write (u, "(A)")  "*   Purpose: resonance history selection"
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

    call pset%init_direct (n_beam = 0, n_in = 2, n_rem = 0, n_vir = 0, n_out = 4, &
         pdg = [1, -1, 1, -2, -3, 4], model = model)

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Prepare resonance history set"
    write (u, "(A)")

    call res_history_set(1)%init ()

    call res_info%init (3, -24, model, 4)
    call res_history%add_resonance (res_info)
    call res_history_set(1)%enter (res_history)
    call res_history%clear ()

    call res_info%init (12, 24, model, 4)
    call res_history%add_resonance (res_info)
    call res_history_set(1)%enter (res_history)
    call res_history%clear ()

    call res_info%init (12, 24, model, 4)
    call res_history%add_resonance (res_info)
    call res_info%init (15, 25, model, 4)
    call res_history%add_resonance (res_info)
    call res_history_set(1)%enter (res_history)
    call res_history%clear ()

    call res_history_set(1)%freeze ()

    write (u, "(A)")  "* Initialize resonance insertion transform"
    write (u, "(A)")

    evt_trivial%next => evt_resonance
    evt_resonance%previous => evt_trivial

    allocate (rng_test_t :: rng)
    call evt_resonance%import_rng (rng)

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)
    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Fill resonance insertion transform"
    write (u, "(A)")

    do i = 1, 6
       write (u, "(A,1x,I0)")  "* Event #", i
       write (u, "(A)")

       call evt_resonance%prepare_new_event (1, 1)
       call evt_resonance%init_selector ([1._default, 2._default, 1._default])
       call evt_resonance%generate_weighted (probability)
       call evt_resonance%make_particle_set (0, .false.)

       call evt_resonance%write (u)
       write (u, "(A)")
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_resonance%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonance_insertion_4"

  end subroutine resonance_insertion_4

  subroutine resonance_insertion_5 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(particle_set_t) :: pset
    type(model_list_t) :: model_list
    type(model_t), pointer :: model
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(evt_trivial_t), target :: evt_trivial
    type(evt_resonance_t), target :: evt_resonance
    type(color_t) :: col
    class(rng_t), allocatable :: rng
    real(default) :: probability
    integer :: i

    write (u, "(A)")  "* Test output: resonance_insertion_5"
    write (u, "(A)")  "*   Purpose: resonance history selection including no resonance"
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

    call pset%init_direct (n_beam = 0, n_in = 2, n_rem = 0, n_vir = 0, n_out = 4, &
         pdg = [1, -1, 1, -2, -3, 4], model = model)

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)

    write (u, "(A)")  "* Prepare resonance history set"
    write (u, "(A)")

    call res_history_set(1)%init ()

    call res_info%init (3, -24, model, 4)
    call res_history%add_resonance (res_info)
    call res_history_set(1)%enter (res_history)
    call res_history%clear ()

    call res_history_set(1)%freeze ()

    write (u, "(A)")  "* Initialize resonance insertion transform"
    write (u, "(A)")

    evt_trivial%next => evt_resonance
    evt_resonance%previous => evt_trivial

    allocate (rng_test_t :: rng)
    call evt_resonance%import_rng (rng)

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)

    write (u, "(A)")  "* Fill resonance insertion transform"
    write (u, "(A)")

    do i = 1, 2
       write (u, "(A,1x,I0)")  "* Event #", i
       write (u, "(A)")

       call evt_resonance%prepare_new_event (1, 1)
       call evt_resonance%init_selector ([1._default, 3._default], offset = -1)
       call evt_resonance%generate_weighted (probability)
       call evt_resonance%make_particle_set (0, .false.)

       call evt_resonance%write (u)
       write (u, "(A)")
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_resonance%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonance_insertion_5"

  end subroutine resonance_insertion_5

  subroutine resonance_insertion_6 (u)
    integer, intent(in) :: u
    type(os_data_t) :: os_data
    type(model_list_t) :: model_list
    type(particle_set_t) :: pset
    type(model_t), pointer :: model
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_history_set_t), dimension(1) :: res_history_set
    type(evt_trivial_t), target :: evt_trivial
    type(evt_resonance_t), target :: evt_resonance
    class(rng_t), allocatable :: rng
    real(default) :: probability

    write (u, "(A)")  "* Test output: resonance_insertion_6"
    write (u, "(A)")  "*   Purpose: resonance insertion with structured beams"
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

    write (u, "(A)")  "* Fill trivial event transform"
    write (u, "(A)")

    call evt_trivial%reset ()
    call evt_trivial%set_particle_set (pset, 1, 1)
    call evt_trivial%write (u)
    call write_separator (u, 2)

    write (u, "(A)")
    write (u, "(A)")  "* Prepare resonance history set"
    write (u, "(A)")

    call res_history_set(1)%init ()

    call res_info%init (3, 23, model, 2)
    call res_history%add_resonance (res_info)
    call res_history_set(1)%enter (res_history)
    call res_history%clear ()

    call res_history_set(1)%freeze ()

    write (u, "(A)")  "* Initialize resonance insertion transform"
    write (u, "(A)")

    evt_trivial%next => evt_resonance
    evt_resonance%previous => evt_trivial

    allocate (rng_test_t :: rng)
    call evt_resonance%import_rng (rng)

    call evt_resonance%set_resonance_data (res_history_set)
    call evt_resonance%select_component (1)
    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Fill resonance insertion transform"
    write (u, "(A)")

    call evt_resonance%prepare_new_event (1, 1)
    call evt_resonance%init_selector ([1._default])
    call evt_resonance%generate_weighted (probability)
    call evt_resonance%make_particle_set (0, .false.)

    call evt_resonance%write (u)

    write (u, "(A)")
    write (u, "(A,1x,F8.5)")  "Event probability =", probability

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call evt_resonance%final ()
    call evt_trivial%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonance_insertion_6"

  end subroutine resonance_insertion_6


end module resonance_insertion_uti

