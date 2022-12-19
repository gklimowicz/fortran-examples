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

module particles_uti

  use kinds, only: default
  use io_units
  use numeric_utils
  use constants, only: one, tiny_07
  use lorentz
  use flavors
  use colors
  use helicities
  use quantum_numbers
  use state_matrices
  use interactions
  use evaluators
  use model_data
  use subevents

  use particles

  implicit none
  private

  public :: particles_1
  public :: particles_2
  public :: particles_3
  public :: particles_4
  public :: particles_5
  public :: particles_6
  public :: particles_7
  public :: particles_8
  public :: particles_9

contains

  subroutine particles_1 (u)
    use os_interface
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t), dimension(3) :: flv
    type(color_t), dimension(3) :: col
    type(helicity_t), dimension(3) :: hel
    type(quantum_numbers_t), dimension(3) :: qn
    type(vector4_t), dimension(3) :: p
    type(interaction_t), target :: int1, int2
    type(quantum_numbers_mask_t) :: qn_mask_conn
    type(evaluator_t), target :: eval
    type(interaction_t) :: int
    type(particle_set_t) :: particle_set1, particle_set2
    type(particle_set_t) :: particle_set3, particle_set4
    type(subevt_t) :: subevt
    logical :: ok
    integer :: unit, iostat

    write (u, "(A)")  "* Test output: Particles"
    write (u, "(A)")  "*   Purpose: test particle_set routines"
    write (u, "(A)")

    write (u, "(A)")  "* Reading model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Initializing production process"

    call int1%basic_init (2, 0, 1, set_relations=.true.)
    call flv%init ([1, -1, 23], model)
    call col%init_col_acl ([0, 0, 0], [0, 0, 0])
    call hel(3)%init (1, 1)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0.25_default, 0._default))
    call hel(3)%init (1,-1)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0._default, 0.25_default))
    call hel(3)%init (-1, 1)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0._default,-0.25_default))
    call hel(3)%init (-1,-1)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0.25_default, 0._default))
    call hel(3)%init (0, 0)
    call qn%init (flv, col, hel)
    call int1%add_state (qn, value=(0.5_default, 0._default))
    call int1%freeze ()
    p(1) = vector4_moving (45._default, 45._default, 3)
    p(2) = vector4_moving (45._default,-45._default, 3)
    p(3) = p(1) + p(2)
    call int1%set_momenta (p)

    write (u, "(A)")
    write (u, "(A)")  "* Setup decay process"

    call int2%basic_init (1, 0, 2, set_relations=.true.)
    call flv%init ([23, 1, -1], model)
    call col%init_col_acl ([0, 501, 0], [0, 0, 501])
    call hel%init ([1, 1, 1], [1, 1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(1._default, 0._default))
    call hel%init ([1, 1, 1], [-1,-1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(0._default, 0.1_default))
    call hel%init ([-1,-1,-1], [1, 1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(0._default,-0.1_default))
    call hel%init ([-1,-1,-1], [-1,-1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(1._default, 0._default))
    call hel%init ([0, 1,-1], [0, 1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(4._default, 0._default))
    call hel%init ([0,-1, 1], [0, 1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(2._default, 0._default))
    call hel%init ([0, 1,-1], [0,-1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(2._default, 0._default))
    call hel%init ([0,-1, 1], [0,-1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(4._default, 0._default))
    call flv%init ([23, 2, -2], model)
    call hel%init ([0, 1,-1], [0, 1,-1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(0.5_default, 0._default))
    call hel%init ([0,-1, 1], [0,-1, 1])
    call qn%init (flv, col, hel)
    call int2%add_state (qn, value=(0.5_default, 0._default))
    call int2%freeze ()
    p(2) = vector4_moving (45._default, 45._default, 2)
    p(3) = vector4_moving (45._default,-45._default, 2)
    call int2%set_momenta (p)
    call int2%set_source_link (1, int1, 3)
    call int1%basic_write (u)
    call int2%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Concatenate production and decay"

    call eval%init_product (int1, int2, qn_mask_conn, &
         connections_are_resonant=.true.)
    call eval%receive_momenta ()
    call eval%evaluate ()
    call eval%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Factorize as subevent (complete, polarized)"
    write (u, "(A)")

    int = eval%interaction_t
    call particle_set1%init &
         (ok, int, int, FM_FACTOR_HELICITY, &
          [0.2_default, 0.2_default], .false., .true.)
    call particle_set1%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Factorize as subevent (in/out only, selected helicity)"
    write (u, "(A)")

    int = eval%interaction_t
    call particle_set2%init &
         (ok, int, int, FM_SELECT_HELICITY, &
          [0.9_default, 0.9_default], .false., .false.)
    call particle_set2%write (u)
    call particle_set2%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Factorize as subevent (complete, selected helicity)"
    write (u, "(A)")

    int = eval%interaction_t
    call particle_set2%init &
         (ok, int, int, FM_SELECT_HELICITY, &
          [0.7_default, 0.7_default], .false., .true.)
    call particle_set2%write (u)

    write (u, "(A)")
    write (u, "(A)")  &
         "* Factorize (complete, polarized, correlated); write and read again"
    write (u, "(A)")

    int = eval%interaction_t
    call particle_set3%init &
         (ok, int, int, FM_FACTOR_HELICITY, &
          [0.7_default, 0.7_default], .true., .true.)
    call particle_set3%write (u)

    unit = free_unit ()
    open (unit, action="readwrite", form="unformatted", status="scratch")
    call particle_set3%write_raw (unit)
    rewind (unit)
    call particle_set4%read_raw (unit, iostat=iostat)
    call particle_set4%set_model (model)
    close (unit)

    write (u, "(A)")
    write (u, "(A)")  "* Result from reading"
    write (u, "(A)")

    call particle_set4%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Transform to a subevt object"
    write (u, "(A)")

    call particle_set4%to_subevt (subevt)
    call subevt%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call particle_set1%final ()
    call particle_set2%final ()
    call particle_set3%final ()
    call particle_set4%final ()
    call eval%final ()
    call int1%final ()
    call int2%final ()

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particles_1"

  end subroutine particles_1

  subroutine particles_2 (u)

    integer, intent(in) :: u
    type(interaction_t) :: int
    type(state_flv_content_t) :: state_flv
    type(particle_set_t) :: pset
    type(flavor_t), dimension(:), allocatable :: flv
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer :: i, j

    write (u, "(A)")  "* Test output: Particles"
    write (u, "(A)")  "*   Purpose: reconstruct simple interaction"
    write (u, "(A)")

    write (u, "(A)")  "* Set up a 2 -> 3 interaction"
    write (u, "(A)")  "    + incoming partons marked as virtual"
    write (u, "(A)")  "    + no quantum numbers"
    write (u, "(A)")

    call reset_interaction_counter ()
    call int%basic_init (0, 2, 3)
    do i = 1, 2
       do j = 3, 5
          call int%relate (i, j)
       end do
    end do

    allocate (qn (5))
    call int%add_state (qn)
    call int%freeze ()

    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually set up a flavor-content record"
    write (u, "(A)")

    call state_flv%init (1, &
         mask = [.false., .false., .true., .true., .true.])
    call state_flv%set_entry (1, &
         pdg = [11, 12, 3, 4, 5], &
         map = [1, 2, 3, 4, 5])

    call state_flv%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually create a matching particle set"
    write (u, "(A)")

    pset%n_beam = 0
    pset%n_in   = 2
    pset%n_vir  = 0
    pset%n_out  = 3
    pset%n_tot  = 5
    allocate (pset%prt (pset%n_tot))
    do i = 1, 2
       call pset%prt(i)%reset_status (PRT_INCOMING)
       call pset%prt(i)%set_children ([3,4,5])
    end do
    do i = 3, 5
       call pset%prt(i)%reset_status (PRT_OUTGOING)
       call pset%prt(i)%set_parents ([1,2])
    end do
    call pset%prt(1)%set_momentum (vector4_at_rest (1._default))
    call pset%prt(2)%set_momentum (vector4_at_rest (2._default))
    call pset%prt(3)%set_momentum (vector4_at_rest (5._default))
    call pset%prt(4)%set_momentum (vector4_at_rest (4._default))
    call pset%prt(5)%set_momentum (vector4_at_rest (3._default))

    allocate (flv (5))
    call flv%init ([11,12,5,4,3])
    do i = 1, 5
       call pset%prt(i)%set_flavor (flv(i))
    end do

    call pset%write (u)

    write (u, "(A)")
    write (u, "(A)")  "*   Fill interaction from particle set"
    write (u, "(A)")

    call pset%fill_interaction (int, 2, state_flv=state_flv)
    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call int%final ()
    call pset%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particles_2"

  end subroutine particles_2

  subroutine particles_3 (u)

    integer, intent(in) :: u
    type(interaction_t) :: int
    type(state_flv_content_t) :: state_flv
    type(particle_set_t) :: pset
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer :: i, j

    write (u, "(A)")  "* Test output: Particles"
    write (u, "(A)")  "*   Purpose: reconstruct simple interaction"
    write (u, "(A)")

    write (u, "(A)")  "* Set up a 2 -> 2 -> 3 interaction with radiation"
    write (u, "(A)")  "    + no quantum numbers"
    write (u, "(A)")

    call reset_interaction_counter ()
    call int%basic_init (0, 6, 3)
    call int%relate (1, 3)
    call int%relate (1, 4)
    call int%relate (2, 5)
    call int%relate (2, 6)
    do i = 4, 6, 2
       do j = 7, 9
          call int%relate (i, j)
       end do
    end do

    allocate (qn (9))
    call int%add_state (qn)
    call int%freeze ()

    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually set up a flavor-content record"
    write (u, "(A)")

    call state_flv%init (1, &
         mask = [.false., .false., .false., .false., .false., .false., &
         .true., .true., .true.])
    call state_flv%set_entry (1, &
         pdg = [2011, 2012, 91, 11, 92, 12, 3, 4, 5], &
         map = [1, 2, 3, 4, 5, 6, 7, 8, 9])

    call state_flv%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually create a matching particle set"
    write (u, "(A)")

    call create_test_particle_set_1 (pset)

    call pset%write (u)

    write (u, "(A)")
    write (u, "(A)")  "*   Fill interaction from particle set"
    write (u, "(A)")

    call pset%fill_interaction (int, 2, state_flv=state_flv)
    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call int%final ()
    call pset%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particles_3"

  end subroutine particles_3

  subroutine particles_4 (u)

    integer, intent(in) :: u
    type(interaction_t) :: int
    type(interaction_t), target :: int_beams
    type(state_flv_content_t) :: state_flv
    type(particle_set_t) :: pset
    type(flavor_t), dimension(:), allocatable :: flv
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer :: i, j

    write (u, "(A)")  "* Test output: Particles"
    write (u, "(A)")  "*   Purpose: reconstruct beams"
    write (u, "(A)")

    call reset_interaction_counter ()

    write (u, "(A)")  "* Set up an interaction that contains beams only"
    write (u, "(A)")

    call int_beams%basic_init (0, 0, 2)
    call int_beams%set_momentum (vector4_at_rest (1._default), 1)
    call int_beams%set_momentum (vector4_at_rest (2._default), 2)
    allocate (qn (2))
    call int_beams%add_state (qn)
    call int_beams%freeze ()

    call int_beams%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Set up a 2 -> 2 -> 3 interaction with radiation"
    write (u, "(A)")  "    + no quantum numbers"
    write (u, "(A)")

    call int%basic_init (0, 6, 3)
    call int%relate (1, 3)
    call int%relate (1, 4)
    call int%relate (2, 5)
    call int%relate (2, 6)
    do i = 4, 6, 2
       do j = 7, 9
          call int%relate (i, j)
       end do
    end do
    do i = 1, 2
       call int%set_source_link (i, int_beams, i)
    end do

    deallocate (qn)
    allocate (qn (9))
    call int%add_state (qn)
    call int%freeze ()

    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually set up a flavor-content record"
    write (u, "(A)")

    call state_flv%init (1, &
         mask = [.false., .false., .false., .false., .false., .false., &
         .true., .true., .true.])
    call state_flv%set_entry (1, &
         pdg = [2011, 2012, 91, 11, 92, 12, 3, 4, 5], &
         map = [1, 2, 3, 4, 5, 6, 7, 8, 9])

    call state_flv%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually create a matching particle set"
    write (u, "(A)")

    pset%n_beam = 0
    pset%n_in   = 2
    pset%n_vir  = 0
    pset%n_out  = 3
    pset%n_tot  = 5

    allocate (pset%prt (pset%n_tot))
    call pset%prt(1)%reset_status (PRT_INCOMING)
    call pset%prt(2)%reset_status (PRT_INCOMING)
    call pset%prt(3)%reset_status (PRT_OUTGOING)
    call pset%prt(4)%reset_status (PRT_OUTGOING)
    call pset%prt(5)%reset_status (PRT_OUTGOING)

    call pset%prt(1)%set_children ([3,4,5])
    call pset%prt(2)%set_children ([3,4,5])

    call pset%prt(3)%set_parents ([1,2])
    call pset%prt(4)%set_parents ([1,2])
    call pset%prt(5)%set_parents ([1,2])

    call pset%prt(1)%set_momentum (vector4_at_rest (6._default))
    call pset%prt(2)%set_momentum (vector4_at_rest (6._default))
    call pset%prt(3)%set_momentum (vector4_at_rest (3._default))
    call pset%prt(4)%set_momentum (vector4_at_rest (4._default))
    call pset%prt(5)%set_momentum (vector4_at_rest (5._default))

    allocate (flv (5))
    call flv%init ([11, 12, 3, 4, 5])
    do i = 1, 5
       call pset%prt(i)%set_flavor (flv(i))
    end do

    call pset%write (u)

    write (u, "(A)")
    write (u, "(A)")  "*   Fill interaction from particle set"
    write (u, "(A)")

    call pset%fill_interaction (int, 2, state_flv=state_flv, &
         recover_beams = .true.)
    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call int%final ()
    call pset%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particles_4"

  end subroutine particles_4

  subroutine particles_5 (u)

    integer, intent(in) :: u
    type(interaction_t) :: int
    type(state_flv_content_t) :: state_flv
    type(particle_set_t) :: pset
    type(flavor_t), dimension(:), allocatable :: flv
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer :: i, j

    write (u, "(A)")  "* Test output: Particles"
    write (u, "(A)")  "*   Purpose: reconstruct event with duplicate entries"
    write (u, "(A)")

    write (u, "(A)")  "* Set up a 2 -> 2 -> 3 interaction with radiation"
    write (u, "(A)")  "    + no quantum numbers"
    write (u, "(A)")

    call reset_interaction_counter ()
    call int%basic_init (0, 6, 3)
    call int%relate (1, 3)
    call int%relate (1, 4)
    call int%relate (2, 5)
    call int%relate (2, 6)
    do i = 4, 6, 2
       do j = 7, 9
          call int%relate (i, j)
       end do
    end do

    allocate (qn (9))
    call int%add_state (qn)
    call int%freeze ()

    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually set up a flavor-content record"
    write (u, "(A)")

    call state_flv%init (1, &
         mask = [.false., .false., .false., .false., .false., .false., &
         .true., .true., .true.])
    call state_flv%set_entry (1, &
         pdg = [2011, 2012, 91, 11, 92, 12, 3, 4, 5], &
         map = [1, 2, 3, 4, 5, 6, 7, 8, 9])

    call state_flv%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually create a matching particle set"
    write (u, "(A)")

    pset%n_beam = 2
    pset%n_in   = 2
    pset%n_vir  = 4
    pset%n_out  = 5
    pset%n_tot  = 13

    allocate (pset%prt (pset%n_tot))
    call pset%prt(1)%reset_status (PRT_BEAM)
    call pset%prt(2)%reset_status (PRT_BEAM)
    call pset%prt(3)%reset_status (PRT_VIRTUAL)
    call pset%prt(4)%reset_status (PRT_VIRTUAL)
    call pset%prt(5)%reset_status (PRT_VIRTUAL)
    call pset%prt(6)%reset_status (PRT_VIRTUAL)
    call pset%prt(7)%reset_status (PRT_INCOMING)
    call pset%prt(8)%reset_status (PRT_INCOMING)
    call pset%prt( 9)%reset_status (PRT_OUTGOING)
    call pset%prt(10)%reset_status (PRT_OUTGOING)
    call pset%prt(11)%reset_status (PRT_OUTGOING)
    call pset%prt(12)%reset_status (PRT_OUTGOING)
    call pset%prt(13)%reset_status (PRT_OUTGOING)

    call pset%prt(1)%set_children ([3,4])
    call pset%prt(2)%set_children ([5,6])
    call pset%prt(3)%set_children ([ 7])
    call pset%prt(4)%set_children ([ 9])
    call pset%prt(5)%set_children ([ 8])
    call pset%prt(6)%set_children ([10])
    call pset%prt(7)%set_children ([11,12,13])
    call pset%prt(8)%set_children ([11,12,13])

    call pset%prt(3)%set_parents ([1])
    call pset%prt(4)%set_parents ([1])
    call pset%prt(5)%set_parents ([2])
    call pset%prt(6)%set_parents ([2])
    call pset%prt( 7)%set_parents ([3])
    call pset%prt( 8)%set_parents ([5])
    call pset%prt( 9)%set_parents ([4])
    call pset%prt(10)%set_parents ([6])
    call pset%prt(11)%set_parents ([7,8])
    call pset%prt(12)%set_parents ([7,8])
    call pset%prt(13)%set_parents ([7,8])

    call pset%prt(1)%set_momentum (vector4_at_rest (1._default))
    call pset%prt(2)%set_momentum (vector4_at_rest (2._default))
    call pset%prt(3)%set_momentum (vector4_at_rest (4._default))
    call pset%prt(4)%set_momentum (vector4_at_rest (3._default))
    call pset%prt(5)%set_momentum (vector4_at_rest (6._default))
    call pset%prt(6)%set_momentum (vector4_at_rest (5._default))
    call pset%prt(7)%set_momentum (vector4_at_rest (4._default))
    call pset%prt(8)%set_momentum (vector4_at_rest (6._default))
    call pset%prt( 9)%set_momentum (vector4_at_rest (3._default))
    call pset%prt(10)%set_momentum (vector4_at_rest (5._default))
    call pset%prt(11)%set_momentum (vector4_at_rest (7._default))
    call pset%prt(12)%set_momentum (vector4_at_rest (8._default))
    call pset%prt(13)%set_momentum (vector4_at_rest (9._default))

    allocate (flv (13))
    call flv%init ([2011, 2012, 11, 91, 12, 92, 11, 12, 91, 92, 3, 4, 5])
    do i = 1, 13
       call pset%prt(i)%set_flavor (flv(i))
    end do

    call pset%write (u)

    write (u, "(A)")
    write (u, "(A)")  "*   Fill interaction from particle set"
    write (u, "(A)")

    call pset%fill_interaction (int, 2, state_flv=state_flv)
    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call int%final ()
    call pset%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particles_5"

  end subroutine particles_5

  subroutine particles_6 (u)

    integer, intent(in) :: u
    type(interaction_t) :: int
    type(state_flv_content_t) :: state_flv
    type(particle_set_t) :: pset
    type(flavor_t), dimension(:), allocatable :: flv
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer :: i, j

    write (u, "(A)")  "* Test output: Particles"
    write (u, "(A)")  "*   Purpose: reconstruct interaction with pair spectrum"
    write (u, "(A)")

    write (u, "(A)")  "* Set up a 2 -> 2 -> 3 interaction with radiation"
    write (u, "(A)")  "    + no quantum numbers"
    write (u, "(A)")

    call reset_interaction_counter ()
    call int%basic_init (0, 6, 3)
    do i = 1, 2
       do j = 3, 6
          call int%relate (i, j)
       end do
    end do
    do i = 5, 6
       do j = 7, 9
          call int%relate (i, j)
       end do
    end do

    allocate (qn (9))
    call int%add_state (qn)
    call int%freeze ()

    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually set up a flavor-content record"
    write (u, "(A)")

    call state_flv%init (1, &
         mask = [.false., .false., .false., .false., .false., .false., &
         .true., .true., .true.])
    call state_flv%set_entry (1, &
         pdg = [1011, 1012, 21, 22, 11, 12, 3, 4, 5], &
         map = [1, 2, 3, 4, 5, 6, 7, 8, 9])

    call state_flv%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually create a matching particle set"
    write (u, "(A)")

    pset%n_beam = 2
    pset%n_in   = 2
    pset%n_vir  = 2
    pset%n_out  = 3
    pset%n_tot  = 9

    allocate (pset%prt (pset%n_tot))
    call pset%prt(1)%reset_status (PRT_BEAM)
    call pset%prt(2)%reset_status (PRT_BEAM)
    call pset%prt(3)%reset_status (PRT_INCOMING)
    call pset%prt(4)%reset_status (PRT_INCOMING)
    call pset%prt(5)%reset_status (PRT_OUTGOING)
    call pset%prt(6)%reset_status (PRT_OUTGOING)
    call pset%prt(7)%reset_status (PRT_OUTGOING)
    call pset%prt(8)%reset_status (PRT_OUTGOING)
    call pset%prt(9)%reset_status (PRT_OUTGOING)

    call pset%prt(1)%set_children ([3,4,5,6])
    call pset%prt(2)%set_children ([3,4,5,6])
    call pset%prt(3)%set_children ([7,8,9])
    call pset%prt(4)%set_children ([7,8,9])

    call pset%prt(3)%set_parents ([1,2])
    call pset%prt(4)%set_parents ([1,2])
    call pset%prt(5)%set_parents ([1,2])
    call pset%prt(6)%set_parents ([1,2])
    call pset%prt(7)%set_parents ([3,4])
    call pset%prt(8)%set_parents ([3,4])
    call pset%prt(9)%set_parents ([3,4])

    call pset%prt(1)%set_momentum (vector4_at_rest (1._default))
    call pset%prt(2)%set_momentum (vector4_at_rest (2._default))
    call pset%prt(3)%set_momentum (vector4_at_rest (5._default))
    call pset%prt(4)%set_momentum (vector4_at_rest (6._default))
    call pset%prt(5)%set_momentum (vector4_at_rest (3._default))
    call pset%prt(6)%set_momentum (vector4_at_rest (4._default))
    call pset%prt(7)%set_momentum (vector4_at_rest (7._default))
    call pset%prt(8)%set_momentum (vector4_at_rest (8._default))
    call pset%prt(9)%set_momentum (vector4_at_rest (9._default))

    allocate (flv (9))
    call flv%init ([1011, 1012, 11, 12, 21, 22, 3, 4, 5])
    do i = 1, 9
       call pset%prt(i)%set_flavor (flv(i))
    end do

    call pset%write (u)

    write (u, "(A)")
    write (u, "(A)")  "*   Fill interaction from particle set"
    write (u, "(A)")

    call pset%fill_interaction (int, 2, state_flv=state_flv)
    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call int%final ()
    call pset%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particles_6"

  end subroutine particles_6

  subroutine particles_7 (u)

    integer, intent(in) :: u
    type(interaction_t) :: int
    type(state_flv_content_t) :: state_flv
    type(particle_set_t) :: pset
    type(flavor_t), dimension(:), allocatable :: flv
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer :: i, j

    write (u, "(A)")  "* Test output: Particles"
    write (u, "(A)")  "*   Purpose: reconstruct decay interaction with reordering"
    write (u, "(A)")

    write (u, "(A)")  "* Set up a 1 -> 3 interaction"
    write (u, "(A)")  "    + no quantum numbers"
    write (u, "(A)")

    call reset_interaction_counter ()
    call int%basic_init (0, 1, 3)
    do j = 2, 4
       call int%relate (1, j)
    end do

    allocate (qn (4))
    call int%add_state (qn)
    call int%freeze ()

    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually set up a flavor-content record"
    write (u, "(A)")  "*   assumed interaction: 6 12 5 -11"
    write (u, "(A)")

    call state_flv%init (1, &
         mask = [.false., .true., .true., .true.])
    call state_flv%set_entry (1, &
         pdg = [6, 5, -11, 12], &
         map = [1, 4, 2, 3])

    call state_flv%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Manually create a matching particle set"
    write (u, "(A)")

    pset%n_beam = 0
    pset%n_in   = 1
    pset%n_vir  = 0
    pset%n_out  = 3
    pset%n_tot  = 4
    allocate (pset%prt (pset%n_tot))
    do i = 1, 1
       call pset%prt(i)%reset_status (PRT_INCOMING)
       call pset%prt(i)%set_children ([2,3,4])
    end do
    do i = 2, 4
       call pset%prt(i)%reset_status (PRT_OUTGOING)
       call pset%prt(i)%set_parents ([1])
    end do
    call pset%prt(1)%set_momentum (vector4_at_rest (1._default))
    call pset%prt(2)%set_momentum (vector4_at_rest (3._default))
    call pset%prt(3)%set_momentum (vector4_at_rest (2._default))
    call pset%prt(4)%set_momentum (vector4_at_rest (4._default))

    allocate (flv (4))
    call flv%init ([6,5,12,-11])
    do i = 1, 4
       call pset%prt(i)%set_flavor (flv(i))
    end do

    call pset%write (u)

    write (u, "(A)")
    write (u, "(A)")  "*   Fill interaction from particle set"
    write (u, "(A)")

    call pset%fill_interaction (int, 1, state_flv=state_flv)
    call int%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call int%final ()
    call pset%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particles_7"

  end subroutine particles_7

  subroutine particles_8 (u)
    integer, intent(in) :: u
    type(particle_set_t) :: particle_set
    type(particle_t), dimension(:), allocatable :: particles
    integer, allocatable, dimension(:) :: children, parents
    integer :: n_particles, i
    write (u, "(A)")  "* Test output: particles_8"
    write (u, "(A)")  "*   Purpose: Test functions on particle sets"
    write (u, "(A)")

    call create_test_particle_set_1 (particle_set)
    call particle_set%write (u)
    call assert_equal (u, particle_set%n_tot, 9)
    call assert_equal (u, particle_set%n_beam, 2)
    allocate (children (particle_set%prt(3)%get_n_children ()))
    children = particle_set%prt(3)%get_children()
    call assert_equal (u, particle_set%prt(children(1))%get_pdg (), 3)
    call assert_equal (u, size (particle_set%prt(1)%get_children ()), 2)
    call assert_equal (u, size (particle_set%prt(2)%get_children ()), 2)

    call particle_set%without_hadronic_remnants &
         (particles, n_particles, 3)
    call particle_set%replace (particles)
    write (u, "(A)")
    call particle_set%write (u)

    call assert_equal (u, n_particles, 7)
    call assert_equal (u, size(particles), 10)
    call assert_equal (u, particle_set%n_tot, 10)
    call assert_equal (u, particle_set%n_beam, 2)
    do i = 3, 4
       if (allocated (children))  deallocate (children)
       allocate (children (particle_set%prt(i)%get_n_children ()))
       children = particle_set%prt(i)%get_children()
       call assert_equal (u, particle_set%prt(children(1))%get_pdg (), 3)
       call assert_equal (u, particle_set%prt(children(2))%get_pdg (), 4)
       call assert_equal (u, particle_set%prt(children(3))%get_pdg (), 5)
    end do
    do i = 5, 7
       if (allocated (parents))  deallocate (parents)
       allocate (parents (particle_set%prt(i)%get_n_parents ()))
       parents = particle_set%prt(i)%get_parents()
       call assert_equal (u, particle_set%prt(parents(1))%get_pdg (), 11)
       call assert_equal (u, particle_set%prt(parents(2))%get_pdg (), 12)
    end do
    call assert_equal (u, size (particle_set%prt(1)%get_children ()), &
         1, "get children of 1")
    call assert_equal (u, size (particle_set%prt(2)%get_children ()), &
         1, "get children of 2")

    call assert_equal (u, particle_set%find_particle &
         (particle_set%prt(1)%get_pdg (), particle_set%prt(1)%p), &
         1, "find 1st particle")
    call assert_equal (u, particle_set%find_particle &
         (particle_set%prt(2)%get_pdg (), particle_set%prt(2)%p * &
          (one + tiny_07), rel_smallness=1.0E-6_default), &
         2, "find 2nd particle fuzzy")

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particles_8"
  end subroutine particles_8

  subroutine particles_9 (u)
    integer, intent(in) :: u
    write (u, "(A)")  "* Test output: particles_9"
    write (u, "(A)")  "*   Purpose: Order into Lund strings, "
    write (u, "(A)")  "*              uncolored beam remnants"
    write (u, "(A)")
  end subroutine particles_9


  subroutine create_test_particle_set_1 (pset)
    type(particle_set_t), intent(out) :: pset
    type(flavor_t), dimension(:), allocatable :: flv
    integer :: i
    pset%n_beam = 2
    pset%n_in   = 2
    pset%n_vir  = 2
    pset%n_out  = 3
    pset%n_tot  = 9

    allocate (pset%prt (pset%n_tot))
    call pset%prt(1)%reset_status (PRT_BEAM)
    call pset%prt(2)%reset_status (PRT_BEAM)
    call pset%prt(3)%reset_status (PRT_INCOMING)
    call pset%prt(4)%reset_status (PRT_INCOMING)
    call pset%prt(5)%reset_status (PRT_BEAM_REMNANT)
    call pset%prt(6)%reset_status (PRT_BEAM_REMNANT)
    call pset%prt(7)%reset_status (PRT_OUTGOING)
    call pset%prt(8)%reset_status (PRT_OUTGOING)
    call pset%prt(9)%reset_status (PRT_OUTGOING)

    call pset%prt(1)%set_children ([3,5])
    call pset%prt(2)%set_children ([4,6])
    call pset%prt(3)%set_children ([7,8,9])
    call pset%prt(4)%set_children ([7,8,9])

    call pset%prt(3)%set_parents ([1])
    call pset%prt(4)%set_parents ([2])
    call pset%prt(5)%set_parents ([1])
    call pset%prt(6)%set_parents ([2])
    call pset%prt(7)%set_parents ([3,4])
    call pset%prt(8)%set_parents ([3,4])
    call pset%prt(9)%set_parents ([3,4])

    call pset%prt(1)%set_momentum (vector4_at_rest (1._default))
    call pset%prt(2)%set_momentum (vector4_at_rest (2._default))
    call pset%prt(3)%set_momentum (vector4_at_rest (4._default))
    call pset%prt(4)%set_momentum (vector4_at_rest (6._default))
    call pset%prt(5)%set_momentum (vector4_at_rest (3._default))
    call pset%prt(6)%set_momentum (vector4_at_rest (5._default))
    call pset%prt(7)%set_momentum (vector4_at_rest (7._default))
    call pset%prt(8)%set_momentum (vector4_at_rest (8._default))
    call pset%prt(9)%set_momentum (vector4_at_rest (9._default))

    allocate (flv (9))
    call flv%init ([2011, 2012, 11, 12, 91, 92, 3, 4, 5])
    do i = 1, 9
       call pset%prt(i)%set_flavor (flv(i))
    end do
  end subroutine create_test_particle_set_1


end module particles_uti
