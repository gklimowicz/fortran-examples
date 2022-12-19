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

module state_matrices_uti

  use kinds, only: default
  use io_units
  use format_defs, only: FMT_19
  use flavors
  use colors
  use helicities
  use quantum_numbers

  use state_matrices

  implicit none
  private

  public :: state_matrix_1
  public :: state_matrix_2
  public :: state_matrix_3
  public :: state_matrix_4
  public :: state_matrix_5
  public :: state_matrix_6
  public :: state_matrix_7

contains

  subroutine state_matrix_1 (u)
    integer, intent(in) :: u
    type(state_matrix_t) :: state1, state2, state3
    type(flavor_t), dimension(3) :: flv
    type(color_t), dimension(3) :: col
    type(quantum_numbers_t), dimension(3) :: qn

    write (u, "(A)")  "* Test output: state_matrix_1"
    write (u, "(A)")  "*   Purpose: create and merge two quantum states"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    write (u, "(A)")  "*  State matrix 1"
    write (u, "(A)")

    call state1%init ()
    call flv%init ([1, 2, 11])
    call qn%init (flv, helicity ([ 1, 1, 1]))
    call state1%add_state (qn)
    call qn%init (flv, helicity ([ 1, 1, 1], [-1, 1, -1]))
    call state1%add_state (qn)
    call state1%freeze ()
    call state1%write (u)

    write (u, "(A)")
    write (u, "(A)")  "*  State matrix 2"
    write (u, "(A)")

    call state2%init ()
    call col(1)%init ([501])
    call col(2)%init ([-501])
    call col(3)%init ([0])
    call qn%init (col, helicity ([-1, -1, 0]))
    call state2%add_state (qn)
    call col(3)%init ([99])
    call qn%init (col, helicity ([-1, -1, 0]))
    call state2%add_state (qn)
    call state2%freeze ()
    call state2%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Merge the state matrices"
    write (u, "(A)")

    call merge_state_matrices (state1, state2, state3)
    call state3%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Collapse the state matrix"
    write (u, "(A)")

    call state3%collapse (quantum_numbers_mask (.false., .false., &
         [.true.,.false.,.false.]))
    call state3%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    call state1%final ()
    call state2%final ()
    call state3%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: state_matrix_1"
    write (u, "(A)")

  end subroutine state_matrix_1

  subroutine state_matrix_2 (u)
    integer, intent(in) :: u
    type(state_matrix_t) :: state
    type(state_matrix_t), dimension(:), allocatable :: single_state
    type(state_matrix_t) :: correlated_state
    integer :: f, h11, h12, h21, h22, i, mode
    type(flavor_t), dimension(2) :: flv
    type(color_t), dimension(2) :: col
    type(helicity_t), dimension(2) :: hel
    type(quantum_numbers_t), dimension(2) :: qn
    logical :: ok

    write (u, "(A)")
    write (u, "(A)")  "* Test output: state_matrix_2"
    write (u, "(A)")  "*   Purpose: factorize correlated 3-particle state"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call state%init ()
    do f = 1, 2
       do h11 = -1, 1, 2
          do h12 = -1, 1, 2
             do h21 = -1, 1, 2
                do h22 = -1, 1, 2
                   call flv%init ([f, -f])
                   call col(1)%init ([1])
                   call col(2)%init ([-1])
                   call hel%init ([h11,h12], [h21, h22])
                   call qn%init (flv, col, hel)
                   call state%add_state (qn)
                end do
             end do
          end do
       end do
    end do
    call state%freeze ()
    call state%write (u)

    write (u, "(A)")
    write (u, "(A,'('," // FMT_19 // ",','," // FMT_19 // ",')')") &
         "* Trace = ", state%trace ()
    write (u, "(A)")

    do mode = 1, 3
       write (u, "(A)")
       write (u, "(A,I1)")  "* Mode = ", mode
       call state%factorize &
            (mode, 0.15_default, ok, single_state, correlated_state)
       do i = 1, size (single_state)
          write (u, "(A)")
          call single_state(i)%write (u)
          write (u, "(A,'('," // FMT_19 // ",','," // FMT_19 // ",')')") &
               "Trace = ", single_state(i)%trace ()
       end do
       write (u, "(A)")
       call correlated_state%write (u)
       write (u, "(A,'('," // FMT_19 // ",','," // FMT_19 // ",')')")  &
            "Trace = ", correlated_state%trace ()
       do i = 1, size(single_state)
          call single_state(i)%final ()
       end do
       call correlated_state%final ()
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call state%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: state_matrix_2"

  end subroutine state_matrix_2

  subroutine state_matrix_3 (u)
    use physics_defs, only: HADRON_REMNANT_TRIPLET, HADRON_REMNANT_OCTET
    integer, intent(in) :: u
    type(state_matrix_t) :: state
    type(flavor_t), dimension(4) :: flv
    type(color_t), dimension(4) :: col
    type(quantum_numbers_t), dimension(4) :: qn

    write (u, "(A)")  "* Test output: state_matrix_3"
    write (u, "(A)")  "*   Purpose: add color connections to colored state"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    call state%init ()
    call flv%init ([ 1, -HADRON_REMNANT_TRIPLET, -1, HADRON_REMNANT_TRIPLET ])
    call col(1)%init ([17])
    call col(2)%init ([-17])
    call col(3)%init ([-19])
    call col(4)%init ([19])
    call qn%init (flv, col)
    call state%add_state (qn)
    call flv%init ([ 1, -HADRON_REMNANT_TRIPLET, 21, HADRON_REMNANT_OCTET ])
    call col(1)%init ([17])
    call col(2)%init ([-17])
    call col(3)%init ([3, -5])
    call col(4)%init ([5, -3])
    call qn%init (flv, col)
    call state%add_state (qn)
    call state%freeze ()

    write (u, "(A)") "* State:"
    write (u, "(A)")

    call state%write (u)
    call state%add_color_contractions ()

    write (u, "(A)") "* State with contractions:"
    write (u, "(A)")

    call state%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call state%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: state_matrx_3"

  end subroutine state_matrix_3

  subroutine state_matrix_4 (u)
    integer, intent(in) :: u
    type(state_matrix_t), allocatable :: state
    integer :: f, h11, h12, h21, h22, i
    type(flavor_t), dimension(2) :: flv
    type(color_t), dimension(2) :: col
    type(helicity_t), dimension(2) :: hel
    type(quantum_numbers_t), dimension(2) :: qn
    integer :: unit, iostat

    write (u, "(A)")
    write (u, "(A)")  "* Test output: state_matrix_4"
    write (u, "(A)")  "*   Purpose: raw I/O for correlated 3-particle state"
    write (u, "(A)")

    write (u, "(A)")  "*  Initialization"
    write (u, "(A)")

    allocate (state)

    call state%init ()
    do f = 1, 2
       do h11 = -1, 1, 2
          do h12 = -1, 1, 2
             do h21 = -1, 1, 2
                do h22 = -1, 1, 2
                   call flv%init ([f, -f])
                   call col(1)%init ([1])
                   call col(2)%init ([-1])
                   call hel%init ([h11, h12], [h21, h22])
                   call qn%init (flv, col, hel)
                   call state%add_state (qn)
                end do
             end do
          end do
       end do
    end do
    call state%freeze ()

    call state%set_norm (3._default)
    do i = 1, state%get_n_leaves ()
       call state%set_matrix_element (i, cmplx (2 * i, 2 * i + 1, default))
    end do

    call state%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Write to file and read again "
    write (u, "(A)")

    unit = free_unit ()
    open (unit, action="readwrite", form="unformatted", status="scratch")
    call state%write_raw (unit)
    call state%final ()
    deallocate (state)

    allocate(state)
    rewind (unit)
    call state%read_raw (unit, iostat=iostat)
    close (unit)

    call state%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call state%final ()
    deallocate (state)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: state_matrix_4"

  end subroutine state_matrix_4

  subroutine state_matrix_5 (u)
    integer, intent(in) :: u
    type(state_matrix_t), allocatable, target :: state
    type(state_iterator_t) :: it
    type(state_flv_content_t), allocatable :: state_flv
    type(flavor_t), dimension(4) :: flv1, flv2, flv3, flv4
    type(color_t), dimension(4) :: col1, col2
    type(helicity_t), dimension(4) :: hel1, hel2, hel3
    type(quantum_numbers_t), dimension(4) :: qn
    logical, dimension(4) :: mask

    write (u, "(A)")  "* Test output: state_matrix_5"
    write (u, "(A)")  "*   Purpose: check flavor-content state"
    write (u, "(A)")

    write (u, "(A)")  "* Set up arbitrary state matrix"
    write (u, "(A)")

    call flv1%init ([1, 4, 2, 7])
    call flv2%init ([1, 3,-3, 8])
    call flv3%init ([5, 6, 3, 7])
    call flv4%init ([6, 3, 5, 8])
    call hel1%init ([0, 1, -1, 0])
    call hel2%init ([0, 1, 1, 1])
    call hel3%init ([1, 0, 0, 0])
    call col1(1)%init ([0])
    call col1(2)%init ([0])
    call col1(3)%init ([0])
    call col1(4)%init ([0])
    call col2(1)%init ([5, -6])
    call col2(2)%init ([0])
    call col2(3)%init ([6, -5])
    call col2(4)%init ([0])

    allocate (state)
    call state%init ()
    call qn%init (flv1, col1, hel1)
    call state%add_state (qn)
    call qn%init (flv1, col1, hel2)
    call state%add_state (qn)
    call qn%init (flv3, col1, hel3)
    call state%add_state (qn)
    call qn%init (flv4, col1, hel3)
    call state%add_state (qn)
    call qn%init (flv1, col2, hel3)
    call state%add_state (qn)
    call qn%init (flv2, col2, hel2)
    call state%add_state (qn)
    call qn%init (flv2, col2, hel1)
    call state%add_state (qn)
    call qn%init (flv2, col1, hel1)
    call state%add_state (qn)
    call qn%init (flv3, col1, hel1)
    call state%add_state (qn)
    call qn%init (flv3, col2, hel3)
    call state%add_state (qn)
    call qn%init (flv1, col1, hel1)
    call state%add_state (qn)

    write (u, "(A)")  "* Quantum number content"
    write (u, "(A)")

    call it%init (state)
    do while (it%is_valid ())
       call quantum_numbers_write (it%get_quantum_numbers (), u)
       write (u, *)
       call it%advance ()
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Extract the flavor content"
    write (u, "(A)")

    mask = [.true., .true., .true., .false.]

    allocate (state_flv)
    call state_flv%fill (state, mask)
    call state_flv%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Match trial sets"
    write (u, "(A)")

    call check ([1, 2, 3, 0])
    call check ([1, 4, 2, 0])
    call check ([4, 2, 1, 0])
    call check ([1, 3, -3, 0])
    call check ([1, -3, 3, 0])
    call check ([6, 3, 5, 0])

    write (u, "(A)")
    write (u, "(A)")  "* Determine the flavor content with mask"
    write (u, "(A)")

    mask = [.false., .true., .true., .false.]

    call state_flv%fill (state, mask)
    call state_flv%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Match trial sets"
    write (u, "(A)")

    call check ([1, 2, 3, 0])
    call check ([1, 4, 2, 0])
    call check ([4, 2, 1, 0])
    call check ([1, 3, -3, 0])
    call check ([1, -3, 3, 0])
    call check ([6, 3, 5, 0])

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    deallocate (state_flv)

    call state%final ()
    deallocate (state)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: state_matrix_5"

  contains

    subroutine check (pdg)
      integer, dimension(4), intent(in) :: pdg
      integer, dimension(4) :: map
      logical :: success
      call state_flv%match (pdg, success, map)
      write (u, "(2x,4(1x,I0),':',1x,L1)", advance="no")  pdg, success
      if (success) then
         write (u, "(2x,'map = (',4(1x,I0),' )')")  map
      else
         write (u, *)
      end if
    end subroutine check

  end subroutine state_matrix_5

  subroutine state_matrix_6 (u)
    integer, intent(in) :: u
    type(state_matrix_t), allocatable :: state_orig, state_reduced
    type(flavor_t), dimension(4) :: flv
    type(helicity_t), dimension(4) :: hel
    type(color_t), dimension(4) :: col
    type(quantum_numbers_t), dimension(4) :: qn
    type(quantum_numbers_mask_t), dimension(4) :: qn_mask
    integer :: h1, h2, h3 , h4
    integer :: n_states = 0

    write (u, "(A)") "* Test output: state_matrix_6"
    write (u, "(A)") "* Purpose: Check state matrix reduction"
    write (u, "(A)")

    write (u, "(A)") "* Set up helicity-diagonal state matrix"
    write (u, "(A)")

    allocate (state_orig)
    call state_orig%init ()

    call flv%init ([11, -11, 1, -1])
    call col(3)%init ([1])
    call col(4)%init ([-1])
    do h1 = -1, 1, 2
       do h2 = -1, 1, 2
          do h3 = -1, 1, 2
             do h4 = -1, 1, 2
                n_states = n_states + 1
                call hel%init ([h1, h2, h3, h4], [h1, h2, h3, h4])
                call qn%init (flv, col, hel)
                call state_orig%add_state (qn)
             end do
          end do
       end do
    end do
    call state_orig%freeze ()

    write (u, "(A)") "* Original state: "
    write (u, "(A)")
    call state_orig%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Setup quantum mask: "

    call qn_mask%init ([.false., .false., .false., .false.], &
                       [.true., .true., .true., .true.], &
                       [.false., .false., .true., .true.])
    call quantum_numbers_mask_write (qn_mask, u)
    write (u, "(A)")
    write (u, "(A)") "* Reducing the state matrix using above mask"
    write (u, "(A)")
    allocate (state_reduced)
    call state_orig%reduce (qn_mask, state_reduced)

    write (u, "(A)") "* Reduced state matrix: "
    call state_reduced%write (u)

    write (u, "(A)") "* Test output end: state_matrix_6"
  end subroutine state_matrix_6

  subroutine state_matrix_7 (u)
    integer, intent(in) :: u
    type(state_matrix_t), allocatable :: state_orig, state_reduced, &
         state_ordered
    type(flavor_t), dimension(4) :: flv
    type(helicity_t), dimension(4) :: hel
    type(color_t), dimension(4) :: col
    type(quantum_numbers_t), dimension(4) :: qn
    type(quantum_numbers_mask_t), dimension(4) :: qn_mask
    integer :: h1, h2, h3 , h4
    integer :: n_states = 0

    write (u, "(A)") "* Test output: state_matrix_7"
    write (u, "(A)") "* Purpose: Check ordered state matrix reduction"
    write (u, "(A)")

    write (u, "(A)") "* Set up helicity-diagonal state matrix"
    write (u, "(A)")

    allocate (state_orig)
    call state_orig%init ()

    call flv%init ([11, -11, 1, -1])
    call col(3)%init ([1])
    call col(4)%init ([-1])
    do h1 = -1, 1, 2
       do h2 = -1, 1, 2
          do h3 = -1, 1, 2
             do h4 = -1, 1, 2
                n_states = n_states + 1
                call hel%init ([h1, h2, h3, h4], [h1, h2, h3, h4])
                call qn%init (flv, col, hel)
                call state_orig%add_state (qn)
             end do
          end do
       end do
    end do
    call state_orig%freeze ()

    write (u, "(A)") "* Original state: "
    write (u, "(A)")
    call state_orig%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Setup quantum mask: "

    call qn_mask%init ([.false., .false., .false., .false.], &
                       [.true., .true., .true., .true.], &
                       [.false., .false., .true., .true.])
    call quantum_numbers_mask_write (qn_mask, u)
    write (u, "(A)")
    write (u, "(A)") "* Reducing the state matrix using above mask and keeping the old indices:"
    write (u, "(A)")
    allocate (state_reduced)
    call state_orig%reduce (qn_mask, state_reduced, keep_me_index = .true.)

    write (u, "(A)") "* Reduced state matrix with kept indices: "
    call state_reduced%write (u)

    write (u, "(A)")
    write (u, "(A)") "* Reordering reduced state matrix:"
    write (u, "(A)")
    allocate (state_ordered)
    call state_reduced%reorder_me (state_ordered)

    write (u, "(A)") "* Reduced and ordered state matrix:"
    call state_ordered%write (u)

    write (u, "(A)") "* Test output end: state_matrix_6"
  end subroutine state_matrix_7


end module state_matrices_uti
