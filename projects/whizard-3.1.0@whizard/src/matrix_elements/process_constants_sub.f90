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

submodule (process_constants) process_constants_s

  use io_units, only: given_output_unit, free_unit
  use format_utils, only: write_integer_array
  use md5, only: md5sum

  implicit none

contains

  elemental module function process_constants_get_n_tot (prc_const) result (n_tot)
    integer :: n_tot
    class(process_constants_t), intent(in) :: prc_const
    n_tot = prc_const%n_in + prc_const%n_out
  end function process_constants_get_n_tot

  module subroutine process_constants_get_flv_state (prc_const, flv_state)
    class(process_constants_t), intent(in) :: prc_const
    integer, dimension(:,:), allocatable, intent(out) :: flv_state
    allocate (flv_state (size (prc_const%flv_state, 1), &
         size (prc_const%flv_state, 2)))
    flv_state = prc_const%flv_state
  end subroutine process_constants_get_flv_state

  module function process_constants_get_n_flv (data) result (n_flv)
    integer :: n_flv
    class(process_constants_t), intent(in) :: data
    n_flv = data%n_flv
  end function process_constants_get_n_flv

  module function process_constants_get_n_hel (data) result (n_hel)
    integer :: n_hel
    class(process_constants_t), intent(in) :: data
    n_hel = data%n_hel
  end function process_constants_get_n_hel

  module subroutine process_constants_get_hel_state (prc_const, hel_state)
    class(process_constants_t), intent(in) :: prc_const
    integer, dimension(:,:), allocatable, intent(out) :: hel_state
    allocate (hel_state (size (prc_const%hel_state, 1), &
         size (prc_const%hel_state, 2)))
    hel_state = prc_const%hel_state
  end subroutine process_constants_get_hel_state

  module subroutine process_constants_get_col_state (prc_const, col_state)
    class(process_constants_t), intent(in) :: prc_const
    integer, dimension(:,:,:), allocatable, intent(out) :: col_state
    allocate (col_state (size (prc_const%col_state, 1), &
         size (prc_const%col_state, 2), size (prc_const%col_state, 3)))
    col_state = prc_const%col_state
  end subroutine process_constants_get_col_state

  module subroutine process_constants_get_ghost_flag (prc_const, ghost_flag)
    class(process_constants_t), intent(in) :: prc_const
    logical, dimension(:,:), allocatable, intent(out) :: ghost_flag
    allocate (ghost_flag (size (prc_const%ghost_flag, 1), &
         size (prc_const%ghost_flag, 2)))
    ghost_flag = prc_const%ghost_flag
  end subroutine process_constants_get_ghost_flag

  module subroutine process_constants_get_color_factors (prc_const, col_facts)
    class(process_constants_t), intent(in) :: prc_const
    complex(default), dimension(:), allocatable, intent(out) :: col_facts
    allocate (col_facts (size (prc_const%color_factors)))
    col_facts = prc_const%color_factors
  end subroutine process_constants_get_color_factors

  module subroutine process_constants_get_cf_index (prc_const, cf_index)
    class(process_constants_t), intent(in) :: prc_const
    integer, intent(out), dimension(:,:), allocatable :: cf_index
    allocate (cf_index (size (prc_const%cf_index, 1), &
         size (prc_const%cf_index, 2)))
    cf_index = prc_const%cf_index
  end subroutine process_constants_get_cf_index

  module subroutine process_constants_set_flv_state (prc_const, flv_state)
    class(process_constants_t), intent(inout) :: prc_const
    integer, intent(in), dimension(:,:), allocatable :: flv_state
    if (allocated (prc_const%flv_state)) deallocate (prc_const%flv_state)
    allocate (prc_const%flv_state (size (flv_state, 1), &
         size (flv_state, 2)))
    prc_const%flv_state = flv_state
    prc_const%n_flv = size (flv_state, 2)
  end subroutine process_constants_set_flv_state

  module subroutine process_constants_set_col_state (prc_const, col_state)
    class(process_constants_t), intent(inout) :: prc_const
    integer, intent(in), dimension(:,:,:), allocatable :: col_state
    allocate (prc_const%col_state (size (col_state, 1), &
         size (col_state, 2), size (col_state, 3)))
    prc_const%col_state = col_state
  end subroutine process_constants_set_col_state

  module subroutine process_constants_set_cf_index (prc_const, cf_index)
    class(process_constants_t), intent(inout) :: prc_const
    integer, dimension(:,:), intent(in), allocatable :: cf_index
    allocate (prc_const%cf_index (size (cf_index, 1), &
         size (cf_index, 2)))
    prc_const%cf_index = cf_index
  end subroutine process_constants_set_cf_index

  module subroutine process_constants_set_color_factors (prc_const, color_factors)
    class(process_constants_t), intent(inout) :: prc_const
    complex(default), dimension(:), intent(in), allocatable :: color_factors
    allocate (prc_const%color_factors (size (color_factors)))
    prc_const%color_factors = color_factors
  end subroutine process_constants_set_color_factors

  module subroutine process_constants_set_ghost_flag (prc_const, ghost_flag)
    class(process_constants_t), intent(inout) :: prc_const
    logical, dimension(:,:), allocatable, intent(in) :: ghost_flag
    allocate (prc_const%ghost_flag (size (ghost_flag, 1), &
         size (ghost_flag, 2)))
    prc_const%ghost_flag = ghost_flag
  end subroutine process_constants_set_ghost_flag
  module function process_constants_get_pdg_in (prc_const) result (pdg_in)
    type(pdg_array_t), dimension(:), allocatable :: pdg_in
    class(process_constants_t), intent(in) :: prc_const
    type(pdg_array_t) :: pdg_tmp
    integer :: i
    allocate (pdg_in (prc_const%n_in))
    do i = 1, prc_const%n_in
       pdg_tmp = prc_const%flv_state(i,:)
       pdg_in(i) = sort_abs (pdg_tmp, unique = .true.)
    end do
  end function process_constants_get_pdg_in

  module subroutine process_constants_compute_md5sum (prc_const, include_id)
    class(process_constants_t), intent(inout) :: prc_const
    logical, intent(in) :: include_id
    integer :: unit
    unit = prc_const%fill_unit_for_md5sum (include_id)
    rewind (unit)
    prc_const%md5sum = md5sum (unit)
    close (unit)
  end subroutine process_constants_compute_md5sum

  module function process_constants_fill_unit_for_md5sum &
       (prc_const, include_id) result (unit)
    integer :: unit
    class(process_constants_t), intent(in) :: prc_const
    logical, intent(in) :: include_id
    integer :: i, j, k
    unit = free_unit ()
    open (unit, status="scratch", action="readwrite")
    if (include_id) write (unit, '(A)') char (prc_const%id)
    write (unit, '(A)') char (prc_const%model_name)
    write (unit, '(L1)') prc_const%openmp_supported
    write (unit, '(I0)') prc_const%n_in
    write (unit, '(I0)') prc_const%n_out
    write (unit, '(I0)') prc_const%n_flv
    write (unit, '(I0)') prc_const%n_hel
    write (unit, '(I0)') prc_const%n_col
    write (unit, '(I0)') prc_const%n_cin
    write (unit, '(I0)') prc_const%n_cf
    do i = 1, size (prc_const%flv_state, dim=1)
       do j = 1, size (prc_const%flv_state, dim=2)
          write (unit, '(I0)') prc_const%flv_state (i, j)
       end do
    end do
    do i = 1, size (prc_const%hel_state, dim=1)
       do j = 1, size (prc_const%hel_state, dim=2)
          write (unit, '(I0)') prc_const%hel_state (i, j)
       end do
    end do
    do i = 1, size (prc_const%col_state, dim=1)
       do j = 1, size (prc_const%col_state, dim=2)
          do k = 1, size (prc_const%col_state, dim=3)
             write (unit, '(I0)') prc_const%col_state (i, j, k)
          end do
      end do
    end do
    do i = 1, size (prc_const%ghost_flag, dim=1)
       do j = 1, size (prc_const%ghost_flag, dim=2)
          write (unit, '(L1)') prc_const%ghost_flag (i, j)
       end do
    end do
    do i = 1, size (prc_const%color_factors)
       write (unit, '(F0.0,F0.0)') real (prc_const%color_factors(i)), &
          aimag (prc_const%color_factors(i))
    end do
    do i = 1, size (prc_const%cf_index, dim=1)
       do j = 1, size (prc_const%cf_index, dim=2)
          write (unit, '(I0)') prc_const%cf_index(i, j)
       end do
    end do
  end function process_constants_fill_unit_for_md5sum

  module subroutine process_constants_write (prc_const, unit)
    class(process_constants_t), intent(in) :: prc_const
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A,A)") "Process data of id: ", char (prc_const%id)
    write (u, "(1x,A,A)") "Associated model: ", char (prc_const%model_name)
    write (u, "(1x,A,I0)") "n_in: ", prc_const%n_in
    write (u, "(1x,A,I0)") "n_out: ", prc_const%n_out
    write (u, "(1x,A,I0)") "n_flv: ", prc_const%n_flv
    write (u, "(1x,A,I0)") "n_hel: ", prc_const%n_hel
    write (u, "(1x,A,I0)") "n_col: ", prc_const%n_col
    write (u, "(1x,A,I0)") "n_cin: ", prc_const%n_cin
    write (u, "(1x,A,I0)") "n_cf: ", prc_const%n_cf
    write (u, "(1x,A)") "Flavors: "
    do i = 1, prc_const%n_flv
       write (u, "(1x,A,I0)") "i_flv: ", i
       call write_integer_array (prc_const%flv_state (:,i))
    end do
  end subroutine process_constants_write


end submodule process_constants_s

