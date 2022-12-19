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

submodule (pdg_arrays) pdg_arrays_s

  use io_units
  use sorting
  use physics_defs

  implicit none

contains

  module subroutine pdg_array_write (aval, unit)
    class(pdg_array_t), intent(in) :: aval
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A)", advance="no")  "PDG("
    if (allocated (aval%pdg)) then
       do i = 1, size (aval%pdg)
          if (i > 1)  write (u, "(A)", advance="no")  ", "
          write (u, "(I0)", advance="no")  aval%pdg(i)
       end do
    end if
    write (u, "(A)", advance="no")  ")"
  end subroutine pdg_array_write

  module subroutine pdg_array_write_set (aval, unit)
    type(pdg_array_t), intent(in), dimension(:) :: aval
    integer, intent(in), optional :: unit
    integer :: i
    do i = 1, size (aval)
       call aval(i)%write (unit)
       print *, ''
    end do
  end subroutine pdg_array_write_set

  module subroutine pdg_array_from_int_array (aval, iarray)
    type(pdg_array_t), intent(out) :: aval
    integer, dimension(:), intent(in) :: iarray
    allocate (aval%pdg (size (iarray)))
    aval%pdg = iarray
  end subroutine pdg_array_from_int_array

  elemental module subroutine pdg_array_from_int (aval, int)
    type(pdg_array_t), intent(out) :: aval
    integer, intent(in) :: int
    allocate (aval%pdg (1))
    aval%pdg = int
  end subroutine pdg_array_from_int

  module subroutine int_array_from_pdg_array (iarray, aval)
    integer, dimension(:), allocatable, intent(out) :: iarray
    type(pdg_array_t), intent(in) :: aval
    if (allocated (aval%pdg)) then
       allocate (iarray (size (aval%pdg)))
       iarray = aval%pdg
    else
       allocate (iarray (0))
    end if
  end subroutine int_array_from_pdg_array

  module subroutine pdg_array_init (aval, n_elements)
    class(pdg_array_t), intent(inout) :: aval
    integer, intent(in) :: n_elements
    allocate(aval%pdg(n_elements))
  end subroutine pdg_array_init

  module subroutine pdg_array_delete (aval)
    class(pdg_array_t), intent(inout) :: aval
    if (allocated (aval%pdg)) deallocate (aval%pdg)
  end subroutine pdg_array_delete

  module subroutine pdg_array_merge (aval1, aval2)
    class(pdg_array_t), intent(inout) :: aval1
    type(pdg_array_t), intent(in) :: aval2
    type(pdg_array_t) :: aval
    if (allocated (aval1%pdg) .and. allocated (aval2%pdg)) then
      if (.not. any (aval1%pdg == aval2%pdg)) aval = aval1 // aval2
    else if (allocated (aval1%pdg)) then
      aval = aval1
    else if (allocated (aval2%pdg)) then
      aval = aval2
    end if
    call pdg_array_delete (aval1)
    call pdg_array_from_int_array (aval1, aval%pdg)
  end subroutine pdg_array_merge

  elemental module function pdg_array_get_length (aval) result (n)
    class(pdg_array_t), intent(in) :: aval
    integer :: n
    if (allocated (aval%pdg)) then
       n = size (aval%pdg)
    else
       n = 0
    end if
  end function pdg_array_get_length

  elemental module function pdg_array_get (aval, i) result (pdg)
    class(pdg_array_t), intent(in) :: aval
    integer, intent(in), optional :: i
    integer :: pdg
    if (present (i)) then
       pdg = aval%pdg(i)
    else
       pdg = aval%pdg(1)
    end if
  end function pdg_array_get

  module subroutine pdg_array_set (aval, i, pdg)
    class(pdg_array_t), intent(inout) :: aval
    integer, intent(in) :: i
    integer, intent(in) :: pdg
    aval%pdg(i) = pdg
  end subroutine pdg_array_set

  module function pdg_array_add (aval, aval_add) result (aval_out)
    type(pdg_array_t) :: aval_out
    class(pdg_array_t), intent(in) :: aval
    type(pdg_array_t), intent(in) :: aval_add
    integer :: n, n_add, i
    n = size (aval%pdg)
    n_add = size (aval_add%pdg)
    allocate (aval_out%pdg (n + n_add))
    aval_out%pdg(1:n) = aval%pdg
    do i = 1, n_add
       aval_out%pdg(n+i) = aval_add%pdg(i)
    end do
  end function pdg_array_add

  module function pdg_array_replace (aval, i, pdg_new) result (aval_new)
    class(pdg_array_t), intent(in) :: aval
    integer, intent(in) :: i
    integer, dimension(:), intent(in) :: pdg_new
    type(pdg_array_t) :: aval_new
    integer :: n, l
    n = size (aval%pdg)
    l = size (pdg_new)
    allocate (aval_new%pdg (n + l - 1))
    aval_new%pdg(:i-1) = aval%pdg(:i-1)
    aval_new%pdg(i:i+l-1) = pdg_new
    aval_new%pdg(i+l:) = aval%pdg(i+1:)
  end function pdg_array_replace

  module function concat_pdg_arrays (aval1, aval2) result (aval)
    type(pdg_array_t) :: aval
    type(pdg_array_t), intent(in) :: aval1, aval2
    integer :: n1, n2
    if (allocated (aval1%pdg) .and. allocated (aval2%pdg)) then
       n1 = size (aval1%pdg)
       n2 = size (aval2%pdg)
       allocate (aval%pdg (n1 + n2))
       aval%pdg(:n1) = aval1%pdg
       aval%pdg(n1+1:) = aval2%pdg
    else if (allocated (aval1%pdg)) then
       aval = aval1
    else if (allocated (aval2%pdg)) then
       aval = aval2
    end if
  end function concat_pdg_arrays

  elemental module function pdg_array_match_integer (aval, pdg) result (flag)
    logical :: flag
    type(pdg_array_t), intent(in) :: aval
    integer, intent(in) :: pdg
    if (allocated (aval%pdg)) then
       flag = pdg == UNDEFINED &
            .or. any (aval%pdg == UNDEFINED) &
            .or. any (aval%pdg == pdg)
    else
       flag = .false.
    end if
  end function pdg_array_match_integer

  elemental module function is_quark (pdg_nr)
    logical :: is_quark
    integer, intent(in) :: pdg_nr
    if (abs (pdg_nr) >= 1 .and. abs (pdg_nr) <= 6) then
       is_quark = .true.
    else
       is_quark = .false.
    end if
  end function is_quark

  elemental module function is_gluon (pdg_nr)
    logical :: is_gluon
    integer, intent(in) :: pdg_nr
    if (pdg_nr == GLUON) then
       is_gluon = .true.
    else
       is_gluon = .false.
    end if
  end function is_gluon

  elemental module function is_photon (pdg_nr)
    logical :: is_photon
    integer, intent(in) :: pdg_nr
    if (pdg_nr == PHOTON) then
       is_photon = .true.
    else
       is_photon = .false.
    end if
  end function is_photon

  elemental module function is_colored (pdg_nr)
    logical :: is_colored
    integer, intent(in) :: pdg_nr
    is_colored = is_quark (pdg_nr) .or. is_gluon (pdg_nr)
  end function is_colored

  elemental module function is_lepton (pdg_nr)
    logical :: is_lepton
    integer, intent(in) :: pdg_nr
    if (abs (pdg_nr) >= ELECTRON .and. &
         abs (pdg_nr) <= TAU_NEUTRINO) then
      is_lepton = .true.
    else
      is_lepton = .false.
    end if
  end function is_lepton

  elemental module function is_charged_lepton (pdg_nr)
    logical :: is_charged_lepton
    integer, intent(in) :: pdg_nr
    if (abs (pdg_nr) == ELECTRON .or. &
         abs (pdg_nr) == MUON .or. &
         abs (pdg_nr) == TAU) then
      is_charged_lepton = .true.
    else
      is_charged_lepton = .false.
    end if
  end function is_charged_lepton
  elemental module function is_fermion (pdg_nr)
    logical :: is_fermion
    integer, intent(in) :: pdg_nr
    is_fermion = is_lepton(pdg_nr) .or. is_quark(pdg_nr)
  end function is_fermion

  elemental module function is_massless_vector (pdg_nr)
    integer, intent(in) :: pdg_nr
    logical :: is_massless_vector
    if (pdg_nr == GLUON .or. pdg_nr == PHOTON) then
      is_massless_vector = .true.
    else
      is_massless_vector = .false.
    end if
  end function is_massless_vector

  elemental module function is_massive_vector (pdg_nr)
    integer, intent(in) :: pdg_nr
    logical :: is_massive_vector
    if (abs (pdg_nr) == Z_BOSON .or. abs (pdg_nr) == W_BOSON) then
      is_massive_vector = .true.
    else
      is_massive_vector = .false.
    end if
  end function is_massive_vector

  elemental module function is_vector (pdg_nr)
    integer, intent(in) :: pdg_nr
    logical :: is_vector
    if (is_massless_vector (pdg_nr) .or. is_massive_vector (pdg_nr)) then
      is_vector = .true.
    else
      is_vector = .false.
    end if
  end function is_vector

  elemental module function is_elementary (pdg_nr)
    integer, intent(in) :: pdg_nr
    logical :: is_elementary
    if (is_vector (pdg_nr) .or. is_fermion (pdg_nr) .or. pdg_nr == 25) then
       is_elementary = .true.
    else
       is_elementary = .false.
    end if
  end function is_elementary

  elemental module function is_ew_boson_scalar (pdg_nr)
    integer, intent(in) :: pdg_nr
    logical :: is_ew_boson_scalar
    if (is_photon (pdg_nr) .or. is_massive_vector (pdg_nr) .or. pdg_nr == 25) then
       is_ew_boson_scalar = .true.
    else
       is_ew_boson_scalar = .false.
    end if
  end function is_ew_boson_scalar

  module function pdg_array_has_colored_particles (pdg) result (colored)
    class(pdg_array_t), intent(in) :: pdg
    logical :: colored
    integer :: i, pdg_nr
    colored = .false.
    do i = 1, size (pdg%pdg)
      pdg_nr = pdg%pdg(i)
      if (is_quark (pdg_nr) .or. is_gluon (pdg_nr)) then
         colored = .true.
         exit
      end if
    end do
  end function pdg_array_has_colored_particles

  module function query_coupling_powers (flv, a_power, as_power) result (valid)
     integer, intent(in), dimension(:) :: flv
     integer, dimension(:, :), allocatable :: power_pair_array
     integer, dimension(2) :: power_pair_ref
     integer, intent(in) :: a_power, as_power
     integer :: i, n_legs, n_gluons, n_quarks, n_gamWZH, n_leptons
     logical, dimension(:), allocatable :: pairs_included
     logical :: valid
     integer :: n_bound
     power_pair_ref = [a_power, as_power]
     n_legs = size (flv)
     allocate (power_pair_array (2, n_legs - 1))
     do i = 1, n_legs - 1
        power_pair_array (1, i) = n_legs - 1 - i
        power_pair_array (2, i) = i - 1
     end do
     allocate (pairs_included (n_legs - 1))
     pairs_included = .true.
     n_gluons = count (is_gluon (flv))
     n_gamWZH = count (is_ew_boson_scalar (flv))
     n_quarks = count (is_quark (flv))
     n_leptons = count (is_lepton (flv))
     if (n_gluons >= 1 .and. n_gluons <= 3) then
        do i = 1, n_gluons
           pairs_included (i) = .false.
        end do
     else if (n_gluons > 2 .and. n_quarks <= 2 .and. n_gluons + n_quarks == n_legs) then
        do i = 1, n_legs - 2
           pairs_included (i) = .false.
        end do
     end if
     n_bound = 0
     if (n_gamWZH + n_leptons == n_legs) then
        n_bound = n_gamWZH + n_leptons - 2
     else if (n_quarks == 2 .and. n_leptons + n_quarks + n_gamWZH == n_legs) then
        n_bound = n_legs - 2
     else if (n_gamWZH + n_leptons > 0) then
        n_bound = n_leptons/2 + n_gamWZH
     end if
     if (n_bound > 0) then
        do i = 1, n_bound
           pairs_included (n_legs - i) = .false.
        end do
     end if
     if (n_quarks == 4 .and. .not. qcd_ew_interferences (flv)) then
        do i = 1, 2
           pairs_included (n_legs - i) = .false.
        end do
     end if
     valid = .false.
     do i = 1, n_legs - 1
        if (all (power_pair_array (:, i) == power_pair_ref) .and. pairs_included (i)) then
           valid = .true.
           exit
        end if
     end do
  end function query_coupling_powers

  module function qcd_ew_interferences (flv) result (valid)
     integer, intent(in), dimension(:) :: flv
     integer :: i, n_pairs
     logical :: valid, qqbar_pair
     n_pairs = 0
     valid = .false.
     qqbar_pair = .false.
     if (count (is_quark (flv)) >= 4) then
        do i = DOWN_Q, TOP_Q
           qqbar_pair = count (abs (flv) == i) >= 2
           if (qqbar_pair) n_pairs = n_pairs + 1
           if (n_pairs > 0) then
              valid = .true.
              exit
           end if
        end do
     end if
  end function qcd_ew_interferences

  module function flv_eqv_expr_class (flv) result (assign_qgA)
     integer, intent(in) :: flv
     logical, dimension(3) :: assign_qgA
     assign_qgA = [is_quark (flv), is_gluon (flv), is_photon (flv)]
  end function flv_eqv_expr_class

  module function pdg_array_match_pdg_array (aval1, aval2) result (flag)
    logical :: flag
    type(pdg_array_t), intent(in) :: aval1, aval2
    if (allocated (aval1%pdg) .and. allocated (aval2%pdg)) then
       flag = any (aval1 .match. aval2%pdg)
    else
       flag = .false.
    end if
  end function pdg_array_match_pdg_array

  elemental module function pdg_array_lt (aval1, aval2) result (flag)
    type(pdg_array_t), intent(in) :: aval1, aval2
    logical :: flag
    integer :: i
    if (size (aval1%pdg) /= size (aval2%pdg)) then
       flag = size (aval1%pdg) < size (aval2%pdg)
    else
       do i = 1, size (aval1%pdg)
          if (abs (aval1%pdg(i)) /= abs (aval2%pdg(i))) then
             flag = abs (aval1%pdg(i)) < abs (aval2%pdg(i))
             return
          end if
       end do
       do i = 1, size (aval1%pdg)
          if (aval1%pdg(i) /= aval2%pdg(i)) then
             flag = aval1%pdg(i) > aval2%pdg(i)
             return
          end if
       end do
       flag = .false.
    end if
  end function pdg_array_lt

  elemental module function pdg_array_gt (aval1, aval2) result (flag)
    type(pdg_array_t), intent(in) :: aval1, aval2
    logical :: flag
    flag = .not. (aval1 < aval2 .or. aval1 == aval2)
  end function pdg_array_gt

  elemental module function pdg_array_le (aval1, aval2) result (flag)
    type(pdg_array_t), intent(in) :: aval1, aval2
    logical :: flag
    flag = aval1 < aval2 .or. aval1 == aval2
  end function pdg_array_le

  elemental module function pdg_array_ge (aval1, aval2) result (flag)
    type(pdg_array_t), intent(in) :: aval1, aval2
    logical :: flag
    flag = .not. (aval1 < aval2)
  end function pdg_array_ge

  elemental module function pdg_array_eq (aval1, aval2) result (flag)
    type(pdg_array_t), intent(in) :: aval1, aval2
    logical :: flag
    if (size (aval1%pdg) /= size (aval2%pdg)) then
       flag = .false.
    else
       flag = all (aval1%pdg == aval2%pdg)
    end if
  end function pdg_array_eq

  elemental module function pdg_array_ne (aval1, aval2) result (flag)
    type(pdg_array_t), intent(in) :: aval1, aval2
    logical :: flag
    flag = .not. (aval1 == aval2)
  end function pdg_array_ne

  elemental module function pdg_array_equivalent (aval1, aval2) result (eq)
    logical :: eq
    type(pdg_array_t), intent(in) :: aval1, aval2
    logical, dimension(:), allocatable :: match1, match2
    integer :: i
    if (allocated (aval1%pdg) .and. allocated (aval2%pdg)) then
       eq = any (aval1%pdg == UNDEFINED) &
            .or. any (aval2%pdg == UNDEFINED)
       if (.not. eq) then
          allocate (match1 (size (aval1%pdg)))
          allocate (match2 (size (aval2%pdg)))
          match1 = .false.
          match2 = .false.
          do i = 1, size (aval1%pdg)
             match2 = match2 .or. aval1%pdg(i) == aval2%pdg
          end do
          do i = 1, size (aval2%pdg)
             match1 = match1 .or. aval2%pdg(i) == aval1%pdg
          end do
          eq = all (match1) .and. all (match2)
       end if
    else
       eq = .false.
    end if
  end function pdg_array_equivalent

  elemental module function pdg_array_inequivalent (aval1, aval2) result (neq)
    logical :: neq
    type(pdg_array_t), intent(in) :: aval1, aval2
    neq = .not. pdg_array_equivalent (aval1, aval2)
  end function pdg_array_inequivalent

  module function pdg_array_sort_abs (aval1, unique) result (aval2)
    class(pdg_array_t), intent(in) :: aval1
    logical, intent(in), optional :: unique
    type(pdg_array_t) :: aval2
    integer, dimension(:), allocatable :: tmp
    logical, dimension(:), allocatable :: mask
    integer :: i, n
    logical :: uni
    uni = .false.;  if (present (unique))  uni = unique
    n = size (aval1%pdg)
    if (uni) then
       allocate (tmp (n), mask(n))
       tmp = sort_abs (aval1%pdg)
       mask(1) = .true.
       do i = 2, n
          mask(i) = tmp(i) /= tmp(i-1)
       end do
       allocate (aval2%pdg (count (mask)))
       aval2%pdg = pack (tmp, mask)
    else
       allocate (aval2%pdg (n))
       aval2%pdg = sort_abs (aval1%pdg)
    end if
  end function pdg_array_sort_abs

  module function pdg_array_intersect (aval1, match) result (aval2)
   class(pdg_array_t), intent(in) :: aval1
   integer, dimension(:) :: match
   type(pdg_array_t) :: aval2
   integer, dimension(:), allocatable :: isec
   integer :: i
   isec = pack (aval1%pdg, [(any(aval1%pdg(i) == match), i=1,size(aval1%pdg))])
   call pdg_array_from_int_array (aval2, isec)
  end function pdg_array_intersect

  elemental module function pdg_array_search_for_particle (pdg, i_part) result (found)
    class(pdg_array_t), intent(in) :: pdg
    integer, intent(in) :: i_part
    logical :: found
    found = any (pdg%pdg == i_part)
  end function pdg_array_search_for_particle

  module function pdg_array_invert (pdg) result (pdg_inverse)
    class(pdg_array_t), intent(in) :: pdg
    type(pdg_array_t) :: pdg_inverse
    integer :: i, n
    n = size (pdg%pdg)
    allocate (pdg_inverse%pdg (n))
    do i = 1, n
       select case (pdg%pdg(i))
       case (GLUON, PHOTON, Z_BOSON, 25)
          pdg_inverse%pdg(i) = pdg%pdg(i)
       case default
          pdg_inverse%pdg(i) = -pdg%pdg(i)
       end select
    end do
  end function pdg_array_invert

  module subroutine pdg_list_write (object, unit)
    class(pdg_list_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    if (allocated (object%a)) then
       do i = 1, size (object%a)
          if (i > 1)  write (u, "(A)", advance="no")  ", "
          call object%a(i)%write (u)
       end do
    end if
  end subroutine pdg_list_write

  module subroutine pdg_list_init_size (pl, n)
    class(pdg_list_t), intent(out) :: pl
    integer, intent(in) :: n
    allocate (pl%a (n))
  end subroutine pdg_list_init_size

  module subroutine pdg_list_init_int_array (pl, pdg)
    class(pdg_list_t), intent(out) :: pl
    integer, dimension(:), intent(in) :: pdg
    integer :: i
    allocate (pl%a (size (pdg)))
    do i = 1, size (pdg)
       call pdg_array_from_int (pl%a(i), pdg(i))
    end do
  end subroutine pdg_list_init_int_array

  module subroutine pdg_list_set_int (pl, i, pdg)
    class(pdg_list_t), intent(inout) :: pl
    integer, intent(in) :: i
    integer, intent(in) :: pdg
    call pdg_array_from_int (pl%a(i), pdg)
  end subroutine pdg_list_set_int

  module subroutine pdg_list_set_int_array (pl, i, pdg)
    class(pdg_list_t), intent(inout) :: pl
    integer, intent(in) :: i
    integer, dimension(:), intent(in) :: pdg
    call pdg_array_from_int_array (pl%a(i), pdg)
  end subroutine pdg_list_set_int_array

  module subroutine pdg_list_set_pdg_array (pl, i, pa)
    class(pdg_list_t), intent(inout) :: pl
    integer, intent(in) :: i
    type(pdg_array_t), intent(in) :: pa
    pl%a(i) = pa
  end subroutine pdg_list_set_pdg_array

  module function pdg_list_get_size (pl) result (n)
    class(pdg_list_t), intent(in) :: pl
    integer :: n
    if (allocated (pl%a)) then
       n = size (pl%a)
    else
       n = 0
    end if
  end function pdg_list_get_size

  module function pdg_list_get (pl, i) result (pa)
    type(pdg_array_t) :: pa
    class(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: i
    pa = pl%a(i)
  end function pdg_list_get

  module function pdg_list_is_regular (pl) result (flag)
    class(pdg_list_t), intent(in) :: pl
    logical :: flag
    integer :: i, j, s
    s = pl%get_size ()
    flag = .true.
    do i = 1, s
       do j = i + 1, s
          if (pl%a(i) .match. pl%a(j)) then
             if (pl%a(i) /= pl%a(j)) then
                flag = .false.
                return
             end if
          end if
       end do
    end do
  end function pdg_list_is_regular

  module function pdg_list_sort_abs (pl, n_in) result (pl_sorted)
    class(pdg_list_t), intent(in) :: pl
    integer, intent(in), optional :: n_in
    type(pdg_list_t) :: pl_sorted
    type(pdg_array_t), dimension(:), allocatable :: pa
    integer, dimension(:), allocatable :: pdg, map
    integer :: i, n0
    call pl_sorted%init (pl%get_size ())
    if (allocated (pl%a)) then
       allocate (pa (size (pl%a)))
       do i = 1, size (pl%a)
          pa(i) = pl%a(i)%sort_abs (unique = .true.)
       end do
       allocate (pdg (size (pa)), source = 0)
       do i = 1, size (pa)
          if (allocated (pa(i)%pdg)) then
             if (size (pa(i)%pdg) > 0) then
                pdg(i) = pa(i)%pdg(1)
             end if
          end if
       end do
       if (present (n_in)) then
          n0 = n_in
       else
          n0 = 0
       end if
       allocate (map (size (pdg)))
       map(:n0) = [(i, i = 1, n0)]
       map(n0+1:) = n0 + order_abs (pdg(n0+1:))
       do i = 1, size (pa)
          call pl_sorted%set (i, pa(map(i)))
       end do
    end if
  end function pdg_list_sort_abs

  module function pdg_list_eq (pl1, pl2) result (flag)
    class(pdg_list_t), intent(in) :: pl1, pl2
    logical :: flag
    integer :: i
    flag = .false.
    if (allocated (pl1%a) .and. allocated (pl2%a)) then
       if (size (pl1%a) == size (pl2%a)) then
          do i = 1, size (pl1%a)
             associate (a1 => pl1%a(i), a2 => pl2%a(i))
               if (allocated (a1%pdg) .and. allocated (a2%pdg)) then
                  if (size (a1%pdg) == size (a2%pdg)) then
                     if (size (a1%pdg) > 0) then
                        if (a1%pdg(1) /= a2%pdg(1)) return
                     end if
                  else
                     return
                  end if
               else
                  return
               end if
             end associate
          end do
          flag = .true.
       end if
    end if
  end function pdg_list_eq

  module function pdg_list_lt (pl1, pl2) result (flag)
    class(pdg_list_t), intent(in) :: pl1, pl2
    logical :: flag
    integer :: i
    flag = .false.
    if (allocated (pl1%a) .and. allocated (pl2%a)) then
       if (size (pl1%a) < size (pl2%a)) then
          flag = .true.;  return
       else if (size (pl1%a) > size (pl2%a)) then
          return
       else
          do i = 1, size (pl1%a)
             associate (a1 => pl1%a(i), a2 => pl2%a(i))
               if (allocated (a1%pdg) .and. allocated (a2%pdg)) then
                  if (size (a1%pdg) < size (a2%pdg)) then
                     flag = .true.;  return
                  else if (size (a1%pdg) > size (a2%pdg)) then
                     return
                  else
                     if (size (a1%pdg) > 0) then
                        if (abs (a1%pdg(1)) < abs (a2%pdg(1))) then
                           flag = .true.;  return
                        else if (abs (a1%pdg(1)) > abs (a2%pdg(1))) then
                           return
                        else if (a1%pdg(1) > 0 .and. a2%pdg(1) < 0) then
                           flag = .true.;  return
                        else if (a1%pdg(1) < 0 .and. a2%pdg(1) > 0) then
                           return
                        end if
                     end if
                  end if
               else
                  return
               end if
             end associate
          end do
          flag = .false.
       end if
    end if
  end function pdg_list_lt

  module function pdg_list_replace (pl, i, pl_insert, n_in) result (pl_out)
    type(pdg_list_t) :: pl_out
    class(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: i
    class(pdg_list_t), intent(in) :: pl_insert
    integer, intent(in), optional :: n_in
    integer :: n, n_insert, n_out, k
    n = pl%get_size ()
    n_insert = pl_insert%get_size ()
    n_out = n + n_insert - 1
    call pl_out%init (n_out)
!    if (allocated (pl%a)) then
       do k = 1, i - 1
          pl_out%a(k) = pl%a(k)
       end do
!    end if
    if (present (n_in)) then
       pl_out%a(i) = pl_insert%a(1)
       do k = i + 1, n_in
          pl_out%a(k) = pl%a(k)
       end do
       do k = 1, n_insert - 1
          pl_out%a(n_in+k) = pl_insert%a(1+k)
       end do
       do k = 1, n - n_in
          pl_out%a(n_in+k+n_insert-1) = pl%a(n_in+k)
       end do
    else
!       if (allocated (pl_insert%a)) then
          do k = 1, n_insert
             pl_out%a(i-1+k) = pl_insert%a(k)
          end do
!       end if
!       if (allocated (pl%a)) then
          do k = 1, n - i
             pl_out%a(i+n_insert-1+k) = pl%a(i+k)
          end do
       end if
!    end if
  end function pdg_list_replace

  module function pdg_list_fusion (pl, pl_insert, i, check_if_existing) result (pl_out)
    type(pdg_list_t) :: pl_out
    class(pdg_list_t), intent(in) :: pl
    type(pdg_list_t), intent(in) :: pl_insert
    integer, intent(in) :: i
    logical, intent(in) :: check_if_existing
    integer :: n, n_insert, k, n_out
    logical :: new_pdg
    n = pl%get_size ()
    n_insert = pl_insert%get_size ()
    new_pdg = .not. check_if_existing .or. &
         (.not. any (pl%search_for_particle (pl_insert%a(1)%pdg)))
    call pl_out%init (n + n_insert - 1)
    do k = 1, n
       if (new_pdg .and. k == i) then
          pl_out%a(k) = pl%a(k)%add (pl_insert%a(1))
       else
          pl_out%a(k) = pl%a(k)
       end if
    end do
    do k = n + 1, n + n_insert - 1
       pl_out%a(k) = pl_insert%a(k-n)
    end do
  end function pdg_list_fusion

  module function pdg_list_get_pdg_sizes (pl) result (i_size)
    integer, dimension(:), allocatable :: i_size
    class(pdg_list_t), intent(in) :: pl
    integer :: i, n
    n = pl%get_size ()
    allocate (i_size (n))
    do i = 1, n
       i_size(i) = size (pl%a(i)%pdg)
    end do
  end function pdg_list_get_pdg_sizes

  module subroutine pdg_list_match_replace (pl, pl_match, success)
    class(pdg_list_t), intent(inout) :: pl
    class(pdg_list_t), intent(in) :: pl_match
    logical, intent(out) :: success
    integer :: i, j
    success = .true.
    SCAN_ENTRIES: do i = 1, size (pl%a)
       do j = 1, size (pl_match%a)
          if (pl%a(i) .match. pl_match%a(j)) then
             pl%a(i) = pl_match%a(j)
             cycle SCAN_ENTRIES
          end if
       end do
       success = .false.
       return
    end do SCAN_ENTRIES
  end subroutine pdg_list_match_replace

  module function pdg_list_match_pdg_array (pl, pa) result (flag)
    class(pdg_list_t), intent(in) :: pl
    type(pdg_array_t), intent(in) :: pa
    logical :: flag
    flag = pl%find_match (pa) /= 0
  end function pdg_list_match_pdg_array

  module function pdg_list_find_match_pdg_array (pl, pa, mask) result (i)
    class(pdg_list_t), intent(in) :: pl
    type(pdg_array_t), intent(in) :: pa
    logical, dimension(:), intent(in), optional :: mask
    integer :: i
    do i = 1, size (pl%a)
       if (present (mask)) then
          if (.not. mask(i))  cycle
       end if
       if (pl%a(i) .match. pa)  return
    end do
    i = 0
  end function pdg_list_find_match_pdg_array

  module subroutine pdg_list_create_pdg_array (pl, pdg)
    class(pdg_list_t), intent(in) :: pl
    type(pdg_array_t), dimension(:), intent(inout), allocatable :: pdg
    integer :: n_elements
    integer :: i
    associate (a => pl%a)
      n_elements = size (a)
      if (allocated (pdg))  deallocate (pdg)
      allocate (pdg (n_elements))
      do i = 1, n_elements
         pdg(i) = a(i)
      end do
    end associate
  end subroutine pdg_list_create_pdg_array

  module subroutine pdg_list_create_antiparticles (pl, pl_anti, n_new_particles)
    class(pdg_list_t), intent(in) :: pl
    type(pdg_list_t), intent(out) :: pl_anti
    integer, intent(out) :: n_new_particles
    type(pdg_list_t) :: pl_inverse
    integer :: i, n
    integer :: n_identical
    logical, dimension(:), allocatable :: collect
    n = pl%get_size (); n_identical = 0
    allocate (collect (n)); collect = .true.
    call pl_inverse%init (n)
    do i = 1, n
       pl_inverse%a(i) = pl%a(i)%invert()
    end do
    do i = 1, n
       if (any (pl_inverse%a(i) == pl%a)) then
          collect(i) = .false.
          n_identical = n_identical + 1
       end if
    end do
    n_new_particles = n - n_identical
    if (n_new_particles > 0) then
       call pl_anti%init (n_new_particles)
       do i = 1, n
          if (collect (i)) pl_anti%a(i) = pl_inverse%a(i)
       end do
   end if
  end subroutine pdg_list_create_antiparticles

  elemental module function pdg_list_search_for_particle (pl, i_part) result (found)
    logical :: found
    class(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: i_part
    integer :: i_pl
    do i_pl = 1, size (pl%a)
       found = pl%a(i_pl)%search_for_particle (i_part)
       if (found)  return
    end do
  end function pdg_list_search_for_particle

  module function pdg_list_contains_colored_particles (pl) result (colored)
    class(pdg_list_t), intent(in) :: pl
    logical :: colored
    integer :: i
    colored = .false.
    do i = 1, size (pl%a)
       if (pl%a(i)%has_colored_particles()) then
          colored = .true.
          exit
       end if
    end do
  end function pdg_list_contains_colored_particles


end submodule pdg_arrays_s

