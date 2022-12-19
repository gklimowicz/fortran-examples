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

submodule (fks_regions) fks_regions_s

  use debug_master, only: debug_on
  use format_utils, only: write_separator
  use numeric_utils
  use string_utils, only: str
  use io_units
  use permutations
  use physics_defs
  use flavors
  use pdg_arrays

  implicit none

contains

  pure module subroutine ftuple_assign (ftuple_out, ftuple_in)
    type(ftuple_t), intent(out) :: ftuple_out
    type(ftuple_t), intent(in) :: ftuple_in
    ftuple_out%ireg = ftuple_in%ireg
    ftuple_out%i_res = ftuple_in%i_res
    ftuple_out%splitting_type = ftuple_in%splitting_type
    ftuple_out%pseudo_isr = ftuple_in%pseudo_isr
    ftuple_out%qcd_split = ftuple_in%qcd_split
  end subroutine ftuple_assign

  elemental module function ftuple_equal (f1, f2) result (value)
    logical :: value
    type(ftuple_t), intent(in) :: f1, f2
    value = all (f1%ireg == f2%ireg) .and. f1%i_res == f2%i_res &
         .and. f1%splitting_type == f2%splitting_type &
         .and. (f1%pseudo_isr .eqv. f2%pseudo_isr) &
         .and. (f1%qcd_split .eqv. f2%qcd_split)
  end function ftuple_equal

  elemental function ftuple_equal_ireg (f1, f2) result (value)
    logical :: value
    type(ftuple_t), intent(in) :: f1, f2
    value = all (f1%ireg == f2%ireg)
  end function ftuple_equal_ireg

  elemental module function ftuple_greater (f1, f2) result (greater)
    logical :: greater
    type(ftuple_t), intent(in) :: f1, f2
    if (f1%ireg(1) == f2%ireg(1)) then
       greater = f1%ireg(2) > f2%ireg(2)
    else
       greater = f1%ireg(1) > f2%ireg(1)
    end if
  end function ftuple_greater

  elemental module function ftuple_less (f1, f2) result (less)
    logical :: less
    type(ftuple_t), intent(in) :: f1, f2
    if (f1%ireg(1) == f2%ireg(1)) then
       less = f1%ireg(2) < f2%ireg(2)
    else
       less = f1%ireg(1) < f2%ireg(1)
    end if
  end function ftuple_less

  subroutine ftuple_sort_array (ftuple_array, equivalences)
    type(ftuple_t), intent(inout), dimension(:), allocatable :: ftuple_array
    logical, intent(inout), dimension(:,:), allocatable :: equivalences
    type(ftuple_t) :: ftuple_tmp
    logical, dimension(:), allocatable :: eq_tmp
    integer :: i1, i2, n
    n = size (ftuple_array)
    allocate (eq_tmp (n))
    do i1 = 2, n
       i2 = i1
       do while (ftuple_array(i2 - 1) > ftuple_array(i2))
          ftuple_tmp = ftuple_array(i2 - 1)
          eq_tmp = equivalences(i2, :)
          ftuple_array(i2 - 1) = ftuple_array(i2)
          ftuple_array(i2) = ftuple_tmp
          equivalences(i2 - 1, :) = equivalences(i2, :)
          equivalences(i2, :) = eq_tmp
          i2 = i2 - 1
          if (i2 == 1) exit
       end do
    end do
  end subroutine ftuple_sort_array

  module subroutine ftuple_write (ftuple, unit, newline)
    class(ftuple_t), intent(in) :: ftuple
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: newline
    integer :: u
    logical :: nl
    u = given_output_unit (unit); if (u < 0) return
    nl = .true.; if (present(newline)) nl = newline
    if (all (ftuple%ireg > -1)) then
       if (ftuple%i_res > 0) then
          if (nl) then
             write (u, "(A1,I1,A1,I1,A1,I1,A1)") &
                 '(', ftuple%ireg(1), ',', ftuple%ireg(2), ';', ftuple%i_res, ')'
          else
             write (u, "(A1,I1,A1,I1,A1,I1,A1)", advance = "no") &
                 '(', ftuple%ireg(1), ',', ftuple%ireg(2), ';', ftuple%i_res, ')'
          end if
       else
          if (nl) then
             write (u, "(A1,I1,A1,I1,A1)") &
                  '(', ftuple%ireg(1), ',', ftuple%ireg(2), ')'
          else
             write (u, "(A1,I1,A1,I1,A1)", advance = "no") &
                  '(', ftuple%ireg(1), ',', ftuple%ireg(2), ')'
          end if
       end if
    else
       write (u, "(A)") "(Empty)"
    end if
  end subroutine ftuple_write

  function ftuple_string (ftuples, latex)
    type(string_t) :: ftuple_string
    type(ftuple_t), intent(in), dimension(:) :: ftuples
    logical, intent(in) :: latex
    integer :: i, nreg
    if (latex) then
       ftuple_string = var_str ("$\left\{")
    else
       ftuple_string = var_str ("{")
    end if
    nreg = size(ftuples)
    do i = 1, nreg
       if (ftuples(i)%i_res == 0) then
          ftuple_string = ftuple_string // var_str ("(") // &
               str (ftuples(i)%ireg(1)) // var_str (",") // &
               str (ftuples(i)%ireg(2)) // var_str (")")
       else
          ftuple_string = ftuple_string // var_str ("(") // &
               str (ftuples(i)%ireg(1)) // var_str (",") // &
               str (ftuples(i)%ireg(2)) // var_str (";") // &
               str (ftuples(i)%i_res) // var_str (")")
       end if
       if (ftuples(i)%pseudo_isr) ftuple_string = ftuple_string // var_str ("*")
       if (i < nreg) ftuple_string = ftuple_string // var_str (",")
    end do
    if (latex) then
       ftuple_string = ftuple_string // var_str ("\right\}$")
    else
       ftuple_string = ftuple_string // var_str ("}")
    end if
  end function ftuple_string

  module subroutine ftuple_get (ftuple, pos1, pos2)
    class(ftuple_t), intent(in) :: ftuple
    integer, intent(out) :: pos1, pos2
    pos1 = ftuple%ireg(1)
    pos2 = ftuple%ireg(2)
  end subroutine ftuple_get

  module subroutine ftuple_set (ftuple, pos1, pos2)
    class(ftuple_t), intent(inout) :: ftuple
    integer, intent(in) ::  pos1, pos2
    ftuple%ireg(1) = pos1
    ftuple%ireg(2) = pos2
  end subroutine ftuple_set

  module subroutine ftuple_determine_splitting_type_fsr (ftuple, flv, i, j)
    class(ftuple_t), intent(inout) :: ftuple
    type(flv_structure_t), intent(in) :: flv
    integer, intent(in) :: i, j
    associate (flst => flv%flst)
       if (is_vector (flst(i)) .and. is_vector (flst(j))) then
          ftuple%splitting_type = V_TO_VV
       else if (flst(i)+flst(j) == 0 &
             .and. is_fermion (flst(i))) then
          ftuple%splitting_type = V_TO_FF
       else if (is_fermion(flst(i)) .and. is_massless_vector (flst(j)) &
             .or. is_fermion(flst(j)) .and. is_massless_vector (flst(i))) then
          ftuple%splitting_type = F_TO_FV
       else
          ftuple%splitting_type = UNDEFINED_SPLITTING
       end if
    end associate
  end subroutine ftuple_determine_splitting_type_fsr

  module subroutine ftuple_determine_splitting_type_isr (ftuple, flv, i, j)
    class(ftuple_t), intent(inout) :: ftuple
    type(flv_structure_t), intent(in) :: flv
    integer, intent(in) :: i, j
    integer :: em
    em = i; if (i == 0) em = 1
    associate (flst => flv%flst)
       if (is_vector (flst(em)) .and. is_vector (flst(j))) then
          ftuple%splitting_type = V_TO_VV
       else if (is_massless_vector(flst(em)) .and. is_fermion(flst(j))) then
          ftuple%splitting_type = F_TO_VF
       else if (is_fermion(flst(em)) .and. is_massless_vector(flst(j))) then
          ftuple%splitting_type = F_TO_FV
       else if (is_fermion(flst(em)) .and. is_fermion(flst(j))) then
          ftuple%splitting_type = V_TO_FF
       else
          ftuple%splitting_type = UNDEFINED_SPLITTING
       end if
    end associate
  end subroutine ftuple_determine_splitting_type_isr

  module subroutine ftuple_determine_sub_correction_type &
       (ftuple, flv_born, flv_real, i, j)
    class(ftuple_t), intent(inout) :: ftuple
    type(flv_structure_t), intent(in) :: flv_born, flv_real
    integer, intent(in) :: i, j
    type(flv_structure_t) :: flv_test_qcd
    integer :: em
    em = i; if (i == 0) em = 1
    select case (ftuple%splitting_type)
    case (V_TO_VV)
       ftuple%qcd_split = is_gluon (flv_real%flst(em)) .and. is_gluon (flv_real%flst(j))
    case (F_TO_VF)
       ftuple%qcd_split = is_gluon (flv_real%flst(em))
    case (F_TO_FV)
       if (i == 0) then
          ftuple%qcd_split = is_gluon (flv_real%flst(j))
       else
          ftuple%qcd_split = is_gluon (flv_real%flst(i)) .or. is_gluon (flv_real%flst(j))
       end if
    case (V_TO_FF)
       if (any ([i, j] <= flv_real%n_in)) then
          flv_test_qcd = flv_real%insert_particle_isr (i, j, GLUON)
       else
          flv_test_qcd = flv_real%insert_particle_fsr (i, j, GLUON)
       end if
       ftuple%qcd_split = flv_test_qcd .equiv. flv_born
    case (UNDEFINED_SPLITTING)
       ftuple%qcd_split = .false.
    end select
  end subroutine ftuple_determine_sub_correction_type

  elemental module function ftuple_has_negative_elements &
       (ftuple) result (value)
    logical :: value
    class(ftuple_t), intent(in) :: ftuple
    value = any (ftuple%ireg < 0)
  end function ftuple_has_negative_elements

  elemental module function ftuple_has_identical_elements &
       (ftuple) result (value)
    logical :: value
    class(ftuple_t), intent(in) :: ftuple
    value = ftuple%ireg(1) == ftuple%ireg(2)
  end function ftuple_has_identical_elements

  module subroutine ftuple_list_write (list, unit, verbose)
    class(ftuple_list_t), intent(in), target :: list
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    type(ftuple_list_t), pointer :: current
    logical :: verb
    integer :: u
    u = given_output_unit (unit); if (u < 0) return
    verb = .false.; if (present (verbose)) verb = verbose
    select type (list)
    type is (ftuple_list_t)
       current => list
       do
          call current%ftuple%write (unit = u, newline = .false.)
          if (verb .and. associated (current%equiv))  &
               write (u, '(A)', advance = "no") "'"
          if (associated (current%next)) then
             current => current%next
          else
             exit
          end if
       end do
       write (u, *) ""
    end select
  end subroutine ftuple_list_write

  module subroutine ftuple_list_append (list, ftuple)
    class(ftuple_list_t), intent(inout), target :: list
    type(ftuple_t), intent(in) :: ftuple
    type(ftuple_list_t), pointer :: current

    select type (list)
    type is (ftuple_list_t)
       if (list%index == 0) then
          nullify (list%next)
          list%index = 1
          list%ftuple = ftuple
       else
          current => list
          do
             if (associated (current%next)) then
                current => current%next
             else
                allocate (current%next)
                nullify (current%next%next)
                nullify (current%next%equiv)
                current%next%prev => current
                current%next%index = current%index + 1
                current%next%ftuple = ftuple
                exit
             end if
          end do
       end if
    end select
  end subroutine ftuple_list_append

  impure elemental module function ftuple_list_get_n_tuples &
       (list) result(n_tuples)
    integer :: n_tuples
    class(ftuple_list_t), intent(in), target :: list
    type(ftuple_list_t), pointer :: current
    n_tuples = 0
    select type (list)
    type is (ftuple_list_t)
       current => list
       if (current%index > 0) then
          n_tuples = 1
          do
             if (associated (current%next)) then
                current => current%next
                n_tuples = n_tuples + 1
             else
                exit
             end if
          end do
       end if
     end select
  end function ftuple_list_get_n_tuples

  module function ftuple_list_get_entry (list, index) result (entry)
    type(ftuple_list_t), pointer :: entry
    class(ftuple_list_t), intent(in), target :: list
    integer, intent(in) :: index
    type(ftuple_list_t), pointer :: current
    integer :: i
    entry => null()
    select type (list)
    type is (ftuple_list_t)
       current => list
       if (index == 1) then
          entry => current
       else
          do i = 1, index - 1
             current => current%next
          end do
          entry => current
       end if
    end select
  end function ftuple_list_get_entry

  module function ftuple_list_get_ftuple (list, index)  result (ftuple)
    type(ftuple_t) :: ftuple
    class(ftuple_list_t), intent(in), target :: list
    integer, intent(in) :: index
    type(ftuple_list_t), pointer :: entry
    entry => list%get_entry (index)
    ftuple = entry%ftuple
  end function ftuple_list_get_ftuple

  module subroutine ftuple_list_set_equiv (list, i1, i2)
    class(ftuple_list_t), intent(in) :: list
    integer, intent(in) :: i1, i2
    type(ftuple_list_t), pointer :: list1, list2 => null ()
    select type (list)
    type is (ftuple_list_t)
       if (list%get_ftuple (i1) > list%get_ftuple (i2)) then
          list1 => list%get_entry (i2)
          list2 => list%get_entry (i1)
       else
          list1 => list%get_entry (i1)
          list2 => list%get_entry (i2)
       end if
       do
          if (associated (list1%equiv)) then
             list1 => list1%equiv
          else
             exit
          end if
       end do
       list1%equiv => list2
    end select
  end subroutine ftuple_list_set_equiv

  module function ftuple_list_check_equiv(list, i1, i2) result (eq)
    class(ftuple_list_t), intent(in) :: list
    integer, intent(in) :: i1, i2
    logical :: eq
    type(ftuple_list_t), pointer :: current
    eq = .false.
    select type (list)
    type is (ftuple_list_t)
       current => list%get_entry (i1)
       do
          if (associated (current%equiv)) then
             current => current%equiv
             if (current%index == i2) then
                eq = .true.
                exit
             end if
          else
             exit
          end if
       end do
    end select
  end function ftuple_list_check_equiv

  module subroutine ftuple_list_to_array &
       (ftuple_list, ftuple_array, equivalences, ordered)
    class(ftuple_list_t), intent(in), target :: ftuple_list
    type(ftuple_t), intent(out), dimension(:), allocatable :: ftuple_array
    logical, intent(out), dimension(:,:), allocatable :: equivalences
    logical, intent(in) :: ordered
    integer :: i_tuple, n
    type(ftuple_list_t), pointer :: current => null ()
    integer :: i1, i2
    type(ftuple_t) :: ftuple_tmp
    logical, dimension(:), allocatable :: eq_tmp
    n = ftuple_list%get_n_tuples ()
    allocate (ftuple_array (n), equivalences (n, n))
    equivalences = .false.
    select type (ftuple_list)
    type is (ftuple_list_t)
       current => ftuple_list
       i_tuple = 1
       do
          ftuple_array(i_tuple) = current%ftuple
          if (associated (current%equiv)) then
             i1 = current%index
             i2 = current%equiv%index
             equivalences (i1, i2) = .true.
          end if
          if (associated (current%next)) then
             current => current%next
             i_tuple = i_tuple + 1
          else
             exit
          end if
       end do
    end select
    if (ordered) call ftuple_sort_array (ftuple_array, equivalences)
  end subroutine ftuple_list_to_array

  subroutine print_equivalence_matrix (ftuple_array, equivalences)
    type(ftuple_t), intent(in), dimension(:) :: ftuple_array
    logical, intent(in), dimension(:,:) :: equivalences
    integer :: i, i1, i2
    print *, 'Equivalence matrix: '
    do i = 1, size (ftuple_array)
       call ftuple_array(i)%get(i1,i2)
       print *, 'i: ', i, '(', i1, i2, '): ', equivalences(i,:)
    end do
  end subroutine print_equivalence_matrix

  module function flv_structure_valid_pair &
     (flv, i, j, flv_ref, model) result (valid)
    logical :: valid
    class(flv_structure_t), intent(in) :: flv
    integer, intent(in) :: i,j
    type(flv_structure_t), intent(in) :: flv_ref
    type(model_t), intent(in) :: model
    integer :: k, n_orig
    type(flv_structure_t) :: flv_test
    integer, dimension(:), allocatable :: flv_orig
    valid = .false.
    if (all ([i, j] <= flv%n_in)) return
    if (i <= flv%n_in .and. is_fermion(flv%flst(i))) then
       call model%match_vertex (-flv%flst(i), flv%flst(j), flv_orig)
    else if (j <= flv%n_in .and. is_fermion(flv%flst(j))) then
       call model%match_vertex (flv%flst(i), -flv%flst(j), flv_orig)
    else
       call model%match_vertex (flv%flst(i), flv%flst(j), flv_orig)
    end if
    n_orig = size (flv_orig)
    if (n_orig == 0) then
       return
    else
      do k = 1, n_orig
         if (any ([i, j] <= flv%n_in)) then
            flv_test = flv%insert_particle_isr (i, j, flv_orig(k))
         else
            flv_test = flv%insert_particle_fsr (i, j, flv_orig(k))
         end if
         valid = flv_ref .equiv. flv_test
         call flv_test%final ()
         if (valid) return
      end do
    end if
    deallocate (flv_orig)
  end function flv_structure_valid_pair

  function flv_structure_equivalent (flv1, flv2, with_tag) result (equiv)
    logical :: equiv
    type(flv_structure_t), intent(in) :: flv1, flv2
    logical, intent(in) :: with_tag
    type(flv_perm_t) :: perm
    integer :: n
    n = size (flv1%flst)
    equiv = .true.
    if (n /= size (flv2%flst)) then
       call msg_fatal &
            ('flv_structure_equivalent: flavor arrays do not have equal lengths')
    else if (flv1%n_in /= flv2%n_in) then
       call msg_fatal &
            ('flv_structure_equivalent: flavor arrays do not have equal n_in')
    else
       call perm%init (flv1, flv2, flv1%n_in, flv1%nlegs, with_tag)
       equiv = perm%eqv (flv2, flv1, with_tag)
       call perm%final ()
    end if
  end function flv_structure_equivalent

  module function flv_structure_equivalent_no_tag (flv1, flv2) result(equiv)
    logical :: equiv
    type(flv_structure_t), intent(in) :: flv1, flv2
    equiv = flv_structure_equivalent (flv1, flv2, .false.)
  end function flv_structure_equivalent_no_tag

  module function flv_structure_equivalent_with_tag (flv1, flv2) result(equiv)
    logical :: equiv
    type(flv_structure_t), intent(in) :: flv1, flv2
    equiv = flv_structure_equivalent (flv1, flv2, .true.)
  end function flv_structure_equivalent_with_tag


  pure module subroutine flv_structure_assign_flv (flv_out, flv_in)
    type(flv_structure_t), intent(out) :: flv_out
    type(flv_structure_t), intent(in) :: flv_in
    flv_out%nlegs = flv_in%nlegs
    flv_out%n_in = flv_in%n_in
    flv_out%prt_symm_fs = flv_in%prt_symm_fs
    if (allocated (flv_in%flst)) then
       allocate (flv_out%flst (size (flv_in%flst)))
       flv_out%flst = flv_in%flst
    end if
    if (allocated (flv_in%tag)) then
       allocate (flv_out%tag (size (flv_in%tag)))
       flv_out%tag = flv_in%tag
    end if
    if (allocated (flv_in%massive)) then
       allocate (flv_out%massive (size (flv_in%massive)))
       flv_out%massive = flv_in%massive
    end if
    if (allocated (flv_in%colored)) then
       allocate (flv_out%colored (size (flv_in%colored)))
       flv_out%colored = flv_in%colored
    end if
  end subroutine flv_structure_assign_flv

  pure module subroutine flv_structure_assign_integer (flv_out, iarray)
    type(flv_structure_t), intent(out) :: flv_out
    integer, intent(in), dimension(:) :: iarray
    integer :: i
    flv_out%nlegs = size (iarray)
    allocate (flv_out%flst (flv_out%nlegs))
    allocate (flv_out%tag (flv_out%nlegs))
    flv_out%flst = iarray
    flv_out%tag = [(i, i = 1, flv_out%nlegs)]
  end subroutine flv_structure_assign_integer

  module function flv_structure_remove_particle (flv, index) result (flv_new)
    type(flv_structure_t) :: flv_new
    class(flv_structure_t), intent(in) :: flv
    integer, intent(in) :: index
    integer :: n1, n2
    integer :: i, removed_tag
    n1 = size (flv%flst); n2 = n1 - 1
    allocate (flv_new%flst (n2), flv_new%tag (n2))
    flv_new%nlegs = n2
    flv_new%n_in = flv%n_in
    removed_tag = flv%tag(index)
    if (index == 1) then
       flv_new%flst(1 : n2) = flv%flst(2 : n1)
       flv_new%tag(1 : n2) = flv%tag(2 : n1)
    else if (index == n1) then
       flv_new%flst(1 : n2) = flv%flst(1 : n2)
       flv_new%tag(1 : n2) = flv%tag(1 : n2)
    else
       flv_new%flst(1 : index - 1) = flv%flst(1 : index - 1)
       flv_new%flst(index : n2) = flv%flst(index + 1 : n1)
       flv_new%tag(1 : index - 1) = flv%tag(1 : index - 1)
       flv_new%tag(index : n2) = flv%tag(index + 1 : n1)
    end if
    do i = 1, n2
       if (flv_new%tag(i) > removed_tag) &
            flv_new%tag(i) = flv_new%tag(i) - 1
    end do
    call flv_new%compute_prt_symm_fs (flv_new%n_in)
  end function flv_structure_remove_particle

  module function flv_structure_insert_particle_fsr &
       (flv, i1, i2, flv_add) result (flv_new)
    type(flv_structure_t) :: flv_new
    class(flv_structure_t), intent(in) :: flv
    integer, intent(in) :: i1, i2, flv_add
    if (flv%flst(i1) + flv_add == 0 .or. flv%flst(i2) + flv_add == 0) then
       flv_new = flv%insert_particle (i1, i2, -flv_add)
    else
       flv_new = flv%insert_particle (i1, i2, flv_add)
    end if
  end function flv_structure_insert_particle_fsr

  module function flv_structure_insert_particle_isr &
       (flv, i_in, i_out, flv_add) result (flv_new)
    type(flv_structure_t) :: flv_new
    class(flv_structure_t), intent(in) :: flv
    integer, intent(in) :: i_in, i_out, flv_add
    if (flv%flst(i_in) + flv_add == 0) then
       flv_new = flv%insert_particle (i_in, i_out, -flv_add)
    else
       flv_new = flv%insert_particle (i_in, i_out, flv_add)
    end if
  end function flv_structure_insert_particle_isr

  module function flv_structure_insert_particle &
       (flv, i1, i2, particle) result (flv_new)
    type(flv_structure_t) :: flv_new
    class(flv_structure_t), intent(in) :: flv
    integer, intent(in) :: i1, i2, particle
    type(flv_structure_t) :: flv_tmp
    integer :: n1, n2
    integer :: new_tag
    n1 = size (flv%flst); n2 = n1 - 1
    allocate (flv_new%flst (n2), flv_new%tag (n2))
    flv_new%nlegs = n2
    flv_new%n_in = flv%n_in
    new_tag = maxval(flv%tag) + 1
    if (i1 < i2) then
       flv_tmp = flv%remove_particle (i1)
       flv_tmp = flv_tmp%remove_particle (i2 - 1)
    else if(i2 < i1) then
       flv_tmp = flv%remove_particle(i2)
       flv_tmp = flv_tmp%remove_particle(i1 - 1)
    else
       call msg_fatal ("flv_structure_insert_particle: Indices are identical!")
    end if
    if (i1 == 1) then
       flv_new%flst(1) = particle
       flv_new%flst(2 : n2) = flv_tmp%flst(1 : n2 - 1)
       flv_new%tag(1) = new_tag
       flv_new%tag(2 : n2) = flv_tmp%tag(1 : n2 - 1)
    else if (i1 == n1 .or. i1 == n2) then
       flv_new%flst(1 : n2 - 1) = flv_tmp%flst(1 : n2 - 1)
       flv_new%flst(n2) = particle
       flv_new%tag(1 : n2 - 1) = flv_tmp%tag(1 : n2 - 1)
       flv_new%tag(n2) = new_tag
    else
       flv_new%flst(1 : i1 - 1) = flv_tmp%flst(1 : i1 - 1)
       flv_new%flst(i1) = particle
       flv_new%flst(i1 + 1 : n2) = flv_tmp%flst(i1 : n2 - 1)
       flv_new%tag(1 : i1 - 1) = flv_tmp%tag(1 : i1 - 1)
       flv_new%tag(i1) = new_tag
       flv_new%tag(i1 + 1 : n2) = flv_tmp%tag(i1 : n2 - 1)
    end if
    call flv_new%compute_prt_symm_fs (flv_new%n_in)
  end function flv_structure_insert_particle

  module function flv_structure_count_particle (flv, part) result (n)
    class(flv_structure_t), intent(in) :: flv
    integer, intent(in) :: part
    integer :: n
    n = count (flv%flst == part)
  end function flv_structure_count_particle

  module subroutine flv_structure_init (flv, aval, n_in, tags)
    class(flv_structure_t), intent(inout) :: flv
    integer, intent(in), dimension(:) :: aval
    integer, intent(in) :: n_in
    integer, intent(in), dimension(:), optional :: tags
    integer :: i, n
    integer, dimension(:), allocatable :: aval_unique
    integer, dimension(:), allocatable :: mult
    n = size (aval)
    allocate (flv%flst (n), flv%tag (n))
    flv%flst = aval
    if (present (tags)) then
       flv%tag = tags
    else
       do i = 1, n
          flv%tag(i) = i
       end do
    end if
    flv%nlegs = n
    flv%n_in = n_in
    call flv%compute_prt_symm_fs (flv%n_in)
  end subroutine flv_structure_init

  module subroutine flv_structure_compute_prt_symm_fs (flv, n_in)
    class(flv_structure_t), intent(inout) :: flv
    integer, intent(in) :: n_in
    integer, dimension(:), allocatable :: flst_unique
    integer, dimension(:), allocatable :: mult
    integer :: i
    flst_unique = remove_duplicates_from_int_array (flv%flst(n_in + 1 :))
    allocate (mult(size (flst_unique)))
    do i = 1, size (flst_unique)
       mult(i) = count (flv%flst(n_in + 1 :) == flst_unique(i))
    end do
    flv%prt_symm_fs = one / product (gamma (real (mult + 1, default)))
  end subroutine flv_structure_compute_prt_symm_fs

  module subroutine flv_structure_write (flv, unit)
    class(flv_structure_t), intent(in) :: flv
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit); if (u < 0) return
    write (u, '(A)') char (flv%to_string ())
  end subroutine flv_structure_write

  module function flv_structure_to_string (flv) result (flv_string)
    type(string_t) :: flv_string
    class(flv_structure_t), intent(in) :: flv
    integer :: i, n
    if (allocated (flv%flst)) then
       flv_string = var_str ("[")
       n = size (flv%flst)
       do i = 1, n - 1
          flv_string = flv_string // str (flv%flst(i)) // var_str(",")
       end do
       flv_string = flv_string // str (flv%flst(n)) // var_str("]")
    else
       flv_string = var_str ("[not allocated]")
    end if
  end function flv_structure_to_string

  module function flv_structure_create_uborn &
       (flv, emitter, nlo_correction_type) result(flv_uborn)
    type(flv_structure_t) :: flv_uborn
    class(flv_structure_t), intent(in) :: flv
    type(string_t), intent(in) :: nlo_correction_type
    integer, intent(in) :: emitter
    integer n_legs
    integer :: f1, f2
    integer :: gauge_boson

    n_legs = size(flv%flst)
    allocate (flv_uborn%flst (n_legs - 1), flv_uborn%tag (n_legs - 1))
    gauge_boson = determine_gauge_boson_to_be_inserted ()

    if (emitter > flv%n_in) then
       f1 = flv%flst(n_legs); f2 = flv%flst(n_legs - 1)
       if (is_massless_vector (f1)) then
          !!! Emitted particle is a gluon or photon => just remove it
          flv_uborn = flv%remove_particle(n_legs)
       else if (is_fermion (f1) .and. is_fermion (f2) .and. f1 + f2 == 0) then
          !!! Emission type is a gauge boson splitting into two fermions
          flv_uborn = flv%insert_particle(n_legs - 1, n_legs, gauge_boson)
       else
          call msg_error ("Create underlying Born: Unsupported splitting type.")
          call msg_error (char (str (flv%flst)))
          call msg_fatal ("FKS - FAIL")
       end if
    else if (emitter > 0) then
       f1 = flv%flst(n_legs); f2 = flv%flst(emitter)
       if (is_massless_vector (f1)) then
          flv_uborn = flv%remove_particle(n_legs)
       else if (is_fermion (f1) .and. is_massless_vector (f2)) then
          flv_uborn = flv%insert_particle (emitter, n_legs, -f1)
       else if (is_fermion (f1) .and. is_fermion (f2) .and. f1 == f2) then
          flv_uborn = flv%insert_particle(emitter, n_legs, gauge_boson)
       end if
    else
       flv_uborn = flv%remove_particle (n_legs)
    end if

  contains
    integer function determine_gauge_boson_to_be_inserted ()
      select case (char (nlo_correction_type))
      case ("QCD")
         determine_gauge_boson_to_be_inserted = GLUON
      case ("EW")
         determine_gauge_boson_to_be_inserted = PHOTON
      case ("Full")
         call msg_fatal ("NLO correction type 'Full' not yet implemented!")
      case default
         call msg_fatal ("Invalid NLO correction type! Valid inputs " // &
              "are: QCD, EW and Full (default: QCD)")
      end select
    end function determine_gauge_boson_to_be_inserted

  end function flv_structure_create_uborn

  module subroutine flv_structure_init_mass_color_and_charge (flv, model)
    class(flv_structure_t), intent(inout) :: flv
    type(model_t), intent(in) :: model
    integer :: i
    type(flavor_t) :: flavor
    allocate (flv%massive (flv%nlegs), flv%colored(flv%nlegs), &
         flv%charge(flv%nlegs))
    do i = 1, flv%nlegs
       call flavor%init (flv%flst(i), model)
       flv%massive(i) = flavor%get_mass () > 0
       flv%colored(i) = &
            is_quark (flv%flst(i)) .or. is_gluon (flv%flst(i))
       flv%charge(i) = flavor%get_charge ()
    end do
  end subroutine flv_structure_init_mass_color_and_charge

  module function flv_structure_get_last_two (flv, n) result (flst_last)
    integer, dimension(2) :: flst_last
    class(flv_structure_t), intent(in) :: flv
    integer, intent(in) :: n
    flst_last = [flv%flst(n - 1), flv%flst(n)]
  end function flv_structure_get_last_two

  module subroutine flv_structure_final (flv)
    class(flv_structure_t), intent(inout) :: flv
    if (allocated (flv%flst)) deallocate (flv%flst)
    if (allocated (flv%tag)) deallocate (flv%tag)
    if (allocated (flv%massive)) deallocate (flv%massive)
    if (allocated (flv%colored)) deallocate (flv%colored)
    if (allocated (flv%charge)) deallocate (flv%charge)
  end subroutine flv_structure_final

  module subroutine flv_perm_init &
       (perm, flv_in, flv_ref, n_first, n_last, with_tag)
    class(flv_perm_t), intent(out) :: perm
    type(flv_structure_t), intent(in) :: flv_in, flv_ref
    integer, intent(in) :: n_first, n_last
    logical, intent(in) :: with_tag
    integer :: flv1, flv2, tmp
    integer :: tag1, tag2
    integer :: i, j, j_min, i_perm
    integer, dimension(:,:), allocatable :: perm_list_tmp
    type(flv_structure_t) :: flv_copy
    logical :: condition
    logical, dimension(:), allocatable :: already_correct
    flv_copy = flv_in
    allocate (perm_list_tmp (factorial (n_last - n_first - 1), 2))
    allocate (already_correct (flv_in%nlegs))
    already_correct = flv_in%flst == flv_ref%flst
    if (with_tag) &
         already_correct = already_correct .and. (flv_in%tag == flv_ref%tag)
    j_min = n_first + 1
    i_perm = 0
    do i = n_first + 1, n_last
       flv1 = flv_ref%flst(i)
       tag1 = flv_ref%tag(i)
       do j = j_min, n_last
          if (already_correct(i) .or. already_correct(j)) cycle
          flv2 = flv_copy%flst(j)
          tag2 = flv_copy%tag(j)
          condition = (flv1 == flv2) .and. i /= j
          if (with_tag) condition = condition .and. (tag1 == tag2)
          if (condition) then
             i_perm = i_perm + 1
             tmp = flv_copy%flst(i)
             flv_copy%flst(i) = flv2
             flv_copy%flst(j) = tmp
             tmp = flv_copy%tag(i)
             flv_copy%tag(i) = tag2
             flv_copy%tag(j) = tmp
             perm_list_tmp (i_perm, 1) = i
             perm_list_tmp (i_perm, 2) = j
             exit
          end if
       end do
       j_min = j_min + 1
    end do
    allocate (perm%perms (i_perm, 2))
    perm%perms = perm_list_tmp (1 : i_perm, :)
    deallocate (perm_list_tmp)
    call flv_copy%final ()
  end subroutine flv_perm_init

  module subroutine flv_perm_write (perm, unit)
    class(flv_perm_t), intent(in) :: perm
    integer, intent(in), optional :: unit
    integer :: i, n, u
    u = given_output_unit (unit); if (u < 0) return
    write (u, "(A)") "Flavor permutation list: "
    n = size (perm%perms, dim = 1)
    if (n > 0) then
       do i = 1, n
          write (u, "(A1,I1,1X,I1,A1)", advance = "no") "[", &
               perm%perms(i,1), perm%perms(i,2), "]"
          if (i < n) write (u, "(A4)", advance = "no") " // "
       end do
       write (u, "(A)") ""
    else
       write (u, "(A)") "[Empty]"
    end if
  end subroutine flv_perm_write

  module subroutine flv_perm_final (perm)
    class(flv_perm_t), intent(inout) :: perm
    if (allocated (perm%perms)) deallocate (perm%perms)
  end subroutine flv_perm_final

  elemental module function flv_perm_apply_flv_structure &
       (perm, flv_in, invert) result (flv_out)
    type(flv_structure_t) :: flv_out
    class(flv_perm_t), intent(in) :: perm
    type(flv_structure_t), intent(in) :: flv_in
    logical, intent(in), optional :: invert
    integer :: i, i1, i2
    integer :: p1, p2, incr
    integer :: flv_tmp, tag_tmp
    logical :: inv
    inv = .false.; if (present(invert)) inv = invert
    flv_out = flv_in
    if (inv) then
       p1 = 1
       p2 = size (perm%perms, dim = 1)
       incr = 1
    else
       p1 = size (perm%perms, dim = 1)
       p2 = 1
       incr = -1
    end if
    do i = p1, p2, incr
       i1 = perm%perms(i,1)
       i2 = perm%perms(i,2)
       flv_tmp = flv_out%flst(i1)
       tag_tmp = flv_out%tag(i1)
       flv_out%flst(i1) = flv_out%flst(i2)
       flv_out%flst(i2) = flv_tmp
       flv_out%tag(i1) = flv_out%tag(i2)
       flv_out%tag(i2) = tag_tmp
    end do
  end function flv_perm_apply_flv_structure

  elemental module function flv_perm_apply_index &
       (perm, i_in) result (i_out)
    integer :: i_out
    class(flv_perm_t), intent(in) :: perm
    integer, intent(in) :: i_in
    integer :: i, i1, i2
    i_out = i_in
    do i = size (perm%perms(:,1)), 1, -1
       i1 = perm%perms(i,1)
       i2 = perm%perms(i,2)
       if (i_out == i1) then
          i_out = i2
       else if (i_out == i2) then
          i_out = i1
       end if
    end do
  end function flv_perm_apply_index

  elemental module function flv_perm_apply_ftuple &
       (perm, f_in) result (f_out)
    type(ftuple_t) :: f_out
    class(flv_perm_t), intent(in) :: perm
    type(ftuple_t), intent(in) :: f_in
    integer :: i, i1, i2
    f_out = f_in
    do i = size (perm%perms, dim = 1), 1, -1
       i1 = perm%perms(i,1)
       i2 = perm%perms(i,2)
       if (f_out%ireg(1) == i1) then
          f_out%ireg(1) = i2
       else if (f_out%ireg(1) == i2) then
          f_out%ireg(1) = i1
       end if
       if (f_out%ireg(2) == i1) then
          f_out%ireg(2) = i2
       else if (f_out%ireg(2) == i2) then
          f_out%ireg(2) = i1
       end if
    end do
    if (f_out%ireg(1) > f_out%ireg(2)) f_out%ireg = f_out%ireg([2,1])
  end function flv_perm_apply_ftuple

  module function flv_perm_eqv &
       (perm, flv1, flv2, with_tag) result (valid)
    logical :: valid
    class(flv_perm_t), intent(in) :: perm
    type(flv_structure_t), intent(in) :: flv1, flv2
    logical, intent(in) :: with_tag
    type(flv_structure_t) :: flv_tmp
    flv_tmp = perm%apply (flv2, invert = .true.)
    valid = all (flv_tmp%flst == flv1%flst)
    if (with_tag) valid = valid .and. all (flv_tmp%tag == flv1%tag)
    call flv_tmp%final ()
  end function flv_perm_eqv

  module subroutine singular_region_init (sregion, alr, mult, i_res, &
       flst_real, flst_uborn, flv_born, emitter, ftuples, equivalences, &
       nlo_correction_type)
    class(singular_region_t), intent(out) :: sregion
    integer, intent(in) :: alr, mult, i_res
    type(flv_structure_t), intent(in) :: flst_real
    type(flv_structure_t), intent(in) :: flst_uborn
    type(flv_structure_t), dimension(:), intent(in) :: flv_born
    integer, intent(in) :: emitter
    type(ftuple_t), intent(inout), dimension(:) :: ftuples
    logical, intent(inout), dimension(:,:) :: equivalences
    type(string_t), intent(in) :: nlo_correction_type
    integer :: i
    call debug_input_values ()
    sregion%alr = alr
    sregion%mult = mult
    sregion%i_res = i_res
    sregion%flst_real = flst_real
    sregion%flst_uborn = flst_uborn
    sregion%emitter = emitter
    sregion%nlo_correction_type = nlo_correction_type
    sregion%nregions = size (ftuples)
    allocate (sregion%ftuples (sregion%nregions))
    sregion%ftuples = ftuples
    do i = 1, size(flv_born)
       if (flv_born (i) .equiv. sregion%flst_uborn) then
          sregion%uborn_index = i
          exit
       end if
    end do
    sregion%sc_required = any (sregion%flst_uborn%flst == GLUON) .or. &
            any (sregion%flst_uborn%flst == PHOTON)
  contains
    subroutine debug_input_values()
      if (debug_on) call msg_debug2 (D_SUBTRACTION, "singular_region_init")
      if (debug2_active (D_SUBTRACTION)) then
         print *, 'alr =    ', alr
         print *, 'mult =    ', mult
         print *, 'i_res =    ', i_res
         call flst_real%write ()
         call flst_uborn%write ()
         print *, 'emitter =    ', emitter
         call print_equivalence_matrix (ftuples, equivalences)
      end if
    end subroutine debug_input_values
  end subroutine singular_region_init

  module subroutine singular_region_write (sregion, unit, maxnregions)
    class(singular_region_t), intent(in) :: sregion
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: maxnregions
    character(len=7), parameter :: flst_format = "(I3,A1)"
    character(len=7), parameter :: ireg_space_format = "(7X,A1)"
    integer :: nreal, nborn, i, u, mr
    integer :: nleft, nright, nreg, nreg_diff
    u = given_output_unit (unit); if (u < 0) return
    mr = sregion%nregions; if (present (maxnregions))  mr = maxnregions
    nreal = size (sregion%flst_real%flst)
    nborn = size (sregion%flst_uborn%flst)
    call write_vline (u)
    write (u, '(A1)', advance = 'no') '['
    do i = 1, nreal - 1
       write (u, flst_format, advance = 'no') sregion%flst_real%flst(i), ','
    end do
    write (u, flst_format, advance = 'no') sregion%flst_real%flst(nreal), ']'
    call write_vline (u)
    write (u, '(I6)', advance = 'no') sregion%real_index
    call write_vline (u)
    write (u, '(I3)', advance = 'no') sregion%emitter
    call write_vline (u)
    write (u, '(I3)', advance = 'no') sregion%mult
    call write_vline (u)
    write (u, '(I4)', advance = 'no') sregion%nregions
    call write_vline (u)
    if (sregion%i_res > 0) then
       write (u, '(I3)', advance = 'no') sregion%i_res
       call write_vline (u)
    end if
    nreg = sregion%nregions
    if (nreg == mr) then
       nleft = 0
       nright = 0
    else
       nreg_diff = mr - nreg
       nleft = nreg_diff / 2
       if (mod(nreg_diff , 2) == 0) then
          nright = nleft
       else
          nright = nleft + 1
       end if
    end if
    if (nleft > 0) then
       do i = 1, nleft
          write(u, ireg_space_format, advance='no') ' '
       end do
    end if
    write (u, '(A)', advance = 'no') char (ftuple_string (sregion%ftuples, .false.))
    call write_vline (u)
    write (u,'(A1)',advance = 'no') '['
    do i = 1, nborn - 1
       write(u, flst_format, advance = 'no') sregion%flst_uborn%flst(i), ','
    end do
    write (u, flst_format, advance = 'no') sregion%flst_uborn%flst(nborn), ']'
    call write_vline (u)
    write (u, '(I7)', advance = 'no') sregion%uborn_index
    call write_vline (u)
    if (sregion%nlo_correction_type == "EW") then
       write (u, '(A3)', advance = 'no') 'ew'
    else if (sregion%nlo_correction_type == "QCD") then
       write (u, '(A4)', advance = 'no') 'qcd'
    else
       write (u, '(A5)', advance = 'no') 'none'
    end if
    write (u, '(A)')
  end subroutine singular_region_write

  module subroutine singular_region_write_latex (region, unit)
    class(singular_region_t), intent(in) :: region
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit); if (u < 0) return
    write (u, "(I2,A3,A,A3,I2,A3,I1,A3,I1,A3,A,A3,I2,A3,A,A3)") &
           region%alr, " & ", char (region%flst_real%to_string ()), &
           " & ", region%real_index, " & ", region%emitter, " & ", &
           region%mult, " & ", char (ftuple_string (region%ftuples, .true.)), &
           " & ", region%uborn_index, " & ",  char (region%flst_uborn%to_string ()), &
           " \\"
  end subroutine singular_region_write_latex

  module subroutine singular_region_set_splitting_info (region, n_in)
    class(singular_region_t), intent(inout) :: region
    integer, intent(in) :: n_in
    integer :: i1, i2
    integer :: reg
    region%double_fsr = .false.
    region%soft_divergence = .false.
    associate (ftuple => region%ftuples)
       do reg = 1, region%nregions
          call ftuple(reg)%get (i1, i2)
          if (i1 /= region%emitter .or. i2 /= region%flst_real%nlegs) then
             cycle
          else
             if (ftuple(reg)%splitting_type == V_TO_VV .or.  &
                 ftuple(reg)%splitting_type == F_TO_FV ) then
                region%soft_divergence = .true.
             end if

             if (i1 == 0) then
                region%coll_divergence = &
                     .not. all (region%flst_real%massive(1:n_in))
             else
                region%coll_divergence = .not. region%flst_real%massive(i1)
             end if

             if (ftuple(reg)%splitting_type == V_TO_VV) then
                if (all (ftuple(reg)%ireg > n_in))  &
                     region%double_fsr = &
                     all (is_gluon (region%flst_real%flst(ftuple(reg)%ireg)))
                exit
             else if (ftuple(reg)%splitting_type == UNDEFINED_SPLITTING) then
                call msg_fatal ("All splittings should be defined!")
             end if
          end if
       end do
       if (.not. region%soft_divergence .and. .not. region%coll_divergence) &
            call msg_fatal ("Singular region defined without divergence!")
     end associate
  end subroutine singular_region_set_splitting_info

  module function singular_region_double_fsr_factor (region, p) result (val)
    class(singular_region_t), intent(in) :: region
    type(vector4_t), intent(in), dimension(:) :: p
    real(default) :: val
    real(default) :: E_rad, E_em
    if (region%double_fsr) then
       E_em = energy (p(region%emitter))
       E_rad = energy (p(region%flst_real%nlegs))
       val = two * E_em / (E_em + E_rad)
    else
       val = one
    end if
  end function singular_region_double_fsr_factor

  module function singular_region_has_soft_divergence (region) result (div)
    logical :: div
    class(singular_region_t), intent(in) :: region
    div = region%soft_divergence
  end function singular_region_has_soft_divergence

  module function singular_region_has_collinear_divergence &
       (region) result (div)
    logical :: div
    class(singular_region_t), intent(in) :: region
    div = region%coll_divergence
  end function singular_region_has_collinear_divergence

  elemental module function singular_region_has_identical_ftuples &
       (sregion) result (value)
    logical :: value
    class(singular_region_t), intent(in) :: sregion
    integer :: alr
    value = .false.
    do alr = 1, sregion%nregions
       value = value .or. (count (sregion%ftuples(alr) == sregion%ftuples) > 1)
    end do
  end function singular_region_has_identical_ftuples

  module subroutine singular_region_assign (reg_out, reg_in)
    type(singular_region_t), intent(out) :: reg_out
    type(singular_region_t), intent(in) :: reg_in
    reg_out%alr = reg_in%alr
    reg_out%i_res = reg_in%i_res
    reg_out%flst_real = reg_in%flst_real
    reg_out%flst_uborn = reg_in%flst_uborn
    reg_out%mult = reg_in%mult
    reg_out%emitter = reg_in%emitter
    reg_out%nregions = reg_in%nregions
    reg_out%real_index = reg_in%real_index
    reg_out%uborn_index = reg_in%uborn_index
    reg_out%double_fsr = reg_in%double_fsr
    reg_out%soft_divergence = reg_in%soft_divergence
    reg_out%coll_divergence = reg_in%coll_divergence
    reg_out%nlo_correction_type = reg_in%nlo_correction_type
    if (allocated (reg_in%ftuples)) then
       allocate (reg_out%ftuples (size (reg_in%ftuples)))
       reg_out%ftuples = reg_in%ftuples
    else
       call msg_bug ("singular_region_assign: Trying to copy a " // &
            "singular region without allocated ftuples!")
    end if
  end subroutine singular_region_assign
  module function singular_region_match (reg1, reg2) result (match)
    type(singular_region_t), intent(in) :: reg1, reg2
    logical :: match
    match = all ([reg1%flst_real%eqv_index, reg2%flst_real%eqv_index] > 0)
    match = match .and. (reg1%flst_real%eqv_index == reg2%flst_real%eqv_index)
    match = match .and. (reg1%flst_uborn%eqv_index == reg2%flst_uborn%eqv_index)
    match = match .and. (reg1%mult == reg2%mult)
    match = match .and. (reg1%emitter == reg2%emitter)
    match = match .and. (reg1%nregions == reg2%nregions)
    match = match .and. (reg1%double_fsr .eqv. reg2%double_fsr)
    match = match .and. (reg1%soft_divergence .eqv. reg2%soft_divergence)
    match = match .and. (reg1%coll_divergence .eqv. reg2%coll_divergence)
    match = match .and. (char (reg1%nlo_correction_type) == char (reg2%nlo_correction_type))
    if (match) match = match .and. (all (reg1%ftuples == reg2%ftuples))
  end function singular_region_match

  module subroutine resonance_mapping_init (res_map, res_hist)
    class(resonance_mapping_t), intent(inout) :: res_map
    type(resonance_history_t), intent(in), dimension(:) :: res_hist
    integer :: n_hist, i_hist1, i_hist2, n_contributors
    n_contributors = 0
    n_hist = size (res_hist)
    allocate (res_map%res_histories (n_hist))
    do i_hist1 = 1, n_hist
       if (i_hist1 + 1 <= n_hist) then
          do i_hist2 = i_hist1 + 1, n_hist
             if (.not. (res_hist(i_hist1) .contains. res_hist(i_hist2))) &
                n_contributors = n_contributors + &
                                    res_hist(i_hist2)%n_resonances
          end do
       else
          n_contributors = n_contributors + res_hist(i_hist1)%n_resonances
       end if
    end do
    allocate (res_map%p_res (n_contributors))
    res_map%res_histories = res_hist
    res_map%p_res = vector4_null
  end subroutine resonance_mapping_init

  module subroutine resonance_mapping_set_alr_to_i_res &
       (res_map, regions, alr_new_to_old)
    class(resonance_mapping_t), intent(inout) :: res_map
    type(singular_region_t), intent(in), dimension(:) :: regions
    integer, intent(out), dimension(:), allocatable :: alr_new_to_old
    integer :: alr, i_res
    integer :: alr_new, n_alr_res
    integer :: k
    if (debug_on)  call msg_debug &
         (D_SUBTRACTION, "resonance_mapping_set_alr_to_i_res")
    n_alr_res = 0
    do alr = 1, size (regions)
       do i_res = 1, size (res_map%res_histories)
          if (res_map%res_histories(i_res)%contains_leg &
               (regions(alr)%emitter)) &
               n_alr_res = n_alr_res + 1
       end do
    end do

    allocate (res_map%alr_to_i_res (n_alr_res))
    allocate (res_map%i_res_to_alr (size (res_map%res_histories), 10))
    res_map%i_res_to_alr = 0
    allocate (alr_new_to_old (n_alr_res))
    alr_new = 1
    do alr = 1, size (regions)
       do i_res = 1, size (res_map%res_histories)
          if (res_map%res_histories(i_res)%contains_leg &
               (regions(alr)%emitter)) then
             res_map%alr_to_i_res (alr_new) = i_res
             alr_new_to_old (alr_new) = alr
             alr_new = alr_new  + 1
          end if
       end do
    end do

    do i_res = 1, size (res_map%res_histories)
       k = 1
       do alr = 1, size (regions)
          if (res_map%res_histories(i_res)%contains_leg &
               (regions(alr)%emitter)) then
             res_map%i_res_to_alr (i_res, k) = alr
             k = k + 1
          end if
       end do
    end do
    if (debug_active (D_SUBTRACTION)) then
       print *, 'i_res_to_alr:'
       do i_res = 1, size(res_map%i_res_to_alr, dim=1)
          print *, res_map%i_res_to_alr (i_res, :)
       end do
       print *, 'alr_new_to_old:', alr_new_to_old
    end if
  end subroutine resonance_mapping_set_alr_to_i_res

  module function resonance_mapping_get_resonance_history &
       (res_map, alr) result (res_hist)
    type(resonance_history_t) :: res_hist
    class(resonance_mapping_t), intent(in) :: res_map
    integer, intent(in) :: alr
    res_hist = res_map%res_histories(res_map%alr_to_i_res (alr))
  end function resonance_mapping_get_resonance_history

  module subroutine resonance_mapping_write (res_map)
    class(resonance_mapping_t), intent(in) :: res_map
    integer :: i_res
    do i_res = 1, size (res_map%res_histories)
       call res_map%res_histories(i_res)%write ()
    end do
  end subroutine resonance_mapping_write

  module function resonance_mapping_get_resonance_value &
       (res_map, i_res, p, i_gluon) result (p_map)
    real(default) :: p_map
    class(resonance_mapping_t), intent(in) :: res_map
    integer, intent(in) :: i_res
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in), optional :: i_gluon
    p_map = res_map%res_histories(i_res)%mapping (p, i_gluon)
  end function resonance_mapping_get_resonance_value

  module function resonance_mapping_get_resonance_all &
       (res_map, alr, p, i_gluon) result (p_map)
    real(default) :: p_map
    class(resonance_mapping_t), intent(in) :: res_map
    integer, intent(in) :: alr
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in), optional :: i_gluon
    integer :: i_res
    p_map = zero
    do i_res = 1, size (res_map%res_histories)
       associate (res => res_map%res_histories(i_res))
          if (any (res_map%i_res_to_alr (i_res, :) == alr)) &
             p_map = p_map + res%mapping (p, i_gluon)
       end associate
    end do
  end function resonance_mapping_get_resonance_all

  module function resonance_mapping_get_weight (res_map, alr, p) result (pfr)
    real(default) :: pfr
    class(resonance_mapping_t), intent(in) :: res_map
    integer, intent(in) :: alr
    type(vector4_t), intent(in), dimension(:) :: p
    real(default) :: sumpfr
    integer :: i_res
    sumpfr = zero
    do i_res = 1, size (res_map%res_histories)
       sumpfr = sumpfr + res_map%get_resonance_value (i_res, p)
    end do
    pfr = res_map%get_resonance_value (res_map%alr_to_i_res (alr), p) / sumpfr
  end function resonance_mapping_get_weight

  module function resonance_mapping_get_resonance_alr &
       (res_map, alr, p, i_gluon) result (p_map)
    real(default) :: p_map
    class(resonance_mapping_t), intent(in) :: res_map
    integer, intent(in) :: alr
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in), optional :: i_gluon
    integer :: i_res
    i_res = res_map%alr_to_i_res (alr)
    p_map = res_map%res_histories(i_res)%mapping (p, i_gluon)
  end function resonance_mapping_get_resonance_alr

  module subroutine resonance_mapping_assign (res_map_out, res_map_in)
    type(resonance_mapping_t), intent(out) :: res_map_out
    type(resonance_mapping_t), intent(in) :: res_map_in
    if (allocated (res_map_in%res_histories)) then
       allocate (res_map_out%res_histories (size (res_map_in%res_histories)))
       res_map_out%res_histories = res_map_in%res_histories
    end if
    if (allocated (res_map_in%alr_to_i_res)) then
       allocate (res_map_out%alr_to_i_res (size (res_map_in%alr_to_i_res)))
       res_map_out%alr_to_i_res = res_map_in%alr_to_i_res
    end if
    if (allocated (res_map_in%i_res_to_alr)) then
       allocate (res_map_out%i_res_to_alr (size (res_map_in%i_res_to_alr, 1), &
            size (res_map_in%i_res_to_alr, 2)))
       res_map_out%i_res_to_alr = res_map_in%i_res_to_alr
    end if
    if (allocated (res_map_in%p_res)) then
       allocate (res_map_out%p_res (size (res_map_in%p_res)))
       res_map_out%p_res = res_map_in%p_res
    end if
  end subroutine resonance_mapping_assign

  module subroutine region_data_init (reg_data, n_in, model, flavor_born, &
       flavor_real, nlo_correction_type, alpha_pow, alphas_pow)
    class(region_data_t), intent(out) :: reg_data
    integer, intent(in) :: n_in, alpha_pow, alphas_pow
    type(model_t), intent(in) :: model
    integer, intent(in), dimension(:,:) :: flavor_born, flavor_real
    type(string_t), intent(in) :: nlo_correction_type
    integer, dimension(:,:), allocatable ::  flv_real_tmp
    type(ftuple_list_t), dimension(:), allocatable :: ftuples
    integer, dimension(:), allocatable :: emitter
    type(flv_structure_t), dimension(:), allocatable :: flst_alr
    integer :: i, n_real
    integer :: n_flv_real_before_check
    reg_data%n_in = n_in
    reg_data%alpha_power = alpha_pow
    reg_data%alphas_power = alphas_pow
    reg_data%n_flv_born = size (flavor_born, dim = 2)
    reg_data%n_legs_born = size (flavor_born, dim = 1)
    reg_data%n_legs_real = reg_data%n_legs_born + 1
    reg_data%nlo_correction_type = nlo_correction_type
    n_flv_real_before_check = size (flavor_real, dim = 2)
    allocate (reg_data%flv_born (reg_data%n_flv_born))
    allocate (flv_real_tmp (reg_data%n_legs_real, n_flv_real_before_check))
    do i = 1, reg_data%n_flv_born
       call reg_data%flv_born(i)%init (flavor_born (:, i), n_in)
    end do
    n_real = 0
    do i = 1, n_flv_real_before_check
       if (nlo_correction_type == "EW") then
          if (.not. (query_coupling_powers &
              (flavor_real (:, i), reg_data%alpha_power + 1, &
               reg_data%alphas_power))) cycle
       end if
       n_real = n_real + 1
       flv_real_tmp (:, n_real) = flavor_real (:, i)
    end do
    allocate (reg_data%flv_real (n_real))
    do i = 1, n_real
       call reg_data%flv_real(i)%init (flv_real_tmp (:, i), n_in)
    end do

    call reg_data%find_regions (model, ftuples, emitter, flst_alr)
    call reg_data%init_singular_regions (ftuples, emitter, &
         flst_alr, nlo_correction_type)
    reg_data%n_flv_real = maxval (reg_data%regions%real_index)
    call reg_data%find_emitters ()
    call reg_data%set_mass_color_and_charge (model)
    call reg_data%set_splitting_info ()
  end subroutine region_data_init

  module subroutine region_data_init_resonance_information (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    call reg_data%enlarge_singular_regions_with_resonances ()
    call reg_data%find_resonances ()
  end subroutine region_data_init_resonance_information

  module subroutine region_data_set_resonance_mappings &
       (reg_data, resonance_histories)
    class(region_data_t), intent(inout) :: reg_data
    type(resonance_history_t), intent(in), dimension(:) :: resonance_histories
    select type (map => reg_data%fks_mapping)
    type is (fks_mapping_resonances_t)
       call map%res_map%init (resonance_histories)
    end select
  end subroutine region_data_set_resonance_mappings

  module subroutine region_data_setup_fks_mappings (reg_data, template, n_in)
    class(region_data_t), intent(inout) :: reg_data
    type(fks_template_t), intent(in) :: template
    integer, intent(in) :: n_in
    call reg_data%allocate_fks_mappings (template%mapping_type)
    select type (map => reg_data%fks_mapping)
    type is (fks_mapping_default_t)
       call map%set_parameter (n_in, template%fks_dij_exp1, &
            template%fks_dij_exp2)
    end select
  end subroutine region_data_setup_fks_mappings

  module subroutine region_data_enlarge_singular_regions_with_resonances &
       (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    integer :: alr
    integer, dimension(:), allocatable :: alr_new_to_old
    integer :: n_alr_new
    type(singular_region_t), dimension(:), allocatable :: save_regions
    if (debug_on)  call msg_debug &
         (D_SUBTRACTION, "region_data_enlarge_singular_regions_with_resonances")
    call debug_input_values ()
    select type (fks_mapping => reg_data%fks_mapping)
    type is (fks_mapping_default_t)
       return
    type is (fks_mapping_resonances_t)
       allocate (save_regions (reg_data%n_regions))
       do alr = 1, reg_data%n_regions
          save_regions(alr) = reg_data%regions(alr)
       end do

       associate (res_map => fks_mapping%res_map)
          call res_map%set_alr_to_i_res (reg_data%regions, alr_new_to_old)
          deallocate (reg_data%regions)
          n_alr_new = size (alr_new_to_old)
          reg_data%n_regions = n_alr_new
          allocate (reg_data%regions (n_alr_new))
          do alr = 1, n_alr_new
             reg_data%regions(alr) = save_regions(alr_new_to_old (alr))
             reg_data%regions(alr)%i_res = res_map%alr_to_i_res (alr)
          end do
       end associate
    end select

  contains

    subroutine debug_input_values ()
      if (debug2_active (D_SUBTRACTION)) then
         call reg_data%write ()
      end if
    end subroutine debug_input_values

  end subroutine region_data_enlarge_singular_regions_with_resonances

  module subroutine region_data_set_isr_pseudo_regions (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    integer :: alr
    integer :: n_alr_new
    !!! Subroutine called for threshold factorization ->
    !!! Size of singular regions at this point is fixed
    type(singular_region_t), dimension(2) :: save_regions
    integer, dimension(4) :: alr_new_to_old
    do alr = 1, reg_data%n_regions
       save_regions(alr) = reg_data%regions(alr)
    end do
    n_alr_new = reg_data%n_regions * 2
    alr_new_to_old = [1, 1, 2, 2]
    deallocate (reg_data%regions)
    allocate (reg_data%regions (n_alr_new))
    reg_data%n_regions = n_alr_new
    do alr = 1, n_alr_new
       reg_data%regions(alr) = save_regions(alr_new_to_old (alr))
       call add_pseudo_emitters (reg_data%regions(alr))
       if (mod (alr, 2) == 0) reg_data%regions(alr)%pseudo_isr = .true.
    end do
  contains
    subroutine add_pseudo_emitters (sregion)
      type(singular_region_t), intent(inout) :: sregion
      type(ftuple_t), dimension(2) :: ftuples_save
      integer :: alr
      do alr = 1, 2
         ftuples_save(alr) = sregion%ftuples(alr)
      end do
      deallocate (sregion%ftuples)
      sregion%nregions = sregion%nregions * 2
      allocate (sregion%ftuples (sregion%nregions))
      do alr = 1, sregion%nregions
         sregion%ftuples(alr) = ftuples_save (alr_new_to_old(alr))
         if (mod (alr, 2) == 0) sregion%ftuples(alr)%pseudo_isr = .true.
      end do
    end subroutine add_pseudo_emitters
  end subroutine region_data_set_isr_pseudo_regions

  module subroutine region_data_split_up_interference_regions_for_threshold &
       (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    integer :: alr, i_ftuple
    integer :: current_emitter
    integer :: i1, i2
    integer :: n_new_reg
    type(ftuple_t), dimension(2) :: ftuples
    do alr = 1, reg_data%n_regions
       associate (region => reg_data%regions(alr))
          current_emitter = region%emitter
          n_new_reg = 0
          do i_ftuple = 1, region%nregions
             call region%ftuples(i_ftuple)%get (i1, i2)
             if (i1 == current_emitter) then
                n_new_reg = n_new_reg + 1
                ftuples(n_new_reg) = region%ftuples(i_ftuple)
             end if
          end do
          deallocate (region%ftuples)
          allocate (region%ftuples(n_new_reg))
          region%ftuples = ftuples (1 : n_new_reg)
          region%nregions = n_new_reg
       end associate
    end do
    reg_data%fks_mapping%normalization_factor = 0.5_default
  end subroutine region_data_split_up_interference_regions_for_threshold

  module subroutine region_data_set_mass_color_and_charge (reg_data, model)
    class(region_data_t), intent(inout) :: reg_data
    type(model_t), intent(in) :: model
    integer :: i
    do i = 1, reg_data%n_regions
       associate (region => reg_data%regions(i))
          call region%flst_uborn%init_mass_color_and_charge (model)
          call region%flst_real%init_mass_color_and_charge (model)
       end associate
    end do
    do i = 1, reg_data%n_flv_born
       call reg_data%flv_born(i)%init_mass_color_and_charge (model)
    end do
    do i = 1, size (reg_data%flv_real)
       call reg_data%flv_real(i)%init_mass_color_and_charge (model)
    end do
  end subroutine region_data_set_mass_color_and_charge

  module function region_data_uses_resonances (reg_data) result (val)
    logical :: val
    class(region_data_t), intent(in) :: reg_data
    select type (fks_mapping => reg_data%fks_mapping)
    type is (fks_mapping_resonances_t)
       val = .true.
    class default
       val = .false.
    end select
  end function region_data_uses_resonances

  pure module function region_data_get_emitter_list (reg_data) result (emitters)
    class(region_data_t), intent(in) :: reg_data
    integer, dimension(:), allocatable :: emitters
    integer :: i
    allocate (emitters (reg_data%n_regions))
    do i = 1, reg_data%n_regions
       emitters(i) = reg_data%regions(i)%emitter
    end do
  end function region_data_get_emitter_list

  module function region_data_get_n_emitters_sc &
       (reg_data) result (n_emitters_sc)
    class(region_data_t), intent(in) :: reg_data
    integer :: n_emitters_sc
    n_emitters_sc = count (reg_data%emitters /= 0)
  end function region_data_get_n_emitters_sc

  module function region_data_get_associated_resonances &
       (reg_data, emitter) result (res)
    integer, dimension(:), allocatable :: res
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: emitter
    integer :: alr, i
    integer :: n_res
    select type (fks_mapping => reg_data%fks_mapping)
    type is (fks_mapping_resonances_t)
       n_res = 0

       do alr = 1, reg_data%n_regions
          if (reg_data%regions(alr)%emitter == emitter) &
             n_res = n_res + 1
       end do

       if (n_res > 0) then
          allocate (res (n_res))
       else
          return
       end if
       i = 1

       do alr = 1, reg_data%n_regions
          if (reg_data%regions(alr)%emitter == emitter) then
             res (i) = fks_mapping%res_map%alr_to_i_res (alr)
             i = i + 1
          end if
       end do
    end select
  end function region_data_get_associated_resonances

  module function region_data_emitter_is_compatible_with_resonance &
       (reg_data, i_res, emitter) result (compatible)
    logical :: compatible
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: i_res, emitter
    integer :: i_res_alr, alr
    compatible = .false.
    select type (fks_mapping => reg_data%fks_mapping)
    type is (fks_mapping_resonances_t)
       do alr = 1, reg_data%n_regions
          i_res_alr = fks_mapping%res_map%alr_to_i_res (alr)
          if (i_res_alr == i_res .and. reg_data%get_emitter(alr) == emitter) then
             compatible = .true.
             exit
          end if
       end do
    end select
  end function region_data_emitter_is_compatible_with_resonance

  module function region_data_emitter_is_in_resonance &
       (reg_data, i_res, emitter) result (exist)
    logical :: exist
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: i_res, emitter
    integer :: i
    exist = .false.
    select type (fks_mapping => reg_data%fks_mapping)
    type is (fks_mapping_resonances_t)
       associate (res_history => fks_mapping%res_map%res_histories(i_res))
          do i = 1, res_history%n_resonances
             exist = exist .or. &
                  any (res_history%resonances(i)%contributors%c == emitter)
          end do
      end associate
    end select
  end function region_data_emitter_is_in_resonance

  module subroutine region_data_get_contributors &
       (reg_data, i_res, emitter, c, success)
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: i_res, emitter
    integer, intent(inout), dimension(:), allocatable :: c
    logical, intent(out) :: success
    integer :: i
    success = .false.
    select type (fks_mapping => reg_data%fks_mapping)
    type is (fks_mapping_resonances_t)
       associate (res_history => fks_mapping%res_map%res_histories (i_res))
          do i = 1, res_history%n_resonances
             if (any (res_history%resonances(i)%contributors%c == emitter)) then
                allocate (c (size (res_history%resonances(i)%contributors%c)))
                c = res_history%resonances(i)%contributors%c
                success = .true.
                exit
             end if
          end do
       end associate
    end select
  end subroutine region_data_get_contributors

  pure module function region_data_get_emitter (reg_data, alr) result (emitter)
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: alr
    integer :: emitter
    emitter = reg_data%regions(alr)%emitter
  end function region_data_get_emitter

  module function region_data_map_real_to_born_index &
       (reg_data, real_index) result (uborn_index)
    integer :: uborn_index
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: real_index
    integer :: alr
    uborn_index = 0
    do alr = 1, size (reg_data%regions)
       if (reg_data%regions(alr)%real_index == real_index) then
          uborn_index = reg_data%regions(alr)%uborn_index
          exit
       end if
    end do
  end function region_data_map_real_to_born_index

  module function region_data_get_flv_states_born_array &
       (reg_data) result (flv_states)
    integer, dimension(:,:), allocatable :: flv_states
    class(region_data_t), intent(in) :: reg_data
    integer :: i_flv
    allocate (flv_states (reg_data%n_legs_born, reg_data%n_flv_born))
    do i_flv = 1, reg_data%n_flv_born
       flv_states (:, i_flv) = reg_data%flv_born(i_flv)%flst
    end do
  end function region_data_get_flv_states_born_array

  module function region_data_get_flv_states_born_single &
       (reg_data, i_flv) result (flv_states)
    integer, dimension(:), allocatable :: flv_states
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: i_flv
    allocate (flv_states (reg_data%n_legs_born))
    flv_states = reg_data%flv_born(i_flv)%flst
  end function region_data_get_flv_states_born_single

  module function region_data_get_flv_states_real_single &
       (reg_data, i_flv) result (flv_states)
    integer, dimension(:), allocatable :: flv_states
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: i_flv
    integer :: i_reg
    allocate (flv_states (reg_data%n_legs_real))
    do i_reg = 1, reg_data%n_regions
       if (i_flv == reg_data%regions(i_reg)%real_index) then
          flv_states = reg_data%regions(i_reg)%flst_real%flst
          exit
       end if
    end do
  end function region_data_get_flv_states_real_single

  module function region_data_get_flv_states_real_array &
       (reg_data) result (flv_states)
    integer, dimension(:,:), allocatable :: flv_states
    class(region_data_t), intent(in) :: reg_data
    integer :: i_flv
    allocate (flv_states (reg_data%n_legs_real, reg_data%n_flv_real))
    do i_flv = 1, reg_data%n_flv_real
       flv_states (:, i_flv) = reg_data%get_flv_states_real (i_flv)
    end do
  end function region_data_get_flv_states_real_array

  module subroutine region_data_get_all_flv_states &
       (reg_data, flv_born, flv_real)
    class(region_data_t), intent(in) :: reg_data
    integer, dimension(:,:), allocatable, intent(out) :: flv_born, flv_real
    allocate (flv_born (reg_data%n_legs_born, reg_data%n_flv_born))
    flv_born = reg_data%get_flv_states_born ()
    allocate (flv_real (reg_data%n_legs_real, reg_data%n_flv_real))
    flv_real = reg_data%get_flv_states_real ()
  end subroutine region_data_get_all_flv_states

  module function region_data_get_n_in (reg_data) result (n_in)
    integer :: n_in
    class(region_data_t), intent(in) :: reg_data
    n_in = reg_data%n_in
  end function region_data_get_n_in

  module function region_data_get_n_legs_real (reg_data) result (n_legs)
    integer :: n_legs
    class(region_data_t), intent(in) :: reg_data
    n_legs = reg_data%n_legs_real
  end function region_data_get_n_legs_real

  module function region_data_get_n_legs_born (reg_data) result (n_legs)
    integer :: n_legs
    class(region_data_t), intent(in) :: reg_data
    n_legs = reg_data%n_legs_born
  end function region_data_get_n_legs_born

  module function region_data_get_n_flv_real (reg_data) result (n_flv)
    integer :: n_flv
    class(region_data_t), intent(in) :: reg_data
    n_flv = reg_data%n_flv_real
  end function region_data_get_n_flv_real

  module function region_data_get_n_flv_born (reg_data) result (n_flv)
    integer :: n_flv
    class(region_data_t), intent(in) :: reg_data
    n_flv = reg_data%n_flv_born
  end function region_data_get_n_flv_born

  module function region_data_get_svalue_ij &
       (reg_data, p_real, alr, i, j, i_res) result (sval)
    class(region_data_t), intent(inout) :: reg_data
    type(vector4_t), intent(in), dimension(:) :: p_real
    integer, intent(in) :: alr, i, j
    integer, intent(in) :: i_res
    real(default) :: sval
    associate (map => reg_data%fks_mapping)
       call map%compute_sumdij (reg_data%regions(alr), p_real)
       select type (map)
       type is (fks_mapping_resonances_t)
          map%i_con = reg_data%alr_to_i_contributor (alr)
       end select
       map%pseudo_isr = reg_data%regions(alr)%pseudo_isr
       sval = map%svalue (p_real, i, j, i_res) * map%normalization_factor
    end associate
  end function region_data_get_svalue_ij

  module function region_data_get_svalue_last_pos &
       (reg_data, p, alr, emitter, i_res) result (sval)
    class(region_data_t), intent(inout) :: reg_data
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in) :: alr, emitter
    integer, intent(in) :: i_res
    real(default) :: sval
    sval = reg_data%get_svalue (p, alr, emitter, reg_data%n_legs_real, i_res)
  end function region_data_get_svalue_last_pos

  module function region_data_get_svalue_soft &
       (reg_data, p_born, p_soft, alr, emitter, i_res) result (sval)
    class(region_data_t), intent(inout) :: reg_data
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(in) :: p_soft
    integer, intent(in) :: alr, emitter, i_res
    real(default) :: sval
    associate (map => reg_data%fks_mapping)
       call map%compute_sumdij_soft (reg_data%regions(alr), p_born, p_soft)
       select type (map)
       type is (fks_mapping_resonances_t)
          map%i_con = reg_data%alr_to_i_contributor (alr)
       end select
       map%pseudo_isr = reg_data%regions(alr)%pseudo_isr
       sval = map%svalue_soft (p_born, p_soft, emitter, i_res) * &
            map%normalization_factor
    end associate
  end function region_data_get_svalue_soft

  module subroutine region_data_find_regions &
       (reg_data, model, ftuples, emitters, flst_alr)
    class(region_data_t), intent(in) :: reg_data
    type(model_t), intent(in) :: model
    type(ftuple_list_t), intent(out), dimension(:), allocatable :: ftuples
    integer, intent(out), dimension(:), allocatable :: emitters
    type(flv_structure_t), intent(out), dimension(:), allocatable :: flst_alr
    type(ftuple_list_t), dimension(:,:), allocatable :: ftuples_tmp
    integer, dimension(:,:), allocatable :: ftuple_index
    integer :: n_born, n_real
    integer :: n_legreal
    integer :: i_born, i_real, i_ftuple
    integer :: last_registered_i_born, last_registered_i_real

    n_born = size (reg_data%flv_born)
    n_real = size (reg_data%flv_real)
    n_legreal = size (reg_data%flv_real(1)%flst)
    allocate (emitters (0))
    allocate (flst_alr (0))
    allocate (ftuples (0))
    i_ftuple = 0
    last_registered_i_born = 0; last_registered_i_real = 0

    do i_real = 1, n_real
       do i_born = 1, n_born
          call setup_flsts_emitters_and_ftuples_fsr &
               (i_real, i_born, i_ftuple, flst_alr, emitters, ftuples)
          call setup_flsts_emitters_and_ftuples_isr &
               (i_real, i_born, i_ftuple, flst_alr, emitters, ftuples)
       end do
    end do

  contains
    function incr_i_ftuple_if_required (i_born, i_real, i_ftuple_in) result (i_ftuple)
      integer :: i_ftuple
      integer, intent(in) :: i_born, i_real, i_ftuple_in
      if (last_registered_i_born /= i_born .or. last_registered_i_real /= i_real) then
         last_registered_i_born = i_born
         last_registered_i_real = i_real
         i_ftuple = i_ftuple_in + 1
      else
         i_ftuple = i_ftuple_in
      end if
    end function incr_i_ftuple_if_required

    subroutine setup_flsts_emitters_and_ftuples_fsr &
           (i_real, i_born, i_ftuple, flst_alr, emitters, ftuples)
      integer, intent(in) :: i_real, i_born
      integer, intent(inout) :: i_ftuple
      type(flv_structure_t), intent(inout), dimension(:), allocatable :: flst_alr
      integer, intent(inout), dimension(:), allocatable :: emitters
      type(ftuple_list_t), intent(inout), dimension(:), allocatable :: ftuples
      type(ftuple_list_t) :: ftuples_tmp
      type(flv_structure_t) :: flst_alr_tmp
      type(ftuple_t) :: current_ftuple
      integer :: leg1, leg2
      logical :: valid
      associate (flv_born => reg_data%flv_born(i_born), &
              flv_real => reg_data%flv_real(i_real))
         do leg1 = reg_data%n_in + 1, n_legreal
            do leg2 = leg1 + 1, n_legreal
               valid = flv_real%valid_pair(leg1, leg2, flv_born, model)
               if (valid) then
                  if (is_vector(flv_real%flst(leg1)) .and. &
                     is_fermion(flv_real%flst(leg2))) then
                     flst_alr_tmp = create_alr (flv_real, &
                          reg_data%n_in, leg2, leg1)
                  else
                     flst_alr_tmp = create_alr (flv_real, &
                          reg_data%n_in, leg1, leg2)
                  end if
                  flst_alr = [flst_alr, flst_alr_tmp]
                  emitters = [emitters, n_legreal - 1]
                  call current_ftuple%set (leg1, leg2)
                  call current_ftuple%determine_splitting_type_fsr &
                       (flv_real, leg1, leg2)
                  call current_ftuple%determine_sub_correction_type &
                       (flv_born, flv_real, leg1, leg2)
                  i_ftuple = incr_i_ftuple_if_required (i_born, i_real, i_ftuple)
                  if (i_ftuple > size (ftuples)) then
                     call ftuples_tmp%append (current_ftuple)
                     ftuples = [ftuples, ftuples_tmp]
                  else
                     call ftuples(i_ftuple)%append (current_ftuple)
                  end if
               end if
            end do
         end do
      end associate
    end subroutine setup_flsts_emitters_and_ftuples_fsr

    subroutine setup_flsts_emitters_and_ftuples_isr &
           (i_real, i_born, i_ftuple, flst_alr, emitters, ftuples)
      integer, intent(in) :: i_real, i_born
      integer, intent(inout) :: i_ftuple
      type(flv_structure_t), intent(inout), dimension(:), allocatable :: flst_alr
      integer, intent(inout), dimension(:), allocatable :: emitters
      type(ftuple_list_t), intent(inout), dimension(:), allocatable :: ftuples
      type(ftuple_list_t) :: ftuples_tmp
      type(flv_structure_t) :: flst_alr_tmp
      type(ftuple_t) :: current_ftuple
      integer :: leg, emitter
      logical :: valid1, valid2
      associate (flv_born => reg_data%flv_born(i_born), &
              flv_real => reg_data%flv_real(i_real))
         do leg = reg_data%n_in + 1, n_legreal
            valid1 = flv_real%valid_pair(1, leg, flv_born, model)
            if (reg_data%n_in > 1) then
               valid2 = flv_real%valid_pair(2, leg, flv_born, model)
            else
               valid2 = .false.
            end if
            if (valid1 .and. valid2) then
               emitter = 0
            else if (valid1 .and. .not. valid2) then
               emitter = 1
            else if (.not. valid1 .and. valid2) then
               emitter = 2
            else
               emitter = -1
            end if
            if (valid1 .or. valid2) then
               flst_alr_tmp = create_alr (flv_real, reg_data%n_in, emitter, leg)
               flst_alr = [flst_alr, flst_alr_tmp]
               emitters = [emitters, emitter]
               call current_ftuple%set(emitter, leg)
               call current_ftuple%determine_splitting_type_isr &
                    (flv_real, emitter, leg)
               call current_ftuple%determine_sub_correction_type &
                    (flv_born, flv_real, emitter, leg)
               i_ftuple = incr_i_ftuple_if_required (i_born, i_real, i_ftuple)
               if (i_ftuple > size (ftuples)) then
                  call ftuples_tmp%append (current_ftuple)
                  ftuples = [ftuples, ftuples_tmp]
               else
                  call ftuples(i_ftuple)%append (current_ftuple)
               end if
            end if
         end do
      end associate
    end subroutine setup_flsts_emitters_and_ftuples_isr

  end subroutine region_data_find_regions

  module subroutine region_data_find_eqv_regions (reg_data, optimize)
    class(region_data_t), intent(inout) :: reg_data
    logical, intent(in) :: optimize
    integer :: n_reg, alr1, alr2
    n_reg = reg_data%n_regions
    if (optimize) then
       do alr1 = 1, n_reg
          reg_data%regions(alr1)%flst_uborn%eqv_index = &
               reg_data%eqv_flv_index_born(reg_data%regions(alr1)%uborn_index)
          reg_data%regions(alr1)%flst_real%eqv_index = &
               reg_data%eqv_flv_index_real(reg_data%regions(alr1)%real_index)
       end do
       do alr1 = 1, n_reg
          do alr2 = 1, alr1
             if (reg_data%regions(alr2) .match. reg_data%regions(alr1)) then
                reg_data%regions(alr1)%eqv_index = alr2
                exit
             end if
          end do
       end do
    else
       do alr1 = 1, n_reg
          reg_data%regions(alr1)%eqv_index = alr1
       end do
    end if
  end subroutine region_data_find_eqv_regions

  module subroutine region_data_init_singular_regions &
       (reg_data, ftuples, emitter, flv_alr, nlo_correction_type)
    class(region_data_t), intent(inout) :: reg_data
    type(ftuple_list_t), intent(inout), dimension(:), allocatable :: ftuples
    type(string_t), intent(in) :: nlo_correction_type
    type(string_t), dimension(:), allocatable :: nlo_correction_type_dyn
    integer, intent(in), dimension(:) :: emitter
    integer :: n_independent_flv
    type(flv_structure_t), intent(in), dimension(:) :: flv_alr
    type(flv_structure_t), dimension(:), allocatable :: flv_uborn, &
         flv_alr_registered
    integer, dimension(:), allocatable :: mult
    integer, dimension(:), allocatable :: flst_emitter
    integer :: n_regions, maxregions
    integer, dimension(:), allocatable :: index
    integer :: i, i_flv, n_legs
    logical :: equiv, valid_fs_splitting, pure_corr_type, &
         corr_type_valid, non_singular_reg
    integer :: i_first, i_reg, i_reg_prev
    integer, dimension(:), allocatable :: region_to_ftuple, alr_limits
    integer, dimension(:), allocatable :: equiv_index

    maxregions = size (emitter)
    n_legs = flv_alr(1)%nlegs

    allocate (flv_uborn (maxregions))
    allocate (flv_alr_registered (maxregions))
    allocate (mult (maxregions))
    mult = 0
    allocate (flst_emitter (maxregions))
    allocate (index (0))
    allocate (region_to_ftuple (maxregions))
    allocate (equiv_index (maxregions))
    allocate (nlo_correction_type_dyn (maxregions))

    call setup_region_mappings (n_independent_flv, alr_limits, region_to_ftuple)
    nlo_correction_type_dyn = nlo_correction_type
    i_first = 1
    i_reg = 1
    SCAN_FLAVORS: do i_flv = 1, n_independent_flv
       SCAN_FTUPLES: do i = i_first, i_first + alr_limits (i_flv) - 1
          equiv = .false.
          corr_type_valid = .true.
          non_singular_reg = .false.
          if (i == i_first) then
             flv_alr_registered(i_reg) = flv_alr(i)
             if (nlo_correction_type == "EW" .and. &
                  reg_data%alphas_power > 0) then
                nlo_correction_type_dyn (i_reg) = &
                     set_dynamic_correction_type (i)
             end if
             flv_uborn(i_reg) = flv_alr(i)%create_uborn (emitter(i), &
                                     nlo_correction_type_dyn (i_reg))
             flst_emitter(i_reg) = emitter(i)
             equiv_index(i_reg) = region_to_ftuple(i)
             if (nlo_correction_type == "EW" .and. reg_data%alphas_power > 0) then
                corr_type_valid = (nlo_correction_type_dyn (i_reg) == "EW" .and. &
                                  query_coupling_powers (flv_uborn(i_reg)%flst, &
                                  reg_data%alpha_power, reg_data%alphas_power)) &
                                  .or. (nlo_correction_type_dyn (i_reg) == "QCD" .and. &
                                  query_coupling_powers (flv_uborn(i_reg)%flst, &
                                  reg_data%alpha_power + 1, reg_data%alphas_power - 1))
                non_singular_reg = .not. corr_type_valid .and. &
                                   qcd_ew_interferences (flv_alr_registered(i_reg)%flst) &
                                   .and. query_coupling_powers (flv_alr_registered(i_reg)%flst, &
                                   reg_data%alpha_power+2, reg_data%alphas_power-1) &
                                   .and. .not. qcd_ew_interferences (flv_uborn(i_reg)%flst) &
                                   .and. nlo_correction_type_dyn (i_reg) == "QCD"
                if (non_singular_reg) nlo_correction_type_dyn (i_reg) = "none"
             end if
             if (corr_type_valid .or. non_singular_reg) then
                mult(i_reg) = mult(i_reg) + 1
                index = [index, region_to_real_index(ftuples, i)]
                i_reg = i_reg + 1
             end if
          else
             !!! Check for equivalent flavor structures
             do i_reg_prev = 1, i_reg - 1
                if (emitter(i) == flst_emitter(i_reg_prev) .and. emitter(i) > reg_data%n_in) then
                   valid_fs_splitting = check_fs_splitting &
                         (flv_alr(i)%get_last_two(n_legs), &
                          flv_alr_registered(i_reg_prev)%get_last_two(n_legs), &
                          flv_alr(i)%tag(n_legs - 1), &
                          flv_alr_registered(i_reg_prev)%tag(n_legs - 1))
                   if (nlo_correction_type == "EW" .and. reg_data%alphas_power > 0) then
                      nlo_correction_type_dyn (i_reg) = set_dynamic_correction_type (i)
                   end if
                   pure_corr_type = nlo_correction_type_dyn (i_reg) &
                                     == nlo_correction_type_dyn (i_reg_prev)
                   if ((flv_alr(i) .equiv. flv_alr_registered(i_reg_prev)) &
                        .and. valid_fs_splitting .and. pure_corr_type) then
                      mult(i_reg_prev) = mult(i_reg_prev) + 1
                      equiv = .true.
                      call ftuples(region_to_real_index(ftuples, i))%set_equiv &
                           (equiv_index(i_reg_prev), region_to_ftuple(i))
                      exit
                   end if
                else if (emitter(i) == flst_emitter(i_reg_prev) .and. emitter(i) <= reg_data%n_in) then
                   if (nlo_correction_type == "EW" .and. reg_data%alphas_power > 0) then
                      nlo_correction_type_dyn (i_reg) = set_dynamic_correction_type (i)
                   end if
                   pure_corr_type = nlo_correction_type_dyn (i_reg) &
                                     == nlo_correction_type_dyn (i_reg_prev)
                   if ((flv_alr(i) .equiv. flv_alr_registered(i_reg_prev)) &
                        .and. pure_corr_type) then
                      mult(i_reg_prev) = mult(i_reg_prev) + 1
                      equiv = .true.
                      call ftuples(region_to_real_index(ftuples, i))%set_equiv &
                           (equiv_index(i_reg_prev), region_to_ftuple(i))
                      exit
                   end if
                end if
             end do
             if (.not. equiv) then
                flv_alr_registered(i_reg) = flv_alr(i)
                if (nlo_correction_type == "EW" .and. reg_data%alphas_power > 0) then
                   nlo_correction_type_dyn (i_reg) = set_dynamic_correction_type (i)
                end if
                flv_uborn(i_reg) = flv_alr(i)%create_uborn (emitter(i), &
                                       nlo_correction_type_dyn (i_reg))
                flst_emitter(i_reg) = emitter(i)
                equiv_index (i_reg) = region_to_ftuple(i)
                if (nlo_correction_type == "EW" .and. reg_data%alphas_power > 0) then
                   corr_type_valid = (nlo_correction_type_dyn (i_reg) == "EW" .and. &
                                query_coupling_powers (flv_uborn(i_reg)%flst, &
                                reg_data%alpha_power, reg_data%alphas_power)) &
                                .or. (nlo_correction_type_dyn (i_reg) == "QCD" .and. &
                                query_coupling_powers (flv_uborn(i_reg)%flst, &
                                reg_data%alpha_power + 1, reg_data%alphas_power - 1))
                   non_singular_reg = .not. corr_type_valid .and. &
                                 qcd_ew_interferences (flv_alr_registered(i_reg)%flst) &
                                 .and. query_coupling_powers (flv_alr_registered(i_reg)%flst, &
                                 reg_data%alpha_power+2, reg_data%alphas_power-1) &
                                 .and. .not. qcd_ew_interferences (flv_uborn(i_reg)%flst) &
                                 .and. nlo_correction_type_dyn (i_reg) == "QCD"
                   if (non_singular_reg) nlo_correction_type_dyn (i_reg) = "none"
                end if
                if (corr_type_valid .or. non_singular_reg) then
                   mult(i_reg) = mult(i_reg) + 1
                   index = [index, region_to_real_index(ftuples, i)]
                   i_reg = i_reg + 1
                end if
             end if
          end if
       end do SCAN_FTUPLES
       i_first = i_first + alr_limits(i_flv)
    end do SCAN_FLAVORS
    n_regions = i_reg - 1

    allocate (reg_data%regions (n_regions))
    reg_data%n_regions = n_regions
    call account_for_regions_from_other_uborns (ftuples)
    call init_regions_with_permuted_flavors ()
    call assign_real_indices ()

    deallocate (flv_uborn)
    deallocate (flv_alr_registered)
    deallocate (mult)
    deallocate (flst_emitter)
    deallocate (index)
    deallocate (region_to_ftuple)
    deallocate (equiv_index)

  contains

    subroutine account_for_regions_from_other_uborns (ftuples)
      type(ftuple_list_t), intent(inout), dimension(:), allocatable :: ftuples
      integer :: alr1, alr2, i
      type(ftuple_t), dimension(:), allocatable :: ftuples_alr1, ftuples_alr2
      logical, dimension(:,:), allocatable :: equivalences
      do alr1 = 1, n_regions
         do alr2 = 1, n_regions
            if (index(alr1) == index(alr2)) cycle
            if (flv_alr_registered(alr1) .equiv. flv_alr_registered(alr2)) then
               call ftuples(index(alr1))%to_array (ftuples_alr1, equivalences, .false.)
               call ftuples(index(alr2))%to_array (ftuples_alr2, equivalences, .false.)
               do i = 1, size (ftuples_alr2)
                  if (.not. any (ftuple_equal_ireg (ftuples_alr1, ftuples_alr2(i)))) then
                     call ftuples(index(alr1))%append (ftuples_alr2(i))
                  end if
               end do
            end if
         end do
      end do
    end subroutine account_for_regions_from_other_uborns

    subroutine setup_region_mappings (n_independent_flv, &
           alr_limits, region_to_ftuple)
       integer, intent(inout) :: n_independent_flv
       integer, intent(inout), dimension(:), allocatable :: alr_limits
       integer, intent(inout), dimension(:), allocatable :: region_to_ftuple
       integer :: i, j, i_flv
       if (any (ftuples%get_n_tuples() == 0)) &
            call msg_fatal ("Inconsistent collection of FKS pairs!")
       n_independent_flv = size (ftuples)
       alr_limits = ftuples%get_n_tuples()
       if (.not. (sum (alr_limits) == maxregions)) &
            call msg_fatal ("Too many regions!")
       j = 1
       do i_flv = 1, n_independent_flv
          do i = 1, alr_limits(i_flv)
             region_to_ftuple(j) = i
             j = j + 1
          end do
       end do
    end subroutine setup_region_mappings

    subroutine check_permutation (perm, flv_perm, flv_orig, i_reg)
      type(flv_perm_t), intent(in) :: perm
      type(flv_structure_t), intent(in) :: flv_perm, flv_orig
      integer, intent(in) :: i_reg
      type(flv_structure_t) :: flv_tmp
      flv_tmp = perm%apply (flv_orig, invert = .true.)
      if (.not. all (flv_tmp%flst == flv_perm%flst)) then
         print *, 'Fail at: ', i_reg
         print *, 'Original flavor structure: ', flv_orig%flst
         call perm%write ()
         print *, 'Permuted flavor: ', flv_perm%flst
         print *, 'Should be: ', flv_tmp%flst
         call msg_fatal ("Permutation does not reproduce original flavor!")
      end if
    end subroutine check_permutation

    subroutine init_regions_with_permuted_flavors ()
       type(flv_perm_t) :: perm_list
       type(ftuple_t), dimension(:), allocatable :: ftuple_array
       logical, dimension(:,:), allocatable :: equivalences
       integer :: i, j
       do j = 1, n_regions
          do i = 1, reg_data%n_flv_born
             if (reg_data%flv_born (i) .equiv. flv_uborn (j)) then
                call perm_list%reset ()
                call perm_list%init (reg_data%flv_born(i), flv_uborn(j), &
                     reg_data%n_in, reg_data%n_legs_born, .true.)
                flv_uborn(j) = perm_list%apply (flv_uborn(j))
                flv_alr_registered(j) = perm_list%apply (flv_alr_registered(j))
                flst_emitter(j) = perm_list%apply (flst_emitter(j))
             end if
          end do
          call ftuples(index(j))%to_array (ftuple_array, equivalences, .false.)
          do i = 1, size (reg_data%flv_real)
             if (reg_data%flv_real(i) .equiv. flv_alr_registered(j)) then
                call perm_list%reset ()
                call perm_list%init (flv_alr_registered(j), reg_data%flv_real(i), &
                     reg_data%n_in, reg_data%n_legs_real, .false.)
                if (debug_active (D_SUBTRACTION)) call check_permutation &
                     (perm_list, reg_data%flv_real(i), flv_alr_registered(j), j)
                ftuple_array = perm_list%apply (ftuple_array)
                call ftuple_sort_array (ftuple_array, equivalences)
             end if
          end do
          call reg_data%regions(j)%init (j, mult(j), 0, flv_alr_registered(j), &
               flv_uborn(j), reg_data%flv_born, flst_emitter(j), ftuple_array, &
               equivalences, nlo_correction_type_dyn (j))
          if (allocated (ftuple_array)) deallocate (ftuple_array)
          if (allocated (equivalences)) deallocate (equivalences)
       end do
    end subroutine init_regions_with_permuted_flavors

    subroutine assign_real_indices ()
      type(flv_structure_t) :: current_flv_real
      type(flv_structure_t), dimension(:), allocatable :: these_flv
      integer :: i_real, current_uborn_index
      integer :: i, j, this_i_real
      allocate (these_flv (size (flv_alr_registered)))
      i_real = 1
      associate (regions => reg_data%regions)
         do i = 1, reg_data%n_regions
            do j = 1, size (these_flv)
               if (.not. allocated (these_flv(j)%flst)) then
                  this_i_real = i_real
                  call these_flv(i_real)%init (flv_alr_registered(i)%flst, reg_data%n_in)
                  i_real = i_real + 1
                  exit
               else if (all (these_flv(j)%flst == flv_alr_registered(i)%flst)) then
                  this_i_real = j
                  exit
               end if
            end do
            regions(i)%real_index = this_i_real
         end do
      end associate
      deallocate (these_flv)
    end subroutine assign_real_indices

    function check_fs_splitting (flv1, flv2, tag1, tag2) result (valid)
      logical :: valid
      integer, intent(in), dimension(2) :: flv1, flv2
      integer, intent(in) :: tag1, tag2
      if (flv1(1) + flv1(2) == 0) then
         valid = abs(flv1(1)) == abs(flv2(1)) .and. abs(flv1(2)) == abs(flv2(2))
      else
         valid = flv1(1) == flv2(1) .and. flv1(2) == flv2(2) .and. tag1 == tag2
      end if
    end function check_fs_splitting

    function set_dynamic_correction_type (i_flv_alr) result (nlo_corr_type_dyn)
       type(string_t) :: nlo_corr_type_dyn
       type(ftuple_t) :: ftuple_tmp
       integer, intent(in) :: i_flv_alr
       ftuple_tmp = ftuples (region_to_real_index(ftuples, i_flv_alr))%get_ftuple &
                         (region_to_ftuple(i_flv_alr))
       if (ftuple_tmp%qcd_split) then
          nlo_corr_type_dyn = var_str ("QCD")
       else
          nlo_corr_type_dyn = var_str ("EW")
       end if
    end function set_dynamic_correction_type
  end subroutine region_data_init_singular_regions

  module subroutine region_data_find_emitters (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    integer :: alr, j, n_em, em
    integer, dimension(:), allocatable :: em_count
    allocate (em_count(reg_data%n_regions))
    em_count = -1
    n_em = 0

    !!!Count the number of different emitters
    do alr = 1, reg_data%n_regions
       em = reg_data%regions(alr)%emitter
       if (.not. any (em_count == em)) then
          n_em = n_em + 1
          em_count(alr) = em
       end if
    end do

    if (n_em < 1) call msg_fatal ("region_data_find_emitters: No emitters found!")
    reg_data%n_emitters = n_em
    allocate (reg_data%emitters (reg_data%n_emitters))
    reg_data%emitters = -1

    j = 1
    do alr = 1, size (reg_data%regions)
       em = reg_data%regions(alr)%emitter
       if (.not. any (reg_data%emitters == em)) then
          reg_data%emitters(j) = em
          j = j + 1
       end if
    end do
  end subroutine region_data_find_emitters

  module subroutine region_data_find_resonances (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    integer :: alr, j, k, n_res, n_contr
    integer :: res
    integer, dimension(10) :: res_count
    type(resonance_contributors_t), dimension(10) :: contributors_count
    type(resonance_contributors_t) :: contributors
    integer :: i_res, emitter
    logical :: share_emitter
    res_count = -1
    n_res = 0; n_contr = 0

    !!! Count the number of different resonances
    do alr = 1, reg_data%n_regions
       select type (fks_mapping => reg_data%fks_mapping)
       type is (fks_mapping_resonances_t)
          res = fks_mapping%res_map%alr_to_i_res (alr)
          if (.not. any (res_count == res)) then
             n_res = n_res + 1
             res_count(alr) = res
          end if
       end select
    end do

    if (n_res > 0) allocate (reg_data%resonances (n_res))

    j = 1
    select type (fks_mapping => reg_data%fks_mapping)
    type is (fks_mapping_resonances_t)
       do alr = 1, size (reg_data%regions)
          res = fks_mapping%res_map%alr_to_i_res (alr)
          if (.not. any (reg_data%resonances == res)) then
             reg_data%resonances(j) = res
             j = j + 1
          end if
       end do

       allocate (reg_data%alr_to_i_contributor (size (reg_data%regions)))
       do alr = 1, size (reg_data%regions)
          i_res = fks_mapping%res_map%alr_to_i_res (alr)
          emitter = reg_data%regions(alr)%emitter
          call reg_data%get_contributors (i_res, emitter, contributors%c, share_emitter)
          if (.not. share_emitter) cycle
          if (.not. any (contributors_count == contributors)) then
             n_contr = n_contr + 1
             contributors_count(alr) = contributors
          end if
          if (allocated (contributors%c)) deallocate (contributors%c)
       end do
       allocate (reg_data%alr_contributors (n_contr))
       j = 1
       do alr = 1, size (reg_data%regions)
          i_res = fks_mapping%res_map%alr_to_i_res (alr)
          emitter = reg_data%regions(alr)%emitter
          call reg_data%get_contributors (i_res, emitter, contributors%c, share_emitter)
          if (.not. share_emitter) cycle
          if (.not. any (reg_data%alr_contributors == contributors)) then
             reg_data%alr_contributors(j) = contributors
             reg_data%alr_to_i_contributor (alr) = j
             j = j + 1
          else
             do k = 1, size (reg_data%alr_contributors)
                if (reg_data%alr_contributors(k) == contributors) exit
             end do
             reg_data%alr_to_i_contributor (alr) = k
          end if
          if (allocated (contributors%c)) deallocate (contributors%c)
       end do
    end select
    call reg_data%extend_ftuples (n_res)
    call reg_data%set_contributors ()

  end subroutine region_data_find_resonances

  module subroutine region_data_set_i_phs_to_i_con (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    integer :: alr
    integer :: i_res, emitter, i_con, i_phs, i_em
    type(phs_identifier_t), dimension(:), allocatable :: phs_id_tmp
    logical :: share_emitter, phs_exist
    type(resonance_contributors_t) :: contributors
    allocate (phs_id_tmp (reg_data%n_phs))
    if (allocated (reg_data%resonances)) then
       allocate (reg_data%i_phs_to_i_con (reg_data%n_phs))
       do i_em = 1, size (reg_data%emitters)
          emitter = reg_data%emitters(i_em)
          do i_res = 1, size (reg_data%resonances)
             if (reg_data%emitter_is_compatible_with_resonance (i_res, emitter)) then
                alr = find_alr (emitter, i_res)
                if (alr == 0) call msg_fatal ("Could not find requested alpha region!")
                i_con = reg_data%alr_to_i_contributor (alr)
                call reg_data%get_contributors (i_res, emitter, contributors%c, share_emitter)
                if (.not. share_emitter) cycle
                call check_for_phs_identifier &
                   (phs_id_tmp, reg_data%n_in, emitter, contributors%c, phs_exist, i_phs)
                if (phs_id_tmp(i_phs)%emitter < 0) then
                   phs_id_tmp(i_phs)%emitter = emitter
                   allocate (phs_id_tmp(i_phs)%contributors (size (contributors%c)))
                   phs_id_tmp(i_phs)%contributors = contributors%c
                end if
                reg_data%i_phs_to_i_con (i_phs) = i_con
             end if
             if (allocated (contributors%c)) deallocate (contributors%c)
          end do
       end do
    end if
  contains
    function find_alr (emitter, i_res) result (alr)
       integer :: alr
       integer, intent(in) :: emitter, i_res
       integer :: i
       do i = 1, reg_data%n_regions
          if (reg_data%regions(i)%emitter == emitter .and. &
              reg_data%regions(i)%i_res == i_res) then
             alr = i
             return
          end if
       end do
       alr = 0
    end function find_alr
  end subroutine region_data_set_i_phs_to_i_con

  module subroutine region_data_set_alr_to_i_phs &
       (reg_data, phs_identifiers, alr_to_i_phs)
    class(region_data_t), intent(inout) :: reg_data
    type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
    integer, intent(out), dimension(:) :: alr_to_i_phs
    integer :: alr, i_phs
    integer :: emitter, i_res
    type(resonance_contributors_t) :: contributors
    logical :: share_emitter, phs_exist
    do alr = 1, reg_data%n_regions
       associate (region => reg_data%regions(alr))
          emitter = region%emitter
          i_res = region%i_res
          if (i_res /= 0) then
             call reg_data%get_contributors (i_res, emitter, &
                contributors%c, share_emitter)
             if (.not. share_emitter) cycle
          end if
          if (allocated (contributors%c)) then
             call check_for_phs_identifier (phs_identifiers, reg_data%n_in, &
                emitter, contributors%c, phs_exist = phs_exist, i_phs = i_phs)
          else
             call check_for_phs_identifier (phs_identifiers, reg_data%n_in, &
                emitter, phs_exist = phs_exist, i_phs = i_phs)
          end if
          if (.not. phs_exist) &
             call msg_fatal ("phs identifiers are not set up correctly!")
          alr_to_i_phs(alr) = i_phs
       end associate
       if (allocated (contributors%c)) deallocate (contributors%c)
    end do
  end subroutine region_data_set_alr_to_i_phs

  module subroutine region_data_set_contributors (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    integer :: alr, i_res, i_reg, i_con
    integer :: i1, i2, i_em
    integer, dimension(:), allocatable :: contributors
    logical :: share_emitter
    do alr = 1, size (reg_data%regions)
       associate (sregion => reg_data%regions(alr))
         allocate (sregion%i_reg_to_i_con (sregion%nregions))
         do i_reg = 1, sregion%nregions
            call sregion%ftuples(i_reg)%get (i1, i2)
            i_em = get_emitter_index (i1, i2, reg_data%n_legs_real)
            i_res = sregion%ftuples(i_reg)%i_res
            call reg_data%get_contributors &
                 (i_res, i_em, contributors, share_emitter)
            !!! Lookup contributor index
            do i_con = 1, size (reg_data%alr_contributors)
               if (all (reg_data%alr_contributors(i_con)%c == contributors)) &
                    then
                  sregion%i_reg_to_i_con (i_reg) = i_con
                  exit
               end if
            end do
            deallocate (contributors)
         end do
       end associate
    end do
  contains
    function get_emitter_index (i1, i2, n) result (i_em)
      integer :: i_em
      integer, intent(in) :: i1, i2, n
      if (i1 == n) then
         i_em = i2
      else
         i_em = i1
      end if
    end function get_emitter_index
  end subroutine region_data_set_contributors

  module subroutine region_data_extend_ftuples (reg_data, n_res)
    class(region_data_t), intent(inout) :: reg_data
    integer, intent(in) :: n_res
    integer :: alr, n_reg_save
    integer :: i_reg, i_res, i_em, k
    type(ftuple_t), dimension(:), allocatable :: ftuple_save
    integer :: n_new
    do alr = 1, size (reg_data%regions)
       associate (sregion => reg_data%regions(alr))
          n_reg_save = sregion%nregions
          allocate (ftuple_save (n_reg_save))
          ftuple_save = sregion%ftuples
          n_new = count_n_new_ftuples (sregion, n_res)
          deallocate (sregion%ftuples)
          sregion%nregions = n_new
          allocate (sregion%ftuples (n_new))
          k = 1
          do i_res = 1, n_res
             do i_reg = 1, n_reg_save
                associate (ftuple_new => sregion%ftuples(k))
                   i_em = ftuple_save(i_reg)%ireg(1)
                   if (reg_data%emitter_is_in_resonance (i_res, i_em)) then
                      call ftuple_new%set (i_em, ftuple_save(i_reg)%ireg(2))
                      ftuple_new%i_res = i_res
                      ftuple_new%splitting_type = ftuple_save(i_reg)%splitting_type
                      k = k + 1
                   end if
                end associate
             end do
          end do
       end associate
       deallocate (ftuple_save)
    end do
  contains
    function count_n_new_ftuples (sregion, n_res) result (n_new)
      integer :: n_new
      type(singular_region_t), intent(in) :: sregion
      integer, intent(in) :: n_res
      integer :: i_reg, i_res, i_em
      n_new = 0
      do i_reg = 1, sregion%nregions
         do i_res = 1, n_res
            i_em = sregion%ftuples(i_reg)%ireg(1)
            if (reg_data%emitter_is_in_resonance (i_res, i_em)) &
               n_new = n_new + 1
         end do
      end do
    end function count_n_new_ftuples
  end subroutine region_data_extend_ftuples

  module function region_data_get_flavor_indices (reg_data, born) result (i_flv)
    integer, dimension(:), allocatable :: i_flv
    class(region_data_t), intent(in) :: reg_data
    logical, intent(in) :: born
    allocate (i_flv (reg_data%n_regions))
    if (born) then
       i_flv = reg_data%regions%uborn_index
    else
       i_flv = reg_data%regions%real_index
    end if
  end function region_data_get_flavor_indices

  module function region_data_get_matrix_element_index &
       (reg_data, i_reg) result (i_me)
    integer :: i_me
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: i_reg
    i_me = reg_data%regions(i_reg)%real_index
  end function region_data_get_matrix_element_index

  module subroutine region_data_compute_number_of_phase_spaces (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    integer :: i_em, i_res, i_phs
    integer :: emitter
    type(resonance_contributors_t) :: contributors
    integer, parameter :: n_max_phs = 10
    type(phs_identifier_t), dimension(n_max_phs) :: phs_id_tmp
    logical :: share_emitter, phs_exist
    if (allocated (reg_data%resonances)) then
       reg_data%n_phs = 0
       do i_em = 1, size (reg_data%emitters)
          emitter = reg_data%emitters(i_em)
          do i_res = 1, size (reg_data%resonances)
             if (reg_data%emitter_is_compatible_with_resonance (i_res, emitter)) then
                call reg_data%get_contributors &
                     (i_res, emitter, contributors%c, share_emitter)
                if (.not. share_emitter) cycle
                call check_for_phs_identifier (phs_id_tmp, reg_data%n_in, &
                     emitter, contributors%c, phs_exist, i_phs)
                if (.not. phs_exist) then
                   reg_data%n_phs = reg_data%n_phs + 1
                   if (reg_data%n_phs > n_max_phs) call msg_fatal &
                      ("Buffer of phase space identifieres: Too much phase spaces!")
                   call phs_id_tmp(i_phs)%init (emitter, contributors%c)
                end if
             end if
             if (allocated (contributors%c)) deallocate (contributors%c)
          end do
       end do
    else
       reg_data%n_phs = size (remove_duplicates_from_int_array &
            (reg_data%emitters))
    end if
  end subroutine region_data_compute_number_of_phase_spaces

  module function region_data_get_n_phs (reg_data) result (n_phs)
    integer :: n_phs
    class(region_data_t), intent(in) :: reg_data
    n_phs = reg_data%n_phs
  end function region_data_get_n_phs

  module subroutine region_data_set_splitting_info (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    integer :: alr
    do alr = 1, reg_data%n_regions
       call reg_data%regions(alr)%set_splitting_info (reg_data%n_in)
    end do
  end subroutine region_data_set_splitting_info

  module subroutine region_data_init_phs_identifiers (reg_data, phs_id)
    class(region_data_t), intent(in) :: reg_data
    type(phs_identifier_t), intent(out), dimension(:), allocatable :: phs_id
    integer :: i_em, i_res, i_phs
    integer :: emitter
    type(resonance_contributors_t) :: contributors
    logical :: share_emitter, phs_exist
    allocate (phs_id (reg_data%n_phs))
    do i_em = 1, size (reg_data%emitters)
       emitter = reg_data%emitters(i_em)
       if (allocated (reg_data%resonances)) then
          do i_res = 1, size (reg_data%resonances)
             call reg_data%get_contributors (i_res, emitter, &
                  contributors%c, share_emitter)
             if (.not. share_emitter) cycle
             call check_for_phs_identifier &
                (phs_id, reg_data%n_in, emitter, contributors%c, &
                phs_exist, i_phs)
             if (.not. phs_exist) &
                call phs_id(i_phs)%init (emitter, contributors%c)
             if (allocated (contributors%c)) deallocate (contributors%c)
          end do
       else
          call check_for_phs_identifier (phs_id, reg_data%n_in, emitter, &
             phs_exist = phs_exist, i_phs = i_phs)
          if (.not. phs_exist) call phs_id(i_phs)%init (emitter)
       end if
    end do
  end subroutine region_data_init_phs_identifiers

  module subroutine region_data_get_all_ftuples (reg_data, ftuples)
    class(region_data_t), intent(in) :: reg_data
    type(ftuple_t), intent(inout), dimension(:), allocatable :: ftuples
    type(ftuple_t), dimension(:), allocatable :: ftuple_tmp
    integer :: i, j, alr, n_fs
    j = 0
    n_fs = reg_data%n_legs_real - reg_data%n_in
    allocate (ftuple_tmp (n_fs * (n_fs - 1) + 3 * n_fs))
    do i = 1, reg_data%n_regions
       associate (region => reg_data%regions(i))
          do alr = 1, region%nregions
             if (.not. any &
                  (ftuple_equal_ireg (region%ftuples(alr), ftuple_tmp))) then
                j = j + 1
                ftuple_tmp(j) = region%ftuples(alr)
             end if
          end do
       end associate
    end do
    allocate (ftuples (j))
    ftuples = ftuple_tmp(1:j)
    deallocate (ftuple_tmp)
  end subroutine region_data_get_all_ftuples

  module subroutine region_data_write_to_file &
       (reg_data, proc_id, latex, os_data)
    class(region_data_t), intent(inout) :: reg_data
    type(string_t), intent(in) :: proc_id
    logical, intent(in) :: latex
    type(os_data_t), intent(in) :: os_data
    type(string_t) :: filename
    integer :: u
    integer :: status

    if (latex) then
       filename = proc_id // "_fks_regions.tex"
    else
       filename = proc_id // "_fks_regions.out"
    end if
    u = free_unit ()
    open (u, file=char(filename), action = "write", status="replace")
    if (latex) then
       call reg_data%write_latex (u)
       close (u)
       call os_data%build_latex_file &
            (proc_id // "_fks_regions", stat_out = status)
       if (status /= 0) &
            call msg_error (char ("Failed to compile " // filename))
    else
       call reg_data%write (u)
       close (u)
    end if
  end subroutine region_data_write_to_file

  module subroutine region_data_write_latex (reg_data, unit)
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in), optional :: unit
    integer :: i, u
    u = given_output_unit (); if (present (unit)) u = unit
    write (u, "(A)") "\documentclass{article}"
    write (u, "(A)") "\begin{document}"
    write (u, "(A)") "%FKS region data, automatically created by WHIZARD"
    write (u, "(A)") "\begin{table}"
    write (u, "(A)") "\begin{center}"
    write (u, "(A)") "\begin{tabular} {|c|c|c|c|c|c|c|c|}"
    write (u, "(A)") "\hline"
    write (u, "(A)") "$\alpha_r$ & $f_r$ & $i_r$ & $\varepsilon$ & " // &
         "$\varsigma$ & $\mathcal{P}_{\rm{FKS}}$ & $i_b$ & $f_b$ \\"
    write (u, "(A)") "\hline"
    do i = 1, reg_data%n_regions
       call reg_data%regions(i)%write_latex (u)
    end do
    write (u, "(A)") "\hline"
    write (u, "(A)") "\end{tabular}"
    write (u, "(A)") "\caption{List of singular regions}"
    write (u, "(A)") "\begin{description}"
    write (u, "(A)") "\item[$\alpha_r$] Index of the singular region"
    write (u, "(A)") "\item[$f_r$] Real flavor structure"
    write (u, "(A)") "\item[$i_r$] Index of the associated real flavor structure"
    write (u, "(A)") "\item[$\varepsilon$] Emitter"
    write (u, "(A)") "\item[$\varsigma$] Multiplicity"
                          !!! The symbol used by 0908.4272 for multiplicities
    write (u, "(A)") "\item[$\mathcal{P}_{\rm{FKS}}$] The set of singular FKS-pairs"
    write (u, "(A)") "\item[$i_b$] Underlying Born index"
    write (u, "(A)") "\item[$f_b$] Underlying Born flavor structure"
    write (u, "(A)") "\end{description}"
    write (u, "(A)") "\end{center}"
    write (u, "(A)") "\end{table}"
    write (u, "(A)") "\end{document}"
  end subroutine region_data_write_latex

  module subroutine region_data_write (reg_data, unit)
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in), optional :: unit
    integer :: j
    integer :: maxnregions, i_reg_max
    type(string_t) :: flst_title, ftuple_title
    integer :: n_res, u
    u = given_output_unit (unit); if (u < 0) return
    maxnregions = 1; i_reg_max = 1
    do j = 1, reg_data%n_regions
       if (size (reg_data%regions(j)%ftuples) > maxnregions) then
          maxnregions = reg_data%regions(j)%nregions
          i_reg_max = j
       end if
    end do
    flst_title = '(A' // flst_title_format(reg_data%n_legs_real) // ')'
    ftuple_title = '(A' // ftuple_title_format() // ')'
    write (u,'(A,1X,I4)') 'Total number of regions: ', size(reg_data%regions)
    write (u, '(A4)', advance = 'no') ' alr'
    call write_vline (u)
    write (u, char (flst_title), advance = 'no') 'flst_real'
    call write_vline (u)
    write (u, '(A6)', advance = 'no') 'i_real'
    call write_vline (u)
    write (u, '(A3)', advance = 'no') 'em'
    call write_vline (u)
    write (u, '(A3)', advance = 'no') 'mult'
    call write_vline (u)
    write (u, '(A4)', advance = 'no') 'nreg'
    call write_vline (u)
    if (allocated (reg_data%fks_mapping)) then
       select type (fks_mapping => reg_data%fks_mapping)
       type is (fks_mapping_resonances_t)
          write (u, '(A3)', advance = 'no') 'res'
          call write_vline (u)
       end select
    end if
    write (u, char (ftuple_title), advance = 'no') 'ftuples'
    call write_vline (u)
    flst_title = '(A' // flst_title_format(reg_data%n_legs_born) // ')'
    write (u, char (flst_title), advance = 'no') 'flst_born'
    call write_vline (u)
    write (u, '(A7)', advance = 'no') 'i_born'
    call write_vline (u)
    write (u, '(A4)') 'corr'
    do j = 1, reg_data%n_regions
       write (u, '(I4)', advance = 'no') j
       call reg_data%regions(j)%write (u, maxnregions)
    end do
    call write_separator (u)
    if (allocated (reg_data%fks_mapping)) then
       select type (fks_mapping => reg_data%fks_mapping)
       type is (fks_mapping_resonances_t)
          write (u, '(A)')
          write (u, '(A)') &
               "The FKS regions are combined with resonance information: "
          n_res = size (fks_mapping%res_map%res_histories)
          write (u, '(A,1X,I1)') "Number of QCD resonance histories: ", n_res
          do j = 1, n_res
             write (u, '(A,1X,I1)') "i_res = ", j
             call fks_mapping%res_map%res_histories(j)%write (u)
             call write_separator (u)
          end do
       end select
    end if

  contains

    function flst_title_format (n) result (frmt)
      integer, intent(in) :: n
      type(string_t) :: frmt
      character(len=2) :: frmt_char
      write (frmt_char, '(I2)') 4 * n + 1
      frmt = var_str (frmt_char)
    end function flst_title_format

    function ftuple_title_format () result (frmt)
      type(string_t) :: frmt
      integer :: n_ftuple_char
      !!! An ftuple (x,x) consists of five characters. In the string, they
      !!! are separated by maxregions - 1 commas. In total these are
      !!! 5 * maxnregions + maxnregions - 1 = 6 * maxnregions - 1 characters.
      !!! The {} brackets at add two additional characters.
      n_ftuple_char = 6 * maxnregions + 1
      !!! If there are resonances, each ftuple with a resonance adds a ";x"
      !!! to the ftuple
      n_ftuple_char = n_ftuple_char + &
           2 * count (reg_data%regions(i_reg_max)%ftuples%i_res > 0)
      !!! Pseudo-ISR regions are denoted with a * at the end
      n_ftuple_char = n_ftuple_char + &
           count (reg_data%regions(i_reg_max)%ftuples%pseudo_isr)
      frmt = str (n_ftuple_char)
    end function ftuple_title_format

  end subroutine region_data_write

  subroutine write_vline (u)
    integer, intent(in) :: u
    character(len=10), parameter :: sep_format = "(1X,A2,1X)"
    write (u, sep_format, advance = 'no') '||'
  end subroutine write_vline

  module subroutine region_data_assign (reg_data_out, reg_data_in)
    type(region_data_t), intent(out) :: reg_data_out
    type(region_data_t), intent(in) :: reg_data_in
    integer :: i
    if (allocated (reg_data_in%regions)) then
       allocate (reg_data_out%regions (size (reg_data_in%regions)))
       do i = 1, size (reg_data_in%regions)
          reg_data_out%regions(i) = reg_data_in%regions(i)
       end do
    else
       call msg_warning &
            ("Copying region data without allocated singular regions!")
    end if
    if (allocated (reg_data_in%flv_born)) then
       allocate (reg_data_out%flv_born (size (reg_data_in%flv_born)))
       do i = 1, size (reg_data_in%flv_born)
          reg_data_out%flv_born(i) = reg_data_in%flv_born(i)
       end do
    else
       call msg_warning &
            ("Copying region data without allocated born flavor structure!")
    end if
    if (allocated (reg_data_in%flv_real)) then
       allocate (reg_data_out%flv_real (size (reg_data_in%flv_real)))
       do i = 1, size (reg_data_in%flv_real)
          reg_data_out%flv_real(i) = reg_data_in%flv_real(i)
       end do
    else
       call msg_warning ("Copying region data without allocated real flavor structure!")
    end if
    if (allocated (reg_data_in%emitters)) then
       allocate (reg_data_out%emitters (size (reg_data_in%emitters)))
       do i = 1, size (reg_data_in%emitters)
          reg_data_out%emitters(i) = reg_data_in%emitters(i)
       end do
    else
       call msg_warning ("Copying region data without allocated emitters!")
    end if
    reg_data_out%n_regions = reg_data_in%n_regions
    reg_data_out%n_emitters = reg_data_in%n_emitters
    reg_data_out%n_flv_born = reg_data_in%n_flv_born
    reg_data_out%n_flv_real = reg_data_in%n_flv_real
    reg_data_out%n_in = reg_data_in%n_in
    reg_data_out%n_legs_born = reg_data_in%n_legs_born
    reg_data_out%n_legs_real = reg_data_in%n_legs_real
    if (allocated (reg_data_in%fks_mapping)) then
       select type (fks_mapping_in => reg_data_in%fks_mapping)
       type is (fks_mapping_default_t)
          allocate (fks_mapping_default_t :: reg_data_out%fks_mapping)
          select type (fks_mapping_out => reg_data_out%fks_mapping)
          type is (fks_mapping_default_t)
             fks_mapping_out = fks_mapping_in
          end select
       type is (fks_mapping_resonances_t)
          allocate (fks_mapping_resonances_t :: reg_data_out%fks_mapping)
          select type (fks_mapping_out => reg_data_out%fks_mapping)
          type is (fks_mapping_resonances_t)
             fks_mapping_out = fks_mapping_in
          end select
       end select
    else
       call msg_warning ("Copying region data without allocated FKS regions!")
    end if
    if (allocated (reg_data_in%resonances)) then
       allocate (reg_data_out%resonances (size (reg_data_in%resonances)))
       reg_data_out%resonances = reg_data_in%resonances
    end if
    reg_data_out%n_phs = reg_data_in%n_phs
    if (allocated (reg_data_in%alr_contributors)) then
       allocate (reg_data_out%alr_contributors (size (reg_data_in%alr_contributors)))
       reg_data_out%alr_contributors = reg_data_in%alr_contributors
    end if
    if (allocated (reg_data_in%alr_to_i_contributor)) then
       allocate (reg_data_out%alr_to_i_contributor &
          (size (reg_data_in%alr_to_i_contributor)))
       reg_data_out%alr_to_i_contributor = reg_data_in%alr_to_i_contributor
    end if
  end subroutine region_data_assign

  function region_to_real_index (list, i) result(index)
    type(ftuple_list_t), intent(in), dimension(:), allocatable :: list
    integer, intent(in) :: i
    integer, dimension(:), allocatable :: nreg
    integer :: index, j
    allocate (nreg (0))
    index = 0
    do j = 1, size (list)
       nreg = [nreg, sum (list(:j)%get_n_tuples ())]
       if (j == 1) then
          if (i <= nreg(j)) then
             index = j
             exit
          end if
       else
          if (i > nreg(j - 1) .and. i <= nreg(j)) then
             index = j
             exit
          end if
       end if
    end do
  end function region_to_real_index

  function create_alr (flv1, n_in, i_em, i_rad) result(flv2)
    type(flv_structure_t), intent(in) :: flv1
    integer, intent(in) :: n_in
    integer, intent(in) :: i_em, i_rad
    type(flv_structure_t) :: flv2
    integer :: n
    n = size (flv1%flst)
    allocate (flv2%flst (n), flv2%tag (n))
    flv2%nlegs = n
    flv2%n_in = n_in
    if (i_em > n_in) then
       flv2%flst(1 : n_in) = flv1%flst(1 : n_in)
       flv2%flst(n - 1) = flv1%flst(i_em)
       flv2%flst(n) = flv1%flst(i_rad)
       flv2%tag(1 : n_in) = flv1%tag(1 : n_in)
       flv2%tag(n - 1) = flv1%tag(i_em)
       flv2%tag(n) = flv1%tag(i_rad)
       call fill_remaining_flavors (n_in, .true.)
    else
       flv2%flst(1 : n_in) = flv1%flst(1 : n_in)
       flv2%flst(n) = flv1%flst(i_rad)
       flv2%tag(1 : n_in) = flv1%tag(1 : n_in)
       flv2%tag(n) = flv1%tag(i_rad)
       call fill_remaining_flavors (n_in, .false.)
    end if
    call flv2%compute_prt_symm_fs (flv2%n_in)
  contains
    subroutine fill_remaining_flavors (n_in, final_final)
      integer, intent(in) :: n_in
      logical, intent(in) :: final_final
      integer :: i, j
      logical :: check
      j = n_in + 1
      do i = n_in + 1, n
         if (final_final) then
            check = (i /= i_em .and. i /= i_rad)
         else
            check = (i /= i_rad)
         end if
         if (check) then
            flv2%flst(j) = flv1%flst(i)
            flv2%tag(j) = flv1%tag(i)
            j = j + 1
         end if
      end do
    end subroutine fill_remaining_flavors
  end function create_alr

  module function region_data_has_pseudo_isr (reg_data) result (flag)
    logical :: flag
    class(region_data_t), intent(in) :: reg_data
    flag = any (reg_data%regions%pseudo_isr)
  end function region_data_has_pseudo_isr

  module subroutine region_data_check_consistency (reg_data, fail_fatal, unit)
    class(region_data_t), intent(in) :: reg_data
    logical, intent(in) :: fail_fatal
    integer, intent(in), optional :: unit
    integer :: u
    integer :: i_reg, alr
    integer :: i1, f1, f2
    logical :: undefined_ftuples, same_ftuple_indices, valid_splitting
    logical, dimension(4) :: no_fail
    u = given_output_unit(unit); if (u < 0) return
    no_fail = .true.
    call msg_message ("Check that no negative ftuple indices occur", unit = u)
    do i_reg = 1, reg_data%n_regions
       if (any (reg_data%regions(i_reg)%ftuples%has_negative_elements ())) then
          !!! This error is so severe that we stop immediately
          call msg_fatal ("Negative ftuple indices!")
       end if
    end do
    call msg_message ("Success!", unit = u)
    call msg_message ("Check that there is no ftuple with identical elements", unit = u)
    do i_reg = 1, reg_data%n_regions
       if (any (reg_data%regions(i_reg)%ftuples%has_identical_elements ())) then
          !!! This error is so severe that we stop immediately
          call msg_fatal ("Identical ftuple indices!")
       end if
    end do
    call msg_message ("Success!", unit = u)
    call msg_message ("Check that there are no duplicate ftuples in a region", unit = u)
    do i_reg = 1, reg_data%n_regions
       if (reg_data%regions(i_reg)%has_identical_ftuples ()) then
          if (no_fail(1)) then
             call msg_error ("FAIL: ", unit = u)
             no_fail(1) = .false.
          end if
          write (u, '(A,1x,I3)') 'i_reg:', i_reg
       end if
    end do
    if (no_fail(1)) call msg_message ("Success!", unit = u)
    call msg_message ("Check that ftuples add up to a valid splitting", unit = u)
    do i_reg = 1, reg_data%n_regions
       do alr = 1, reg_data%regions(i_reg)%nregions
          associate (region => reg_data%regions(i_reg))
             i1 = region%ftuples(alr)%ireg(1)
             if (i1 == 0) i1 = 1 !!! Gluon emission from both initial-state particles
             f1 = region%flst_real%flst(i1)
             f2 = region%flst_real%flst(region%ftuples(alr)%ireg(2))
             ! Flip PDG sign of IS fermions to allow a q -> g q splitting
             ! in which the ftuple has the flavors (q,q).
             if (i1 <= reg_data%n_in .and. is_fermion(f1)) then
                f1 = -f1
             end if
             valid_splitting = f1 + f2 == 0 &
                  .or. (is_gluon(f1) .and. is_gluon(f2)) &
                  .or. (is_massive_vector(f1) .and. is_photon(f2)) &
                  .or. is_fermion_vector_splitting (f1, f2)
             if (.not. valid_splitting) then
                if (no_fail(2)) then
                   call msg_error ("FAIL: ", unit = u)
                   no_fail(2) = .false.
                end if
                write (u, '(A,1x,I3)') 'i_reg:', i_reg
                exit
             end if
          end associate
       end do
    end do
    if (no_fail(2)) call msg_message ("Success!", unit = u)
    call msg_message ("Check that at least one ftuple contains the emitter", unit = u)
    do i_reg = 1, reg_data%n_regions
       associate (region => reg_data%regions(i_reg))
          if (.not. any (region%emitter == region%ftuples%ireg(1))) then
             if (no_fail(3)) then
                call msg_error ("FAIL: ", unit = u)
                no_fail(3) = .false.
             end if
             write (u, '(A,1x,I3)') 'i_reg:', i_reg
          end if
       end associate
    end do
    if (no_fail(3)) call msg_message ("Success!", unit = u)
    call msg_message ("Check that each region has at least one ftuple &
         &with index n + 1", unit = u)
    do i_reg = 1, reg_data%n_regions
       if (.not. any (reg_data%regions(i_reg)%ftuples%ireg(2) == reg_data%n_legs_real)) then
          if (no_fail(4)) then
             call msg_error ("FAIL: ", unit = u)
             no_fail(4) = .false.
          end if
          write (u, '(A,1x,I3)') 'i_reg:', i_reg
       end if
    end do
    if (no_fail(4)) call msg_message ("Success!", unit = u)
    if (.not. all (no_fail)) &
         call abort_with_message ("Stop due to inconsistent region data!")

  contains
    subroutine abort_with_message (msg)
      character(len=*), intent(in) :: msg
      if (fail_fatal) then
         call msg_fatal (msg)
      else
         call msg_error (msg, unit = u)
      end if
    end subroutine abort_with_message

    function is_fermion_vector_splitting (pdg_1, pdg_2) result (value)
      logical :: value
      integer, intent(in) :: pdg_1, pdg_2
      value = (is_fermion (pdg_1) .and. is_massless_vector (pdg_2)) .or. &
           (is_fermion (pdg_2) .and. is_massless_vector (pdg_1))
    end function
  end subroutine region_data_check_consistency

  module function region_data_requires_spin_correlations &
       (reg_data) result (flag)
    class(region_data_t), intent(in) :: reg_data
    logical :: flag
    integer :: alr
    flag = .false.
    do alr = 1, reg_data%n_regions
       flag = reg_data%regions(alr)%sc_required
       if (flag) return
    end do
  end function region_data_requires_spin_correlations

  module function region_data_born_to_real_symm_factor_fs &
       (reg_data, alr) result (factor)
    class(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: alr
    real(default) :: factor
    associate (flv_real => reg_data%regions(alr)%flst_real, &
            flv_uborn => reg_data%regions(alr)%flst_uborn)
       factor = flv_real%prt_symm_fs / flv_uborn%prt_symm_fs
    end associate
  end function region_data_born_to_real_symm_factor_fs

  module subroutine region_data_final (reg_data)
    class(region_data_t), intent(inout) :: reg_data
    if (allocated (reg_data%regions)) deallocate (reg_data%regions)
    if (allocated (reg_data%flv_born)) deallocate (reg_data%flv_born)
    if (allocated (reg_data%flv_real)) deallocate (reg_data%flv_real)
    if (allocated (reg_data%emitters)) deallocate (reg_data%emitters)
    if (allocated (reg_data%fks_mapping)) deallocate (reg_data%fks_mapping)
    if (allocated (reg_data%resonances)) deallocate (reg_data%resonances)
    if (allocated (reg_data%alr_contributors)) &
         deallocate (reg_data%alr_contributors)
    if (allocated (reg_data%alr_to_i_contributor)) &
         deallocate (reg_data%alr_to_i_contributor)
  end subroutine region_data_final

  module subroutine fks_mapping_default_set_parameter &
       (map, n_in, dij_exp1, dij_exp2)
    class(fks_mapping_default_t), intent(inout) :: map
    integer, intent(in) :: n_in
    real(default), intent(in) :: dij_exp1, dij_exp2
    map%n_in = n_in
    map%exp_1 = dij_exp1
    map%exp_2 = dij_exp2
  end subroutine fks_mapping_default_set_parameter

  module function fks_mapping_default_dij (map, p, i, j, i_con) result (d)
    real(default) :: d
    class(fks_mapping_default_t), intent(in) :: map
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in) :: i, j
    integer, intent(in), optional :: i_con
    d = zero
    if (map%pseudo_isr) then
       d = dij_threshold_gluon_from_top (i, j, p, map%exp_1)
    else if (i > map%n_in .and. j > map%n_in) then
       d = dij_fsr (p(i), p(j), map%exp_1)
    else
       d = dij_isr (map%n_in, i, j, p, map%exp_2)
    end if
  contains

    function dij_fsr (p1, p2, expo) result (d_ij)
      real(default) :: d_ij
      type(vector4_t), intent(in) :: p1, p2
      real(default), intent(in) :: expo
      real(default) :: E1, E2
      E1 = p1%p(0); E2 = p2%p(0)
      d_ij = (two * p1 * p2 * E1 * E2 / (E1 + E2)**2)**expo
    end function dij_fsr

    function dij_threshold_gluon_from_top (i, j, p, expo) result (d_ij)
      real(default) :: d_ij
      integer, intent(in) :: i, j
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in) :: expo
      type(vector4_t) :: p_top
      if (i == THR_POS_B) then
         p_top = p(THR_POS_WP) + p(THR_POS_B)
      else
         p_top = p(THR_POS_WM) + p(THR_POS_BBAR)
      end if
      d_ij = dij_fsr (p_top, p(j), expo)
    end function dij_threshold_gluon_from_top

    function dij_isr (n_in, i, j, p, expo) result (d_ij)
      real(default) :: d_ij
      integer, intent(in) :: n_in, i, j
      type(vector4_t), intent(in), dimension(:) :: p
      real(default), intent(in) :: expo
      real(default) :: E, y
      select case (n_in)
      case (1)
         call get_emitter_variables (1, i, j, p, E, y)
         d_ij = (E**2 * (one - y**2))**expo
      case (2)
         if ((i == 0 .and. j > 2) .or. (j == 0 .and. i > 2)) then
            call get_emitter_variables (0, i, j, p, E, y)
            d_ij = (E**2 * (one - y**2))**expo
         else if ((i == 1 .and. j > 2) .or. (j == 1 .and. i > 2)) then
            call get_emitter_variables (1, i, j, p, E, y)
            d_ij = (two * E**2 * (one - y))**expo
         else if ((i == 2 .and. j > 2) .or. (j == 2 .and. i > 2)) then
            call get_emitter_variables (2, i, j, p, E, y)
            d_ij = (two * E**2 * (one + y))**expo
         end if
      end select
    end function dij_isr

    subroutine get_emitter_variables (i_check, i, j, p, E, y)
       integer, intent(in) :: i_check, i, j
       type(vector4_t), intent(in), dimension(:) :: p
       real(default), intent(out) :: E, y
       if (j == i_check) then
           E = energy (p(i))
           y = polar_angle_ct (p(i))
       else
           E = energy (p(j))
           y = polar_angle_ct(p(j))
       end if
    end subroutine get_emitter_variables

  end function fks_mapping_default_dij

  module subroutine fks_mapping_default_compute_sumdij (map, sregion, p_real)
    class(fks_mapping_default_t), intent(inout) :: map
    type(singular_region_t), intent(in) :: sregion
    type(vector4_t), intent(in), dimension(:) :: p_real
    real(default) :: d
    integer :: i_ftuple, i, j

    associate (ftuples => sregion%ftuples)
      d = zero
      do i_ftuple = 1, sregion%nregions
         call ftuples(i_ftuple)%get (i, j)
         map%pseudo_isr = ftuples(i_ftuple)%pseudo_isr
         d = d + one / map%dij (p_real, i, j)
      end do
    end associate
    map%sumdij = d
  end subroutine fks_mapping_default_compute_sumdij

  module function fks_mapping_default_svalue &
       (map, p, i, j, i_res) result (value)
    real(default) :: value
    class(fks_mapping_default_t), intent(in) :: map
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in) :: i, j
    integer, intent(in), optional :: i_res
    value = one / (map%dij (p, i, j) * map%sumdij)
  end function fks_mapping_default_svalue

  module function fks_mapping_default_dij_soft &
       (map, p_born, p_soft, em, i_con) result (d)
    real(default) :: d
    class(fks_mapping_default_t), intent(in) :: map
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(in) :: p_soft
    integer, intent(in) :: em
    integer, intent(in), optional :: i_con
    if (map%pseudo_isr) then
       d = dij_soft_threshold_gluon_from_top (em, p_born, p_soft, map%exp_1)
    else if (em <= map%n_in) then
       d = dij_soft_isr (map%n_in, p_soft, map%exp_2)
    else
       d = dij_soft_fsr (p_born(em), p_soft, map%exp_1)
    end if
  contains

    function dij_soft_threshold_gluon_from_top &
         (em, p_born, p_soft, expo) result (dij_soft)
      real(default) :: dij_soft
      integer, intent(in) :: em
      type(vector4_t), intent(in), dimension(:) :: p_born
      type(vector4_t), intent(in) :: p_soft
      real(default), intent(in) :: expo
      type(vector4_t) :: p_top
      if (em == THR_POS_B) then
         p_top = p_born(THR_POS_WP) + p_born(THR_POS_B)
      else
         p_top = p_born(THR_POS_WM) + p_born(THR_POS_BBAR)
      end if
      dij_soft = dij_soft_fsr (p_top, p_soft, expo)
    end function dij_soft_threshold_gluon_from_top

    function dij_soft_fsr (p_em, p_soft, expo) result (dij_soft)
      real(default) :: dij_soft
      type(vector4_t), intent(in) :: p_em, p_soft
      real(default), intent(in) :: expo
      dij_soft = (two * p_em * p_soft / p_em%p(0))**expo
    end function dij_soft_fsr

    function dij_soft_isr (n_in, p_soft, expo) result (dij_soft)
       real(default) :: dij_soft
       integer, intent(in) :: n_in
       type(vector4_t), intent(in) :: p_soft
       real(default), intent(in) :: expo
       real(default) :: y
       y = polar_angle_ct (p_soft)
       select case (n_in)
       case (1)
          dij_soft = one - y**2
       case (2)
          select case (em)
          case (0)
             dij_soft = one - y**2
          case (1)
             dij_soft = two * (one - y)
          case (2)
             dij_soft = two * (one + y)
          case default
             dij_soft = zero
             call msg_fatal ("fks_mappings_default_dij_soft: n_in > 2")
          end select
       case default
          dij_soft = zero
          call msg_fatal ("fks_mappings_default_dij_soft: n_in > 2")
       end select
       dij_soft = dij_soft**expo
    end function dij_soft_isr
  end function fks_mapping_default_dij_soft

  module subroutine fks_mapping_default_compute_sumdij_soft &
       (map, sregion, p_born, p_soft)
    class(fks_mapping_default_t), intent(inout) :: map
    type(singular_region_t), intent(in) :: sregion
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(in) :: p_soft
    real(default) :: d
    integer :: i_ftuple, i, j
    integer :: nlegs
    d = zero
    nlegs = size (sregion%flst_real%flst)
    associate (ftuples => sregion%ftuples)
      do i_ftuple = 1, sregion%nregions
         call ftuples(i_ftuple)%get (i ,j)
         if (j == nlegs) then
            map%pseudo_isr = ftuples(i_ftuple)%pseudo_isr
            d = d + one / map%dij_soft (p_born, p_soft, i)
         end if
      end do
    end associate
    map%sumdij_soft = d
  end subroutine fks_mapping_default_compute_sumdij_soft

  module function fks_mapping_default_svalue_soft &
       (map, p_born, p_soft, em, i_res) result (value)
    real(default) :: value
    class(fks_mapping_default_t), intent(in) :: map
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(in) :: p_soft
    integer, intent(in) :: em
    integer, intent(in), optional :: i_res
    value = one / (map%sumdij_soft * map%dij_soft (p_born, p_soft, em))
  end function fks_mapping_default_svalue_soft

  module subroutine fks_mapping_default_assign (fks_map_out, fks_map_in)
    type(fks_mapping_default_t), intent(out) :: fks_map_out
    type(fks_mapping_default_t), intent(in) :: fks_map_in
    fks_map_out%exp_1 = fks_map_in%exp_1
    fks_map_out%exp_2 = fks_map_in%exp_2
    fks_map_out%n_in = fks_map_in%n_in
  end subroutine fks_mapping_default_assign

  module function fks_mapping_resonances_dij (map, p, i, j, i_con) result (d)
    real(default) :: d
    class(fks_mapping_resonances_t), intent(in) :: map
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in) :: i, j
    integer, intent(in), optional :: i_con
    real(default) :: E1, E2
    integer :: ii_con
    if (present (i_con)) then
       ii_con = i_con
    else
       call msg_fatal ("Resonance mappings require resonance index as input!")
    end if
    d = 0
    if (i /= j) then
       if (i > 2 .and. j > 2) then
          associate (p_res => map%res_map%p_res (ii_con))
             E1 = p(i) * p_res
             E2 = p(j) * p_res
             d = two * p(i) * p(j) * E1 * E2 / (E1 + E2)**2
          end associate
       else
          call msg_fatal ("Resonance mappings are not implemented for ISR")
       end if
    end if
  end function fks_mapping_resonances_dij

  module subroutine fks_mapping_resonances_compute_sumdij (map, sregion, p_real)
    class(fks_mapping_resonances_t), intent(inout) :: map
    type(singular_region_t), intent(in) :: sregion
    type(vector4_t), intent(in), dimension(:) :: p_real
    real(default) :: d, pfr
    integer :: i_res, i_reg, i, j, i_con
    integer :: nlegreal

    nlegreal = size (p_real)
    d = zero
    do i_reg = 1, sregion%nregions
       associate (ftuple => sregion%ftuples(i_reg))
          call ftuple%get (i, j)
          i_res = ftuple%i_res
       end associate
       pfr = map%res_map%get_resonance_value (i_res, p_real, nlegreal)
       i_con = sregion%i_reg_to_i_con (i_reg)
       d = d + pfr / map%dij (p_real, i, j, i_con)
    end do
    map%sumdij = d
  end subroutine fks_mapping_resonances_compute_sumdij

  module function fks_mapping_resonances_svalue &
       (map, p, i, j, i_res) result (value)
    real(default) :: value
    class(fks_mapping_resonances_t), intent(in) :: map
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in) :: i, j
    integer, intent(in), optional :: i_res
    real(default) :: pfr
    integer :: i_gluon
    i_gluon = size (p)
    pfr = map%res_map%get_resonance_value (i_res, p, i_gluon)
    value = pfr / (map%dij (p, i, j, map%i_con) * map%sumdij)
  end function fks_mapping_resonances_svalue

  module function fks_mapping_resonances_get_resonance_weight &
       (map, alr, p) result (pfr)
    real(default) :: pfr
    class(fks_mapping_resonances_t), intent(in) :: map
    integer, intent(in) :: alr
    type(vector4_t), intent(in), dimension(:) :: p
    pfr = map%res_map%get_weight (alr, p)
  end function fks_mapping_resonances_get_resonance_weight

  module function fks_mapping_resonances_dij_soft &
       (map, p_born, p_soft, em, i_con) result (d)
    real(default) :: d
    class(fks_mapping_resonances_t), intent(in) :: map
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(in) :: p_soft
    integer, intent(in) :: em
    integer, intent(in), optional :: i_con
    real(default) :: E1, E2
    integer :: ii_con
    type(vector4_t) :: pb
    if (present (i_con)) then
       ii_con = i_con
    else
       call msg_fatal ("fks_mapping_resonances requires resonance index")
    end if
    associate (p_res => map%res_map%p_res(ii_con))
      pb = p_born(em)
      E1 = pb * p_res
      E2 = p_soft * p_res
      d = two * pb * p_soft * E1 * E2 / E1**2
    end associate
  end function fks_mapping_resonances_dij_soft

  module subroutine fks_mapping_resonances_compute_sumdij_soft &
       (map, sregion, p_born, p_soft)
    class(fks_mapping_resonances_t), intent(inout) :: map
    type(singular_region_t), intent(in) :: sregion
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(in) :: p_soft
    real(default) :: d
    real(default) :: pfr
    integer :: i_res, i, j, i_reg, i_con
    integer :: nlegs

    d = zero
    nlegs = size (sregion%flst_real%flst)
    do i_reg = 1, sregion%nregions
       associate (ftuple => sregion%ftuples(i_reg))
          call ftuple%get(i, j)
          i_res = ftuple%i_res
       end associate
       pfr = map%res_map%get_resonance_value (i_res, p_born)
       i_con = sregion%i_reg_to_i_con (i_reg)
       if (j == nlegs) d = d + pfr / map%dij_soft (p_born, p_soft, i, i_con)
    end do
    map%sumdij_soft = d
  end subroutine fks_mapping_resonances_compute_sumdij_soft

  module function fks_mapping_resonances_svalue_soft &
       (map, p_born, p_soft, em, i_res) result (value)
    real(default) :: value
    class(fks_mapping_resonances_t), intent(in) :: map
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(in) :: p_soft
    integer, intent(in) :: em
    integer, intent(in), optional :: i_res
    real(default) :: pfr
    pfr = map%res_map%get_resonance_value (i_res, p_born)
    value = pfr / (map%sumdij_soft * map%dij_soft &
         (p_born, p_soft, em, map%i_con))
  end function fks_mapping_resonances_svalue_soft

  module subroutine fks_mapping_resonances_set_resonance_momentum (map, p)
    class(fks_mapping_resonances_t), intent(inout) :: map
    type(vector4_t), intent(in) :: p
    map%res_map%p_res = p
  end subroutine fks_mapping_resonances_set_resonance_momentum

  module subroutine fks_mapping_resonances_set_resonance_momenta (map, p)
    class(fks_mapping_resonances_t), intent(inout) :: map
    type(vector4_t), intent(in), dimension(:) :: p
    map%res_map%p_res = p
  end subroutine fks_mapping_resonances_set_resonance_momenta

  module subroutine fks_mapping_resonances_assign (fks_map_out, fks_map_in)
    type(fks_mapping_resonances_t), intent(out) :: fks_map_out
    type(fks_mapping_resonances_t), intent(in) :: fks_map_in
    fks_map_out%exp_1 = fks_map_in%exp_1
    fks_map_out%exp_2 = fks_map_in%exp_2
    fks_map_out%res_map = fks_map_in%res_map
  end subroutine fks_mapping_resonances_assign

  module function create_resonance_histories_for_threshold &
       () result (res_history)
    type(resonance_history_t) :: res_history
    res_history%n_resonances = 2
    allocate (res_history%resonances (2))
    allocate (res_history%resonances(1)%contributors%c(2))
    allocate (res_history%resonances(2)%contributors%c(2))
    res_history%resonances(1)%contributors%c = [THR_POS_WP, THR_POS_B]
    res_history%resonances(2)%contributors%c = [THR_POS_WM, THR_POS_BBAR]
  end function create_resonance_histories_for_threshold

  module subroutine setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, &
       nlo_corr_type, alpha_pow, alphas_pow)
    integer, intent(in) :: n_in, alpha_pow, alphas_pow
    integer, intent(in), dimension(:,:) :: flv_born, flv_real
    type(string_t), intent(in) :: nlo_corr_type
    type(region_data_t), intent(out) :: reg_data
    type(model_t), pointer :: test_model => null ()
    call create_test_model (var_str ("SM"), test_model)
    call test_model%set_real (var_str ("me"), 0._default)
    call test_model%set_real (var_str ("mmu"), 0._default)
    call test_model%set_real (var_str ("mtau"), 0._default)
    call test_model%set_real (var_str ("ms"), 0._default)
    call test_model%set_real (var_str ("mc"), 0._default)
    call test_model%set_real (var_str ("mb"), 0._default)
    call reg_data%init (n_in, test_model, flv_born, flv_real, nlo_corr_type, alpha_pow, &
         alphas_pow)
  end subroutine setup_region_data_for_test


end submodule fks_regions_s

