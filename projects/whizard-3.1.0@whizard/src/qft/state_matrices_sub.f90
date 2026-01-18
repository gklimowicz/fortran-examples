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

submodule (state_matrices) state_matrices_s

  use constants, only: zero
  use format_utils, only: pac_fmt
  use format_defs, only: FMT_17, FMT_19
  use io_units
  use diagnostics
  use sorting

  implicit none

contains

  pure recursive subroutine node_delete_offspring (node)
    type(node_t), pointer :: node
    type(node_t), pointer :: child
    child => node%child_first
    do while (associated (child))
       node%child_first => node%child_first%next
       call node_delete_offspring (child)
       deallocate (child)
       child => node%child_first
    end do
    node%child_last => null ()
  end subroutine node_delete_offspring

  pure subroutine node_delete (node)
    type(node_t), pointer :: node
    call node_delete_offspring (node)
    if (associated (node%previous)) then
       node%previous%next => node%next
    else if (associated (node%parent)) then
       node%parent%child_first => node%next
    end if
    if (associated (node%next)) then
       node%next%previous => node%previous
    else if (associated (node%parent)) then
       node%parent%child_last => node%previous
    end if
    deallocate (node)
  end subroutine node_delete

  subroutine node_append_child (node, child)
    type(node_t), target, intent(inout) :: node
    type(node_t), pointer :: child
    allocate (child)
    if (associated (node%child_last)) then
       node%child_last%next => child
       child%previous => node%child_last
    else
       node%child_first => child
    end if
    node%child_last => child
    child%parent => node
  end subroutine node_append_child

  subroutine node_write (node, me_array, verbose, unit, col_verbose, testflag)
    type(node_t), intent(in) :: node
    complex(default), dimension(:), intent(in), optional :: me_array
    logical, intent(in), optional :: verbose, col_verbose, testflag
    integer, intent(in), optional :: unit
    logical :: verb
    integer :: u
    character(len=7) :: fmt
    call pac_fmt (fmt, FMT_19, FMT_17, testflag)
    verb = .false.;  if (present (verbose)) verb = verbose
    u = given_output_unit (unit);  if (u < 0)  return
    call node%qn%write (u, col_verbose)
    if (node%me_index /= 0) then
       write (u, "(A,I0,A)", advance="no")  " => ME(", node%me_index, ")"
       if (present (me_array)) then
          write (u, "(A)", advance="no")  " = "
          write (u, "('('," // fmt // ",','," // fmt // ",')')", &
               advance="no") pacify_complex (me_array(node%me_index))
       end if
    end if
    write (u, *)
    if (verb) then
       call ptr_write ("parent     ", node%parent)
       call ptr_write ("child_first", node%child_first)
       call ptr_write ("child_last ", node%child_last)
       call ptr_write ("next       ", node%next)
       call ptr_write ("previous   ", node%previous)
    end if
  contains
    subroutine ptr_write (label, node)
      character(*), intent(in) :: label
      type(node_t), pointer :: node
      if (associated (node)) then
         write (u, "(10x,A,1x,'->',1x)", advance="no") label
         call node%qn%write (u, col_verbose)
         write (u, *)
      end if
    end subroutine ptr_write
  end subroutine node_write

  recursive subroutine node_write_rec (node, me_array, verbose, &
        indent, unit, col_verbose, testflag)
    type(node_t), intent(in), target :: node
    complex(default), dimension(:), intent(in), optional :: me_array
    logical, intent(in), optional :: verbose, col_verbose, testflag
    integer, intent(in), optional :: indent
    integer, intent(in), optional :: unit
    type(node_t), pointer :: current
    logical :: verb
    integer :: i, u
    verb = .false.;  if (present (verbose))  verb = verbose
    i = 0;  if (present (indent)) i = indent
    u = given_output_unit (unit);  if (u < 0)  return
    current => node%child_first
    do while (associated (current))
       write (u, "(A)", advance="no")  repeat (" ", i)
       call node_write (current, me_array, verbose = verb, &
            unit = u, col_verbose = col_verbose, testflag = testflag)
       call node_write_rec (current, me_array, verbose = verb, &
            indent = i + 2, unit = u, col_verbose = col_verbose, testflag = testflag)
       current => current%next
    end do
  end subroutine node_write_rec

  recursive subroutine node_write_raw_rec (node, u)
    type(node_t), intent(in), target :: node
    integer, intent(in) :: u
    logical :: associated_child_first, associated_next
    call node%qn%write_raw (u)
    associated_child_first = associated (node%child_first)
    write (u) associated_child_first
    associated_next = associated (node%next)
    write (u) associated_next
    if (associated_child_first) then
       call node_write_raw_rec (node%child_first, u)
    else
       write (u)  node%me_index
       write (u)  node%me
    end if
    if (associated_next) then
       call node_write_raw_rec (node%next, u)
    end if
  end subroutine node_write_raw_rec

  recursive subroutine node_read_raw_rec (node, u, parent, iostat)
    type(node_t), intent(out), target :: node
    integer, intent(in) :: u
    type(node_t), intent(in), optional, target :: parent
    integer, intent(out), optional :: iostat
    logical :: associated_child_first, associated_next
    type(node_t), pointer :: child
    call node%qn%read_raw (u, iostat=iostat)
    read (u, iostat=iostat) associated_child_first
    read (u, iostat=iostat) associated_next
    if (present (parent))  node%parent => parent
    if (associated_child_first) then
       allocate (child)
       node%child_first => child
       node%child_last => null ()
       call node_read_raw_rec (child, u, node, iostat=iostat)
       do while (associated (child))
          child%previous => node%child_last
          node%child_last => child
          child => child%next
       end do
    else
       read (u, iostat=iostat)  node%me_index
       read (u, iostat=iostat)  node%me
    end if
    if (associated_next) then
       allocate (node%next)
       call node_read_raw_rec (node%next, u, parent, iostat=iostat)
    end if
  end subroutine node_read_raw_rec

  module subroutine state_matrix_init (state, store_values, n_counters)
    class(state_matrix_t), intent(out) :: state
    logical, intent(in), optional :: store_values
    integer, intent(in), optional :: n_counters
    allocate (state%root)
    if (present (store_values)) &
       state%leaf_nodes_store_values = store_values
    if (present (n_counters)) state%n_counters = n_counters
  end subroutine state_matrix_init

  module subroutine state_matrix_final (state)
    class(state_matrix_t), intent(inout) :: state
    if (allocated (state%me))  deallocate (state%me)
    if (associated (state%root))  call node_delete (state%root)
    state%depth = 0
    state%n_matrix_elements = 0
  end subroutine state_matrix_final

  module subroutine state_matrix_write (state, unit, write_value_list, &
         verbose, col_verbose, testflag)
    class(state_matrix_t), intent(in) :: state
    logical, intent(in), optional :: write_value_list, verbose, col_verbose
    logical, intent(in), optional :: testflag
    integer, intent(in), optional :: unit
    complex(default) :: me_dum
    character(len=7) :: fmt
    integer :: u
    integer :: i
    call pac_fmt (fmt, FMT_19, FMT_17, testflag)
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A," // fmt // ")") "State matrix:  norm = ", state%norm
    if (associated (state%root)) then
       if (allocated (state%me)) then
          call node_write_rec (state%root, state%me, verbose = verbose, &
               indent = 1, unit = u, col_verbose = col_verbose, &
               testflag = testflag)
       else
          call node_write_rec (state%root, verbose = verbose, indent = 1, &
               unit = u, col_verbose = col_verbose, testflag = testflag)
       end if
    end if
    if (present (write_value_list)) then
       if (write_value_list .and. allocated (state%me)) then
          do i = 1, size (state%me)
             write (u, "(1x,I0,A)", advance="no")  i, ":"
             me_dum = state%me(i)
             if (real(state%me(i)) == -real(state%me(i))) then
                me_dum = &
                     cmplx (0._default, aimag(me_dum), kind=default)
             end if
             if (aimag(me_dum) == -aimag(me_dum)) then
                me_dum = &
                     cmplx (real(me_dum), 0._default, kind=default)
             end if
             write (u, "('('," // fmt // ",','," // fmt // &
                  ",')')")  me_dum
          end do
       end if
    end if
  end subroutine state_matrix_write

  module subroutine state_matrix_write_raw (state, u)
    class(state_matrix_t), intent(in), target :: state
    integer, intent(in) :: u
    logical :: is_defined
    integer :: depth, j
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    is_defined = state%is_defined ()
    write (u)  is_defined
    if (is_defined) then
       write (u)  state%get_norm ()
       write (u)  state%get_n_leaves ()
       depth = state%get_depth ()
       write (u)  depth
       allocate (qn (depth))
       call it%init (state)
       do while (it%is_valid ())
          qn = it%get_quantum_numbers ()
          do j = 1, depth
             call qn(j)%write_raw (u)
          end do
          write (u)  it%get_me_index ()
          write (u)  it%get_matrix_element ()
          call it%advance ()
       end do
    end if
  end subroutine state_matrix_write_raw

  module subroutine state_matrix_read_raw (state, u, iostat)
    class(state_matrix_t), intent(out) :: state
    integer, intent(in) :: u
    integer, intent(out) :: iostat
    logical :: is_defined
    real(default) :: norm
    integer :: n_leaves, depth, i, j
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer :: me_index
    complex(default) :: me
    read (u, iostat=iostat)  is_defined
    if (iostat /= 0)  goto 1
    if (is_defined) then
       call state%init (store_values = .true.)
       read (u, iostat=iostat)  norm
       if (iostat /= 0)  goto 1
       call state_matrix_set_norm (state, norm)
       read (u)  n_leaves
       if (iostat /= 0)  goto 1
       read (u)  depth
       if (iostat /= 0)  goto 1
       allocate (qn (depth))
       do i = 1, n_leaves
          do j = 1, depth
             call qn(j)%read_raw (u, iostat=iostat)
             if (iostat /= 0)  goto 1
          end do
          read (u, iostat=iostat)  me_index
          if (iostat /= 0)  goto 1
          read (u, iostat=iostat)  me
          if (iostat /= 0)  goto 1
          call state%add_state (qn, index = me_index, value = me)
       end do
       call state_matrix_freeze (state)
    end if
    return

    ! Clean up on error
1   continue
    call state%final ()
  end subroutine state_matrix_read_raw

  module subroutine state_matrix_set_model (state, model)
    class(state_matrix_t), intent(inout), target :: state
    class(model_data_t), intent(in), target :: model
    type(state_iterator_t) :: it
    call it%init (state)
    do while (it%is_valid ())
       call it%set_model (model)
       call it%advance ()
    end do
  end subroutine state_matrix_set_model

  module subroutine state_matrix_tag_hard_process (state, tagged_state, tag)
    class(state_matrix_t), intent(in), target :: state
    type(state_matrix_t), intent(out) :: tagged_state
    integer, dimension(:), intent(in), optional :: tag
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    complex(default) :: value
    integer :: i
    call tagged_state%init (store_values = .true.)
    call it%init (state)
    do while (it%is_valid ())
       qn = it%get_quantum_numbers ()
       value = it%get_matrix_element ()
       if (present (tag)) then
          do i = 1, size (tag)
             call qn(tag(i))%tag_hard_process ()
          end do
       else
          call qn%tag_hard_process ()
       end if
       call tagged_state%add_state (qn, index = it%get_me_index (), value = value)
       call it%advance ()
    end do
    call tagged_state%freeze ()
  end subroutine state_matrix_tag_hard_process

  elemental module function state_matrix_is_defined (state) result (defined)
    logical :: defined
    class(state_matrix_t), intent(in) :: state
    defined = associated (state%root)
  end function state_matrix_is_defined

  elemental module function state_matrix_is_empty (state) result (flag)
    logical :: flag
    class(state_matrix_t), intent(in) :: state
    flag = state%depth == 0
  end function state_matrix_is_empty

  pure module function state_matrix_get_n_matrix_elements_all (state) result (n)
    integer :: n
    class(state_matrix_t), intent(in) :: state
    n = state%n_matrix_elements
  end function state_matrix_get_n_matrix_elements_all

  module function state_matrix_get_n_matrix_elements_mask (state, qn_mask) result (n)
    integer :: n
    class(state_matrix_t), intent(in) :: state
    type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(size(qn_mask)) :: qn
    type(state_matrix_t) :: state_tmp
    call state_tmp%init ()
    call it%init (state)
    do while (it%is_valid ())
       qn = it%get_quantum_numbers ()
       call qn%undefine (qn_mask)
       call state_tmp%add_state (qn)
       call it%advance ()
    end do
    n = state_tmp%n_matrix_elements
    call state_tmp%final ()
  end function state_matrix_get_n_matrix_elements_mask

  pure module function state_matrix_get_me_size (state) result (n)
    integer :: n
    class(state_matrix_t), intent(in) :: state
    if (allocated (state%me)) then
       n = size (state%me)
    else
       n = 0
    end if
  end function state_matrix_get_me_size

  module function state_matrix_compute_n_sub (state) result (n_sub)
    integer :: n_sub
    class(state_matrix_t), intent(in) :: state
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(state%depth) :: qn
    integer :: sub, sub_pos
    n_sub = 0
    call it%init (state)
    do while (it%is_valid ())
       qn = it%get_quantum_numbers ()
       sub = 0
       sub_pos = qn_array_sub_pos ()
       if (sub_pos > 0)  sub = qn(sub_pos)%get_sub ()
       if (sub > n_sub)  n_sub = sub
       call it%advance ()
    end do
  contains
    function qn_array_sub_pos () result (pos)
      integer :: pos
      integer :: i
      pos = 0
      do i = 1, state%depth
         if (qn(i)%get_sub () > 0) then
            pos = i
            exit
         end if
      end do
    end function qn_array_sub_pos
  end function state_matrix_compute_n_sub

  module subroutine state_matrix_set_n_sub (state)
    class(state_matrix_t), intent(inout) :: state
    state%n_sub = state%compute_n_sub ()
  end subroutine state_matrix_set_n_sub

  module function state_matrix_get_n_sub (state) result (n_sub)
    integer :: n_sub
    class(state_matrix_t), intent(in) :: state
    if (state%n_sub < 0) then
       call msg_bug ("[state_matrix_get_n_sub] number of subtractions not set.")
    end if
    n_sub = state%n_sub
  end function state_matrix_get_n_sub

  module function state_matrix_get_n_leaves (state) result (n)
    integer :: n
    class(state_matrix_t), intent(in) :: state
    type(state_iterator_t) :: it
    n = 0
    call it%init (state)
    do while (it%is_valid ())
       n = n + 1
       call it%advance ()
    end do
  end function state_matrix_get_n_leaves

  pure module function state_matrix_get_depth (state) result (depth)
    integer :: depth
    class(state_matrix_t), intent(in) :: state
    depth = state%depth
  end function state_matrix_get_depth

  pure module function state_matrix_get_norm (state) result (norm)
    real(default) :: norm
    class(state_matrix_t), intent(in) :: state
    norm = state%norm
  end function state_matrix_get_norm

  module function state_matrix_get_quantum_number (state, i, by_me_index) result (qn)
    class(state_matrix_t), intent(in), target :: state
    integer, intent(in) :: i
    logical, intent(in), optional :: by_me_index
    logical :: opt_by_me_index
    type(quantum_numbers_t), dimension(state%depth) :: qn
    type(state_iterator_t) :: it
    integer :: k
    opt_by_me_index = .false.
    if (present (by_me_index)) opt_by_me_index = by_me_index
    k = 0
    call it%init (state)
    do while (it%is_valid ())
       if (opt_by_me_index) then
          k = it%get_me_index ()
       else
          k = k + 1
       end if
       if (k == i) then
          qn = it%get_quantum_numbers ()
          exit
       end if
       call it%advance ()
    end do
  end function state_matrix_get_quantum_number

  module subroutine state_matrix_get_quantum_numbers_all (state, qn)
    class(state_matrix_t), intent(in), target :: state
    type(quantum_numbers_t), intent(out), dimension(:,:), allocatable :: qn
    integer :: i
    allocate (qn (state%get_n_matrix_elements (), &
       state%get_depth()))
    do i = 1, state%get_n_matrix_elements ()
       qn (i, :) = state%get_quantum_number (i)
    end do
  end subroutine state_matrix_get_quantum_numbers_all

  module subroutine state_matrix_get_quantum_numbers_mask (state, qn_mask, qn)
    class(state_matrix_t), intent(in), target :: state
    type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask
    type(quantum_numbers_t), intent(out), dimension(:,:), allocatable :: qn
    type(quantum_numbers_t), dimension(:), allocatable :: qn_tmp
    type(state_matrix_t) :: state_tmp
    type(state_iterator_t) :: it
    integer :: i, n
    n = state%get_n_matrix_elements (qn_mask)
    allocate (qn (n, state%get_depth ()))
    allocate (qn_tmp (state%get_depth ()))
    call it%init (state)
    call state_tmp%init ()
    do while (it%is_valid ())
       qn_tmp = it%get_quantum_numbers ()
       call qn_tmp%undefine (qn_mask)
       call state_tmp%add_state (qn_tmp)
       call it%advance ()
    end do
    do i = 1, n
       qn (i, :) = state_tmp%get_quantum_number (i)
    end do
    call state_tmp%final ()
  end subroutine state_matrix_get_quantum_numbers_mask

  module subroutine state_matrix_get_flavors (state, only_elementary, qn_mask, flv)
    class(state_matrix_t), intent(in), target :: state
    logical, intent(in) :: only_elementary
    type(quantum_numbers_mask_t), intent(in), dimension(:), optional :: qn_mask
    integer, intent(out), dimension(:,:), allocatable :: flv
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn
    integer :: i_flv, n_partons
    type(flavor_t), dimension(:), allocatable :: flv_flv
    if (present (qn_mask)) then
       call state%get_quantum_numbers (qn_mask, qn)
    else
       call state%get_quantum_numbers (qn)
    end if
    allocate (flv_flv (size (qn, dim=2)))
    if (only_elementary) then
       flv_flv = qn(1, :)%get_flavor ()
       n_partons = count (is_elementary (flv_flv%get_pdg ()))
    end if
    allocate (flv (n_partons, size (qn, dim=1)))
    associate (n_flv => size (qn, dim=1))
      do i_flv = 1, size (qn, dim=1)
         flv_flv = qn(i_flv, :)%get_flavor ()
         flv(:, i_flv) = pack (flv_flv%get_pdg (), is_elementary(flv_flv%get_pdg()))
      end do
    end associate
  contains
    elemental function is_elementary (pdg)
      logical :: is_elementary
      integer, intent(in) :: pdg
      is_elementary = abs(pdg) /= 2212 .and. abs(pdg) /= 92 .and. abs(pdg) /= 93
    end function is_elementary
  end subroutine state_matrix_get_flavors

  elemental module function state_matrix_get_matrix_element_single (state, i) result (me)
    complex(default) :: me
    class(state_matrix_t), intent(in) :: state
    integer, intent(in) :: i
    if (allocated (state%me)) then
       me = state%me(i)
    else
       me = 0
    end if
  end function state_matrix_get_matrix_element_single

  module function state_matrix_get_matrix_element_array (state) result (me)
     complex(default), dimension(:), allocatable :: me
     class(state_matrix_t), intent(in) :: state
     if (allocated (state%me)) then
        allocate (me (size (state%me)))
        me = state%me
     else
        me = 0
     end if
  end function state_matrix_get_matrix_element_array

  module function state_matrix_get_max_color_value (state) result (cmax)
    integer :: cmax
    class(state_matrix_t), intent(in) :: state
    if (associated (state%root)) then
       cmax = node_get_max_color_value (state%root)
    else
       cmax = 0
    end if
  contains
    recursive function node_get_max_color_value (node) result (cmax)
      integer :: cmax
      type(node_t), intent(in), target :: node
      type(node_t), pointer :: current
      cmax = quantum_numbers_get_max_color_value (node%qn)
      current => node%child_first
      do while (associated (current))
         cmax = max (cmax, node_get_max_color_value (current))
         current => current%next
      end do
    end function node_get_max_color_value
  end function state_matrix_get_max_color_value

  module subroutine state_matrix_add_state (state, qn, index, value, &
         sum_values, counter_index, ignore_sub_for_qn, me_index)
    class(state_matrix_t), intent(inout) :: state
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    integer, intent(in), optional :: index
    complex(default), intent(in), optional :: value
    logical, intent(in), optional :: sum_values
    integer, intent(in), optional :: counter_index
    logical, intent(in), optional :: ignore_sub_for_qn
    integer, intent(out), optional :: me_index
    logical :: set_index, get_index, add
    set_index = present (index)
    get_index = present (me_index)
    add = .false.;  if (present (sum_values))  add = sum_values
    if (state%depth == 0) then
       state%depth = size (qn)
    else if (state%depth /= size (qn)) then
       call state%write ()
       call msg_bug ("State matrix: depth mismatch")
    end if
    if (size (qn) > 0)  call node_make_branch (state%root, qn)
  contains
     recursive subroutine node_make_branch (parent, qn)
       type(node_t), pointer :: parent
       type(quantum_numbers_t), dimension(:), intent(in) :: qn
       type(node_t), pointer :: child
       logical :: match
       match = .false.
       child => parent%child_first
       SCAN_CHILDREN: do while (associated (child))
          if (present (ignore_sub_for_qn)) then
             if (ignore_sub_for_qn) then
                match = quantum_numbers_eq_wo_sub (child%qn, qn(1))
             else
                match = child%qn == qn(1)
             end if
          else
             match = child%qn == qn(1)
          end if
          if (match)  exit SCAN_CHILDREN
          child => child%next
       end do SCAN_CHILDREN
       if (.not. match) then
          call node_append_child (parent, child)
          child%qn = qn(1)
       end if
       select case (size (qn))
       case (1)
          if (.not. match) then
             state%n_matrix_elements = state%n_matrix_elements + 1
             child%me_index = state%n_matrix_elements
          end if
          if (set_index) then
             child%me_index = index
          end if
          if (get_index) then
             me_index = child%me_index
          end if
          if (present (counter_index)) then
             if (.not. allocated (child%me_count)) then
                allocate (child%me_count (state%n_counters))
                child%me_count = 0
             end if
             child%me_count(counter_index) = child%me_count(counter_index) + 1
          end if
          if (present (value)) then
             if (add) then
                child%me = child%me + value
             else
                child%me = value
             end if
          end if
       case (2:)
          call node_make_branch (child, qn(2:))
       end select
     end subroutine node_make_branch
   end subroutine state_matrix_add_state

  module subroutine state_matrix_collapse (state, mask)
    class(state_matrix_t), intent(inout) :: state
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
    type(state_matrix_t) :: red_state
    if (state%is_defined ()) then
       call state%reduce (mask, red_state)
       call state%final ()
       state = red_state
    end if
   end subroutine state_matrix_collapse

  module subroutine state_matrix_reduce (state, mask, red_state, keep_me_index)
    class(state_matrix_t), intent(in), target :: state
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
    type(state_matrix_t), intent(out) :: red_state
    logical, optional, intent(in)  :: keep_me_index
    logical :: opt_keep_me_index
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(size(mask)) :: qn
    opt_keep_me_index = .false.
    if (present (keep_me_index)) opt_keep_me_index = keep_me_index
    call red_state%init ()
    call it%init (state)
    do while (it%is_valid ())
       qn = it%get_quantum_numbers ()
       call qn%undefine (mask)
       if (opt_keep_me_index) then
          call red_state%add_state (qn, index = it%get_me_index ())
       else
          call red_state%add_state (qn)
       end if
       call it%advance ()
    end do
    if (.not. opt_keep_me_index) then
       call red_state%freeze ()
    end if
  end subroutine state_matrix_reduce

  module subroutine state_matrix_reorder_me (state, ordered_state)
    class(state_matrix_t), intent(in), target :: state
    type(state_matrix_t), intent(out) :: ordered_state
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(state%depth) :: qn
    integer,  dimension(:), allocatable :: me_index
    integer :: i
    call ordered_state%init ()
    call get_me_index_sorted (state, me_index)
    i = 1; call it%init (state)
    do while (it%is_valid ())
       qn = it%get_quantum_numbers ()
       call ordered_state%add_state (qn, index = me_index(i))
       i = i + 1; call it%advance ()
    end do
    call ordered_state%freeze ()
  contains
    subroutine get_me_index_sorted (state, me_index)
      class(state_matrix_t), intent(in), target :: state
      integer, dimension(:), allocatable, intent(out) :: me_index
      type(state_iterator_t) :: it
      integer :: i, j
      integer, dimension(:), allocatable :: me_index_unsorted, me_index_sorted
      associate (n_matrix_elements => state%get_n_matrix_elements ())
        allocate (me_index(n_matrix_elements), source = 0)
        allocate (me_index_sorted(n_matrix_elements), source = 0)
        allocate (me_index_unsorted(n_matrix_elements), source = 0)
        i = 1; call it%init (state)
        do while (it%is_valid ())
           me_index_unsorted(i) = it%get_me_index ()
           i = i + 1
           call it%advance ()
        end do
        me_index_sorted = sort (me_index_unsorted)
        ! We do not care about efficiency at this point.
        UNSORTED: do i = 1, n_matrix_elements
           SORTED: do j = 1, n_matrix_elements
              if (me_index_unsorted(i) == me_index_sorted(j)) then
                 me_index(i) = j
                 cycle UNSORTED
              end if
           end do SORTED
        end do UNSORTED
      end associate
    end subroutine get_me_index_sorted
  end subroutine state_matrix_reorder_me

  module subroutine state_matrix_set_duplicate_flv_zero (state)
    class(state_matrix_t), intent(inout), target :: state
    type(quantum_numbers_t), dimension(state%depth) :: qn
    type(flavor_t) :: flv
    type(state_flv_content_t), allocatable :: state_flv
    logical, dimension(:), allocatable :: hard_mask, sort_mask, duplicate_mask
    integer :: i, j, n_in, n_flvs
    n_flvs = state%get_depth ()
    n_in = 2
    !!! TODO (PS-28-07-21) n_in should not be hard coded to work for decays
    !!! This assumes that the positions of the non-hard flavors are the same for all flavor structures.
    qn = state%get_quantum_number(1)
    allocate (hard_mask(n_flvs))
    do i = 1, n_flvs
       flv = qn(i)%get_flavor()
       hard_mask(i) = flv%is_hard_process ()
    end do
    allocate (sort_mask(n_flvs))
    sort_mask = hard_mask
    j = 0
    do i = 1, n_flvs
       if (j == n_in) exit
       if (sort_mask(i)) then
          sort_mask(i) = .false.
          j = j + 1
       end if
    end do
    allocate (state_flv)
    call state_flv%fill (state, sort_mask)
    call state_flv%find_duplicates (duplicate_mask)
    do i = 1, state%get_n_matrix_elements ()
       if (duplicate_mask(i)) then
          call state%set_matrix_element_single(i, cmplx(zero, zero, default))
       end if
    end do
  end subroutine state_matrix_set_duplicate_flv_zero

  module subroutine state_matrix_freeze (state)
    class(state_matrix_t), intent(inout), target :: state
    type(state_iterator_t) :: it
    if (associated (state%root)) then
       if (allocated (state%me))  deallocate (state%me)
       allocate (state%me (state%n_matrix_elements))
       state%me = 0
       call state%set_n_sub ()
    end if
    if (state%leaf_nodes_store_values) then
       call it%init (state)
       do while (it%is_valid ())
          state%me(it%get_me_index ()) = it%get_matrix_element ()
          call it%advance ()
       end do
       state%leaf_nodes_store_values = .false.
    end if
  end subroutine state_matrix_freeze

  module subroutine state_matrix_set_matrix_element_qn (state, qn, value)
    class(state_matrix_t), intent(inout), target :: state
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    complex(default), intent(in) :: value
    type(state_iterator_t) :: it
    if (.not. allocated (state%me)) then
       allocate (state%me (size(qn)))
    end if
    call it%init (state)
    call it%go_to_qn (qn)
    call it%set_matrix_element (value)
  end subroutine state_matrix_set_matrix_element_qn

  module subroutine state_matrix_set_matrix_element_all (state, value)
    class(state_matrix_t), intent(inout) :: state
    complex(default), intent(in) :: value
    if (.not. allocated (state%me)) then
       allocate (state%me (state%n_matrix_elements))
    end if
    state%me = value
  end subroutine state_matrix_set_matrix_element_all

  module subroutine state_matrix_set_matrix_element_array (state, value, range)
    class(state_matrix_t), intent(inout) :: state
    complex(default), intent(in), dimension(:) :: value
    integer, intent(in), dimension(:), optional :: range
    if (present (range)) then
       state%me(range) = value
    else
       if (.not. allocated (state%me)) &
            allocate (state%me (size (value)))
       state%me(:) = value
    end if
  end subroutine state_matrix_set_matrix_element_array

  pure module subroutine state_matrix_set_matrix_element_single (state, i, value)
    class(state_matrix_t), intent(inout) :: state
    integer, intent(in) :: i
    complex(default), intent(in) :: value
    if (.not. allocated (state%me)) then
       allocate (state%me (state%n_matrix_elements))
    end if
    state%me(i) = value
  end subroutine state_matrix_set_matrix_element_single

  module subroutine state_matrix_set_matrix_element_clone (state, state1)
    class(state_matrix_t), intent(inout) :: state
    type(state_matrix_t), intent(in) :: state1
    if (.not. allocated (state1%me)) return
    if (.not. allocated (state%me)) allocate (state%me (size (state1%me)))
    state%me = state1%me
  end subroutine state_matrix_set_matrix_element_clone

  module subroutine state_matrix_add_to_matrix_element (state, qn, value, match_only_flavor)
    class(state_matrix_t), intent(inout), target :: state
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    complex(default), intent(in) :: value
    logical, intent(in), optional :: match_only_flavor
    type(state_iterator_t) :: it
    call it%init (state)
    call it%go_to_qn (qn, match_only_flavor)
    if (it%is_valid ()) then
       call it%add_to_matrix_element (value)
    else
       call msg_fatal ("Cannot add to matrix element - it%node not allocated")
    end if
  end subroutine state_matrix_add_to_matrix_element

  module subroutine state_iterator_init (it, state)
    class(state_iterator_t), intent(out) :: it
    type(state_matrix_t), intent(in), target :: state
    it%state => state
    it%depth = state%depth
    if (state%is_defined ()) then
       it%node => state%root
       do while (associated (it%node%child_first))
          it%node => it%node%child_first
       end do
    else
       it%node => null ()
    end if
  end subroutine state_iterator_init

  module subroutine state_iterator_advance (it)
    class(state_iterator_t), intent(inout) :: it
    call find_next (it%node)
  contains
    recursive subroutine find_next (node_in)
      type(node_t), intent(in), target :: node_in
      type(node_t), pointer :: node
      node => node_in
      if (associated (node%next)) then
         node => node%next
         do while (associated (node%child_first))
            node => node%child_first
         end do
         it%node => node
      else if (associated (node%parent)) then
         call find_next (node%parent)
      else
         it%node => null ()
      end if
    end subroutine find_next
  end subroutine state_iterator_advance

  module function state_iterator_is_valid (it) result (defined)
    logical :: defined
    class(state_iterator_t), intent(in) :: it
    defined = associated (it%node)
  end function state_iterator_is_valid

  module function state_iterator_get_me_index (it) result (n)
    integer :: n
    class(state_iterator_t), intent(in) :: it
    n = it%node%me_index
  end function state_iterator_get_me_index

  module function state_iterator_get_me_count (it) result (n)
    integer, dimension(:), allocatable :: n
    class(state_iterator_t), intent(in) :: it
    if (allocated (it%node%me_count)) then
       allocate (n (size (it%node%me_count)))
       n = it%node%me_count
    else
       allocate (n (0))
    end if
  end function state_iterator_get_me_count

  pure module function state_iterator_get_depth (state_iterator) result (depth)
    integer :: depth
    class(state_iterator_t), intent(in) :: state_iterator
    depth = state_iterator%depth
  end function state_iterator_get_depth

  module subroutine state_iterator_go_to_qn (it, qn, match_only_flavor)
    class(state_iterator_t), intent(inout) :: it
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    logical, intent(in), optional :: match_only_flavor
    type(quantum_numbers_t), dimension(:), allocatable :: qn_hard, qn_tmp
    logical :: match_flv
    match_flv = .false.; if (present (match_only_flavor)) match_flv = .true.
    do while (it%is_valid ())
       if (match_flv) then
          qn_tmp = it%get_quantum_numbers ()
          qn_hard = pack (qn_tmp, qn_tmp%are_hard_process ())
          if (all (qn .fmatch. qn_hard)) then
             return
          else
             call it%advance ()
          end if
       else
          if (all (qn == it%get_quantum_numbers ())) then
             return
          else
             call it%advance ()
          end if
       end if
    end do
  end subroutine state_iterator_go_to_qn

  module function state_iterator_get_qn_multi (it) result (qn)
    class(state_iterator_t), intent(in) :: it
    type(quantum_numbers_t), dimension(it%depth) :: qn
    type(node_t), pointer :: node
    integer :: i
    node => it%node
    do i = it%depth, 1, -1
       qn(i) = node%qn
       node => node%parent
    end do
  end function state_iterator_get_qn_multi

  module function state_iterator_get_flv_multi (it) result (flv)
    class(state_iterator_t), intent(in) :: it
    type(flavor_t), dimension(it%depth) :: flv
    flv = quantum_numbers_get_flavor &
         (it%get_quantum_numbers ())
  end function state_iterator_get_flv_multi

  module function state_iterator_get_col_multi (it) result (col)
    class(state_iterator_t), intent(in) :: it
    type(color_t), dimension(it%depth) :: col
    col = quantum_numbers_get_color &
         (it%get_quantum_numbers ())
  end function state_iterator_get_col_multi

  module function state_iterator_get_hel_multi (it) result (hel)
    class(state_iterator_t), intent(in) :: it
    type(helicity_t), dimension(it%depth) :: hel
    hel = quantum_numbers_get_helicity &
         (it%get_quantum_numbers ())
  end function state_iterator_get_hel_multi

  module function state_iterator_get_qn_slice (it, index) result (qn)
    class(state_iterator_t), intent(in) :: it
    integer, dimension(:), intent(in) :: index
    type(quantum_numbers_t), dimension(size(index)) :: qn
    type(quantum_numbers_t), dimension(it%depth) :: qn_tmp
    qn_tmp = state_iterator_get_qn_multi (it)
    qn = qn_tmp(index)
  end function state_iterator_get_qn_slice

  module function state_iterator_get_flv_slice (it, index) result (flv)
    class(state_iterator_t), intent(in) :: it
    integer, dimension(:), intent(in) :: index
    type(flavor_t), dimension(size(index)) :: flv
    flv = quantum_numbers_get_flavor &
         (it%get_quantum_numbers (index))
  end function state_iterator_get_flv_slice

  module function state_iterator_get_col_slice (it, index) result (col)
    class(state_iterator_t), intent(in) :: it
    integer, dimension(:), intent(in) :: index
    type(color_t), dimension(size(index)) :: col
    col = quantum_numbers_get_color &
         (it%get_quantum_numbers (index))
  end function state_iterator_get_col_slice

  module function state_iterator_get_hel_slice (it, index) result (hel)
    class(state_iterator_t), intent(in) :: it
    integer, dimension(:), intent(in) :: index
    type(helicity_t), dimension(size(index)) :: hel
    hel = quantum_numbers_get_helicity &
         (it%get_quantum_numbers (index))
  end function state_iterator_get_hel_slice

  module function state_iterator_get_qn_range (it, k1, k2) result (qn)
    class(state_iterator_t), intent(in) :: it
    integer, intent(in) :: k1, k2
    type(quantum_numbers_t), dimension(k2-k1+1) :: qn
    type(node_t), pointer :: node
    integer :: i
    node => it%node
    SCAN: do i = it%depth, 1, -1
       if (k1 <= i .and. i <= k2) then
          qn(i-k1+1) = node%qn
       else
          node => node%parent
       end if
    end do SCAN
  end function state_iterator_get_qn_range

  module function state_iterator_get_flv_range (it, k1, k2) result (flv)
    class(state_iterator_t), intent(in) :: it
    integer, intent(in) :: k1, k2
    type(flavor_t), dimension(k2-k1+1) :: flv
    flv = quantum_numbers_get_flavor &
         (it%get_quantum_numbers (k1, k2))
  end function state_iterator_get_flv_range

  module function state_iterator_get_col_range (it, k1, k2) result (col)
    class(state_iterator_t), intent(in) :: it
    integer, intent(in) :: k1, k2
    type(color_t), dimension(k2-k1+1) :: col
    col = quantum_numbers_get_color &
         (it%get_quantum_numbers (k1, k2))
  end function state_iterator_get_col_range

  module function state_iterator_get_hel_range (it, k1, k2) result (hel)
    class(state_iterator_t), intent(in) :: it
    integer, intent(in) :: k1, k2
    type(helicity_t), dimension(k2-k1+1) :: hel
    hel = quantum_numbers_get_helicity &
         (it%get_quantum_numbers (k1, k2))
  end function state_iterator_get_hel_range

  module function state_iterator_get_qn_single (it, k) result (qn)
    class(state_iterator_t), intent(in) :: it
    integer, intent(in) :: k
    type(quantum_numbers_t) :: qn
    type(node_t), pointer :: node
    integer :: i
    node => it%node
    SCAN: do i = it%depth, 1, -1
       if (i == k) then
          qn = node%qn
          exit SCAN
       else
          node => node%parent
       end if
    end do SCAN
  end function state_iterator_get_qn_single

  module function state_iterator_get_flv_single (it, k) result (flv)
    class(state_iterator_t), intent(in) :: it
    integer, intent(in) :: k
    type(flavor_t) :: flv
    flv = quantum_numbers_get_flavor &
         (it%get_quantum_numbers (k))
  end function state_iterator_get_flv_single

  module function state_iterator_get_col_single (it, k) result (col)
    class(state_iterator_t), intent(in) :: it
    integer, intent(in) :: k
    type(color_t) :: col
    col = quantum_numbers_get_color &
         (it%get_quantum_numbers (k))
  end function state_iterator_get_col_single

  module function state_iterator_get_hel_single (it, k) result (hel)
    class(state_iterator_t), intent(in) :: it
    integer, intent(in) :: k
    type(helicity_t) :: hel
    hel = quantum_numbers_get_helicity &
         (it%get_quantum_numbers (k))
  end function state_iterator_get_hel_single

  module subroutine state_iterator_set_model (it, model)
    class(state_iterator_t), intent(inout) :: it
    class(model_data_t), intent(in), target :: model
    type(node_t), pointer :: node
    integer :: i
    node => it%node
    do i = it%depth, 1, -1
       call node%qn%set_model (model)
       node => node%parent
    end do
  end subroutine state_iterator_set_model

  module subroutine state_iterator_retag_hard_process (it, i, hard)
    class(state_iterator_t), intent(inout) :: it
    integer, intent(in) :: i
    logical, intent(in) :: hard
    type(node_t), pointer :: node
    integer :: j
    node => it%node
    do j = 1, it%depth-i
       node => node%parent
    end do
    call node%qn%tag_hard_process (hard)
  end subroutine state_iterator_retag_hard_process

  module function state_iterator_get_matrix_element (it) result (me)
    complex(default) :: me
    class(state_iterator_t), intent(in) :: it
    if (it%state%leaf_nodes_store_values) then
       me = it%node%me
    else if (it%node%me_index /= 0) then
       me = it%state%me(it%node%me_index)
    else
       me = 0
    end if
  end function state_iterator_get_matrix_element

  module subroutine state_iterator_set_matrix_element (it, value)
    class(state_iterator_t), intent(inout) :: it
    complex(default), intent(in) :: value
    if (it%node%me_index /= 0) it%state%me(it%node%me_index) = value
  end subroutine state_iterator_set_matrix_element

  module subroutine state_iterator_add_to_matrix_element (it, value)
    class(state_iterator_t), intent(inout) :: it
    complex(default), intent(in) :: value
    if (it%node%me_index /= 0) &
         it%state%me(it%node%me_index) = it%state%me(it%node%me_index) + value
  end subroutine state_iterator_add_to_matrix_element

  module subroutine state_matrix_assign (state_out, state_in)
    type(state_matrix_t), intent(out) :: state_out
    type(state_matrix_t), intent(in), target :: state_in
    type(state_iterator_t) :: it
    if (.not. state_in%is_defined ())  return
    call state_out%init ()
    call it%init (state_in)
    do while (it%is_valid ())
       call state_out%add_state (it%get_quantum_numbers (), &
            it%get_me_index ())
       call it%advance ()
    end do
    if (allocated (state_in%me)) then
       allocate (state_out%me (size (state_in%me)))
       state_out%me = state_in%me
    end if
    state_out%n_sub = state_in%n_sub
  end subroutine state_matrix_assign

  module subroutine state_matrix_get_diagonal_entries (state, i)
    class(state_matrix_t), intent(in) :: state
    integer, dimension(:), allocatable, intent(out) :: i
    integer, dimension(state%n_matrix_elements) :: tmp
    integer :: n
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    n = 0
    call it%init (state)
    allocate (qn (it%depth))
    do while (it%is_valid ())
       qn = it%get_quantum_numbers ()
       if (all (qn%are_diagonal ())) then
          n = n + 1
          tmp(n) = it%get_me_index ()
       end if
       call it%advance ()
    end do
    allocate (i(n))
    if (n > 0) i = tmp(:n)
  end subroutine state_matrix_get_diagonal_entries

  module subroutine state_matrix_renormalize (state, factor)
    class(state_matrix_t), intent(inout) :: state
    complex(default), intent(in) :: factor
    state%me = state%me * factor
  end subroutine state_matrix_renormalize

  module subroutine state_matrix_normalize_by_trace (state)
    class(state_matrix_t), intent(inout) :: state
    real(default) :: trace
    trace = state%trace ()
    if (trace /= 0) then
       state%me = state%me / trace
       state%norm = state%norm * trace
    end if
  end subroutine state_matrix_normalize_by_trace

  module subroutine state_matrix_normalize_by_max (state)
    class(state_matrix_t), intent(inout) :: state
    real(default) :: m
    m = maxval (abs (state%me))
    if (m /= 0) then
       state%me = state%me / m
       state%norm = state%norm * m
    end if
  end subroutine state_matrix_normalize_by_max

  module subroutine state_matrix_set_norm (state, norm)
    class(state_matrix_t), intent(inout) :: state
    real(default), intent(in) :: norm
    state%norm = norm
  end subroutine state_matrix_set_norm

  pure module function state_matrix_sum (state) result (value)
    complex(default) :: value
    class(state_matrix_t), intent(in) :: state
    value = sum (state%me)
  end function state_matrix_sum

  module function state_matrix_trace (state, qn_in) result (trace)
    complex(default) :: trace
    class(state_matrix_t), intent(in), target :: state
    type(quantum_numbers_t), dimension(:), intent(in), optional :: qn_in
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    type(state_iterator_t) :: it
    allocate (qn (state%get_depth ()))
    trace = 0
    call it%init (state)
    do while (it%is_valid ())
       qn = it%get_quantum_numbers ()
       if (present (qn_in)) then
          if (.not. all (qn .fhmatch. qn_in)) then
             call it%advance ();  cycle
          end if
       end if
       if (all (qn%are_diagonal ())) then
          trace = trace + it%get_matrix_element ()
       end if
       call it%advance ()
    end do
  end function state_matrix_trace

  module subroutine state_matrix_add_color_contractions (state)
    class(state_matrix_t), intent(inout), target :: state
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn_con
    integer, dimension(:), allocatable :: me_index
    integer :: depth, n_me, i, j
    depth = state%get_depth ()
    n_me = state%get_n_matrix_elements ()
    allocate (qn (depth, n_me))
    allocate (me_index (n_me))
    i = 0
    call it%init (state)
    do while (it%is_valid ())
       i = i + 1
       qn(:,i) = it%get_quantum_numbers ()
       me_index(i) = it%get_me_index ()
       call it%advance ()
    end do
    do i = 1, n_me
       call quantum_number_array_make_color_contractions (qn(:,i), qn_con)
       do j = 1, size (qn_con, 2)
          call state%add_state (qn_con(:,j), index = me_index(i))
       end do
    end do
  end subroutine state_matrix_add_color_contractions

  module subroutine merge_state_matrices (state1, state2, state3)
    type(state_matrix_t), intent(in), target :: state1, state2
    type(state_matrix_t), intent(out) :: state3
    type(state_iterator_t) :: it1, it2
    type(quantum_numbers_t), dimension(state1%depth) :: qn1, qn2
    if (state1%depth /= state2%depth) then
       call state1%write ()
       call state2%write ()
       call msg_bug ("State matrices merge impossible: incompatible depths")
    end if
    call state3%init ()
    call it1%init (state1)
    do while (it1%is_valid ())
       qn1 = it1%get_quantum_numbers ()
       call it2%init (state2)
       do while (it2%is_valid ())
          qn2 = it2%get_quantum_numbers ()
          call state3%add_state (qn1 .merge. qn2)
          call it2%advance ()
       end do
       call it1%advance ()
    end do
    call state3%freeze ()
  end subroutine merge_state_matrices

  pure module subroutine state_matrix_evaluate_product &
       (state, i, state1, state2, index1, index2)
    class(state_matrix_t), intent(inout) :: state
    integer, intent(in) :: i
    type(state_matrix_t), intent(in) :: state1, state2
    integer, dimension(:), intent(in) :: index1, index2
    state%me(i) = &
         dot_product (conjg (state1%me(index1)), state2%me(index2))
    state%norm = state1%norm * state2%norm
  end subroutine state_matrix_evaluate_product

  pure module subroutine state_matrix_evaluate_product_cf &
       (state, i, state1, state2, index1, index2, factor)
    class(state_matrix_t), intent(inout) :: state
    integer, intent(in) :: i
    type(state_matrix_t), intent(in) :: state1, state2
    integer, dimension(:), intent(in) :: index1, index2
    complex(default), dimension(:), intent(in) :: factor
    state%me(i) = &
         dot_product (state1%me(index1), factor * state2%me(index2))
    state%norm = state1%norm * state2%norm
  end subroutine state_matrix_evaluate_product_cf

  pure module subroutine state_matrix_evaluate_square_c (state, i, state1, index1)
    class(state_matrix_t), intent(inout) :: state
    integer, intent(in) :: i
    type(state_matrix_t), intent(in) :: state1
    integer, dimension(:), intent(in) :: index1
    state%me(i) = &
         dot_product (state1%me(index1), state1%me(index1))
    state%norm = abs (state1%norm) ** 2
  end subroutine state_matrix_evaluate_square_c

  pure module subroutine state_matrix_evaluate_sum (state, i, state1, index1)
    class(state_matrix_t), intent(inout) :: state
    integer, intent(in) :: i
    type(state_matrix_t), intent(in) :: state1
    integer, dimension(:), intent(in) :: index1
    state%me(i) = &
         sum (state1%me(index1)) * state1%norm
  end subroutine state_matrix_evaluate_sum

  pure module subroutine state_matrix_evaluate_me_sum (state, i, state1, index1)
    class(state_matrix_t), intent(inout) :: state
    integer, intent(in) :: i
    type(state_matrix_t), intent(in) :: state1
    integer, dimension(:), intent(in) :: index1
    state%me(i) = sum (state1%me(index1))
  end subroutine state_matrix_evaluate_me_sum

  module subroutine outer_multiply_pair (state1, state2, state3)
    type(state_matrix_t), intent(in), target :: state1, state2
    type(state_matrix_t), intent(out) :: state3
    type(state_iterator_t) :: it1, it2
    type(quantum_numbers_t), dimension(state1%depth) :: qn1
    type(quantum_numbers_t), dimension(state2%depth) :: qn2
    type(quantum_numbers_t), dimension(state1%depth+state2%depth) :: qn3
    complex(default) :: val1, val2
    call state3%init (store_values = .true.)
    call it1%init (state1)
    do while (it1%is_valid ())
       qn1 = it1%get_quantum_numbers ()
       val1 = it1%get_matrix_element ()
       call it2%init (state2)
       do while (it2%is_valid ())
          qn2 = it2%get_quantum_numbers ()
          val2 = it2%get_matrix_element ()
          qn3(:state1%depth) = qn1
          qn3(state1%depth+1:) = qn2
          call state3%add_state (qn3, value=val1 * val2)
          call it2%advance ()
       end do
       call it1%advance ()
    end do
    call state3%freeze ()
  end subroutine outer_multiply_pair

  module subroutine outer_multiply_array (state_in, state_out)
    type(state_matrix_t), dimension(:), intent(in), target :: state_in
    type(state_matrix_t), intent(out) :: state_out
    type(state_matrix_t), dimension(:), allocatable, target :: state_tmp
    integer :: i, n
    n = size (state_in)
    select case (n)
    case (0)
       call state_out%init ()
    case (1)
       state_out = state_in(1)
    case (2)
       call outer_multiply_pair (state_in(1), state_in(2), state_out)
    case default
       allocate (state_tmp (n-2))
       call outer_multiply_pair (state_in(1), state_in(2), state_tmp(1))
       do i = 2, n - 2
          call outer_multiply_pair (state_tmp(i-1), state_in(i+1), state_tmp(i))
       end do
       call outer_multiply_pair (state_tmp(n-2), state_in(n), state_out)
       do i = 1, size(state_tmp)
          call state_tmp(i)%final ()
       end do
    end select
  end subroutine outer_multiply_array

  module subroutine state_matrix_factorize &
       (state, mode, x, ok, single_state, correlated_state, qn_in)
    class(state_matrix_t), intent(in), target :: state
    integer, intent(in) :: mode
    real(default), intent(in) :: x
    logical, intent(out) :: ok
    type(state_matrix_t), &
         dimension(:), allocatable, intent(out) :: single_state
    type(state_matrix_t), intent(out), optional :: correlated_state
    type(quantum_numbers_t), dimension(:), intent(in), optional :: qn_in
    type(state_iterator_t) :: it
    real(default) :: s, xt
    complex(default) :: value
    integer :: i, depth
    type(quantum_numbers_t), dimension(:), allocatable :: qn, qn1
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask
    logical, dimension(:), allocatable :: diagonal
    logical, dimension(:,:), allocatable :: mask
    ok = .true.
    if (x /= 0) then
       xt = x * abs (state%trace (qn_in))
    else
       xt = 0
    end if
    s = 0
    depth = state%get_depth ()
    allocate (qn (depth), qn1 (depth), diagonal (depth))
    call it%init (state)
    do while (it%is_valid ())
       qn = it%get_quantum_numbers ()
       if (present (qn_in)) then
          if (.not. all (qn .fhmatch. qn_in)) then
             call it%advance (); cycle
          end if
       end if
       if (all (qn%are_diagonal ())) then
          value = abs (it%get_matrix_element ())
          s = s + value
          if (s > xt)  exit
       end if
       call it%advance ()
    end do
    if (.not. it%is_valid ()) then
       if (s == 0)  ok = .false.
       call it%init (state)
    end if
    allocate (single_state (depth))
    do i = 1, depth
       call single_state(i)%init (store_values = .true.)
    end do
    if (present (correlated_state)) &
         call correlated_state%init (store_values = .true.)
    qn = it%get_quantum_numbers ()
    select case (mode)
    case (FM_SELECT_HELICITY)  ! single branch selected; shortcut
       do i = 1, depth
          call single_state(i)%add_state ([qn(i)], value=value)
       end do
       if (.not. present (correlated_state)) then
          do i = 1, size(single_state)
             call single_state(i)%freeze ()
          end do
          return
       end if
    end select
    allocate (qn_mask (depth))
    call qn_mask%init (.false., .false., .false., .true.)
    call qn%undefine (qn_mask)
    select case (mode)
    case (FM_FACTOR_HELICITY)
       allocate (mask (depth, depth))
       mask = .false.
       forall (i = 1:depth)  mask(i,i) = .true.
    end select
    call it%init (state)
    do while (it%is_valid ())
       qn1 = it%get_quantum_numbers ()
       if (all (qn .match. qn1)) then
          diagonal = qn1%are_diagonal ()
          value = it%get_matrix_element ()
          select case (mode)
          case (FM_IGNORE_HELICITY, FM_CORRELATED_HELICITY)
             !!! trace over diagonal states that match qn
             if (all (diagonal)) then
                do i = 1, depth
                   call single_state(i)%add_state &
                        ([qn(i)], value=value, sum_values=.true.)
                end do
             end if
          case (FM_FACTOR_HELICITY)  !!! trace over all other particles
             do i = 1, depth
                if (all (diagonal .or. mask(:,i))) then
                   call single_state(i)%add_state &
                        ([qn1(i)], value=value, sum_values=.true.)
                end if
             end do
          end select
          if (present (correlated_state)) &
               call correlated_state%add_state (qn1, value=value)
       end if
       call it%advance ()
    end do
    do i = 1, depth
       call single_state(i)%freeze ()
    end do
    if (present (correlated_state)) &
         call correlated_state%freeze ()
  end subroutine state_matrix_factorize

  module function state_matrix_get_polarization_density_matrix &
       (state) result (pol_matrix)
    real(default), dimension(:,:), allocatable :: pol_matrix
    class(state_matrix_t), intent(in) :: state
    type(node_t), pointer :: current => null ()
    !!! What's the generic way to allocate the matrix?
    allocate (pol_matrix (4,4)); pol_matrix = 0
    if (associated (state%root%child_first)) then
       current => state%root%child_first
       do while (associated (current))
         call current%qn%write ()
         current => current%next
       end do
    else
       call msg_fatal ("Polarization state not allocated!")
    end if
  end function state_matrix_get_polarization_density_matrix

  module subroutine state_flv_content_write (state_flv, unit)
    class(state_flv_content_t), intent(in), target :: state_flv
    integer, intent(in), optional :: unit
    integer :: u, n, d, i, j
    u = given_output_unit (unit)
    d = size (state_flv%pdg, 1)
    n = size (state_flv%pdg, 2)
    do i = 1, n
       write (u, "(2x,'PDG =')", advance="no")
       do j = 1, d
          write (u, "(1x,I0)", advance="no")  state_flv%pdg(j,i)
       end do
       write (u, "(' :: map = (')", advance="no")
       do j = 1, d
          write (u, "(1x,I0)", advance="no")  state_flv%map(j,i)
       end do
       write (u, "(' )')")
    end do
  end subroutine state_flv_content_write

  module subroutine state_flv_content_init (state_flv, n, mask)
    class(state_flv_content_t), intent(out) :: state_flv
    integer, intent(in) :: n
    logical, dimension(:), intent(in) :: mask
    integer :: d, i
    d = size (mask)
    allocate (state_flv%pdg (d, n), source = 0)
    allocate (state_flv%map (d, n), source = spread ([(i, i = 1, d)], 2, n))
    allocate (state_flv%mask (d), source = mask)
  end subroutine state_flv_content_init

  module subroutine state_flv_content_set_entry (state_flv, i, pdg, map)
    class(state_flv_content_t), intent(inout) :: state_flv
    integer, intent(in) :: i
    integer, dimension(:), intent(in) :: pdg, map
    state_flv%pdg(:,i) = pdg
    where (map /= 0)
       state_flv%map(:,i) = map
    end where
  end subroutine state_flv_content_set_entry

  module subroutine state_flv_content_fill &
       (state_flv, state_full, mask)
    class(state_flv_content_t), intent(out) :: state_flv
    type(state_matrix_t), intent(in), target :: state_full
    logical, dimension(:), intent(in) :: mask
    type(state_matrix_t), target :: state_tmp
    type(state_iterator_t) :: it
    type(flavor_t), dimension(:), allocatable :: flv
    integer, dimension(:), allocatable :: pdg, pdg_subset
    integer, dimension(:), allocatable :: idx, map_subset, idx_subset, map
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer :: n, d, c, i
    call state_tmp%init ()
    d = state_full%get_depth ()
    allocate (flv (d), qn (d), pdg (d), idx (d), map (d))
    idx = [(i, i = 1, d)]
    c = count (mask)
    allocate (pdg_subset (c), map_subset (c), idx_subset (c))
    call it%init (state_full)
    do while (it%is_valid ())
       flv = it%get_flavor ()
       call qn%init (flv)
       call state_tmp%add_state (qn)
       call it%advance ()
    end do
    n = state_tmp%get_n_leaves ()
    call state_flv%init (n, mask)
    i = 0
    call it%init (state_tmp)
    do while (it%is_valid ())
       i = i + 1
       flv = it%get_flavor ()
       pdg = flv%get_pdg ()
       idx_subset = pack (idx, mask)
       pdg_subset = pack (pdg, mask)
       map_subset = order_abs (pdg_subset)
       map = unpack (idx_subset (map_subset), mask, idx)
       call state_flv%set_entry (i, &
            unpack (pdg_subset(map_subset), mask, pdg), &
            order (map))
       call it%advance ()
    end do
    call state_tmp%final ()
  end subroutine state_flv_content_fill

  module subroutine state_flv_content_match (state_flv, pdg, success, map)
    class(state_flv_content_t), intent(in) :: state_flv
    integer, dimension(:), intent(in) :: pdg
    logical, intent(out) :: success
    integer, dimension(:), intent(out) :: map
    integer, dimension(:), allocatable :: pdg_subset, pdg_sorted, map1, map2
    integer, dimension(:), allocatable :: idx, map_subset, idx_subset
    integer :: i, n, c, d
    c = count (state_flv%mask)
    d = size (state_flv%pdg, 1)
    n = size (state_flv%pdg, 2)
    allocate (idx (d), source = [(i, i = 1, d)])
    allocate (idx_subset (c), pdg_subset (c), map_subset (c))
    allocate (pdg_sorted (d), map1 (d), map2 (d))
    idx_subset = pack (idx, state_flv%mask)
    pdg_subset = pack (pdg, state_flv%mask)
    map_subset = order_abs (pdg_subset)
    pdg_sorted = unpack (pdg_subset(map_subset), state_flv%mask, pdg)
    success = .false.
    do i = 1, n
       if (all (pdg_sorted == state_flv%pdg(:,i) &
            .or. pdg_sorted == 0)) then
          success = .true.
          exit
       end if
    end do
    if (success) then
       map1 = state_flv%map(:,i)
       map2 = unpack (idx_subset(map_subset), state_flv%mask, idx)
       map = map2(map1)
       where (pdg == 0)  map = 0
    end if
  end subroutine state_flv_content_match

  module function state_flv_content_contains (state_flv, pdg) result (success)
    class(state_flv_content_t), intent(in) :: state_flv
    integer, intent(in) :: pdg
    logical :: success
    success = any (state_flv%pdg == pdg)
  end function state_flv_content_contains

  elemental function pacify_complex (c_in) result (c_pac)
    complex(default), intent(in) :: c_in
    complex(default) :: c_pac
    c_pac = c_in
    if (real(c_pac) == -real(c_pac)) then
       c_pac = &
            cmplx (0._default, aimag(c_pac), kind=default)
    end if
    if (aimag(c_pac) == -aimag(c_pac)) then
       c_pac = &
            cmplx (real(c_pac), 0._default, kind=default)
    end if
  end function pacify_complex

  module subroutine state_flv_content_find_duplicates (state_flv, duplicate_mask)
    class(state_flv_content_t), intent(in) :: state_flv
    logical, dimension(:), allocatable, intent(out) :: duplicate_mask
    integer, dimension(:), allocatable :: flvst
    integer :: i1, i2, n_flvsts
    logical :: found_once
    n_flvsts = size (state_flv%pdg, 2)
    allocate (duplicate_mask (n_flvsts))
    duplicate_mask = .false.
    do i1 = 1, n_flvsts
       found_once = .false.
       flvst = state_flv%pdg(:,i1)
       do i2 = 1, i1
          if (all(flvst == state_flv%pdg(:,i2))) then
             if (found_once) then
                duplicate_mask(i1) = .true.
                exit
             else
                found_once = .true.
             end if
          end if
       end do
    end do
  end subroutine state_flv_content_find_duplicates


end submodule state_matrices_s

