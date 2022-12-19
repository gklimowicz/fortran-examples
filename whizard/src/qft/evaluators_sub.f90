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

submodule (evaluators) evaluators_s

  use io_units
  use format_defs, only: FMT_19
  use physics_defs, only: n_beams_rescaled
  use diagnostics
  use lorentz

  implicit none

contains

  elemental module subroutine pairing_array_init (pa, n, has_i2, has_factor)
    type(pairing_array_t), intent(out) :: pa
    integer, intent(in) :: n
    logical, intent(in) :: has_i2, has_factor
    allocate (pa%i1 (n))
    if (has_i2)  allocate (pa%i2 (n))
    if (has_factor)  allocate (pa%factor (n))
  end subroutine pairing_array_init

  module subroutine pairing_array_write (pa, unit)
    type(pairing_array_t), intent(in) :: pa
    integer, intent(in), optional :: unit
    integer :: i, u
    u = given_output_unit (unit); if (u < 0) return
    write (u, "(A)", advance = "no") "["
    if (allocated (pa%i1)) then
       write (u, "(I0,A)", advance = "no") pa%i1, ","
    else
       write (u, "(A)", advance = "no") "x,"
    end if
    if (allocated (pa%i2)) then
       write (u, "(I0,A)", advance = "no") pa%i1, ","
    else
       write (u, "(A)", advance = "no") "x,"
    end if
    write (u, "(A)", advance = "no") "]"
    if (allocated (pa%factor)) then
       write (u, "(A,F5.4,A,F5.4,A)") ";(", &
            real(pa%factor), ",", aimag(pa%factor), ")]"
    else
       write (u, "(A)") ""
    end if
  end subroutine pairing_array_write

  module subroutine evaluator_write (eval, unit, &
       verbose, show_momentum_sum, show_mass, show_state, show_table, &
       col_verbose, testflag)
    class(evaluator_t), intent(in) :: eval
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, show_momentum_sum, show_mass
    logical, intent(in), optional :: show_state, show_table, col_verbose
    logical, intent(in), optional :: testflag
    logical :: conjugate, square, show_tab
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    show_tab = .true.;  if (present (show_table))  show_tab = .false.
    call eval%basic_write &
         (unit, verbose, show_momentum_sum, show_mass, &
            show_state, col_verbose, testflag)
    if (show_tab) then
       write (u, "(1x,A)")  "Matrix-element multiplication"
       write (u, "(2x,A)", advance="no")  "Input interaction 1:"
       if (associated (eval%int_in1)) then
          write (u, "(1x,I0)")  eval%int_in1%get_tag ()
       else
          write (u, "(A)")  " [undefined]"
       end if
       write (u, "(2x,A)", advance="no")  "Input interaction 2:"
       if (associated (eval%int_in2)) then
          write (u, "(1x,I0)")  eval%int_in2%get_tag ()
       else
          write (u, "(A)")  " [undefined]"
       end if
       select case (eval%type)
       case (EVAL_SQUARED_FLOWS, EVAL_SQUARE_WITH_COLOR_FACTORS)
          conjugate = .true.
          square = .true.
       case (EVAL_IDENTITY)
          write (u, "(1X,A)") "Identity evaluator, pairing array unused"
          return
       case default
          conjugate = .false.
          square = .false.
       end select
       call eval%write_pairing_array (conjugate, square, u)
    end if
  end subroutine evaluator_write

  module subroutine evaluator_write_pairing_array (eval, conjugate, square, unit)
    class(evaluator_t), intent(in) :: eval
    logical, intent(in) :: conjugate, square
    integer, intent(in), optional :: unit
    integer :: u, i, j
    u = given_output_unit (unit);  if (u < 0)  return
    if (allocated (eval%pairing_array)) then
       do i = 1, size (eval%pairing_array)
          write (u, "(2x,A,I0,A)")  "ME(", i, ") = "
          do j = 1, size (eval%pairing_array(i)%i1)
             write (u, "(4x,A)", advance="no")  "+"
             if (allocated (eval%pairing_array(i)%i2)) then
                write (u, "(1x,A,I0,A)", advance="no")  &
                     "ME1(", eval%pairing_array(i)%i1(j), ")"
                if (conjugate) then
                   write (u, "(A)", advance="no")  "* x"
                else
                   write (u, "(A)", advance="no")  " x"
                end if
                write (u, "(1x,A,I0,A)", advance="no")  &
                     "ME2(", eval%pairing_array(i)%i2(j), ")"
             else if (square) then
                write (u, "(1x,A)", advance="no")  "|"
                write (u, "(A,I0,A)", advance="no")  &
                     "ME1(", eval%pairing_array(i)%i1(j), ")"
                write (u, "(A)", advance="no")  "|^2"
             else
                write (u, "(1x,A,I0,A)", advance="no")  &
                     "ME1(", eval%pairing_array(i)%i1(j), ")"
             end if
             if (allocated (eval%pairing_array(i)%factor)) then
                write (u, "(1x,A)", advance="no")  "x"
                write (u, "(1x,'('," // FMT_19 // ",','," // FMT_19 // &
                     ",')')") eval%pairing_array(i)%factor(j)
             else
                write (u, *)
             end if
          end do
       end do
    end if
  end subroutine evaluator_write_pairing_array

  module subroutine evaluator_assign (eval_out, eval_in)
    type(evaluator_t), intent(out) :: eval_out
    type(evaluator_t), intent(in) :: eval_in
    eval_out%type = eval_in%type
    eval_out%int_in1 => eval_in%int_in1
    eval_out%int_in2 => eval_in%int_in2
    eval_out%interaction_t = eval_in%interaction_t
    if (allocated (eval_in%pairing_array)) then
       allocate (eval_out%pairing_array (size (eval_in%pairing_array)))
       eval_out%pairing_array = eval_in%pairing_array
    end if
  end subroutine evaluator_assign

  elemental module subroutine index_map_init (map, n)
    type(index_map_t), intent(out) :: map
    integer, intent(in) :: n
    allocate (map%entry (n))
    map%entry = 0
  end subroutine index_map_init

  function index_map_exists (map) result (flag)
    logical :: flag
    type(index_map_t), intent(in) :: map
    flag = allocated (map%entry)
  end function index_map_exists

  module function index_map_size (map) result (s)
    integer :: s
    type(index_map_t), intent(in) :: map
    if (allocated (map%entry)) then
       s = size (map%entry)
    else
       s = 0
    end if
  end function index_map_size

  elemental module subroutine index_map_assign_int (map, ival)
    type(index_map_t), intent(inout) :: map
    integer, intent(in) :: ival
    map%entry = ival
  end subroutine index_map_assign_int

  module subroutine index_map_assign_array (map, array)
    type(index_map_t), intent(inout) :: map
    integer, dimension(:), intent(in) :: array
    map%entry = array
  end subroutine index_map_assign_array

  elemental module subroutine index_map_set_entry (map, i, ival)
    type(index_map_t), intent(inout) :: map
    integer, intent(in) :: i
    integer, intent(in) :: ival
    map%entry(i) = ival
  end subroutine index_map_set_entry

  elemental module function index_map_get_entry (map, i) result (ival)
    integer :: ival
    type(index_map_t), intent(in) :: map
    integer, intent(in) :: i
    ival = map%entry(i)
  end function index_map_get_entry

  elemental subroutine index_map2_init (map, n)
    type(index_map2_t), intent(out) :: map
    integer, intent(in) :: n
    map%s = n
    allocate (map%entry (n, n))
  end subroutine index_map2_init

  function index_map2_exists (map) result (flag)
    logical :: flag
    type(index_map2_t), intent(in) :: map
    flag = allocated (map%entry)
  end function index_map2_exists

  module function index_map2_size (map) result (s)
    integer :: s
    type(index_map2_t), intent(in) :: map
    s = map%s
  end function index_map2_size

  elemental module subroutine index_map2_assign_int (map, ival)
    type(index_map2_t), intent(inout) :: map
    integer, intent(in) :: ival
    map%entry = ival
  end subroutine index_map2_assign_int

  elemental subroutine index_map2_set_entry (map, i, j, ival)
    type(index_map2_t), intent(inout) :: map
    integer, intent(in) :: i, j
    integer, intent(in) :: ival
    map%entry(i,j) = ival
  end subroutine index_map2_set_entry

  elemental function index_map2_get_entry (map, i, j) result (ival)
    integer :: ival
    type(index_map2_t), intent(in) :: map
    integer, intent(in) :: i, j
    ival = map%entry(i,j)
  end function index_map2_get_entry

  subroutine prt_mask_init (mask, n)
    type(prt_mask_t), intent(out) :: mask
    integer, intent(in) :: n
    allocate (mask%entry (n))
  end subroutine prt_mask_init

  module function prt_mask_size (mask) result (s)
    integer :: s
    type(prt_mask_t), intent(in) :: mask
    s = size (mask%entry)
  end function prt_mask_size

  subroutine connection_entry_init &
      (entry, n_count, n_map, qn_conn, count, n_rest)
    type(connection_entry_t), intent(out) :: entry
    integer, intent(in) :: n_count, n_map
    type(quantum_numbers_t), dimension(:), intent(in) :: qn_conn
    integer, dimension(n_count), intent(in) :: count
    integer, dimension(n_count), intent(in) :: n_rest
    integer :: i
    allocate (entry%qn_conn (size (qn_conn)))
    allocate (entry%n_index (n_count))
    allocate (entry%count (n_count))
    allocate (entry%index_in (n_map))
    allocate (entry%qn_in_list (n_count))
    entry%qn_conn = qn_conn
    entry%n_index = count
    entry%count = 0
    if (size (entry%index_in) == size (count)) then
       call index_map_init (entry%index_in, count)
    else
       call index_map_init (entry%index_in, count(1))
    end if
    do i = 1, n_count
      allocate (entry%qn_in_list(i)%qn (n_rest(i), count(i)))
    end do
  end subroutine connection_entry_init

  subroutine connection_entry_write (entry, unit)
    type(connection_entry_t), intent(in) :: entry
    integer, intent(in), optional :: unit
    integer :: i, j
    integer :: u
    u = given_output_unit (unit)
    call quantum_numbers_write (entry%qn_conn, unit)
    write (u, *)
    do i = 1, size (entry%n_index)
       write (u, *)  "Input interaction", i
       do j = 1, entry%n_index(i)
          if (size (entry%n_index) == size (entry%index_in)) then
             write (u, "(2x,I0,4x,I0,2x)", advance = "no") &
                  j, index_map_get_entry (entry%index_in(i), j)
          else
             write (u, "(2x,I0,4x,I0,2x,I0,2x)", advance = "no") &
                  j, index_map_get_entry (entry%index_in(1), j), &
                     index_map_get_entry (entry%index_in(2), j)
          end if
          call quantum_numbers_write (entry%qn_in_list(i)%qn(:,j), unit)
          write (u, *)
       end do
    end do
  end subroutine connection_entry_write

  subroutine color_table_init (color_table, state, n_tot)
    type(color_table_t), intent(out) :: color_table
    type(state_matrix_t), intent(in) :: state
    integer, intent(in) :: n_tot
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    type(state_matrix_t) :: state_col
    integer :: index, n_col_state
    allocate (color_table%index (state%get_n_matrix_elements ()))
    color_table%index = 0
    allocate (qn (n_tot))
    call state_col%init ()
    call it%init (state)
    do while (it%is_valid ())
       index = it%get_me_index ()
       call qn%init (col = it%get_color ())
       call state_col%add_state (qn, me_index = color_table%index(index))
       call it%advance ()
    end do
    n_col_state = state_col%get_n_matrix_elements ()
    allocate (color_table%col (n_tot, n_col_state))
    call it%init (state_col)
    do while (it%is_valid ())
       index = it%get_me_index ()
       color_table%col(:,index) = it%get_color ()
       call it%advance ()
    end do
    call state_col%final ()
    allocate (color_table%factor_is_known (n_col_state, n_col_state))
    allocate (color_table%factor (n_col_state, n_col_state))
    color_table%factor_is_known = .false.
  end subroutine color_table_init

  subroutine color_table_write (color_table, unit)
    type(color_table_t), intent(in) :: color_table
    integer, intent(in), optional :: unit
    integer :: i, j
    integer :: u
    u = given_output_unit (unit)
    write (u, *) "Color table:"
    if (allocated (color_table%index)) then
       write (u, *) "  Index mapping state => color table:"
       do i = 1, size (color_table%index)
          write (u, "(3x,I0,2x,I0,2x)")  i, color_table%index(i)
       end do
       write (u, *) "  Color table:"
       do i = 1, size (color_table%col, 2)
          write (u, "(3x,I0,2x)", advance = "no")  i
          call color_write (color_table%col(:,i), unit)
          write (u, *)
       end do
       write (u, *) "  Defined color factors:"
       do i = 1, size (color_table%factor, 1)
          do j = 1, size (color_table%factor, 2)
             if (color_table%factor_is_known(i,j)) then
                write (u, *)  i, j, color_table%factor(i,j)
             end if
          end do
       end do
    end if
  end subroutine color_table_write

  subroutine color_table_set_color_factors (color_table, &
       col_flow_index, col_factor, col_index_hi)
    type(color_table_t), intent(inout) :: color_table
    integer, dimension(:,:), intent(in) :: col_flow_index
    complex(default), dimension(:), intent(in) :: col_factor
    integer, dimension(:), intent(in) :: col_index_hi
    integer, dimension(:), allocatable :: hi_to_ct
    integer :: n_cflow
    integer :: hi_index, me_index, ct_index, cf_index
    integer, dimension(2) :: hi_index_pair, ct_index_pair
    n_cflow = size (col_index_hi)
    if (size (color_table%index) /= n_cflow) &
         call msg_bug ("Mismatch between hard matrix element and color table")
    allocate (hi_to_ct (n_cflow))
    do me_index = 1, size (color_table%index)
       ct_index = color_table%index(me_index)
       hi_index = col_index_hi(me_index)
       hi_to_ct(hi_index) = ct_index
    end do
    do cf_index = 1, size (col_flow_index, 2)
       hi_index_pair = col_flow_index(:,cf_index)
       ct_index_pair = hi_to_ct(hi_index_pair)
       color_table%factor(ct_index_pair(1), ct_index_pair(2)) = &
            col_factor(cf_index)
       color_table%factor_is_known(ct_index_pair(1), ct_index_pair(2)) = .true.
    end do
  end subroutine color_table_set_color_factors

  function color_table_get_color_factor (color_table, index1, index2, nc) &
      result (factor)
    real(default) :: factor
    type(color_table_t), intent(inout) :: color_table
    integer, intent(in) :: index1, index2
    integer, intent(in), optional :: nc
    integer :: i1, i2
    i1 = color_table%index(index1)
    i2 = color_table%index(index2)
    if (color_table%factor_is_known(i1,i2)) then
       factor = real(color_table%factor(i1,i2), kind=default)
    else
       factor = compute_color_factor &
            (color_table%col(:,i1), color_table%col(:,i2), nc)
       color_table%factor(i1,i2) = factor
       color_table%factor_is_known(i1,i2) = .true.
    end if
  end function color_table_get_color_factor

  module subroutine evaluator_init_product &
       (eval, int_in1, int_in2, qn_mask_conn, qn_filter_conn, qn_mask_rest, &
        connections_are_resonant, ignore_sub_for_qn)

    class(evaluator_t), intent(out), target :: eval
    class(interaction_t), intent(in), target :: int_in1, int_in2
    type(quantum_numbers_mask_t), intent(in) :: qn_mask_conn
    type(quantum_numbers_t), intent(in), optional :: qn_filter_conn
    type(quantum_numbers_mask_t), intent(in), optional :: qn_mask_rest
    logical, intent(in), optional :: connections_are_resonant
    logical, intent(in), optional :: ignore_sub_for_qn

    type(qn_mask_array_t), dimension(2) :: qn_mask_in
    type(state_matrix_t), pointer :: state_in1, state_in2
    type(connection_table_t) :: connection_table

    integer :: n_in, n_vir, n_out, n_tot
    integer, dimension(2) :: n_rest
    integer :: n_conn

    integer, dimension(:,:), allocatable :: connection_index
    type(index_map_t), dimension(2) :: prt_map_in
    type(index_map_t) :: prt_map_conn
    type(prt_mask_t), dimension(2) :: prt_is_connected
    type(quantum_numbers_mask_t), dimension(:), allocatable :: &
         qn_mask_conn_initial, int_in1_mask, int_in2_mask

    integer :: i

    eval%type = EVAL_PRODUCT
    eval%int_in1 => int_in1
    eval%int_in2 => int_in2

    state_in1 => int_in1%get_state_matrix_ptr ()
    state_in2 => int_in2%get_state_matrix_ptr ()

    call find_connections (int_in1, int_in2, n_conn, connection_index)
    if (n_conn == 0) then
       call msg_message ("First interaction:")
       call int_in1%basic_write (col_verbose=.true.)
       call msg_message ("Second interaction:")
       call int_in2%basic_write (col_verbose=.true.)
       call msg_fatal ("Evaluator product: no connections found between factors")
    end if
    call compute_index_bounds_and_mappings &
         (int_in1, int_in2, n_conn, &
          n_in, n_vir, n_out, n_tot, &
          n_rest, prt_map_in, prt_map_conn)

    call prt_mask_init (prt_is_connected(1), int_in1%get_n_tot ())
    call prt_mask_init (prt_is_connected(2), int_in2%get_n_tot ())
    do i = 1, 2
       prt_is_connected(i)%entry = .true.
       prt_is_connected(i)%entry(connection_index(:,i)) = .false.
    end do
    allocate (qn_mask_conn_initial (n_conn), &
         int_in1_mask (n_conn), int_in2_mask (n_conn))
    int_in1_mask = int_in1%get_mask (connection_index(:,1))
    int_in2_mask = int_in2%get_mask (connection_index(:,2))
    do i = 1, n_conn
       qn_mask_conn_initial(i) = int_in1_mask(i) .or. int_in2_mask(i)
    end do
    allocate (qn_mask_in(1)%mask (int_in1%get_n_tot ()))
    allocate (qn_mask_in(2)%mask (int_in2%get_n_tot ()))
    qn_mask_in(1)%mask = int_in1%get_mask ()
    qn_mask_in(2)%mask = int_in2%get_mask ()
    call connection_table_init (connection_table, &
         state_in1, state_in2, &
         qn_mask_conn_initial,  &
         n_conn, connection_index, n_rest, &
         qn_filter_conn, ignore_sub_for_qn)
    call connection_table_fill (connection_table, &
         state_in1, state_in2, &
         connection_index, prt_is_connected)
    call make_product_interaction (eval%interaction_t, &
         n_in, n_vir, n_out, &
         connection_table, &
         prt_map_in, prt_is_connected, &
         qn_mask_in, qn_mask_conn_initial, &
         qn_mask_conn, qn_filter_conn, qn_mask_rest)
    call make_pairing_array (eval%pairing_array, &
         eval%get_n_matrix_elements (), &
         connection_table)
    call record_links (eval%interaction_t, &
         int_in1, int_in2, connection_index, prt_map_in, prt_map_conn, &
         prt_is_connected, connections_are_resonant)
    call connection_table_final (connection_table)

    if (eval%get_n_matrix_elements () == 0) then
       print *, "Evaluator product"
       print *, "First interaction"
       call int_in1%basic_write (col_verbose=.true.)
       print *
       print *, "Second interaction"
       call int_in2%basic_write (col_verbose=.true.)
       print *
       call msg_fatal ("Product of density matrices is empty", &
           [var_str ("   --------------------------------------------"), &
            var_str ("This happens when two density matrices are convoluted "), &
            var_str ("but the processes they belong to (e.g., production "), &
            var_str ("and decay) do not match. This could happen if the "), &
            var_str ("beam specification does not match the hard "), &
            var_str ("process. Or it may indicate a WHIZARD bug.")])
    end if

  contains

    subroutine compute_index_bounds_and_mappings &
         (int1, int2, n_conn, &
          n_in, n_vir, n_out, n_tot, &
          n_rest, prt_map_in, prt_map_conn)
      class(interaction_t), intent(in) :: int1, int2
      integer, intent(in) :: n_conn
      integer, intent(out) :: n_in, n_vir, n_out, n_tot
      integer, dimension(2), intent(out) :: n_rest
      type(index_map_t), dimension(2), intent(out) :: prt_map_in
      type(index_map_t), intent(out) :: prt_map_conn
      integer, dimension(:), allocatable :: index
      integer :: n_in1, n_vir1, n_out1
      integer :: n_in2, n_vir2, n_out2
      integer :: k
      n_in1  = int1%get_n_in  ()
      n_vir1 = int1%get_n_vir ()
      n_out1 = int1%get_n_out () - n_conn
      n_rest(1) = n_in1 + n_vir1 + n_out1
      n_in2  = int2%get_n_in  () - n_conn
      n_vir2 = int2%get_n_vir ()
      n_out2 = int2%get_n_out ()
      n_rest(2) = n_in2 + n_vir2 + n_out2
      n_in  = n_in1  + n_in2
      n_vir = n_vir1 + n_vir2 + n_conn
      n_out = n_out1 + n_out2
      n_tot = n_in + n_vir + n_out
      call index_map_init (prt_map_in, n_rest)
      call index_map_init (prt_map_conn, n_conn)
      allocate (index (n_tot))
      index = [ (i, i = 1, n_tot) ]
      prt_map_in(1)%entry(1 : n_in1) = index(  1 :   n_in1)
      k =     n_in1
      prt_map_in(2)%entry(1 : n_in2) = index(k + 1 : k + n_in2)
      k = k + n_in2
      prt_map_in(1)%entry(n_in1 + 1 : n_in1 + n_vir1) = index(k + 1 : k + n_vir1)
      k = k + n_vir1
      prt_map_in(2)%entry(n_in2 + 1 : n_in2 + n_vir2) = index(k + 1 : k + n_vir2)
      k = k + n_vir2
      prt_map_conn%entry = index(k + 1 : k + n_conn)
      k = k + n_conn
      prt_map_in(1)%entry(n_in1 + n_vir1 + 1 : n_rest(1)) = index(k + 1 : k + n_out1)
      k = k + n_out1
      prt_map_in(2)%entry(n_in2 + n_vir2 + 1 : n_rest(2)) = index(k + 1 : k + n_out2)
    end subroutine compute_index_bounds_and_mappings

    subroutine connection_table_init &
        (connection_table, state_in1, state_in2, qn_mask_conn, &
         n_conn, connection_index, n_rest, &
         qn_filter_conn, ignore_sub_for_qn_in)
      type(connection_table_t), intent(out) :: connection_table
      type(state_matrix_t), intent(in), target :: state_in1, state_in2
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask_conn
      integer, intent(in) :: n_conn
      integer, dimension(:,:), intent(in) :: connection_index
      integer, dimension(2), intent(in) :: n_rest
      type(quantum_numbers_t), intent(in), optional :: qn_filter_conn
      logical, intent(in), optional :: ignore_sub_for_qn_in
      integer, dimension(2) :: n_me_in
      type(state_iterator_t) :: it
      type(quantum_numbers_t), dimension(n_conn) :: qn
      integer :: i, me_index_in, me_index_conn, n_me_conn
      integer, dimension(2) :: me_count
      logical :: ignore_sub_for_qn, has_sub_qn
      integer :: i_beam_sub
      connection_table%n_conn = n_conn
      connection_table%n_rest = n_rest
      n_me_in(1) = state_in1%get_n_matrix_elements ()
      n_me_in(2) = state_in2%get_n_matrix_elements ()
      allocate (connection_table%index_conn (2))
      call index_map_init (connection_table%index_conn, n_me_in)
      call connection_table%state%init (n_counters = 2)
      do i = 1, 2
         select case (i)
         case (1);  call it%init (state_in1)
         case (2);  call it%init (state_in2)
         end select
         do while (it%is_valid ())
            qn = it%get_quantum_numbers (connection_index(:,i))
            call qn%undefine (qn_mask_conn)
            if (present (qn_filter_conn)) then
               if (.not. all (qn .match. qn_filter_conn)) then
                  call it%advance ();  cycle
               end if
            end if
            call quantum_numbers_canonicalize_color (qn)
            me_index_in = it%get_me_index ()
            ignore_sub_for_qn = .false.; if (present (ignore_sub_for_qn_in)) ignore_sub_for_qn = ignore_sub_for_qn_in
            has_sub_qn = .false.
            do i_beam_sub = 1, n_beams_rescaled
               has_sub_qn = has_sub_qn .or. any (qn%get_sub () == i_beam_sub)
            end do
            call connection_table%state%add_state (qn, &
                 counter_index = i, &
                 ignore_sub_for_qn = .not. (ignore_sub_for_qn .and. has_sub_qn), &
                 me_index = me_index_conn)
            call index_map_set_entry (connection_table%index_conn(i), &
                 me_index_in, me_index_conn)
            call it%advance ()
         end do
      end do
      n_me_conn = connection_table%state%get_n_matrix_elements ()
      connection_table%n_me_conn = n_me_conn
      allocate (connection_table%entry (n_me_conn))
      call it%init (connection_table%state)
      do while (it%is_valid ())
         i = it%get_me_index ()
         me_count = it%get_me_count ()
         call connection_entry_init (connection_table%entry(i), 2, 2, &
              it%get_quantum_numbers (), me_count, n_rest)
         call it%advance ()
      end do
    end subroutine connection_table_init

    subroutine connection_table_final (connection_table)
      type(connection_table_t), intent(inout) :: connection_table
      call connection_table%state%final ()
    end subroutine connection_table_final

    subroutine connection_table_write (connection_table, unit)
      type(connection_table_t), intent(in) :: connection_table
      integer, intent(in), optional :: unit
      integer :: i, j
      integer :: u
      u = given_output_unit (unit)
      write (u, *) "Connection table:"
      call connection_table%state%write (unit)
      if (allocated (connection_table%index_conn)) then
         write (u, *) "  Index mapping input => connection table:"
         do i = 1, size (connection_table%index_conn)
            write (u, *) "  Input state", i
            do j = 1, size (connection_table%index_conn(i))
               write (u, *)  j, &
                    index_map_get_entry (connection_table%index_conn(i), j)
            end do
         end do
      end if
      if (allocated (connection_table%entry)) then
         write (u, *) "  Connection table contents:"
         do i = 1, size (connection_table%entry)
            call connection_entry_write (connection_table%entry(i), unit)
         end do
      end if
      if (index_map_exists (connection_table%index_result)) then
         write (u, *) "  Index mapping connection table => output:"
         do i = 1, size (connection_table%index_result)
            write (u, *)  i, &
                 index_map_get_entry (connection_table%index_result, i)
         end do
      end if
    end subroutine connection_table_write

    subroutine connection_table_fill &
        (connection_table, state_in1, state_in2, &
         connection_index, prt_is_connected)
      type(connection_table_t), intent(inout) :: connection_table
      type(state_matrix_t), intent(in), target :: state_in1, state_in2
      integer, dimension(:,:), intent(in) :: connection_index
      type(prt_mask_t), dimension(2), intent(in) :: prt_is_connected
      type(state_iterator_t) :: it
      integer :: index_in, index_conn
      integer :: color_offset
      integer :: n_result_entries
      integer :: i, k
      color_offset = connection_table%state%get_max_color_value ()
      do i = 1, 2
         select case (i)
         case (1);  call it%init (state_in1)
         case (2);  call it%init (state_in2)
         end select
         do while (it%is_valid ())
            index_in = it%get_me_index ()
            index_conn = index_map_get_entry &
                              (connection_table%index_conn(i), index_in)
            if (index_conn /= 0) then
               call connection_entry_add_state &
                    (connection_table%entry(index_conn), i, &
                    index_in, it%get_quantum_numbers (), &
                    connection_index(:,i), prt_is_connected(i), &
                    color_offset)
            end if
            call it%advance ()
         end do
         color_offset = color_offset + state_in1%get_max_color_value ()
      end do
      n_result_entries = 0
      do k = 1, size (connection_table%entry)
         n_result_entries = &
              n_result_entries + product (connection_table%entry(k)%n_index)
      end do
      call index_map_init (connection_table%index_result, n_result_entries)
    end subroutine connection_table_fill

    subroutine connection_entry_add_state &
        (entry, i, index_in, qn_in, connection_index, prt_is_connected, &
         color_offset)
      type(connection_entry_t), intent(inout) :: entry
      integer, intent(in) :: i
      integer, intent(in) :: index_in
      type(quantum_numbers_t), dimension(:), intent(in) :: qn_in
      integer, dimension(:), intent(in) :: connection_index
      type(prt_mask_t), intent(in) :: prt_is_connected
      integer, intent(in) :: color_offset
      integer :: c
      integer, dimension(:,:), allocatable :: color_map
      entry%count(i) = entry%count(i) + 1
      c = entry%count(i)
      call make_color_map (color_map, &
           qn_in(connection_index), entry%qn_conn)
      call index_map_set_entry (entry%index_in(i), c, index_in)
      entry%qn_in_list(i)%qn(:,c) = pack (qn_in, prt_is_connected%entry)
      call quantum_numbers_translate_color &
           (entry%qn_in_list(i)%qn(:,c), color_map, color_offset)
    end subroutine connection_entry_add_state

    subroutine make_product_interaction (int, &
         n_in, n_vir, n_out, &
         connection_table, &
         prt_map_in, prt_is_connected, &
         qn_mask_in, qn_mask_conn_initial, &
         qn_mask_conn, qn_filter_conn, qn_mask_rest)
      type(interaction_t), intent(out), target :: int
      integer, intent(in) :: n_in, n_vir, n_out
      type(connection_table_t), intent(inout), target :: connection_table
      type(index_map_t), dimension(2), intent(in) :: prt_map_in
      type(prt_mask_t), dimension(2), intent(in) :: prt_is_connected
      type(qn_mask_array_t), dimension(2), intent(in) :: qn_mask_in
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: &
           qn_mask_conn_initial
      type(quantum_numbers_mask_t), intent(in) :: qn_mask_conn
      type(quantum_numbers_t), intent(in), optional :: qn_filter_conn
      type(quantum_numbers_mask_t), intent(in), optional :: qn_mask_rest
      type(index_map_t), dimension(2) :: prt_index_in
      type(index_map_t) :: prt_index_conn
      integer :: n_tot, n_conn
      integer, dimension(2) :: n_rest
      integer :: i, j, k, m
      type(quantum_numbers_t), dimension(:), allocatable :: qn
      type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask
      type(connection_entry_t), pointer :: entry
      integer :: result_index
      n_conn = connection_table%n_conn
      n_rest = connection_table%n_rest
      n_tot = sum (n_rest) + n_conn
      allocate (qn (n_tot), qn_mask (n_tot))
      do i = 1, 2
         call index_map_init (prt_index_in(i), n_rest(i))
         prt_index_in(i)%entry = &
              prt_map_in(i)%entry ([ (j, j = 1, n_rest(i)) ])
      end do
      call index_map_init (prt_index_conn, n_conn)
      prt_index_conn%entry = prt_map_conn%entry ([ (j, j = 1, n_conn) ])
      do i = 1, 2
         if (present (qn_mask_rest)) then
            qn_mask(prt_index_in(i)%entry) = &
                 pack (qn_mask_in(i)%mask, prt_is_connected(i)%entry) &
                 .or. qn_mask_rest
         else
            qn_mask(prt_index_in(i)%entry) = &
                 pack (qn_mask_in(i)%mask, prt_is_connected(i)%entry)
         end if
      end do
      qn_mask(prt_index_conn%entry) = qn_mask_conn_initial .or. qn_mask_conn
      call eval%interaction_t%basic_init (n_in, n_vir, n_out, mask = qn_mask)
      m = 1
      do i = 1, connection_table%n_me_conn
         entry => connection_table%entry(i)
         qn(prt_index_conn%entry) = &
              quantum_numbers_undefined (entry%qn_conn, qn_mask_conn)
         if (present (qn_filter_conn)) then
            if (.not. all (qn(prt_index_conn%entry) .match. qn_filter_conn)) &
                 cycle
         end if
         do j = 1, entry%n_index(1)
            qn(prt_index_in(1)%entry) = entry%qn_in_list(1)%qn(:,j)
            do k = 1, entry%n_index(2)
               qn(prt_index_in(2)%entry) = entry%qn_in_list(2)%qn(:,k)
               call int%add_state (qn, me_index = result_index)
               call index_map_set_entry &
                    (connection_table%index_result, m, result_index)
               m = m + 1
            end do
         end do
      end do
      call int%freeze ()
    end subroutine make_product_interaction

    subroutine make_pairing_array (pa, n_matrix_elements, connection_table)
      type(pairing_array_t), dimension(:), intent(out), allocatable :: pa
      integer, intent(in) :: n_matrix_elements
      type(connection_table_t), intent(in), target :: connection_table
      type(connection_entry_t), pointer :: entry
      integer, dimension(:), allocatable :: n_entries
      integer :: i, j, k, m, r
      allocate (pa (n_matrix_elements))
      allocate (n_entries (n_matrix_elements))
      n_entries = 0
      do m = 1, size (connection_table%index_result)
         r = index_map_get_entry (connection_table%index_result, m)
         n_entries(r) = n_entries(r) + 1
      end do
      call pairing_array_init &
           (pa, n_entries, has_i2=.true., has_factor=.false.)
      m = 1
      n_entries = 0
      do i = 1, connection_table%n_me_conn
         entry => connection_table%entry(i)
         do j = 1, entry%n_index(1)
            do k = 1, entry%n_index(2)
               r = index_map_get_entry (connection_table%index_result, m)
               n_entries(r) = n_entries(r) + 1
               pa(r)%i1(n_entries(r)) = &
                    index_map_get_entry (entry%index_in(1), j)
               pa(r)%i2(n_entries(r)) = &
                    index_map_get_entry (entry%index_in(2), k)
               m = m + 1
            end do
         end do
      end do
    end subroutine make_pairing_array

    subroutine record_links (int, &
         int_in1, int_in2, connection_index, prt_map_in, prt_map_conn, &
         prt_is_connected, connections_are_resonant)
      class(interaction_t), intent(inout) :: int
      class(interaction_t), intent(in), target :: int_in1, int_in2
      integer, dimension(:,:), intent(in) :: connection_index
      type(index_map_t), dimension(2), intent(in) :: prt_map_in
      type(index_map_t), intent(in) :: prt_map_conn
      type(prt_mask_t), dimension(2), intent(in) :: prt_is_connected
      logical, intent(in), optional :: connections_are_resonant
      type(index_map_t), dimension(2) :: prt_map_all
      integer :: i, j, k, ival
      call index_map_init (prt_map_all(1), size (prt_is_connected(1)))
      k = 0
      j = 0
      do i = 1, size (prt_is_connected(1))
         if (prt_is_connected(1)%entry(i)) then
            j = j + 1
            ival = index_map_get_entry (prt_map_in(1), j)
            call index_map_set_entry (prt_map_all(1), i, ival)
         else
            k = k + 1
            ival = index_map_get_entry (prt_map_conn, k)
            call index_map_set_entry (prt_map_all(1), i, ival)
         end if
         call int%set_source_link (ival, int_in1, i)
      end do
      call int_in1%transfer_relations (int, prt_map_all(1)%entry)
      call index_map_init (prt_map_all(2), size (prt_is_connected(2)))
      j = 0
      do i = 1, size (prt_is_connected(2))
         if (prt_is_connected(2)%entry(i)) then
            j = j + 1
            ival = index_map_get_entry (prt_map_in(2), j)
            call index_map_set_entry (prt_map_all(2), i, ival)
            call int%set_source_link (ival, int_in2, i)
         else
            call index_map_set_entry (prt_map_all(2), i, 0)
         end if
      end do
      call int_in2%transfer_relations (int, prt_map_all(2)%entry)
      call int%relate_connections &
           (int_in2, connection_index(:,2), prt_map_all(2)%entry, &
           prt_map_conn%entry, connections_are_resonant)
    end subroutine record_links

  end subroutine evaluator_init_product

  module subroutine evaluator_init_square (eval, int_in, qn_mask, &
       col_flow_index, col_factor, col_index_hi, expand_color_flows, nc)
    class(evaluator_t), intent(out), target :: eval
    class(interaction_t), intent(in), target :: int_in
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
    integer, dimension(:,:), intent(in), optional :: col_flow_index
    complex(default), dimension(:), intent(in), optional :: col_factor
    integer, dimension(:), intent(in), optional :: col_index_hi
    logical, intent(in), optional :: expand_color_flows
    integer, intent(in), optional :: nc
    if (all (qn_mask%diagonal_helicity ())) then
       call eval%init_square_diag (int_in, qn_mask, &
            col_flow_index, col_factor, col_index_hi, expand_color_flows, nc)
    else
       call eval%init_square_nondiag (int_in, qn_mask, &
            col_flow_index, col_factor, col_index_hi, expand_color_flows, nc)
    end if
  end subroutine evaluator_init_square

  module subroutine evaluator_init_square_diag (eval, int_in, qn_mask, &
       col_flow_index, col_factor, col_index_hi, expand_color_flows, nc)

    class(evaluator_t), intent(out), target :: eval
    class(interaction_t), intent(in), target :: int_in
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
    integer, dimension(:,:), intent(in), optional :: col_flow_index
    complex(default), dimension(:), intent(in), optional :: col_factor
    integer, dimension(:), intent(in), optional :: col_index_hi
    logical, intent(in), optional :: expand_color_flows
    integer, intent(in), optional :: nc

    integer :: n_in, n_vir, n_out, n_tot
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask_initial
    type(state_matrix_t), pointer :: state_in
    type(connection_table_diag_t) :: connection_table

    logical :: sum_colors
    type(color_table_t) :: color_table

    if (present (expand_color_flows)) then
       sum_colors = .not. expand_color_flows
    else
       sum_colors = .true.
    end if

    if (sum_colors) then
       eval%type = EVAL_SQUARE_WITH_COLOR_FACTORS
    else
       eval%type = EVAL_SQUARED_FLOWS
    end if
    eval%int_in1 => int_in

    n_in  = int_in%get_n_in  ()
    n_vir = int_in%get_n_vir ()
    n_out = int_in%get_n_out ()
    n_tot = int_in%get_n_tot ()

    state_in => int_in%get_state_matrix_ptr ()

    allocate (qn_mask_initial (n_tot))
    qn_mask_initial = int_in%get_mask ()
    call qn_mask_initial%set_color (sum_colors, mask_cg=.false.)
    if (sum_colors) then
       call color_table_init (color_table, state_in, n_tot)
       if (present (col_flow_index) .and. present (col_factor) &
           .and. present (col_index_hi)) then
          call color_table_set_color_factors &
               (color_table, col_flow_index, col_factor, col_index_hi)
       end if
    end if

    call connection_table_init (connection_table, state_in, &
         qn_mask_initial, qn_mask, n_tot)
    call connection_table_fill (connection_table, state_in)
    call make_squared_interaction (eval%interaction_t, &
         n_in, n_vir, n_out, n_tot, &
         connection_table, sum_colors, qn_mask_initial .or. qn_mask)
    call make_pairing_array (eval%pairing_array, &
         eval%get_n_matrix_elements (), &
         connection_table, sum_colors, color_table, n_in, n_tot, nc)
    call record_links (eval, int_in, n_tot)
    call connection_table_final (connection_table)

  contains

    subroutine connection_table_init &
         (connection_table, state_in, qn_mask_in, qn_mask, n_tot)
      type(connection_table_diag_t), intent(out) :: connection_table
      type(state_matrix_t), intent(in), target :: state_in
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask_in
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
      integer, intent(in) :: n_tot
      type(quantum_numbers_t), dimension(n_tot) :: qn
      type(state_iterator_t) :: it
      integer :: i, n_me_in, me_index_in
      integer :: me_index_conn, n_me_conn
      integer, dimension(1) :: me_count
      logical :: qn_passed
      connection_table%n_tot = n_tot
      n_me_in = state_in%get_n_matrix_elements ()
      call index_map_init (connection_table%index_conn, n_me_in)
      call connection_table%state%init (n_counters=1)
      call it%init (state_in)
      do while (it%is_valid ())
         qn = it%get_quantum_numbers ()
         if (all (quantum_numbers_are_physical (qn, qn_mask))) then
            call qn%undefine (qn_mask_in)
            qn_passed = .true.
            if (qn_passed) then
               me_index_in = it%get_me_index ()
               call connection_table%state%add_state (qn, &
                    counter_index = 1, me_index = me_index_conn)
               call index_map_set_entry (connection_table%index_conn, &
                    me_index_in, me_index_conn)
            end if
         end if
         call it%advance ()
      end do
      n_me_conn = connection_table%state%get_n_matrix_elements ()
      connection_table%n_me_conn = n_me_conn
      allocate (connection_table%entry (n_me_conn))
      call it%init (connection_table%state)
      do while (it%is_valid ())
         i = it%get_me_index ()
         me_count = it%get_me_count ()
         call connection_entry_init (connection_table%entry(i), 1, 2, &
              it%get_quantum_numbers (), me_count, [n_tot])
         call it%advance ()
      end do
    end subroutine connection_table_init

    subroutine connection_table_final (connection_table)
      type(connection_table_diag_t), intent(inout) :: connection_table
      call connection_table%state%final ()
    end subroutine connection_table_final

    subroutine connection_table_write (connection_table, unit)
      type(connection_table_diag_t), intent(in) :: connection_table
      integer, intent(in), optional :: unit
      integer :: i
      integer :: u
      u = given_output_unit (unit)
      write (u, *) "Connection table:"
      call connection_table%state%write (unit)
      if (index_map_exists (connection_table%index_conn)) then
         write (u, *) "  Index mapping input => connection table:"
         do i = 1, size (connection_table%index_conn)
            write (u, *)  i, &
                   index_map_get_entry (connection_table%index_conn, i)
         end do
      end if
      if (allocated (connection_table%entry)) then
         write (u, *) "  Connection table contents"
         do i = 1, size (connection_table%entry)
            call connection_entry_write (connection_table%entry(i), unit)
         end do
      end if
      if (index_map_exists (connection_table%index_result)) then
         write (u, *) "  Index mapping connection table => output"
         do i = 1, size (connection_table%index_result)
            write (u, *)  i, &
                  index_map_get_entry (connection_table%index_result, i)
         end do
      end if
    end subroutine connection_table_write

    subroutine connection_table_fill (connection_table, state)
      type(connection_table_diag_t), intent(inout) :: connection_table
      type(state_matrix_t), intent(in), target :: state
      integer :: index_in, index_conn, n_result_entries
      type(state_iterator_t) :: it
      integer :: k
      call it%init (state)
      do while (it%is_valid ())
         index_in = it%get_me_index ()
         index_conn = &
              index_map_get_entry (connection_table%index_conn, index_in)
         if (index_conn /= 0) then
            call connection_entry_add_state &
                 (connection_table%entry(index_conn), &
                 index_in, it%get_quantum_numbers ())
         end if
         call it%advance ()
      end do
      n_result_entries = 0
      do k = 1, size (connection_table%entry)
         n_result_entries = &
              n_result_entries + connection_table%entry(k)%n_index(1) ** 2
      end do
      call index_map_init (connection_table%index_result, n_result_entries)
    end subroutine connection_table_fill

    subroutine connection_entry_add_state (entry, index_in, qn_in)
      type(connection_entry_t), intent(inout) :: entry
      integer, intent(in) :: index_in
      type(quantum_numbers_t), dimension(:), intent(in) :: qn_in
      integer :: c
      entry%count = entry%count + 1
      c = entry%count(1)
      call index_map_set_entry (entry%index_in(1), c, index_in)
      entry%qn_in_list(1)%qn(:,c) = qn_in
    end subroutine connection_entry_add_state

    subroutine make_squared_interaction (int, &
         n_in, n_vir, n_out, n_tot, &
         connection_table, sum_colors, qn_mask)
      type(interaction_t), intent(out), target :: int
      integer, intent(in) :: n_in, n_vir, n_out, n_tot
      type(connection_table_diag_t), intent(inout), target :: connection_table
      logical, intent(in) :: sum_colors
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
      type(connection_entry_t), pointer :: entry
      integer :: result_index, n_contrib
      integer :: i, m
      type(quantum_numbers_t), dimension(n_tot) :: qn
      call eval%interaction_t%basic_init (n_in, n_vir, n_out, mask=qn_mask)
      m = 0
      do i = 1, connection_table%n_me_conn
         entry => connection_table%entry(i)
         qn = quantum_numbers_undefined (entry%qn_conn, qn_mask)
         if (.not. sum_colors)   call qn(1:n_in)%invert_color ()
         call int%add_state (qn, me_index = result_index)
         n_contrib = entry%n_index(1) ** 2
         connection_table%index_result%entry(m+1:m+n_contrib) = result_index
         m = m + n_contrib
      end do
      call int%freeze ()
    end subroutine make_squared_interaction

    subroutine make_pairing_array (pa, &
         n_matrix_elements, connection_table, sum_colors, color_table, &
         n_in, n_tot, nc)
      type(pairing_array_t), dimension(:), intent(out), allocatable :: pa
      integer, intent(in) :: n_matrix_elements
      type(connection_table_diag_t), intent(in), target :: connection_table
      logical, intent(in) :: sum_colors
      type(color_table_t), intent(inout) :: color_table
      type(connection_entry_t), pointer :: entry
      integer, intent(in) :: n_in, n_tot
      integer, intent(in), optional :: nc
      integer, dimension(:), allocatable :: n_entries
      integer :: i, k, l, ks, ls, m, r
      integer :: color_multiplicity_in
      allocate (pa (n_matrix_elements))
      allocate (n_entries (n_matrix_elements))
      n_entries = 0
      do m = 1, size (connection_table%index_result)
         r = index_map_get_entry (connection_table%index_result, m)
         n_entries(r) = n_entries(r) + 1
      end do
      call pairing_array_init &
           (pa, n_entries, has_i2 = sum_colors, has_factor = sum_colors)
      m = 1
      n_entries = 0
      do i = 1, connection_table%n_me_conn
         entry => connection_table%entry(i)
         do k = 1, entry%n_index(1)
            if (sum_colors) then
               color_multiplicity_in = product (abs &
                    (entry%qn_in_list(1)%qn(:n_in, k)%get_color_type ()))
               do l = 1, entry%n_index(1)
                  r = index_map_get_entry (connection_table%index_result, m)
                  n_entries(r) = n_entries(r) + 1
                  ks = index_map_get_entry (entry%index_in(1), k)
                  ls = index_map_get_entry (entry%index_in(1), l)
                  pa(r)%i1(n_entries(r)) = ks
                  pa(r)%i2(n_entries(r)) = ls
                  pa(r)%factor(n_entries(r)) = &
                       color_table_get_color_factor (color_table, ks, ls, nc) &
                       / color_multiplicity_in
                  m = m + 1
               end do
            else
               r = index_map_get_entry (connection_table%index_result, m)
               n_entries(r) = n_entries(r) + 1
               ks = index_map_get_entry (entry%index_in(1), k)
               pa(r)%i1(n_entries(r)) = ks
               m = m + 1
            end if
         end do
      end do
    end subroutine make_pairing_array

    subroutine record_links (int, int_in, n_tot)
      class(interaction_t), intent(inout) :: int
      class(interaction_t), intent(in), target :: int_in
      integer, intent(in) :: n_tot
      integer, dimension(n_tot) :: map
      integer :: i
      do i = 1, n_tot
         call int%set_source_link (i, int_in, i)
      end do
      map = [ (i, i = 1, n_tot) ]
      call int_in%transfer_relations (int, map)
    end subroutine record_links

  end subroutine evaluator_init_square_diag

  module subroutine evaluator_init_square_nondiag (eval, int_in, qn_mask, &
       col_flow_index, col_factor, col_index_hi, expand_color_flows, nc)

    class(evaluator_t), intent(out), target :: eval
    class(interaction_t), intent(in), target :: int_in
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
    integer, dimension(:,:), intent(in), optional :: col_flow_index
    complex(default), dimension(:), intent(in), optional :: col_factor
    integer, dimension(:), intent(in), optional :: col_index_hi
    logical, intent(in), optional :: expand_color_flows
    integer, intent(in), optional :: nc

    integer :: n_in, n_vir, n_out, n_tot
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask_initial
    type(state_matrix_t), pointer :: state_in
    type(connection_table_nondiag_t) :: connection_table

    logical :: sum_colors
    type(color_table_t) :: color_table

    if (present (expand_color_flows)) then
       sum_colors = .not. expand_color_flows
    else
       sum_colors = .true.
    end if

    if (sum_colors) then
       eval%type = EVAL_SQUARE_WITH_COLOR_FACTORS
    else
       eval%type = EVAL_SQUARED_FLOWS
    end if
    eval%int_in1 => int_in

    n_in  = int_in%get_n_in  ()
    n_vir = int_in%get_n_vir ()
    n_out = int_in%get_n_out ()
    n_tot = int_in%get_n_tot ()

    state_in => int_in%get_state_matrix_ptr ()

    allocate (qn_mask_initial (n_tot))
    qn_mask_initial = int_in%get_mask ()
    call qn_mask_initial%set_color (sum_colors, mask_cg=.false.)
    if (sum_colors) then
       call color_table_init (color_table, state_in, n_tot)
       if (present (col_flow_index) .and. present (col_factor) &
           .and. present (col_index_hi)) then
          call color_table_set_color_factors &
               (color_table, col_flow_index, col_factor, col_index_hi)
       end if
    end if

    call connection_table_init (connection_table, state_in, &
         qn_mask_initial, qn_mask, n_tot)
    call connection_table_fill (connection_table, state_in)
    call make_squared_interaction (eval%interaction_t, &
         n_in, n_vir, n_out, n_tot, &
         connection_table, sum_colors, qn_mask_initial .or. qn_mask)
    call make_pairing_array (eval%pairing_array, &
         eval%get_n_matrix_elements (), &
         connection_table, sum_colors, color_table, n_in, n_tot, nc)
    call record_links (eval, int_in, n_tot)
    call connection_table_final (connection_table)

  contains

    subroutine connection_table_init &
         (connection_table, state_in, qn_mask_in, qn_mask, n_tot)
      type(connection_table_nondiag_t), intent(out) :: connection_table
      type(state_matrix_t), intent(in), target :: state_in
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask_in
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
      integer, intent(in) :: n_tot
      type(quantum_numbers_t), dimension(n_tot) :: qn1, qn2, qn
      type(state_iterator_t) :: it1, it2, it
      integer :: i, n_me_in, me_index_in1, me_index_in2
      integer :: me_index_conn, n_me_conn
      integer, dimension(1) :: me_count
      logical :: qn_passed
      connection_table%n_tot = n_tot
      n_me_in = state_in%get_n_matrix_elements ()
      call index_map2_init (connection_table%index_conn, n_me_in)
      connection_table%index_conn = 0
      call connection_table%state%init (n_counters=1)
      call it1%init (state_in)
      do while (it1%is_valid ())
         qn1 = it1%get_quantum_numbers ()
         me_index_in1 = it1%get_me_index ()
         call it2%init (state_in)
         do while (it2%is_valid ())
            qn2 = it2%get_quantum_numbers ()
            if (all (quantum_numbers_are_compatible (qn1, qn2, qn_mask))) then
               qn = qn1 .merge. qn2
               call qn%undefine (qn_mask_in)
               qn_passed = .true.
               if (qn_passed) then
                  me_index_in2 = it2%get_me_index ()
                  call connection_table%state%add_state (qn, &
                       counter_index = 1, me_index = me_index_conn)
                  call index_map2_set_entry (connection_table%index_conn, &
                       me_index_in1, me_index_in2, me_index_conn)
               end if
            end if
            call it2%advance ()
         end do
         call it1%advance ()
      end do
      n_me_conn = connection_table%state%get_n_matrix_elements ()
      connection_table%n_me_conn = n_me_conn
      allocate (connection_table%entry (n_me_conn))
      call it%init (connection_table%state)
      do while (it%is_valid ())
         i = it%get_me_index ()
         me_count = it%get_me_count ()
         call connection_entry_init (connection_table%entry(i), 1, 2, &
              it%get_quantum_numbers (), me_count, [n_tot])
         call it%advance ()
      end do
    end subroutine connection_table_init

    subroutine connection_table_final (connection_table)
      type(connection_table_nondiag_t), intent(inout) :: connection_table
      call connection_table%state%final ()
    end subroutine connection_table_final

    subroutine connection_table_write (connection_table, unit)
      type(connection_table_nondiag_t), intent(in) :: connection_table
      integer, intent(in), optional :: unit
      integer :: i, j
      integer :: u
      u = given_output_unit (unit)
      write (u, *) "Connection table:"
      call connection_table%state%write (unit)
      if (index_map2_exists (connection_table%index_conn)) then
         write (u, *) "  Index mapping input => connection table:"
         do i = 1, size (connection_table%index_conn)
            do j = 1, size (connection_table%index_conn)
               write (u, *)  i, j, &
                    index_map2_get_entry (connection_table%index_conn, i, j)
            end do
         end do
      end if
      if (allocated (connection_table%entry)) then
         write (u, *) "  Connection table contents"
         do i = 1, size (connection_table%entry)
            call connection_entry_write (connection_table%entry(i), unit)
         end do
      end if
      if (index_map_exists (connection_table%index_result)) then
         write (u, *) "  Index mapping connection table => output"
         do i = 1, size (connection_table%index_result)
            write (u, *)  i, &
                 index_map_get_entry (connection_table%index_result, i)
         end do
      end if
    end subroutine connection_table_write

    subroutine connection_table_fill (connection_table, state)
      type(connection_table_nondiag_t), intent(inout), target :: connection_table
      type(state_matrix_t), intent(in), target :: state
      integer :: index1_in, index2_in, index_conn, n_result_entries
      type(state_iterator_t) :: it1, it2
      integer :: k
      call it1%init (state)
      do while (it1%is_valid ())
         index1_in = it1%get_me_index ()
         call it2%init (state)
         do while (it2%is_valid ())
            index2_in = it2%get_me_index ()
            index_conn = index_map2_get_entry &
                            (connection_table%index_conn, index1_in, index2_in)
            if (index_conn /= 0) then
               call connection_entry_add_state &
                    (connection_table%entry(index_conn), &
                     index1_in, index2_in, &
                     it1%get_quantum_numbers () &
                     .merge. &
                     it2%get_quantum_numbers ())
            end if
            call it2%advance ()
         end do
         call it1%advance ()
      end do
      n_result_entries = 0
      do k = 1, size (connection_table%entry)
         n_result_entries = &
              n_result_entries + connection_table%entry(k)%n_index(1)
      end do
      call index_map_init (connection_table%index_result, n_result_entries)
    end subroutine connection_table_fill

    subroutine connection_entry_add_state (entry, index1_in, index2_in, qn_in)
      type(connection_entry_t), intent(inout) :: entry
      integer, intent(in) :: index1_in, index2_in
      type(quantum_numbers_t), dimension(:), intent(in) :: qn_in
      integer :: c
      entry%count = entry%count + 1
      c = entry%count(1)
      call index_map_set_entry (entry%index_in(1), c, index1_in)
      call index_map_set_entry (entry%index_in(2), c, index2_in)
      entry%qn_in_list(1)%qn(:,c) = qn_in
    end subroutine connection_entry_add_state

    subroutine make_squared_interaction (int, &
         n_in, n_vir, n_out, n_tot, &
         connection_table, sum_colors, qn_mask)
      type(interaction_t), intent(out), target :: int
      integer, intent(in) :: n_in, n_vir, n_out, n_tot
      type(connection_table_nondiag_t), intent(inout), target :: connection_table
      logical, intent(in) :: sum_colors
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
      type(connection_entry_t), pointer :: entry
      integer :: result_index
      integer :: i, k, m
      type(quantum_numbers_t), dimension(n_tot) :: qn
      call eval%interaction_t%basic_init (n_in, n_vir, n_out, mask=qn_mask)
      m = 0
      do i = 1, connection_table%n_me_conn
         entry => connection_table%entry(i)
         do k = 1, size (entry%qn_in_list(1)%qn, 2)
            qn = quantum_numbers_undefined &
                    (entry%qn_in_list(1)%qn(:,k), qn_mask)
            if (.not. sum_colors)  call qn(1:n_in)%invert_color ()
            call int%add_state (qn, me_index = result_index)
            call index_map_set_entry (connection_table%index_result, m + 1, &
                 result_index)
            m = m + 1
         end do
      end do
      call int%freeze ()
    end subroutine make_squared_interaction

    subroutine make_pairing_array (pa, &
         n_matrix_elements, connection_table, sum_colors, color_table, &
         n_in, n_tot, nc)
      type(pairing_array_t), dimension(:), intent(out), allocatable :: pa
      integer, intent(in) :: n_matrix_elements
      type(connection_table_nondiag_t), intent(in), target :: connection_table
      logical, intent(in) :: sum_colors
      type(color_table_t), intent(inout) :: color_table
      type(connection_entry_t), pointer :: entry
      integer, intent(in) :: n_in, n_tot
      integer, intent(in), optional :: nc
      integer, dimension(:), allocatable :: n_entries
      integer :: i, k, k1s, k2s, m, r
      integer :: color_multiplicity_in
      allocate (pa (n_matrix_elements))
      allocate (n_entries (n_matrix_elements))
      n_entries = 0
      do m = 1, size (connection_table%index_result)
         r = index_map_get_entry (connection_table%index_result, m)
         n_entries(r) = n_entries(r) + 1
      end do
      call pairing_array_init &
           (pa, n_entries, has_i2 = sum_colors, has_factor = sum_colors)
      m = 1
      n_entries = 0
      do i = 1, connection_table%n_me_conn
         entry => connection_table%entry(i)
         do k = 1, entry%n_index(1)
            r = index_map_get_entry (connection_table%index_result, m)
            n_entries(r) = n_entries(r) + 1
            if (sum_colors) then
               k1s = index_map_get_entry (entry%index_in(1), k)
               k2s = index_map_get_entry (entry%index_in(2), k)
               pa(r)%i1(n_entries(r)) = k1s
               pa(r)%i2(n_entries(r)) = k2s
               color_multiplicity_in = product (abs &
                    (entry%qn_in_list(1)%qn(:n_in, k)%get_color_type ()))
               pa(r)%factor(n_entries(r)) = &
                    color_table_get_color_factor (color_table, k1s, k2s, nc) &
                    / color_multiplicity_in
            else
               k1s = index_map_get_entry (entry%index_in(1), k)
               pa(r)%i1(n_entries(r)) = k1s
            end if
            m = m + 1
         end do
      end do
    end subroutine make_pairing_array

    subroutine record_links (int, int_in, n_tot)
      class(interaction_t), intent(inout) :: int
      class(interaction_t), intent(in), target :: int_in
      integer, intent(in) :: n_tot
      integer, dimension(n_tot) :: map
      integer :: i
      do i = 1, n_tot
         call int%set_source_link (i, int_in, i)
      end do
      map = [ (i, i = 1, n_tot) ]
      call int_in%transfer_relations (int, map)
    end subroutine record_links

  end subroutine evaluator_init_square_nondiag

  module subroutine evaluator_init_color_contractions (eval, int_in)
    class(evaluator_t), intent(out), target :: eval
    type(interaction_t), intent(in), target :: int_in
    integer :: n_in, n_vir, n_out, n_tot
    type(state_matrix_t) :: state_with_contractions
    integer, dimension(:), allocatable :: me_index
    integer, dimension(:), allocatable :: result_index
    eval%type = EVAL_COLOR_CONTRACTION
    eval%int_in1 => int_in
    n_in  = int_in%get_n_in  ()
    n_vir = int_in%get_n_vir ()
    n_out = int_in%get_n_out ()
    n_tot = int_in%get_n_tot ()
    state_with_contractions = int_in%get_state_matrix_ptr ()
    call state_with_contractions%add_color_contractions ()
    call make_contracted_interaction (eval%interaction_t, &
         me_index, result_index, &
         n_in, n_vir, n_out, n_tot, &
         state_with_contractions, int_in%get_mask ())
    call make_pairing_array (eval%pairing_array, me_index, result_index)
    call record_links (eval, int_in, n_tot)
    call state_with_contractions%final ()

  contains

    subroutine make_contracted_interaction (int, &
         me_index, result_index, &
         n_in, n_vir, n_out, n_tot, state, qn_mask)
      type(interaction_t), intent(out), target :: int
      integer, dimension(:), intent(out), allocatable :: me_index
      integer, dimension(:), intent(out), allocatable :: result_index
      integer, intent(in) :: n_in, n_vir, n_out, n_tot
      type(state_matrix_t), intent(in) :: state
      type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
      type(state_iterator_t) :: it
      integer :: n_me, i
      type(quantum_numbers_t), dimension(n_tot) :: qn
      call int%basic_init (n_in, n_vir, n_out, mask=qn_mask)
      n_me = state%get_n_leaves ()
      allocate (me_index (n_me))
      allocate (result_index (n_me))
      call it%init (state)
      i = 0
      do while (it%is_valid ())
         i = i + 1
         me_index(i) = it%get_me_index ()
         qn = it%get_quantum_numbers ()
         call int%add_state (qn, me_index = result_index(i))
         call it%advance ()
      end do
      call int%freeze ()
    end subroutine make_contracted_interaction

    subroutine make_pairing_array (pa, me_index, result_index)
      type(pairing_array_t), dimension(:), intent(out), allocatable :: pa
      integer, dimension(:), intent(in) :: me_index, result_index
      integer, dimension(:), allocatable :: n_entries
      integer :: n_matrix_elements, r, i, k
      !!! The result indices of the appended color contracted states
      !!! start counting from 1 again. For the pairing array, we currently
      !!! only take the first part of ascending indices into account
      !!! excluding the color contracted states.
      n_matrix_elements = size (me_index)
      k = 0
      do i = 1, n_matrix_elements
         r = result_index(i)
         if (r < i) exit
         k = r
      end do
      allocate (pa (k))
      allocate (n_entries (k))
      n_entries = 1
      call pairing_array_init &
           (pa, n_entries, has_i2=.false., has_factor=.false.)
      do i = 1, k
         r = result_index(i)
         pa(r)%i1(1) = me_index(i)
      end do
    end subroutine make_pairing_array

    subroutine record_links (int, int_in, n_tot)
      class(interaction_t), intent(inout) :: int
      class(interaction_t), intent(in), target :: int_in
      integer, intent(in) :: n_tot
      integer, dimension(n_tot) :: map
      integer :: i
      do i = 1, n_tot
         call int%set_source_link (i, int_in, i)
      end do
      map = [ (i, i = 1, n_tot) ]
      call int_in%transfer_relations (int, map)
    end subroutine record_links

  end subroutine evaluator_init_color_contractions

  module subroutine evaluator_reassign_links_eval (eval, eval_src, eval_target)
    type(evaluator_t), intent(inout) :: eval
    type(evaluator_t), intent(in) :: eval_src
    type(evaluator_t), intent(in), target :: eval_target
    if (associated (eval%int_in1)) then
       if (eval%int_in1%get_tag () == eval_src%get_tag ()) then
          eval%int_in1 => eval_target%interaction_t
       end if
    end if
    if (associated (eval%int_in2)) then
       if (eval%int_in2%get_tag () == eval_src%get_tag ()) then
          eval%int_in2 => eval_target%interaction_t
       end if
    end if
    call interaction_reassign_links &
         (eval%interaction_t, eval_src%interaction_t, &
         eval_target%interaction_t)
  end subroutine evaluator_reassign_links_eval

  module subroutine evaluator_reassign_links_int (eval, int_src, int_target)
    type(evaluator_t), intent(inout) :: eval
    type(interaction_t), intent(in) :: int_src
    type(interaction_t), intent(in), target :: int_target
    if (associated (eval%int_in1)) then
       if (eval%int_in1%get_tag () == int_src%get_tag ()) then
          eval%int_in1 => int_target
       end if
    end if
    if (associated (eval%int_in2)) then
       if (eval%int_in2%get_tag () == int_src%get_tag ()) then
          eval%int_in2 => int_target
       end if
    end if
    call interaction_reassign_links (eval%interaction_t, int_src, int_target)
  end subroutine evaluator_reassign_links_int

  module function evaluator_get_int_in_ptr (eval, i) result (int_in)
    class(interaction_t), pointer :: int_in
    type(evaluator_t), intent(in), target :: eval
    integer, intent(in) :: i
    if (i == 1) then
       int_in => eval%int_in1
    else if (i == 2) then
       int_in => eval%int_in2
    else
       int_in => null ()
    end if
  end function evaluator_get_int_in_ptr

  module subroutine evaluator_init_identity (eval, int)
    class(evaluator_t), intent(out), target :: eval
    class(interaction_t), intent(in), target :: int
    integer :: n_in, n_out, n_vir, n_tot
    integer :: i
    integer, dimension(:), allocatable :: map
    type(state_matrix_t), pointer :: state
    type(state_iterator_t) :: it
    eval%type = EVAL_IDENTITY
    eval%int_in1 => int
    nullify (eval%int_in2)
    n_in = int%get_n_in ()
    n_out = int%get_n_out ()
    n_vir = int%get_n_vir ()
    n_tot = int%get_n_tot ()
    call eval%interaction_t%basic_init (n_in, n_vir, n_out, &
       mask = int%get_mask (), &
       resonant = int%get_resonance_flags ())
    do i = 1, n_tot
       call eval%set_source_link (i, int, i)
    end do
    allocate (map(n_tot))
    map = [(i, i = 1, n_tot)]
    call int%transfer_relations (eval, map)
    state => int%get_state_matrix_ptr ()
    call it%init (state)
    do while (it%is_valid ())
       call eval%add_state (it%get_quantum_numbers (), &
            it%get_me_index ())
       call it%advance ()
    end do
    call eval%freeze ()

  end subroutine evaluator_init_identity

  module subroutine evaluator_init_qn_sum (eval, int, qn_mask, drop)
    class(evaluator_t), intent(out), target :: eval
    class(interaction_t), target, intent(in) :: int
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: qn_mask
    logical, intent(in), optional, dimension(:) :: drop
    type(state_iterator_t) :: it_old, it_new
    integer, dimension(:), allocatable :: pairing_size, pairing_target, i_new
    integer, dimension(:), allocatable :: map
    integer :: n_in, n_out, n_vir, n_tot, n_me_old, n_me_new
    integer :: i, j
    type(state_matrix_t), pointer :: state_new, state_old
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    logical :: matched
    logical, dimension(size (qn_mask)) :: dropped
    integer :: ndropped
    integer, dimension(:), allocatable :: inotdropped
    type(quantum_numbers_mask_t), dimension(:), allocatable :: mask
    logical, dimension(:), allocatable :: resonant

    eval%type = EVAL_QN_SUM
    eval%int_in1 => int
    nullify (eval%int_in2)
    if (present (drop)) then
       dropped = drop
    else
       dropped = .false.
    end if
    ndropped = count (dropped)

    n_in = int%get_n_in ()
    n_out = int%get_n_out () - ndropped
    n_vir = int%get_n_vir ()
    n_tot = int%get_n_tot () - ndropped

    allocate (inotdropped (n_tot))
    i = 1
    do j = 1, n_tot + ndropped
       if (dropped (j)) cycle
       inotdropped(i) = j
       i = i + 1
    end do

    allocate (mask(n_tot + ndropped))
    mask = int%get_mask ()
    allocate (resonant(n_tot + ndropped))
    resonant = int%get_resonance_flags ()
    call eval%interaction_t%basic_init (n_in, n_vir, n_out, &
         mask = mask(inotdropped) .or. qn_mask(inotdropped), &
         resonant = resonant(inotdropped))
    i = 1
    do j = 1, n_tot + ndropped
       if (dropped(j)) cycle
       call eval%set_source_link (i, int, j)
       i = i + 1
    end do
    allocate (map(n_tot + ndropped))
    i = 1
    do j = 1, n_tot + ndropped
       if (dropped (j)) then
          map(j) = 0
       else
          map(j) = i
          i = i + 1
       end if
    end do
    call int%transfer_relations (eval, map)

    n_me_old = int%get_n_matrix_elements ()
    allocate (pairing_size (n_me_old), source = 0)
    allocate (pairing_target (n_me_old), source = 0)
    pairing_size = 0
    state_old => int%get_state_matrix_ptr ()
    state_new => eval%get_state_matrix_ptr ()
    call it_old%init (state_old)
    allocate (qn(n_tot + ndropped))
    do while (it_old%is_valid ())
       qn = it_old%get_quantum_numbers ()
       if (.not. all (qn%are_diagonal ())) then
          call it_old%advance ()
          cycle
       end if
       matched = .false.
       call it_new%init (state_new)
       if (eval%get_n_matrix_elements () > 0) then
          do while (it_new%is_valid ())
             if (all (qn(inotdropped) .match. &
                it_new%get_quantum_numbers ())) &
             then
                matched = .true.
                i = it_new%get_me_index ()
                exit
             end if
             call it_new%advance ()
          end do
       end if
       if (.not. matched) then
          call eval%add_state (qn(inotdropped))
          i = eval%get_n_matrix_elements ()
       end if
       pairing_size(i) = pairing_size(i) + 1
       pairing_target(it_old%get_me_index ()) = i
       call it_old%advance ()
    end do
    call eval%freeze ()

    n_me_new = eval%get_n_matrix_elements ()
    allocate (eval%pairing_array (n_me_new))
    do i = 1, n_me_new
       call pairing_array_init (eval%pairing_array(i), &
            pairing_size(i), .false., .false.)
    end do

    allocate (i_new (n_me_new), source = 0)
    do i = 1, n_me_old
       j = pairing_target(i)
       if (j > 0) then
          i_new(j) = i_new(j) + 1
          eval%pairing_array(j)%i1(i_new(j)) = i
       end if
    end do

  end subroutine evaluator_init_qn_sum

  module subroutine evaluator_evaluate (eval)
    class(evaluator_t), intent(inout), target :: eval
    integer :: i
    select case (eval%type)
    case (EVAL_PRODUCT)
       do i = 1, size(eval%pairing_array)
          call eval%evaluate_product (i, &
               eval%int_in1, eval%int_in2, &
               eval%pairing_array(i)%i1, eval%pairing_array(i)%i2)
          if (debug2_active (D_QFT)) then
             print *, 'eval%pairing_array(i)%i1, eval%pairing_array(i)%i2 =    ', &
                  eval%pairing_array(i)%i1, eval%pairing_array(i)%i2
             print *, 'MEs =    ', &
                  eval%int_in1%get_matrix_element (eval%pairing_array(i)%i1), &
                  eval%int_in2%get_matrix_element (eval%pairing_array(i)%i2)
          end if
       end do
    case (EVAL_SQUARE_WITH_COLOR_FACTORS)
       do i = 1, size(eval%pairing_array)
          call eval%evaluate_product_cf (i, &
               eval%int_in1, eval%int_in1, &
               eval%pairing_array(i)%i1, eval%pairing_array(i)%i2, &
               eval%pairing_array(i)%factor)
       end do
    case (EVAL_SQUARED_FLOWS)
       do i = 1, size(eval%pairing_array)
          call eval%evaluate_square_c (i, &
               eval%int_in1, &
               eval%pairing_array(i)%i1)
       end do
    case (EVAL_COLOR_CONTRACTION)
       do i = 1, size(eval%pairing_array)
          call eval%evaluate_sum (i, &
               eval%int_in1, &
               eval%pairing_array(i)%i1)
       end do
    case (EVAL_IDENTITY)
       call eval%set_matrix_element (eval%int_in1)
    case (EVAL_QN_SUM)
       do i = 1, size (eval%pairing_array)
          call eval%evaluate_me_sum (i, &
             eval%int_in1, eval%pairing_array(i)%i1)
          call eval%set_norm (eval%int_in1%get_norm ())
       end do
    end select
  end subroutine evaluator_evaluate


end submodule evaluators_s

