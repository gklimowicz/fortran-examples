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

submodule (auto_components) auto_components_s

  use io_units
  use diagnostics
  use physics_defs, only: PHOTON, GLUON, Z_BOSON, W_BOSON
  use numeric_utils, only: extend_integer_array

  implicit none

  integer, parameter :: PROC_UNDEFINED = 0
  integer, parameter :: PROC_DECAY = 1
  integer, parameter :: PROC_SCATTER = 2


contains

  module subroutine split_constraint_check_before_split (c, table, pl, k, passed)
    class(split_constraint_t), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: k
    logical, intent(out) :: passed
    passed = .true.
  end subroutine split_constraint_check_before_split

  module subroutine split_constraint_check_before_insert (c, table, pa, pl, passed)
    class(split_constraint_t), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_array_t), intent(in) :: pa
    type(pdg_list_t), intent(inout) :: pl
    logical, intent(out) :: passed
    passed = .true.
  end subroutine split_constraint_check_before_insert

  module subroutine split_constraint_check_before_record (c, table, pl, n_loop, passed)
    class(split_constraint_t), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop
    logical, intent(out) :: passed
    passed = .true.
  end subroutine split_constraint_check_before_record

  module subroutine split_constraints_set (constraints, i, c)
    class(split_constraints_t), intent(inout) :: constraints
    integer, intent(in) :: i
    class(split_constraint_t), intent(in) :: c
    allocate (constraints%cc(i)%c, source = c)
  end subroutine split_constraints_set

  module subroutine split_constraints_check_before_split &
       (constraints, table, pl, k, passed)
    class(split_constraints_t), intent(in) :: constraints
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: k
    logical, intent(out) :: passed
    integer :: i
    passed = .true.
    do i = 1, size (constraints%cc)
       call constraints%cc(i)%c%check_before_split (table, pl, k, passed)
       if (.not. passed)  return
    end do
  end subroutine split_constraints_check_before_split

  module subroutine split_constraints_check_before_insert &
       (constraints, table, pa, pl, passed)
    class(split_constraints_t), intent(in) :: constraints
    class(ps_table_t), intent(in) :: table
    type(pdg_array_t), intent(in) :: pa
    type(pdg_list_t), intent(inout) :: pl
    logical, intent(out) :: passed
    integer :: i
    passed = .true.
    do i = 1, size (constraints%cc)
       call constraints%cc(i)%c%check_before_insert (table, pa, pl, passed)
       if (.not. passed)  return
    end do
  end subroutine split_constraints_check_before_insert

  module subroutine split_constraints_check_before_record &
       (constraints, table, pl, n_loop, passed)
    class(split_constraints_t), intent(in) :: constraints
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop
    logical, intent(out) :: passed
    integer :: i
    passed = .true.
    do i = 1, size (constraints%cc)
       call constraints%cc(i)%c%check_before_record (table, pl, n_loop, passed)
       if (.not. passed)  return
    end do
  end subroutine split_constraints_check_before_record

  module function constrain_n_tot (n_max) result (c)
    integer, intent(in) :: n_max
    type(constraint_n_tot) :: c
    c%n_max = n_max
  end function constrain_n_tot

  module subroutine constraint_n_tot_check_before_split (c, table, pl, k, passed)
    class(constraint_n_tot), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: k
    logical, intent(out) :: passed
    passed = pl%get_size () < c%n_max
  end subroutine constraint_n_tot_check_before_split

  module subroutine constraint_n_tot_check_before_record (c, table, pl, n_loop, passed)
    class(constraint_n_tot), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop
    logical, intent(out) :: passed
    passed = pl%get_size () + n_loop <= c%n_max
  end subroutine constraint_n_tot_check_before_record

  module function constrain_n_loop (n_loop_max) result (c)
    integer, intent(in) :: n_loop_max
    type(constraint_n_loop) :: c
    c%n_loop_max = n_loop_max
  end function constrain_n_loop

  module subroutine constraint_n_loop_check_before_record &
       (c, table, pl, n_loop, passed)
    class(constraint_n_loop), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop
    logical, intent(out) :: passed
    passed = n_loop <= c%n_loop_max
  end subroutine constraint_n_loop_check_before_record

  module function constrain_splittings &
       (pl_match, pl_excluded_gauge_splittings) result (c)
    type(pdg_list_t), intent(in) :: pl_match
    type(pdg_list_t), intent(in) :: pl_excluded_gauge_splittings
    type(constraint_splittings) :: c
    c%pl_match = pl_match
    c%pl_excluded_gauge_splittings = pl_excluded_gauge_splittings
  end function constrain_splittings

  module subroutine constraint_splittings_check_before_insert &
       (c, table, pa, pl, passed)
    class(constraint_splittings), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_array_t), intent(in) :: pa
    type(pdg_list_t), intent(inout) :: pl
    logical, intent(out) :: passed
    logical :: has_massless_vector
    integer :: i
    has_massless_vector = .false.
    do i = 1, pa%get_length ()
       if (is_massless_vector(pa%get(i))) then
          has_massless_vector = .true.
          exit
       end if
    end do
    passed = .false.
    if (has_massless_vector .and. count (is_fermion(pl%a%get ())) == 2) then
       do i = 1, c%pl_excluded_gauge_splittings%get_size ()
          if (pl .match. c%pl_excluded_gauge_splittings%a(i)) return
       end do
       call pl%match_replace (c%pl_match, passed)
       passed = .true.
    else
       call pl%match_replace (c%pl_match, passed)
    end if
  end subroutine constraint_splittings_check_before_insert

  module function constrain_insert (pl_match) result (c)
    type(pdg_list_t), intent(in) :: pl_match
    type(constraint_insert) :: c
    c%pl_match = pl_match
  end function constrain_insert

  module subroutine constraint_insert_check_before_insert (c, table, pa, pl, passed)
    class(constraint_insert), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_array_t), intent(in) :: pa
    type(pdg_list_t), intent(inout) :: pl
    logical, intent(out) :: passed
    call pl%match_replace (c%pl_match, passed)
  end subroutine constraint_insert_check_before_insert

  module function constrain_require (pl) result (c)
    type(pdg_list_t), intent(in) :: pl
    type(constraint_require) :: c
    c%pl = pl
  end function constrain_require

  module subroutine constraint_require_check_before_record &
       (c, table, pl, n_loop, passed)
    class(constraint_require), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop
    logical, intent(out) :: passed
    logical, dimension(:), allocatable :: mask
    integer :: i, k, n_in
    select type (table)
    type is (if_table_t)
       if (table%proc_type > 0) then
          select case (table%proc_type)
          case (PROC_DECAY)
             n_in = 1
          case (PROC_SCATTER)
             n_in = 2
          end select
       else
          call msg_fatal ("Neither a decay nor a scattering process")
       end if
    class default
       n_in = 0
    end select
    allocate (mask (c%pl%get_size ()), source = .true.)
    do i = n_in + 1, pl%get_size ()
       k = c%pl%find_match (pl%get (i), mask)
       if (k /= 0)  mask(k) = .false.
    end do
    passed = .not. any (mask)
  end subroutine constraint_require_check_before_record

  module function constrain_radiation () result (c)
    type(constraint_radiation) :: c
  end function constrain_radiation

  module subroutine constraint_radiation_check_before_insert &
       (c, table, pa, pl, passed)
    class(constraint_radiation), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_array_t), intent(in) :: pa
    type(pdg_list_t), intent(inout) :: pl
    logical, intent(out) :: passed
    passed = .not. (pl .match. pa)
  end subroutine constraint_radiation_check_before_insert

  module function constrain_mass_sum (mass_limit, margin) result (c)
    real(default), intent(in) :: mass_limit
    real(default), intent(in), optional :: margin
    type(constraint_mass_sum) :: c
    c%mass_limit = mass_limit
    if (present (margin)) then
       c%strictly_less = .true.
       c%margin = margin
    end if
  end function constrain_mass_sum

  module subroutine constraint_mass_sum_check_before_record &
       (c, table, pl, n_loop, passed)
    class(constraint_mass_sum), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop
    logical, intent(out) :: passed
    real(default) :: limit
    if (c%strictly_less) then
       limit = c%mass_limit - c%margin
       select type (table)
       type is (if_table_t)
          passed = mass_sum (pl, 1, 2, table%model) < limit  &
               .and. mass_sum (pl, 3, pl%get_size (), table%model) < limit
       class default
          passed = mass_sum (pl, 1, pl%get_size (), table%model) < limit
       end select
    else
       limit = c%mass_limit
       select type (table)
       type is (if_table_t)
          passed = mass_sum (pl, 1, 2, table%model) <= limit &
               .and. mass_sum (pl, 3, pl%get_size (), table%model) <= limit
       class default
          passed = mass_sum (pl, 1, pl%get_size (), table%model) <= limit
       end select
    end if
  end subroutine constraint_mass_sum_check_before_record

  module function constrain_in_state (pl) result (c)
    type(pdg_list_t), intent(in) :: pl
    type(constraint_in_state) :: c
    c%pl = pl
  end function constrain_in_state

  module subroutine constraint_in_state_check_before_record &
       (c, table, pl, n_loop, passed)
    class(constraint_in_state), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop
    logical, intent(out) :: passed
    integer :: i
    select type (table)
    type is (if_table_t)
       passed = .false.
       do i = 1, 2
          if (.not. (c%pl .match. pl%get (i)))  return
       end do
    end select
    passed = .true.
  end subroutine constraint_in_state_check_before_record

  module function constrain_photon_induced_processes (n_in) result (c)
    integer, intent(in) :: n_in
    type(constraint_photon_induced_processes) :: c
    c%n_in = n_in
  end function constrain_photon_induced_processes

  module subroutine constraint_photon_induced_processes_check_before_record &
       (c, table, pl, n_loop, passed)
    class(constraint_photon_induced_processes), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop
    logical, intent(out) :: passed
    integer :: i
    select type (table)
    type is (if_table_t)
       passed = .false.
       do i = 1, c%n_in
          if (pl%a(i)%get () == 22)  return
       end do
    end select
    passed = .true.
  end subroutine constraint_photon_induced_processes_check_before_record

  module function constrain_couplings (qcd, qed, n_nlo_correction_types) result (c)
    type(constraint_coupling_t) :: c
    logical, intent(in) :: qcd, qed
    integer, intent(in) :: n_nlo_correction_types
    c%qcd = qcd; c%qed = qed
    c%n_nlo_correction_types = n_nlo_correction_types
  end function constrain_couplings

  module subroutine constraint_coupling_check_before_insert (c, table, pa, pl, passed)
    class(constraint_coupling_t), intent(in) :: c
    class(ps_table_t), intent(in) :: table
    type(pdg_array_t), intent(in) :: pa
    type(pdg_list_t), intent(inout) :: pl
    logical, intent(out) :: passed
    type(pdg_list_t) :: pl_vertex
    type(pdg_array_t) :: pdg_gluon, pdg_photon, pdg_W_Z, pdg_gauge_bosons
    integer :: i, j
    pdg_gluon = GLUON; pdg_photon = PHOTON
    pdg_W_Z = [W_BOSON,-W_BOSON, Z_BOSON]
    if (c%qcd) pdg_gauge_bosons = pdg_gauge_bosons // pdg_gluon
    if (c%qed) pdg_gauge_bosons = pdg_gauge_bosons // pdg_photon
    if (c%ew) pdg_gauge_bosons = pdg_gauge_bosons // pdg_W_Z
    do j = 1, pa%get_length ()
       call pl_vertex%init (pl%get_size () + 1)
       call pl_vertex%set (1, pa%get(j))
       do i = 1, pl%get_size ()
          call pl_vertex%set (i + 1, pl%get(i))
       end do
       if (pl_vertex%get_size () > 3) then
          passed = .false.
          cycle
       end if
       if (is_massless_vector(pa%get(j))) then
          if (.not. table%model%check_vertex &
               (pl_vertex%a(1)%get (), pl_vertex%a(2)%get (), pl_vertex%a(3)%get ())) then
             passed = .false.
             cycle
          end if
       else if (.not. table%model%check_vertex &
            (- pl_vertex%a(1)%get (), pl_vertex%a(2)%get (), pl_vertex%a(3)%get ())) then
          passed = .false.
          cycle
       end if
       if (.not. (pl_vertex .match. pdg_gauge_bosons)) then
          passed = .false.
          cycle
       end if
       passed = .true.
       exit
    end do
  end subroutine constraint_coupling_check_before_insert

  module subroutine ps_table_final (object)
    class(ps_table_t), intent(inout) :: object
    type(ps_entry_t), pointer :: current
    do while (associated (object%first))
       current => object%first
       object%first => current%next
       deallocate (current)
    end do
    nullify (object%last)
  end subroutine ps_table_final

  module subroutine ps_table_base_write (object, unit, n_in)
    class(ps_table_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: n_in
    integer, dimension(:), allocatable :: pdg
    type(ps_entry_t), pointer :: entry
    type(field_data_t), pointer :: prt
    integer :: u, i, j, n0
    u = given_output_unit (unit)
    entry => object%first
    do while (associated (entry))
       write (u, "(2x)", advance = "no")
       if (present (n_in)) then
          do i = 1, n_in
             write (u, "(1x)", advance = "no")
             pdg = entry%get (i)
             do j = 1, size (pdg)
                prt => object%model%get_field_ptr (pdg(j))
                if (j > 1)  write (u, "(':')", advance = "no")
                write (u, "(A)", advance = "no") &
                     char (prt%get_name (pdg(j) >= 0))
             end do
          end do
          write (u, "(1x,A)", advance = "no")  "=>"
          n0 = n_in + 1
       else
          n0 = 1
       end if
       do i = n0, entry%get_size ()
          write (u, "(1x)", advance = "no")
          pdg = entry%get (i)
          do j = 1, size (pdg)
             prt => object%model%get_field_ptr (pdg(j))
             if (j > 1)  write (u, "(':')", advance = "no")
             write (u, "(A)", advance = "no") &
                  char (prt%get_name (pdg(j) < 0))
          end do
       end do
       if (object%loops) then
          write (u, "(2x,'[',I0,',',I0,']')")  entry%n_loop, entry%n_rad
       else
          write (u, "(A)")
       end if
       entry => entry%next
    end do
  end subroutine ps_table_base_write

  module subroutine ds_table_write (object, unit)
    class(ds_table_t), intent(in) :: object
    integer, intent(in), optional :: unit
    type(field_data_t), pointer :: prt
    integer :: u
    u = given_output_unit (unit)
    prt => object%model%get_field_ptr (object%pdg_in)
    write (u, "(1x,A,1x,A)")  "Decays for particle:", &
         char (prt%get_name (object%pdg_in < 0))
    call object%base_write (u)
  end subroutine ds_table_write

  module subroutine fs_table_write (object, unit)
    class(fs_table_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Table of final states:"
    call object%base_write (u)
  end subroutine fs_table_write

  module subroutine if_table_write (object, unit)
    class(if_table_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Table of in/out states:"
    select case (object%proc_type)
    case (PROC_DECAY)
       call object%base_write (u, n_in = 1)
    case (PROC_SCATTER)
       call object%base_write (u, n_in = 2)
     end select
  end subroutine if_table_write

  module subroutine ps_table_get_particle_string (object, index, prt_in, prt_out)
    class(ps_table_t), intent(in) :: object
    integer, intent(in) :: index
    type(string_t), intent(out), dimension(:), allocatable :: prt_in, prt_out
    integer :: n_in
    type(field_data_t), pointer :: prt
    type(ps_entry_t), pointer :: entry
    integer, dimension(:), allocatable :: pdg
    integer :: n0
    integer :: i, j
    entry => object%first
    i = 1
    do while (i < index)
      if (associated (entry%next)) then
         entry => entry%next
         i = i + 1
      else
         call msg_fatal ("ps_table: entry with requested index does not exist!")
      end if
    end do

    if (object%proc_type > 0) then
       select case (object%proc_type)
       case (PROC_DECAY)
          n_in = 1
       case (PROC_SCATTER)
          n_in = 2
       end select
    else
       call msg_fatal ("Neither decay nor scattering process")
    end if

    n0 = n_in + 1
    allocate (prt_in (n_in), prt_out (entry%get_size () - n_in))
    do i = 1, n_in
       prt_in(i) = ""
       pdg = entry%get(i)
       do j = 1, size (pdg)
          prt => object%model%get_field_ptr (pdg(j))
          prt_in(i) = prt_in(i) // prt%get_name (pdg(j) >= 0)
          if (j /= size (pdg))  prt_in(i) = prt_in(i) // ":"
       end do
    end do
    do i = n0, entry%get_size ()
       prt_out(i-n_in) = ""
       pdg = entry%get(i)
       do j = 1, size (pdg)
          prt => object%model%get_field_ptr (pdg(j))
          prt_out(i-n_in) = prt_out(i-n_in) // prt%get_name (pdg(j) < 0)
          if (j /= size (pdg))  prt_out(i-n_in) = prt_out(i-n_in) // ":"
       end do
    end do
  end subroutine ps_table_get_particle_string

  module subroutine ps_table_init &
       (table, model, pl, constraints, n_in, do_not_check_regular)
    class(ps_table_t), intent(out) :: table
    class(model_data_t), intent(in), target :: model
    type(pdg_list_t), dimension(:), intent(in) :: pl
    type(split_constraints_t), intent(in) :: constraints
    integer, intent(in), optional :: n_in
    logical, intent(in), optional :: do_not_check_regular
    logical :: passed
    integer :: i
    table%model => model

    if (present (n_in)) then
       select case (n_in)
       case (1)
          table%proc_type = PROC_DECAY
       case (2)
          table%proc_type = PROC_SCATTER
       case default
          table%proc_type = PROC_UNDEFINED
       end select
    else
       table%proc_type = PROC_UNDEFINED
    end if

    do i = 1, size (pl)
       call table%record (pl(i), 0, 0, constraints, &
            do_not_check_regular, passed)
       if (.not. passed) then
          call msg_fatal ("ps_table: Registering process components failed")
       end if
    end do
  end subroutine ps_table_init

  module subroutine if_table_init (table, model, pl_in, pl_out, constraints)
    class(if_table_t), intent(out) :: table
    class(model_data_t), intent(in), target :: model
    type(pdg_list_t), dimension(:), intent(in) :: pl_in, pl_out
    type(split_constraints_t), intent(in) :: constraints
    integer :: i, j, k, p, n_in, n_out
    type(pdg_array_t), dimension(:), allocatable :: pa_in
    type(pdg_list_t), dimension(:), allocatable :: pl
    allocate (pl (size (pl_in) * size (pl_out)))
    k = 0
    do i = 1, size (pl_in)
       n_in = pl_in(i)%get_size ()
       allocate (pa_in (n_in))
       do p = 1, n_in
          pa_in(p) = pl_in(i)%get (p)
       end do
       do j = 1, size (pl_out)
          n_out = pl_out(j)%get_size ()
          k = k + 1
          call pl(k)%init (n_in + n_out)
          do p = 1, n_in
             call pl(k)%set (p, invert_pdg_array (pa_in(p), model))
          end do
          do p = 1, n_out
             call pl(k)%set (n_in + p, pl_out(j)%get (p))
          end do
       end do
       deallocate (pa_in)
    end do
    n_in = size (pl_in(1)%a)
    call table%init (model, pl, constraints, n_in, do_not_check_regular = .true.)
  end subroutine if_table_init

  module subroutine ps_table_enable_loops (table)
    class(ps_table_t), intent(inout) :: table
    table%loops = .true.
  end subroutine ps_table_enable_loops

  module subroutine ds_table_make (table, model, pdg_in, constraints)
    class(ds_table_t), intent(out) :: table
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: pdg_in
    type(split_constraints_t), intent(in) :: constraints
    type(pdg_list_t) :: pl_in
    type(pdg_list_t), dimension(0) :: pl
    call table%init (model, pl, constraints)
    table%pdg_in = pdg_in
    call pl_in%init (1)
    call pl_in%set (1, [pdg_in])
    call table%split (pl_in, 0, constraints)
  end subroutine ds_table_make

  module subroutine fs_table_radiate (table, constraints, do_not_check_regular)
    class(fs_table_t), intent(inout) :: table
    type(split_constraints_t) :: constraints
    logical, intent(in), optional :: do_not_check_regular
    type(ps_entry_t), pointer :: current
    current => table%first
    do while (associated (current))
       call table%split (current, 0, constraints, record = .true., &
            do_not_check_regular = do_not_check_regular)
       current => current%next
    end do
  end subroutine fs_table_radiate

  recursive module subroutine ps_table_split (table, pl, n_rad, constraints, &
        record, do_not_check_regular)
    class(ps_table_t), intent(inout) :: table
    class(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_rad
    type(split_constraints_t), intent(in) :: constraints
    logical, intent(in), optional :: record, do_not_check_regular
    integer :: n_loop, i
    logical :: passed, save_pdg_index
    type(vertex_iterator_t) :: vit
    integer, dimension(:), allocatable :: pdg1
    integer, dimension(:), allocatable :: pdg2
    if (present (record)) then
       if (record) then
          n_loop = 0
          INCR_LOOPS: do
             call table%record_sorted (pl, n_loop, n_rad, constraints, &
                  do_not_check_regular, passed)
             if (.not. passed)  exit INCR_LOOPS
             if (.not. table%loops)  exit INCR_LOOPS
             n_loop = n_loop + 1
          end do INCR_LOOPS
       end if
    end if
    select type (table)
    type is (if_table_t)
       save_pdg_index = .true.
    class default
       save_pdg_index = .false.
    end select
    do i = 1, pl%get_size ()
       call constraints%check_before_split (table, pl, i, passed)
       if (passed) then
          pdg1 = pl%get (i)
          call vit%init (table%model, pdg1, save_pdg_index)
          SCAN_VERTICES: do
             call vit%get_next_match (pdg2)
             if (allocated (pdg2)) then
                call table%insert (pl, n_rad, i, pdg2, constraints, &
                     do_not_check_regular = do_not_check_regular)
             else
                exit SCAN_VERTICES
             end if
          end do SCAN_VERTICES
       end if
    end do
  end subroutine ps_table_split

  recursive module subroutine ps_table_insert &
       (table, pl, n_rad, i, pdg, constraints, n_in, do_not_check_regular)
    class(ps_table_t), intent(inout) :: table
    class(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_rad, i
    integer, dimension(:), intent(in) :: pdg
    type(split_constraints_t), intent(in) :: constraints
    integer, intent(in), optional :: n_in
    logical, intent(in), optional :: do_not_check_regular
    type(pdg_list_t) :: pl_insert
    logical :: passed
    integer :: k, s
    s = size (pdg)
    call pl_insert%init (s)
    do k = 1, s
       call pl_insert%set (k, pdg(k))
    end do
    call constraints%check_before_insert (table, pl%get (i), pl_insert, passed)
    if (passed) then
       if (.not. is_colored_isr ()) return
       call table%split (pl%replace (i, pl_insert, n_in), n_rad + s - 1, &
            constraints, record = .true., do_not_check_regular = .true.)
    end if
  contains
    logical function is_colored_isr () result (ok)
      type(pdg_list_t) :: pl_replaced
      ok = .true.
      if (present (n_in)) then
         if (i <= n_in) then
            ok = pl_insert%contains_colored_particles ()
            if (.not. ok) then
               pl_replaced = pl%replace (i, pl_insert, n_in)
               associate (size_replaced => pl_replaced%get_pdg_sizes (), &
                    size => pl%get_pdg_sizes ())
                  ok = all (size_replaced(:n_in) == size(:n_in))
               end associate
            end if
         end if
      end if
    end function is_colored_isr
 end subroutine ps_table_insert

  recursive module subroutine if_table_insert  &
       (table, pl, n_rad, i, pdg, constraints, n_in, do_not_check_regular)
    class(if_table_t), intent(inout) :: table
    class(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_rad, i
    integer, dimension(:), intent(in) :: pdg
    type(split_constraints_t), intent(in) :: constraints
    integer, intent(in), optional :: n_in
    logical, intent(in), optional :: do_not_check_regular
    integer, dimension(:), allocatable :: pdg_work
    integer :: p
    if (i > 2) then
       call ps_table_insert (table, pl, n_rad, i, pdg, constraints, &
            do_not_check_regular = do_not_check_regular)
    else
       allocate (pdg_work (size (pdg)))
       do p = 1, size (pdg)
          pdg_work(1) = pdg(p)
          pdg_work(2:p) = pdg(1:p-1)
          pdg_work(p+1:) = pdg(p+1:)
          select case (table%proc_type)
          case (PROC_DECAY)
             call ps_table_insert (table, &
                  pl, n_rad, i, pdg_work, constraints, n_in = 1, &
                  do_not_check_regular = do_not_check_regular)
          case (PROC_SCATTER)
             call ps_table_insert (table, &
                  pl, n_rad, i, pdg_work, constraints, n_in = 2, &
                  do_not_check_regular = do_not_check_regular)
          end select
       end do
    end if
  end subroutine if_table_insert

  module subroutine ps_table_record_sorted &
       (table, pl, n_loop, n_rad, constraints, do_not_check_regular, passed)
    class(ps_table_t), intent(inout) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop, n_rad
    type(split_constraints_t), intent(in) :: constraints
    logical, intent(in), optional :: do_not_check_regular
    logical, intent(out) :: passed
    call table%record (pl%sort_abs (), n_loop, n_rad, constraints, &
         do_not_check_regular, passed)
  end subroutine ps_table_record_sorted

  module subroutine if_table_record_sorted &
       (table, pl, n_loop, n_rad, constraints, do_not_check_regular, passed)
    class(if_table_t), intent(inout) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop, n_rad
    type(split_constraints_t), intent(in) :: constraints
    logical, intent(in), optional :: do_not_check_regular
    logical, intent(out) :: passed
    call table%record (pl%sort_abs (2), n_loop, n_rad, constraints, &
         do_not_check_regular, passed)
  end subroutine if_table_record_sorted

  module subroutine ps_table_record (table, pl, n_loop, n_rad, constraints, &
       do_not_check_regular, passed)
    class(ps_table_t), intent(inout) :: table
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n_loop, n_rad
    type(split_constraints_t), intent(in) :: constraints
    logical, intent(in), optional :: do_not_check_regular
    logical, intent(out) :: passed
    type(ps_entry_t), pointer :: current
    logical :: needs_check
    passed = .false.
    needs_check = .true.
    if (present (do_not_check_regular))  needs_check = .not. do_not_check_regular
    if (needs_check .and. .not. pl%is_regular ()) then
       call msg_warning ("Record ps_table entry: Irregular pdg-list encountered!")
       return
    end if
    call constraints%check_before_record (table, pl, n_loop, passed)
    if (.not. passed)  then
       return
    end if
    current => table%first
    do while (associated (current))
       if (pl == current) then
          if (n_loop == current%n_loop)  return
       else if (pl < current) then
          call record_insert ()
          return
       end if
       current => current%next
    end do
    call record_insert ()
  contains
    subroutine record_insert ()
      type(ps_entry_t), pointer :: entry
      allocate (entry)
      entry%pdg_list_t = pl
      entry%n_loop = n_loop
      entry%n_rad = n_rad
      if (associated (current)) then
         if (associated (current%previous)) then
            current%previous%next => entry
            entry%previous => current%previous
         else
            table%first => entry
         end if
         entry%next => current
         current%previous => entry
      else
         if (associated (table%last)) then
            table%last%next => entry
            entry%previous => table%last
         else
            table%first => entry
         end if
         table%last => entry
      end if
    end subroutine record_insert
  end subroutine ps_table_record

  function mass_sum (pl, n1, n2, model) result (m)
    type(pdg_list_t), intent(in) :: pl
    integer, intent(in) :: n1, n2
    class(model_data_t), intent(in), target :: model
    integer, dimension(:), allocatable :: pdg
    real(default) :: m
    type(field_data_t), pointer :: prt
    integer :: i
    m = 0
    do i = n1, n2
       pdg = pl%get (i)
       prt => model%get_field_ptr (pdg(1))
       m = m + prt%get_mass ()
    end do
  end function mass_sum

  function invert_pdg_array (pa, model) result (pa_inv)
    type(pdg_array_t), intent(in) :: pa
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t) :: pa_inv
    type(field_data_t), pointer :: prt
    integer :: i, pdg
    pa_inv = pa
    do i = 1, pa_inv%get_length ()
       pdg = pa_inv%get (i)
       prt => model%get_field_ptr (pdg)
       if (prt%has_antiparticle ())  call pa_inv%set (i, -pdg)
    end do
  end function invert_pdg_array

  module function ps_table_get_length (ps_table) result (n)
    class(ps_table_t), intent(in) :: ps_table
    integer :: n
    type(ps_entry_t), pointer :: entry
    n = 0
    entry => ps_table%first
    do while (associated (entry))
       n = n + 1
       entry => entry%next
    end do
  end function ps_table_get_length

  module subroutine ps_table_get_emitters (table, constraints, emitters)
    class(ps_table_t), intent(in) :: table
    type(split_constraints_t), intent(in) :: constraints
    integer, dimension(:), allocatable, intent(out) :: emitters
    class(pdg_list_t), pointer :: pl
    integer :: i
    logical :: passed
    type(vertex_iterator_t) :: vit
    integer, dimension(:), allocatable :: pdg1, pdg2
    integer :: n_emitters
    integer, dimension(:), allocatable :: emitters_tmp
    integer, parameter :: buf0 = 6
    n_emitters = 0
    pl => table%first
    allocate (emitters_tmp (buf0))
    do i = 1, pl%get_size ()
       call constraints%check_before_split (table, pl, i, passed)
       if (passed) then
          pdg1 = pl%get(i)
          call vit%init (table%model, pdg1, .false.)
          do
             call vit%get_next_match(pdg2)
             if (allocated (pdg2)) then
                if (n_emitters + 1 > size (emitters_tmp)) &
                     call extend_integer_array (emitters_tmp, 10)
                emitters_tmp (n_emitters + 1) = pdg1(1)
                n_emitters = n_emitters + 1
             else
                exit
             end if
          end do
       end if
    end do
    allocate (emitters (n_emitters))
    emitters = emitters_tmp (1:n_emitters)
    deallocate (emitters_tmp)
  end subroutine ps_table_get_emitters

  module subroutine ps_table_get_pdg_out (ps_table, i, pa_out, n_loop, n_rad)
    class(ps_table_t), intent(in) :: ps_table
    integer, intent(in) :: i
    type(pdg_array_t), dimension(:), allocatable, intent(out) :: pa_out
    integer, intent(out), optional :: n_loop, n_rad
    type(ps_entry_t), pointer :: entry
    integer :: n, j
    n = 0
    entry => ps_table%first
    FIND_ENTRY: do while (associated (entry))
       n = n + 1
       if (n == i) then
          allocate (pa_out (entry%get_size ()))
          do j = 1, entry%get_size ()
             pa_out(j) = entry%get (j)
             if (present (n_loop))  n_loop = entry%n_loop
             if (present (n_rad))  n_rad = entry%n_rad
          end do
          exit FIND_ENTRY
       end if
       entry => entry%next
    end do FIND_ENTRY
  end subroutine ps_table_get_pdg_out


end submodule auto_components_s

