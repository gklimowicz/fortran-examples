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

submodule (interactions) interactions_s

  use io_units
  use diagnostics
  use sorting

  implicit none

contains

  module subroutine qn_index_map_init_trivial (self, int)
    class(qn_index_map_t), intent(out) :: self
    class(interaction_t), intent(in) :: int
    integer :: qn
    self%n_flv = int%get_n_matrix_elements ()
    self%n_hel = 1
    self%n_sub = 0
    allocate (self%index(self%n_flv, self%n_hel, 0:self%n_sub), source = 0)
    do qn = 1, self%n_flv
       self%index(qn, 1, 0) = qn
    end do
  end subroutine qn_index_map_init_trivial

  module subroutine qn_index_map_init_involved (self, int, qn_flv, n_sub, qn_hel)
    class(qn_index_map_t), intent(out) :: self
    type(interaction_t), intent(in) :: int
    type(quantum_numbers_t), dimension(:, :), intent(in) :: qn_flv
    integer, intent(in) :: n_sub
    type(quantum_numbers_t), dimension(:, :), intent(in), optional :: qn_hel
    type(quantum_numbers_t), dimension(:), allocatable :: qn, qn_int
    integer :: i, i_flv, i_hel, i_sub
    self%qn_flv = qn_flv
    self%n_flv = size (qn_flv, dim=2)
    self%n_sub = n_sub
    if (present (qn_hel)) then
       if (size (qn_flv, dim=1) /= size (qn_hel, dim=1)) then
          call msg_bug ("[qn_index_map_init] number of particles does not match.")
       end if
       self%qn_hel = qn_hel
       self%n_hel = size (qn_hel, dim=2)
    else
       self%n_hel = 1
    end if
    allocate (self%index (self%n_flv, self%n_hel, 0:self%n_sub), source=0)
    associate (n_me => int%get_n_matrix_elements ())
       do i = 1, n_me
          qn_int = int%get_quantum_numbers (i, by_me_index = .true.)
          qn = pack (qn_int, qn_int%are_hard_process ())
          i_flv = find_flv_index (self, qn)
          i_hel = 1; if (allocated (self%qn_hel)) &
               i_hel = find_hel_index (self, qn)
          i_sub = find_sub_index (self, qn)
          self%index(i_flv, i_hel, i_sub) = i
       end do
    end associate
  contains
    integer function find_flv_index (self, qn) result (i_flv)
      type(qn_index_map_t), intent(in) :: self
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      integer :: j
      i_flv = 0
      do j = 1, self%n_flv
         if (.not. all (qn .fmatch. self%qn_flv(:, j))) cycle
         i_flv = j
         exit
      end do
      if (i_flv < 1) then
         call msg_message ("QN:")
         call quantum_numbers_write (qn)
         call msg_message ("")
         call msg_message ("QN_FLV:")
         do j = 1, self%n_flv
            call quantum_numbers_write (self%qn_flv(:, j))
            call msg_message ("")
         end do
         call msg_bug ("[find_flv_index] could not find flv in qn_flv.")
      end if
    end function find_flv_index

    integer function find_hel_index (self, qn) result (i_hel)
      type(qn_index_map_t), intent(in) :: self
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      integer :: j
      i_hel = 0
      do j = 1, self%n_hel
         if (.not. all (qn .hmatch. self%qn_hel(:, j))) cycle
         i_hel = j
         exit
      end do
      if (i_hel < 1) then
         call msg_message ("QN:")
         call quantum_numbers_write (qn)
         call msg_message ("")
         call msg_message ("QN_HEL:")
         do j = 1, self%n_hel
            call quantum_numbers_write (self%qn_hel(:, j))
            call msg_message ("")
         end do
         call msg_bug ("[find_hel_index] could not find hel in qn_hel.")
      end if
    end function find_hel_index

    integer function find_sub_index (self, qn) result (i_sub)
      type(qn_index_map_t), intent(in) :: self
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      integer :: s
      i_sub = -1
      do s = 0, self%n_sub
         if ((all (pack(qn%get_sub (), qn%get_sub () > 0) == s)) &
              .or. (all (qn%get_sub () == 0) .and. s == 0)) then
            i_sub = s
            exit
         end if
      end do
      if (i_sub < 0) then
         call msg_message ("QN:")
         call quantum_numbers_write (qn)
         call msg_bug ("[find_sub_index] could not find sub in qn.")
      end if
    end function find_sub_index
  end subroutine qn_index_map_init_involved

  module subroutine qn_index_map_init_sf (self, int, qn_flv, n_flv_born, n_flv_real)
    class(qn_index_map_t), intent(out) :: self
    type(interaction_t), intent(in) :: int
    integer, intent(in) :: n_flv_born, n_flv_real
    type(quantum_numbers_t), dimension(:,:), intent(in) :: qn_flv
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn_int
    type(quantum_numbers_t), dimension(:), allocatable :: qn_int_tmp
    integer :: i, i_sub, n_flv, n_hard
    n_flv = int%get_n_matrix_elements ()
    qn_int_tmp = int%get_quantum_numbers (1, by_me_index = .true.)
    n_hard = count (qn_int_tmp%are_hard_process ())
    allocate (qn_int(n_hard, n_flv))
    do i = 1, n_flv
       qn_int_tmp = int%get_quantum_numbers (i, by_me_index = .true.)
       qn_int(:, i) = pack (qn_int_tmp, qn_int_tmp%are_hard_process ())
    end do
    call self%init (int, qn_int, int%get_n_sub ())
    allocate (self%sf_index_born(n_flv_born, 0:self%n_sub))
    allocate (self%sf_index_real(n_flv_real, 0:self%n_sub))
    do i_sub = 0, self%n_sub
       do i = 1, n_flv_born
          self%sf_index_born(i, i_sub) = self%get_index_by_qn (qn_flv(:,i), i_sub)
       end do
       do i = 1, n_flv_real
          self%sf_index_real(i, i_sub) = &
               self%get_index_by_qn (qn_flv(:,n_flv_born + i), i_sub)
       end do
    end do
    deallocate (self%index)
  end subroutine qn_index_map_init_sf

  module subroutine qn_index_map_write (self, unit)
    class(qn_index_map_t), intent(in) :: self
    integer, intent(in), optional :: unit
    integer :: u, i_flv, i_hel, i_sub
    u = given_output_unit (unit); if (u < 0) return
    write (u, *) "flip_hel: ", self%flip_hel
    do i_flv = 1, self%n_flv
       if (allocated (self%qn_flv)) &
            call quantum_numbers_write (self%qn_flv(:, i_flv))
       write (u, *)
       do i_hel = 1, self%n_hel
          if (allocated (self%qn_hel)) then
             call quantum_numbers_write (self%qn_hel(:, i_hel))
             write (u, *)
          end if
          do i_sub = 0, self%n_sub
             write (u, *) &
                  "(", i_flv, ",", i_hel, ",", i_sub, ") => ", self%index(i_flv, i_hel, i_sub)
          end do
       end do
    end do
  end subroutine qn_index_map_write

  module subroutine qn_index_map_set_helicity_flip (self, yorn)
    class(qn_index_map_t), intent(inout) :: self
    logical, intent(in) :: yorn
    integer :: i, i_flv, i_hel, i_hel_new
    type(quantum_numbers_t), dimension(:, :), allocatable :: qn_hel_flip
    integer, dimension(:, :, :), allocatable :: index
    if (.not. allocated (self%qn_hel)) then
       call msg_bug ("[qn_index_map_set_helicity_flip] &
            &cannot flip not-given helicity.")
    end if
    allocate (index (self%n_flv, self%n_hel, 0:self%n_sub), &
         source=self%index)
    self%flip_hel = yorn
    if (self%flip_hel) then
       do i_flv = 1, self%n_flv
          qn_hel_flip = self%qn_hel
          do i_hel = 1, self%n_hel
             do i = 1, size (self%qn_flv, dim=1)
                if (is_anti_particle (self%qn_flv(i, i_flv))) then
                   call qn_hel_flip(i, i_hel)%flip_helicity ()
                end if
             end do
          end do
          do i_hel = 1, self%n_hel
             i_hel_new = find_hel_index (qn_hel_flip, self%qn_hel(:, i_hel))
             self%index(i_flv, i_hel_new, :) = index(i_flv, i_hel, :)
          end do
       end do
    end if
  contains
    logical function is_anti_particle (qn) result (yorn)
      type(quantum_numbers_t), intent(in) :: qn
      type(flavor_t) :: flv
      flv = qn%get_flavor ()
      yorn = flv%get_pdg () < 0
    end function is_anti_particle

    integer function find_hel_index (qn_sort, qn) result (i_hel)
      type(quantum_numbers_t), dimension(:, :), intent(in) :: qn_sort
      type(quantum_numbers_t), dimension(:), intent(in) :: qn
      integer :: j
      do j = 1, size(qn_sort, dim=2)
         if (.not. all (qn .hmatch. qn_sort(:, j))) cycle
         i_hel = j
         exit
      end do
    end function find_hel_index
  end subroutine qn_index_map_set_helicity_flip

  module function qn_index_map_get_index (self, i_flv, i_hel, i_sub) result (index)
    class(qn_index_map_t), intent(in) :: self
    integer :: index
    integer, intent(in) :: i_flv
    integer, intent(in), optional :: i_hel
    integer, intent(in), optional :: i_sub
    integer :: i_sub_opt, i_hel_opt
    i_sub_opt = 0; if (present (i_sub)) &
         i_sub_opt = i_sub
    i_hel_opt = 1; if (present (i_hel)) &
         i_hel_opt = i_hel
    index = 0
    if (.not. allocated (self%index)) then
       call msg_bug ("[qn_index_map_get_index] The index map is not allocated.")
    end if
    index = self%index(i_flv, i_hel_opt, i_sub_opt)
    if (index <= 0) then
       call self%write ()
       call msg_bug ("[qn_index_map_get_index] The index for the given quantum numbers could not be retrieved.")
    end if
  end function qn_index_map_get_index

  module function qn_index_map_get_n_flv (self) result (n_flv)
    class(qn_index_map_t), intent(in) :: self
    integer :: n_flv
    n_flv = self%n_flv
  end function qn_index_map_get_n_flv

  module function qn_index_map_get_n_hel (self) result (n_hel)
    class(qn_index_map_t), intent(in) :: self
    integer :: n_hel
    n_hel = self%n_hel
  end function qn_index_map_get_n_hel

  module function qn_index_map_get_n_sub (self) result (n_sub)
    class(qn_index_map_t), intent(in) :: self
    integer :: n_sub
    n_sub = self%n_sub
  end function qn_index_map_get_n_sub

  module function qn_index_map_get_index_by_qn (self, qn, i_sub) result (index)
    class(qn_index_map_t), intent(in) :: self
    integer :: index
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    integer, intent(in), optional :: i_sub
    integer :: i_qn
    if (size (qn) /= size (self%qn_flv, dim = 1)) &
         call msg_bug ("[qn_index_map_get_index_by_qn] number of particles does not match.")
    do i_qn = 1, self%n_flv
       if (all (qn .fmatch. self%qn_flv(:, i_qn))) then
          index = self%get_index (i_qn, i_sub = i_sub)
          return
       end if
    end do
    call self%write ()
    call msg_bug ("[qn_index_map_get_index_by_qn] The index for the given quantum &
         & numbers could not be retrieved.")
  end function qn_index_map_get_index_by_qn

  module function qn_index_map_get_sf_index_born (self, i_born, i_sub) result (index)
    class(qn_index_map_t), intent(in) :: self
    integer, intent(in) :: i_born, i_sub
    integer :: index
    index = self%sf_index_born(i_born, i_sub)
  end function qn_index_map_get_sf_index_born

  module function qn_index_map_get_sf_index_real (self, i_real, i_sub) result (index)
    class(qn_index_map_t), intent(in) :: self
    integer, intent(in) :: i_real, i_sub
    integer :: index
    index = self%sf_index_real(i_real, i_sub)
  end function qn_index_map_get_sf_index_real

  module subroutine external_link_set (link, int, i)
    type(external_link_t), intent(out) :: link
    type(interaction_t), target, intent(in) :: int
    integer, intent(in) :: i
    if (i /= 0) then
       link%int => int
       link%i = i
    end if
  end subroutine external_link_set

  module subroutine external_link_reassign (link, int_src, int_target)
    type(external_link_t), intent(inout) :: link
    type(interaction_t), intent(in) :: int_src
    type(interaction_t), intent(in), target :: int_target
    if (associated (link%int)) then
       if (link%int%tag == int_src%tag)  link%int => int_target
    end if
  end subroutine external_link_reassign

  module function external_link_is_set (link) result (flag)
    logical :: flag
    type(external_link_t), intent(in) :: link
    flag = associated (link%int)
  end function external_link_is_set

  module function external_link_get_ptr (link) result (int)
    type(interaction_t), pointer :: int
    type(external_link_t), intent(in) :: link
    int => link%int
  end function external_link_get_ptr

  module function external_link_get_index (link) result (i)
    integer :: i
    type(external_link_t), intent(in) :: link
    i = link%i
  end function external_link_get_index

  module function external_link_get_momentum_ptr (link) result (p)
    type(vector4_t), pointer :: p
    type(external_link_t), intent(in) :: link
    if (associated (link%int)) then
       p => link%int%p(link%i)
    else
       p => null ()
    end if
  end function external_link_get_momentum_ptr

  module subroutine internal_link_list_write (object, unit)
    class(internal_link_list_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    do i = 1, object%length
       write (u, "(1x,I0)", advance="no")  object%link(i)
    end do
  end subroutine internal_link_list_write

  module subroutine internal_link_list_append (link_list, link)
    class(internal_link_list_t), intent(inout) :: link_list
    integer, intent(in) :: link
    integer :: l, j
    integer, dimension(:), allocatable :: tmp
    l = link_list%length
    if (allocated (link_list%link)) then
       if (l == size (link_list%link)) then
          allocate (tmp (2 * l))
          tmp(:l) = link_list%link
          call move_alloc (from = tmp, to = link_list%link)
       end if
    else
       allocate (link_list%link (2))
    end if
    link_list%link(l+1) = link
    SHIFT_LINK_IN_PLACE: do j = l, 1, -1
       if (link >= link_list%link(j)) then
          exit SHIFT_LINK_IN_PLACE
       else
          link_list%link(j+1) = link_list%link(j)
          link_list%link(j) = link
       end if
    end do SHIFT_LINK_IN_PLACE
    link_list%length = l + 1
  end subroutine internal_link_list_append

  module function internal_link_list_has_entries (link_list) result (flag)
    class(internal_link_list_t), intent(in) :: link_list
    logical :: flag
    flag = link_list%length > 0
  end function internal_link_list_has_entries

  module function internal_link_list_get_length (link_list) result (length)
    class(internal_link_list_t), intent(in) :: link_list
    integer :: length
    length = link_list%length
  end function internal_link_list_get_length

  module function internal_link_list_get_link (link_list, i) result (link)
    class(internal_link_list_t), intent(in) :: link_list
    integer, intent(in) :: i
    integer :: link
    if (i <= link_list%length) then
       link = link_list%link(i)
    else
       call msg_bug ("Internal link list: out of bounds")
    end if
  end function internal_link_list_get_link

  module subroutine interaction_init &
       (int, n_in, n_vir, n_out, &
        tag, resonant, mask, hel_lock, set_relations, store_values)
    class(interaction_t), intent(out) :: int
    integer, intent(in) :: n_in, n_vir, n_out
    integer, intent(in), optional :: tag
    logical, dimension(:), intent(in), optional :: resonant
    type(quantum_numbers_mask_t), dimension(:), intent(in), optional :: mask
    integer, dimension(:), intent(in), optional :: hel_lock
    logical, intent(in), optional :: set_relations, store_values
    logical :: set_rel
    integer :: i, j
    set_rel = .false.;  if (present (set_relations))  set_rel = set_relations
    call interaction_set_tag (int, tag)
    call int%state_matrix%init (store_values)
    int%n_in = n_in
    int%n_vir = n_vir
    int%n_out = n_out
    int%n_tot = n_in + n_vir + n_out
    allocate (int%p_is_known (int%n_tot))
    int%p_is_known = .false.
    allocate (int%p (int%n_tot))
    allocate (int%source (int%n_tot))
    allocate (int%parents (int%n_tot))
    allocate (int%children (int%n_tot))
    allocate (int%resonant (int%n_tot))
    if (present (resonant)) then
       int%resonant = resonant
    else
       int%resonant = .false.
    end if
    allocate (int%mask (int%n_tot))
    allocate (int%hel_lock (int%n_tot))
    if (present (mask)) then
       int%mask = mask
    end if
    if (present (hel_lock)) then
       int%hel_lock = hel_lock
    else
       int%hel_lock = 0
    end if
    int%update_state_matrix = .false.
    int%update_values = .true.
    if (set_rel) then
       do i = 1, n_in
          do j = 1, n_out
             call int%relate (i, n_in + j)
          end do
       end do
    end if
  end subroutine interaction_init

  module subroutine interaction_init_qn_index_trivial (int)
    class(interaction_t), intent(inout) :: int
    call int%qn_index%init (int)
  end subroutine interaction_init_qn_index_trivial

  module subroutine interaction_init_qn_index_involved (int, qn_flv, n_sub, qn_hel)
    class(interaction_t), intent(inout) :: int
    type(quantum_numbers_t), dimension(:, :), intent(in) :: qn_flv
    integer, intent(in) :: n_sub
    type(quantum_numbers_t), dimension(:, :), intent(in), optional :: qn_hel
    call int%qn_index%init (int, qn_flv, n_sub, qn_hel)
  end subroutine interaction_init_qn_index_involved

  module subroutine interaction_init_qn_index_sf (int, qn_flv, n_flv_born, n_flv_real)
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: n_flv_born, n_flv_real
    type(quantum_numbers_t), dimension(:,:), intent(in) :: qn_flv
    call int%qn_index%init (int, qn_flv, n_flv_born, n_flv_real)
  end subroutine interaction_init_qn_index_sf

  module subroutine interaction_set_qn_index_helicity_flip (int, yorn)
    class(interaction_t), intent(inout) :: int
    logical, intent(in) :: yorn
    call int%qn_index%set_helicity_flip (yorn)
  end subroutine interaction_set_qn_index_helicity_flip

  module function interaction_get_qn_index (int, i_flv, i_hel, i_sub) result (index)
    class(interaction_t), intent(in) :: int
    integer :: index
    integer, intent(in) :: i_flv
    integer, intent(in), optional :: i_hel
    integer, intent(in), optional :: i_sub
    index = int%qn_index%get_index (i_flv, i_hel, i_sub)
  end function interaction_get_qn_index

  module function interaction_get_sf_qn_index_born (int, i_born, i_sub) result (index)
    class(interaction_t), intent(in) :: int
    integer :: index
    integer, intent(in) :: i_born, i_sub
    index = int%qn_index%get_sf_index_born (i_born, i_sub)
  end function interaction_get_sf_qn_index_born

  module function interaction_get_sf_qn_index_real (int, i_real, i_sub) result (index)
    class(interaction_t), intent(in) :: int
    integer :: index
    integer, intent(in) :: i_real, i_sub
    index = int%qn_index%get_sf_index_real (i_real, i_sub)
  end function interaction_get_sf_qn_index_real

  module function interaction_get_qn_index_n_flv (int) result (index)
    class(interaction_t), intent(in) :: int
    integer :: index
    index = int%qn_index%get_n_flv ()
  end function interaction_get_qn_index_n_flv

  module function interaction_get_qn_index_n_hel (int) result (index)
    class(interaction_t), intent(in) :: int
    integer :: index
    index = int%qn_index%get_n_hel ()
  end function interaction_get_qn_index_n_hel

  module function interaction_get_qn_index_n_sub (int) result (index)
    class(interaction_t), intent(in) :: int
    integer :: index
    index = int%qn_index%get_n_sub ()
  end function interaction_get_qn_index_n_sub

  module subroutine interaction_set_tag (int, tag)
    type(interaction_t), intent(inout), optional :: int
    integer, intent(in), optional :: tag
    integer, save :: stored_tag = 1
    if (present (int)) then
       if (present (tag)) then
          int%tag = tag
       else
          int%tag = stored_tag
          stored_tag = stored_tag + 1
       end if
    else if (present (tag)) then
       stored_tag = tag
    else
       stored_tag = 1
    end if
  end subroutine interaction_set_tag

  module subroutine reset_interaction_counter (tag)
    integer, intent(in), optional :: tag
    call interaction_set_tag (tag=tag)
  end subroutine reset_interaction_counter

  module subroutine interaction_final (object)
    class(interaction_t), intent(inout) :: object
    call object%state_matrix%final ()
  end subroutine interaction_final

  module subroutine interaction_write &
       (int, unit, verbose, show_momentum_sum, show_mass, show_state, &
       col_verbose, testflag)
    class(interaction_t), intent(in) :: int
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, show_momentum_sum, show_mass
    logical, intent(in), optional :: show_state, col_verbose, testflag
    integer :: u
    integer :: i, index_link
    type(interaction_t), pointer :: int_link
    logical :: show_st
    u = given_output_unit (unit);  if (u < 0)  return
    show_st = .true.;  if (present (show_state))  show_st = show_state
    if (int%tag /= 0) then
       write (u, "(1x,A,I0)")  "Interaction: ", int%tag
       do i = 1, int%n_tot
          if (i == 1 .and. int%n_in > 0) then
             write (u, "(1x,A)") "Incoming:"
          else if (i == int%n_in + 1 .and. int%n_vir > 0) then
             write (u, "(1x,A)") "Virtual:"
          else if (i == int%n_in + int%n_vir + 1 .and. int%n_out > 0) then
             write (u, "(1x,A)") "Outgoing:"
          end if
          write (u, "(1x,A,1x,I0)", advance="no") "Particle", i
          if (allocated (int%resonant)) then
             if (int%resonant(i)) then
                write (u, "(A)") "[r]"
             else
                write (u, *)
             end if
          else
             write (u, *)
          end if
          if (allocated (int%p)) then
             if (int%p_is_known(i)) then
                call vector4_write (int%p(i), u, show_mass, testflag)
             else
                write (u, "(A)")  "  [momentum undefined]"
             end if
          else
             write (u, "(A)") "  [momentum not allocated]"
          end if
          if (allocated (int%mask)) then
             write (u, "(1x,A)", advance="no")  "mask [fch] = "
             call int%mask(i)%write (u)
             write (u, *)
          end if
          if (int%parents(i)%has_entries () &
               .or. int%children(i)%has_entries ()) then
             write (u, "(1x,A)", advance="no") "internal links:"
             call int%parents(i)%write (u)
             if (int%parents(i)%has_entries ()) &
                  write (u, "(1x,A)", advance="no") "=>"
             write (u, "(1x,A)", advance="no") "X"
             if (int%children(i)%has_entries ()) &
                  write (u, "(1x,A)", advance="no") "=>"
             call int%children(i)%write (u)
             write (u, *)
          end if
          if (allocated (int%hel_lock)) then
             if (int%hel_lock(i) /= 0) then
                write (u, "(1x,A,1x,I0)")  "helicity lock:", int%hel_lock(i)
             end if
          end if
          if (external_link_is_set (int%source(i))) then
             write (u, "(1x,A)", advance="no") "source:"
             int_link => external_link_get_ptr (int%source(i))
             index_link = external_link_get_index (int%source(i))
             write (u, "(1x,'(',I0,')',I0)", advance="no") &
                  int_link%tag, index_link
             write (u, *)
          end if
       end do
       if (present (show_momentum_sum)) then
          if (allocated (int%p) .and. show_momentum_sum) then
             write (u, "(1x,A)") "Incoming particles (sum):"
             call vector4_write &
                  (sum (int%p(1 : int%n_in)), u, show_mass = show_mass)
             write (u, "(1x,A)") "Outgoing particles (sum):"
             call vector4_write &
                  (sum (int%p(int%n_in + int%n_vir + 1 : )), &
                   u, show_mass = show_mass)
             write (u, *)
          end if
       end if
       if (show_st) then
          call int%write_state_matrix (write_value_list = verbose, &
             verbose = verbose, unit = unit, col_verbose = col_verbose, &
             testflag = testflag)
       end if
    else
       write (u, "(1x,A)") "Interaction: [empty]"
    end if
  end subroutine interaction_write

  module subroutine interaction_write_state_matrix (int, unit, write_value_list, &
     verbose, col_verbose, testflag)
    class(interaction_t), intent(in) :: int
    logical, intent(in), optional :: write_value_list, verbose, col_verbose
    logical, intent(in), optional :: testflag
    integer, intent(in), optional :: unit
    call int%state_matrix%write (write_value_list = verbose, &
       verbose = verbose, unit = unit, col_verbose = col_verbose, &
       testflag = testflag)
  end subroutine interaction_write_state_matrix

  module subroutine interaction_reduce_state_matrix (int, qn_mask, keep_order)
    class(interaction_t), intent(inout) :: int
    type(quantum_numbers_mask_t), intent(in), dimension(:) :: qn_mask
    logical, optional, intent(in) :: keep_order
    type(state_matrix_t) :: state
    logical :: opt_keep_order
    opt_keep_order = .false.
    if (present (keep_order)) opt_keep_order = keep_order
    call int%state_matrix%reduce (qn_mask, state, keep_me_index = keep_order)
    int%state_matrix = state
    if (opt_keep_order) then
       call int%state_matrix%reorder_me (state)
       int%state_matrix = state
    end if
  end subroutine interaction_reduce_state_matrix

  module subroutine interaction_assign (int_out, int_in)
    type(interaction_t), intent(out) :: int_out
    type(interaction_t), intent(in), target :: int_in
    call interaction_set_tag (int_out)
    int_out%state_matrix = int_in%state_matrix
    int_out%n_in  = int_in%n_in
    int_out%n_out = int_in%n_out
    int_out%n_vir = int_in%n_vir
    int_out%n_tot = int_in%n_tot
    if (allocated (int_in%p_is_known)) then
       allocate (int_out%p_is_known (size (int_in%p_is_known)))
       int_out%p_is_known = int_in%p_is_known
    end if
    if (allocated (int_in%p)) then
       allocate (int_out%p (size (int_in%p)))
       int_out%p = int_in%p
    end if
    if (allocated (int_in%source)) then
       allocate (int_out%source (size (int_in%source)))
       int_out%source = int_in%source
    end if
    if (allocated (int_in%parents)) then
       allocate (int_out%parents (size (int_in%parents)))
       int_out%parents = int_in%parents
    end if
    if (allocated (int_in%children)) then
       allocate (int_out%children (size (int_in%children)))
       int_out%children = int_in%children
    end if
    if (allocated (int_in%resonant)) then
       allocate (int_out%resonant (size (int_in%resonant)))
       int_out%resonant = int_in%resonant
    end if
    if (allocated (int_in%mask)) then
       allocate (int_out%mask (size (int_in%mask)))
       int_out%mask = int_in%mask
    end if
    if (allocated (int_in%hel_lock)) then
       allocate (int_out%hel_lock (size (int_in%hel_lock)))
       int_out%hel_lock = int_in%hel_lock
    end if
    int_out%update_state_matrix = int_in%update_state_matrix
    int_out%update_values = int_in%update_values
  end subroutine interaction_assign

  module subroutine interaction_add_state &
       (int, qn, index, value, sum_values, counter_index, ignore_sub_for_qn, me_index)
    class(interaction_t), intent(inout) :: int
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    integer, intent(in), optional :: index
    complex(default), intent(in), optional :: value
    logical, intent(in), optional :: sum_values
    integer, intent(in), optional :: counter_index
    logical, intent(in), optional :: ignore_sub_for_qn
    integer, intent(out), optional :: me_index
    type(quantum_numbers_t), dimension(size(qn)) :: qn_tmp
    qn_tmp = qn
    call qn_tmp%undefine (int%mask)
    call int%state_matrix%add_state (qn_tmp, index, value, sum_values, &
         counter_index, ignore_sub_for_qn, me_index)
    int%update_values = .true.
  end subroutine interaction_add_state

  module subroutine interaction_set_duplicate_flv_zero (int)
    class(interaction_t), intent(inout) :: int
    call int%state_matrix%set_duplicate_flv_zero ()
  end subroutine interaction_set_duplicate_flv_zero

  module subroutine interaction_freeze (int)
    class(interaction_t), intent(inout) :: int
    if (int%update_state_matrix) then
       call int%state_matrix%collapse (int%mask)
       int%update_state_matrix = .false.
       int%update_values = .true.
    end if
    if (int%update_values) then
       call int%state_matrix%freeze ()
       int%update_values = .false.
    end if
  end subroutine interaction_freeze

  pure module function interaction_is_empty (int) result (flag)
    logical :: flag
    class(interaction_t), intent(in) :: int
    flag = int%state_matrix%is_empty ()
  end function interaction_is_empty

  pure module function interaction_get_n_matrix_elements (int) result (n)
    integer :: n
    class(interaction_t), intent(in) :: int
    n = int%state_matrix%get_n_matrix_elements ()
  end function interaction_get_n_matrix_elements

  module function interaction_get_state_depth (int) result (n)
    integer :: n
    class(interaction_t), intent(in) :: int
    n = int%state_matrix%get_depth ()
  end function interaction_get_state_depth

  module function interaction_get_n_in_helicities (int) result (n_hel)
    integer :: n_hel
    class(interaction_t), intent(in) :: int
    type(interaction_t) :: int_copy
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn
    integer :: i
    allocate (qn_mask (int%n_tot))
    do i = 1, int%n_tot
       if (i <= int%n_in) then
          call qn_mask(i)%init (.true., .true., .false.)
       else
          call qn_mask(i)%init (.true., .true., .true.)
       end if
    end do
    int_copy = int
    call int_copy%set_mask (qn_mask)
    call int_copy%freeze ()
    allocate (qn (int_copy%state_matrix%get_n_matrix_elements (), &
         int_copy%state_matrix%get_depth ()))
    qn = int_copy%get_quantum_numbers ()
    n_hel = 0
    do i = 1, size (qn, dim=1)
       if (all (qn(:, i)%get_subtraction_index () == 0)) n_hel = n_hel + 1
    end do
    call int_copy%final ()
    deallocate (qn_mask)
    deallocate (qn)
  end function interaction_get_n_in_helicities

  pure module function interaction_get_me_size (int) result (n)
    integer :: n
    class(interaction_t), intent(in) :: int
    n = int%state_matrix%get_me_size ()
  end function interaction_get_me_size

  pure module function interaction_get_norm (int) result (norm)
    real(default) :: norm
    class(interaction_t), intent(in) :: int
    norm = int%state_matrix%get_norm ()
  end function interaction_get_norm

  module function interaction_get_n_sub (int) result (n_sub)
    integer :: n_sub
    class(interaction_t), intent(in) :: int
    n_sub = int%state_matrix%get_n_sub ()
  end function interaction_get_n_sub

  module function interaction_get_quantum_numbers_single (int, i, by_me_index) result (qn)
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    class(interaction_t), intent(in), target :: int
    integer, intent(in) :: i
    logical, intent(in), optional :: by_me_index
    allocate (qn (int%state_matrix%get_depth ()))
    qn = int%state_matrix%get_quantum_number (i, by_me_index)
  end function interaction_get_quantum_numbers_single

  module function interaction_get_quantum_numbers_all (int) result (qn)
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn
    class(interaction_t), intent(in), target :: int
    integer :: i
    allocate (qn (int%state_matrix%get_depth(), &
         int%state_matrix%get_n_matrix_elements ()))
    do i = 1, int%state_matrix%get_n_matrix_elements ()
       qn (:, i) = int%state_matrix%get_quantum_number (i)
    end do
  end function interaction_get_quantum_numbers_all

  module function interaction_get_quantum_numbers_all_qn_mask (int, qn_mask) &
     result (qn)
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn
    class(interaction_t), intent(in) :: int
    type(quantum_numbers_mask_t), intent(in) :: qn_mask
    integer :: n_redundant, n_all, n_me
    integer :: i
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn_all
    call int%state_matrix%get_quantum_numbers (qn_all)
    n_redundant = count (qn_all%are_redundant (qn_mask))
    n_all = size (qn_all)
    !!! Number of matrix elements = survivors / n_particles
    n_me = (n_all - n_redundant) / int%state_matrix%get_depth ()
    allocate (qn (int%state_matrix%get_depth(), n_me))
    do i = 1, n_me
       if (.not. any (qn_all(i, :)%are_redundant (qn_mask))) &
          qn (:, i) = qn_all (i, :)
    end do
  end function interaction_get_quantum_numbers_all_qn_mask

  module subroutine interaction_get_quantum_numbers_all_sub (int, qn)
    class(interaction_t), intent(in) :: int
    type(quantum_numbers_t), dimension(:,:), allocatable, intent(out) :: qn
    integer :: i
    allocate (qn (int%state_matrix%get_depth(), &
         int%state_matrix%get_n_matrix_elements ()))
    do i = 1, int%state_matrix%get_n_matrix_elements ()
       qn (:, i) = int%state_matrix%get_quantum_number (i)
    end do
  end subroutine interaction_get_quantum_numbers_all_sub

  module subroutine interaction_get_flavors (int, only_elementary, qn_mask, flv)
    class(interaction_t), intent(in), target :: int
    logical, intent(in) :: only_elementary
    type(quantum_numbers_mask_t), intent(in), dimension(:), optional :: qn_mask
    integer, intent(out), dimension(:,:), allocatable :: flv
    call int%state_matrix%get_flavors (only_elementary, qn_mask, flv)
  end subroutine interaction_get_flavors

  module subroutine interaction_get_quantum_numbers_mask (int, qn_mask, qn)
    class(interaction_t), intent(in) :: int
    type(quantum_numbers_mask_t), intent(in) :: qn_mask
    type(quantum_numbers_t), dimension(:,:), allocatable, intent(out) :: qn
    integer :: n_redundant, n_all, n_me
    integer :: i
    type(quantum_numbers_t), dimension(:,:), allocatable :: qn_all
    call int%state_matrix%get_quantum_numbers (qn_all)
    n_redundant = count (qn_all%are_redundant (qn_mask))
    n_all = size (qn_all)
    !!! Number of matrix elements = survivors / n_particles
    n_me = (n_all - n_redundant) / int%state_matrix%get_depth ()
    allocate (qn (int%state_matrix%get_depth(), n_me))
    do i = 1, n_me
       if (.not. any (qn_all(i, :)%are_redundant (qn_mask))) &
          qn (:, i) = qn_all (i, :)
    end do
  end subroutine interaction_get_quantum_numbers_mask

  elemental module function interaction_get_matrix_element_single (int, i) result (me)
    complex(default) :: me
    class(interaction_t), intent(in) :: int
    integer, intent(in) :: i
    me = int%state_matrix%get_matrix_element (i)
  end function interaction_get_matrix_element_single

  module function interaction_get_matrix_element_array (int) result (me)
    complex(default), dimension(:), allocatable :: me
    class(interaction_t), intent(in) :: int
    allocate (me (int%get_n_matrix_elements ()))
    me = int%state_matrix%get_matrix_element ()
  end function interaction_get_matrix_element_array

  module subroutine interaction_set_matrix_element_qn (int, qn, val)
    class(interaction_t), intent(inout) :: int
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    complex(default), intent(in) :: val
    call int%state_matrix%set_matrix_element (qn, val)
  end subroutine interaction_set_matrix_element_qn

  module subroutine interaction_set_matrix_element_all (int, value)
    class(interaction_t), intent(inout) :: int
    complex(default), intent(in) :: value
    call int%state_matrix%set_matrix_element (value)
  end subroutine interaction_set_matrix_element_all

  module subroutine interaction_set_matrix_element_array (int, value, range)
    class(interaction_t), intent(inout) :: int
    complex(default), intent(in), dimension(:) :: value
    integer, intent(in), dimension(:), optional :: range
    call int%state_matrix%set_matrix_element (value, range)
  end subroutine interaction_set_matrix_element_array

  pure module subroutine interaction_set_matrix_element_single (int, i, value)
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: i
    complex(default), intent(in) :: value
    call int%state_matrix%set_matrix_element (i, value)
  end subroutine interaction_set_matrix_element_single

  module subroutine interaction_set_matrix_element_clone (int, int1)
    class(interaction_t), intent(inout) :: int
    class(interaction_t), intent(in) :: int1
    call int%state_matrix%set_matrix_element (int1%state_matrix)
  end subroutine interaction_set_matrix_element_clone

  module subroutine interaction_set_only_matrix_element (int, i, value)
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: i
    complex(default), intent(in) :: value
    call int%set_matrix_element (cmplx (0, 0, default))
    call int%set_matrix_element (i, value)
  end subroutine interaction_set_only_matrix_element

  module subroutine interaction_add_to_matrix_element (int, qn, value, match_only_flavor)
    class(interaction_t), intent(inout) :: int
    type(quantum_numbers_t), dimension(:), intent(in) :: qn
    complex(default), intent(in) :: value
    logical, intent(in), optional :: match_only_flavor
    call int%state_matrix%add_to_matrix_element (qn, value, match_only_flavor)
  end subroutine interaction_add_to_matrix_element

  module subroutine interaction_get_diagonal_entries (int, i)
    class(interaction_t), intent(in) :: int
    integer, dimension(:), allocatable, intent(out) :: i
    call int%state_matrix%get_diagonal_entries (i)
  end subroutine interaction_get_diagonal_entries

  module subroutine interaction_normalize_by_trace (int)
    class(interaction_t), intent(inout) :: int
    call int%state_matrix%normalize_by_trace ()
  end subroutine interaction_normalize_by_trace

  module subroutine interaction_normalize_by_max (int)
    class(interaction_t), intent(inout) :: int
    call int%state_matrix%normalize_by_max ()
  end subroutine interaction_normalize_by_max

  module subroutine interaction_set_norm (int, norm)
    class(interaction_t), intent(inout) :: int
    real(default), intent(in) :: norm
    call int%state_matrix%set_norm (norm)
  end subroutine interaction_set_norm

  module subroutine interaction_set_state_matrix (int, state)
    class(interaction_t), intent(inout) :: int
    type(state_matrix_t), intent(in) :: state
    int%state_matrix = state
  end subroutine interaction_set_state_matrix

  module function interaction_get_max_color_value (int) result (cmax)
    class(interaction_t), intent(in) :: int
    integer :: cmax
    cmax = int%state_matrix%get_max_color_value ()
  end function interaction_get_max_color_value

  module subroutine interaction_factorize &
       (int, mode, x, ok, single_state, correlated_state, qn_in)
    class(interaction_t), intent(in), target :: int
    integer, intent(in) :: mode
    real(default), intent(in) :: x
    logical, intent(out) :: ok
    type(state_matrix_t), &
         dimension(:), allocatable, intent(out) :: single_state
    type(state_matrix_t), intent(out), optional :: correlated_state
    type(quantum_numbers_t), dimension(:), intent(in), optional :: qn_in
    call int%state_matrix%factorize &
         (mode, x, ok, single_state, correlated_state, qn_in)
  end subroutine interaction_factorize

  module function interaction_sum (int) result (value)
    class(interaction_t), intent(in) :: int
    complex(default) :: value
    value = int%state_matrix%sum ()
  end function interaction_sum

  module subroutine interaction_add_color_contractions (int)
    class(interaction_t), intent(inout) :: int
    call int%state_matrix%add_color_contractions ()
  end subroutine interaction_add_color_contractions

  pure module subroutine interaction_evaluate_product &
       (int, i, int1, int2, index1, index2)
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: i
    type(interaction_t), intent(in) :: int1, int2
    integer, dimension(:), intent(in) :: index1, index2
    call int%state_matrix%evaluate_product &
         (i, int1%state_matrix, int2%state_matrix, &
          index1, index2)
  end subroutine interaction_evaluate_product

  pure module subroutine interaction_evaluate_product_cf &
       (int, i, int1, int2, index1, index2, factor)
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: i
    type(interaction_t), intent(in) :: int1, int2
    integer, dimension(:), intent(in) :: index1, index2
    complex(default), dimension(:), intent(in) :: factor
    call int%state_matrix%evaluate_product_cf &
         (i, int1%state_matrix, int2%state_matrix, &
          index1, index2, factor)
  end subroutine interaction_evaluate_product_cf

  pure module subroutine interaction_evaluate_square_c (int, i, int1, index1)
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: i
    type(interaction_t), intent(in) :: int1
    integer, dimension(:), intent(in) :: index1
    call int%state_matrix%evaluate_square_c (i, int1%state_matrix, index1)
  end subroutine interaction_evaluate_square_c

  pure module subroutine interaction_evaluate_sum (int, i, int1, index1)
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: i
    type(interaction_t), intent(in) :: int1
    integer, dimension(:), intent(in) :: index1
    call int%state_matrix%evaluate_sum (i, int1%state_matrix, index1)
  end subroutine interaction_evaluate_sum

  pure module subroutine interaction_evaluate_me_sum (int, i, int1, index1)
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: i
    type(interaction_t), intent(in) :: int1
    integer, dimension(:), intent(in) :: index1
    call int%state_matrix%evaluate_me_sum (i, int1%state_matrix, index1)
  end subroutine interaction_evaluate_me_sum

  module subroutine interaction_tag_hard_process (int, tag)
    class(interaction_t), intent(inout) :: int
    integer, dimension(:), intent(in), optional :: tag
    type(state_matrix_t) :: state
    call int%state_matrix%tag_hard_process (state, tag)
    call int%state_matrix%final ()
    int%state_matrix = state
  end subroutine interaction_tag_hard_process

  module subroutine interaction_retag_hard_process (int, i, hard)
    class(interaction_t), intent(inout), target :: int
    integer, intent(in) :: i
    logical, intent(in) :: hard
    type(state_iterator_t) :: it
    call it%init (int%get_state_matrix_ptr ())
    do while (it%is_valid ())
       call it%retag_hard_process (i, hard)
       call it%advance ()
    end do
  end subroutine interaction_retag_hard_process

  module function interaction_get_tag (int) result (tag)
    class(interaction_t), intent(in) :: int
    integer :: tag
    tag = int%tag
  end function interaction_get_tag

  pure module function interaction_get_n_tot (object) result (n_tot)
    class(interaction_t), intent(in) :: object
    integer :: n_tot
    n_tot = object%n_tot
  end function interaction_get_n_tot

  pure module function interaction_get_n_in (object) result (n_in)
    class(interaction_t), intent(in) :: object
    integer :: n_in
    n_in = object%n_in
  end function interaction_get_n_in

  pure module function interaction_get_n_vir (object) result (n_vir)
    class(interaction_t), intent(in) :: object
    integer :: n_vir
    n_vir = object%n_vir
  end function interaction_get_n_vir

  pure module function interaction_get_n_out (object) result (n_out)
    class(interaction_t), intent(in) :: object
    integer :: n_out
    n_out = object%n_out
  end function interaction_get_n_out

  module function idx (int, i, outgoing)
    integer :: idx
    type(interaction_t), intent(in) :: int
    integer, intent(in) :: i
    logical, intent(in), optional :: outgoing
    logical :: in, vir, out
    if (present (outgoing)) then
       in  = .not. outgoing
       vir = .false.
       out = outgoing
    else
       in = .true.
       vir = .true.
       out = .true.
    end if
    idx = 0
    if (in) then
       if (vir) then
          if (out) then
             if (i <= int%n_tot)  idx = i
          else
             if (i <= int%n_in + int%n_vir)  idx = i
          end if
       else if (out) then
          if (i <= int%n_in) then
             idx = i
          else if (i <= int%n_in + int%n_out) then
             idx = int%n_vir + i
          end if
       else
          if (i <= int%n_in)  idx = i
       end if
    else if (vir) then
       if (out) then
          if (i <= int%n_vir + int%n_out)  idx = int%n_in + i
       else
          if (i <= int%n_vir)  idx = int%n_in + i
       end if
    else if (out) then
       if (i <= int%n_out)  idx = int%n_in + int%n_vir + i
    end if
    if (idx == 0) then
       call int%basic_write ()
       print *, i, in, vir, out
       call msg_bug (" Momentum index is out of range for this interaction")
    end if
  end function idx

  module function interaction_get_momenta_all (int, outgoing) result (p)
    class(interaction_t), intent(in) :: int
    type(vector4_t), dimension(:), allocatable :: p
    logical, intent(in), optional :: outgoing
    integer :: i
    if (present (outgoing)) then
       if (outgoing) then
          allocate (p (int%n_out))
       else
          allocate (p (int%n_in))
       end if
    else
       allocate (p (int%n_tot))
    end if
    do i = 1, size (p)
       p(i) = int%p(idx (int, i, outgoing))
    end do
  end function interaction_get_momenta_all

  module function interaction_get_momenta_idx (int, jj) result (p)
    class(interaction_t), intent(in) :: int
    type(vector4_t), dimension(:), allocatable :: p
    integer, dimension(:), intent(in) :: jj
    allocate (p (size (jj)))
    p = int%p(jj)
  end function interaction_get_momenta_idx

  module function interaction_get_momentum (int, i, outgoing) result (p)
    class(interaction_t), intent(in) :: int
    type(vector4_t) :: p
    integer, intent(in) :: i
    logical, intent(in), optional :: outgoing
    p = int%p(idx (int, i, outgoing))
  end function interaction_get_momentum

  module function interaction_get_state_matrix_ptr (int) result (state)
    class(interaction_t), intent(in), target :: int
    type(state_matrix_t), pointer :: state
    state => int%state_matrix
  end function interaction_get_state_matrix_ptr

  module function interaction_get_resonance_flags (int) result (resonant)
    class(interaction_t), intent(in) :: int
    logical, dimension(size(int%resonant)) :: resonant
    resonant = int%resonant
  end function interaction_get_resonance_flags

  module function interaction_get_mask_all (int) result (mask)
    class(interaction_t), intent(in) :: int
    type(quantum_numbers_mask_t), dimension(size(int%mask)) :: mask
    mask = int%mask
  end function interaction_get_mask_all

  module function interaction_get_mask_slice (int, index) result (mask)
    class(interaction_t), intent(in) :: int
    integer, dimension(:), intent(in) :: index
    type(quantum_numbers_mask_t), dimension(size(index)) :: mask
    mask = int%mask(index)
  end function interaction_get_mask_slice

  module function interaction_get_s (int) result (s)
    real(default) :: s
    class(interaction_t), intent(in) :: int
    if (int%n_in /= 0) then
       s = sum (int%p(:int%n_in)) ** 2
    else
       s = sum (int%p(int%n_vir + 1 : )) ** 2
    end if
  end function interaction_get_s

  module function interaction_get_cm_transformation (int) result (lt)
    type(lorentz_transformation_t) :: lt
    class(interaction_t), intent(in) :: int
    type(vector4_t) :: p_cm
    real(default) :: s
    if (int%n_in /= 0) then
       p_cm = sum (int%p(:int%n_in))
    else
       p_cm = sum (int%p(int%n_vir+1:))
    end if
    s = p_cm ** 2
    if (s > 0) then
       lt = boost (p_cm, sqrt (s))
    else
       lt = identity
    end if
  end function interaction_get_cm_transformation

  module subroutine interaction_get_unstable_particle (int, flv, p, i)
    class(interaction_t), intent(in), target :: int
    type(flavor_t), intent(out) :: flv
    type(vector4_t), intent(out) :: p
    integer, intent(out) :: i
    type(state_iterator_t) :: it
    type(flavor_t), dimension(int%n_tot) :: flv_array
    call it%init (int%state_matrix)
    flv_array = it%get_flavor ()
    do i = int%n_in + int%n_vir + 1, int%n_tot
       if (.not. flv_array(i)%is_stable ()) then
          flv = flv_array(i)
          p = int%p(i)
          return
       end if
    end do
  end subroutine interaction_get_unstable_particle

  module subroutine interaction_get_flv_out (int, flv)
    class(interaction_t), intent(in), target :: int
    type(flavor_t), dimension(:,:), allocatable, intent(out) :: flv
    type(state_iterator_t) :: it
    type(flavor_t), dimension(:), allocatable :: flv_state
    integer :: n_in, n_vir, n_out, n_tot, n_state, i
    n_in = int%get_n_in ()
    n_vir = int%get_n_vir ()
    n_out = int%get_n_out ()
    n_tot = int%get_n_tot ()
    n_state = int%get_n_matrix_elements ()
    allocate (flv (n_out, n_state))
    allocate (flv_state (n_tot))
    i = 1
    call it%init (int%get_state_matrix_ptr ())
    do while (it%is_valid ())
       flv_state = it%get_flavor ()
       flv(:,i) = flv_state(n_in + n_vir + 1 : )
       i = i + 1
       call it%advance ()
    end do
  end subroutine interaction_get_flv_out

  module subroutine interaction_get_flv_content (int, state_flv, n_out_hard)
    class(interaction_t), intent(in), target :: int
    type(state_flv_content_t), intent(out) :: state_flv
    integer, intent(in) :: n_out_hard
    logical, dimension(:), allocatable :: mask
    integer :: n_tot
    n_tot = int%get_n_tot ()
    allocate (mask (n_tot), source = .false.)
    mask(n_tot-n_out_hard + 1 : ) = .true.
    call state_flv%fill (int%get_state_matrix_ptr (), mask)
  end subroutine interaction_get_flv_content

  module subroutine interaction_set_mask (int, mask)
    class(interaction_t), intent(inout) :: int
    type(quantum_numbers_mask_t), dimension(:), intent(in) :: mask
    if (size (int%mask) /= size (mask)) &
       call msg_fatal ("Attempting to set mask with unfitting size!")
    int%mask = mask
    int%update_state_matrix = .true.
  end subroutine interaction_set_mask

  subroutine interaction_merge_mask_entry (int, i, mask)
    type(interaction_t), intent(inout) :: int
    integer, intent(in) :: i
    type(quantum_numbers_mask_t), intent(in) :: mask
    type(quantum_numbers_mask_t) :: mask_tmp
    integer :: ii
    ii = idx (int, i)
    if (int%mask(ii) .neqv. mask) then
       int%mask(ii) = int%mask(ii) .or. mask
       if (int%hel_lock(ii) /= 0) then
          call mask_tmp%assign (mask, helicity=.true.)
          int%mask(int%hel_lock(ii)) = int%mask(int%hel_lock(ii)) .or. mask_tmp
       end if
    end if
    int%update_state_matrix = .true.
  end subroutine interaction_merge_mask_entry

  module subroutine interaction_reset_momenta (int)
    class(interaction_t), intent(inout) :: int
    int%p = vector4_null
    int%p_is_known = .true.
  end subroutine interaction_reset_momenta

  module subroutine interaction_set_momenta (int, p, outgoing)
    class(interaction_t), intent(inout) :: int
    type(vector4_t), dimension(:), intent(in) :: p
    logical, intent(in), optional :: outgoing
    integer :: i, index
    do i = 1, size (p)
       index = idx (int, i, outgoing)
       int%p(index) = p(i)
       int%p_is_known(index) = .true.
    end do
  end subroutine interaction_set_momenta

  module subroutine interaction_set_momentum (int, p, i, outgoing)
    class(interaction_t), intent(inout) :: int
    type(vector4_t), intent(in) :: p
    integer, intent(in) :: i
    logical, intent(in), optional :: outgoing
    integer :: index
    index = idx (int, i, outgoing)
    int%p(index) = p
    int%p_is_known(index) = .true.
  end subroutine interaction_set_momentum

  module subroutine interaction_set_flavored_values (int, value, flv_in, pos)
    class(interaction_t), intent(inout) :: int
    complex(default), dimension(:), intent(in) :: value
    type(flavor_t), dimension(:), intent(in) :: flv_in
    integer, intent(in) :: pos
    type(state_iterator_t) :: it
    type(flavor_t) :: flv
    integer :: i
    if (size (value) == 1) then
       call int%set_matrix_element (value(1))
    else
       call it%init (int%state_matrix)
       do while (it%is_valid ())
          flv = it%get_flavor (pos)
          SCAN_FLV: do i = 1, size (value)
             if (flv == flv_in(i)) then
                call it%set_matrix_element (value(i))
                exit SCAN_FLV
             end if
          end do SCAN_FLV
          call it%advance ()
       end do
    end if
  end subroutine interaction_set_flavored_values

  module subroutine interaction_relate (int, i1, i2)
    class(interaction_t), intent(inout), target :: int
    integer, intent(in) :: i1, i2
    if (i1 /= 0 .and. i2 /= 0) then
       call int%children(i1)%append (i2)
       call int%parents(i2)%append (i1)
    end if
  end subroutine interaction_relate

  module subroutine interaction_transfer_relations (int1, int2, map)
    class(interaction_t), intent(in) :: int1
    class(interaction_t), intent(inout), target :: int2
    integer, dimension(:), intent(in) :: map
    integer :: i, j, k
    do i = 1, size (map)
       do j = 1, int1%parents(i)%get_length ()
          k = int1%parents(i)%get_link (j)
          call int2%relate (map(k), map(i))
       end do
       if (map(i) /= 0) then
          int2%resonant(map(i)) = int1%resonant(i)
       end if
    end do
  end subroutine interaction_transfer_relations

  module subroutine interaction_relate_connections &
       (int, int_in, connection_index, &
        map, map_connections, resonant)
    class(interaction_t), intent(inout), target :: int
    class(interaction_t), intent(in) :: int_in
    integer, dimension(:), intent(in) :: connection_index
    integer, dimension(:), intent(in) :: map, map_connections
    logical, intent(in), optional :: resonant
    logical :: reson
    integer :: i, j, i2, k2
    reson = .false.;  if (present (resonant))  reson = resonant
    do i = 1, size (map_connections)
       k2 = connection_index(i)
       do j = 1, int_in%children(k2)%get_length ()
          i2 = int_in%children(k2)%get_link (j)
          call int%relate (map_connections(i), map(i2))
          if (reson)  call int%retag_hard_process (map(i2), .false.)
       end do
       int%resonant(map_connections(i)) = reson
    end do
  end subroutine interaction_relate_connections

  module function interaction_get_n_children (int, i) result (n)
    integer :: n
    class(interaction_t), intent(in) :: int
    integer, intent(in) :: i
    n = int%children(i)%get_length ()
  end function interaction_get_n_children

  module function interaction_get_n_parents (int, i) result (n)
    integer :: n
    class(interaction_t), intent(in) :: int
    integer, intent(in) :: i
    n = int%parents(i)%get_length ()
  end function interaction_get_n_parents

  module function interaction_get_children (int, i) result (idx)
    integer, dimension(:), allocatable :: idx
    class(interaction_t), intent(in) :: int
    integer, intent(in) :: i
    integer :: k, l
    l = int%children(i)%get_length ()
    allocate (idx (l))
    do k = 1, l
       idx(k) = int%children(i)%get_link (k)
    end do
  end function interaction_get_children

  module function interaction_get_parents (int, i) result (idx)
    integer, dimension(:), allocatable :: idx
    class(interaction_t), intent(in) :: int
    integer, intent(in) :: i
    integer :: k, l
    l = int%parents(i)%get_length ()
    allocate (idx (l))
    do k = 1, l
       idx(k) = int%parents(i)%get_link (k)
    end do
  end function interaction_get_parents

  module subroutine interaction_set_source_link (int, i, int1, i1)
    class(interaction_t), intent(inout) :: int
    integer, intent(in) :: i
    class(interaction_t), intent(in), target :: int1
    integer, intent(in) :: i1
    if (i /= 0)  call external_link_set (int%source(i), int1, i1)
  end subroutine interaction_set_source_link

  module subroutine interaction_reassign_links (int, int_src, int_target)
    type(interaction_t), intent(inout) :: int
    type(interaction_t), intent(in) :: int_src
    type(interaction_t), intent(in), target :: int_target
    integer :: i
    if (allocated (int%source)) then
       do i = 1, size (int%source)
          call external_link_reassign (int%source(i), int_src, int_target)
       end do
    end if
  end subroutine interaction_reassign_links

  module function interaction_find_link (int, int1, i1) result (i)
    integer :: i
    type(interaction_t), intent(in) :: int, int1
    integer, intent(in) :: i1
    type(interaction_t), pointer :: int_tmp
    do i = 1, int%n_tot
       int_tmp => external_link_get_ptr (int%source(i))
       if (int_tmp%tag == int1%tag) then
          if (external_link_get_index (int%source(i)) == i1)  return
       end if
    end do
    i = 0
  end function interaction_find_link

  module subroutine interaction_find_source (int, i, int1, i1)
    class(interaction_t), intent(in) :: int
    integer, intent(in) :: i
    type(interaction_t), intent(out), pointer :: int1
    integer, intent(out) :: i1
    type(external_link_t) :: link
    link = interaction_get_ultimate_source (int, i)
    int1 => external_link_get_ptr (link)
    i1 = external_link_get_index (link)
  end subroutine interaction_find_source

  module function interaction_get_ultimate_source (int, i) result (link)
    type(external_link_t) :: link
    type(interaction_t), intent(in) :: int
    integer, intent(in) :: i
    type(interaction_t), pointer :: int_src
    integer :: i_src
    link = int%source(i)
    if (external_link_is_set (link)) then
       do
          int_src => external_link_get_ptr (link)
          i_src = external_link_get_index (link)
          if (external_link_is_set (int_src%source(i_src))) then
             link = int_src%source(i_src)
          else
             exit
          end if
       end do
    end if
  end function interaction_get_ultimate_source

  module subroutine interaction_exchange_mask (int)
    class(interaction_t), intent(inout) :: int
    integer :: i, index_link
    type(interaction_t), pointer :: int_link
    do i = 1, int%n_tot
       if (external_link_is_set (int%source(i))) then
          int_link => external_link_get_ptr (int%source(i))
          index_link = external_link_get_index (int%source(i))
          call interaction_merge_mask_entry &
               (int, i, int_link%mask(index_link))
          call interaction_merge_mask_entry &
               (int_link, index_link, int%mask(i))
       end if
    end do
    call int%freeze ()
  end subroutine interaction_exchange_mask

  module subroutine interaction_receive_momenta (int)
    class(interaction_t), intent(inout) :: int
    integer :: i, index_link
    type(interaction_t), pointer :: int_link
    do i = 1, int%n_tot
       if (external_link_is_set (int%source(i))) then
          int_link => external_link_get_ptr (int%source(i))
          index_link = external_link_get_index (int%source(i))
          call int%set_momentum (int_link%p(index_link), i)
       end if
    end do
  end subroutine interaction_receive_momenta

  module subroutine interaction_send_momenta (int)
    class(interaction_t), intent(in) :: int
    integer :: i, index_link
    type(interaction_t), pointer :: int_link
    do i = 1, int%n_tot
       if (external_link_is_set (int%source(i))) then
          int_link => external_link_get_ptr (int%source(i))
          index_link = external_link_get_index (int%source(i))
          call int_link%set_momentum (int%p(i), index_link)
       end if
    end do
  end subroutine interaction_send_momenta

  module subroutine interaction_pacify_momenta (int, acc)
    class(interaction_t), intent(inout) :: int
    real(default), intent(in) :: acc
    integer :: i
    do i = 1, int%n_tot
       call pacify (int%p(i), acc)
    end do
  end subroutine interaction_pacify_momenta

  module subroutine interaction_declare_subtraction (int, n_sub)
    class(interaction_t), intent(inout), target :: int
    integer, intent(in) :: n_sub
    integer :: i_sub
    type(state_iterator_t) :: it
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    type(state_matrix_t) :: state_matrix
    call state_matrix%init (store_values = .true.)
    allocate (qn (int%get_state_depth ()))
    do i_sub = 0, n_sub
       call it%init (int%state_matrix)
       do while (it%is_valid ())
          qn = it%get_quantum_numbers ()
          call qn%set_subtraction_index (i_sub)
          call state_matrix%add_state (qn, value = it%get_matrix_element ())
          call it%advance ()
       end do
    end do
    call state_matrix%freeze ()
    call state_matrix%set_n_sub ()
    call int%state_matrix%final ()
    int%state_matrix = state_matrix
  end subroutine interaction_declare_subtraction

  module subroutine find_connections (int1, int2, n, connection_index)
    class(interaction_t), intent(in) :: int1, int2
    integer, intent(out) :: n
    integer, dimension(:,:), intent(out), allocatable :: connection_index
    integer, dimension(:,:), allocatable :: conn_index_tmp
    integer, dimension(:), allocatable :: ordering
    integer :: i, j, k
    type(external_link_t) :: link1, link2
    type(interaction_t), pointer :: int_link1, int_link2
    n = 0
    do i = 1, size (int2%source)
       link2 = interaction_get_ultimate_source (int2, i)
       if (external_link_is_set (link2)) then
          int_link2 => external_link_get_ptr (link2)
          if (int_link2%tag == int1%tag) then
             n = n + 1
          else
             k = external_link_get_index (link2)
             do j = 1, size (int1%source)
                link1 = interaction_get_ultimate_source (int1, j)
                if (external_link_is_set (link1)) then
                   int_link1 => external_link_get_ptr (link1)
                   if (int_link1%tag == int_link2%tag) then
                      if (external_link_get_index (link1) == k) &
                           n = n + 1
                   end if
                end if
             end do
          end if
       end if
    end do
    allocate (conn_index_tmp (n, 2))
    n = 0
    do i = 1, size (int2%source)
       link2 = interaction_get_ultimate_source (int2, i)
       if (external_link_is_set (link2)) then
          int_link2 => external_link_get_ptr (link2)
          if (int_link2%tag == int1%tag) then
             n = n + 1
             conn_index_tmp(n,1) = external_link_get_index (int2%source(i))
             conn_index_tmp(n,2) = i
          else
             k = external_link_get_index (link2)
             do j = 1, size (int1%source)
                link1 = interaction_get_ultimate_source (int1, j)
                if (external_link_is_set (link1)) then
                   int_link1 => external_link_get_ptr (link1)
                   if (int_link1%tag == int_link2%tag) then
                      if (external_link_get_index (link1) == k) then
                         n = n + 1
                         conn_index_tmp(n,1) = j
                         conn_index_tmp(n,2) = i
                      end if
                   end if
                end if
             end do
          end if
       end if
    end do
    allocate (connection_index (n, 2))
    if (n > 1) then
       allocate (ordering (n))
       ordering = order (conn_index_tmp(:,1))
       connection_index = conn_index_tmp(ordering,:)
    else
       connection_index = conn_index_tmp
    end if
  end subroutine find_connections


end submodule interactions_s

