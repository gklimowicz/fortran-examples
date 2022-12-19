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

submodule (resonances) resonances_s

  use debug_master, only: debug_on
  use string_utils, only: str
  use format_utils, only: write_indent
  use constants, only: one
  use io_units
  use diagnostics

  implicit none

contains

  elemental module function resonance_contributors_equal &
       (c1, c2) result (equal)
    logical :: equal
    class(resonance_contributors_t), intent(in) :: c1, c2
    equal = allocated (c1%c) .and. allocated (c2%c)
    if (equal) equal = size (c1%c) == size (c2%c)
    if (equal) equal = all (c1%c == c2%c)
  end function resonance_contributors_equal

  pure module subroutine resonance_contributors_assign &
       (contributors_out, contributors_in)
    class(resonance_contributors_t), intent(inout) :: contributors_out
    class(resonance_contributors_t), intent(in) :: contributors_in
    if (allocated (contributors_out%c))  deallocate (contributors_out%c)
    if (allocated (contributors_in%c)) then
       allocate (contributors_out%c (size (contributors_in%c)))
       contributors_out%c = contributors_in%c
    end if
  end subroutine resonance_contributors_assign

  module subroutine resonance_info_copy (resonance_in, resonance_out)
    class(resonance_info_t), intent(in) :: resonance_in
    type(resonance_info_t), intent(out) :: resonance_out
    resonance_out%flavor = resonance_in%flavor
    if (allocated (resonance_in%contributors%c)) then
       associate (c => resonance_in%contributors%c)
          allocate (resonance_out%contributors%c (size (c)))
          resonance_out%contributors%c = c
       end associate
    end if
  end subroutine resonance_info_copy

  module subroutine resonance_info_write (resonance, unit, verbose)
    class(resonance_info_t), intent(in) :: resonance
    integer, optional, intent(in) :: unit
    logical, optional, intent(in) :: verbose
    integer :: u, i
    logical :: verb
    u = given_output_unit (unit);  if (u < 0)  return
    verb = .true.;  if (present (verbose))  verb = verbose
    if (verb) then
       write (u, '(A)', advance='no') "Resonance contributors: "
    else
       write (u, '(1x)', advance="no")
    end if
    if (allocated (resonance%contributors%c)) then
       do i = 1, size(resonance%contributors%c)
          write (u, '(I0,1X)', advance='no') resonance%contributors%c(i)
       end do
    else if (verb) then
       write (u, "(A)", advance="no")  "[not allocated]"
    end if
    if (resonance%flavor%is_defined ()) call resonance%flavor%write (u)
    write (u, '(A)')
  end subroutine resonance_info_write

  module subroutine resonance_info_init_pdg &
       (resonance, mom_id, pdg, model, n_out)
    class(resonance_info_t), intent(out) :: resonance
    integer, intent(in) :: mom_id
    integer, intent(in) :: pdg, n_out
    class(model_data_t), intent(in), target :: model
    type(flavor_t) :: flv
    if (debug_on) call msg_debug (D_PHASESPACE, "resonance_info_init_pdg")
    call flv%init (pdg, model)
    call resonance%init (mom_id, flv, n_out)
  end subroutine resonance_info_init_pdg

  module subroutine resonance_info_init_flv (resonance, mom_id, flv, n_out)
    class(resonance_info_t), intent(out) :: resonance
    integer, intent(in) :: mom_id
    type(flavor_t), intent(in) :: flv
    integer, intent(in) :: n_out
    integer :: i
    logical, dimension(n_out) :: contrib
    integer, dimension(n_out) :: tmp
    if (debug_on) call msg_debug (D_PHASESPACE, "resonance_info_init_flv")
    resonance%flavor = flv
    do i = 1, n_out
       tmp(i) = i
    end do
    contrib = btest (mom_id, tmp - 1)
    allocate (resonance%contributors%c (count (contrib)))
    resonance%contributors%c = pack (tmp, contrib)
  end subroutine resonance_info_init_flv

  elemental module function resonance_info_equal (r1, r2) result (equal)
    logical :: equal
    class(resonance_info_t), intent(in) :: r1, r2
    equal = r1%flavor == r2%flavor .and. r1%contributors == r2%contributors
  end function resonance_info_equal

  module function resonance_info_mapping (resonance, s) result (bw)
    real(default) :: bw
    class(resonance_info_t), intent(in) :: resonance
    real(default), intent(in) :: s
    real(default) :: m, gamma
    if (resonance%flavor%is_defined ()) then
       m = resonance%flavor%get_mass ()
       gamma = resonance%flavor%get_width ()
       bw = m**4 / ((s - m**2)**2 + gamma**2 * m**2)
    else
       bw = one
    end if
  end function resonance_info_mapping

  elemental module function resonance_info_get_n_contributors &
       (resonance) result (n)
    class(resonance_info_t), intent(in) :: resonance
    integer :: n
    if (allocated (resonance%contributors%c)) then
       n = size (resonance%contributors%c)
    else
       n = 0
    end if
  end function resonance_info_get_n_contributors

  elemental module function resonance_info_contains &
       (resonance, c) result (flag)
    class(resonance_info_t), intent(in) :: resonance
    integer, intent(in) :: c
    logical :: flag
    if (allocated (resonance%contributors%c)) then
       flag = any (resonance%contributors%c == c)
    else
       flag = .false.
    end if
  end function resonance_info_contains

  module subroutine resonance_history_clear (res_hist)
    class(resonance_history_t), intent(out) :: res_hist
  end subroutine resonance_history_clear

  module subroutine resonance_history_copy (res_hist_in, res_hist_out)
    class(resonance_history_t), intent(in) :: res_hist_in
    type(resonance_history_t), intent(out) :: res_hist_out
    integer :: i
    res_hist_out%n_resonances = res_hist_in%n_resonances
    allocate (res_hist_out%resonances (size (res_hist_in%resonances)))
    do i = 1, size (res_hist_in%resonances)
       call res_hist_in%resonances(i)%copy (res_hist_out%resonances(i))
    end do
  end subroutine resonance_history_copy

  module subroutine resonance_history_write (res_hist, unit, verbose, indent)
    class(resonance_history_t), intent(in) :: res_hist
    integer, optional, intent(in) :: unit
    logical, optional, intent(in) :: verbose
    integer, optional, intent(in) :: indent
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    call write_indent (u, indent)
    write(u, '(A,I0,A)') "Resonance history with ", &
         res_hist%n_resonances, " resonances:"
    do i = 1, res_hist%n_resonances
       call write_indent (u, indent)
       write (u, "(2x)", advance="no")
       call res_hist%resonances(i)%write (u, verbose)
    end do
  end subroutine resonance_history_write

  module subroutine resonance_history_assign (res_hist_out, res_hist_in)
    class(resonance_history_t), intent(out) :: res_hist_out
    class(resonance_history_t), intent(in) :: res_hist_in
    if (allocated (res_hist_in%resonances)) then
       res_hist_out%resonances = res_hist_in%resonances
       res_hist_out%n_resonances = res_hist_in%n_resonances
    end if
  end subroutine resonance_history_assign

  elemental module function resonance_history_equal (rh1, rh2) result (equal)
    logical :: equal
    class(resonance_history_t), intent(in) :: rh1, rh2
    integer :: i
    equal = .false.
    if (rh1%n_resonances == rh2%n_resonances) then
       do i = 1, rh1%n_resonances
          if (.not. rh1%resonances(i) == rh2%resonances(i)) then
             return
          end if
       end do
       equal = .true.
    end if
  end function resonance_history_equal

  elemental module function resonance_history_contains &
       (rh1, rh2) result (flag)
    logical :: flag
    class(resonance_history_t), intent(in) :: rh1, rh2
    integer :: i
    if (rh1%n_resonances > rh2%n_resonances) then
       flag = .true.
       do i = 1, rh2%n_resonances
          flag = flag .and. any (rh1%resonances == rh2%resonances(i))
       end do
    else
       flag = .false.
    end if
  end function resonance_history_contains

  module subroutine resonance_history_add_resonance (res_hist, resonance)
    class(resonance_history_t), intent(inout) :: res_hist
    type(resonance_info_t), intent(in) :: resonance
    type(resonance_info_t), dimension(:), allocatable :: tmp
    integer :: n, i
    if (debug_on)  call msg_debug &
         (D_PHASESPACE, "resonance_history_add_resonance")
    if (.not. allocated (res_hist%resonances)) then
       n = 0
       allocate (res_hist%resonances (1))
    else
       n = res_hist%n_resonances
       allocate (tmp (n))
       do i = 1, n
          call res_hist%resonances(i)%copy (tmp(i))
       end do
       deallocate (res_hist%resonances)
       allocate (res_hist%resonances (n+1))
       do i = 1, n
          call tmp(i)%copy (res_hist%resonances(i))
       end do
       deallocate (tmp)
    end if
    call resonance%copy (res_hist%resonances(n+1))
    res_hist%n_resonances = n + 1
    if (debug_on) call msg_debug &
         (D_PHASESPACE, "res_hist%n_resonances", res_hist%n_resonances)
  end subroutine resonance_history_add_resonance

  module subroutine resonance_history_remove_resonance (res_hist, i_res)
    class(resonance_history_t), intent(inout) :: res_hist
    integer, intent(in) :: i_res
    type(resonance_info_t), dimension(:), allocatable :: tmp_1, tmp_2
    integer :: i, j, n
    n = res_hist%n_resonances
    res_hist%n_resonances = n - 1
    if (res_hist%n_resonances == 0) then
       deallocate (res_hist%resonances)
    else
       if (i_res > 1) allocate (tmp_1(1:i_res-1))
       if (i_res < n) allocate (tmp_2(i_res+1:n))
       if (allocated (tmp_1)) then
          do i = 1, i_res - 1
             call res_hist%resonances(i)%copy (tmp_1(i))
          end do
       end if
       if (allocated (tmp_2)) then
          do i = i_res + 1, n
             call res_hist%resonances(i)%copy (tmp_2(i))
          end do
       end if
       deallocate (res_hist%resonances)
       allocate (res_hist%resonances (res_hist%n_resonances))
       j = 1
       if (allocated (tmp_1)) then
          do i = 1, i_res - 1
             call tmp_1(i)%copy (res_hist%resonances(j))
             j = j + 1
          end do
          deallocate (tmp_1)
       end if
       if (allocated (tmp_2)) then
          do i = i_res + 1, n
             call tmp_2(i)%copy (res_hist%resonances(j))
             j = j + 1
          end do
          deallocate (tmp_2)
       end if
    end if
  end subroutine resonance_history_remove_resonance

  module subroutine resonance_history_add_offset (res_hist, n)
    class(resonance_history_t), intent(inout) :: res_hist
    integer, intent(in) :: n
    integer :: i_res
    do i_res = 1, res_hist%n_resonances
       associate (contributors => res_hist%resonances(i_res)%contributors%c)
          contributors = contributors + n
       end associate
    end do
  end subroutine resonance_history_add_offset

  module function resonance_history_contains_leg &
       (res_hist, i_leg) result (val)
    logical :: val
    class(resonance_history_t), intent(in) :: res_hist
    integer, intent(in) :: i_leg
    integer :: i_res
    val = .false.
    do i_res = 1, res_hist%n_resonances
       if (any (res_hist%resonances(i_res)%contributors%c == i_leg)) then
          val = .true.
          exit
       end if
    end do
  end function resonance_history_contains_leg

  module function resonance_history_mapping &
       (res_hist, p, i_gluon) result (p_map)
    real(default) :: p_map
    class(resonance_history_t), intent(in) :: res_hist
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in), optional :: i_gluon
    integer :: i_res
    real(default) :: s
    p_map = one
    do i_res = 1, res_hist%n_resonances
       associate (res => res_hist%resonances(i_res))
          s = compute_resonance_mass (p, res%contributors%c, i_gluon)**2
          p_map = p_map * res%mapping (s)
       end associate
    end do
  end function resonance_history_mapping

  module function resonance_history_only_has_n_contributors &
       (res_hist, n) result (value)
    logical :: value
    class(resonance_history_t), intent(in) :: res_hist
    integer, intent(in) :: n
    integer :: i_res
    value = .true.
    do i_res = 1, res_hist%n_resonances
       associate (res => res_hist%resonances(i_res))
          value = value .and. size (res%contributors%c) == n
       end associate
    end do
  end function resonance_history_only_has_n_contributors

  module function resonance_history_has_flavor &
       (res_hist, flv) result (has_flv)
    logical :: has_flv
    class(resonance_history_t), intent(in) :: res_hist
    type(flavor_t), intent(in) :: flv
    integer :: i
    has_flv = .false.
    do i = 1, res_hist%n_resonances
       has_flv = has_flv .or. res_hist%resonances(i)%flavor == flv
    end do
  end function resonance_history_has_flavor

  module subroutine resonance_info_evaluate_distance (res_info, p, dist)
    class(resonance_info_t), intent(in) :: res_info
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(out) :: dist
    real(default) :: m, w
    type(vector4_t) :: q
    m = res_info%flavor%get_mass ()
    w = res_info%flavor%get_width ()
    q = sum (p(res_info%contributors%c))
    dist = abs (q**2 - m**2) / (m * w)
  end subroutine resonance_info_evaluate_distance

  module subroutine resonance_history_evaluate_distances (res_hist, p, dist)
    class(resonance_history_t), intent(in) :: res_hist
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), dimension(:), intent(out) :: dist
    integer :: i
    do i = 1, res_hist%n_resonances
       call res_hist%resonances(i)%evaluate_distance (p, dist(i))
    end do
  end subroutine resonance_history_evaluate_distances

  module function resonance_info_evaluate_gaussian &
       (res_info, p, gw) result (factor)
    class(resonance_info_t), intent(in) :: res_info
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in) :: gw
    real(default) :: factor
    real(default) :: dist, w
    if (gw > 0) then
       w = res_info%flavor%get_width ()
       call res_info%evaluate_distance (p, dist)
       factor = exp (- (dist / (gw * w)) **2)
    else
       factor = 1
    end if
  end function resonance_info_evaluate_gaussian

  module function resonance_history_evaluate_gaussian &
       (res_hist, p, gw) result (factor)
    class(resonance_history_t), intent(in) :: res_hist
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in) :: gw
    real(default), dimension(:), allocatable :: dist
    real(default) :: factor
    integer :: i
    factor = 1
    do i = 1, res_hist%n_resonances
       factor = factor * res_hist%resonances(i)%evaluate_gaussian (p, gw)
    end do
  end function resonance_history_evaluate_gaussian

  module function resonance_info_is_on_shell (res_info, p, on_shell_limit) &
       result (flag)
    class(resonance_info_t), intent(in) :: res_info
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in) :: on_shell_limit
    logical :: flag
    real(default) :: dist
    call res_info%evaluate_distance (p, dist)
    flag = dist < on_shell_limit
  end function resonance_info_is_on_shell

  module function resonance_history_is_on_shell &
       (res_hist, p, on_shell_limit) result (flag)
    class(resonance_history_t), intent(in) :: res_hist
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in) :: on_shell_limit
    logical :: flag
    integer :: i
    flag = .true.
    do i = 1, res_hist%n_resonances
       flag = flag .and. res_hist%resonances(i)%is_on_shell (p, on_shell_limit)
    end do
  end function resonance_history_is_on_shell

  module function resonance_info_as_omega_string &
       (res_info, n_in) result (string)
    class(resonance_info_t), intent(in) :: res_info
    integer, intent(in) :: n_in
    type(string_t) :: string
    integer :: i
    string = ""
    if (allocated (res_info%contributors%c)) then
       do i = 1, size (res_info%contributors%c)
          if (i > 1)  string = string // "+"
          string = string // str (res_info%contributors%c(i) + n_in)
       end do
       string = string // "~" // res_info%flavor%get_name ()
    end if
  end function resonance_info_as_omega_string

  module function resonance_history_as_omega_string &
       (res_hist, n_in) result (string)
    class(resonance_history_t), intent(in) :: res_hist
    integer, intent(in) :: n_in
    type(string_t) :: string
    integer :: i
    string = ""
    do i = 1, res_hist%n_resonances
       if (i > 1)  string = string // " && "
       string = string // res_hist%resonances(i)%as_omega_string (n_in)
    end do
  end function resonance_history_as_omega_string

  module subroutine resonance_tree_write (tree, unit, indent)
    class(resonance_tree_t), intent(in) :: tree
    integer, intent(in), optional :: unit, indent
    integer :: u, b, c
    u = given_output_unit (unit)
    call write_indent (u, indent)
    write (u, "(A)", advance="no")  "Resonance tree:"
    if (tree%n > 0) then
       write (u, *)
       do b = 1, tree%n
          call write_indent (u, indent)
          write (u, "(2x,'r',I0,':',1x)", advance="no")  b
          associate (branch => tree%branch(b))
            call branch%flv%write (u)
            write (u, "(1x,'=>')", advance="no")
            if (allocated (branch%r_child)) then
               do c = 1, size (branch%r_child)
                  write (u, "(1x,'r',I0)", advance="no")  branch%r_child(c)
               end do
            end if
            if (allocated (branch%o_child)) then
               do c = 1, size (branch%o_child)
                  write (u, "(1x,I0)", advance="no")  branch%o_child(c)
               end do
            end if
            write (u, *)
          end associate
       end do
    else
       write (u, "(1x,A)")  "[empty]"
    end if
  end subroutine resonance_tree_write

  module function resonance_tree_get_n_resonances (tree) result (n)
    class(resonance_tree_t), intent(in) :: tree
    integer :: n
    n = tree%n
  end function resonance_tree_get_n_resonances

  module function resonance_tree_get_flv (tree, i) result (flv)
    class(resonance_tree_t), intent(in) :: tree
    integer, intent(in) :: i
    type(flavor_t) :: flv
    flv = tree%branch(i)%flv
  end function resonance_tree_get_flv

  module function resonance_tree_get_children (tree, i, offset_r, offset_o) &
       result (child)
    class(resonance_tree_t), intent(in) :: tree
    integer, intent(in) :: i, offset_r, offset_o
    integer, dimension(:), allocatable :: child
    integer :: nr, no
    associate (branch => tree%branch(i))
      nr = size (branch%r_child)
      no = size (branch%o_child)
      allocate (child (nr + no))
      child(1:nr) = branch%r_child + offset_r
      child(nr+1:nr+no) = branch%o_child + offset_o
    end associate
  end function resonance_tree_get_children

  module subroutine resonance_history_to_tree (res_hist, tree)
    class(resonance_history_t), intent(in) :: res_hist
    type(resonance_tree_t), intent(out) :: tree
    integer :: nr
    integer, dimension(:), allocatable :: r_branch, r_source

    nr = res_hist%n_resonances
    tree%n = nr
    allocate (tree%branch (tree%n), r_branch (tree%n), r_source (tree%n))
    if (tree%n > 0) then
       call find_branch_ordering ()
       call set_flavors ()
       call set_child_resonances ()
       call set_child_outgoing ()
    end if

  contains

    subroutine find_branch_ordering ()
      integer, dimension(:), allocatable :: nc_array
      integer :: r, ir, nc
      allocate (nc_array (tree%n))
      nc_array(:) = res_hist%resonances%get_n_contributors ()
      ir = 0
      do nc = maxval (nc_array), minval (nc_array), -1
         do r = 1, nr
            if (nc_array(r) == nc) then
               ir = ir + 1
               r_branch(r) = ir
               r_source(ir) = r
            end if
         end do
      end do
    end subroutine find_branch_ordering

    subroutine set_flavors ()
      integer :: r
      do r = 1, nr
         tree%branch(r_branch(r))%flv = res_hist%resonances(r)%flavor
      end do
    end subroutine set_flavors

    subroutine set_child_resonances ()
      integer, dimension(:), allocatable :: r_child, r_parent
      integer :: r, ir, pr
      allocate (r_parent (nr), source = 0)
      SCAN_RES: do r = 1, nr
         associate (this_res => res_hist%resonances(r))
           SCAN_PARENT: do ir = 1, nr
              pr = r_source(ir)
              if (pr == r)  cycle SCAN_PARENT
              if (all (res_hist%resonances(pr)%contains &
                   (this_res%contributors%c))) then
                 r_parent (r) = pr
              end if
           end do SCAN_PARENT
         end associate
      end do SCAN_RES
      allocate (r_child (nr), source = [(r, r = 1, nr)])
      do r = 1, nr
         ir = r_branch(r)
         tree%branch(ir)%r_child = r_branch (pack (r_child, r_parent == r))
      end do
    end subroutine set_child_resonances

    subroutine set_child_outgoing ()
      integer, dimension(:), allocatable :: o_child, o_parent
      integer :: o_max, r, o, ir
      o_max = 0
      do r = 1, nr
         associate (this_res => res_hist%resonances(r))
           o_max = max (o_max, maxval (this_res%contributors%c))
         end associate
      end do
      allocate (o_parent (o_max), source=0)
      SCAN_OUT: do o = 1, o_max
         SCAN_PARENT: do ir = 1, nr
            r = r_source(ir)
            associate (this_res => res_hist%resonances(r))
              if (this_res%contains (o))  o_parent(o) = r
            end associate
         end do SCAN_PARENT
      end do SCAN_OUT
      allocate (o_child (o_max), source = [(o, o = 1, o_max)])
      do r = 1, nr
         ir = r_branch(r)
         tree%branch(ir)%o_child = pack (o_child, o_parent == r)
      end do
    end subroutine set_child_outgoing

  end subroutine resonance_history_to_tree

  module subroutine resonance_history_set_write &
       (res_set, unit, indent, show_trees)
    class(resonance_history_set_t), intent(in) :: res_set
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: indent
    logical, intent(in), optional :: show_trees
    logical :: s_trees
    integer :: u, i, j, ind
    u = given_output_unit (unit)
    s_trees = .false.;  if (present (show_trees))  s_trees = show_trees
    ind = 0;  if (present (indent))  ind = indent
    call write_indent (u, indent)
    write (u, "(A)", advance="no")  "Resonance history set:"
    if (res_set%complete) then
       write (u, *)
    else
       write (u, "(1x,A)")  "[incomplete]"
    end if
    do i = 1, res_set%last
       write (u, "(1x,I0,1x)", advance="no")  i
       call res_set%history(i)%write (u, verbose=.false., indent=indent)
       if (allocated (res_set%contains_this)) then
          call write_indent (u, indent)
          write (u, "(3x,A)", advance="no")  "contained in ("
           do j = 1, size (res_set%contains_this(i)%i)
             if (j>1)  write (u, "(',')", advance="no")
             write (u, "(I0)", advance="no")  res_set%contains_this(i)%i(j)
          end do
          write (u, "(A)")  ")"
       end if
       if (s_trees .and. allocated (res_set%tree)) then
          call res_set%tree(i)%write (u, ind + 1)
       end if
    end do
  end subroutine resonance_history_set_write

  module subroutine resonance_history_set_init &
       (res_set, n_filter, initial_size)
    class(resonance_history_set_t), intent(out) :: res_set
    integer, intent(in), optional :: n_filter
    integer, intent(in), optional :: initial_size
    if (present (n_filter))  res_set%n_filter = n_filter
    if (present (initial_size)) then
       allocate (res_set%history (initial_size))
    else
       allocate (res_set%history (resonance_history_set_initial_size))
    end if
  end subroutine resonance_history_set_init

  module subroutine resonance_history_set_enter &
       (res_set, res_history, trivial)
    class(resonance_history_set_t), intent(inout) :: res_set
    type(resonance_history_t), intent(in) :: res_history
    logical, intent(in), optional :: trivial
    integer :: i, new
    if (res_history%n_resonances == 0) then
       if (present (trivial)) then
          if (.not. trivial) return
       else
          return
       end if
    end if
    if (res_set%n_filter > 0) then
       if (.not. res_history%only_has_n_contributors (res_set%n_filter))  return
    end if
    do i = 1, res_set%last
       if (res_set%history(i) == res_history)  return
    end do
    new = res_set%last + 1
    if (new > size (res_set%history))  call res_set%expand ()
    res_set%history(new) = res_history
    res_set%last = new
  end subroutine resonance_history_set_enter

  module subroutine resonance_history_set_freeze (res_set)
    class(resonance_history_set_t), intent(inout) :: res_set
    integer :: i, n, c
    logical, dimension(:), allocatable :: contains_this
    integer, dimension(:), allocatable :: index_array
    n = res_set%last
    allocate (contains_this (n))
    allocate (index_array (n), source = [(i, i=1, n)])
    allocate (res_set%contains_this (n))
    do i = 1, n
       contains_this = resonance_history_contains &
            (res_set%history(1:n), res_set%history(i))
       c = count (contains_this)
       allocate (res_set%contains_this(i)%i (c))
       res_set%contains_this(i)%i = pack (index_array, contains_this)
    end do
    allocate (res_set%tree (n))
    do i = 1, n
       call res_set%history(i)%to_tree (res_set%tree(i))
    end do
    res_set%complete = .true.
  end subroutine resonance_history_set_freeze

  module subroutine resonance_history_set_determine_on_shell_histories &
       (res_set, p, on_shell_limit, index_array)
    class(resonance_history_set_t), intent(in) :: res_set
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in) :: on_shell_limit
    integer, dimension(:), allocatable, intent(out) :: index_array
    integer :: n, i
    integer, dimension(:), allocatable :: i_array
    if (res_set%complete) then
       n = res_set%last
       allocate (i_array (n), source=0)
       do i = 1, n
          if (res_set%history(i)%is_on_shell (p, on_shell_limit))  i_array(i) = i
       end do
       do i = 1, n
          if (any (i_array(res_set%contains_this(i)%i) /= 0)) then
             i_array(i) = 0
          end if
       end do
       allocate (index_array (count (i_array /= 0)))
       index_array(:) = pack (i_array, i_array /= 0)
    end if
  end subroutine resonance_history_set_determine_on_shell_histories

  module function resonance_history_set_evaluate_gaussian &
       (res_set, p, gw, i) result (factor)
    class(resonance_history_set_t), intent(in) :: res_set
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in) :: gw
    integer, intent(in) :: i
    real(default) :: factor
    factor = res_set%history(i)%evaluate_gaussian (p, gw)
  end function resonance_history_set_evaluate_gaussian

  module function resonance_history_set_get_n_history (res_set) result (n)
    class(resonance_history_set_t), intent(in) :: res_set
    integer :: n
    if (res_set%complete) then
       n = res_set%last
    else
       n = 0
    end if
  end function resonance_history_set_get_n_history

  module function resonance_history_set_get_history &
       (res_set, i) result (res_history)
    class(resonance_history_set_t), intent(in) :: res_set
    integer, intent(in) :: i
    type(resonance_history_t) :: res_history
    if (res_set%complete .and. i <= res_set%last) then
       res_history = res_set%history(i)
    end if
  end function resonance_history_set_get_history

  module subroutine resonance_history_set_to_array (res_set, res_history)
    class(resonance_history_set_t), intent(in) :: res_set
    type(resonance_history_t), dimension(:), allocatable, intent(out) :: &
         res_history
    if (res_set%complete) then
       allocate (res_history (res_set%last))
       res_history(:) = res_set%history(1:res_set%last)
    end if
  end subroutine resonance_history_set_to_array

  module subroutine resonance_history_set_get_tree (res_set, i, res_tree)
    class(resonance_history_set_t), intent(in) :: res_set
    integer, intent(in) :: i
    type(resonance_tree_t), intent(out) :: res_tree
    if (res_set%complete) then
       res_tree = res_set%tree(i)
    end if
  end subroutine resonance_history_set_get_tree

  module subroutine resonance_history_set_expand (res_set)
    class(resonance_history_set_t), intent(inout) :: res_set
    type(resonance_history_t), dimension(:), allocatable :: history_new
    integer :: s
    s = size (res_set%history)
    allocate (history_new (2 * s))
    history_new(1:s) = res_set%history(1:s)
    call move_alloc (history_new, res_set%history)
  end subroutine resonance_history_set_expand


end submodule resonances_s

