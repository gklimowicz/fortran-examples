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

submodule (phs_trees) phs_trees_s

  use io_units
  use constants, only: twopi, twopi2, twopi5
  use format_defs, only: FMT_19
  use numeric_utils, only: vanishes
  use diagnostics

  implicit none

contains

  elemental module subroutine phs_prt_set_defined (prt)
    class(phs_prt_t), intent(inout) :: prt
    prt%defined = .true.
  end subroutine phs_prt_set_defined

  elemental module subroutine phs_prt_set_undefined (prt)
    class(phs_prt_t), intent(inout) :: prt
    prt%defined = .false.
  end subroutine phs_prt_set_undefined

  elemental module subroutine phs_prt_set_momentum (prt, p)
    class(phs_prt_t), intent(inout) :: prt
    type(vector4_t), intent(in) :: p
    prt%p = p
  end subroutine phs_prt_set_momentum

  elemental module subroutine phs_prt_set_msq (prt, p2)
    class(phs_prt_t), intent(inout) :: prt
    real(default), intent(in) :: p2
    prt%p2 = p2
  end subroutine phs_prt_set_msq

  elemental module function phs_prt_is_defined (prt) result (defined)
    logical :: defined
    class(phs_prt_t), intent(in) :: prt
    defined = prt%defined
  end function phs_prt_is_defined

  elemental module function phs_prt_get_momentum (prt) result (p)
    type(vector4_t) :: p
    class(phs_prt_t), intent(in) :: prt
    p = prt%p
  end function phs_prt_get_momentum

  elemental module function phs_prt_get_msq (prt) result (p2)
    real(default) :: p2
    class(phs_prt_t), intent(in) :: prt
    p2 = prt%p2
  end function phs_prt_get_msq

  elemental module subroutine phs_prt_combine (prt, prt1, prt2)
    class(phs_prt_t), intent(inout) :: prt
    type(phs_prt_t), intent(in) :: prt1, prt2
    prt%defined = .true.
    prt%p = prt1%p + prt2%p
    prt%p2 = prt%p ** 2
    call phs_prt_check (prt)
  end subroutine phs_prt_combine

  module subroutine phs_prt_write (prt, unit)
    class(phs_prt_t), intent(in) :: prt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    if (prt%defined) then
       call vector4_write (prt%p, u)
       write (u, "(1x,A,1x," // FMT_19 // ")") "T = ", prt%p2
    else
       write (u, "(3x,A)") "[undefined]"
    end if
  end subroutine phs_prt_write

  elemental module subroutine phs_prt_check (prt)
    class(phs_prt_t), intent(inout) :: prt
    if (prt%p2 < 0._default) then
       prt%p2 = 0._default
    end if
  end subroutine phs_prt_check

  elemental module subroutine phs_tree_init &
       (tree, n_in, n_out, n_masses, n_angles)
    class(phs_tree_t), intent(inout) :: tree
    integer, intent(in) :: n_in, n_out, n_masses, n_angles
    integer(TC) :: i
    tree%n_externals = n_in + n_out
    tree%n_branches_tot = 2**(n_in+n_out) - 1
    tree%n_branches_out = 2**n_out - 1
    tree%mask = 0
    do i = 0, n_in + n_out - 1
       tree%mask = ibset (tree%mask, i)
    end do
    tree%n_in = n_in
    tree%mask_in = 0
    do i = n_out, n_in + n_out - 1
       tree%mask_in = ibset (tree%mask_in, i)
    end do
    tree%mask_out = ieor (tree%mask, tree%mask_in)
    tree%n_msq = n_masses
    tree%n_angles = n_angles
    allocate (tree%branch (tree%n_branches_tot))
    tree%n_branches  = 0
    allocate (tree%mapping (tree%n_branches_out))
    allocate (tree%mass_sum (tree%n_branches_out))
    allocate (tree%effective_mass (tree%n_branches_out))
    allocate (tree%effective_width (tree%n_branches_out))
  end subroutine phs_tree_init

  elemental module subroutine phs_tree_final (tree)
    class(phs_tree_t), intent(inout) :: tree
    deallocate (tree%branch)
    deallocate (tree%mapping)
    deallocate (tree%mass_sum)
    deallocate (tree%effective_mass)
    deallocate (tree%effective_width)
  end subroutine phs_tree_final

  module subroutine phs_tree_write (tree, unit)
    class(phs_tree_t), intent(in) :: tree
    integer, intent(in), optional :: unit
    integer :: u
    integer(TC) :: k
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, '(3X,A,1x,I0,5X,A,I3)') &
         'External:', tree%n_externals, 'Mask:', tree%mask
    write (u, '(3X,A,1x,I0,5X,A,I3)') &
         'Incoming:', tree%n_in, 'Mask:', tree%mask_in
    write (u, '(3X,A,1x,I0,5X,A,I3)') &
         'Branches:', tree%n_branches
    do k = size (tree%branch), 1, -1
       if (tree%branch(k)%set) &
            call phs_branch_write (tree%branch(k), unit=unit, kval=k)
    end do
    do k = 1, size (tree%mapping)
       call tree%mapping (k)%write (unit, verbose=.true.)
    end do
    write (u, "(3x,A)") "Arrays: mass_sum, effective_mass, effective_width"
    do k = 1, size (tree%mass_sum)
       if (tree%branch(k)%set) then
          write (u, "(5x,I0,3(2x," // FMT_19 // "))") k, tree%mass_sum(k), &
               tree%effective_mass(k), tree%effective_width(k)
       end if
    end do
  end subroutine phs_tree_write

  subroutine phs_branch_write (b, unit, kval)
    type(phs_branch_t), intent(in) :: b
    integer, intent(in), optional :: unit
    integer(TC), intent(in), optional :: kval
    integer :: u
    integer(TC) :: k
    character(len=6) :: tmp
    character(len=1) :: firstborn(2), sign_decay, sign_axis
    integer :: i
    u = given_output_unit (unit);  if (u < 0)  return
    k = 0;  if (present (kval))  k = kval
    if (b%origin /= 0) then
       write(tmp, '(A,I4,A)') '(', b%origin, ')'
    else
       tmp = ' '
    end if
    do i=1, 2
       if (b%firstborn == i) then
          firstborn(i) = "*"
       else
          firstborn(i) = " "
       end if
    end do
    if (b%inverted_decay) then
       sign_decay = "-"
    else
       sign_decay = "+"
    end if
    if (b%inverted_axis) then
       sign_axis = "-"
    else
       sign_axis = "+"
    end if
    if (b%has_children) then
       if (b%has_friend) then
          write(u,'(4X,A1,I0,3x,A,1X,A,I0,A1,1x,I0,A1,1X,A1,1X,A,1x,I0)') &
               &   '*', k, tmp, &
               &   'Daughters: ', &
               &   b%daughter(1), firstborn(1), &
               &   b%daughter(2), firstborn(2), sign_decay, &
               &   'Friend:    ', b%friend
       else
          write(u,'(4X,A1,I0,3x,A,1X,A,I0,A1,1x,I0,A1,1X,A1,1X,A)') &
               &   '*', k, tmp, &
               &   'Daughters: ', &
               &   b%daughter(1), firstborn(1), &
               &   b%daughter(2), firstborn(2), sign_decay, &
               &   '(axis '//sign_axis//')'
       end if
    else
       write(u,'(5X,I0)') k
    end if
  end subroutine phs_branch_write

  module subroutine phs_tree_from_array (tree, a)
    class(phs_tree_t), intent(inout) :: tree
    integer(TC), dimension(:), intent(in) :: a
    integer :: i
    integer(TC) :: k
    do i=1, size(a)
       k = a(i)
       if (iand(k, tree%mask_in) == tree%mask_in)  k = ieor(tree%mask, k)
       tree%branch(k)%set = .true.
       tree%n_branches = tree%n_branches+1
    end do
    do i=0, tree%n_externals-1
       k = ibset(0,i)
       if (iand(k, tree%mask_in) == tree%mask_in)  k = ieor(tree%mask, k)
       if (tree%branch(ieor(tree%mask, k))%set) then
          tree%branch(ieor(tree%mask, k))%set = .false.
          tree%branch(k)%set = .true.
       else if (.not.tree%branch(k)%set) then
          tree%branch(k)%set = .true.
          tree%n_branches = tree%n_branches+1
       end if
    end do
    if (tree%n_branches /= tree%n_externals*2-3) then
       call phs_tree_write (tree)
       call msg_bug &
            & (" Wrong number of branches set in phase space tree")
    end if
    do k=1, size (tree%branch)
       if (tree%branch(k)%set .and. tc_decay_level (k) /= 1) then
          call branch_set_relatives(k)
       end if
    end do
  contains
    subroutine branch_set_relatives (k)
      integer(TC), intent(in) :: k
      integer(TC) :: m,n
      do m=1, k-1
         if (iand(k,m)==m) then
            n = ieor(k,m)
            if ( tree%branch(m)%set .and. tree%branch(n)%set ) then
               tree%branch(k)%daughter(1) = m;  tree%branch(k)%daughter(2) = n
               tree%branch(m)%mother      = k;  tree%branch(n)%mother      = k
               tree%branch(m)%sibling     = n;  tree%branch(n)%sibling     = m
               tree%branch(k)%has_children = .true.
               return
            end if
         end if
      end do
      call phs_tree_write (tree)
      call msg_bug &
           & (" Missing daughter branch(es) in phase space tree")
    end subroutine branch_set_relatives

  end subroutine phs_tree_from_array

  module subroutine phs_tree_flip_t_to_s_channel (tree)
    class(phs_tree_t), intent(inout) :: tree
    integer(TC) :: k, f, m, n, d, s
    if (tree%n_in == 2) then
       FLIP: do k=3, tree%mask-1
          if (.not. tree%branch(k)%set) cycle FLIP
          f = iand(k,tree%mask_in)
          if (f==0 .or. f==k) cycle FLIP
          m = tree%branch(k)%mother
          s = tree%branch(k)%sibling
          if (s==0) call find_orphan(s)
          d = tree%branch(k)%daughter(1)
          n = ior(d,s)
          tree%branch(k)%set = .false.
          tree%branch(n)%set = .true.
          tree%branch(n)%origin = k
          tree%branch(n)%daughter(1) = d; tree%branch(d)%mother  = n
          tree%branch(n)%daughter(2) = s; tree%branch(s)%mother  = n
          tree%branch(n)%has_children = .true.
          tree%branch(d)%sibling = s;  tree%branch(s)%sibling = d
          tree%branch(n)%sibling = f;  tree%branch(f)%sibling = n
          tree%branch(n)%mother      = m
          tree%branch(f)%mother      = m
          if (m/=0) then
             tree%branch(m)%daughter(1) = n
             tree%branch(m)%daughter(2) = f
          end if
          tree%branch(n)%friend = f
          tree%branch(n)%has_friend = .true.
          tree%branch(n)%firstborn = 2
       end do FLIP
    end if
  contains
    subroutine find_orphan(s)
      integer(TC) :: s
      do s=1, tree%mask_out
         if (tree%branch(s)%set .and. tree%branch(s)%mother==0) return
      end do
      call phs_tree_write (tree)
      call msg_bug (" Can't flip phase space tree to channel")
    end subroutine find_orphan
  end subroutine phs_tree_flip_t_to_s_channel

  function tc_flipped (tree, kt) result (ks)
    type(phs_tree_t), intent(in) :: tree
    integer(TC), intent(in) :: kt
    integer(TC) :: ks
    if (iand (kt, tree%mask_in) == 0) then
       ks = kt
    else
       ks = tree%branch(iand (kt, tree%mask_out))%mother
    end if
  end function tc_flipped

  module subroutine phs_tree_canonicalize (tree)
    class(phs_tree_t), intent(inout) :: tree
    integer :: n_out
    integer(TC) :: k_out
    call branch_canonicalize (tree%branch(tree%mask_out))
    n_out = tree%n_externals - tree%n_in
    k_out = tree%mask_out
    if (tree%branch(k_out)%has_friend &
         & .and. tree%branch(k_out)%friend == ibset (0, n_out)) then
       tree%branch(k_out)%inverted_axis = .not.tree%branch(k_out)%inverted_axis
    end if
    tree%branch(k_out)%has_friend = .false.
    tree%branch(k_out)%friend = 0
  contains
    recursive subroutine branch_canonicalize (b)
      type(phs_branch_t), intent(inout) :: b
      integer(TC) :: d1, d2
      if (b%has_children) then
         d1 = b%daughter(1)
         d2 = b%daughter(2)
         if (d1 > d2) then
            b%daughter(1) = d2
            b%daughter(2) = d1
            b%inverted_decay = .not.b%inverted_decay
            if (b%firstborn /= 0)  b%firstborn = 3 - b%firstborn
         end if
         call branch_canonicalize (tree%branch(b%daughter(1)))
         call branch_canonicalize (tree%branch(b%daughter(2)))
      end if
    end subroutine branch_canonicalize
  end subroutine phs_tree_canonicalize

  module subroutine phs_tree_init_mapping (tree, k, type, pdg, model)
    class(phs_tree_t), intent(inout) :: tree
    integer(TC), intent(in) :: k
    type(string_t), intent(in) :: type
    integer, intent(in) :: pdg
    class(model_data_t), intent(in), target :: model
    integer(TC) :: kk
    kk = tc_flipped (tree, k)
    call tree%mapping(kk)%init (kk, type, pdg, model)
  end subroutine phs_tree_init_mapping

  module subroutine phs_tree_set_mapping_parameters &
       (tree, mapping_defaults, variable_limits)
    class(phs_tree_t), intent(inout) :: tree
    type(mapping_defaults_t), intent(in) :: mapping_defaults
    logical, intent(in) :: variable_limits
    integer(TC) :: k
    do k = 1, tree%n_branches_out
       call tree%mapping(k)%set_parameters (mapping_defaults, variable_limits)
    end do
  end subroutine phs_tree_set_mapping_parameters

  module subroutine phs_tree_assign_s_mapping (tree, mapping)
    class(phs_tree_t), intent(in) :: tree
    type(mapping_t), intent(out) :: mapping
    mapping = tree%mapping(tree%mask_out)
  end subroutine phs_tree_assign_s_mapping

  module subroutine phs_tree_set_mass_sum (tree, flv)
    class(phs_tree_t), intent(inout) :: tree
    type(flavor_t), dimension(:), intent(in) :: flv
    integer(TC) :: k
    integer :: i
    tree%mass_sum = 0
    do k = 1, tree%n_branches_out
       do i = 0, size (flv) - 1
          if (btest(k,i)) then
             if (ibclr(k,i) == 0) then
                tree%mass_sum(k) = flv(i+1)%get_mass ()
             else
                tree%mass_sum(k) = &
                     tree%mass_sum(ibclr(k,i)) + tree%mass_sum(ibset(0,i))
             end if
             exit
          end if
       end do
    end do
  end subroutine phs_tree_set_mass_sum

  module subroutine phs_tree_set_effective_masses (tree)
    class(phs_tree_t), intent(inout) :: tree
    tree%effective_mass = 0
    tree%effective_width = 0
    call set_masses_x (tree%mask_out)
  contains
    recursive subroutine set_masses_x (k)
      integer(TC), intent(in) :: k
      integer(TC) :: k1, k2
      if (tree%branch(k)%has_children) then
         k1 = tree%branch(k)%daughter(1)
         k2 = tree%branch(k)%daughter(2)
         call set_masses_x (k1)
         call set_masses_x (k2)
         if (tree%mapping(k)%is_s_channel ()) then
            tree%effective_mass(k) = tree%mapping(k)%get_mass ()
            tree%effective_width(k) = tree%mapping(k)%get_width ()
         else
            tree%effective_mass(k) = &
                 tree%effective_mass(k1) + tree%effective_mass(k2)
            tree%effective_width(k) = &
                 tree%effective_width(k1) + tree%effective_width(k2)
         end if
      else
         tree%effective_mass(k) = tree%mass_sum(k)
      end if
    end subroutine set_masses_x
  end subroutine phs_tree_set_effective_masses

  module subroutine phs_tree_set_step_mappings &
       (tree, exp_type, variable_limits)
    class(phs_tree_t), intent(inout) :: tree
    logical, intent(in) :: exp_type
    logical, intent(in) :: variable_limits
    type(string_t) :: map_str
    integer(TC) :: k
    if (exp_type) then
       map_str = "step_exp"
    else
       map_str = "step_hyp"
    end if
    k = tree%mask_out
    call set_step_mappings_x (k, 0._default, 0._default)
  contains
    recursive subroutine set_step_mappings_x (k, m_limit, w_limit)
      integer(TC), intent(in) :: k
      real(default), intent(in) :: m_limit, w_limit
      integer(TC), dimension(2) :: kk
      real(default), dimension(2) :: m, w
      if (tree%branch(k)%has_children) then
         if (m_limit > 0) then
            if (.not. tree%mapping(k)%is_set ()) then
               call tree%mapping(k)%init (k, map_str)
               call tree%mapping(k)%set_step_mapping_parameters (m_limit, &
                    w_limit, variable_limits)
            end if
         end if
         kk = tree%branch(k)%daughter
         m = tree%effective_mass(kk)
         w = tree%effective_width(kk)
         if (tree%mapping(k)%is_s_channel ()) then
            call set_step_mappings_x (kk(1), &
                 tree%mapping(k)%get_mass () - m(2), &
                 tree%mapping(k)%get_width () + w(2))
            call set_step_mappings_x (kk(2), &
                 tree%mapping(k)%get_mass () - m(1), &
                 tree%mapping(k)%get_width () + w(1))
         else if (m_limit > 0) then
            call set_step_mappings_x (kk(1), &
                 m_limit - m(2), &
                 w_limit + w(2))
            call set_step_mappings_x (kk(2), &
                 m_limit - m(1), &
                 w_limit + w(1))
         else
            call set_step_mappings_x (kk(1), &
                 - m(2), &
                 + w(2))
            call set_step_mappings_x (kk(2), &
                 - m(1), &
                 + w(1))
         end if
      end if
    end subroutine set_step_mappings_x
  end subroutine phs_tree_set_step_mappings

  module subroutine phs_tree_extract_resonance_history (tree, res_history)
    class(phs_tree_t), intent(in) :: tree
    type(resonance_history_t), intent(out) :: res_history
    type(resonance_info_t) :: res_info
    integer :: i
    if (allocated (tree%mapping)) then
       do i = 1, size (tree%mapping)
          associate (mapping => tree%mapping(i))
             if (mapping%is_s_channel ()) then
                call res_info%init (mapping%get_bincode (), mapping%get_flv (), &
                     n_out = tree%n_externals - tree%n_in)
                call res_history%add_resonance (res_info)
             end if
          end associate
       end do
    end if
  end subroutine phs_tree_extract_resonance_history

  module function phs_tree_equivalent (t1, t2, perm) result (is_equal)
    type(phs_tree_t), intent(in) :: t1, t2
    type(permutation_t), intent(in) :: perm
    logical :: equal, is_equal
    integer(TC) :: k1, k2, mask_in
    k1 = t1%mask_out
    k2 = t2%mask_out
    mask_in = t1%mask_in
    equal = .true.
    call check (t1%branch(k1), t2%branch(k2), k1, k2)
    is_equal = equal
  contains
    recursive subroutine check (b1, b2, k1, k2)
      type(phs_branch_t), intent(in) :: b1, b2
      integer(TC), intent(in) :: k1, k2
      integer(TC), dimension(2) :: d1, d2, pd2
      integer :: i
      if (.not.b1%has_friend .and. .not.b2%has_friend) then
         equal = .true.
      else if (b1%has_friend .and. b2%has_friend) then
         equal = (b1%friend == tc_permute (b2%friend, perm, mask_in))
      end if
      if (equal) then
         if (b1%has_children .and. b2%has_children) then
            d1 = b1%daughter
            d2 = b2%daughter
            do i=1, 2
               pd2(i) = tc_permute (d2(i), perm, mask_in)
            end do
            if (d1(1)==pd2(1) .and. d1(2)==pd2(2)) then
               equal = (b1%firstborn == b2%firstborn)
               if (equal) call check &
                    &     (t1%branch(d1(1)), t2%branch(d2(1)), d1(1), d2(1))
               if (equal) call check &
                    &     (t1%branch(d1(2)), t2%branch(d2(2)), d1(2), d2(2))
            else if (d1(1)==pd2(2) .and. d1(2)==pd2(1)) then
               equal = ( (b1%firstborn == 0 .and. b2%firstborn == 0) &
                    &   .or. (b1%firstborn == 3 - b2%firstborn) )
               if (equal) call check &
                    &     (t1%branch(d1(1)), t2%branch(d2(2)), d1(1), d2(2))
               if (equal) call check &
                    &     (t1%branch(d1(2)), t2%branch(d2(1)), d1(2), d2(1))
            else
               equal = .false.
            end if
         end if
      end if
      if (equal) then
         equal = (t1%mapping(k1) == t2%mapping(k2))
      end if
    end subroutine check
  end function phs_tree_equivalent

  module subroutine phs_tree_find_msq_permutation &
       (tree1, tree2, perm2, msq_perm)
    type(phs_tree_t), intent(in) :: tree1, tree2
    type(permutation_t), intent(in) :: perm2
    type(permutation_t), intent(out) :: msq_perm
    type(permutation_t) :: perm1
    integer(TC) :: mask_in, root
    integer(TC), dimension(:), allocatable :: index1, index2
    integer :: i
    allocate (index1 (tree1%n_msq), index2 (tree2%n_msq))
    call permutation_init (perm1, permutation_size (perm2))
    mask_in = tree1%mask_in
    root = tree1%mask_out
    i = 0
    call tree_scan (tree1, root, perm1, index1)
    i = 0
    call tree_scan (tree2, root, perm2, index2)
    call permutation_find (msq_perm, index1, index2)
  contains
    recursive subroutine tree_scan (tree, k, perm, index)
      type(phs_tree_t), intent(in) :: tree
      integer(TC), intent(in) :: k
      type(permutation_t), intent(in) :: perm
      integer, dimension(:), intent(inout) :: index
      if (tree%branch(k)%has_children) then
         call tree_scan (tree, tree%branch(k)%daughter(1), perm, index)
         call tree_scan (tree, tree%branch(k)%daughter(2), perm, index)
         i = i + 1
         if (i <= size (index))  index(i) = tc_permute (k, perm, mask_in)
      end if
    end subroutine tree_scan
  end subroutine phs_tree_find_msq_permutation

  module subroutine phs_tree_find_angle_permutation &
       (tree1, tree2, perm2, angle_perm, sig2)
    type(phs_tree_t), intent(in) :: tree1, tree2
    type(permutation_t), intent(in) :: perm2
    type(permutation_t), intent(out) :: angle_perm
    logical, dimension(:), allocatable, intent(out) :: sig2
    type(permutation_t) :: perm1
    integer(TC) :: mask_in, root
    integer(TC), dimension(:), allocatable :: index1, index2
    logical, dimension(:), allocatable :: sig1
    integer :: i
    allocate (index1 (tree1%n_angles), index2 (tree2%n_angles))
    allocate (sig1 (tree1%n_angles), sig2 (tree2%n_angles))
    call permutation_init (perm1, permutation_size (perm2))
    mask_in = tree1%mask_in
    root = tree1%mask_out
    i = 0
    call tree_scan (tree1, root, perm1, index1, sig1)
    i = 0
    call tree_scan (tree2, root, perm2, index2, sig2)
    call permutation_find (angle_perm, index1, index2)
  contains
    recursive subroutine tree_scan (tree, k, perm, index, sig)
      type(phs_tree_t), intent(in) :: tree
      integer(TC), intent(in) :: k
      type(permutation_t), intent(in) :: perm
      integer, dimension(:), intent(inout) :: index
      logical, dimension(:), intent(inout) :: sig
      integer(TC) :: k1, k2, kp
      logical :: s
      if (tree%branch(k)%has_children) then
         k1 = tree%branch(k)%daughter(1)
         k2 = tree%branch(k)%daughter(2)
         s = (tc_permute(k1, perm, mask_in) < tc_permute(k2, perm, mask_in))
         kp = tc_permute (k, perm, mask_in)
         i = i + 1
         index(i) = kp
         sig(i) = s
         i = i + 1
         index(i) = - kp
         sig(i) = s
         call tree_scan (tree, k1, perm, index, sig)
         call tree_scan (tree, k2, perm, index, sig)
      end if
    end subroutine tree_scan
  end subroutine phs_tree_find_angle_permutation

  module subroutine phs_tree_compute_volume (tree, sqrts, volume)
    class(phs_tree_t), intent(in) :: tree
    real(default), intent(in) :: sqrts
    real(default), intent(out) :: volume
    integer(TC) :: k
    k  = tree%mask_out
    if (tree%branch(k)%has_children) then
       call compute_volume_x (tree%branch(k), k, volume, .true.)
    else
       volume = 1
    end if
  contains
    recursive subroutine compute_volume_x (b, k, volume, initial)
      type(phs_branch_t), intent(in) :: b
      integer(TC), intent(in) :: k
      real(default), intent(out) :: volume
      logical, intent(in) :: initial
      integer(TC) :: k1, k2
      real(default) :: v1, v2
      k1 = b%daughter(1);  k2 = b%daughter(2)
      if (tree%branch(k1)%has_children) then
         call compute_volume_x (tree%branch(k1), k1, v1, .false.)
      else
         v1 = 1
      end if
      if (tree%branch(k2)%has_children) then
         call compute_volume_x (tree%branch(k2), k2, v2, .false.)
      else
         v2 = 1
      end if
      if (initial) then
         volume = v1 * v2 / (4 * twopi5)
      else
         volume = v1 * v2 * sqrts**2 / (4 * twopi2)
      end if
    end subroutine compute_volume_x
  end subroutine phs_tree_compute_volume

  module subroutine phs_tree_compute_momenta_from_x &
       (tree, prt, factor, volume, sqrts, x, ok)
    class(phs_tree_t), intent(inout) :: tree
    type(phs_prt_t), dimension(:), intent(inout) :: prt
    real(default), intent(out) :: factor, volume
    real(default), intent(in) :: sqrts
    real(default), dimension(:), intent(in) :: x
    logical, intent(out) :: ok
    real(default), dimension(tree%mask_out) :: decay_p
    integer :: n1, n2
    integer :: n_out
    if (tree%real_phsp) then
      n_out = tree%n_externals - tree%n_in - 1
      n1 = max (n_out-2, 0)
      n2 = n1 + max (2*n_out, 0)
    else
      n1 = tree%n_msq
      n2 = n1 + tree%n_angles
    end if
    call phs_tree_set_msq &
         (tree, prt, factor, volume, decay_p, sqrts, x(1:n1), ok)
    if (ok) call phs_tree_set_angles &
         (tree, prt, factor, decay_p, sqrts, x(n1+1:n2))
  end subroutine phs_tree_compute_momenta_from_x

  subroutine phs_tree_set_msq &
       (tree, prt, factor, volume, decay_p, sqrts, x, ok)
    type(phs_tree_t), intent(inout) :: tree
    type(phs_prt_t), dimension(:), intent(inout) :: prt
    real(default), intent(out) :: factor, volume
    real(default), dimension(:), intent(out) :: decay_p
    real(default), intent(in) :: sqrts
    real(default), dimension(:), intent(in) :: x
    logical, intent(out) :: ok
    integer :: ix
    integer(TC) :: k
    real(default) :: m_tot
    ok =.true.
    ix = 1
    k  = tree%mask_out
    m_tot = tree%mass_sum(k)
    decay_p(k) = 0.
    if (m_tot < sqrts .or. k == 1) then
       if (tree%branch(k)%has_children) then
          call set_msq_x (tree%branch(k), k, factor, volume, .true.)
       else
          factor = 1
          volume = 1
       end if
    else
       ok = .false.
    end if
  contains
    recursive subroutine set_msq_x (b, k, factor, volume, initial)
      type(phs_branch_t), intent(in) :: b
      integer(TC), intent(in) :: k
      real(default), intent(out) :: factor, volume
      logical, intent(in) :: initial
      real(default) :: msq, m, m_min, m_max, m1, m2, msq1, msq2, lda, rlda
      integer(TC) :: k1, k2
      real(default) :: f1, f2, v1, v2
      k1 = b%daughter(1);  k2 = b%daughter(2)
      if (tree%branch(k1)%has_children) then
         call set_msq_x (tree%branch(k1), k1, f1, v1, .false.)
         if (.not.ok) return
      else
         f1 = 1;  v1 = 1
      end if
      if (tree%branch(k2)%has_children) then
         call set_msq_x (tree%branch(k2), k2, f2, v2, .false.)
         if (.not.ok) return
      else
         f2 = 1;  v2 = 1
      end if
      m_min = tree%mass_sum(k)
      if (initial) then
         msq = sqrts**2
         m = sqrts
         m_max = sqrts
         factor = f1 * f2
         volume = v1 * v2 / (4 * twopi5)
      else
         m_max = sqrts - m_tot + m_min
         call tree%mapping(k)%compute_msq_from_x (sqrts**2, m_min**2, &
              m_max**2, msq, factor, x(ix)); ix = ix + 1
         if (msq >= 0) then
            m = sqrt (msq)
            factor = f1 * f2 * factor
            volume = v1 * v2 * sqrts**2 / (4 * twopi2)
            call prt(k)%set_msq (msq)
            call prt(k)%set_defined ()
         else
            ok = .false.
         end if
      end if
      if (ok) then
         msq1 = prt(k1)%get_msq ();  m1 = sqrt (msq1)
         msq2 = prt(k2)%get_msq ();  m2 = sqrt (msq2)
         lda = lambda (msq, msq1, msq2)
         if (lda > 0 .and. m > m1 + m2 .and. m <= m_max) then
            rlda = sqrt (lda)
            decay_p(k1) = rlda / (2*m)
            decay_p(k2) = - decay_p(k1)
            factor = rlda / msq * factor
         else
            ok = .false.
         end if
      end if
    end subroutine set_msq_x

  end subroutine phs_tree_set_msq

  subroutine phs_tree_set_angles (tree, prt, factor, decay_p, sqrts, x)
    type(phs_tree_t), intent(inout) :: tree
    type(phs_prt_t), dimension(:), intent(inout) :: prt
    real(default), intent(inout) :: factor
    real(default), dimension(:), intent(in) :: decay_p
    real(default), intent(in) :: sqrts
    real(default), dimension(:), intent(in) :: x
    integer :: ix
    integer(TC) :: k
    ix = 1
    k  = tree%mask_out
    call set_angles_x (tree%branch(k), k)
  contains
    recursive subroutine set_angles_x (b, k, L0)
      type(phs_branch_t), intent(in) :: b
      integer(TC), intent(in) :: k
      type(lorentz_transformation_t), intent(in), optional :: L0
      real(default) :: m, msq, ct, st, phi, f, E, p, bg
      type(lorentz_transformation_t) :: L, LL
      integer(TC) :: k1, k2
      type(vector3_t) :: axis
      p = decay_p(k)
      msq = prt(k)%get_msq ();  m = sqrt (msq)
      E = sqrt (msq + p**2)
      if (present (L0)) then
         call prt(k)%set_momentum (L0 * vector4_moving (E,p,3))
      else
         call prt(k)%set_momentum (vector4_moving (E,p,3))
      end if
      call prt(k)%set_defined ()
      if (b%has_children) then
         k1 = b%daughter(1)
         k2 = b%daughter(2)
         if (m > 0) then
            bg = p / m
         else
            bg = 0
         end if
         phi = x(ix) * twopi;  ix = ix + 1
         call tree%mapping(k)%compute_ct_from_x (sqrts**2, ct, st, f, &
              x(ix));  ix = ix + 1
         factor = factor * f
         if (.not. b%has_friend) then
            L = LT_compose_r2_r3_b3 (ct, st, cos(phi), sin(phi), bg)
            !!! The function above is equivalent to:
            ! L = boost (bg,3) * rotation (phi,3) * rotation (ct,st,2)
         else
            LL = boost (-bg,3);  if (present (L0))  LL = LL * inverse(L0)
            axis = space_part ( &
                 LL * prt(tree%branch(k)%friend)%get_momentum () )
            L = boost(bg,3) * rotation_to_2nd (vector3_canonical(3), axis) &
                 * LT_compose_r2_r3_b3 (ct, st, cos(phi), sin(phi), 0._default)
         end if
         if (present (L0))  L = L0 * L
         call set_angles_x (tree%branch(k1), k1, L)
         call set_angles_x (tree%branch(k2), k2, L)
      end if
    end subroutine set_angles_x

  end subroutine phs_tree_set_angles

  module subroutine phs_tree_compute_x_from_momenta &
       (tree, prt, factor, sqrts, x)
    class(phs_tree_t), intent(inout) :: tree
    type(phs_prt_t), dimension(:), intent(in) :: prt
    real(default), intent(out) :: factor
    real(default), intent(in) :: sqrts
    real(default), dimension(:), intent(inout) :: x
    real(default), dimension(tree%mask_out) :: decay_p
    integer :: n1, n2
    n1 = tree%n_msq
    n2 = n1 + tree%n_angles
    call phs_tree_get_msq &
         (tree, prt, factor, decay_p, sqrts, x(1:n1))
    call phs_tree_get_angles &
         (tree, prt, factor, decay_p, sqrts, x(n1+1:n2))
  end subroutine phs_tree_compute_x_from_momenta

  subroutine phs_tree_get_msq (tree, prt, factor, decay_p, sqrts, x)
    type(phs_tree_t), intent(inout) :: tree
    type(phs_prt_t), dimension(:), intent(in) :: prt
    real(default), intent(out) :: factor
    real(default), dimension(:), intent(out) :: decay_p
    real(default), intent(in) :: sqrts
    real(default), dimension(:), intent(inout) :: x
    integer :: ix
    integer(TC) :: k
    real(default) :: m_tot
    ix = 1
    k  = tree%mask_out
    m_tot = tree%mass_sum(k)
    decay_p(k) = 0.
    if (tree%branch(k)%has_children) then
       call get_msq_x (tree%branch(k), k, factor, .true.)
    else
       factor = 1
    end if
  contains
    recursive subroutine get_msq_x (b, k, factor, initial)
      type(phs_branch_t), intent(in) :: b
      integer(TC), intent(in) :: k
      real(default), intent(out) :: factor
      logical, intent(in) :: initial
      real(default) :: msq, m, m_min, m_max, msq1, msq2, lda, rlda
      integer(TC) :: k1, k2
      real(default) :: f1, f2
      k1 = b%daughter(1);  k2 = b%daughter(2)
      if (tree%branch(k1)%has_children) then
         call get_msq_x (tree%branch(k1), k1, f1, .false.)
      else
         f1 = 1
      end if
      if (tree%branch(k2)%has_children) then
         call get_msq_x (tree%branch(k2), k2, f2, .false.)
      else
         f2 = 1
      end if
      m_min = tree%mass_sum(k)
      m_max = sqrts - m_tot + m_min
      msq = prt(k)%get_msq ();  m = sqrt (msq)
      if (initial) then
         factor = f1 * f2
      else
         call tree%mapping(k)%compute_x_from_msq (sqrts**2, m_min**2, &
              m_max**2, msq, factor, x(ix));  ix = ix + 1
         factor = f1 * f2 * factor
      end if
      msq1 = prt(k1)%get_msq ()
      msq2 = prt(k2)%get_msq ()
      lda = lambda (msq, msq1, msq2)
      if (lda > 0) then
         rlda = sqrt (lda)
         decay_p(k1) = rlda / (2 * m)
         decay_p(k2) = - decay_p(k1)
         factor = rlda / msq * factor
      else
         decay_p(k1) = 0
         decay_p(k2) = 0
         factor = 0
      end if
    end subroutine get_msq_x

  end subroutine phs_tree_get_msq

  subroutine phs_tree_get_angles (tree, prt, factor, decay_p, sqrts, x)
    type(phs_tree_t), intent(inout) :: tree
    type(phs_prt_t), dimension(:), intent(in) :: prt
    real(default), intent(inout) :: factor
    real(default), dimension(:), intent(in) :: decay_p
    real(default), intent(in) :: sqrts
    real(default), dimension(:), intent(out) :: x
    integer :: ix
    integer(TC) :: k
    ix = 1
    k  = tree%mask_out
    if (tree%branch(k)%has_children) then
       call get_angles_x (tree%branch(k), k)
    end if
  contains
    recursive subroutine get_angles_x (b, k, ct0, st0, phi0, L0)
      type(phs_branch_t), intent(in) :: b
      integer(TC), intent(in) :: k
      real(default), intent(in), optional :: ct0, st0, phi0
      type(lorentz_transformation_t), intent(in), optional :: L0
      real(default) :: cp0, sp0, m, msq, ct, st, phi, bg, f
      type(lorentz_transformation_t) :: L, LL
      type(vector4_t) :: p1, pf
      type(vector3_t) :: n, axis
      integer(TC) :: k1, k2, kf
      logical :: has_friend, need_L
      k1 = b%daughter(1)
      k2 = b%daughter(2)
      kf = b%friend
      has_friend = b%has_friend
      if (present(L0)) then
         p1 = L0 * prt(k1)%get_momentum ()
         if (has_friend)  pf = L0 * prt(kf)%get_momentum ()
      else
         p1 = prt(k1)%get_momentum ()
         if (has_friend)  pf = prt(kf)%get_momentum ()
      end if
      if (present(phi0)) then
         cp0 = cos (phi0)
         sp0 = sin (phi0)
      end if
      msq = prt(k)%get_msq ();  m = sqrt (msq)
      if (m > 0) then
         bg = decay_p(k) / m
      else
         bg = 0
      end if
      if (has_friend) then
         if (present (phi0)) then
            axis = axis_from_p_r3_r2_b3 (pf, cp0, -sp0, ct0, -st0, -bg)
            LL = rotation_to_2nd (axis, vector3_canonical (3)) &
                 * LT_compose_r3_r2_b3 (cp0, -sp0, ct0, -st0, -bg)
         else
            axis = axis_from_p_b3 (pf, -bg)
            LL = rotation_to_2nd (axis, vector3_canonical(3))
            if (.not. vanishes (bg))  LL = LL * boost(-bg, 3)
         end if
         n = space_part (LL * p1)
      else if (present (phi0)) then
         n = axis_from_p_r3_r2_b3 (p1, cp0, -sp0, ct0, -st0, -bg)
      else
         n = axis_from_p_b3 (p1, -bg)
      end if
      phi = azimuthal_angle (n)
      x(ix) = phi / twopi;  ix = ix + 1
      ct = polar_angle_ct (n)
      st = sqrt (1 - ct**2)
      call tree%mapping(k)%compute_x_from_ct (sqrts**2, ct, f, &
           x(ix)); ix = ix + 1
      factor = factor * f
      if (tree%branch(k1)%has_children .or. tree%branch(k2)%has_children) then
         need_L = .true.
         if (has_friend) then
            if (present (L0)) then
               L = LL * L0
            else
               L = LL
            end if
         else if (present (L0)) then
            L = LT_compose_r3_r2_b3 (cp0, -sp0, ct0, -st0, -bg) * L0
         else if (present (phi0)) then
            L = LT_compose_r3_r2_b3 (cp0, -sp0, ct0, -st0, -bg)
         else if (bg /= 0) then
            L = boost(-bg, 3)
         else
            need_L = .false.
         end if
         if (need_L) then
            if (tree%branch(k1)%has_children) &
                 call get_angles_x (tree%branch(k1), k1, ct, st, phi, L)
            if (tree%branch(k2)%has_children) &
                 call get_angles_x (tree%branch(k2), k2, ct, st, phi, L)
         else
            if (tree%branch(k1)%has_children) &
                 call get_angles_x (tree%branch(k1), k1, ct, st, phi)
            if (tree%branch(k2)%has_children) &
                 call get_angles_x (tree%branch(k2), k2, ct, st, phi)
         end if
      end if
    end subroutine get_angles_x
  end subroutine phs_tree_get_angles

  module subroutine phs_tree_combine_particles (tree, prt)
    type(phs_tree_t), intent(in) :: tree
    type(phs_prt_t), dimension(:), intent(inout) :: prt
    call combine_particles_x (tree%mask_out)
  contains
    recursive subroutine combine_particles_x (k)
      integer(TC), intent(in) :: k
      integer :: k1, k2
      if (tree%branch(k)%has_children) then
         k1 = tree%branch(k)%daughter(1);  k2 = tree%branch(k)%daughter(2)
         call combine_particles_x (k1)
         call combine_particles_x (k2)
         if (.not. prt(k)%defined) then
            call prt(k)%combine (prt(k1), prt(k2))
         end if
      end if
    end subroutine combine_particles_x
  end subroutine phs_tree_combine_particles

  module subroutine phs_tree_setup_prt_combinations (tree, comb)
    type(phs_tree_t), intent(in) :: tree
    integer, dimension(:,:), intent(out) :: comb
    comb = 0
    call setup_prt_combinations_x (tree%mask_out)
  contains
    recursive subroutine setup_prt_combinations_x (k)
      integer(TC), intent(in) :: k
      integer, dimension(2) :: kk
      if (tree%branch(k)%has_children) then
         kk = tree%branch(k)%daughter
         call setup_prt_combinations_x (kk(1))
         call setup_prt_combinations_x (kk(2))
         comb(:,k) = kk
      end if
    end subroutine setup_prt_combinations_x
  end subroutine phs_tree_setup_prt_combinations

  module subroutine phs_tree_reshuffle_mappings (tree)
   class(phs_tree_t), intent(inout) :: tree
   integer(TC) :: k0, k_old, k_new, k2
   integer :: i
   type(mapping_t) :: mapping_tmp
   real(default) :: mass_tmp
   do i = 1, size (tree%momentum_link)
     if (i /= tree%momentum_link (i)) then
       k_old = 2**(i-tree%n_in-1)
       k_new = 2**(tree%momentum_link(i)-tree%n_in-1)
       k0 = tree%branch(k_old)%mother
       k2 = k_new + tree%branch(k_old)%sibling
       mapping_tmp = tree%mapping(k0)
       mass_tmp = tree%mass_sum(k0)
       tree%mapping(k0) = tree%mapping(k2)
       tree%mapping(k2) = mapping_tmp
       tree%mass_sum(k0) = tree%mass_sum(k2)
       tree%mass_sum(k2) = mass_tmp
     end if
   end do
  end subroutine phs_tree_reshuffle_mappings

  module subroutine phs_tree_set_momentum_links (tree, list)
    type(phs_tree_t), intent(inout) :: tree
    integer, dimension(:), allocatable :: list
    tree%momentum_link = list
  end subroutine phs_tree_set_momentum_links


end submodule phs_trees_s

