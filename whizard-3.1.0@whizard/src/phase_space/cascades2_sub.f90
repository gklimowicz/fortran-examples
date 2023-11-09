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

submodule (cascades2) cascades2_s

  use sorting
  use io_units
  use physics_defs, only: SCALAR, SPINOR, VECTOR, VECTORSPINOR, TENSOR
  use hashes
  use cascades, only: phase_space_vanishes, MAX_WARN_RESONANCE

  implicit none

contains

  module subroutine part_prop_final (part)
    class(part_prop_t), intent(inout) :: part
    part%anti => null ()
  end subroutine part_prop_final

  module subroutine tree_final (tree)
    class(tree_t), intent(inout) :: tree
    if (allocated (tree%bc)) deallocate (tree%bc)
    if (allocated (tree%pdg)) deallocate (tree%pdg)
    if (allocated (tree%mapping)) deallocate (tree%mapping)
  end subroutine tree_final

  module subroutine tree_assign (tree1, tree2)
    type(tree_t), intent(inout) :: tree1
    type(tree_t), intent(in) :: tree2
    if (allocated (tree2%bc)) then
       allocate (tree1%bc(size(tree2%bc)))
       tree1%bc = tree2%bc
    end if
    if (allocated (tree2%pdg)) then
       allocate (tree1%pdg(size(tree2%pdg)))
       tree1%pdg = tree2%pdg
    end if
    if (allocated (tree2%mapping)) then
       allocate (tree1%mapping(size(tree2%mapping)))
       tree1%mapping = tree2%mapping
    end if
    tree1%n_entries = tree2%n_entries
    tree1%keep = tree2%keep
    tree1%empty = tree2%empty
  end subroutine tree_assign

  module subroutine tree_add_entry_from_numbers (tree, bincode, pdg, mapping)
    class(tree_t), intent(inout) :: tree
    integer(TC), intent(in) :: bincode
    integer, intent(in) :: pdg
    integer, intent(in) :: mapping
    integer :: pos
    if (tree%empty) then
       allocate (tree%bc(1))
       allocate (tree%pdg(1))
       allocate (tree%mapping(1))
       pos = tree%n_entries + 1
       tree%bc(pos) = bincode
       tree%pdg(pos) = pdg
       tree%mapping(pos) = mapping
       tree%n_entries = pos
       tree%empty = .false.
    end if
  end subroutine tree_add_entry_from_numbers

  subroutine tree_merge (tree, tree1, tree2, bc, pdg, mapping)
    class(tree_t), intent(inout) :: tree
    type(tree_t), intent(in) :: tree1, tree2
    integer(TC), intent(in) :: bc
    integer, intent(in) :: pdg, mapping
    integer :: tree_size
    integer :: i1, i2
    if (tree%empty) then
       i1 = tree1%n_entries
       i2 = tree1%n_entries + tree2%n_entries
       !! Proof: tree_size > 0 (always)
       tree_size = tree1%n_entries + tree2%n_entries + 1
       allocate (tree%bc (tree_size))
       allocate (tree%pdg (tree_size))
       allocate (tree%mapping (tree_size))
       if (.not. tree1%empty) then
          tree%bc(:i1) = tree1%bc
          tree%pdg(:i1) = tree1%pdg
          tree%mapping(:i1) = tree1%mapping
       end if
       if (.not. tree2%empty) then
          tree%bc(i1+1:i2) = tree2%bc
          tree%pdg(i1+1:i2) = tree2%pdg
          tree%mapping(i1+1:i2) = tree2%mapping
       end if
       tree%bc(tree_size) = bc
       tree%pdg(tree_size) = pdg
       tree%mapping(tree_size) = mapping
       tree%n_entries = tree_size
       tree%empty = .false.
    end if
  end subroutine tree_merge

  module subroutine tree_add_entry_from_node (tree, node)
    class(tree_t), intent(inout) :: tree
    type(k_node_t), intent(in) :: node
    integer :: pdg
    if (node%t_line) then
       pdg = abs (node%particle%pdg)
    else
       pdg = node%particle%pdg
    end if
    if (associated (node%daughter1) .and. &
         associated (node%daughter2)) then
       call tree_merge (tree, node%daughter1%subtree, &
            node%daughter2%subtree, node%bincode, &
            node%particle%pdg, node%mapping)
    else
       call tree_add_entry_from_numbers (tree, node%bincode, &
            node%particle%pdg, node%mapping)
    end if
    call tree%sort ()
  end subroutine tree_add_entry_from_node

  module subroutine tree_sort (tree)
    class(tree_t), intent(inout) :: tree
    integer(TC), dimension(size(tree%bc)) :: bc_tmp
    integer, dimension(size(tree%pdg)) :: pdg_tmp, mapping_tmp
    integer, dimension(1) :: pos
    integer :: i
    bc_tmp = tree%bc
    pdg_tmp = tree%pdg
    mapping_tmp = tree%mapping
    do i = size(tree%bc),1,-1
       pos = maxloc (bc_tmp)
       tree%bc(i) = bc_tmp (pos(1))
       tree%pdg(i) = pdg_tmp (pos(1))
       tree%mapping(i) = mapping_tmp (pos(1))
       bc_tmp(pos(1)) = 0
    end do
  end subroutine tree_sort

  module subroutine feyngraph_final (graph)
    class(feyngraph_t), intent(inout) :: graph
    type(kingraph_t), pointer :: current
    graph%root => null ()
    graph%kin_last => null ()
    do while (associated (graph%kin_first))
       current => graph%kin_first
       graph%kin_first => graph%kin_first%next
       call current%final ()
       deallocate (current)
    end do
  end subroutine feyngraph_final

  module subroutine kingraph_final (graph)
    class(kingraph_t), intent(inout) :: graph
    graph%root => null ()
    graph%next => null ()
    graph%grove_next => null ()
    call graph%tree%final ()
  end subroutine kingraph_final

  module subroutine k_node_entry_final (entry)
    class(k_node_entry_t), intent(inout) :: entry
    if (associated (entry%node)) then
       call entry%node%final
       deallocate (entry%node)
    end if
    entry%next => null ()
  end subroutine k_node_entry_final

  module subroutine k_node_entry_write (k_node_entry, u)
    class(k_node_entry_t), intent(in) :: k_node_entry
    integer, intent(in) :: u
  end subroutine k_node_entry_write

  module subroutine k_node_list_final (list)
    class(k_node_list_t), intent(inout) :: list
    type(k_node_entry_t), pointer :: current
    do while (associated (list%first))
       current => list%first
       list%first => list%first%next
       if (list%observer) current%node => null ()
       call current%final ()
       deallocate (current)
    end do
  end subroutine k_node_list_final

  recursive module subroutine f_node_final (node)
    class(f_node_t), intent(inout) :: node
    call node%k_node_list%final ()
    node%daughter1 => null ()
    node%daughter2 => null ()
  end subroutine f_node_final

  module subroutine f_node_entry_final (entry)
    class(f_node_entry_t), intent(inout) :: entry
    if (associated (entry%node)) then
       call entry%node%final ()
       deallocate (entry%node)
    end if
    entry%next => null ()
  end subroutine f_node_entry_final

  module subroutine f_node_set_index (f_node)
    class(f_node_t), intent(inout) :: f_node
    integer, save :: counter = 0
    if (f_node%index == 0) then
       counter = counter + 1
       f_node%index = counter
    end if
  end subroutine f_node_set_index

  module subroutine f_node_ptr_final (f_node_ptr)
    class(f_node_ptr_t), intent(inout) :: f_node_ptr
    f_node_ptr%node => null ()
  end subroutine f_node_ptr_final

  module subroutine f_node_ptr_assign (ptr1, ptr2)
    type(f_node_ptr_t), intent(out) :: ptr1
    type(f_node_ptr_t), intent(in) :: ptr2
    ptr1%node => ptr2%node
  end subroutine f_node_ptr_assign

  module subroutine k_node_assign (k_node1, k_node2)
    type(k_node_t), intent(inout) :: k_node1
    type(k_node_t), intent(in) :: k_node2
    k_node1%f_node => k_node2%f_node
    k_node1%particle => k_node2%particle
    k_node1%incoming = k_node2%incoming
    k_node1%t_line = k_node2%t_line
    k_node1%keep = k_node2%keep
    k_node1%n_subtree_nodes = k_node2%n_subtree_nodes
    k_node1%ext_mass_sum = k_node2%ext_mass_sum
    k_node1%effective_mass = k_node2%effective_mass
    k_node1%resonant = k_node2%resonant
    k_node1%on_shell = k_node2%on_shell
    k_node1%log_enhanced = k_node2%log_enhanced
    k_node1%mapping = k_node2%mapping
    k_node1%bincode = k_node2%bincode
    k_node1%mapping_assigned = k_node2%mapping_assigned
    k_node1%is_nonresonant_copy = k_node2%is_nonresonant_copy
    k_node1%n_off_shell = k_node2%n_off_shell
    k_node1%n_log_enhanced = k_node2%n_log_enhanced
    k_node1%n_resonances = k_node2%n_resonances
    k_node1%multiplicity = k_node2%multiplicity
    k_node1%n_t_channel = k_node2%n_t_channel
    k_node1%f_node_index = k_node2%f_node_index
  end subroutine k_node_assign

  recursive module subroutine k_node_final (k_node)
    class(k_node_t), intent(inout) :: k_node
    k_node%daughter1 => null ()
    k_node%daughter2 => null ()
    k_node%inverse_daughter1 => null ()
    k_node%inverse_daughter2 => null ()
    k_node%f_node => null ()
  end subroutine k_node_final

  module subroutine k_node_set_index (k_node)
    class(k_node_t), intent(inout) :: k_node
    integer, save :: counter = 0
    if (k_node%index == 0) then
       counter = counter + 1
       k_node%index = counter
    end if
  end subroutine k_node_set_index

  module subroutine f_node_entry_write (f_node_entry, u)
    class(f_node_entry_t), intent(in) :: f_node_entry
    integer, intent(in) :: u
    write (unit=u, fmt='(A)') trim(f_node_entry%subtree_string)
  end subroutine f_node_entry_write

  module subroutine f_node_entry_assign (entry1, entry2)
    type(f_node_entry_t), intent(out) :: entry1
    type(f_node_entry_t), intent(in) :: entry2
    entry1%node => entry2%node
    entry1%subtree_string = entry2%subtree_string
    entry1%string_len = entry2%string_len
    entry1%subtree_size = entry2%subtree_size
  end subroutine f_node_entry_assign

  module subroutine f_node_list_add_entry (list, subtree_string, &
       ptr_to_node, recycle, subtree_size)
    class(f_node_list_t), intent(inout) :: list
    character(len=*), intent(in) :: subtree_string
    type(f_node_t), pointer, intent(out) :: ptr_to_node
    logical, intent(in) :: recycle
    integer, intent(in), optional :: subtree_size
    type(f_node_entry_t), pointer :: current
    type(f_node_entry_t), pointer :: second
    integer :: subtree_len
    ptr_to_node => null ()
    if (recycle) then
       subtree_len = len_trim (subtree_string)
       current => list%first
       do while (associated (current))
          if (present (subtree_size)) then
             if (current%subtree_size /= subtree_size) exit
          end if
          if (current%string_len == subtree_len) then
             if (trim (current%subtree_string) == trim (subtree_string)) then
                ptr_to_node => current%node
                exit
             end if
          end if
          current => current%next
       end do
    end if
    if (.not. associated (ptr_to_node)) then
       if (list%n_entries == 0) then
          allocate (list%first)
          list%last => list%first
       else
          second => list%first
          list%first => null ()
          allocate (list%first)
          list%first%next => second
       end if
       list%n_entries = list%n_entries + 1
       list%first%subtree_string = trim(subtree_string)
       list%first%string_len = subtree_len
       if (present (subtree_size)) list%first%subtree_size = subtree_size
       allocate (list%first%node)
       call list%first%node%set_index ()
       ptr_to_node => list%first%node
    end if
  end subroutine f_node_list_add_entry

  module subroutine f_node_list_write (f_node_list, u)
    class(f_node_list_t), intent(in) :: f_node_list
    integer, intent(in) :: u
    type(f_node_entry_t), pointer :: current
    integer :: pos = 0
    current => f_node_list%first
    do while (associated (current))
       pos = pos + 1
       write (unit=u, fmt='(A,I10)') 'entry #: ', pos
       call current%write (u)
       write (unit=u, fmt=*)
       current => current%next
    end do
  end subroutine f_node_list_write

  module subroutine k_node_entry_assign (entry1, entry2)
    type(k_node_entry_t), intent(out) :: entry1
    type(k_node_entry_t), intent(in) :: entry2
    entry1%node => entry2%node
    entry1%recycle = entry2%recycle
  end subroutine k_node_entry_assign

  recursive module subroutine k_node_list_add_entry &
       (list, ptr_to_node, recycle)
    class(k_node_list_t), intent(inout) :: list
    type(k_node_t), pointer, intent(out) :: ptr_to_node
    logical, intent(in) :: recycle
    if (list%n_entries == 0) then
       allocate (list%first)
       list%last => list%first
    else
       allocate (list%last%next)
       list%last => list%last%next
    end if
    list%n_entries = list%n_entries + 1
    list%last%recycle = recycle
    allocate (list%last%node)
    call list%last%node%set_index ()
    ptr_to_node => list%last%node
  end subroutine k_node_list_add_entry

  module subroutine k_node_list_add_pointer (list, ptr_to_node, recycle)
    class(k_node_list_t), intent(inout) :: list
    type(k_node_t), pointer, intent(in) :: ptr_to_node
    logical, optional, intent(in) :: recycle
    logical :: rec
    if (present (recycle)) then
       rec = recycle
    else
       rec = .false.
    end if
    if (list%n_entries == 0) then
       allocate (list%first)
       list%last => list%first
    else
       allocate (list%last%next)
       list%last => list%last%next
    end if
    list%n_entries = list%n_entries + 1
    list%last%recycle = rec
    list%last%node => ptr_to_node
  end subroutine k_node_list_add_pointer

  module subroutine k_node_list_check_subtree_equivalences (list, model)
    class(k_node_list_t), intent(inout) :: list
    type(model_data_t), intent(in) :: model
    type(k_node_ptr_t), dimension (:), allocatable :: set
    type(k_node_entry_t), pointer :: current
    integer :: pos
    integer :: i,j
    if (list%n_entries == 0) return
    allocate (set (list%n_entries))
    current => list%first
    pos = 0
    do while (associated (current))
       pos = pos + 1
       set(pos)%node => current%node
       current => current%next
    end do
    do i=1, list%n_entries
       if (set(i)%node%keep) then
          do j=i+1, list%n_entries
             if (set(j)%node%keep) then
                if (set(i)%node%bincode == set(j)%node%bincode) then
                   call subtree_select (set(i)%node%subtree,set(j)%node%subtree, model)
                   if (.not. set(i)%node%subtree%keep) then
                      set(i)%node%keep = .false.
                      exit
                   else if (.not. set(j)%node%subtree%keep) then
                      set(j)%node%keep = .false.
                   end if
                end if
             end if
          end do
       end if
    end do
    deallocate (set)
  end subroutine k_node_list_check_subtree_equivalences

  module subroutine k_node_list_get_nodes (list, nodes)
    class(k_node_list_t), intent(inout) :: list
    type(k_node_ptr_t), dimension(:), allocatable, intent(out) :: nodes
    integer :: n_nodes
    integer :: pos
    type(k_node_entry_t), pointer :: current, garbage
    n_nodes = 0
    current => list%first
    do while (associated (current))
       if (current%recycle .and. current%node%keep) n_nodes = n_nodes + 1
       current => current%next
    end do
    if (n_nodes /= 0) then
       pos = 1
       allocate (nodes (n_nodes))
       do while (associated (list%first) .and. .not. list%first%node%keep)
          garbage => list%first
          list%first => list%first%next
          call garbage%final ()
          deallocate (garbage)
       end do
       current => list%first
       do while (associated (current))
          do while (associated (current%next))
             if (.not. current%next%node%keep) then
                garbage => current%next
                current%next => current%next%next
                call garbage%final
                deallocate (garbage)
             else
                exit
             end if
          end do
          if (current%recycle .and. current%node%keep) then
             nodes(pos)%node => current%node
             pos = pos + 1
          end if
          current => current%next
       end do
    end if
  end subroutine k_node_list_get_nodes

  module subroutine compare_tree_final (ctree)
    class(compare_tree_t), intent(inout) :: ctree
    integer :: i
    if (associated (ctree%entry)) then
       do i=1, size (ctree%entry)
          call ctree%entry(i)%final ()
          deallocate (ctree%entry)
       end do
    end if
  end subroutine compare_tree_final

  recursive module subroutine compare_tree_entry_final (ct_entry)
    class(compare_tree_entry_t), intent(inout) :: ct_entry
    integer :: i
    if (associated (ct_entry%entry)) then
       do i=1, size (ct_entry%entry)
          call ct_entry%entry(i)%final ()
       end do
       deallocate (ct_entry%entry)
    else
       deallocate (ct_entry%graph_entry)
    end if
  end subroutine compare_tree_entry_final

  module subroutine compare_tree_check_kingraph &
       (ctree, kingraph, model, preliminary)
    class(compare_tree_t), intent(inout) :: ctree
    type(kingraph_t), intent(inout), pointer :: kingraph
    type(model_data_t), intent(in) :: model
    logical, intent(in) :: preliminary
    integer :: i
    integer :: pos
    integer(TC) :: sz
    integer(TC), dimension(:), allocatable :: identifier
    if (.not. associated (ctree%entry)) then
       sz = 0_TC
       do i = size(kingraph%tree%bc), 1, -1
          sz = ior (sz, kingraph%tree%bc(i))
       end do
       if (sz > 0) then
          allocate (ctree%entry (sz))
       else
          call msg_bug ("Compare tree could not be created")
       end if
    end if
    allocate (identifier (ctree%depth))
    pos = 0
    do i = size(kingraph%tree%bc), 1, -1
       if (popcnt (kingraph%tree%bc(i)) /= 1) then
          pos = pos + 1
          identifier(pos) = kingraph%tree%bc(i)
          if (pos == ctree%depth) exit
       end if
    end do
    if (size (identifier) > 1) then
       call ctree%entry(identifier(1))%check_kingraph (kingraph, model, &
            preliminary, identifier(1), identifier(2:))
    else if (size (identifier) == 1) then
       call ctree%entry(identifier(1))%check_kingraph &
            (kingraph, model, preliminary)
    end if
    deallocate (identifier)
  end subroutine compare_tree_check_kingraph

  recursive module subroutine compare_tree_entry_check_kingraph (ct_entry, &
       kingraph, model, preliminary, subtree_size, identifier)
    class(compare_tree_entry_t), intent(inout) :: ct_entry
    type(kingraph_t), pointer, intent(inout) :: kingraph
    type(model_data_t), intent(in) :: model
    logical, intent(in) :: preliminary
    integer, intent(in), optional :: subtree_size
    integer, dimension (:), intent(in), optional :: identifier
    if (present (identifier)) then
       if (.not. associated (ct_entry%entry)) &
            allocate (ct_entry%entry(subtree_size))
       if (size (identifier) > 1) then
          call ct_entry%entry(identifier(1))%check_kingraph (kingraph, &
               model, preliminary, identifier(1), identifier(2:))
       else if (size (identifier) == 1) then
          call ct_entry%entry(identifier(1))%check_kingraph (kingraph, &
               model, preliminary)
       end if
    else
       if (allocated (ct_entry%graph_entry)) then
          call perform_check
       else
          allocate (ct_entry%graph_entry(1))
          ct_entry%graph_entry(1)%graph => kingraph
       end if
    end if

    contains

      subroutine perform_check
        integer :: i
        logical :: rebuild
        rebuild = .true.
        do i=1, size(ct_entry%graph_entry)
           if (ct_entry%graph_entry(i)%graph%keep) then
              if (preliminary .or. &
                   ct_entry%graph_entry(i)%graph%prc_component /= &
                   kingraph%prc_component) then
                 call kingraph_select (ct_entry%graph_entry(i)%graph, &
                      kingraph, model, preliminary)
                 if (.not. kingraph%keep) then
                    return
                 else if (rebuild .and. .not. &
                      ct_entry%graph_entry(i)%graph%keep) then
                    ct_entry%graph_entry(i)%graph => kingraph
                    rebuild = .false.
                 end if
              end if
           end if
        end do
        if (rebuild) call rebuild_graph_entry
      end subroutine perform_check

      subroutine rebuild_graph_entry
        type(kingraph_ptr_t), dimension(:), allocatable :: tmp_ptr
        integer :: i
        integer :: pos
        allocate (tmp_ptr(size(ct_entry%graph_entry)+1))
        pos = 0
        do i=1, size(ct_entry%graph_entry)
           pos = pos + 1
           tmp_ptr(pos)%graph => ct_entry%graph_entry(i)%graph
        end do
        pos = pos + 1
        tmp_ptr(pos)%graph => kingraph
        deallocate (ct_entry%graph_entry)
        allocate (ct_entry%graph_entry (pos))
        do i=1, pos
           ct_entry%graph_entry(i)%graph => tmp_ptr(i)%graph
        end do
        deallocate (tmp_ptr)
      end subroutine rebuild_graph_entry
  end subroutine compare_tree_entry_check_kingraph

  module subroutine grove_final (grove)
    class(grove_t), intent(inout) :: grove
    grove%first => null ()
    grove%last  => null ()
    grove%next => null ()
  end subroutine grove_final

  module subroutine feyngraph_set_build (feyngraph_set, u_in)
    class(feyngraph_set_t), intent(inout) :: feyngraph_set
    integer, intent(in) :: u_in
    integer :: stat = 0
    character(len=FEYNGRAPH_LEN) :: omega_feyngraph_output
    type(feyngraph_t), pointer :: current_graph
    type(feyngraph_t), pointer :: compare_graph
    logical :: present
    if (feyngraph_set%use_dag) then
       allocate (feyngraph_set%dag)
       if (.not. associated (feyngraph_set%first)) then
          call feyngraph_set%dag%read_string (u_in, feyngraph_set%flv(:,1))
          call feyngraph_set%dag%construct (feyngraph_set)
          call feyngraph_set%dag%make_feyngraphs (feyngraph_set)
       end if
    else
       if (.not. associated (feyngraph_set%first)) then
          read (unit=u_in, fmt='(A)', iostat=stat, advance='yes') &
               omega_feyngraph_output
          if (omega_feyngraph_output(1:1) == '(') then
             allocate (feyngraph_set%first)
             feyngraph_set%first%omega_feyngraph_output = &
                  trim(omega_feyngraph_output)
             feyngraph_set%last => feyngraph_set%first
             feyngraph_set%n_graphs = feyngraph_set%n_graphs + 1
          else
             call msg_fatal ("Invalid input file")
          end if
          read (unit=u_in, fmt='(A)', iostat=stat, advance='yes') &
               omega_feyngraph_output
          do while (stat == 0)
             if (omega_feyngraph_output(1:1) == '(') then
                compare_graph => feyngraph_set%first
                present = .false.
                do while (associated (compare_graph))
                   if (len_trim(compare_graph%omega_feyngraph_output) &
                        == len_trim(omega_feyngraph_output)) then
                      if (compare_graph%omega_feyngraph_output == &
                           omega_feyngraph_output) then
                         present = .true.
                         exit
                      end if
                   end if
                   compare_graph => compare_graph%next
                end do
                if (.not. present) then
                   allocate (feyngraph_set%last%next)
                   feyngraph_set%last => feyngraph_set%last%next
                   feyngraph_set%last%omega_feyngraph_output = &
                        trim(omega_feyngraph_output)
                   feyngraph_set%n_graphs = feyngraph_set%n_graphs + 1
                end if
                read (unit=u_in, fmt='(A)', iostat=stat, advance='yes') &
                     omega_feyngraph_output
             else
                exit
             end if
          end do
          current_graph => feyngraph_set%first
          do while (associated (current_graph))
             call feyngraph_construct (feyngraph_set, current_graph)
             current_graph => current_graph%next
          end do
          feyngraph_set%f_node_list%max_tree_size = feyngraph_set%first%n_nodes
       end if
    end if
  end subroutine feyngraph_set_build

  module subroutine dag_read_string (dag, u_in, flv)
    class(dag_t), intent(inout) :: dag
    integer, intent(in) :: u_in
    type(flavor_t), dimension(:), intent(in) :: flv
    character(len=BUFFER_LEN) :: process_string
    logical :: process_found
    logical :: rewound
    !!! Find process string in file
    process_found = .false.
    rewound = .false.
    do while (.not. process_found)
       process_string = ""
       read (unit=u_in, fmt='(A)') process_string
       if (len_trim(process_string) /= 0) then
          if (index (process_string, "::") > 0) then
             process_found = process_string_match (trim (process_string), flv)
          end if
       else if (.not. rewound) then
          rewind (u_in)
          rewound = .true.
       else
          call msg_bug ("Process string not found in O'Mega input file.")
       end if
    end do
    call fds_file_get_line (u_in, dag%string)
    call dag%string%clean ()
    if (.not. allocated (dag%string%t) .or. dag%string%char_len == 0) &
         call msg_bug ("Process string not found in O'Mega input file.")
  end subroutine dag_read_string

  subroutine fds_file_get_line (u, string)
    integer, intent(in) :: u
    type(dag_string_t), intent(out) :: string
    type(dag_chain_t) :: chain
    integer :: string_size, current_len
    character(len=BUFFER_LEN) :: buffer
    integer :: fragment_len
    integer :: stat
    current_len = 0
    stat = 0
    string_size = 0
    do while (stat == 0)
       read (unit=u, fmt='(A)', iostat=stat) buffer
       if (stat /= 0) exit
       fragment_len = len_trim (buffer)
       if (fragment_len == 0) then
          exit
       else if (buffer (fragment_len:fragment_len) == BACKSLASH_CHAR) then
          fragment_len = fragment_len - 1
       end if
       call chain%append (buffer(:fragment_len))
       if (buffer(fragment_len+1:fragment_len+1) /= BACKSLASH_CHAR) exit
    end do
    if (associated (chain%first)) then
       call chain%compress ()
       string = chain%first
       call chain%final ()
    end if
  end subroutine fds_file_get_line

  function process_string_match (string, flv) result (match)
    character(len=*), intent(in) :: string
    type(flavor_t), dimension(:), intent(in) :: flv
    logical :: match
    integer :: pos
    integer :: occurence
    integer :: i
    pos = 1
    match = .false.
    do i=1, size (flv)
       occurence = index (string(pos:), char(flv(i)%get_name()))
       if (occurence > 0) then
          pos = pos + occurence
          match = .true.
       else
          match = .false.
          exit
       end if
    end do
  end function process_string_match

  module subroutine init_sm_full_test (model)
    class(model_data_t), intent(out) :: model
    type(field_data_t), pointer :: field
    integer, parameter :: n_real = 17
    integer, parameter :: n_field = 21
    integer, parameter :: n_vtx = 56
    integer :: i
    call model%init (var_str ("SM_vertex_test"), &
         n_real, 0, n_field, n_vtx)
    call model%init_par (1, var_str ("mZ"), 91.1882_default)
    call model%init_par (2, var_str ("mW"), 80.419_default)
    call model%init_par (3, var_str ("mH"), 125._default)
    call model%init_par (4, var_str ("me"), 0.000510997_default)
    call model%init_par (5, var_str ("mmu"), 0.105658389_default)
    call model%init_par (6, var_str ("mtau"), 1.77705_default)
    call model%init_par (7, var_str ("ms"), 0.095_default)
    call model%init_par (8, var_str ("mc"), 1.2_default)
    call model%init_par (9, var_str ("mb"), 4.2_default)
    call model%init_par (10, var_str ("mtop"), 173.1_default)
    call model%init_par (11, var_str ("wtop"), 1.523_default)
    call model%init_par (12, var_str ("wZ"), 2.443_default)
    call model%init_par (13, var_str ("wW"), 2.049_default)
    call model%init_par (14, var_str ("wH"), 0.004143_default)
    call model%init_par (15, var_str ("ee"), 0.3079561542961_default)
    call model%init_par (16, var_str ("cw"), 8.819013863636E-01_default)
    call model%init_par (17, var_str ("sw"), 4.714339240339E-01_default)
    i = 0
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("D_QUARK"), 1)
    call field%set (spin_type=2, color_type=3, charge_type=-2, isospin_type=-2)
    call field%set (name = [var_str ("d")], anti = [var_str ("dbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("U_QUARK"), 2)
    call field%set (spin_type=2, color_type=3, charge_type=3, isospin_type=2)
    call field%set (name = [var_str ("u")], anti = [var_str ("ubar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("S_QUARK"), 3)
    call field%set (spin_type=2, color_type=3, charge_type=-2, isospin_type=-2)
    call field%set (mass_data=model%get_par_real_ptr (7))
    call field%set (name = [var_str ("s")], anti = [var_str ("sbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("C_QUARK"), 4)
    call field%set (spin_type=2, color_type=3, charge_type=3, isospin_type=2)
    call field%set (mass_data=model%get_par_real_ptr (8))
    call field%set (name = [var_str ("c")], anti = [var_str ("cbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("B_QUARK"), 5)
    call field%set (spin_type=2, color_type=3, charge_type=-2, isospin_type=-2)
    call field%set (mass_data=model%get_par_real_ptr (9))
    call field%set (name = [var_str ("b")], anti = [var_str ("bbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("T_QUARK"), 6)
    call field%set (spin_type=2, color_type=3, charge_type=3, isospin_type=2)
    call field%set (mass_data=model%get_par_real_ptr (10))
    call field%set (width_data=model%get_par_real_ptr (11))
    call field%set (name = [var_str ("t")], anti = [var_str ("tbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("E_LEPTON"), 11)
    call field%set (spin_type=2)
    call field%set (mass_data=model%get_par_real_ptr (4))
    call field%set (name = [var_str ("e-")], anti = [var_str ("e+")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("E_NEUTRINO"), 12)
    call field%set (spin_type=2, is_left_handed=.true.)
    call field%set (name = [var_str ("nue")], anti = [var_str ("nuebar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("MU_LEPTON"), 13)
    call field%set (spin_type=2)
    call field%set (mass_data=model%get_par_real_ptr (5))
    call field%set (name = [var_str ("mu-")], anti = [var_str ("mu+")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("MU_NEUTRINO"), 14)
    call field%set (spin_type=2, is_left_handed=.true.)
    call field%set (name = [var_str ("numu")], anti = [var_str ("numubar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("TAU_LEPTON"), 15)
    call field%set (spin_type=2)
    call field%set (mass_data=model%get_par_real_ptr (6))
    call field%set (name = [var_str ("tau-")], anti = [var_str ("tau+")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("TAU_NEUTRINO"), 16)
    call field%set (spin_type=2, is_left_handed=.true.)
    call field%set (name = [var_str ("nutau")], anti = [var_str ("nutaubar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("GLUON"), 21)
    call field%set (spin_type=3, color_type=8)
    call field%set (name = [var_str ("gl")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("PHOTON"), 22)
    call field%set (spin_type=3)
    call field%set (name = [var_str ("A")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("Z_BOSON"), 23)
    call field%set (spin_type=3)
    call field%set (mass_data=model%get_par_real_ptr (1))
    call field%set (width_data=model%get_par_real_ptr (12))
    call field%set (name = [var_str ("Z")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("W_BOSON"), 24)
    call field%set (spin_type=3)
    call field%set (mass_data=model%get_par_real_ptr (2))
    call field%set (width_data=model%get_par_real_ptr (13))
    call field%set (name = [var_str ("W+")], anti = [var_str ("W-")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("HIGGS"), 25)
    call field%set (spin_type=1)
    call field%set (mass_data=model%get_par_real_ptr (3))
    call field%set (width_data=model%get_par_real_ptr (14))
    call field%set (name = [var_str ("H")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("PROTON"), 2212)
    call field%set (spin_type=2)
    call field%set (name = [var_str ("p")], anti = [var_str ("pbar")])
!    call field%set (mass_data=model%get_par_real_ptr (12))
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("HADRON_REMNANT_SINGLET"), 91)
    call field%set (color_type=1)
    call field%set (name = [var_str ("hr1")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("HADRON_REMNANT_TRIPLET"), 92)
    call field%set (color_type=3)
    call field%set (name = [var_str ("hr3")], anti = [var_str ("hr3bar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("HADRON_REMNANT_OCTET"), 93)
    call field%set (color_type=8)
    call field%set (name = [var_str ("hr8")])
    call model%freeze_fields ()
    i = 0
    i = i + 1
!!! QED
    call model%set_vertex (i, [var_str ("dbar"), var_str ("d"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("ubar"), var_str ("u"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("sbar"), var_str ("s"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("cbar"), var_str ("c"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("bbar"), var_str ("b"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("tbar"), var_str ("t"), var_str ("A")])
    i = i + 1
!!!
    call model%set_vertex (i, [var_str ("e+"), var_str ("e-"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("mu+"), var_str ("mu-"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("tau+"), var_str ("tau-"), var_str ("A")])
    i = i + 1
!!! QCD
    call model%set_vertex (i, [var_str ("gl"), var_str ("gl"), var_str ("gl")])
    i = i + 1
    call model%set_vertex (i, [var_str ("gl"), var_str ("gl"), &
         var_str ("gl"), var_str ("gl")])
    i = i + 1
!!!
    call model%set_vertex (i, [var_str ("dbar"), var_str ("d"), var_str ("gl")])
    i = i + 1
    call model%set_vertex (i, [var_str ("ubar"), var_str ("u"), var_str ("gl")])
    i = i + 1
    call model%set_vertex (i, [var_str ("sbar"), var_str ("s"), var_str ("gl")])
    i = i + 1
    call model%set_vertex (i, [var_str ("cbar"), var_str ("c"), var_str ("gl")])
    i = i + 1
    call model%set_vertex (i, [var_str ("bbar"), var_str ("b"), var_str ("gl")])
    i = i + 1
    call model%set_vertex (i, [var_str ("tbar"), var_str ("t"), var_str ("gl")])
    i = i + 1
!!! Neutral currents
    call model%set_vertex (i, [var_str ("dbar"), var_str ("d"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("ubar"), var_str ("u"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("sbar"), var_str ("s"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("cbar"), var_str ("c"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("bbar"), var_str ("b"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("tbar"), var_str ("t"), var_str ("Z")])
    i = i + 1
!!!
    call model%set_vertex (i, [var_str ("e+"), var_str ("e-"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("mu+"), var_str ("muu-"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("tau+"), var_str ("tau-"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("nuebar"), var_str ("nue"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("numubar"), var_str ("numu"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("nutaubar"), var_str ("nutau"), &
         var_str ("Z")])
    i = i + 1
!!! Charged currents
    call model%set_vertex (i, [var_str ("ubar"), var_str ("d"), var_str ("W+")])
    i = i + 1
    call model%set_vertex (i, [var_str ("cbar"), var_str ("s"), var_str ("W+")])
    i = i + 1
    call model%set_vertex (i, [var_str ("tbar"), var_str ("b"), var_str ("W+")])
    i = i + 1
    call model%set_vertex (i, [var_str ("dbar"), var_str ("u"), var_str ("W-")])
    i = i + 1
    call model%set_vertex (i, [var_str ("sbar"), var_str ("c"), var_str ("W-")])
    i = i + 1
    call model%set_vertex (i, [var_str ("bbar"), var_str ("t"), var_str ("W-")])
    i = i + 1
!!!
    call model%set_vertex (i, [var_str ("nuebar"), var_str ("e-"), var_str ("W+")])
    i = i + 1
    call model%set_vertex (i, [var_str ("numubar"), var_str ("mu-"), var_str ("W+")])
    i = i + 1
    call model%set_vertex (i, [var_str ("nutaubar"), var_str ("tau-"), var_str ("W+")])
    i = i + 1
    call model%set_vertex (i, [var_str ("e+"), var_str ("nue"), var_str ("W-")])
    i = i + 1
    call model%set_vertex (i, [var_str ("mu+"), var_str ("numu"), var_str ("W-")])
    i = i + 1
    call model%set_vertex (i, [var_str ("tau+"), var_str ("nutau"), var_str ("W-")])
    i = i + 1
!!! Yukawa
!!! keeping only 3rd generation for the moment
    ! call model%set_vertex (i, [var_str ("sbar"), var_str ("s"), var_str ("H")])
    ! i = i + 1
    ! call model%set_vertex (i, [var_str ("cbar"), var_str ("c"), var_str ("H")])
    ! i = i + 1
    call model%set_vertex (i, [var_str ("bbar"), var_str ("b"), var_str ("H")])
    i = i + 1
    call model%set_vertex (i, [var_str ("tbar"), var_str ("t"), var_str ("H")])
    i = i + 1
    ! call model%set_vertex (i, [var_str ("mubar"), var_str ("mu"), var_str ("H")])
    ! i = i + 1
    call model%set_vertex (i, [var_str ("taubar"), var_str ("tau"), var_str ("H")])
    i = i + 1
!!! Vector-boson self-interactions
    call model%set_vertex (i, [var_str ("W+"), var_str ("W-"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("W+"), var_str ("W-"), var_str ("Z")])
    i = i + 1
!!!
    call model%set_vertex (i, [var_str ("W+"), var_str ("W-"), var_str ("Z"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("W+"), var_str ("W+"), var_str ("W-"), var_str ("W-")])
    i = i + 1
    call model%set_vertex (i, [var_str ("W+"), var_str ("W-"), var_str ("Z"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("W+"), var_str ("W-"), var_str ("A"), var_str ("A")])
    i = i + 1
!!! Higgs - vector boson
    ! call model%set_vertex (i, [var_str ("H"), var_str ("Z"), var_str ("A")])
    ! i = i + 1
    ! call model%set_vertex (i, [var_str ("H"), var_str ("A"), var_str ("A")])
    ! i = i + 1
    ! call model%set_vertex (i, [var_str ("H"), var_str ("gl"), var_str ("gl")])
    ! i = i + 1
!!!
    call model%set_vertex (i, [var_str ("H"), var_str ("W+"), var_str ("W-")])
    i = i + 1
    call model%set_vertex (i, [var_str ("H"), var_str ("Z"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("H"), var_str ("H"), var_str ("W+"), var_str ("W-")])
    i = i + 1
    call model%set_vertex (i, [var_str ("H"), var_str ("H"), var_str ("Z"), var_str ("Z")])
    i = i + 1
!!! Higgs self-interactions
    call model%set_vertex (i, [var_str ("H"), var_str ("H"), var_str ("H")])
    i = i + 1
    call model%set_vertex (i, [var_str ("H"), var_str ("H"), var_str ("H"), var_str ("H")])
    i = i + 1
    call model%freeze_vertices ()
  end subroutine init_sm_full_test

  recursive module subroutine part_prop_init &
       (part_prop, feyngraph_set, particle_label)
    class(part_prop_t), intent(out), target :: part_prop
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    character(len=*), intent(in) :: particle_label
    type(flavor_t) :: flv, anti
    type(string_t) :: name
    integer :: i
    name = particle_label
    call flv%init (name, feyngraph_set%model)
    part_prop%particle_label = particle_label
    part_prop%pdg = flv%get_pdg ()
    part_prop%mass = flv%get_mass ()
    part_prop%width = flv%get_width()
    part_prop%spin_type = flv%get_spin_type ()
    part_prop%is_vector = flv%get_spin_type () == VECTOR
    part_prop%empty = .false.
    part_prop%tex_name = flv%get_tex_name ()
    anti = flv%anti ()
    if (flv%get_pdg() == anti%get_pdg()) then
       select type (part_prop)
       type is (part_prop_t)
          part_prop%anti => part_prop
       end select
    else
       do i=1, size (feyngraph_set%particle)
          if (feyngraph_set%particle(i)%pdg == (- part_prop%pdg)) then
             part_prop%anti => feyngraph_set%particle(i)
             exit
          else if (feyngraph_set%particle(i)%empty) then
             part_prop%anti => feyngraph_set%particle(i)
             call feyngraph_set%particle(i)%init &
                  (feyngraph_set, char(anti%get_name()))
             exit
          end if
       end do
    end if
  end subroutine part_prop_init

  module subroutine f_node_assign_particle_properties (node, feyngraph_set)
    class(f_node_t), intent(inout ) :: node
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    character(len=LABEL_LEN) :: particle_label
    integer :: i
    particle_label = node%particle_label(1:index (node%particle_label, '[')-1)
    if (.not. associated (feyngraph_set%particle)) then
       allocate (feyngraph_set%particle (PRT_ARRAY_SIZE))
    end if
    do i = 1, size (feyngraph_set%particle)
       if (particle_label == feyngraph_set%particle(i)%particle_label) then
          node%particle => feyngraph_set%particle(i)
          exit
       else if (feyngraph_set%particle(i)%empty) then
          call feyngraph_set%particle(i)%init (feyngraph_set, particle_label)
          node%particle => feyngraph_set%particle(i)
          exit
       end if
    end do
    !!! Since the O'Mega output uses the anti-particles instead of the
    !!! particles specified in the process definition, we revert this
    !!! here. An exception is the first particle in the parsable DAG output
    node%particle => node%particle%anti
  end subroutine f_node_assign_particle_properties

  function get_n_daughters (subtree_string, pos_first_colon) &
       result (n_daughters)
    character(len=*), intent(in) :: subtree_string
    integer, intent(in) :: pos_first_colon
    integer :: n_daughters
    integer :: n_open_par
    integer :: i
    n_open_par = 1
    n_daughters = 0
    if (len_trim(subtree_string) > 0) then
       if (pos_first_colon > 0) then
          do i=pos_first_colon, len_trim(subtree_string)
             if (subtree_string(i:i) == ',') then
                if (n_open_par == 1) n_daughters = n_daughters + 1
             else if (subtree_string(i:i) == '(') then
                n_open_par = n_open_par + 1
             else if (subtree_string(i:i) == ')') then
                n_open_par = n_open_par - 1
             end if
          end do
          if (n_open_par == 0) then
             n_daughters = n_daughters + 1
          end if
       end if
    end if
  end function get_n_daughters

  recursive subroutine node_construct_subtree_rec (feyngraph_set, &
       feyngraph, subtree_string, mother_node)
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(feyngraph_t), intent(inout) :: feyngraph
    character(len=*), intent(in) :: subtree_string
    type(f_node_t), pointer, intent(inout) :: mother_node
    integer :: n_daughters
    integer :: pos_first_colon
    integer :: current_daughter
    integer :: pos_subtree_begin, pos_subtree_end
    integer :: i
    integer :: n_open_par
    if (.not. associated (mother_node)) then
       call feyngraph_set%f_node_list%add_entry (subtree_string, mother_node, .true.)
       current_daughter = 1
       n_open_par = 1
       pos_first_colon = index (subtree_string, ':')
       n_daughters = get_n_daughters (subtree_string, pos_first_colon)
       if (pos_first_colon == 0) then
          mother_node%particle_label = subtree_string
       else
          mother_node%particle_label = subtree_string(2:pos_first_colon-1)
       end if
       if (.not. associated (mother_node%particle)) then
          call mother_node%assign_particle_properties (feyngraph_set)
       end if
       if (n_daughters /= 2 .and. n_daughters /= 0) then
          mother_node%keep = .false.
          feyngraph%keep = .false.
          return
       end if
       pos_subtree_begin = pos_first_colon + 1
       do i = pos_first_colon + 1, len(trim(subtree_string))
          if (current_daughter == 2) then
             pos_subtree_end = len(trim(subtree_string)) - 1
             call node_construct_subtree_rec (feyngraph_set, feyngraph, &
                  subtree_string(pos_subtree_begin:pos_subtree_end), &
                  mother_node%daughter2)
             exit
          else if (subtree_string(i:i) == ',') then
             if (n_open_par == 1) then
                pos_subtree_end = i - 1
                call node_construct_subtree_rec (feyngraph_set, feyngraph, &
                     subtree_string(pos_subtree_begin:pos_subtree_end), &
                     mother_node%daughter1)
                current_daughter = 2
                pos_subtree_begin = i + 1
             end if
          else if (subtree_string(i:i) == '(') then
             n_open_par = n_open_par + 1
          else if (subtree_string(i:i) == ')') then
             n_open_par = n_open_par - 1
          end if
       end do
    end if
    if (associated (mother_node%daughter1)) then
       if (.not. mother_node%daughter1%keep) then
          mother_node%keep = .false.
       end if
    end if
    if (associated (mother_node%daughter2)) then
       if (.not. mother_node%daughter2%keep) then
          mother_node%keep = .false.
       end if
    end if
    if (associated (mother_node%daughter1) .and. &
         associated (mother_node%daughter2)) then
       mother_node%n_subtree_nodes = &
            mother_node%daughter1%n_subtree_nodes &
            + mother_node%daughter2%n_subtree_nodes + 1
    end if
    if (.not. mother_node%keep) then
       feyngraph%keep = .false.
    end if
  end subroutine node_construct_subtree_rec

  subroutine feyngraph_construct (feyngraph_set, feyngraph)
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(feyngraph_t), pointer, intent(inout) :: feyngraph
    call node_construct_subtree_rec (feyngraph_set, feyngraph, &
         char(feyngraph%omega_feyngraph_output), feyngraph%root)
    feyngraph%n_nodes = feyngraph%root%n_subtree_nodes
  end subroutine feyngraph_construct

  module subroutine dag_node_final (dag_node)
    class(dag_node_t), intent(inout) :: dag_node
    integer :: i
    call dag_node%string%final ()
    if (allocated (dag_node%f_node)) then
       do i=1, size (dag_node%f_node)
          if (associated (dag_node%f_node(i)%node)) then
             call dag_node%f_node(i)%node%final ()
             deallocate (dag_node%f_node(i)%node)
          end if
       end do
       deallocate (dag_node%f_node)
    end if
  end subroutine dag_node_final

  module subroutine dag_options_final (dag_options)
    class(dag_options_t), intent(inout) :: dag_options
    integer :: i
    call dag_options%string%final ()
    if (allocated (dag_options%f_node_ptr1)) then
       do i=1, size (dag_options%f_node_ptr1)
          dag_options%f_node_ptr1(i)%node => null ()
       end do
       deallocate (dag_options%f_node_ptr1)
    end if
        if (allocated (dag_options%f_node_ptr2)) then
       do i=1, size (dag_options%f_node_ptr2)
          dag_options%f_node_ptr2(i)%node => null ()
       end do
       deallocate (dag_options%f_node_ptr2)
    end if
  end subroutine dag_options_final

  module subroutine dag_combination_final (dag_combination)
    class(dag_combination_t), intent(inout) :: dag_combination
    integer :: i
    call dag_combination%string%final ()
    if (allocated (dag_combination%f_node_ptr1)) then
       do i=1, size (dag_combination%f_node_ptr1)
          dag_combination%f_node_ptr1(i)%node => null ()
       end do
       deallocate (dag_combination%f_node_ptr1)
    end if
    if (allocated (dag_combination%f_node_ptr2)) then
       do i=1, size (dag_combination%f_node_ptr2)
          dag_combination%f_node_ptr2(i)%node => null ()
       end do
       deallocate (dag_combination%f_node_ptr2)
    end if
  end subroutine dag_combination_final

  module subroutine dag_final (dag)
    class(dag_t), intent(inout) :: dag
    integer :: i
    call dag%string%final ()
    if (allocated (dag%node)) then
       do i=1, size (dag%node)
          call dag%node(i)%final ()
       end do
       deallocate (dag%node)
    end if
    if (allocated (dag%options)) then
       do i=1, size (dag%options)
          call dag%options(i)%final ()
       end do
       deallocate (dag%options)
    end if
    if (allocated (dag%combination)) then
       do i=1, size (dag%combination)
          call dag%combination(i)%final ()
       end do
       deallocate (dag%combination)
    end if
  end subroutine dag_final

  module subroutine dag_construct (dag, feyngraph_set)
    class(dag_t), intent(inout) :: dag
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    integer :: n_nodes
    integer :: n_options
    integer :: n_combinations
    logical :: continue_loop
    integer :: subtree_size
    integer :: i,j
    subtree_size = 1
    call dag%get_nodes_and_combinations (leaves = .true.)
    do i=1, dag%n_nodes
       call dag%node(i)%make_f_nodes (feyngraph_set, dag)
    end do
    continue_loop = .true.
    subtree_size = subtree_size + 2
    do while (continue_loop)
       n_nodes = dag%n_nodes
       n_options = dag%n_options
       n_combinations = dag%n_combinations
       call dag%get_nodes_and_combinations (leaves = .false.)
       if (n_nodes /= dag%n_nodes) then
          dag%node(n_nodes+1:dag%n_nodes)%subtree_size = subtree_size
          do i = n_nodes+1, dag%n_nodes
             call dag%node(i)%make_f_nodes (feyngraph_set, dag)
          end do
          subtree_size = subtree_size + 2
       end if
       if (n_combinations /= dag%n_combinations) then
          !$OMP PARALLEL DO
          do i = n_combinations+1, dag%n_combinations
             call dag%combination(i)%make_f_nodes (feyngraph_set, dag)
          end do
          !$OMP END PARALLEL DO
       end if
       call dag%get_options ()
       if (n_options /= dag%n_options) then
          !$OMP PARALLEL DO
          do i = n_options+1, dag%n_options
             call dag%options(i)%make_f_nodes (feyngraph_set, dag)
          end do
          !$OMP END PARALLEL DO
       end if
       if (n_nodes == dag%n_nodes .and. n_options == dag%n_options &
            .and. n_combinations == dag%n_combinations) then
          continue_loop = .false.
       end if
    end do
!!! add root node to dag
    call dag%add_node (dag%string%t, leaf = .false.)
    dag%node(dag%n_nodes)%subtree_size = subtree_size
    call dag%node(dag%n_nodes)%make_f_nodes (feyngraph_set, dag)
    if (debug2_active (D_PHASESPACE)) then
       call dag%write (output_unit)
    end if
!!! set indices for all f_nodes
    do i=1, dag%n_nodes
       if (allocated (dag%node(i)%f_node)) then
          do j=1, size (dag%node(i)%f_node)
             if (associated (dag%node(i)%f_node(j)%node)) &
                  call dag%node(i)%f_node(j)%node%set_index ()
          end do
       end if
    end do
  end subroutine dag_construct

  module subroutine dag_get_nodes_and_combinations (dag, leaves)
    class(dag_t), intent(inout) :: dag
    logical, intent(in) :: leaves
    type(dag_string_t) :: new_string
    integer :: i, j, k
    integer :: i_node
    integer :: new_size
    integer :: first_colon
    logical :: combination
    !!! Create nodes also for external particles, except for the incoming one
    !!! which appears as the root of the tree. These can easily be identified
    !!! by their bincodes, since they should contain only one bit which is set.
    if (leaves) then
       first_colon = &
            minloc (dag%string%t%type, 1, dag%string%t%type == COLON_TK)
       do i = first_colon + 1, size (dag%string%t)
          if (dag%string%t(i)%type == NODE_TK) then
             if (popcnt(dag%string%t(i)%bincode) == 1) then
                call dag%add_node (dag%string%t(i:i), .true., i_node)
                call dag%string%t(i)%init_dag_object_token (DAG_NODE_TK, i_node)
             end if
          end if
       end do
       call dag%string%update_char_len ()
    else
    !!! Create a node or combination for every closed pair of parentheses
    !!! which do not contain any other parentheses or curly braces.
    !!! A node (not outgoing) contains a colon. This is not the case
    !!! for combinations, which we use as the criteria to distinguish
    !!! between both.
       allocate (new_string%t (size (dag%string%t)))
       i = 1
       new_size = 0
       do while (i <= size(dag%string%t))
          if (dag%string%t(i)%type == OPEN_PAR_TK) then
             combination = .true.
             do j = i+1, size (dag%string%t)
                select case (dag%string%t(j)%type)
                case (CLOSED_PAR_TK)
                   new_size = new_size + 1
                   if (combination) then
                      call dag%add_combination (dag%string%t(i:j), i_node)
                      call new_string%t(new_size)%init_dag_object_token &
                           (DAG_COMBINATION_TK, i_node)
                   else
                      call dag%add_node (dag%string%t(i:j), leaves, i_node)
                      call new_string%t(new_size)%init_dag_object_token &
                           (DAG_NODE_TK, i_node)
                   end if
                   i = j + 1
                   exit
                case (OPEN_PAR_TK, OPEN_CURLY_TK, CLOSED_CURLY_TK)
                   new_size = new_size + 1
                   new_string%t(new_size) = dag%string%t(i)
                   i = i + 1
                   exit
                case (COLON_TK)
                   combination = .false.
                end select
             end do
          else
             new_size = new_size + 1
             new_string%t(new_size) = dag%string%t(i)
             i = i + 1
          end if
       end do
       dag%string = new_string%t(:new_size)
       call dag%string%update_char_len ()
    end if
  end subroutine dag_get_nodes_and_combinations

  module subroutine dag_get_options (dag)
    class(dag_t), intent(inout) :: dag
    type(dag_string_t) :: new_string
    integer :: i, j, k
    integer :: new_size
    integer :: i_options
    character(len=10) :: index_char
    integer :: index_start, index_end
    !!! Create a node or combination for every closed pair of parentheses
    !!! which do not contain any other parentheses or curly braces.
    !!! A node (not outgoing) contains a colon. This is not the case
    !!! for combinations, which we use as the criteria to distinguish
    !!! between both.
    allocate (new_string%t (size (dag%string%t)))
    i = 1
    new_size = 0
    do while (i <= size(dag%string%t))
       if (dag%string%t(i)%type == OPEN_CURLY_TK) then
          do j = i+1, size (dag%string%t)
             select case (dag%string%t(j)%type)
             case (CLOSED_CURLY_TK)
                new_size = new_size + 1
                call dag%add_options (dag%string%t(i:j), i_options)
                call new_string%t(new_size)%init_dag_object_token (DAG_OPTIONS_TK, i_options)
                i = j + 1
                exit
             case (OPEN_PAR_TK, CLOSED_PAR_TK, OPEN_CURLY_TK)
                new_size = new_size + 1
                new_string%t(new_size) = dag%string%t(i)
                i = i + 1
                exit
             end select
          end do
       else
          new_size = new_size + 1
          new_string%t(new_size) = dag%string%t(i)
          i = i + 1
       end if
    end do
    dag%string = new_string%t(:new_size)
    call dag%string%update_char_len ()
  end subroutine dag_get_options

  module subroutine dag_add_node (dag, string, leaf, i_node)
    class(dag_t), intent(inout) :: dag
    type(dag_token_t), dimension (:), intent(in) :: string
    logical, intent(in) :: leaf
    integer, intent(out), optional :: i_node
    type(dag_node_t), dimension (:), allocatable :: tmp_node
    integer :: string_len
    integer :: i
    string_len = sum (string%char_len)
    if (.not. allocated (dag%node)) then
        allocate (dag%node (DAG_STACK_SIZE))
     else if (dag%n_nodes == size (dag%node)) then
        allocate (tmp_node (dag%n_nodes))
        tmp_node = dag%node
        deallocate (dag%node)
        allocate (dag%node (dag%n_nodes+DAG_STACK_SIZE))
        dag%node(:dag%n_nodes) = tmp_node
        deallocate (tmp_node)
     end if
     do i = 1, dag%n_nodes
        if (dag%node(i)%string_len == string_len) then
           if (size (dag%node(i)%string%t) == size (string)) then
              if (all(dag%node(i)%string%t == string)) then
                 if (present (i_node)) i_node = i
                 return
              end if
           end if
        end if
     end do
     dag%n_nodes = dag%n_nodes + 1
     dag%node(dag%n_nodes)%string = string
     dag%node(dag%n_nodes)%string_len = string_len
     if (present (i_node)) i_node = dag%n_nodes
     dag%node(dag%n_nodes)%leaf = leaf
  end subroutine dag_add_node

  module subroutine dag_add_options (dag, string, i_options)
    class(dag_t), intent(inout) :: dag
    type(dag_token_t), dimension (:), intent(in) :: string
    integer, intent(out), optional :: i_options
    type(dag_options_t), dimension (:), allocatable :: tmp_options
    integer :: string_len
    integer :: i
    string_len = sum (string%char_len)
    if (.not. allocated (dag%options)) then
        allocate (dag%options (DAG_STACK_SIZE))
     else if (dag%n_options == size (dag%options)) then
        allocate (tmp_options (dag%n_options))
        tmp_options = dag%options
        deallocate (dag%options)
        allocate (dag%options (dag%n_options+DAG_STACK_SIZE))
        dag%options(:dag%n_options) = tmp_options
        deallocate (tmp_options)
     end if
     do i = 1, dag%n_options
        if (dag%options(i)%string_len == string_len) then
           if (size (dag%options(i)%string%t) == size (string)) then
              if (all(dag%options(i)%string%t == string)) then
                 if (present (i_options)) i_options = i
                 return
              end if
           end if
        end if
     end do
     dag%n_options = dag%n_options + 1
     dag%options(dag%n_options)%string = string
     dag%options(dag%n_options)%string_len = string_len
     if (present (i_options)) i_options = dag%n_options
  end subroutine dag_add_options

  module subroutine dag_add_combination (dag, string, i_combination)
    class(dag_t), intent(inout) :: dag
    type(dag_token_t), dimension (:), intent(in) :: string
    integer, intent(out), optional :: i_combination
    type(dag_combination_t), dimension (:), allocatable :: tmp_combination
    integer :: string_len
    integer :: i
    string_len = sum (string%char_len)
    if (.not. allocated (dag%combination)) then
        allocate (dag%combination (DAG_STACK_SIZE))
     else if (dag%n_combinations == size (dag%combination)) then
        allocate (tmp_combination (dag%n_combinations))
        tmp_combination = dag%combination
        deallocate (dag%combination)
        allocate (dag%combination (dag%n_combinations+DAG_STACK_SIZE))
        dag%combination(:dag%n_combinations) = tmp_combination
        deallocate (tmp_combination)
     end if
     do i = 1, dag%n_combinations
        if (dag%combination(i)%string_len == string_len) then
           if (size (dag%combination(i)%string%t) == size (string)) then
              if (all(dag%combination(i)%string%t == string)) then
                 i_combination = i
                 return
              end if
           end if
        end if
     end do
     dag%n_combinations = dag%n_combinations + 1
     dag%combination(dag%n_combinations)%string = string
     dag%combination(dag%n_combinations)%string_len = string_len
     if (present (i_combination)) i_combination = dag%n_combinations
  end subroutine dag_add_combination

  module subroutine dag_node_make_f_nodes (dag_node, feyngraph_set, dag)
    class(dag_node_t), intent(inout) :: dag_node
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(dag_t), intent(inout) :: dag
    character(len=LABEL_LEN) :: particle_label
    integer :: i, j
    integer, dimension (2) :: obj
    integer, dimension (2) :: i_obj
    integer :: n_obj
    integer :: pos
    integer :: new_size, size1, size2
    integer, dimension(:), allocatable :: match
    if (allocated (dag_node%f_node)) return
    pos = minloc (dag_node%string%t%type, 1,dag_node%string%t%type == NODE_TK)
    particle_label = char (dag_node%string%t(pos))
    if (dag_node%leaf) then
!!! construct subtree with procedure similar to the one for the old output
       allocate (dag_node%f_node(1))
       allocate (dag_node%f_node(1)%node)
       dag_node%f_node(1)%node%particle_label = particle_label
       call dag_node%f_node(1)%node%assign_particle_properties (feyngraph_set)
       if (.not. dag_node%f_node(1)%node%keep) then
          deallocate (dag_node%f_node)
          return
       end if
    else
       n_obj = 0
       do i = 1, size (dag_node%string%t)
          select case (dag_node%string%t(i)%type)
          case (DAG_NODE_TK, DAG_OPTIONS_TK, DAG_COMBINATION_TK)
             n_obj = n_obj + 1
             if (n_obj > 2) return
             obj(n_obj) = dag_node%string%t(i)%type
             i_obj(n_obj) = dag_node%string%t(i)%index
          end select
       end do
       if (n_obj == 1) then
          if (obj(1) == DAG_OPTIONS_TK) then
             if (allocated (dag%options(i_obj(1))%f_node_ptr1)) then
                size1 = size(dag%options(i_obj(1))%f_node_ptr1)
                allocate (dag_node%f_node(size1))
                do i=1, size1
                   allocate (dag_node%f_node(i)%node)
                   dag_node%f_node(i)%node%particle_label = particle_label
                   call dag_node%f_node(i)%node%assign_particle_properties (feyngraph_set)
                   dag_node%f_node(i)%node%daughter1 => dag%options(i_obj(1))%f_node_ptr1(i)%node
                   dag_node%f_node(i)%node%daughter2 => dag%options(i_obj(1))%f_node_ptr2(i)%node
                   dag_node%f_node(i)%node%n_subtree_nodes = &
                        dag%options(i_obj(1))%f_node_ptr1(i)%node%n_subtree_nodes &
                        + dag%options(i_obj(1))%f_node_ptr2(i)%node%n_subtree_nodes + 1
                end do
             end if
          else if (obj(1) == DAG_COMBINATION_TK) then
             if (allocated (dag%combination(i_obj(1))%f_node_ptr1)) then
                size1 = size(dag%combination(i_obj(1))%f_node_ptr1)
                allocate (dag_node%f_node(size1))
                do i=1, size1
                   allocate (dag_node%f_node(i)%node)
                   dag_node%f_node(i)%node%particle_label = particle_label
                   call dag_node%f_node(i)%node%assign_particle_properties (feyngraph_set)
                   dag_node%f_node(i)%node%daughter1 => dag%combination(i_obj(1))%f_node_ptr1(i)%node
                   dag_node%f_node(i)%node%daughter2 => dag%combination(i_obj(1))%f_node_ptr2(i)%node
                   dag_node%f_node(i)%node%n_subtree_nodes = &
                        dag%combination(i_obj(1))%f_node_ptr1(i)%node%n_subtree_nodes &
                     + dag%combination(i_obj(1))%f_node_ptr2(i)%node%n_subtree_nodes + 1
                end do
             end if
          end if
!!! simply set daughter pointers, daughters are already combined correctly
       else if (n_obj == 2) then
          size1 = 0
          size2 = 0
          if (obj(1) == DAG_NODE_TK) then
             if (allocated (dag%node(i_obj(1))%f_node)) then
                do i=1, size (dag%node(i_obj(1))%f_node)
                   if (dag%node(i_obj(1))%f_node(i)%node%keep) size1 = size1 + 1
                end do
             end if
          else if (obj(1) == DAG_OPTIONS_TK) then
             if (allocated (dag%options(i_obj(1))%f_node_ptr1)) then
                do i=1, size (dag%options(i_obj(1))%f_node_ptr1)
                   if (dag%options(i_obj(1))%f_node_ptr1(i)%node%keep) size1 = size1 + 1
                end do
             end if
          end if
          if (obj(2) == DAG_NODE_TK) then
             if (allocated (dag%node(i_obj(2))%f_node)) then
                do i=1, size (dag%node(i_obj(2))%f_node)
                   if (dag%node(i_obj(2))%f_node(i)%node%keep) size2 = size2 + 1
                end do
             end if
          else if (obj(2) == DAG_OPTIONS_TK) then
             if (allocated (dag%options(i_obj(2))%f_node_ptr1)) then
                do i=1, size (dag%options(i_obj(2))%f_node_ptr1)
                   if (dag%options(i_obj(2))%f_node_ptr1(i)%node%keep) size2 = size2 + 1
                end do
             end if
          end if
!!! make all combinations of daughters
          select case (obj(1))
          case (DAG_NODE_TK)
             select case (obj(2))
             case (DAG_NODE_TK)
                call combine_all_daughters(dag%node(i_obj(1))%f_node, &
                     dag%node(i_obj(2))%f_node)
             case (DAG_OPTIONS_TK)
                call combine_all_daughters(dag%node(i_obj(1))%f_node, &
                     dag%options(i_obj(2))%f_node_ptr1)
             end select
          case (DAG_OPTIONS_TK)
             select case (obj(2))
             case (DAG_NODE_TK)
                call combine_all_daughters(dag%options(i_obj(1))%f_node_ptr1, &
                     dag%node(i_obj(2))%f_node)
             case (DAG_OPTIONS_TK)
                call combine_all_daughters(dag%options(i_obj(1))%f_node_ptr1, &
                     dag%options(i_obj(2))%f_node_ptr1)
             end select
          end select
       end if
    end if

  contains

    subroutine combine_all_daughters (daughter1_ptr, daughter2_ptr)
      type(f_node_ptr_t), dimension (:), intent(in) :: daughter1_ptr
      type(f_node_ptr_t), dimension (:), intent(in) :: daughter2_ptr
      integer :: i, j
      integer :: pos
      new_size = size1*size2
      allocate (dag_node%f_node(new_size))
      pos = 0
      do i = 1, size (daughter1_ptr)
         if (daughter1_ptr(i)%node%keep) then
            do j = 1, size (daughter2_ptr)
               if (daughter2_ptr(j)%node%keep) then
                  pos = pos + 1
                  allocate (dag_node%f_node(pos)%node)
                  dag_node%f_node(pos)%node%particle_label = particle_label
                  call dag_node%f_node(pos)%node%assign_particle_properties (feyngraph_set)
                  dag_node%f_node(pos)%node%daughter1 => daughter1_ptr(i)%node
                  dag_node%f_node(pos)%node%daughter2 => daughter2_ptr(j)%node
                  dag_node%f_node(pos)%node%n_subtree_nodes = daughter1_ptr(i)%node%n_subtree_nodes &
                       + daughter2_ptr(j)%node%n_subtree_nodes + 1
                  call feyngraph_set%model%match_vertex (daughter1_ptr(i)%node%particle%pdg, &
                       daughter2_ptr(j)%node%particle%pdg, match)
                  if (allocated (match)) then
                     if (any (abs(match) == abs(dag_node%f_node(pos)%node%particle%pdg))) then
                        dag_node%f_node(pos)%node%keep = .true.
                     else
                        dag_node%f_node(pos)%node%keep = .false.
                     end if
                     deallocate (match)
                  else
                     dag_node%f_node(pos)%node%keep = .false.
                  end if
               end if
            end do
         end if
      end do
    end subroutine combine_all_daughters
  end subroutine dag_node_make_f_nodes

  module subroutine dag_options_make_f_nodes (dag_options, &
       feyngraph_set, dag)
    class(dag_options_t), intent(inout) :: dag_options
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(dag_t), intent(inout) :: dag
    integer, dimension (:), allocatable :: obj, i_obj
    integer :: n_obj
    integer :: i
    integer :: pos
!!! read options
    if (allocated (dag_options%f_node_ptr1)) return
    n_obj = count ((dag_options%string%t%type == DAG_NODE_TK) .or. &
         (dag_options%string%t%type == DAG_OPTIONS_TK) .or. &
         (dag_options%string%t%type == DAG_COMBINATION_TK), 1)
    allocate (obj(n_obj)); allocate (i_obj(n_obj))
    pos = 0
    do i = 1, size (dag_options%string%t)
       select case (dag_options%string%t(i)%type)
       case (DAG_NODE_TK, DAG_OPTIONS_TK, DAG_COMBINATION_TK)
          pos = pos + 1
          obj(pos) = dag_options%string%t(i)%type
          i_obj(pos) = dag_options%string%t(i)%index
       end select
    end do
    if (any (dag_options%string%t%type == DAG_NODE_TK)) then
       call dag_options_make_f_nodes_single
    else if (any (dag_options%string%t%type == DAG_COMBINATION_TK)) then
       call dag_options_make_f_nodes_pair
    end if
    deallocate (obj, i_obj)

  contains

    subroutine dag_options_make_f_nodes_single
      integer :: i_start, i_end
      integer :: n_nodes
      n_nodes = 0
      do i=1, n_obj
         if (allocated (dag%node(i_obj(i))%f_node)) then
            n_nodes = n_nodes + size (dag%node(i_obj(i))%f_node)
         end if
      end do
      if (n_nodes /= 0) then
         allocate (dag_options%f_node_ptr1 (n_nodes))
         i_end = 0
         do i = 1, n_obj
            if (allocated (dag%node(i_obj(i))%f_node)) then
               i_start = i_end + 1
               i_end = i_end + size (dag%node(i_obj(i))%f_node)
               dag_options%f_node_ptr1(i_start:i_end) = dag%node(i_obj(i))%f_node
            end if
         end do
      end if
    end subroutine dag_options_make_f_nodes_single

    subroutine dag_options_make_f_nodes_pair
      integer :: i_start, i_end
      integer :: n_nodes
!!! get f_nodes from each combination
      n_nodes = 0
      do i=1, n_obj
         if (allocated (dag%combination(i_obj(i))%f_node_ptr1)) then
            n_nodes = n_nodes + size (dag%combination(i_obj(i))%f_node_ptr1)
         end if
      end do
      if (n_nodes /= 0) then
         allocate (dag_options%f_node_ptr1 (n_nodes))
         allocate (dag_options%f_node_ptr2 (n_nodes))
         i_end = 0
         do i=1, n_obj
            if (allocated (dag%combination(i_obj(i))%f_node_ptr1)) then
               i_start = i_end + 1
               i_end = i_end + size (dag%combination(i_obj(i))%f_node_ptr1)
               dag_options%f_node_ptr1(i_start:i_end) = dag%combination(i_obj(i))%f_node_ptr1
               dag_options%f_node_ptr2(i_start:i_end) = dag%combination(i_obj(i))%f_node_ptr2
            end if
         end do
      end if
    end subroutine dag_options_make_f_nodes_pair
  end subroutine dag_options_make_f_nodes

  module subroutine dag_combination_make_f_nodes (dag_combination, &
       feyngraph_set, dag)
    class(dag_combination_t), intent(inout) :: dag_combination
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(dag_t), intent(inout) :: dag
    integer, dimension (2) :: obj, i_obj
    integer :: n_obj
    integer :: new_size, size1, size2
    integer :: i, j, pos
    if (allocated (dag_combination%f_node_ptr1)) return
    n_obj = 0
    do i = 1, size (dag_combination%string%t)
       select case (dag_combination%string%t(i)%type)
       case (DAG_NODE_TK, DAG_OPTIONS_TK, DAG_COMBINATION_TK)
          n_obj = n_obj + 1
          if (n_obj > 2) return
          obj(n_obj) = dag_combination%string%t(i)%type
          i_obj(n_obj) = dag_combination%string%t(i)%index
       end select
    end do
    size1 = 0
    size2 = 0
    if (obj(1) == DAG_NODE_TK) then
       if (allocated (dag%node(i_obj(1))%f_node)) &
            size1 = size (dag%node(i_obj(1))%f_node)
    else if (obj(1) == DAG_OPTIONS_TK) then
       if (allocated (dag%options(i_obj(1))%f_node_ptr1)) &
            size1 = size (dag%options(i_obj(1))%f_node_ptr1)
    end if
    if (obj(2) == DAG_NODE_TK) then
       if (allocated (dag%node(i_obj(2))%f_node)) &
            size2 = size (dag%node(i_obj(2))%f_node)
    else if (obj(2) == DAG_OPTIONS_TK) then
       if (allocated (dag%options(i_obj(2))%f_node_ptr1)) &
            size2 = size (dag%options(i_obj(2))%f_node_ptr1)
    end if
!!! combine the 2 arrays of f_nodes
    new_size = size1*size2
    if (new_size /= 0) then
       allocate (dag_combination%f_node_ptr1 (new_size))
       allocate (dag_combination%f_node_ptr2 (new_size))
       pos = 0
       select case (obj(1))
       case (DAG_NODE_TK)
          select case (obj(2))
          case (DAG_NODE_TK)
             do i = 1, size1
                do j = 1, size2
                   pos = pos + 1
                   dag_combination%f_node_ptr1(pos) = &
                        dag%node(i_obj(1))%f_node(i)
                   dag_combination%f_node_ptr2(pos) = &
                        dag%node(i_obj(2))%f_node(j)
                end do
             end do
          case (DAG_OPTIONS_TK)
             do i = 1, size1
                do j = 1, size2
                   pos = pos + 1
                   dag_combination%f_node_ptr1(pos) = &
                        dag%node(i_obj(1))%f_node(i)
                   dag_combination%f_node_ptr2(pos) = &
                        dag%options(i_obj(2))%f_node_ptr1(j)
                end do
             end do
          end select
       case (DAG_OPTIONS_TK)
          select case (obj(2))
          case (DAG_NODE_TK)
             do i = 1, size1
                do j = 1, size2
                   pos = pos + 1
                   dag_combination%f_node_ptr1(pos) = &
                        dag%options(i_obj(1))%f_node_ptr1(i)
                   dag_combination%f_node_ptr2(pos) = &
                        dag%node(i_obj(2))%f_node(j)
                end do
             end do
          case (DAG_OPTIONS_TK)
             do i = 1, size1
                do j = 1, size2
                   pos = pos + 1
                   dag_combination%f_node_ptr1(pos) = &
                        dag%options(i_obj(1))%f_node_ptr1(i)
                   dag_combination%f_node_ptr2(pos) = &
                        dag%options(i_obj(2))%f_node_ptr1(j)
                end do
             end do
          end select
       end select
    end if
  end subroutine dag_combination_make_f_nodes

  module subroutine dag_make_feyngraphs (dag, feyngraph_set)
    class(dag_t), intent(inout) :: dag
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    integer :: i
    integer :: max_subtree_size
    max_subtree_size = dag%node(dag%n_nodes)%subtree_size
    if (allocated (dag%node(dag%n_nodes)%f_node)) then
       do i = 1, size (dag%node(dag%n_nodes)%f_node)
          if (.not. associated (feyngraph_set%first)) then
             allocate (feyngraph_set%last)
             feyngraph_set%first => feyngraph_set%last
          else
             allocate (feyngraph_set%last%next)
             feyngraph_set%last => feyngraph_set%last%next
          end if
          feyngraph_set%last%root => dag%node(dag%n_nodes)%f_node(i)%node
          !!! The first particle was correct in the O'Mega parsable DAG output.
          !!! It was however changed to its anti-particle in
          !!! f_node_assign_particle_properties, which we revert here.
          feyngraph_set%last%root%particle => &
               feyngraph_set%last%root%particle%anti
          feyngraph_set%last%n_nodes = feyngraph_set%last%root%n_subtree_nodes
          feyngraph_set%n_graphs = feyngraph_set%n_graphs + 1
       end do
       feyngraph_set%f_node_list%max_tree_size = feyngraph_set%first%n_nodes
    end if
  end subroutine dag_make_feyngraphs

  module subroutine dag_write (dag, u)
    class(dag_t), intent(in) :: dag
    integer, intent(in) :: u
    integer :: i
    write (u,fmt='(A)') 'nodes'
    do i=1, dag%n_nodes
       write (u,fmt='(I5,3X,A)') i, char (dag%node(i)%string)
    end do
    write (u,fmt='(A)') 'options'
    do i=1, dag%n_options
       write (u,fmt='(I5,3X,A)') i, char (dag%options(i)%string)
    end do
    write (u,fmt='(A)') 'combination'
    do i=1, dag%n_combinations
       write (u,fmt='(I5,3X,A)') i, char (dag%combination(i)%string)
    end do
  end subroutine dag_write

  subroutine k_node_make_nonresonant_copy (k_node)
    type(k_node_t), intent(in) :: k_node
    type(k_node_t), pointer :: copy
    call k_node%f_node%k_node_list%add_entry (copy, recycle=.true.)
    copy%daughter1 => k_node%daughter1
    copy%daughter2 => k_node%daughter2
    copy = k_node
    copy%mapping = NONRESONANT
    copy%resonant = .false.
    copy%on_shell = .false.
    copy%mapping_assigned = .true.
    copy%is_nonresonant_copy = .true.
  end subroutine k_node_make_nonresonant_copy

  module subroutine feyngraph_make_kingraphs (feyngraph, feyngraph_set)
    class(feyngraph_t), intent(inout) :: feyngraph
    type(feyngraph_set_t), intent(in) :: feyngraph_set
    type(k_node_ptr_t), dimension (:), allocatable :: kingraph_root
    integer :: i
    if (.not. associated (feyngraph%kin_first)) then
       call k_node_init_from_f_node (feyngraph%root, &
            kingraph_root, feyngraph_set)
       if (.not. feyngraph%root%keep) return
       if (feyngraph_set%process_type == SCATTERING) then
          call split_up_t_lines (kingraph_root)
       end if
       do i=1, size (kingraph_root)
          if (associated (feyngraph%kin_last)) then
             allocate (feyngraph%kin_last%next)
             feyngraph%kin_last => feyngraph%kin_last%next
          else
             allocate (feyngraph%kin_last)
             feyngraph%kin_first => feyngraph%kin_last
          end if
          feyngraph%kin_last%root => kingraph_root(i)%node
          feyngraph%kin_last%n_nodes = feyngraph%n_nodes
          feyngraph%kin_last%keep = feyngraph%keep
          if (feyngraph_set%process_type == SCATTERING) then
             feyngraph%kin_last%root%bincode = &
                  f_node_get_external_bincode (feyngraph_set, feyngraph%root)
          end if
       end do
       deallocate (kingraph_root)
    end if
  end subroutine feyngraph_make_kingraphs

  recursive subroutine k_node_init_from_f_node (f_node, k_node_ptr, feyngraph_set)
    type(f_node_t), target, intent(inout) :: f_node
    type(k_node_ptr_t), allocatable, dimension (:), intent(out) :: k_node_ptr
    type(feyngraph_set_t), intent(in) :: feyngraph_set
    type(k_node_ptr_t), allocatable, dimension(:) :: daughter_ptr1, daughter_ptr2
    integer :: n_nodes
    integer :: i, j
    integer :: pos
    integer, save :: counter = 0
    if (.not. (f_node%incoming .or. f_node%t_line)) then
       call f_node%k_node_list%get_nodes (k_node_ptr)
       if (.not. allocated (k_node_ptr) .and. f_node%k_node_list%n_entries > 0) then
          f_node%keep = .false.
          return
       end if
    end if
    if (.not. allocated (k_node_ptr)) then
       if (associated (f_node%daughter1) .and. associated (f_node%daughter2)) then
          call k_node_init_from_f_node (f_node%daughter1, daughter_ptr1, &
               feyngraph_set)
          call k_node_init_from_f_node (f_node%daughter2, daughter_ptr2, &
               feyngraph_set)
          if (.not. (f_node%daughter1%keep .and. f_node%daughter2%keep)) then
             f_node%keep = .false.
             return
          end if
          n_nodes = size (daughter_ptr1) * size (daughter_ptr2)
          allocate (k_node_ptr (n_nodes))
          pos = 1
          do i=1, size (daughter_ptr1)
             do j=1, size (daughter_ptr2)
                if (f_node%incoming .or. f_node%t_line) then
                   call f_node%k_node_list%add_entry (k_node_ptr(pos)%node, recycle = .false.)
                else
                   call f_node%k_node_list%add_entry (k_node_ptr(pos)%node, recycle = .true.)
                end if
                k_node_ptr(pos)%node%f_node => f_node
                k_node_ptr(pos)%node%daughter1 => daughter_ptr1(i)%node
                k_node_ptr(pos)%node%daughter2 => daughter_ptr2(j)%node
                k_node_ptr(pos)%node%f_node_index = f_node%index
                k_node_ptr(pos)%node%incoming = f_node%incoming
                k_node_ptr(pos)%node%t_line = f_node%t_line
                k_node_ptr(pos)%node%particle => f_node%particle
                pos = pos + 1
             end do
          end do
          deallocate (daughter_ptr1, daughter_ptr2)
       else
          allocate (k_node_ptr(1))
          if (f_node%incoming .or. f_node%t_line) then
             call f_node%k_node_list%add_entry (k_node_ptr(1)%node, recycle=.false.)
          else
             call f_node%k_node_list%add_entry (k_node_ptr(1)%node, recycle=.true.)
          end if
          k_node_ptr(1)%node%f_node => f_node
          k_node_ptr(1)%node%f_node_index = f_node%index
          k_node_ptr(1)%node%incoming = f_node%incoming
          k_node_ptr(1)%node%t_line = f_node%t_line
          k_node_ptr(1)%node%particle => f_node%particle
          k_node_ptr(1)%node%bincode = f_node_get_external_bincode (feyngraph_set, &
               f_node)
       end if
    end if
  end subroutine k_node_init_from_f_node

  recursive subroutine split_up_t_lines (t_node)
    type(k_node_ptr_t), dimension(:), intent(inout) :: t_node
    type(k_node_t), pointer :: ref_node => null ()
    type(k_node_t), pointer :: ref_daughter => null ()
    type(k_node_t), pointer :: new_daughter => null ()
    type(k_node_ptr_t), dimension(:), allocatable :: t_daughter
    integer :: ref_daughter_index
    integer :: i, j
    allocate (t_daughter (size (t_node)))
    do i=1, size (t_node)
       ref_node => t_node(i)%node
       if (associated (ref_node%daughter1) .and. associated (ref_node%daughter2)) then
          ref_daughter => null ()
          if (ref_node%daughter1%incoming .or. ref_node%daughter1%t_line) then
             ref_daughter => ref_node%daughter1
             ref_daughter_index = 1
          else if (ref_node%daughter2%incoming .or. ref_node%daughter2%t_line) then
             ref_daughter => ref_node%daughter2
             ref_daughter_index = 2
          end if
          do j=1, size (t_daughter)
             if (.not. associated (t_daughter(j)%node)) then
                t_daughter(j)%node => ref_daughter
                exit
             else if (t_daughter(j)%node%index == ref_daughter%index) then
                new_daughter => null ()
                call ref_daughter%f_node%k_node_list%add_entry (new_daughter, recycle=.false.)
                new_daughter = ref_daughter
                new_daughter%daughter1 => ref_daughter%daughter1
                new_daughter%daughter2 => ref_daughter%daughter2
                if (ref_daughter_index == 1) then
                   ref_node%daughter1 => new_daughter
                else if (ref_daughter_index == 2) then
                   ref_node%daughter2 => new_daughter
                end if
                ref_daughter => new_daughter
             end if
          end do
       else
          return
       end if
    end do
    call split_up_t_lines (t_daughter)
    deallocate (t_daughter)
  end subroutine split_up_t_lines

  subroutine kingraph_set_inverse_daughters (kingraph)
    type(kingraph_t), intent(inout) :: kingraph
    type(k_node_t), pointer :: mother
    type(k_node_t), pointer :: t_daughter
    type(k_node_t), pointer :: s_daughter
    mother => kingraph%root
    do while (associated (mother))
       if (associated (mother%daughter1) .and. &
            associated (mother%daughter2)) then
          if (mother%daughter1%t_line .or. mother%daughter1%incoming) then
             t_daughter => mother%daughter1; s_daughter => mother%daughter2
          else if (mother%daughter2%t_line .or. mother%daughter2%incoming) then
             t_daughter => mother%daughter2; s_daughter => mother%daughter1
          else
             exit
          end if
          t_daughter%inverse_daughter1 => mother
          t_daughter%inverse_daughter2 => s_daughter
          mother => t_daughter
       else
          exit
       end if
    end do
  end subroutine kingraph_set_inverse_daughters

  function f_node_get_external_bincode (feyngraph_set, f_node) result (bincode)
    type(feyngraph_set_t), intent(in) :: feyngraph_set
    type(f_node_t), intent(in) :: f_node
    integer (TC) :: bincode
    character(len=LABEL_LEN) :: particle_label
    integer :: start_pos, end_pos, n_out_decay
    integer :: n_prt ! for DAG
    integer :: i
    bincode = 0
    if (feyngraph_set%process_type == DECAY) then
       n_out_decay = feyngraph_set%n_out
    else
       n_out_decay = feyngraph_set%n_out + 1
    end if
    particle_label = f_node%particle_label
    start_pos = index (particle_label, '[') + 1
    end_pos = index (particle_label, ']') - 1
    particle_label = particle_label(start_pos:end_pos)
!!! n_out_decay is the number of outgoing particles in the
!!! O'Mega output, which is always represented as a decay
    if (feyngraph_set%use_dag) then
       n_prt = 1
       do i=1, len(particle_label)
          if (particle_label(i:i) == '/') n_prt = n_prt + 1
       end do
    else
       n_prt = end_pos - start_pos + 1
    end if
    if (n_prt == 1) then
       bincode = calculate_external_bincode (particle_label, &
            feyngraph_set%process_type, n_out_decay)
    else if (n_prt == n_out_decay) then
       bincode = ibset (0, n_out_decay)
    end if
  end function f_node_get_external_bincode

  subroutine node_assign_bincode (node)
    type(k_node_t), intent(inout) :: node
    if (associated (node%daughter1) .and. associated (node%daughter2) &
         .and. .not. node%incoming) then
       node%bincode = ior(node%daughter1%bincode, node%daughter2%bincode)
    end if
  end subroutine node_assign_bincode

  function calculate_external_bincode (label_number_string, process_type, n_out_decay) result (bincode)
    character(len=*), intent(in) :: label_number_string
    integer, intent(in) :: process_type
    integer, intent(in) :: n_out_decay
    character :: number_char
    integer :: number_int
    integer (kind=TC) :: bincode
    bincode = 0
    read (label_number_string, fmt='(A)') number_char
!!! check if the character is a letter (A,B,C,...) or a number (1...9)
!!! numbers 1 and 2 are special cases
    select case (number_char)
    case ('1')
       if (process_type == SCATTERING) then
          number_int = n_out_decay + 3
       else
          number_int = n_out_decay + 2
       end if
    case ('2')
       if (process_type == SCATTERING) then
          number_int = n_out_decay + 2
       else
          number_int = 2
       end if
    case ('A')
       number_int = 10
    case ('B')
       number_int = 11
    case ('C')
       number_int = 12
    case ('D')
       number_int = 13
    case default
       read (number_char, fmt='(I1)') number_int
    end select
    bincode = ibset (bincode, number_int - process_type - 1)
  end function calculate_external_bincode

  subroutine node_assign_mapping_s (feyngraph, node, feyngraph_set)
    type(feyngraph_t), intent(inout) :: feyngraph
    type(k_node_t), intent(inout) :: node
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    real(default) :: eff_mass_sum
    logical :: keep
    if (.not. node%mapping_assigned) then
       if (node%particle%mass > feyngraph_set%phs_par%m_threshold_s) then
          node%effective_mass = node%particle%mass
       end if
       if (associated (node%daughter1) .and. associated (node%daughter2)) then
          if (.not. (node%daughter1%keep .and. node%daughter2%keep)) then
             node%keep = .false.; return
          end if
          node%ext_mass_sum = node%daughter1%ext_mass_sum &
               + node%daughter2%ext_mass_sum
          keep = .false.
!!! Potentially resonant cases [sqrts = m_rea for on-shell decay]
          if (node%particle%mass > node%ext_mass_sum &
               .and. node%particle%mass <= feyngraph_set%phs_par%sqrts) then
             if (node%particle%width /= 0) then
                if (node%daughter1%on_shell .or. node%daughter2%on_shell) then
                   keep = .true.
                   node%mapping = S_CHANNEL
                   node%resonant = .true.
                end if
             else
                call warn_decay (node%particle)
             end if
!!! Collinear and IR singular cases
          else if (node%particle%mass < feyngraph_set%phs_par%sqrts) then
!!! Massless splitting
             if (node%daughter1%effective_mass == 0 &
                  .and. node%daughter2%effective_mass == 0 &
                  .and. .not. associated (node%daughter1%daughter1) &
                  .and. .not. associated (node%daughter1%daughter2) &
                  .and. .not. associated (node%daughter2%daughter1) &
                  .and. .not. associated (node%daughter2%daughter2)) then
                keep = .true.
                node%log_enhanced = .true.
                if (node%particle%is_vector) then
                   if (node%daughter1%particle%is_vector &
                        .and. node%daughter2%particle%is_vector) then
                      node%mapping = COLLINEAR   !!! three-vector-splitting
                   else
                      node%mapping = INFRARED    !!! vector spliiting into matter
                   end if
                else
                   if (node%daughter1%particle%is_vector &
                        .or. node%daughter2%particle%is_vector) then
                      node%mapping = COLLINEAR   !!! vector radiation off matter
                   else
                      node%mapping = INFRARED    !!! scalar radiation/splitting
                   end if
                end if
!!! IR radiation off massive particle [cascades]
             else if (node%effective_mass > 0 .and. &
                  node%daughter1%effective_mass > 0 .and. &
                  node%daughter2%effective_mass == 0 .and. &
                  (node%daughter1%on_shell .or. &
                  node%daughter1%mapping == RADIATION) .and. &
                  abs (node%effective_mass - &
                  node%daughter1%effective_mass) < feyngraph_set%phs_par%m_threshold_s) &
                  then
                keep = .true.
                node%log_enhanced = .true.
                node%mapping = RADIATION
             else if (node%effective_mass > 0 .and. &
                  node%daughter2%effective_mass > 0 .and. &
                  node%daughter1%effective_mass == 0 .and. &
                  (node%daughter2%on_shell .or. &
                  node%daughter2%mapping == RADIATION) .and. &
                  abs (node%effective_mass - &
                  node%daughter2%effective_mass) < feyngraph_set%phs_par%m_threshold_s) &
                  then
                keep = .true.
                node%log_enhanced = .true.
                node%mapping = RADIATION
             end if
          end if
!!! Non-singular cases, including failed resonances [from cascades]
          if (.not. keep) then
!!! Two on-shell particles from a virtual mother [from cascades, here eventually more than 2]
             if (node%daughter1%on_shell .or. node%daughter2%on_shell) then
                keep = .true.
                eff_mass_sum = node%daughter1%effective_mass &
                     + node%daughter2%effective_mass
                node%effective_mass = max (node%ext_mass_sum, eff_mass_sum)
                if (node%effective_mass < feyngraph_set%phs_par%m_threshold_s) then
                   node%effective_mass = 0
                end if
             end if
          end if
!!! Complete and register feyngraph (make copy in case of resonance)
          if (keep) then
             node%on_shell = node%resonant .or. node%log_enhanced
             if (node%resonant) then
                if (feyngraph_set%phs_par%keep_nonresonant) then
                   call k_node_make_nonresonant_copy (node)
                end if
                node%ext_mass_sum = node%particle%mass
             end if
          end if
          node%mapping_assigned = .true.
          call node_assign_bincode (node)
          call node%subtree%add_entry (node)
       else !!! external (outgoing) particle
          node%ext_mass_sum = node%particle%mass
          node%mapping = EXTERNAL_PRT
          node%multiplicity = 1
          node%mapping_assigned = .true.
          call node%subtree%add_entry (node)
          node%on_shell = .true.
          if (node%particle%mass >= feyngraph_set%phs_par%m_threshold_s) then
             node%effective_mass = node%particle%mass
          end if
       end if
    else if (node%is_nonresonant_copy) then
       call node_assign_bincode (node)
       call node%subtree%add_entry (node)
       node%is_nonresonant_copy = .false.
    end if
    call node_count_specific_properties (node)
    if (node%n_off_shell > feyngraph_set%phs_par%off_shell) then
       node%keep = .false.
    end if
  contains
    subroutine warn_decay (particle)
      type(part_prop_t), intent(in) :: particle
      integer :: i
      integer, dimension(MAX_WARN_RESONANCE), save :: warned_code = 0
      LOOP_WARNED: do i = 1, MAX_WARN_RESONANCE
         if (warned_code(i) == 0) then
            warned_code(i) = particle%pdg
            write (msg_buffer, "(A)") &
                 & " Intermediate decay of zero-width particle " &
                 & // trim(particle%particle_label) &
                 & // " may be possible."
            call msg_warning
            exit LOOP_WARNED
         else if (warned_code(i) == particle%pdg) then
            exit LOOP_WARNED
         end if
      end do LOOP_WARNED
    end subroutine warn_decay
  end subroutine node_assign_mapping_s

  subroutine node_count_specific_properties (node)
    type(k_node_t), intent(inout) :: node
    if (associated (node%daughter1) .and. associated(node%daughter2)) then
       if (node%resonant) then
          node%multiplicity = 1
          node%n_resonances &
               = node%daughter1%n_resonances &
               + node%daughter2%n_resonances + 1
       else
          node%multiplicity &
               = node%daughter1%multiplicity &
               + node%daughter2%multiplicity
          node%n_resonances &
               = node%daughter1%n_resonances &
               + node%daughter2%n_resonances
       end if
       if (node%log_enhanced) then
          node%n_log_enhanced &
               = node%daughter1%n_log_enhanced &
               + node%daughter2%n_log_enhanced + 1
       else
          node%n_log_enhanced &
               = node%daughter1%n_log_enhanced &
               + node%daughter2%n_log_enhanced
       end if
       if (node%resonant) then
          node%n_off_shell = 0
       else if (node%log_enhanced) then
          node%n_off_shell &
               = node%daughter1%n_off_shell &
               + node%daughter2%n_off_shell
       else
          node%n_off_shell &
               = node%daughter1%n_off_shell &
               + node%daughter2%n_off_shell + 1
       end if
       if (node%t_line) then
          if (node%daughter1%t_line .or. node%daughter1%incoming) then
             node%n_t_channel = node%daughter1%n_t_channel + 1
          else if (node%daughter2%t_line .or. node%daughter2%incoming) then
             node%n_t_channel = node%daughter2%n_t_channel + 1
          end if
       end if
    end if
  end subroutine node_count_specific_properties

  subroutine kingraph_assign_mappings_s (feyngraph, kingraph, feyngraph_set)
    type(feyngraph_t), intent(inout) :: feyngraph
    type(kingraph_t), pointer, intent(inout) :: kingraph
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    if (.not. (kingraph%root%daughter1%keep .and. kingraph%root%daughter2%keep)) then
       kingraph%keep = .false.
       call kingraph%tree%final ()
    end if
    if (kingraph%keep) then
       kingraph%root%on_shell = .true.
       kingraph%root%mapping = EXTERNAL_PRT
       kingraph%root%mapping_assigned = .true.
       call node_assign_bincode (kingraph%root)
       kingraph%root%ext_mass_sum = &
            kingraph%root%daughter1%ext_mass_sum + &
            kingraph%root%daughter2%ext_mass_sum
       if (kingraph%root%ext_mass_sum >= feyngraph_set%phs_par%sqrts) then
          kingraph%root%keep = .false.
          kingraph%keep = .false.; call kingraph%tree%final (); return
       end if
       call kingraph%root%subtree%add_entry (kingraph%root)
       kingraph%root%multiplicity &
            = kingraph%root%daughter1%multiplicity &
            + kingraph%root%daughter2%multiplicity
       kingraph%root%n_resonances &
            = kingraph%root%daughter1%n_resonances &
            + kingraph%root%daughter2%n_resonances
       kingraph%root%n_off_shell &
            = kingraph%root%daughter1%n_off_shell &
            + kingraph%root%daughter2%n_off_shell
       kingraph%root%n_log_enhanced &
            = kingraph%root%daughter1%n_log_enhanced &
            + kingraph%root%daughter2%n_log_enhanced
       if (kingraph%root%n_off_shell > feyngraph_set%phs_par%off_shell) then
          kingraph%root%keep = .false.
          kingraph%keep = .false.; call kingraph%tree%final (); return
       else
          kingraph%grove_prop%multiplicity = &
               kingraph%root%multiplicity
          kingraph%grove_prop%n_resonances = &
               kingraph%root%n_resonances
          kingraph%grove_prop%n_off_shell = &
               kingraph%root%n_off_shell
          kingraph%grove_prop%n_log_enhanced = &
               kingraph%root%n_log_enhanced
       end if
       kingraph%tree = kingraph%root%subtree
    end if
  end subroutine kingraph_assign_mappings_s

  subroutine kingraph_compute_mappings_t_line (feyngraph, kingraph, feyngraph_set)
    type(feyngraph_t), intent(inout) :: feyngraph
    type(kingraph_t), pointer, intent(inout) :: kingraph
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    call node_compute_t_line (feyngraph, kingraph, kingraph%root, feyngraph_set)
    if (.not. kingraph%root%keep) then
       kingraph%keep = .false.
       call kingraph%tree%final ()
    end if
    if (kingraph%keep) kingraph%tree = kingraph%root%subtree
  end subroutine kingraph_compute_mappings_t_line

  recursive subroutine node_compute_t_line (feyngraph, kingraph, node, feyngraph_set)
    type(feyngraph_t), intent(inout) :: feyngraph
    type(kingraph_t), intent(inout) :: kingraph
    type(k_node_t), intent(inout) :: node
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(k_node_t), pointer :: s_node
    type(k_node_t), pointer :: t_node
    type(k_node_t), pointer :: new_s_node
    if (.not. (node%daughter1%keep .and. node%daughter2%keep)) then
       node%keep = .false.
       return
    end if
    s_node => null ()
    t_node => null ()
    new_s_node => null ()
    if (associated (node%daughter1) .and. associated (node%daughter2)) then
       if (node%daughter1%t_line .or. node%daughter1%incoming) then
          t_node => node%daughter1; s_node => node%daughter2
       else if (node%daughter2%t_line .or. node%daughter2%incoming) then
          t_node => node%daughter2; s_node => node%daughter1
       end if
       if (t_node%t_line) then
          call node_compute_t_line (feyngraph, kingraph, t_node, feyngraph_set)
          if (.not. t_node%keep) then
             node%keep = .false.
             return
          end if
       else if (t_node%incoming) then
          t_node%mapping = EXTERNAL_PRT
          t_node%on_shell = .true.
          t_node%ext_mass_sum = t_node%particle%mass
          if (t_node%particle%mass >= feyngraph_set%phs_par%m_threshold_t) then
             t_node%effective_mass = t_node%particle%mass
          end if
          call t_node%subtree%add_entry (t_node)
       end if
!!! root:
       if (.not. node%incoming) then
          if (t_node%incoming) then
             node%ext_mass_sum = s_node%ext_mass_sum
          else
             node%ext_mass_sum &
                  = node%daughter1%ext_mass_sum &
                  + node%daughter2%ext_mass_sum
          end if
          if (node%particle%mass > feyngraph_set%phs_par%m_threshold_t) then
             node%effective_mass = max (node%particle%mass, &
                  s_node%effective_mass)
          else if (s_node%effective_mass > feyngraph_set%phs_par%m_threshold_t) then
             node%effective_mass = s_node%effective_mass
          else
             node%effective_mass = 0
          end if
!!! Allowed decay of beam particle
          if (t_node%incoming &
               .and. t_node%particle%mass > s_node%particle%mass &
               + node%particle%mass) then
             call beam_decay (feyngraph_set%fatal_beam_decay)
!!! Massless splitting
          else if (t_node%effective_mass == 0 &
               .and. s_node%effective_mass < feyngraph_set%phs_par%m_threshold_t &
               .and. node%effective_mass == 0) then
             node%mapping = U_CHANNEL
             node%log_enhanced = .true.
!!! IR radiation off massive particle
          else if (t_node%effective_mass /= 0 &
               .and. s_node%effective_mass == 0 &
               .and. node%effective_mass /= 0 &
               .and. (t_node%on_shell &
               .or. t_node%mapping == RADIATION) &
               .and. abs (t_node%effective_mass - node%effective_mass) &
               < feyngraph_set%phs_par%m_threshold_t) then
             node%log_enhanced = .true.
             node%mapping = RADIATION
          end if
          node%mapping_assigned = .true.
          call node_assign_bincode (node)
          call node%subtree%add_entry (node)
          call node_count_specific_properties (node)
          if (node%n_off_shell > feyngraph_set%phs_par%off_shell) then
             node%keep = .false.
             kingraph%keep = .false.; call kingraph%tree%final (); return
          else if (node%n_t_channel > feyngraph_set%phs_par%t_channel) then
             node%keep = .false.;
             kingraph%keep = .false.; call kingraph%tree%final (); return
          end if
       else
          node%mapping = EXTERNAL_PRT
          node%on_shell = .true.
          node%ext_mass_sum &
               = t_node%ext_mass_sum &
               + s_node%ext_mass_sum
          node%effective_mass = node%particle%mass
          if (.not. (node%ext_mass_sum < feyngraph_set%phs_par%sqrts)) then
             node%keep = .false.
             kingraph%keep = .false.; call kingraph%tree%final (); return
          end if
          if (kingraph%keep) then
             if (t_node%incoming .and. s_node%log_enhanced) then
                call s_node%f_node%k_node_list%add_entry (new_s_node, recycle=.false.)
                new_s_node = s_node
                new_s_node%daughter1 => s_node%daughter1
                new_s_node%daughter2 => s_node%daughter2
                if (s_node%index == node%daughter1%index) then
                   node%daughter1 => new_s_node
                else if (s_node%index ==  node%daughter2%index) then
                   node%daughter2 => new_s_node
                end if
                new_s_node%subtree = s_node%subtree
                new_s_node%mapping = NO_MAPPING
                new_s_node%log_enhanced = .false.
                new_s_node%n_log_enhanced &
                     = new_s_node%n_log_enhanced - 1
                new_s_node%log_enhanced = .false.
                where (new_s_node%subtree%bc == new_s_node%bincode)
                   new_s_node%subtree%mapping = NO_MAPPING
                endwhere
             else if ((t_node%t_line .or. t_node%incoming) .and. &
                  t_node%mapping == U_CHANNEL) then
                t_node%mapping = T_CHANNEL
                where (t_node%subtree%bc == t_node%bincode)
                   t_node%subtree%mapping = T_CHANNEL
                endwhere
             else if (t_node%incoming .and. &
                  .not. associated (s_node%daughter1) .and. &
                  .not. associated (s_node%daughter2)) then
                call s_node%f_node%k_node_list%add_entry (new_s_node, recycle=.false.)
                new_s_node = s_node
                new_s_node%mapping = ON_SHELL
                new_s_node%daughter1 => s_node%daughter1
                new_s_node%daughter2 => s_node%daughter2
                new_s_node%subtree = s_node%subtree
                if (s_node%index == node%daughter1%index) then
                   node%daughter1 => new_s_node
                else if (s_node%index == node%daughter2%index) then
                   node%daughter2 => new_s_node
                end if
                where (new_s_node%subtree%bc == new_s_node%bincode)
                   new_s_node%subtree%mapping = ON_SHELL
                endwhere
             end if
          end if
          call node%subtree%add_entry (node)
          node%multiplicity &
               = node%daughter1%multiplicity &
               + node%daughter2%multiplicity
          node%n_resonances &
               = node%daughter1%n_resonances &
               + node%daughter2%n_resonances
          node%n_off_shell &
               = node%daughter1%n_off_shell &
               + node%daughter2%n_off_shell
          node%n_log_enhanced &
               = node%daughter1%n_log_enhanced &
               + node%daughter2%n_log_enhanced
          node%n_t_channel &
               = node%daughter1%n_t_channel &
               + node%daughter2%n_t_channel
          if (node%n_off_shell > feyngraph_set%phs_par%off_shell) then
             node%keep = .false.
             kingraph%keep = .false.; call kingraph%tree%final (); return
          else if (node%n_t_channel > feyngraph_set%phs_par%t_channel) then
             node%keep = .false.
             kingraph%keep = .false.; call kingraph%tree%final (); return
          else
             kingraph%grove_prop%multiplicity = node%multiplicity
             kingraph%grove_prop%n_resonances = node%n_resonances
             kingraph%grove_prop%n_off_shell = node%n_off_shell
             kingraph%grove_prop%n_log_enhanced = node%n_log_enhanced
             kingraph%grove_prop%n_t_channel = node%n_t_channel
          end if
       end if
    end if
  contains
    subroutine beam_decay (fatal_beam_decay)
      logical, intent(in) :: fatal_beam_decay
      write (msg_buffer, "(1x,A,1x,'->',1x,A,1x,A)") &
           t_node%particle%particle_label, &
           node%particle%particle_label, &
           s_node%particle%particle_label
      call msg_message
      write (msg_buffer, "(1x,'mass(',A,') =',1x,E17.10)") &
           t_node%particle%particle_label, t_node%particle%mass
      call msg_message
      write (msg_buffer, "(1x,'mass(',A,') =',1x,E17.10)") &
           node%particle%particle_label, node%particle%mass
      call msg_message
      write (msg_buffer, "(1x,'mass(',A,') =',1x,E17.10)") &
           s_node%particle%particle_label, s_node%particle%mass
      call msg_message
      if (fatal_beam_decay) then
         call msg_fatal (" Phase space: Initial beam particle can decay")
      else
         call msg_warning (" Phase space: Initial beam particle can decay")
      end if
    end subroutine beam_decay
  end subroutine node_compute_t_line

  module subroutine feyngraph_make_inverse_kingraphs (feyngraph)
    class(feyngraph_t), intent(inout) :: feyngraph
    type(kingraph_t), pointer :: current
    current => feyngraph%kin_first
    do while (associated (current))
       if (current%inverse) exit
       call current%make_inverse_copy (feyngraph)
       current => current%next
    end do
  end subroutine feyngraph_make_inverse_kingraphs

  module subroutine feyngraph_compute_mappings (feyngraph, feyngraph_set)
    class(feyngraph_t), intent(inout) :: feyngraph
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(kingraph_t), pointer :: current
    current => feyngraph%kin_first
    do while (associated (current))
       if (feyngraph_set%process_type == DECAY) then
          call kingraph_assign_mappings_s (feyngraph, current, feyngraph_set)
       else if (feyngraph_set%process_type == SCATTERING) then
          call kingraph_compute_mappings_t_line &
               (feyngraph, current, feyngraph_set)
       end if
       current => current%next
    end do
  end subroutine feyngraph_compute_mappings

  subroutine f_node_list_compute_mappings_s (feyngraph_set)
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(f_node_ptr_t), dimension(:), allocatable :: set
    type(k_node_ptr_t), dimension(:), allocatable :: k_set
    type(k_node_entry_t), pointer :: k_entry
    type(f_node_entry_t), pointer :: current
    type(k_node_list_t), allocatable :: compare_list
    integer :: n_entries
    integer :: pos
    integer :: i, j, k
    do i = 1, feyngraph_set%f_node_list%max_tree_size - 2, 2
!!! Counter number of f_nodes with subtree size i for s channel calculations
       n_entries = 0
       if (feyngraph_set%use_dag) then
          do j=1, feyngraph_set%dag%n_nodes
             if (allocated (feyngraph_set%dag%node(j)%f_node)) then
                do k=1, size(feyngraph_set%dag%node(j)%f_node)
                   if (associated (feyngraph_set%dag%node(j)%f_node(k)%node)) then
                      if (.not. (feyngraph_set%dag%node(j)%f_node(k)%node%incoming &
                           .or. feyngraph_set%dag%node(j)%f_node(k)%node%t_line) &
                           .and. feyngraph_set%dag%node(j)%f_node(k)%node%n_subtree_nodes == i) then
                         n_entries = n_entries + 1
                      end if
                   end if
                end do
             end if
          end do
       else
          current => feyngraph_set%f_node_list%first
          do while (associated (current))
             if (.not. (current%node%incoming .or. current%node%t_line) &
                  .and. current%node%n_subtree_nodes == i) then
                n_entries = n_entries + 1
             end if
             current => current%next
          end do
       end if
       if (n_entries == 0) exit
!!! Create a temporary k node list for comparison
       allocate (set(n_entries))
       pos = 0
       if (feyngraph_set%use_dag) then
          do j=1, feyngraph_set%dag%n_nodes
             if (allocated (feyngraph_set%dag%node(j)%f_node)) then
                do k=1, size(feyngraph_set%dag%node(j)%f_node)
                   if (associated (feyngraph_set%dag%node(j)%f_node(k)%node)) then
                      if (.not. (feyngraph_set%dag%node(j)%f_node(k)%node%incoming &
                           .or. feyngraph_set%dag%node(j)%f_node(k)%node%t_line) &
                           .and. feyngraph_set%dag%node(j)%f_node(k)%node%n_subtree_nodes == i) then
                         pos = pos + 1
                         set(pos)%node => feyngraph_set%dag%node(j)%f_node(k)%node
                      end if
                   end if
                end do
             end if
          end do
       else
          current => feyngraph_set%f_node_list%first
          do while (associated (current))
             if (.not. (current%node%incoming .or. current%node%t_line) &
                  .and. current%node%n_subtree_nodes == i) then
                pos = pos + 1
                set(pos)%node => current%node
             end if
             current => current%next
          end do
       end if
       allocate (compare_list)
       compare_list%observer = .true.
       do j = 1, n_entries
          call k_node_init_from_f_node (set(j)%node, k_set, &
               feyngraph_set)
          if (allocated (k_set)) deallocate (k_set)
       end do
       !$OMP PARALLEL DO PRIVATE (k_entry)
       do j = 1, n_entries
          k_entry => set(j)%node%k_node_list%first
          do while (associated (k_entry))
             call node_assign_mapping_s(feyngraph_set%first, k_entry%node, feyngraph_set)
             k_entry => k_entry%next
          end do
       end do
       !$OMP END PARALLEL DO
       do j = 1, size (set)
          k_entry => set(j)%node%k_node_list%first
          do while (associated (k_entry))
             if (k_entry%node%keep) then
                if (k_entry%node%mapping == NO_MAPPING .or. k_entry%node%mapping == NONRESONANT) then
                   call compare_list%add_pointer (k_entry%node)
                end if
             end if
             k_entry => k_entry%next
          end do
       end do
       deallocate (set)
       call compare_list%check_subtree_equivalences(feyngraph_set%model)
       call compare_list%final
       deallocate (compare_list)
    end do
  end subroutine f_node_list_compute_mappings_s

  module subroutine grove_list_get_grove (grove_list, kingraph, &
       return_grove, preliminary)
    class(grove_list_t), intent(inout) :: grove_list
    type(kingraph_t), intent(in), pointer :: kingraph
    type(grove_t), intent(inout), pointer :: return_grove
    logical, intent(in) :: preliminary
    type(grove_t), pointer :: current_grove
    return_grove => null ()
    if (.not. associated(grove_list%first)) then
       allocate (grove_list%first)
       grove_list%first%grove_prop = kingraph%grove_prop
       return_grove => grove_list%first
       return
    end if
    current_grove => grove_list%first
    do while (associated (current_grove))
       if ((preliminary .and. &
            (current_grove%grove_prop .match. kingraph%grove_prop)) .or. &
            (.not. preliminary .and. &
            current_grove%grove_prop == kingraph%grove_prop)) then
          return_grove => current_grove
          exit
       else if (.not. associated (current_grove%next)) then
          allocate (current_grove%next)
          current_grove%next%grove_prop = kingraph%grove_prop
          if (size (kingraph%tree%bc) < 9) &
               current_grove%compare_tree%depth = 1
          return_grove => current_grove%next
          exit
       end if
       if (associated (current_grove%next)) then
          current_grove => current_grove%next
       end if
    end do
  end subroutine grove_list_get_grove

  module subroutine grove_list_add_kingraph (grove_list, kingraph, &
       preliminary, check, model)
    class(grove_list_t), intent(inout) :: grove_list
    type(kingraph_t), pointer, intent(inout) :: kingraph
    logical, intent(in) :: preliminary
    logical, intent(in) :: check
    type(model_data_t), optional, intent(in) :: model
    type(grove_t), pointer :: grove
    type(kingraph_t), pointer :: current
    integer, save :: index = 0
    grove => null ()
    current => null ()
    if (preliminary) then
       if (kingraph%index == 0) then
          index = index + 1
          kingraph%index = index
       end if
    end if
    call grove_list%get_grove (kingraph, grove, preliminary)
    if (check) then
       call grove%compare_tree%check_kingraph (kingraph, model, preliminary)
    end if
    if (kingraph%keep) then
       if (associated (grove%first)) then
          grove%last%grove_next => kingraph
          grove%last => kingraph
       else
          grove%first => kingraph
          grove%last => kingraph
       end if
    end if
  end subroutine grove_list_add_kingraph

  module subroutine grove_list_add_feyngraph (grove_list, feyngraph, model)
    class(grove_list_t), intent(inout) :: grove_list
    type(feyngraph_t), intent(inout) :: feyngraph
    type(model_data_t), intent(in) :: model
    type(kingraph_t), pointer :: current_kingraph, add_kingraph
    do while (associated (feyngraph%kin_first))
       if (feyngraph%kin_first%keep) then
          add_kingraph => feyngraph%kin_first
          feyngraph%kin_first => feyngraph%kin_first%next
          add_kingraph%next => null ()
          call grove_list%add_kingraph (kingraph=add_kingraph, &
               preliminary=.true., check=.true., model=model)
       else
          exit
       end if
    end do
    if (associated (feyngraph%kin_first)) then
       current_kingraph => feyngraph%kin_first
       do while (associated (current_kingraph%next))
          if (current_kingraph%next%keep) then
             add_kingraph => current_kingraph%next
             current_kingraph%next => current_kingraph%next%next
             add_kingraph%next => null ()
             call grove_list%add_kingraph (kingraph=add_kingraph, &
                  preliminary=.true., check=.true., model=model)
          else
             current_kingraph => current_kingraph%next
          end if
       end do
    end if
  end subroutine grove_list_add_feyngraph

  module function grove_prop_match (grove_prop1, grove_prop2) result (gp_match)
    type(grove_prop_t), intent(in) :: grove_prop1
    type(grove_prop_t), intent(in) :: grove_prop2
    logical :: gp_match
    gp_match = (grove_prop1%n_resonances == grove_prop2%n_resonances) &
         .and. (grove_prop1%n_log_enhanced == grove_prop2%n_log_enhanced) &
         .and. (grove_prop1%n_t_channel == grove_prop2%n_t_channel)
  end function grove_prop_match

  module function grove_prop_equal (grove_prop1, grove_prop2) result (gp_equal)
    type(grove_prop_t), intent(in) :: grove_prop1
    type(grove_prop_t), intent(in) :: grove_prop2
    logical :: gp_equal
    gp_equal = (grove_prop1%res_hash == grove_prop2%res_hash) &
         .and. (grove_prop1%n_resonances == grove_prop2%n_resonances) &
         .and. (grove_prop1%n_log_enhanced == grove_prop2%n_log_enhanced) &
         .and. (grove_prop1%n_off_shell == grove_prop2%n_off_shell) &
         .and. (grove_prop1%multiplicity == grove_prop2%multiplicity) &
         .and. (grove_prop1%n_t_channel == grove_prop2%n_t_channel)
  end function grove_prop_equal

  function kingraph_eqv (kingraph1, kingraph2) result (eqv)
    type(kingraph_t), intent(in) :: kingraph1
    type(kingraph_t), intent(inout) :: kingraph2
    logical :: eqv
    integer :: i
    logical :: equal
    eqv = .false.
    do i = kingraph1%tree%n_entries, 1, -1
       if (kingraph1%tree%bc(i) /= kingraph2%tree%bc(i)) return
    end do
    do i = kingraph1%tree%n_entries, 1, -1
       if ( .not. (kingraph1%tree%mapping(i) == kingraph2%tree%mapping(i) &
            .or. ((kingraph1%tree%mapping(i) == NO_MAPPING .or. &
            kingraph1%tree%mapping(i) == NONRESONANT) .and. &
            (kingraph2%tree%mapping(i) == NO_MAPPING .or. &
            kingraph2%tree%mapping(i) == NONRESONANT)))) return
    end do
    equal = .true.
    do i = kingraph1%tree%n_entries, 1, -1
       if (abs(kingraph1%tree%pdg(i)) /= abs(kingraph2%tree%pdg(i))) then
          equal = .false.;
          select case (kingraph1%tree%mapping(i))
          case (S_CHANNEL, RADIATION)
             select case (kingraph2%tree%mapping(i))
             case (S_CHANNEL, RADIATION)
                return
             end select
          end select
       end if
    end do
    if (equal) then
       kingraph2%keep = .false.
       call kingraph2%tree%final ()
    else
       eqv = .true.
    end if
  end function kingraph_eqv

  subroutine kingraph_select (kingraph1, kingraph2, model, preliminary)
    type(kingraph_t), intent(inout) :: kingraph1
    type(kingraph_t), intent(inout) :: kingraph2
    type(model_data_t), intent(in) :: model
    logical, intent(in) :: preliminary
    integer(TC), dimension(:), allocatable :: tmp_bc, daughter_bc
    integer, dimension(:), allocatable :: tmp_pdg, daughter_pdg
    integer, dimension (:), allocatable :: pdg_match
    integer :: i, j
    integer :: n_ext1, n_ext2
    if (kingraph_eqv (kingraph1, kingraph2)) then
       if (.not. preliminary) then
          kingraph2%keep = .false.; call kingraph2%tree%final ()
          return
       end if
       do i=1, size (kingraph1%tree%bc)
          if (abs(kingraph1%tree%pdg(i)) /= abs(kingraph2%tree%pdg(i))) then
             if (kingraph1%tree%mapping(i) /= EXTERNAL_PRT) then
                n_ext1 = popcnt (kingraph1%tree%bc(i))
                n_ext2 = n_ext1
                do j=i+1, size (kingraph1%tree%bc)
                   if (abs(kingraph1%tree%pdg(j)) /= abs(kingraph2%tree%pdg(j))) then
                      n_ext2 = popcnt (kingraph1%tree%bc(j))
                      if (n_ext2 < n_ext1) exit
                   end if
                end do
                if (n_ext2 < n_ext1) cycle
                allocate (tmp_bc(i-1))
                tmp_bc = kingraph1%tree%bc(:i-1)
                allocate (tmp_pdg(i-1))
                tmp_pdg = kingraph1%tree%pdg(:i-1)
                do j=i-1, 1, - 1
                   where (iand (tmp_bc(:j-1),tmp_bc(j)) /= 0 &
                        .or. iand(tmp_bc(:j-1),kingraph1%tree%bc(i)) == 0)
                      tmp_bc(:j-1) = 0
                      tmp_pdg(:j-1) = 0
                   endwhere
                end do
                allocate (daughter_bc(size(pack(tmp_bc, tmp_bc /= 0))))
                daughter_bc = pack (tmp_bc, tmp_bc /= 0)
                allocate (daughter_pdg(size(pack(tmp_pdg, tmp_pdg /= 0))))
                daughter_pdg = pack (tmp_pdg, tmp_pdg /= 0)
                if (size (daughter_pdg) == 2) then
                   call model%match_vertex(daughter_pdg(1), daughter_pdg(2), pdg_match)
                end if
                do j=1, size (pdg_match)
                   if (abs(pdg_match(j)) == abs(kingraph1%tree%pdg(i))) then
                      kingraph2%keep = .false.; call kingraph2%tree%final ()
                      exit
                   else if (abs(pdg_match(j)) == abs(kingraph2%tree%pdg(i))) then
                      kingraph1%keep = .false.; call kingraph1%tree%final ()
                      exit
                   end if
                end do
                deallocate (tmp_bc, tmp_pdg, daughter_bc, daughter_pdg, pdg_match)
                if (.not. (kingraph1%keep .and. kingraph2%keep)) exit
             end if
          end if
       end do
    end if
  end subroutine kingraph_select

  module subroutine grove_list_merge (target_list, grove_list, model, &
       prc_component)
    class(grove_list_t), intent(inout) :: target_list
    type(grove_list_t), intent(inout) :: grove_list
    type(model_data_t), intent(in) :: model
    integer, intent(in) :: prc_component
    type(grove_t), pointer :: current_grove
    type(kingraph_t), pointer :: current_graph
    current_grove => grove_list%first
    do while (associated (current_grove))
       do while (associated (current_grove%first))
          current_graph => current_grove%first
          current_grove%first => current_grove%first%grove_next
          current_graph%grove_next => null ()
          if (current_graph%keep) then
             current_graph%prc_component = prc_component
             call target_list%add_kingraph(kingraph=current_graph, &
                  preliminary=.false., check=.true., model=model)
          else
             call current_graph%final ()
             deallocate (current_graph)
          end if
       end do
       current_grove => current_grove%next
    end do
  end subroutine grove_list_merge

  module subroutine grove_list_rebuild (grove_list)
    class(grove_list_t), intent(inout) :: grove_list
    type(grove_list_t) :: tmp_list
    type(grove_t), pointer :: current_grove
    type(grove_t), pointer :: remove_grove
    type(kingraph_t), pointer :: current_graph
    type(kingraph_t), pointer :: next_graph
    tmp_list%first => grove_list%first
    grove_list%first => null ()
    current_grove => tmp_list%first
    do while (associated (current_grove))
       current_graph => current_grove%first
       do while (associated (current_graph))
          call current_graph%assign_resonance_hash ()
          next_graph => current_graph%grove_next
          current_graph%grove_next => null ()
          if (current_graph%keep) then
             call grove_list%add_kingraph (kingraph=current_graph, &
                  preliminary=.false., check=.false.)
          end if
          current_graph => next_graph
       end do
       current_grove => current_grove%next
    end do
    call tmp_list%final
  end subroutine grove_list_rebuild

  module subroutine feyngraph_set_write_file_format (feyngraph_set, u)
    type(feyngraph_set_t), intent(in) :: feyngraph_set
    integer, intent(in) :: u
    type(grove_t), pointer :: grove
    integer :: channel_number
    integer :: grove_number
    channel_number = 0
    grove_number = 0
    grove => feyngraph_set%grove_list%first
    do while (associated (grove))
       grove_number = grove_number + 1
       call grove%write_file_format &
            (feyngraph_set, grove_number, channel_number, u)
       grove => grove%next
    end do
  end subroutine feyngraph_set_write_file_format

  recursive module subroutine grove_write_file_format &
       (grove, feyngraph_set, gr_number, ch_number, u)
    class(grove_t), intent(in) :: grove
    type(feyngraph_set_t), intent(in) :: feyngraph_set
    integer, intent(in) :: u
    integer, intent(inout) :: gr_number
    integer, intent(inout) :: ch_number
    type(kingraph_t), pointer :: current
1   format(3x,A,1x,40(1x,I4))
    write (u, "(A)")
    write (u, "(1x,'!',1x,A,1x,I0,A)", advance='no') &
         'Multiplicity =', grove%grove_prop%multiplicity, ","
    select case (grove%grove_prop%n_resonances)
    case (0)
       write (u, '(1x,A)', advance='no') 'no resonances, '
    case (1)
       write (u, '(1x,A)', advance='no') '1 resonance,  '
    case default
       write (u, '(1x,I0,1x,A)', advance='no') &
            grove%grove_prop%n_resonances, 'resonances, '
    end select
    write (u, '(1x,I0,1x,A)', advance='no') &
         grove%grove_prop%n_log_enhanced, 'logs, '
    write (u, '(1x,I0,1x,A)', advance='no') &
         grove%grove_prop%n_off_shell, 'off-shell, '
    select case (grove%grove_prop%n_t_channel)
    case (0);  write (u, '(1x,A)') 's-channel graph'
    case (1);  write (u, '(1x,A)') '1 t-channel line'
    case default
       write(u,'(1x,I0,1x,A)') &
            grove%grove_prop%n_t_channel, 't-channel lines'
    end select
    write (u, '(1x,A,I0)') 'grove #', gr_number
    current => grove%first
    do while (associated (current))
       if (current%keep) then
          ch_number = ch_number + 1
          call current%write_file_format (feyngraph_set, ch_number, u)
       end if
       current => current%grove_next
    end do
  end subroutine grove_write_file_format

  module subroutine kingraph_write_file_format &
       (kingraph, feyngraph_set, ch_number, u)
    class(kingraph_t), intent(in) :: kingraph
    type(feyngraph_set_t), intent(in) :: feyngraph_set
    integer, intent(in) :: ch_number
    integer, intent(in) :: u
    integer :: i
    integer(TC) :: bincode_incoming
2   format(3X,'map',1X,I3,1X,A,1X,I9,1X,'!',1X,A)
    !!! determine bincode of incoming particle from tree
    bincode_incoming = maxval (kingraph%tree%bc)
    write (unit=u, fmt='(1X,A,I0)') '! Channel #', ch_number
    write (unit=u, fmt='(3X,A,1X)', advance='no') 'tree'
    do i=1, size (kingraph%tree%bc)
       if (kingraph%tree%mapping(i) >=0 &
            .or. kingraph%tree%mapping(i) == NONRESONANT &
            .or. (kingraph%tree%bc(i) == bincode_incoming &
            .and. feyngraph_set%process_type == DECAY)) then
          write (unit=u, fmt='(1X,I0)', advance='no') kingraph%tree%bc(i)
       end if
    end do
    write (unit=u, fmt='(A)', advance='yes')
    do i=1, size(kingraph%tree%bc)
       select case (kingraph%tree%mapping(i))
       case (NO_MAPPING, NONRESONANT, EXTERNAL_PRT)
       case (S_CHANNEL)
          write (unit=u, fmt=2) kingraph%tree%bc(i), 's_channel', &
               kingraph%tree%pdg(i), &
               trim(get_particle_name (feyngraph_set, kingraph%tree%pdg(i)))
       case (T_CHANNEL)
          write (unit=u, fmt=2) kingraph%tree%bc(i), 't_channel', &
               abs (kingraph%tree%pdg(i)), &
               trim(get_particle_name (feyngraph_set, abs(kingraph%tree%pdg(i))))
       case (U_CHANNEL)
          write (unit=u, fmt=2) kingraph%tree%bc(i), 'u_channel', &
               abs (kingraph%tree%pdg(i)), &
               trim(get_particle_name (feyngraph_set, abs(kingraph%tree%pdg(i))))
       case (RADIATION)
          write (unit=u, fmt=2) kingraph%tree%bc(i), 'radiation', &
               kingraph%tree%pdg(i), &
               trim(get_particle_name (feyngraph_set, kingraph%tree%pdg(i)))
       case (COLLINEAR)
          write (unit=u, fmt=2) kingraph%tree%bc(i), 'collinear', &
               kingraph%tree%pdg(i), &
               trim(get_particle_name (feyngraph_set, kingraph%tree%pdg(i)))
       case (INFRARED)
          write (unit=u, fmt=2) kingraph%tree%bc(i), 'infrared ', &
               kingraph%tree%pdg(i), &
               trim(get_particle_name (feyngraph_set, kingraph%tree%pdg(i)))
       case (ON_SHELL)
          write (unit=u, fmt=2) kingraph%tree%bc(i), 'on_shell ', &
               kingraph%tree%pdg(i), &
               trim(get_particle_name (feyngraph_set, kingraph%tree%pdg(i)))
       case default
          call msg_bug (" Impossible mapping mode encountered")
       end select
    end do
  end subroutine kingraph_write_file_format

   function get_particle_name (feyngraph_set, pdg) result (particle_name)
     type(feyngraph_set_t), intent(in) :: feyngraph_set
     integer, intent(in) :: pdg
     character(len=LABEL_LEN) :: particle_name
     integer :: i
     do i=1, size (feyngraph_set%particle)
        if (feyngraph_set%particle(i)%pdg == pdg) then
           particle_name = feyngraph_set%particle(i)%particle_label
           exit
        end if
     end do
   end function get_particle_name

  module subroutine feyngraph_make_invertible (feyngraph)
    class(feyngraph_t), intent(inout) :: feyngraph
    logical :: t_line_found
    feyngraph%root%incoming = .true.
    t_line_found = .false.
    if (associated (feyngraph%root%daughter1)) then
       call f_node_t_line_check (feyngraph%root%daughter1, t_line_found)
       if (.not. t_line_found) then
          if (associated (feyngraph%root%daughter2)) then
             call f_node_t_line_check (feyngraph%root%daughter2, t_line_found)
          end if
       end if
    end if

  contains

  recursive subroutine f_node_t_line_check (node, t_line_found)
    type(f_node_t), target, intent(inout) :: node
    integer :: pos
    logical, intent(inout) :: t_line_found
    if (associated (node%daughter1)) then
       call f_node_t_line_check (node%daughter1, t_line_found)
       if (node%daughter1%incoming .or. node%daughter1%t_line) then
          node%t_line = .true.
       else if (associated (node%daughter2)) then
          call f_node_t_line_check (node%daughter2, t_line_found)
          if (node%daughter2%incoming .or. node%daughter2%t_line) then
             node%t_line = .true.
          end if
       end if
    else
       pos = index (node%particle_label, '[') + 1
       if (node%particle_label(pos:pos) == '2') then
          node%incoming = .true.
          t_line_found = .true.
       end if
    end if
  end subroutine f_node_t_line_check

  end subroutine feyngraph_make_invertible

  module subroutine kingraph_make_inverse_copy (original_kingraph, feyngraph)
    class(kingraph_t), intent(inout) :: original_kingraph
    type(feyngraph_t), intent(inout) :: feyngraph
    type(kingraph_t), pointer :: kingraph_copy
    type(k_node_t), pointer :: potential_root
    allocate (kingraph_copy)
    if (associated (feyngraph%kin_last)) then
       allocate (feyngraph%kin_last%next)
       feyngraph%kin_last => feyngraph%kin_last%next
    else
       allocate(feyngraph%kin_first)
       feyngraph%kin_last => feyngraph%kin_first
    end if
    kingraph_copy => feyngraph%kin_last
    call kingraph_set_inverse_daughters (original_kingraph)
    kingraph_copy%inverse = .true.
    kingraph_copy%n_nodes = original_kingraph%n_nodes
    kingraph_copy%keep = original_kingraph%keep
    potential_root => original_kingraph%root
    do while (.not. potential_root%incoming .or. &
         (associated (potential_root%daughter1) .and. &
          associated (potential_root%daughter2)))
       if (potential_root%daughter1%incoming .or. &
           potential_root%daughter1%t_line) then
          potential_root => potential_root%daughter1
       else if (potential_root%daughter2%incoming .or. &
            potential_root%daughter2%t_line) then
          potential_root => potential_root%daughter2
       end if
    end do
    call node_inverse_deep_copy (potential_root, kingraph_copy%root)
  end subroutine kingraph_make_inverse_copy

  recursive subroutine node_inverse_deep_copy (original_node, node_copy)
    type(k_node_t), intent(in) :: original_node
    type(k_node_t), pointer, intent(out) :: node_copy
    call original_node%f_node%k_node_list%add_entry(node_copy, recycle=.false.)
    node_copy = original_node
    if (node_copy%t_line .or. node_copy%incoming) then
       node_copy%particle => original_node%particle%anti
    else
       node_copy%particle => original_node%particle
    end if
    if (associated (original_node%inverse_daughter1) .and. associated (original_node%inverse_daughter2)) then
       if (original_node%inverse_daughter1%incoming .or. original_node%inverse_daughter1%t_line) then
          node_copy%daughter2 => original_node%inverse_daughter2
          call node_inverse_deep_copy (original_node%inverse_daughter1, &
               node_copy%daughter1)
       else if (original_node%inverse_daughter2%incoming .or. original_node%inverse_daughter2%t_line) then
          node_copy%daughter1 => original_node%inverse_daughter1
          call node_inverse_deep_copy (original_node%inverse_daughter2, &
               node_copy%daughter2)
       end if
    end if
  end subroutine node_inverse_deep_copy

  module subroutine feyngraph_set_generate_single (feyngraph_set, model, &
       n_in, n_out, phs_par, fatal_beam_decay, u_in)
    type(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(model_data_t), target, intent(in) :: model
    integer, intent(in) :: n_in, n_out
    type(phs_parameters_t), intent(in) :: phs_par
    logical, intent(in) :: fatal_beam_decay
    integer, intent(in) :: u_in
    feyngraph_set%n_in = n_in
    feyngraph_set%n_out = n_out
    feyngraph_set%process_type = n_in
    feyngraph_set%phs_par = phs_par
    feyngraph_set%model => model
    if (debug_on)  call msg_debug &
         (D_PHASESPACE, "Construct relevant Feynman diagrams from Omega output")
    call feyngraph_set%build (u_in)
    if (debug_on)  call msg_debug &
         (D_PHASESPACE, "Find phase-space parametrizations")
    call feyngraph_set_find_phs_parametrizations(feyngraph_set)
  end subroutine feyngraph_set_generate_single

  subroutine feyngraph_set_find_phs_parametrizations (feyngraph_set)
    class(feyngraph_set_t), intent(inout) :: feyngraph_set
    type(feyngraph_t), pointer :: current => null ()
    type(feyngraph_ptr_t), dimension (:), allocatable :: set
    integer :: pos
    integer :: i
    allocate (set (feyngraph_set%n_graphs))
    pos = 0
    current => feyngraph_set%first
    do while (associated (current))
       pos = pos + 1
       set(pos)%graph => current
       current => current%next
    end do
    if (feyngraph_set%process_type == SCATTERING) then
       !$OMP PARALLEL DO
       do i=1, feyngraph_set%n_graphs
          if (set(i)%graph%keep) then
             call set(i)%graph%make_invertible ()
          end if
       end do
       !$OMP END PARALLEL DO
    end if
    call f_node_list_compute_mappings_s (feyngraph_set)
    do i=1, feyngraph_set%n_graphs
       if (set(i)%graph%keep) then
          call set(i)%graph%make_kingraphs (feyngraph_set)
       end if
    end do
    if (feyngraph_set%process_type == SCATTERING) then
       do i=1, feyngraph_set%n_graphs
          if (set(i)%graph%keep) then
             call set(i)%graph%make_inverse_kingraphs ()
          end if
       end do
    end if
    do i=1, feyngraph_set%n_graphs
       if (set(i)%graph%keep) then
          call set(i)%graph%compute_mappings (feyngraph_set)
       end if
    end do
    do i=1, feyngraph_set%n_graphs
       if (set(i)%graph%keep) then
          call feyngraph_set%grove_list%add_feyngraph (set(i)%graph, &
               feyngraph_set%model)
       end if
    end do
  end subroutine feyngraph_set_find_phs_parametrizations

  elemental module function tree_equal (tree1, tree2) result (flag)
    type(tree_t), intent(in) :: tree1, tree2
    logical :: flag
    if (tree1%n_entries == tree2%n_entries) then
       if (tree1%bc(size(tree1%bc)) == tree2%bc(size(tree2%bc))) then
          flag = all (tree1%mapping == tree2%mapping) .and. &
               all (tree1%bc == tree2%bc) .and. &
               all (abs(tree1%pdg) == abs(tree2%pdg))
       else
          flag = .false.
       end if
    else
       flag = .false.
    end if
  end function tree_equal

  pure module function subtree_eqv (subtree1, subtree2) result (eqv)
    type(tree_t), intent(in) :: subtree1, subtree2
    logical :: eqv
    integer :: root_pos
    integer :: i
    logical :: equal
    eqv = .false.
    if (subtree1%n_entries /= subtree2%n_entries) return
    root_pos = subtree1%n_entries
    if (subtree1%mapping(root_pos) == NONRESONANT .or. &
         subtree2%mapping(root_pos) == NONRESONANT .or. &
         (subtree1%mapping(root_pos) == NO_MAPPING .and. &
         subtree2%mapping(root_pos) == NO_MAPPING .and. &
         abs(subtree1%pdg(root_pos)) == abs(subtree2%pdg(root_pos)))) then
       do i = subtree1%n_entries, 1, -1
          if (subtree1%bc(i) /= subtree2%bc(i)) return
       end do
       equal = .true.
       do i = subtree1%n_entries, 1, -1
          if (abs(subtree1%pdg(i)) /= abs (subtree2%pdg(i))) then
             select case (subtree1%mapping(i))
             case (NO_MAPPING, NONRESONANT)
                select case (subtree2%mapping(i))
                case (NO_MAPPING, NONRESONANT)
                   equal = .false.
                case default
                   return
                end select
             case default
                return
             end select
          end if
       end do
       do i = subtree1%n_entries, 1, -1
          if (subtree1%mapping(i) /= subtree2%mapping(i)) then
             select case (subtree1%mapping(i))
             case (NO_MAPPING, NONRESONANT)
                select case (subtree2%mapping(i))
                case (NO_MAPPING, NONRESONANT)
                case default
                   return
                end select
             case default
                return
             end select
          end if
       end do
       if (.not. equal) eqv = .true.
    end if
  end function subtree_eqv

  subroutine subtree_select (subtree1, subtree2, model)
    type(tree_t), intent(inout) :: subtree1, subtree2
    type(model_data_t), intent(in) :: model
    integer :: j, k
    integer(TC), dimension(:), allocatable :: tmp_bc, daughter_bc
    integer, dimension(:), allocatable :: tmp_pdg, daughter_pdg
    integer, dimension (:), allocatable :: pdg_match
    if (subtree1 .eqv. subtree2) then
       do j=1, subtree1%n_entries
          if (abs(subtree1%pdg(j)) /= abs(subtree2%pdg(j))) then
             tmp_bc = subtree1%bc(:j-1); tmp_pdg = subtree1%pdg(:j-1)
             do k=j-1, 1, - 1
                where (iand (tmp_bc(:k-1),tmp_bc(k)) /= 0 &
                     .or. iand(tmp_bc(:k-1),subtree1%bc(j)) == 0)
                   tmp_bc(:k-1) = 0
                   tmp_pdg(:k-1) = 0
                endwhere
             end do
             daughter_bc = pack (tmp_bc, tmp_bc /= 0)
             daughter_pdg = pack (tmp_pdg, tmp_pdg /= 0)
             if (size (daughter_pdg) == 2) then
                call model%match_vertex(daughter_pdg(1), daughter_pdg(2), pdg_match)
                if (.not. allocated (pdg_match)) then
!!! Relevant if tree contains only abs (pdg). In this case, changing the
!!! sign of one of the pdg codes should give a result.
                   call model%match_vertex(-daughter_pdg(1), daughter_pdg(2), pdg_match)
                end if
             end if
             do k=1, size (pdg_match)
                if (abs(pdg_match(k)) == abs(subtree1%pdg(j))) then
                   if (subtree1%keep) subtree2%keep = .false.
                   exit
                else if (abs(pdg_match(k)) == abs(subtree2%pdg(j))) then
                   if (subtree2%keep) subtree1%keep = .false.
                   exit
                end if
             end do
             deallocate (tmp_bc, tmp_pdg, daughter_bc, daughter_pdg, pdg_match)
             if (.not. (subtree1%keep .and. subtree2%keep)) exit
          end if
       end do
    end if
  end subroutine subtree_select

  module subroutine kingraph_assign_resonance_hash (kingraph)
    class(kingraph_t), intent(inout) :: kingraph
    logical, dimension (:), allocatable :: tree_resonant
    integer(i8), dimension(1) :: mold
    allocate (tree_resonant (kingraph%tree%n_entries))
    tree_resonant = (kingraph%tree%mapping == S_CHANNEL)
    kingraph%grove_prop%res_hash = hash (transfer &
         ([sort (pack (kingraph%tree%pdg, tree_resonant)), &
           sort (pack (abs (kingraph%tree%pdg), &
           kingraph%tree%mapping == T_CHANNEL .or. &
           kingraph%tree%mapping == U_CHANNEL))], mold))
    deallocate (tree_resonant)
  end subroutine kingraph_assign_resonance_hash

  module subroutine feyngraph_set_write_process_bincode_format &
       (feyngraph_set, unit)
    type(feyngraph_set_t), intent(in), target :: feyngraph_set
    integer, intent(in), optional :: unit
    integer, dimension(:), allocatable :: bincode, field_width
    integer :: n_in, n_out, n_tot, n_flv
    integer :: u, f, i, bc
    character(20) :: str
    type(string_t) :: fmt_head
    type(string_t), dimension(:), allocatable :: fmt_proc
    u = given_output_unit (unit);  if (u < 0)  return
    if (.not. allocated (feyngraph_set%flv)) return
    write (u, "('!',1x,A)")  "List of subprocesses with particle bincodes:"
    n_in  = feyngraph_set%n_in
    n_out = feyngraph_set%n_out
    n_tot = n_in + n_out
    n_flv = size (feyngraph_set%flv, 2)
    allocate (bincode (n_tot), field_width (n_tot), fmt_proc (n_tot))
    bc = 1
    do i = 1, n_out
       bincode(n_in + i) = bc
       bc = 2 * bc
    end do
    do i = n_in, 1, -1
       bincode(i) = bc
       bc = 2 * bc
    end do
    do i = 1, n_tot
       write (str, "(I0)")  bincode(i)
       field_width(i) = len_trim (str)
       do f = 1, n_flv
          field_width(i) = max (field_width(i), &
               len (feyngraph_set%flv(i,f)%get_name ()))
       end do
    end do
    fmt_head = "('!'"
    do i = 1, n_tot
       fmt_head = fmt_head // ",1x,"
       fmt_proc(i) = "(1x,"
       write (str, "(I0)")  field_width(i)
       fmt_head = fmt_head // "I" // trim(str)
       fmt_proc(i) = fmt_proc(i) // "A" // trim(str)
       if (i == n_in) then
          fmt_head = fmt_head // ",1x,'  '"
       end if
    end do
    do i = 1, n_tot
       fmt_proc(i) = fmt_proc(i) // ")"
    end do
    fmt_head = fmt_head // ")"
    write (u, char (fmt_head))  bincode
    do f = 1, n_flv
       write (u, "('!')", advance="no")
       do i = 1, n_tot
          write (u, char (fmt_proc(i)), advance="no") &
               char (feyngraph_set%flv(i,f)%get_name ())
          if (i == n_in)  write (u, "(1x,'=>')", advance="no")
       end do
       write (u, *)
    end do
    write (u, char (fmt_head))  bincode
  end subroutine feyngraph_set_write_process_bincode_format

  module subroutine feyngraph_set_write_graph_format &
       (feyngraph_set, filename, process_id, unit)
    type(feyngraph_set_t), intent(in), target :: feyngraph_set
    type(string_t), intent(in) :: filename, process_id
    integer, intent(in), optional :: unit
    type(kingraph_t), pointer :: kingraph
    type(grove_t), pointer :: grove
    integer :: u, n_grove, count, pgcount
    logical :: first_in_grove
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, '(A)') "\documentclass[10pt]{article}"
    write (u, '(A)') "\usepackage{amsmath}"
    write (u, '(A)') "\usepackage{feynmp}"
    write (u, '(A)') "\usepackage{url}"
    write (u, '(A)') "\usepackage{color}"
    write (u, *)
    write (u, '(A)') "\textwidth 18.5cm"
    write (u, '(A)') "\evensidemargin -1.5cm"
    write (u, '(A)') "\oddsidemargin -1.5cm"
    write (u, *)
    write (u, '(A)') "\newcommand{\blue}{\color{blue}}"
    write (u, '(A)') "\newcommand{\green}{\color{green}}"
    write (u, '(A)') "\newcommand{\red}{\color{red}}"
    write (u, '(A)') "\newcommand{\magenta}{\color{magenta}}"
    write (u, '(A)') "\newcommand{\cyan}{\color{cyan}}"
    write (u, '(A)') "\newcommand{\sm}{\footnotesize}"
    write (u, '(A)') "\setlength{\parindent}{0pt}"
    write (u, '(A)') "\setlength{\parsep}{20pt}"
    write (u, *)
    write (u, '(A)') "\begin{document}"
    write (u, '(A)') "\begin{fmffile}{" // char (filename) // "}"
    write (u, '(A)') "\fmfcmd{color magenta; magenta = red + blue;}"
    write (u, '(A)') "\fmfcmd{color cyan; cyan = green + blue;}"
    write (u, '(A)') "\begin{fmfshrink}{0.5}"
    write (u, '(A)') "\begin{flushleft}"
    write (u, *)
    write (u, '(A)') "\noindent" // &
         & "\textbf{\large\texttt{WHIZARD} phase space channels}" // &
         & "\hfill\today"
    write (u, *)
    write (u, '(A)') "\vspace{10pt}"
    write (u, '(A)') "\noindent" // &
         & "\textbf{Process:} \url{" // char (process_id) // "}"
    call feyngraph_set_write_process_tex_format (feyngraph_set, u)
    write (u, *)
    write (u, '(A)') "\noindent" // &
         & "\textbf{Note:} These are pseudo Feynman graphs that "
    write (u, '(A)') "visualize phase-space parameterizations " // &
         & "(``integration channels'').  "
    write (u, '(A)') "They do \emph{not} indicate Feynman graphs used for the " // &
         & "matrix element."
    write (u, *)
    write (u, '(A)') "\textbf{Color code:} " // &
         & "{\blue resonance,} " // &
         & "{\cyan t-channel,} " // &
         & "{\green radiation,} "
    write (u, '(A)') "{\red infrared,} " // &
         & "{\magenta collinear,} " // &
         & "external/off-shell"
    write (u, *)
    write (u, '(A)') "\noindent" // &
         & "\textbf{Black square:} Keystone, indicates ordering of " // &
         & "phase space parameters."
    write (u, *)
    write (u, '(A)') "\vspace{-20pt}"
    count = 0
    pgcount = 0
    n_grove = 0
    grove => feyngraph_set%grove_list%first
    do while (associated (grove))
       n_grove = n_grove + 1
       write (u, *)
       write (u, '(A)') "\vspace{20pt}"
       write (u, '(A)') "\begin{tabular}{l}"
       write (u, '(A,I5,A)') &
            & "\fbox{\bf Grove \boldmath$", n_grove, "$} \\[10pt]"
       write (u, '(A,I1,A)') "Multiplicity: ", &
            grove%grove_prop%multiplicity, "\\"
       write (u, '(A,I1,A)') "Resonances:   ", &
            grove%grove_prop%n_resonances, "\\"
       write (u, '(A,I1,A)') "Log-enhanced: ", &
            grove%grove_prop%n_log_enhanced, "\\"
       write (u, '(A,I1,A)') "Off-shell:    ", &
            grove%grove_prop%n_off_shell, "\\"
       write (u, '(A,I1,A)') "t-channel:    ", &
            grove%grove_prop%n_t_channel, ""
       write (u, '(A)') "\end{tabular}"
       kingraph => grove%first
       do while (associated (kingraph))
          count = count + 1
          call kingraph_write_graph_format (kingraph, count, unit)
          kingraph => kingraph%grove_next
       end do
       grove => grove%next
    end do
    write (u, '(A)') "\end{flushleft}"
    write (u, '(A)') "\end{fmfshrink}"
    write (u, '(A)') "\end{fmffile}"
    write (u, '(A)') "\end{document}"
  end subroutine feyngraph_set_write_graph_format

  subroutine feyngraph_set_write_process_tex_format (feyngraph_set, unit)
    type(feyngraph_set_t), intent(in), target :: feyngraph_set
    integer, intent(in), optional :: unit
    integer :: n_tot
    integer :: u, f, i
    n_tot = feyngraph_set%n_in + feyngraph_set%n_out
    u = given_output_unit (unit);  if (u < 0)  return
    if (.not. allocated (feyngraph_set%flv)) return
    write (u, "(A)")  "\begin{align*}"
    do f = 1, size (feyngraph_set%flv, 2)
       do i = 1, feyngraph_set%n_in
          if (i > 1)  write (u, "(A)", advance="no") "\quad "
          write (u, "(A)", advance="no") &
               char (feyngraph_set%flv(i,f)%get_tex_name ())
       end do
       write (u, "(A)", advance="no")  "\quad &\to\quad "
       do i = feyngraph_set%n_in + 1, n_tot
          if (i > feyngraph_set%n_in + 1)  write (u, "(A)", advance="no") "\quad "
          write (u, "(A)", advance="no") &
               char (feyngraph_set%flv(i,f)%get_tex_name ())
       end do
       if (f < size (feyngraph_set%flv, 2)) then
          write (u, "(A)")  "\\"
       else
          write (u, "(A)")  ""
       end if
    end do
    write (u, "(A)")  "\end{align*}"
  end subroutine feyngraph_set_write_process_tex_format

  subroutine kingraph_write_graph_format (kingraph, count, unit)
    type(kingraph_t), intent(in) :: kingraph
    integer, intent(in) :: count
    integer, intent(in), optional :: unit
    integer :: u
    type(string_t) :: left_str, right_str
    u = given_output_unit (unit);  if (u < 0)  return
    left_str = ""
    right_str = ""
    write (u, '(A)') "\begin{minipage}{105pt}"
    write (u, '(A)') "\vspace{30pt}"
    write (u, '(A)') "\begin{center}"
    write (u, '(A)') "\begin{fmfgraph*}(55,55)"
    call graph_write_node (kingraph%root)
    write (u, '(A)') "\fmfleft{" // char (extract (left_str, 2)) // "}"
    write (u, '(A)') "\fmfright{" // char (extract (right_str, 2)) // "}"
    write (u, '(A)') "\end{fmfgraph*}\\"
    write (u, '(A,I5,A)') "\fbox{$", count, "$}"
    write (u, '(A)') "\end{center}"
    write (u, '(A)') "\end{minipage}"
    write (u, '(A)') "%"
  contains
    recursive subroutine graph_write_node (node)
      type(k_node_t), intent(in) :: node
      if (associated (node%daughter1) .or. associated (node%daughter2)) then
         if (node%daughter2%t_line .or. node%daughter2%incoming) then
            call vertex_write (node, node%daughter2)
            call vertex_write (node, node%daughter1)
         else
            call vertex_write (node, node%daughter1)
            call vertex_write (node, node%daughter2)
         end if
         if (node%mapping == EXTERNAL_PRT) then
            call line_write (node%bincode, 0, node%particle)
            call external_write (node%bincode, node%particle%tex_name, &
                 left_str)
            write (u, '(A,I0,A)') "\fmfv{d.shape=square}{v0}"
         end if
      else
         if (node%incoming) then
            call external_write (node%bincode, node%particle%anti%tex_name, &
                 left_str)
         else
            call external_write (node%bincode, node%particle%tex_name, &
                 right_str)
         end if
      end if
    end subroutine graph_write_node
    recursive subroutine vertex_write (node, daughter)
      type(k_node_t), intent(in) :: node, daughter
      integer :: bincode
      if (associated (node%daughter1) .and. associated (node%daughter2) &
           .and. node%mapping == EXTERNAL_PRT) then
         bincode = 0
      else
         bincode = node%bincode
      end if
      call graph_write_node (daughter)
      if (associated (node%daughter1) .or. associated (node%daughter2)) then
         call line_write (bincode, daughter%bincode, daughter%particle, &
              mapping=daughter%mapping)
      else
         call line_write (bincode, daughter%bincode, daughter%particle)
      end if
    end subroutine vertex_write
    subroutine line_write (i1, i2, particle, mapping)
      integer(TC), intent(in) :: i1, i2
      type(part_prop_t), intent(in) :: particle
      integer, intent(in), optional :: mapping
      integer :: k1, k2
      type(string_t) :: prt_type
      select case (particle%spin_type)
      case (SCALAR);       prt_type = "plain"
      case (SPINOR);       prt_type = "fermion"
      case (VECTOR);       prt_type = "boson"
      case (VECTORSPINOR); prt_type = "fermion"
      case (TENSOR);       prt_type = "dbl_wiggly"
      case default;        prt_type = "dashes"
      end select
      if (particle%pdg < 0) then
!!! anti-particle
         k1 = i2;  k2 = i1
      else
         k1 = i1;  k2 = i2
      end if
      if (present (mapping)) then
         select case (mapping)
         case (S_CHANNEL)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=blue,lab=\sm\blue$" // &
                 & char (particle%tex_name) // "$}" // &
                 & "{v", k1, ",v", k2, "}"
         case (T_CHANNEL, U_CHANNEL)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=cyan,lab=\sm\cyan$" // &
                 & char (particle%tex_name) // "$}" // &
                 & "{v", k1, ",v", k2, "}"
         case (RADIATION)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=green,lab=\sm\green$" // &
                 & char (particle%tex_name) // "$}" // &
                 & "{v", k1, ",v", k2, "}"
         case (COLLINEAR)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=magenta,lab=\sm\magenta$" // &
                 & char (particle%tex_name) // "$}" // &
                 & "{v", k1, ",v", k2, "}"
         case (INFRARED)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=red,lab=\sm\red$" // &
                 & char (particle%tex_name) // "$}" // &
                 & "{v", k1, ",v", k2, "}"
         case default
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=black}" // &
                 & "{v", k1, ",v", k2, "}"
         end select
      else
         write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & "}" // &
                 & "{v", k1, ",v", k2, "}"
      end if
    end subroutine line_write
    subroutine external_write (bincode, name, ext_str)
      integer(TC), intent(in) :: bincode
      type(string_t), intent(in) :: name
      type(string_t), intent(inout) :: ext_str
      character(len=20) :: str
      write (str, '(A2,I0)') ",v", bincode
      ext_str = ext_str // trim (str)
      write (u, '(A,I0,A,I0,A)') "\fmflabel{\sm$" &
        // char (name) &
        // "\,(", bincode, ")" &
        // "$}{v", bincode, "}"
    end subroutine external_write
  end subroutine kingraph_write_graph_format

  module subroutine feyngraph_set_generate &
    (feyngraph_set, model, n_in, n_out, flv, phs_par, fatal_beam_decay, &
    u_in, vis_channels, use_dag)
    type(feyngraph_set_t), intent(out) :: feyngraph_set
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: n_in, n_out
    type(flavor_t), dimension(:,:), intent(in) :: flv
    type(phs_parameters_t), intent(in) :: phs_par
    logical, intent(in) :: fatal_beam_decay
    integer, intent(in) :: u_in
    logical, intent(in) :: vis_channels
    logical, optional, intent(in) :: use_dag
    type(grove_t), pointer :: grove
    integer :: i, j
    type(kingraph_t), pointer :: kingraph
    if (phase_space_vanishes (phs_par%sqrts, n_in, flv))  return
    if (present (use_dag)) feyngraph_set%use_dag = use_dag
    feyngraph_set%process_type = n_in
    feyngraph_set%n_in = n_in
    feyngraph_set%n_out = n_out
    allocate (feyngraph_set%flv (size (flv, 1), size (flv, 2)))
    do i = 1, size (flv, 2)
       do j = 1, size (flv, 1)
          call feyngraph_set%flv(j,i)%init (flv(j,i)%get_pdg (), model)
       end do
    end do
    allocate (feyngraph_set%particle (PRT_ARRAY_SIZE))
    allocate (feyngraph_set%grove_list)
    allocate (feyngraph_set%fset (size (flv, 2)))
    do i = 1, size (feyngraph_set%fset)
       feyngraph_set%fset(i)%use_dag = feyngraph_set%use_dag
       allocate (feyngraph_set%fset(i)%flv(size (flv,1),1))
       feyngraph_set%fset(i)%flv(:,1) = flv(:,i)
       feyngraph_set%fset(i)%particle => feyngraph_set%particle
       allocate (feyngraph_set%fset(i)%grove_list)
       call feyngraph_set_generate_single (feyngraph_set%fset(i), &
            model, n_in, n_out, phs_par, fatal_beam_decay, u_in)
       call feyngraph_set%grove_list%merge &
            (feyngraph_set%fset(i)%grove_list, model, i)
       if (.not. vis_channels) call feyngraph_set%fset(i)%final()
    end do
    call feyngraph_set%grove_list%rebuild ()
  end subroutine feyngraph_set_generate

  module function feyngraph_set_is_valid (feyngraph_set) result (flag)
    class(feyngraph_set_t), intent(in) :: feyngraph_set
    type(kingraph_t), pointer :: kingraph
    type(grove_t), pointer :: grove
    logical :: flag
    flag = .false.
    if (associated (feyngraph_set%grove_list)) then
       grove => feyngraph_set%grove_list%first
       do while (associated (grove))
          kingraph => grove%first
          do while (associated (kingraph))
             if (kingraph%keep) then
                flag = .true.
                return
             end if
             kingraph => kingraph%next
          end do
          grove => grove%next
       end do
    end if
  end function feyngraph_set_is_valid

  module subroutine kingraph_extract_resonance_history &
       (kingraph, res_hist, model, n_out)
    class(kingraph_t), intent(in), target :: kingraph
    type(resonance_history_t), intent(out) :: res_hist
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: n_out
    type(resonance_info_t) :: resonance
    integer :: i, mom_id, pdg
    if (debug_on)  call msg_debug2 &
         (D_PHASESPACE, "kingraph_extract_resonance_history")
    if (kingraph%grove_prop%n_resonances > 0) then
       if (associated (kingraph%root%daughter1) .or. &
            associated (kingraph%root%daughter2)) then
          if (debug_on)  call msg_debug2 &
               (D_PHASESPACE, "kingraph has resonances, root has children")
          do i = 1, kingraph%tree%n_entries
             if (kingraph%tree%mapping(i) == S_CHANNEL) then
                mom_id = kingraph%tree%bc (i)
                pdg = kingraph%tree%pdg (i)
                call resonance%init (mom_id, pdg, model, n_out)
                if (debug2_active (D_PHASESPACE)) then
                   print *, 'D: Adding resonance'
                   call resonance%write ()
                end if
                call res_hist%add_resonance (resonance)
             end if
          end do
       end if
    end if
  end subroutine kingraph_extract_resonance_history

  module function grove_list_get_n_trees (grove_list) result (n)
    class(grove_list_t), intent(in) :: grove_list
    integer :: n
    type(kingraph_t), pointer :: kingraph
    type(grove_t), pointer :: grove
    if (debug_on) call msg_debug (D_PHASESPACE, "grove_list_get_n_trees")
    n = 0
    grove => grove_list%first
    do while (associated (grove))
       kingraph => grove%first
       do while (associated (kingraph))
          if (kingraph%keep) n = n + 1
          kingraph => kingraph%grove_next
       end do
       grove => grove%next
    end do
    if (debug_on) call msg_debug (D_PHASESPACE, "n", n)
  end function grove_list_get_n_trees

  module subroutine feyngraph_set_get_resonance_histories &
       (feyngraph_set, n_filter, res_hists)
    type(feyngraph_set_t), intent(in), target :: feyngraph_set
    integer, intent(in), optional :: n_filter
    type(resonance_history_t), dimension(:), allocatable, intent(out) :: &
         res_hists
    type(kingraph_t), pointer :: kingraph
    type(grove_t), pointer :: grove
    type(resonance_history_t) :: res_hist
    type(resonance_history_set_t) :: res_hist_set
    integer :: i_grove
    if (debug_on)  call msg_debug &
         (D_PHASESPACE, "grove_list_get_resonance_histories")
    call res_hist_set%init (n_filter = n_filter)
    grove => feyngraph_set%grove_list%first
    i_grove = 0
    do while (associated (grove))
       i_grove = i_grove + 1
       kingraph => grove%first
       do while (associated (kingraph))
          if (kingraph%keep) then
             if (debug_on) call msg_debug2 (D_PHASESPACE, "grove", i_grove)
             call kingraph%extract_resonance_history &
                  (res_hist, feyngraph_set%model, feyngraph_set%n_out)
             call res_hist_set%enter (res_hist)
          end if
          kingraph => kingraph%grove_next
       end do
    end do
    call res_hist_set%freeze ()
    call res_hist_set%to_array (res_hists)
  end subroutine feyngraph_set_get_resonance_histories


end submodule cascades2_s

