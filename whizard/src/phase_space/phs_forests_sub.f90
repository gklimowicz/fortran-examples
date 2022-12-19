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

submodule (phs_forests) phs_forests_s

  use io_units
  use format_defs, only: FMT_19
  use diagnostics
  use numeric_utils
  use ifiles
  use lexers

  implicit none

contains

  module subroutine phs_parameters_write (phs_par, unit)
    class(phs_parameters_t), intent(in) :: phs_par
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A," // FMT_19 // ")") "sqrts         = ", phs_par%sqrts
    write (u, "(3x,A," // FMT_19 // ")") "m_threshold_s = ", phs_par%m_threshold_s
    write (u, "(3x,A," // FMT_19 // ")") "m_threshold_t = ", phs_par%m_threshold_t
    write (u, "(3x,A,I0)") "off_shell = ", phs_par%off_shell
    write (u, "(3x,A,I0)") "t_channel = ", phs_par%t_channel
    write (u, "(3x,A,L1)") "keep_nonresonant = ", phs_par%keep_nonresonant
  end subroutine phs_parameters_write

  module subroutine phs_parameters_read (phs_par, unit)
    class(phs_parameters_t), intent(out) :: phs_par
    integer, intent(in) :: unit
    character(20) :: dummy
    character :: equals
    read (unit, *)  dummy, equals, phs_par%sqrts
    read (unit, *)  dummy, equals, phs_par%m_threshold_s
    read (unit, *)  dummy, equals, phs_par%m_threshold_t
    read (unit, *)  dummy, equals, phs_par%off_shell
    read (unit, *)  dummy, equals, phs_par%t_channel
    read (unit, *)  dummy, equals, phs_par%keep_nonresonant
  end subroutine phs_parameters_read

  module function phs_parameters_eq (phs_par1, phs_par2) result (equal)
    logical :: equal
    type(phs_parameters_t), intent(in) :: phs_par1, phs_par2
    equal = phs_par1%sqrts == phs_par2%sqrts &
         .and. phs_par1%m_threshold_s == phs_par2%m_threshold_s &
         .and. phs_par1%m_threshold_t == phs_par2%m_threshold_t &
         .and. phs_par1%off_shell == phs_par2%off_shell &
         .and. phs_par1%t_channel == phs_par2%t_channel &
         .and.(phs_par1%keep_nonresonant .eqv. phs_par2%keep_nonresonant)
  end function phs_parameters_eq

  module function phs_parameters_ne (phs_par1, phs_par2) result (ne)
    logical :: ne
    type(phs_parameters_t), intent(in) :: phs_par1, phs_par2
    ne = phs_par1%sqrts /= phs_par2%sqrts &
         .or. phs_par1%m_threshold_s /= phs_par2%m_threshold_s &
         .or. phs_par1%m_threshold_t /= phs_par2%m_threshold_t &
         .or. phs_par1%off_shell /= phs_par2%off_shell &
         .or. phs_par1%t_channel /= phs_par2%t_channel &
         .or.(phs_par1%keep_nonresonant .neqv. phs_par2%keep_nonresonant)
  end function phs_parameters_ne

  subroutine equivalence_list_add (eql, left, right, perm)
    type(equivalence_list_t), intent(inout) :: eql
    integer, intent(in) :: left, right
    type(permutation_t), intent(in) :: perm
    type(equivalence_t), pointer :: eq
    allocate (eq)
    eq%left = left
    eq%right = right
    eq%perm = perm
    if (associated (eql%last)) then
       eql%last%next => eq
    else
       eql%first => eq
    end if
    eql%last => eq
    eql%length = eql%length + 1
  end subroutine equivalence_list_add

  pure subroutine equivalence_list_final (eql)
    type(equivalence_list_t), intent(inout) :: eql
    type(equivalence_t), pointer :: eq
    do while (associated (eql%first))
       eq => eql%first
       eql%first => eql%first%next
       deallocate (eq)
    end do
    eql%last => null ()
    eql%length = 0
  end subroutine equivalence_list_final

  elemental function equivalence_list_length (eql) result (length)
    integer :: length
    type(equivalence_list_t), intent(in) :: eql
    length = eql%length
  end function equivalence_list_length

  subroutine equivalence_list_write (eql, unit)
    type(equivalence_list_t), intent(in) :: eql
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    if (associated (eql%first)) then
       call equivalence_write_rec (eql%first, u)
    else
       write (u, *) " [empty]"
    end if
  contains
    recursive subroutine equivalence_write_rec (eq, u)
      type(equivalence_t), intent(in) :: eq
      integer, intent(in) :: u
      integer :: i
      write (u, "(3x,A,1x,I0,1x,I0,2x,A)", advance="no") &
           "Equivalence:", eq%left, eq%right, "Final state permutation:"
      call permutation_write (eq%perm, u)
      write (u, "(1x,12x,1x,A,1x)", advance="no") &
           "       msq permutation:  "
      call permutation_write (eq%msq_perm, u)
      write (u, "(1x,12x,1x,A,1x)", advance="no") &
           "       angle permutation:"
      call permutation_write (eq%angle_perm, u)
      write (u, "(1x,12x,1x,26x)", advance="no")
      do i = 1, size (eq%angle_sig)
         if (eq%angle_sig(i)) then
            write (u, "(1x,A)", advance="no") "+"
         else
            write (u, "(1x,A)", advance="no") "-"
         end if
      end do
      write (u, *)
      if (associated (eq%next))  call equivalence_write_rec (eq%next, u)
    end subroutine equivalence_write_rec
  end subroutine equivalence_list_write

  elemental subroutine phs_grove_init &
       (grove, n_trees, n_in, n_out, n_masses, n_angles)
    type(phs_grove_t), intent(inout) :: grove
    integer, intent(in) :: n_trees, n_in, n_out, n_masses, n_angles
    grove%tree_count_offset = 0
    allocate (grove%tree (n_trees))
    call grove%tree%init (n_in, n_out, n_masses, n_angles)
  end subroutine phs_grove_init

  elemental subroutine phs_grove_final (grove)
    type(phs_grove_t), intent(inout) :: grove
    deallocate (grove%tree)
    call equivalence_list_final (grove%equivalence_list)
  end subroutine phs_grove_final

  subroutine phs_grove_assign_s_mappings (grove, mapping)
    type(phs_grove_t), intent(in) :: grove
    type(mapping_t), dimension(:), intent(out) :: mapping
    integer :: i
    if (size (mapping) == size (grove%tree)) then
       do i = 1, size (mapping)
          call grove%tree(i)%assign_s_mapping (mapping(i))
       end do
    else
       call msg_bug ("phs_grove_assign_s_mappings: array size mismatch")
    end if
  end subroutine phs_grove_assign_s_mappings

  module subroutine phs_forest_init (forest, n_tree, n_in, n_out)
    class(phs_forest_t), intent(inout) :: forest
    integer, dimension(:), intent(in) :: n_tree
    integer, intent(in) :: n_in, n_out
    integer :: g, count, k_root
    forest%n_in = n_in
    forest%n_out = n_out
    forest%n_tot = n_in + n_out
    forest%n_masses = max (n_out - 2, 0)
    forest%n_angles = max (2*n_out - 2, 0)
    forest%n_dimensions = forest%n_masses + forest%n_angles
    forest%n_trees = sum (n_tree)
    forest%n_equivalences = 0
    allocate (forest%grove (size (n_tree)))
    call phs_grove_init &
         (forest%grove, n_tree, n_in, n_out, forest%n_masses, &
          forest%n_angles)
    allocate (forest%grove_lookup (forest%n_trees))
    count = 0
    do g = 1, size (forest%grove)
       forest%grove(g)%tree_count_offset = count
       forest%grove_lookup (count+1:count+n_tree(g)) = g
       count = count + n_tree(g)
    end do
    allocate (forest%prt_in  (n_in))
    allocate (forest%prt_out (forest%n_out))
    k_root = 2**forest%n_tot - 1
    allocate (forest%prt (k_root))
    allocate (forest%prt_combination (2, k_root))
    allocate (forest%s_mapping (forest%n_trees))
  end subroutine phs_forest_init

  module subroutine phs_forest_set_s_mappings (forest)
    class(phs_forest_t), intent(inout) :: forest
    integer :: g, i0, i1, n
    do g = 1, size (forest%grove)
       call forest%get_grove_bounds (g, i0, i1, n)
       call phs_grove_assign_s_mappings &
            (forest%grove(g), forest%s_mapping(i0:i1))
    end do
  end subroutine phs_forest_set_s_mappings

  module subroutine phs_forest_final (forest)
    class(phs_forest_t), intent(inout) :: forest
    if (allocated (forest%grove)) then
       call phs_grove_final (forest%grove)
       deallocate (forest%grove)
    end if
    if (allocated (forest%grove_lookup))  deallocate (forest%grove_lookup)
    if (allocated (forest%prt))  deallocate (forest%prt)
    if (allocated (forest%s_mapping))  deallocate (forest%s_mapping)
  end subroutine phs_forest_final

  module subroutine phs_forest_write (forest, unit)
    class(phs_forest_t), intent(in) :: forest
    integer, intent(in), optional :: unit
    integer :: u
    integer :: i, g, k
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)") "Phase space forest:"
    write (u, "(3x,A,I0)") "n_in  = ", forest%n_in
    write (u, "(3x,A,I0)") "n_out = ", forest%n_out
    write (u, "(3x,A,I0)") "n_tot = ", forest%n_tot
    write (u, "(3x,A,I0)") "n_masses = ", forest%n_masses
    write (u, "(3x,A,I0)") "n_angles = ", forest%n_angles
    write (u, "(3x,A,I0)") "n_dim    = ", forest%n_dimensions
    write (u, "(3x,A,I0)") "n_trees  = ", forest%n_trees
    write (u, "(3x,A,I0)") "n_equiv  = ", forest%n_equivalences
    write (u, "(3x,A)", advance="no") "flavors  ="
    if (allocated (forest%flv)) then
       do i = 1, size (forest%flv)
          write (u, "(1x,I0)", advance="no")  forest%flv(i)%get_pdg ()
       end do
       write (u, "(A)")
    else
       write (u, "(1x,A)") "[empty]"
    end if
    write (u, "(1x,A)") "Particle combinations:"
    if (allocated (forest%prt_combination)) then
       do k = 1, size (forest%prt_combination, 2)
          if (forest%prt_combination(1, k) /= 0) then
             write (u, "(3x,I0,1x,'<=',1x,I0,1x,'+',1x,I0)") &
                  k, forest%prt_combination(:,k)
          end if
       end do
    else
       write (u, "(3x,A)") "  [empty]"
    end if
    write (u, "(1x,A)") "Groves and trees:"
    if (allocated (forest%grove)) then
       do g = 1, size (forest%grove)
          write (u, "(3x,A,1x,I0)") "Grove    ", g
          call phs_grove_write (forest%grove(g), unit)
       end do
    else
       write (u, "(3x,A)") "  [empty]"
    end if
    write (u, "(1x,A,I0)") "Total number of equivalences: ", &
         forest%n_equivalences
    write (u, "(A)")
    write (u, "(1x,A)") "Global s-channel mappings:"
    if (allocated (forest%s_mapping)) then
       do i = 1, size (forest%s_mapping)
          associate (mapping => forest%s_mapping(i))
            if (mapping%is_s_channel () .or. mapping%is_on_shell ()) then
               write (u, "(1x,I0,':',1x)", advance="no")  i
               call forest%s_mapping(i)%write (unit)
            end if
          end associate
       end do
    else
       write (u, "(3x,A)") "  [empty]"
    end if
    write (u, "(A)")
    write (u, "(1x,A)") "Incoming particles:"
    if (allocated (forest%prt_in)) then
       if (any (forest%prt_in%is_defined ())) then
          do i = 1, size (forest%prt_in)
             if (forest%prt_in(i)%is_defined ()) then
                write (u, "(1x,A,1x,I0)")  "Particle", i
                call forest%prt_in(i)%write (u)
             end if
          end do
       else
          write (u, "(3x,A)")  "[all undefined]"
       end if
    else
       write (u, "(3x,A)")  "  [empty]"
    end if
    write (u, "(A)")
    write (u, "(1x,A)") "Outgoing particles:"
    if (allocated (forest%prt_out)) then
       if (any (forest%prt_out%is_defined ())) then
          do i = 1, size (forest%prt_out)
             if (forest%prt_out(i)%is_defined ()) then
                write (u, "(1x,A,1x,I0)")  "Particle", i
                call forest%prt_out(i)%write (u)
             end if
          end do
       else
          write (u, "(3x,A)")  "[all undefined]"
       end if
    else
       write (u, "(1x,A)")  "  [empty]"
    end if
    write (u, "(A)")
    write (u, "(1x,A)") "Tree particles:"
    if (allocated (forest%prt)) then
       if (any (forest%prt%is_defined ())) then
          do i = 1, size (forest%prt)
             if (forest%prt(i)%is_defined ()) then
                write (u, "(1x,A,1x,I0)")  "Particle", i
                call forest%prt(i)%write (u)
             end if
          end do
       else
          write (u, "(3x,A)")  "[all undefined]"
       end if
    else
       write (u, "(3x,A)")  "  [empty]"
    end if
  end subroutine phs_forest_write

  subroutine phs_grove_write (grove, unit)
    type(phs_grove_t), intent(in) :: grove
    integer, intent(in), optional :: unit
    integer :: u
    integer :: t
    u = given_output_unit (unit);  if (u < 0)  return
    do t = 1, size (grove%tree)
       write (u, "(3x,A,I0)") "Tree      ", t
       call grove%tree(t)%write (unit)
    end do
    write (u, "(1x,A)") "Equivalence list:"
    call equivalence_list_write (grove%equivalence_list, unit)
  end subroutine phs_grove_write

  module subroutine phs_forest_assign (forest_out, forest_in)
    type(phs_forest_t), intent(out) :: forest_out
    type(phs_forest_t), intent(in) :: forest_in
    forest_out%n_in  = forest_in%n_in
    forest_out%n_out = forest_in%n_out
    forest_out%n_tot = forest_in%n_tot
    forest_out%n_masses = forest_in%n_masses
    forest_out%n_angles = forest_in%n_angles
    forest_out%n_dimensions  = forest_in%n_dimensions
    forest_out%n_trees  = forest_in%n_trees
    forest_out%n_equivalences  = forest_in%n_equivalences
    if (allocated (forest_in%flv)) then
       allocate (forest_out%flv (size (forest_in%flv)))
       forest_out%flv = forest_in%flv
    end if
    if (allocated (forest_in%grove)) then
       allocate (forest_out%grove (size (forest_in%grove)))
       forest_out%grove = forest_in%grove
    end if
    if (allocated (forest_in%grove_lookup)) then
       allocate (forest_out%grove_lookup (size (forest_in%grove_lookup)))
       forest_out%grove_lookup = forest_in%grove_lookup
    end if
    if (allocated (forest_in%prt_in)) then
       allocate (forest_out%prt_in (size (forest_in%prt_in)))
       forest_out%prt_in = forest_in%prt_in
    end if
    if (allocated (forest_in%prt_out)) then
       allocate (forest_out%prt_out (size (forest_in%prt_out)))
       forest_out%prt_out = forest_in%prt_out
    end if
    if (allocated (forest_in%prt)) then
       allocate (forest_out%prt (size (forest_in%prt)))
       forest_out%prt = forest_in%prt
    end if
    if (allocated (forest_in%s_mapping)) then
       allocate (forest_out%s_mapping (size (forest_in%s_mapping)))
       forest_out%s_mapping = forest_in%s_mapping
    end if
    if (allocated (forest_in%prt_combination)) then
       allocate (forest_out%prt_combination &
            (2, size (forest_in%prt_combination, 2)))
       forest_out%prt_combination = forest_in%prt_combination
    end if
  end subroutine phs_forest_assign

  module function phs_forest_get_n_parameters (forest) result (n)
    integer :: n
    class(phs_forest_t), intent(in) :: forest
    n = forest%n_dimensions
  end function phs_forest_get_n_parameters

  module function phs_forest_get_n_channels (forest) result (n)
    integer :: n
    class(phs_forest_t), intent(in) :: forest
    n = forest%n_trees
  end function phs_forest_get_n_channels

  module function phs_forest_get_n_groves (forest) result (n)
    integer :: n
    class(phs_forest_t), intent(in) :: forest
    n = size (forest%grove)
  end function phs_forest_get_n_groves

  module subroutine phs_forest_get_grove_bounds (forest, g, i0, i1, n)
    class(phs_forest_t), intent(in) :: forest
    integer, intent(in) :: g
    integer, intent(out) :: i0, i1, n
    n = size (forest%grove(g)%tree)
    i0 = forest%grove(g)%tree_count_offset + 1
    i1 = forest%grove(g)%tree_count_offset + n
  end subroutine phs_forest_get_grove_bounds

  module function phs_forest_get_n_equivalences (forest) result (n)
    integer :: n
    class(phs_forest_t), intent(in) :: forest
    n = forest%n_equivalences
  end function phs_forest_get_n_equivalences

  module subroutine phs_forest_get_s_mapping &
       (forest, channel, flag, mass, width)
    class(phs_forest_t), intent(in) :: forest
    integer, intent(in) :: channel
    logical, intent(out) :: flag
    real(default), intent(out) :: mass, width
    flag = forest%s_mapping(channel)%is_s_channel ()
    if (flag) then
       mass = forest%s_mapping(channel)%get_mass ()
       width = forest%s_mapping(channel)%get_width ()
    else
       mass = 0
       width = 0
    end if
  end subroutine phs_forest_get_s_mapping

  module subroutine phs_forest_get_on_shell (forest, channel, flag, mass)
    class(phs_forest_t), intent(in) :: forest
    integer, intent(in) :: channel
    logical, intent(out) :: flag
    real(default), intent(out) :: mass
    flag = forest%s_mapping(channel)%is_on_shell ()
    if (flag) then
       mass = forest%s_mapping(channel)%get_mass ()
    else
       mass = 0
    end if
  end subroutine phs_forest_get_on_shell

  module subroutine phs_forest_extract_resonance_history_set &
       (forest, res_set, include_trivial)
    class(phs_forest_t), intent(in) :: forest
    type(resonance_history_set_t), intent(out) :: res_set
    logical, intent(in), optional :: include_trivial
    type(resonance_history_t) :: rh
    integer :: g, t
    logical :: triv
    triv = .false.;  if (present (include_trivial))  triv = include_trivial
    call res_set%init ()
    do g = 1, size (forest%grove)
       associate (grove => forest%grove(g))
          do t = 1, size (grove%tree)
             call grove%tree(t)%extract_resonance_history (rh)
             call res_set%enter (rh, include_trivial)
          end do
       end associate
    end do
    call res_set%freeze ()
  end subroutine phs_forest_extract_resonance_history_set

  subroutine define_phs_forest_syntax (ifile)
    type(ifile_t) :: ifile
    call ifile_append (ifile, "SEQ phase_space_list = process_phase_space*")
    call ifile_append (ifile, "SEQ process_phase_space = " &
         // "process_def process_header phase_space")
    call ifile_append (ifile, "SEQ process_def = process process_list")
    call ifile_append (ifile, "KEY process")
    call ifile_append (ifile, "LIS process_list = process_tag*")
    call ifile_append (ifile, "IDE process_tag")
    call ifile_append (ifile, "SEQ process_header = " &
         // "md5sum_process = md5sum " &
         // "md5sum_model_par = md5sum " &
         // "md5sum_phs_config = md5sum " &
         // "sqrts = real " &
         // "m_threshold_s = real " &
         // "m_threshold_t = real " &
         // "off_shell = integer " &
         // "t_channel = integer " &
         // "keep_nonresonant = logical")
    call ifile_append (ifile, "KEY '='")
    call ifile_append (ifile, "KEY '-'")
    call ifile_append (ifile, "KEY md5sum_process")
    call ifile_append (ifile, "KEY md5sum_model_par")
    call ifile_append (ifile, "KEY md5sum_phs_config")
    call ifile_append (ifile, "KEY sqrts")
    call ifile_append (ifile, "KEY m_threshold_s")
    call ifile_append (ifile, "KEY m_threshold_t")
    call ifile_append (ifile, "KEY off_shell")
    call ifile_append (ifile, "KEY t_channel")
    call ifile_append (ifile, "KEY keep_nonresonant")
    call ifile_append (ifile, "QUO md5sum = '""' ... '""'")
    call ifile_append (ifile, "REA real")
    call ifile_append (ifile, "INT integer")
    call ifile_append (ifile, "IDE logical")
    call ifile_append (ifile, "SEQ phase_space = grove_def+")
    call ifile_append (ifile, "SEQ grove_def = grove tree_def+")
    call ifile_append (ifile, "KEY grove")
    call ifile_append (ifile, "SEQ tree_def = tree bincodes mapping*")
    call ifile_append (ifile, "KEY tree")
    call ifile_append (ifile, "SEQ bincodes = bincode*")
    call ifile_append (ifile, "INT bincode")
    call ifile_append (ifile, "SEQ mapping = map bincode channel signed_pdg")
    call ifile_append (ifile, "KEY map")
    call ifile_append (ifile, "ALT channel = &
         &s_channel | t_channel | u_channel | &
         &collinear | infrared | radiation | on_shell")
    call ifile_append (ifile, "KEY s_channel")
    ! call ifile_append (ifile, "KEY t_channel")   !!! Key already exists
    call ifile_append (ifile, "KEY u_channel")
    call ifile_append (ifile, "KEY collinear")
    call ifile_append (ifile, "KEY infrared")
    call ifile_append (ifile, "KEY radiation")
    call ifile_append (ifile, "KEY on_shell")
    call ifile_append (ifile, "ALT signed_pdg = &
         &pdg | negative_pdg")
    call ifile_append (ifile, "SEQ negative_pdg = '-' pdg")
    call ifile_append (ifile, "INT pdg")
  end subroutine define_phs_forest_syntax

  module subroutine syntax_phs_forest_init ()
    type(ifile_t) :: ifile
    call define_phs_forest_syntax (ifile)
    call syntax_init (syntax_phs_forest, ifile)
    call ifile_final (ifile)
  end subroutine syntax_phs_forest_init

  subroutine lexer_init_phs_forest (lexer)
    type(lexer_t), intent(out) :: lexer
    call lexer_init (lexer, &
         comment_chars = "#!", &
         quote_chars = '"', &
         quote_match = '"', &
         single_chars = "-", &
         special_class = ["="] , &
         keyword_list = syntax_get_keyword_list_ptr (syntax_phs_forest))
  end subroutine lexer_init_phs_forest

  module subroutine syntax_phs_forest_final ()
    call syntax_final (syntax_phs_forest)
  end subroutine syntax_phs_forest_final

  module subroutine syntax_phs_forest_write (unit)
    integer, intent(in), optional :: unit
    call syntax_write (syntax_phs_forest, unit)
  end subroutine syntax_phs_forest_write

  module subroutine phs_forest_read_file &
       (forest, filename, process_id, n_in, n_out, model, found, &
        md5sum_process, md5sum_model_par, &
        md5sum_phs_config, phs_par, match)
    class(phs_forest_t), intent(out) :: forest
    type(string_t), intent(in) :: filename
    type(string_t), intent(in) :: process_id
    integer, intent(in) :: n_in, n_out
    class(model_data_t), intent(in), target :: model
    logical, intent(out) :: found
    character(32), intent(in), optional :: &
         md5sum_process, md5sum_model_par, md5sum_phs_config
    type(phs_parameters_t), intent(in), optional :: phs_par
    logical, intent(out), optional :: match
    type(parse_tree_t), target :: parse_tree
    type(stream_t), target :: stream
    type(lexer_t) :: lexer
    call lexer_init_phs_forest (lexer)
    call stream_init (stream, char (filename))
    call lexer_assign_stream (lexer, stream)
    call parse_tree_init (parse_tree, syntax_phs_forest, lexer)
    call phs_forest_read_parse_tree (forest, parse_tree, &
         process_id, n_in, n_out, model, found, &
         md5sum_process, md5sum_model_par, md5sum_phs_config, phs_par, match)
    call stream_final (stream)
    call lexer_final (lexer)
    call parse_tree_final (parse_tree)
  end subroutine phs_forest_read_file

  module subroutine phs_forest_read_unit &
       (forest, unit, process_id, n_in, n_out, model, found, &
        md5sum_process, md5sum_model_par, md5sum_phs_config, &
        phs_par, match)
    class(phs_forest_t), intent(out) :: forest
    integer, intent(in) :: unit
    type(string_t), intent(in) :: process_id
    integer, intent(in) :: n_in, n_out
    class(model_data_t), intent(in), target :: model
    logical, intent(out) :: found
    character(32), intent(in), optional :: &
         md5sum_process, md5sum_model_par, md5sum_phs_config
    type(phs_parameters_t), intent(in), optional :: phs_par
    logical, intent(out), optional :: match
    type(parse_tree_t), target :: parse_tree
    type(stream_t), target :: stream
    type(lexer_t) :: lexer
    call lexer_init_phs_forest (lexer)
    call stream_init (stream, unit)
    call lexer_assign_stream (lexer, stream)
    call parse_tree_init (parse_tree, syntax_phs_forest, lexer)
    call phs_forest_read_parse_tree (forest, parse_tree, &
         process_id, n_in, n_out, model, found, &
         md5sum_process, md5sum_model_par, md5sum_phs_config, &
         phs_par, match)
    call stream_final (stream)
    call lexer_final (lexer)
    call parse_tree_final (parse_tree)
  end subroutine phs_forest_read_unit

  module subroutine phs_forest_read_parse_tree &
       (forest, parse_tree, process_id, n_in, n_out, model, found, &
        md5sum_process, md5sum_model_par, md5sum_phs_config, &
        phs_par, match)
    class(phs_forest_t), intent(out) :: forest
    type(parse_tree_t), intent(in), target :: parse_tree
    type(string_t), intent(in) :: process_id
    integer, intent(in) :: n_in, n_out
    class(model_data_t), intent(in), target :: model
    logical, intent(out) :: found
    character(32), intent(in), optional :: &
         md5sum_process, md5sum_model_par, md5sum_phs_config
    type(phs_parameters_t), intent(in), optional :: phs_par
    logical, intent(out), optional :: match
    type(parse_node_t), pointer :: node_header, node_phs, node_grove
    integer :: n_grove, g
    integer, dimension(:), allocatable :: n_tree
    integer :: t
    node_header => parse_tree_get_process_ptr (parse_tree, process_id)
    found = associated (node_header);  if (.not. found)  return
    if (present (match)) then
       call phs_forest_check_input (node_header, &
            md5sum_process, md5sum_model_par, md5sum_phs_config, phs_par, match)
       if (.not. match)  return
    end if
    node_phs => parse_node_get_next_ptr (node_header)
    n_grove = parse_node_get_n_sub (node_phs)
    allocate (n_tree (n_grove))
    do g = 1, n_grove
       node_grove => parse_node_get_sub_ptr (node_phs, g)
       n_tree(g) = parse_node_get_n_sub (node_grove) - 1
    end do
    call forest%init (n_tree, n_in, n_out)
    do g = 1, n_grove
       node_grove => parse_node_get_sub_ptr (node_phs, g)
       do t = 1, n_tree(g)
          call phs_tree_set (forest%grove(g)%tree(t), &
               parse_node_get_sub_ptr (node_grove, t+1), model)
       end do
    end do
  end subroutine phs_forest_read_parse_tree

  subroutine phs_forest_check_input (pn_header, &
       md5sum_process, md5sum_model_par, md5sum_phs_config, phs_par, match)
    type(parse_node_t), intent(in), target :: pn_header
    character(32), intent(in) :: &
         md5sum_process, md5sum_model_par, md5sum_phs_config
    type(phs_parameters_t), intent(in), optional :: phs_par
    logical, intent(out) :: match
    type(parse_node_t), pointer :: pn_md5sum, pn_rval, pn_ival, pn_lval
    character(32) :: md5sum
    type(phs_parameters_t) :: phs_par_old
    character(1) :: lstr
    pn_md5sum => parse_node_get_sub_ptr (pn_header, 3)
    md5sum = parse_node_get_string (pn_md5sum)
    if (md5sum /= "" .and. md5sum /= md5sum_process) then
       call msg_message ("Phase space: discarding old configuration &
            &(process changed)")
       match = .false.;  return
    end if
    pn_md5sum => parse_node_get_next_ptr (pn_md5sum, 3)
    md5sum = parse_node_get_string (pn_md5sum)
    if (md5sum /= "" .and. md5sum /= md5sum_model_par) then
       call msg_message ("Phase space: discarding old configuration &
            &(model parameters changed)")
       match = .false.;  return
    end if
    pn_md5sum => parse_node_get_next_ptr (pn_md5sum, 3)
    md5sum = parse_node_get_string (pn_md5sum)
    if (md5sum /= "" .and. md5sum /= md5sum_phs_config) then
       call msg_message ("Phase space: discarding old configuration &
            &(configuration parameters changed)")
       match = .false.;  return
    end if
    if (present (phs_par)) then
       pn_rval => parse_node_get_next_ptr (pn_md5sum, 3)
       phs_par_old%sqrts = parse_node_get_real (pn_rval)
       pn_rval => parse_node_get_next_ptr (pn_rval, 3)
       phs_par_old%m_threshold_s = parse_node_get_real (pn_rval)
       pn_rval => parse_node_get_next_ptr (pn_rval, 3)
       phs_par_old%m_threshold_t = parse_node_get_real (pn_rval)
       pn_ival => parse_node_get_next_ptr (pn_rval, 3)
       phs_par_old%off_shell = parse_node_get_integer (pn_ival)
       pn_ival => parse_node_get_next_ptr (pn_ival, 3)
       phs_par_old%t_channel = parse_node_get_integer (pn_ival)
       pn_lval => parse_node_get_next_ptr (pn_ival, 3)
       lstr = parse_node_get_string (pn_lval)
       read (lstr, "(L1)")  phs_par_old%keep_nonresonant
       if (phs_par_old /= phs_par) then
          call msg_message &
               ("Phase space: discarding old configuration &
               &(configuration parameters changed)")
          match = .false.;  return
       end if
    end if
    match = .true.
  end subroutine phs_forest_check_input

  subroutine phs_tree_set (tree, node, model)
    type(phs_tree_t), intent(inout) :: tree
    type(parse_node_t), intent(in), target :: node
    class(model_data_t), intent(in), target :: model
    type(parse_node_t), pointer :: node_bincodes, node_mapping, pn_pdg
    integer :: n_bincodes, offset
    integer(TC), dimension(:), allocatable :: bincode
    integer :: b, n_mappings, m
    integer(TC) :: k
    type(string_t) :: type
    integer :: pdg
    node_bincodes => parse_node_get_sub_ptr (node, 2)
    if (associated (node_bincodes)) then
       select case (char (parse_node_get_rule_key (node_bincodes)))
       case ("bincodes")
          n_bincodes = parse_node_get_n_sub (node_bincodes)
          offset = 2
       case default
          n_bincodes = 0
          offset = 1
       end select
    else
       n_bincodes = 0
       offset = 2
    end if
    allocate (bincode (n_bincodes))
    do b = 1, n_bincodes
       bincode(b) = parse_node_get_integer &
            (parse_node_get_sub_ptr (node_bincodes, b))
    end do
    call phs_tree_from_array (tree, bincode)
    call tree%flip_t_to_s_channel ()
    call tree%canonicalize ()
    n_mappings = parse_node_get_n_sub (node) - offset
    do m = 1, n_mappings
       node_mapping => parse_node_get_sub_ptr (node, m + offset)
       k = parse_node_get_integer &
            (parse_node_get_sub_ptr (node_mapping, 2))
       type = parse_node_get_key &
            (parse_node_get_sub_ptr (node_mapping, 3))
       pn_pdg => parse_node_get_sub_ptr (node_mapping, 4)
       select case (char (pn_pdg%get_rule_key ()))
       case ("pdg")
          pdg = pn_pdg%get_integer ()
       case ("negative_pdg")
          pdg = - parse_node_get_integer (pn_pdg%get_sub_ptr (2))
       end select
       call tree%init_mapping (k, type, pdg, model)
    end do
  end subroutine phs_tree_set

  module subroutine phs_forest_set_flavors (forest, flv, reshuffle, flv_extra)
    class(phs_forest_t), intent(inout) :: forest
    type(flavor_t), dimension(:), intent(in) :: flv
    integer, intent(in), dimension(:), allocatable, optional :: reshuffle
    type(flavor_t), intent(in), optional :: flv_extra
    integer :: i, n_flv0
    if (present (reshuffle) .and. present (flv_extra)) then
       n_flv0 = size (flv)
       do i = 1, n_flv0
          if (reshuffle(i) <= n_flv0) then
             forest%flv(i) = flv (reshuffle(i))
          else
             forest%flv(i) = flv_extra
          end if
       end do
    else
       allocate (forest%flv (size (flv)))
       forest%flv = flv
    end if
  end subroutine phs_forest_set_flavors

  module subroutine phs_forest_set_momentum_links (forest, list)
    class(phs_forest_t), intent(inout) :: forest
    integer, intent(in), dimension(:), allocatable :: list
    integer :: g, t
    do g = 1, size (forest%grove)
      do t = 1, size (forest%grove(g)%tree)
        associate (tree => forest%grove(g)%tree(t))
          call phs_tree_set_momentum_links (tree, list)
          !!! call tree%reshuffle_mappings ()
        end associate
      end do
    end do
  end subroutine phs_forest_set_momentum_links

  module subroutine phs_forest_set_parameters &
       (forest, mapping_defaults, variable_limits)
    class(phs_forest_t), intent(inout) :: forest
    type(mapping_defaults_t), intent(in) :: mapping_defaults
    logical, intent(in) :: variable_limits
    integer :: g, t
    do g = 1, size (forest%grove)
       do t = 1, size (forest%grove(g)%tree)
          call forest%grove(g)%tree(t)%set_mass_sum (forest%flv(forest%n_in+1:))
          call forest%grove(g)%tree(t)%set_mapping_parameters &
               (mapping_defaults, variable_limits)
          call forest%grove(g)%tree(t)%set_effective_masses ()
          if (mapping_defaults%step_mapping) then
             call forest%grove(g)%tree(t)%set_step_mappings &
                  (mapping_defaults%step_mapping_exp, variable_limits)
          end if
       end do
    end do
  end subroutine phs_forest_set_parameters

  module subroutine phs_forest_setup_prt_combinations (forest)
    class(phs_forest_t), intent(inout) :: forest
    integer :: g, t
    integer, dimension(:,:), allocatable :: tree_prt_combination
    forest%prt_combination = 0
    allocate (tree_prt_combination (2, size (forest%prt_combination, 2)))
    do g = 1, size (forest%grove)
       do t = 1, size (forest%grove(g)%tree)
          call phs_tree_setup_prt_combinations &
               (forest%grove(g)%tree(t), tree_prt_combination)
          where (tree_prt_combination /= 0 .and. forest%prt_combination == 0)
             forest%prt_combination = tree_prt_combination
          end where
       end do
    end do
  end subroutine phs_forest_setup_prt_combinations

  module subroutine phs_forest_set_prt_in_int (forest, int, lt_cm_to_lab)
    class(phs_forest_t), intent(inout) :: forest
    type(interaction_t), intent(in) :: int
    type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    if (present (lt_cm_to_lab)) then
       call forest%prt_in%set_momentum (inverse (lt_cm_to_lab) * &
            int%get_momenta (outgoing=.false.))
    else
       call forest%prt_in%set_momentum (int%get_momenta (outgoing=.false.))
    end if
    associate (m_in => forest%flv(:forest%n_in)%get_mass ())
      call forest%prt_in%set_msq (m_in ** 2)
    end associate
    call forest%prt_in%set_defined ()
  end subroutine phs_forest_set_prt_in_int

  module subroutine phs_forest_set_prt_in_mom (forest, mom, lt_cm_to_lab)
    class(phs_forest_t), intent(inout) :: forest
    type(vector4_t), dimension(size (forest%prt_in)), intent(in) :: mom
    type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    if (present (lt_cm_to_lab)) then
       call forest%prt_in%set_momentum (inverse (lt_cm_to_lab) * mom)
    else
       call forest%prt_in%set_momentum (mom)
    end if
    associate (m_in => forest%flv(:forest%n_in)%get_mass ())
      call forest%prt_in%set_msq (m_in ** 2)
    end associate
    call forest%prt_in%set_defined ()
  end subroutine phs_forest_set_prt_in_mom

  module subroutine phs_forest_set_prt_out_int (forest, int, lt_cm_to_lab)
    class(phs_forest_t), intent(inout) :: forest
    type(interaction_t), intent(in) :: int
    type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    if (present (lt_cm_to_lab)) then
       call forest%prt_out%set_momentum (inverse (lt_cm_to_lab) * &
            int%get_momenta (outgoing=.true.))
    else
       call forest%prt_out%set_momentum (int%get_momenta (outgoing=.true.))
    end if
    associate (m_out => forest%flv(forest%n_in+1:)%get_mass ())
      call forest%prt_out%set_msq (m_out ** 2)
    end associate
    call forest%prt_out%set_defined ()
  end subroutine phs_forest_set_prt_out_int

  module subroutine phs_forest_set_prt_out_mom (forest, mom, lt_cm_to_lab)
    class(phs_forest_t), intent(inout) :: forest
    type(vector4_t), dimension(size (forest%prt_out)), intent(in) :: mom
    type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    if (present (lt_cm_to_lab)) then
       call forest%prt_out%set_momentum (inverse (lt_cm_to_lab) * mom)
    else
       call forest%prt_out%set_momentum (mom)
    end if
    associate (m_out => forest%flv(forest%n_in+1:)%get_mass ())
      call forest%prt_out%set_msq (m_out ** 2)
    end associate
    call forest%prt_out%set_defined ()
  end subroutine phs_forest_set_prt_out_mom

  module subroutine phs_forest_combine_particles (forest)
    class(phs_forest_t), intent(inout) :: forest
    integer :: k
    integer, dimension(2) :: kk
    do k = 1, size (forest%prt_combination, 2)
       kk = forest%prt_combination(:,k)
       if (kk(1) /= 0) then
          call forest%prt(k)%combine (forest%prt(kk(1)), forest%prt(kk(2)))
       end if
    end do
  end subroutine phs_forest_combine_particles

  module subroutine phs_forest_get_prt_out (forest, int, lt_cm_to_lab)
    class(phs_forest_t), intent(in) :: forest
    type(interaction_t), intent(inout) :: int
    type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    if (present (lt_cm_to_lab)) then
       call int%set_momenta (lt_cm_to_lab * &
            forest%prt_out%get_momentum (), outgoing=.true.)
    else
       call int%set_momenta (forest%prt_out%get_momentum (), &
            outgoing=.true.)
    end if
  end subroutine phs_forest_get_prt_out

  module function phs_forest_get_momenta_out (forest, lt_cm_to_lab) result (p)
    class(phs_forest_t), intent(in) :: forest
    type(lorentz_transformation_t), intent(in), optional :: lt_cm_to_lab
    type(vector4_t), dimension(size (forest%prt_out)) :: p
    p = forest%prt_out%get_momentum ()
    if (present (lt_cm_to_lab)) p = p * lt_cm_to_lab
  end function phs_forest_get_momenta_out

  subroutine phs_grove_set_equivalences (grove, perm_array)
    type(phs_grove_t), intent(inout) :: grove
    type(permutation_t), dimension(:), intent(in) :: perm_array
    type(equivalence_t), pointer :: eq
    integer :: t1, t2, i
    do t1 = 1, size (grove%tree)
       do t2 = 1, size (grove%tree)
          SCAN_PERM: do i = 1, size (perm_array)
             if (phs_tree_equivalent &
                  (grove%tree(t1), grove%tree(t2), perm_array(i))) then
                call equivalence_list_add &
                     (grove%equivalence_list, t1, t2, perm_array(i))
                eq => grove%equivalence_list%last
                call phs_tree_find_msq_permutation &
                     (grove%tree(t1), grove%tree(t2), eq%perm, &
                      eq%msq_perm)
                call phs_tree_find_angle_permutation &
                     (grove%tree(t1), grove%tree(t2), eq%perm, &
                      eq%angle_perm, eq%angle_sig)
             end if
          end do SCAN_PERM
       end do
    end do
  end subroutine phs_grove_set_equivalences

  module subroutine phs_forest_set_equivalences (forest)
    class(phs_forest_t), intent(inout) :: forest
    type(permutation_t), dimension(:), allocatable :: perm_array
    integer :: i
    call permutation_array_make &
         (perm_array, forest%flv(forest%n_in+1:)%get_pdg ())
    do i = 1, size (forest%grove)
       call phs_grove_set_equivalences (forest%grove(i), perm_array)
    end do
    forest%n_equivalences = sum (forest%grove%equivalence_list%length)
  end subroutine phs_forest_set_equivalences

  module subroutine phs_forest_get_equivalences &
       (forest, channel, azimuthal_dependence)
    class(phs_forest_t), intent(in) :: forest
    type(phs_channel_t), dimension(:), intent(out) :: channel
    logical, intent(in) :: azimuthal_dependence
    integer :: n_masses, n_angles
    integer :: mode_azimuthal_angle
    integer, dimension(:), allocatable :: n_eq
    type(equivalence_t), pointer :: eq
    integer, dimension(:), allocatable :: perm, mode
    integer :: g, c, j, left, right
    n_masses = forest%n_masses
    n_angles = forest%n_angles
    allocate (n_eq (forest%n_trees), source = 0)
    allocate (perm (forest%n_dimensions))
    allocate (mode (forest%n_dimensions), source = EQ_IDENTITY)
    do g = 1, size (forest%grove)
       eq => forest%grove(g)%equivalence_list%first
       do while (associated (eq))
          left = eq%left + forest%grove(g)%tree_count_offset
          n_eq(left) = n_eq(left) + 1
          eq => eq%next
       end do
    end do
    do c = 1, size (channel)
       allocate (channel(c)%eq (n_eq(c)))
       do j = 1, n_eq(c)
          call channel(c)%eq(j)%init (forest%n_dimensions)
       end do
    end do
    n_eq = 0
    if (azimuthal_dependence) then
       mode_azimuthal_angle = EQ_IDENTITY
    else
       mode_azimuthal_angle = EQ_INVARIANT
    end if
    do g = 1, size (forest%grove)
       eq => forest%grove(g)%equivalence_list%first
       do while (associated (eq))
          left = eq%left + forest%grove(g)%tree_count_offset
          right = eq%right + forest%grove(g)%tree_count_offset
          do j = 1, n_masses
             perm(j) = permute (j, eq%msq_perm)
             mode(j) = EQ_IDENTITY
          end do
          do j = 1, n_angles
             perm(n_masses+j) = n_masses + permute (j, eq%angle_perm)
             if (j == 1) then
                mode(n_masses+j) = mode_azimuthal_angle   ! first az. angle
             else if (mod(j,2) == 1) then
                mode(n_masses+j) = EQ_SYMMETRIC          ! other az. angles
             else if (eq%angle_sig(j)) then
                mode(n_masses+j) = EQ_IDENTITY           ! polar angle +
             else
                mode(n_masses+j) = EQ_INVERT             ! polar angle -
             end if
          end do
          n_eq(left) = n_eq(left) + 1
          associate (eq_cur => channel(left)%eq(n_eq(left)))
            eq_cur%c = right
            eq_cur%perm = perm
            eq_cur%mode = mode
          end associate
          eq => eq%next
       end do
    end do
  end subroutine phs_forest_get_equivalences

  module subroutine phs_forest_evaluate_selected_channel &
       (forest, channel, active, sqrts, x, phs_factor, volume, ok)
    class(phs_forest_t), intent(inout) :: forest
    integer, intent(in) :: channel
    logical, dimension(:), intent(in) :: active
    real(default), intent(in) :: sqrts
    real(default), dimension(:,:), intent(inout) :: x
    real(default), dimension(:), intent(out) :: phs_factor
    real(default), intent(out) :: volume
    logical, intent(out) :: ok
    integer :: g, t
    integer(TC) :: k, k_root, k_in

    g = forest%grove_lookup (channel)
    t = channel - forest%grove(g)%tree_count_offset
    call forest%prt%set_undefined ()
    call forest%prt_out%set_undefined ()
    k_in = forest%n_tot

    do k = 1,forest%n_in
       forest%prt(ibset(0,k_in-k)) = forest%prt_in(k)
    end do

    do k = 1, forest%n_out
       call forest%prt(ibset(0,k-1))%set_msq &
            (forest%flv(forest%n_in+k)%get_mass () ** 2)
    end do


    k_root = 2**forest%n_out - 1
    select case (forest%n_in)
    case (1)
       forest%prt(k_root) = forest%prt_in(1)
    case (2)
       call forest%prt(k_root)%combine (forest%prt_in(1), forest%prt_in(2))
    end select
    call forest%grove(g)%tree(t)%compute_momenta_from_x (forest%prt,  &
         phs_factor(channel), volume, sqrts, x(:,channel), ok)
    if (ok) then
       do k = 1, forest%n_out
          forest%prt_out(k) = forest%prt(ibset(0,k-1))
       end do
    end if
  end subroutine phs_forest_evaluate_selected_channel

  module subroutine phs_forest_evaluate_other_channels &
       (forest, channel, active, sqrts, x, phs_factor, combine)
    class(phs_forest_t), intent(inout) :: forest
    integer, intent(in) :: channel
    logical, dimension(:), intent(in) :: active
    real(default), intent(in) :: sqrts
    real(default), dimension(:,:), intent(inout) :: x
    real(default), dimension(:), intent(inout) :: phs_factor
    logical, intent(in) :: combine
    integer :: g, t, ch, n_channel

    g = forest%grove_lookup (channel)
    t = channel - forest%grove(g)%tree_count_offset

    n_channel = forest%n_trees
    if (combine) then
       do ch = 1, n_channel
          if (ch == channel)  cycle
          if (active(ch)) then
             g = forest%grove_lookup(ch)
             t = ch - forest%grove(g)%tree_count_offset
             call phs_tree_combine_particles &
                  (forest%grove(g)%tree(t), forest%prt)
          end if
       end do
    end if

    !OMP PARALLEL PRIVATE (g,t,ch) SHARED(active,forest,sqrts,x,channel)
    !OMP DO SCHEDULE(STATIC)
    do ch = 1, n_channel
       if (ch == channel)  cycle
       if (active(ch)) then
          g = forest%grove_lookup(ch)
          t = ch - forest%grove(g)%tree_count_offset
          call forest%grove(g)%tree(t)%compute_x_from_momenta (forest%prt, &
               phs_factor(ch), sqrts, x(:,ch))
       end if
    end do
    !OMP END DO
    !OMP END PARALLEL

  end subroutine phs_forest_evaluate_other_channels

  module subroutine phs_forest_recover_channel &
       (forest, channel, sqrts, x, phs_factor, volume)
    class(phs_forest_t), intent(inout) :: forest
    integer, intent(in) :: channel
    real(default), intent(in) :: sqrts
    real(default), dimension(:,:), intent(inout) :: x
    real(default), dimension(:), intent(inout) :: phs_factor
    real(default), intent(out) :: volume
    integer :: g, t
    integer(TC) :: k, k_in
    g = forest%grove_lookup (channel)
    t = channel - forest%grove(g)%tree_count_offset
    call forest%prt%set_undefined ()
    k_in = forest%n_tot
    forall (k = 1:forest%n_in)
       forest%prt(ibset(0,k_in-k)) = forest%prt_in(k)
    end forall
    forall (k = 1:forest%n_out)
       forest%prt(ibset(0,k-1)) = forest%prt_out(k)
    end forall
    call forest%combine_particles ()
    call forest%grove(g)%tree(t)%compute_volume (sqrts, volume)
    call forest%grove(g)%tree(t)%compute_x_from_momenta (forest%prt, &
         phs_factor(channel), sqrts, x(:,channel))
  end subroutine phs_forest_recover_channel


end submodule phs_forests_s

