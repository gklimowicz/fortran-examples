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

submodule (cascades) cascades_s

  use io_units
  use constants, only: one
  use format_defs, only: FMT_12, FMT_19
  use numeric_utils
  use diagnostics
  use hashes
  use sorting
  use lorentz

  implicit none

contains

  subroutine cascade_init (cascade, depth)
    type(cascade_t), intent(out) :: cascade
    integer, intent(in) :: depth
    integer, save :: index = 0
    index = cascade_index ()
    cascade%index = index
    cascade%depth = depth
    cascade%active = .true.
    allocate (cascade%tree (depth))
    allocate (cascade%tree_pdg (depth))
    allocate (cascade%tree_mapping (depth))
    allocate (cascade%tree_resonant (depth))
  end subroutine cascade_init
  function cascade_index (seed) result (index)
    integer :: index
    integer, intent(in), optional :: seed
    integer, save :: i = 0
    if (present (seed))  i = seed
    i = i + 1
    index = i
  end function cascade_index

  subroutine cascade_write_file_format (cascade, model, unit)
    type(cascade_t), intent(in) :: cascade
    class(model_data_t), intent(in), target :: model
    integer, intent(in), optional :: unit
    type(flavor_t) :: flv
    integer :: u, i
2   format(3x,A,1x,I3,1x,A,1x,I9,1x,'!',1x,A)
    u = given_output_unit (unit);  if (u < 0)  return
    call write_reduced (cascade%tree, u)
    write (u, "(A)")
    do i = 1, cascade%depth
       call flv%init (cascade%tree_pdg(i), model)
       select case (cascade%tree_mapping(i))
       case (NO_MAPPING, EXTERNAL_PRT)
       case (S_CHANNEL)
          write(u,2) 'map', &
               cascade%tree(i), 's_channel', cascade%tree_pdg(i), &
               char (flv%get_name ())
       case (T_CHANNEL)
          write(u,2) 'map', &
               cascade%tree(i), 't_channel', abs (cascade%tree_pdg(i)), &
               char (flv%get_name ())
       case (U_CHANNEL)
          write(u,2) 'map', &
               cascade%tree(i), 'u_channel', abs (cascade%tree_pdg(i)), &
               char (flv%get_name ())
       case (RADIATION)
          write(u,2) 'map', &
               cascade%tree(i), 'radiation', cascade%tree_pdg(i), &
               char (flv%get_name ())
       case (COLLINEAR)
          write(u,2) 'map', &
               cascade%tree(i), 'collinear', cascade%tree_pdg(i), &
               char (flv%get_name ())
       case (INFRARED)
          write(u,2) 'map', &
               cascade%tree(i), 'infrared ', cascade%tree_pdg(i), &
               char (flv%get_name ())
       case (ON_SHELL)
          write(u,2) 'map', &
               cascade%tree(i), 'on_shell ', cascade%tree_pdg(i), &
               char (flv%get_name ())
       case default
          call msg_bug (" Impossible mapping mode encountered")
       end select
    end do
  contains
    subroutine write_reduced (array, unit)
      integer(TC), dimension(:), intent(in) :: array
      integer, intent(in) :: unit
      integer :: i
      write (u, "(3x,A,1x)", advance="no")  "tree"
      do i = 1, size (array)
         if (decay_level (array(i)) > 1) then
            write (u, "(1x,I0)", advance="no")  array(i)
         end if
      end do
    end subroutine write_reduced

    elemental function decay_level (k) result (l)
      integer(TC), intent(in) :: k
      integer :: l
      integer :: i
      l = 0
      do i = 0, bit_size(k) - 1
         if (btest(k,i)) l = l + 1
      end do
    end function decay_level
    subroutine start_comment (u)
      integer, intent(in) :: u
      write(u, '(1x,A)', advance='no') '!'
    end subroutine start_comment
  end subroutine cascade_write_file_format

  subroutine cascade_write_graph_format (cascade, count, unit)
    type(cascade_t), intent(in) :: cascade
    integer, intent(in) :: count
    integer, intent(in), optional :: unit
    integer :: u
    integer(TC) :: mask
    type(string_t) :: left_str, right_str
    u = given_output_unit (unit);  if (u < 0)  return
    mask = 2**((cascade%depth+3)/2) - 1
    left_str = ""
    right_str = ""
    write (u, '(A)') "\begin{minipage}{105pt}"
    write (u, '(A)') "\vspace{30pt}"
    write (u, '(A)') "\begin{center}"
    write (u, '(A)') "\begin{fmfgraph*}(55,55)"
    call graph_write (cascade, mask)
    write (u, '(A)') "\fmfleft{" // char (extract (left_str, 2)) // "}"
    write (u, '(A)') "\fmfright{" // char (extract (right_str, 2)) // "}"
    write (u, '(A)') "\end{fmfgraph*}\\"
    write (u, '(A,I5,A)') "\fbox{$", count, "$}"
    write (u, '(A)') "\end{center}"
    write (u, '(A)') "\end{minipage}"
    write (u, '(A)') "%"
  contains
    recursive subroutine graph_write (cascade, mask, reverse)
      type(cascade_t), intent(in) :: cascade
      integer(TC), intent(in) :: mask
      logical, intent(in), optional :: reverse
      type(flavor_t) :: anti
      logical :: rev
      rev = .false.;  if (present(reverse))  rev = reverse
      if (cascade%has_children) then
         if (.not.rev) then
            call vertex_write (cascade, cascade%daughter1, mask)
            call vertex_write (cascade, cascade%daughter2, mask)
         else
            call vertex_write (cascade, cascade%daughter2, mask, .true.)
            call vertex_write (cascade, cascade%daughter1, mask, .true.)
         end if
         if (cascade%complete) then
            call vertex_write (cascade, cascade%mother, mask, .true.)
            write (u, '(A,I0,A)') "\fmfv{d.shape=square}{v0}"
         end if
      else
         if (cascade%incoming) then
            anti = cascade%flv%anti ()
            call external_write (cascade%bincode, anti%get_tex_name (), &
                 left_str)
         else
            call external_write (cascade%bincode, cascade%flv%get_tex_name (), &
                 right_str)
         end if
      end if
    end subroutine graph_write
    recursive subroutine vertex_write (cascade, daughter, mask, reverse)
      type(cascade_t), intent(in) :: cascade, daughter
      integer(TC), intent(in) :: mask
      logical, intent(in), optional :: reverse
      integer :: bincode
      if (cascade%complete) then
         bincode = 0
      else
         bincode = cascade%bincode
      end if
      call graph_write (daughter, mask, reverse)
      if (daughter%has_children) then
         call line_write (bincode, daughter%bincode, daughter%flv, &
              mapping=daughter%mapping)
      else
         call line_write (bincode, daughter%bincode, daughter%flv)
      end if
    end subroutine vertex_write
    subroutine line_write (i1, i2, flv, mapping)
      integer(TC), intent(in) :: i1, i2
      type(flavor_t), intent(in) :: flv
      integer, intent(in), optional :: mapping
      integer :: k1, k2
      type(string_t) :: prt_type
      select case (flv%get_spin_type ())
      case (SCALAR);       prt_type = "plain"
      case (SPINOR);       prt_type = "fermion"
      case (VECTOR);       prt_type = "boson"
      case (VECTORSPINOR); prt_type = "fermion"
      case (TENSOR);       prt_type = "dbl_wiggly"
      case default;        prt_type = "dashes"
      end select
      if (flv%is_antiparticle ()) then
         k1 = i2;  k2 = i1
      else
         k1 = i1;  k2 = i2
      end if
      if (present (mapping)) then
         select case (mapping)
         case (S_CHANNEL)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=blue,lab=\sm\blue$" // &
                 & char (flv%get_tex_name ()) // "$}" // &
                 & "{v", k1, ",v", k2, "}"
         case (T_CHANNEL, U_CHANNEL)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=cyan,lab=\sm\cyan$" // &
                 & char (flv%get_tex_name ()) // "$}" // &
                 & "{v", k1, ",v", k2, "}"
         case (RADIATION)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=green,lab=\sm\green$" // &
                 & char (flv%get_tex_name ()) // "$}" // &
                 & "{v", k1, ",v", k2, "}"
         case (COLLINEAR)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=magenta,lab=\sm\magenta$" // &
                 & char (flv%get_tex_name ()) // "$}" // &
                 & "{v", k1, ",v", k2, "}"
         case (INFRARED)
            write (u, '(A,I0,A,I0,A)') "\fmf{" // char (prt_type) // &
                 & ",f=red,lab=\sm\red$" // &
                 & char (flv%get_tex_name ()) // "$}" // &
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
  end subroutine cascade_write_graph_format

  subroutine cascade_write (cascade, unit)
    type(cascade_t), intent(in) :: cascade
    integer, intent(in), optional :: unit
    integer :: u
    character(9) :: depth
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A,(1x,I7))") 'Cascade #', cascade%index
    write (u, "(A,(1x,I7))") '  Grove:       #', cascade%grove
    write (u, "(A,3(1x,L1))") '  act/cmp/inc:  ', &
         cascade%active, cascade%complete, cascade%incoming
    write (u, "(A,I0)") '  Bincode:      ', cascade%bincode
    write (u, "(A)", advance="no") '  Flavor:       '
    call cascade%flv%write (unit)
    write (u, "(A,I9)") '  Active flavor:', cascade%pdg
    write (u, "(A,L1)") '  Is vector:    ', cascade%is_vector
    write (u, "(A,3(1x," // FMT_19 // "))") '  Mass (m/r/e): ', &
         cascade%m_min, cascade%m_rea, cascade%m_eff
    write (u, "(A,I1)") '  Mapping:      ', cascade%mapping
    write (u, "(A,3(1x,L1))") '  res/log/tch:  ', &
         cascade%resonant, cascade%log_enhanced, cascade%t_channel
    write (u, "(A,(1x,I7))") '  Multiplicity: ', cascade%multiplicity
    write (u, "(A,2(1x,I7))") '  n intern/off: ', &
         cascade%internal, cascade%n_off_shell
    write (u, "(A,3(1x,I7))") '  n res/log/tch:', &
         cascade%n_resonances, cascade%n_log_enhanced, cascade%n_t_channel
    write (u, "(A,I7)") '  Depth:        ', cascade%depth
    write (depth, "(I7)") cascade%depth
    write (u, "(A," // depth // "(1x,I7))") &
       '  Tree:         ', cascade%tree
    write (u, "(A," // depth // "(1x,I7))") &
       '  Tree(PDG):    ', cascade%tree_pdg
    write (u, "(A," // depth // "(1x,I7))") &
       '  Tree(mapping):', cascade%tree_mapping
    write (u, "(A," // depth // "(1x,L1))") &
       '  Tree(res):    ', cascade%tree_resonant
    if (cascade%has_children) then
       write (u, "(A,I7,1x,I7)") '  Daughter1/2:  ', &
            cascade%daughter1%index, cascade%daughter2%index
    end if
    if (associated (cascade%mother)) then
       write (u, "(A,I7)") '  Mother:       ', cascade%mother%index
    end if
  end subroutine cascade_write

  subroutine cascade_init_outgoing (cascade, flv, pos, m_thr)
    type(cascade_t), intent(out) :: cascade
    type(flavor_t), intent(in) :: flv
    integer, intent(in) :: pos
    real(default), intent(in) :: m_thr
    call cascade_init (cascade, 1)
    cascade%bincode = ibset (0_TC, pos-1)
    cascade%flv = flv
    cascade%pdg = cascade%flv%get_pdg ()
    cascade%is_vector = flv%get_spin_type () == VECTOR
    cascade%m_min = flv%get_mass ()
    cascade%m_rea = cascade%m_min
    if (cascade%m_rea >= m_thr) then
       cascade%m_eff = cascade%m_rea
    end if
    cascade%on_shell = .true.
    cascade%multiplicity = 1
    cascade%tree(1) = cascade%bincode
    cascade%tree_pdg(1) = cascade%pdg
    cascade%tree_mapping(1) = EXTERNAL_PRT
    cascade%tree_resonant(1) = .false.
  end subroutine cascade_init_outgoing

  subroutine cascade_init_incoming (cascade, flv, pos, m_thr)
    type(cascade_t), intent(out) :: cascade
    type(flavor_t), intent(in) :: flv
    integer, intent(in) :: pos
    real(default), intent(in) :: m_thr
    call cascade_init (cascade, 1)
    cascade%incoming = .true.
    cascade%bincode = ibset (0_TC, pos-1)
    cascade%flv = flv%anti ()
    cascade%pdg = cascade%flv%get_pdg ()
    cascade%is_vector = flv%get_spin_type () == VECTOR
    cascade%m_min = flv%get_mass ()
    cascade%m_rea = cascade%m_min
    if (cascade%m_rea >= m_thr) then
       cascade%m_eff = cascade%m_rea
    end if
    cascade%on_shell = .true.
    cascade%n_t_channel = 0
    cascade%n_off_shell = 0
    cascade%tree(1) = cascade%bincode
    cascade%tree_pdg(1) = cascade%pdg
    cascade%tree_mapping(1) = EXTERNAL_PRT
    cascade%tree_resonant(1) = .false.
  end subroutine cascade_init_incoming

  module function cascade_disjunct (cascade1, cascade2) result (flag)
    logical :: flag
    type(cascade_t), intent(in) :: cascade1, cascade2
    flag = iand (cascade1%bincode, cascade2%bincode) == 0
  end function cascade_disjunct

  subroutine cascade_assign_resonance_hash (cascade)
    type(cascade_t), intent(inout) :: cascade
    integer(i8), dimension(1) :: mold
    cascade%res_hash = hash (transfer &
         ([sort (pack (cascade%tree_pdg, &
                 cascade%tree_resonant)), &
           sort (pack (abs (cascade%tree_pdg), &
                 cascade%tree_mapping == T_CHANNEL .or. &
                 cascade%tree_mapping == U_CHANNEL))], &
          mold))
  end subroutine cascade_assign_resonance_hash

  module subroutine hash_entry_init (entry, entry_in)
    type(hash_entry_t), intent(out) :: entry
    type(hash_entry_t), intent(in) :: entry_in
    type(cascade_p), pointer :: casc_iter, casc_copy
    entry%hashval = entry_in%hashval
    entry%key = entry_in%key
    casc_iter => entry_in%first
    do while (associated (casc_iter))
       allocate (casc_copy)
       casc_copy = casc_iter
       casc_copy%next => null ()
       if (associated (entry%first)) then
          entry%last%next => casc_copy
       else
          entry%first => casc_copy
       end if
       entry%last => casc_copy
       casc_iter => casc_iter%next
    end do
  end subroutine hash_entry_init

  subroutine hash_entry_final (hash_entry)
    type(hash_entry_t), intent(inout) :: hash_entry
    type(cascade_p), pointer :: current
    do while (associated (hash_entry%first))
       current => hash_entry%first
       hash_entry%first => current%next
       deallocate (current)
    end do
  end subroutine hash_entry_final

  subroutine hash_entry_write (hash_entry, unit)
    type(hash_entry_t), intent(in) :: hash_entry
    integer, intent(in), optional :: unit
    type(cascade_p), pointer :: current
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)", advance="no")  "Entry:"
    do i = 1, size (hash_entry%key)
       write (u, "(1x,I0)", advance="no")  hash_entry%key(i)
    end do
    write (u, "(1x,A)", advance="no")  "->"
    current => hash_entry%first
    do while (associated (current))
       write (u, "(1x,I7)", advance="no") current%cascade%index
       current => current%next
    end do
    write (u, *)
  end subroutine hash_entry_write

  subroutine hash_entry_add_cascade_ptr (hash_entry, cascade, ok, cascade_ptr)
    type(hash_entry_t), intent(inout) :: hash_entry
    type(cascade_t), intent(in), target :: cascade
    logical, intent(out), optional :: ok
    type(cascade_t), optional, pointer :: cascade_ptr
    type(cascade_p), pointer :: current
    if (present (ok)) then
       call hash_entry_check_cascade (hash_entry, cascade, ok, cascade_ptr)
       if (.not. ok)  return
    end if
    allocate (current)
    current%cascade => cascade
    if (associated (hash_entry%last)) then
       hash_entry%last%next => current
    else
       hash_entry%first => current
    end if
    hash_entry%last => current
  end subroutine hash_entry_add_cascade_ptr

  subroutine hash_entry_check_cascade (hash_entry, cascade, ok, cascade_ptr)
    type(hash_entry_t), intent(in), target :: hash_entry
    type(cascade_t), intent(in), target :: cascade
    logical, intent(out) :: ok
    type(cascade_t), optional, pointer :: cascade_ptr
    type(cascade_p), pointer :: current
    integer, dimension(:), allocatable :: tree_pdg
    ok = .true.
    allocate (tree_pdg (size (cascade%tree_pdg)))
    if (cascade%complete) then
       where (cascade%tree_mapping == INFRARED .or. &
            cascade%tree_mapping == COLLINEAR .or. &
            cascade%tree_mapping == T_CHANNEL .or. &
            cascade%tree_mapping == U_CHANNEL)
          tree_pdg = 0
       elsewhere
          tree_pdg = cascade%tree_pdg
       end where
    else
       tree_pdg = cascade%tree_pdg
    end if
    current => hash_entry%first
    do while (associated (current))
       if (current%cascade%depth == cascade%depth) then
          if (all (current%cascade%tree == cascade%tree)) then
             if (all (current%cascade%tree_mapping == cascade%tree_mapping)) &
                  then
                if (all (current%cascade%tree_pdg .match. tree_pdg)) then
                   if (present (cascade_ptr))  cascade_ptr => current%cascade
                   ok = .false.;  return
                end if
             end if
          end if
       end if
       current => current%next
    end do
    if (present (cascade_ptr))  cascade_ptr => cascade
  end subroutine hash_entry_check_cascade

  elemental module function pdg_match (pdg1, pdg2) result (flag)
    logical :: flag
    integer(TC), intent(in) :: pdg1, pdg2
    select case (pdg1)
    case (0)
       flag = .true.
    case default
       select case (pdg2)
       case (0)
          flag = .true.
       case default
          flag = pdg1 == pdg2
       end select
    end select
  end function pdg_match

  module subroutine cascade_set_init_from_cascade &
       (cascade_set, cascade_set_in)
    type(cascade_set_t), intent(out) :: cascade_set
    type(cascade_set_t), intent(in), target :: cascade_set_in
    type(cascade_t), pointer :: casc_iter, casc_copy
    cascade_set%model => cascade_set_in%model
    cascade_set%n_in = cascade_set_in%n_in
    cascade_set%n_out = cascade_set_in%n_out
    cascade_set%n_tot = cascade_set_in%n_tot
    cascade_set%flv = cascade_set_in%flv
    cascade_set%depth_out = cascade_set_in%depth_out
    cascade_set%depth_tot = cascade_set_in%depth_tot
    cascade_set%sqrts = cascade_set_in%sqrts
    cascade_set%m_threshold_s = cascade_set_in%m_threshold_s
    cascade_set%m_threshold_t = cascade_set_in%m_threshold_t
    cascade_set%off_shell = cascade_set_in%off_shell
    cascade_set%t_channel = cascade_set_in%t_channel
    cascade_set%keep_nonresonant = cascade_set_in%keep_nonresonant
    cascade_set%n_groves = cascade_set_in%n_groves

    casc_iter => cascade_set_in%first
    do while (associated (casc_iter))
       allocate (casc_copy)
       casc_copy = casc_iter
       casc_copy%next => null ()
       if (associated (cascade_set%first)) then
          cascade_set%last%next => casc_copy
       else
          cascade_set%first => casc_copy
       end if
       cascade_set%last => casc_copy
       casc_iter => casc_iter%next
    end do

    cascade_set%n_entries = cascade_set_in%n_entries
    cascade_set%fill_ratio = cascade_set_in%fill_ratio
    cascade_set%n_entries_max = cascade_set_in%n_entries_max
    cascade_set%mask = cascade_set_in%mask
    cascade_set%fatal_beam_decay = cascade_set_in%fatal_beam_decay
    allocate (cascade_set%entry (0:cascade_set%mask))
    cascade_set%entry = cascade_set_in%entry
  end subroutine cascade_set_init_from_cascade

  module function cascade_set_is_valid (cascade_set) result (flag)
    logical :: flag
    type(cascade_set_t), intent(in) :: cascade_set
    type(cascade_t), pointer :: cascade
    flag = .false.
    cascade => cascade_set%first_k
    do while (associated (cascade))
       if (cascade%active .and. cascade%complete) then
          flag = .true.
          return
       end if
       cascade => cascade%next
    end do
  end function cascade_set_is_valid

  module subroutine cascade_set_init_base (cascade_set, model, &
       n_in, n_out, phs_par, fatal_beam_decay, flv)
    type(cascade_set_t), intent(out) :: cascade_set
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: n_in, n_out
    type(phs_parameters_t), intent(in) :: phs_par
    logical, intent(in) :: fatal_beam_decay
    type(flavor_t), dimension(:,:), intent(in), optional :: flv
    integer :: size_guess
    integer :: i, j
    cascade_set%model => model
    cascade_set%n_in = n_in
    cascade_set%n_out = n_out
    cascade_set%n_tot = n_in + n_out
    if (present (flv)) then
       allocate (cascade_set%flv (size (flv, 1), size (flv, 2)))
       do i = 1, size (flv, 2)
          do j = 1, size (flv, 1)
             call cascade_set%flv(j,i)%init (flv(j,i)%get_pdg (), model)
          end do
       end do
    end if
    select case (n_in)
    case (1);  cascade_set%depth_out = 2 * n_out - 3
    case (2);  cascade_set%depth_out = 2 * n_out - 1
    end select
    cascade_set%depth_tot = 2 * cascade_set%n_tot - 3
    cascade_set%sqrts = phs_par%sqrts
    cascade_set%m_threshold_s = phs_par%m_threshold_s
    cascade_set%m_threshold_t = phs_par%m_threshold_t
    cascade_set%off_shell = phs_par%off_shell
    cascade_set%t_channel = phs_par%t_channel
    cascade_set%keep_nonresonant = phs_par%keep_nonresonant
    cascade_set%fill_ratio = CASCADE_SET_FILL_RATIO
    size_guess = ishft (256, min (2 * (cascade_set%n_tot - 3), 22))
    cascade_set%n_entries_max = size_guess * cascade_set%fill_ratio
    cascade_set%mask = size_guess - 1
    allocate (cascade_set%entry (0:cascade_set%mask))
    cascade_set%fatal_beam_decay = fatal_beam_decay
  end subroutine cascade_set_init_base

  module subroutine cascade_set_final (cascade_set)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), pointer :: current
    integer :: i
    if (allocated (cascade_set%entry)) then
       do i = 0, cascade_set%mask
          call hash_entry_final (cascade_set%entry(i))
       end do
       deallocate (cascade_set%entry)
    end if
    do while (associated (cascade_set%first))
       current => cascade_set%first
       cascade_set%first => cascade_set%first%next
       deallocate (current)
    end do
  end subroutine cascade_set_final

  module subroutine cascade_set_write_process_bincode_format &
       (cascade_set, unit)
    type(cascade_set_t), intent(in), target :: cascade_set
    integer, intent(in), optional :: unit
    integer, dimension(:), allocatable :: bincode, field_width
    integer :: n_in, n_out, n_tot, n_flv
    integer :: u, f, i, bc
    character(20) :: str
    type(string_t) :: fmt_head
    type(string_t), dimension(:), allocatable :: fmt_proc
    u = given_output_unit (unit);  if (u < 0)  return
    if (.not. allocated (cascade_set%flv)) return
    write (u, "('!',1x,A)")  "List of subprocesses with particle bincodes:"
    n_in  = cascade_set%n_in
    n_out = cascade_set%n_out
    n_tot = cascade_set%n_tot
    n_flv = size (cascade_set%flv, 2)
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
               len (cascade_set%flv(i,f)%get_name ()))
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
               char (cascade_set%flv(i,f)%get_name ())
          if (i == n_in)  write (u, "(1x,'=>')", advance="no")
       end do
       write (u, *)
    end do
    write (u, char (fmt_head))  bincode
  end subroutine cascade_set_write_process_bincode_format

  subroutine cascade_set_write_process_tex_format (cascade_set, unit)
    type(cascade_set_t), intent(in), target :: cascade_set
    integer, intent(in), optional :: unit
    integer :: u, f, i
    u = given_output_unit (unit);  if (u < 0)  return
    if (.not. allocated (cascade_set%flv)) return
    write (u, "(A)")  "\begin{align*}"
    do f = 1, size (cascade_set%flv, 2)
       do i = 1, cascade_set%n_in
          if (i > 1)  write (u, "(A)", advance="no") "\quad "
          write (u, "(A)", advance="no") &
               char (cascade_set%flv(i,f)%get_tex_name ())
       end do
       write (u, "(A)", advance="no")  "\quad &\to\quad "
       do i = cascade_set%n_in + 1, cascade_set%n_tot
          if (i > cascade_set%n_in + 1)  write (u, "(A)", advance="no") "\quad "
          write (u, "(A)", advance="no") &
               char (cascade_set%flv(i,f)%get_tex_name ())
       end do
       if (f < size (cascade_set%flv, 2)) then
          write (u, "(A)")  "\\"
       else
          write (u, "(A)")  ""
       end if
    end do
    write (u, "(A)")  "\end{align*}"
  end subroutine cascade_set_write_process_tex_format

  module subroutine cascade_set_write_file_format (cascade_set, unit)
    type(cascade_set_t), intent(in), target :: cascade_set
    integer, intent(in), optional :: unit
    type(cascade_t), pointer :: cascade
    integer :: u, grove, count
    logical :: first_in_grove
    u = given_output_unit (unit);  if (u < 0)  return
    count = 0
    do grove = 1, cascade_set%n_groves
       first_in_grove = .true.
       cascade => cascade_set%first_k
       do while (associated (cascade))
          if (cascade%active .and. cascade%complete) then
             if (cascade%grove == grove) then
                if (first_in_grove) then
                   first_in_grove = .false.
                   write (u, "(A)")
                   write (u, "(1x,'!',1x,A,1x,I0,A)", advance='no') &
                      'Multiplicity =', cascade%multiplicity, ","
                   select case (cascade%n_resonances)
                   case (0)
                      write (u, '(1x,A)', advance='no') 'no resonances, '
                   case (1)
                      write (u, '(1x,A)', advance='no') '1 resonance,  '
                   case default
                      write (u, '(1x,I0,1x,A)', advance='no') &
                           cascade%n_resonances, 'resonances, '
                   end select
                   write (u, '(1x,I0,1x,A)', advance='no') &
                        cascade%n_log_enhanced, 'logs, '
                   write (u, '(1x,I0,1x,A)', advance='no') &
                        cascade%n_off_shell, 'off-shell, '
                   select case (cascade%n_t_channel)
                   case (0);  write (u, '(1x,A)') 's-channel graph'
                   case (1);  write (u, '(1x,A)') '1 t-channel line'
                   case default
                      write(u,'(1x,I0,1x,A)') &
                           cascade%n_t_channel, 't-channel lines'
                   end select
                   write (u, '(1x,A,I0)') 'grove #', grove
                end if
                count = count + 1
                write (u, "(1x,'!',1x,A,I0)")  "Channel #", count
                call cascade_write_file_format (cascade, cascade_set%model, u)
             end if
          end if
          cascade => cascade%next
       end do
    end do
  end subroutine cascade_set_write_file_format

  module subroutine cascade_set_write_graph_format &
      (cascade_set, filename, process_id, unit)
    type(cascade_set_t), intent(in), target :: cascade_set
    type(string_t), intent(in) :: filename, process_id
    integer, intent(in), optional :: unit
    type(cascade_t), pointer :: cascade
    integer :: u, grove, count, pgcount
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
    call cascade_set_write_process_tex_format (cascade_set, u)
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
    do grove = 1, cascade_set%n_groves
       first_in_grove = .true.
       cascade => cascade_set%first
       do while (associated (cascade))
          if (cascade%active .and. cascade%complete) then
             if (cascade%grove == grove) then
                if (first_in_grove) then
                   first_in_grove = .false.
                   write (u, *)
                   write (u, '(A)') "\vspace{20pt}"
                   write (u, '(A)') "\begin{tabular}{l}"
                   write (u, '(A,I5,A)') &
                        & "\fbox{\bf Grove \boldmath$", grove, "$} \\[10pt]"
                   write (u, '(A,I1,A)') "Multiplicity: ", &
                        cascade%multiplicity, "\\"
                   write (u, '(A,I1,A)') "Resonances:   ", &
                        cascade%n_resonances, "\\"
                   write (u, '(A,I1,A)') "Log-enhanced: ", &
                        cascade%n_log_enhanced, "\\"
                   write (u, '(A,I1,A)') "Off-shell:    ", &
                        cascade%n_off_shell, "\\"
                   write (u, '(A,I1,A)') "t-channel:    ", &
                        cascade%n_t_channel, ""
                   write (u, '(A)') "\end{tabular}"
                end if
                count = count + 1
                call cascade_write_graph_format (cascade, count, unit)
                if (pgcount >= 250) then
                   write (u, '(A)') "\clearpage"
                   pgcount = 0
                end if
             end if
          end if
          cascade => cascade%next
       end do
    end do
    write (u, '(A)') "\end{flushleft}"
    write (u, '(A)') "\end{fmfshrink}"
    write (u, '(A)') "\end{fmffile}"
    write (u, '(A)') "\end{document}"
 end subroutine cascade_set_write_graph_format

  module subroutine cascade_set_write &
       (cascade_set, unit, active_only, complete_only)
    type(cascade_set_t), intent(in), target :: cascade_set
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: active_only, complete_only
    logical :: active, complete
    type(cascade_t), pointer :: cascade
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    active = .true.;  if (present (active_only))  active = active_only
    complete = .false.;  if (present (complete_only))  complete = complete_only
    write (u, "(A)") "Cascade set:"
    write (u, "(3x,A)", advance="no")  "Model:"
    if (associated (cascade_set%model)) then
       write (u, "(1x,A)") char (cascade_set%model%get_name ())
    else
       write (u, "(1x,A)") "[none]"
    end if
    write (u, "(3x,A)", advance="no")  "n_in/out/tot  ="
    write (u, "(3(1x,I7))")  &
         cascade_set%n_in, cascade_set%n_out, cascade_set%n_tot
    write (u, "(3x,A)", advance="no")  "depth_out/tot ="
    write (u, "(2(1x,I7))")  cascade_set%depth_out, cascade_set%depth_tot
    write (u, "(3x,A)", advance="no")  "mass thr(s/t) ="
    write (u, "(2(1x," // FMT_19 // "))")  &
         cascade_set%m_threshold_s, cascade_set%m_threshold_t
    write (u, "(3x,A)", advance="no")  "off shell     ="
    write (u, "(1x,I7)")  cascade_set%off_shell
    write (u, "(3x,A)", advance="no")  "keep_nonreson ="
    write (u, "(1x,L1)")  cascade_set%keep_nonresonant
    write (u, "(3x,A)", advance="no")  "n_groves      ="
    write (u, "(1x,I7)")  cascade_set%n_groves
    write (u, "(A)")
    write (u, "(A)") "Cascade list:"
    if (associated (cascade_set%first)) then
       cascade => cascade_set%first
       do while (associated (cascade))
          if (active .and. .not. cascade%active)  cycle
          if (complete .and. .not. cascade%complete)  cycle
          call cascade_write (cascade, unit)
          cascade => cascade%next
       end do
    else
       write (u, "(A)") "[empty]"
    end if
    write (u, "(A)") "Hash array"
    write (u, "(3x,A)", advance="no")  "n_entries     ="
    write (u, "(1x,I7)")  cascade_set%n_entries
    write (u, "(3x,A)", advance="no")  "fill_ratio    ="
    write (u, "(1x," // FMT_12 // ")")  cascade_set%fill_ratio
    write (u, "(3x,A)", advance="no")  "n_entries_max ="
    write (u, "(1x,I7)")  cascade_set%n_entries_max
    write (u, "(3x,A)", advance="no")  "mask          ="
    write (u, "(1x,I0)")  cascade_set%mask
    do i = 0, ubound (cascade_set%entry, 1)
       if (allocated (cascade_set%entry(i)%key)) then
          write (u, "(1x,I7)") i
          call hash_entry_write (cascade_set%entry(i), u)
       end if
    end do
  end subroutine cascade_set_write

  recursive subroutine cascade_set_add_copy &
       (cascade_set, cascade_in, cascade_ptr)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), intent(in) :: cascade_in
    type(cascade_t), optional, pointer :: cascade_ptr
    type(cascade_t), pointer :: cascade
    logical :: ok
    allocate (cascade)
    cascade = cascade_in
    if (associated (cascade_in%daughter1))  call cascade_set_add_copy &
         (cascade_set, cascade_in%daughter1, cascade%daughter1)
    if (associated (cascade_in%daughter2))  call cascade_set_add_copy &
         (cascade_set, cascade_in%daughter2, cascade%daughter2)
    if (associated (cascade_in%mother))  call cascade_set_add_copy &
         (cascade_set, cascade_in%mother, cascade%mother)
    cascade%next => null ()
    call cascade_set_add (cascade_set, cascade, ok, cascade_ptr)
    if (.not. ok)  deallocate (cascade)
  end subroutine cascade_set_add_copy

  subroutine cascade_set_add (cascade_set, cascade, ok, cascade_ptr)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), intent(in), target :: cascade
    logical, intent(out) :: ok
    type(cascade_t), optional, pointer :: cascade_ptr
    integer(i8), dimension(1) :: mold
    call cascade_set_hash_insert &
         (cascade_set, transfer (cascade%tree, mold), cascade, ok, cascade_ptr)
    if (ok)  call cascade_set_list_add (cascade_set, cascade)
  end subroutine cascade_set_add

  subroutine cascade_set_list_add (cascade_set, cascade)
    type(cascade_set_t), intent(inout) :: cascade_set
    type(cascade_t), intent(in), target :: cascade
    if (associated (cascade_set%last)) then
       cascade_set%last%next => cascade
    else
       cascade_set%first => cascade
    end if
    cascade_set%last => cascade
  end subroutine cascade_set_list_add

  subroutine cascade_set_hash_insert &
       (cascade_set, key, cascade, ok, cascade_ptr)
    type(cascade_set_t), intent(inout), target :: cascade_set
    integer(i8), dimension(:), intent(in) :: key
    type(cascade_t), intent(in), target :: cascade
    logical, intent(out) :: ok
    type(cascade_t), optional, pointer :: cascade_ptr
    integer(i32) :: h
    if (cascade_set%n_entries >= cascade_set%n_entries_max) &
         call cascade_set_hash_expand (cascade_set)
    h = hash (key)
    call cascade_set_hash_insert_rec &
         (cascade_set, h, h, key, cascade, ok, cascade_ptr)
  end subroutine cascade_set_hash_insert

  subroutine cascade_set_hash_expand (cascade_set)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(hash_entry_t), dimension(:), allocatable, target :: table_tmp
    type(cascade_p), pointer :: current
    integer :: i, s
    allocate (table_tmp (0:cascade_set%mask))
    table_tmp = cascade_set%entry
    deallocate (cascade_set%entry)
    s = 2 * size (table_tmp)
    cascade_set%n_entries = 0
    cascade_set%n_entries_max = s * cascade_set%fill_ratio
    cascade_set%mask = s - 1
    allocate (cascade_set%entry (0:cascade_set%mask))
    do i = 0, ubound (table_tmp, 1)
       current => table_tmp(i)%first
       do while (associated (current))
          call cascade_set_hash_insert_rec &
               (cascade_set, table_tmp(i)%hashval, table_tmp(i)%hashval, &
                table_tmp(i)%key, current%cascade)
          current => current%next
       end do
    end do
  end subroutine cascade_set_hash_expand

  recursive subroutine cascade_set_hash_insert_rec &
       (cascade_set, h, hashval, key, cascade, ok, cascade_ptr)
    type(cascade_set_t), intent(inout) :: cascade_set
    integer(i32), intent(in) :: h, hashval
    integer(i8), dimension(:), intent(in) :: key
    type(cascade_t), intent(in), target :: cascade
    logical, intent(out), optional :: ok
    type(cascade_t), optional, pointer :: cascade_ptr
    integer(i32) :: i
    i = iand (h, cascade_set%mask)
    if (allocated (cascade_set%entry(i)%key)) then
       if (size (cascade_set%entry(i)%key) /= size (key)) then
          call cascade_set_hash_insert_rec &
               (cascade_set, h + 1, hashval, key, cascade, ok, cascade_ptr)
       else if (any (cascade_set%entry(i)%key /= key)) then
          call cascade_set_hash_insert_rec &
               (cascade_set, h + 1, hashval, key, cascade, ok, cascade_ptr)
       else
          call hash_entry_add_cascade_ptr &
               (cascade_set%entry(i), cascade, ok, cascade_ptr)
       end if
    else
       cascade_set%entry(i)%hashval = hashval
       allocate (cascade_set%entry(i)%key (size (key)))
       cascade_set%entry(i)%key = key
       call hash_entry_add_cascade_ptr &
            (cascade_set%entry(i), cascade, ok, cascade_ptr)
       cascade_set%n_entries = cascade_set%n_entries + 1
    end if
  end subroutine cascade_set_hash_insert_rec

  module subroutine cascade_set_add_outgoing2 (cascade_set, flv)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(flavor_t), dimension(:,:), intent(in) :: flv
    integer :: pos, prc, n_out, n_prc
    type(cascade_t), pointer :: cascade
    logical :: ok
    n_out = size (flv, dim=1)
    n_prc = size (flv, dim=2)
    do prc = 1, n_prc
       do pos = 1, n_out
          allocate (cascade)
          call cascade_init_outgoing &
               (cascade, flv(pos,prc), pos, cascade_set%m_threshold_s)
          call cascade_set_add (cascade_set, cascade, ok)
          if (.not. ok) then
             deallocate (cascade)
          end if
       end do
    end do
  end subroutine cascade_set_add_outgoing2

  module subroutine cascade_set_add_outgoing1 (cascade_set, flv)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(flavor_t), dimension(:), intent(in) :: flv
    integer :: pos, n_out
    type(cascade_t), pointer :: cascade
    logical :: ok
    n_out = size (flv, dim=1)
    do pos = 1, n_out
       allocate (cascade)
       call cascade_init_outgoing &
            (cascade, flv(pos), pos, cascade_set%m_threshold_s)
       call cascade_set_add (cascade_set, cascade, ok)
       if (.not. ok) then
          deallocate (cascade)
       end if
    end do
  end subroutine cascade_set_add_outgoing1

  module subroutine cascade_set_add_incoming1 (cascade_set, n1, n2, pos, flv)
    type(cascade_set_t), intent(inout), target :: cascade_set
    integer, intent(out) :: n1, n2
    integer, intent(in) :: pos
    type(flavor_t), dimension(:), intent(in) :: flv
    integer :: prc, n_prc
    type(cascade_t), pointer :: cascade
    logical :: ok
    n1 = 0
    n2 = 0
    n_prc = size (flv)
    do prc = 1, n_prc
       allocate (cascade)
       call cascade_init_incoming &
            (cascade, flv(prc), pos, cascade_set%m_threshold_t)
       call cascade_set_add (cascade_set, cascade, ok)
       if (ok) then
          if (n1 == 0)  n1 = cascade%index
          n2 = cascade%index
          if (.not. associated (cascade_set%first_t)) then
             cascade_set%first_t => cascade
          end if
       else
          deallocate (cascade)
       end if
    end do
  end subroutine cascade_set_add_incoming1

  module subroutine cascade_set_add_incoming0 (cascade_set, n1, n2, pos, flv)
    type(cascade_set_t), intent(inout), target :: cascade_set
    integer, intent(out) :: n1, n2
    integer, intent(in) :: pos
    type(flavor_t), intent(in) :: flv
    type(cascade_t), pointer :: cascade
    logical :: ok
    n1 = 0
    n2 = 0
    allocate (cascade)
    call cascade_init_incoming &
         (cascade, flv, pos, cascade_set%m_threshold_t)
    call cascade_set_add (cascade_set, cascade, ok)
    if (ok) then
       if (n1 == 0)  n1 = cascade%index
       n2 = cascade%index
       if (.not. associated (cascade_set%first_t)) then
          cascade_set%first_t => cascade
       end if
    else
       deallocate (cascade)
    end if
  end subroutine cascade_set_add_incoming0

  subroutine cascade_match_pair (cascade_set, cascade1, cascade2, s_channel)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), intent(in), target :: cascade1, cascade2
    logical, intent(in) :: s_channel
    integer, dimension(:), allocatable :: pdg3
    integer :: i, depth_max
    type(flavor_t) :: flv
    if (s_channel) then
       depth_max = cascade_set%depth_out
    else
       depth_max = cascade_set%depth_tot
    end if
    if (cascade1%depth + cascade2%depth < depth_max) then
       call cascade_set%model%match_vertex ( &
            cascade1%flv%get_pdg (), &
            cascade2%flv%get_pdg (), &
            pdg3)
       do i = 1, size (pdg3)
          call flv%init (pdg3(i), cascade_set%model)
          if (s_channel) then
             call cascade_combine_s (cascade_set, cascade1, cascade2, flv)
          else
             call cascade_combine_t (cascade_set, cascade1, cascade2, flv)
          end if
       end do
       deallocate (pdg3)
    end if
  end subroutine cascade_match_pair

  subroutine cascade_match_triplet &
       (cascade_set, cascade1, cascade2, cascade3, s_channel)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), intent(in), target :: cascade1, cascade2, cascade3
    logical, intent(in) :: s_channel
    integer :: depth_max
    depth_max = cascade_set%depth_tot
    if (cascade1%depth + cascade2%depth + cascade3%depth == depth_max) then
       if (cascade_set%model%check_vertex ( &
            cascade1%flv%get_pdg (), &
            cascade2%flv%get_pdg (), &
            cascade3%flv%get_pdg ())) then
          call cascade_combine_keystone &
               (cascade_set, cascade1, cascade2, cascade3, s_channel)
       end if
    end if
  end subroutine cascade_match_triplet

  subroutine cascade_combine_s (cascade_set, cascade1, cascade2, flv)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), intent(in), target :: cascade1, cascade2
    type(flavor_t), intent(in) :: flv
    type(cascade_t), pointer :: cascade3, cascade4
    logical :: keep
    keep = .false.
    allocate (cascade3)
    call cascade_init (cascade3, cascade1%depth + cascade2%depth + 1)
    cascade3%bincode = ior (cascade1%bincode, cascade2%bincode)
    cascade3%flv = flv%anti ()
    cascade3%pdg = cascade3%flv%get_pdg ()
    cascade3%is_vector = flv%get_spin_type () == VECTOR
    cascade3%m_min = cascade1%m_min + cascade2%m_min
    cascade3%m_rea = flv%get_mass ()
    if (cascade3%m_rea > cascade_set%m_threshold_s) then
       cascade3%m_eff = cascade3%m_rea
    end if
    ! Potentially resonant cases [sqrts = m_rea for on-shell decay]
    if (cascade3%m_rea > cascade3%m_min &
         .and. cascade3%m_rea <= cascade_set%sqrts) then
       if (flv%get_width () /= 0) then
          if (cascade1%on_shell .or. cascade2%on_shell) then
             keep = .true.
             cascade3%mapping = S_CHANNEL
             cascade3%resonant = .true.
          end if
       else
          call warn_decay (flv)
       end if
    ! Collinear and IR singular cases
    else if (cascade3%m_rea < cascade_set%sqrts) then
       ! Massless splitting
       if (cascade1%m_eff == 0 .and. cascade2%m_eff == 0 &
            .and. cascade3%depth <= 3) then
          keep = .true.
          cascade3%log_enhanced = .true.
          if (cascade3%is_vector) then
             if (cascade1%is_vector .and. cascade2%is_vector) then
                cascade3%mapping = COLLINEAR   ! three-vector-vertex
             else
                cascade3%mapping = INFRARED    ! vector splitting into matter
             end if
          else
             if (cascade1%is_vector .or. cascade2%is_vector) then
                cascade3%mapping = COLLINEAR   ! vector radiation off matter
             else
                cascade3%mapping = INFRARED    ! scalar radiation/splitting
             end if
          end if
       ! IR radiation off massive particle
       else if (cascade3%m_eff > 0 .and. cascade1%m_eff > 0 &
            .and. cascade2%m_eff == 0 &
            .and. (cascade1%on_shell .or. cascade1%mapping == RADIATION) &
            .and. abs (cascade3%m_eff - cascade1%m_eff) &
                       < cascade_set%m_threshold_s) &
            then
          keep = .true.
          cascade3%log_enhanced = .true.
          cascade3%mapping = RADIATION
       else if (cascade3%m_eff > 0 .and. cascade2%m_eff > 0 &
            .and. cascade1%m_eff == 0 &
            .and. (cascade2%on_shell .or. cascade2%mapping == RADIATION) &
            .and. abs (cascade3%m_eff - cascade2%m_eff) &
                  < cascade_set%m_threshold_s) &
            then
          keep = .true.
          cascade3%log_enhanced = .true.
          cascade3%mapping = RADIATION
       end if
    end if
    ! Non-singular cases, including failed resonances
    if (.not. keep) then
       ! Two on-shell particles from a virtual mother
       if (cascade1%on_shell .or. cascade2%on_shell) then
          keep = .true.
          cascade3%m_eff = max (cascade3%m_min, &
                                cascade1%m_eff + cascade2%m_eff)
          if (cascade3%m_eff < cascade_set%m_threshold_s) then
             cascade3%m_eff = 0
          end if
       end if
    end if
    ! Complete and register the cascade (two in case of resonance)
    if (keep) then
       cascade3%on_shell = cascade3%resonant .or. cascade3%log_enhanced
       if (cascade3%resonant) then
          cascade3%pdg = cascade3%flv%get_pdg ()
          if (cascade_set%keep_nonresonant) then
             allocate (cascade4)
             cascade4 = cascade3
             cascade4%index = cascade_index ()
             cascade4%pdg = UNDEFINED
             cascade4%mapping = NO_MAPPING
             cascade4%resonant = .false.
             cascade4%on_shell = .false.
          end if
          cascade3%m_min = cascade3%m_rea
          call cascade_fusion (cascade_set, cascade1, cascade2, cascade3)
          if (cascade_set%keep_nonresonant) then
             call cascade_fusion (cascade_set, cascade1, cascade2, cascade4)
          end if
       else
          call cascade_fusion (cascade_set, cascade1, cascade2, cascade3)
       end if
    else
       deallocate (cascade3)
    end if
  contains
    subroutine warn_decay (flv)
      type(flavor_t), intent(in) :: flv
      integer :: i
      integer, dimension(MAX_WARN_RESONANCE), save :: warned_code = 0
      LOOP_WARNED: do i = 1, MAX_WARN_RESONANCE
         if (warned_code(i) == 0) then
            warned_code(i) = flv%get_pdg ()
            write (msg_buffer, "(A)") &
                 & " Intermediate decay of zero-width particle " &
                 & // char (flv%get_name ()) &
                 & // " may be possible."
            call msg_warning
            exit LOOP_WARNED
         else if (warned_code(i) == flv%get_pdg ()) then
            exit LOOP_WARNED
         end if
      end do LOOP_WARNED
    end subroutine warn_decay
  end subroutine cascade_combine_s

  subroutine cascade_combine_t (cascade_set, cascade1, cascade2, flv)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), intent(in), target :: cascade1, cascade2
    type(flavor_t), intent(in) :: flv
    type(cascade_t), pointer :: cascade3
    allocate (cascade3)
    call cascade_init (cascade3, cascade1%depth + cascade2%depth + 1)
    cascade3%bincode = ior (cascade1%bincode, cascade2%bincode)
    cascade3%flv = flv%anti ()
    cascade3%pdg = abs (cascade3%flv%get_pdg ())
    cascade3%is_vector = flv%get_spin_type () == VECTOR
    if (cascade1%incoming) then
       cascade3%m_min = cascade2%m_min
    else
       cascade3%m_min = cascade1%m_min + cascade2%m_min
    end if
    cascade3%m_rea = flv%get_mass ()
    if (cascade3%m_rea > cascade_set%m_threshold_t) then
       cascade3%m_eff = max (cascade3%m_rea, cascade2%m_eff)
    else if (cascade2%m_eff > cascade_set%m_threshold_t) then
       cascade3%m_eff = cascade2%m_eff
    else
       cascade3%m_eff = 0
    end if
    ! Allowed decay of beam particle
    if (cascade1%incoming &
         .and. cascade1%m_rea > cascade2%m_rea + cascade3%m_rea) then
         call beam_decay (cascade_set%fatal_beam_decay)
    ! Massless splitting
    else if (cascade1%m_eff == 0 &
         .and. cascade2%m_eff < cascade_set%m_threshold_t &
         .and. cascade3%m_eff == 0) then
       cascade3%mapping = U_CHANNEL
       cascade3%log_enhanced = .true.
    ! IR radiation off massive particle
    else if (cascade1%m_eff /= 0 .and. cascade2%m_eff == 0 &
         .and. cascade3%m_eff /= 0 &
         .and. (cascade1%on_shell .or. cascade1%mapping == RADIATION) &
         .and. abs (cascade1%m_eff - cascade3%m_eff) &
               < cascade_set%m_threshold_t) &
         then
       cascade3%pdg = flv%get_pdg ()
       cascade3%log_enhanced = .true.
       cascade3%mapping = RADIATION
    end if
    cascade3%t_channel = .true.
    call cascade_fusion (cascade_set, cascade1, cascade2, cascade3)
  contains
    subroutine beam_decay (fatal_beam_decay)
      logical, intent(in) :: fatal_beam_decay
      write (msg_buffer, "(1x,A,1x,'->',1x,A,1x,A)") &
           char (cascade1%flv%get_name ()), &
           char (cascade3%flv%get_name ()), &
           char (cascade2%flv%get_name ())
      call msg_message
      write (msg_buffer, "(1x,'mass(',A,') =',1x,E17.10)") &
           char (cascade1%flv%get_name ()), cascade1%m_rea
      call msg_message
      write (msg_buffer, "(1x,'mass(',A,') =',1x,E17.10)") &
           char (cascade3%flv%get_name ()), cascade3%m_rea
      call msg_message
      write (msg_buffer, "(1x,'mass(',A,') =',1x,E17.10)") &
           char (cascade2%flv%get_name ()), cascade2%m_rea
      call msg_message
      if (fatal_beam_decay) then
         call msg_fatal (" Phase space: Initial beam particle can decay")
      else
         call msg_warning (" Phase space: Initial beam particle can decay")
      end if
    end subroutine beam_decay
  end subroutine cascade_combine_t

  subroutine cascade_combine_keystone &
       (cascade_set, cascade1, cascade2, cascade3, s_channel)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), intent(in), target :: cascade1, cascade2, cascade3
    logical, intent(in) :: s_channel
    type(cascade_t), pointer :: cascade4, cascade0
    logical :: keep, ok
    keep = .false.
    allocate (cascade4)
    call cascade_init &
         (cascade4, cascade1%depth + cascade2%depth + cascade3%depth)
    cascade4%complete = .true.
    if (s_channel) then
       cascade4%bincode = ior (cascade1%bincode, cascade2%bincode)
    else
       cascade4%bincode = cascade3%bincode
    end if
    cascade4%flv = cascade3%flv
    cascade4%pdg = cascade3%pdg
    cascade4%mapping = EXTERNAL_PRT
    cascade4%is_vector = cascade3%is_vector
    cascade4%m_min = cascade1%m_min + cascade2%m_min
    cascade4%m_rea = cascade3%m_rea
    cascade4%m_eff = cascade3%m_rea
    if (cascade4%m_min < cascade_set%sqrts) then
       keep = .true.
    end if
    if (keep) then
       if (cascade1%incoming .and. cascade2%log_enhanced) then
          allocate (cascade0)
          cascade0 = cascade2
          cascade0%next => null ()
          cascade0%index = cascade_index ()
          cascade0%mapping = NO_MAPPING
          cascade0%log_enhanced = .false.
          cascade0%n_log_enhanced = cascade0%n_log_enhanced - 1
          cascade0%tree_mapping(cascade0%depth) = NO_MAPPING
          call cascade_keystone &
               (cascade_set, cascade1, cascade0, cascade3, cascade4, ok)
          if (ok) then
             call cascade_set_add (cascade_set, cascade0, ok)
          else
             deallocate (cascade0)
          end if
       else if (cascade1%t_channel .and. cascade1%mapping == U_CHANNEL) then
          allocate (cascade0)
          cascade0 = cascade1
          cascade0%next => null ()
          cascade0%index = cascade_index ()
          cascade0%mapping = T_CHANNEL
          cascade0%tree_mapping(cascade0%depth) = T_CHANNEL
          call cascade_keystone &
               (cascade_set, cascade0, cascade2, cascade3, cascade4, ok)
          if (ok) then
             call cascade_set_add (cascade_set, cascade0, ok)
          else
             deallocate (cascade0)
          end if
       else if (cascade1%incoming .and. cascade2%depth == 1) then
          allocate (cascade0)
          cascade0 = cascade2
          cascade0%next => null ()
          cascade0%index = cascade_index ()
          cascade0%mapping = ON_SHELL
          cascade0%tree_mapping(cascade0%depth) = ON_SHELL
          call cascade_keystone &
               (cascade_set, cascade1, cascade0, cascade3, cascade4, ok)
          if (ok) then
             call cascade_set_add (cascade_set, cascade0, ok)
          else
             deallocate (cascade0)
          end if
       else
          call cascade_keystone &
               (cascade_set, cascade1, cascade2, cascade3, cascade4, ok)
       end if
    else
       deallocate (cascade4)
    end if
  end subroutine cascade_combine_keystone

  subroutine cascade_fusion (cascade_set, cascade1, cascade2, cascade3)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), intent(in), target :: cascade1, cascade2
    type(cascade_t), pointer :: cascade3
    integer :: i1, i2, i3, i4
    logical :: ok
    cascade3%internal = (cascade3%depth - 3) / 2
    if (cascade3%resonant) then
       cascade3%multiplicity = 1
       cascade3%n_resonances = &
            cascade1%n_resonances + cascade2%n_resonances + 1
    else
       cascade3%multiplicity = cascade1%multiplicity + cascade2%multiplicity
       cascade3%n_resonances = cascade1%n_resonances + cascade2%n_resonances
    end if
    if (cascade3%log_enhanced) then
       cascade3%n_log_enhanced = &
            cascade1%n_log_enhanced + cascade2%n_log_enhanced + 1
    else
       cascade3%n_log_enhanced = &
            cascade1%n_log_enhanced + cascade2%n_log_enhanced
    end if
    if (cascade3%resonant) then
       cascade3%n_off_shell = 0
    else if (cascade3%log_enhanced) then
       cascade3%n_off_shell = cascade1%n_off_shell + cascade2%n_off_shell
    else
       cascade3%n_off_shell = cascade1%n_off_shell + cascade2%n_off_shell + 1
    end if
    if (cascade3%t_channel) then
       cascade3%n_t_channel = cascade1%n_t_channel + 1
    end if
    if (cascade3%n_off_shell > cascade_set%off_shell) then
       deallocate (cascade3)
    else if (cascade3%n_t_channel > cascade_set%t_channel) then
       deallocate (cascade3)
    else
       i1 = cascade1%depth
       i2 = i1 + 1
       i3 = i1 + cascade2%depth
       i4 = cascade3%depth
       cascade3%tree(:i1) = cascade1%tree
       where (cascade1%tree_mapping > NO_MAPPING)
          cascade3%tree_pdg(:i1) = cascade1%tree_pdg
       elsewhere
          cascade3%tree_pdg(:i1) = UNDEFINED
       end where
       cascade3%tree_mapping(:i1) = cascade1%tree_mapping
       cascade3%tree_resonant(:i1) = cascade1%tree_resonant
       cascade3%tree(i2:i3) = cascade2%tree
       where (cascade2%tree_mapping > NO_MAPPING)
          cascade3%tree_pdg(i2:i3) = cascade2%tree_pdg
       elsewhere
          cascade3%tree_pdg(i2:i3) = UNDEFINED
       end where
       cascade3%tree_mapping(i2:i3) = cascade2%tree_mapping
       cascade3%tree_resonant(i2:i3) = cascade2%tree_resonant
       cascade3%tree(i4) = cascade3%bincode
       cascade3%tree_pdg(i4) = cascade3%pdg
       cascade3%tree_mapping(i4) = cascade3%mapping
       cascade3%tree_resonant(i4) = cascade3%resonant
       call tree_sort (cascade3%tree, &
            cascade3%tree_pdg, cascade3%tree_mapping, cascade3%tree_resonant)
       cascade3%has_children = .true.
       cascade3%daughter1 => cascade1
       cascade3%daughter2 => cascade2
       call cascade_set_add (cascade_set, cascade3, ok)
       if (.not. ok)  deallocate (cascade3)
    end if
  end subroutine cascade_fusion

  subroutine cascade_keystone &
       (cascade_set, cascade1, cascade2, cascade3, cascade4, ok)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), intent(in), target :: cascade1, cascade2, cascade3
    type(cascade_t), pointer :: cascade4
    logical, intent(out) :: ok
    integer :: i1, i2, i3, i4
    cascade4%internal = (cascade4%depth - 3) / 2
    cascade4%multiplicity = cascade1%multiplicity + cascade2%multiplicity
    cascade4%n_resonances = cascade1%n_resonances + cascade2%n_resonances
    cascade4%n_off_shell = cascade1%n_off_shell + cascade2%n_off_shell
    cascade4%n_log_enhanced = &
            cascade1%n_log_enhanced + cascade2%n_log_enhanced
    cascade4%n_t_channel = cascade1%n_t_channel + cascade2%n_t_channel
    if (cascade4%n_off_shell > cascade_set%off_shell) then
       deallocate (cascade4)
       ok = .false.
    else if (cascade4%n_t_channel > cascade_set%t_channel) then
       deallocate (cascade4)
       ok = .false.
    else
       i1 = cascade1%depth
       i2 = i1 + 1
       i3 = i1 + cascade2%depth
       i4 = cascade4%depth
       cascade4%tree(:i1) = cascade1%tree
       where (cascade1%tree_mapping > NO_MAPPING)
          cascade4%tree_pdg(:i1) = cascade1%tree_pdg
       elsewhere
          cascade4%tree_pdg(:i1) = UNDEFINED
       end where
       cascade4%tree_mapping(:i1) = cascade1%tree_mapping
       cascade4%tree_resonant(:i1) = cascade1%tree_resonant
       cascade4%tree(i2:i3) = cascade2%tree
       where (cascade2%tree_mapping > NO_MAPPING)
          cascade4%tree_pdg(i2:i3) = cascade2%tree_pdg
       elsewhere
          cascade4%tree_pdg(i2:i3) = UNDEFINED
       end where
       cascade4%tree_mapping(i2:i3) = cascade2%tree_mapping
       cascade4%tree_resonant(i2:i3) = cascade2%tree_resonant
       cascade4%tree(i4) = cascade4%bincode
       cascade4%tree_pdg(i4) = UNDEFINED
       cascade4%tree_mapping(i4) = cascade4%mapping
       cascade4%tree_resonant(i4) = .false.
       call tree_sort (cascade4%tree, &
            cascade4%tree_pdg, cascade4%tree_mapping, cascade4%tree_resonant)
       cascade4%has_children = .true.
       cascade4%daughter1 => cascade1
       cascade4%daughter2 => cascade2
       cascade4%mother => cascade3
       call cascade_set_add (cascade_set, cascade4, ok)
       if (ok) then
          if (.not. associated (cascade_set%first_k)) then
             cascade_set%first_k => cascade4
          end if
       else
          deallocate (cascade4)
       end if
    end if
  end subroutine cascade_keystone

  subroutine tree_sort (tree, pdg, mapping, resonant)
    integer(TC), dimension(:), intent(inout) :: tree
    integer, dimension(:), intent(inout) :: pdg, mapping
    logical, dimension(:), intent(inout) :: resonant
    integer(TC), dimension(size(tree)) :: tree_tmp
    integer, dimension(size(pdg)) :: pdg_tmp, mapping_tmp
    logical, dimension(size(resonant)) :: resonant_tmp
    integer, dimension(1) :: pos
    integer :: i
    tree_tmp = tree
    pdg_tmp = pdg
    mapping_tmp = mapping
    resonant_tmp = resonant
    do i = size(tree),1,-1
       pos = maxloc (tree_tmp)
       tree(i) = tree_tmp (pos(1))
       pdg(i) = pdg_tmp (pos(1))
       mapping(i) = mapping_tmp (pos(1))
       resonant(i) = resonant_tmp (pos(1))
       tree_tmp(pos(1)) = 0
    end do
  end subroutine tree_sort

  subroutine cascade_set_generate_s (cascade_set)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), pointer :: cascade1, cascade2
    cascade1 => cascade_set%first
    LOOP1: do while (associated (cascade1))
       cascade2 => cascade_set%first
       LOOP2: do while (associated (cascade2))
          if (cascade2%index >= cascade1%index)  exit LOOP2
          if (cascade1 .disjunct. cascade2) then
             call cascade_match_pair (cascade_set, cascade1, cascade2, .true.)
          end if
          call terminate_now_if_signal ()
          cascade2 => cascade2%next
       end do LOOP2
       cascade1 => cascade1%next
    end do LOOP1
  end subroutine cascade_set_generate_s

  subroutine cascade_set_generate_t (cascade_set, pos_seed, pos_target)
    type(cascade_set_t), intent(inout), target :: cascade_set
    integer, intent(in) :: pos_seed, pos_target
    type(cascade_t), pointer :: cascade_seed, cascade_target
    type(cascade_t), pointer :: cascade1, cascade2
    integer(TC) :: bc_seed, bc_target
    bc_seed = ibset (0_TC, pos_seed-1)
    bc_target = ibset (0_TC, pos_target-1)
    cascade_seed => cascade_set%first_t
    LOOP_SEED: do while (associated (cascade_seed))
       if (cascade_seed%bincode == bc_seed) then
          cascade_target => cascade_set%first_t
          LOOP_TARGET: do while (associated (cascade_target))
             if (cascade_target%bincode == bc_target) then
                cascade1 => cascade_set%first_t
                LOOP_T: do while (associated (cascade1))
                   if ((cascade1 .disjunct. cascade_target) &
                        .and. .not. (cascade1 .disjunct. cascade_seed)) then
                      cascade2 => cascade_set%first
                      LOOP_S: do while (associated (cascade2))
                         if ((cascade2 .disjunct. cascade_target) &
                              .and. (cascade2 .disjunct. cascade1)) then
                            call cascade_match_pair &
                                 (cascade_set, cascade1, cascade2, .false.)
                         end if
                         call terminate_now_if_signal ()
                         cascade2 => cascade2%next
                      end do LOOP_S
                   end if
                   call terminate_now_if_signal ()
                   cascade1 => cascade1%next
                end do LOOP_T
             end if
             call terminate_now_if_signal ()
             cascade_target => cascade_target%next
          end do LOOP_TARGET
       end if
       call terminate_now_if_signal ()
       cascade_seed => cascade_seed%next
    end do LOOP_SEED
  end subroutine cascade_set_generate_t

  subroutine cascade_set_generate_decay (cascade_set)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), pointer :: cascade1, cascade2
    type(cascade_t), pointer :: cascade_in
    cascade_in => cascade_set%first_t
    cascade1 => cascade_set%first
    do while (associated (cascade1))
       if (cascade1 .disjunct. cascade_in) then
          cascade2 => cascade1%next
          do while (associated (cascade2))
             if ((cascade2 .disjunct. cascade1) &
                  .and. (cascade2 .disjunct. cascade_in)) then
                call cascade_match_triplet (cascade_set, &
                     cascade1, cascade2, cascade_in, .true.)
             end if
             call terminate_now_if_signal ()
             cascade2 => cascade2%next
          end do
       end if
       call terminate_now_if_signal ()
       cascade1 => cascade1%next
    end do
  end subroutine cascade_set_generate_decay

  subroutine cascade_set_generate_scattering &
       (cascade_set, ns1, ns2, nt1, nt2, pos_seed, pos_target)
    type(cascade_set_t), intent(inout), target :: cascade_set
    integer, intent(in) :: pos_seed, pos_target
    integer, intent(in) :: ns1, ns2, nt1, nt2
    type(cascade_t), pointer :: cascade_seed, cascade_target
    type(cascade_t), pointer :: cascade1, cascade2
    integer(TC) :: bc_seed, bc_target
    bc_seed = ibset (0_TC, pos_seed-1)
    bc_target = ibset (0_TC, pos_target-1)
    cascade_seed => cascade_set%first_t
    LOOP_SEED: do while (associated (cascade_seed))
       if (cascade_seed%index < ns1) then
          cascade_seed => cascade_seed%next
          cycle LOOP_SEED
       else if (cascade_seed%index > ns2) then
          exit LOOP_SEED
       else if (cascade_seed%bincode == bc_seed) then
          cascade_target => cascade_set%first_t
          LOOP_TARGET: do while (associated (cascade_target))
             if (cascade_target%index < nt1) then
                cascade_target => cascade_target%next
                cycle LOOP_TARGET
             else if (cascade_target%index > nt2) then
                exit LOOP_TARGET
             else if (cascade_target%bincode == bc_target) then
                cascade1 => cascade_set%first_t
                LOOP_T: do while (associated (cascade1))
                   if ((cascade1 .disjunct. cascade_target) &
                        .and. .not. (cascade1 .disjunct. cascade_seed)) then
                      cascade2 => cascade_set%first
                      LOOP_S: do while (associated (cascade2))
                         if ((cascade2 .disjunct. cascade_target) &
                              .and. (cascade2 .disjunct. cascade1)) then
                            call cascade_match_triplet (cascade_set, &
                                 cascade1, cascade2, cascade_target, .false.)
                         end if
                         call terminate_now_if_signal ()
                         cascade2 => cascade2%next
                      end do LOOP_S
                   end if
                   call terminate_now_if_signal ()
                   cascade1 => cascade1%next
                end do LOOP_T
             end if
             call terminate_now_if_signal ()
             cascade_target => cascade_target%next
          end do LOOP_TARGET
       end if
       call terminate_now_if_signal ()
       cascade_seed => cascade_seed%next
    end do LOOP_SEED
  end subroutine cascade_set_generate_scattering

  subroutine cascade_set_assign_resonance_hash (cascade_set)
    type(cascade_set_t), intent(inout) :: cascade_set
    type(cascade_t), pointer :: cascade
    cascade => cascade_set%first_k
    do while (associated (cascade))
       call cascade_assign_resonance_hash (cascade)
       cascade => cascade%next
    end do
  end subroutine cascade_set_assign_resonance_hash

  subroutine cascade_set_assign_groves (cascade_set)
    type(cascade_set_t), intent(inout), target :: cascade_set
    type(cascade_t), pointer :: cascade1, cascade2
    integer :: multiplicity
    integer :: n_resonances, n_log_enhanced, n_t_channel, n_off_shell
    integer :: res_hash
    integer :: grove
    grove = 0
    cascade1 => cascade_set%first_k
    do while (associated (cascade1))
       if (cascade1%active .and. cascade1%complete &
            .and. cascade1%grove == 0) then
          grove = grove + 1
          cascade1%grove = grove
          multiplicity = cascade1%multiplicity
          n_resonances = cascade1%n_resonances
          n_log_enhanced = cascade1%n_log_enhanced
          n_off_shell = cascade1%n_off_shell
          n_t_channel = cascade1%n_t_channel
          res_hash = cascade1%res_hash
          cascade2 => cascade1%next
          do while (associated (cascade2))
             if (cascade2%grove == 0) then
                if (cascade2%multiplicity == multiplicity &
                     .and. cascade2%n_resonances == n_resonances &
                     .and. cascade2%n_log_enhanced == n_log_enhanced &
                     .and. cascade2%n_off_shell == n_off_shell &
                     .and. cascade2%n_t_channel == n_t_channel &
                     .and. cascade2%res_hash == res_hash) then
                   cascade2%grove = grove
                end if
             end if
             call terminate_now_if_signal ()
             cascade2 => cascade2%next
          end do
       end if
       call terminate_now_if_signal ()
       cascade1 => cascade1%next
    end do
    cascade_set%n_groves = grove
  end subroutine cascade_set_assign_groves

  module subroutine cascade_set_generate &
       (cascade_set, model, n_in, n_out, flv, phs_par, fatal_beam_decay)
    type(cascade_set_t), intent(out) :: cascade_set
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: n_in, n_out
    type(flavor_t), dimension(:,:), intent(in) :: flv
    type(phs_parameters_t), intent(in) :: phs_par
    logical, intent(in) :: fatal_beam_decay
    type(cascade_set_t), dimension(:), allocatable :: cset
    type(cascade_t), pointer :: cascade
    integer :: i
    if (phase_space_vanishes (phs_par%sqrts, n_in, flv))  return
    call cascade_set_init (cascade_set, model, n_in, n_out, phs_par, &
       fatal_beam_decay, flv)
    allocate (cset (size (flv, 2)))
    do i = 1, size (cset)
       call cascade_set_generate_single (cset(i), &
            model, n_in, n_out, flv(:,i), phs_par, fatal_beam_decay)
       cascade => cset(i)%first_k
       do while (associated (cascade))
          if (cascade%active .and. cascade%complete) then
             call cascade_set_add_copy (cascade_set, cascade)
          end if
          cascade => cascade%next
       end do
       call cascade_set_final (cset(i))
    end do
    cascade_set%first_k => cascade_set%first
    call cascade_set_assign_resonance_hash (cascade_set)
    call cascade_set_assign_groves (cascade_set)
  end subroutine cascade_set_generate

  subroutine cascade_set_generate_single (cascade_set, &
      model, n_in, n_out, flv, phs_par, fatal_beam_decay)
    type(cascade_set_t), intent(out) :: cascade_set
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: n_in, n_out
    type(flavor_t), dimension(:), intent(in) :: flv
    type(phs_parameters_t), intent(in) :: phs_par
    logical, intent(in) :: fatal_beam_decay
    integer :: n11, n12, n21, n22
    call cascade_set_init (cascade_set, model, n_in, n_out, phs_par, &
       fatal_beam_decay)
    call cascade_set_add_outgoing (cascade_set, flv(n_in+1:))
    call cascade_set_generate_s (cascade_set)
    select case (n_in)
    case(1)
       call cascade_set_add_incoming &
            (cascade_set, n11, n12, n_out + 1, flv(1))
       call cascade_set_generate_decay (cascade_set)
    case(2)
       call cascade_set_add_incoming &
            (cascade_set, n11, n12, n_out + 1, flv(2))
       call cascade_set_add_incoming &
            (cascade_set, n21, n22, n_out + 2, flv(1))
       call cascade_set_generate_t (cascade_set, n_out + 1, n_out + 2)
       call cascade_set_generate_t (cascade_set, n_out + 2, n_out + 1)
       call cascade_set_generate_scattering &
            (cascade_set, n11, n12, n21, n22, n_out + 1, n_out + 2)
       call cascade_set_generate_scattering &
            (cascade_set, n21, n22, n11, n12, n_out + 2, n_out + 1)
    end select
  end subroutine cascade_set_generate_single

  module function phase_space_vanishes (sqrts, n_in, flv) result (flag)
    logical :: flag
    real(default), intent(in) :: sqrts
    integer, intent(in) :: n_in
    type(flavor_t), dimension(:,:), intent(in) :: flv
    real(default), dimension(:,:), allocatable :: mass
    real(default), dimension(:), allocatable :: mass_in, mass_out
    integer :: n_prt, n_flv, i, j
    flag = .false.
    if (sqrts <= 0) then
       call msg_error ("Phase space vanishes (sqrts must be positive)")
       flag = .true.;  return
    end if
    n_prt = size (flv, 1)
    n_flv = size (flv, 2)
    allocate (mass (n_prt, n_flv), mass_in (n_flv), mass_out (n_flv))
    mass = flv%get_mass ()
    mass_in = sum (mass(:n_in,:), 1)
    mass_out = sum (mass(n_in+1:,:), 1)
    if (any (mass_in > sqrts)) then
       call msg_error ("Mass sum of incoming particles " &
            // "is more than available energy")
       flag = .true.;  return
    end if
    if (any (mass_out > sqrts)) then
       call msg_error ("Mass sum of outgoing particles " &
            // "is more than available energy")
       flag = .true.;  return
    end if
  end function phase_space_vanishes

  module subroutine cascade_extract_resonance_history &
       (cascade, res_hist, model, n_out)
    class(cascade_t), intent(in), target :: cascade
    type(resonance_history_t), intent(out) :: res_hist
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: n_out
    type(resonance_info_t) :: resonance
    integer :: i, mom_id, pdg
    if (debug_on) call msg_debug2 (D_PHASESPACE, "cascade_extract_resonance_history")
    if (cascade%n_resonances > 0) then
       if (cascade%has_children) then
          if (debug_on) call msg_debug2 (D_PHASESPACE, "cascade has resonances and children")
          do i = 1, size(cascade%tree_resonant)
             if (cascade%tree_resonant (i)) then
                mom_id = cascade%tree (i)
                pdg = cascade%tree_pdg (i)
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
  end subroutine cascade_extract_resonance_history

  module function cascade_set_get_n_trees (cascade_set) result (n)
    type(cascade_set_t), intent(in), target :: cascade_set
    integer :: n
    type(cascade_t), pointer :: cascade
    integer :: grove
    if (debug_on) call msg_debug (D_PHASESPACE, "cascade_set_get_n_trees")
    n = 0
    do grove = 1, cascade_set%n_groves
       cascade => cascade_set%first_k
       do while (associated (cascade))
          if (cascade%active .and. cascade%complete) then
             if (cascade%grove == grove) then
                n = n + 1
             end if
          end if
          cascade => cascade%next
       end do
    end do
    if (debug_on) call msg_debug (D_PHASESPACE, "n", n)
  end function cascade_set_get_n_trees

  module subroutine cascade_set_get_resonance_histories &
       (cascade_set, n_filter, res_hists)
    type(cascade_set_t), intent(in), target :: cascade_set
    integer, intent(in), optional :: n_filter
    type(resonance_history_t), dimension(:), allocatable, intent(out) :: &
         res_hists
    type(resonance_history_t), dimension(:), allocatable :: tmp
    type(cascade_t), pointer :: cascade
    type(resonance_history_t) :: res_hist
    type(resonance_history_set_t) :: res_hist_set
    integer :: grove, i, n_hists
    logical :: included, add_to_list
    if (debug_on)  call msg_debug &
         (D_PHASESPACE, "cascade_set_get_resonance_histories")
    call res_hist_set%init (n_filter = n_filter)
    do grove = 1, cascade_set%n_groves
       cascade => cascade_set%first_k
       do while (associated (cascade))
          if (cascade%active .and. cascade%complete) then
             if (cascade%grove == grove) then
                if (debug_on) call msg_debug2 (D_PHASESPACE, "grove", grove)
                call cascade%extract_resonance_history &
                     (res_hist, cascade_set%model, cascade_set%n_out)
                call res_hist_set%enter (res_hist)
             end if
          end if
          cascade => cascade%next
       end do
    end do
    call res_hist_set%freeze ()
    call res_hist_set%to_array (res_hists)
  end subroutine cascade_set_get_resonance_histories


end submodule cascades_s

