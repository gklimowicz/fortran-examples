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

submodule (subevents) subevents_s

  use io_units
  use format_defs, only: FMT_14, FMT_19
  use format_utils, only: pac_fmt
  use physics_defs
  use sorting

  implicit none

contains

  subroutine prt_init_beam (prt, pdg, p, p2, src)
    type(prt_t), intent(out) :: prt
    integer, intent(in) :: pdg
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: p2
    integer, dimension(:), intent(in) :: src
    prt%type = PRT_BEAM
    call prt_set (prt, pdg, - p, p2, src)
  end subroutine prt_init_beam

  subroutine prt_init_incoming (prt, pdg, p, p2, src)
    type(prt_t), intent(out) :: prt
    integer, intent(in) :: pdg
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: p2
    integer, dimension(:), intent(in) :: src
    prt%type = PRT_INCOMING
    call prt_set (prt, pdg, - p, p2, src)
  end subroutine prt_init_incoming

  subroutine prt_init_outgoing (prt, pdg, p, p2, src)
    type(prt_t), intent(out) :: prt
    integer, intent(in) :: pdg
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: p2
    integer, dimension(:), intent(in) :: src
    prt%type = PRT_OUTGOING
    call prt_set (prt, pdg, p, p2, src)
  end subroutine prt_init_outgoing

  subroutine prt_init_composite (prt, p, src)
    type(prt_t), intent(out) :: prt
    type(vector4_t), intent(in) :: p
    integer, dimension(:), intent(in) :: src
    prt%type = PRT_COMPOSITE
    call prt_set (prt, 0, p, p**2, src)
  end subroutine prt_init_composite

  module subroutine prt_init_combine (prt, prt1, prt2)
    type(prt_t), intent(out) :: prt
    type(prt_t), intent(in) :: prt1, prt2
    type(vector4_t) :: p
    integer, dimension(0) :: src
    prt%type = PRT_COMPOSITE
    p = prt1%p + prt2%p
    call prt_set (prt, 0, p, p**2, src)
  end subroutine prt_init_combine

  subroutine prt_init_pseudojet (prt, jet, src, pdg, is_b_jet, is_c_jet)
    type(prt_t), intent(out) :: prt
    type(pseudojet_t), intent(in) :: jet
    integer, dimension(:), intent(in) :: src
    integer, intent(in) :: pdg
    logical, intent(in) :: is_b_jet, is_c_jet
    type(vector4_t) :: p
    prt%type = PRT_COMPOSITE
    p = vector4_moving (jet%e(), &
         vector3_moving ([jet%px(), jet%py(), jet%pz()]))
    call prt_set (prt, pdg, p, p**2, src)
    prt%is_b_jet = is_b_jet
    prt%is_c_jet = is_c_jet
    prt%clustered = .true.
  end subroutine prt_init_pseudojet

  elemental module function prt_get_pdg (prt) result (pdg)
    integer :: pdg
    type(prt_t), intent(in) :: prt
    pdg = prt%pdg
  end function prt_get_pdg

  elemental module function prt_get_momentum (prt) result (p)
    type(vector4_t) :: p
    type(prt_t), intent(in) :: prt
    p = prt%p
  end function prt_get_momentum

  elemental module function prt_get_msq (prt) result (msq)
    real(default) :: msq
    type(prt_t), intent(in) :: prt
    msq = prt%p2
  end function prt_get_msq

  elemental module function prt_is_polarized (prt) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt
    flag = prt%polarized
  end function prt_is_polarized

  elemental module function prt_get_helicity (prt) result (h)
    integer :: h
    type(prt_t), intent(in) :: prt
    h = prt%h
  end function prt_get_helicity

  elemental module function prt_is_colorized (prt) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt
    flag = prt%colorized
  end function prt_is_colorized

  elemental module function prt_is_clustered (prt) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt
    flag = prt%clustered
  end function prt_is_clustered

  elemental module function prt_is_recombinable (prt) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt
    flag = prt_is_parton (prt) .or. &
           abs(prt%pdg) == TOP_Q .or. &
           prt_is_lepton (prt) .or. &
           prt_is_photon (prt)
  end function prt_is_recombinable

  elemental module function prt_is_photon (prt) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt
    flag = prt%pdg == PHOTON
  end function prt_is_photon

  elemental module function prt_is_parton (prt) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt
    flag = abs(prt%pdg) == DOWN_Q .or. &
           abs(prt%pdg) == UP_Q .or. &
           abs(prt%pdg) == STRANGE_Q .or. &
           abs(prt%pdg) == CHARM_Q .or. &
           abs(prt%pdg) == BOTTOM_Q .or. &
           prt%pdg == GLUON
  end function prt_is_parton

  elemental module function prt_is_lepton (prt) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt
    flag = abs(prt%pdg) == ELECTRON .or. &
           abs(prt%pdg) == MUON     .or. &
           abs(prt%pdg) == TAU
  end function prt_is_lepton

  elemental module function prt_is_b_jet (prt) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt
    flag = prt%is_b_jet
  end function prt_is_b_jet

  elemental module function prt_is_c_jet (prt) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt
    flag = prt%is_c_jet
  end function prt_is_c_jet

  elemental module function prt_get_n_col (prt) result (n)
    integer :: n
    type(prt_t), intent(in) :: prt
    integer, dimension(:), allocatable :: col, acl
    integer :: i
    n = 0
    if (prt%colorized) then
       do i = 1, size (prt%col)
          if (all (prt%col(i) /= prt%acl)) n = n + 1
       end do
    end if
  end function prt_get_n_col

  elemental module function prt_get_n_acl (prt) result (n)
    integer :: n
    type(prt_t), intent(in) :: prt
    integer, dimension(:), allocatable :: col, acl
    integer :: i
    n = 0
    if (prt%colorized) then
       do i = 1, size (prt%acl)
          if (all (prt%acl(i) /= prt%col)) n = n + 1
       end do
    end if
  end function prt_get_n_acl

  module subroutine prt_get_color_indices (prt, col, acl)
    type(prt_t), intent(in) :: prt
    integer, dimension(:), allocatable, intent(out) :: col, acl
    if (prt%colorized) then
       col = prt%col
       acl = prt%acl
    else
       col = [integer::]
       acl = [integer::]
    end if
  end subroutine prt_get_color_indices

  subroutine prt_set (prt, pdg, p, p2, src)
    type(prt_t), intent(inout) :: prt
    integer, intent(in) :: pdg
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: p2
    integer, dimension(:), intent(in) :: src
    prt%pdg = pdg
    prt%p = p
    prt%p2 = p2
    if (allocated (prt%src)) then
       if (size (src) /= size (prt%src)) then
          deallocate (prt%src)
          allocate (prt%src (size (src)))
       end if
    else
       allocate (prt%src (size (src)))
    end if
    prt%src = src
  end subroutine prt_set

  elemental subroutine prt_set_pdg (prt, pdg)
    type(prt_t), intent(inout) :: prt
    integer, intent(in) :: pdg
    prt%pdg = pdg
  end subroutine prt_set_pdg

  elemental subroutine prt_set_p (prt, p)
    type(prt_t), intent(inout) :: prt
    type(vector4_t), intent(in) :: p
    prt%p = p
  end subroutine prt_set_p

  elemental subroutine prt_set_p2 (prt, p2)
    type(prt_t), intent(inout) :: prt
    real(default), intent(in) :: p2
    prt%p2 = p2
  end subroutine prt_set_p2

  subroutine prt_polarize (prt, h)
    type(prt_t), intent(inout) :: prt
    integer, intent(in) :: h
    prt%polarized = .true.
    prt%h = h
  end subroutine prt_polarize

  subroutine prt_colorize (prt, col, acl)
    type(prt_t), intent(inout) :: prt
    integer, dimension(:), intent(in) :: col, acl
    prt%colorized = .true.
    prt%col = col
    prt%acl = acl
  end subroutine prt_colorize

  elemental module function c_prt_from_prt (prt) result (c_prt)
    type(c_prt_t) :: c_prt
    type(prt_t), intent(in) :: prt
    c_prt = prt%p
    c_prt%type = prt%type
    c_prt%pdg = prt%pdg
    if (prt%polarized) then
       c_prt%polarized = 1
    else
       c_prt%polarized = 0
    end if
    c_prt%h = prt%h
  end function c_prt_from_prt

  module subroutine prt_write (prt, unit, testflag)
    type(prt_t), intent(in) :: prt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    logical :: pacified
    type(prt_t) :: tmp
    character(len=7) :: fmt
    integer :: u, i
    call pac_fmt (fmt, FMT_19, FMT_14, testflag)
    u = given_output_unit (unit);  if (u < 0)  return
    pacified = .false. ; if (present (testflag))  pacified = testflag
    tmp = prt
    if (pacified) call pacify (tmp)
    write (u, "(1x,A)", advance="no")  "prt("
    select case (prt%type)
    case (PRT_UNDEFINED);    write (u, "('?')", advance="no")
    case (PRT_BEAM);         write (u, "('b:')", advance="no")
    case (PRT_INCOMING);     write (u, "('i:')", advance="no")
    case (PRT_OUTGOING);     write (u, "('o:')", advance="no")
    case (PRT_COMPOSITE);    write (u, "('c:')", advance="no")
    end select
    select case (prt%type)
    case (PRT_BEAM, PRT_INCOMING, PRT_OUTGOING)
       if (prt%polarized) then
          write (u, "(I0,'/',I0,'|')", advance="no")  prt%pdg, prt%h
       else
          write (u, "(I0,'|')", advance="no") prt%pdg
       end if
    end select
    select case (prt%type)
    case (PRT_BEAM, PRT_INCOMING, PRT_OUTGOING, PRT_COMPOSITE)
       if (prt%colorized) then
          write (u, "(*(I0,:,','))", advance="no")  prt%col
          write (u, "('/')", advance="no")
          write (u, "(*(I0,:,','))", advance="no")  prt%acl
          write (u, "('|')", advance="no")
       end if
    end select
    select case (prt%type)
    case (PRT_BEAM, PRT_INCOMING, PRT_OUTGOING, PRT_COMPOSITE)
       write (u, "(" // FMT_14 // ",';'," // FMT_14 // ",','," // &
            FMT_14 // ",','," // FMT_14 // ")", advance="no") tmp%p
       write (u, "('|'," // fmt // ")", advance="no") tmp%p2
    end select
    if (allocated (prt%src)) then
       write (u, "('|')", advance="no")
       do i = 1, size (prt%src)
          write (u, "(1x,I0)", advance="no")  prt%src(i)
       end do
    end if
    if (prt%is_b_jet) then
       write (u, "('|b jet')", advance="no")
    end if
    if (prt%is_c_jet) then
       write (u, "('|c jet')", advance="no")
    end if
    write (u, "(A)")  ")"
  end subroutine prt_write

  elemental module function prt_match (prt1, prt2) result (match)
    logical :: match
    type(prt_t), intent(in) :: prt1, prt2
    if (size (prt1%src) == size (prt2%src)) then
       match = all (prt1%src == prt2%src)
    else
       match = .false.
    end if
  end function prt_match

  subroutine prt_combine (prt, prt_in1, prt_in2, ok)
    type(prt_t), intent(inout) :: prt
    type(prt_t), intent(in) :: prt_in1, prt_in2
    logical :: ok
    integer, dimension(:), allocatable :: src
    integer, dimension(:), allocatable :: col1, acl1, col2, acl2
    call combine_index_lists (src, prt_in1%src, prt_in2%src)
    ok = allocated (src)
    if (ok) then
       call prt_init_composite (prt, prt_in1%p + prt_in2%p, src)
       if (prt_in1%colorized .or. prt_in2%colorized) then
          select case (prt_in1%type)
          case default
             call prt_get_color_indices (prt_in1, col1, acl1)
          case (PRT_BEAM, PRT_INCOMING)
             call prt_get_color_indices (prt_in1, acl1, col1)
          end select
          select case (prt_in2%type)
          case default
             call prt_get_color_indices (prt_in2, col2, acl2)
          case (PRT_BEAM, PRT_INCOMING)
             call prt_get_color_indices (prt_in2, acl2, col2)
          end select
          call prt_colorize (prt, [col1, col2], [acl1, acl2])
       end if
    end if
  end subroutine prt_combine

  module function are_disjoint (prt_in1, prt_in2) result (flag)
    logical :: flag
    type(prt_t), intent(in) :: prt_in1, prt_in2
    flag = index_lists_are_disjoint (prt_in1%src, prt_in2%src)
  end function are_disjoint

  subroutine combine_index_lists (res, src1, src2)
    integer, dimension(:), intent(in) :: src1, src2
    integer, dimension(:), allocatable :: res
    integer :: i1, i2, i
    allocate (res (size (src1) + size (src2)))
    if (size (src1) == 0) then
       res = src2
       return
    else if (size (src2) == 0) then
       res = src1
       return
    end if
    i1 = 1
    i2 = 1
    LOOP: do i = 1, size (res)
       if (src1(i1) < src2(i2)) then
          res(i) = src1(i1);  i1 = i1 + 1
          if (i1 > size (src1)) then
             res(i+1:) = src2(i2:)
             exit LOOP
          end if
       else if (src1(i1) > src2(i2)) then
          res(i) = src2(i2);  i2 = i2 + 1
          if (i2 > size (src2)) then
             res(i+1:) = src1(i1:)
             exit LOOP
          end if
       else
          deallocate (res)
          exit LOOP
       end if
    end do LOOP
  end subroutine combine_index_lists

  function index_lists_are_disjoint (src1, src2) result (flag)
    logical :: flag
    integer, dimension(:), intent(in) :: src1, src2
    integer :: i1, i2, i
    flag = .true.
    i1 = 1
    i2 = 1
    LOOP: do i = 1, size (src1) + size (src2)
       if (src1(i1) < src2(i2)) then
          i1 = i1 + 1
          if (i1 > size (src1)) then
             exit LOOP
          end if
       else if (src1(i1) > src2(i2)) then
          i2 = i2 + 1
          if (i2 > size (src2)) then
             exit LOOP
          end if
       else
          flag = .false.
          exit LOOP
       end if
    end do LOOP
  end function index_lists_are_disjoint

  module subroutine subevt_init (subevt, n_active)
    type(subevt_t), intent(out) :: subevt
    integer, intent(in), optional :: n_active
    if (present (n_active))  subevt%n_active = n_active
    allocate (subevt%prt (subevt%n_active))
  end subroutine subevt_init

  module subroutine subevt_reset (subevt, n_active)
    class(subevt_t), intent(inout) :: subevt
    integer, intent(in) :: n_active
    subevt%n_active = n_active
    if (subevt%n_active > size (subevt%prt)) then
       deallocate (subevt%prt)
       allocate (subevt%prt (subevt%n_active))
    end if
  end subroutine subevt_reset

  module subroutine subevt_write (object, unit, prefix, pacified)
    class(subevt_t), intent(in) :: object
    integer, intent(in), optional :: unit
    character(*), intent(in), optional :: prefix
    logical, intent(in), optional :: pacified
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)") "subevent:"
    do i = 1, object%n_active
       if (present (prefix))  write (u, "(A)", advance="no") prefix
       write (u, "(1x,I0)", advance="no")  i
       call prt_write (object%prt(i), unit = unit, testflag = pacified)
    end do
  end subroutine subevt_write

  module subroutine subevt_assign (subevt, subevt_in)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: subevt_in
    if (.not. allocated (subevt%prt)) then
       call subevt_init (subevt, subevt_in%n_active)
    else
       call subevt%reset (subevt_in%n_active)
    end if
    subevt%prt(:subevt%n_active) = subevt_in%prt(:subevt%n_active)
  end subroutine subevt_assign

  module subroutine subevt_set_beam (subevt, i, pdg, p, p2, src)
    class(subevt_t), intent(inout) :: subevt
    integer, intent(in) :: i
    integer, intent(in) :: pdg
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: p2
    integer, dimension(:), intent(in), optional :: src
    if (present (src)) then
       call prt_init_beam (subevt%prt(i), pdg, p, p2, src)
    else
       call prt_init_beam (subevt%prt(i), pdg, p, p2, [i])
    end if
  end subroutine subevt_set_beam

  module subroutine subevt_set_incoming (subevt, i, pdg, p, p2, src)
    class(subevt_t), intent(inout) :: subevt
    integer, intent(in) :: i
    integer, intent(in) :: pdg
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: p2
    integer, dimension(:), intent(in), optional :: src
    if (present (src)) then
       call prt_init_incoming (subevt%prt(i), pdg, p, p2, src)
    else
       call prt_init_incoming (subevt%prt(i), pdg, p, p2, [i])
    end if
  end subroutine subevt_set_incoming

  module subroutine subevt_set_outgoing (subevt, i, pdg, p, p2, src)
    class(subevt_t), intent(inout) :: subevt
    integer, intent(in) :: i
    integer, intent(in) :: pdg
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: p2
    integer, dimension(:), intent(in), optional :: src
    if (present (src)) then
       call prt_init_outgoing (subevt%prt(i), pdg, p, p2, src)
    else
       call prt_init_outgoing (subevt%prt(i), pdg, p, p2, [i])
    end if
  end subroutine subevt_set_outgoing

  module subroutine subevt_set_composite (subevt, i, p, src)
    class(subevt_t), intent(inout) :: subevt
    integer, intent(in) :: i
    type(vector4_t), intent(in) :: p
    integer, dimension(:), intent(in) :: src
    call prt_init_composite (subevt%prt(i), p, src)
  end subroutine subevt_set_composite

  module subroutine subevt_set_pdg_beam (subevt, pdg)
    class(subevt_t), intent(inout) :: subevt
    integer, dimension(:), intent(in) :: pdg
    integer :: i, j
    j = 1
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_BEAM) then
          call prt_set_pdg (subevt%prt(i), pdg(j))
          j = j + 1
          if (j > size (pdg))  exit
       end if
    end do
  end subroutine subevt_set_pdg_beam

  module subroutine subevt_set_pdg_incoming (subevt, pdg)
    class(subevt_t), intent(inout) :: subevt
    integer, dimension(:), intent(in) :: pdg
    integer :: i, j
    j = 1
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_INCOMING) then
          call prt_set_pdg (subevt%prt(i), pdg(j))
          j = j + 1
          if (j > size (pdg))  exit
       end if
    end do
  end subroutine subevt_set_pdg_incoming

  module subroutine subevt_set_pdg_outgoing (subevt, pdg)
    class(subevt_t), intent(inout) :: subevt
    integer, dimension(:), intent(in) :: pdg
    integer :: i, j
    j = 1
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_OUTGOING) then
          call prt_set_pdg (subevt%prt(i), pdg(j))
          j = j + 1
          if (j > size (pdg))  exit
       end if
    end do
  end subroutine subevt_set_pdg_outgoing

  module subroutine subevt_set_p_beam (subevt, p)
    class(subevt_t), intent(inout) :: subevt
    type(vector4_t), dimension(:), intent(in) :: p
    integer :: i, j
    j = 1
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_BEAM) then
          call prt_set_p (subevt%prt(i), p(j))
          j = j + 1
          if (j > size (p))  exit
       end if
    end do
  end subroutine subevt_set_p_beam

  module subroutine subevt_set_p_incoming (subevt, p)
    class(subevt_t), intent(inout) :: subevt
    type(vector4_t), dimension(:), intent(in) :: p
    integer :: i, j
    j = 1
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_INCOMING) then
          call prt_set_p (subevt%prt(i), p(j))
          j = j + 1
          if (j > size (p))  exit
       end if
    end do
  end subroutine subevt_set_p_incoming

  module subroutine subevt_set_p_outgoing (subevt, p)
    class(subevt_t), intent(inout) :: subevt
    type(vector4_t), dimension(:), intent(in) :: p
    integer :: i, j
    j = 1
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_OUTGOING) then
          call prt_set_p (subevt%prt(i), p(j))
          j = j + 1
          if (j > size (p))  exit
       end if
    end do
  end subroutine subevt_set_p_outgoing

  module subroutine subevt_set_p2_beam (subevt, p2)
    class(subevt_t), intent(inout) :: subevt
    real(default), dimension(:), intent(in) :: p2
    integer :: i, j
    j = 1
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_BEAM) then
          call prt_set_p2 (subevt%prt(i), p2(j))
          j = j + 1
          if (j > size (p2))  exit
       end if
    end do
  end subroutine subevt_set_p2_beam

  module subroutine subevt_set_p2_incoming (subevt, p2)
    class(subevt_t), intent(inout) :: subevt
    real(default), dimension(:), intent(in) :: p2
    integer :: i, j
    j = 1
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_INCOMING) then
          call prt_set_p2 (subevt%prt(i), p2(j))
          j = j + 1
          if (j > size (p2))  exit
       end if
    end do
  end subroutine subevt_set_p2_incoming

  module subroutine subevt_set_p2_outgoing (subevt, p2)
    class(subevt_t), intent(inout) :: subevt
    real(default), dimension(:), intent(in) :: p2
    integer :: i, j
    j = 1
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_OUTGOING) then
          call prt_set_p2 (subevt%prt(i), p2(j))
          j = j + 1
          if (j > size (p2))  exit
       end if
    end do
  end subroutine subevt_set_p2_outgoing

  module subroutine subevt_polarize (subevt, i, h)
    type(subevt_t), intent(inout) :: subevt
    integer, intent(in) :: i, h
    call prt_polarize (subevt%prt(i), h)
  end subroutine subevt_polarize

  module subroutine subevt_colorize (subevt, i, col, acl)
    type(subevt_t), intent(inout) :: subevt
    integer, intent(in) :: i, col, acl
    if (col > 0 .and. acl > 0) then
       call prt_colorize (subevt%prt(i), [col], [acl])
    else if (col > 0) then
       call prt_colorize (subevt%prt(i), [col], [integer ::])
    else if (acl > 0) then
       call prt_colorize (subevt%prt(i), [integer ::], [acl])
    else
       call prt_colorize (subevt%prt(i), [integer ::], [integer ::])
    end if
  end subroutine subevt_colorize

  module function subevt_is_nonempty (subevt) result (flag)
    logical :: flag
    class(subevt_t), intent(in) :: subevt
    flag = subevt%n_active /= 0
  end function subevt_is_nonempty

  module function subevt_get_length (subevt) result (length)
    integer :: length
    class(subevt_t), intent(in) :: subevt
    length = subevt%n_active
  end function subevt_get_length

  module function subevt_get_prt (subevt, i) result (prt)
    type(prt_t) :: prt
    class(subevt_t), intent(in) :: subevt
    integer, intent(in) :: i
    prt = subevt%prt(i)
  end function subevt_get_prt

  module function subevt_get_sqrts_hat (subevt) result (sqrts_hat)
    class(subevt_t), intent(in) :: subevt
    real(default) :: sqrts_hat
    type(vector4_t) :: p
    integer :: i
    do i = 1, subevt%n_active
       if (subevt%prt(i)%type == PRT_INCOMING) then
          p = p + prt_get_momentum (subevt%prt(i))
       end if
    end do
    sqrts_hat = p ** 1
  end function subevt_get_sqrts_hat

  module function subevt_get_n_in (subevt) result (n_in)
    class(subevt_t), intent(in) :: subevt
    integer :: n_in
    n_in = count (subevt%prt(:subevt%n_active)%type == PRT_INCOMING)
  end function subevt_get_n_in

  module function subevt_get_n_out (subevt) result (n_out)
    class(subevt_t), intent(in) :: subevt
    integer :: n_out
    n_out = count (subevt%prt(:subevt%n_active)%type == PRT_OUTGOING)
  end function subevt_get_n_out

  module function c_prt_from_subevt (subevt, i) result (c_prt)
    type(c_prt_t) :: c_prt
    type(subevt_t), intent(in) :: subevt
    integer, intent(in) :: i
    c_prt = c_prt_from_prt (subevt%prt(i))
  end function c_prt_from_subevt

  module function c_prt_array_from_subevt (subevt) result (c_prt_array)
    type(subevt_t), intent(in) :: subevt
    type(c_prt_t), dimension(subevt%n_active) :: c_prt_array
    c_prt_array = c_prt_from_prt (subevt%prt(1:subevt%n_active))
  end function c_prt_array_from_subevt

  module subroutine subevt_join (subevt, pl1, pl2, mask2)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl1, pl2
    logical, dimension(:), intent(in), optional :: mask2
    integer :: n1, n2, i, n
    n1 = pl1%n_active
    n2 = pl2%n_active
    call subevt%reset (n1 + n2)
    subevt%prt(:n1)   = pl1%prt(:n1)
    n = n1
    if (present (mask2)) then
       do i = 1, pl2%n_active
          if (mask2(i)) then
             if (disjoint (i)) then
                n = n + 1
                subevt%prt(n) = pl2%prt(i)
             end if
          end if
       end do
    else
       do i = 1, pl2%n_active
          if (disjoint (i)) then
             n = n + 1
             subevt%prt(n) = pl2%prt(i)
          end if
       end do
    end if
    subevt%n_active = n
  contains
    function disjoint (i) result (flag)
      integer, intent(in) :: i
      logical :: flag
      integer :: j
      do j = 1, pl1%n_active
         if (.not. are_disjoint (pl1%prt(j), pl2%prt(i))) then
            flag = .false.
            return
         end if
      end do
      flag = .true.
    end function disjoint
  end subroutine subevt_join

  module subroutine subevt_combine (subevt, pl1, pl2, mask12)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl1, pl2
    logical, dimension(:,:), intent(in), optional :: mask12
    integer :: n1, n2, i1, i2, n, j
    logical :: ok
    n1 = pl1%n_active
    n2 = pl2%n_active
    call subevt%reset (n1 * n2)
    n = 1
    do i1 = 1, n1
       do i2 = 1, n2
          if (present (mask12)) then
             ok = mask12(i1,i2)
          else
             ok = .true.
          end if
          if (ok)  call prt_combine &
               (subevt%prt(n), pl1%prt(i1), pl2%prt(i2), ok)
          if (ok) then
             CHECK_DOUBLES: do j = 1, n - 1
                if (subevt%prt(n) .match. subevt%prt(j)) then
                   ok = .false.;  exit CHECK_DOUBLES
                end if
             end do CHECK_DOUBLES
             if (ok)  n = n + 1
          end if
       end do
    end do
    subevt%n_active = n - 1
  end subroutine subevt_combine

  module subroutine subevt_collect (subevt, pl1, mask1)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl1
    logical, dimension(:), intent(in) :: mask1
    type(prt_t) :: prt
    integer :: i
    logical :: ok
    call subevt%reset (1)
    subevt%n_active = 0
    do i = 1, pl1%n_active
       if (mask1(i)) then
          if (subevt%n_active == 0) then
             subevt%n_active = 1
             subevt%prt(1) = pl1%prt(i)
          else
             call prt_combine (prt, subevt%prt(1), pl1%prt(i), ok)
             if (ok)  subevt%prt(1) = prt
          end if
       end if
    end do
  end subroutine subevt_collect

  module subroutine subevt_cluster (subevt, pl1, dcut, mask1, jet_def, &
       keep_jets, exclusive)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl1
    real(default), intent(in) :: dcut
    logical, dimension(:), intent(in) :: mask1
    type(jet_definition_t), intent(in) :: jet_def
    logical, intent(in) :: keep_jets, exclusive
    integer, dimension(:), allocatable :: map, jet_index
    type(pseudojet_t), dimension(:), allocatable :: jet_in, jet_out
    type(pseudojet_vector_t) :: jv_in, jv_out
    type(cluster_sequence_t) :: cs
    integer :: i, n_src, n_active
    call map_prt_index (pl1, mask1, n_src, map)
    n_active = count (map /= 0)
    allocate (jet_in (n_active))
    allocate (jet_index (n_active))
    do i = 1, n_active
       call jet_in(i)%init (prt_get_momentum (pl1%prt(map(i))))
    end do
    call jv_in%init (jet_in)
    call cs%init (jv_in, jet_def)
    if (exclusive) then
       jv_out = cs%exclusive_jets (dcut)
    else
       jv_out = cs%inclusive_jets ()
    end if
    call cs%assign_jet_indices (jv_out, jet_index)
    allocate (jet_out (jv_out%size ()))
    jet_out = jv_out
    call fill_pseudojet (subevt, pl1, jet_out, jet_index, n_src, map)
    do i = 1, size (jet_out)
       call jet_out(i)%final ()
    end do
    call jv_out%final ()
    call cs%final ()
    call jv_in%final ()
    do i = 1, size (jet_in)
       call jet_in(i)%final ()
    end do
  contains
    ! Uniquely combine sources and add map those new indices to the old ones
    subroutine map_prt_index (pl1, mask1, n_src, map)
      type(subevt_t), intent(in) :: pl1
      logical, dimension(:), intent(in) :: mask1
      integer, intent(out) :: n_src
      integer, dimension(:), allocatable, intent(out) :: map
      integer, dimension(:), allocatable :: src, src_tmp
      integer :: i
      allocate (src(0))
      allocate (map (pl1%n_active), source = 0)
      n_active = 0
      do i = 1, pl1%n_active
         if (.not. mask1(i)) cycle
         call combine_index_lists (src_tmp, src, pl1%prt(i)%src)
         if (.not. allocated (src_tmp)) cycle
         call move_alloc (from=src_tmp, to=src)
         n_active = n_active + 1
         map(n_active) = i
      end do
      n_src = size (src)
    end subroutine map_prt_index
    ! Retrieve source(s) of a jet and fill corresponding subevent
    subroutine fill_pseudojet (subevt, pl1, jet_out, jet_index, n_src, map)
      type(subevt_t), intent(inout) :: subevt
      type(subevt_t), intent(in) :: pl1
      type(pseudojet_t), dimension(:), intent(in) :: jet_out
      integer, dimension(:), intent(in) :: jet_index
      integer, dimension(:), intent(in) :: map
      integer, intent(in) :: n_src
      integer, dimension(n_src) :: src_fill
      integer :: i, jet, k, combined_pdg, pdg, n_quarks, n_src_fill
      logical :: is_b, is_c
      call subevt%reset (size (jet_out))
      do jet = 1, size (jet_out)
         pdg = 0; src_fill = 0; n_src_fill = 0; combined_pdg = 0; n_quarks = 0
         is_b = .false.; is_c = .false.
         PARTICLE: do i = 1, size (jet_index)
            if (jet_index(i) /= jet) cycle PARTICLE
            associate (prt => pl1%prt(map(i)), n_src_prt => size(pl1%prt(map(i))%src))
              do k = 1, n_src_prt
                 src_fill(n_src_fill + k) = prt%src(k)
              end do
              n_src_fill = n_src_fill + n_src_prt
              if (is_quark (prt%pdg)) then
                 n_quarks = n_quarks + 1
                 if (.not. is_b) then
                    if (abs (prt%pdg) == 5) then
                       is_b = .true.
                       is_c = .false.
                    else if (abs (prt%pdg) == 4) then
                       is_c = .true.
                    end if
                 end if
                 if (combined_pdg == 0) combined_pdg = prt%pdg
              end if
            end associate
         end do PARTICLE
         if (keep_jets .and. n_quarks == 1) pdg = combined_pdg
         call prt_init_pseudojet (subevt%prt(jet), jet_out(jet), &
              src_fill(:n_src_fill), pdg, is_b, is_c)
      end do
    end subroutine fill_pseudojet
  end subroutine subevt_cluster

  module subroutine subevt_recombine (subevt, pl, mask1, reco_r0, keep_flv)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl
    type(prt_t), dimension(:), allocatable :: prt_rec
    logical, dimension(:), intent(in) :: mask1
    logical, intent(in) :: keep_flv
    real(default), intent(in) :: reco_r0
    real(default), dimension(:), allocatable :: del_rij
    integer, dimension(:), allocatable :: i_sortr
    type(prt_t) :: prt_gam, prt_comb
    logical :: recombine, ok
    integer :: i, n, i_gam, n_gam, n_rec, pdg_orig
    n = pl%get_length ()
    n_gam = 0
    FIND_FIRST_PHOTON: do i = 1, n
       if (prt_is_photon (pl%prt (i))) then
          n_gam = n_gam + 1
          prt_gam = pl%prt (i)
          i_gam = i
          exit FIND_FIRST_PHOTON
       end if
    end do FIND_FIRST_PHOTON
    n_rec = n - n_gam
    if (n_gam == 0) then
       subevt = pl
    else
       if (n_rec > 0) then
          allocate (prt_rec (n_rec))
          do i = 1, n_rec
             if (i == i_gam)  cycle
             if (i < i_gam) then
                prt_rec(i) = pl%prt(i)
             else
                prt_rec(i) = pl%prt(i+n_gam)
             end if
          end do
          allocate (del_rij (n_rec), i_sortr (n_rec))
          del_rij(1:n_rec) = eta_phi_distance(prt_get_momentum (prt_gam), &
               prt_get_momentum (prt_rec(1:n_rec)))
          i_sortr = order (del_rij)
          recombine = del_rij (i_sortr (1)) <= reco_r0 .and. mask1(i_gam)
          if (recombine) then
             call subevt%reset (pl%n_active-n_gam)
             do i = 1, n_rec
                if (i == i_sortr(1)) then
                   pdg_orig = prt_get_pdg (prt_rec(i_sortr (1)))
                   call prt_combine (prt_comb, prt_gam, prt_rec(i_sortr (1)), ok)
                   if (ok) then
                      subevt%prt(i_sortr (1)) = prt_comb
                      if (keep_flv)  call prt_set_pdg &
                           (subevt%prt(i_sortr (1)), pdg_orig)
                   end if
                else
                   subevt%prt(i) = prt_rec(i)
                end if
             end do
          else
             subevt = pl
          end if
       else
          subevt = pl
       end if
    end if
  end subroutine subevt_recombine

  module subroutine subevt_select (subevt, pl, mask1)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl
    logical, dimension(:), intent(in) :: mask1
    integer :: i, n
    call subevt%reset (pl%n_active)
    n = 0
    do i = 1, pl%n_active
       if (mask1(i)) then
          n = n + 1
          subevt%prt(n) = pl%prt(i)
       end if
    end do
    subevt%n_active = n
  end subroutine subevt_select

  module subroutine subevt_extract (subevt, pl, index)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl
    integer, intent(in) :: index
    if (index > 0) then
       if (index <= pl%n_active) then
          call subevt%reset (1)
          subevt%prt(1) = pl%prt(index)
       else
          call subevt%reset (0)
       end if
    else if (index < 0) then
       if (abs (index) <= pl%n_active) then
          call subevt%reset (1)
          subevt%prt(1) = pl%prt(pl%n_active + 1 + index)
       else
          call subevt%reset (0)
       end if
    else
       call subevt%reset (0)
    end if
  end subroutine subevt_extract

  module subroutine subevt_sort_pdg (subevt, pl)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl
    integer :: n
    n = subevt%n_active
    call subevt_sort_int (subevt, pl, abs (3 * subevt%prt(:n)%pdg - 1))
  end subroutine subevt_sort_pdg

  module subroutine subevt_sort_int (subevt, pl, ival)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl
    integer, dimension(:), intent(in) :: ival
    call subevt%reset (pl%n_active)
    subevt%n_active = pl%n_active
    subevt%prt = pl%prt( order (ival) )
  end subroutine subevt_sort_int

  module subroutine subevt_sort_real (subevt, pl, rval)
    type(subevt_t), intent(inout) :: subevt
    type(subevt_t), intent(in) :: pl
    real(default), dimension(:), intent(in) :: rval
    integer :: i
    integer, dimension(size(rval)) :: idx
    call subevt%reset (pl%n_active)
    subevt%n_active = pl%n_active
    if (allocated (subevt%prt))  deallocate (subevt%prt)
    allocate (subevt%prt (size(pl%prt)))
    idx = order (rval)
    do i = 1, size (idx)
       subevt%prt(i) = pl%prt (idx(i))
    end do
  end subroutine subevt_sort_real

  module subroutine subevt_select_pdg_code (subevt, aval, subevt_in, prt_type)
    type(subevt_t), intent(inout) :: subevt
    type(pdg_array_t), intent(in) :: aval
    type(subevt_t), intent(in) :: subevt_in
    integer, intent(in), optional :: prt_type
    integer :: n_active, n_match
    logical, dimension(:), allocatable :: mask
    integer :: i, j
    n_active = subevt_in%n_active
    allocate (mask (n_active))
    forall (i = 1:n_active) &
         mask(i) = aval .match. subevt_in%prt(i)%pdg
    if (present (prt_type)) &
         mask = mask .and. subevt_in%prt(:n_active)%type == prt_type
    n_match = count (mask)
    call subevt%reset (n_match)
    j = 0
    do i = 1, n_active
       if (mask(i)) then
          j = j + 1
          subevt%prt(j) = subevt_in%prt(i)
       end if
    end do
  end subroutine subevt_select_pdg_code

  module subroutine pacify_prt (prt)
    class(prt_t), intent(inout) :: prt
    real(default) :: e
    e = max (1E-10_default * energy (prt%p), 1E-13_default)
    call pacify (prt%p, e)
    call pacify (prt%p2, 1E3_default * e)
  end subroutine pacify_prt

  module subroutine pacify_subevt (subevt)
    class(subevt_t), intent(inout) :: subevt
    integer :: i
    do i = 1, subevt%n_active
       call pacify (subevt%prt(i))
    end do
  end subroutine pacify_subevt


end submodule subevents_s

