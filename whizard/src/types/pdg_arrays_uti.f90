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

module pdg_arrays_uti

  use pdg_arrays

  implicit none
  private

  public :: pdg_arrays_1
  public :: pdg_arrays_2
  public :: pdg_arrays_3
  public :: pdg_arrays_4
  public :: pdg_arrays_5

contains

  subroutine pdg_arrays_1 (u)
    integer, intent(in) :: u

    type(pdg_array_t) :: pa, pa1, pa2, pa3, pa4, pa5, pa6
    integer, dimension(:), allocatable :: pdg

    write (u, "(A)")  "* Test output: pdg_arrays_1"
    write (u, "(A)")  "*   Purpose: create and sort PDG arrays"
    write (u, "(A)")

    write (u, "(A)")  "* Assignment"
    write (u, "(A)")

    call pa%write (u)
    write (u, *)
    write (u, "(A,I0)")  "length = ", pa%get_length ()
    pdg = pa
    write (u, "(A,3(1x,I0))")  "contents = ", pdg

    write (u, *)
    pa = 1
    call pa%write (u)
    write (u, *)
    write (u, "(A,I0)")  "length = ", pa%get_length ()
    pdg = pa
    write (u, "(A,3(1x,I0))")  "contents = ", pdg

    write (u, *)
    pa = [1, 2, 3]
    call pa%write (u)
    write (u, *)
    write (u, "(A,I0)")  "length = ", pa%get_length ()
    pdg = pa
    write (u, "(A,3(1x,I0))")  "contents = ", pdg
    write (u, "(A,I0)")  "element #2 = ", pa%get (2)

    write (u, *)
    write (u, "(A)")  "* Replace"
    write (u, *)

    pa = pa%replace (2, [-5, 5, -7])
    call pa%write (u)
    write (u, *)

    write (u, *)
    write (u, "(A)")  "* Sort"
    write (u, *)

    pa = [1, -7, 3, -5, 5, 3]
    call pa%write (u)
    write (u, *)
    pa1 = pa%sort_abs ()
    pa2 = pa%sort_abs (unique = .true.)
    call pa1%write (u)
    write (u, *)
    call pa2%write (u)
    write (u, *)

    write (u, *)
    write (u, "(A)")  "* Compare"
    write (u, *)

    pa1 = [1, 3]
    pa2 = [1, 2, -2]
    pa3 = [1, 2, 4]
    pa4 = [1, 2, 4]
    pa5 = [1, 2, -4]
    pa6 = [1, 2, -3]

    write (u, "(A,6(1x,L1))")  "< ", &
         pa1 < pa2, pa2 < pa3, pa3 < pa4, pa4 < pa5, pa5 < pa6, pa6 < pa1
    write (u, "(A,6(1x,L1))")  "> ", &
         pa1 > pa2, pa2 > pa3, pa3 > pa4, pa4 > pa5, pa5 > pa6, pa6 > pa1
    write (u, "(A,6(1x,L1))")  "<=", &
         pa1 <= pa2, pa2 <= pa3, pa3 <= pa4, pa4 <= pa5, pa5 <= pa6, pa6 <= pa1
    write (u, "(A,6(1x,L1))")  ">=", &
         pa1 >= pa2, pa2 >= pa3, pa3 >= pa4, pa4 >= pa5, pa5 >= pa6, pa6 >= pa1
    write (u, "(A,6(1x,L1))")  "==", &
         pa1 == pa2, pa2 == pa3, pa3 == pa4, pa4 == pa5, pa5 == pa6, pa6 == pa1
    write (u, "(A,6(1x,L1))")  "/=", &
         pa1 /= pa2, pa2 /= pa3, pa3 /= pa4, pa4 /= pa5, pa5 /= pa6, pa6 /= pa1

    write (u, *)
    pa1 = [0]
    pa2 = [1, 2]
    pa3 = [1, -2]

    write (u, "(A,6(1x,L1))")  "eqv ", &
         pa1 .eqv. pa1, pa1 .eqv. pa2, &
         pa2 .eqv. pa2, pa2 .eqv. pa3

    write (u, "(A,6(1x,L1))")  "neqv", &
         pa1 .neqv. pa1, pa1 .neqv. pa2, &
         pa2 .neqv. pa2, pa2 .neqv. pa3


    write (u, *)
    write (u, "(A,6(1x,L1))")  "match", &
         pa1 .match. 0, pa1 .match. 1, &
         pa2 .match. 0, pa2 .match. 1, pa2 .match. 3

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: pdg_arrays_1"

  end subroutine pdg_arrays_1

  subroutine pdg_arrays_2 (u)
    integer, intent(in) :: u

    type(pdg_array_t) :: pa
    type(pdg_list_t) :: pl, pl1

    write (u, "(A)")  "* Test output: pdg_arrays_2"
    write (u, "(A)")  "*   Purpose: create and sort PDG lists"
    write (u, "(A)")

    write (u, "(A)")  "* Assignment"
    write (u, "(A)")

    call pl%init (3)
    call pl%set (1, 42)
    call pl%set (2, [3, 2])
    pa = [5, -5]
    call pl%set (3, pa)
    call pl%write (u)
    write (u, *)
    write (u, "(A,I0)")  "size = ", pl%get_size ()

    write (u, "(A)")
    write (u, "(A)")  "* Sort"
    write (u, "(A)")

    pl = pl%sort_abs ()
    call pl%write (u)
    write (u, *)

    write (u, "(A)")
    write (u, "(A)")  "* Extract item #3"
    write (u, "(A)")

    pa = pl%get (3)
    call pa%write (u)
    write (u, *)

    write (u, "(A)")
    write (u, "(A)")  "* Replace item #3"
    write (u, "(A)")

    call pl1%init (2)
    call pl1%set (1, [2, 4])
    call pl1%set (2, -7)

    pl = pl%replace (3, pl1)
    call pl%write (u)
    write (u, *)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: pdg_arrays_2"

  end subroutine pdg_arrays_2

  subroutine pdg_arrays_3 (u)
    integer, intent(in) :: u

    type(pdg_list_t) :: pl

    write (u, "(A)")  "* Test output: pdg_arrays_3"
    write (u, "(A)")  "*   Purpose: check for regular PDG lists"
    write (u, "(A)")

    write (u, "(A)")  "* Regular list"
    write (u, "(A)")

    call pl%init (4)
    call pl%set (1, [1, 2])
    call pl%set (2, [1, 2])
    call pl%set (3, [5, -5])
    call pl%set (4, 42)
    call pl%write (u)
    write (u, *)
    write (u, "(L1)") pl%is_regular ()

    write (u, "(A)")
    write (u, "(A)")  "* Irregular list"
    write (u, "(A)")

    call pl%init (4)
    call pl%set (1, [1, 2])
    call pl%set (2, [1, 2])
    call pl%set (3, [2, 5, -5])
    call pl%set (4, 42)
    call pl%write (u)
    write (u, *)
    write (u, "(L1)") pl%is_regular ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: pdg_arrays_3"

  end subroutine pdg_arrays_3

  subroutine pdg_arrays_4 (u)
    integer, intent(in) :: u

    type(pdg_list_t) :: pl1, pl2, pl3

    write (u, "(A)")  "* Test output: pdg_arrays_4"
    write (u, "(A)")  "*   Purpose: check for regular PDG lists"
    write (u, "(A)")

    write (u, "(A)")  "* Create lists"
    write (u, "(A)")

    call pl1%init (4)
    call pl1%set (1, [1, 2])
    call pl1%set (2, [1, 2])
    call pl1%set (3, [5, -5])
    call pl1%set (4, 42)
    write (u, "(I1,1x)", advance = "no")  1
    call pl1%write (u)
    write (u, *)

    call pl2%init (2)
    call pl2%set (1, 3)
    call pl2%set (2, [5, -5])
    write (u, "(I1,1x)", advance = "no")  2
    call pl2%write (u)
    write (u, *)

    call pl3%init (2)
    call pl3%set (1, 4)
    call pl3%set (2, [5, -5])
    write (u, "(I1,1x)", advance = "no")  3
    call pl3%write (u)
    write (u, *)

    write (u, "(A)")
    write (u, "(A)")  "* a == b"
    write (u, "(A)")

    write (u, "(2x,A)")  "123"
    write (u, *)
    write (u, "(I1,1x,4L1)")  1, pl1 == pl1, pl1 == pl2, pl1 == pl3
    write (u, "(I1,1x,4L1)")  2, pl2 == pl1, pl2 == pl2, pl2 == pl3
    write (u, "(I1,1x,4L1)")  3, pl3 == pl1, pl3 == pl2, pl3 == pl3

    write (u, "(A)")
    write (u, "(A)")  "* a < b"
    write (u, "(A)")

    write (u, "(2x,A)")  "123"
    write (u, *)
    write (u, "(I1,1x,4L1)")  1, pl1 < pl1, pl1 < pl2, pl1 < pl3
    write (u, "(I1,1x,4L1)")  2, pl2 < pl1, pl2 < pl2, pl2 < pl3
    write (u, "(I1,1x,4L1)")  3, pl3 < pl1, pl3 < pl2, pl3 < pl3

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: pdg_arrays_4"

  end subroutine pdg_arrays_4

  subroutine pdg_arrays_5 (u)
    integer, intent(in) :: u

    type(pdg_list_t) :: pl1, pl2, pl3
    logical :: success

    write (u, "(A)")  "* Test output: pdg_arrays_5"
    write (u, "(A)")  "*   Purpose: match-replace"
    write (u, "(A)")

    write (u, "(A)")  "* Create lists"
    write (u, "(A)")

    call pl1%init (2)
    call pl1%set (1, [1, 2])
    call pl1%set (2, 42)
    call pl1%write (u)
    write (u, *)
    call pl3%init (2)
    call pl3%set (1, [42, -42])
    call pl3%set (2, [1, 2, 3, 4])
    call pl1%match_replace (pl3, success)
    call pl3%write (u)
    write (u, "(1x,A,1x,L1,':',1x)", advance="no")  "=>", success
    call pl1%write (u)
    write (u, *)

    write (u, *)

    call pl2%init (2)
    call pl2%set (1, 9)
    call pl2%set (2, 42)
    call pl2%write (u)
    write (u, *)
    call pl2%match_replace (pl3, success)
    call pl3%write (u)
    write (u, "(1x,A,1x,L1,':',1x)", advance="no")  "=>", success
    call pl2%write (u)
    write (u, *)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: pdg_arrays_5"

  end subroutine pdg_arrays_5


end module pdg_arrays_uti
