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

module beam_structures_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  use beam_structures

  implicit none
  private

  public :: beam_structures_1
  public :: beam_structures_2
  public :: beam_structures_3
  public :: beam_structures_4
  public :: beam_structures_5
  public :: beam_structures_6

contains

  subroutine beam_structures_1 (u)
    integer, intent(in) :: u
    type(beam_structure_t) :: beam_structure

    write (u, "(A)")  "* Test output: beam_structures_1"
    write (u, "(A)")  "*   Purpose: display empty beam structure record"
    write (u, "(A)")

    call beam_structure%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: beam_structures_1"

  end subroutine beam_structures_1

  subroutine beam_structures_2 (u)
    integer, intent(in) :: u
    type(beam_structure_t) :: beam_structure
    integer, dimension(0) :: empty_array
    type(string_t) :: s

    write (u, "(A)")  "* Test output: beam_structures_2"
    write (u, "(A)")  "*   Purpose: setup beam structure records"
    write (u, "(A)")

    s = "s"

    call beam_structure%init_sf ([s], empty_array)
    call beam_structure%write (u)

    write (u, "(A)")

    call beam_structure%init_sf ([s, s], [1])
    call beam_structure%set_sf (1, 1, var_str ("a"))
    call beam_structure%write (u)

    write (u, "(A)")

    call beam_structure%init_sf ([s, s], [2])
    call beam_structure%set_sf (1, 1, var_str ("a"))
    call beam_structure%set_sf (1, 2, var_str ("b"))
    call beam_structure%write (u)

    write (u, "(A)")

    call beam_structure%init_sf ([s, s], [2, 1])
    call beam_structure%set_sf (1, 1, var_str ("a"))
    call beam_structure%set_sf (1, 2, var_str ("b"))
    call beam_structure%set_sf (2, 1, var_str ("c"))
    call beam_structure%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: beam_structures_2"

  end subroutine beam_structures_2

  subroutine beam_structures_3 (u)
    integer, intent(in) :: u
    type(beam_structure_t) :: beam_structure
    type(string_t) :: s

    write (u, "(A)")  "* Test output: beam_structures_3"
    write (u, "(A)")  "*   Purpose: expand beam structure records"
    write (u, "(A)")

    s = "s"

    write (u, "(A)")  "* Pair spectrum (keep as-is)"
    write (u, "(A)")

    call beam_structure%init_sf ([s, s], [1])
    call beam_structure%set_sf (1, 1, var_str ("a"))
    call beam_structure%write (u)

    write (u, "(A)")

    call beam_structure%expand (test_strfun_mode)
    call beam_structure%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Structure function pair (expand)"
    write (u, "(A)")

    call beam_structure%init_sf ([s, s], [2])
    call beam_structure%set_sf (1, 1, var_str ("b"))
    call beam_structure%set_sf (1, 2, var_str ("b"))
    call beam_structure%write (u)

    write (u, "(A)")

    call beam_structure%expand (test_strfun_mode)
    call beam_structure%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Structure function (separate and expand)"
    write (u, "(A)")

    call beam_structure%init_sf ([s, s], [1])
    call beam_structure%set_sf (1, 1, var_str ("b"))
    call beam_structure%write (u)

    write (u, "(A)")

    call beam_structure%expand (test_strfun_mode)
    call beam_structure%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Combination"
    write (u, "(A)")

    call beam_structure%init_sf ([s, s], [1, 1])
    call beam_structure%set_sf (1, 1, var_str ("a"))
    call beam_structure%set_sf (2, 1, var_str ("b"))
    call beam_structure%write (u)

    write (u, "(A)")

    call beam_structure%expand (test_strfun_mode)
    call beam_structure%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: beam_structures_3"

  end subroutine beam_structures_3

  subroutine beam_structures_4 (u)
    integer, intent(in) :: u
    type(beam_structure_t) :: beam_structure
    type(string_t) :: s
    type(string_t), dimension(2) :: prt
    integer :: i

    write (u, "(A)")  "* Test output: beam_structures_4"
    write (u, "(A)")  "*   Purpose: check the API"
    write (u, "(A)")

    s = "s"

    write (u, "(A)")  "* Structure-function combination"
    write (u, "(A)")

    call beam_structure%init_sf ([s, s], [1, 2, 2])
    call beam_structure%set_sf (1, 1, var_str ("a"))
    call beam_structure%set_sf (2, 1, var_str ("b"))
    call beam_structure%set_sf (3, 2, var_str ("c"))
    call beam_structure%write (u)

    write (u, *)
    write (u, "(1x,A,I0)")  "n_beam = ", beam_structure%get_n_beam ()
    prt = beam_structure%get_prt ()
    write (u, "(1x,A,2(1x,A))")  "prt =", char (prt(1)), char (prt(2))

    write (u, *)
    write (u, "(1x,A,I0)")  "n_record = ", beam_structure%get_n_record ()

    do i = 1, 3
       write (u, "(A)")
       write (u, "(1x,A,I0,A,A)")  "name(", i, ") = ", &
            char (beam_structure%get_name (i))
       write (u, "(1x,A,I0,A,2(1x,I0))")  "i_entry(", i, ") =", &
            beam_structure%get_i_entry (i)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: beam_structures_4"

  end subroutine beam_structures_4

  subroutine beam_structures_5 (u)
    integer, intent(in) :: u
    type(beam_structure_t) :: beam_structure
    integer, dimension(0) :: empty_array
    type(string_t) :: s

    write (u, "(A)")  "* Test output: beam_structures_5"
    write (u, "(A)")  "*   Purpose: setup polarization in beam structure records"
    write (u, "(A)")

    s = "s"

    call beam_structure%init_sf ([s], empty_array)
    call beam_structure%init_pol (1)
    call beam_structure%init_smatrix (1, 1)
    call beam_structure%set_sentry (1, 1, [0,0], (1._default, 0._default))
    call beam_structure%set_pol_f ([0.5_default])
    call beam_structure%write (u)


    write (u, "(A)")
    call beam_structure%final_sf ()
    call beam_structure%final_pol ()

    call beam_structure%init_sf ([s, s], [1])
    call beam_structure%set_sf (1, 1, var_str ("a"))
    call beam_structure%init_pol (2)
    call beam_structure%init_smatrix (1, 2)
    call beam_structure%set_sentry (1, 1, [-1,1], (0.5_default,-0.5_default))
    call beam_structure%set_sentry (1, 2, [ 1,1], (1._default, 0._default))
    call beam_structure%init_smatrix (2, 0)
    call beam_structure%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: beam_structures_5"

  end subroutine beam_structures_5

  subroutine beam_structures_6 (u)
    integer, intent(in) :: u
    type(beam_structure_t) :: beam_structure
    integer, dimension(0) :: empty_array
    type(string_t) :: s

    write (u, "(A)")  "* Test output: beam_structures_6"
    write (u, "(A)")  "*   Purpose: setup momenta in beam structure records"
    write (u, "(A)")

    s = "s"

    call beam_structure%init_sf ([s], empty_array)
    call beam_structure%set_momentum ([500._default])
    call beam_structure%write (u)


    write (u, "(A)")
    call beam_structure%final_sf ()
    call beam_structure%final_mom ()

    call beam_structure%init_sf ([s, s], [1])
    call beam_structure%set_momentum ([500._default, 700._default])
    call beam_structure%set_theta ([0._default, 0.1_default])
    call beam_structure%set_phi ([0._default, 1.51_default])
    call beam_structure%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: beam_structures_6"

  end subroutine beam_structures_6


  function test_strfun_mode (name) result (n)
    type(string_t), intent(in) :: name
    integer :: n
    select case (char (name))
    case ("a");  n = 2
    case ("b");  n = 1
    case default;  n = 0
    end select
  end function test_strfun_mode


end module beam_structures_uti
