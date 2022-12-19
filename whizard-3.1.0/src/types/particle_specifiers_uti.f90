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

module particle_specifiers_uti

  use iso_varying_string, string_t => varying_string

  use particle_specifiers

  implicit none
  private

  public :: particle_specifiers_1
  public :: particle_specifiers_2

contains

  subroutine particle_specifiers_1 (u)
    integer, intent(in) :: u
    type(prt_spec_t), dimension(:), allocatable :: prt_spec
    type(string_t), dimension(:), allocatable :: decay
    type(string_t), dimension(0) :: no_decay
    integer :: i, j

    write (u, "(A)")  "* Test output: particle_specifiers_1"
    write (u, "(A)")  "*   Purpose: Read and write a particle specifier array"
    write (u, "(A)")

    allocate (prt_spec (5))
    prt_spec = [ &
         new_prt_spec (var_str ("a")), &
         new_prt_spec (var_str ("b"), .true.), &
         new_prt_spec (var_str ("c"), [var_str ("dec1")]), &
         new_prt_spec (var_str ("d"), [var_str ("dec1"), var_str ("dec2")]), &
         new_prt_spec (var_str ("e"), no_decay) &
         ]
    do i = 1, size (prt_spec)
       write (u, "(A)")  char (prt_spec(i)%to_string ())
    end do
    write (u, "(A)")

    call prt_spec_read (prt_spec, &
         var_str (" a, b( *), c( dec1), d (dec1 + dec2 ), e()"))
    call prt_spec_write (prt_spec, u)

    do i = 1, size (prt_spec)
       write (u, "(A)")
       write (u, "(A,A)")  char (prt_spec(i)%get_name ()), ":"
       write (u, "(A,L1)")  "polarized = ", prt_spec(i)%is_polarized ()
       write (u, "(A,L1)")  "unstable  = ", prt_spec(i)%is_unstable ()
       write (u, "(A,I0)")  "n_decays  = ", prt_spec(i)%get_n_decays ()
       call prt_spec(i)%get_decays (decay)
       write (u, "(A)", advance="no") "decays    ="
       do j = 1, size (decay)
          write (u, "(1x,A)", advance="no") char (decay(j))
       end do
       write (u, "(A)")
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particle_specifiers_1"
  end subroutine particle_specifiers_1

  subroutine particle_specifiers_2 (u)
    integer, intent(in) :: u
    type(prt_spec_t) :: a, b, c, d, e, f
    type(prt_expr_t) :: pe1, pe2, pe3
    type(prt_expr_t) :: pe4, pe5, pe6, pe7, pe8, pe9
    integer :: i
    type(prt_spec_t), dimension(:), allocatable :: pa

    write (u, "(A)")  "* Test output: particle_specifiers_2"
    write (u, "(A)")  "*   Purpose: Create and display particle expressions"
    write (u, "(A)")

    write (u, "(A)")  "* Basic expressions"
    write (u, *)

    a = new_prt_spec (var_str ("a"))
    b = new_prt_spec (var_str ("b"))
    c = new_prt_spec (var_str ("c"))
    d = new_prt_spec (var_str ("d"))
    e = new_prt_spec (var_str ("e"))
    f = new_prt_spec (var_str ("f"))

    call pe1%init_spec (a)
    write (u, "(A)")  char (pe1%to_string ())

    call pe2%init_sum (2)
    select type (x => pe2%x)
    type is (prt_spec_sum_t)
       call x%expr(1)%init_spec (a)
       call x%expr(2)%init_spec (b)
    end select
    write (u, "(A)")  char (pe2%to_string ())

    call pe3%init_list (2)
    select type (x => pe3%x)
    type is (prt_spec_list_t)
       call x%expr(1)%init_spec (a)
       call x%expr(2)%init_spec (b)
    end select
    write (u, "(A)")  char (pe3%to_string ())

    write (u, *)
    write (u, "(A)")  "* Nested expressions"
    write (u, *)

    call pe4%init_list (2)
    select type (x => pe4%x)
    type is (prt_spec_list_t)
       call x%expr(1)%init_sum (2)
       select type (y => x%expr(1)%x)
       type is (prt_spec_sum_t)
          call y%expr(1)%init_spec (a)
          call y%expr(2)%init_spec (b)
       end select
       call x%expr(2)%init_spec (c)
    end select
    write (u, "(A)")  char (pe4%to_string ())

    call pe5%init_list (2)
    select type (x => pe5%x)
    type is (prt_spec_list_t)
       call x%expr(1)%init_list (2)
       select type (y => x%expr(1)%x)
       type is (prt_spec_list_t)
          call y%expr(1)%init_spec (a)
          call y%expr(2)%init_spec (b)
       end select
       call x%expr(2)%init_spec (c)
    end select
    write (u, "(A)")  char (pe5%to_string ())

    call pe6%init_sum (2)
    select type (x => pe6%x)
    type is (prt_spec_sum_t)
       call x%expr(1)%init_spec (a)
       call x%expr(2)%init_sum (2)
       select type (y => x%expr(2)%x)
       type is (prt_spec_sum_t)
          call y%expr(1)%init_spec (b)
          call y%expr(2)%init_spec (c)
       end select
    end select
    write (u, "(A)")  char (pe6%to_string ())

    call pe7%init_list (2)
    select type (x => pe7%x)
    type is (prt_spec_list_t)
       call x%expr(1)%init_sum (2)
       select type (y => x%expr(1)%x)
       type is (prt_spec_sum_t)
          call y%expr(1)%init_spec (a)
          call y%expr(2)%init_list (2)
          select type (z => y%expr(2)%x)
          type is (prt_spec_list_t)
             call z%expr(1)%init_spec (b)
             call z%expr(2)%init_spec (c)
          end select
       end select
       call x%expr(2)%init_spec (d)
    end select
    write (u, "(A)")  char (pe7%to_string ())

    call pe8%init_sum (2)
    select type (x => pe8%x)
    type is (prt_spec_sum_t)
       call x%expr(1)%init_list (2)
       select type (y => x%expr(1)%x)
       type is (prt_spec_list_t)
          call y%expr(1)%init_spec (a)
          call y%expr(2)%init_spec (b)
       end select
       call x%expr(2)%init_list (2)
       select type (y => x%expr(2)%x)
       type is (prt_spec_list_t)
          call y%expr(1)%init_spec (c)
          call y%expr(2)%init_spec (d)
       end select
    end select
    write (u, "(A)")  char (pe8%to_string ())

    call pe9%init_list (3)
    select type (x => pe9%x)
    type is (prt_spec_list_t)
       call x%expr(1)%init_sum (2)
       select type (y => x%expr(1)%x)
       type is (prt_spec_sum_t)
          call y%expr(1)%init_spec (a)
          call y%expr(2)%init_spec (b)
       end select
       call x%expr(2)%init_spec (c)
       call x%expr(3)%init_sum (3)
       select type (y => x%expr(3)%x)
       type is (prt_spec_sum_t)
          call y%expr(1)%init_spec (d)
          call y%expr(2)%init_spec (e)
          call y%expr(3)%init_spec (f)
       end select
    end select
    write (u, "(A)")  char (pe9%to_string ())

    write (u, *)
    write (u, "(A)")  "* Expand as sum"
    write (u, *)

    call pe1%expand ()
    write (u, "(A)")  char (pe1%to_string ())

    call pe4%expand ()
    write (u, "(A)")  char (pe4%to_string ())

    call pe5%expand ()
    write (u, "(A)")  char (pe5%to_string ())

    call pe6%expand ()
    write (u, "(A)")  char (pe6%to_string ())

    call pe7%expand ()
    write (u, "(A)")  char (pe7%to_string ())

    call pe8%expand ()
    write (u, "(A)")  char (pe8%to_string ())

    call pe9%expand ()
    write (u, "(A)")  char (pe9%to_string ())

    write (u, *)
    write (u, "(A)")  "* Transform to arrays:"

    write (u, "(A)")  "* Atomic specifier"
    do i = 1, pe1%get_n_terms ()
       call pe1%term_to_array (pa, i)
       call prt_spec_write (pa, u)
    end do

    write (u, *)
    write (u, "(A)")  "* List"
    do i = 1, pe5%get_n_terms ()
       call pe5%term_to_array (pa, i)
       call prt_spec_write (pa, u)
    end do

    write (u, *)
    write (u, "(A)")  "* Sum of atoms"
    do i = 1, pe6%get_n_terms ()
       call pe6%term_to_array (pa, i)
       call prt_spec_write (pa, u)
    end do

    write (u, *)
    write (u, "(A)")  "* Sum of lists"
    do i = 1, pe9%get_n_terms ()
       call pe9%term_to_array (pa, i)
       call prt_spec_write (pa, u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: particle_specifiers_2"
  end subroutine particle_specifiers_2


end module particle_specifiers_uti
