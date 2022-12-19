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

module interactions_uti

  use kinds, only: default
  use lorentz
  use flavors
  use colors
  use helicities
  use quantum_numbers
  use state_matrices

  use interactions

  implicit none
  private

  public :: interaction_1

contains

  subroutine interaction_1 (u)
    integer, intent(in) :: u
    type(interaction_t), target :: int, rad
    type(vector4_t), dimension(3) :: p
    type(quantum_numbers_mask_t), dimension(3) :: mask
    p(2) = vector4_moving (500._default, 500._default, 1)
    p(3) = vector4_moving (500._default,-500._default, 1)
    p(1) = p(2) + p(3)

    write (u, "(A)")  "* Test output: interaction"
    write (u, "(A)")  "*   Purpose: check routines for interactions"
    write (u, "(A)")

    call int%basic_init (1, 0, 2, set_relations=.true., &
         store_values = .true. )
    call int_set (int, 1, -1, 1, 1, &
         cmplx (0.3_default, 0.1_default, kind=default))
    call int_set (int, 1, -1,-1, 1, &
         cmplx (0.5_default,-0.7_default, kind=default))
    call int_set (int, 1, 1, 1, 1, &
         cmplx (0.1_default, 0._default, kind=default))
    call int_set (int, -1, 1, -1, 2, &
         cmplx (0.4_default, -0.1_default, kind=default))
    call int_set (int, 1, 1, 1, 2, &
         cmplx (0.2_default, 0._default, kind=default))
    call int%freeze ()
    call int%set_momenta (p)
    mask = quantum_numbers_mask (.false.,.false., [.true.,.true.,.true.])
    call rad%basic_init (1, 0, 2, &
         mask=mask, set_relations=.true., store_values = .true.)
    call rad_set (1)
    call rad_set (2)
    call rad%set_source_link (1, int, 2)
    call rad%exchange_mask ()
    call rad%receive_momenta ()
    p(1) = rad%get_momentum (1)
    p(2) = 0.4_default * p(1)
    p(3) = p(1) - p(2)
    call rad%set_momenta (p(2:3), outgoing=.true.)
    call int%freeze ()
    call rad%freeze ()
    call rad%set_matrix_element &
         (cmplx (0._default, 0._default, kind=default))
    call int%basic_write (u)
    write (u, "(A)")
    call rad%basic_write (u)
    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"
    call int%final ()
    call rad%final ()
    write (u, "(A)")
    write (u, "(A)")  "* Test interaction_1: successful."
  contains
    subroutine int_set (int, h1, h2, hq, q, val)
      type(interaction_t), target, intent(inout) :: int
      integer, intent(in) :: h1, h2, hq, q
      type(flavor_t), dimension(3) :: flv
      type(color_t), dimension(3) :: col
      type(helicity_t), dimension(3) :: hel
      type(quantum_numbers_t), dimension(3) :: qn
      complex(default), intent(in) :: val
      call flv%init ([21, q, -q])
      call col(2)%init_col_acl (5, 0)
      call col(3)%init_col_acl (0, 5)
      call hel%init ([h1, hq, -hq], [h2, hq, -hq])
      call qn%init (flv, col, hel)
      call int%add_state (qn)
      call int%set_matrix_element (val)
    end subroutine int_set
    subroutine rad_set (q)
      integer, intent(in) :: q
      type(flavor_t), dimension(3) :: flv
      type(quantum_numbers_t), dimension(3) :: qn
      call flv%init ([ q, q, 21 ])
      call qn%init (flv)
      call rad%add_state (qn)
    end subroutine rad_set
  end subroutine interaction_1


end module interactions_uti
