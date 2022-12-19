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

module sm_qed_uti

  use kinds, only: default
  use physics_defs, only: ME_REF

  use sm_qed

  implicit none
  private

  public :: sm_qed_1

contains

  subroutine sm_qed_1 (u)
    integer, intent(in) :: u
    type(qed_t) :: qed

    write (u, "(A)")  "* Test output: sm_qed_1"
    write (u, "(A)")  "*   Purpose: compute running alpha"
    write (u, "(A)")

    write (u, "(A)")  "* Fixed:"
    write (u, "(A)")

    allocate (alpha_qed_fixed_t :: qed%alpha)
    call qed%compute_alpha_md5sum ()

    call qed%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha (me)     =", &
         qed%alpha%get (ME_REF)
    write (u, "(1x,A,F10.7)")  "alpha (10 GeV) =", &
         qed%alpha%get (10._default)
    write (u, "(1x,A,F10.7)")  "alpha (1 TeV)  =", &
         qed%alpha%get (1000._default)
    write (u, *)
    deallocate (qed%alpha)

    write (u, "(A)")  "* Running from me (LO):"
    write (u, "(A)")

    allocate (alpha_qed_from_scale_t :: qed%alpha)
    call qed%compute_alpha_md5sum ()

    call qed%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha (me)     =", &
         qed%alpha%get (ME_REF)
    write (u, "(1x,A,F10.7)")  "alpha (10 GeV) =", &
         qed%alpha%get (10._default)
    write (u, "(1x,A,F10.7)")  "alpha (1 TeV)  =", &
         qed%alpha%get (1000._default)
    write (u, *)

    write (u, "(A)")  "* Running from me (NLO, analytic):"
    write (u, "(A)")

    select type (alpha => qed%alpha)
    type is (alpha_qed_from_scale_t)
       alpha%order = 1
    end select
    call qed%compute_alpha_md5sum ()

    call qed%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha (me)     =", &
         qed%alpha%get (ME_REF)
    write (u, "(1x,A,F10.7)")  "alpha (10 GeV) =", &
         qed%alpha%get (10._default)
    write (u, "(1x,A,F10.7)")  "alpha (1 TeV)  =", &
         qed%alpha%get (1000._default)
    write (u, *)

    write (u, "(A)")  "* Running from me (NLO, numeric):"
    write (u, "(A)")

    select type (alpha => qed%alpha)
    type is (alpha_qed_from_scale_t)
       alpha%order = 1
       alpha%analytic = .false.
    end select
    call qed%compute_alpha_md5sum ()

    call qed%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha (me)     =", &
         qed%alpha%get (ME_REF)
    write (u, "(1x,A,F10.7)")  "alpha (10 GeV) =", &
         qed%alpha%get (10._default)
    write (u, "(1x,A,F10.7)")  "alpha (1 TeV)  =", &
         qed%alpha%get (1000._default)
    write (u, *)
    deallocate (qed%alpha)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sm_qed_1"

  end subroutine sm_qed_1


end module sm_qed_uti
