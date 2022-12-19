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

module sm_qcd_uti

  use kinds, only: default
  use physics_defs, only: MZ_REF

  use sm_qcd

  implicit none
  private

  public :: sm_qcd_1

contains

  subroutine sm_qcd_1 (u)
    integer, intent(in) :: u
    type(qcd_t) :: qcd

    write (u, "(A)")  "* Test output: sm_qcd_1"
    write (u, "(A)")  "*   Purpose: compute running alpha_s"
    write (u, "(A)")

    write (u, "(A)")  "* Fixed:"
    write (u, "(A)")

    allocate (alpha_qcd_fixed_t :: qcd%alpha)
    call qcd%compute_alphas_md5sum ()

    call qcd%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha_s (mz)    =", &
         qcd%alpha%get (MZ_REF)
    write (u, "(1x,A,F10.7)")  "alpha_s (1 TeV) =", &
         qcd%alpha%get (1000._default)
    write (u, *)
    deallocate (qcd%alpha)

    write (u, "(A)")  "* Running from MZ (LO):"
    write (u, "(A)")

    allocate (alpha_qcd_from_scale_t :: qcd%alpha)
    call qcd%compute_alphas_md5sum ()

    call qcd%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha_s (mz)    =", &
         qcd%alpha%get (MZ_REF)
    write (u, "(1x,A,F10.7)")  "alpha_s (1 TeV) =", &
         qcd%alpha%get (1000._default)
    write (u, *)

    write (u, "(A)")  "* Running from MZ (NLO):"
    write (u, "(A)")

    select type (alpha => qcd%alpha)
    type is (alpha_qcd_from_scale_t)
       alpha%order = 1
    end select
    call qcd%compute_alphas_md5sum ()

    call qcd%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha_s (mz)    =", &
         qcd%alpha%get (MZ_REF)
    write (u, "(1x,A,F10.7)")  "alpha_s (1 TeV) =", &
         qcd%alpha%get (1000._default)
    write (u, *)

    write (u, "(A)")  "* Running from MZ (NNLO):"
    write (u, "(A)")

    select type (alpha => qcd%alpha)
    type is (alpha_qcd_from_scale_t)
       alpha%order = 2
    end select
    call qcd%compute_alphas_md5sum ()

    call qcd%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha_s (mz)    =", &
         qcd%alpha%get (MZ_REF)
    write (u, "(1x,A,F10.7)")  "alpha_s (1 TeV) =", &
         qcd%alpha%get (1000._default)
    write (u, *)

    deallocate (qcd%alpha)
    write (u, "(A)")  "* Running from Lambda_QCD (LO):"
    write (u, "(A)")

    allocate (alpha_qcd_from_lambda_t :: qcd%alpha)
    call qcd%compute_alphas_md5sum ()

    call qcd%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha_s (mz)    =", &
         qcd%alpha%get (MZ_REF)
    write (u, "(1x,A,F10.7)")  "alpha_s (1 TeV) =", &
         qcd%alpha%get (1000._default)
    write (u, *)

    write (u, "(A)")  "* Running from Lambda_QCD (NLO):"
    write (u, "(A)")

    select type (alpha => qcd%alpha)
    type is (alpha_qcd_from_lambda_t)
       alpha%order = 1
    end select
    call qcd%compute_alphas_md5sum ()

    call qcd%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha_s (mz)    =", &
         qcd%alpha%get (MZ_REF)
    write (u, "(1x,A,F10.7)")  "alpha_s (1 TeV) =", &
         qcd%alpha%get (1000._default)
    write (u, *)

    write (u, "(A)")  "* Running from Lambda_QCD (NNLO):"
    write (u, "(A)")

    select type (alpha => qcd%alpha)
    type is (alpha_qcd_from_lambda_t)
       alpha%order = 2
    end select
    call qcd%compute_alphas_md5sum ()

    call qcd%write (u)
    write (u, *)
    write (u, "(1x,A,F10.7)")  "alpha_s (mz)    =", &
         qcd%alpha%get (MZ_REF)
    write (u, "(1x,A,F10.7)")  "alpha_s (1 TeV) =", &
         qcd%alpha%get (1000._default)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sm_qcd_1"

  end subroutine sm_qcd_1


end module sm_qcd_uti
