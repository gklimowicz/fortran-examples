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

module sm_qcd

  use kinds, only: default
  use physics_defs

  implicit none
  private

  public :: alpha_qcd_t
  public :: alpha_qcd_fixed_t
  public :: alpha_qcd_from_scale_t
  public :: alpha_qcd_from_lambda_t
  public :: qcd_t

  type, abstract :: alpha_qcd_t
   contains
     procedure (alpha_qcd_write), deferred :: write
     procedure (alpha_qcd_get), deferred :: get
  end type alpha_qcd_t

  type, extends (alpha_qcd_t) :: alpha_qcd_fixed_t
     real(default) :: val = ALPHA_QCD_MZ_REF
   contains
     procedure :: write => alpha_qcd_fixed_write
     procedure :: get => alpha_qcd_fixed_get
  end type alpha_qcd_fixed_t

  type, extends (alpha_qcd_t) :: alpha_qcd_from_scale_t
     real(default) :: mu_ref = MZ_REF
     real(default) :: ref = ALPHA_QCD_MZ_REF
     integer :: order = 0
     integer :: nf = 5
   contains
     procedure :: write => alpha_qcd_from_scale_write
     procedure :: get => alpha_qcd_from_scale_get
  end type alpha_qcd_from_scale_t

  type, extends (alpha_qcd_t) :: alpha_qcd_from_lambda_t
     real(default) :: lambda = LAMBDA_QCD_REF
     integer :: order = 0
     integer :: nf = 5
   contains
     procedure :: write => alpha_qcd_from_lambda_write
     procedure :: get => alpha_qcd_from_lambda_get
  end type alpha_qcd_from_lambda_t

  type :: qcd_t
     class(alpha_qcd_t), allocatable :: alpha
     character(32) :: md5sum = ""
     integer :: n_f = -1
   contains
     procedure :: write => qcd_write
     procedure :: compute_alphas_md5sum => qcd_compute_alphas_md5sum
     procedure :: get_md5sum => qcd_get_md5sum
  end type qcd_t


  abstract interface
     subroutine alpha_qcd_write (object, unit)
       import
       class(alpha_qcd_t), intent(in) :: object
       integer, intent(in), optional :: unit
     end subroutine alpha_qcd_write
  end interface

  abstract interface
     function alpha_qcd_get (alpha_qcd, scale) result (alpha)
       import
       class(alpha_qcd_t), intent(in) :: alpha_qcd
       real(default), intent(in) :: scale
       real(default) :: alpha
     end function alpha_qcd_get
  end interface


  interface
    module subroutine alpha_qcd_fixed_write (object, unit)
      class(alpha_qcd_fixed_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine alpha_qcd_fixed_write
    module function alpha_qcd_fixed_get (alpha_qcd, scale) result (alpha)
      class(alpha_qcd_fixed_t), intent(in) :: alpha_qcd
      real(default), intent(in) :: scale
      real(default) :: alpha
    end function alpha_qcd_fixed_get
    module subroutine alpha_qcd_from_scale_write (object, unit)
      class(alpha_qcd_from_scale_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine alpha_qcd_from_scale_write
    module function alpha_qcd_from_scale_get (alpha_qcd, scale) result (alpha)
      class(alpha_qcd_from_scale_t), intent(in) :: alpha_qcd
      real(default), intent(in) :: scale
      real(default) :: alpha
    end function alpha_qcd_from_scale_get
    module subroutine alpha_qcd_from_lambda_write (object, unit)
      class(alpha_qcd_from_lambda_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine alpha_qcd_from_lambda_write
    module function alpha_qcd_from_lambda_get (alpha_qcd, scale) result (alpha)
      class(alpha_qcd_from_lambda_t), intent(in) :: alpha_qcd
      real(default), intent(in) :: scale
      real(default) :: alpha
    end function alpha_qcd_from_lambda_get
    module subroutine qcd_write (qcd, unit, show_md5sum)
      class(qcd_t), intent(in) :: qcd
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_md5sum
    end subroutine qcd_write
    module subroutine qcd_compute_alphas_md5sum (qcd)
      class(qcd_t), intent(inout) :: qcd
      integer :: unit
    end subroutine qcd_compute_alphas_md5sum
    module function qcd_get_md5sum (qcd) result (md5sum)
      character(32) :: md5sum
      class(qcd_t), intent(inout) :: qcd
    end function qcd_get_md5sum
  end interface

end module sm_qcd
