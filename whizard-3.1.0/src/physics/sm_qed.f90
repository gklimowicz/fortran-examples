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

module sm_qed

  use kinds, only: default
  use physics_defs

  implicit none
  private

  public :: alpha_qed_t
  public :: alpha_qed_fixed_t
  public :: alpha_qed_from_scale_t
  public :: qed_t

  type, abstract :: alpha_qed_t
   contains
     procedure (alpha_qed_write), deferred :: write
     procedure (alpha_qed_get), deferred :: get
  end type alpha_qed_t

  type, extends (alpha_qed_t) :: alpha_qed_fixed_t
     real(default) :: val = ALPHA_QED_ME_REF
   contains
     procedure :: write => alpha_qed_fixed_write
     procedure :: get => alpha_qed_fixed_get
  end type alpha_qed_fixed_t

  type, extends (alpha_qed_t) :: alpha_qed_from_scale_t
     real(default) :: mu_ref = ME_REF
     real(default) :: ref = ALPHA_QED_ME_REF
     integer :: order = 0
     integer :: nf = 5
     integer :: nlep = 1
     logical :: analytic = .true.
   contains
     procedure :: write => alpha_qed_from_scale_write
     procedure :: get => alpha_qed_from_scale_get
  end type alpha_qed_from_scale_t

  type :: qed_t
     class(alpha_qed_t), allocatable :: alpha
     character(32) :: md5sum = ""
     integer :: n_f = -1
     integer :: n_lep = -1
   contains
     procedure :: write => qed_write
     procedure :: compute_alpha_md5sum => qed_compute_alpha_md5sum
     procedure :: get_md5sum => qed_get_md5sum
  end type qed_t


  abstract interface
     subroutine alpha_qed_write (object, unit)
       import
       class(alpha_qed_t), intent(in) :: object
       integer, intent(in), optional :: unit
     end subroutine alpha_qed_write
  end interface

  abstract interface
     function alpha_qed_get (alpha_qed, scale) result (alpha)
       import
       class(alpha_qed_t), intent(in) :: alpha_qed
       real(default), intent(in) :: scale
       real(default) :: alpha
     end function alpha_qed_get
  end interface


  interface
    module subroutine alpha_qed_fixed_write (object, unit)
      class(alpha_qed_fixed_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine alpha_qed_fixed_write
    module function alpha_qed_fixed_get (alpha_qed, scale) result (alpha)
      class(alpha_qed_fixed_t), intent(in) :: alpha_qed
      real(default), intent(in) :: scale
      real(default) :: alpha
    end function alpha_qed_fixed_get
    module subroutine alpha_qed_from_scale_write (object, unit)
      class(alpha_qed_from_scale_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine alpha_qed_from_scale_write
    module function alpha_qed_from_scale_get (alpha_qed, scale) result (alpha)
      class(alpha_qed_from_scale_t), intent(in) :: alpha_qed
      real(default), intent(in) :: scale
      real(default) :: alpha
    end function alpha_qed_from_scale_get
    module subroutine qed_write (qed, unit, show_md5sum)
      class(qed_t), intent(in) :: qed
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: show_md5sum
    end subroutine qed_write
    module subroutine qed_compute_alpha_md5sum (qed)
      class(qed_t), intent(inout) :: qed
      integer :: unit
    end subroutine qed_compute_alpha_md5sum
    module function qed_get_md5sum (qed) result (md5sum)
      character(32) :: md5sum
      class(qed_t), intent(inout) :: qed
    end function qed_get_md5sum
  end interface

end module sm_qed
