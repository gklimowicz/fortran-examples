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

submodule (sm_qed) sm_qed_s

  use io_units
  use format_defs, only: FMT_12
  use md5
  use sm_physics

  implicit none

contains

  module subroutine alpha_qed_fixed_write (object, unit)
    class(alpha_qed_fixed_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3x,A)")  "QED parameters (fixed coupling):"
    write (u, "(5x,A," // FMT_12 // ")")  "alpha = ", object%val
  end subroutine alpha_qed_fixed_write

  module function alpha_qed_fixed_get (alpha_qed, scale) result (alpha)
    class(alpha_qed_fixed_t), intent(in) :: alpha_qed
    real(default), intent(in) :: scale
    real(default) :: alpha
    alpha = alpha_qed%val
  end function alpha_qed_fixed_get

  module subroutine alpha_qed_from_scale_write (object, unit)
    class(alpha_qed_from_scale_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3x,A)")  "QED parameters (running coupling):"
    write (u, "(5x,A," // FMT_12 // ")")  "Scale mu  = ", object%mu_ref
    write (u, "(5x,A," // FMT_12 // ")")  "alpha(mu) = ", object%ref
    write (u, "(5x,A,I0)")      "LL order  = ", object%order
    write (u, "(5x,A,I0)")      "N(flv)    = ", object%nf
    write (u, "(5x,A,I0)")      "N(lep)    = ", object%nlep
    write (u, "(5x,A,L1)")      "analytic  = ", object%analytic
  end subroutine alpha_qed_from_scale_write

  module function alpha_qed_from_scale_get (alpha_qed, scale) result (alpha)
    class(alpha_qed_from_scale_t), intent(in) :: alpha_qed
    real(default), intent(in) :: scale
    real(default) :: alpha
    if (alpha_qed%analytic) then
       alpha = running_alpha (scale, alpha_qed%ref, alpha_qed%mu_ref, &
            alpha_qed%order, alpha_qed%nf, alpha_qed%nlep)
    else
       alpha = running_alpha_num (scale, alpha_qed%ref, alpha_qed%mu_ref, &
            alpha_qed%order, alpha_qed%nf, alpha_qed%nlep)
    end if
  end function alpha_qed_from_scale_get

  module subroutine qed_write (qed, unit, show_md5sum)
    class(qed_t), intent(in) :: qed
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_md5sum
    logical :: show_md5
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    show_md5 = .true.;  if (present (show_md5sum))  show_md5 = show_md5sum
    if (allocated (qed%alpha)) then
       call qed%alpha%write (u)
    else
       write (u, "(3x,A)")  "QED parameters (coupling undefined)"
    end if
    if (show_md5 .and. qed%md5sum /= "") &
         write (u, "(5x,A,A,A)") "md5sum = '", qed%md5sum, "'"
  end subroutine qed_write

  module subroutine qed_compute_alpha_md5sum (qed)
    class(qed_t), intent(inout) :: qed
    integer :: unit
    if (allocated (qed%alpha)) then
       unit = free_unit ()
       open (unit, status="scratch", action="readwrite")
       call qed%alpha%write (unit)
       rewind (unit)
       qed%md5sum = md5sum (unit)
       close (unit)
    end if
  end subroutine qed_compute_alpha_md5sum

  module function qed_get_md5sum (qed) result (md5sum)
    character(32) :: md5sum
    class(qed_t), intent(inout) :: qed
    md5sum = qed%md5sum
  end function qed_get_md5sum


end submodule sm_qed_s

