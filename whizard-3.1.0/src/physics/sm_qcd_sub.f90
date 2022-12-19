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

submodule (sm_qcd) sm_qcd_s

  use io_units
  use format_defs, only: FMT_12
  use numeric_utils
  use diagnostics
  use md5
  use sm_physics

  implicit none

contains

  module subroutine alpha_qcd_fixed_write (object, unit)
    class(alpha_qcd_fixed_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3x,A)")  "QCD parameters (fixed coupling):"
    write (u, "(5x,A," // FMT_12 // ")")  "alpha = ", object%val
  end subroutine alpha_qcd_fixed_write

  module function alpha_qcd_fixed_get (alpha_qcd, scale) result (alpha)
    class(alpha_qcd_fixed_t), intent(in) :: alpha_qcd
    real(default), intent(in) :: scale
    real(default) :: alpha
    alpha = alpha_qcd%val
  end function alpha_qcd_fixed_get

  module subroutine alpha_qcd_from_scale_write (object, unit)
    class(alpha_qcd_from_scale_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3x,A)")  "QCD parameters (running coupling):"
    write (u, "(5x,A," // FMT_12 // ")")  "Scale mu  = ", object%mu_ref
    write (u, "(5x,A," // FMT_12 // ")")  "alpha(mu) = ", object%ref
    write (u, "(5x,A,I0)")      "LL order  = ", object%order
    write (u, "(5x,A,I0)")      "N(flv)    = ", object%nf
  end subroutine alpha_qcd_from_scale_write

  module function alpha_qcd_from_scale_get (alpha_qcd, scale) result (alpha)
    class(alpha_qcd_from_scale_t), intent(in) :: alpha_qcd
    real(default), intent(in) :: scale
    real(default) :: alpha
    alpha = running_as (scale, alpha_qcd%ref, alpha_qcd%mu_ref, &
         alpha_qcd%order, real (alpha_qcd%nf, kind=default))
  end function alpha_qcd_from_scale_get

  module subroutine alpha_qcd_from_lambda_write (object, unit)
    class(alpha_qcd_from_lambda_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3x,A)")  "QCD parameters (Lambda_QCD as input):"
    write (u, "(5x,A," // FMT_12 // ")")  "Lambda_QCD = ", object%lambda
    write (u, "(5x,A,I0)")      "LL order   = ", object%order
    write (u, "(5x,A,I0)")      "N(flv)     = ", object%nf
  end subroutine alpha_qcd_from_lambda_write

  module function alpha_qcd_from_lambda_get (alpha_qcd, scale) result (alpha)
    class(alpha_qcd_from_lambda_t), intent(in) :: alpha_qcd
    real(default), intent(in) :: scale
    real(default) :: alpha
    alpha = running_as_lam (real (alpha_qcd%nf, kind=default), scale, &
         alpha_qcd%lambda, alpha_qcd%order)
  end function alpha_qcd_from_lambda_get

  module subroutine qcd_write (qcd, unit, show_md5sum)
    class(qcd_t), intent(in) :: qcd
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_md5sum
    logical :: show_md5
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    show_md5 = .true.;  if (present (show_md5sum))  show_md5 = show_md5sum
    if (allocated (qcd%alpha)) then
       call qcd%alpha%write (u)
    else
       write (u, "(3x,A)")  "QCD parameters (coupling undefined)"
    end if
    if (show_md5 .and. qcd%md5sum /= "") &
         write (u, "(5x,A,A,A)") "md5sum = '", qcd%md5sum, "'"
  end subroutine qcd_write

  module subroutine qcd_compute_alphas_md5sum (qcd)
    class(qcd_t), intent(inout) :: qcd
    integer :: unit
    if (allocated (qcd%alpha)) then
       unit = free_unit ()
       open (unit, status="scratch", action="readwrite")
       call qcd%alpha%write (unit)
       rewind (unit)
       qcd%md5sum = md5sum (unit)
       close (unit)
    end if
  end subroutine qcd_compute_alphas_md5sum

  module function qcd_get_md5sum (qcd) result (md5sum)
    character(32) :: md5sum
    class(qcd_t), intent(inout) :: qcd
    md5sum = qcd%md5sum
  end function qcd_get_md5sum


end submodule sm_qcd_s

