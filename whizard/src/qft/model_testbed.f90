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

module model_testbed

  use iso_varying_string, string_t => varying_string
  use model_data
  use var_base

  implicit none
  private

  public :: prepare_model
  public :: cleanup_model

  procedure (prepare_model_proc), pointer :: prepare_model => null ()
  procedure (cleanup_model_proc), pointer :: cleanup_model => null ()

  abstract interface
     subroutine prepare_model_proc (model, name, vars)
       import
       class(model_data_t), intent(inout), pointer :: model
       type(string_t), intent(in) :: name
       class(vars_t), pointer, intent(out), optional :: vars
     end subroutine prepare_model_proc
  end interface

  abstract interface
     subroutine cleanup_model_proc (model)
       import
       class(model_data_t), intent(inout), target :: model
     end subroutine cleanup_model_proc
  end interface


end module model_testbed
