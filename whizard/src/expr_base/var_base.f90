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

module var_base

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: vars_t

  type, abstract :: vars_t
   contains
     procedure (vars_link), deferred :: link
     procedure (vars_final), deferred :: final
     procedure (vars_get_lval), deferred :: contains
     procedure (vars_get_lval), deferred :: is_known
     procedure (vars_get_ival), deferred :: get_ival
     procedure (vars_get_rval), deferred :: get_rval
     procedure (vars_get_cval), deferred :: get_cval
     procedure (vars_get_lval), deferred :: get_lval
     procedure (vars_get_sval), deferred :: get_sval
     procedure (vars_unset), deferred :: unset
     procedure (vars_set_ival), deferred :: set_ival
     procedure (vars_set_rval), deferred :: set_rval
     procedure (vars_set_cval), deferred :: set_cval
     procedure (vars_set_lval), deferred :: set_lval
     procedure (vars_set_sval), deferred :: set_sval
  end type vars_t


  abstract interface
     subroutine vars_link (vars, target_vars)
       import
       class(vars_t), intent(inout) :: vars
       class(vars_t), intent(in), target :: target_vars
     end subroutine vars_link
  end interface

  abstract interface
     subroutine vars_final (vars, follow_link)
       import
       class(vars_t), intent(inout) :: vars
       logical, intent(in), optional :: follow_link
     end subroutine vars_final
  end interface

  abstract interface
     function vars_get_lval (vars, name, follow_link) result (lval)
       import
       logical :: lval
       class(vars_t), intent(in) :: vars
       type(string_t), intent(in) :: name
       logical, intent(in), optional :: follow_link
     end function vars_get_lval
  end interface

  abstract interface
     function vars_get_ival (vars, name, follow_link) result (ival)
       import
       integer :: ival
       class(vars_t), intent(in) :: vars
       type(string_t), intent(in) :: name
       logical, intent(in), optional :: follow_link
     end function vars_get_ival
  end interface

  abstract interface
     function vars_get_rval (vars, name, follow_link) result (rval)
       import
       real(default) :: rval
       class(vars_t), intent(in) :: vars
       type(string_t), intent(in) :: name
       logical, intent(in), optional :: follow_link
     end function vars_get_rval
  end interface

  abstract interface
     function vars_get_cval (vars, name, follow_link) result (cval)
       import
       complex(default) :: cval
       class(vars_t), intent(in) :: vars
       type(string_t), intent(in) :: name
       logical, intent(in), optional :: follow_link
     end function vars_get_cval
  end interface

  abstract interface
     function vars_get_sval (vars, name, follow_link) result (sval)
       import
       type(string_t) :: sval
       class(vars_t), intent(in) :: vars
       type(string_t), intent(in) :: name
       logical, intent(in), optional :: follow_link
     end function vars_get_sval
  end interface

  abstract interface
     subroutine vars_unset (vars, name, follow_link)
       import
       class(vars_t), intent(inout) :: vars
       type(string_t), intent(in) :: name
       logical, intent(in), optional :: follow_link
     end subroutine vars_unset
  end interface

  abstract interface
     subroutine vars_set_ival (vars, name, ival, follow_link)
       import
       class(vars_t), intent(inout) :: vars
       type(string_t), intent(in) :: name
       integer, intent(in) :: ival
       logical, intent(in), optional :: follow_link
     end subroutine vars_set_ival
  end interface

  abstract interface
     subroutine vars_set_rval (vars, name, rval, follow_link)
       import
       class(vars_t), intent(inout) :: vars
       type(string_t), intent(in) :: name
       real(default), intent(in) :: rval
       logical, intent(in), optional :: follow_link
     end subroutine vars_set_rval
  end interface

  abstract interface
     subroutine vars_set_cval (vars, name, cval, follow_link)
       import
       class(vars_t), intent(inout) :: vars
       type(string_t), intent(in) :: name
       complex(default), intent(in) :: cval
       logical, intent(in), optional :: follow_link
     end subroutine vars_set_cval
  end interface

  abstract interface
     subroutine vars_set_lval (vars, name, lval, follow_link)
       import
       class(vars_t), intent(inout) :: vars
       type(string_t), intent(in) :: name
       logical, intent(in) :: lval
       logical, intent(in), optional :: follow_link
     end subroutine vars_set_lval
  end interface

  abstract interface
     subroutine vars_set_sval (vars, name, sval, follow_link)
       import
       class(vars_t), intent(inout) :: vars
       type(string_t), intent(in) :: name
       type(string_t), intent(in) :: sval
       logical, intent(in), optional :: follow_link
     end subroutine vars_set_sval
  end interface


end module var_base
