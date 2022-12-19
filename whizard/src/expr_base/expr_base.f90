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

module expr_base

  use kinds, only: default
  use var_base

  implicit none
  private

  public :: expr_t
  public :: expr_factory_t

  type, abstract :: expr_t
   contains
    procedure(expr_final), deferred :: final
    procedure(expr_write), deferred :: write
    procedure(expr_setup_expr), deferred :: setup_expr
    procedure(expr_setup_expr), deferred :: setup_lexpr
    procedure(expr_evaluate), deferred :: evaluate
    procedure(expr_is_known), deferred :: is_known
    procedure(expr_get_log), deferred :: get_log
    procedure(expr_get_real), deferred :: get_real
  end type expr_t

  type, abstract :: expr_factory_t
   contains
    procedure(expr_factory_write), deferred :: write
    procedure(expr_factory_build), deferred :: build
  end type expr_factory_t


  abstract interface
     subroutine expr_final (expr)
       import
       class(expr_t), intent(inout) :: expr
     end subroutine expr_final
  end interface

  abstract interface
     subroutine expr_write (expr, unit, write_vars)
       import
       class(expr_t), intent(in) :: expr
       integer, intent(in), optional :: unit
       logical, intent(in), optional :: write_vars
     end subroutine expr_write
  end interface

  abstract interface
     subroutine expr_setup_expr (expr, vars)
       import
       class(expr_t), intent(inout), target :: expr
       class(vars_t), intent(in), target :: vars
     end subroutine expr_setup_expr
  end interface

  abstract interface
     subroutine expr_evaluate (expr)
       import
       class(expr_t), intent(inout) :: expr
     end subroutine expr_evaluate
  end interface

  abstract interface
     function expr_is_known (expr) result (flag)
       import
       class(expr_t), intent(in) :: expr
       logical :: flag
     end function expr_is_known
  end interface

  abstract interface
     function expr_get_log (expr) result (lval)
       import
       class(expr_t), intent(in) :: expr
       logical :: lval
     end function expr_get_log
  end interface

  abstract interface
     function expr_get_real (expr) result (rval)
       import
       class(expr_t), intent(in) :: expr
       real(default) :: rval
     end function expr_get_real
  end interface

  abstract interface
     subroutine expr_factory_write (expr_factory, unit)
       import
       class(expr_factory_t), intent(in) :: expr_factory
       integer, intent(in), optional :: unit
     end subroutine expr_factory_write
  end interface

  abstract interface
     subroutine expr_factory_build (expr_factory, expr)
       import
       class(expr_factory_t), intent(in) :: expr_factory
       class(expr_t), intent(out), allocatable :: expr
     end subroutine expr_factory_build
  end interface


end module expr_base
