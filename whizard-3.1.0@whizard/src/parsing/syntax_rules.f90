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

module syntax_rules

  use iso_varying_string, string_t => varying_string
  use ifiles, only: ifile_t
  use lexers

  implicit none
  private

  public :: S_UNKNOWN
  public :: S_LOGICAL, S_INTEGER, S_REAL, S_COMPLEX, S_QUOTED
  public :: S_IDENTIFIER, S_KEYWORD
  public :: S_SEQUENCE, S_LIST, S_GROUP, S_ARGS
  public :: S_ALTERNATIVE
  public :: S_IGNORE
  public :: syntax_rule_get_type
  public :: syntax_rule_get_key
  public :: syntax_rule_get_separator
  public :: syntax_rule_get_delimiter
  public :: syntax_rule_get_n_sub
  public :: syntax_rule_get_sub_ptr
  public :: syntax_rule_last_optional
  public :: syntax_rule_last_repetitive
  public :: syntax_rule_is_atomic
  public :: syntax_rule_write
  public :: syntax_t
  public :: syntax_init
  public :: syntax_final
  public :: syntax_get_rule_ptr
  public :: syntax_get_top_rule_ptr
  public :: syntax_get_keyword_list_ptr
  public :: syntax_write

  integer, parameter :: &
       S_UNKNOWN = 0, &
       S_LOGICAL = 1, S_INTEGER = 2, S_REAL = 3, S_COMPLEX = 4, &
       S_QUOTED = 5, S_IDENTIFIER = 6, S_KEYWORD = 7, &
       S_SEQUENCE = 8, S_LIST = 9, S_GROUP = 10, S_ARGS = 11, &
       S_ALTERNATIVE = 12, &
       S_IGNORE = 99

  type :: rule_p
     private
     type(syntax_rule_t), pointer :: p => null ()
  end type rule_p

  public :: syntax_rule_t
  type :: syntax_rule_t
     private
     integer :: type = S_UNKNOWN
     logical :: used = .false.
     type(string_t) :: keyword
     type(string_t) :: separator
     type(string_t), dimension(2) :: delimiter
     type(rule_p), dimension(:), allocatable :: child
     character(1) :: modifier = ""
     logical :: opt = .false., rep = .false.
   contains
     procedure :: get_key => syntax_rule_get_key
     procedure :: write => syntax_rule_write
  end type syntax_rule_t

  type :: syntax_t
     private
     type(syntax_rule_t), dimension(:), allocatable :: rule
     type(keyword_list_t) :: keyword_list
  end type syntax_t


  interface syntax_init
     module procedure syntax_init_from_ifile
  end interface


  interface
    module function syntax_rule_get_type (rule) result (type)
      integer :: type
      type(syntax_rule_t), intent(in) :: rule
    end function syntax_rule_get_type
    module function syntax_rule_get_key (rule) result (key)
      class(syntax_rule_t), intent(in) :: rule
      type(string_t) :: key
    end function syntax_rule_get_key
    module function syntax_rule_get_separator (rule) result (separator)
      type(string_t) :: separator
      type(syntax_rule_t), intent(in) :: rule
    end function syntax_rule_get_separator
    module function syntax_rule_get_delimiter (rule) result (delimiter)
      type(string_t), dimension(2) :: delimiter
      type(syntax_rule_t), intent(in) :: rule
    end function syntax_rule_get_delimiter
    module function syntax_rule_get_sub_ptr (rule, i) result (sub)
      type(syntax_rule_t), pointer :: sub
      type(syntax_rule_t), intent(in), target :: rule
      integer, intent(in) :: i
    end function syntax_rule_get_sub_ptr
    module function syntax_rule_get_n_sub (rule) result (n)
      integer :: n
      type(syntax_rule_t), intent(in) :: rule
    end function syntax_rule_get_n_sub
    module function syntax_rule_last_optional (rule) result (opt)
      logical :: opt
      type(syntax_rule_t), intent(in) :: rule
    end function syntax_rule_last_optional
    module function syntax_rule_last_repetitive (rule) result (rep)
      logical :: rep
      type(syntax_rule_t), intent(in) :: rule
    end function syntax_rule_last_repetitive
    module function syntax_rule_is_atomic (rule) result (atomic)
      logical :: atomic
      type(syntax_rule_t), intent(in) :: rule
    end function syntax_rule_is_atomic
    module subroutine syntax_rule_write (rule, unit, short, key_only, advance)
      class(syntax_rule_t), intent(in) :: rule
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: short, key_only, advance
    end subroutine syntax_rule_write
    module subroutine syntax_init_from_ifile (syntax, ifile)
      type(syntax_t), intent(out), target :: syntax
      type(ifile_t), intent(in) :: ifile
    end subroutine syntax_init_from_ifile
    module function syntax_get_rule_ptr (syntax, key) result (rule)
      type(syntax_rule_t), pointer :: rule
      type(syntax_t), intent(in), target :: syntax
      type(string_t), intent(in) :: key
    end function syntax_get_rule_ptr
    module subroutine syntax_final (syntax)
      type(syntax_t), intent(inout) :: syntax
    end subroutine syntax_final
    module function syntax_get_top_rule_ptr (syntax) result (rule)
      type(syntax_rule_t), pointer :: rule
      type(syntax_t), intent(in), target :: syntax
    end function syntax_get_top_rule_ptr
    module function syntax_get_keyword_list_ptr (syntax) result (keyword_list)
      type(keyword_list_t), pointer :: keyword_list
      type(syntax_t), intent(in), target :: syntax
    end function syntax_get_keyword_list_ptr
    module subroutine syntax_write (syntax, unit)
      type(syntax_t), intent(in) :: syntax
      integer, intent(in), optional :: unit
    end subroutine syntax_write
  end interface

end module syntax_rules
