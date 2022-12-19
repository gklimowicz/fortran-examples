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

module lexers

  use iso_varying_string, string_t => varying_string
  use ifiles, only: ifile_t
  use ifiles, only: line_p

  implicit none
  private

  public :: stream_t
  public :: stream_init
  public :: stream_final
  public :: stream_get_record
  public :: stream_get_source_info_string
  public :: stream_get_record_info_string
  public :: keyword_list_t
  public :: keyword_list_add
  public :: keyword_list_contains
  public :: keyword_list_write
  public :: keyword_list_final
  public :: T_KEYWORD, T_IDENTIFIER, T_QUOTED, T_NUMERIC
  public :: lexeme_t
  public :: lexeme_write
  public :: lexeme_get_string
  public :: lexeme_get_contents
  public :: lexeme_get_delimiters
  public :: lexeme_get_type
  public :: lexeme_is_break
  public :: lexeme_is_eof
  public :: lexer_t
  public :: lexer_init
  public :: lexer_clear
  public :: lexer_final
  public :: lexer_assign_stream
  public :: lex
  public :: lexer_put_back
  public :: lexer_write_setup
  public :: lexer_show_location

  integer, parameter :: T_KEYWORD = 1
  integer, parameter :: T_IDENTIFIER = 2, T_QUOTED = 3, T_NUMERIC = 4
  integer, parameter :: EMPTY = 0, WHITESPACE = 10
  integer, parameter :: NO_MATCH = 11, IO_ERROR = 12, OVERFLOW = 13
  integer, parameter :: UNMATCHED_QUOTE = 14
  integer, parameter :: CASE_KEEP = 0, CASE_UP = 1, CASE_DOWN = 2


  type :: stream_t
     type(string_t), pointer :: filename => null ()
     integer, pointer :: unit => null ()
     type(string_t), pointer :: string => null ()
     type(ifile_t), pointer :: ifile => null ()
     type(line_p), pointer :: line => null ()
     integer :: record = 0
     logical :: eof = .false.
   contains
     generic :: init => &
          stream_init_filename, &
          stream_init_unit, &
          stream_init_string, &
          stream_init_ifile, &
          stream_init_line
     procedure, private :: stream_init_filename
     procedure, private :: stream_init_unit
     procedure, private :: stream_init_string
     procedure, private :: stream_init_ifile
     procedure, private :: stream_init_line
     procedure :: final => stream_final
  end type stream_t

  type :: keyword_entry_t
     private
     type(string_t) :: string
     type(keyword_entry_t), pointer :: next => null ()
  end type keyword_entry_t

  type :: keyword_list_t
     private
     type(keyword_entry_t), pointer :: first => null ()
     type(keyword_entry_t), pointer :: last => null ()
  end type keyword_list_t

  type :: template_t
     !!! !!! Compiler bug in ifort 20/21/22: no structure constants for 
     !!! !!! types with private components in submodules possible
     !!! private
     integer :: type
     character(256) :: charset1, charset2
     integer :: len1, len2
  end type template_t

  type :: lexer_setup_t
     private
     type(template_t), dimension(:), allocatable :: tt
     integer, dimension(:), allocatable :: type
     integer :: keyword_case = CASE_KEEP
     type(keyword_list_t), pointer :: keyword_list => null ()
  end type lexer_setup_t

  type :: lexeme_t
     private
     integer :: type = EMPTY
     type(string_t) :: s
     integer :: b = 0, e = 0
  end type lexeme_t

  type :: lexer_t
     private
     type(lexer_setup_t) :: setup
     type(stream_t), pointer :: stream => null ()
     type(lexeme_t) :: lexeme
     type(string_t) :: previous_line2
     type(string_t) :: previous_line1
     type(string_t) :: current_line
     integer :: lines_read = 0
     integer :: current_column = 0
     integer :: previous_column = 0
     type(string_t) :: buffer
     type(lexer_t), pointer :: parent => null ()
   contains
     procedure :: init => lexer_init
     procedure :: clear => lexer_clear
     procedure :: final => lexer_final
     procedure :: assign_stream => lexer_assign_stream
  end type lexer_t

  interface stream_init
     module procedure stream_init_filename
     module procedure stream_init_unit
     module procedure stream_init_string
     module procedure stream_init_ifile
     module procedure stream_init_line
  end interface
  interface keyword_list_write
     module procedure keyword_list_write_unit
  end interface

  interface
    module subroutine stream_init_filename (stream, filename)
      class(stream_t), intent(out) :: stream
      character(*), intent(in) :: filename
    end subroutine stream_init_filename
    module subroutine stream_init_unit (stream, unit)
      class(stream_t), intent(out) :: stream
      integer, intent(in) :: unit
    end subroutine stream_init_unit
    module subroutine stream_init_string (stream, string)
      class(stream_t), intent(out) :: stream
      type(string_t), intent(in) :: string
    end subroutine stream_init_string
    module subroutine stream_init_ifile (stream, ifile)
      class(stream_t), intent(out) :: stream
      type(ifile_t), intent(in) :: ifile
    end subroutine stream_init_ifile
    module subroutine stream_init_line (stream, line)
      class(stream_t), intent(out) :: stream
      type(line_p), intent(in) :: line
    end subroutine stream_init_line
    module subroutine stream_final (stream)
      class(stream_t), intent(inout) :: stream
    end subroutine stream_final
    module subroutine stream_get_record (stream, string, iostat)
      type(stream_t), intent(inout) :: stream
      type(string_t), intent(out) :: string
      integer, intent(out) :: iostat
    end subroutine stream_get_record
    module function stream_get_source_info_string (stream) result (string)
      type(string_t) :: string
      type(stream_t), intent(in) :: stream
    end function stream_get_source_info_string
    module function stream_get_record_info_string (stream) result (string)
      type(string_t) :: string
      type(stream_t), intent(in) :: stream
    end function stream_get_record_info_string
    module subroutine keyword_list_add (keylist, string)
      type(keyword_list_t), intent(inout) :: keylist
      type(string_t), intent(in) :: string
    end subroutine keyword_list_add
    module function keyword_list_contains (keylist, string) result (found)
      type(keyword_list_t), intent(in) :: keylist
      type(string_t), intent(in) :: string
      logical :: found
    end function keyword_list_contains
    module subroutine keyword_list_write_unit (keylist, unit)
      type(keyword_list_t), intent(in) :: keylist
      integer, intent(in) :: unit
    end subroutine keyword_list_write_unit
    module subroutine keyword_list_final (keylist)
      type(keyword_list_t), intent(inout) :: keylist
    end subroutine keyword_list_final
    module subroutine lexeme_write (t, unit)
      type(lexeme_t), intent(in) :: t
      integer, intent(in) :: unit
    end subroutine lexeme_write
    module function lexeme_get_string (t) result (s)
      type(string_t) :: s
      type(lexeme_t), intent(in) :: t
    end function lexeme_get_string
    module function lexeme_get_contents (t) result (s)
      type(string_t) :: s
      type(lexeme_t), intent(in) :: t
    end function lexeme_get_contents
    module function lexeme_get_delimiters (t) result (del)
      type(string_t), dimension(2) :: del
      type(lexeme_t), intent(in) :: t
    end function lexeme_get_delimiters
    module function lexeme_get_type (t) result (type)
      integer :: type
      type(lexeme_t), intent(in) :: t
    end function lexeme_get_type
    module function lexeme_is_break (t) result (break)
      logical :: break
      type(lexeme_t), intent(in) :: t
    end function lexeme_is_break
    module function lexeme_is_eof (t) result (ok)
      logical :: ok
      type(lexeme_t), intent(in) :: t
    end function lexeme_is_eof
    module subroutine lexer_init (lexer, &
         comment_chars, quote_chars, quote_match, &
         single_chars, special_class, &
         keyword_list, upper_case_keywords, &
         parent)
      class(lexer_t), intent(inout) :: lexer
      character(*), intent(in) :: comment_chars
      character(*), intent(in) :: quote_chars, quote_match
      character(*), intent(in) :: single_chars
      character(*), dimension(:), intent(in) :: special_class
      type(keyword_list_t), pointer :: keyword_list
      logical, intent(in), optional :: upper_case_keywords
      type(lexer_t), target, intent(in), optional :: parent
    end subroutine lexer_init
    module subroutine lexer_clear (lexer)
      class(lexer_t), intent(inout) :: lexer
    end subroutine lexer_clear
    module subroutine lexer_final (lexer)
      class(lexer_t), intent(inout) :: lexer
    end subroutine lexer_final
    module subroutine lexer_assign_stream (lexer, stream)
      class(lexer_t), intent(inout) :: lexer
      type(stream_t), intent(in), target :: stream
    end subroutine lexer_assign_stream
    module subroutine lex (lexeme, lexer)
      type(lexeme_t), intent(out) :: lexeme
      type(lexer_t), intent(inout) :: lexer
    end subroutine lex
    module subroutine lexer_put_back (lexer, lexeme)
      type(lexer_t), intent(inout) :: lexer
      type(lexeme_t), intent(in) :: lexeme
    end subroutine lexer_put_back
    module subroutine lexer_write_setup (lexer, unit)
      type(lexer_t), intent(in) :: lexer
      integer, intent(in), optional :: unit
    end subroutine lexer_write_setup
    module subroutine lexer_show_location (lexer)
      type(lexer_t), intent(in) :: lexer
    end subroutine lexer_show_location
  end interface

end module lexers
