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

submodule (lexers) lexers_s

  use io_units
  use string_utils
  use system_defs, only: EOF, EOR
  use system_defs, only: LF
  use system_defs, only: WHITESPACE_CHARS, LCLETTERS, UCLETTERS, DIGIT_CHARS
  use ifiles, only: line_get_string_advance
  use ifiles, only: line_is_associated, line_init, line_final
  use diagnostics

  implicit none

contains

  module subroutine stream_init_filename (stream, filename)
    class(stream_t), intent(out) :: stream
    character(*), intent(in) :: filename
    integer :: unit
    unit = free_unit ()
    open (unit=unit, file=filename, status="old", action="read")
    call stream_init_unit (stream, unit)
    allocate (stream%filename)
    stream%filename = filename
  end subroutine stream_init_filename

  module subroutine stream_init_unit (stream, unit)
    class(stream_t), intent(out) :: stream
    integer, intent(in) :: unit
    allocate (stream%unit)
    stream%unit = unit
    stream%eof = .false.
  end subroutine stream_init_unit

  module subroutine stream_init_string (stream, string)
    class(stream_t), intent(out) :: stream
    type(string_t), intent(in) :: string
    allocate (stream%string)
    stream%string = string
  end subroutine stream_init_string

  module subroutine stream_init_ifile (stream, ifile)
    class(stream_t), intent(out) :: stream
    type(ifile_t), intent(in) :: ifile
    type(line_p) :: line
    call line_init (line, ifile)
    call stream_init_line (stream, line)
    allocate (stream%ifile)
    stream%ifile = ifile
  end subroutine stream_init_ifile

  module subroutine stream_init_line (stream, line)
    class(stream_t), intent(out) :: stream
    type(line_p), intent(in) :: line
    allocate (stream%line)
    stream%line = line
  end subroutine stream_init_line

  module subroutine stream_final (stream)
    class(stream_t), intent(inout) :: stream
    if (associated (stream%filename)) then
       close (stream%unit)
       deallocate (stream%unit)
       deallocate (stream%filename)
    else if (associated (stream%unit)) then
       deallocate (stream%unit)
    else if (associated (stream%string)) then
       deallocate (stream%string)
    else if (associated (stream%ifile)) then
       call line_final (stream%line)
       deallocate (stream%line)
       deallocate (stream%ifile)
    else if (associated (stream%line)) then
       call line_final (stream%line)
       deallocate (stream%line)
    end if
  end subroutine stream_final

  module subroutine stream_get_record (stream, string, iostat)
    type(stream_t), intent(inout) :: stream
    type(string_t), intent(out) :: string
    integer, intent(out) :: iostat
    if (associated (stream%unit)) then
       if (stream%eof) then
          iostat = EOF
       else
          call get (stream%unit, string, iostat=iostat)
          if (iostat == EOR) then
             iostat = 0
             stream%record = stream%record + 1
          end if
          if (iostat == EOF) then
             iostat = 0
             stream%eof = .true.
             if (len (string) /= 0) stream%record = stream%record + 1
          end if
       end if
    else if (associated (stream%string)) then
       if (len (stream%string) /= 0) then
          string = stream%string
          stream%string = ""
          iostat = 0
          stream%record = stream%record + 1
       else
          string = ""
          iostat = EOF
       end if
    else if (associated (stream%line)) then
       if (line_is_associated (stream%line)) then
          string = line_get_string_advance (stream%line)
          iostat = 0
          stream%record = stream%record + 1
       else
          string = ""
          iostat = EOF
       end if
    else
       call msg_bug (" Attempt to read from uninitialized input stream")
    end if
  end subroutine stream_get_record

  module function stream_get_source_info_string (stream) result (string)
    type(string_t) :: string
    type(stream_t), intent(in) :: stream
    character(20) :: buffer
    if (associated (stream%filename)) then
       string = "File '" // stream%filename // "' (unit = "
       write (buffer, "(I0)")  stream%unit
       string = string // trim (buffer) // ")"
    else if (associated (stream%unit)) then
       write (buffer, "(I0)")  stream%unit
       string = "Unit " // trim (buffer)
    else if (associated (stream%string)) then
       string = "Input string"
    else if (associated (stream%ifile) .or. associated (stream%line)) then
       string = "Internal file"
    else
       string = ""
    end if
  end function stream_get_source_info_string

  module function stream_get_record_info_string (stream) result (string)
    type(string_t) :: string
    type(stream_t), intent(in) :: stream
    character(20) :: buffer
    string = stream_get_source_info_string (stream)
    if (string /= "")  string = string // ", "
    write (buffer, "(I0)")  stream%record
    string = string // "line " // trim (buffer)
  end function stream_get_record_info_string

  module subroutine keyword_list_add (keylist, string)
    type(keyword_list_t), intent(inout) :: keylist
    type(string_t), intent(in) :: string
    type(keyword_entry_t), pointer :: k_entry_new
    if (.not. keyword_list_contains (keylist, string)) then
       allocate (k_entry_new)
       k_entry_new%string = string
       if (associated (keylist%first)) then
          keylist%last%next => k_entry_new
       else
          keylist%first => k_entry_new
       end if
       keylist%last => k_entry_new
    end if
  end subroutine keyword_list_add

  module function keyword_list_contains (keylist, string) result (found)
    type(keyword_list_t), intent(in) :: keylist
    type(string_t), intent(in) :: string
    logical :: found
    found = .false.
    call check_rec (keylist%first)
  contains
    recursive subroutine check_rec (k_entry)
      type(keyword_entry_t), pointer :: k_entry
      if (associated (k_entry)) then
         if (k_entry%string /= string) then
            call check_rec (k_entry%next)
         else
            found = .true.
         end if
      end if
    end subroutine check_rec
  end function keyword_list_contains

  module subroutine keyword_list_write_unit (keylist, unit)
    type(keyword_list_t), intent(in) :: keylist
    integer, intent(in) :: unit
    write (unit, "(A)") "Keyword list:"
    if (associated (keylist%first)) then
       call keyword_write_rec (keylist%first)
       write (unit, *)
    else
       write (unit, "(1x,A)") "[empty]"
    end if
  contains
    recursive subroutine keyword_write_rec (k_entry)
      type(keyword_entry_t), intent(in), pointer :: k_entry
      if (associated (k_entry)) then
         write (unit, "(1x,A)", advance="no")  char (k_entry%string)
         call keyword_write_rec (k_entry%next)
      end if
    end subroutine keyword_write_rec
  end subroutine keyword_list_write_unit

  module subroutine keyword_list_final (keylist)
    type(keyword_list_t), intent(inout) :: keylist
    call keyword_destroy_rec (keylist%first)
    nullify (keylist%last)
  contains
    recursive subroutine keyword_destroy_rec (k_entry)
      type(keyword_entry_t), pointer :: k_entry
      if (associated (k_entry)) then
         call keyword_destroy_rec (k_entry%next)
         deallocate (k_entry)
      end if
    end subroutine keyword_destroy_rec
  end subroutine keyword_list_final

  subroutine lexeme_type_write (type, unit)
    integer, intent(in) :: type
    integer, intent(in) :: unit
    select case (type)
    case (EMPTY);       write(unit,"(A)",advance="no") " EMPTY      "
    case (WHITESPACE);  write(unit,"(A)",advance="no") " WHITESPACE "
    case (T_IDENTIFIER);write(unit,"(A)",advance="no") " IDENTIFIER "
    case (T_QUOTED);    write(unit,"(A)",advance="no") " QUOTED     "
    case (T_NUMERIC);   write(unit,"(A)",advance="no") " NUMERIC    "
    case (IO_ERROR);    write(unit,"(A)",advance="no") " IO_ERROR   "
    case (OVERFLOW);    write(unit,"(A)",advance="no") " OVERFLOW   "
    case (UNMATCHED_QUOTE);    write(unit,"(A)",advance="no") " UNMATCHEDQ "
    case (NO_MATCH);    write(unit,"(A)",advance="no") " NO_MATCH   "
    case (EOF);         write(unit,"(A)",advance="no") " EOF        "
    case default;       write(unit,"(A)",advance="no") " [illegal]  "
    end select
  end subroutine lexeme_type_write

  subroutine template_write (tt, unit)
    type(template_t), intent(in) :: tt
    integer, intent(in) :: unit
    call lexeme_type_write (tt%type, unit)
    write (unit, "(A)", advance="no") "'" // tt%charset1(1:tt%len1) // "'"
    write (unit, "(A)", advance="no") " '" // tt%charset2(1:tt%len2) // "'"
  end subroutine template_write

  pure function template_whitespace (chars) result (tt)
    character(*), intent(in) :: chars
    type(template_t) :: tt
    tt = template_t (WHITESPACE, chars, "", len (chars), 0)
  end function template_whitespace

  subroutine match_whitespace (tt, s, n)
    type(template_t), intent(in) :: tt
    character(*), intent(in) :: s
    integer, intent(out) :: n
    n = verify (s, tt%charset1(1:tt%len1)) - 1
    if (n < 0)  n = len (s)
  end subroutine match_whitespace

  pure function template_identifier (chars1, chars2) result (tt)
    character(*), intent(in) :: chars1, chars2
    type(template_t) :: tt
    tt = template_t (T_IDENTIFIER, chars1, chars2, len(chars1), len(chars2))
  end function template_identifier

  subroutine match_identifier (tt, s, n)
    type(template_t), intent(in) :: tt
    character(*), intent(in) :: s
    integer, intent(out) :: n
    if (verify (s(1:1), tt%charset1(1:tt%len1)) == 0) then
       n = verify (s(2:), tt%charset2(1:tt%len2))
       if (n == 0)  n = len (s)
    else
       n = 0
    end if
  end subroutine match_identifier

  pure function template_quoted (chars1, chars2) result (tt)
    character(*), intent(in) :: chars1, chars2
    type(template_t) :: tt
    tt = template_t (T_QUOTED, chars1, chars2, len (chars1), len (chars2))
  end function template_quoted

  subroutine match_quoted (tt, s, n, range)
    type(template_t), intent(in) :: tt
    character(*), intent(in) :: s
    integer, intent(out) :: n
    integer, dimension(2), intent(out) :: range
    character(tt%len1) :: ch1
    character(tt%len2) :: ch2
    integer :: i
    ch1 = tt%charset1
    if (s(1:tt%len1) == ch1) then
       ch2 = tt%charset2
       do i = tt%len1 + 1, len (s) - tt%len2 + 1
          if (s(i:i+tt%len2-1) == ch2) then
             n = i + tt%len2 - 1
             range(1) = tt%len1 + 1
             range(2) = i - 1
             return
          end if
       end do
       n = -1
       range = 0
    else
       n = 0
       range = 0
    end if
  end subroutine match_quoted

  pure function template_numeric (chars) result (tt)
    character(*), intent(in) :: chars
    type(template_t) :: tt
    tt = template_t (T_NUMERIC, chars, "", len (chars), 0)
  end function template_numeric

  subroutine match_numeric (tt, s, n)
    type(template_t), intent(in) :: tt
    character(*), intent(in) :: s
    integer, intent(out) :: n
    integer :: i, n0
    character(10), parameter :: digits = "0123456789"
    character(2), parameter :: signs = "-+"
    n = verify (s, digits) - 1
    if (n < 0) then
       n = 0
       return
    else if (s(n+1:n+1) == ".") then
       i = verify (s(n+2:), digits) - 1
       if (i < 0) then
          n = len (s)
          return
       else if (i > 0 .or. n > 0) then
          n = n + 1 + i
       end if
    end if
    n0 = n
    if (n > 0) then
       if (verify (s(n+1:n+1), tt%charset1(1:tt%len1)) == 0) then
          n = n + 1
          if (verify (s(n+1:n+1), signs) == 0)  n = n + 1
          i = verify (s(n+1:), digits) - 1
          if (i < 0) then
             n = len (s)
          else if (i == 0) then
             n = n0
          else
             n = n + i
          end if
       end if
    end if
  end subroutine match_numeric

  subroutine match_template (tt, s, n, range)
    type(template_t), intent(in) :: tt
    character(*), intent(in) :: s
    integer, intent(out) :: n
    integer, dimension(2), intent(out) :: range
    select case (tt%type)
    case (WHITESPACE)
       call match_whitespace (tt, s, n)
       range = 0
    case (T_IDENTIFIER)
       call match_identifier (tt, s, n)
       range(1) = 1
       range(2) = len_trim (s)
    case (T_QUOTED)
       call match_quoted (tt, s, n, range)
    case (T_NUMERIC)
       call match_numeric (tt, s, n)
       range(1) = 1
       range(2) = len_trim (s)
    case default
       call msg_bug ("Invalid lexeme template encountered")
    end select
  end subroutine match_template

  subroutine match (tt, s, n, range, ii)
    type(template_t), dimension(:), intent(in) :: tt
    character(*), intent(in) :: s
    integer, intent(out) :: n
    integer, dimension(2), intent(out) :: range
    integer, intent(out) :: ii
    integer :: i
    do i = 1, size (tt)
       call match_template (tt(i), s, n, range)
       if (n /= 0) then
          ii = i
          return
       end if
    end do
    n = 0
    ii = 0
  end subroutine match

  subroutine lexer_setup_init (setup, &
       comment_chars, quote_chars, quote_match, &
       single_chars, special_class, &
       keyword_list, upper_case_keywords)
    type(lexer_setup_t), intent(inout) :: setup
    character(*), intent(in) :: comment_chars
    character(*), intent(in) :: quote_chars, quote_match
    character(*), intent(in) :: single_chars
    character(*), dimension(:), intent(in) :: special_class
    type(keyword_list_t), pointer :: keyword_list
    logical, intent(in), optional :: upper_case_keywords
    integer :: n, i
    if (present (upper_case_keywords)) then
       if (upper_case_keywords) then
          setup%keyword_case = CASE_UP
       else
          setup%keyword_case = CASE_DOWN
       end if
    else
       setup%keyword_case = CASE_KEEP
    end if
    n = 1 + len (comment_chars) + len (quote_chars) + 1 &
         + len (single_chars) + size (special_class) + 1
    allocate (setup%tt(n))
    allocate (setup%type(0:n))
    n = 0
    setup%type(n) = NO_MATCH
    n = n + 1
    setup%tt(n) = template_whitespace (WHITESPACE_CHARS)
    setup%type(n) = EMPTY
    forall (i = 1:len(comment_chars))
       setup%tt(n+i) = template_quoted (comment_chars(i:i), LF)
       setup%type(n+i) = EMPTY
    end forall
    n = n + len (comment_chars)
    forall (i = 1:len(quote_chars))
       setup%tt(n+i) = template_quoted (quote_chars(i:i), quote_match(i:i))
       setup%type(n+i) = T_QUOTED
    end forall
    n = n + len (quote_chars)
    setup%tt(n+1) = template_numeric ("EeDd")
    setup%type(n+1) = T_NUMERIC
    n = n + 1
    forall (i = 1:len (single_chars))
       setup%tt(n+i) = template_identifier (single_chars(i:i), "")
       setup%type(n+i) = T_IDENTIFIER
    end forall
    n = n + len (single_chars)
    forall (i = 1:size (special_class))
       setup%tt(n+i) = template_identifier &
            (trim (special_class(i)), trim (special_class(i)))
       setup%type(n+i) = T_IDENTIFIER
    end forall
    n = n + size (special_class)
    setup%tt(n+1) = template_identifier &
         (LCLETTERS//UCLETTERS, LCLETTERS//DIGIT_CHARS//"_"//UCLETTERS)
    setup%type(n+1) = T_IDENTIFIER
    n = n + 1
    if (n /= size (setup%tt)) &
         call msg_bug ("Size mismatch in lexer setup")
    setup%keyword_list => keyword_list
  end subroutine lexer_setup_init

  subroutine lexer_setup_final (setup)
    type(lexer_setup_t), intent(inout) :: setup
    deallocate (setup%tt, setup%type)
    setup%keyword_list => null ()
  end subroutine lexer_setup_final

  subroutine lexer_setup_write (setup, unit)
    type(lexer_setup_t), intent(in) :: setup
    integer, intent(in) :: unit
    integer :: i
    write (unit, "(A)") "Lexer setup:"
    if (allocated (setup%tt)) then
       do i = 1, size (setup%tt)
          call template_write (setup%tt(i), unit)
          write (unit, '(A)', advance = "no")  " -> "
          call lexeme_type_write (setup%type(i), unit)
          write (unit, *)
       end do
    else
       write (unit, *) "[empty]"
    end if
    if (associated (setup%keyword_list)) then
       call keyword_list_write (setup%keyword_list, unit)
    end if
   end subroutine lexer_setup_write

  module subroutine lexeme_write (t, unit)
    type(lexeme_t), intent(in) :: t
    integer, intent(in) :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    select case (t%type)
    case (T_KEYWORD)
       write (u, *) "KEYWORD:    '" // char (t%s) // "'"
    case (T_IDENTIFIER)
       write (u, *) "IDENTIFIER: '" // char (t%s) // "'"
    case (T_QUOTED)
       write (u, *) "QUOTED:     '" // char (t%s) // "'"
    case (T_NUMERIC)
       write (u, *) "NUMERIC:    '" // char (t%s) // "'"
    case (UNMATCHED_QUOTE)
       write (u, *) "Unmatched quote: "// char (t%s)
    case (OVERFLOW); write (u, *) "Overflow: "// char (t%s)
    case (EMPTY);    write (u, *) "Empty lexeme"
    case (NO_MATCH); write (u, *) "No match"
    case (IO_ERROR); write (u, *) "IO error"
    case (EOF);      write (u, *) "EOF"
    case default
       write (u, *) "Error"
    end select
  end subroutine lexeme_write

  subroutine lexeme_set (t, keyword_list, s, range, type, keyword_case)
    type(lexeme_t), intent(out) :: t
    type(keyword_list_t), pointer :: keyword_list
    type(string_t), intent(in) :: s
    type(string_t) :: keyword
    integer, dimension(2), intent(in) :: range
    integer, intent(in) :: type
    integer, intent(in), optional :: keyword_case
    t%type = type
    if (present (keyword_case)) then
       select case (keyword_case)
       case (CASE_KEEP);   keyword = s
       case (CASE_UP);     keyword = upper_case (s)
       case (CASE_DOWN);   keyword = lower_case (s)
       end select
    else
       keyword = s
    end if
    if (type == T_IDENTIFIER) then
       if (associated (keyword_list)) then
          if (keyword_list_contains (keyword_list, keyword)) &
               t%type = T_KEYWORD
       end if
    end if
    select case (t%type)
    case (T_KEYWORD);  t%s = keyword
    case default;      t%s = s
    end select
    t%b = range(1)
    t%e = range(2)
  end subroutine lexeme_set

  subroutine lexeme_clear (t)
    type(lexeme_t), intent(out) :: t
    t%type = EMPTY
    t%s = ""
  end subroutine lexeme_clear

  module function lexeme_get_string (t) result (s)
    type(string_t) :: s
    type(lexeme_t), intent(in) :: t
    s = t%s
  end function lexeme_get_string

  module function lexeme_get_contents (t) result (s)
    type(string_t) :: s
    type(lexeme_t), intent(in) :: t
    s = extract (t%s, t%b, t%e)
  end function lexeme_get_contents

  module function lexeme_get_delimiters (t) result (del)
    type(string_t), dimension(2) :: del
    type(lexeme_t), intent(in) :: t
    del(1) = extract (t%s, finish = t%b-1)
    del(2) = extract (t%s, start = t%e+1)
  end function lexeme_get_delimiters

  module function lexeme_get_type (t) result (type)
    integer :: type
    type(lexeme_t), intent(in) :: t
    type = t%type
  end function lexeme_get_type

  module function lexeme_is_break (t) result (break)
    logical :: break
    type(lexeme_t), intent(in) :: t
    select case (t%type)
    case (EOF, IO_ERROR, OVERFLOW, NO_MATCH)
       break = .true.
    case default
       break = .false.
    end select
  end function lexeme_is_break

  module function lexeme_is_eof (t) result (ok)
    logical :: ok
    type(lexeme_t), intent(in) :: t
    ok = t%type == EOF
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
    call lexer_setup_init (lexer%setup, &
         comment_chars = comment_chars, &
         quote_chars = quote_chars, &
         quote_match = quote_match, &
         single_chars = single_chars, &
         special_class = special_class, &
         keyword_list = keyword_list, &
         upper_case_keywords = upper_case_keywords)
    if (present (parent))  lexer%parent => parent
    call lexer_clear (lexer)
  end subroutine lexer_init

  module subroutine lexer_clear (lexer)
    class(lexer_t), intent(inout) :: lexer
    call lexeme_clear (lexer%lexeme)
    lexer%previous_line2 = ""
    lexer%previous_line1 = ""
    lexer%current_line = ""
    lexer%lines_read = 0
    lexer%current_column = 0
    lexer%previous_column = 0
    lexer%buffer = ""
  end subroutine lexer_clear

  module subroutine lexer_final (lexer)
    class(lexer_t), intent(inout) :: lexer
    call lexer%clear ()
    call lexer_setup_final (lexer%setup)
  end subroutine lexer_final

  module subroutine lexer_assign_stream (lexer, stream)
    class(lexer_t), intent(inout) :: lexer
    type(stream_t), intent(in), target :: stream
    lexer%stream => stream
  end subroutine lexer_assign_stream

  module subroutine lex (lexeme, lexer)
    type(lexeme_t), intent(out) :: lexeme
    type(lexer_t), intent(inout) :: lexer
    integer :: iostat1, iostat2
    integer :: pos
    integer, dimension(2) :: range
    integer :: template_index, type
    if (.not. associated (lexer%stream)) &
        call msg_bug ("Lexer called without assigned stream")
    GET_LEXEME: do while (lexeme_get_type (lexer%lexeme) == EMPTY)
       if (len (lexer%buffer) /= 0) then
          iostat1 = 0
       else
          call lexer_read_line (lexer, iostat1)
       end if
       select case (iostat1)
       case (0)
          MATCH_BUFFER: do
             call match (lexer%setup%tt, char (lexer%buffer), &
                         pos, range, template_index)
             if (pos >= 0) then
                type = lexer%setup%type(template_index)
                exit MATCH_BUFFER
             else
                pos = 0
                call lexer_read_line (lexer, iostat2)
                select case (iostat2)
                case (EOF); type = UNMATCHED_QUOTE; exit MATCH_BUFFER
                case (1);   type = IO_ERROR;        exit MATCH_BUFFER
                case (2);   type = OVERFLOW;        exit MATCH_BUFFER
                end select
             end if
          end do MATCH_BUFFER
       case (EOF); type = EOF
       case (1);   type = IO_ERROR
       case (2);   type = OVERFLOW
       end select
       call lexeme_set (lexer%lexeme, lexer%setup%keyword_list, &
            extract (lexer%buffer, finish=pos), range, type, &
            lexer%setup%keyword_case)
       lexer%buffer = remove (lexer%buffer, finish=pos)
       lexer%previous_column = lexer%current_column
       lexer%current_column = lexer%current_column + pos
    end do GET_LEXEME
    lexeme = lexer%lexeme
    call lexeme_clear (lexer%lexeme)
  end subroutine lex

  subroutine lexer_read_line (lexer, iostat)
    type(lexer_t), intent(inout) :: lexer
    integer, intent(out) :: iostat
    type(string_t) :: current_line
    current_line = lexer%current_line
    call stream_get_record (lexer%stream, lexer%current_line, iostat)
    if (iostat == 0) then
       lexer%lines_read = lexer%lines_read + 1
       lexer%previous_line2 = lexer%previous_line1
       lexer%previous_line1 = current_line
       lexer%buffer = lexer%buffer // lexer%current_line // LF
       lexer%previous_column = 0
       lexer%current_column = 0
    end if
  end subroutine lexer_read_line

  module subroutine lexer_put_back (lexer, lexeme)
    type(lexer_t), intent(inout) :: lexer
    type(lexeme_t), intent(in) :: lexeme
    if (lexeme_get_type (lexer%lexeme) == EMPTY) then
       lexer%lexeme = lexeme
    else
       call msg_bug (" Lexer: lex_back fails; probably called twice")
    end if
  end subroutine lexer_put_back

  module subroutine lexer_write_setup (lexer, unit)
    type(lexer_t), intent(in) :: lexer
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    call lexer_setup_write (lexer%setup, u)
  end subroutine lexer_write_setup

  module subroutine lexer_show_location (lexer)
    type(lexer_t), intent(in) :: lexer
    type(string_t) :: loc_str
    if (associated (lexer%parent)) then
       call lexer_show_source (lexer%parent)
       call msg_message ("[includes]")
    else
       call msg_message ()
    end if
    if (associated (lexer%stream)) then
       call msg_message &
            (char (stream_get_record_info_string (lexer%stream)) // ":")
    end if
    if (lexer%lines_read >= 4)  call msg_result ("[...]")
    if (lexer%lines_read >= 3)  call msg_result (char (lexer%previous_line2))
    if (lexer%lines_read >= 2)  call msg_result (char (lexer%previous_line1))
    if (lexer%lines_read >= 1) then
       call msg_result (char (lexer%current_line))
       loc_str = repeat (" ", lexer%previous_column)
       loc_str = loc_str // "^"
       if (lexer%current_column > lexer%previous_column) then
          loc_str = loc_str &
               // repeat ("-", max (lexer%current_column &
                                    - lexer%previous_column - 1, 0)) &
               // "^"
       end if
       call msg_result (char (loc_str))
    end if
  end subroutine lexer_show_location

  recursive subroutine lexer_show_source (lexer)
    type(lexer_t), intent(in) :: lexer
    if (associated (lexer%parent)) then
       call lexer_show_source (lexer%parent)
       call msg_message ("[includes]")
    else
       call msg_message ()
    end if
    if (associated (lexer%stream)) then
       call msg_message &
            (char (stream_get_source_info_string (lexer%stream)) // ":")
    end if
  end subroutine lexer_show_source


end submodule lexers_s

