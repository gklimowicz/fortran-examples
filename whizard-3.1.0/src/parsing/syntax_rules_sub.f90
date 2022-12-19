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

submodule (syntax_rules) syntax_rules_s

  use system_defs, only: LCLETTERS, UCLETTERS, DIGIT_CHARS
  use io_units
  use diagnostics
  use ifiles, only: ifile_get_length
  use ifiles, only: line_p, line_init, line_get_string_advance, line_final

  implicit none

  character(*), parameter :: &
           UNQUOTED = "(),|_"//LCLETTERS//UCLETTERS//DIGIT_CHARS

contains

  elemental function rule_is_associated (rp) result (ok)
    logical :: ok
    type (rule_p), intent(in) :: rp
    ok = associated (rp%p)
  end function rule_is_associated

  subroutine syntax_rule_init (rule, key, type)
    type(syntax_rule_t), intent(inout) :: rule
    type(string_t), intent(in) :: key
    integer, intent(in) :: type
    rule%keyword = key
    rule%type = type
    select case (rule%type)
    case (S_GROUP)
       call syntax_rule_set_delimiter (rule)
    case (S_LIST)
       call syntax_rule_set_separator (rule)
    case (S_ARGS)
       call syntax_rule_set_delimiter (rule)
       call syntax_rule_set_separator (rule)
    end select
  end subroutine syntax_rule_init

  module function syntax_rule_get_type (rule) result (type)
    integer :: type
    type(syntax_rule_t), intent(in) :: rule
    type = rule%type
  end function syntax_rule_get_type

  module function syntax_rule_get_key (rule) result (key)
    class(syntax_rule_t), intent(in) :: rule
    type(string_t) :: key
    key = rule%keyword
  end function syntax_rule_get_key

  module function syntax_rule_get_separator (rule) result (separator)
    type(string_t) :: separator
    type(syntax_rule_t), intent(in) :: rule
    separator = rule%separator
  end function syntax_rule_get_separator

  module function syntax_rule_get_delimiter (rule) result (delimiter)
    type(string_t), dimension(2) :: delimiter
    type(syntax_rule_t), intent(in) :: rule
    delimiter = rule%delimiter
  end function syntax_rule_get_delimiter

  module function syntax_rule_get_n_sub (rule) result (n)
    integer :: n
    type(syntax_rule_t), intent(in) :: rule
    if (allocated (rule%child)) then
       n = size (rule%child)
    else
       n = 0
    end if
  end function syntax_rule_get_n_sub

  module function syntax_rule_get_sub_ptr (rule, i) result (sub)
    type(syntax_rule_t), pointer :: sub
    type(syntax_rule_t), intent(in), target :: rule
    integer, intent(in) :: i
    sub => rule%child(i)%p
  end function syntax_rule_get_sub_ptr

  subroutine syntax_rule_set_sub (rule, i, sub)
    type(syntax_rule_t), intent(inout) :: rule
    integer, intent(in) :: i
    type(syntax_rule_t), intent(in), target :: sub
    rule%child(i)%p => sub
  end subroutine syntax_rule_set_sub

  module function syntax_rule_last_optional (rule) result (opt)
    logical :: opt
    type(syntax_rule_t), intent(in) :: rule
    opt = rule%opt
  end function syntax_rule_last_optional
  module function syntax_rule_last_repetitive (rule) result (rep)
    logical :: rep
    type(syntax_rule_t), intent(in) :: rule
    rep = rule%rep
  end function syntax_rule_last_repetitive

  module function syntax_rule_is_atomic (rule) result (atomic)
    logical :: atomic
    type(syntax_rule_t), intent(in) :: rule
    select case (rule%type)
    case (S_LOGICAL, S_INTEGER, S_REAL, S_COMPLEX, S_IDENTIFIER, &
             S_KEYWORD, S_QUOTED)
       atomic = .true.
    case default
       atomic = .false.
    end select
  end function syntax_rule_is_atomic

  module subroutine syntax_rule_write (rule, unit, short, key_only, advance)
    class(syntax_rule_t), intent(in) :: rule
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: short, key_only, advance
    logical :: typ, def, adv
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    typ = .true.;  if (present (short))    typ = .not. short
    def = .true.;  if (present (key_only)) def = .not. key_only
    adv = .true.;  if (present (advance))  adv = advance
    select case (rule%type)
    case (S_UNKNOWN);    call write_atom ("???", typ)
    case (S_IGNORE);     call write_atom ("IGNORE", typ)
    case (S_LOGICAL);    call write_atom ("LOGICAL", typ)
    case (S_INTEGER);    call write_atom ("INTEGER", typ)
    case (S_REAL);       call write_atom ("REAL", typ)
    case (S_COMPLEX);    call write_atom ("COMPLEX", typ)
    case (S_IDENTIFIER); call write_atom ("IDENTIFIER", typ)
    case (S_KEYWORD);    call write_atom ("KEYWORD", typ)
    case (S_QUOTED)
       call write_quotes (typ, def, &
            del = rule%delimiter)
    case (S_SEQUENCE)
       call write_sequence ("SEQUENCE", typ, def, size (rule%child))
    case (S_GROUP)
       call write_sequence ("GROUP", typ, def, size (rule%child), &
            del = rule%delimiter)
    case (S_LIST)
       call write_sequence ("LIST", typ, def, size (rule%child), &
            sep = rule%separator)
    case (S_ARGS)
       call write_sequence ("ARGUMENTS", typ, def, size (rule%child), &
            del = rule%delimiter, &
            sep = rule%separator)
    case (S_ALTERNATIVE)
       call write_sequence ("ALTERNATIVE", typ, def, size (rule%child), &
            sep = var_str ("|"))
    end select
    if (adv)  write (u, *)
  contains
    subroutine write_type (type)
      character(*), intent(in) :: type
      character(11) :: str
      str = type
      write (u, "(1x,A)", advance="no")  str
    end subroutine write_type
    subroutine write_key
      write (u, "(1x,A)", advance="no")  char (wkey (rule))
    end subroutine write_key
    subroutine write_atom (type, typ)
      character(*), intent(in) :: type
      logical, intent(in) :: typ
      if (typ)  call write_type (type)
      call write_key
    end subroutine write_atom
    subroutine write_maybe_quoted (string)
      character(*), intent(in) :: string
      character, parameter :: q = "'"
      character, parameter :: qq = '"'
      if (verify (string, UNQUOTED) == 0) then
         write (u, "(1x,A)", advance = "no") trim (string)
      else if (verify (string, q) == 0) then
         write (u, "(1x,A)", advance = "no") qq // trim (string) // qq
      else
         write (u, "(1x,A)", advance = "no") q // trim (string) // q
      end if
    end subroutine write_maybe_quoted
    subroutine write_quotes (typ, def, del)
      logical, intent(in) :: typ, def
      type(string_t), dimension(2), intent(in) :: del
      if (typ)  call write_type ("QUOTED")
      call write_key
      if (def) then
         write (u, "(1x,'=')", advance="no")
         call write_maybe_quoted (char (del(1)))
         write (u, "(1x,A)", advance="no") "..."
         call write_maybe_quoted (char (del(2)))
      end if
    end subroutine write_quotes
    subroutine write_sequence (type, typ, def, n, del, sep)
      character(*), intent(in) :: type
      logical, intent(in) :: typ, def
      integer, intent(in) :: n
      type(string_t), dimension(2), intent(in), optional :: del
      type(string_t), intent(in), optional :: sep
      integer :: i
      if (typ)  call write_type (type)
      call write_key
      if (def) then
         write (u, "(1x,'=')", advance="no")
         if (present (del))  call write_maybe_quoted (char (del(1)))
         do i = 1, n
            if (i > 1 .and. present (sep)) &
                 call write_maybe_quoted (char (sep))
            write (u, "(1x,A)", advance="no")  &
                 char (wkey (syntax_rule_get_sub_ptr(rule, i)))
            if (i == n)  write (u, "(A)", advance="no")  trim (rule%modifier)
         end do
         if (present (del))  call write_maybe_quoted (char (del(2)))
      end if
    end subroutine write_sequence
  end subroutine syntax_rule_write

  function wkey (rule) result (string)
    type(string_t) :: string
    type(syntax_rule_t), intent(in) :: rule
    select case (rule%type)
    case (S_KEYWORD)
       if (verify (rule%keyword, UNQUOTED) == 0) then
          string = rule%keyword
       else if (scan (rule%keyword, "'") == 0) then
          string = "'" // rule%keyword // "'"
       else
          string = '"' // rule%keyword // '"'
       end if
    case default
       string = "<" // rule%keyword // ">"
    end select
  end function wkey

  subroutine syntax_rule_set_separator (rule, separator)
    type(syntax_rule_t), intent(inout) :: rule
    type(string_t), intent(in), optional :: separator
    if (present (separator)) then
       rule%separator = separator
    else
       rule%separator = ","
    end if
  end subroutine syntax_rule_set_separator

  subroutine syntax_rule_set_delimiter (rule, delimiter)
    type(syntax_rule_t), intent(inout) :: rule
    type(string_t), dimension(2), intent(in), optional :: delimiter
    if (present (delimiter)) then
       rule%delimiter = delimiter
    else
       rule%delimiter(1) = "("
       rule%delimiter(2) = ")"
    end if
  end subroutine syntax_rule_set_delimiter

  function is_modifier (string) result (ok)
    logical :: ok
    type(string_t), intent(in) :: string
    select case (char (string))
    case (" ", "?", "*", "+");  ok = .true.
    case default;               ok = .false.
    end select
  end function is_modifier

  subroutine syntax_rule_set_modifier (rule, modifier)
    type(syntax_rule_t), intent(inout) :: rule
    type(string_t), intent(in) :: modifier
    rule%modifier = char (modifier)
    select case (rule%modifier)
    case (" ")
    case ("?");  rule%opt = .true.
    case ("*");  rule%opt = .true.;  rule%rep = .true.
    case ("+");  rule%rep = .true.
    case default
       call msg_bug (" Syntax: sequence modifier '" // rule%modifier &
            // "' is not one of '+' '*' '?'")
    end select
  end subroutine syntax_rule_set_modifier

  subroutine syntax_rule_check (rule)
    type(syntax_rule_t), intent(in) :: rule
    if (rule%keyword == "")  call msg_bug ("Rule key not set")
    select case (rule%type)
    case (S_UNKNOWN);  call bug (" Undefined rule")
    case (S_IGNORE, S_LOGICAL, S_INTEGER, S_REAL, S_COMPLEX, &
             S_IDENTIFIER, S_KEYWORD)
    case (S_QUOTED)
       if (rule%delimiter(1) == "" .or. rule%delimiter(2) == "") &
            call bug (" Missing quote character(s)")
    case (S_SEQUENCE)
    case (S_GROUP)
       if (rule%delimiter(1) == "" .or. rule%delimiter(2) == "") &
            call bug (" Missing delimiter(s)")
    case (S_LIST)
       if (rule%separator == "") call bug (" Missing separator")
    case (S_ARGS)
       if (rule%delimiter(1) == "" .or. rule%delimiter(2) == "") &
            call bug (" Missing delimiter(s)")
       if (rule%separator == "") call bug (" Missing separator")
    case (S_ALTERNATIVE)
    case default
       call bug (" Undefined syntax code")
    end select
    select case (rule%type)
    case (S_SEQUENCE, S_GROUP, S_LIST, S_ARGS, S_ALTERNATIVE)
       if (allocated (rule%child)) then
          if (.not.all (rule_is_associated (rule%child))) &
               call bug (" Child rules not all associated")
       else
          call bug (" Parent rule without children")
       end if
    case default
       if (allocated (rule%child))  call bug (" Non-parent rule with children")
    end select
  contains
    subroutine bug (string)
      character(*), intent(in) :: string
      call msg_bug (" Syntax table: Rule " // char (rule%keyword) // ": " &
           // string)
    end subroutine bug
  end subroutine syntax_rule_check

  module subroutine syntax_init_from_ifile (syntax, ifile)
    type(syntax_t), intent(out), target :: syntax
    type(ifile_t), intent(in) :: ifile
    type(lexer_t) :: lexer
    type(line_p) :: line
    type(string_t) :: string
    integer :: n_token
    integer :: i
    call lexer_init (lexer, &
       comment_chars = "", &
       quote_chars = "<'""", &
       quote_match = ">'""", &
       single_chars = "?*+|=,()", &
       special_class = ["."], &
       keyword_list = null ())
    allocate (syntax%rule (ifile_get_length (ifile)))
    call line_init (line, ifile)
    do i = 1, size (syntax%rule)
       string = line_get_string_advance (line)
       call set_rule_type_and_key (syntax%rule(i), string, lexer)
    end do
    call line_init (line, ifile)
    do i = 1, size (syntax%rule)
       string = line_get_string_advance (line)
       select case (syntax%rule(i)%type)
       case (S_QUOTED, S_SEQUENCE, S_GROUP, S_LIST, S_ARGS, S_ALTERNATIVE)
          n_token = get_n_token (string, lexer)
          call set_rule_contents &
               (syntax%rule(i), syntax, n_token, string, lexer)
       end select
    end do
    call line_final (line)
    call lexer_final (lexer)
    call syntax_make_keyword_list (syntax)
    if (.not. all (syntax%rule%used)) then
       do i = 1, size (syntax%rule)
          if (.not. syntax%rule(i)%used) then
             call syntax_rule_write (syntax%rule(i), 6)
          end if
       end do
       call msg_bug (" Syntax table: unused rules")
    end if
  end subroutine syntax_init_from_ifile

  subroutine set_rule_type_and_key (rule, string, lexer)
    type(syntax_rule_t), intent(inout) :: rule
    type(string_t), intent(in) :: string
    type(lexer_t), intent(inout) :: lexer
    type(stream_t), target :: stream
    type(lexeme_t) :: lexeme
    type(string_t) :: key
    character(2) :: type
    call lexer_clear (lexer)
    call stream_init (stream, string)
    call lexer_assign_stream (lexer, stream)
    call lex (lexeme, lexer)
    type = lexeme_get_string (lexeme)
    call lex (lexeme, lexer)
    key = lexeme_get_contents (lexeme)
    call stream_final (stream)
    if (trim (key) /= "") then
       select case (type)
       case ("IG");  call syntax_rule_init (rule, key, S_IGNORE)
       case ("LO");  call syntax_rule_init (rule, key, S_LOGICAL)
       case ("IN");  call syntax_rule_init (rule, key, S_INTEGER)
       case ("RE");  call syntax_rule_init (rule, key, S_REAL)
       case ("CO");  call syntax_rule_init (rule, key, S_COMPLEX)
       case ("ID");  call syntax_rule_init (rule, key, S_IDENTIFIER)
       case ("KE");  call syntax_rule_init (rule, key, S_KEYWORD)
       case ("QU");  call syntax_rule_init (rule, key, S_QUOTED)
       case ("SE");  call syntax_rule_init (rule, key, S_SEQUENCE)
       case ("GR");  call syntax_rule_init (rule, key, S_GROUP)
       case ("LI");  call syntax_rule_init (rule, key, S_LIST)
       case ("AR");  call syntax_rule_init (rule, key, S_ARGS)
       case ("AL");  call syntax_rule_init (rule, key, S_ALTERNATIVE)
       case default
          call lexer_show_location (lexer)
          call msg_bug (" Syntax definition: unknown type '" // type // "'")
       end select
    else
       print *, char (string)
       call msg_bug (" Syntax definition: empty rule key")
    end if
  end subroutine set_rule_type_and_key

  function get_n_token (string, lexer) result (n)
    integer :: n
    type(string_t), intent(in) :: string
    type(lexer_t), intent(inout) :: lexer
    type(stream_t), target :: stream
    type(lexeme_t) :: lexeme
    integer :: i
    call lexer_clear (lexer)
    call stream_init (stream, string)
    call lexer_assign_stream (lexer, stream)
    i = 0
    do
       call lex (lexeme, lexer)
       if (lexeme_is_break (lexeme))  exit
       i = i + 1
    end do
    n = i
    call stream_final (stream)
  end function get_n_token

  module function syntax_get_rule_ptr (syntax, key) result (rule)
    type(syntax_rule_t), pointer :: rule
    type(syntax_t), intent(in), target :: syntax
    type(string_t), intent(in) :: key
    integer :: i
    do i = 1, size (syntax%rule)
       if (syntax%rule(i)%keyword == key) then
          rule => syntax%rule(i)
          return
       end if
    end do
    call msg_bug (" Syntax table: Rule " // char (key) // " not found")
  end function syntax_get_rule_ptr

  subroutine set_rule_contents (rule, syntax, n_token, string, lexer)
    type(syntax_rule_t), intent(inout) :: rule
    type(syntax_t), intent(in), target :: syntax
    integer, intent(in) :: n_token
    type(string_t), intent(in) :: string
    type(lexer_t), intent(inout) :: lexer
    type(stream_t), target :: stream
    type(lexeme_t), dimension(n_token) :: lexeme
    integer :: i, n_children
    call lexer_clear (lexer)
    call stream_init (stream, string)
    call lexer_assign_stream (lexer, stream)
    do i = 1, n_token
       call lex (lexeme(i), lexer)
    end do
    call stream_final (stream)
    n_children = get_n_children ()
    call set_delimiters
    if (n_children > 1)  call set_separator
    if (n_children > 0)  call set_children
  contains
    function get_n_children () result (n_children)
      integer :: n_children
      select case (rule%type)
      case (S_QUOTED)
         if (n_token /= 6)  call broken_rule (rule)
         n_children = 0
      case (S_GROUP)
         if (n_token /= 6)  call broken_rule (rule)
         n_children = 1
      case (S_SEQUENCE)
         if (is_modifier (lexeme_get_string (lexeme(n_token)))) then
            if (n_token <= 4)  call broken_rule (rule)
            call syntax_rule_set_modifier &
                 (rule, lexeme_get_string (lexeme(n_token)))
            n_children = n_token - 4
         else
            if (n_token <= 3) call broken_rule (rule)
            n_children = n_token - 3
         end if
      case (S_LIST)
         if (is_modifier (lexeme_get_string (lexeme(n_token)))) then
            if (n_token <= 4 .or. mod (n_token, 2) /= 1) &
                 call broken_rule (rule)
            call syntax_rule_set_modifier &
                 (rule, lexeme_get_string (lexeme(n_token)))
         else if (n_token <= 3 .or. mod (n_token, 2) /= 0) then
            call broken_rule (rule)
         end if
         n_children = (n_token - 2) / 2
      case (S_ARGS)
         if (is_modifier (lexeme_get_string (lexeme(n_token-1)))) then
            if (n_token <= 6 .or. mod (n_token, 2) /= 1) &
                 call broken_rule (rule)
            call syntax_rule_set_modifier &
                 (rule, lexeme_get_string (lexeme(n_token-1)))
         else if (n_token <= 5 .or. mod (n_token, 2) /= 0) then
            call broken_rule (rule)
         end if
         n_children = (n_token - 4) / 2
      case (S_ALTERNATIVE)
         if (n_token <= 3 .or. mod (n_token, 2) /= 0)  call broken_rule (rule)
         n_children = (n_token - 2) / 2
      end select
    end function get_n_children
    subroutine set_delimiters
      type(string_t), dimension(2) :: delimiter
      select case (rule%type)
      case (S_QUOTED, S_GROUP, S_ARGS)
         delimiter(1) = lexeme_get_contents (lexeme(4))
         delimiter(2) = lexeme_get_contents (lexeme(n_token))
         call syntax_rule_set_delimiter (rule, delimiter)
      end select
    end subroutine set_delimiters
    subroutine set_separator
      type(string_t) :: separator
      select case (rule%type)
      case (S_LIST)
         separator = lexeme_get_contents (lexeme(5))
         call syntax_rule_set_separator (rule, separator)
      case (S_ARGS)
         separator = lexeme_get_contents (lexeme(6))
         call syntax_rule_set_separator (rule, separator)
      end select
    end subroutine set_separator
    subroutine set_children
      allocate (rule%child(n_children))
      select case (rule%type)
      case (S_GROUP)
         call syntax_rule_set_sub (rule, 1, syntax_get_rule_ptr (syntax, &
              lexeme_get_contents (lexeme(5))))
      case (S_SEQUENCE)
         do i = 1, n_children
            call syntax_rule_set_sub (rule, i, syntax_get_rule_ptr (syntax, &
                 lexeme_get_contents (lexeme(i+3))))
         end do
      case (S_LIST, S_ALTERNATIVE)
         do i = 1, n_children
            call syntax_rule_set_sub (rule, i, syntax_get_rule_ptr (syntax, &
                 lexeme_get_contents (lexeme(2*i+2))))
         end do
      case (S_ARGS)
         do i = 1, n_children
            call syntax_rule_set_sub (rule, i, syntax_get_rule_ptr (syntax, &
                 lexeme_get_contents (lexeme(2*i+3))))
         end do
      end select
    end subroutine set_children
    subroutine broken_rule (rule)
      type(syntax_rule_t), intent(in) :: rule
      call lexer_show_location (lexer)
      call msg_bug (" Syntax definition: broken rule '" &
           // char (wkey (rule)) // "'")
    end subroutine broken_rule
  end subroutine set_rule_contents

  subroutine syntax_make_keyword_list (syntax)
    type(syntax_t), intent(inout), target :: syntax
    type(syntax_rule_t), pointer :: rule
    rule => syntax%rule(1)
    call rule_scan_rec (rule, syntax%keyword_list)
  contains
    recursive subroutine rule_scan_rec (rule, keyword_list)
      type(syntax_rule_t), pointer :: rule
      type(keyword_list_t), intent(inout) :: keyword_list
      integer :: i
      if (rule%used)  return
      rule%used = .true.
      select case (rule%type)
      case (S_UNKNOWN)
         call msg_bug (" Syntax: rule tree contains undefined rule")
      case (S_KEYWORD)
         call keyword_list_add (keyword_list, rule%keyword)
      end select
      select case (rule%type)
      case (S_LIST, S_ARGS)
         call keyword_list_add (keyword_list, rule%separator)
      end select
      select case (rule%type)
      case (S_GROUP, S_ARGS)
         call keyword_list_add (keyword_list, rule%delimiter(1))
         call keyword_list_add (keyword_list, rule%delimiter(2))
      end select
      select case (rule%type)
      case (S_SEQUENCE, S_GROUP, S_LIST, S_ARGS, S_ALTERNATIVE)
         if (.not. allocated (rule%child)) &
              call msg_bug (" Syntax: Non-terminal rule without children")
      case default
         if (allocated (rule%child)) &
              call msg_bug (" Syntax: Terminal rule with children")
      end select
      if (allocated (rule%child)) then
         do i = 1, size (rule%child)
            call rule_scan_rec (rule%child(i)%p, keyword_list)
         end do
      end if
    end subroutine rule_scan_rec
  end subroutine syntax_make_keyword_list

  module subroutine syntax_final (syntax)
    type(syntax_t), intent(inout) :: syntax
    if (allocated (syntax%rule))  deallocate (syntax%rule)
    call keyword_list_final (syntax%keyword_list)
  end subroutine syntax_final

  module function syntax_get_top_rule_ptr (syntax) result (rule)
    type(syntax_rule_t), pointer :: rule
    type(syntax_t), intent(in), target :: syntax
    if (allocated (syntax%rule)) then
       rule => syntax%rule(1)
    else
       rule => null ()
    end if
  end function syntax_get_top_rule_ptr

  module function syntax_get_keyword_list_ptr (syntax) result (keyword_list)
    type(keyword_list_t), pointer :: keyword_list
    type(syntax_t), intent(in), target :: syntax
    keyword_list => syntax%keyword_list
  end function syntax_get_keyword_list_ptr

  module subroutine syntax_write (syntax, unit)
    type(syntax_t), intent(in) :: syntax
    integer, intent(in), optional :: unit
    integer :: u
    integer :: i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A)") "Syntax table:"
    if (allocated (syntax%rule)) then
       do i = 1, size (syntax%rule)
          call syntax_rule_write (syntax%rule(i), u)
       end do
    else
       write (u, "(1x,A)") "[not allocated]"
    end if
    call keyword_list_write (syntax%keyword_list, u)
  end subroutine syntax_write


end submodule syntax_rules_s

