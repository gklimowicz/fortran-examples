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

submodule (parser) parser_s

  use io_units
  use diagnostics
  use md5
  use system_defs, only: DIGIT_CHARS
  use format_defs, only: FMT_19

  implicit none

contains

  subroutine token_init (token, lexeme, requested_type, key)
    type(token_t), intent(out) :: token
    type(lexeme_t), intent(in) :: lexeme
    integer, intent(in) :: requested_type
    type(string_t), intent(in) :: key
    integer :: type
    type = lexeme_get_type (lexeme)
    token%type = S_UNKNOWN
    select case (requested_type)
    case (S_LOGICAL)
       if (type == T_IDENTIFIER)  call read_logical &
            (char (lexeme_get_string (lexeme)))
    case (S_INTEGER)
       if (type == T_NUMERIC)  call read_integer &
            (char (lexeme_get_string (lexeme)))
    case (S_REAL)
       if (type == T_NUMERIC)  call read_real &
            (char (lexeme_get_string (lexeme)))
    case (S_COMPLEX)
       if (type == T_NUMERIC)  call read_complex &
            (char (lexeme_get_string (lexeme)))
    case (S_IDENTIFIER)
       if (type == T_IDENTIFIER)  call read_identifier &
            (lexeme_get_string (lexeme))
    case (S_KEYWORD)
       if (type == T_KEYWORD)  call check_keyword &
            (lexeme_get_string (lexeme), key)
    case (S_QUOTED)
       if (type == T_QUOTED)   call read_quoted &
            (lexeme_get_contents (lexeme), lexeme_get_delimiters (lexeme))
    case default
       print *, requested_type
       call msg_bug (" Invalid token type code requested by the parser")
    end select
    if (token%type /= S_UNKNOWN) then
       allocate (token%kval)
       token%kval = key
    end if
  contains
    subroutine read_logical (s)
      character(*), intent(in) :: s
      select case (s)
      case ("t", "T", "true", "TRUE", "y", "Y", "yes", "YES")
         allocate (token%lval)
         token%lval = .true.
         token%type = S_LOGICAL
      case ("f", "F", "false", "FALSE", "n", "N", "no", "NO")
         allocate (token%lval)
         token%lval = .false.
         token%type = S_LOGICAL
      end select
    end subroutine read_logical
    subroutine read_integer (s)
      character(*), intent(in) :: s
      integer :: tmp, iostat
      if (verify (s, DIGIT_CHARS) == 0) then
         read (s, *, iostat=iostat) tmp
         if (iostat == 0) then
            allocate (token%ival)
            token%ival = tmp
            token%type = S_INTEGER
         end if
      end if
    end subroutine read_integer
    subroutine read_real (s)
      character(*), intent(in) :: s
      real(default) :: tmp
      integer :: iostat
      read (s, *, iostat=iostat) tmp
      if (iostat == 0) then
         allocate (token%rval)
         token%rval = tmp
         token%type = S_REAL
      end if
    end subroutine read_real
    subroutine read_complex (s)
      character(*), intent(in) :: s
      complex(default) :: tmp
      integer :: iostat
      read (s, *, iostat=iostat) tmp
      if (iostat == 0) then
         allocate (token%cval)
         token%cval = tmp
         token%type = S_COMPLEX
      end if
    end subroutine read_complex
    subroutine read_identifier (s)
      type(string_t), intent(in) :: s
      allocate (token%sval)
      token%sval = s
      token%type = S_IDENTIFIER
    end subroutine read_identifier
    subroutine check_keyword (s, key)
      type(string_t), intent(in) :: s
      type(string_t), intent(in) :: key
      if (key == s)  token%type = S_KEYWORD
    end subroutine check_keyword
    subroutine read_quoted (s, del)
      type(string_t), intent(in) :: s
      type(string_t), dimension(2), intent(in) :: del
      allocate (token%sval, token%quote(2))
      token%sval = s
      token%quote(1) = del(1)
      token%quote(2) = del(2)
      token%type = S_QUOTED
    end subroutine read_quoted
  end subroutine token_init

  subroutine token_init_key (token, key)
    type(token_t), intent(out) :: token
    type(string_t), intent(in) :: key
    token%type = S_KEYWORD
    allocate (token%kval)
    token%kval = key
  end subroutine token_init_key

  subroutine token_final (token)
    type(token_t), intent(inout) :: token
    token%type = S_UNKNOWN
    if (associated (token%lval))  deallocate (token%lval)
    if (associated (token%ival))  deallocate (token%ival)
    if (associated (token%rval))  deallocate (token%rval)
    if (associated (token%sval))  deallocate (token%sval)
    if (associated (token%kval))  deallocate (token%kval)
    if (associated (token%quote))  deallocate (token%quote)
  end subroutine token_final

  function token_is_valid (token) result (valid)
    logical :: valid
    type(token_t), intent(in) :: token
    valid = token%type /= S_UNKNOWN
  end function token_is_valid

  subroutine token_write (token, unit)
    type(token_t), intent(in) :: token
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    select case (token%type)
    case (S_LOGICAL)
       write (u, "(L1)") token%lval
    case (S_INTEGER)
       write (u, "(I0)") token%ival
    case (S_REAL)
       write (u, "(" // FMT_19 // ")") token%rval
    case (S_COMPLEX)
       write (u, "('('," // FMT_19 // ",','," // FMT_19 // ",')')") token%cval
    case (S_IDENTIFIER)
       write (u, "(A)") char (token%sval)
    case (S_KEYWORD)
       write (u, "(A,A)") '[keyword] ' // char (token%kval)
    case (S_QUOTED)
       write (u, "(A)") &
            char (token%quote(1)) // char (token%sval) // char (token%quote(2))
    case default
       write (u, "(A)") '[empty]'
    end select
  end subroutine token_write

  module subroutine token_assign (token, token_in)
    type(token_t), intent(out) :: token
    type(token_t), intent(in) :: token_in
    token%type = token_in%type
    select case (token%type)
    case (S_LOGICAL);    allocate (token%lval);  token%lval = token_in%lval
    case (S_INTEGER);    allocate (token%ival);  token%ival = token_in%ival
    case (S_REAL);       allocate (token%rval);  token%rval = token_in%rval
    case (S_COMPLEX); allocate (token%cval);  token%cval = token_in%cval
    case (S_IDENTIFIER); allocate (token%sval);  token%sval = token_in%sval
    case (S_QUOTED);     allocate (token%sval);  token%sval = token_in%sval
       allocate (token%quote(2));  token%quote = token_in%quote
    end select
    if (token%type /= S_UNKNOWN) then
       allocate (token%kval);  token%kval = token_in%kval
    end if
  end subroutine token_assign

  module subroutine token_assign_integer (token, ival)
    type(token_t), intent(out) :: token
    integer, intent(in) :: ival
    token%type = S_INTEGER
    allocate (token%ival)
    token%ival = ival
  end subroutine token_assign_integer

  module subroutine token_assign_real (token, rval)
    type(token_t), intent(out) :: token
    real(default), intent(in) :: rval
    token%type = S_REAL
    allocate (token%rval)
    token%rval = rval
  end subroutine token_assign_real

  module subroutine token_assign_complex (token, cval)
    type(token_t), intent(out) :: token
    complex(default), intent(in) :: cval
    token%type = S_COMPLEX
    allocate (token%cval)
    token%cval = cval
  end subroutine token_assign_complex

  module subroutine token_assign_logical (token, lval)
    type(token_t), intent(out) :: token
    logical, intent(in) :: lval
    token%type = S_LOGICAL
    allocate (token%lval)
    token%lval = lval
  end subroutine token_assign_logical

  module subroutine token_assign_string (token, sval)
    type(token_t), intent(out) :: token
    type(string_t), intent(in) :: sval
    token%type = S_QUOTED
    allocate (token%sval)
    token%sval = sval
    allocate (token%quote(2));  token%quote = '"'
  end subroutine token_assign_string

  function token_get_logical (token) result (lval)
    logical :: lval
    type(token_t), intent(in) :: token
    if (associated (token%lval)) then
       lval = token%lval
    else
       call token_mismatch (token, "logical")
    end if
  end function token_get_logical

  function token_get_integer (token) result (ival)
    integer :: ival
    type(token_t), intent(in) :: token
    if (associated (token%ival)) then
       ival = token%ival
    else
       call token_mismatch (token, "integer")
    end if
  end function token_get_integer

  function token_get_real (token) result (rval)
    real(default) :: rval
    type(token_t), intent(in) :: token
    if (associated (token%rval)) then
       rval = token%rval
    else
       call token_mismatch (token, "real")
    end if
  end function token_get_real

  function token_get_cmplx (token) result (cval)
    complex(default) :: cval
    type(token_t), intent(in) :: token
    if (associated (token%cval)) then
       cval = token%cval
    else
       call token_mismatch (token, "complex")
    end if
  end function token_get_cmplx

  function token_get_string (token) result (sval)
    type(string_t) :: sval
    type(token_t), intent(in) :: token
    if (associated (token%sval)) then
       sval = token%sval
    else
       call token_mismatch (token, "string")
    end if
  end function token_get_string

  function token_get_key (token) result (kval)
    type(string_t) :: kval
    type(token_t), intent(in) :: token
    if (associated (token%kval)) then
       kval = token%kval
    else
       call token_mismatch (token, "keyword")
    end if
  end function token_get_key

  function token_get_quote (token) result (quote)
    type(string_t), dimension(2) :: quote
    type(token_t), intent(in) :: token
    if (associated (token%quote)) then
       quote = token%quote
    else
       call token_mismatch (token, "quote")
    end if
  end function token_get_quote

  subroutine token_mismatch (token, type)
    type(token_t), intent(in) :: token
    character(*), intent(in) :: type
    write (6, "(A)", advance="no")  "Token: "
    call token_write (token)
    call msg_bug (" Token type mismatch; value required as " // type)
  end subroutine token_mismatch

  recursive module subroutine parse_node_write_rec (node, unit, short, depth)
    class(parse_node_t), intent(in), target :: node
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: short
    integer, intent(in), optional :: depth
    integer :: u, d
    type(parse_node_t), pointer :: current
    u = given_output_unit (unit);  if (u < 0)  return
    d = 0;  if (present (depth))  d = depth
    call parse_node_write (node, u, short=short)
    current => node%sub_first
    do while (associated (current))
       write (u, "(A)", advance = "no")  repeat ("|  ", d)
       call parse_node_write_rec (current, unit, short, d+1)
       current => current%next
    end do
  end subroutine parse_node_write_rec

  module subroutine parse_node_write (node, unit, short)
    class(parse_node_t), intent(in) :: node
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: short
    integer :: u
    type(parse_node_t), pointer :: current
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "('+ ')", advance = "no")
    if (associated (node%rule)) then
       call syntax_rule_write (node%rule, u, &
            short=short, key_only=.true., advance=.false.)
       if (token_is_valid (node%token)) then
          write (u, "('  = ')", advance="no")
          call token_write (node%token, u)
       else if (associated (node%sub_first)) then
          write (u, "('  = ')", advance="no")
          current => node%sub_first
          do while (associated (current))
             call syntax_rule_write (current%rule, u, &
                  short=.true., key_only=.true., advance=.false.)
             current => current%next
          end do
          write (u, *)
       else
          write (u, *)
       end if
    else
       write (u, *) "[empty]"
    end if
  end subroutine parse_node_write

  recursive module subroutine parse_node_final (node, recursive)
    type(parse_node_t), intent(inout) :: node
    type(parse_node_t), pointer :: current
    logical, intent(in), optional :: recursive
    logical :: rec
    rec = .true.;  if (present (recursive))  rec = recursive
    call token_final (node%token)
    if (rec) then
       do while (associated (node%sub_first))
          current => node%sub_first
          node%sub_first => node%sub_first%next
          call parse_node_final (current)
          deallocate (current)
       end do
    end if
  end subroutine parse_node_final

  subroutine parse_node_create_leaf (node, rule, lexeme)
    type(parse_node_t), pointer :: node
    type(syntax_rule_t), intent(in), target :: rule
    type(lexeme_t), intent(in) :: lexeme
    allocate (node)
    node%rule => rule
    call token_init (node%token, lexeme, &
         syntax_rule_get_type (rule), syntax_rule_get_key (rule))
    if (.not. token_is_valid (node%token))  deallocate (node)
  end subroutine parse_node_create_leaf

  module subroutine parse_node_create_key (node, rule)
    type(parse_node_t), intent(out) :: node
    type(syntax_rule_t), intent(in), target :: rule
    node%rule => rule
    call token_init_key (node%token, syntax_rule_get_key (rule))
  end subroutine parse_node_create_key

  module subroutine parse_node_create_value (node, rule, ival, rval, cval, sval, lval)
    type(parse_node_t), intent(out) :: node
    type(syntax_rule_t), intent(in), target :: rule
    integer, intent(in), optional :: ival
    real(default), intent(in), optional :: rval
    complex(default), intent(in), optional :: cval
    type(string_t), intent(in), optional :: sval
    logical, intent(in), optional :: lval
    node%rule => rule
    call parse_node_set_value (node, ival, rval, cval, sval, lval)
  end subroutine parse_node_create_value

  module subroutine parse_node_set_value (node, ival, rval, cval, sval, lval)
    type(parse_node_t), intent(inout) :: node
    integer, intent(in), optional :: ival
    real(default), intent(in), optional :: rval
    complex(default), intent(in), optional :: cval
    type(string_t), intent(in), optional :: sval
    logical, intent(in), optional :: lval
    if (present (ival)) then
       node%token = ival
    else if (present (rval)) then
       node%token = rval
    else if (present (cval)) then
       node%token = cval
    else if (present (lval)) then
       node%token = lval
    else if (present (sval)) then
       node%token = sval
    end if
  end subroutine parse_node_set_value

  module subroutine parse_node_create_branch (node, rule)
    type(parse_node_t), pointer :: node
    type(syntax_rule_t), intent(in), target :: rule
    allocate (node)
    node%rule => rule
  end subroutine parse_node_create_branch

  module subroutine parse_node_copy (node, copy)
    class(parse_node_t), intent(in) :: node
    type(parse_node_t), pointer, intent(out) :: copy
    allocate (copy)
    select type (node)
    type is (parse_node_t)
       copy = node
    end select
    copy%next => null ()
  end subroutine parse_node_copy

  module subroutine parse_node_append_sub (node, sub)
    class(parse_node_t), intent(inout) :: node
    type(parse_node_t), pointer :: sub
    if (associated (sub)) then
       if (associated (node%sub_last)) then
          node%sub_last%next => sub
       else
          node%sub_first => sub
       end if
       node%sub_last => sub
    end if
  end subroutine parse_node_append_sub

  module subroutine parse_node_freeze_branch (node)
    type(parse_node_t), pointer :: node
    type(parse_node_t), pointer :: current
    node%n_sub = 0
    current => node%sub_first
    do while (associated (current))
       node%n_sub = node%n_sub + 1
       current => current%next
    end do
    if (node%n_sub == 0)  deallocate (node)
  end subroutine parse_node_freeze_branch

  module subroutine parse_node_replace_rule (node, rule)
    type(parse_node_t), pointer :: node
    type(syntax_rule_t), intent(in), target :: rule
    node%rule => rule
  end subroutine parse_node_replace_rule

  module subroutine parse_node_replace_last_sub (node, pn_target)
    type(parse_node_t), intent(inout), target :: node
    type(parse_node_t), intent(in), target :: pn_target
    type(parse_node_t), pointer :: current, current_copy, previous
    integer :: i
    select case (node%n_sub)
    case (1)
       allocate (current_copy)
       current_copy = pn_target
       node%sub_first => current_copy
    case (2:)
       current => node%sub_first
       allocate (current_copy)
       current_copy = current
       node%sub_first => current_copy
       previous => current_copy
       do i = 1, node%n_sub - 2
          current => current%next
          allocate (current_copy)
          current_copy = current
          previous%next => current_copy
          previous => current_copy
       end do
       allocate (current_copy)
       current_copy = pn_target
       previous%next => current_copy
    case default
       call parse_node_write (node)
       call msg_bug ("'replace_last_sub' called for non-branch parse node")
    end select
    current_copy%next => null ()
    node%sub_last => current_copy
  end subroutine parse_node_replace_last_sub

  module function parse_node_get_rule_ptr (node) result (rule)
    class(parse_node_t), intent(in) :: node
    type(syntax_rule_t), pointer :: rule
    if (associated (node%rule)) then
       rule => node%rule
    else
       rule => null ()
       call parse_node_undefined (node, "rule")
    end if
  end function parse_node_get_rule_ptr

  module function parse_node_get_n_sub (node) result (n)
    class(parse_node_t), intent(in) :: node
    integer :: n
    n = node%n_sub
  end function parse_node_get_n_sub

  module function parse_node_get_sub_ptr (node, n, tag, required) result (sub)
    class(parse_node_t), intent(in), target :: node
    type(parse_node_t), pointer :: sub
    integer, intent(in), optional :: n
    character(*), intent(in), optional :: tag
    logical, intent(in), optional :: required
    integer :: i
    sub => node%sub_first
    if (present (n)) then
       do i = 2, n
          if (associated (sub)) then
             sub => sub%next
          else
             return
          end if
       end do
    end if
    call parse_node_check (sub, tag, required)
  end function parse_node_get_sub_ptr

  module function parse_node_get_next_ptr (sub, n, tag, required) result (next)
    class(parse_node_t), intent(in), target :: sub
    type(parse_node_t), pointer :: next
    integer, intent(in), optional :: n
    character(*), intent(in), optional :: tag
    logical, intent(in), optional :: required
    integer :: i
    next => sub%next
    if (present (n)) then
       do i = 2, n
          if (associated (next)) then
             next => next%next
          else
             exit
          end if
       end do
    end if
    call parse_node_check (next, tag, required)
  end function parse_node_get_next_ptr

  module function parse_node_get_last_sub_ptr (node, tag, required) result (sub)
    type(parse_node_t), pointer :: sub
    type(parse_node_t), intent(in), target :: node
    character(*), intent(in), optional :: tag
    logical, intent(in), optional :: required
    sub => node%sub_last
    call parse_node_check (sub, tag, required)
  end function parse_node_get_last_sub_ptr

  subroutine parse_node_undefined (node, obj)
    type(parse_node_t), intent(in) :: node
    character(*), intent(in) :: obj
    call parse_node_write (node, 6)
    call msg_bug (" Parse-tree node: " // obj // " requested, but undefined")
  end subroutine parse_node_undefined

  module subroutine parse_node_check (node, tag, required)
    type(parse_node_t), pointer :: node
    character(*), intent(in), optional :: tag
    logical, intent(in), optional :: required
    if (associated (node)) then
       if (present (tag)) then
          if (parse_node_get_rule_key (node) /= tag) &
               call parse_node_mismatch (tag, node)
       end if
    else
       if (present (required)) then
          if (required) &
               call msg_bug (" Missing node, expected <" // tag // ">")
       end if
    end if
  end subroutine parse_node_check

  module subroutine parse_node_mismatch (string, parse_node)
    character(*), intent(in) :: string
    type(parse_node_t), intent(in) :: parse_node
    call parse_node_write (parse_node)
    call msg_bug (" Syntax mismatch, expected <" // string // ">.")
  end subroutine parse_node_mismatch

  module function parse_node_get_logical (node) result (lval)
    class(parse_node_t), intent(in), target :: node
    logical :: lval
    lval = token_get_logical (parse_node_get_token_ptr (node))
  end function parse_node_get_logical

  module function parse_node_get_integer (node) result (ival)
    class(parse_node_t), intent(in), target :: node
    integer :: ival
    ival = token_get_integer (parse_node_get_token_ptr (node))
  end function parse_node_get_integer

  module function parse_node_get_real (node) result (rval)
    class(parse_node_t), intent(in), target :: node
    real(default) :: rval
    rval = token_get_real (parse_node_get_token_ptr (node))
  end function parse_node_get_real

  module function parse_node_get_cmplx (node) result (cval)
    class(parse_node_t), intent(in), target :: node
    complex(default) :: cval
    cval = token_get_cmplx (parse_node_get_token_ptr (node))
  end function parse_node_get_cmplx

  module function parse_node_get_string (node) result (sval)
    class(parse_node_t), intent(in), target :: node
    type(string_t) :: sval
    sval = token_get_string (parse_node_get_token_ptr (node))
  end function parse_node_get_string

  module function parse_node_get_key (node) result (kval)
    class(parse_node_t), intent(in), target :: node
    type(string_t) :: kval
    kval = token_get_key (parse_node_get_token_ptr (node))
  end function parse_node_get_key

  module function parse_node_get_rule_key (node) result (kval)
    class(parse_node_t), intent(in), target :: node
    type(string_t) :: kval
    kval = syntax_rule_get_key (parse_node_get_rule_ptr (node))
  end function parse_node_get_rule_key

  module function parse_node_get_token_ptr (node) result (token)
    type(token_t), pointer :: token
    type(parse_node_t), intent(in), target :: node
    if (token_is_valid (node%token)) then
       token => node%token
    else
       call parse_node_undefined (node, "token")
    end if
  end function parse_node_get_token_ptr

  module function parse_node_get_md5sum (pn) result (md5sum_pn)
    character(32) :: md5sum_pn
    type(parse_node_t), intent(in) :: pn
    integer :: u
    u = free_unit ()
    open (unit = u, status = "scratch", action = "readwrite")
    call parse_node_write_rec (pn, unit=u)
    rewind (u)
    md5sum_pn = md5sum (u)
    close (u)
  end function parse_node_get_md5sum

  module subroutine parse_tree_init &
       (parse_tree, syntax, lexer, key, check_eof)
    class(parse_tree_t), intent(inout) :: parse_tree
    type(lexer_t), intent(inout) :: lexer
    type(syntax_t), intent(in), target :: syntax
    type(string_t), intent(in), optional :: key
    logical, intent(in), optional :: check_eof
    type(syntax_rule_t), pointer :: rule
    type(lexeme_t) :: lexeme
    type(parse_node_t), pointer :: node
    logical :: ok, check
    check = .true.;  if (present (check_eof)) check = check_eof
    call lexer_clear (lexer)
    if (present (key)) then
       rule => syntax_get_rule_ptr (syntax, key)
    else
       rule => syntax_get_top_rule_ptr (syntax)
    end if
    if (associated (rule)) then
       call parse_node_match_rule (node, rule, ok)
       if (ok) then
          parse_tree%root_node => node
       else
          call parse_error (rule, lexeme)
       end if
       if (check) then
          call lex (lexeme, lexer)
          if (.not. lexeme_is_eof (lexeme)) then
             call lexer_show_location (lexer)
             call msg_fatal (" Syntax error " &
                  // "(at or before the location indicated above)")
          end if
       end if
    else
       call msg_bug (" Parser failed because syntax is empty")
    end if
  contains
    recursive subroutine parse_node_match_rule (node, rule, ok)
      type(parse_node_t), pointer :: node
      type(syntax_rule_t), intent(in), target :: rule
      logical, intent(out) :: ok
      logical, parameter :: debug = .false.
      integer :: type
      if (debug)  write (6, "(A)", advance="no") "Parsing rule: "
      if (debug)  call syntax_rule_write (rule, 6)
      node => null ()
      type = syntax_rule_get_type (rule)
      if (syntax_rule_is_atomic (rule)) then
         call lex (lexeme, lexer)
         if (debug)  write (6, "(A)", advance="no") "Token: "
         if (debug)  call lexeme_write (lexeme, 6)
         call parse_node_create_leaf (node, rule, lexeme)
         ok = associated (node)
         if (.not. ok)  call lexer_put_back (lexer, lexeme)
      else
         select case (type)
         case (S_ALTERNATIVE);  call parse_alternative (node, rule, ok)
         case (S_GROUP);        call parse_group (node, rule, ok)
         case (S_SEQUENCE);     call parse_sequence (node, rule, .false., ok)
         case (S_LIST);         call parse_sequence (node, rule, .true., ok)
         case (S_ARGS);         call parse_args (node, rule, ok)
         case (S_IGNORE);       call parse_ignore (node, ok)
         end select
      end if
      if (debug) then
         if (ok) then
            write (6, "(A)", advance="no") "Matched rule: "
         else
            write (6, "(A)", advance="no") "Failed rule: "
         end if
         call syntax_rule_write (rule)
         if (associated (node)) call parse_node_write (node)
      end if
    end subroutine parse_node_match_rule
    recursive subroutine parse_alternative (node, rule, ok)
      type(parse_node_t), pointer :: node
      type(syntax_rule_t), intent(in), target :: rule
      logical, intent(out) :: ok
      integer :: i
      do i = 1, syntax_rule_get_n_sub (rule)
         call parse_node_match_rule (node, syntax_rule_get_sub_ptr (rule, i), ok)
         if (ok)  return
      end do
      ok = .false.
    end subroutine parse_alternative
    recursive subroutine parse_group (node, rule, ok)
      type(parse_node_t), pointer :: node
      type(syntax_rule_t), intent(in), target :: rule
      logical, intent(out) :: ok
      type(string_t), dimension(2) :: delimiter
      delimiter = syntax_rule_get_delimiter (rule)
      call lex (lexeme, lexer)
      if (lexeme_get_string (lexeme) == delimiter(1)) then
         call parse_node_match_rule (node, syntax_rule_get_sub_ptr (rule, 1), ok)
         if (ok) then
            call lex (lexeme, lexer)
            if (lexeme_get_string (lexeme) == delimiter(2)) then
               ok = .true.
            else
               call parse_error (rule, lexeme)
            end if
         else
            call parse_error (rule, lexeme)
         end if
      else
         call lexer_put_back (lexer, lexeme)
         ok = .false.
      end if
    end subroutine parse_group
    recursive subroutine parse_sequence (node, rule, sep, ok)
      type(parse_node_t), pointer :: node
      type(syntax_rule_t), intent(in), target :: rule
      logical, intent(in) :: sep
      logical, intent(out) :: ok
      type(parse_node_t), pointer :: current
      integer :: i, n
      logical :: opt, rep, cont
      type(string_t) :: separator
      call parse_node_create_branch (node, rule)
      if (sep)  separator = syntax_rule_get_separator (rule)
      n = syntax_rule_get_n_sub (rule)
      opt = syntax_rule_last_optional (rule)
      rep = syntax_rule_last_repetitive (rule)
      ok = .true.
      cont = .true.
      SCAN_RULE: do i = 1, n
         call parse_node_match_rule &
              (current, syntax_rule_get_sub_ptr (rule, i), cont)
         if (cont) then
            call parse_node_append_sub (node, current)
            if (sep .and. (i<n .or. rep)) then
               call lex (lexeme, lexer)
               if (lexeme_get_string (lexeme) /= separator) then
                  call lexer_put_back (lexer, lexeme)
                  cont = .false.
                  exit SCAN_RULE
               end if
            end if
         else
            if (i == n .and. opt) then
               exit SCAN_RULE
            else if (i == 1) then
               ok = .false.
               exit SCAN_RULE
            else
               call parse_error (rule, lexeme)
            end if
         end if
      end do SCAN_RULE
      if (rep) then
         do while (cont)
            call parse_node_match_rule &
                 (current, syntax_rule_get_sub_ptr (rule, n), cont)
            if (cont) then
               call parse_node_append_sub (node, current)
               if (sep) then
                  call lex (lexeme, lexer)
                  if (lexeme_get_string (lexeme) /= separator) then
                     call lexer_put_back (lexer, lexeme)
                     cont = .false.
                  end if
               end if
            else
               if (sep)  call parse_error (rule, lexeme)
            end if
         end do
      end if
      call parse_node_freeze_branch (node)
    end subroutine parse_sequence
    recursive subroutine parse_args (node, rule, ok)
      type(parse_node_t), pointer :: node
      type(syntax_rule_t), intent(in), target :: rule
      logical, intent(out) :: ok
      type(string_t), dimension(2) :: delimiter
      delimiter = syntax_rule_get_delimiter (rule)
      call lex (lexeme, lexer)
      if (lexeme_get_string (lexeme) == delimiter(1)) then
         call parse_sequence (node, rule, .true., ok)
         if (ok) then
            call lex (lexeme, lexer)
            if (lexeme_get_string (lexeme) == delimiter(2)) then
               ok = .true.
            else
               call parse_error (rule, lexeme)
            end if
         else
            call parse_error (rule, lexeme)
         end if
      else
         call lexer_put_back (lexer, lexeme)
         ok = .false.
      end if
    end subroutine parse_args
    subroutine parse_ignore (node, ok)
      type(parse_node_t), pointer :: node
      logical, intent(out) :: ok
      call lex (lexeme, lexer)
      select case (lexeme_get_type (lexeme))
      case (T_NUMERIC, T_IDENTIFIER, T_QUOTED)
         ok = .true.
      case default
         ok = .false.
      end select
      node => null ()
    end subroutine parse_ignore
    subroutine parse_error (rule, lexeme)
      type(syntax_rule_t), intent(in) :: rule
      type(lexeme_t), intent(in) :: lexeme
      character(80) :: buffer
      integer :: u, iostat
      call lexer_show_location (lexer)
      u = free_unit ()
      open (u, status = "scratch")
      write (u, "(A)", advance="no")  "Expected syntax:"
      call syntax_rule_write (rule, u)
      write (u, "(A)", advance="no")  "Found token:"
      call lexeme_write (lexeme, u)
      rewind (u)
      do
         read (u, "(A)", iostat=iostat)  buffer
         if (iostat /= 0)  exit
         call msg_message (trim (buffer))
      end do
      call msg_fatal (" Syntax error " &
           // "(at or before the location indicated above)")
    end subroutine parse_error
  end subroutine parse_tree_init

  module subroutine parse_tree_final (parse_tree)
    class(parse_tree_t), intent(inout) :: parse_tree
    if (associated (parse_tree%root_node)) then
       call parse_node_final (parse_tree%root_node)
       deallocate (parse_tree%root_node)
    end if
  end subroutine parse_tree_final

  module subroutine parse_tree_write (parse_tree, unit, verbose)
    class(parse_tree_t), intent(in) :: parse_tree
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    logical :: short
    u = given_output_unit (unit);  if (u < 0)  return
    short = .true.;  if (present (verbose))  short = .not. verbose
    write (u, "(A)") "Parse tree:"
    if (associated (parse_tree%root_node)) then
       call parse_node_write_rec (parse_tree%root_node, unit, short, 1)
    else
       write (u, *) "[empty]"
    end if
  end subroutine parse_tree_write

  module subroutine parse_tree_bug (node, keys)
    type(parse_node_t), intent(in) :: node
    character(*), intent(in) :: keys
    call parse_node_write (node)
    call msg_bug (" Inconsistency in parse tree: expected " // keys)
  end subroutine parse_tree_bug

  module function parse_tree_get_root_ptr (parse_tree) result (node)
    class(parse_tree_t), intent(in) :: parse_tree
    type(parse_node_t), pointer :: node
    node => parse_tree%root_node
  end function parse_tree_get_root_ptr

  module subroutine parse_tree_reduce (parse_tree, rule_key)
    type(parse_tree_t), intent(inout) :: parse_tree
    type(string_t), dimension(:), intent(in) :: rule_key
    type(parse_node_t), pointer :: pn
    pn => parse_tree%root_node
    if (associated (pn)) then
       call parse_node_reduce (pn, null(), null())
    end if
  contains
    recursive subroutine parse_node_reduce (pn, pn_prev, pn_parent)
      type(parse_node_t), intent(inout), pointer :: pn
      type(parse_node_t), intent(in), pointer :: pn_prev, pn_parent
      type(parse_node_t), pointer :: pn_sub, pn_sub_prev, pn_tmp
      pn_sub_prev => null ()
      pn_sub => pn%sub_first
      do while (associated (pn_sub))
         call parse_node_reduce (pn_sub, pn_sub_prev, pn)
         pn_sub_prev => pn_sub
         pn_sub => pn_sub%next
      end do
      if (parse_node_get_n_sub (pn) == 1) then
         if (matches (parse_node_get_rule_key (pn), rule_key)) then
            pn_tmp => pn
            pn => pn%sub_first
            if (associated (pn_prev)) then
               pn_prev%next => pn
            else if (associated (pn_parent)) then
               pn_parent%sub_first => pn
            else
               parse_tree%root_node => pn
            end if
            if (associated (pn_tmp%next)) then
               pn%next => pn_tmp%next
            else if (associated (pn_parent)) then
               pn_parent%sub_last => pn
            end if
            call parse_node_final (pn_tmp, recursive=.false.)
            deallocate (pn_tmp)
         end if
      end if
    end subroutine parse_node_reduce
    function matches (key, key_list) result (flag)
      logical :: flag
      type(string_t), intent(in) :: key
      type(string_t), dimension(:), intent(in) :: key_list
      integer :: i
      flag = .true.
      do i = 1, size (key_list)
         if (key == key_list(i))  return
      end do
      flag = .false.
    end function matches
  end subroutine parse_tree_reduce

  module function parse_tree_get_process_ptr (parse_tree, process) result (node)
    type(parse_node_t), pointer :: node
    type(parse_tree_t), intent(in), target :: parse_tree
    type(string_t), intent(in) :: process
    type(parse_node_t), pointer :: node_root, node_process_def
    type(parse_node_t), pointer :: node_process_phs, node_process_list
    integer :: j
    node_root => parse_tree%get_root_ptr ()
    if (associated (node_root)) then
       node_process_phs => parse_node_get_sub_ptr (node_root)
       SCAN_FILE: do while (associated (node_process_phs))
          node_process_def => parse_node_get_sub_ptr (node_process_phs)
          node_process_list => parse_node_get_sub_ptr (node_process_def, 2)
          do j = 1, parse_node_get_n_sub (node_process_list)
             if (parse_node_get_string &
                  (parse_node_get_sub_ptr (node_process_list, j)) &
                  == process) then
                node => parse_node_get_next_ptr (node_process_def)
                return
             end if
          end do
          node_process_phs => parse_node_get_next_ptr (node_process_phs)
       end do SCAN_FILE
       node => null ()
    else
       node => null ()
    end if
  end function parse_tree_get_process_ptr


end submodule parser_s

