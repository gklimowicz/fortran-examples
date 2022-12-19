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

module parser_uti

  use syntax_rules

  use parser

  implicit none
  private

  public :: parse_1

contains

  subroutine parse_1 (u)
    use ifiles
    use lexers
    integer, intent(in) :: u

    type(ifile_t) :: ifile
    type(syntax_t), target :: syntax
    type(lexer_t) :: lexer
    type(stream_t), target :: stream
    type(parse_tree_t), target :: parse_tree

    write (u, "(A)")  "* Test output: Parsing"
    write (u, "(A)")  "*   Purpose: test parse routines"
    write (u, "(A)")

    call ifile_append (ifile, "SEQ expr = term addition*")
    call ifile_append (ifile, "SEQ addition = plus_or_minus term")
    call ifile_append (ifile, "SEQ term = factor multiplication*")
    call ifile_append (ifile, "SEQ multiplication = times_or_over factor")
    call ifile_append (ifile, "SEQ factor = atom exponentiation*")
    call ifile_append (ifile, "SEQ exponentiation = '^' atom")
    call ifile_append (ifile, "ALT atom = real | delimited_expr")
    call ifile_append (ifile, "GRO delimited_expr = ( expr )")
    call ifile_append (ifile, "ALT plus_or_minus = '+' | '-'")
    call ifile_append (ifile, "ALT times_or_over = '*' | '/'")
    call ifile_append (ifile, "KEY '+'")
    call ifile_append (ifile, "KEY '-'")
    call ifile_append (ifile, "KEY '*'")
    call ifile_append (ifile, "KEY '/'")
    call ifile_append (ifile, "KEY '^'")
    call ifile_append (ifile, "REA real")

    write (u, "(A)")  "* File contents (syntax definition):"
    call ifile_write (ifile, u)
    write (u, "(A)")  "EOF"
    write (u, "(A)")

    call syntax_init (syntax, ifile)
    call ifile_final (ifile)
    call syntax_write (syntax, u)
    write (u, "(A)")

    call lexer_init (lexer, &
         comment_chars = "", &
         quote_chars = "'", &
         quote_match = "'", &
         single_chars = "+-*/^()", &
         special_class = [""] , &
         keyword_list = syntax_get_keyword_list_ptr (syntax))
    call lexer_write_setup (lexer, u)
    write (u, "(A)")

    call ifile_append (ifile, "(27+8^3-2/3)*(4+7)^2*99")
    write (u, "(A)")  "* File contents (input file):"
    call ifile_write (ifile, u)
    write (u, "(A)")  "EOF"
    print *

    call stream_init (stream, ifile)
    call lexer_assign_stream (lexer, stream)
    call parse_tree_init (parse_tree, syntax, lexer)
    call stream_final (stream)
    call parse_tree_write (parse_tree, u, .true.)
    print *

    write (u, "(A)")  "* Cleanup, everything should now be empty:"
    write (u, "(A)")

    call parse_tree_final (parse_tree)
    call parse_tree_write (parse_tree, u, .true.)
    write (u, "(A)")

    call lexer_final (lexer)
    call lexer_write_setup (lexer, u)
    write (u, "(A)")

    call ifile_final (ifile)
    write (u, "(A)")  "* File contents:"
    call ifile_write (ifile, u)
    write (u, "(A)")  "EOF"
    write (u, "(A)")

    call syntax_final (syntax)
    call syntax_write (syntax, u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: parser_1"

  end subroutine parse_1

end module parser_uti
