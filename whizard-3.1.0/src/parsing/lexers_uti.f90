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

module lexers_uti

  use iso_varying_string, string_t => varying_string

  use lexers

  implicit none
  private

  public :: lexer_1

contains

  subroutine lexer_1 (u)
    integer, intent(in) :: u
    type(lexer_t), target :: lexer
    type(stream_t), target :: stream
    type(string_t) :: string
    type(lexeme_t) :: lexeme
    string = "abcdefghij"
    call lexer_init (lexer, &
       comment_chars = "", &
       quote_chars = "<'""", &
       quote_match = ">'""", &
       single_chars = "?*+|=,()", &
       special_class = ["."], &
       keyword_list = null ())
    call stream_init (stream, string)
    call lexer_assign_stream (lexer, stream)
    do
       call lex (lexeme, lexer)
       call lexeme_write (lexeme, u)
       if (lexeme_is_break (lexeme))  exit
    end do
    call stream_final (stream)
    call lexer_final (lexer)
  end subroutine lexer_1


end module lexers_uti
