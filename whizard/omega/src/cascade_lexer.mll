(* cascade_lexer.mll --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       Christian Speckner <cnspeckn@googlemail.com>

   WHIZARD is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   WHIZARD is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  *)

{
open Cascade_parser
let unquote s =
  String.sub s 1 (String.length s - 2)
}

let digit = ['0'-'9']
let upper = ['A'-'Z']
let lower = ['a'-'z']
let char = upper | lower
let white = [' ' '\t' '\n']

(* We use a very liberal definition of strings for flavor names. *)
rule token = parse
    white      { token lexbuf }     (* skip blanks *)
  | '%' [^'\n']* '\n'
               { token lexbuf }     (* skip comments *)
  | digit+     { INT (int_of_string (Lexing.lexeme lexbuf)) }
  | '+'        { PLUS }
  | ':'        { COLON }
  | '~'        { OFFSHELL }
  | '='        { ONSHELL }
  | '#'        { GAUSS }
  | '!'        { NOT }
  | '&' '&'?   { AND }
  | '('        { LPAREN }
  | ')'        { RPAREN }
  | '^'        { HAT }
  | ','        { COMMA }
  | '['        { LBRACKET }
  | ']'        { RBRACKET }
  | char [^ ' ' '\t' '\n' '&' '(' ')' '[' ']' ':' ',' ]*
               { STRING (Lexing.lexeme lexbuf) }
  | '"' [^ '"']* '"'
               { STRING (unquote (Lexing.lexeme lexbuf)) }
  | eof        { END }
