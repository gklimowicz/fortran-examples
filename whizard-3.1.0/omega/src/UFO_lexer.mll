(* vertex_lexer.mll --

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
open Lexing
open UFO_parser

let string_of_char c =
  String.make 1 c

let init_position fname lexbuf =
  let curr_p = lexbuf.lex_curr_p in
  lexbuf.lex_curr_p <-
    { curr_p with
      pos_fname = fname;
      pos_lnum = 1;
      pos_bol = curr_p.pos_cnum };
  lexbuf

}

let digit = ['0'-'9']
let upper = ['A'-'Z']
let lower = ['a'-'z']
let char = upper | lower
let word = char | digit | '_'
let white = [' ' '\t']
let esc = ['\'' '"' '\\']
let crlf = ['\r' '\n']
let not_crlf = [^'\r' '\n']

rule token = parse
    white             { token lexbuf }     (* skip blanks *)
  | '#' not_crlf*      { token lexbuf }     (* skip comments *)
  | crlf               { new_line lexbuf; token lexbuf }
  | "from" not_crlf*   { token lexbuf }     (* skip imports *)
  | "import" not_crlf* { token lexbuf }     (* skip imports (for now) *)
  | "try:" not_crlf*   { token lexbuf }     (* skip imports (for now) *)
  | "except" not_crlf* { token lexbuf }     (* skip imports (for now) *)
  | "pass"            { token lexbuf }     (* skip imports (for now) *)
  | '('        	      { LPAREN }
  | ')'        	      { RPAREN }
  | '{'        	      { LBRACE }
  | '}'        	      { RBRACE }
  | '['        	      { LBRACKET }
  | ']'        	      { RBRACKET }
  | '='        	      { EQUAL }
  | '+'        	      { PLUS }
  | '-'        	      { MINUS }
  | '/'        	      { DIV }
  | '.'        	      { DOT }
  | ','        	      { COMMA }
  | ':'        	      { COLON }
  | '-'? ( digit+ '.' digit* | digit* '.' digit+ )
         ( ['E''e'] '-'? digit+ )? as x
                      { FLOAT (float_of_string x) }
  | '-'? digit+ as i  { INT (int_of_string i) }
  | char word* as s   { ID s }
  | '\\' '[' (word+ as stem) ']' (word* as suffix)
                      { ID (UFO_tools.mathematica_symbol stem suffix) }
  | '\''              { let sbuf = Buffer.create 20 in
                        STRING (string1 sbuf lexbuf) }
  | '"'               { let sbuf = Buffer.create 20 in
                        STRING (string2 sbuf lexbuf) }
  | _ as c            { raise (UFO_tools.Lexical_Error
                                 ("invalid character `" ^ string_of_char c ^ "'",
                                  lexbuf.lex_start_p, lexbuf.lex_curr_p)) }
  | eof               { END }
and string1 sbuf = parse
    '\''              { Buffer.contents sbuf }
  | '\\' (esc as c)   { Buffer.add_char sbuf c; string1 sbuf lexbuf }
  | eof               { raise End_of_file }
  | '\\' '[' (word+ as stem) ']' (word* as suffix)
                      { Buffer.add_string
                          sbuf (UFO_tools.mathematica_symbol stem suffix);
                        string1 sbuf lexbuf }
  | _ as c            { Buffer.add_char sbuf c; string1 sbuf lexbuf }
and string2 sbuf = parse
    '"'               { Buffer.contents sbuf }
  | '\\' (esc as c)   { Buffer.add_char sbuf c; string2 sbuf lexbuf }
  | eof               { raise End_of_file }
  | '\\' '[' (word+ as stem) ']' (word* as suffix)
                      { Buffer.add_string
                          sbuf (UFO_tools.mathematica_symbol stem suffix);
                        string2 sbuf lexbuf }
  | _ as c            { Buffer.add_char sbuf c; string2 sbuf lexbuf }
