/* cascade_parser.mly --

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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

%{
open Cascade_syntax
let parse_error msg =
  raise (Syntax_Error (msg, symbol_start (), symbol_end ()))
%}

%token < string > STRING
%token < int > INT
%token LPAREN RPAREN LBRACKET RBRACKET
%token AND PLUS COLON COMMA NOT HAT
%token ONSHELL OFFSHELL GAUSS
%token END
%left AND
%left PLUS COLON COMMA
%left NOT HAT

%start main
%type < (string, int list, string) Cascade_syntax.t > main

%%

main:
    END                             { mk_true () }
  | cascades END                    { $1 }
;

cascades:
    exclusion                       { $1 }
  | vertex                          { $1 }
  | cascade                         { $1 }
  | LPAREN cascades RPAREN          { $2 }
  | cascades AND cascades           { mk_and $1 $3 }
;

exclusion:
    NOT string_list                 { mk_x_flavor $2 }
;

vertex:
    HAT string_list                 { mk_x_vertex $2 [] }
  | HAT string_list LBRACKET RBRACKET
                                    { mk_x_vertex $2 [] }
  | HAT LBRACKET string_lists RBRACKET
                                    { mk_x_vertex [] $3 }
  | HAT string_list LBRACKET string_lists RBRACKET
                                    { mk_x_vertex $2 $4 }
;

cascade:
    momentum_list                   { mk_any_flavor $1 }
  | momentum_list ONSHELL string_list
                                    { mk_on_shell $3 $1 }
  | momentum_list ONSHELL NOT string_list
                                    { mk_on_shell_not $4 $1 }
  | momentum_list OFFSHELL string_list
                                    { mk_off_shell $3 $1 }
  | momentum_list OFFSHELL NOT string_list
                                    { mk_off_shell_not $4 $1 }
  | momentum_list GAUSS string_list { mk_gauss $3 $1 }
  | momentum_list GAUSS NOT string_list
                                    { mk_gauss_not $4 $1 }
;

momentum_list:
  | momentum                        { [$1] }
  | momentum_list PLUS momentum     { $3 :: $1 }
;

momentum:
    INT                             { $1 }
;

string_list:
    STRING                          { [$1] }
  | string_list COLON STRING        { $3 :: $1 }
;

string_lists:
    string_list                     { [$1] }
  | string_lists COMMA string_list  { $3 :: $1 }
;

