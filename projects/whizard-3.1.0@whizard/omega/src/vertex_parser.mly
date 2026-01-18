/* vertex_parser.mly --

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

/* Right recursion is more convenient for constructing
   the value.  Since the lists will always be short,
   there is no performace or stack size reason for
   prefering left recursion. */

%{
module T = Vertex_syntax.Token
module E = Vertex_syntax.Expr
module P = Vertex_syntax.Particle
module V = Vertex_syntax.Parameter
module I = Vertex_syntax.Index
module X = Vertex_syntax.Tensor
module F = Vertex_syntax.File_Tree

let parse_error msg =
  raise (Vertex_syntax.Syntax_Error
	   (msg, symbol_start_pos (), symbol_end_pos ()))

let invalid_parameter_attr () =
  parse_error "invalid parameter attribute"

%}

%token < int > DIGIT
%token < string > CHAR
%token < string > PREFIX TOKEN
%token SUPER SUB PRIME LBRACE RBRACE LBRACKET RBRACKET
%token LPAREN RPAREN
%token COMMA
%token PLUS MINUS TIMES DIV EQUAL

%token < string > INCLUDE
%token END

%token NEUTRAL CHARGED
%token ANTI ALIAS TEX FORTRAN SPIN COLOR CHARGE MASS WIDTH
%token PARAMETER DERIVED
%token TENSOR INDEX FLAVOR LORENTZ
%token VERTEX

%left PLUS MINUS
%nonassoc NEG UPLUS
%left TIMES DIV

%start file
%type < Vertex_syntax.File_Tree.t > file

%%

file:
 | declarations END { $1 }
;

declarations:
 |                                 { [] }
 | declaration declarations        { $1 :: $2 }
;

declaration:
 | particle           { F.Particle $1 }
 | parameter          { F.Parameter $1 }
 | index              { F.Index $1 }
 | tensor             { F.Tensor $1 }
 | vertex             { let e, t = $1 in
			F.Vertex (e, t) }
 | INCLUDE            { F.Include $1 }
;

particle:
 | NEUTRAL token_arg particle_attributes
     { { P.name = P.Neutral $2; P.attr = $3 } }
 | CHARGED token_arg_pair particle_attributes
     { let p, ap = $2 in
       { P.name = P.Charged (p, ap); P.attr = $3 } }
;

expr_arg:
 | LBRACKET expr RBRACKET { $2 }
 | LBRACKET expr RBRACE   { parse_error "expected `]', found `}'" }
 | LBRACKET expr END      { parse_error "missing `]'" }
;

token_arg:
 | LBRACE scripted_token RBRACE   { $2 }
 | LBRACE scripted_token END      { parse_error "missing `}'" }
;

token_arg_pair:
 | token_arg token_arg { ($1, $2) }
;

token_list_arg:
 | LBRACE token_list RBRACE   { $2 }
 | LBRACE token_list END      { parse_error "missing `}'" }
/* This results in a reduce/reduce conflict:\hfil\goodbreak
\verb+ | LBRACE token_list RBRACKET { parse_error "expected `}', found `]'" }+ */
;

token_list_opt_arg:
 | LBRACKET token_list RBRACKET   { $2 }
 | LBRACKET token_list END        { parse_error "missing `}'" }
;

particle_attributes:
 |                                        { [ ] }
 | particle_attribute particle_attributes { $1 :: $2 }
;

particle_attribute:
 |      ALIAS      token_list_arg   		       { P.Alias $2 }
 | ANTI ALIAS      token_list_arg   		       { P.Alias $3 }
 |      TEX        token_list_arg   		       { P.TeX $2 }
 | ANTI TEX        token_list_arg   		       { P.TeX_Anti $3 }
 |      FORTRAN    token_list_arg   		       { P.Fortran $2 }
 | ANTI FORTRAN    token_list_arg   		       { P.Fortran_Anti $3 }
 |      SPIN       arg              		       { P.Spin $2 }
 |      COLOR                         token_list_arg   { P.Color ([], $2) }
 |      COLOR      token_list_opt_arg token_list_arg   { P.Color ($2, $3) }
 |      CHARGE     arg              		       { P.Charge $2 }
 |      MASS       token_list_arg   		       { P.Mass $2 }
 |      WIDTH      token_list_arg   		       { P.Width $2 }
;

parameter:
 | PARAMETER token_arg arg parameter_attributes
     { V.Parameter { V.name = $2; V.value = $3; V.attr = $4 } }
 | DERIVED   token_arg arg parameter_attributes
     { V.Derived { V.name = $2; V.value = $3; V.attr = $4 } }
;

parameter_attributes:
 |                                          { [ ] }
 | parameter_attribute parameter_attributes { $1 :: $2 }
;

parameter_attribute:
 | ALIAS   token_list_arg { V.Alias $2 }
 | TEX     token_list_arg { V.TeX $2 }
 | FORTRAN token_list_arg { V.Fortran $2 }
 | ANTI                   { invalid_parameter_attr () }
 | SPIN                   { invalid_parameter_attr () }
 | COLOR                  { invalid_parameter_attr () }
 | CHARGE                 { invalid_parameter_attr () }
 | MASS                   { invalid_parameter_attr () }
 | WIDTH                  { invalid_parameter_attr () }
;

index:
 | INDEX token_arg index_attributes { { I.name = $2; I.attr = $3 } }
;

index_attributes:
 |                                  { [ ] }
 | index_attribute index_attributes { $1 :: $2 }
;

index_attribute:
 | COLOR                      token_list_arg { I.Color ([], $2) }
 | COLOR  token_list_opt_arg  token_list_arg { I.Color ($2, $3) }
 | FLAVOR                     token_list_arg { I.Flavor ([], $2) }
 | FLAVOR token_list_opt_arg  token_list_arg { I.Flavor ($2, $3) }
 | LORENTZ                    token_list_arg { I.Lorentz $2 }
;

tensor:
 | TENSOR token_arg tensor_attributes { { X.name = $2; X.attr = $3 } }
;

tensor_attributes:
 |                                    { [ ] }
 | tensor_attribute tensor_attributes { $1 :: $2 }
;

tensor_attribute:
 | COLOR                      token_list_arg { X.Color ([], $2) }
 | COLOR  token_list_opt_arg  token_list_arg { X.Color ($2, $3) }
 | FLAVOR                     token_list_arg { X.Flavor ([], $2) }
 | FLAVOR token_list_opt_arg  token_list_arg { X.Flavor ($2, $3) }
 | LORENTZ                    token_list_arg { X.Lorentz $2 }
;

vertex:
 | VERTEX token_list_arg           { (E.integer 1, T.list $2) }
 | VERTEX expr_arg token_list_arg  { ($2, T.list $3) }
 | VERTEX expr_arg LBRACE RBRACE   { ($2, T.list []) }
 | VERTEX expr_arg LBRACE END      { parse_error "missing `}'" }
 | VERTEX not_arg_or_token_list    { parse_error "expected `[' or `{'" }
/* This results in a shift/reduce conflict:\hfil\goodbreak
\verb+ | VERTEX expr_arg LBRACE RBRACKET { parse_error "expected `}', found `]'" }+ */
;

expr:
 | integer                 	{ E.integer $1 }
 | LPAREN expr RPAREN      	{ $2 }
 | LPAREN expr RBRACKET      	{ parse_error "expected `)', found `]'" }
 | LPAREN expr RBRACE      	{ parse_error "expected `)', found `}'" }
 | LPAREN expr END      	{ parse_error "missing `)'" }
 | expr PLUS expr          	{ E.add $1 $3 }
 | expr MINUS expr         	{ E.sub $1 $3 }
 | expr TIMES expr         	{ E.mult $1 $3 }
 | expr DIV expr           	{ E.div $1 $3 }
 | bare_scripted_token arg_list { E.apply $1 $2 }
/* Making `\verb+*+' optional introduces \emph{many}
   shift/reduce and reduce/reduce conflicts:\hfil\goodbreak
\verb+ | expr expr { E.mult $1 $2 }+ */
;

arg_list:
 |                         { [] }
 | arg arg_list            { $1 :: $2 }
;

arg:
 | LBRACE expr RBRACE   { $2 }
 | LBRACE expr RBRACKET { parse_error "expected `}', found `]'" }
 | LBRACE expr END      { parse_error "missing `}'" }
;

integer:
 | DIGIT           { $1 }
 | integer DIGIT   { 10 * $1 + $2 }
;

token:
 | bare_token                              { $1 }
 | LBRACE scripted_token RBRACE            { $2 }
 | LBRACE scripted_token END               { parse_error "missing `}'" }
 | LBRACE scripted_token token_list RBRACE { T.list ($2 :: $3) }
 | LBRACE scripted_token token_list END    { parse_error "missing `}'" }
/* This results in a shift/reduce conflict because
   RBRACKET is a bare token:\hfil\goodbreak
\verb+ | LBRACE scripted_token RBRACKET     { parse_error "expected `}', found `]'" }+ */
;

token_list:
 | scripted_token            { [$1] }
 | scripted_token token_list { $1 :: $2 }
;

scripted_token:
 | prefixes token optional_scripts { T.scripted $1 $2 $3 }
;

bare_scripted_token:
 | prefixes name optional_scripts  { T.scripted $1 $2 $3 }
;

optional_scripts:
 |                { (None, None) }
 | super          { ($1, None) }
 | sub            { (None, $1) }
 | super sub      { ($1, $2) }
 | sub super      { ($2, $1) }
 | primes         { ($1, None) }
 | primes sub     { ($1, $2) }
 | sub primes     { ($2, $1) }
;

super:
 | SUPER token  { Some $2 }
 | SUPER RBRACE { parse_error "superscript can't start with `}'" }
/* This results in many reduce/reduce conflicts:\hfil\goodbreak
\verb+ | SUPER RBRACKET { parse_error "superscript can't start with `]'" }+ */
;

sub:
 | SUB token    { Some $2 }
 | SUB RBRACE   { parse_error "subscript can't start with `}'" }
/* This results in many reduce/reduce conflicts:\hfil\goodbreak
\verb+ | SUB RBRACKET { parse_error "subscript can't start with `]'" }+ */
;

prefixes:
 |                  { [] }
 | PREFIX prefixes  { $1 :: $2 }
;

primes:
 | prime_list   { Some (T.list $1) }
;

prime_list:
 | PRIME            { [T.token "\\prime"] }
 | PRIME prime_list { T.token "\\prime" :: $2 }
;

name:
 | CHAR     { T.token $1 }
 | TOKEN    { T.token $1 }
;

bare_token:
 | DIGIT    { T.digit $1 }
 | CHAR     { T.token $1 }
 | TOKEN    { T.token $1 }
 | PLUS     { T.token "+" }
 | MINUS    { T.token "-" }
 | TIMES    { T.token "*" }
 | DIV      { T.token "/" }
 | COMMA    { T.token "," }
 | LPAREN   { T.token "(" }
 | RPAREN   { T.token ")" }
;

not_arg_or_token_list:
 | DIGIT    { () }
 | CHAR     { () }
 | TOKEN    { () }
 | PLUS     { () }
 | MINUS    { () }
 | TIMES    { () }
 | DIV      { () }
 | COMMA    { () }
 | RPAREN   { () }
 | RBRACKET { () }
 | RBRACE   { () }
;
