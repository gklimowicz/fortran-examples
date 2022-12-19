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
module X = UFOx_syntax

let parse_error msg =
  raise (UFOx_syntax.Syntax_Error
	   (msg, symbol_start_pos (), symbol_end_pos ()))

let invalid_parameter_attr () =
  parse_error "invalid parameter attribute"

%}

%token < int > INT
%token < float > FLOAT
%token < string > ID QUOTED
%token PLUS MINUS TIMES POWER DIV
%token LPAREN RPAREN COMMA DOT

%token END

%left PLUS MINUS
%left TIMES DIV
%left POWER
%nonassoc UNARY

%start input
%type < UFOx_syntax.expr > input

%%

input:
 | expr END { $1 }
;

expr:
 | MINUS INT %prec UNARY  { X.integer (- $2) }
 | MINUS FLOAT %prec UNARY{ X.float (-. $2) }
 | INT             	  { X.integer $1 }
 | FLOAT           	  { X.float $1 }
 | ID              	  { X.variable $1 }
 | QUOTED             	  { X.quoted $1 }
 | expr PLUS expr  	  { X.add $1 $3 }
 | expr MINUS expr 	  { X.subtract $1 $3 }
 | expr TIMES expr 	  { X.multiply $1 $3 }
 | expr DIV expr   	  { X.divide $1 $3 }
 | PLUS expr  %prec UNARY { $2 }
 | MINUS expr %prec UNARY { X.multiply (X.integer (-1)) $2 }
 | expr POWER expr  	  { X.power $1 $3 }
 | LPAREN expr RPAREN     { $2 }
 | ID LPAREN RPAREN       { X.apply $1 [] }
 | ID LPAREN args RPAREN  { X.apply $1 $3 }
;

args:
 | expr            { [$1] }
 | expr COMMA args { $1 :: $3 }
;
