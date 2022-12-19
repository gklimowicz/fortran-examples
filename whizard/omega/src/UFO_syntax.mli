(* vertex_syntax.mli --

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

(* \thocwmodulesection{Abstract Syntax} *)

exception Syntax_Error of string * Lexing.position * Lexing.position

type name = string list

type string_atom =
  | Macro of name
  | Literal of string

type value =
  | Name of name
  | Integer of int
  | Float of float
  | Fraction of int * int
  | String of string
  | String_Expr of string_atom list
  | Empty_List
  | Name_List of name list
  | Integer_List of int list
  | String_List of string list
  | Order_Dictionary of (string * int) list
  | Coupling_Dictionary of (int * int * name) list
  | Decay_Dictionary of (name list * string) list

type attrib =
  { a_name : string;
    a_value : value }
  
type declaration =
  { name : string;
    kind : name;
    attribs : attrib list }

type t = declaration list

(* A macro expansion is encoded as a special [declaration], with
   [kind = "$"] and a single attribute.  There should not never
   be the risk of a name clash.  *)
val macro : string -> value -> declaration

val to_strings : t -> string list
