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

type expr =
  | Integer of int
  | Float of float
  | Variable of string
  | Quoted of string
  | Sum of expr * expr
  | Difference of expr * expr
  | Product of expr * expr
  | Quotient of expr * expr
  | Power of expr * expr
  | Application of string * expr list

val integer : int -> expr
val float : float -> expr
val variable : string -> expr
val quoted : string -> expr
val add : expr -> expr -> expr
val subtract : expr -> expr -> expr
val multiply : expr -> expr -> expr
val divide : expr -> expr -> expr
val power : expr -> expr -> expr
val apply : string -> expr list -> expr

(* Return the sets of variable and function names referenced
   in the expression. *)
val variables : expr -> Sets.String_Caseless.t
val functions : expr -> Sets.String_Caseless.t
