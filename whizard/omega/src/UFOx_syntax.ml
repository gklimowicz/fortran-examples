(* vertex_syntax.ml --

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

let integer i =
  Integer i

let float x =
  Float x

let variable s =
  Variable s

let quoted s =
  Quoted s

let add e1 e2 =
  Sum (e1, e2)
    
let subtract e1 e2 =
  Difference (e1, e2)
    
let multiply e1 e2 =
  Product (e1, e2)
    
let divide e1 e2 =
  Quotient (e1, e2)
    
let power e p =
  Power (e, p)

let apply f args =
  Application (f, args)

module CSet = Sets.String_Caseless

let rec variables = function
  | Integer _ | Float _ | Quoted _ -> CSet.empty
  | Variable name -> CSet.singleton name
  | Sum (e1, e2) | Difference (e1, e2)
  | Product (e1, e2) | Quotient (e1, e2)
  | Power (e1, e2) -> CSet.union (variables e1) (variables e2)
  | Application (_, elist) ->
     List.fold_left CSet.union CSet.empty (List.map variables elist)

let rec functions = function
  | Integer _ | Float _ | Variable _ | Quoted _ -> CSet.empty
  | Sum (e1, e2) | Difference (e1, e2)
  | Product (e1, e2) | Quotient (e1, e2)
  | Power (e1, e2) -> CSet.union (functions e1) (functions e2)
  | Application (f, elist) ->
     List.fold_left CSet.union (CSet.singleton f) (List.map functions elist)
