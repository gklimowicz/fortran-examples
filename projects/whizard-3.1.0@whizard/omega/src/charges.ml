(* charges.ml --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>

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

module type T = 
  sig
    type t
    val add : t -> t -> t
    val sum : t list -> t
    val is_null : t -> bool
  end

module Null : T with type t = unit = 
  struct
    type t = unit
    let add () () = ()
    let sum _ = ()
    let is_null _ = true
  end

module Z : T with type t = int = 
  struct
    type t = int
    let add = ( + )
    let sum = List.fold_left add 0
    let is_null n = (n = 0)
  end

module ZZ : T with type t = int list = 
  struct
    type t = int list
    let add = List.map2 ( + )
    let sum = function 
      | [] -> []
      | [charges] -> charges
      | charges :: rest -> List.fold_left add charges rest
    let is_null = List.for_all (fun n -> n = 0)
  end

module Rat = Algebra.Small_Rational

module Q : T with type t = Rat.t = 
  struct
    type t = Rat.t
    let add = Rat.add
    let sum = List.fold_left Rat.add Rat.null
    let is_null = Rat.is_null
  end

module QQ : T with type t = Rat.t list = 
  struct
    type t = Rat.t list
    let add = List.map2 Rat.add
    let sum = function 
      | [] -> []
      | [charges] -> charges
      | charges :: rest -> List.fold_left add charges rest
    let is_null = List.for_all Rat.is_null
  end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
