(* powSet.mli --

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

(* Manipulate the power set, i.\,e.~the set of all subsets, of an
   set [Ordered_Type].  The concrete order is actually irrelevant, we just
   need it to construct [Set.S]s in the implementation.
   In fact, what we are implementating is the \textit{free semilattice}
   generated from the set of subsets of [Ordered_Type], where the join
   operation is the set union.

   The non trivial operation is [basis], which takes a set of subsets
   and returns the smallest set of disjoint subsets from which the argument
   can be reconstructed by forming unions.
   It is used in O'Mega for finding coarsest partitions of sets of
   partiticles.

   \begin{dubious}
     Eventually, this could be generalized from \textit{power set} or
     \textit{semi lattice} to \textit{lattice} with a notion of subtraction.
   \end{dubious} *)

module type Ordered_Type =
  sig
    type t
    val compare : t -> t -> int

    (* Debugging \ldots *)
    val to_string : t -> string
  end

module type T =
  sig
    type elt
    type t

    val empty : t
    val is_empty : t -> bool

    (* Set union (a.\,k.\,a.~join).  *)
    val union : t list -> t

    (* Construct the abstract type from a list of subsets represented as
       lists and the inverse operation. *)
    val of_lists : elt list list -> t
    val to_lists : t -> elt list list

    (* The smallest set of disjoint subsets that generates the given subset. *)
    val basis : t -> t

    (* Debugging \ldots *)
    val to_string : t -> string
  end

module Make (E : Ordered_Type) : T with type elt = E.t


(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
