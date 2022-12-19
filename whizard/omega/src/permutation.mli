(* permutation.mli --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       cf. main AUTHORS file

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

    (* The argument list $\lbrack p_1;\ldots;p_n\rbrack$ must
       contain every integer from~$0$ to~$n-1$ exactly once. *)
    val of_list : int list -> t
    val of_array : int array -> t

    (* [list (of_lists l l') l = l'] *)
    val of_lists : 'a list -> 'a list -> t

    val inverse : t -> t
    val compose : t -> t -> t

    (* [compose_inv p q = compose p (inverse q)], but more efficient. *)
    val compose_inv : t -> t -> t

    (* If [p] is [of_list] $\lbrack p_1;\ldots;p_n\rbrack$, then
       [list p] $\lbrack a_1;\ldots;a_n\rbrack$ reorders the list
       $\lbrack a_1;\ldots;a_n\rbrack$
       in the sequence given by $\lbrack p_1;\ldots;p_n\rbrack$.
       Thus the $\lbrack p_1;\ldots;p_n\rbrack$ are
       \emph{not} used as a map of the indices reshuffling an array.
       Instead they denote the new positions of the elements
       of $\lbrack a_1;\ldots;a_n\rbrack$.
       However [list (inverse p)] $\lbrack a_1;\ldots;a_n\rbrack$
       is $\lbrack a_{p_1};\ldots;a_{p_n}\rbrack$, by duality. *)

    val list : t -> 'a list -> 'a list
    val array : t -> 'a array -> 'a array

    val all : int -> t list
    val even : int -> t list
    val odd : int -> t list
    val cyclic : int -> t list
    val signed : int -> (int * t) list

    (* Assuming fewer than 10 elements! *)
    val to_string : t -> string

  end

module Using_Lists : T
module Using_Arrays : T

module Default : T

module Test : functor (P : T) ->
  sig val suite : OUnit.test val time : unit -> unit end
