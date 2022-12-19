(* thoList.mli --

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

(* [splitn n l = (hdn l, tln l)], but more efficient. *)
val hdn : int -> 'a list -> 'a list
val tln : int -> 'a list -> 'a list
val splitn : int -> 'a list -> 'a list * 'a list

(* [split_last (l @ [a]) = (l, a)] *)
val split_last : 'a list -> 'a list * 'a

(* [chop n l] chops [l] into pieces of size [n] (except for the last
   one, which contains th remainder).  *)
val chopn : int -> 'a list -> 'a list list

(* [cycle_until a l] finds a member [a] in the list [l] and returns the
   cyclically permuted list with [a] as head.  Raises [Not_found] if
   [a] is not in [l]. *)
val cycle_until : 'a -> 'a list -> 'a list

(* [cycle n l] cyclically permute the list [l] by [n >= 0]
   positions. Raises [Not_found] [List.length l > n].
   NB: [cycle n l = tln n l @ hdn n l], but more efficient. *)
val cycle : int -> 'a list -> 'a list

(* [of_subarray n m a] is $[\ocwlowerid{a.}(\ocwlowerid{n});
   \ocwlowerid{a.}(\ocwlowerid{n}+1);\ldots;
   \ocwlowerid{a.}(\ocwlowerid{m})]$.  Values of~[n] and~[m]
   out of bounds are silently shifted towards these bounds.  *)
val of_subarray : int -> int -> 'a array -> 'a list

(* [range s n m] is $[\ocwlowerid{n}; \ocwlowerid{n}+\ocwlowerid{s};
   \ocwlowerid{n}+2\ocwlowerid{s};\ldots;
   \ocwlowerid{m} - ((\ocwlowerid{m}-\ocwlowerid{n})\mod s)]$ *)
val range : ?stride:int -> int -> int -> int list

(* [enumerate s n [a1;a2;...] is [(n,a1); (n+s,a2); ...] *)
val enumerate : ?stride:int -> int -> 'a list -> (int * 'a) list

(* [alist_of_list ~predicate ~offset list] takes the elements of
   [list] that satisfy [predicate] and forms a list of pairs of
   an offset into the original [list] and the element with the
   offsets starting from [offset].  NB: the order of the returned
   alist is not specified! *)
val alist_of_list :
  ?predicate:('a -> bool) -> ?offset:int -> 'a list -> (int * 'a) list

(* Compress identical elements in a sorted list.  Identity
   is determined using the polymorphic equality function
   [Pervasives.(=)]. *)
val uniq : 'a list -> 'a list

(* Test if all members of a list are structurally identical
   (actually [homogeneous l] and [List.length (uniq l) <= 1]
   are equivalent, but the former is more efficient if a mismatch
   comes early). *)
val homogeneous : 'a list -> bool

(* If all elements of the list [l] appear exactly twice,
   [pairs l] returns a sorted list with these elements appearing
   once.  Otherwise [Invalid_argument] is raised. *)
val pairs : 'a list -> 'a list

(* [compare cmp l1 l2] compare two lists [l1] and [l2] according to
   [cmp].  [cmp] defaults to the polymorphic [Pervasives.compare].  *)
val compare : ?cmp:('a -> 'a -> int) -> 'a list -> 'a list -> int

(* Collect and count identical elements in a list.  Identity
   is determined using the polymorphic equality function
   [Pervasives.(=)].  [classify] does not assume that the list
   is sorted.  However, it is~$O(n)$ for sorted lists and~$O(n^2)$
   in the worst case.  *)
val classify : 'a list -> (int * 'a) list

(* Collect the second factors with a common first factor in lists.
   \label{ThoList.factorize} *)
val factorize : ('a * 'b) list -> ('a * 'b list) list

(* [flatmap f] is equivalent to $\ocwlowerid{flatten} \circ
   (\ocwlowerid{map}\;\ocwlowerid{f})$, but more efficient,
   because no intermediate lists are built.  Unfortunately, it is
   not tail recursive. *)
val flatmap : ('a -> 'b list) -> 'a list -> 'b list

(* [rev_flatmap f] is equivalent to $\ocwlowerid{flatten} \circ
   (\ocwlowerid{rev\_map}\;(\ocwlowerid{rev}\circ\ocwlowerid{f}))
   = \ocwlowerid{rev}\circ(\ocwlowerid{flatmap}\;\ocwlowerid{f})$,
   but more efficient, because no intermediate lists are built.
   It is tail recursive. *)
val rev_flatmap : ('a -> 'b list) -> 'a list -> 'b list

(* [clone a n] builds a list from [n] copies of the element [a]. *)
val clone : 'a -> int -> 'a list

(* [multiply n l] concatenates [n] copies of the list [l]. *)
val multiply : int -> 'a list -> 'a list

(* [filtermap f l] applies [f] to each element of [l] and drops
   the results [None]. *)
val filtermap : ('a -> 'b option) -> 'a list -> 'b list

(* [power a_list] computes the list of all sublists of [a_list],
   i.\,e.~the power set.  The elements of the sublists are \emph{not}
   required to have been sequential in [a_list]. *)
val power : 'a list -> 'a list list

(* \begin{dubious}
     Invent other names to avoid confusions with [List.fold_left2]
     and [List.fold_right2].
   \end{dubious} *)
val fold_right2 : ('a -> 'b -> 'b) -> 'a list list -> 'b -> 'b
val fold_left2 : ('b -> 'a -> 'b) -> 'b -> 'a list list -> 'b

(* [iteri f n [a;b;c]] evaluates [f n a], [f (n+1) b] and [f (n+2) c]. *)
val iteri : (int -> 'a -> unit) -> int -> 'a list -> unit
val mapi : (int -> 'a -> 'b) -> int -> 'a list -> 'b list

(* [iteri2 f n m [[aa;ab];[ba;bb]]] evaluates [f n m aa], [f n (m+1) ab],
   [f (n+1) m ba] and [f (n+1) (m+1) bb].
   NB: the nested lists need not be rectangular. *)
val iteri2 : (int -> int -> 'a -> unit) -> int -> int -> 'a list list -> unit

(* Just like [List.map3]: *)
val map3 : ('a -> 'b -> 'c -> 'd) -> 'a list -> 'b list -> 'c list -> 'd list

(* Transpose a \emph{rectangular} list of lists like a matrix.  *)
val transpose : 'a list list -> 'a list list

(* [interleave f list] walks through [list] and inserts the result
   of [f] applied to the reversed list of elements before and the
   list of elements after.  The empty lists at the beginning and
   end are included! *)
val interleave : ('a list -> 'a list -> 'a list) -> 'a list -> 'a list

(* [interleave_nearest f list] is like [interleave f list], but
   [f] looks only at the nearest neighbors. *)
val interleave_nearest : ('a -> 'a -> 'a list) -> 'a list -> 'a list

(* [partitioned_sort cmp index_sets list] sorts the sublists of [list] specified
   by the [index_sets] and the complement of their union.  \textbf{NB:} the sorting
   follows to order in the lists in [index_sets].  \textbf{NB:} the indices are
   0-based. *)
val partitioned_sort : ('a -> 'a -> int) -> int list list -> 'a list -> 'a list
exception Overlapping_indices
exception Out_of_bounds

(* [ariadne_sort cmp list] sorts [list] according to [cmp]
   (default [Pervasives.compare]) keeping track of the original order
   by a 0-based list of indices. *)
val ariadne_sort : ?cmp:('a -> 'a -> int) -> 'a list -> 'a list * int list

(* [ariadne_unsort (ariadne_sort cmp list)] returns [list]. *)
val ariadne_unsort : 'a list * int list -> 'a list

(* [lexicographic cmp list1 list2] compares [list1] and [list2]
   lexicographically. *)
val lexicographic : ?cmp:('a -> 'a -> int) -> 'a list -> 'a list -> int

(* [common l1 l2] returns the elements common to the lists [l1] and [l2].
   The lists are not required to be ordered and the result will also
   not be ordered. *)
val common : 'a list -> 'a list -> 'a list

(* [complement l1 l2] returns the list [l1] with elements of list [l2]
   removed. The lists are not required to be ordered.  Raises
   [Invalid_argument "ThoList.complement"], if a member of [l1] is not
   in [l1]. *)
val complement : 'a list -> 'a list -> 'a list

val to_string : ('a -> string) -> 'a list -> string

module Test : sig val suite : OUnit.test end
