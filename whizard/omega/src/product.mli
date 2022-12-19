(* product.mli --

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

(* \thocwmodulesection{Lists}
   Since April 2001, we preserve lexicographic ordering.  *)

val fold2 : ('a -> 'b -> 'c -> 'c) -> 'a list -> 'b list -> 'c -> 'c
val fold3 : ('a -> 'b -> 'c -> 'd -> 'd) -> 'a list -> 'b list -> 'c list -> 'd -> 'd
val fold : ('a list -> 'b -> 'b) -> 'a list list -> 'b -> 'b

val list2 : ('a -> 'b -> 'c) -> 'a list -> 'b list -> 'c list
val list3 : ('a -> 'b -> 'c -> 'd) -> 'a list -> 'b list -> 'c list -> 'd list
val list : ('a list -> 'b) -> 'a list list -> 'b list

(* Suppress all [None] in the results. *)
val list2_opt :
  ('a -> 'b -> 'c option) -> 'a list -> 'b list -> 'c list
val list3_opt :
  ('a -> 'b -> 'c -> 'd option) -> 'a list -> 'b list -> 'c list -> 'd list
val list_opt :
  ('a list -> 'b option) -> 'a list list -> 'b list

val power : int -> 'a list -> 'a list list

val thread : 'a list list -> 'a list list

(* \thocwmodulesection{Sets} *)

(* ['a_set] is actually ['a set] for a suitable [set], but this
   relation can not be expressed polymorphically (in [set]) in O'Caml.
   The two sets can be of different type, but we provide a symmetric
   version as syntactic sugar. *)

type 'a set

type ('a, 'a_set, 'b) fold = ('a -> 'b -> 'b) -> 'a_set -> 'b -> 'b
type ('a, 'a_set, 'b, 'b_set, 'c) fold2 =
    ('a -> 'b -> 'c -> 'c) -> 'a_set -> 'b_set -> 'c -> 'c

val outer : ('a, 'a_set, 'c) fold -> ('b, 'b_set, 'c) fold ->
  ('a, 'a_set, 'b, 'b_set, 'c) fold2
val outer_self : ('a, 'a_set, 'b) fold -> ('a, 'a_set, 'a, 'a_set, 'b) fold2

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
