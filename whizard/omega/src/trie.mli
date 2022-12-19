(* trie.mli --

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

(* \thocwmodulesection{Monomorphically} *)

module type T =
  sig

    type key
    type (+'a) t
    val empty : 'a t
    val is_empty : 'a t -> bool

(* Standard trie interface: *)

    val add : key -> 'a -> 'a t -> 'a t
    val find : key -> 'a t -> 'a

(* Functionals: *)

    val remove : key -> 'a t -> 'a t
    val mem : key -> 'a t -> bool
    val map : ('a -> 'b) -> 'a t -> 'b t
    val mapi : (key -> 'a -> 'b) -> 'a t -> 'b t
    val iter : (key -> 'a -> unit) -> 'a t -> unit
    val fold : (key -> 'a -> 'b -> 'b) -> 'a t -> 'b -> 'b

(* Try to match a longest prefix and return the unmatched rest. *)

    val longest : key -> 'a t -> 'a option * key

(* Try to match a shortest prefix and return the unmatched rest. *)

    val shortest : key -> 'a t -> 'a option * key

(* \thocwmodulesection{New in O'Caml 3.08} *)

    val compare : ('a -> 'a -> int) -> 'a t -> 'a t -> int
    val equal : ('a -> 'a -> bool) -> 'a t -> 'a t -> bool

(* \thocwmodulesection{O'Mega customization}
   [export f_open f_close f_descend f_match trie] allows us to export the
   trie [trie] as source code to another programming language. *)

    val export : (int -> unit) -> (int -> unit) ->
      (int -> key -> unit) -> (int -> key -> 'a -> unit) -> 'a t -> unit

  end

(* O'Caml's [Map.S] prior to Version 3.12: *)

module type Map_S =
  sig
    type key
    type (+'a) t
    val empty: 'a t
    val is_empty: 'a t -> bool
    val add: key -> 'a -> 'a t -> 'a t
    val find: key -> 'a t -> 'a
    val remove: key -> 'a t -> 'a t
    val mem: key -> 'a t -> bool
    val iter: (key -> 'a -> unit) -> 'a t -> unit
    val map: ('a -> 'b) -> 'a t -> 'b t
    val mapi: (key -> 'a -> 'b) -> 'a t -> 'b t
    val fold: (key -> 'a -> 'b -> 'b) -> 'a t -> 'b -> 'b
    val compare: ('a -> 'a -> int) -> 'a t -> 'a t -> int
    val equal: ('a -> 'a -> bool) -> 'a t -> 'a t -> bool
  end

module Make (M : Map_S) : T with type key = M.key list
module MakeMap (M : Map_S) : Map_S with type key = M.key list

(* \thocwmodulesection{Polymorphically} *)

module type Poly =
  sig

    type ('a, 'b) t
    val empty : ('a, 'b) t

(* Standard trie interface: *)

    val add : ('a -> 'a -> int) -> 'a list -> 'b -> ('a, 'b) t -> ('a, 'b) t
    val find : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> 'b

(* Functionals: *)

    val remove : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> ('a, 'b) t
    val mem : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> bool
    val map : ('b -> 'c) -> ('a, 'b) t -> ('a, 'c) t
    val mapi : ('a list -> 'b -> 'c) -> ('a, 'b) t -> ('a, 'c) t
    val iter : ('a list -> 'b -> unit) -> ('a, 'b) t -> unit
    val fold : ('a list -> 'b -> 'c -> 'c) -> ('a, 'b) t -> 'c -> 'c

(* Try to match a longest prefix and return the unmatched rest. *)

    val longest : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> 'b option * 'a list

(* Try to match a shortest prefix and return the unmatched rest. *)

    val shortest : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> 'b option * 'a list

(* \thocwmodulesection{O'Mega customization}
   [export f_open f_close f_descend f_match trie] allows us to export the
   trie [trie] as source code to another programming language. *)

    val export : (int -> unit) -> (int -> unit) ->
      (int -> 'a list -> unit) -> (int -> 'a list -> 'b -> unit) -> ('a, 'b) t -> unit

  end

module MakePoly (M : Pmap.T) : Poly

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
