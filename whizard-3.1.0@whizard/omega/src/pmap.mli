(* pmap.mli --

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

(* Module [Pmap]: association tables over a polymorphic
   type\footnote{Extension of code \textcopyright~1996 by Xavier Leroy}. *)

module type T =
  sig
    type ('key, 'a) t
    val empty : ('key, 'a) t
    val is_empty : ('key, 'a) t -> bool
    val singleton : 'key -> 'a -> ('key, 'a) t
    val add : ('key -> 'key -> int) -> 'key -> 'a -> ('key, 'a) t -> ('key, 'a) t
    val update : ('key -> 'key -> int) -> ('a -> 'a -> 'a) ->
      'key -> 'a -> ('key, 'a) t -> ('key, 'a) t
    val cons : ('key -> 'key -> int) -> ('a -> 'a -> 'a option) ->
      'key -> 'a -> ('key, 'a) t -> ('key, 'a) t
    val find : ('key -> 'key -> int) -> 'key -> ('key, 'a) t -> 'a
    val find_opt : ('key -> 'key -> int) -> 'key -> ('key, 'a) t -> 'a option
    val choose : ('key, 'a) t -> 'key * 'a
    val choose_opt : ('key, 'a) t -> ('key * 'a) option
    val uncons : ('key, 'a) t -> 'key * 'a * ('key, 'a) t
    val uncons_opt : ('key, 'a) t -> ('key * 'a * ('key, 'a) t) option
    val elements : ('key, 'a) t -> ('key * 'a) list
    val mem :  ('key -> 'key -> int) -> 'key -> ('key, 'a) t -> bool
    val remove : ('key -> 'key -> int) -> 'key -> ('key, 'a) t -> ('key, 'a) t
    val union : ('key -> 'key -> int) -> ('a -> 'a -> 'a) ->
      ('key, 'a) t -> ('key, 'a) t -> ('key, 'a) t
    val compose : ('key -> 'key -> int) -> ('a -> 'a -> 'a option) ->
      ('key, 'a) t -> ('key, 'a) t -> ('key, 'a) t
    val iter : ('key -> 'a -> unit) -> ('key, 'a) t -> unit
    val map : ('a -> 'b) -> ('key, 'a) t -> ('key, 'b) t
    val mapi : ('key -> 'a -> 'b) -> ('key, 'a) t -> ('key, 'b) t
    val fold : ('key -> 'a -> 'b -> 'b) -> ('key, 'a) t -> 'b -> 'b
    val compare : ('key -> 'key -> int) -> ('a -> 'a -> int) ->
      ('key, 'a) t -> ('key, 'a) t -> int
    val canonicalize : ('key -> 'key -> int) -> ('key, 'a) t -> ('key, 'a) t
  end

(* Balanced trees: logarithmic access, but representation not unique. *)

module Tree : T

(* Sorted lists: representation unique, but linear access. *)

module List : T

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
