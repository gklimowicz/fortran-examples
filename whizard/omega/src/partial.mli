(* partial.mli --

   Copyright (C) 1999-2015 by

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

(* Partial maps that are constructed from assoc lists. *)

module type T =
  sig

    (* The domain of the map.
       It needs to be compatible with [Map.OrderedType.t] *)
    type domain

    (* The codomain ['a] can be anything we want. *)
    type 'a t

    (* A list of argument-value pairs is mapped to a partial map.
       If an argument appears twice, the later value takes
       precedence. *)
    val of_list : (domain * 'a) list -> 'a t

    (* Two lists of arguments and values (both must have the
       same length) are mapped to a partial map.  Again the
       later value takes precedence. *)
    val of_lists : domain list -> 'a list -> 'a t

    (* If domain and codomain disagree, we must raise an exception
       or provide a fallback. *)
    exception Undefined of domain
    val apply : 'a t -> domain -> 'a
    val apply_with_fallback : (domain -> 'a) -> 'a t -> domain -> 'a

    (* Iff domain and codomain of the map agree, we can
       fall back to the identity map. *)
    val auto : domain t -> domain -> domain

  end

module Make : functor (D : Map.OrderedType) -> T with type domain = D.t
module Test : sig val suite : OUnit.test end
