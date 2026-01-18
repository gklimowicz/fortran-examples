(* phasespace.mli --

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

module type T =
  sig
    type momentum

    type 'a t
    type 'a decay

(* Sort individual decays and complete phasespaces in a canonical order
   to determine topological equivalence classes.  *)
    val sort : ('a -> 'a -> int) -> 'a t -> 'a t
    val sort_decay : ('a -> 'a -> int) -> 'a decay -> 'a decay

(* Functionals: *)
    val map : ('a -> 'b) -> 'a t -> 'b t
    val map_decay : ('a -> 'b) -> 'a decay -> 'b decay

    val eval : ('a -> 'b) -> ('a -> 'b) -> ('a -> 'b -> 'b -> 'b) -> 'a t -> 'b t
    val eval_decay : ('a -> 'b) -> ('a -> 'b -> 'b -> 'b) -> 'a decay -> 'b decay

(* [of_momenta f1 f2 plist] constructs the phasespace parameterization
   for a process $f_1 f_2 \to X$ with flavor decoration from pairs 
   of outgoing momenta and flavors [plist] and initial flavors~$f1$
   and~$f2$ *)
    val of_momenta : 'a -> 'a -> (momentum * 'a) list -> (momentum * 'a) t
    val decay_of_momenta : (momentum * 'a) list -> (momentum * 'a) decay

    exception Duplicate of momentum
    exception Unordered of momentum
    exception Incomplete of momentum

  end

module Make (M : Momentum.T) : T with type momentum = M.t

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
