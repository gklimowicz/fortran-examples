(* cascade.mli --

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

    type constant
    type flavor
    type p

    type t
    val of_string_list : int -> string list -> t
    val to_string : t -> string

(* An opaque type that describes the set of all constraints on an amplitude
   and how to construct it from a cascade description. *)
    type selectors
    val to_selectors : t -> selectors

(* Don't throw anything away: *)
    val no_cascades : selectors

(* [select_wf s is_timelike f p ps] returns [true] iff either
   \begin{itemize}
     \item the flavor [f] and momentum [p] match the selection [s] or
     \item \emph{all} combinations of the momenta in [ps]
       are compatible, i.\,e.~$\pm\sum p_i\leq q$.
    \end{itemize}
    The latter test is only required in theories with quartic
    or higher vertices, where [ps] will be the list of all
    incoming momenta in a fusion.  [is_timelike] is required
    to determine, whether particles and anti-particles should
    be distinct. *)
    val select_wf : selectors -> (p -> bool) -> flavor -> p -> p list -> bool

(* [select_p s p ps] same as [select_wf s f p ps], but ignores the flavor [f] *)
    val select_p : selectors -> p -> p list -> bool

(* [on_shell s p] *)
    val on_shell : selectors -> flavor -> p -> bool

(* [is_gauss s p] *)
    val is_gauss : selectors -> flavor -> p -> bool

    val select_vtx : selectors -> constant Coupling.t ->
      flavor -> flavor list -> bool

(* [partition s] returns a partition of the external particles that can not
   be reordered without violating the cascade constraints. *)
    val partition : selectors -> int list list

(* Diagnostics: *)
    val description : selectors -> string option

  end

module Make (M : Model.T) (P : Momentum.T) :
    T with type flavor = M.flavor
       and type constant = M.constant
       and type p = P.t

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)

