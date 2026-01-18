(* omega.mli --

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
    val main : unit -> unit

(* \begin{dubious}
     This used to be only intended for debugging O'Giga,
     but might live longer \ldots
   \end{dubious} *)
    type flavor
    val diagrams : flavor -> flavor -> flavor list ->
      ((flavor * Momentum.Default.t) *
         (flavor * Momentum.Default.t,
          flavor * Momentum.Default.t) Tree.t) list
  end


(* Wrap the two instances of [Fusion.Maker] for
   amplitudes and phase space into a single functor to
   make sure that the Dirac and Majorana versions match.
   Don't export the slightly unsafe
   [module Make (FM : Fusion.Maker) (PM : Fusion.Maker)
    (TM : Target.Maker) (M : Model.T) : T with type flavor = M.flavor]. *)

module Binary (TM : Target.Maker) (M : Model.T) : T with type flavor = M.flavor
module Binary_Majorana (TM : Target.Maker) (M : Model.T) : T with type flavor = M.flavor
   
module Mixed23 (TM : Target.Maker) (M : Model.T) : T with type flavor = M.flavor
module Mixed23_Majorana (TM : Target.Maker) (M : Model.T) : T with type flavor = M.flavor
module Mixed23_Majorana_vintage (TM : Target.Maker) (M : Model.T) : T with type flavor = M.flavor

module Nary (TM : Target.Maker) (M : Model.T) : T with type flavor = M.flavor
module Nary_Majorana (TM : Target.Maker) (M : Model.T) : T with type flavor = M.flavor
