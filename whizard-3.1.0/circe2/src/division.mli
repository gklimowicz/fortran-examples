(* circe2/division.mli --  *)
(* Copyright (C) 2001-2022 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
   Circe2 is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by 
   the Free Software Foundation; either version 2, or (at your option)
   any later version.
   Circe2 is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  *)  

(* We have divisions ([Mono]) and divisions of divisions ([Poly]).
   Except for creation, they share the same interface ([T]), which
   can be used as a signature for functor arguments.  In particular,
   both kinds of divisions can be used with the [Grid.Make] functor. *)

module type T =
  sig

    type t

    (* Copy a division, allocating fresh arrays with identical contents. *)
    val copy : t -> t

    (* Using~$\{x_0,x_1,\ldots,x_n\}$, find~$i$, such that~$x_i\le x<x_{i+1}$.
       We need to export this, if we want to maintain additional histograms
       in user modules.  *)
    val find : t -> float -> int

    (* [record d x f] records the value~$f$ at coordinate~$x$.  NB: this
       function modifies~[d]. *)
    val record : t -> float -> float -> unit

    (* VEGAS style rebinning.  The default values for [power] and both
        [fixed_min], [fixed_max] are~$1.5$ and~[false] respectively.  *)
    val rebin : ?power:float -> ?fixed_min:bool -> ?fixed_max:bool -> t -> t

    (* $J^*(y)$
       \begin{dubious}
         Should this include the $1/\Delta y$?
       \end{dubious} *)
    val caj : t -> float -> float

    val n_bins : t -> int
    val bins : t -> float array
    val to_channel : out_channel -> t -> unit

  end

exception Above_max of float * (float * float) * int
exception Below_min of float * (float * float) * int
exception Out_of_range of float * (float * float)
exception Rebinning_failure of string

(* \subsubsection{Primary Divisions} *)

module type Mono =
  sig
    include T

    (* [create bias n x_min x_max] creates a division with~$n$
       equidistant bins spanning $[x_{\min},x_{\max}]$. The [bias]
       is a function that is multiplied with the weights for
       VEGAS/VAMP rebinning.  It can be used to highlight the
       regions of phasespace that are expected to be most relevant
       in applications.  The default is [fun x -> 1.0], of course. *)
    val create : ?bias:(float -> float) -> int -> float -> float -> t

  end

module Mono : Mono

(* \subsubsection{Polydivisions} *)

module type Poly =
  sig

    module M : Diffmaps.Real

    include T

    (* [create n x_min x_max intervals] creates a polydivision of the
       interval from [x_min] to [x_max] described by the list of [intervals],
       filling the gaps among intervals and between the intervals and the
       outer borders with an unmapped divisions with [n] bins each.  *)
    val create : ?bias:(float -> float) ->
      (int * M.t) list -> int -> float -> float -> t

  end

module Make_Poly (M : Diffmaps.Real) : Poly with module M = M

module Poly : Poly

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)


