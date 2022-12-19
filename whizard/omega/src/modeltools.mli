(* modeltools.mli --

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

(* \thocwmodulesection{Compilation} *)

module type Flavor =
  sig
    type f
    type c
    val compare : f -> f -> int
    val conjugate : f -> f
  end

module type Fusions =
  sig
    type t
    type f
    type c
    val fuse2 : t -> f -> f -> (f * c Coupling.t) list
    val fuse3 : t -> f -> f -> f -> (f * c Coupling.t) list
    val fuse : t -> f list -> (f * c Coupling.t) list
    val of_vertices :
        (((f * f * f) * c Coupling.vertex3 * c) list
           * ((f * f * f * f) * c Coupling.vertex4 * c) list
           * (f list * c Coupling.vertexn * c) list) -> t
  end

module Fusions : functor (F : Flavor) ->
  Fusions with type f = F.f and type c = F.c

(* \thocwmodulesection{Coupling Constants} *)

(* There is no [Model.constant_of_string] function, but we can
   construct one by inverting [Model.constant_symbol] on the set
   of all coupling constants appearing in the vertices. *)

module type Constant =
  sig
    type t
    val of_string : string -> t
  end

module Constant : functor (M : Model.T) -> Constant with type t = M.constant

(* \thocwmodulesection{Mutable Models} *)

module Mutable : functor (FGC : sig type f and g and c end) ->
  Model.Mutable with type flavor = FGC.f and type gauge = FGC.g 
  and type constant = FGC.c

module Static (M : Model.T) : Model.Mutable

(* \thocwmodulesection{Topology Only} *)

module Topology (M : Model.T) : Model.T
  with type flavor = M.flavor
  and type gauge = M.gauge
  and type constant = M.constant

module Topology3 (M : Model.T) : Model.T
  with type flavor = M.flavor
  and type gauge = M.gauge
  and type constant = M.constant
