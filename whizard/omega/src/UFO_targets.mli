(* UFO_targets.mli --

   Copyright (C) 1999-2017 by

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

(* \thocwmodulesection{Generating Code for UFO Lorentz Structures} *)

module type T =
  sig

    (* [lorentz ff name spins lorentz] writes the Fortran code
       implementing the fusion corresponding to the Lorentz
       structure [lorentz] to [ff].
       NB: The [spins : int list] element of [UFO.Lorentz.t]
       from the UFO file is \emph{not} sufficient
       to determine the domain and codomain of the function.  We
       had to inspect the flavors, where the Lorentz structure
       is referenced to heuristically compute the [spins]
       as a [Coupling.lorentz array] . *)
    val lorentz :
      Format_Fortran.formatter -> string ->
      Coupling.lorentz array -> UFO_Lorentz.t -> unit

    val propagator :
      Format_Fortran.formatter -> string -> string -> string list ->
      Coupling.lorentz * Coupling.lorentz ->
      UFO_Lorentz.t -> UFO_Lorentz.t -> unit

    (* [fusion_name name perm cc_list] forms a name for the fusion
       [name] with the permutations [perm] and charge conjugations
       applied to the fermion lines [cc_list]. *)
    val fusion_name :
      string -> Permutation.Default.t -> Coupling.fermion_lines -> string

    (* [fuse c v s fl g wfs ps fusion]
       fuses the wavefunctions named [wfs] with momenta named [ps]
       using the vertex named [v] with legs reordered according to [fusion].
       The overall coupling constant named [g] is multiplied by the rational
       coefficient [c].  The list of spins [s] and the fermion
       lines [fl] are used for selecting the appropriately
       transformed version of the vertex [v]. *)
    val fuse :
      Algebra.QC.t -> string ->
      Coupling.lorentzn -> Coupling.fermion_lines ->
      string -> string list -> string list -> Coupling.fusen -> unit

    val eps4_g4_g44_decl : Format_Fortran.formatter -> unit -> unit
    val eps4_g4_g44_init : Format_Fortran.formatter -> unit -> unit
    val inner_product_functions : Format_Fortran.formatter -> unit -> unit

    module type Test =
      sig
        val suite : OUnit.test
      end

    module Test : Test

  end

module Fortran : T


