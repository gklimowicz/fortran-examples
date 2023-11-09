(* vertex.mli --

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

val parse_string : string -> UFO_syntax.t
val parse_file : string -> UFO_syntax.t

(* These are the contents of the Python files after lexical
   analysis as context-free variable declarations, before
   any semantic interpretation. *)

module type Files =
  sig
    
    type t = private
      { particles : UFO_syntax.t;
	couplings : UFO_syntax.t;
	coupling_orders : UFO_syntax.t;
	vertices : UFO_syntax.t;
	lorentz : UFO_syntax.t;
	parameters : UFO_syntax.t;
	propagators : UFO_syntax.t;
	decays : UFO_syntax.t }

    val parse_directory : string -> t

  end

type t

exception Unhandled of string

module Model : Model.T

val parse_directory : string -> t

module type Fortran_Target =
  sig

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

    val lorentz_module :
      ?only:Sets.String.t -> ?name:string ->
      ?fortran_module:string -> ?parameter_module:string ->
      Format_Fortran.formatter -> unit -> unit

  end

module Targets :
  sig
    module Fortran : Fortran_Target
  end

(* Export some functions for testing: *)

module Propagator_UFO :
  sig
    type t = (* private *)
      { name : string;
	numerator : UFOx.Lorentz.t;
	denominator : UFOx.Lorentz.t }
  end

module Propagator :
  sig
    type t = (* private *)
      { name : string;
        spins : Coupling.lorentz * Coupling.lorentz;
	numerator : UFO_Lorentz.t;
	denominator : UFO_Lorentz.t;
        variables : string list }
    val of_propagator_UFO : ?majorana:bool -> Propagator_UFO.t -> t
    val transpose : t -> t
  end

module type Test =
  sig
    val suite : OUnit.test
  end

module Test : Test
