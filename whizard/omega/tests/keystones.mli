(* keystones.mli --

   Copyright (C) 2019-2022 by

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

(* A [field] has a Lorentz representation, which will be translated
   to a Fortran type, i.\,e.~the corresponding mnemonic,
   and a position index. *)
type field = Coupling.lorentz * int

(* The different kind of arguments of the fusions. *)
type argument =
  | G of int (* complex coupling *)
  | N of int (* negative of complex coupling *)
  | M of int (* real mass (or width) *)
  | P of int (* momentum *)
  | F of field (* field *)
  | V of string (* verbatim *)

(* A [keystone] is translated to the Fortran expression
   \texttt{bra * name (args)}. *)
type keystone =
  { bra : field;
    name : string;
    args : argument list }

(* A vertex has a unique name [tag] used for the Fortran
   routine and a list of keystones that must \emph{all} produce
   the same result within a reasonable numerical accuracy. *)
type vertex =
  { tag : string;
    keystones : keystone list }

val generate :
  ?reps:int -> ?threshold:float ->
  ?program:string -> ?omega_module:string ->  ?modules:string list ->
  vertex list -> unit

(* In the case of UFO Lorentz structures, we can generate the
   [keystone list] automatically: *)
type ufo_vertex =
  { v_tag : string;
    v_spins : Coupling.lorentz array;
    v_tensor : UFO_Lorentz.t;
    v_flines : Coupling.fermion_lines }

type ufo_propagator =
  { p_tag : string;
    p_omega : string;
    p_spins : Coupling.lorentz * Coupling.lorentz;
    p_propagator : UFO.Propagator.t }

(* Almost always, there is more than one way to write the
   \emph{same} Lorentz structure.  Produce the corresponding
   [ufo_vertex list].  NB: despite the name there is no checking
   for equivalences done here. *)
val equivalent_tensors :
  ?fermion_lines:Coupling.fermion_lines ->
  Coupling.lorentz array -> (string * string) list -> ufo_vertex list

val transpose : ufo_propagator -> ufo_propagator

val generate_ufo :
  ?program:string -> ?omega_module:string -> ?reps:int -> ?threshold:float ->
  ?only_fusions:ufo_vertex list ->
  string -> (ufo_vertex list * vertex) list -> ufo_propagator list -> unit

(* We need different tests for the Majorana spinor permutations,
   transpositions and conjugations. *)
val generate_ufo_bispinors :
  ?program:string -> ?omega_module:string -> ?reps:int -> ?threshold:float ->
  ?only_fusions:ufo_vertex list ->
  string -> (ufo_vertex list * vertex) list -> ufo_propagator list -> unit
