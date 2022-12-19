(* UFO_Lorentz.mli --

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

(* \thocwmodulesection{Processed UFO Lorentz Structures} *)

(* Just like [UFOx.Lorentz_Atom.dirac], but without the Dirac matrix indices. *)
type dirac = (* [private] *)
  | Gamma5
  | ProjM
  | ProjP
  | Gamma of int
  | Sigma of int * int
  | C
  | Minus

(* A sandwich of a string of $\gamma$-matrices. [bra] and [ket] are
   positions of fields in the vertex, \emph{not} spinor indices. *)
type dirac_string = (* [private] *)
  { bra : int;
    ket : int;
    conjugated : bool;
    gammas : dirac list }

(* In the case of Majorana spinors, we have to insert charge conjugation
   matrices. *)

(* $\Gamma\to - \Gamma$: *)
val minus : dirac_string -> dirac_string

(* $\Gamma\to C\Gamma$: *)
val cc_times : dirac_string -> dirac_string

(* $\Gamma\to - \Gamma C$: *)
val times_minus_cc : dirac_string -> dirac_string

(* $\Gamma\to \Gamma^T$: *)
val transpose : dirac_string -> dirac_string

(* $\Gamma\to C\Gamma C^{-1}$: *)
val conjugate : dirac_string -> dirac_string

(* $\Gamma\to C\Gamma^T C^{-1}$, i.\,e.~the composition of [conjugate]
   and [transpose]: *)
val conjugate_transpose : dirac_string -> dirac_string

(* The Lorentz indices appearing in a term are either negative
   internal summation indices or positive external polarization
   indices.  Note that the external
   indices are not really indices, but denote the position
   of the particle in the vertex. *)
type 'a term =  (* [private] *)
  { indices : int list;
    atom : 'a }

(* Split the list of indices into summation and polarization indices. *)
val classify_indices : int list -> int list * int list

(* Replace the atom keeping the associated indices. *)
val map_atom : ('a -> 'b) -> 'a term -> 'b term

(* A contraction consists of a (possibly empty) product of
   Dirac strings and a (possibly empty) product of Lorentz
   tensors with a rational coefficient.
   The [denominator] is required for the poorly documented
   propagator extensions.  The type [atom linear] is
   a [list] and an empty list is interpreted as~$1$. *)
(* \begin{dubious}
     The [denominator] is a [contraction list] to allow code reuse,
     though a [(A.scalar list * A.scalar list * QC.t) list] would
     suffice.
   \end{dubious} *)
type contraction = (* [private] *)
  { coeff : Algebra.QC.t;
    dirac : dirac_string term list;
    vector : UFOx.Lorentz_Atom.vector term list;
    scalar : UFOx.Lorentz_Atom.scalar list;
    inverse : UFOx.Lorentz_Atom.scalar list;
    denominator : contraction list }

(* A sum of [contraction]s. *)
type t = contraction list

(* Fermion line connections. *)
val fermion_lines : t -> Coupling.fermion_lines

(* $\Gamma\to C\Gamma C^{-1}$ *)
val charge_conjugate : int * int -> t -> t

(* [parse spins lorentz] uses the [spins] to parse the
   UFO [lorentz] structure as a list of [contraction]s. *)
val parse : ?allow_denominator:bool -> Coupling.lorentz list -> UFOx.Lorentz.t -> t

(* [map_indices f lorentz] applies the map [f] to the free
   indices in [lorentz]. *)
val map_indices : (int -> int) -> t -> t
val map_fermion_lines :
  (int -> int) -> Coupling.fermion_lines -> Coupling.fermion_lines

(* Create a readable representation for debugging and
   documenting generated code. *)
val to_string : t -> string
val fermion_lines_to_string : Coupling.fermion_lines -> string

(* Punting \ldots *)
val dummy : t

(* More debugging and documenting. *)
val dirac_string_to_string : dirac_string -> string

(* [dirac_string_to_matrix substitute ds] take a string
   of $\gamma$-matrices [ds], applies [substitute] to
   the indices and returns the product as a matrix. *)
val dirac_string_to_matrix : (int -> int) -> dirac_string -> Dirac.Chiral.t

module type Test =
  sig
    val suite : OUnit.test
  end

module Test : Test
