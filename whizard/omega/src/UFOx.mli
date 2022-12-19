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

module Expr :
  sig
    type t
    val of_string : string -> t
    val of_strings : string list -> t
    val substitute : string -> t -> t -> t
    val rename : (string * string) list -> t -> t
    val half : string -> t
    val variables : t -> Sets.String_Caseless.t
    val functions : t -> Sets.String_Caseless.t
  end

module Value :
  sig
    type t
    val of_expr : Expr.t -> t
    val to_string : t -> string
    val to_coupling : (string -> 'b) -> t -> 'b Coupling.expr
  end

    (* \begin{dubious}
         UFO represents rank-2 indices $(i,j)$ as $1000\cdot j + i$.
         This should be replaced by a proper union type eventually.
         Unfortunately, this requires many changes in the [Atom]s in
         [UFOx].  Therefore, we try a quick'n'dirty proof of principle
         first.
       \end{dubious} *)
module type Index =
  sig
    type t = int
               
    val position : t -> int
    val factor : t -> int
    val unpack : t -> int * int
    val pack : int -> int -> t
    val map_position : (int -> int) -> t -> t
    val to_string : t -> string
    val list_to_string : t list -> string


    (* Indices are represented by a pair [int * 'r], where
       ['r] denotes the representation the index belongs to.  *)

    (* [free indices] returns all free indices in the
       list [indices], i.\,e.~all positive indices. *)
    val free : (t * 'r) list -> (t * 'r) list

    (* [summation indices] returns all summation indices in the
       list [indices], i.\,e.~all negative indices.  *)
    val summation : (t * 'r) list -> (t * 'r) list

    val classes_to_string : ('r -> string) -> (t * 'r) list -> string

    (* Generate summation indices, starting from~$-1001$.
       TODO: check that there are no clashes with explicitely
       named indices. *)
    val fresh_summation : unit -> t
    val named_summation : string -> unit -> t

  end

module Index : Index

module type Tensor =
  sig

    type atom

    (* A tensor is a linear combination of products of [atom]s
       with rational coefficients. The following could be refined
       by introducing [scalar] atoms and restricting the denominators
       to [(scalar list * Algebra.QC.t) list].  At the moment, this
       restriction is implemented dynamically by [of_expr] and not
       statically in the type system.
       Polymorphic variants appear to be the right tool, either
       directly or as phantom types.
       However, this is certainly only \textit{nice-to-have}
       and is not essential. *)
    type 'a linear = ('a list * Algebra.QC.t) list
    type t =
      | Linear of atom linear
      | Ratios of (atom linear * atom linear) list

    (* We might need to replace atoms if the syntax is not
       context free. *)
    val map_atoms : (atom -> atom) -> t -> t

    (* We need to rename indices to implement permutations \ldots *)
    val map_indices : (int -> int) -> t -> t

    (* \ldots{} but in order to to clean up inconsistencies
       in the syntax of \texttt{lorentz.py} and
       \texttt{propagators.py} we also need to rename indices
       without touching the second argument of \texttt{P}, the
       argument of \texttt{Mass} etc. *)
    val rename_indices : (int -> int) -> t -> t

    (* We need scale coefficients. *)
    val map_coeff : (Algebra.QC.t -> Algebra.QC.t) -> t -> t

    (* Try to contract adjacent pairs of [atoms] as allowed
       but [Atom.contract_pair].  This is not exhaustive, but
       helps a lot with invariant squares of momenta in
       applications of [Lorentz]. *)
    val contract_pairs : t -> t

    (* The list of variable referenced in the tensor expression,
       that will need to be imported by the numerical code.  *)
    val variables : t -> string list

    (* Parsing and unparsing.  Lists of [string]s are
       interpreted as sums. *)
    val of_expr : UFOx_syntax.expr -> t
    val of_string : string -> t
    val of_strings : string list -> t
    val to_string : t -> string

    (* The supported representations. *)
    type r
    val classify_indices : t -> (int * r) list 
    val rep_to_string : r -> string
    val rep_to_string_whizard : r -> string
    val rep_of_int : bool -> int -> r
    val rep_conjugate : r -> r
    val rep_trivial : r -> bool

    (* There is not a 1-to-1 mapping between the representations
       in the model files and the representations used by O'Mega,
       e.\,g.~in [Coupling.lorentz].  We might need to use heuristics. *)
    type r_omega
    val omega : r -> r_omega

  end

module type Atom =
  sig
    type t
    val map_indices : (int -> int) -> t -> t
    val rename_indices : (int -> int) -> t -> t
    val contract_pair : t -> t -> t option
    val variable : t -> string option
    val scalar : t -> bool
    val is_unit : t -> bool
    val invertible : t -> bool
    val invert : t -> t
    val of_expr : string -> UFOx_syntax.expr list -> t list
    val to_string : t -> string
    type r
    val classify_indices : t list -> (int * r) list
    val disambiguate_indices : t list -> t list
    val rep_to_string : r -> string
    val rep_to_string_whizard : r -> string
    val rep_of_int : bool -> int -> r
    val rep_conjugate : r -> r
    val rep_trivial : r -> bool
    type r_omega
    val omega : r -> r_omega
  end

module type Lorentz_Atom =
  sig

    type dirac = private
      | C of int * int
      | Gamma of int * int * int
      | Gamma5 of int * int
      | Identity of int * int
      | ProjP of int * int
      | ProjM of int * int
      | Sigma of int * int * int * int

    type vector = (* private *)
      | Epsilon of int * int * int * int
      | Metric of int * int
      | P of int * int

    type scalar = (* private *)
      | Mass of int
      | Width of int
      | P2 of int
      | P12 of int * int
      | Variable of string
      | Coeff of Value.t

    type t = (* private *)
      | Dirac of dirac
      | Vector of vector
      | Scalar of scalar
      | Inverse of scalar

    val map_indices_scalar : (int -> int) -> scalar -> scalar
    val map_indices_vector : (int -> int) -> vector -> vector
    val rename_indices_vector : (int -> int) -> vector -> vector

  end

module Lorentz_Atom : Lorentz_Atom

module Lorentz : Tensor
  with type atom = Lorentz_Atom.t and type r_omega = Coupling.lorentz

module type Color_Atom =
  sig
    type t = (* private *)
      | Identity of int * int
      | Identity8 of int * int
      | T of int * int * int
      | F of int * int * int
      | D of int * int * int
      | Epsilon of int * int * int
      | EpsilonBar of int * int * int
      | T6 of int * int * int
      | K6 of int * int * int
      | K6Bar of int * int * int
  end

module Color_Atom : Color_Atom

module Color : Tensor
  with type atom = Color_Atom.t and type r_omega = Color.t

module type Test =
  sig
    val example : unit -> unit
    val suite : OUnit.test
  end
