(* algebra.mli --

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

module type Test =
  sig
    val suite : OUnit.test
  end

(* \thocwmodulesection{Coefficients} *)

(* For our algebra, we need coefficient rings. *)

module type CRing =
  sig
    type t
    val null : t
    val unit : t
    val mul : t -> t -> t
    val add : t -> t -> t
    val sub : t -> t -> t
    val neg : t -> t
    val to_string : t -> string
  end

(* And rational numbers provide a particularly important example: *)

module type Rational =
  sig
    include CRing
    val is_null : t -> bool
    val is_unit : t -> bool
    val is_positive : t -> bool
    val is_negative : t -> bool
    val is_integer : t -> bool
    val make : int -> int -> t
    val abs : t -> t
    val inv : t -> t
    val div : t -> t -> t
    val pow : t -> int -> t
    val sum : t list -> t
    val to_ratio : t -> int * int
    val to_float : t -> float
    val to_integer : t -> int
    module Test : Test
  end

(* \thocwmodulesection{Naive Rational Arithmetic} *)

(* \begin{dubious}
     This \emph{is} dangerous and will overflow even for simple
     applications.  The production code will have to be linked to
     a library for large integer arithmetic.
   \end{dubious} *)

module Small_Rational : Rational
module Q : Rational

(* \thocwmodulesection{Rational Complex Numbers} *)

module type QComplex =
  sig

    type q
    type t

    val make : q -> q -> t
    val null : t
    val unit : t

    val real : t -> q
    val imag : t -> q

    val conj : t -> t
    val neg : t -> t

    val add : t -> t -> t
    val sub : t -> t -> t
    val mul : t -> t -> t
    val inv : t -> t
    val div : t -> t -> t

    val pow : t -> int -> t
    val sum : t list -> t

    val is_null : t -> bool
    val is_unit : t -> bool
    val is_positive : t -> bool
    val is_negative : t -> bool
    val is_integer : t -> bool
    val is_real : t -> bool

    val to_string : t -> string

    module Test : Test

  end

module QComplex : functor (Q' : Rational) -> QComplex with type q = Q'.t
module QC : QComplex with type q = Q.t

(* \thocwmodulesection{Laurent Polynomials} *)

(* Polynomials, including negative powers, in one variable.
   In our applications, the variable~$x$ will often be~$N_C$,
   the number of colors
   \begin{equation}
     \sum_n c_n N_C^n
   \end{equation} *)
module type Laurent =
  sig

    (* The type of coefficients.  In the implementation below,
       it is [QComplex.t]: complex numbers with rational real
       and imaginary parts. *)
    type c
    type t

    (* Elementary constructors *)
    val null : t
    val is_null : t -> bool
    val unit : t

    (* [atom c n] constructs a term $c x^n$, where $x$ denotes
       the variable. *)
    val atom : c -> int -> t

    (* Shortcut: [const c = atom c 0] *)
    val const : c -> t

    (* Elementary arithmetic *)
    val scale : c -> t -> t
    val neg : t -> t
    val add : t -> t -> t
    val diff : t -> t -> t
    val sum : t list -> t
    val mul : t -> t -> t
    val product : t list -> t
    val pow : int -> t -> t

    (* [eval c p] evaluates the polynomial [p] by substituting
       the constant [c] for the variable. *)
    val eval : c -> t -> c

    (* A total ordering.  Does not correspond to any mathematical order. *)
    val compare : t -> t -> int

    (* Logging, debugging and toplevel integration. *)
    val to_string : string -> t -> string
    val pp : Format.formatter -> t -> unit
    module Test : Test
  end

(* \begin{dubious}
     Could (should?) be functorialized over [QComplex].
     We had to wait until we upgraded our O'Caml requirements to 4.02,
     but that has been done.
   \end{dubious} *)

module Laurent : Laurent with type c = QC.t

(* \thocwmodulesection{Expressions: Terms, Rings and Linear Combinations} *)

(* The tensor algebra will be spanned by an abelian monoid: *)

module type Term =
  sig
    type 'a t
    val unit : unit -> 'a t
    val is_unit : 'a t -> bool
    val atom : 'a -> 'a t
    val power : int -> 'a t -> 'a t
    val mul : 'a t -> 'a t -> 'a t
    val map : ('a -> 'b) -> 'a t -> 'b t
    val to_string : ('a -> string) -> 'a t -> string

    (* The derivative of a term is \emph{not} a term,
       but a sum of terms instead:
       \begin{equation}
           D (f_1^{p_1}f_2^{p_2}\cdots f_n^{p_n}) =
             \sum_i (Df_i) p_i f_1^{p_1}f_2^{p_2}\cdots f_i^{p_i-1} \cdots f_n^{p_n}
       \end{equation}
       The function returns the sum as a list of triples
       $(Df_i,p_i, f_1^{p_1}f_2^{p_2}\cdots f_i^{p_i-1} \cdots f_n^{p_n})$.
       Summing the terms is left to the calling module and the $Df_i$ are
       \emph{not} guaranteed to be different.
       NB: The function implementating the inner derivative, is supposed to
       return~[Some]~$Df_i$ and [None], iff~$Df_i$ vanishes. *)
    val derive : ('a -> 'b option) -> 'a t -> ('b * int * 'a t) list

    (* convenience function *)
    val product : 'a t list -> 'a t
    val atoms : 'a t -> 'a list

  end

module type Ring =
  sig
    module C : Rational
    type 'a t
    val null : unit -> 'a t
    val unit : unit -> 'a t
    val is_null : 'a t -> bool
    val is_unit : 'a t -> bool
    val atom : 'a -> 'a t
    val scale : C.t -> 'a t -> 'a t
    val add : 'a t -> 'a t -> 'a t
    val sub : 'a t -> 'a t -> 'a t
    val mul : 'a t -> 'a t -> 'a t
    val neg : 'a t -> 'a t

    (* Again
       \begin{equation}
           D (f_1^{p_1}f_2^{p_2}\cdots f_n^{p_n}) =
             \sum_i (Df_i) p_i f_1^{p_1}f_2^{p_2}\cdots f_i^{p_i-1} \cdots f_n^{p_n}
       \end{equation}
       but, iff~$Df_i$ can be identified with a~$f'$, we know how to perform
       the sum. *)

    val derive_inner : ('a -> 'a t) -> 'a t -> 'a t (* this? *)
    val derive_inner' : ('a -> 'a t option) -> 'a t -> 'a t (* or that? *)

(* Below, we will need partial derivatives that lead out of the ring:
   [derive_outer derive_atom term] returns a list of partial derivatives
   ['b] with non-zero coefficients ['a t]: *)
    val derive_outer : ('a -> 'b option) -> 'a t -> ('b * 'a t) list

    (* convenience functions *)
    val sum : 'a t list -> 'a t
    val product : 'a t list -> 'a t

(* The list of all generators appearing in an expression: *)
    val atoms : 'a t -> 'a list

    val to_string : ('a -> string) -> 'a t -> string

  end

module type Linear =
  sig
    module C : Ring
    type ('a, 'c) t
    val null : unit -> ('a, 'c) t
    val atom : 'a -> ('a, 'c) t
    val singleton : 'c C.t -> 'a -> ('a, 'c) t
    val scale : 'c C.t -> ('a, 'c) t -> ('a, 'c) t
    val add : ('a, 'c) t -> ('a, 'c) t -> ('a, 'c) t
    val sub : ('a, 'c) t -> ('a, 'c) t -> ('a, 'c) t

(* A partial derivative w.\,r.\,t.~a vector maps from a coefficient ring to
   the dual vector space.  *)
    val partial : ('c -> ('a, 'c) t) -> 'c C.t -> ('a, 'c) t

(* A linear combination of vectors
   \begin{equation}
     \text{[linear]} \lbrack (v_1, c_1); (v_2, c_2); \ldots; (v_n, c_n)\rbrack
        = \sum_{i=1}^{n} c_i\cdot v_i
   \end{equation} *)
    val linear : (('a, 'c) t * 'c C.t) list -> ('a, 'c) t

(* Some convenience functions *)
    val map : ('a -> 'c C.t -> ('b, 'd) t) -> ('a, 'c) t ->  ('b, 'd) t
    val sum : ('a, 'c) t list -> ('a, 'c) t

(* The list of all generators and the list of all generators of coefficients
   appearing in an expression: *)
    val atoms : ('a, 'c) t -> 'a list * 'c list

    val to_string : ('a -> string) -> ('c -> string) -> ('a, 'c) t -> string

  end

module Term : Term

module Make_Ring (C : Rational) (T : Term) : Ring
module Make_Linear (C : Ring) : Linear with module C = C
