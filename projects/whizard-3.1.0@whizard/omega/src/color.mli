(* color.mli --

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
    val suite_long : OUnit.test
  end

(* \thocwmodulesection{Quantum Numbers} *)

(* Color is not necessarily the~$\textrm{SU}(3)$ of QCD.  Conceptually,
   it can be any \emph{unbroken} symmetry (\emph{broken} symmetries correspond
   to [Model.flavor]).  In order to keep the group theory simple, we confine
   ourselves to the fundamental and adjoint representation
   of a single~$\textrm{SU}(N_C)$ for the moment.  Therefore,
   particles are either color singlets or live in the defining
   representation of $\textrm{SU}(N_C)$: [SUN]$(|N_C|)$, its conjugate
   [SUN]$(-|N_C|)$ or in the adjoint representation of
   $\textrm{SU}(N_C)$: [AdjSUN]$(N_C)$. *)

type t = Singlet | SUN of int | AdjSUN of int

val conjugate : t -> t
val compare : t -> t -> int

(* \thocwmodulesection{Color Flows} *)

(* This computes the color flow as used by WHIZARD: *)

module type Flow =
  sig

    type color
    type t = color list * color list
    val rank : t -> int

    val of_list : int list -> color
    val ghost : unit -> color
    val to_lists : t -> int list list
    val in_to_lists : t -> int list list
    val out_to_lists : t -> int list list
    val ghost_flags : t -> bool list
    val in_ghost_flags : t -> bool list
    val out_ghost_flags : t -> bool list

(* A factor is a list of powers
   \begin{equation}
     \sum_{i}
        \left( \frac{\ocwlowerid{num}_i}{\ocwlowerid{den}_i}
                  \right)^{\ocwlowerid{power}_i}
   \end{equation} *)
    type power = { num : int; den : int; power : int }
    type factor = power list

    val factor : t -> t -> factor
    val zero : factor

    module Test : Test

  end

module Flow : Flow

(* \thocwmodulesection{Vertex Color Flows} *)

(* \begin{dubious}
     The following is (still work-in-progress) infrastructure for
     translating UFO style color factors into color flows.
   \end{dubious} *)

(* \begin{dubious}
     It might be beneficial, to use the color flow representation
     here.  This will simplify the colorizer at the price of
     some complexity in [UFO] or here.
   \end{dubious} *)

(* The datatypes [Arrow.free] and [Arrow.factor] will be used as
   building blocks for [Birdtracks.t] below. *)
module type Arrow =
  sig

    (* For fundamental and adjoint representations, the endpoints
       of arrows are uniquely specified by a vertex (which will
       be represented by a number).  For representations with more
       than one outgoing or incoming arrow, we need an additional index.
       This is abrcated in the [endpoint] type. *)
    type endpoint

    (* Endpoints can be the the tip or tail of an arrow or a ghost.
       Currently, we use the types for illustration only, but we
       might eventually try to make them abstract for additional
       safety.. *)
    type tip = endpoint
    type tail = endpoint
    type ghost = endpoint

    (* The position of the endpoint is encoded as an integer, which
       can be mapped, if necessary. *)
    val position : endpoint -> int
    val relocate : (int -> int) -> endpoint -> endpoint

    (* An [Arrow.t] is either a genuine arrow or a ghost \ldots *)
    type ('tail, 'tip, 'ghost) t =
      | Arrow of 'tail * 'tip
      | Ghost of 'ghost
      | Epsilon of 'tip list
      | Epsilon_bar of 'tail list

    (* {}\ldots and we distuish [free] arrows that must not contain
       summation indices from [factor]s that may.  Indices are
       opaque.  [('tail, 'tip, 'ghost) t] is polymorphic so that
       we can use richer ['tail], ['tip] and ['ghost] in [factor]. *)
    type free = (tail, tip, ghost) t
    type factor

    (* For debugging, logging, etc. *)
    val free_to_string : free -> string
    val factor_to_string : factor -> string

    (* Change the [endpoint]s in a [free] arrow. *)
    val map : (endpoint -> endpoint) -> free -> free

    (* Turn the [endpoint]s satisfying the predicate into a
       left or right hand side summation index.  Left and right
       refer to the two factors in a product and
       we must only match arrows with [endpoint]s in both
       factors, not double lines on either side.
       Typically, the predicate will be set up to select only the
       summation indices that appear on both sides.*)
    
    val to_left_factor : (endpoint -> bool) -> free -> factor
    val to_right_factor : (endpoint -> bool) -> free -> factor

    (* The incomplete inverse [of_factor] raises an exception
       if there are remaining summation indices.  [is_free] can
       be used to check first. *)
    val of_factor : factor -> free
    val is_free : factor -> bool

    (* Return all the endpoints of the arrow that have a [position]
       encoded as a negative integer.  These are treated as summation
       indices in our applications. *)
    val negatives : free -> endpoint list

    (* We will need to test whether an arrow represents a ghost. *)
    val is_ghost : free -> bool

    (* An arrow looping back to itself. *)
    val is_tadpole : factor -> bool

    (* An $\epsilon$ or an $\bar\epsilon$ *)
    val is_epsilon : factor -> bool

(* If [arrow] is an~$\epsilon$ (or $\bar\epsilon$) and [arrows] contains
   an~$\bar\epsilon$ (or $\epsilon$), use
   \begin{equation}
      \forall n, N \in\mathbf{N}, 2\le n \le N:\;
      \epsilon_{i_1i_2\cdots i_n} \bar\epsilon^{j_1j_2\cdots j_n}
        = \sum_{\sigma\in S_n} (-1)^{\varepsilon(\sigma)}
            \delta_{i_1}^{\sigma(j_1)} 
            \delta_{i_2}^{\sigma(j_2)} 
            \cdots
            \delta_{i_n}^{\sigma(j_n)}\,,
   \end{equation}
   where~$N=\delta_i^i$ is the dimension, to expand the pair into two lists of
   list of arrows: the first corresponding to the even permutations, the
   second to the odd ones.  In addition, return the remaining arrows. *)
    val match_epsilon : factor -> factor list -> (factor list list * factor list list * factor list) option

    (* Merging two arrows can give a variety of results.
       NB: $\epsilon$-$\bar\epsilon$ pairs are assumed to have been
       already expanded by [match_epsilon]. *)
    type merge =
      | Match of factor  (* a tip fits the other's tail: make one arrow out of two *)
      | Ghost_Match (* two matching ghosts *)
      | Loop_Match (* both tips fit both tails: drop the arrows *)
      | Mismatch (* ghost meets arrow: error *)
      | No_Match (* nothing to be done *)
    val merge : factor -> factor -> merge

(* Break up an arrow [tee a (i => j) -> [i => a; a => j]], i.\,e.~insert
   a gluon. Returns an empty list for a ghost and raises an exception
   for~$\epsilon$ and~$\bar\epsilon$. *)
    val tee : int -> free -> free list

(* [dir i j arrow] returns the direction of the arrow relative to [j => i].
   Returns 0 for a ghost and raises an exception for~$\epsilon$
   and~$\bar\epsilon$. *)
    val dir : int -> int -> free -> int

(* It's intuitive to use infix operators to construct the lines. *)
    val single : endpoint -> endpoint -> free
    val double : endpoint -> endpoint -> free list
    val ghost : endpoint -> free

    module Infix : sig

      (* [single i j] or [i => j] creates a single line from [i] to [j] and
         [i ==> j] is a shorthard for [[i => j]]. *)
      val (=>) : int -> int -> free
      val (==>) : int -> int -> free list

      (* [double i j] or [i <=> j] creates a double line from [i] to [j]
         and back. *)
      val (<=>) : int -> int -> free list

      (* Single lines with subindices at the tip and/or tail *)
      val (>=>) : int * int -> int -> free
      val (=>>) : int -> int * int -> free
      val (>=>>) : int * int -> int * int -> free

      (* [?? i] creates a ghost at [i]. *)
      val (??) : int -> free

      (* NB: I wanted to use [~~] instead of [??], but ocamlweb can't handle
         operators starting with [~] in the index properly. *)

    end

    val epsilon : int list -> free
    val epsilon_bar : int list -> free

    (* [chain [1;2;3]] is a shorthand for [[1 => 2; 2 => 3]] and
       [cycle [1;2;3]] for [[1 => 2; 2 => 3; 3 => 1]].  Other lists
       and edge cases are handled in the natural way. *)
    val chain : int list -> free list
    val cycle : int list -> free list

    module Test : Test

    (* Pretty printer for the toplevel. *)
    val pp_free : Format.formatter -> free -> unit
    val pp_factor : Format.formatter -> factor -> unit

  end

module Arrow : Arrow

(* Possible color flows for a single propagator, as currently
   supported by WHIZARD. *)
module type Propagator =
  sig
    type cf_in = int
    type cf_out = int
    type t = W | I of cf_in | O of cf_out | IO of cf_in * cf_out | G
    val to_string : t -> string
  end

module Propagator : Propagator

(* Implement birdtracks operations as generally as possible.
   Below, the signature will be extended with group specific
   generators for $\mathrm{SU}(N_C)$ and $\mathrm{U}(N_C)$ and
   even $N_C=3$. *)
module type Birdtracks =
  sig
    type t

    (* Strip out redundancies. *)
    val canonicalize : t -> t

    (* Debugging, logging, etc. *)
    val to_string : t -> string

    (* Test for trivial color flows that are just a number. *)
    val trivial : t -> bool

    (* Test for vanishing coefficients. *)
    val is_null : t -> bool

    (* Purely numeric factors, implemented as Laurent polynomials
       (cf.~[Algebra.Laurent] in~$N_C$ with complex rational
       coefficients. *)
    val const : Algebra.Laurent.t -> t
    val null : t (* $0$ *)
    val one : t (* $1$ *)
    val two : t (* $2$ *)
    val half : t (* $1/2$ *)
    val third : t (* $1/3$ *)
    val minus : t (* $-1$ *)
    val int : int -> t (* $n$ *)
    val fraction : int -> t (* $1/n$ *)
    val nc : t (* $N_C$ *)
    val over_nc : t (* $1/N_C$ *)
    val imag : t (* $\ii$ *)

    (* Shorthand: $\{(c_i,p_i)\}_i\to \sum_i c_i (N_C)^{p_i}$*)
    val ints : (int * int) list -> t

    val scale : Algebra.QC.t -> t -> t

    val sum : t list -> t
    val diff : t -> t -> t
    val times : t -> t -> t
    val multiply : t list -> t

    (* For convenience, here are infix versions of the above operations. *)
    module Infix : sig
      val ( +++ ) : t -> t -> t
      val ( --- ) : t -> t -> t
      val ( *** ) : t -> t -> t
    end

   (* We can compute the $f_{abc}$ and $d_{abc}$ invariant tensors
      from the generators of an arbitrary representation:
      \begin{subequations}
      \begin{align}
       f_{a_1a_2a_3} &=
        - \ii \tr\left(T_{a_1}\left\lbrack T_{a_2},T_{a_3}\right\rbrack_-\right)
          = - \ii \tr\left(T_{a_1}T_{a_2}T_{a_3}\right)
            + \ii \tr\left(T_{a_1}T_{a_3}T_{a_2}\right) \\
       d_{a_1a_2a_3} &=
         \tr\left(T_{a_1}\left\lbrack T_{a_2},T_{a_3}\right\rbrack_+\right)
          =   \tr\left(T_{a_1}T_{a_2}T_{a_3}\right)
            + \tr\left(T_{a_1}T_{a_3}T_{a_2}\right)\,
      \end{align}
      \end{subequations}
      assuming the normalization $ \tr(T_aT_b) = \delta_{ab}$.

      NB: this uses the summation indices $-1$, $-2$ and $-3$.  Therefore
      it \emph{must not} appear unevaluated more than once in a product! *)
    val f_of_rep : (int -> int -> int -> t) -> int -> int -> int -> t
    val d_of_rep : (int -> int -> int -> t) -> int -> int -> int -> t

    (* Rename the indices of endpoints in a birdtrack. *)
    val relocate : (int -> int) -> t -> t

    (* [fuse nc vertex children] use the color flows in the [vertex]
       to combine the color flows in the incoming [children] and return
       the color flows for outgoing particle together with their weights. *)
    val fuse : int -> t -> Propagator.t list -> (Algebra.QC.t * Propagator.t) list

    module Test : Test

    (* Pretty printer for the toplevel. *)
    val pp : Format.formatter -> t -> unit
  end

module Birdtracks : Birdtracks

module type SU3 =
  sig
    include Birdtracks
    val delta3 : int -> int -> t
    val delta8 : int -> int -> t
    val delta8_loop : int -> int -> t
    val gluon : int -> int -> t
    val delta6 : int -> int -> t
    val delta10 : int -> int -> t
    val t : int -> int -> int -> t
    val f : int -> int -> int -> t
    val d : int -> int -> int -> t
    val epsilon : int list -> t
    val epsilon_bar : int list -> t
    val t8 : int -> int -> int -> t
    val t6 : int -> int -> int -> t
    val t10 : int -> int -> int -> t
    val k6 : int -> int -> int -> t
    val k6bar : int -> int -> int -> t
    val delta_of_tableau : int Young.tableau -> int -> int -> t
    val t_of_tableau : int Young.tableau -> int -> int -> int -> t
  end

module SU3 : SU3
module Vertex : SU3

(* \begin{dubious}
     This must not be used, because it has not yet been updated
     to the correctly symmetrized version!
   \end{dubious} *)
module U3 : SU3

