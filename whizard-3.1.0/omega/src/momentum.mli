(* momentum.mli --

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

(* Model the finite combinations
   \begin{equation}
     p = \sum_{n=1}^k c_k \bar p_n,\qquad \text{(with $c_k\in\{0,1\}$)}
   \end{equation}
   of~$n_{\text{in}}$ incoming and~$k-n_{\text{in}}$ outgoing momenta~$p_n$
   \begin{equation}
     \bar p_n =
       \begin{cases}
          - p_n & \text{for $1\le n \le n_{\text{in}}$} \\
            p_n & \text{for $n_{\text{in}}+1\le n\le k$}
       \end{cases}
   \end{equation}
   where momentum is conserved
   \begin{equation}
     \sum_{n=1}^k \bar p_n = 0
   \end{equation}
   below, we need the notion of `rank' and `dimension':
   \begin{subequations}
   \begin{align}
      \text{\ocwlowerid{dim}} (p) &= k \\
      \text{\ocwlowerid{rank}} (p) &= \sum_{n=1}^{k} c_k
   \end{align}
   \end{subequations}
   where `dimension' is \emph{not} the dimension of the
   underlying space-time, of course. *)

module type T =
  sig
    type t

(* Constructor: $(k,N)\to p = \sum_{n\in N} \bar p_n$ and
   $k=\text{\ocwlowerid{dim}}(p)$ is the \emph{overall} number
   of independent momenta, while $\text{\ocwlowerid{rank}}(p)=|N|$
   is the number of momenta in~$p$. It would be possible to
   fix~[dim] as a functor argument instead.  This might
   be slightly faster and allow a few more compile time checks,
   but would be much more tedious to use, since the number
   of particles will be chosen at runtime. *) 
    val of_ints : int -> int list -> t

(* No two indices may be the same.  Implementions of [of_ints] can
   either raise the exception [Duplicate] or ignore the duplicate,
   but implementations of [add] are required to raise [Duplicate]. *)
    exception Duplicate of int

(* Raise [Range] iff $n>k$: *)
    exception Range of int

(* Binary oparations require that both momenta have the same dimension.
   [Mismatch] is raised if this condition is violated.  *)
    exception Mismatch of string * t * t

(* [Negative] is raised if the result of [sub] is undefined.  *)
    exception Negative

(* The inverses of the constructor (we have
   [rank p = List.length (to_ints p)], but [rank] might be more efficient): *)
    val to_ints : t -> int list
    val dim : t -> int
    val rank : t -> int

(* Shortcuts: [singleton d p = of_ints d [p]] and [zero d = of_ints d []]: *)
    val singleton : int -> int -> t
    val zero : int -> t

(* An arbitrary total order, with the condition
   $\text{\ocwlowerid{rank}}(p_1)<\text{\ocwlowerid{rank}}(p_2)
    \Rightarrow p_1<p_2$.  *)
    val compare : t -> t -> int

(* Use momentum conservation to construct the negative momentum with
   positive coefficients: *)
    val neg : t -> t

(* Return the momentum or its negative, whichever has the lower rank.
   NB: the present implementation does \emph{not} guarantee that
   \begin{equation}
     \text{abs} p = \text{abs} q \Longleftrightarrow p = p \lor p = - q
   \end{equation}
   for momenta with $\text{rank} = \text{dim}/2$. *)
    val abs : t -> t

(* Add and subtract momenta.  This can fail, since the coefficients~$c_k$ must
   me either~$0$ or~$1$. *)
    val add : t -> t -> t
    val sub : t -> t -> t

(* Once more, but not raising exceptions this time: *)
    val try_add : t -> t -> t option
    val try_sub : t -> t -> t option

(* \emph{Not} the total order provided by [compare], but set inclusion of
   non-zero coefficients instead: *)
    val less : t -> t -> bool
    val lesseq : t -> t -> bool

(* $p_1 + (\pm p_2) + (\pm p_3) = 0$ *)
    val try_fusion : t -> t -> t -> (bool * bool) option

(* A textual representation for debugging: *)
    val to_string : t -> string

(* [split i n p] splits~$\bar p_i$ into~$n$ momenta~$\bar p_i \to
   \bar p_i + \bar p_{i+1} + \ldots + \bar p_{i+n-1}$ and makes room
   via~$\bar p_{j>i} \to \bar p_{j+n-1}$.  This is used for implementating
   cascade decays, like combining
   \begin{subequations}
   \begin{align}
     \mathrm{e}^+(p_1) \mathrm{e}^-(p_2) \to
       &\mathrm{W}^-(p_3) \nu_{\mathrm{e}}(p_4) \mathrm{e}^+(p_5)\\
       &\mathrm{W}^-(p_3)\to \mathrm{d}(p_3') \bar{\mathrm{u}}(p_4')
   \end{align}
   \end{subequations}
   to
   \begin{equation}
     \mathrm{e}^+(p_1) \mathrm{e}^-(p_2) \to
       \mathrm{d}(p_3) \bar{\mathrm{u}}(p_4)
        \nu_{\mathrm{e}}(p_5) \mathrm{e}^+(p_6)
   \end{equation}
   in narrow width approximation for the~$\mathrm{W}^-$. *)
    val split : int -> int -> t -> t

(* \thocwmodulesection{Scattering Kinematics}
   From here on, we assume scattering kinematics $\{1,2\}\to\{3,4,\ldots\}$,
   i.\,e.~$n_{\text{in}}=2$.
   \begin{dubious}
     Since functions like [timelike] can be used for decays as well (in which
     case they must \emph{always} return [true], the representation---and
     consequently the constructors---should be extended by a flag discriminating
     between the two cases!
   \end{dubious} *)

    module Scattering :
        sig

(* Test if the momentum is an incoming one: $p=\bar p_1\lor p=\bar p_2$ *)
          val incoming : t -> bool

(* $p=\bar p_3\lor p=\bar p_4\lor \ldots$ *)
          val outgoing : t -> bool

(* $p^2 \ge 0$.  NB: \textit{par abus de langange}, we report the incoming
   individual momenta as spacelike, instead as timelike.  This will be useful
   for phasespace constructions below. *)
          val timelike : t -> bool

(* $p^2 \le 0$.  NB: the simple algebraic criterion can be violated for heavy
   initial state particles. *)
          val spacelike : t -> bool

(* $p = \bar p_1 + \bar p_2$ *)
          val s_channel_in : t -> bool

(* $p = \bar p_3 + \bar p_4 + \ldots + \bar p_n$ *)
          val s_channel_out : t -> bool

(* $p = \bar p_1 + \bar p_2 \lor p = \bar p_3 + \bar p_4 + \ldots + \bar p_n$ *)
          val s_channel : t -> bool

(* $ \bar p_1 + \bar p_2 \to \bar p_3 + \bar p_4 + \ldots + \bar p_n$ *)
          val flip_s_channel_in : t -> t
      end

(* \thocwmodulesection{Decay Kinematics} *)
    module Decay :
        sig

(* Test if the momentum is an incoming one: $p=\bar p_1$ *)
          val incoming : t -> bool

(* $p=\bar p_2\lor p=\bar p_3\lor \ldots$ *)
          val outgoing : t -> bool

(* $p^2 \ge 0$.  NB: here, we report the incoming
   individual momenta as timelike. *)
          val timelike : t -> bool

(* $p^2 \le 0$. *)
          val spacelike : t -> bool

        end

  end

module Lists : T
module Bits : T
module Default : T

(* Wolfgang's funny tree codes:
   \begin{equation}
   (2^n, 2^{n-1}) \to (1, 2, 4, \ldots, 2^{n-2})
   \end{equation} *)

module type Whizard =
  sig
    type t
    val of_momentum : t -> int
    val to_momentum : int -> int -> t
  end

module ListsW : Whizard with type t = Lists.t
module BitsW : Whizard with type t = Bits.t
module DefaultW : Whizard with type t = Default.t

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
