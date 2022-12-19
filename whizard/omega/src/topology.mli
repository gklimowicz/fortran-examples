(* topology.mli --

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

module type T =
  sig

(* [partition] is a collection of integers, with arity one larger than
   the arity of ['a children] below.  These arities can one fixed number
   corresponding to homogeneous tuples or a collection of tupes or
   lists. *)
    type partition

(* [partitions n] returns the union of
   all~$\lbrack n_1; n_2; \ldots; n_d\rbrack$
   with~$1\le n_1\le n_2\le\ldots\le n_d\le \lfloor n/2\rfloor$ and
   \begin{equation}
     \sum_{i=1}^d n_i = n
   \end{equation}
   for~[d] from~$3$ to~$d_{\max}$, where $d_{\max}$ is a fixed number
   for each module implementating [T].  In particular, if
   [type partition = int * int * int], then [partitions n] returns
   all~$(n_1,n_2,n_3)$ with~$n_1\le n_2\le n_3$ and~$n_1+n_2+n_3=n$. *)
    val partitions : int -> partition list

(* A (poly)tuple as implemented by the modules in [Tuple]: *)
    type 'a children

(* [keystones externals] returns all keystones for the amplitude with
   external states [externals] in the vanilla scalar theory with a
   \begin{equation}
     \sum_{3\le k\le d_{\max}} \lambda_k\phi^k
   \end{equation}
   interaction.  One factor of the products is factorized.  In particular, if
   \begin{quote}
     [type 'a children = 'a Tuple.Binary.t = 'a * 'a],
   \end{quote}
   then [keystones externals] returns all keystones for the amplitude with
   external states [externals] in the vanilla scalar
   $\lambda\phi^3$-theory. *)
    val keystones : 'a list -> ('a list * 'a list children list) list

(* The maximal depth of subtrees for a given number of external lines.  *)
    val max_subtree : int -> int

(* Only for diagnostics: *)
    val inspect_partition : partition -> int list
  end

module Binary : T with type 'a children = 'a Tuple.Binary.t
module Ternary : T with type 'a children = 'a Tuple.Ternary.t
module Mixed23 : T with type 'a children = 'a Tuple.Mixed23.t
module Nary : functor (B : Tuple.Bound) ->
  (T with type 'a children = 'a Tuple.Nary(B).t)

(* \thocwmodulesection{%
     Diagnostics: Counting Diagrams and Factorizations for $\sum_n\lambda_n\phi^n$}
   The number of diagrams for many particles can easily exceed the range of native
   integers.  Even if we can not calculate the corresponding amplitudes, we want
   to check combinatorical factors.  Therefore we code a functor that can use
   arbitray implementations of integers. *)

module type Integer =
  sig
    type t
    val zero : t
    val one : t
    val ( + ) : t -> t -> t
    val ( - ) : t -> t -> t
    val ( * ) : t -> t -> t
    val ( / ) : t -> t -> t
    val pred : t -> t
    val succ : t -> t
    val ( = ) : t -> t -> bool
    val ( <> ) : t -> t -> bool
    val ( < ) : t -> t -> bool
    val ( <= ) : t -> t -> bool
    val ( > ) : t -> t -> bool
    val ( >= ) : t -> t -> bool
    val of_int : int -> t
    val to_int : t -> int
    val to_string : t -> string
    val compare : t -> t -> int
    val factorial : t -> t
  end

(* Of course, native integers will provide the fastest implementation: *)
module Int : Integer

module type Count =
  sig
    type integer

(* [diagrams f d n] returns the number of tree diagrams contributing
   to the $n$-point amplitude in vanilla scalar theory with
   \begin{equation}
     \sum_{3\le k\le d \land f(k)} \lambda_k\phi^k
   \end{equation}
   interaction.  The default value of~[f] returns [true] for all
   arguments.  *)
    val diagrams : ?f:(integer -> bool) -> integer -> integer -> integer
    val diagrams_via_keystones : integer -> integer -> integer

(* \begin{equation}
     \frac{1}{S(n_k,n-n_k)} \frac{1}{S(n_1,n_2,\ldots,n_k)}
       \binom{n_1+n_2+\ldots+n_k}{n_1,n_2,\ldots,n_k}
   \end{equation} *)
    val keystones : integer list -> integer

(* [diagrams_via_keystones d n] must produce the same
   results as [diagrams d n].  This is shown explicitely in
   tables~\ref{tab:keystone-check}, \ref{tab:keystone-check4} and
   \ref{tab:keystone-check6} for small values of~[d] and~[n].
   The test program in appendix~\ref{sec:count} can be used to
   verify this relation for larger values. *)
    val diagrams_per_keystone : integer -> integer list -> integer

  end

module Count : functor (I : Integer) -> Count with type integer = I.t

(* \thocwmodulesection{Emulating HELAC} *)

(* We can also proceed \'a la~\cite{HELAC:2000}.  *)
module Helac : functor (B : Tuple.Bound) ->
  (T with type 'a children = 'a Tuple.Nary(B).t)

(* \begin{dubious}
     The following has never been tested, but it is no rocket science and
     should work anyway \ldots
   \end{dubious} *)
module Helac_Binary : T with type 'a children = 'a Tuple.Binary.t

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)

