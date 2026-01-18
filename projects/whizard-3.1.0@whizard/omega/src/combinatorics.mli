(* combinatorics.mli --

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

(* This type is defined just for documentation.  Below, most functions will 
   construct a (possibly nested) [list] of partitions or permutations of
   a ['a seq].  *)
type 'a seq = 'a list

(* \thocwmodulesection{Simple Combinatorial Functions} *)

(* The functions
   \begin{subequations}
   \begin{align}
     \ocwlowerid{factorial}:\;& n \to n! \\
     \ocwlowerid{binomial}:\; & (n, k) \to
        \binom{n}{k} = \frac{n!}{k!(n-k)!} \\
     \ocwlowerid{multinomial}:\; & \lbrack n_1; n_2; \ldots; n_k \rbrack \to
        \binom{n_1+n_2+\ldots+n_k}{n_1,n_2,\ldots,n_k} =
        \frac{(n_1+n_2+\ldots+n_k)!}{n_1!n_2!\cdots n_k!}
   \end{align}
   \end{subequations}
   have not been optimized. They can quickly run out of the range of
   native integers. *)
val factorial : int -> int
val binomial : int -> int -> int
val multinomial : int list -> int

(* [symmetry l] returns the size of the symmetric group on~[l],
   i.\,e.~the product of the factorials of the numbers of identical
   elements. *)
val symmetry : 'a list -> int

(* \thocwmodulesection{Partitions} *)

(* $\ocwlowerid{partitions}\,
    \lbrack n_1;n_2;\ldots;n_k \rbrack\, \lbrack x_1;x_2;\ldots;x_n\rbrack$,
   where $n=n_1+n_2+\ldots+n_k$, returns all inequivalent partitions of
   $\lbrack x_1;x_2;\ldots;x_n\rbrack$ into parts of size $n_1$, $n_2$, \ldots,
   $n_k$.  The order of the $n_i$ is not respected.  There are
   \begin{equation}
     \frac{1}{S(n_1,n_2,\ldots,n_k)}
       \binom{n_1+n_2+\ldots+n_k}{n_1,n_2,\ldots,n_k}
   \end{equation}
   such partitions, where the symmetry factor~$S(n_1,n_2,\ldots,n_k)$ is
   the size of the permutation group of~$\lbrack n_1;n_2;\ldots;n_k \rbrack$
   as determined by the function [symmetry]. *)
val partitions : int list -> 'a seq -> 'a seq list list

(* [ordered_partitions] is identical to [partitions], except that the
   order of the $n_i$ is respected.  There are
   \begin{equation}
       \binom{n_1+n_2+\ldots+n_k}{n_1,n_2,\ldots,n_k}
   \end{equation}
   such partitions. *)
val ordered_partitions : int list -> 'a seq -> 'a seq list list

(* [keystones m l] is equivalent to [partitions m l], except for the
   special case when the length of~[l] is even and~[m] contains a part
   that has exactly half the length of~[l].  In this case only the half
   of the partitions is created that has the head of~[l] in the longest
   part.  *)
val keystones : int list -> 'a seq -> 'a seq list list

(* It can be beneficial to factorize a common part in the partitions and
   keystones: *)
val factorized_partitions : int list -> 'a seq -> ('a seq * 'a seq list list) list
val factorized_keystones : int list -> 'a seq -> ('a seq * 'a seq list list) list

(* \thocwmodulesubsection{Special Cases} *)

(* [partitions] is built from components that can be convenient by themselves,
   even thepugh they are just special cases of [partitions].

   [split k l] returns the list of all inequivalent splits of the list~[l] into
   one part of length~[k] and the rest. There are
   \begin{equation}
     \frac{1}{S(|l|-k,k)} \binom{|l|}{k}
   \end{equation}
   such splits.   After replacing the pairs by two-element lists,
   [split k l] is equivalent to [partitions [k; length l - k] l].*)

val split : int -> 'a seq -> ('a seq * 'a seq) list

(* Create both equipartitions of lists of even length.  There are
   \begin{equation}
     \binom{|l|}{k}
   \end{equation}
   such splits.  After replacing the pairs by two-element lists,
   the result of [ordered_split k l] is equivalent to
   [ordered_partitions [k; length l - k] l].*)

val ordered_split : int -> 'a seq -> ('a seq * 'a seq) list

(* [multi_split n k l] returns the list of all inequivalent splits of the list~[l]
   into~[n] parts of length~[k] and the rest.  *)

val multi_split : int -> int -> 'a seq -> ('a seq list * 'a seq) list
val ordered_multi_split : int -> int -> 'a seq -> ('a seq list * 'a seq) list

(* \thocwmodulesection{Choices} *)

(* $\ocwlowerid{choose}\,n\,\lbrack x_1;x_2;\ldots;x_n\rbrack$
   returns the list of all $n$-element subsets
   of~$\lbrack x_1;x_2;\ldots;x_n\rbrack$.
   [choose n] is equivalent to $(\ocwlowerid{map}\,\ocwlowerid{fst})\circ
   (\ocwlowerid{ordered\_split}\,\ocwlowerid{n})$.  *)

val choose : int -> 'a seq -> 'a seq list

(* [multi_choose n k] is equivalent to $(\ocwlowerid{map}\,\ocwlowerid{fst})\circ
   (\ocwlowerid{multi\_split}\,\ocwlowerid{n}\,\ocwlowerid{k})$.  *)

val multi_choose : int -> int -> 'a seq -> 'a seq list list
val ordered_multi_choose : int -> int -> 'a seq -> 'a seq list list

(* \thocwmodulesection{Permutations} *)

val permute : 'a seq -> 'a seq list

(* \thocwmodulesubsection{Graded Permutations} *)

val permute_signed : 'a seq -> (int * 'a seq) list
val permute_even : 'a seq -> 'a seq list
val permute_odd : 'a seq -> 'a seq list
val permute_cyclic : 'a seq -> 'a seq list
val permute_cyclic_signed : 'a seq -> (int * 'a seq) list

(* \thocwmodulesubsection{Tensor Products of Permutations} *)

(* In other words: permutations which respect compartmentalization. *)
val permute_tensor : 'a seq list -> 'a seq list list
val permute_tensor_signed : 'a seq list -> (int * 'a seq list) list
val permute_tensor_even : 'a seq list -> 'a seq list list
val permute_tensor_odd : 'a seq list -> 'a seq list list

val sign : ?cmp:('a -> 'a -> int) -> 'a seq -> int

(* \thocwmodulesubsection{Sorting} *)

val sort_signed : ?cmp:('a -> 'a -> int) -> 'a seq -> int * 'a seq

(* \thocwmodulesubsection{Unit Tests} *)

module Test : sig val suite : OUnit.test end

(*i
 *  Local Variables:
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
