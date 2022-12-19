(* tuple.mli --

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

(* The [Tuple.Poly] interface abstracts the notion of tuples with variable
   arity.  Simple cases are binary polytuples, which are simply pairs and
   indefinite polytuples, which are nothing but lists.  Another example is
   the union of pairs and triples.  The interface is very
   similar to [List] from the O'Caml standard library, but the [Tuple.Poly]
   signature allows a more fine grained control of arities. The latter
   provides typesafe linking of models, targets and topologies.  *)

module type Mono =
  sig
    type 'a t

    (* The size of the tuple, i.\,e.~[arity (a1,a2,a3) = 3]. *)
    val arity : 'a t -> int

    (* The maximum size of tuples supported by the module.
       A negative value means that there is no limit.  In this
       case the functions [power] and [power_fold] may raise
       the exception [No_termination]. *)
    val max_arity : unit -> int

    val compare : ('a -> 'a -> int) -> 'a t -> 'a t -> int

    val for_all : ('a -> bool) -> 'a t -> bool

    val map : ('a -> 'b) -> 'a t -> 'b t
    val iter : ('a -> unit) -> 'a t -> unit
    val fold_left : ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a
    val fold_right : ('a -> 'b -> 'b) -> 'a t -> 'b -> 'b

(* We have applications, where no sensible intial value can be defined: *)
    val fold_left_internal : ('a -> 'a -> 'a) -> 'a t -> 'a
    val fold_right_internal : ('a -> 'a -> 'a) -> 'a t -> 'a

    val map2 : ('a -> 'b -> 'c) -> 'a t -> 'b t -> 'c t

    val split : ('a * 'b) t -> 'a t * 'b t

(* The distributive tensor product expands a tuple of lists into
   list of tuples, e.\,g.~for binary tuples:
   \begin{equation}
     \ocwlowerid{product}\, (\lbrack x_1;x_2\rbrack,\lbrack y_1;y_2\rbrack)
       = \lbrack (x_1,y_1);(x_1,y_2);(x_2,y_1);(x_2,y_2)\rbrack
   \end{equation}
   NB: [product_fold] is usually much more memory efficient than
   the combination of [product] and [List.fold_right] for large sets.  *) 
    val product : 'a list t -> 'a t list
    val product_fold : ('a t -> 'b -> 'b) -> 'a list t -> 'b -> 'b

(* For homogeneous tuples the [power] function could trivially be built from
   [product], e.\,g.:
   \begin{equation}
     \ocwlowerid{power}\,\lbrack x_1;x_2\rbrack
       = \ocwlowerid{product}\,(\lbrack x_1;x_2\rbrack,\lbrack x_1;x_2\rbrack)
       = \lbrack (x_1,x_1);(x_1,x_2);(x_2,x_1);(x_2,x_2)\rbrack
   \end{equation}
   but it is also well defined for polytuples, e.\,g.~for pairs and triples
   \begin{equation}
     \ocwlowerid{power}\,\lbrack x_1;x_2\rbrack
       = \ocwlowerid{product}\,(\lbrack x_1;x_2\rbrack,\lbrack x_1;x_2\rbrack)
         \cup \ocwlowerid{product}\,
           (\lbrack x_1;x_2\rbrack,\lbrack x_1;x_2\rbrack,\lbrack x_1;x_2\rbrack)
   \end{equation}
   For tuples and polytuples with bounded arity, the [power]
   and [power_fold] functions terminate. In polytuples with unbounded arity, the
   the [power] function raises [No_termination] unless a limit is given
   by [?truncate].  [power_fold]
   also raises [No_termination], but could be changed to run until the
   argument function raises an exception.  However, if we need this behaviour,
   we should probably implement [power_iter] instead. *)
    val power : ?truncate:int -> 'a list -> 'a t list
    val power_fold : ?truncate:int -> ('a t -> 'b -> 'b) -> 'a list -> 'b -> 'b

(* We can also identify all (poly)tuples with permuted elements and return
   only one representative, e.\,g.:
   \begin{equation}
     \ocwlowerid{sym\_power}\,\lbrack x_1;x_2\rbrack
       = \lbrack (x_1,x_1);(x_1,x_2);(x_2,x_2)\rbrack
   \end{equation}
   NB: this function has not yet been implemented, because O'Mega only needs
   the more efficient special case [graded_sym_power].  *)

(* If a set $X$ is graded (i.\,e.~there is a map $\phi:X\to\mathbf{N}$,
   called [rank] below), the results of [power] or [sym_power] can
   canonically be filtered by requiring that the sum of the ranks in
   each (poly)tuple has one chosen value.  Implementing such a function
   directly is much more efficient than constructing and subsequently
   disregarding many (poly)tuples.  The elements of rank $n$ are at offset
   $(n-1)$ in the array.  The array is assumed to be \emph{immutable}, even
   if O'Caml doesn't support immutable arrays.  NB: [graded_power] has not
   yet been implemented, because O'Mega only needs [graded_sym_power]. *)
    type 'a graded = 'a list array
    val graded_sym_power : int -> 'a graded -> 'a t list
    val graded_sym_power_fold : int -> ('a t -> 'b -> 'b) -> 'a graded ->
      'b -> 'b

(* \begin{dubious}
     We hope to be able to avoid the next one in the long run, because it mildly
     breaks typesafety for arities.  Unfortunately, we're still working on it \ldots
   \end{dubious} *)
    val to_list : 'a t -> 'a list

(* \begin{dubious}
     The next one is only used for Fermi statistics in the obsolescent
     [Fusion_vintage] module below, but can not
     be implemented if there are no binary tuples.  It must be retired
     as soon as possible.
   \end{dubious} *)
    val of2_kludge : 'a -> 'a -> 'a t

  end

module type Poly =
  sig
    include Mono
    exception Mismatched_arity
    exception No_termination
  end

module type Binary =
    sig
      include Poly (* should become [Mono]! *)
      val of2 : 'a -> 'a -> 'a t
    end
module Binary : Binary

module type Ternary =
    sig
      include Mono 
      val of3 : 'a -> 'a -> 'a -> 'a t
    end
module Ternary : Ternary

type 'a pair_or_triple = T2 of 'a * 'a | T3 of 'a * 'a *'a

module type Mixed23 =
    sig
      include Poly
      val of2 : 'a -> 'a -> 'a t
      val of3 : 'a -> 'a -> 'a -> 'a t
    end
module Mixed23 : Mixed23

module type Nary =
    sig
      include Poly
      val of2 : 'a -> 'a -> 'a t
      val of3 : 'a -> 'a -> 'a -> 'a t
      val of_list : 'a list -> 'a t
    end
module Unbounded_Nary : Nary

(* \begin{dubious}
     It seemed like a good idea, but hardcoding [max_arity] here prevents
     optimizations for processes with fewer external particles than
     [max_arity].  For [max_arity >= 8] things become bad!
     Need to implement a truncating version of [power] and [power_fold].
   \end{dubious} *)

module type Bound = sig val max_arity : unit -> int end
module Nary (B: Bound) : Nary

(* \begin{dubious}
     For compleneteness sake, we could add most of the [List] signature
     \begin{itemize}
       \item{} [val length : 'a t -> int]
       \item{} [val hd : 'a t -> 'a]
       \item{} [val nth : 'a t -> int -> 'a]
       \item{} [val rev : 'a t -> 'a t]
       \item{} [val rev_map : ('a -> 'b) -> 'a t -> 'b t]
       \item{} [val iter2 : ('a -> 'b -> unit) -> 'a t -> 'b t -> unit]
       \item{} [val rev_map2 : ('a -> 'b -> 'c) -> 'a t -> 'b t -> 'c t]
       \item{} [val fold_left2 : ('a -> 'b -> 'c -> 'a) -> 'a -> 'b t -> 'c t -> 'a]
       \item{} [val fold_right2 : ('a -> 'b -> 'c -> 'c) -> 'a t -> 'b t -> 'c -> 'c]
       \item{} [val exists : ('a -> bool) -> 'a t -> bool]
       \item{} [val for_all2 : ('a -> 'b -> bool) -> 'a t -> 'b t -> bool]
       \item{} [val exists2 : ('a -> 'b -> bool) -> 'a t -> 'b t -> bool]
       \item{} [val mem : 'a -> 'a t -> bool]
       \item{} [val memq : 'a -> 'a t -> bool]
       \item{} [val find : ('a -> bool) -> 'a t -> 'a]
       \item{} [val find_all : ('a -> bool) -> 'a t -> 'a list]
       \item{} [val assoc : 'a -> ('a * 'b) t -> 'b]
       \item{} [val assq : 'a -> ('a * 'b) t -> 'b]
       \item{} [val mem_assoc : 'a -> ('a * 'b) t -> bool]
       \item{} [val mem_assq : 'a -> ('a * 'b) t -> bool]
       \item{} [val combine : 'a t -> 'b t -> ('a * 'b) t]
       \item{} [val sort : ('a -> 'a -> int) -> 'a t -> 'a t]
       \item{} [val stable_sort : ('a -> 'a -> int) -> 'a t -> 'a t]
     \end{itemize}
   \end{dubious}
   but only if we ever have too much time on our hand \ldots *)

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
