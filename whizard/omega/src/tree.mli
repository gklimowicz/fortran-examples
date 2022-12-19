(* tree.mli --

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

(* This module provides utilities for generic decorated trees, such as
   FeynMF output. *)

(* \thocwmodulesection{Abstract Data Type} *)
type ('n, 'l) t

(* [leaf n l] returns a tree consisting of a single leaf node
   of type [n] with a label [l]. *)
val leaf : 'n -> 'l -> ('n, 'l) t

(* [cons n ch] returns a tree node. *)
val cons : 'n -> ('n, 'l) t list -> ('n, 'l) t

(* Note that [cons node []] constructs a terminal node, but
   \emph{not} a leaf, since the latter \emph{must} have a label!
   \begin{dubious}
     \label{Tree.Leaf}
     This approach was probably tailored to Feynman diagrams,
     where we have external propagators as nodes with additional
     labels (cf.~the function [to_feynmf] on page~\pageref{Tree.to_feynmf}
     below). I'm not so sure anymore that this was a good choice.
   \end{dubious} *)

(* [node t] returns the top node of the tree [t]. *)
val node : ('n, 'l) t -> 'n

(* [leafs t] returns a list of all leaf labels \textit{in order}. *)
val leafs : ('n, 'l) t -> 'l list

(* [nodes t] returns a list of all nodes that are not leafs
   in post-order. This guarantees
   that the root node can be stripped from the result by [List.tl]. *)
val nodes :  ('n, 'l) t -> 'n list

(* [fuse conjg root contains_root trees] joins the [trees], using
   the leaf [root] in one of the trees as root of the new tree.
   [contains_root] guides the search for the subtree containing [root]
   as a leaf. [fun t -> List.mem root (leafs t)] is acceptable, but more
   efficient solutions could be available in special circumstances.  *)
val fuse : ('n -> 'n) -> 'l -> (('n, 'l) t -> bool) -> ('n, 'l) t list -> ('n, 'l) t

(* [sort lesseq t] return a sorted copy of the tree~[t]: node
   labels are ignored and nodes are according to the supremum of the
   leaf labels in the corresponding subtree. *)
val sort : ('l -> 'l -> bool) -> ('n, 'l) t -> ('n, 'l) t
val canonicalize : ('n, 'l) t -> ('n, 'l) t

(* \thocwmodulesection{Homomorphisms} *)
val map : ('n1 -> 'n2) -> ('l1 -> 'l2) -> ('n1, 'l1) t -> ('n2, 'l2) t
val fold : ('n -> 'l -> 'a) -> ('n -> 'a list -> 'a) -> ('n, 'l) t -> 'a
val fan : ('n -> 'l -> 'a list) -> ('n -> 'a list -> 'a list) ->
  ('n, 'l) t -> 'a list

(* \thocwmodulesection{Output} *)
val to_string : (string, string) t -> string

(* \thocwmodulesubsection{Feynmf} *)
(* \begin{dubious}
      [style : (string * string) option] should be replaced by
      [style : string option; tex_label : string option]
   \end{dubious} *)
type feynmf =
    { style : (string * string) option;
      rev : bool;
      label : string option;
      tension : float option } 
val vanilla : feynmf
val sty : (string * string) * bool * string -> feynmf

(* [to_feynmf file to_string incoming t] write the trees in the
   list~[t] to the file named~[file].  The leaves~[incoming] are
   used as incoming particles and~[to_string] is use to convert
   leaf labels to \LaTeX-strings. *)
(* \label{Tree.to_feynmf} *)

type 'l feynmf_set =
  { header : string;
    incoming : 'l list;
    diagrams : (feynmf, 'l) t list }

type ('l, 'm) feynmf_sets =
  { outer : 'l feynmf_set;
    inner : 'm feynmf_set list }

val feynmf_sets_plain : bool -> int -> string ->
  ('l -> string) -> ('l -> string) ->
  ('m -> string) -> ('m -> string) -> ('l, 'm) feynmf_sets list -> unit

val feynmf_sets_wrapped : bool -> string ->
  ('l -> string) -> ('l -> string) ->
  ('m -> string) -> ('m -> string) -> ('l, 'm) feynmf_sets list -> unit

(* If the diagrams at all levels are of the same type,
   we can recurse to arbitrary depth. *)

type 'l feynmf_levels =
  { this : 'l feynmf_set;
    lower : 'l feynmf_levels list }

(* [to_feynmf_levels_plain sections level file wf_to_TeX p_to_TeX levels]
   \ldots *)

val feynmf_levels_plain : bool -> int -> string ->
  ('l -> string) -> ('l -> string) -> 'l feynmf_levels list -> unit

(* [to_feynmf_levels_wrapped file wf_to_TeX p_to_TeX levels]
   \ldots *)

val feynmf_levels_wrapped : string ->
  ('l -> string) -> ('l -> string) -> 'l feynmf_levels list -> unit

(* \thocwmodulesubsection{Least Squares Layout} *)

(* A general graph with edges of type~['e], internal nodes of type~['n],
   and external nodes of type ['ext].  *)
type ('e, 'n, 'ext) graph
val graph_of_tree : ('n -> 'n -> 'e) -> ('n -> 'n) ->
  'n -> ('n, 'n) t -> ('e, 'n, 'n) graph

(* A general graph with the layout of the external nodes fixed.  *)
type ('e, 'n, 'ext) ext_layout
val left_to_right : int -> ('e, 'n, 'ext) graph -> ('e, 'n, 'ext) ext_layout

(* A general graph with the layout of all nodes fixed.  *)
type ('e, 'n, 'ext) layout
val layout : ('e, 'n, 'ext) ext_layout -> ('e, 'n, 'ext) layout

val dump : ('e, 'n, 'ext) layout -> unit
val iter_edges : ('e -> float * float -> float * float -> unit) ->
  ('e, 'n, 'ext) layout -> unit
val iter_internal : (float * float -> unit) ->
  ('e, 'n, 'ext) layout -> unit
val iter_incoming : ('ext * float * float -> unit) ->
  ('e, 'n, 'ext) layout -> unit
val iter_outgoing : ('ext * float * float -> unit) ->
  ('e, 'n, 'ext) layout -> unit

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
