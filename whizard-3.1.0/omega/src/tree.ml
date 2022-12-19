(* tree.ml --

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

(* \thocwmodulesection{Abstract Data Type} *)

type ('n, 'l) t =
  | Leaf of 'n * 'l
  | Node of 'n * ('n, 'l) t list

let leaf n l = Leaf (n, l)

let cons n children = Node (n, children)

(* Presenting the leafs \textit{in order} comes naturally, but will be
   useful below. *)
let rec leafs = function
  | Leaf (_, l) -> [l]
  | Node (_, ch) -> ThoList.flatmap leafs ch

let node = function
  | Leaf (n, _) -> n
  | Node (n, _) -> n

(* This guarantees that the root node can be stripped from the result
   by [List.tl]. *)
let rec nodes = function
  | Leaf _ -> []
  | Node (n, ch) -> n :: ThoList.flatmap nodes ch

(* [first_match p list] returns [(x,list')], where [x] is the first element
   of [list] for which [p x = true] and [list'] is [list] sans [x]. *)
let first_match p list =
  let rec first_match' no_match = function
    | [] -> invalid_arg "Tree.fuse: prospective root not found"
    | t :: rest when p t -> (t, List.rev_append no_match rest)
    | t :: rest -> first_match' (t :: no_match) rest in
  first_match' [] list

(* One recursion step in [fuse'] rotates the topmost tree node, moving
   the prospective root up:
   \begin{equation}
   \label{eq:tree-rotation}
     \parbox{46\unitlength}{%
       \fmfframe(0,0)(0,4){%
         \begin{fmfgraph*}(45,30)
           \fmfstraight
           \fmftop{r}
           \fmfbottom{l11,l12,l1x,l1n,db1,l21,l22,l2x,l2n,db2,db3,db4,db5,db6,%
                      lx1,lx2,lxx,lxn,db7,ln1,ln2,lnx,lnn}
           \fmf{plain,tension=4}{r,vr1}
           \fmf{plain,tension=4,lab=$p$,lab.side=left}{r,vr2}
           \fmf{dots,tension=4}{r,vrx}
           \fmf{plain,tension=4}{r,vrn}
           \fmf{plain}{vr1,l11}\fmf{plain}{vr1,l12}
             \fmf{dots}{vr1,l1x}\fmf{plain}{vr1,l1n}
           \fmf{plain}{vr2,l21}\fmf{plain}{vr2,l22}
             \fmf{dots}{vr2,l2x}\fmf{plain}{vr2,l2n}
           \fmf{dots}{vrx,lx1}\fmf{dots}{vrx,lx2}
             \fmf{dots}{vrx,lxx}\fmf{dots}{vrx,lxn}
           \fmf{plain}{vrn,ln1}\fmf{plain}{vrn,ln2}
             \fmf{dots}{vrn,lnx}\fmf{plain}{vrn,lnn}
           \fmfv{l=$r$,l.ang=-90}{l22}
           \fmfv{d.shape=circle,d.filled=empty,d.size=7thick,%
                 back=.8white}{r,vr1,vrx,vrn}
           \fmfv{d.shape=circle,d.filled=empty,d.size=7thick,%
                 lab=$R$,lab.dist=0}{vr2}
         \end{fmfgraph*}}}
     \to
     \parbox{61\unitlength}{%
       \fmfframe(0,0)(0,4){%
         \begin{fmfgraph*}(60,30)
           \fmfstraight
           \fmftop{r}
           \fmfbottom{l21,d1,d2,l22,d3,d4,l2x,d5,d6,l2n,d7,d8,db2,%
                      l11,l12,l1x,l1n,db1,db2,db3,lx1,lx2,lxx,lxn,db4,%
                      ln1,ln2,lnx,lnn}
           \fmf{plain}{r,vr1}\fmf{phantom}{vr1,l21}
           \fmf{plain}{r,vr2}\fmf{phantom}{vr2,l22}
           \fmf{dots}{r,vrx}\fmf{phantom}{vrx,l2x}
           \fmf{plain}{r,vr3}\fmf{phantom}{vr3,l2n}
           \fmf{plain,tension=12,lab=$-p$,lab.side=left}{r,vrn}
           \fmf{plain,tension=4}{vrn,vvr1}
           \fmf{dots,tension=4}{vrn,vvrx}
           \fmf{plain,tension=4}{vrn,vvrn}
           \fmf{plain}{vvr1,l11}\fmf{plain}{vvr1,l12}
             \fmf{dots}{vvr1,l1x}\fmf{plain}{vvr1,l1n}
           \fmf{dots}{vvrx,lx1}\fmf{dots}{vvrx,lx2}
             \fmf{dots}{vvrx,lxx}\fmf{dots}{vvrx,lxn}
           \fmf{plain}{vvrn,ln1}\fmf{plain}{vvrn,ln2}
             \fmf{dots}{vvrn,lnx}\fmf{plain}{vvrn,lnn}
           \fmfv{l=$r$,l.ang=-90}{vr2}
           \fmfv{d.shape=circle,d.filled=empty,d.size=7thick,%
                 back=.8white}{vrn,vvr1,vvrx,vvrn}
           \fmfv{d.shape=circle,d.filled=empty,d.size=7thick,%
                 lab=$R$,lab.dist=0}{r}
         \end{fmfgraph*}}}
   \end{equation} *)

let fuse conjg root contains_root trees =
  let rec fuse' subtrees =
    match first_match contains_root subtrees with

(* If the prospective root is contained in a leaf, we have either found
   the root---in which case we're done---or have failed catastrophically: *)
    | Leaf (n, l), children ->
        if l = root then
          Node (conjg n, children)
        else
          invalid_arg "Tree.fuse: root predicate inconsistent"

(* Otherwise, we perform a rotation as in~(\ref{eq:tree-rotation}) and
   connect all nodes that do not contain the root to a new node.
   For efficiency, we append the new node at the end and prevent
   [first_match] from searching for the root in it in vain again.
   Since [root_children] is probably rather short, this should be
   a good strategy. *)
    | Node (n, root_children), other_children ->
        fuse' (root_children @ [Node (conjg n, other_children)]) in
  fuse' trees

(* Sorting is also straightforward, we only have to keep track of the
   suprema of the subtrees: *)

type ('a, 'b) with_supremum = { sup : 'a; data : 'b }

(* Since the lists are rather short, [List.sort] could be replaced by
   an optimized version, but we're not (yet) dealing with the most
   important speed bottleneck here: *)

let rec sort' lesseq = function
  | Leaf (_, l) as e -> { sup = l; data = e }
  | Node (n, ch) ->
      let ch' = List.sort
          (fun x y -> compare x.sup y.sup) (List.map (sort' lesseq) ch) in
      { sup = (List.hd (List.rev ch')).sup;
        data = Node (n, List.map (fun x -> x.data) ch') }

(* finally, throw away the overall supremum: *)

let sort lesseq t = (sort' lesseq t).data

let rec canonicalize = function
  | Leaf (_, _) as l -> l
  | Node (n, ch) ->
    Node (n, List.sort compare (List.map canonicalize ch))

(* \thocwmodulesection{Homomorphisms} *)

(* Isomophisms are simple: *)

let rec map fn fl = function
  | Leaf (n, l) -> Leaf (fn n, fl l)
  | Node (n, ch) -> Node (fn n, List.map (map fn fl) ch)

(* homomorphisms are not more complicated: *)

let rec fold leaf node = function
  | Leaf (n, l) -> leaf n l
  | Node (n, ch) -> node n (List.map (fold leaf node) ch)

(* and tensor products are fun: *)

let rec fan leaf node = function
  | Leaf (n, l) -> leaf n l
  | Node (n, ch) -> Product.fold
        (fun ch' t -> node n ch' @ t) (List.map (fan leaf node) ch) []

(* \thocwmodulesection{Output} *)

let leaf_to_string n l =
  if n = "" then
    l
  else if l = "" then
    n
  else
    n ^ "(" ^ l ^ ")"

let node_to_string n ch =
  "(" ^ (if n = "" then "" else n ^ ":") ^ (String.concat "," ch) ^ ")"

let to_string t =
  fold leaf_to_string node_to_string t

(* \thocwmodulesubsection{Feynmf}
   Add a value that is greater than all suprema *)

type 'a supremum_or_infinity = Infinity | Sup of 'a

type ('a, 'b) with_supremum_or_infinity =
    { sup : 'a supremum_or_infinity; data : 'b }

let with_infinity cmp x y =
  match x.sup, y.sup with
  | Infinity, _ -> 1
  | _, Infinity -> -1
  | Sup x', Sup y' -> cmp x' y'

(* Using this, we can sort the tree in another way that guarantees that
   a particular leaf ([i2]) is moved as far to the end as possible.  We
   can then flip this leaf from outgoing to incoming without introducing
   a crossing: *)

let rec sort_2i' lesseq i2 = function
  | Leaf (_, l) as e ->
      { sup = if l = i2 then Infinity else Sup l; data = e }
  | Node (n, ch) ->
      let ch' = List.sort (with_infinity compare)
          (List.map (sort_2i' lesseq i2) ch) in
      { sup = (List.hd (List.rev ch')).sup;
        data = Node (n, List.map (fun x -> x.data) ch') }

(* again, throw away the overall supremum: *)

let sort_2i lesseq i2 t = (sort_2i' lesseq i2 t).data

type feynmf =
    { style : (string * string) option;
      rev : bool;
      label : string option;
      tension : float option } 

open Printf

let style prop =
  match prop.style with
  | None -> ("plain","")
  | Some s -> s

let species prop = fst (style prop)
let tex_lbl prop = snd (style prop)

let leaf_label tex io leaf lab = function
  | None -> fprintf tex "    \\fmflabel{${%s}$}{%s%s}\n" lab io leaf 
  | Some s ->
      fprintf tex "    \\fmflabel{${%s{}^{(%s)}}$}{%s%s}\n" s lab io leaf

let leaf_label tex io leaf lab label =
  ()

(* We try to draw diagrams more symmetrically by reducing the tension
   on the outgoing external lines.
   \begin{dubious}
   \index{shortcomings!algorithmical}
      This is insufficient for asymmetrical cascade decays.
   \end{dubious} *)

let rec leaf_node tex to_label i2 n prop leaf =
  let io, tension, rev =
    if leaf = i2 then
      ("i", "", not prop.rev)
    else
      ("o", ",tension=0.5", prop.rev) in
  leaf_label tex io (to_label leaf) (tex_lbl prop) prop.label ;
  fprintf tex "    \\fmfdot{v%d}\n"  n;
  if rev then 
    fprintf tex "    \\fmf{%s%s}{%s%s,v%d}\n"
      (species prop) tension io (to_label leaf) n
  else
    fprintf tex "    \\fmf{%s%s}{v%d,%s%s}\n"
      (species prop) tension n io (to_label leaf)

and int_node tex to_label i2 n n' prop t =
  if prop.rev then
    fprintf tex 
      "    \\fmf{%s,label=\\begin{scriptsize}${%s}$\\end{scriptsize}}{v%d,v%d}\n" 
      (species prop) (tex_lbl prop) n' n
  else
    fprintf tex 
      "    \\fmf{%s,label=\\begin{scriptsize}${%s}$\\end{scriptsize}}{v%d,v%d}\n" 
      (species prop) (tex_lbl prop) n n';
  fprintf tex "    \\fmfdot{v%d,v%d}\n" n n';
  edges_feynmf' tex to_label i2 n' t

and leaf_or_int_node tex to_label i2 n n' = function
  | Leaf (prop, l) -> leaf_node tex to_label i2 n prop l
  | Node (prop, _) as t -> int_node tex to_label i2 n n' prop t

and edges_feynmf' tex to_label i2 n = function
  | Leaf (prop, l) -> leaf_node tex to_label i2 n prop l
  | Node (_, ch) ->
      ignore (List.fold_right
                (fun t' n' ->
                  leaf_or_int_node tex to_label i2 n n' t';
                  succ n') ch (4*n))

let edges_feynmf tex to_label i1 i2 t =
  let n = 1 in
  begin match t with
  | Leaf _ -> ()
  | Node (prop, _) ->
      leaf_label tex "i" "1" (tex_lbl prop) prop.label;
      if prop.rev then
        fprintf tex "    \\fmf{%s}{v%d,i%s}\n" (species prop) n (to_label i1)
      else
        fprintf tex "    \\fmf{%s}{i%s,v%d}\n" (species prop) (to_label i1) n
  end;
  fprintf tex "    \\fmfdot{v%d}\n" n;
  edges_feynmf' tex to_label i2 n t

let to_feynmf_channel tex to_TeX to_label incoming t =
  match incoming with
  | i1 :: i2 :: _ ->
      let t' = sort_2i (<=) i2 t in
      let out = List.filter (fun a -> i2 <> a) (leafs t') in
      fprintf tex "\\fmfframe(8,7)(8,6){%%\n";
      fprintf tex "  \\begin{fmfgraph*}(35,30)\n";
      fprintf tex "   \\fmfpen{thin}\n";
      fprintf tex "   \\fmfset{arrow_len}{2mm}\n";
      fprintf tex "    \\fmfleft{i%s,i%s}\n" (to_label i1) (to_label i2);
      fprintf tex "    \\fmfright{o%s}\n"
        (String.concat ",o" (List.map to_label out));
      List.iter
        (fun s ->
          fprintf tex "    \\fmflabel{${%s}$}{i%s}\n"
            (to_TeX s) (to_label s))
        [i1; i2];
      List.iter
        (fun s ->
          fprintf tex "    \\fmflabel{${%s}$}{o%s}\n"
            (to_TeX s) (to_label s))
        out;
      edges_feynmf tex to_label i1 i2 t';
      fprintf tex "  \\end{fmfgraph*}}\\hfil\\allowbreak\n"
  | _ -> ()

(* \begin{figure}
   \fmfframe(3,5)(3,5){%
     \begin{fmfgraph*}(30,30)
       \fmfleft{i1,i2}
       \fmfright{o3,o4,o5,o6}
       \fmflabel{$1$}{i1}
       \fmflabel{$2$}{i2}
       \fmflabel{$3$}{o3}
       \fmflabel{$4$}{o4}
       \fmflabel{$5$}{o5}
       \fmflabel{$6$}{o6}
     \fmf{plain}{i1,v1}
     \fmf{plain}{v1,v3}
     \fmf{plain,tension=0.5}{v3,o3}
     \fmf{plain}{v3,v9}
     \fmf{plain,tension=0.5}{v9,o4}
     \fmf{plain}{v9,v27}
     \fmf{plain,tension=0.5}{v27,o5}
     \fmf{plain,tension=0.5}{v27,o6}
     \fmf{plain}{v1,i2}
     \end{fmfgraph*}}
   \fmfframe(3,5)(3,5){%
     \begin{fmfgraph*}(30,30)
       \fmfleft{i1,i2}
       \fmfright{o3,o4,o6,o5}
       \fmflabel{$1$}{i1}
       \fmflabel{$2$}{i2}
       \fmflabel{$3$}{o3}
       \fmflabel{$4$}{o4}
       \fmflabel{$6$}{o6}
       \fmflabel{$5$}{o5}
     \fmf{plain}{i1,v1}
     \fmf{plain}{v1,v3}
     \fmf{plain,tension=0.5}{v3,o3}
     \fmf{plain}{v3,v9}
     \fmf{plain}{v9,v27}
     \fmf{plain,tension=0.5}{v27,o4}
     \fmf{plain,tension=0.5}{v27,o6}
     \fmf{plain,tension=0.5}{v9,o5}
     \fmf{plain}{v1,i2}
     \end{fmfgraph*}}
   \fmfframe(3,5)(3,5){%
     \begin{fmfgraph*}(30,30)
       \fmfleft{i1,i2}
       \fmfright{o3,o4,o5,o6}
       \fmflabel{$1$}{i1}
       \fmflabel{$2$}{i2}
       \fmflabel{$3$}{o3}
       \fmflabel{$4$}{o4}
       \fmflabel{$5$}{o5}
       \fmflabel{$6$}{o6}
     \fmf{plain}{i1,v1}
     \fmf{plain}{v1,v3}
     \fmf{plain}{v3,v9}
     \fmf{plain,tension=0.5}{v9,o3}
     \fmf{plain,tension=0.5}{v9,o4}
     \fmf{plain}{v3,v10}
     \fmf{plain,tension=0.5}{v10,o5}
     \fmf{plain,tension=0.5}{v10,o6}
     \fmf{plain}{v1,i2}
     \end{fmfgraph*}}
     \caption{\label{fig:to_feynmf}%
       Note that this is subtly different \ldots}
   \end{figure} *)

let vanilla = { style = None; rev = false; label = None; tension = None }

let sty (s, r, l) = { vanilla with style = Some s; rev = r; label = Some l }

type 'l feynmf_set =
  { header : string;
    incoming : 'l list;
    diagrams : (feynmf, 'l) t list }

type ('l, 'm) feynmf_sets =
  { outer : 'l feynmf_set;
    inner : 'm feynmf_set list }

type 'l feynmf_levels =
  { this : 'l feynmf_set;
    lower : 'l feynmf_levels list }

let latex_section = function
  | level when level < 0 -> "part"
  | 0 -> "chapter"
  | 1 -> "section"
  | 2 -> "subsection"
  | 3 -> "subsubsection"
  | 4 -> "paragraph"
  | _ -> "subparagraph"

let rec feynmf_set tex sections level to_TeX to_label set =
  fprintf tex "%s\\%s{%s}\n"
    (if sections then "" else "%%% ")
    (latex_section level)
    set.header;
  List.iter
    (to_feynmf_channel tex to_TeX to_label set.incoming)
    set.diagrams

let feynmf_sets tex sections level
    to_TeX_outer to_label_outer to_TeX_inner to_label_inner set =
  feynmf_set tex sections level to_TeX_outer to_label_outer set.outer;
  List.iter
    (feynmf_set tex sections (succ level) to_TeX_inner to_label_inner)
    set.inner

let feynmf_sets_plain sections level file
    to_TeX_outer to_label_outer to_TeX_inner to_label_inner sets =
  let tex = open_out (file ^ ".tex") in
  List.iter
    (feynmf_sets tex sections level
       to_TeX_outer to_label_outer to_TeX_inner to_label_inner)
    sets;
  close_out tex

let feynmf_header tex file =
  fprintf tex "\\documentclass[10pt]{article}\n";
  fprintf tex "\\usepackage{ifpdf}\n";
  fprintf tex "\\usepackage[colorlinks]{hyperref}\n";
  fprintf tex "\\usepackage[a4paper,margin=1cm]{geometry}\n";
  fprintf tex "\\usepackage{feynmp}\n";
  fprintf tex "\\ifpdf\n";
  fprintf tex "   \\DeclareGraphicsRule{*}{mps}{*}{}\n";
  fprintf tex "\\else\n";
  fprintf tex "   \\DeclareGraphicsRule{*}{eps}{*}{}\n";
  fprintf tex "\\fi\n";
  fprintf tex "\\setlength{\\unitlength}{1mm}\n";
  fprintf tex "\\setlength{\\parindent}{0pt}\n";
  fprintf tex
    "\\renewcommand{\\mathstrut}{\\protect\\vphantom{\\hat{0123456789}}}\n";
  fprintf tex "\\begin{document}\n";
  fprintf tex "\\tableofcontents\n";
  fprintf tex "\\begin{fmffile}{%s-fmf}\n\n" file

let feynmf_footer tex =
  fprintf tex "\n";   
  fprintf tex "\\end{fmffile} \n";
  fprintf tex "\\end{document} \n"

let feynmf_sets_wrapped latex file 
    to_TeX_outer to_label_outer to_TeX_inner to_label_inner sets =
  let tex = open_out (file ^ ".tex") in
  if latex then feynmf_header tex file;
  List.iter
    (feynmf_sets tex latex 1
       to_TeX_outer to_label_outer to_TeX_inner to_label_inner)
    sets;
  if latex then feynmf_footer tex;
  close_out tex

let rec feynmf_levels tex sections level to_TeX to_label set =
  fprintf tex "%s\\%s{%s}\n"
    (if sections then "" else "%%% ")
    (latex_section level)
    set.this.header;
  List.iter
    (to_feynmf_channel tex to_TeX to_label set.this.incoming)
    set.this.diagrams;
  List.iter (feynmf_levels tex sections (succ level) to_TeX to_label) set.lower

let feynmf_levels_plain sections level file to_TeX to_label sets =
  let tex = open_out (file ^ ".tex") in
  List.iter (feynmf_levels tex sections level to_TeX to_label) sets;
  close_out tex
    
let feynmf_levels_wrapped file to_TeX to_label sets =
  let tex = open_out (file ^ ".tex") in
  feynmf_header tex file;
  List.iter (feynmf_levels tex true 1 to_TeX to_label) sets;
  feynmf_footer tex;
  close_out tex
    
(* \thocwmodulesection{Least Squares Layout}
   \begin{equation}
     L = \frac{1}{2} \sum_{i\not=i'} T_{ii'} \left(x_i-x_{i'}\right)^2
       + \frac{1}{2} \sum_{i,j} T'_{ij} \left(x_i-e_j\right)^2
   \end{equation}
   and thus
   \begin{equation}
     0 = \frac{\partial L}{\partial x_i}
       = \sum_{i'\not=i} T_{ii'} \left(x_i-x_{i'}\right)
       + \sum_{j} T'_{ij} \left(x_i-e_j\right)
   \end{equation}
   or
   \begin{equation}
   \label{eq:layout}
         \left(\sum_{i'\not=i} T_{ii'} + \sum_{j} T'_{ij}\right) x_i
       - \sum_{i'\not=i} T_{ii'} x_{i'}
       = \sum_{j} T'_{ij} e_j
   \end{equation}
   where we can assume that
   \begin{subequations}
   \begin{align}
      T_{ii'} &= T_{i'i} \\
      T_{ii} &= 0
   \end{align}
   \end{subequations} *)
type 'a node_with_tension = { node : 'a; tension : float }

let unit_tension t =
  map (fun n -> { node = n; tension = 1.0 }) (fun l -> l) t

let leafs_and_nodes i2 t =
  let t' = sort_2i (<=) i2 t in
  match nodes t' with
  | [] -> failwith "Tree.nodes_and_leafs: impossible"
  | i1 :: _ as n -> (i1, i2, List.filter (fun l -> l <> i2) (leafs t'), n)

(* Not tail recursive, but they're unlikely to meet any deep trees: *)
let rec internal_edges_from n = function
  | Leaf _ -> []
  | Node (n', ch) -> (n', n) :: (ThoList.flatmap (internal_edges_from n') ch)

(* The root node of the tree represents a vertex (node) and an
   external line (leaf) of the Feynman diagram simultaneously.  Thus
   it requires special treatment: *)
let internal_edges = function
  | Leaf _ -> []
  | Node (n, ch) -> ThoList.flatmap (internal_edges_from n) ch

let rec external_edges_from n = function
  | Leaf (n', _) -> [(n', n)]
  | Node (n', ch) -> ThoList.flatmap (external_edges_from n') ch

let external_edges = function
  | Leaf (n, _) -> [(n, n)]
  | Node (n, ch) -> (n, n) :: ThoList.flatmap (external_edges_from n) ch

type ('edge, 'node, 'ext) graph =
    { int_nodes : 'node array;
      ext_nodes : 'ext array;
      int_edges : ('edge * int * int) list;
      ext_edges : ('edge * int * int) list } 

module M = Pmap.Tree

(* Invert an array, viewed as a map from non-negative integers
   into a set. The result is a map from the set to the integers:
   [val invert_array : 'a array -> ('a, int) M.t] *)

let invert_array_unsafe a =
  fst (Array.fold_left (fun (m, i) a_i ->
    (M.add compare a_i i m, succ i)) (M.empty, 0) a)

exception Not_invertible

let add_unique key data map =
  if M.mem compare key map then
    raise Not_invertible
  else
    M.add compare key data map

let invert_array a =
  fst (Array.fold_left (fun (m, i) a_i ->
    (add_unique a_i i m, succ i)) (M.empty, 0) a)

let graph_of_tree nodes2edge conjugate i2 t =
  let i1, i2, out, vertices = leafs_and_nodes i2 t in
  let int_nodes = Array.of_list vertices
  and ext_nodes = Array.of_list (conjugate i1 :: i2 :: out) in
  let int_nodes_index_table = invert_array int_nodes
  and ext_nodes_index_table = invert_array ext_nodes in
  let int_nodes_index n = M.find compare n int_nodes_index_table
  and ext_nodes_index n = M.find compare n ext_nodes_index_table in
  { int_nodes = int_nodes;
    ext_nodes = ext_nodes;
    int_edges = List.map
      (fun (n1, n2) ->
        (nodes2edge n1 n2, int_nodes_index n1, int_nodes_index n2))
      (internal_edges t);
    ext_edges = List.map
      (fun (e, n) ->
        let e' = 
          if e = i1 then
            conjugate e
          else
            e in
        (nodes2edge e' n, ext_nodes_index e', int_nodes_index n))
      (external_edges t) } 

let int_incidence f null g =
  let n = Array.length g.int_nodes in
  let incidence = Array.make_matrix n n null in
  List.iter (fun (edge, n1, n2) ->
    if n1 <> n2 then begin
      let edge' = f edge g.int_nodes.(n1) g.int_nodes.(n2) in
      incidence.(n1).(n2) <- edge';
      incidence.(n2).(n1) <- edge'
    end)
    g.int_edges;
  incidence

let ext_incidence f null g =
  let n_int = Array.length g.int_nodes
  and n_ext = Array.length g.ext_nodes in
  let incidence = Array.make_matrix n_int n_ext null in
  List.iter (fun (edge, e, n) ->
    incidence.(n).(e) <- f edge g.ext_nodes.(e) g.int_nodes.(n))
    g.ext_edges;
  incidence

let division n =
  if n < 0 then
    []
  else if n = 1 then
    [0.5]
  else
    let n' = pred n in
    let d = 1.0 /. (float n') in
    let rec division' i acc =
      if i < 0 then
        acc
      else
        division' (pred i) (float i *. d :: acc) in
    division' n' []

type ('e, 'n, 'ext) ext_layout = ('e, 'n, 'ext * float * float) graph
type ('e, 'n, 'ext) layout = ('e, 'n * float * float, 'ext) ext_layout

let left_to_right num_in g =
  if num_in < 1 then
    invalid_arg "left_to_right"
  else
    let num_out = Array.length g.ext_nodes - num_in in
    if num_out < 1 then
      invalid_arg "left_to_right"
    else
      let incoming =
        List.map2 (fun e y -> (e, 0.0, y))
          (Array.to_list (Array.sub g.ext_nodes 0 num_in))
          (division num_in)
      and outgoing =
        List.map2 (fun e y -> (e, 1.0, y))
          (Array.to_list (Array.sub g.ext_nodes num_in num_out))
          (division num_out) in
      { g with ext_nodes = Array.of_list (incoming @ outgoing) }

(* Reformulating~(\ref{eq:layout})
   \begin{subequations}
   \begin{align}
     Ax &= b_x \\
     Ay &= b_y
   \end{align}
   \end{subequations}
   with
   \begin{subequations}
   \begin{align}
      A_{ii'} &=
        \left( \sum_{i''\not=i} T_{ii''}
                 + \sum_j T'_{ij} \right) \delta_{ii'} - T_{ii'} \\
      (b_{x/y})_i &= \sum_j T'_{ij} (e_{x/y})_j
   \end{align}
   \end{subequations} *)
let sum a = Array.fold_left (+.) 0.0 a
  
let tension_to_equation t t' e =
  let xe, ye = List.split e in
  let bx = Linalg.matmulv t' (Array.of_list xe)
  and by = Linalg.matmulv t' (Array.of_list ye)
  and a = Array.init (Array.length t)
      (fun i ->
        let a_i = Array.map (~-.) t.(i) in
        a_i.(i) <- a_i.(i) +. sum t.(i) +. sum t'.(i);
        a_i) in
  (a, bx, by)
  
let layout g =
  let ext_nodes =
    List.map (fun (_, x, y) -> (x, y)) (Array.to_list g.ext_nodes) in
  let a, bx, by =
    tension_to_equation
      (int_incidence (fun _ _ _ -> 1.0) 0.0 g)
      (ext_incidence (fun _ _ _ -> 1.0) 0.0 g) ext_nodes in
  match Linalg.solve_many a [bx; by] with
  | [x; y] -> { g with int_nodes = Array.mapi
                  (fun i n -> (n, x.(i), y.(i))) g.int_nodes }
  | _ -> failwith "impossible"

let iter_edges f g =
  List.iter (fun (edge, n1, n2) ->
    let _, x1, y1 = g.int_nodes.(n1)
    and _, x2, y2 = g.int_nodes.(n2) in
    f edge (x1, y1) (x2, y2)) g.int_edges;
  List.iter (fun (edge, e, n) ->
    let _, x1, y1 = g.ext_nodes.(e)
    and _, x2, y2 = g.int_nodes.(n) in
    f edge (x1, y1) (x2, y2)) g.ext_edges
  
let iter_internal f g =
  Array.iter (fun (node, x, y) -> f (x, y)) g.int_nodes
  
let iter_incoming f g =
  f g.ext_nodes.(0);
  f g.ext_nodes.(1)

let iter_outgoing f g =
  for i = 2 to pred (Array.length g.ext_nodes) do
    f g.ext_nodes.(i)
  done
  
let dump g =
  Array.iter (fun (_, x, y) -> Printf.eprintf "(%g,%g) " x y) g.ext_nodes;
  Printf.eprintf "\n => ";
  Array.iter (fun (_, x, y) -> Printf.eprintf "(%g,%g) " x y) g.int_nodes;
  Printf.eprintf "\n"

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
