(* powSet.ml --

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

module type Ordered_Type =
  sig
    type t
    val compare : t -> t -> int
    val to_string : t -> string
  end

module type T =
  sig
    type elt
    type t
    val empty : t
    val is_empty : t -> bool
    val union : t list -> t
    val of_lists : elt list list -> t
    val to_lists : t -> elt list list
    val basis : t -> t
    val to_string : t -> string
  end

module Make (E : Ordered_Type) =
  struct

    type elt = E.t

    module ESet = Set.Make (E)
    type set = ESet.t

    module EPowSet = Set.Make (ESet)
    type t = EPowSet.t

    let empty = EPowSet.empty
    let is_empty = EPowSet.is_empty
(*i let elements = EPowSet.elements i*)

    let union s_list =
      List.fold_right EPowSet.union s_list EPowSet.empty

    let set_to_string set =
      "{" ^ String.concat "," (List.map E.to_string (ESet.elements set)) ^ "}"

    let to_string powset =
      "{" ^ String.concat "," (List.map set_to_string (EPowSet.elements powset)) ^ "}"

    let set_of_list = ESet.of_list

    let of_lists lists =
      List.fold_right
        (fun list acc -> EPowSet.add (ESet.of_list list) acc)
        lists EPowSet.empty

    let to_lists ps =
      List.map ESet.elements (EPowSet.elements ps)

(* [product] $(s_1,s_2) = s_1 \circ s_2 =
     \{s_1\setminus s_2, s_1 \cap s_2, s_2\setminus s_1\} \setminus \{\emptyset\}$ *)
    let product s1 s2 =
      List.fold_left
        (fun pset set -> if ESet.is_empty set then pset else EPowSet.add set pset)
        EPowSet.empty [ESet.diff s1 s2; ESet.inter s1 s2; ESet.diff s2 s1]

(*i let product s1 s2 =
      Printf.eprintf "product %s %s" (set_to_string s1) (set_to_string s2);
      flush stderr;
      let result = product s1 s2 in
      Printf.eprintf " => %s\n" (to_string result);
      flush stderr;
      result i*)

    let disjoint s1 s2 =
      ESet.is_empty (ESet.inter s1 s2)

(* In [augment_basis_overlapping] $(s, \{s_i\}_i)$, we are guaranteed
   that
   \begin{subequations}
   \begin{align}
   \label{eq:powset:overlap}
     \forall_i        :\;& s  \cap s_i\not=\emptyset\\
   \label{eq:powset:disjoint}
     \forall_{i\not=j}:\;& s_i\cap s_j    =\emptyset\,.
   \end{align}
   \end{subequations}
   Therefore from~(\ref{eq:powset:disjoint}) 
   \begin{subequations}
   \begin{align}
     \forall_{i\not=j}:\;& (s  \cap      s_i) \cap (s  \cap      s_j)
                       = s \cap (s_i \cap s_j) = s \cap \emptyset = \emptyset\\
     \forall_{i\not=j}:\;& (s_i\setminus s  ) \cap (s_j\setminus s  )
                       \subset s_i \cap s_j = \emptyset\\
     \forall_{i\not=j}:& (s  \setminus s_i) \cap (s_j\setminus s  )
                       \subset s \cap \bar s = \emptyset\\
     \forall_{i\not=j}:& (s  \cap      s_i) \cap (s_j\setminus s  )
                       \subset s \cap \bar s = \emptyset\,,
   \end{align}
   \end{subequations}
   but in general
   \begin{subequations}
   \begin{align}
     \exists_{i\not=j} :& (s  \setminus s_i) \cap (s  \setminus s_j) \not=\emptyset\\
     \exists_{i\not=j}:& (s  \setminus s_i) \cap (s  \cap      s_j) \not=\emptyset\,,
   \end{align}
   \end{subequations}
   because, e.\,g., for $s_i=\{i\}$ and $s=\{1,2,3\}$
   \begin{subequations}
   \begin{align}
     (s  \setminus s_1) \cap (s  \setminus s_2) &= \{2,3\} \cap \{1,3\} = \{3\} \\
     (s  \setminus s_1) \cap (s  \cap      s_2) &= \{2,3\} \cap \{2\}   = \{2\}\,.
   \end{align}
   \end{subequations}
   Summarizing:
   \begin{center}
     \begin{tabular}{c||c|c|c}
        $\forall_{i\not=j}:\;A_i\cap A_j$&$s_j\setminus s  $&$s  \cap s_j   $&$s  \setminus s_j$\\
          \hline\hline
                       $s_i\setminus s  $&$\emptyset       $&$\emptyset     $&$\emptyset       $\\
          \hline
                       $s  \cap      s_i$&$\emptyset       $&$\emptyset     $&$\not=\emptyset  $\\
          \hline
                       $s  \setminus s_i$&$\emptyset       $&$\not=\emptyset$&$\not=\emptyset  $
     \end{tabular}
   \end{center}
   Fortunately, we also know from~(\ref{eq:powset:overlap}) that
   \begin{subequations}
   \begin{align}
     \forall_i:\;& |s  \setminus s_i| < |s| \\
     \forall_i:\;& |s  \cap      s_i| < \min(|s|,|s_i|) \\
     \forall_i:\;& |s_i\setminus s  | < |s_i|
   \end{align}
   \end{subequations}
   and can call [basis] recursively without risking non-termination. *)

    let rec basis ps =
      EPowSet.fold augment_basis ps EPowSet.empty

    and augment_basis s ps =
      if EPowSet.mem s ps then
        ps
      else
        let no_overlaps, overlaps = EPowSet.partition (disjoint s) ps in
        if EPowSet.is_empty overlaps then
          EPowSet.add s ps
        else
          EPowSet.union no_overlaps (augment_basis_overlapping s overlaps)

    and augment_basis_overlapping s ps =
      basis (EPowSet.fold (fun s' -> EPowSet.union (product s s')) ps EPowSet.empty)

  end

(*i

module EPowSet =
  Make (struct type t = int let compare = compare let to_string = string_of_int end)

let test lists =
  let ps = EPowSet.of_lists lists in
  let basis = EPowSet.basis ps in
  Printf.eprintf "basis %s -> %s\n" (EPowSet.to_string ps) (EPowSet.to_string basis);
  flush stderr

let _ = List.iter test
    [ [[1;3];[2;4];[3;4];[5;6]];
      [[1;2];[3;4];[5;6]];
      [[1;2;3;4];[3;4];[5;6]];
      [[1;2];[1;3;4];[1;4;5]];
      [[1;3;4];[1;3;4];[1;3;4]]
    ]

i*)

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
