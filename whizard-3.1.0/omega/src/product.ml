(* product.ml --

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

(* \thocwmodulesection{Lists} *)

(* We use the tail recursive [List.fold_left] over [List.fold_right]
   for efficiency, but revert the argument lists in order to preserve
   lexicographic ordering.   The argument lists are much shorter than
   the results, so the cost of the [List.rev] is negligible. *)

let fold2_rev f l1 l2 acc =
  List.fold_left (fun acc1 x1 ->
    List.fold_left (fun acc2 x2 -> f x1 x2 acc2) acc1 l2) acc l1

let fold2 f l1 l2 acc =
  fold2_rev f (List.rev l1) (List.rev l2) acc

let fold3_rev f l1 l2 l3 acc =
  List.fold_left (fun acc1 x1 -> fold2 (f x1) l2 l3 acc1) acc l1

let fold3 f l1 l2 l3 acc =
  fold3_rev f (List.rev l1) (List.rev l2) (List.rev l3) acc

(* If all lists have the same type, there's also *)

let rec fold_rev f ll acc =
  match ll with
  | [] -> acc
  | [l] ->  List.fold_left (fun acc' x -> f [x] acc') acc l
  | l :: rest ->
      List.fold_left (fun acc' x -> fold_rev (fun xr -> f (x::xr)) rest acc') acc l

let fold f ll acc = fold_rev f (List.map List.rev ll) acc

let list2 op l1 l2 =
  fold2 (fun x1 x2 c -> op x1 x2 :: c) l1 l2 []

let list3 op l1 l2 l3 =
  fold3 (fun x1 x2 x3 c -> op x1 x2 x3 :: c) l1 l2 l3 []

let list op ll =
  fold (fun l c -> op l :: c) ll []

let list2_opt op l1 l2 =
  fold2
    (fun x1 x2 c ->
      match op x1 x2 with
      | None -> c
      | Some op_x1_x2 -> op_x1_x2 :: c)
    l1 l2 []

let list3_opt op l1 l2 l3 =
  fold3
    (fun x1 x2 x3 c ->
      match op x1 x2 x3 with
      | None -> c
      | Some op_x1_x2_x3 -> op_x1_x2_x3 :: c)
    l1 l2 l3 []

let list_opt op ll =
  fold
    (fun l c ->
      match op l with
      | None -> c
      | Some op_l -> op_l :: c)
    ll []

let power n l =
  list (fun x -> x) (ThoList.clone l n)

(* Reshuffling lists: 
   \begin{equation}
    \lbrack 
       \lbrack a_1;\ldots;a_k \rbrack;
       \lbrack b_1;\ldots;b_k \rbrack;
       \lbrack c_1;\ldots;c_k \rbrack; 
    \ldots\rbrack \rightarrow  
    \lbrack
       \lbrack a_1;b_1;c_1;\ldots\rbrack;
       \lbrack a_2;b_2;c_2;\ldots\rbrack; 
     \ldots\rbrack
   \end{equation}  
*)

(*i JR/WK
let thread l =
  List.map List.rev 
    (List.fold_left (fun i acc -> List.map2 (fun a b -> b::a) i acc)
       (List.map (fun i -> [i]) (List.hd l)) (List.tl l))
i*)

(* \begin{dubious}
     [tho:] Is this really an optimal implementation?
   \end{dubious} *)

let thread = function
  | head :: tail ->
      List.map List.rev 
        (List.fold_left (fun i acc -> List.map2 (fun a b -> b::a) i acc)
           (List.map (fun i -> [i]) head) tail)
  | [] -> []

(* \thocwmodulesection{Sets} *)

(* The implementation is amazingly simple: *)

type 'a set

type ('a, 'a_set, 'b) fold = ('a -> 'b -> 'b) -> 'a_set -> 'b -> 'b
type ('a, 'a_set, 'b, 'b_set, 'c) fold2 =
    ('a -> 'b -> 'c -> 'c) -> 'a_set -> 'b_set -> 'c -> 'c

let outer fold1 fold2 f l1 l2 = fold1 (fun x1 -> fold2 (f x1) l2) l1
let outer_self fold f l1 l2 = fold (fun x1 -> fold (f x1) l2) l1

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)





