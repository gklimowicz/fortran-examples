(* young.ml --

   Copyright (C) 2022- by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>

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

(* Avoid refering to [Pervasives.compare], because [Pervasives] will
   become [Stdlib.Pervasives] in O'Caml 4.07 and [Stdlib] in O'Caml 4.08. *)
let pcompare = compare

type diagram = int list
type 'a tableau = 'a list list

(* Not exposed.  Just for documentation. *)
type 'a table = 'a option array array

(* The following three are candidates for [ThoList]. *)
let rec sum = function
  | [] -> 0
  | n :: rest -> n + sum rest

let rec product = function
  | [] -> 1
  | n :: rest -> n * product rest

(* Test a predicate for each pair of consecutive elements of a list.
   Trivially true for empty and one-element lists. *)
let rec for_all_pairs predicate = function
  | [] | [_] -> true
  | a1 :: (a2 :: _ as a_list) ->
     if not (predicate a1 a2) then
       false
     else
       for_all_pairs predicate a_list

let decreasing l = for_all_pairs (fun a1 a2 -> pcompare a1 a2 > 0) l
let increasing l = for_all_pairs (fun a1 a2 -> pcompare a1 a2 < 0) l
let non_increasing l = for_all_pairs (fun a1 a2 -> pcompare a1 a2 >= 0) l
let non_decreasing l = for_all_pairs (fun a1 a2 -> pcompare a1 a2 <= 0) l

let valid_diagram = non_increasing

let diagram_rows d =
  List.length d

let diagram_columns = function
  | [] -> 0
  | nc :: _ -> nc

let take_column d =
  let rec take_column' len acc = function
    | [] -> (len, List.rev acc)
    | cols :: rest ->
       if cols <= 1 then
         take_column' (succ len) acc rest
       else
         take_column' (succ len) (pred cols :: acc) rest in
  take_column' 0 [] d

let conjugate_diagram_new d =
  let rec conjugate_diagram' rows =
    match take_column rows with
    | n, [] -> [n]
    | n, rest -> n :: conjugate_diagram' rest in
  conjugate_diagram' d

let tableau_rows t =
  List.length t

let tableau_columns = function
  | [] -> 0
  | row :: _ -> List.length row

let num_cells_diagram d =
  sum d

let cells_tableau t =
  List.flatten t

let num_cells_tableau t =
  List.fold_left (fun acc row -> acc + List.length row) 0 t

let diagram_of_tableau t =
  List.map List.length t

let tableau_of_diagram cell d =
  List.map (ThoList.clone cell) d

(* Note that the first index counts the rows and the second the columns! *)
let array_of_tableau t =
  let nr = tableau_rows t
  and nc = tableau_columns t in
  let a = Array.make_matrix nr nc None in
  List.iteri
    (fun ir -> List.iteri (fun ic cell -> a.(ir).(ic) <- Some cell))
    t;
  a

let transpose_array a =
  let nr = Array.length a in
  if nr <= 0 then
    invalid_arg "Young.transpose_array"
  else
    let nc = Array.length a.(0) in
    let a' = Array.make_matrix nc nr None in
    for ic = 0 to pred nc do
      for ir = 0 to pred nr do
        a'.(ic).(ir) <- a.(ir).(ic)
      done
    done;
    a'
         
let list_of_array_row a =
  let n = Array.length a in
  let rec list_of_array_row' ic =
    if ic >= n then
      []
    else
      match a.(ic) with
      | None -> []
      | Some cell -> cell :: list_of_array_row' (succ ic) in
  list_of_array_row' 0

let tableau_of_array a =
  Array.fold_right (fun row acc -> list_of_array_row row :: acc) a []

let conjugate_tableau t =
  array_of_tableau t |> transpose_array |> tableau_of_array

let conjugate_diagram d =
  tableau_of_diagram () d |> conjugate_tableau |> diagram_of_tableau

let valid_tableau t =
  valid_diagram (diagram_of_tableau t)

let semistandard_tableau t =
  let rows = t
  and columns = conjugate_tableau t in
  valid_tableau t
  && List.for_all non_decreasing rows
  && List.for_all increasing columns

let standard_tableau ?offset t =
  match List.sort pcompare (cells_tableau t) with
  | [] -> true
  | cell :: _ as cell_list ->
     (match offset with None -> true | Some o -> cell = o)
     && for_all_pairs (fun c1 c2 -> c2 = c1 + 1) cell_list
     && semistandard_tableau t

let hook_lengths_table d =
  let nr = diagram_rows d
  and nc = diagram_columns d in
  if min nr nc <= 0 then
    invalid_arg "Young.hook_lengths_table"
  else
    let a = array_of_tableau (tableau_of_diagram 0 d) in
    let cols = Array.of_list d
    and rows = transpose_array a |> tableau_of_array
               |> diagram_of_tableau |> Array.of_list in
    for ir = 0 to pred nr do
      for ic = 0 to pred cols.(ir) do
        a.(ir).(ic) <- Some (rows.(ic) - ir + cols.(ir) - ic - 1)
      done
    done;
    a

(* \begin{dubious}
     The following products and factorials can easily overflow,
     even if the final ratio is a smallish number.  We can avoid
     this by representing them as lists of factors (or maps from
     factors to powers).  The ratio can be computed by first
     cancelling all common factors and multiplying the remaining
     factors at the very end.
   \end{dubious} *)

let hook_lengths_product d =
  let nr = diagram_rows d
  and nc = diagram_columns d in
  if min nr nc <= 0 then
    0
  else
    let cols = Array.of_list d
    and rows = Array.of_list (conjugate_diagram d) in
    let n = ref 1 in
    for ir = 0 to pred nr do
      for ic = 0 to pred cols.(ir) do
        n := !n * (rows.(ic) - ir + cols.(ir) - ic - 1)
      done
    done;
    !n

let num_standard_tableaux d =
  let num = Combinatorics.factorial (num_cells_diagram d)
  and den = hook_lengths_product d in
  if num mod den <> 0 then
    failwith "Young.num_standard_tableaux"
  else
    num / den

(* Note that [hook_lengths_product] calls [conjugate_diagram]
   and this calls it again.
   This is wasteful, but probably no big deal for our applications. *)
let normalization d =
  let num =
    product (List.map Combinatorics.factorial (d @ conjugate_diagram d))
  and den = hook_lengths_product d in
  (num, den)

module type Test =
  sig
    val suite : OUnit.test
    val suite_long : OUnit.test
  end

module Test =
  struct
    open OUnit

    let random_int ratio =
      truncate (Random.float ratio +. 0.5)

    let random_diagram ?(ratio=1.0) rows =
      let rec random_diagram' acc row cols =
        if row >= rows then
          acc
        else
          let cols' = cols + random_int ratio in
          random_diagram' (cols' :: acc) (succ row) cols' in
      random_diagram' [] 0 (1 + random_int ratio)

    let suite_hook_lengths_product =
      "hook_lengths_product" >:::

        [ "[4;3;2]" >::
	    (fun () -> assert_equal 2160 (hook_lengths_product [4; 3; 2])) ]

    let suite_num_standard_tableaux =
      "num_standard_tableaux" >:::

        [ "[4;3;2]" >::
	    (fun () -> assert_equal 168 (num_standard_tableaux [4; 3; 2])) ]

    let suite_normalization =
      "normalization" >:::

        [ "[2;1]" >::
	    (fun () -> assert_equal (4, 3) (normalization [2; 1])) ]

    let suite =
      "Young" >:::
	[suite_hook_lengths_product;
         suite_num_standard_tableaux;
         suite_normalization]

    let suite_long =
      "Young long" >:::
	[]

  end
