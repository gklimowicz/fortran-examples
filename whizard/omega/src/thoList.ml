(* thoList.ml --

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

(* Avoid refering to [Pervasives.compare], because [Pervasives] will
   become [Stdlib.Pervasives] in O'Caml 4.07 and [Stdlib] in O'Caml 4.08. *)
let pcompare = compare

let rec hdn n l =
  if n <= 0 then
    []
  else
    match l with
    | x :: rest -> x :: hdn (pred n) rest
    | [] -> invalid_arg "ThoList.hdn"

let rec tln n l =
  if n <= 0 then
    l
  else
    match l with
    | _ :: rest -> tln (pred n) rest
    | [] -> invalid_arg "ThoList.tln"

let rec splitn' n l1_rev l2 =
  if n <= 0 then
    (List.rev l1_rev, l2)
  else
    match l2 with
    | x :: l2' -> splitn' (pred n) (x :: l1_rev) l2'
    | [] -> invalid_arg "ThoList.splitn n > len"

let splitn n l =
  if n < 0 then
    invalid_arg "ThoList.splitn n < 0"
  else
    splitn' n [] l

let split_last l =
  match List.rev l with
  | [] -> invalid_arg "ThoList.split_last []"
  | ln :: l12_rev -> (List.rev l12_rev, ln)

(* This is [splitn'] all over again, but without the exception. *)
let rec chopn'' n l1_rev l2 =
  if n <= 0 then
    (List.rev l1_rev, l2)
  else
    match l2 with
    | x :: l2' -> chopn'' (pred n) (x :: l1_rev) l2'
    | [] -> (List.rev l1_rev, [])
  
let rec chopn' n ll_rev = function
  | [] -> List.rev ll_rev
  | l ->
      begin match chopn'' n [] l with
      | [], [] -> List.rev ll_rev
      | l1, [] -> List.rev (l1 :: ll_rev)
      | l1, l2 -> chopn' n (l1 :: ll_rev) l2
      end

let chopn n l =
  if n <= 0 then
    invalid_arg "ThoList.chopn n <= 0"
  else
    chopn' n [] l

(* Find a member [a] in the list [l] and return the
   cyclically permuted list with [a] as head. *)
let cycle_until a l =
  let rec cycle_until' acc = function
    | [] -> raise Not_found
    | a' :: l' as al' ->
       if a' = a then
         al' @ List.rev acc
       else
         cycle_until' (a' :: acc) l' in
  cycle_until' [] l

let rec cycle' i acc l =
  if i <= 0 then
    l @ List.rev acc
  else
    match l with
    | [] -> invalid_arg "ThoList.cycle"
    | a' :: l' ->
       cycle' (pred i) (a' :: acc) l'

let cycle n l =
  if n < 0 then
    invalid_arg "ThoList.cycle"
  else
    cycle' n [] l

let of_subarray n1 n2 a =
  let rec of_subarray' n1 n2 =
    if n1 > n2 then
      []
    else
      a.(n1) :: of_subarray' (succ n1) n2 in
  of_subarray' (max 0 n1) (min n2 (pred (Array.length a)))

let range ?(stride=1) n1 n2 =
  if stride <= 0 then
    invalid_arg "ThoList.range: stride <= 0"
  else
    let rec range' n =
      if n > n2 then
        []
      else
        n :: range' (n + stride) in
    range' n1

(* Tail recursive: *)
let enumerate ?(stride=1) n l =
  let _, l_rev =
    List.fold_left
      (fun (i, acc) a -> (i + stride, (i, a) :: acc))
      (n, []) l in
  List.rev l_rev

(* Take the elements of [list] that satisfy [predicate] and
   form a list of pairs of an offset into the original list
   and the element with the offsets
   starting from [offset].  NB: the order of the returned alist
   is not specified! *)
let alist_of_list ?(predicate=(fun _ -> true)) ?(offset=0) list =
  let _, alist =
    List.fold_left
      (fun (n, acc) x ->
	(succ n, if predicate x then (n, x) :: acc else acc))
      (offset, []) list in
  alist

(* This is \emph{not} tail recursive! *)
let rec flatmap f = function
  | [] -> []
  | x :: rest -> f x @ flatmap f rest

(* This is! *)
let rev_flatmap f l =
  let rec rev_flatmap' acc f = function
    | [] -> acc
    | x :: rest -> rev_flatmap' (List.rev_append (f x) acc) f rest in
  rev_flatmap' [] f l

let rec power = function
  | [] -> [[]]
  | a :: a_list ->
     let power_a_list = power a_list in
     power_a_list @ List.map (fun a_list -> a :: a_list) power_a_list
             
let fold_left2 f acc lists =
  List.fold_left (List.fold_left f) acc lists

let fold_right2 f lists acc =
  List.fold_right (List.fold_right f) lists acc

let iteri f start list =
  ignore (List.fold_left (fun i a -> f i a; succ i) start list)

let iteri2 f start_outer star_inner lists =
  iteri (fun j -> iteri (f j) star_inner) start_outer lists

let mapi f start list =
  let next, list' =
    List.fold_left (fun (i, acc) a -> (succ i, f i a :: acc)) (start, []) list in
  List.rev list'

let rec map3 f l1 l2 l3 =
  match l1, l2, l3 with
  | [], [], [] -> []
  | a1 :: l1, a2 :: l2, a3 :: l3 ->
     let fa123 = f a1 a2 a3 in
     fa123 :: map3 f l1 l2 l3
  | _, _, _ -> invalid_arg "ThoList.map3"

(* Is there a more efficient implementation? *)
let transpose lists =
  let rec transpose' rest =
    if List.for_all ((=) []) rest then
      []
    else
      List.map List.hd rest :: transpose' (List.map List.tl rest) in
  try
    transpose' lists
  with
  | Failure s ->
     if s = "tl" then
       invalid_arg "ThoList.transpose: not rectangular"
     else
       failwith ("ThoList.transpose: unexpected Failure(" ^ s ^ ")")

let compare ?(cmp=pcompare) l1 l2 =
  let rec compare' l1' l2' =
    match l1', l2' with
    | [], [] -> 0
    | [], _ -> -1
    | _, [] -> 1
    | n1 :: r1, n2 :: r2 ->
        let c = cmp n1 n2 in
        if c <> 0 then
          c
        else
          compare' r1 r2
  in
  compare' l1 l2

let rec uniq' x = function
  | [] -> []
  | x' :: rest ->
      if x' = x then
        uniq' x rest
      else
        x' :: uniq' x' rest

let uniq = function
  | [] -> []
  | x :: rest -> x :: uniq' x rest

let rec homogeneous = function
  | [] | [_] -> true
  | a1 :: (a2 :: _ as rest) ->
      if a1 <> a2 then
        false
      else
        homogeneous rest
          
let rec pairs' acc = function
  | [] -> acc
  | [x] -> invalid_arg "pairs: odd number of elements"
  | x :: y :: indices ->
     if x <> y then
       invalid_arg "pairs: not in pairs"
     else
       begin match acc with
       | [] -> pairs' [x] indices
       | x' :: _ ->
          if x = x' then
            invalid_arg "pairs: more than twice"
          else
            pairs' (x :: acc) indices
       end

let pairs l =
  pairs' [] (List.sort pcompare l)

(* If we needed it, we could use a polymorphic version of [Set] to
   speed things up from~$O(n^2)$ to~$O(n\ln n)$.  But not before it
   matters somewhere \ldots *)
let classify l =
  let rec add_to_class a = function
    | [] -> [1, a]
    | (n, a') :: rest ->
        if a = a' then
          (succ n, a) :: rest
        else
          (n, a') :: add_to_class a rest
  in
  let rec classify' cl = function
    | [] -> cl
    | a :: rest -> classify' (add_to_class a cl) rest
  in
  classify' [] l

let rec factorize l =
  let rec add_to_class x y = function
    | [] -> [(x, [y])]
    | (x', ys) :: rest ->
        if x = x' then
          (x, y :: ys) :: rest
        else
          (x', ys) :: add_to_class x y rest
  in
  let rec factorize' fl = function
    | [] -> fl
    | (x, y) :: rest -> factorize' (add_to_class x y fl) rest
  in
  List.map (fun (x, ys) -> (x, List.rev ys)) (factorize' [] l)
    
let rec clone x n =
  if n < 0 then
    invalid_arg "ThoList.clone"
  else if n = 0 then
    []
  else
    x :: clone x (pred n)

let interleave f list =
  let rec interleave' rev_head tail =
    let rev_head' = List.rev_append (f rev_head tail) rev_head in
    match tail with
    | [] -> List.rev rev_head'
    | x :: tail' -> interleave' (x :: rev_head') tail'
  in
  interleave' [] list

let interleave_nearest f list =
  interleave
    (fun head tail ->
      match head, tail with
      | h :: _, t :: _ -> f h t
      | _ -> [])
    list

let rec rev_multiply n rl l =
  if n < 0 then
    invalid_arg "ThoList.multiply"
  else if n = 0 then
    []
  else
    List.rev_append rl (rev_multiply (pred n) rl l)

let multiply n l = rev_multiply n (List.rev l) l

let filtermap f l =
  let rec rev_filtermap acc = function
    | [] -> List.rev acc
    | a :: a_list ->
       match f a with
       | None -> rev_filtermap acc a_list
       | Some fa -> rev_filtermap (fa :: acc) a_list
  in
  rev_filtermap [] l
  
exception Overlapping_indices
exception Out_of_bounds

let iset_list_union list =
  List.fold_right Sets.Int.union list Sets.Int.empty

let complement_index_sets n index_set_lists =
  let index_sets = List.map Sets.Int.of_list index_set_lists in
  let index_set = iset_list_union index_sets in
  let size_index_sets =
    List.fold_left (fun acc s -> Sets.Int.cardinal s + acc) 0 index_sets in
  if size_index_sets <> Sets.Int.cardinal index_set then
    raise Overlapping_indices
  else if Sets.Int.exists (fun i -> i < 0 || i >= n) index_set then
    raise Overlapping_indices
  else
    match Sets.Int.elements
            (Sets.Int.diff (Sets.Int.of_list (range 0 (pred n))) index_set) with
    | [] -> index_set_lists
    | complement -> complement :: index_set_lists

let sort_section cmp array index_set =
  List.iter2
    (Array.set array)
    index_set (List.sort cmp (List.map (Array.get array) index_set))

let partitioned_sort cmp index_sets list =
  let array = Array.of_list list in
  List.fold_left
    (fun () -> sort_section cmp array)
    () (complement_index_sets (List.length list) index_sets);
  Array.to_list array

let ariadne_sort ?(cmp=pcompare) list =
  let sorted =
    List.sort (fun (n1, a1) (n2, a2) -> cmp a1 a2) (enumerate 0 list) in
  (List.map snd sorted, List.map fst sorted)

let ariadne_unsort (sorted, indices) =
  List.map snd
    (List.sort
       (fun (n1, a1) (n2, a2) -> pcompare n1 n2)
       (List.map2 (fun n a -> (n, a)) indices sorted))

let lexicographic ?(cmp=pcompare) l1 l2 =
  let rec lexicographic' = function
    | [], [] -> 0
    | [], _ -> -1
    | _, [] -> 1
    | x1 :: rest1, x2 :: rest2 ->
       let res = cmp x1 x2 in
       if res <> 0 then
	 res
       else
	 lexicographic' (rest1, rest2) in
  lexicographic' (l1, l2)

(* If there was a polymorphic [Set], we could also say
   [Set.elements (Set.union (Set.of_list l1) (Set.of_list l2))]. *)
let common l1 l2 =
  List.fold_left
    (fun acc x1 ->
      if List.mem x1 l2 then
	x1 :: acc
      else
	acc)
    [] l1

let complement l1 = function
  | [] -> l1
  | l2 ->
     if List.for_all (fun x -> List.mem x l1) l2 then
       List.filter (fun x -> not (List.mem x l2)) l1
     else
       invalid_arg "ThoList.complement"


let to_string a2s alist =
  "[" ^ String.concat "; " (List.map a2s alist) ^ "]"

let random_int_list imax n =
  let imax_plus = succ imax in
  Array.to_list (Array.init n (fun _ -> Random.int imax_plus))

module Test =
  struct

    let int_list2_to_string l2 =
      to_string (to_string string_of_int) l2

    (* Inefficient, must only be used for unit tests. *)
    let compare_lists_by_size l1 l2 =
      let lengths = pcompare (List.length l1) (List.length l2) in
      if lengths = 0 then
        pcompare l1 l2
      else
        lengths

    open OUnit

    let suite_filtermap =
      "filtermap" >:::
        [ "filtermap Some []" >::
            (fun () ->
              assert_equal ~printer:(to_string string_of_int)
                [] (filtermap (fun x -> Some x) []));

          "filtermap None []" >::
            (fun () ->
              assert_equal ~printer:(to_string string_of_int)
                [] (filtermap (fun x -> None) []));

          "filtermap even_neg []" >::
            (fun () ->
              assert_equal ~printer:(to_string string_of_int)
                [0; -2; -4]
                (filtermap
                   (fun n -> if n mod 2 = 0 then Some (-n) else None)
                   (range 0 5)));

          "filtermap odd_neg []" >::
            (fun () ->
              assert_equal ~printer:(to_string string_of_int)
                [-1; -3; -5]
                (filtermap
                   (fun n -> if n mod 2 <> 0 then Some (-n) else None)
                   (range 0 5))) ]
          
    let assert_power power_a_list a_list =
      assert_equal ~printer:int_list2_to_string
        power_a_list
        (List.sort compare_lists_by_size (power a_list))

    let suite_power =
      "power" >:::
        [ "power []" >::
            (fun () ->
              assert_power [[]] []);

          "power [1]" >::
            (fun () ->
              assert_power [[]; [1]] [1]);

          "power [1;2]" >::
            (fun () ->
              assert_power [[]; [1]; [2]; [1;2]] [1;2]);

          "power [1;2;3]" >::
            (fun () ->
              assert_power
                [[];
                 [1]; [2]; [3];
                 [1;2]; [1;3]; [2;3];
                 [1;2;3]]
                [1;2;3]);

          "power [1;2;3;4]" >::
            (fun () ->
              assert_power
                [[];
                 [1]; [2]; [3]; [4];
                 [1;2]; [1;3]; [1;4]; [2;3]; [2;4]; [3;4];
                 [1;2;3]; [1;2;4]; [1;3;4]; [2;3;4];
                 [1;2;3;4]]
                [1;2;3;4]) ]

    let suite_split =
      "split*" >:::
	[ "split_last []" >::
	    (fun () ->
	      assert_raises
                (Invalid_argument "ThoList.split_last []")
                (fun () -> split_last []));
          "split_last [1]" >::
	    (fun () ->
	      assert_equal
                ([], 1)
                (split_last [1]));
          "split_last [2;3;1;4]" >::
	    (fun () ->
	      assert_equal
                ([2;3;1], 4)
                (split_last [2;3;1;4])) ]

    let test_list = random_int_list 1000 100

    let assert_equal_int_list =
      assert_equal ~printer:(to_string string_of_int)

    let suite_cycle =
      "cycle_until" >:::
	[ "cycle (-1) [1;2;3]" >::
	    (fun () ->
	      assert_raises
                (Invalid_argument "ThoList.cycle")
                (fun () -> cycle 4 [1;2;3]));
          "cycle 4 [1;2;3]" >::
	    (fun () ->
	      assert_raises
                (Invalid_argument "ThoList.cycle")
                (fun () -> cycle 4 [1;2;3]));
          "cycle 42 [...]" >::
	    (fun () ->
              let n = 42 in
	      assert_equal_int_list
                (tln n test_list @ hdn n test_list)
                (cycle n test_list));
          "cycle_until 1 []" >::
	    (fun () ->
	      assert_raises
                (Not_found)
                (fun () -> cycle_until 1 []));
          "cycle_until 1 [2;3;4]" >::
	    (fun () ->
	      assert_raises
                (Not_found)
                (fun () -> cycle_until 1 [2;3;4]));
          "cycle_until 1 [1;2;3;4]" >::
	    (fun () ->
	      assert_equal
                [1;2;3;4]
                (cycle_until 1 [1;2;3;4]));
          "cycle_until 3 [1;2;3;4]" >::
	    (fun () ->
	      assert_equal
                [3;4;1;2]
                (cycle_until 3 [3;4;1;2]));
          "cycle_until 4 [1;2;3;4]" >::
	    (fun () ->
	      assert_equal
                [4;1;2;3]
                (cycle_until 4 [4;1;2;3])) ]

    let suite_alist_of_list =
      "alist_of_list" >:::
	[ "simple" >::
	    (fun () ->
	      assert_equal
                [(46, 4); (44, 2); (42, 0)]
                (alist_of_list
                   ~predicate:(fun n -> n mod 2 = 0) ~offset:42 [0;1;2;3;4;5])) ]

    let suite_complement =
      "complement" >:::
	[ "simple" >::
	    (fun () ->
	      assert_equal [2;4] (complement [1;2;3;4] [1; 3]));
          "empty" >::
	    (fun () ->
	      assert_equal [1;2;3;4] (complement [1;2;3;4] []));
          "failure" >::
	    (fun () ->
              assert_raises
                (Invalid_argument ("ThoList.complement"))
	        (fun () -> complement (complement [1;2;3;4] [5]))) ]

    let suite =
      "ThoList" >:::
	[suite_filtermap;
         suite_power;
         suite_split;
         suite_cycle;
         suite_alist_of_list;
         suite_complement]

  end
