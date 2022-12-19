(* permutation.ml --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       cf. main AUTHORS file

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
    type t
    val of_list : int list -> t
    val of_array : int array -> t
    val of_lists : 'a list -> 'a list -> t
    val inverse : t -> t
    val compose : t -> t -> t
    val compose_inv : t -> t -> t
    val list : t -> 'a list -> 'a list
    val array : t -> 'a array -> 'a array
    val all : int -> t list
    val even : int -> t list
    val odd : int -> t list
    val cyclic : int -> t list
    val signed : int -> (int * t) list
    val to_string : t -> string
  end

let same_elements l1 l2 =
  List.sort compare l1 = List.sort compare l2

module PM = Pmap.Tree

let offset_map l =
  let _, offsets =
    List.fold_left
      (fun (i, map) a -> (succ i, PM.add compare a i map))
      (0, PM.empty) l in
  offsets

(* TODO: this algorithm fails if the lists contain duplicate elements. *)
let of_lists_list l l' =
  if same_elements l l' then
    let offsets' = offset_map l' in
    let _, p_rev =
      List.fold_left
        (fun (i, acc) a -> (succ i, PM.find compare a offsets' :: acc))
        (0, []) l in
    List.rev p_rev
  else
    invalid_arg "Permutation.of_lists: incompatible lists"

module Using_Lists : T =
  struct

    type t = int list

    let of_list p =
      if List.sort compare p <> (ThoList.range 0 (List.length p - 1)) then
	invalid_arg "Permutation.of_list"
      else
	p

    let of_array p =
      try
	of_list (Array.to_list p)
      with 
      | Invalid_argument s ->
         if s = "Permutation.of_list" then
	   invalid_arg "Permutation.of_array"
         else
           failwith ("Permutation.of_array: unexpected Invalid_argument(" ^
                       s ^ ")")

    let of_lists = of_lists_list

    let inverse p = snd (ThoList.ariadne_sort p)

    let list p l =
      List.map snd
	(List.sort (fun (i, _) (j, _) -> compare i j)
	   (try
	      List.rev_map2 (fun i x -> (i, x)) p l
	    with
	    | Invalid_argument s ->
               if s = "List.rev_map2" then
	         invalid_arg "Permutation.list: length mismatch"
               else
                 failwith ("Permutation.list: unexpected Invalid_argument(" ^
                             s ^ ")")))

    let array p a =
      try
	Array.of_list (list p (Array.to_list a))
      with 
      | Invalid_argument s ->
         if s = "Permutation.list: length mismatch" then
	   invalid_arg "Permutation.array: length mismatch"
         else
           failwith ("Permutation.array: unexpected Invalid_argument(" ^ s ^ ")")

    let compose_inv p q =
      list q p

(* Probably not optimal (or really inefficient), but correct by
   associativity. *)

    let compose p q =
      list (inverse q) p

    let all n =
      List.map of_list (Combinatorics.permute (ThoList.range 0 (pred n)))

    let even n =
      List.map of_list (Combinatorics.permute_even (ThoList.range 0 (pred n)))

    let odd n =
      List.map of_list (Combinatorics.permute_odd (ThoList.range 0 (pred n)))

    let cyclic n =
      List.map of_list (Combinatorics.permute_cyclic (ThoList.range 0 (pred n)))

    let signed n =
      List.map
        (fun (eps, l) -> (eps, of_list l))
        (Combinatorics.permute_signed (ThoList.range 0 (pred n)))

    let to_string p =
      String.concat "" (List.map string_of_int p)

  end

module Using_Arrays : T =
  struct

    type t = int array

    let of_list p =
      if List.sort compare p <> (ThoList.range 0 (List.length p - 1)) then
	invalid_arg "Permutation.of_list"
      else
	Array.of_list p

    let of_array p =
      try
	of_list (Array.to_list p)
      with 
      | Invalid_argument s ->
         if s = "Permutation.of_list" then
	   invalid_arg "Permutation.of_array"
         else
           failwith ("Permutation.of_array: unexpected Invalid_argument(" ^
                       s ^ ")")

    let of_lists l l' =
      Array.of_list (of_lists_list l l')

    let inverse p =
      let len_p = Array.length p in
      let p' = Array.make len_p p.(0) in
      for i = 0 to pred len_p do
	p'.(p.(i)) <- i
      done;
      p'

    let array p a =
      let len_a = Array.length a
      and len_p = Array.length p in
      if len_a <> len_p then
	invalid_arg "Permutation.array: length mismatch";
      let a' = Array.make len_a a.(0) in
      for i = 0 to pred len_a do
	a'.(p.(i)) <- a.(i)
      done;
      a'

    let list p l =
      try
	Array.to_list (array p (Array.of_list l))
      with 
      | Invalid_argument s ->
         if s = "Permutation.array: length mismatch" then
	   invalid_arg "Permutation.list: length mismatch"
         else
           failwith ("Permutation.list: unexpected Invalid_argument(" ^ s ^ ")")

    let compose_inv p q =
      array q p

    let compose p q =
      array (inverse q) p

    let all n =
      List.map of_list (Combinatorics.permute (ThoList.range 0 (pred n)))

    let even n =
      List.map of_list (Combinatorics.permute_even (ThoList.range 0 (pred n)))

    let odd n =
      List.map of_list (Combinatorics.permute_odd (ThoList.range 0 (pred n)))

    let cyclic n =
      List.map of_list (Combinatorics.permute_cyclic (ThoList.range 0 (pred n)))

    let signed n =
      List.map
        (fun (eps, l) -> (eps, of_list l))
        (Combinatorics.permute_signed (ThoList.range 0 (pred n)))

    let to_string p =
      String.concat "" (List.map string_of_int (Array.to_list p))

  end

module Default = Using_Arrays

(*
  This is the Fisher-Yates shuffle, cf. D. Knuth, {\em Seminumerical
  algorithms.  The Art of Computer Programming. 2}. Reading, MA:
  Addison–Wesley. pp. 139-140.
 *)

(*i
  To shuffle an array a of n elements (indices 0..n-1):

     for i from n − 1 downto 1 do
          j ← random integer with 0 ≤ j ≤ i
          exchange a[j] and a[i]

   To initialize an array a of n elements to a randomly shuffled copy
   of source, both 0-based: 

     a[0] ← source[0]
     for i from 1 to n − 1 do
         j ← random integer with 0 ≤ j ≤ i
         a[i] ← a[j]
         a[j] ← source[i]
i*)

let shuffle l =
  let a = Array.of_list l in
  for n = Array.length a - 1 downto 1 do
    let k = Random.int (succ n) in
    if k <> n then
      let tmp  = Array.get a n in
      Array.set a n (Array.get a k);
      Array.set a k tmp
  done;
  Array.to_list a

let time f x =
  let start = Sys.time () in
  let f_x = f x in
  let stop = Sys.time () in
  (f_x, stop -. start)
  
let print_time msg f x =
  let f_x, seconds = time f x in
  Printf.printf "%s took %10.2f ms\n" msg (seconds *. 1000.);
  f_x
  
let random_int_list imax n =
  let imax_plus = succ imax in
  Array.to_list (Array.init n (fun _ -> Random.int imax_plus))

module Test (P : T) : sig val suite : OUnit.test val time : unit -> unit end =
  struct

    open OUnit
    open P

    let of_list_overlap =
      "overlap" >::
	(fun () ->
	  assert_raises (Invalid_argument "Permutation.of_list")
	    (fun () ->
	      of_list [0;1;2;2]))
	
    let of_list_gap =
      "gap" >::
	(fun () ->
	  assert_raises (Invalid_argument "Permutation.of_list")
	    (fun () ->
	      of_list [0;1;2;4;5]))

    let of_list_ok =
      "ok" >::
	(fun () ->
	  let l = ThoList.range 0 10 in
	  assert_equal (of_list l) (of_list l))

    let suite_of_list =
      "of_list" >:::
	[of_list_overlap;
	 of_list_gap;
	 of_list_ok]

    let suite_of_lists =
      "of_lists" >:::
	[ "ok" >::
	    (fun () ->
              for i = 1 to 10 do
	        let l = random_int_list 1000000 100 in
                let l' = shuffle l in
	        assert_equal
                  ~printer:(ThoList.to_string string_of_int)
                  l' (list (of_lists l l') l)
              done) ]

    let apply_invalid_lengths =
      "invalid/lengths" >::
	(fun () ->
	  assert_raises
	    (Invalid_argument "Permutation.list: length mismatch")
	    (fun () ->
	      list (of_list [0;1;2;3;4]) [0;1;2;3]))

    let apply_ok =
      "ok" >::
	(fun () ->
	  assert_equal [2;0;1;3;5;4]
	    (list (of_list [1;2;0;3;5;4]) [0;1;2;3;4;5]))

    let suite_apply =
      "apply" >:::
	[apply_invalid_lengths;
	 apply_ok]

    let inverse_ok =
      "ok" >::
	(fun () ->
	  let l = shuffle (ThoList.range 0 1000) in
	  let p = of_list (shuffle l) in
	  assert_equal l (list (inverse p) (list p l)))

    let suite_inverse =
      "inverse" >:::
	[inverse_ok]

    let compose_ok =
      "ok" >::
	(fun () ->
	  let id = ThoList.range 0 1000 in
	  let p = of_list (shuffle id)
	  and q = of_list (shuffle id)
	  and l = id in
	  assert_equal (list p (list q l)) (list (compose p q) l))
		
    let compose_inverse_ok =
      "inverse/ok" >::
	(fun () ->
	  let id = ThoList.range 0 1000 in
	  let p = of_list (shuffle id)
	  and q = of_list (shuffle id) in
	  assert_equal
	    (compose (inverse p) (inverse q))
	    (inverse (compose q p)))
		
    let suite_compose =
      "compose" >:::
	[compose_ok;
	 compose_inverse_ok]

    let suite =
      "Permutations" >:::
	[suite_of_list;
	 suite_of_lists;
	 suite_apply;
	 suite_inverse;
	 suite_compose]

    let repeat repetitions size =
      let id = ThoList.range 0 size in
      let p = of_list (shuffle id)
      and l = shuffle (List.map string_of_int id) in
      print_time (Printf.sprintf "reps=%d, len=%d" repetitions size)
	(fun () ->
	  for i = 1 to repetitions do
	    ignore (P.list p l)
	  done)
	()
      
    let time () =
      repeat 100000 10;
      repeat 10000 100;
      repeat 1000 1000;
      repeat 100 10000;
      repeat 10 100000;
      ()

  end

