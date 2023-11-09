(* thoArray.ml --

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

type 'a compressed = 
    { uniq : 'a array;
      embedding: int array }

let uniq a = a.uniq
let embedding a = a.embedding

type 'a compressed2 = 
    { uniq2 : 'a array array;
      embedding1: int array;
      embedding2: int array }

let uniq2 a = a.uniq2
let embedding1 a = a.embedding1
let embedding2 a = a.embedding2

module PMap = Pmap.Tree

let compress a =
  let last = Array.length a - 1 in
  let embedding = Array.make (succ last) (-1) in
  let rec scan num_uniq uniq elements n =
    if n > last then
      { uniq = Array.of_list (List.rev elements);
        embedding = embedding }
    else
      match PMap.find_opt compare a.(n) uniq with
      | Some n' ->
          embedding.(n) <- n';
          scan num_uniq uniq elements (succ n)
      | None ->
          embedding.(n) <- num_uniq;
          scan
            (succ num_uniq)
            (PMap.add compare a.(n) num_uniq uniq)
            (a.(n) :: elements)
            (succ n) in
  scan 0 PMap.empty [] 0
  
let uncompress a =
  Array.map (Array.get a.uniq) a.embedding

(* \begin{dubious}
     Using [transpose] simplifies the algorithms, but can be inefficient.
     If this turns out to be the case, we should add special treatments
     for symmetric matrices.
   \end{dubious} *)

let transpose a =
  let dim1 = Array.length a
  and dim2 = Array.length a.(0) in
  let a' = Array.make_matrix dim2 dim1 a.(0).(0) in
  for i1 = 0 to pred dim1 do
    for i2 = 0 to pred dim2 do
      a'.(i2).(i1) <- a.(i1).(i2)
    done
  done;
  a'

let compress2 a =
  let c2 = compress a in
  let c12_transposed = compress (transpose c2.uniq) in
  { uniq2 = transpose c12_transposed.uniq;
    embedding1 = c12_transposed.embedding;
    embedding2 = c2.embedding }

let uncompress2 a =
  let a2 = uncompress { uniq = a.uniq2; embedding = a.embedding2 } in
  transpose (uncompress { uniq = transpose a2; embedding = a.embedding1 })

(* FIXME: not tail recursive! *)
let compare ?(cmp=pcompare) a1 a2 =
  let l1 = Array.length a1
  and l2 = Array.length a2 in
  if l1 < l2 then
    -1
  else if l1 > l2 then
    1
  else
    let rec scan i =
      if i = l1 then
        0
      else
        let c = cmp a1.(i) a2.(i) in
        if c < 0 then
          -1
        else if c > 0 then
          1
        else
          scan (succ i) in
    scan 0

let find_first f a =
  let l = Array.length a in
  let rec find_first' i =
    if i >= l then
      raise Not_found
    else if f (a.(i)) then
      i
    else
      find_first' (succ i)
  in
  find_first' 0

let match_first x a =
  find_first (fun x' -> x = x') a

let find_all f a =
  let matches = ref [] in
  for i = Array.length a - 1 downto 0 do
    if f (a.(i)) then
      matches := i :: !matches
  done;
  !matches

let match_all x a =
  find_all (fun x' -> x = x') a

let num_rows a =
  Array.length a

let num_columns a =
  match ThoList.classify (List.map Array.length (Array.to_list a)) with
  | [ (_, n) ] -> n
  | _ -> invalid_arg "ThoArray.num_columns: inhomogeneous array"

module Test =
  struct

    open OUnit

    let test_compare_empty =
      "empty" >::
	(fun () -> assert_equal   0  (compare [| |] [| |]))
        
    let test_compare_shorter =
      "shorter" >::
	(fun () -> assert_equal (-1) (compare [|0|] [|0; 1|]))
        
    let test_compare_longer =
      "longer" >::
	(fun () -> assert_equal ( 1) (compare [|0; 1|] [|0|]))
        
    let test_compare_less =
      "longer" >::
	(fun () -> assert_equal (-1) (compare [|0; 1|] [|0; 2|]))
        
    let test_compare_equal =
      "equal" >::
	(fun () -> assert_equal ( 0) (compare [|0; 1|] [|0; 1|]))
        
    let test_compare_more =
      "more" >::
	(fun () -> assert_equal ( 1) (compare [|0; 2|] [|0; 1|]))
        
    let suite_compare =
      "compare" >:::
        [test_compare_empty;
         test_compare_shorter;
         test_compare_longer;
         test_compare_less;
         test_compare_equal;
         test_compare_more]

    let test_find_first_not_found =
      "not found" >::
	(fun () ->
	  assert_raises Not_found
            (fun () -> find_first (fun n -> n mod 2 = 0) [|1;3;5|]))
        
    let test_find_first_first =
      "first" >::
	(fun () ->
	  assert_equal 0
            (find_first (fun n -> n mod 2 = 0) [|2;3;4;5|]))
        
    let test_find_first_not_last =
      "last" >::
	(fun () ->
	  assert_equal 1
            (find_first (fun n -> n mod 2 = 0) [|1;2;3;4|]))
        
    let test_find_first_last =
      "not last" >::
	(fun () ->
	  assert_equal 1
            (find_first (fun n -> n mod 2 = 0) [|1;2|]))
        
    let suite_find_first =
      "find_first" >:::
	[test_find_first_not_found;
         test_find_first_first;
         test_find_first_not_last;
         test_find_first_last]

    let test_find_all_empty =
      "empty" >::
	(fun () ->
	  assert_equal []
            (find_all (fun n -> n mod 2 = 0) [|1;3;5|]))
        
    let test_find_all_first =
      "first" >::
	(fun () ->
	  assert_equal [0;2]
            (find_all (fun n -> n mod 2 = 0) [|2;3;4;5|]))
        
    let test_find_all_not_last =
      "last" >::
	(fun () ->
	  assert_equal [1;3]
            (find_all (fun n -> n mod 2 = 0) [|1;2;3;4;5|]))
        
    let test_find_all_last =
      "not last" >::
	(fun () ->
	  assert_equal [1;3]
            (find_all (fun n -> n mod 2 = 0) [|1;2;3;4|]))
        
    let suite_find_all =
      "find_all" >:::
	[test_find_all_empty;
         test_find_all_first;
         test_find_all_last;
         test_find_all_not_last]

    let test_num_columns_ok2 =
      "ok/2" >::
	(fun () ->
	  assert_equal 2
            (num_columns [| [| 11; 12 |];
                            [| 21; 22 |];
                            [| 31; 32 |] |]))

    let test_num_columns_ok0 =
      "ok/0" >::
	(fun () ->
	  assert_equal 0
            (num_columns [| [| |];
                            [| |];
                            [| |] |]))

    let test_num_columns_not_ok =
      "not_ok" >::
	(fun () ->
	  assert_raises (Invalid_argument
                           "ThoArray.num_columns: inhomogeneous array")
            (fun () -> num_columns [| [| 11; 12 |];
                                      [| 21 |];
                                      [| 31; 32 |] |]))

    let suite_num_columns =
      "num_columns" >:::
	[test_num_columns_ok2;
         test_num_columns_ok0;
         test_num_columns_not_ok]

    let suite =
      "ThoArrays" >:::
	[suite_compare;
         suite_find_first;
         suite_find_all;
         suite_num_columns]

  end

(*i
 *  Local Variables:
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)





