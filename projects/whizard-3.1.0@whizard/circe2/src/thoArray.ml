(* circe2/thoArray.ml --  *)
(* Copyright (C) 2001-2022 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
   Circe2 is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by 
   the Free Software Foundation; either version 2, or (at your option)
   any later version.
   Circe2 is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  *)  

exception Out_of_bounds of int * int 

let decode_limit i a =
  let n = Array.length a in
  if i >= n then
    raise (Out_of_bounds (i, n))
  else if i >= 0 then
    i
  else if i >= -n then
    n + i
  else
    raise (Out_of_bounds (i, n))

let decode_inf ?inf a =
  match inf with
  | None -> 0
  | Some i -> decode_limit i a

let decode_sup ?sup a =
  match sup with
  | None -> Array.length a - 1
  | Some i -> decode_limit i a

let decode_limit_suite =
  let ten = Array.init 10 (fun i -> i) in
  let open OUnit in
  "decode_limit" >:::
  ["0" >:: (fun () -> assert_equal 0 (decode_limit 0 ten));
   "9" >:: (fun () -> assert_equal 9 (decode_limit 9 ten));
   "10" >::
   (fun () ->
     assert_raises (Out_of_bounds (10, 10))
       (fun () -> decode_limit 10 ten));
   "-1" >:: (fun () -> assert_equal 9 (decode_limit (-1) ten));
   "-10" >:: (fun () -> assert_equal 0 (decode_limit (-10) ten));
   "-11" >::
   (fun () ->
     assert_raises (Out_of_bounds (-11, 10))
       (fun () -> decode_limit (-11) ten))]


let map ?inf ?sup f a =
  let n = decode_inf ?inf a in
  Array.init (decode_sup ?sup a - n + 1) (fun i -> f a.(n+i))

let copy ?inf ?sup a =
  map ?inf ?sup (fun x -> x) a

let map_suite =
  let five = Array.init 5 succ in
  let twice n = 2 * n in
  let open OUnit in
  "map" >:::
  ["2 * .." >:: (fun () ->
    assert_equal [|2;4;6;8;10|] (map twice five));
   "2 * 1.." >:: (fun () ->
     assert_equal [|4;6;8;10|] (map twice ~inf:1 five));
   "2 * ..-2" >:: (fun () ->
     assert_equal [|2;4;6;8|] (map twice ~sup:(-2) five));
   "2 * 1..-2" >:: (fun () ->
     assert_equal [|4;6;8|] (map twice ~inf:1 ~sup:(-2) five));
   "2 * 1..2" >:: (fun () ->
     assert_equal [|4;6|] (map twice ~inf:1 ~sup:2 five))]


let copy_suite =
  let five = Array.init 5 succ in
  let open OUnit in
  "copy" >:::
  [".." >:: (fun () -> assert_equal five (copy five));
   "1.." >:: (fun () -> assert_equal [|2;3;4;5|] (copy ~inf:1 five));
   "..-2" >:: (fun () -> assert_equal [|1;2;3;4|] (copy ~sup:(-2) five));
   "1..-2" >:: (fun () -> assert_equal [|2;3;4|] (copy ~inf:1 ~sup:(-2) five));
   "1..2" >:: (fun () -> assert_equal [|2;3|] (copy ~inf:1 ~sup:2 five))]


let fold_left ?inf ?sup f x a =
  let acc = ref x in
  try
    for i = decode_inf ?inf a to decode_sup ?sup a do
      acc := f !acc a.(i)
    done;
    !acc
  with
  | Out_of_bounds (_, _) -> x

let iter ?inf ?sup f a =
  fold_left ?inf ?sup (fun () x -> f x) () a

let iter ?inf ?sup f a =
  try
    for i = decode_inf ?inf a to decode_sup ?sup a do
      f a.(i)
    done
  with
  | Out_of_bounds (_, _) -> ()

let sum_float ?inf ?sup a =
  fold_left ?inf ?sup (+.) 0.0 a

let sum_float_suite =
  let ten = Array.init 10 (fun i -> float i +. 1.0) in
  let open OUnit in
  "sum_float" >:::
  [".." >:: (fun () -> assert_equal 55.0 (sum_float ten));
   "1.." >:: (fun () -> assert_equal 54.0 (sum_float ~inf:1 ten));
   "..-2" >:: (fun () -> assert_equal 45.0 (sum_float ~sup:(-2) ten));
   "1..-2" >:: (fun () -> assert_equal 44.0 (sum_float ~inf:1 ~sup:(-2) ten));
   "1..2" >:: (fun () -> assert_equal 5.0 (sum_float ~inf:1 ~sup:2 ten))]

let suite =
  let open OUnit in
  "Array" >:::
  [decode_limit_suite;
   map_suite;
   copy_suite;
   sum_float_suite]

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
