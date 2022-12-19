(* circe2/filter.ml --  *)
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

(* We will assume [left.(0) = center = right.(0)] and use only [center]. *)

type t' = 
    { left' : float array;
      center' : float;
      right' : float array }

type t = 
    { left : float array;
      center : float;
      right : float array;
      norm : float array array }

let unit =
   { left = [| 1.0 |];
     center = 1.0;
     right = [| 1.0 |];
     norm = [| [| 1.0 |] |] }

let normalize f =
  let left_sum = ThoArray.sum_float ~inf:1 f.left'
  and right_sum = ThoArray.sum_float ~inf:1 f.right' in
  let norm = f.center' +. left_sum +. right_sum in
  let left = Array.map (fun x -> x /. norm) f.left'
  and center = f.center' /. norm
  and right = Array.map (fun x -> x /. norm) f.right' in
  let norm =
    Array.make_matrix (Array.length left) (Array.length right) center in
  for i = 1 to Array.length left - 1 do
    norm.(i).(0) <- norm.(pred i).(0) +. left.(i)
  done;
  for i = 0 to Array.length left - 1 do
    for j = 1 to Array.length right - 1 do
      norm.(i).(j) <- norm.(i).(pred j) +. right.(j)
    done
  done;
  { left; center; right; norm }

let upper x =
  truncate (ceil x)

let gaussian width =
  let n = upper (width *. sqrt (2. *. log 1e6)) in
  let weights =
    Array.init (succ n) (fun i -> exp (-. 0.5 *. (float i /. width) ** 2.)) in
  normalize
    { left' = weights;
      center' = 1.0;
      right' = weights }

(* Idea: avoid bleeding into empty regions by treating their
   edges like boundaries. *)

let apply ?inf ?sup f a =
  let inf = ThoArray.decode_inf ?inf a
  and sup = ThoArray.decode_sup ?sup a in
  let n_left = Array.length f.left
  and n_right = Array.length f.right
  and a' = Array.copy a in
  for i = inf to sup do
    let num_left = min (pred n_left) (i - inf)
    and num_right = min (pred n_right) (sup - i) in
    let sum = ref (f.center *. a.(i)) in
    for j = 1 to num_left do
      sum := !sum +. f.left.(j) *. a.(i-j)
    done;
    for j = 1 to num_right do
      sum := !sum +. f.right.(j) *. a.(i+j)
    done;
    a'.(i) <- !sum /. f.norm.(num_left).(num_right)
  done;
  a'

module Real = 
  struct 
    type t = float
    let compare = compare
    let compare x y =
      if abs_float (x -. y) <=
        Float.Double.epsilon *. (max (abs_float x) (abs_float y)) then
        0
      else if x < y then
        -1
      else
        1
    let pp_printer = Format.pp_print_float
    let pp_print_sep = OUnitDiff.pp_comma_separator
  end

module Reals = OUnitDiff.ListSimpleMake (Real)

let array_assert_equal a1 a2 =
  Reals.assert_equal (Array.to_list a1) (Array.to_list a2)

let limits_suite =
  let fence = Array.init 10 (fun i -> if i = 0 || i = 9 then 1.0 else 0.0) in
  let open OUnit in
  "limits" >:::
  ["1..-2" >::
   (fun () ->
     array_assert_equal fence
       (apply ~inf:1 ~sup:(-2) (gaussian 10.0) fence))]

let norm_suite =
  let flat = Array.make 10 1.0 in
  let open OUnit in
  "norm" >:::
  ["gausian 1" >::
   (fun () ->
     array_assert_equal flat (apply (gaussian 1.0) flat));
   "gausian 5" >::
   (fun () ->
     array_assert_equal flat (apply (gaussian 5.0) flat));
   "gausian 10" >::
   (fun () ->
     array_assert_equal flat (apply (gaussian 10.0) flat))]

let apply_suite =
  let open OUnit in
  "apply" >:::
  [limits_suite;
   norm_suite]

let array_map ?inf ?sup f a =
  let a' = Array.copy a in
  for i = ThoArray.decode_inf ?inf a to ThoArray.decode_sup ?sup a do
    a'.(i) <- f a.(i)
  done;
  a'

let array_map_suite =
  let five = Array.init 5 (fun i -> float (succ i)) in
  let open OUnit in
  "array_map" >:::
  ["..-2" >::
   (fun () ->
     array_assert_equal [| 2.0; 4.0; 6.0; 8.0; 5.0 |]
       (array_map ~sup:(-2) (fun x -> 2.0 *. x) five));
   "2.." >::
   (fun () ->
     array_assert_equal [| 1.0; 2.0; 6.0; 8.0; 10.0 |]
       (array_map ~inf:2 (fun x -> 2.0 *. x) five));
   "1..-2" >::
   (fun () ->
     array_assert_equal [| 1.0; 4.0; 6.0; 8.0; 5.0 |]
       (array_map ~inf:1 ~sup:(-2) (fun x -> 2.0 *. x) five))]

let apply1 ?inf1 ?sup1 ?inf2 ?sup2 f a =
  ThoMatrix.transpose
    (array_map ?inf:inf2 ?sup:sup2 
       (apply ?inf:inf1 ?sup:sup1 f)
       (ThoMatrix.transpose a))

let apply2 ?inf1 ?sup1 ?inf2 ?sup2 f a =
  array_map ?inf:inf1 ?sup:sup1
    (apply ?inf:inf2 ?sup:sup2 f) a

let apply12 ?inf1 ?sup1 ?inf2 ?sup2 f1 f2 a =
  array_map ?inf:inf1 ?sup:sup1
    (apply ?inf:inf2 ?sup:sup2 f2)
    (ThoMatrix.transpose
       (array_map ?inf:inf2 ?sup:sup2 
          (apply ?inf:inf1 ?sup:sup1 f1)
          (ThoMatrix.transpose a)))

let apply12_suite =
  let open OUnit in
  "apply12" >:::
  []

let suite =
  let open OUnit in
  "Filter" >:::
  [apply_suite;
   array_map_suite;
   apply12_suite]


(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
