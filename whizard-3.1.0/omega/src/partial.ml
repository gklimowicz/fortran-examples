(* partial.ml --

   Copyright (C) 1999-2015 by

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

module type T =
  sig
    type domain
    type 'a t
    val of_list : (domain * 'a) list -> 'a t
    val of_lists : domain list -> 'a list -> 'a t
    exception Undefined of domain
    val apply : 'a t -> domain -> 'a
    val apply_with_fallback : (domain -> 'a) -> 'a t -> domain -> 'a
    val auto : domain t -> domain -> domain
  end

module Make (D : Map.OrderedType) : T with type domain = D.t =
  struct

    module M = Map.Make (D)

    type domain = D.t
    type 'a t = 'a M.t

    let of_list l =
      List.fold_left (fun m (d, v) -> M.add d v m) M.empty l

    let of_lists domain values =
      of_list
	(try
	   List.map2 (fun d v -> (d, v)) domain values
	 with
	 | Invalid_argument _ (* ["List.map2"] *) ->
	    invalid_arg "Partial.of_lists: length mismatch")

    let auto partial d =
      try
	M.find d partial
      with
      | Not_found -> d

    exception Undefined of domain

    let apply partial d =
      try
	M.find d partial
      with
      | Not_found -> raise (Undefined d)

    let apply_with_fallback fallback partial d =
      try
	M.find d partial
      with
      | Not_found -> fallback d

  end

(* \thocwmodulesection{Unit Tests} *)

module Test : sig val suite : OUnit.test end =
  struct

    open OUnit

    module P = Make (struct type t = int let compare = compare end)

    let apply_ok =
      "apply/ok" >::
	(fun () ->
	  let p = P.of_list [ (0,"a"); (1,"b"); (2,"c") ]
	  and l = [ 0; 1; 2 ] in
	  assert_equal [ "a"; "b"; "c" ] (List.map (P.apply p) l))
	
    let apply_ok2 =
      "apply/ok2" >::
	(fun () ->
	  let p = P.of_lists [0; 1; 2] ["a"; "b"; "c"]
	  and l = [ 0; 1; 2 ] in
	  assert_equal [ "a"; "b"; "c" ] (List.map (P.apply p) l))

    let apply_shadowed =
      "apply/shadowed" >::
	(fun () ->
	  let p = P.of_list [ (0,"a"); (1,"b"); (2,"c"); (1,"d") ]
	  and l = [ 0; 1; 2 ] in
	  assert_equal [ "a"; "d"; "c" ] (List.map (P.apply p) l))

    let apply_shadowed2 =
      "apply/shadowed2" >::
	(fun () ->
	  let p = P.of_lists [0; 1; 2; 1] ["a"; "b"; "c"; "d"]
	  and l = [ 0; 1; 2 ] in
	  assert_equal [ "a"; "d"; "c" ] (List.map (P.apply p) l))

    let apply_mismatch =
      "apply/mismatch" >::
	(fun () ->
	  assert_raises
	    (Invalid_argument "Partial.of_lists: length mismatch")
	    (fun () -> P.of_lists [0; 1; 2] ["a"; "b"; "c"; "d"]))

    let suite_apply =
      "apply" >:::
	[apply_ok;
	 apply_ok2;
	 apply_shadowed;
	 apply_shadowed2;
	 apply_mismatch]

    let auto_ok =
      "auto/ok" >::
	(fun () ->
	  let p = P.of_list [ (0,10); (1,11)]
	  and l = [ 0; 1; 2 ] in
	  assert_equal [ 10; 11; 2 ] (List.map (P.auto p) l))

    let suite_auto =
      "auto" >:::
	[auto_ok]

    let apply_with_fallback_ok =
      "apply_with_fallback/ok" >::
	(fun () ->
	  let p = P.of_list [ (0,10); (1,11)]
	  and l = [ 0; 1; 2 ] in
	  assert_equal
	    [ 10; 11; -2 ] (List.map (P.apply_with_fallback (fun n -> - n) p) l))

    let suite_apply_with_fallback =
      "apply_with_fallback" >:::
	[apply_with_fallback_ok]

    let suite =
      "Partial" >:::
	[suite_apply;
	 suite_auto;
	 suite_apply_with_fallback]

    let time () =
      ()

  end
