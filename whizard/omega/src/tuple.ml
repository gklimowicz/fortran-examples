(* tuple.ml --

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

module type Mono =
  sig
    type 'a t
    val arity : 'a t -> int
    val max_arity : unit -> int
    val compare : ('a -> 'a -> int) -> 'a t -> 'a t -> int
    val for_all : ('a -> bool) -> 'a t -> bool
    val map : ('a -> 'b) -> 'a t -> 'b t
    val iter : ('a -> unit) -> 'a t -> unit
    val fold_left : ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a
    val fold_right : ('a -> 'b -> 'b) -> 'a t -> 'b -> 'b
    val fold_left_internal : ('a -> 'a -> 'a) -> 'a t -> 'a
    val fold_right_internal : ('a -> 'a -> 'a) -> 'a t -> 'a
    val map2 : ('a -> 'b -> 'c) -> 'a t -> 'b t -> 'c t
    val split : ('a * 'b) t -> 'a t * 'b t
    val product : 'a list t -> 'a t list
    val product_fold : ('a t -> 'b -> 'b) -> 'a list t -> 'b -> 'b
    val power : ?truncate:int -> 'a list -> 'a t list
    val power_fold : ?truncate:int -> ('a t -> 'b -> 'b) -> 'a list -> 'b -> 'b
    type 'a graded = 'a list array
    val graded_sym_power : int -> 'a graded -> 'a t list
    val graded_sym_power_fold : int -> ('a t -> 'b -> 'b) -> 'a graded ->
      'b -> 'b
    val to_list : 'a t -> 'a list
    val of2_kludge : 'a -> 'a -> 'a t
  end

module type Poly =
  sig
    include Mono
    exception Mismatched_arity
    exception No_termination
  end

(* \thocwmodulesection{Typesafe Combinatorics} *)

(* Wrap the combinatorical functions with varying arities into typesafe functions
   with fixed arities.  We could provide specialized implementations, but since
   we \emph{know} that [Impossible] is \emph{never} raised, the present approach
   is just as good (except for a tiny inefficiency).  *)

exception Impossible of string
let impossible name = raise (Impossible name)

let choose2 set =
  List.map (function [x; y] -> (x, y) | _ -> impossible "choose2")
      (Combinatorics.choose 2 set)

let choose3 set =
  List.map (function [x; y; z] -> (x, y, z) | _ -> impossible "choose3")
    (Combinatorics.choose 3 set)

(* \thocwmodulesection{Pairs} *)

module type Binary =
    sig
      include Poly (* should become [Mono]! *)
      val of2 : 'a -> 'a -> 'a t
    end

module Binary =
  struct

    type 'a t = 'a * 'a

    let arity _ = 2
    let max_arity () = 2

    let of2 x y = (x, y)

    let compare cmp (x1, y1) (x2, y2) =
      let cx = cmp x1 x2 in
      if cx <> 0 then
        cx
      else
        cmp y1 y2

    let for_all p (x, y) = p x && p y

    let map f (x, y) = (f x, f y)
    let iter f (x, y) = f x; f y
    let fold_left f init (x, y) = f (f init x) y
    let fold_right f (x, y) init = f x (f y init)
    let fold_left_internal f (x, y) = f x y
    let fold_right_internal f (x, y) = f x y

    exception Mismatched_arity
    let map2 f (x1, y1) (x2, y2) = (f x1 x2, f y1 y2)

    let split ((x1, x2), (y1, y2)) = ((x1, y1), (x2, y2))

    let product (lx, ly) =
      Product.list2 (fun x y -> (x, y)) lx ly
    let product_fold f (lx, ly) init =
      Product.fold2 (fun x y -> f (x, y)) lx ly init

    let power ?truncate l =
      match truncate with
      | None -> product (l, l)
      | Some n ->
	  if n >= 2 then
	    product (l, l)
	  else
	    invalid_arg "Tuple.Binary.power: truncate < 2"

    let power_fold ?truncate f l =
      match truncate with
      | None -> product_fold f (l, l)
      | Some n ->
	  if n >= 2 then
	    product_fold f (l, l)
	  else
	    invalid_arg "Tuple.Binary.power_fold: truncate < 2"

(* In the special case of binary fusions, the implementation is very concise. *)
    type 'a graded = 'a list array

    let fuse2 f set (i, j) acc =
      if i = j then
        List.fold_right (fun (x, y) -> f x y) (choose2 set.(pred i)) acc
      else
        Product.fold2 f set.(pred i) set.(pred j) acc

    let graded_sym_power_fold rank f set acc =
      let max_rank = Array.length set in
      List.fold_right (fuse2 (fun x y -> f (of2 x y)) set)
        (Partition.pairs rank 1 max_rank) acc

    let graded_sym_power rank set =
      graded_sym_power_fold rank (fun pair acc -> pair :: acc) set []

    let to_list (x, y) = [x; y]

    let of2_kludge  = of2

    exception No_termination
  end

(* \thocwmodulesection{Triples} *)

module type Ternary =
    sig
      include Mono 
      val of3 : 'a -> 'a -> 'a -> 'a t
    end

module Ternary =
  struct

    type 'a t = 'a * 'a * 'a

    let arity _ = 3
    let max_arity () = 3

    let of3 x y z = (x, y, z)

    let compare cmp (x1, y1, z1) (x2, y2, z2) =
      let cx = cmp x1 x2 in
      if cx <> 0 then
        cx
      else
        let cy = cmp y1 y2 in
        if cy <> 0 then
          cy
        else
          cmp z1 z2

    let for_all p (x, y, z) = p x && p y && p z

    let map f (x, y, z) = (f x, f y, f z)
    let iter f (x, y, z) = f x; f y; f z
    let fold_left f init (x, y, z) = f (f (f init x) y) z
    let fold_right f (x, y, z) init = f x (f y (f z init))
    let fold_left_internal f (x, y, z) = f (f x y) z
    let fold_right_internal f (x, y, z) = f x (f y z)

    exception Mismatched_arity
    let map2 f (x1, y1, z1) (x2, y2, z2) = (f x1 x2, f y1 y2, f z1 z2)

    let split ((x1, x2), (y1, y2), (z1, z2)) = ((x1, y1, z1), (x2, y2, z2))

    let product (lx,ly,lz) =
      Product.list3 (fun x y z -> (x, y, z)) lx ly lz
    let product_fold f (lx, ly, lz) init =
      Product.fold3 (fun x y z -> f (x, y, z)) lx ly lz init

    let power ?truncate l =
      match truncate with
      | None -> product (l, l, l)
      | Some n ->
	  if n >= 3 then
	    product (l, l, l)
	  else
	    invalid_arg "Tuple.Ternary.power: truncate < 3"

    let power_fold ?truncate f l =
      match truncate with
      | None -> product_fold f (l, l, l)
      | Some n ->
	  if n >= 3 then
	    product_fold f (l, l, l)
	  else
	    invalid_arg "Tuple.Ternary.power_fold: truncate < 3"

    type 'a graded = 'a list array

    let fuse3 f set (i, j, k) acc =
      if i = j then begin
        if j = k then
          List.fold_right (fun (x, y, z) -> f x y z) (choose3 set.(pred i)) acc
        else
          Product.fold2 (fun (x, y) z -> f x y z)
            (choose2 set.(pred i)) set.(pred k) acc
      end else begin
        if j = k then
          Product.fold2 (fun x (y, z) -> f x y z)
            set.(pred i) (choose2 set.(pred j)) acc
        else
          Product.fold3 (fun x y z -> f x y z)
            set.(pred i) set.(pred j) set.(pred k) acc
      end

    let graded_sym_power_fold rank f set acc =
      let max_rank = Array.length set in
      List.fold_right (fuse3 (fun x y z -> f (of3 x y z)) set)
        (Partition.triples rank 1 max_rank) acc

    let graded_sym_power rank set =
      graded_sym_power_fold rank (fun pair acc -> pair :: acc) set []

    let to_list (x, y, z) = [x; y; z]

    let of2_kludge _ = failwith "Tuple.Ternary.of2_kludge"

  end

(* \thocwmodulesection{Pairs and Triples} *)

type 'a pair_or_triple = T2 of 'a * 'a | T3 of 'a * 'a *'a

module type Mixed23 =
    sig
      include Poly
      val of2 : 'a -> 'a -> 'a t
      val of3 : 'a -> 'a -> 'a -> 'a t
    end

module Mixed23 =
  struct

    type 'a t = 'a pair_or_triple

    let arity = function
      | T2 _ -> 2
      | T3 _ -> 3
    let max_arity () = 3

    let of2 x y = T2 (x, y)
    let of3 x y z = T3 (x, y, z)

    let compare cmp m1 m2 =
      match m1, m2 with
      | T2 _, T3 _ -> -1
      | T3 _, T2 _ -> 1
      | T2 (x1, y1), T2 (x2, y2) ->
          let cx = cmp x1 x2 in
          if cx <> 0 then
            cx
          else
            cmp y1 y2
      | T3 (x1, y1, z1), T3 (x2, y2, z2) ->
          let cx = cmp x1 x2 in
          if cx <> 0 then
            cx
          else
            let cy = cmp y1 y2 in
            if cy <> 0 then
              cy
            else
              cmp z1 z2

    let for_all p = function
      | T2 (x, y) -> p x && p y
      | T3 (x, y, z) -> p x && p y && p z

    let map f = function
      | T2 (x, y) -> T2 (f x, f y)
      | T3 (x, y, z) -> T3 (f x, f y, f z)

    let iter f = function
      | T2 (x, y) -> f x; f y 
      | T3 (x, y, z) -> f x; f y; f z

    let fold_left f init = function
      | T2 (x, y) -> f (f init x) y
      | T3 (x, y, z) -> f (f (f init x) y) z

    let fold_right f m init =
      match m with
      | T2 (x, y) -> f x (f y init)
      | T3 (x, y, z) -> f x (f y (f z init))

    let fold_left_internal f m =
      match m with
      | T2 (x, y) -> f x y
      | T3 (x, y, z) -> f (f x y) z

    let fold_right_internal f m =
      match m with
      | T2 (x, y) -> f x y
      | T3 (x, y, z) -> f x (f y z)

    exception Mismatched_arity
    let map2 f m1 m2 =
      match m1, m2 with
      | T2 (x1, y1), T2 (x2, y2) -> T2 (f x1 x2, f y1 y2)
      | T3 (x1, y1, z1), T3 (x2, y2, z2) -> T3 (f x1 x2, f y1 y2, f z1 z2)
      | T2 _, T3 _ | T3 _, T2 _ -> raise Mismatched_arity

    let split = function
      | T2 ((x1, x2), (y1, y2)) -> (T2 (x1, y1), T2 (x2, y2))
      | T3 ((x1, x2), (y1, y2), (z1, z2)) -> (T3 (x1, y1, z1), T3 (x2, y2, z2))

    let product = function
      | T2 (lx, ly) -> Product.list2 (fun x y -> T2 (x, y)) lx ly
      | T3 (lx, ly, lz) -> Product.list3 (fun x y z -> T3 (x, y, z)) lx ly lz
    let product_fold f m init =
      match m with
      | T2 (lx, ly) -> Product.fold2 (fun x y -> f (T2 (x, y))) lx ly init
      | T3 (lx, ly, lz) ->
          Product.fold3 (fun x y z -> f (T3 (x, y, z))) lx ly lz init

    exception No_termination

    let power_fold23 f l init =
      product_fold f (T2 (l, l)) (product_fold f (T3 (l, l, l)) init)

    let power_fold2 f l init =
      product_fold f (T2 (l, l)) init

    let power_fold ?truncate f l init =
      match truncate with
      | None -> power_fold23 f l init
      | Some n ->
	  if n >= 3 then
	    power_fold23 f l init
	  else if n = 2 then
	    power_fold2 f l init
	  else
	    invalid_arg "Tuple.Mixed23.power_fold: truncate < 2"

    let power ?truncate l =
      power_fold ?truncate (fun m acc -> m :: acc) l []

    type 'a graded = 'a list array

    let graded_sym_power_fold rank f set acc =
      let max_rank = Array.length set in
      List.fold_right (Binary.fuse2 (fun x y -> f (of2 x y)) set)
        (Partition.pairs rank 1 max_rank)
        (List.fold_right (Ternary.fuse3 (fun x y z -> f (of3 x y z)) set)
           (Partition.triples rank 1 max_rank) acc)

    let graded_sym_power rank set =
      graded_sym_power_fold rank (fun pair acc -> pair :: acc) set []

    let to_list = function
      | T2 (x, y) -> [x; y]
      | T3 (x, y, z) -> [x; y; z]

    let of2_kludge = of2

  end

(* \thocwmodulesection{\ldots{} and All The Rest} *)

module type Nary =
    sig
      include Poly
      val of2 : 'a -> 'a -> 'a t
      val of3 : 'a -> 'a -> 'a -> 'a t
      val of_list : 'a list -> 'a t
    end

module Nary (A : sig val max_arity : unit -> int end) =
  struct

    type 'a t = 'a * 'a list

    let arity (_, y) = succ (List.length y)

    let max_arity () =
      try A.max_arity () with _ -> -1

    let of2 x y = (x, [y])
    let of3 x y z = (x, [y; z])

    let of_list = function
      | x :: y -> (x, y)
      | [] -> invalid_arg "Tuple.Nary.of_list: empty"

    let compare cmp (x1, y1) (x2, y2) =
      let c = cmp x1 x2 in
      if c <> 0 then
        c
      else
        ThoList.compare ~cmp y1 y2

    let for_all p (x, y) = p x && List.for_all p y

    let map f (x, y) = (f x, List.map f y)
    let iter f (x, y) = f x; List.iter f y
    let fold_left f init (x, y) = List.fold_left f (f init x) y
    let fold_right f (x, y) init = f x (List.fold_right f y init)
    let fold_left_internal f (x, y) = List.fold_left f x y
    let fold_right_internal f (x, y) =
      match List.rev y with
      | [] -> x
      | y0 :: y_sans_y0 ->
          f x (List.fold_right f (List.rev y_sans_y0) y0)

    exception Mismatched_arity
    let map2 f (x1, y1) (x2, y2) =
      try (f x1 x2, List.map2 f y1 y2) with
      | Invalid_argument _ -> raise Mismatched_arity

    let split ((x1, x2), y12) =
      let y1, y2 = List.split y12 in
      ((x1, y1), (x2, y2))

    let product (xl, yl) =
      Product.list (function
        | x :: y -> (x, y)
        | [] -> failwith "Tuple.Nary.product") (xl :: yl)

    let product_fold f (xl, yl) init =
      Product.fold (function
        | x :: y -> f (x, y)
        | [] -> failwith "Tuple.Nary.product_fold") (xl :: yl) init

    exception No_termination

    let truncated_arity ?truncate () =
      let ma = max_arity () in
      match truncate with
      | None -> ma
      | Some n ->
	  if n < 2 then
	    invalid_arg "Tuple.Nary.power: truncate < 2"
	  else if ma >= 2 then
	    min n ma
	  else
	    n

    let power_fold ?truncate f l init =
      let ma = truncated_arity ?truncate () in
      if ma > 0 then
        List.fold_right
          (fun n -> product_fold f (l, ThoList.clone l (pred n)))
          (ThoList.range 2 ma) init
      else
        raise No_termination

    let power ?truncate l =
      power_fold ?truncate (fun t acc -> t :: acc) l []

    type 'a graded = 'a list array

    let fuse_n f set partition acc =
      let choose (n, r) = 
        Printf.printf "chose: n=%d r=%d len=%d\n"
          n r (List.length set.(pred r));
        Combinatorics.choose n set.(pred r) in
      Product.fold (fun wfs -> f (List.concat wfs))
        (List.map choose (ThoList.classify partition)) acc

    let fuse_n f set partition acc =
      let choose (n, r) = Combinatorics.choose n set.(pred r) in
      Product.fold (fun wfs -> f (List.concat wfs))
        (List.map choose (ThoList.classify partition)) acc

(* \begin{dubious}
     [graded_sym_power_fold] is well defined for unbounded arities as well: derive
     a reasonable replacement from [set].  The length of the flattened [set] is
     an upper limit, of course, but too pessimistic in most cases.
   \end{dubious} *)

    let graded_sym_power_fold rank f set acc =
      let max_rank = Array.length set in
      let degrees = ThoList.range 2 (max_arity ()) in
      let partitions =
        ThoList.flatmap
          (fun deg -> Partition.tuples deg rank 1 max_rank) degrees in
      List.fold_right (fuse_n (fun wfs -> f (of_list wfs)) set) partitions acc

    let graded_sym_power rank set =
      graded_sym_power_fold rank (fun pair acc -> pair :: acc) set []

    let to_list (x, y) = x :: y

    let of2_kludge = of2

  end

module type Bound = sig val max_arity : unit -> int end
module Unbounded_Nary = Nary (struct let max_arity () = -1 end)

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
