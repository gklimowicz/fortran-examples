(* pmap.ml --

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

module type T =
  sig
    type ('key, 'a) t
    val empty : ('key, 'a) t
    val is_empty : ('key, 'a) t -> bool
    val singleton : 'key -> 'a -> ('key, 'a) t
    val add : ('key -> 'key -> int) -> 'key -> 'a -> ('key, 'a) t -> ('key, 'a) t
    val update : ('key -> 'key -> int) -> ('a -> 'a -> 'a) ->
      'key -> 'a -> ('key, 'a) t -> ('key, 'a) t
    val cons : ('key -> 'key -> int) -> ('a -> 'a -> 'a option) ->
      'key -> 'a -> ('key, 'a) t -> ('key, 'a) t
    val find : ('key -> 'key -> int) -> 'key -> ('key, 'a) t -> 'a
    val find_opt : ('key -> 'key -> int) -> 'key -> ('key, 'a) t -> 'a option
    val choose : ('key, 'a) t -> 'key * 'a
    val choose_opt : ('key, 'a) t -> ('key * 'a) option
    val uncons : ('key, 'a) t -> 'key * 'a * ('key, 'a) t
    val uncons_opt : ('key, 'a) t -> ('key * 'a * ('key, 'a) t) option
    val elements : ('key, 'a) t -> ('key * 'a) list
    val mem :  ('key -> 'key -> int) -> 'key -> ('key, 'a) t -> bool
    val remove : ('key -> 'key -> int) -> 'key -> ('key, 'a) t -> ('key, 'a) t
    val union : ('key -> 'key -> int) -> ('a -> 'a -> 'a) ->
      ('key, 'a) t -> ('key, 'a) t -> ('key, 'a) t
    val compose : ('key -> 'key -> int) -> ('a -> 'a -> 'a option) ->
      ('key, 'a) t -> ('key, 'a) t -> ('key, 'a) t
    val iter : ('key -> 'a -> unit) -> ('key, 'a) t -> unit
    val map : ('a -> 'b) -> ('key, 'a) t -> ('key, 'b) t
    val mapi : ('key -> 'a -> 'b) -> ('key, 'a) t -> ('key, 'b) t
    val fold : ('key -> 'a -> 'b -> 'b) -> ('key, 'a) t -> 'b -> 'b
    val compare : ('key -> 'key -> int) -> ('a -> 'a -> int) ->
      ('key, 'a) t -> ('key, 'a) t -> int
    val canonicalize : ('key -> 'key -> int) -> ('key, 'a) t -> ('key, 'a) t
  end

module Tree  =
  struct
    type ('key, 'a) t =
      | Empty
      | Node of ('key, 'a) t * 'key * 'a * ('key, 'a) t * int

    let empty = Empty

    let is_empty = function
      | Empty -> true
      | _ -> false

    let singleton k d =
      Node (Empty, k, d, Empty, 1)

    let height = function
      | Empty -> 0
      | Node (_,_,_,_,h) -> h

    let create l x d r =
      let hl = height l and hr = height r in
      Node (l, x, d, r, (if hl >= hr then hl + 1 else hr + 1))

    let bal l x d r =
      let hl = match l with Empty -> 0 | Node (_,_,_,_,h) -> h in
      let hr = match r with Empty -> 0 | Node (_,_,_,_,h) -> h in
      if hl > hr + 2 then begin
        match l with
        | Empty -> invalid_arg "Map.bal"
        | Node (ll, lv, ld, lr, _) ->
            if height ll >= height lr then
              create ll lv ld (create lr x d r)
            else begin
              match lr with
              | Empty -> invalid_arg "Map.bal"
              | Node (lrl, lrv, lrd, lrr, _)->
                  create (create ll lv ld lrl) lrv lrd (create lrr x d r)
            end
      end else if hr > hl + 2 then begin
        match r with
        | Empty -> invalid_arg "Map.bal"
        | Node (rl, rv, rd, rr, _) ->
            if height rr >= height rl then
              create (create l x d rl) rv rd rr
            else begin
              match rl with
              | Empty -> invalid_arg "Map.bal"
              | Node (rll, rlv, rld, rlr, _) ->
                  create (create l x d rll) rlv rld (create rlr rv rd rr)
            end
      end else
        Node (l, x, d, r, (if hl >= hr then hl + 1 else hr + 1))

    let rec join l x d r =
      match bal l x d r with
      | Empty -> invalid_arg "Pmap.join"
      | Node (l', x', d', r', _) as t' ->
          let d = height l' - height r' in
          if d < -2 || d > 2 then
            join l' x' d' r'
          else
            t'

(* Merge two trees [t1] and [t2] into one. All elements of [t1] must
   precede the elements of [t2]. Assumes [height t1 - height t2 <= 2]. *)

    let rec merge t1 t2 =
      match t1, t2 with
      | Empty, t -> t
      | t, Empty -> t
      | Node (l1, v1, d1, r1, h1), Node (l2, v2, d2, r2, h2) ->
          bal l1 v1 d1 (bal (merge r1 l2) v2 d2 r2)

(* Same as merge, but does not assume anything about [t1] and [t2]. *)

    let rec concat t1 t2 =
      match t1, t2 with
      | Empty, t -> t
      | t, Empty -> t
      | Node (l1, v1, d1, r1, h1), Node (l2, v2, d2, r2, h2) ->
          join l1 v1 d1 (join (concat r1 l2) v2 d2 r2)
 
(* Splitting *)

    let rec split cmp x = function
      | Empty -> (Empty, None, Empty)
      | Node (l, v, d, r, _) ->
          let c = cmp x v in
          if c = 0 then
            (l, Some d, r)
          else if c < 0 then
            let ll, vl, rl = split cmp x l in
            (ll, vl, join rl v d r)
          else (* [if c > 0 then] *)
            let lr, vr, rr = split cmp x r in
            (join l v d lr, vr, rr)

    let rec find cmp x = function
      | Empty -> raise Not_found
      | Node (l, v, d, r, _) ->
          let c = cmp x v in
          if c = 0 then
            d
          else if c < 0 then
            find cmp x l
          else (* [if c > 0] *)
            find cmp x r

    let rec find_opt cmp x = function
      | Empty -> None
      | Node (l, v, d, r, _) ->
          let c = cmp x v in
          if c = 0 then
            Some d
          else if c < 0 then
            find_opt cmp x l
          else (* [if c > 0] *)
            find_opt cmp x r

    let rec mem cmp x = function
      | Empty -> false
      | Node (l, v, d, r, _) ->
          let c = cmp x v in
          if c = 0 then
            true
          else if c < 0 then
            mem cmp x l
          else (* [if c > 0] *)
            mem cmp x r

    let choose = function
      | Empty -> raise Not_found
      | Node (l, v, d, r, _) -> (v, d)

    let choose_opt = function
      | Empty -> None
      | Node (l, v, d, r, _) -> Some (v, d)

    let uncons = function
      | Empty -> raise Not_found
      | Node (l, v, d, r, h) -> (v, d, merge l r)

    let uncons_opt = function
      | Empty -> None
      | Node (l, v, d, r, h) -> Some (v, d, merge l r)

    let rec remove cmp x = function
      | Empty -> Empty
      | Node (l, v, d, r, h) ->
          let c = cmp x v in
          if c = 0 then
            merge l r
          else if c < 0 then
            bal (remove cmp x l) v d r
          else (* [if c > 0] *)
            bal l v d (remove cmp x r)

    let rec cons cmp resolve x data' = function
      | Empty -> Node (Empty, x, data', Empty, 1)
      | Node (l, v, data, r, h) ->
          let c = cmp x v in
          if c = 0 then
            match resolve data' data with
            | Some data'' -> Node (l, x, data'', r, h)
            | None -> merge l r
          else if c < 0 then
            bal (cons cmp resolve x data' l) v data r
          else (* [if c > 0] *)
            bal l v data (cons cmp resolve x data' r)

    let rec update cmp resolve x data' = function
      | Empty -> Node (Empty, x, data', Empty, 1)
      | Node (l, v, data, r, h) ->
          let c = cmp x v in
          if c = 0 then
            Node (l, x, resolve data' data, r, h)
          else if c < 0 then
            bal (update cmp resolve x data' l) v data r
          else (* [if c > 0] *)
            bal l v data (update cmp resolve x data' r)

    let add cmp x data = update cmp (fun n o -> n) x data

    let rec compose cmp resolve s1 s2 =
      match s1, s2 with
      | Empty, t2 -> t2
      | t1, Empty -> t1
      | Node (l1, v1, d1, r1, h1), Node (l2, v2, d2, r2, h2) ->
          if h1 >= h2 then
            if h2 = 1 then
              cons cmp (fun o n -> resolve n o) v2 d2 s1
            else begin
              match split cmp v1 s2 with
              | l2', None, r2' ->
                  join (compose cmp resolve l1 l2') v1 d1
                    (compose cmp resolve r1 r2')
              | l2', Some d, r2' ->
                  begin match resolve d1 d with
                  | None ->
                      concat (compose cmp resolve l1 l2')
                        (compose cmp resolve r1 r2')
                  | Some d ->
                      join (compose cmp resolve l1 l2') v1 d
                        (compose cmp resolve r1 r2')
                  end
            end
          else
            if h1 = 1 then
              cons cmp resolve v1 d1 s2
            else begin
              match split cmp v2 s1 with
              | l1', None, r1' ->
                  join (compose cmp resolve l1' l2) v2 d2
                    (compose cmp resolve r1' r2)
              | l1', Some d, r1' ->
                  begin match resolve d d2 with
                  | None ->
                      concat (compose cmp resolve l1' l2)
                        (compose cmp resolve r1' r2)
                  | Some d ->
                      join (compose cmp resolve l1' l2) v2 d
                        (compose cmp resolve r1' r2)
                  end
            end

    let rec union cmp resolve s1 s2 =
      match s1, s2 with
      | Empty, t2 -> t2
      | t1, Empty -> t1
      | Node (l1, v1, d1, r1, h1), Node (l2, v2, d2, r2, h2) ->

          if h1 >= h2 then
            if h2 = 1 then
              update cmp (fun o n -> resolve n o) v2 d2 s1
            else begin
              match split cmp v1 s2 with
              | l2', None, r2' ->
                  join (union cmp resolve l1 l2') v1 d1
                    (union cmp resolve r1 r2')
              | l2', Some d, r2' ->
                  join (union cmp resolve l1 l2') v1 (resolve d1 d)
                    (union cmp resolve r1 r2')
            end
          else
            if h1 = 1 then
              update cmp resolve v1 d1 s2
            else begin
              match split cmp v2 s1 with
              | l1', None, r1' ->
                  join (union cmp resolve l1' l2) v2 d2
                    (union cmp resolve r1' r2)
              | l1', Some d, r1' ->
                  join (union cmp resolve l1' l2) v2 (resolve d d2)
                    (union cmp resolve r1' r2)
            end

    let rec iter f = function
      | Empty -> ()
      | Node (l, v, d, r, _) -> iter f l; f v d; iter f r

    let rec map f = function
      | Empty -> Empty
      | Node (l, v, d, r, h) -> Node (map f l, v, f d, map f r, h)

    let rec mapi f = function
      | Empty -> Empty
      | Node(l, v, d, r, h) -> Node (mapi f l, v, f v d, mapi f r, h)

    let rec fold f m accu =
      match m with
      | Empty -> accu
      | Node (l, v, d, r, _) -> fold f l (f v d (fold f r accu))

    let rec compare' cmp_k cmp_d l1 l2 =
      match l1, l2 with
      | [], [] -> 0
      | [], _ -> -1
      | _, [] -> 1
      | Empty :: t1, Empty :: t2 -> compare' cmp_k cmp_d t1 t2
      | Node (Empty, v1, d1, r1, _) :: t1,
          Node (Empty, v2, d2, r2, _) :: t2 ->
          let cv = cmp_k v1 v2 in
          if cv <> 0 then begin
            cv
          end else begin
            let cd = cmp_d d1 d2 in
            if cd <> 0 then
              cd
            else
              compare' cmp_k cmp_d (r1::t1) (r2::t2)
          end
      | Node (l1, v1, d1, r1, _) :: t1, t2 ->
          compare' cmp_k cmp_d (l1 :: Node (Empty, v1, d1, r1, 0) :: t1) t2
      | t1, Node (l2, v2, d2, r2, _) :: t2 ->
          compare' cmp_k cmp_d t1 (l2 :: Node (Empty, v2, d2, r2, 0) :: t2)

    let compare cmp_k cmp_d m1 m2 = compare' cmp_k cmp_d [m1] [m2]

    let rec elements' accu = function
      | Empty -> accu
      | Node (l, v, d, r, _) -> elements' ((v, d) :: elements' accu r) l

    let elements s =
      elements' [] s

    let canonicalize cmp m =
      fold (add cmp) m empty
      
  end
    
module List  =
  struct
    type ('key, 'a) t = ('key * 'a) list

    let empty = []

    let is_empty = function
      | [] -> true
      | _ -> false

    let singleton k d = [(k, d)]

    let rec cons cmp resolve k' d' = function
      | [] -> [(k', d')]
      | ((k, d) as kd :: rest) as list ->
          let c = cmp k' k in
          if c = 0 then
            match resolve d' d with
            | None -> rest
            | Some d'' -> (k', d'') :: rest
          else if c < 0 then (* [k' < k] *)
            (k', d') :: list
          else (* [if c > 0], i.\,e.~[k < k'] *)
            kd :: cons cmp resolve k' d' rest

    let rec update cmp resolve k' d' = function
      | [] -> [(k', d')]
      | ((k, d) as kd :: rest) as list ->
          let c = cmp k' k in
          if c = 0 then
            (k', resolve d' d) :: rest
          else if c < 0 then (* [k' < k] *)
            (k', d') :: list
          else (* [if c > 0], i.\,e.~[k < k'] *)
            kd :: update cmp resolve k' d' rest

    let add cmp k' d' list =
      update cmp (fun n o -> n) k' d' list

    let rec find cmp k' = function
      | [] -> raise Not_found
      | (k, d) :: rest ->
          let c = cmp k' k in
          if c = 0 then
            d
          else if c < 0 then (* [k' < k] *)
            raise Not_found
          else (* [if c > 0], i.\,e.~[k < k'] *)
            find cmp k' rest

    let rec find_opt cmp k' = function
      | [] -> None
      | (k, d) :: rest ->
          let c = cmp k' k in
          if c = 0 then
            Some d
          else if c < 0 then (* [k' < k] *)
            None
          else (* [if c > 0], i.\,e.~[k < k'] *)
            find_opt cmp k' rest

    let choose = function
      | [] -> raise Not_found
      | kd :: _ -> kd

    let rec choose_opt = function
      | [] -> None
      | kd :: _ -> Some kd

    let uncons = function
      | [] -> raise Not_found
      | (k, d) :: rest -> (k, d, rest)

    let uncons_opt = function
      | [] -> None
      | (k, d) :: rest -> Some (k, d, rest)

    let elements list = list

    let rec mem cmp k' = function
      | [] -> false
      | (k, d) :: rest ->
          let c = cmp k' k in
          if c = 0 then
            true
          else if c < 0 then (* [k' < k] *)
            false
          else (* [if c > 0], i.\,e.~[k < k'] *)
            mem cmp k' rest

    let rec remove cmp k' = function
      | [] -> []
      | ((k, d) as kd :: rest) as list ->
          let c = cmp k' k in
          if c = 0 then
            rest
          else if c < 0 then (* [k' < k] *)
            list
          else (* [if c > 0], i.\,e.~[k < k'] *)
            kd :: remove cmp k' rest

    let rec compare cmp_k cmp_d m1 m2 =
      match m1, m2 with
      | [], [] -> 0
      | [], _ -> -1
      | _, [] -> 1
      | (k1, d1) :: rest1, (k2, d2) :: rest2 ->
          let c = cmp_k k1 k2 in
          if c = 0 then begin
            let c' = cmp_d d1 d2 in
            if c' = 0 then
              compare cmp_k cmp_d rest1 rest2
            else
              c'
          end else
            c

    let rec iter f = function
      | [] -> ()
      | (k, d) :: rest -> f k d; iter f rest

    let rec map f = function
      | [] -> []
      | (k, d) :: rest -> (k, f d) :: map f rest

    let rec mapi f = function
      | [] -> []
      | (k, d) :: rest -> (k, f k d) :: mapi f rest

    let rec fold f m accu =
      match m with
      | [] -> accu
      | (k, d) :: rest -> fold f rest (f k d accu)

    let rec compose cmp resolve m1 m2 =
      match m1, m2 with
      | [], [] -> []
      | [], m -> m
      | m, [] -> m
      | ((k1, d1) as kd1 :: rest1), ((k2, d2) as kd2 :: rest2) ->
          let c = cmp k1 k2 in
          if c = 0 then
            match resolve d1 d2 with
            | None -> compose cmp resolve rest1 rest2
            | Some d -> (k1, d) :: compose cmp resolve rest1 rest2
          else if c < 0 then (* [k1 < k2] *)
            kd1 :: compose cmp resolve rest1 m2
          else (* [if c > 0], i.\,e.~[k2 < k1] *)
            kd2 :: compose cmp resolve m1 rest2

    let rec union cmp resolve m1 m2 =
      match m1, m2 with
      | [], [] -> []
      | [], m -> m
      | m, [] -> m
      | ((k1, d1) as kd1 :: rest1), ((k2, d2) as kd2 :: rest2) ->
          let c = cmp k1 k2 in
          if c = 0 then
            (k1, resolve d1 d2) :: union cmp resolve rest1 rest2
          else if c < 0 then (* [k1 < k2] *)
            kd1 :: union cmp resolve rest1 m2
          else (* [if c > 0], i.\,e.~[k2 < k1] *)
            kd2 :: union cmp resolve m1 rest2

    let canonicalize cmp x = x
      
  end

(*i
   Local Variables:
   mode:caml
   indent-tabs-mode:nil
   page-delimiter:"^(\\* .*\n"
   End:
i*)
