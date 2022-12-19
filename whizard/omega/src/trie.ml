(* trie.ml --

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

(* \thocwmodulesection{Monomorphically} *)

module type T =
  sig
    type key
    type (+'a) t
    val empty : 'a t
    val is_empty : 'a t -> bool
    val add : key -> 'a -> 'a t -> 'a t
    val find : key -> 'a t -> 'a
    val remove : key -> 'a t -> 'a t
    val mem : key -> 'a t -> bool
    val map : ('a -> 'b) -> 'a t -> 'b t
    val mapi : (key -> 'a -> 'b) -> 'a t -> 'b t
    val iter : (key -> 'a -> unit) -> 'a t -> unit
    val fold : (key -> 'a -> 'b -> 'b) -> 'a t -> 'b -> 'b
    val longest : key -> 'a t -> 'a option * key
    val shortest : key -> 'a t -> 'a option * key
    val compare : ('a -> 'a -> int) -> 'a t -> 'a t -> int
    val equal : ('a -> 'a -> bool) -> 'a t -> 'a t -> bool
    val export : (int -> unit) -> (int -> unit) ->
      (int -> key -> unit) -> (int -> key -> 'a -> unit) -> 'a t -> unit

  end

(* O'Caml's [Map.S] prior to Version 3.12: *)

module type Map_S =
  sig
    type key
    type (+'a) t
    val empty: 'a t
    val is_empty: 'a t -> bool
    val add: key -> 'a -> 'a t -> 'a t
    val find: key -> 'a t -> 'a
    val remove: key -> 'a t -> 'a t
    val mem: key -> 'a t -> bool
    val iter: (key -> 'a -> unit) -> 'a t -> unit
    val map: ('a -> 'b) -> 'a t -> 'b t
    val mapi: (key -> 'a -> 'b) -> 'a t -> 'b t
    val fold: (key -> 'a -> 'b -> 'b) -> 'a t -> 'b -> 'b
    val compare: ('a -> 'a -> int) -> 'a t -> 'a t -> int
    val equal: ('a -> 'a -> bool) -> 'a t -> 'a t -> bool
  end

module Make (M : Map_S) : (T with type key = M.key list) =
  struct

(* Derived from SML code by Chris Okasaki~\cite{Okasaki:1998:book}.  *)

    type key = M.key list

    type 'a t = Trie of 'a option * 'a t M.t

    let empty = Trie (None, M.empty)

    let is_empty = function
      | Trie (None, m) -> M.is_empty m
      | _ -> false

    let rec add key data trie =
      match key, trie with
      | [], Trie (_, children) -> Trie (Some data, children)
      | k :: rest, Trie (node, children) ->
          let t = try M.find k children with Not_found -> empty in
          Trie (node, M.add k (add rest data t) children)

    let rec find key trie =
      match key, trie with
      | [], Trie (None, _) -> raise Not_found
      | [], Trie (Some data, _) -> data
      | k :: rest, Trie (_, children) -> find rest (M.find k children)

(* The rest is my own fault \ldots{} *)

    let find1 k children =
      try Some (M.find k children) with Not_found -> None

    let add_non_empty k t children =
      if t = empty then
        M.remove k children
      else
        M.add k t children

    let rec remove key trie =
      match key, trie with
      | [], Trie (_, children) -> Trie (None, children)
      | k :: rest, (Trie (node, children) as orig) ->
          match find1 k children with
          | None -> orig
          | Some t -> Trie (node, add_non_empty k (remove rest t) children)

    let rec mem key trie =
      match key, trie with
      | [], Trie (None, _) -> false
      | [], Trie (Some data, _) -> true
      | k :: rest, Trie (_, children) ->
          match find1 k children with
          | None -> false
          | Some t -> mem rest t

    let rec map f = function
      | Trie (Some data, children) ->
          Trie (Some (f data), M.map (map f) children)
      | Trie (None, children) -> Trie (None, M.map (map f) children)

    let rec mapi' key f = function
      | Trie (Some data, children) ->
          Trie (Some (f key data), descend key f children)
      | Trie (None, children) -> Trie (None, descend key f children)
    and descend key f = M.mapi (fun k -> mapi' (key @ [k]) f)
    let mapi f = mapi' [] f

    let rec iter' key f = function
      | Trie (Some data, children) -> f key data; descend key f children
      | Trie (None, children) -> descend key f children
    and descend key f = M.iter (fun k -> iter' (key @ [k]) f)
    let iter f = iter' [] f

    let rec fold' key f t acc =
      match t with
      | Trie (Some data, children) -> descend key f children (f key data acc)
      | Trie (None, children) -> descend key f children acc
    and descend key f = M.fold (fun k -> fold' (key @ [k]) f)
    let fold f t acc = fold' [] f t acc

    let rec longest' partial partial_rest key trie =
      match key, trie with
      | [], Trie (data, _) -> (data, [])
      | k :: rest, Trie (data, children) ->
          match data, find1 k children with
          | None, None -> (partial, partial_rest)
          | Some _, None -> (data, key)
          | _, Some t -> longest' partial partial_rest rest t
    let longest key = longest' None key key

    let rec shortest' partial partial_rest key trie =
      match key, trie with
      | [], Trie (data, _) -> (data, [])
      | k :: rest, Trie (Some _ as data, children) -> (data, key)
      | k :: rest, Trie (None, children) ->
          match find1 k children with
          | None -> (partial, partial_rest)
          | Some t -> shortest' partial partial_rest rest t
    let shortest key = shortest' None key key

(* \thocwmodulesection{O'Mega customization} *)

    let rec export' n key f_open f_close f_descend f_match = function
      | Trie (Some data, children) ->
          f_match n key data;
          if children <> M.empty then
            descend n key f_open f_close f_descend f_match children
      | Trie (None, children) ->
          if children <> M.empty then begin
            f_descend n key;
            descend n key f_open f_close f_descend f_match children
          end
    and descend n key f_open f_close f_descend f_match children =
      f_open n;
      M.iter (fun k ->
        export' (succ n) (k :: key) f_open f_close f_descend f_match) children;
      f_close n

    let export f_open f_close f_descend f_match =
      export' 0 [] f_open f_close f_descend f_match

    let compare _ _ _ =
      failwith "incomplete"

(*i
    let compare cmp m1 m2 =
      let rec compare_aux e1 e2 =
        match (e1, e2) with
        | (End, End) -> 0
        | (End, _)  -> -1
        | (_, End) -> 1
        | (More(v1, d1, r1, e1), More(v2, d2, r2, e2)) ->
            let c = Ord.compare v1 v2 in
            if c <> 0 then c else
            let c = cmp d1 d2 in
            if c <> 0 then c else
            compare_aux (cons_enum r1 e1) (cons_enum r2 e2) in
      compare_aux (cons_enum m1 End) (cons_enum m2 End)
i*)

    let equal _ _ _ =
      failwith "incomplete"

(*i
    let equal cmp m1 m2 =
      let rec equal_aux e1 e2 =
        match (e1, e2) with
        | (End, End) -> true
        | (End, _)  -> false
        | (_, End) -> false
        | (More(v1, d1, r1, e1), More(v2, d2, r2, e2)) ->
            Ord.compare v1 v2 = 0 && cmp d1 d2 &&
            equal_aux (cons_enum r1 e1) (cons_enum r2 e2) in
      equal_aux (cons_enum m1 End) (cons_enum m2 End)
i*)

  end

module MakeMap (M : Map_S) : (Map_S with type key = M.key list) = Make(M)

(* \thocwmodulesection{Polymorphically} *)

module type Poly =
  sig
    type ('a, 'b) t
    val empty : ('a, 'b) t
    val add : ('a -> 'a -> int) -> 'a list -> 'b -> ('a, 'b) t -> ('a, 'b) t
    val find : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> 'b
    val remove : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> ('a, 'b) t
    val mem : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> bool
    val map : ('b -> 'c) -> ('a, 'b) t -> ('a, 'c) t
    val mapi : ('a list -> 'b -> 'c) -> ('a, 'b) t -> ('a, 'c) t
    val iter : ('a list -> 'b -> unit) -> ('a, 'b) t -> unit
    val fold : ('a list -> 'b -> 'c -> 'c) -> ('a, 'b) t -> 'c -> 'c
    val longest : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> 'b option * 'a list
    val shortest : ('a -> 'a -> int) -> 'a list -> ('a, 'b) t -> 'b option * 'a list
    val export : (int -> unit) -> (int -> unit) ->
      (int -> 'a list -> unit) -> (int -> 'a list -> 'b -> unit) -> ('a, 'b) t -> unit
  end

module MakePoly (M : Pmap.T) : Poly =
  struct

(* Derived from SML code by Chris Okasaki~\cite{Okasaki:1998:book}.  *)


    type ('a, 'b) t = Trie of 'b option * ('a, ('a, 'b) t) M.t

    let empty = Trie (None, M.empty)

    let rec add cmp key data trie =
      match key, trie with
      | [], Trie (_, children) -> Trie (Some data, children)
      | k :: rest, Trie (node, children) ->
          let t = try M.find cmp k children with Not_found -> empty in
          Trie (node, M.add cmp k (add cmp rest data t) children)

    let rec find cmp key trie =
      match key, trie with
      | [], Trie (None, _) -> raise Not_found
      | [], Trie (Some data, _) -> data
      | k :: rest, Trie (_, children) -> find cmp rest (M.find cmp k children)

(* The rest is my own fault \ldots{} *)

    let find1 cmp k children =
      try Some (M.find cmp k children) with Not_found -> None

    let add_non_empty cmp k t children =
      if t = empty then
        M.remove cmp k children
      else
        M.add cmp k t children

    let rec remove cmp key trie =
      match key, trie with
      | [], Trie (_, children) -> Trie (None, children)
      | k :: rest, (Trie (node, children) as orig) ->
          match find1 cmp k children with
          | None -> orig
          | Some t -> Trie (node, add_non_empty cmp k (remove cmp rest t) children)

    let rec mem cmp key trie =
      match key, trie with
      | [], Trie (None, _) -> false
      | [], Trie (Some data, _) -> true
      | k :: rest, Trie (_, children) ->
          match find1 cmp k children with
          | None -> false
          | Some t -> mem cmp rest t

    let rec map f = function
      | Trie (Some data, children) ->
          Trie (Some (f data), M.map (map f) children)
      | Trie (None, children) -> Trie (None, M.map (map f) children)

    let rec mapi' key f = function
      | Trie (Some data, children) ->
          Trie (Some (f key data), descend key f children)
      | Trie (None, children) -> Trie (None, descend key f children)
    and descend key f = M.mapi (fun k -> mapi' (key @ [k]) f)
    let mapi f = mapi' [] f

    let rec iter' key f = function
      | Trie (Some data, children) -> f key data; descend key f children
      | Trie (None, children) -> descend key f children
    and descend key f = M.iter (fun k -> iter' (key @ [k]) f)
    let iter f = iter' [] f

    let rec fold' key f t acc =
      match t with
      | Trie (Some data, children) -> descend key f children (f key data acc)
      | Trie (None, children) -> descend key f children acc
    and descend key f = M.fold (fun k -> fold' (key @ [k]) f)
    let fold f t acc = fold' [] f t acc

    let rec longest' cmp partial partial_rest key trie =
      match key, trie with
      | [], Trie (data, _) -> (data, [])
      | k :: rest, Trie (data, children) ->
          match data, find1 cmp k children with
          | None, None -> (partial, partial_rest)
          | Some _, None -> (data, key)
          | _, Some t -> longest' cmp partial partial_rest rest t
    let longest cmp key = longest' cmp None key key

    let rec shortest' cmp partial partial_rest key trie =
      match key, trie with
      | [], Trie (data, _) -> (data, [])
      | k :: rest, Trie (Some _ as data, children) -> (data, key)
      | k :: rest, Trie (None, children) ->
          match find1 cmp k children with
          | None -> (partial, partial_rest)
          | Some t -> shortest' cmp partial partial_rest rest t
    let shortest cmp key = shortest' cmp None key key

(* \thocwmodulesection{O'Mega customization} *)

    let rec export' n key f_open f_close f_descend f_match = function
      | Trie (Some data, children) ->
          f_match n key data;
          if children <> M.empty then
            descend n key f_open f_close f_descend f_match children
      | Trie (None, children) ->
          if children <> M.empty then begin
            f_descend n key;
            descend n key f_open f_close f_descend f_match children
          end
    and descend n key f_open f_close f_descend f_match children =
      f_open n;
      M.iter (fun k ->
        export' (succ n) (k :: key) f_open f_close f_descend f_match) children;
      f_close n

    let export f_open f_close f_descend f_match =
      export' 0 [] f_open f_close f_descend f_match

  end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
