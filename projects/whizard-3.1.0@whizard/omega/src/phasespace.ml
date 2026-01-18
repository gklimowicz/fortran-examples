(* phasespace.ml --

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

(* \thocwmodulesection{Tools} *)

(* These are candidates for [ThoList] and not specific to phase space. *)

let rec first_match' mismatch f = function
  | [] -> None
  | x :: rest ->
      if f x then
        Some (x, List.rev_append mismatch rest)
      else
        first_match' (x :: mismatch) f rest

(* Returns $(x,X\setminus\{x\})$ if $\exists x\in X: f(x)$. *)

let first_match f l = first_match' [] f l

let rec first_pair' mismatch1 f l1 l2 =
  match l1 with
  | [] -> None
  | x1 :: rest1 ->
      begin match first_match (f x1) l2 with
      | None -> first_pair' (x1 :: mismatch1) f rest1 l2
      | Some (x2, rest2) ->
          Some ((x1, x2), (List.rev_append mismatch1 rest1, rest2))
      end

(* Returns $((x,y),(X\setminus\{x\},Y\setminus\{y\}))$ if
   $\exists x\in X: \exists y\in Y: f(x,y)$. *)

let first_pair f l1 l2 = first_pair' [] f l1 l2

(* \thocwmodulesection{Phase Space Parameterization Trees} *)

module type T =
  sig
    type momentum
    type 'a t
    type 'a decay
    val sort : ('a -> 'a -> int) -> 'a t -> 'a t
    val sort_decay : ('a -> 'a -> int) -> 'a decay -> 'a decay
    val map : ('a -> 'b) -> 'a t -> 'b t
    val map_decay : ('a -> 'b) -> 'a decay -> 'b decay
    val eval : ('a -> 'b) -> ('a -> 'b) -> ('a -> 'b -> 'b -> 'b) -> 'a t -> 'b t
    val eval_decay : ('a -> 'b) -> ('a -> 'b -> 'b -> 'b) -> 'a decay -> 'b decay
    val of_momenta : 'a -> 'a -> (momentum * 'a) list -> (momentum * 'a) t
    val decay_of_momenta : (momentum * 'a) list -> (momentum * 'a) decay
    exception Duplicate of momentum
    exception Unordered of momentum
    exception Incomplete of momentum
  end

module Make (M : Momentum.T) =
  struct

    type momentum = M.t

(* \begin{dubious}
     Finally, we came back to binary trees \ldots
   \end{dubious} *)

(* \thocwmodulesubsection{Cascade Decays} *)

    type 'a decay =
      | Leaf of 'a
      | Branch of 'a * 'a decay * 'a decay

(* \begin{dubious}
     Trees of type [(momentum * 'a option) decay] can be build easily and
     mapped to [(momentum * 'a) decay] later, once all the ['a] slots are
     filled.  A more elegant functor operating on ['b decay] directly (with
     [Momentum] style functions defined for ['b]) would not allow holes in
     the ['b decay] during the construction.
   \end{dubious} *)

    let label = function
      | Leaf p -> p
      | Branch (p, _, _) -> p

    let rec sort_decay cmp = function
      | Leaf _ as l -> l
      | Branch (p, d1, d2) ->
          let d1' = sort_decay cmp d1
          and d2' = sort_decay cmp d2 in
          if cmp (label d1') (label d2') <= 0 then
            Branch (p, d1', d2')
          else
            Branch (p, d2', d1')

    let rec map_decay f = function
      | Leaf p -> Leaf (f p)
      | Branch (p, d1, d2) -> Branch (f p, map_decay f d1, map_decay f d2)

    let rec eval_decay fl fb = function
      | Leaf p -> Leaf (fl p)
      | Branch (p, d1, d2) ->
          let d1' = eval_decay fl fb d1
          and d2' = eval_decay fl fb d2 in
          Branch (fb p (label d1') (label d2'), d1', d2')

(* Assuming that $p>p_D \lor p=p_D \lor p<p_D$, where~$p_D$ is the overall
   momentum of a decay tree~$D$, we can add $p$ to $D$ at the top or somewhere
   in the middle.  Note that `$<$' is not a total ordering and the operation
   can fail (raise exceptions) if the set of momenta does not correspond to
   a tree.  Also note that a momentum can already be present without flavor
   as a complement in a branching entered earlier.  *)

    exception Duplicate of momentum
    exception Unordered of momentum

    let rec embed_in_decay (p, f as pf) = function
      | Leaf (p', f' as pf') as d' ->
          if M.less p' p then
            Branch ((p, Some f), d', Leaf (M.sub p p', None))
          else if M.less p p' then
            Branch (pf', Leaf (p, Some f), Leaf (M.sub p' p, None))
          else if p = p' then
            begin match f' with
            | None -> Leaf (p, Some f)
            | Some _ -> raise (Duplicate p)
            end
          else
            raise (Unordered p)
      | Branch ((p', f' as pf'), d1, d2) as d' ->
          let p1, _ = label d1
          and p2, _ = label d2 in
          if M.less p' p then
            Branch ((p, Some f), d', Leaf (M.sub p p', None))
          else if M.lesseq p p1 then
            Branch (pf', embed_in_decay pf d1, d2) 
          else if M.lesseq p p2 then
            Branch (pf', d1, embed_in_decay pf d2) 
          else if p = p' then
            begin match f' with
            | None -> Branch ((p, Some f), d1, d2)
            | Some _ -> raise (Duplicate p)
            end
          else
            raise (Unordered p)

(* \begin{dubious}
     Note that both [embed_in_decay] and [embed_in_decays] below do
     \emph{not} commute, and should process `bigger' momenta first,
     because disjoint sub-momenta will create disjoint subtrees in
     the latter and raise exceptions in the former.
   \end{dubious} *)

    exception Incomplete of momentum

    let finalize1 = function
      | p, Some f -> (p, f)
      | p, None -> raise (Incomplete p)

    let finalize_decay t = map_decay finalize1 t

(* Process the momenta starting in with the highest [M.rank]: *)

    let sort_momenta plist =
      List.sort (fun (p1, _) (p2, _) -> M.compare p1 p2) plist

    let decay_of_momenta plist =
      match sort_momenta plist with
      | (p, f) :: rest ->
          finalize_decay (List.fold_right embed_in_decay rest (Leaf (p, Some f)))
      | [] -> invalid_arg "Phasespace.decay_of_momenta: empty"

(* \thocwmodulesubsection{$2\to n$ Scattering } *)

(* \begin{figure}
     \begin{center}
       \begin{fmfgraph*}(80,50)
         %%%\fmfstraight
         \fmftopn{i}{2}
         \fmfbottomn{o}{20}
         \fmf{plain,label=$p_1$}{i1,v1}
         \fmf{plain,label=$p_2$}{i2,v2}
         \fmf{phantom}{o1,v1,w1,w2,w3,w4,w5,v2,o20}
         \fmfdot{v1,v2}
         \fmfdot{w2,w4}
         \fmffreeze
         \fmfshift{(0,.2h)}{w1,w3,w5}
         \fmflabel{$t_1$}{w1}
         \fmflabel{$t_2$}{w3}
	 %%% Workaround for MetaPost 1.504 bug
         \fmfcmd{pair fubara, fubarb, fubarc; fubara = vloc(__v1); fubarb = vloc(__w2); fubarc = vloc(__w4);}
         \fmfi{plain}{fubara...{right}vloc(__w1){right}...vloc(__w2)}
         \fmfi{plain}{fubarb...{right}vloc(__w3){right}...vloc(__w4)}
         \fmfi{dashes}{fubarc...{right}vloc(__w5){right}...vloc(__v2)}
         \fmf{plain,tension=2,label=$s_1$}{v1,p1}
         \fmf{plain}{o1,p1,q1,o4}
         \fmf{plain,tension=0}{q1,o3}
         \fmf{plain,tension=2,label=$s_2$}{w2,p2}
         \fmf{plain}{o6,p2,q2,o9}
         \fmf{plain,tension=0}{q2,o8}
         \fmf{plain,tension=2,label=$s_3$}{w4,p3}
         \fmf{plain}{o12,q3,p3,o15}
         \fmf{plain,tension=0}{q3,o13}
         \fmf{plain,tension=2,label=$s_4$}{v2,p4}
         \fmf{plain}{o17,q4,p4,o20}
         \fmf{plain,tension=0}{q4,o18}
         \fmfdotn{p}{4}
         \fmfdotn{q}{4}
       \end{fmfgraph*}
     \end{center}
     \caption{\label{fig:phasespace}%
       Phasespace parameterization for~$2\to n$ scattering by a sequence
       of cascade decays.}
   \end{figure}
   A general $2\to n$ scattering process can be parameterized by a sequence
   of cascade decays.  The most symmetric representation is a little bit
   redundant and enters each $t$-channel momentum twice.  *)

    type 'a t = ('a * 'a decay * 'a) list

(* \begin{dubious}
     [let topology = map snd] has type [(momentum * 'a) t -> 'a t]
     and can be used to define topological equivalence classes ``up to
     permutations of momenta,'' which are useful for calculating Whizard
     ``groves''\footnote{Not to be confused with gauge invariant classes
     of Feynman diagrams~\cite{Boos/Ohl:groves}.}~\cite{Kilian:WHIZARD}.
   \end{dubious} *)

    let sort cmp = List.map (fun (l, d, r) -> (l, sort_decay cmp d, r))
    let map f = List.map (fun (l, d, r) -> (f l, map_decay f d, f r))
    let eval ft fl fb = List.map (fun (l, d, r) -> (ft l, eval_decay fl fb d, ft r))

(* Find a tree with a defined ordering relation with respect to~$p$ or create
   a new one at the end of the list.  *) 

    let rec embed_in_decays (p, f as pf) = function
      | [] -> [Leaf (p, Some f)]
      | d' :: rest ->
          let p', _ = label d' in
          if M.lesseq p' p || M.less p p' then
            embed_in_decay pf d' :: rest
          else
            d' :: embed_in_decays pf rest

(* \thocwmodulesubsection{Collecting Ingredients} *)

    type 'a unfinished_decays =
        { n : int;
          t_channel : (momentum * 'a option) list;
          decays : (momentum * 'a option) decay list }

    let empty n = { n = n; t_channel = []; decays = [] }

    let insert_in_unfinished_decays (p, f as pf) d =
      if M.Scattering.spacelike p then
        { d with t_channel = (p, Some f) :: d.t_channel }
      else
        { d with decays = embed_in_decays pf d.decays }

    let flip_incoming plist =
      List.map (fun (p', f') -> (M.Scattering.flip_s_channel_in p', f')) plist

    let unfinished_decays_of_momenta n f2 p =
      List.fold_right insert_in_unfinished_decays
        (sort_momenta (flip_incoming ((M.of_ints n [2], f2) :: p))) (empty n)

(* \thocwmodulesubsection{Assembling Ingredients} *)

    let sort3 compare x y z =
      let a = [| x; y; z |] in
      Array.sort compare a;
      (a.(0), a.(1), a.(2))

(* Take advantage of the fact that sorting with [M.compare]
   sorts with \emph{rising} values of [M.rank]: *)

    let allows_momentum_fusion (p, _) (p1, _) (p2, _) =
      let p2', p1', p' = sort3 M.compare p p1 p2 in
      match M.try_fusion p' p1' p2' with
      | Some _ -> true
      | None -> false

    let allows_fusion p1 p2 d = allows_momentum_fusion (label d) p1 p2

    let rec thread_unfinished_decays' p acc tlist dlist =
      match first_pair (allows_fusion p) tlist dlist with
      | None -> (p, acc, tlist, dlist)
      | Some ((t, _ as td), (tlist', dlist')) ->
          thread_unfinished_decays' t (td :: acc) tlist' dlist'

    let thread_unfinished_decays p c =
      match thread_unfinished_decays' p [] c.t_channel c.decays with
      | _, pairs, [], [] -> pairs
      | _ -> failwith "thread_unfinished_decays"

    let rec combine_decays = function 
      | [] -> []
      | ((t, f as tf), d) :: rest ->
          let p, _ = label d in
          begin match M.try_sub t p with
          | Some p' -> (tf, d, (p', f)) :: combine_decays rest
          | None -> (tf, d, (M.sub (M.neg t) p, f)) :: combine_decays rest
          end

    let finalize t = map finalize1 t

    let of_momenta f1 f2 = function
      | (p, _) :: _ as l ->
          let n = M.dim p in
          finalize (combine_decays
                      (thread_unfinished_decays (M.of_ints n [1], Some f1)
                         (unfinished_decays_of_momenta n f2 l)))
      | [] -> []

(* \thocwmodulesubsection{Diagnostics} *)

    let p_to_string p =
      String.concat "" (List.map string_of_int (M.to_ints (M.abs p)))

    let rec to_string1 = function
      | Leaf p -> "(" ^ p_to_string p ^ ")"
      | Branch (_, d1, d2) -> "(" ^ to_string1 d1 ^ to_string1 d2 ^ ")"

    let to_string ps =
      String.concat "/"
        (List.map (fun (p1, d, p2) ->
          p_to_string p1 ^ to_string1 d ^ p_to_string p2) ps)

(* \thocwmodulesubsection{Examples} *)

    let try_thread_unfinished_decays p c =
      thread_unfinished_decays' p [] c.t_channel c.decays

    let try_of_momenta f  = function
      | (p, _) :: _ as l ->
          let n = M.dim p in
          try_thread_unfinished_decays
            (M.of_ints n [1], None) (unfinished_decays_of_momenta n f l)
      | [] -> invalid_arg "try_of_momenta"

  end

(*i 
   module M = Momentum.Lists
   module PS = Phasespace.Make (M)
   open PS
   let u n = List.map (fun p -> (M.of_ints n p, ()))
   let four_t = u 6 [[3;4]; [1;3;4]; [5;6]]
   let four_s = u 6 [[3;4;5;6]; [3;4]; [5;6]]
   let six_mp_1 = u 8 [[3;4]; [1;3;4]; [5;6]; [1;3;4;5;6]; [7;8]]
   let six_mp_2 = u 8 [[3;4]; [1;3;4]; [5;6]; [2;7;8]; [7;8]]
   let f = map (fun (p, ()) -> M.to_ints p)
   let four_t' = f (of_momenta () () four_t)
   let four_s' = f (of_momenta () () four_s)
   let six_mp_1' = f (of_momenta () () six_mp_1)
   let six_mp_2' = f (of_momenta () () six_mp_2)
i*)

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
