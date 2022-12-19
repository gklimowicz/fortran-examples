(* color.ml --

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

module type Test =
  sig
    val suite : OUnit.test
    val suite_long : OUnit.test
  end

(* \thocwmodulesection{Quantum Numbers} *)

type t =
  | Singlet
  | SUN of int
  | AdjSUN of int

let conjugate = function
  | Singlet -> Singlet
  | SUN n -> SUN (-n)
  | AdjSUN n -> AdjSUN n

let compare c1 c2 =
  match c1, c2 with
  | Singlet, Singlet -> 0
  | Singlet, _ -> -1
  | _, Singlet -> 1
  | SUN n, SUN n' -> compare n n'
  | SUN _, AdjSUN _ -> -1
  | AdjSUN _, SUN _ -> 1
  | AdjSUN n, AdjSUN n' -> compare n n'

module type Line =
  sig
    type t
    val conj : t -> t
    val equal : t -> t -> bool
    val to_string : t -> string
  end

module type Cycles =
  sig

    type line
    type t = (line * line) list

(* Contract the graph by connecting lines and return the number of
   cycles together with the contracted graph.
   \begin{dubious}
     The semantics of the contracted graph is not yet 100\%ly fixed.
   \end{dubious} *)
    val contract : t -> int * t

(* The same as [contract], but returns only the number of cycles
   and raises [Open_line] when not all lines are closed. *)
    val count : t -> int
    exception Open_line

    (* Mainly for debugging \ldots *)
    val to_string : t -> string

  end

module Cycles (L : Line) : Cycles with type line = L.t =
  struct

    type line = L.t
    type t = (line * line) list

    exception Open_line

(* NB: The following algorithm for counting the cycles is quadratic since it
   performs nested scans of the lists.  If this was a serious problem one could
   replace the lists of pairs by a [Map] and replace one power by a logarithm. *)

    let rec find_fst c_final c1 disc seen = function
      | [] -> ((L.conj c_final, c1) :: disc, List.rev seen)
      | (c1', c2') as c12' :: rest ->
          if L.equal c1 c1' then
            find_snd c_final (L.conj c2') disc [] (List.rev_append seen rest)
          else
            find_fst c_final c1 disc (c12' :: seen) rest

    and find_snd c_final c2 disc seen = function
      | [] -> ((L.conj c_final, L.conj c2) :: disc, List.rev seen)
      | (c1', c2') as c12' :: rest->
          if L.equal c2' c2 then begin
            if L.equal c1' c_final then
              (disc, List.rev_append seen rest)
            else
              find_fst c_final (L.conj c1') disc [] (List.rev_append seen rest)
          end else
            find_snd c_final c2 disc (c12' :: seen) rest

    let consume = function
      | [] -> ([], [])
      | (c1, c2) :: rest -> find_snd (L.conj c1) (L.conj c2) [] [] rest

    let contract lines =
      let rec contract' acc disc = function
        | [] -> (acc, List.rev disc)
        | rest ->
            begin match consume rest with
            | [], rest' -> contract' (succ acc) disc rest'
            | disc', rest' -> contract' acc (List.rev_append disc' disc) rest'
            end in
      contract' 0 [] lines

    let count lines =
      match contract lines with
      | n, [] -> n
      | n, _ -> raise Open_line

    let to_string lines =
      String.concat ""
        (List.map
           (fun (c1, c2) -> "[" ^ L.to_string c1 ^ "," ^ L.to_string c2 ^ "]")
           lines)

  end

(* \thocwmodulesection{Color Flows} *)

module type Flow =
  sig
    type color
    type t = color list * color list
    val rank : t -> int
    val of_list : int list -> color
    val ghost : unit -> color
    val to_lists : t -> int list list
    val in_to_lists : t -> int list list
    val out_to_lists : t -> int list list
    val ghost_flags : t -> bool list
    val in_ghost_flags : t -> bool list
    val out_ghost_flags : t -> bool list
    type power = { num : int; den : int; power : int }
    type factor = power list
    val factor : t -> t -> factor
    val zero : factor
    module Test : Test
  end

module Flow : Flow = 
  struct

    (* All [int]s are non-zero! *)
    type color =
      | N of int
      | N_bar of int
      | SUN of int * int
      | Singlet
      | Ghost

    (* Incoming and outgoing, since we need to cross the incoming states. *)
    type t = color list * color list

    let rank cflow =
      2

(* \thocwmodulesubsection{Constructors} *)

    let ghost () =
      Ghost

    let of_list = function
      | [0; 0] -> Singlet
      | [c; 0] -> N c
      | [0; c] -> N_bar c
      | [c1; c2] -> SUN (c1, c2)
      | _ -> invalid_arg "Color.Flow.of_list: num_lines != 2"

    let to_list = function
      | N c -> [c; 0]
      | N_bar c -> [0; c]
      | SUN (c1, c2) -> [c1; c2]
      | Singlet -> [0; 0]
      | Ghost -> [0; 0]

    let to_lists (cfin, cfout) =
      (List.map to_list cfin) @ (List.map to_list cfout)

    let in_to_lists (cfin, _) =
      List.map to_list cfin

    let out_to_lists (_, cfout) =
      List.map to_list cfout

    let ghost_flag = function
      | N _ | N_bar _ | SUN (_, _) | Singlet -> false
      | Ghost -> true

    let ghost_flags (cfin, cfout) =
      (List.map ghost_flag cfin) @ (List.map ghost_flag cfout)

    let in_ghost_flags (cfin, _) =
      List.map ghost_flag cfin

    let out_ghost_flags (_, cfout) =
      List.map ghost_flag cfout

(* \thocwmodulesubsection{Evaluation} *)

    type power = { num : int; den : int; power : int }
    type factor = power list
    let zero = []

    let count_ghosts1 colors =
      List.fold_left
        (fun acc -> function Ghost -> succ acc | _ -> acc)
        0 colors

    let count_ghosts (fin, fout) =
      count_ghosts1 fin + count_ghosts1 fout

    type 'a square =
      | Square of 'a
      | Mismatch

    let conjugate = function
      | N c -> N_bar (-c)
      | N_bar c -> N (-c)
      | SUN (c1, c2) -> SUN (-c2, -c1)
      | Singlet -> Singlet
      | Ghost -> Ghost

    let cross_in (cin, cout) =
      cin @ (List.map conjugate cout)

    let cross_out (cin, cout) =
      (List.map conjugate cin) @ cout
      
    module C = Cycles (struct
      type t = int
      let conj = (~-)
      let equal = (=)
      let to_string = string_of_int
    end)

(* Match lines in the color flows [f1] and [f2] after crossing the
   incoming states.  This will be used to compute squared diagrams
   in [square] and [square2] below. *)

    let match_lines match1 match2 f1 f2 =
      let rec match_lines' acc f1' f2' =
        match f1', f2' with

        (* If we encounter an empty list, we're done --- unless the
           lengths don't match (which should never happen!): *)
        | [], [] -> Square (List.rev acc)
        | _ :: _, [] | [], _ :: _ -> Mismatch

        (* Handle matching \ldots *)
        | Ghost :: rest1, Ghost :: rest2
        | Singlet :: rest1, Singlet :: rest2 ->
           match_lines' acc rest1 rest2

        (* \ldots{} and mismatched ghosts and singlet gluons: *)
        | Ghost :: _, Singlet :: _
        | Singlet :: _, Ghost :: _ ->
           Mismatch

        (* Ghosts and singlet gluons can't match anything else *)
        | (Ghost | Singlet) :: _, (N _ | N_bar _ | SUN (_, _)) :: _
        | (N _ | N_bar _ | SUN (_, _)) :: _, (Ghost | Singlet) :: _ ->
           Mismatch

        (* Handle matching \ldots *)
        | N_bar c1 :: rest1, N_bar c2 :: rest2
        | N c1 :: rest1, N c2 :: rest2 ->
           match_lines' (match1 c1 c2 acc) rest1 rest2

        (* \ldots{} and mismatched $N$ or $\bar N$ states: *)
        | N _ :: _, N_bar _ :: _
        | N_bar _ :: _, N _ :: _ ->
           Mismatch

        (* The $N$ and $\bar N$ don't match non-singlet gluons: *)
        | (N _ | N_bar _) :: _, SUN (_, _) :: _
        | SUN (_, _) :: _, (N _ | N_bar _) :: _ ->
           Mismatch

        (* Now we're down to non-singlet gluons: *)
        | SUN (c1, c1') :: rest1, SUN (c2, c2') :: rest2 ->
           match_lines' (match2 c1 c1' c2 c2' acc) rest1 rest2 in

      match_lines' [] (cross_out f1) (cross_out f2)

(* NB: in WHIZARD versions before 3.0, the code for [match_lines]
   contained a bug in the pattern matching of [Singlet], [N], [N_bar]
   and [SUN] states, because they all were represented as
   [SUN (c1, c2)], only distinguished by the numeric conditions
   [c1 = 0] and/or [c2 = 0].
   This prevented the use of exhaustiveness checking and introduced a
   subtle dependence on the pattern order. *)

    let square f1 f2 =
      match_lines
        (fun c1 c2 pairs -> (c1, c2) :: pairs)
        (fun c1 c1' c2 c2' pairs -> (c1', c2') :: (c1, c2) :: pairs)
        f1 f2

(*i
    let square f1 f2 =
      let ll2s ll =
        String.concat "; "
          (List.map (ThoList.to_string string_of_int) ll)
      and lp2s lp =
        String.concat "; "
          (List.map
             (fun (c1, c2) ->
               string_of_int c1 ^ ", " ^ string_of_int c2)
             lp) in
      Printf.eprintf
        "square ([%s], [%s]) ([%s], [%s]) = "
        (ll2s (in_to_lists f1)) (ll2s (out_to_lists f1))
        (ll2s (in_to_lists f2)) (ll2s (out_to_lists f2));
      let res = square f1 f2 in
      begin match res with
      | Mismatch -> Printf.eprintf "Mismatch!\n"
      | Square f12 -> Printf.eprintf "Square [%s]\n" (lp2s f12)
      end;
      res
i*)

(* In addition to counting closed color loops, we also need to count closed
   gluon loops.  Fortunately, we can use the same algorithm on a different
   data type, provided it doesn't require all lines to be closed. *)

    module C2 = Cycles (struct
      type t = int * int
      let conj (c1, c2) = (- c2, - c1)
      let equal (c1, c2) (c1', c2') = c1 = c1' && c2 = c2'
      let to_string (c1, c2) = "(" ^ string_of_int c1 ^ "," ^ string_of_int c2 ^ ")"
    end)

    let square2 f1 f2 =
      match_lines
        (fun c1 c2 pairs -> pairs)
        (fun c1 c1' c2 c2' pairs -> ((c1, c1'), (c2, c2')) :: pairs)
        f1 f2

(* $\ocwlowerid{int\_power}: n\, p \to n^p$
   for integers is missing from [Pervasives]! *)

    let int_power n p =
      let rec int_power' acc i =
        if i < 0 then
          invalid_arg "int_power"
        else if i = 0 then
          acc
        else
          int_power' (n * acc) (pred i) in
      int_power' 1 p

(* Instead of implementing a full fledged algebraic evaluator, let's
   simply expand the binomial by hand:
   \begin{equation}
    \left(\frac{N_C^2-2}{N_C^2}\right)^n =
      \sum_{i=0}^n \binom{n}{i} (-2)^i N_C^{-2i}
   \end{equation} *)

(* NB: Any result of [square] other than [Mismatch] guarantees
   [count_ghosts f1 = count_ghosts f2]. *)

    let factor f1 f2 =
      match square f1 f2, square2 f1 f2 with
      | Mismatch, _ | _, Mismatch -> []
      | Square f12, Square f12' ->
          let num_cycles = C.count f12
          and num_cycles2, disc = C2.contract f12'
          and num_ghosts = count_ghosts f1 in
(*i       Printf.eprintf "f12  = %s -> #loops = %d\n"
            (C.to_string f12) num_cycles;
          Printf.eprintf "f12' = %s -> #loops = %d, disc = %s\n"
            (C2.to_string f12') num_cycles2 (C2.to_string disc);
          flush stderr; i*)
          List.map
            (fun i ->
              let parity = if num_ghosts mod 2 = 0 then 1 else -1
              and power = num_cycles - num_ghosts in
              let coeff = int_power (-2) i * Combinatorics.binomial num_cycles2 i
              and power2 = - 2 * i in
              { num = parity * coeff;
                den = 1;
                power = power + power2 })
            (ThoList.range 0 num_cycles2)

    module Test : Test =
      struct

        open OUnit

(* Here and elsewhere, we have to resist the temptation to define
   these tests as functions with an additional argument [()] in the
   hope to avoid having to package them into an explicit thunk
   [fun () -> eq v1 v2] in order to delay
   evaluation. It turns out that the runtime would then sometimes
   evaluate the argument [v1] or [v2] even \emph{before} the test
   is run.  For pure functions, there is no difference, but the
   compiler appears to treat explicit thunks specially.
   \begin{dubious}
     I haven't yet managed to construct a small demonstrator to find
     out in which circumstances the premature evaluation happens.
   \end{dubious} *)

        let suite_square =
          "square" >:::

            [ "square ([], []) ([], [])" >::
                (fun () ->
	          assert_equal (Square []) (square ([], []) ([], [])));

              "square ([3], [3; 0]) ([3], [3; 0])" >::
                (fun () ->
	          assert_equal
                    (Square [(-1, -1); (1, 1)])
                    (square
                       ([N 1], [N 1; Singlet])
                       ([N 1], [N 1; Singlet])));

              "square ([0], [3; -3]) ([0], [3; -3])" >::
                (fun () ->
	          assert_equal
                    (Square [(1, 1); (-1, -1)])
                    (square
                       ([Singlet], [N 1; N_bar (-1)])
                       ([Singlet], [N 1; N_bar (-1)])));
 
              "square ([3], [3; 0]) ([0], [3; -3])" >::
                (fun () ->
	          assert_equal
                    Mismatch
                    (square
                       ([N 1], [N 1; Singlet])
                       ([Singlet], [N 1; N_bar (-1)])));

              "square ([3; 8], [3]) ([3; 8], [3])" >::
                (fun () ->
	          assert_equal
                    (Square [-1, -1; 1, 1; -2, -2; 2, 2])
                    (square
                       ([N 1; SUN (2, -1)], [N 2])
                       ([N 1; SUN (2, -1)], [N 2]))) ]

        let suite =
          "Color.Flow" >:::
	    [suite_square]

        let suite_long =
          "Color.Flow long" >:::
	    []

      end
  end

(* later: *)

module General_Flow = 
  struct

    type color =
      | Lines of int list
      | Ghost of int

    type t = color list * color list

    let rank_default = 2 (* Standard model *)

    let rank cflow =
      try
        begin match List.hd cflow with
        | Lines lines -> List.length lines
        | Ghost n_lines -> n_lines
        end
      with
      | _ -> rank_default
  end

(* \thocwmodulesection{Vertex Color Flows} *)

(* \newcommand{\setupFourAmp}{%
     \fmfleft{i1,i2}
     \fmfright{o1,o2}
     \fmf{phantom}{i1,v1,i2}
     \fmf{phantom}{o2,v2,o1}
     \fmf{phantom}{v1,v2}
     \fmffreeze}
   \fmfcmd{%
     numeric joindiameter;
     joindiameter := 7thick;}
   \fmfcmd{%
     vardef sideways_at (expr d, p, frac) =
       save len; len = length p;
       (point frac*len of p) shifted ((d,0) rotated (90 + angle direction frac*len of p))
     enddef;
     secondarydef p sideways d =
       for frac = 0 step 0.01 until 0.99:
         sideways_at (d, p, frac) ..
       endfor
       sideways_at (d, p, 1)
     enddef;
     secondarydef p choptail d =
      subpath (ypart (fullcircle scaled d shifted (point 0 of p) intersectiontimes p), infinity) of p
     enddef;
     secondarydef p choptip d =
      reverse ((reverse p) choptail d)
     enddef;
     secondarydef p pointtail d =
       fullcircle scaled d shifted (point 0 of p) intersectionpoint p
     enddef;
     secondarydef p pointtip d =
       (reverse p) pointtail d
     enddef;
     secondarydef pa join pb =
       pa choptip joindiameter .. pb choptail joindiameter
     enddef;
     vardef cyclejoin (expr p) =
       subpath (0.5*length p, infinity) of p join subpath (0, 0.5*length p) of p .. cycle
     enddef;}
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \fmfcmd{%
     style_def double_line_arrow expr p =
       save pi, po; 
       path pi, po;
       pi = reverse (p sideways thick);
       po = p sideways -thick;
       cdraw pi;
       cdraw po;
       cfill (arrow (subpath (0, 0.9 length pi) of pi));
       cfill (arrow (subpath (0, 0.9 length po) of po));
     enddef;}
   \fmfcmd{%
     style_def double_line_arrow_beg expr p =
       save pi, po, pc; 
       path pi, po, pc;
       pc = p choptail 7thick;
       pi = reverse (pc sideways thick);
       po = pc sideways -thick;
       cdraw pi .. p pointtail 5thick .. po;
       cfill (arrow pi);
       cfill (arrow po);
     enddef;}
   \fmfcmd{%
     style_def double_line_arrow_end expr p =
       save pi, po, pc; 
       path pi, po, pc;
       pc = p choptip 7thick;
       pi = reverse (pc sideways thick);
       po = pc sideways -thick;
       cdraw po .. p pointtip 5thick .. pi;
       cfill (arrow pi);
       cfill (arrow po);
     enddef;}
   \fmfcmd{%
     style_def double_line_arrow_both expr p =
       save pi, po, pc; 
       path pi, po, pc;
       pc = p choptip 7thick choptail 7thick;
       pi = reverse (pc sideways thick);
       po = pc sideways -thick;
       cdraw po .. p pointtip 5thick .. pi .. p pointtail 5thick .. cycle;
       cfill (arrow pi);
       cfill (arrow po);
     enddef;}
   \fmfcmd{%
     style_def double_arrow_parallel expr p =
       save pi, po; 
       path pi, po;
       pi = p sideways thick;
       po = p sideways -thick;
       save li, lo;
       li = length pi;
       lo = length po;
       cdraw pi;
       cdraw po;
       cfill (arrow pi);
       cfill (arrow po);
     enddef;}
   \fmfcmd{%
     style_def double_arrow_crossed_beg expr p =
       save lp;  lp = length p;
       save pi, po; 
       path pi, po;
       pi = p sideways thick;
       po = p sideways -thick;
       save li, lo;
       li = length pi;
       lo = length po;
       cdraw subpath (0, 0.1 li) of pi .. subpath (0.3 lo, lo) of po;
       cdraw subpath (0, 0.1 lo) of po .. subpath (0.3 li, li) of pi;
       cfill (arrow pi);
       cfill (arrow po);
     enddef;}
   \fmfcmd{%
     style_def double_arrow_crossed_end expr p =
       save lp;  lp = length p;
       save pi, po; 
       path pi, po;
       pi = p sideways thick;
       po = p sideways -thick;
       save li, lo;
       li = length pi;
       lo = length po;
       cdraw subpath (0, 0.7 li) of pi .. subpath (0.9 lo, lo) of po;
       cdraw subpath (0, 0.7 lo) of po .. subpath (0.9 li, li) of pi;
       cfill (arrow pi);
       cfill (arrow po);
     enddef;} *)

module Q = Algebra.Q
module QC = Algebra.QC

module type Arrow =
  sig
    type endpoint
    type tip = endpoint
    type tail = endpoint
    type ghost = endpoint
    val position : endpoint -> int
    val relocate : (int -> int) -> endpoint -> endpoint
    type ('tail, 'tip, 'ghost) t =
      | Arrow of 'tail * 'tip
      | Ghost of 'ghost
      | Epsilon of 'tip list
      | Epsilon_bar of 'tail list
    type free = (tail, tip, ghost) t
    type factor
    val free_to_string : free -> string
    val factor_to_string : factor -> string
    val map : (endpoint -> endpoint) -> free -> free
    val to_left_factor : (endpoint -> bool) -> free -> factor
    val to_right_factor : (endpoint -> bool) -> free -> factor
    val of_factor : factor -> free
    val is_free : factor -> bool
    val negatives : free -> endpoint list
    val is_ghost : free -> bool
    val is_tadpole : factor -> bool
    val is_epsilon : factor -> bool
    val match_epsilon : factor -> factor list -> (factor list list * factor list list * factor list) option
    type merge =
      | Match of factor
      | Ghost_Match
      | Loop_Match
      | Mismatch
      | No_Match
    val merge : factor -> factor -> merge
    val tee : int -> free -> free list
    val dir : int -> int -> free -> int
    val single : endpoint -> endpoint -> free
    val double : endpoint -> endpoint -> free list
    val ghost : endpoint -> free
    module Infix : sig
      val (=>) : int -> int -> free
      val (==>) : int -> int -> free list
      val (<=>) : int -> int -> free list
      val (>=>) : int * int -> int -> free
      val (=>>) : int -> int * int -> free
      val (>=>>) : int * int -> int * int -> free
      val (??) : int -> free
    end
    val epsilon : int list -> free
    val epsilon_bar : int list -> free
    val chain : int list -> free list
    val cycle : int list -> free list
    module Test : Test
    val pp_free : Format.formatter -> free -> unit
    val pp_factor : Format.formatter -> factor -> unit
  end

module Arrow : Arrow =
  struct

    type endpoint =
      | I of int
      | M of int * int

    let position = function
      | I i -> i
      | M (i, _) -> i

    let relocate f = function
      | I i -> I (f i)
      | M (i, n) -> M (f i, n)

    type tip = endpoint
    type tail = endpoint
    type ghost = endpoint

    (* Note that in the case of double lines for the adjoint
       representation the \emph{same} [endpoint] appears twice:
       once as a [tip] and once as a [tail].  If we want to
       multiply two factors by merging arrows with matching
       [tip] and [tail], we must make sure that the [tip] is from
       one factor and the [tail] from the other factor. *)
               
    (* The [Free] variant contains positive indices
       as well as negative indices that don't appear on both sides
       and will be summed in a later product.  [SumL] and [SumR]
       indices appear on both sides. *)
    type 'a index =
      | Free of 'a
      | SumL of 'a
      | SumR of 'a

    let is_free_index = function
      | Free _ -> true
      | SumL _ | SumR _ -> false

    type ('tail, 'tip, 'ghost) t =
      | Arrow of 'tail * 'tip
      | Ghost of 'ghost
      | Epsilon of 'tip list
      | Epsilon_bar of 'tail list

    type free = (tail, tip, ghost) t
    type factor = (tail index, tip index, ghost index) t

    let endpoint_to_string = function
      | I i -> string_of_int i
      | M (i, n) -> Printf.sprintf "%d.%d" i n

    let index_to_string = function
      | Free i -> endpoint_to_string i
      | SumL i -> endpoint_to_string i ^ "L"
      | SumR i -> endpoint_to_string i ^ "R"

    let to_string i2s = function
      | Arrow (tail, tip) -> Printf.sprintf "%s>%s" (i2s tail) (i2s tip)
      | Ghost ghost -> Printf.sprintf "{%s}" (i2s ghost)
      | Epsilon tips -> Printf.sprintf ">>>%s" (ThoList.to_string i2s tips)
      | Epsilon_bar tails -> Printf.sprintf "<<<%s" (ThoList.to_string i2s tails)

    let free_to_string = to_string endpoint_to_string

    let factor_to_string = to_string index_to_string

    let index_matches i1 i2 =
      match i1, i2 with
      | SumL i1, SumR i2 | SumR i1, SumL i2 -> i1 = i2
      | _ -> false

    let map f = function
      | Arrow (tail, tip) -> Arrow (f tail, f tip)
      | Ghost ghost -> Ghost (f ghost)
      | Epsilon tips -> Epsilon (List.map f tips)
      | Epsilon_bar tails -> Epsilon_bar (List.map f tails)

    let free_index = function
      | Free i -> i
      | SumL i -> invalid_arg "Color.Arrow.free_index: leftover LHS summation"
      | SumR i -> invalid_arg "Color.Arrow.free_index: leftover RHS summation"

    let to_left_index is_sum i =
      if is_sum i then
        SumL i
      else
        Free i

    let to_right_index is_sum i =
      if is_sum i then
        SumR i
      else
        Free i

    let to_left_factor is_sum = map (to_left_index is_sum)
    let to_right_factor is_sum = map (to_right_index is_sum)
    let of_factor = map free_index

    let negatives = function
      | Arrow (tail, tip) ->
         if position tail < 0 then
           if position tip < 0 then
             [tail; tip]
           else
             [tail]
         else if position tip < 0 then
           [tip]
         else
           []
      | Ghost ghost ->
         if position ghost < 0 then
           [ghost]
         else
           []
      | Epsilon tips -> List.filter (fun tip -> position tip < 0) tips
      | Epsilon_bar tails -> List.filter (fun tail -> position tail < 0) tails

    let is_free = function
      | Arrow (Free _, Free _) | Ghost (Free _) -> true
      | Arrow (_, _) | Ghost _ -> false
      | Epsilon tips -> List.for_all is_free_index tips
      | Epsilon_bar tails -> List.for_all is_free_index tails

    let is_ghost = function
      | Ghost _ -> true
      | Arrow _ -> false
      | Epsilon _ | Epsilon_bar _ -> false
                 
    let is_epsilon = function
      | Epsilon _ | Epsilon_bar _ -> true
      | Ghost _ | Arrow _ -> false

    let single tail tip =
      Arrow (tail, tip)

    let double a b =
      if a = b then
        [single a b]
      else
        [single a b; single b a]

    let ghost g =
      Ghost g

    module Infix =
      struct
        let ( => ) i j = single (I i) (I j)
        let ( ==> ) i j = [i => j]
        let ( <=> ) i j = double (I i) (I j)
        let ( >=> ) (i, n) j = single (M (i, n)) (I j)
        let ( =>> ) i (j, m) = single (I i) (M (j, m))
        let ( >=>> ) (i, n) (j, m) = single (M (i, n)) (M (j, m))
        let ( ?? ) i = ghost (I i)
      end

    open Infix

(* Split [a_list] at the first element equal to [a] according
   to [eq].  Return the reversed first part and the rest as a
   pair and wrap it in [Some]. Return [None] if there is no match.  *)
    let take_first_match_opt ?(eq=(=)) a a_list =
      let rec take_first_match_opt' rev_head = function
        | [] -> None
        | elt :: tail ->
           if eq elt a then
             Some (rev_head, tail)
           else
             take_first_match_opt' (elt :: rev_head) tail in
      take_first_match_opt' [] a_list

(* Split [a_list] and [b_list] at the first element equal according
   to [eq].  Return the reversed first part and the rest of each
   as a pair of pairs wrap it in [Some].
   Return [None] if there is no match.  *)
    let take_first_matching_pair_opt ?(eq=(=)) a_list b_list =
      let rec take_first_matching_pair_opt' rev_a_head = function
        | [] -> None
        | a :: a_tail ->
           begin match take_first_match_opt ~eq a b_list with
           | Some (rev_b_head, b_tail) ->
              Some ((rev_a_head, a_tail), (rev_b_head, b_tail))
           | None ->
              take_first_matching_pair_opt' (a :: rev_a_head) a_tail
           end in
      take_first_matching_pair_opt' [] a_list

(* Replace the first occurence of an element equal to [a] according
   to [eq] in [a_list] by [a'] and wrap the new list in [Some].
   Return [None] if there is no match.  *)
    let replace_first_opt ?(eq=(=)) a a' a_list =
      match take_first_match_opt ~eq a a_list with
      | Some (rev_head, tail) -> Some (List.rev_append rev_head (a' :: tail))
      | None -> None

    let tee a = function
      | Arrow (tail, tip) -> [Arrow (tail, I a); Arrow (I a, tip)]
      | Ghost _ as g -> [g]
      | Epsilon _ -> invalid_arg "Arrow.tee not defined for Epsilon"
      | Epsilon_bar _ -> invalid_arg "Arrow.tee not defined for Epsilon_bar"

    let dir i j = function
      | Arrow (tail, tip) ->
         let tail = position tail
         and tip = position tip in
         if tip = i && tail = j then
            1
         else if tip = j && tail = i then
           -1
         else
           invalid_arg "Arrow.dir"
      | Ghost _ | Epsilon _ | Epsilon_bar _ -> 0

    type merge =
      | Match of factor
      | Ghost_Match
      | Loop_Match
      | Mismatch
      | No_Match

(* When computing
    \begin{equation}
    \label{eq:epsilon*epsilonbar-single-contraction-N=n}
      \epsilon_{ki_1i_2\cdots i_n} \bar\epsilon^{kj_1j_2\cdots j_n}
        = \sum_{\sigma\in S_n} (-1)^{\epsilon(\sigma)}
            \delta_{i_1}^{\sigma(j_1)} 
            \delta_{i_2}^{\sigma(j_2)} 
            \cdots
            \delta_{i_n}^{\sigma(j_n)}\,,
    \end{equation}
    we must keep track of the position of summation indices.
    We can use the fact that cyclic permutations are even for
    $\epsilon$-tensors with an odd number of indices, corresponding
    to $n$ even and odd otherwise. *)

    let fuse_epsilons1 tails tips =
      match take_first_matching_pair_opt ~eq:index_matches tails tips with
      | None -> None
      | Some ((rev_tails_head, tails_tail), (rev_tips_head, tips_tail)) ->
         let tails = tails_tail @ List.rev rev_tails_head
         and tips = tips_tail @ List.rev rev_tips_head  in
         let num_tails = List.length tails
         and num_tips = List.length tips in
         if num_tails <> num_tips then
           invalid_arg
             (Printf.sprintf
                "Color.Arrow.fuse_epsilons1: length mismatch %d <> %d"
                (succ num_tails) (succ num_tips))
         else
           let is_odd n = n mod 2 <> 0 in
           let flip =
             is_odd num_tips &&
               is_odd (List.length rev_tails_head - List.length rev_tips_head) in
           let even_tips = Combinatorics.permute_even tips
           and odd_tips = Combinatorics.permute_odd tips in
           let even = List.rev_map (List.rev_map2 single tails) even_tips
           and odd = List.rev_map (List.rev_map2 single tails) odd_tips in
           if flip then
             Some (odd, even)
           else
             Some (even, odd)

(* We can also use the following (slightly less efficient)
   version that does not need to keep track of signs by itself and
   is more general, since it does not depend on~$N=n$ and works
   for~$N\ge n$.

   Starting with the case of matching dimension~$N$ and rank
   of~$\epsilon$ and $\bar\epsilon$, there is the well known formula
   \begin{equation}
      \forall k, n = N \in \mathbf{N}, 0\le k \le n \ge 2:\;
      \epsilon_{i_1\cdots i_n}
      \bar\epsilon^{i_1\cdots i_kj_{k+1}\cdots j_n}
        = k! \sum_{\sigma\in S_{n-k}} (-1)^{\varepsilon(\sigma)}
            \delta_{i_{k+1}}^{\sigma(j_{k+1})} 
            \delta_{i_{k+2}}^{\sigma(j_{k+2})} 
            \cdots
            \delta_{i_n}^{\sigma(j_n)}\,.
   \end{equation}
   In the general case, we have from anti-symmetry alone
   \begin{equation}
      \forall n, N \in\mathbf{N}, 2\le n \le N:\;
      \epsilon_{i_1i_2\cdots i_n} \bar\epsilon^{j_1j_2\cdots j_n}
        = \sum_{\sigma\in S_n} (-1)^{\varepsilon(\sigma)}
            \delta_{i_1}^{\sigma(j_1)} 
            \delta_{i_2}^{\sigma(j_2)} 
            \cdots
            \delta_{i_n}^{\sigma(j_n)}\,,
   \end{equation}
   where~$N=\delta_i^i$ is the dimension. *)

    let fuse_epsilons tails tips =
      if List.length tails = List.length tips then
        List.fold_left
          (fun (even, odd) (eps, tips) ->
            if eps > 0 then
              (List.rev_map2 single tails tips :: even, odd)
            else if eps < 0 then
              (even, List.rev_map2 single tails tips :: odd)
            else
              failwith "Color.Arrow.fuse_epsilons: Combinatorics.permute_signed returned garbage")
          ([], []) (Combinatorics.permute_signed tips)
      else
        invalid_arg "Color.Arrow.fuse_epsilons: |tails| <> |tips|"

(* From this, we derive
   \begin{multline}
   \label{eq:epsilon*epsilonbar-single-contraction}
      \forall n, N \in\mathbf{N}, 2\le n \le N:\;
      \epsilon_{ki_2\cdots i_n} \bar\epsilon^{kj_2\cdots j_n}
        = \sum_{\sigma\in S_n} (-1)^{\varepsilon(\sigma)}
            \delta_{k}^{\sigma(k)} 
            \delta_{i_2}^{\sigma(j_2)} 
            \cdots
            \delta_{i_n}^{\sigma(j_n)} \\
        = (N-n+1)
            \sum_{\sigma\in S_{n-1}} (-1)^{\varepsilon(\sigma)}
            \delta_{i_2}^{\sigma(j_2)} 
            \cdots
            \delta_{i_n}^{\sigma(j_n)}\,,
   \end{multline}
   where the~$N=\delta_k^k$ comes from the permutations with~$\sigma(k)=k$
   that correspond to a loop in the color flow and the~$n-1$ from the
   permutations with~$\sigma(k)\in\{i_2,\ldots,i_n\}$ that do not
   lead to a loop.  Note that~$N-n+1=1$ in the special case~$N=n$ when
   rank and dimension match.

   By induction
   \begin{multline}
      \forall k, n, N \in \mathbf{N}, 2\le n \le N \land 1\le k \le n:\;
      \epsilon_{i_1\cdots i_n}
      \bar\epsilon^{i_1\cdots i_kj_{k+1}\cdots j_n}\\
        = \frac{(N-n+k)!}{(N-n)!}
            \sum_{\sigma\in S_{n-k}} (-1)^{\varepsilon(\sigma)}
            \delta_{i_{k+1}}^{\sigma(j_{k+1})} 
            \delta_{i_{k+2}}^{\sigma(j_{k+2})} 
            \cdots
            \delta_{i_n}^{\sigma(j_n)}\,,
   \end{multline}
   where
   \begin{equation}
     \frac{(N-n+k)!}{(N-n)!} = (N-n+1)(N-n+2)\cdots(N-n+k)
   \end{equation}
   and in the special case~$N=n$
   \begin{equation}
     \frac{(N-n+k)!}{(N-n)!} = k!\,.
   \end{equation}
   In the case~$k=1$ we
   get~\eqref{eq:epsilon*epsilonbar-single-contraction}, which reduces
   to~\eqref{eq:epsilon*epsilonbar-single-contraction-N=n} for $N=n$,
   of course. *)
 
(* \begin{dubious}
      We also need to handle disconnected pairs of~$\epsilon$
      and~$\bar\epsilon$.  These never appear in [merge']
      below, because [merge] skips all cases without matching
      summation indices.  Handling them in [merge] doesn't work
      yet.
   \end{dubious} *)

    let merge' arrow1 arrow2 =
      match arrow1, arrow2 with
      | Ghost g1, Ghost g2 ->
         if index_matches g1 g2 then
           Ghost_Match
         else
           No_Match
      | Arrow (tail, tip), Ghost g
      | Ghost g, Arrow (tail, tip) ->
         if index_matches g tail || index_matches g tip then
           Mismatch
         else
           No_Match
      | Arrow (tail, tip), Arrow (tail', tip') ->
         if index_matches tip tail' then
           if index_matches tip' tail then
             Loop_Match
           else
             Match (Arrow (tail, tip'))
         else if index_matches tip' tail then
           Match (Arrow (tail', tip))
         else
           No_Match
      | Arrow (tail, tip), Epsilon tips | Epsilon tips, Arrow (tail, tip) ->
         begin match replace_first_opt ~eq:index_matches tail tip tips with
         | None -> No_Match
         | Some tips -> Match (Epsilon tips)
         end
      | Arrow (tail, tip), Epsilon_bar tails | Epsilon_bar tails, Arrow (tail, tip) ->
         begin match replace_first_opt ~eq:index_matches tip tail tails with
         | None -> No_Match
         | Some tails -> Match (Epsilon_bar tails)
         end
      | Epsilon tips, Ghost g | Ghost g, Epsilon tips ->
         if List.exists (index_matches g) tips then
           Mismatch
         else
           No_Match
      | Epsilon_bar tails, Ghost g | Ghost g, Epsilon_bar tails ->
         if List.exists (index_matches g) tails then
           Mismatch
         else
           No_Match
      | Epsilon _, Epsilon _ | Epsilon_bar _, Epsilon_bar _ ->
         No_Match
      | Epsilon tips, Epsilon_bar tails | Epsilon_bar tails, Epsilon tips ->
         failwith "Color.Arrow.merge': impossible"

(* As an optimization, don't attempt to merge if neither of the arrows
   contains a summation index and return immediately. *)

    let merge arrow1 arrow2 =
      if is_free arrow1 || is_free arrow2 then
        No_Match
      else
        merge' arrow1 arrow2

    let merge_to_string = function
      | Match factor ->
         Printf.sprintf "Match (%s)" (factor_to_string factor)
      | Ghost_Match -> "Ghost"
      | Loop_Match -> "Loop"
      | Mismatch -> "Mismatch"
      | No_Match -> "No_Match"

    let logging_merge arrow1 arrow2 =
      let result = merge arrow1 arrow2 in
      Printf.eprintf
        "merge %s with %s ==> %s\n"
        (factor_to_string arrow1)
        (factor_to_string arrow2)
        (merge_to_string result);
      result

    let is_tadpole = function
      | Arrow (tail, tip) ->
         index_matches tail tip
      | _ -> false

    let merge_epsilon_pair arrow1 arrow2 =
      match arrow1, arrow2 with
      | Epsilon tips, Epsilon_bar tails | Epsilon_bar tails, Epsilon tips ->
         Some (fuse_epsilons tails tips)
      | _ -> None

    let match_epsilon arrow arrows =
      let rec match_epsilon' seen = function
        | [] -> None
        | arrow' :: arrows' ->
           begin match merge_epsilon_pair arrow arrow' with
           | Some (even, odd) ->
              Some (even, odd, List.rev_append seen arrows')
           | None -> match_epsilon' (arrow' :: seen) arrows'
           end in
      match_epsilon' [] arrows

(*i
    let merge = logging_merge
i*)
    let epsilon = function
      | [] -> invalid_arg "Color.Arrow.epsilon []"
      | [_] -> invalid_arg "Color.Arrow.epsilon lone index"
      | tips ->
         Epsilon (List.map (fun tip -> I tip) tips)

    let epsilon_bar = function
      | [] -> invalid_arg "Color.Arrow.epsilon []"
      | [_] -> invalid_arg "Color.Arrow.epsilon lone index"
      | tails ->
         Epsilon_bar (List.map (fun tail -> I tail) tails)

    (* Composite Arrows. *)

    let rec chain = function
      | [] -> []
      | [a] -> [a => a]
      | [a; b] -> [a => b]
      | a :: (b :: _ as rest) -> (a => b) :: chain rest

    let rec cycle' a = function
      | [] -> [a => a]
      | [b] -> [b => a]
      | b :: (c :: _ as rest) -> (b => c) :: cycle' a rest

    let cycle = function
      | [] -> []
      | a :: _ as a_list -> cycle' a a_list

    module Test : Test =
      struct

        open OUnit

        let suite_chain =
          "chain" >:::
            [ "[]" >:: (fun () -> assert_equal [] (chain []));
              "[1]" >:: (fun () -> assert_equal [1 => 1] (chain [1]));
              "[1;2]" >:: (fun () -> assert_equal [1 => 2] (chain [1; 2]));
              "[1;2;3]" >:: (fun () -> assert_equal [1 => 2; 2 => 3] (chain [1; 2; 3]));
              "[1;2;3;4]" >:: (fun () -> assert_equal [1 => 2; 2 => 3; 3 => 4] (chain [1; 2; 3; 4])) ]

        let suite_cycle =
          "cycle" >:::
            [ "[]" >:: (fun () -> assert_equal [] (cycle []));
              "[1]" >:: (fun () -> assert_equal [1 => 1] (cycle [1]));
              "[1;2]" >:: (fun () -> assert_equal [1 => 2; 2 => 1] (cycle [1; 2]));
              "[1;2;3]" >:: (fun () -> assert_equal [1 => 2; 2 => 3; 3 => 1] (cycle [1; 2; 3]));

              "[1;2;3;4]" >:: (fun () -> assert_equal [1 => 2; 2 => 3; 3 => 4; 4 => 1] (cycle [1; 2; 3; 4])) ]

        let suite_take =
          "take" >:::
            [ "1 []" >:: (fun () -> assert_equal None (take_first_match_opt 1 []));
              "1 [1]" >:: (fun () -> assert_equal (Some ([], [])) (take_first_match_opt 1 [1]));
              "1 [2;3;4]" >:: (fun () -> assert_equal None (take_first_match_opt 1 [2;3;4]));
              "1 [1;2;3]" >:: (fun () -> assert_equal (Some ([], [2;3])) (take_first_match_opt 1 [1;2;3]));
              "2 [1;2;3]" >:: (fun () -> assert_equal (Some ([1], [3])) (take_first_match_opt 2 [1;2;3]));
              "3 [1;2;3]" >:: (fun () -> assert_equal (Some ([2;1], [])) (take_first_match_opt 3 [1;2;3])) ]

        let suite_take2 =
          "take2" >:::
            [ "[] []" >::
	        (fun () -> assert_equal None (take_first_matching_pair_opt [] []));

              "[] [1;2;3]" >::
	        (fun () -> assert_equal None (take_first_matching_pair_opt [] [1;2;3]));

              "[1] [2;3;4]" >::
	        (fun () -> assert_equal None (take_first_matching_pair_opt [1] [2;3;4]));

              "[2;3;4] [1]" >::
	        (fun () -> assert_equal None (take_first_matching_pair_opt [2;3;4] [1]));

              "[1;2;3] [4;5;6;7]" >::
	        (fun () -> assert_equal None (take_first_matching_pair_opt [1;2;3] [4;5;6;7]));

              "[1] [1;2;3]" >::
	        (fun () ->
                  assert_equal
                    (Some (([],[]), ([],[2;3])))
                    (take_first_matching_pair_opt [1] [1;2;3]));

              "[1;2;3] [1;20;30]" >::
	        (fun () ->
                  assert_equal
                    (Some (([],[2;3]), ([],[20;30])))
                    (take_first_matching_pair_opt [1;2;3] [1;20;30]));

              "[1;2;3;4;5;6] [10;20;4;30;40]" >::
	        (fun () ->
                  assert_equal
                    (Some (([3;2;1],[5;6]), ([20;10],[30;40])))
                    (take_first_matching_pair_opt [1;2;3;4;5;6] [10;20;4;30;40])) ]

        let suite_replace =
          "replace" >:::
            [ "1 10 []" >:: (fun () -> assert_equal None (replace_first_opt 1 2 []));
              "1 10 [1]" >:: (fun () -> assert_equal (Some [10]) (replace_first_opt 1 10 [1]));
              "1 [2;3;4]" >:: (fun () -> assert_equal None (replace_first_opt 1 10 [2;3;4]));
              "1 [1;2;3]" >:: (fun () -> assert_equal (Some [10;2;3]) (replace_first_opt 1 10 [1;2;3]));
              "2 [1;2;3]" >:: (fun () -> assert_equal (Some [1;10;3]) (replace_first_opt 2 10 [1;2;3]));
              "3 [1;2;3]" >:: (fun () -> assert_equal (Some [1;2;10]) (replace_first_opt 3 10 [1;2;3])) ]

(*i     let determinant_to_string = function
          | Determinant (even, odd) ->
             Printf.sprintf
               "Determinant (even = %s, odd = %s)"
               (ThoList.to_string (ThoList.to_string factor_to_string) even)
               (ThoList.to_string (ThoList.to_string factor_to_string) odd)
          | _ -> "not a Determinant"

        let make_determinant even odd =
          let make_free_single (tail, tip) =
            single (Free (I tail)) (Free (I tip)) in
          Determinant
            (List.map (List.map make_free_single) even,
             List.map (List.map make_free_single) odd)

        let canonicalize_determinant = function
          | Determinant (even, odd) ->
             Determinant
               (List.sort pcompare (List.map (List.sort pcompare) even),
                List.sort pcompare (List.map (List.sort pcompare) odd))
          | other -> other

        let merge_epsilon_pair eps eps_bar =
          merge
            (to_left_factor (fun i -> position i < 0) (epsilon eps))
            (to_right_factor (fun i -> position i < 0) (epsilon_bar eps_bar))

        let make_even tails tips =
          List.rev_map
            (List.rev_map2 (fun tail tip -> (tail, tip)) tails)
            (Combinatorics.permute_even tips)

        let make_odd tails tips =
          List.rev_map
            (List.rev_map2 (fun tail tip -> (tail, tip)) tails)
            (Combinatorics.permute_odd tips)

        let assert_eps_aux even odd eps eps_bar =
          assert_equal ~printer:determinant_to_string
            (canonicalize_determinant (make_determinant even odd))
            (canonicalize_determinant (merge_epsilon_pair eps eps_bar))

        let assert_eps unit eps eps_bar =
          let tips, tails = List.split unit in
          let even = make_even tails tips
          and odd = make_odd tails tips in
          assert_eps_aux even odd eps eps_bar

        (* A single arrow needs special treatment to get the
           sign to the proper place.*)
        let assert_eps1 odd (tip, tail) eps eps_bar =
          if odd then
            assert_eps_aux [] [[(tail,tip)]] eps eps_bar
          else
            assert_eps_aux [[(tail,tip)]] [] eps eps_bar

        let suite_fuse_epsilons1 =
          "fuse_epsilons1" >:::

            [ "1a*2a" >:: (fun () -> assert_eps1 false (1,2) [1;-9] [2;-9]);
              "a1*a2" >:: (fun () -> assert_eps1 false (1,2) [-9;1] [-9;2]);
              "1a*a2" >:: (fun () -> assert_eps1 true (1,2) [1;-9] [-9;2]);
              "a1*2a" >:: (fun () -> assert_eps1 true (1,2) [-9;1] [2;-9]);

              "13a*24a" >:: (fun () -> assert_eps [(1,2);(3,4)] [1;3;-9] [2;4;-9]);
              "1a3*24a" >:: (fun () -> assert_eps [(1,4);(3,2)] [1;-9;3] [2;4;-9]);
              "a13*2a4" >:: (fun () -> assert_eps [(1,4);(3,2)] [-9;1;3] [2;-9;4]);
              "1a3*2a4" >:: (fun () -> assert_eps [(1,2);(3,4)] [1;-9;3] [2;-9;4]);

              "135a*246a" >:: (fun () -> assert_eps [(1,2);(3,4);(5,6)] [1;3;5;-9] [2;4;6;-9]);
              "315a*246a" >:: (fun () -> assert_eps [(3,2);(1,4);(5,6)] [3;1;5;-9] [2;4;6;-9]);
              "5a13*246a" >:: (fun () -> assert_eps [(1,2);(3,4);(5,6)] [5;-9;1;3] [2;4;6;-9]);
              "a135*2a46" >:: (fun () -> assert_eps [(3,2);(1,4);(5,6)] [-9;1;3;5] [2;-9;4;6]) ] i*)

        let suite =
          "Color.Arrow" >:::
	    [suite_chain;
             suite_cycle;
             suite_take;
             suite_take2;
             suite_replace]

        let suite_long =
          "Color.Arrow long" >:::
	    []

      end

    let pp_free fmt f =
      Format.fprintf fmt "%s" (free_to_string f)

    let pp_factor fmt f =
      Format.fprintf fmt "%s" (factor_to_string f)

  end

module type Propagator =
  sig
    type cf_in = int
    type cf_out = int
    type t = W | I of cf_in | O of cf_out | IO of cf_in * cf_out | G
    val to_string : t -> string
  end

module Propagator : Propagator =
  struct
    type cf_in = int
    type cf_out = int
    type t = W | I of cf_in | O of cf_out | IO of cf_in * cf_out | G
    let to_string = function
      | W -> "W"
      | I cf -> Printf.sprintf "I(%d)" cf
      | O cf' -> Printf.sprintf "O(%d)" cf'
      | IO (cf, cf') -> Printf.sprintf "IO(%d,%d)" cf cf'
      | G -> "G"
  end

module type LP =
  sig
    val rationals : (Algebra.Q.t * int) list -> Algebra.Laurent.t
    val ints : (int * int) list -> Algebra.Laurent.t

    val rational : Algebra.Q.t -> Algebra.Laurent.t
    val int : int -> Algebra.Laurent.t
    val fraction : int -> Algebra.Laurent.t
    val imag : int -> Algebra.Laurent.t
    val nc : int -> Algebra.Laurent.t
    val over_nc : int -> Algebra.Laurent.t
  end

module LP : LP =
  struct
    module L = Algebra.Laurent

    (* Rationals from integers. *)
    let q_int n = Q.make n 1
    let q_fraction n = Q.make 1 n

    (* Complex rationals: *)
    let qc_rational q = QC.make q Q.null
    let qc_int n = qc_rational (q_int n)
    let qc_fraction n = qc_rational (q_fraction n)
    let qc_imag n = QC.make Q.null (q_int n)

    (* Laurent polynomials: *)
    let of_pairs f pairs =
      L.sum (List.map (fun (coeff, power) -> L.atom (f coeff) power) pairs)

    let rationals = of_pairs qc_rational
    let ints = of_pairs qc_int

    let rational q = rationals [(q, 0)]
    let int n = ints [(n, 0)]
    let fraction n = L.const (qc_fraction n)
    let imag n = L.const (qc_imag n)
    let nc n = ints [(n, 1)]
    let over_nc n = ints [(n, -1)]

  end

module type Birdtracks =
  sig
    type t
    val canonicalize : t -> t
    val to_string : t -> string
    val trivial : t -> bool
    val is_null : t -> bool
    val const : Algebra.Laurent.t -> t
    val null : t
    val one : t
    val two : t
    val half : t
    val third : t
    val minus : t
    val int : int -> t
    val fraction : int -> t
    val nc : t
    val over_nc : t
    val imag : t
    val ints : (int * int) list -> t
    val scale : QC.t -> t -> t
    val sum : t list -> t
    val diff : t -> t -> t
    val times : t -> t -> t
    val multiply : t list -> t
    module Infix : sig
      val ( +++ ) : t -> t -> t
      val ( --- ) : t -> t -> t
      val ( *** ) : t -> t -> t
    end
    val f_of_rep : (int -> int -> int -> t) -> int -> int -> int -> t
    val d_of_rep : (int -> int -> int -> t) -> int -> int -> int -> t
    val relocate : (int -> int) -> t -> t
    val fuse : int -> t -> Propagator.t list -> (QC.t * Propagator.t) list
    module Test : Test
    val pp : Format.formatter -> t -> unit
  end

module Birdtracks =
  struct

    module A = Arrow
    open A.Infix
    module P = Propagator
    module L = Algebra.Laurent

    type connection = L.t * A.free list
    type t = connection list

    let trivial = function
      | [] -> true
      | [(coeff, [])] -> coeff = L.unit
      | _ -> false

    (* Rationals from integers. *)
    let q_int n = Q.make n 1
    let q_fraction n = Q.make 1 n

    (* Complex rationals: *)
    let qc_rational q = QC.make q Q.null
    let qc_int n = qc_rational (q_int n)
    let qc_fraction n = qc_rational (q_fraction n)
    let qc_imag n = QC.make Q.null (q_int n)

    (* Laurent polynomials: *)
    let laurent_of_pairs f pairs =
      L.sum (List.map (fun (coeff, power) -> L.atom (f coeff) power) pairs)

    let l_rationals = laurent_of_pairs qc_rational
    let l_ints = laurent_of_pairs qc_int

    let l_rational q = l_rationals [(q, 0)]
    let l_int n = l_ints [(n, 0)]
    let l_fraction n = L.const (qc_fraction n)
    let l_imag n = L.const (qc_imag n)
    let l_nc n = l_ints [(n, 1)]
    let l_over_nc n = l_ints [(n, -1)]

    (* Expressions *)
    let const c = [c, []]
    let ints pairs = const (LP.ints pairs)
    let null = const L.null
    let half = const (LP.fraction 2)
    let third = const (LP.fraction 3)
    let fraction n = const (LP.fraction n)
    let one = const (LP.int 1)
    let two = const (LP.int 2)
    let minus = const (LP.int (-1))
    let int n = const (LP.int n)
    let nc = const (LP.nc 1)
    let over_nc = const (LP.ints [(1, -1)])
    let imag = const (LP.imag 1)

    module AMap = Pmap.Tree

    let find_arrows_opt arrows map =
      try Some (AMap.find pcompare arrows map) with Not_found -> None

    let canonicalize1 (coeff, io_list) =
      (coeff, List.sort pcompare io_list)

    let canonicalize terms =
      let map =
        List.fold_left
          (fun acc term ->
            let coeff, arrows = canonicalize1 term in
            if L.is_null coeff then
              acc
            else
              match find_arrows_opt arrows acc with
              | None -> AMap.add pcompare arrows coeff acc
              | Some coeff' ->
                 let coeff'' = L.add coeff coeff' in
                 if L.is_null coeff'' then
                   AMap.remove pcompare arrows acc
                 else
                   AMap.add pcompare arrows coeff'' acc)
          AMap.empty terms in
      if AMap.is_empty map then
        null
      else
        AMap.fold (fun arrows coeff acc -> (coeff, arrows) :: acc) map []

    let arrows_to_string_aux f arrows =
      ThoList.to_string f arrows

    let to_string1_aux f (coeff, arrows) =
      Printf.sprintf
        "(%s) * %s"
        (L.to_string "N" coeff) (arrows_to_string_aux f arrows)

    let to_string_raw_aux f v =
      ThoList.to_string (to_string1_aux f) v

    let to_string_aux f v =
      to_string_raw_aux f (canonicalize v)

    let factor_arrows_to_string = arrows_to_string_aux A.factor_to_string
    let factor_to_string1 = to_string1_aux A.factor_to_string
    let factor_to_string_raw = to_string_raw_aux A.factor_to_string
    let factor_to_string = to_string_aux A.factor_to_string

    let arrows_to_string = arrows_to_string_aux A.free_to_string
    let to_string1 = to_string1_aux A.free_to_string
    let to_string_raw = to_string_raw_aux A.free_to_string
    let to_string = to_string_aux A.free_to_string

    let pp fmt v =
      Format.fprintf fmt "%s" (to_string v)

    let is_null v =
      List.for_all (fun (c, _) -> L.is_null c) (canonicalize v)

    let is_white = function
      | P.W -> true
      | _ -> false

    let relocate1 f (c, v) =
      (c, List.map (A.map (A.relocate f)) v)

    let relocate f = List.map (relocate1 f)

    (* Only for documentiation: a [term] is a list of arrows with
       a coefficient. *)
    type term = L.t * A.factor list

    (* \begin{dubious}
         New version: there can be $\epsilon$ or $\bar\epsilon$, but
         not both at the same time.
       \end{dubious} *)
    type term_new =
      | Only_Deltas of L.t * A.factor list
      | Epsilons of L.t * A.factor list * A.endpoint list list
      | Epsilon_Bars of L.t * A.factor list * A.endpoint list list

    (* Avoid the recursion, if there is no summation index in [arrow].
       If [arrow] loops back to itself, replace it by a factor of~$N_C$. *)
    let rec add_arrow : A.factor -> term -> term option =
      fun arrow (coeff, arrows) ->
      if A.is_free arrow then
        Some (coeff, arrow :: arrows)
      else if A.is_tadpole arrow then
        Some (L.mul (LP.nc 1) coeff, arrows)
      else
        add_arrow' coeff [] arrow arrows

    (* Add one [arrow] to a list of arrows, updating [coeff]
       if necessary. Accumulate already processed arrows in [seen].
       Returns an empty list if there is a mismatch (a gluon meeting
       a ghost) and a list of pairs consisting of a coefficient and a
       list of arrows otherwise.   There can be more than one pair,
       because matching $\epsilon$ and $\bar\epsilon$ results
       in a sum over permutations. *)

    and add_arrow' : L.t -> A.factor list -> A.factor -> A.factor list -> term option =
      fun coeff seen arrow -> function
      | [] -> (* visited all [arrows]: no opportunities for further matches *)
         Some (coeff, arrow :: seen)
      | arrow' :: arrows' ->
         begin match A.merge arrow arrow' with
         | A.Mismatch ->
            None
         | A.Ghost_Match -> (* replace matching ghosts by $-1/N_C$ *)
            Some (L.mul (LP.over_nc (-1)) coeff, List.rev_append seen arrows')
         | A.Loop_Match -> (* replace a loop by $N_C$ *)
            Some (L.mul (LP.nc 1) coeff, List.rev_append seen arrows')
         | A.Match arrow'' -> (* two arrows have been merged into one *)
            if A.is_free arrow'' then (* no opportunities for further matches *)
              Some (coeff, arrow'' :: List.rev_append seen arrows')
            else (* the new [arrow''] ist not yet saturated, try again: *)
              add_arrow' coeff seen arrow'' arrows'
         | A.No_Match -> (* recurse to the remaining arrows *)
            add_arrow' coeff (arrow' :: seen)  arrow arrows'
         end

    and add_determinant : A.factor list -> A.factor list list -> A.factor list list -> term -> term list =
      fun seen even odd (coeff, arrows as term) ->
      distribute seen even term (distribute seen odd (L.neg coeff, arrows) [])

    and distribute : A.factor list -> A.factor list list-> term -> term list -> term list =
      fun seen permutations term terms ->
      List.fold_left
        (fun acc permutation -> splice_arrows seen permutation term :: acc)
        terms permutations

    and splice_arrows : A.factor list -> A.factor list -> term -> term =
      fun seen arrows term  ->
      let coeff', arrows' = add_arrow_list arrows term in
      (coeff', List.rev_append seen arrows')

    (* \begin{dubious}
         Here we would like to use the type system to prove
         that the two failing cases can't happen.
         In real life they can't happen, because [arrow] is
         never [A.Epsilon].  Can we use the version of GADTs
         that are available in O'Caml 4.05?
       \end{dubious} *)
    and add_arrow_list : A.factor list -> term -> term =
      fun arrows term ->
      match arrows with
      | [] -> term
      | arrow :: rest ->
         begin match add_arrow arrow term with
         | Some term -> add_arrow_list rest term
         | None -> failwith "add_arrow_list: unexpected None"
         end

(*i   and add_arrow_list arrows (coeff, acc as term) =
      let result = add_arrow_list_raw arrows term in
      Printf.eprintf
        "add_arrow_list (%s) * %s %s ==> %s\n"
        (factor_arrows_to_string arrows)
        (L.to_string "N" coeff)
        (factor_arrows_to_string acc)
        (factor_to_string1 result);
      result
i*)

(*i
    let logging_add_arrow arrow (coeff, arrows) =
      let result = add_arrow arrow (coeff, arrows) in
      Printf.eprintf
        "add_arrow %s to %s ==> %s\n"
        (A.factor_to_string arrow)
        (factor_to_string1 (coeff, arrows))
        (ThoList.to_string factor_to_string1 result);
      result

    let add_arrow = logging_add_arrow
i*)

(* The return type is [term list], because adding an~$\epsilon$
   (or~$\bar\epsilon$) will turn a term to a sum of terms iff
   the term contains a $\bar\epsilon$ (or~$\epsilon$), since
   $\epsilon$-$\bar\epsilon$ pairs will
   be expanded by [add_determinant]. *)

    let add_arrow_or_epsilon : A.factor -> term -> term list =
      fun arrow (coeff, arrows as term) ->
      if A.is_epsilon arrow then
        match A.match_epsilon arrow arrows with
        | None ->
           begin match add_arrow' coeff [] arrow arrows with
           | None -> []
           | Some term -> [term]
           end
        | Some (even, odd, arrows') ->
           add_determinant [] even odd (coeff, arrows')
      else
        match add_arrow arrow term with
        | None -> []
        | Some term -> [term]

    (* [add_arrows_or_epsilons arrows term] add the [arrows] to [term] by
       calling [add_arrow_or_epsilon] for each one.  Return an empty list
       if there are leftover summation indices in the end. *)

    (* NB: we can reject the contributions with unsaturated summation indices
       from Ghost contributions to~$T_a$ only \emph{after} adding all
       arrows that might saturate an open index. *)

    let rec add_arrows_or_epsilons : A.factor list -> term -> term list =
      fun arrows (_, acc_arrows as term) ->
      match arrows with
      | [] ->
         if List.for_all A.is_free acc_arrows then
           [term]
         else
           []
      | arrow :: rest ->
         ThoList.flatmap (add_arrows_or_epsilons rest) (add_arrow_or_epsilon arrow term)


    let logging_add_arrows_or_epsilons arrows term =
      let result = add_arrows_or_epsilons arrows term in
      Printf.eprintf
        "add_arrows_or_epsilons %s to %s ==> %s\n"
        (factor_to_string1 term)
        (factor_arrows_to_string arrows)
        (ThoList.to_string factor_to_string1 result);
      result

(*i
    let add_arrows_or_epsilons = logging_add_arrows_or_epsilons
i*)
    (* Note that a negative index might be summed only
       later in a sequence of binary products and must
       therefore be treated as free in this product.  Therefore,
       we have to classify the indices as summation indices
       \emph{not only} based on their sign, but in addition based on
       whether they appear in both factors. Only then can we reject
       surviving ghosts. *)

    module ESet =
      Set.Make
        (struct
          type t = A.endpoint
          let compare = pcompare
        end)

    let negatives arrows =
      List.fold_left
        (fun acc arrow ->
          List.fold_left
            (fun acc' i -> ESet.add i acc')
            acc (A.negatives arrow))
        ESet.empty arrows

    let times1 (coeff1, arrows1) (coeff2, arrows2) =
      let summations = ESet.inter (negatives arrows1) (negatives arrows2) in
      let is_sum i = ESet.mem i summations in
      let arrows1' = List.map (A.to_left_factor is_sum) arrows1
      and arrows2' = List.map (A.to_right_factor is_sum) arrows2 in
      List.map
        (fun (coeff1, arrows) ->
          (L.mul coeff1 coeff2, List.map A.of_factor arrows))
        (add_arrows_or_epsilons arrows2' (coeff1, arrows1'))

    let logging_times1 factor1 factor2 =
      let result = times1 factor1 factor2 in
      Printf.eprintf
        "%s times1 %s ==> %s\n"
        (to_string1 factor1)
        (to_string1 factor2)
        (ThoList.to_string to_string1 result);
      result

    let sum terms =
      canonicalize (List.concat terms)

    let times term term' =
      canonicalize
        (Product.fold2
           (fun x y -> List.rev_append (times1 x y))
           term term' [])

    (* \begin{dubious}
         Is that more efficient than the following implementation?
       \end{dubious} *)

    let rec multiply1' acc = function
      | [] -> [acc]
      | factor :: factors ->
         List.fold_right multiply1' (times1 acc factor) factors

    let multiply1 = function
      | [] -> [(L.unit, [])]
      | [factor] -> [factor]
      | factor :: factors -> multiply1' factor factors

    let multiply terms =
      canonicalize
        (Product.fold (fun x -> List.rev_append (multiply1 x)) terms [])

    (* \begin{dubious}
         Isn't that the more straightforward implementation?
       \end{dubious} *)

    let multiply = function
      | [] -> []
      | term :: terms ->
         canonicalize (List.fold_left times term terms)

    let scale1 q (coeff, arrows) =
      (L.scale q coeff, arrows)
    let scale q = List.map (scale1 q)

    let diff term1 term2 =
      canonicalize (List.rev_append term1 (scale (qc_int (-1)) term2))

    module Infix =
      struct
        let ( +++ ) term term' = sum [term; term']
        let ( --- ) = diff
        let ( *** ) = times
      end

    open Infix

    (* Compute $ \tr(r(T_a) r(T_b) r(T_c)) $.  NB: this uses the
       summation indices $-1$, $-2$ and $-3$.  Therefore
       it \emph{must not} appear unevaluated more than once in a product! *)
    let trace3 r a b c =
      r a (-1) (-2) *** r b (-2) (-3) *** r c (-3) (-1)

    let f_of_rep r a b c =
      minus *** imag *** (trace3 r a b c --- trace3 r a c b)

    (* $ d_{abc} = \tr(r(T_a) [r(T_b), r(T_c)]_+) $ *)
    let d_of_rep r a b c =
      trace3 r a b c +++ trace3 r a c b

(* \thocwmodulesubsection{Feynman Rules} *)
    module IMap =
      Map.Make (struct type t = int let compare = pcompare end)

    let line_map lines =
      let _, map =
        List.fold_left
          (fun (i, acc) line ->
            (succ i,
             match line with
             | P.W -> acc
             | _ -> IMap.add i line acc))
          (1, IMap.empty)
          lines in
      map

(*i Redundant since ocaml 4.05
     let find_opt i map =
      try Some (IMap.find i map) with Not_found -> None
i*)

    let lines_to_string lines =
      match IMap.bindings lines with
      | [] -> "W"
      | lines ->
         String.concat
           " "
           (List.map
              (fun (i, c) -> Printf.sprintf "%s@%d" (P.to_string c) i)
              lines)

    let clear = IMap.remove

    let add_in i cf lines =
      match IMap.find_opt i lines with
      | Some (P.O cf') -> IMap.add i (P.IO (cf, cf')) lines
      | _ -> IMap.add i (P.I cf) lines

    let add_out i cf' lines =
      match IMap.find_opt i lines with
      | Some (P.I cf) -> IMap.add i (P.IO (cf, cf')) lines
      | _ -> IMap.add i (P.O cf') lines

    let add_ghost i lines =
      IMap.add i P.G lines

    let connect1 n arrow lines =
      match arrow with
      | A.Ghost g ->
         let g = A.position g in
         if g = n then
           Some (add_ghost n lines)
         else
           begin match IMap.find_opt g lines with
           | Some P.G -> Some (clear g lines)
           | _ -> None
           end
      | A.Arrow (i, o) ->
         let i = A.position i
         and o = A.position o in
         if o = n then
           begin match IMap.find_opt i lines with
           | Some (P.I cfi) -> Some (add_in o cfi (clear i lines))
           | Some (P.IO (cfi, cfi')) -> Some (add_in o cfi (add_out i cfi' lines))
           | _ -> None
           end
         else if i = n then
           begin match IMap.find_opt o lines with
           | Some (P.O cfo') -> Some (add_out i cfo' (clear o lines))
           | Some (P.IO (cfo, cfo')) -> Some (add_out i cfo' (add_in o cfo lines))
           | _ -> None
           end
         else
           begin match IMap.find_opt i lines, IMap.find_opt o lines with
           | Some (P.I cfi), Some (P.O cfo') when cfi = cfo' ->
              Some (clear o (clear i lines))
           | Some (P.I cfi), Some (P.IO (cfo, cfo')) when cfi = cfo'->
              Some (add_in o cfo (clear i lines))
           | Some (P.IO (cfi, cfi')), Some (P.O cfo') when cfi = cfo' ->
              Some (add_out i cfi' (clear o lines))
           | Some (P.IO (cfi, cfi')), Some (P.IO (cfo, cfo')) when cfi = cfo' ->
              Some (add_in o cfo (add_out i cfi' lines))
           | _ -> None
           end
      | A.Epsilon _  ->
        failwith "Birdtracks.connect not yet defined for Epsilon"
      | A.Epsilon_bar _ ->
         failwith "Birdtracks.connect not yet defined for Epsilon_bar"
        
    let connect connections lines =
      let n = succ (List.length lines)
      and lines = line_map lines in
      let rec connect' acc = function
        | arrow :: arrows ->
           begin match connect1 n arrow acc with
           | None -> None
           | Some acc -> connect' acc arrows
           end
        | [] -> Some acc in
      match connect' lines connections with
      | None -> None
      | Some acc ->
         begin match IMap.bindings acc with
         | [] -> Some P.W
         | [(i, cf)] when i = n -> Some cf
         | _ -> None
         end

    let fuse1 nc lines (c, vertex) =
      match connect vertex lines with
      | None -> []
      | Some cf -> [(L.eval (qc_int nc) c, cf)]
             
    let fuse nc vertex lines =
      match vertex with
      | [] ->
         if List.for_all is_white lines then
           [(QC.unit, P.W)]
         else
           []
      | vertex ->
         ThoList.flatmap (fuse1 nc lines) vertex

    module Test : Test =
      struct
        open OUnit

        let vertices_equal v1 v2 =
          (canonicalize v1) = (canonicalize v2)

        let eq v1 v2 =
          assert_equal ~printer:(ThoList.to_string to_string1) ~cmp:vertices_equal v1 v2

        let suite_times1 =
          "times1" >:::

            [ "merge two" >::
	        (fun () ->
	          eq
                    [(L.unit, 1 ==> 2)]
                    (times1 (L.unit,  1 ==> -1) (L.unit, -1 ==>  2)));

              "merge two exchanged" >::
	        (fun () ->
	          eq
                    [(L.unit, 1 ==> 2)]
                    (times1 (L.unit, -1 ==>  2) (L.unit,  1 ==> -1)));

              "ghost1" >::
	        (fun () ->
	          eq
                    [(l_over_nc (-1), 1 ==> 2)]
                    (times1
                       (L.unit, [-1 =>  2; ?? (-3)])
                       (L.unit, [ 1 => -1; ?? (-3)])));

              "ghost2" >::
	        (fun () ->
	          eq
                    []
                    (times1
                       (L.unit, [ 1 => -1; ?? (-3)])
                       (L.unit, [-1 =>  2; -3 => -4; -4 => -3])));

              "ghost2 exchanged" >::
	        (fun () ->
	          eq
                    []
                    (times1
                       (L.unit, [-1 =>  2; -3 => -4; -4 => -3])
                       (L.unit, [ 1 => -1; ?? (-3)]))) ]

        let suite_canonicalize =
          "canonicalize" >:::

            [ ]

        let line_option_to_string = function
          | None -> "no match"
          | Some line -> P.to_string line

        let test_connect_msg vertex formatter (expected, result) =
          Format.fprintf
            formatter
            "[%s]: expected %s, got %s"
            (arrows_to_string vertex)
            (line_option_to_string expected)
            (line_option_to_string result)

        let test_connect expected lines vertex =
	  assert_equal
            ~printer:line_option_to_string
            expected (connect vertex lines)

        let test_connect_permutations expected lines vertex =
          List.iter
            (fun v ->
	      assert_equal
                ~pp_diff:(test_connect_msg v)
                expected (connect v lines))
            (Combinatorics.permute vertex)

        let suite_connect =
          "connect" >:::

            [ "delta" >::
	        (fun () ->
                  test_connect_permutations
                    (Some (P.I 1))
                    [ P.I 1; P.W ]
                    ( 1 ==> 3 ));

              "f: 1->3->2->1" >::
                (fun () ->
                  test_connect_permutations
                    (Some (P.IO (1, 3)))
                    [P.IO (1, 2); P.IO (2, 3)]
                    (A.cycle [1; 3; 2]));

              "f: 1->2->3->1" >::
                (fun () ->
                  test_connect_permutations
                    (Some (P.IO (1, 2)))
                    [P.IO (3, 2); P.IO (1, 3)]
                    (A.cycle [1; 2; 3])) ]

        let suite =
          "Color.Birdtracks" >:::
	    [suite_times1;
             suite_canonicalize;
             suite_connect]

        let suite_long =
          "Color.Birdtracks long" >:::
	    []

      end

    let vertices_equal v1 v2 =
      is_null (v1 --- v2)

    let assert_zero_vertex v =
      OUnit.assert_equal ~printer:to_string ~cmp:vertices_equal null v

    (* As an extra protection agains vacuous tests, we make
       sure that the LHS does not vanish.  *)
    let eq v1 v2 =
      OUnit.assert_bool "LHS = 0" (not (is_null v1));
      OUnit.assert_equal ~printer:to_string ~cmp:vertices_equal v1 v2

  end
    
(* \thocwmodulesection{$\mathrm{SU}(N_C)$}
   We're computing with a general $N_C$, but [epsilon] and [epsilonbar]
   make only sense for $N_C=3$.  Also some of the terminology alludes
   to $N_C=3$: triplet, sextet, octet. *)

(* Using the normalization~$\tr(T_{a}T_{b}) = \delta_{ab}$, we can
   check the selfconsistency of the completeness relation
   \begin{equation}
       T_{a}^{i_1j_1} T_{a}^{i_2j_2} =
         \left(                 \delta^{i_1j_2} \delta^{i_2j_1}
                - \frac{1}{N_C} \delta^{i_1j_1} \delta^{j_1j_2}\right)
   \end{equation}
   as
   \begin{multline}
     T_{a}^{i_1j_1} T_{a}^{i_2j_2}
       = \tr\left(T_{a_1}T_{a_2}\right) T_{a_1}^{i_1j_1} T_{a_2}^{i_2j_2}
       = T_{a_1}^{l_1l_2} T_{a_2}^{l_2l_1}
         T_{a_1}^{i_1j_1} T_{a_2}^{i_2j_2} \\
       = \left(                 \delta^{l_1j_1} \delta^{i_1l_2}
                - \frac{1}{N_C} \delta^{l_1l_2} \delta^{i_1j_1}\right)
         \left(                 \delta^{l_2j_2} \delta^{i_2l_1}
                - \frac{1}{N_C} \delta^{l_2l_1} \delta^{i_2j_2}\right)
       = \left(                 \delta^{i_1j_2} \delta^{i_2j_1}
                - \frac{1}{N_C} \delta^{i_1i_2} \delta^{j_2j_1}\right)
   \end{multline}
   With
   \begin{equation}
   \label{eq:f=tr(TTT)'}
     \ii f_{a_1a_2a_3}
       = \tr\left(T_{a_1}\left\lbrack T_{a_2},T_{a_3}\right\rbrack\right)
       = \tr\left(T_{a_1}T_{a_2}T_{a_3}\right)
       - \tr\left(T_{a_1}T_{a_3}T_{a_2}\right)
   \end{equation}
   and
   \begin{multline}
     \tr\left(T_{a_1}T_{a_2}T_{a_3}\right)
         T_{a_1}^{i_1j_1} T_{a_2}^{i_2j_2} T_{a_3}^{i_3j_3}
       = T_{a_1}^{l_1l_2} T_{a_2}^{l_2l_3} T_{a_3}^{l_3l_1}
         T_{a_1}^{i_1j_1} T_{a_2}^{i_2j_2} T_{a_3}^{i_3j_3} = \\
         \left(                 \delta^{l_1j_1} \delta^{i_1l_2}
                - \frac{1}{N_C} \delta^{l_1l_2} \delta^{i_1j_1}\right)
         \left(                 \delta^{l_2j_2} \delta^{i_2l_3}
                - \frac{1}{N_C} \delta^{l_2l_3} \delta^{i_2j_2}\right)
         \left(                 \delta^{l_3j_3} \delta^{i_3l_1}
                - \frac{1}{N_C} \delta^{l_3l_1} \delta^{i_3j_3}\right)
   \end{multline}
   we find the decomposition
   \begin{equation}
   \label{eq:fTTT'}
       \ii f_{a_1a_2a_3} T_{a_1}^{i_1j_1}T_{a_2}^{i_2j_2}T_{a_3}^{i_3j_3}
     = \delta^{i_1j_2}\delta^{i_2j_3}\delta^{i_3j_1}
     - \delta^{i_1j_3}\delta^{i_3j_2}\delta^{i_2j_1}\,.
   \end{equation} *)

(*  Indeed,
\begin{verbatim}
symbol nc;
Dimension nc;
vector i1, i2, i3, j1, j2, j3;
index l1, l2, l3;

local [TT] =
        ( j1(l1) * i1(l2) - d_(l1,l2) * i1.j1 / nc )
      * ( j2(l2) * i2(l1) - d_(l2,l1) * i2.j2 / nc );

#procedure TTT(sign)
local [TTT`sign'] =
        ( j1(l1) * i1(l2) - d_(l1,l2) * i1.j1 / nc )
      * ( j2(l2) * i2(l3) - d_(l2,l3) * i2.j2 / nc )
      * ( j3(l3) * i3(l1) - d_(l3,l1) * i3.j3 / nc )
 `sign' ( j1(l1) * i1(l2) - d_(l1,l2) * i1.j1 / nc )
      * ( j3(l2) * i3(l3) - d_(l2,l3) * i3.j3 / nc )
      * ( j2(l3) * i2(l1) - d_(l3,l1) * i2.j2 / nc );
#endprocedure

#call TTT(-)
#call TTT(+)

bracket nc;
print;
.sort
.end
\end{verbatim}
gives
\begin{verbatim}
   [TT] =
       + nc^-1 * (  - i1.j1*i2.j2 )
       + i1.j2*i2.j1;

   [TTT-] =
       + i1.j2*i2.j3*i3.j1 - i1.j3*i2.j1*i3.j2;

   [TTT+] =
       + nc^-2 * (    4*i1.j1*i2.j2*i3.j3 )
       + nc^-1 * (  - 2*i1.j1*i2.j3*i3.j2
                    - 2*i1.j2*i2.j1*i3.j3
                    - 2*i1.j3*i2.j2*i3.j1 )
       + i1.j2*i2.j3*i3.j1 + i1.j3*i2.j1*i3.j2;
\end{verbatim}
*)

module type SU3 =
  sig
    include Birdtracks
    val delta3 : int -> int -> t
    val delta8 : int -> int -> t
    val delta8_loop : int -> int -> t
    val gluon : int -> int -> t
    val delta6 : int -> int -> t
    val delta10 : int -> int -> t
    val t : int -> int -> int -> t
    val f : int -> int -> int -> t
    val d : int -> int -> int -> t
    val epsilon : int list -> t
    val epsilon_bar : int list -> t
    val t8 : int -> int -> int -> t
    val t6 : int -> int -> int -> t
    val t10 : int -> int -> int -> t
    val k6 : int -> int -> int -> t
    val k6bar : int -> int -> int -> t
    val delta_of_tableau : int Young.tableau -> int -> int -> t
    val t_of_tableau : int Young.tableau -> int -> int -> int -> t
  end

module SU3 : SU3 =
  struct

    module A = Arrow
    open Arrow.Infix

    module B = Birdtracks
    type t = B.t
    let canonicalize = B.canonicalize
    let to_string = B.to_string
    let pp = B.pp
    let trivial = B.trivial
    let is_null = B.is_null
    let null = B.null
    let const = B.const
    let one = B.one
    let two = B.two
    let int = B.int
    let half = B.half
    let third = B.third
    let fraction = B.fraction
    let nc = B.nc
    let over_nc = B.over_nc
    let minus = B.minus
    let imag = B.imag
    let ints = B.ints
    let sum = B.sum
    let diff = B.diff
    let scale = B.scale
    let times = B.times
    let multiply = B.multiply
    let relocate = B.relocate
    let fuse = B.fuse
    let f_of_rep = B.f_of_rep
    let d_of_rep = B.d_of_rep
    module Infix = B.Infix

(* \thocwmodulesubsection{Fundamental and Adjoint Representation} *)

    let delta3 i j =
      [(LP.int 1, j ==> i)]

    let delta8 a b =
      [(LP.int 1, a <=> b)]

    (* If the~$\delta_{ab}$ originates from
       a~$\tr(T_aT_b)$, like an effective~$gg\to H$
       coupling, it makes a difference in the color
       flow basis and we must write the full expression~(6.2)
       from~\cite{Kilian:2012pz} including the ghosts instead.
       Note that the sign for the terms with one ghost
       has not been spelled out in that reference. *)

    let delta8_loop a b =
      [(LP.int 1, a <=> b);
       (LP.int (-1), [a => a; ?? b]);
       (LP.int (-1), [?? a; b => b]);
       (LP.nc 1, [?? a; ?? b])]

    (* The following can be used for computing polarization sums
       (eventually, this could make the [Flow] module redundant).
       Note that we have $-N_C$ instead of $-1/N_C$ in the ghost
       contribution here, because
       two factors of $-1/N_C$ will be produced by [add_arrow]
       below, when contracting two ghost indices.
       Indeed, with this definition we can maintain
       [multiply [delta8 1 (-1); gluon (-1) (-2); delta8 (-2) 2]
        = delta8 1 2]. *)

    let ghost a b =
      [ (LP.nc (-1), [?? a; ?? b])]

    let gluon a b =
      delta8 a b @ ghost a b

    (* Note that the arrow is directed from the second to the first
       index, opposite to our color flow paper~\cite{Kilian:2012pz}.
       Fortunately, this is just a matter of conventions.
\begin{subequations}
\begin{align}
\parbox{28\unitlength}{%
  \fmfframe(4,4)(4,4){%
  \begin{fmfgraph*}(20,20)
    \fmfleft{f1,f2}
    \fmfright{g}
    \fmfv{label=$i$}{f2}
    \fmfv{label=$j$}{f1}
    \fmfv{label=$a$}{g}
    \fmf{fermion}{f1,v}
    \fmf{fermion}{v,f2}
    \fmf{gluon}{v,g}
  \end{fmfgraph*}}} &\Longrightarrow
\parbox{28\unitlength}{%
  \fmfframe(4,4)(4,4){%
  \begin{fmfgraph*}(20,20)
    \fmfleft{f1,f2}
    \fmfright{g}
    \fmfv{label=$i$}{f2}
    \fmfv{label=$j$}{f1}
    \fmfv{label=$a$}{g}
    \fmf{phantom}{f1,v}
    \fmf{phantom}{v,f2}
    \fmf{phantom}{v,g}
    \fmffreeze
    \fmfi{phantom_arrow}{vpath (__v, __g) sideways -thick}
    \fmfi{phantom_arrow}{(reverse vpath (__v, __g)) sideways -thick}
    \fmfi{phantom_arrow}{vpath (__f1, __v)}
    \fmfi{phantom_arrow}{vpath (__v, __f2)}
    \fmfi{plain}{%
      (vpath (__f1, __v) join (vpath (__v, __g)) sideways -thick)}
    \fmfi{plain}{%
      ((reverse vpath (__g, __v) sideways -thick) join vpath (__v, __f2))}
  \end{fmfgraph*}}}
\parbox{28\unitlength}{%
  \fmfframe(4,4)(4,4){%
  \begin{fmfgraph*}(20,20)
    \fmfleft{f1,f2}
    \fmfright{g}
    \fmfv{label=$i$}{f1}
    \fmfv{label=$j$}{f2}
    \fmfv{label=$a$}{g}
    \fmf{fermion}{f1,v}
    \fmf{fermion}{v,f2}
    \fmf{dots}{v,g}
  \end{fmfgraph*}}}\\
  T_a^{ij} \qquad\quad
    &\Longrightarrow \qquad\quad \delta^{ia}\delta^{aj}
       \qquad\qquad\qquad - \delta^{ij}
\end{align}
\end{subequations} *)

    let t a i j =
      [ (LP.int 1, [j => a; a => i]);
        (LP.int (-1), [j => i; ?? a]) ]

(* Note that while we expect $\tr(T_a)=T_a^{ii}=0$,
   the evaluation of the expression [t 1 (-1) (-1)] will stop
   at [ [ -1 => 1; 1 => -1 ] --- [ -1 => -1; ?? 1 ] ], because the
   summation index appears in a single term.
   However, a naive further evaluation would get stuck at
   [ [ 1 => 1 ] --- nc *** [ ?? 1 ] ].
   Fortunately, traces of single generators are never needed in our
   applications.  We just have to resist the temptation to use them
   in unit tests. *)

(*
\begin{equation}
\parbox{29\unitlength}{%
  \fmfframe(2,2)(2,2){%
  \begin{fmfgraph*}(25,25)
    \fmfleft{g1,g2}
    \fmfright{g3}
    \fmfv{label=$a$}{g1}
    \fmfv{label=$b$}{g2}
    \fmfv{label=$c$}{g3}
    \fmf{gluon}{g1,v}
    \fmf{gluon}{g2,v}
    \fmf{gluon}{g3,v}
  \end{fmfgraph*}}}
\qquad\Longrightarrow
\parbox{29\unitlength}{%
  \fmfframe(2,2)(2,2){%
  \begin{fmfgraph*}(25,25)
    \fmfleft{g1,g2}
    \fmfright{g3}
    \fmfv{label=$a$}{g1}
    \fmfv{label=$b$}{g2}
    \fmfv{label=$c$}{g3}
    \fmf{phantom}{g1,v}
    \fmf{phantom}{g2,v}
    \fmf{phantom}{g3,v}
    \fmffreeze
    \fmfi{plain}{(vpath(__g1,__v) join (reverse vpath(__g2,__v))) 
                 sideways thick}
    \fmfi{plain}{(vpath(__g2,__v) join (reverse vpath(__g3,__v)))
                 sideways thick}
    \fmfi{plain}{(vpath(__g3,__v) join (reverse vpath(__g1,__v)))
                 sideways thick}
    \fmfi{phantom_arrow}{vpath (__g1, __v) sideways thick}
    \fmfi{phantom_arrow}{vpath (__g2, __v) sideways thick}
    \fmfi{phantom_arrow}{vpath (__g3, __v) sideways thick}
    \fmfi{phantom_arrow}{(reverse vpath (__g1, __v)) sideways thick}
    \fmfi{phantom_arrow}{(reverse vpath (__g2, __v)) sideways thick}
    \fmfi{phantom_arrow}{(reverse vpath (__g3, __v)) sideways thick}
  \end{fmfgraph*}}}
\qquad
\parbox{29\unitlength}{%
  \fmfframe(2,2)(2,2){%
  \begin{fmfgraph*}(25,25)
    \fmfleft{g1,g2}
    \fmfright{g3}
    \fmfv{label=$a$}{g1}
    \fmfv{label=$b$}{g2}
    \fmfv{label=$c$}{g3}
    \fmf{phantom}{g1,v}
    \fmf{phantom}{g2,v}
    \fmf{phantom}{g3,v}
    \fmffreeze
    \fmfi{plain}{(vpath(__g1,__v) join (reverse vpath(__g3,__v))) 
                 sideways thick}
    \fmfi{plain}{(vpath(__g2,__v) join (reverse vpath(__g1,__v)))
                 sideways thick}
    \fmfi{plain}{(vpath(__g3,__v) join (reverse vpath(__g2,__v)))
                 sideways thick}
    \fmfi{phantom_arrow}{vpath (__g1, __v) sideways thick}
    \fmfi{phantom_arrow}{vpath (__g2, __v) sideways thick}
    \fmfi{phantom_arrow}{vpath (__g3, __v) sideways thick}
    \fmfi{phantom_arrow}{(reverse vpath (__g1, __v)) sideways thick}
    \fmfi{phantom_arrow}{(reverse vpath (__g2, __v)) sideways thick}
    \fmfi{phantom_arrow}{(reverse vpath (__g3, __v)) sideways thick}
  \end{fmfgraph*}}}
\end{equation} *)

    let f a b c =
      [ (LP.imag ( 1), A.cycle [a; b; c]);
        (LP.imag (-1), A.cycle [a; c; b]) ]

(* The generator in the adjoint representation $T_a^{bc}=-\ii f_{abc}$: *)
    let t8 a b c =
      Birdtracks.Infix.( minus *** imag *** f a b c )

(* This $d_{abc}$ is now compatible with~(6.11) in our color
   flow paper~\cite{Kilian:2012pz}.  The signs had been wrong
   in earlier versions of the code to match the missing
   sign in the ghost contribution to the generator~$T_a^{ij}$
   above. *)

    let d a b c =
      [ (LP.int 1, A.cycle [a; b; c]);
        (LP.int 1, A.cycle [a; c; b]);
        (LP.int (-2), (a <=> b) @ [?? c]);
        (LP.int (-2), (b <=> c) @ [?? a]);
        (LP.int (-2), (c <=> a) @ [?? b]);
        (LP.int 2, [a => a; ?? b; ?? c]);
        (LP.int 2, [?? a; b => b; ?? c]);
        (LP.int 2, [?? a; ?? b; c => c]);
        (LP.nc (-2), [?? a; ?? b; ?? c]) ]

(* \thocwmodulesubsection{Decomposed Tensor Product Representations} *)

    let pass_through m n incoming outgoing =
      List.rev_map2 (fun i o -> (m, i) >=>> (n, o)) incoming outgoing

    let delta_of_permutations n permutations k l =
      let incoming = ThoList.range 0 (pred n)
      and normalization = List.length permutations in
      List.rev_map
        (fun (eps, outgoing) ->
          (LP.fraction (eps * normalization),
           pass_through l k incoming outgoing))
        permutations

    let totally_symmetric n =
      List.map
        (fun p -> (1, p))
        (Combinatorics.permute (ThoList.range 0 (pred n)))

    let totally_antisymmetric n =
        (Combinatorics.permute_signed (ThoList.range 0 (pred n)))

    let delta_S n k l =
      delta_of_permutations n (totally_symmetric n) k l

    let delta_A n k l =
      delta_of_permutations n (totally_antisymmetric n) k l

    let delta6 = delta_S 2
    let delta10 = delta_S 3
    let delta15 = delta_S 4

    let delta3bar = delta_A 2

    (* Mixed symmetries, as in section 9.4 of the birdtracks book. *)

    module IM = Partial.Make (struct type t = int let compare = pcompare end)
    module P = Permutation.Default

(* Map the elements of [original] to [permuted] in [all], with [all]
   a list of $n$ integers from $0$ to $n-1$ in order, and use the resulting
   list to define a permutation.
   E.\,g.~[permute_partial [1;3] [3;1] [0;1;2;3;4]] will define a
   permutation that transposes the second and fourth element in
   a 5 element list. *)
    let permute_partial original permuted all =
      P.of_list (List.map (IM.auto (IM.of_lists original permuted)) all)
                         
    let apply1 (sign, indices) (eps, p) =
      (eps * sign, P.list p indices)

    let apply signed_permutations signed_indices =
      List.rev_map (apply1 signed_indices) signed_permutations

    let apply_list signed_permutations signed_indices =
      ThoList.flatmap (apply signed_permutations) signed_indices

    let symmetrizer_of_permutations n original signed_permutations =
      let incoming = ThoList.range 0 (pred n) in
      List.rev_map
        (fun (eps, permuted) ->
          (eps, permute_partial original permuted incoming))
        signed_permutations

    let symmetrizer n indices =
      symmetrizer_of_permutations
        n indices
        (List.rev_map (fun p -> (1, p)) (Combinatorics.permute indices))

    let anti_symmetrizer n indices =
      symmetrizer_of_permutations
        n indices
        (Combinatorics.permute_signed indices)

    let symmetrize n elements indices =
      apply_list (symmetrizer n elements) indices

    let anti_symmetrize n elements indices =
      apply_list (anti_symmetrizer n elements) indices
      
    let id n =
      [(1, ThoList.range 0 (pred n))]

    (* \begin{dubious}
         We can avoid the recursion here, if we use
         [Combinatorics.permute_tensor_signed] in
         [symmetrizer] above.
       \end{dubious} *)
    let rec apply_tableau f n tableau indices =
      match tableau with
      | [] | [_] :: _ -> indices
      | cells :: rest ->
         apply_tableau f n rest (f n cells indices)

(* \begin{dubious}
     Here we should at a sanity test for [tableau]: all integers should
     be consecutive starting from 0 with no duplicates.  In additions
     the rows must not grow in length.
   \end{dubious} *)

    let delta_of_tableau tableau i j =
      let n = Young.num_cells_tableau tableau
      and num, den = Young.normalization (Young.diagram_of_tableau tableau)
      and rows = tableau
      and cols = Young.conjugate_tableau tableau in
      let permutations =
        apply_tableau symmetrize n rows (apply_tableau anti_symmetrize n cols (id n)) in
      Birdtracks.Infix.( int num *** fraction den *** delta_of_permutations n permutations i j )

    let incomplete tensor =
      failwith ("Color.Vertex: " ^ tensor ^ " not supported yet!")

    let experimental tensor =
      Printf.eprintf
        "Color.Vertex: %s support still experimental and untested!\n"
        tensor

    let distinct integers =
      let rec distinct' seen = function
        | [] -> true
        | i :: rest ->
           if Sets.Int.mem i seen then
             false
           else
             distinct' (Sets.Int.add i seen) rest in
      distinct' Sets.Int.empty integers
      
    (* All lines start here: they point towards the vertex. *)
    let epsilon tips =
      if distinct tips then
        [(LP.int 1, [Arrow.epsilon tips])]
      else
        null

    (* All lines end here: they point away from the vertex. *)
    let epsilon_bar tails =
      if distinct tails then
        [(LP.int 1, [Arrow.epsilon_bar tails])]
      else
        null


(* In order to get the correct $N_C$ dependence of
   quadratic Casimir operators, the arrows in the vertex must
   have the same permutation symmetry as the propagator.  This
   is demonstrated by the unit tests involving Casimir operators
   on page \pageref{pg:casimir-tests} below.  These tests also
   provide a check of our normalization.

   The implementation takes a propagator and uses [Arrow.tee] to
   replace one arrow by the pair of arrows corresponding to the
   insertion of a gluon.  This is repeated for each arrow.
   The normalization remains unchanged from the propagator.
   A minus sign is added for antiparallel arrows, since the
   conjugate representation is~$-T^*_a$.

   To this, we add the diagrams with a gluon connected to one arrow.
   Since these are identical, only one diagram multiplied by the
   difference of the number of parallel and antiparallel arrows
   is added. *)

    let insert_gluon a k l (norm, arrows) =
      let rec insert_gluon' acc left = function
        | [] -> acc
        | arrow :: right ->
           insert_gluon'
             ((Algebra.Laurent.mul (LP.int (A.dir k l arrow)) norm,
               List.rev_append left ((A.tee a arrow) @ right)) :: acc)
             (arrow :: left)
             right in
      insert_gluon' [] [] arrows

    let t_of_delta delta a k l =
      match delta k l with
      | [] -> []
      | (_, arrows) :: _ as delta_kl ->
         let n =
           List.fold_left
             (fun acc arrow -> acc + A.dir k l arrow)
             0 arrows in
         let ghosts =
           List.rev_map
             (fun (norm, arrows) ->
               (Algebra.Laurent.mul (LP.int (-n)) norm, ?? a :: arrows))
             delta_kl in
         List.fold_left
           (fun acc arrows -> insert_gluon a k l arrows @ acc)
           ghosts delta_kl

    let t_of_delta delta a k l =
      canonicalize (t_of_delta delta a k l)

    let t_S n a k l =
      t_of_delta (delta_S n) a k l

    let t_A n a k l =
      t_of_delta (delta_A n) a k l

    let t6 = t_S 2
    let t10 = t_S 3
    let t15 = t_S 4
    let t3bar = t_A 2

(* Equivalent definition: *)
    let t8' a b c =
      t_of_delta delta8 a b c

    let t_of_tableau tableau a k l =
      t_of_delta (delta_of_tableau tableau) a k l

(* \begin{dubious}
     Check the following for a real live UFO file!
   \end{dubious} *)

(* In the UFO paper, the Clebsh-Gordan is defined
   as~$K^{(6),ij}_{\hphantom{(6),ij}m}$.  Therefore, keeping
   our convention for the generators~$T_{a\hphantom{(6),j}i}^{(6),j}$,
   the must arrows \emph{end} at~$m$. *)
    let k6 m i j =
      experimental "k6";
      [ (LP.int 1, [i =>> (m, 0); j =>> (m, 1)]);
        (LP.int 1, [i =>> (m, 1); j =>> (m, 0)]) ]

(* The arrow are reversed for~$\bar K^{(6),m}_{\hphantom{(6),m}ij}$
   and \emph{start} at~$m$. *)
    let k6bar m i j =
      experimental "k6bar";
      [ (LP.int 1, [(m, 0) >=> i; (m, 1) >=> j]);
        (LP.int 1, [(m, 1) >=> i; (m, 0) >=> j]) ]

    (* \thocwmodulesubsection{Unit Tests} *)

    module Test : Test =
      struct

        open OUnit
        module L = Algebra.Laurent

        module B = Birdtracks

        open Birdtracks
        open Birdtracks.Infix

        let exorcise vertex =
          List.filter
            (fun (_, arrows) -> not (List.exists A.is_ghost arrows))
            vertex

        let eqx v1 v2 =
          eq (exorcise v1) (exorcise v2)

(* \thocwmodulesubsection{Trivia} *)

        let suite_sum =
          "sum" >:::

            [ "atoms" >::
                (fun () ->
                  eq
                    (two *** delta3 1 2)
                    (delta3 1 2 +++ delta3 1 2)) ]

        let suite_diff =
          "diff" >:::

            [ "atoms" >::
                (fun () ->
                  eq
                    (delta3 3 4)
                    (delta3 1 2 +++ delta3 3 4 --- delta3 1 2)) ]


(* \begin{equation}
      \prod_{k=i}^j f(k)
   \end{equation} *)
        let rec product f i j =
          if i > j then
            null
          else if i = j then
            f i
          else
            f i *** product f (succ i) j

(* In particular
   \begin{multline}
      \text{[product (nc_minus_n_plus n) i j]}\, \mapsto \\
         \prod_{k=i}^j (N_C-n+k)
          = \frac{(N_C-n+j)!}{(N_C-n+i-1)!}
          = (N_C-n+j)(N_C-n+j-1)\cdots(N_C-n+i)
   \end{multline} *)
        let nc_minus_n_plus n k =
          const (LP.ints [ (1, 1); (-n + k, 0) ])

        let contractions rank k =
          product (nc_minus_n_plus rank) 1 k

        let suite_times =
          "times" >:::

            [ "reorder components t1*t2" >:: (* trivial $T_a^{ik}T_a^{kj}=T_a^{kj}T_a^{ik}$ *)
	        (fun () ->
                  let t1 = t (-1) 1 (-2)
                  and t2 = t (-1) (-2) 2 in
	          eq (t1 *** t2) (t2 *** t1));

              "reorder components tr(t1*t2)" >:: (* trivial $T_a^{ij}T_a^{ji}=T_a^{ji}T_a^{ij}$ *)
	        (fun () ->
                  let t1 = t 1 (-1) (-2)
                  and t2 = t 2 (-2) (-1) in
	          eq (t1 *** t2) (t2 *** t1));

              "reorderings" >::
	        (fun () ->
                  let v1 = [(L.unit, [ 1 => -2; -2 => -1; -1 =>  1])]
                  and v2 = [(L.unit, [-1 =>  2;  2 => -2; -2 => -1])]
                  and v' = [(L.unit, [ 1 =>  1;  2 =>  2])] in
	          eq v' (v1 *** v2));
 
              "eps*epsbar" >::
	        (fun () ->
	          eq
                    (delta3 1 2 *** delta3 3 4 --- delta3 1 4 *** delta3 3 2)
                    (epsilon [1; 3] *** epsilon_bar [2; 4]));
 
              "eps*epsbar -" >::
	        (fun () ->
	          eq
                    (delta3 1 4 *** delta3 3 2 --- delta3 1 2 *** delta3 3 4)
                    (epsilon [1; 3] *** epsilon_bar [4; 2]));
 
              "eps*epsbar 1" >::
	        (fun () ->
	          eq (* $N_C-3+1=(N_C-2)$, for $NC=3$: $1$ *)
                    (contractions 3 1 ***
                       (delta3 1 2 *** delta3 3 4 --- delta3 1 4 *** delta3 3 2))
                    (epsilon [-1; 1; 3] *** epsilon_bar [-1; 2; 4]));
 
              "eps*epsbar cyclic 1" >::
	        (fun () ->
	          eq (* $N_C-3+1=(N_C-2)$, for $NC=3$: $1$ *)
                    (contractions 3 1 ***
                       (delta3 1 2 *** delta3 3 4 --- delta3 1 4 *** delta3 3 2))
                    (epsilon [3; -1; 1] *** epsilon_bar [-1; 2; 4]));
 
              "eps*epsbar cyclic 2" >::
	        (fun () ->
	          eq (* $N_C-3+1=(N_C-2)$, for $NC=3$: $1$ *)
                    (contractions 3 1 ***
                       (delta3 1 2 *** delta3 3 4 --- delta3 1 4 *** delta3 3 2))
                    (epsilon [-1; 1; 3] *** epsilon_bar [4; -1; 2]));
 
              "eps*epsbar 2" >::
	        (fun () ->
	          eq (* $(N_C-3+2)(N_C-3+1)=(N_C-1)(N_C-2)$, for $NC=3$: $2$ *)
                    (contractions 3 2 *** delta3 1 2)
                    (epsilon [-1; -2; 1] *** epsilon_bar [-1; -2; 2]));
 
              "eps*epsbar 3" >::
	        (fun () ->
	          eq (* $(N_C-3+3)(N_C-3+2)(N_C-3+1)=N_C(N_C-1)(N_C-2)$, for $NC=3$: $3!$ *)
                    (contractions 3 3)
                    (epsilon [-1; -2; -3] *** epsilon_bar [-1; -2; -3]));
 
              "eps*epsbar big" >::
	        (fun () ->
	          eq (* $(N_C-5+3)(N_C-5+2)(N_C-5+1)=(N_C-2)(N_C-3)(N_C-4)$, for $NC=5$: $3!$ *)
                    (contractions 5 3 ***
                       (epsilon [4; 5] *** epsilon_bar [6; 7]))
                    (epsilon [-1; -2; -3; 4; 5] *** epsilon_bar [-1; -2; -3; 6; 7]));
 
              "eps*epsbar big -" >::
	        (fun () ->
	          eq (* $(N_C-5+3)(N_C-5+2)(N_C-5+1)=(N_C-2)(N_C-3)(N_C-4)$, for $NC=5$: $3!$ *)
                    (contractions 5 3 ***
                       (epsilon [5; 4] *** epsilon_bar [6; 7]))
                    (epsilon [-1; 4; -3; -2; 5] *** epsilon_bar [-1; -2; -3; 6; 7])) ]

(* \thocwmodulesubsection{Propagators} *)

(* Verify the normalization of the propagators by making sure
   that $D^{ij}D^{jk}=D^{ik}$ *)
        let projection_id rep_d =
	  eq (rep_d 1 2) (rep_d 1 (-1) *** rep_d (-1) 2)

        let orthogonality d d' =
          assert_zero_vertex (d 1 (-1) *** d' (-1) 2)

(* Pass every arrow straight through, without (anti-)symmetrization. *)
        let delta_unsymmetrized n k l =
          delta_of_permutations n [(1, ThoList.range 0 (pred n))] k l

        let completeness n tableaux =
          eq
            (delta_unsymmetrized n 1 2)
            (sum (List.map (fun t -> delta_of_tableau t 1 2) tableaux))

(* The following names are of historical origin. From the time,
   when we didn't have full support for Young tableaux and
   implemented figure 9.1 from the birdtrack book.
   \ytableausetup{centertableaux,smalltableaux}
   \begin{equation}
     \ytableaushort{01,2}
   \end{equation} *)

        let delta_SAS i j =
          delta_of_tableau [[0;1];[2]] i j

(* \begin{equation}
     \ytableaushort{02,1}
   \end{equation} *)

        let delta_ASA i j =
          delta_of_tableau [[0;2];[1]] i j

        let suite_propagators =
          "propagators" >:::
            [ "D*D=D" >:: (fun () -> projection_id delta3);
              "D8*D8=D8" >:: (fun () -> projection_id delta8);
              "G*G=G" >:: (fun () -> projection_id gluon);
              "D6*D6=D6" >:: (fun () -> projection_id delta6);
              "D10*D10=D10" >:: (fun () -> projection_id delta10);
              "D15*D15=D15" >:: (fun () -> projection_id delta15);
              "D3bar*D3bar=D3bar" >:: (fun () -> projection_id delta3bar);
              "D6*D3bar=0" >:: (fun () -> orthogonality delta6 delta3bar);
              "D_A3*D_A3=D_A3" >:: (fun () -> projection_id (delta_A 3));
              "D10*D_A3=0" >:: (fun () -> orthogonality delta10 (delta_A 3));
              "D_SAS*D_SAS=D_SAS" >:: (fun () -> projection_id delta_SAS);
              "D_ASA*D_ASA=D_ASA" >:: (fun () -> projection_id delta_ASA);
              "D_SAS*D_S3=0" >:: (fun () -> orthogonality delta_SAS (delta_S 3));
              "D_SAS*D_A3=0" >:: (fun () -> orthogonality delta_SAS (delta_A 3));
              "D_SAS*D_ASA=0" >:: (fun () -> orthogonality delta_SAS delta_ASA);
              "D_ASA*D_SAS=0" >:: (fun () -> orthogonality delta_ASA delta_SAS);
              "D_ASA*D_S3=0" >:: (fun () -> orthogonality delta_ASA (delta_S 3));
              "D_ASA*D_A3=0" >:: (fun () -> orthogonality delta_ASA (delta_A 3));
              "DU*DU=DU" >:: (fun () -> projection_id (delta_unsymmetrized 3));

              "S3=[0123]" >::
                (fun () ->
                  eq (delta_S 4 1 2) (delta_of_tableau [[0;1;2;3]] 1 2));

              "A3=[0,1,2,3]" >::
                (fun () ->
                  eq (delta_A 4 1 2) (delta_of_tableau [[0];[1];[2];[3]] 1 2));

              "[0123]*[012,3]=0" >::
                (fun () ->
                  orthogonality
                    (delta_of_tableau [[0;1;2;3]])
                    (delta_of_tableau [[0;1;2];[3]]));

              "[0123]*[01,23]=0" >::
                (fun () ->
                  orthogonality
                    (delta_of_tableau [[0;1;2;3]])
                    (delta_of_tableau [[0;1];[2;3]]));

              "[012,3]*[012,3]=[012,3]" >::
                (fun () -> projection_id (delta_of_tableau [[0;1;2];[3]]));

(* \ytableausetup{centertableaux,smalltableaux}
   \begin{equation}
                       \ytableaushort{01}
     +                 \ytableaushort{0,1}
   \end{equation} *)

              "completeness 2" >:: (fun () -> completeness 2 [ [[0;1]]; [[0];[1]] ]) ;

              "completeness 2'" >::
                (fun () ->
                  eq
                    (delta_unsymmetrized 2 1 2)
                    (delta_S 2 1 2 +++ delta_A 2 1 2));

(* The normalization factors are written for illustration.  They are
   added by [delta_of_tableau] automatically.
   \ytableausetup{centertableaux,smalltableaux}
   \begin{equation}
                       \ytableaushort{012}
     + \frac{4}{3}\cdot\ytableaushort{01,2}
     + \frac{4}{3}\cdot\ytableaushort{02,1}
     +                 \ytableaushort{0,1,2}
   \end{equation} *)

              "completeness 3" >::
                (fun () -> completeness 3 [ [[0;1;2]]; [[0;1];[2]]; [[0;2];[1]]; [[0];[1];[2]] ]);

              "completeness 3'" >::
                (fun () ->
                  eq
                    (delta_unsymmetrized 3 1 2)
                    (delta_S 3 1 2 +++ delta_SAS 1 2 +++ delta_ASA 1 2 +++ delta_A 3 1 2));

(* \ytableausetup{centertableaux,smalltableaux}
   \begin{equation}
                       \ytableaushort{0123}
     + \frac{3}{2}\cdot\ytableaushort{012,3}
     + \frac{3}{2}\cdot\ytableaushort{013,2}
     + \frac{3}{2}\cdot\ytableaushort{023,1}
     + \frac{4}{3}\cdot\ytableaushort{01,23}
     + \frac{4}{3}\cdot\ytableaushort{02,13}
     + \frac{3}{2}\cdot\ytableaushort{01,2,3}
     + \frac{3}{2}\cdot\ytableaushort{02,1,3}
     + \frac{3}{2}\cdot\ytableaushort{03,1,2}
     +                 \ytableaushort{0,1,2,3}
   \end{equation} *)

              "completeness 4" >::
                (fun () ->
                  completeness 4
                    [ [[0;1;2;3]];
                      [[0;1;2];[3]]; [[0;1;3];[2]]; [[0;2;3];[1]];
                      [[0;1];[2;3]]; [[0;2];[1;3]];
                      [[0;1];[2];[3]]; [[0;2];[1];[3]]; [[0;3];[1];[2]];
                      [[0];[1];[2];[3]] ]) ]

(* \thocwmodulesubsection{Normalization} *)

        let suite_normalization =
          "normalization" >:::

            [ "tr(t*t)" >:: (* $\tr(T_aT_b)=\delta_{ab} + \text{ghosts}$ *)
	        (fun () ->
	          eq
                    (delta8_loop 1 2)
                    (t 1 (-1) (-2) *** t 2 (-2) (-1)));

              "tr(t*t) sans ghosts" >:: (* $\tr(T_aT_b)=\delta_{ab}$ *)
	        (fun () ->
	          eqx
                    (delta8 1 2)
                    (t 1 (-1) (-2) *** t 2 (-2) (-1)));

(* The additional ghostly terms were unexpected, but 
   arises like~(6.2) in our color flow paper~\cite{Kilian:2012pz}. *)
              "t*t*t" >:: (* $T_aT_bT_a=-T_b/N_C + \ldots$ *)
	        (fun () ->
	          eq
                    (minus *** over_nc *** t 1 2 3
                     +++ [(LP.int 1, [1 => 1; 3 => 2]);
                          (LP.nc (-1), [3 => 2; ?? 1])])
                    (t (-1) 2 (-2) *** t 1 (-2) (-3) *** t (-1) (-3) 3));

(* As expected, these ghostly terms cancel in the summed squares
   \begin{equation}
     \tr(T_aT_bT_aT_cT_bT_c)
       = \tr(T_bT_b)/N_C^2
       = \delta_{bb}/N_C^2
       = (N_C^2-1) / N_C^2
       = 1 - 1 / N_C^2
   \end{equation} *)
              "sum((t*t*t)^2)" >:: 
	        (fun () ->
	          eq
                    (ints [(1, 0); (-1, -2)])
                    (t (-1) (-11) (-12) *** t (-2) (-12) (-13) *** t (-1) (-13) (-14)
                     *** t (-3) (-14) (-15) *** t (-2) (-15) (-16) *** t (-3) (-16) (-11)));

              "d*d" >::
                (fun () ->
                  eqx
                    [ (LP.ints [(2, 1); (-8,-1)], 1 <=> 2);
                      (LP.ints [(2, 0); ( 4,-2)], [1=>1; 2=>2]) ]
                    (d 1 (-1) (-2) *** d 2 (-2) (-1))) ]


(* As proposed in our color flow paper~\cite{Kilian:2012pz},
   we can get the correct (anti-)symmetrized generators
   by sandwiching the following unsymmetrized generators
   between the corresponding (anti-)symmetrized projectors.
   Therefore, the unsymmetrized generators work as long as
   they're used in Feynman diagrams, where they are connected
   by propagators that contain (anti-)symmetrized projectors.
   They even work in the Lie algebra relations and give the
   correct normalization there.

   They fail however for more general color algebra expressions
   that can appear in UFO files.
   In particular, the Casimir operators come out really wrong. *)
        let t_unsymmetrized n k l =
          t_of_delta (delta_unsymmetrized n) k l

(* The following trivial vertices are \emph{not} used anymore,
   since they don't get the normalization of the Ward identities
   right.  For the quadratic casimir operators, they always produce a
   result proportional to~$C_F=C_2(S_1)$.  This can be understood because
   they correspond to a fundamental representation with spectators.

   (Anti-)symmetrizing by sandwiching with projectors almost works,
   but they must be multiplied by hand by the number of arrows to get the
   normalization right.
   They're here just for documenting what doesn't work. *)
        let t_trivial n a k l =
          let sterile =
            List.map (fun i -> (l, i) >=>> (k, i)) (ThoList.range 1 (pred n)) in
          [ (LP.int ( 1), ((l, 0) >=> a) :: (a =>> (k, 0)) :: sterile);
            (LP.int (-1), (?? a) :: ((l, 0) >=>> (k, 0)) :: sterile) ]

        let t6_trivial = t_trivial 2
        let t10_trivial = t_trivial 3
        let t15_trivial = t_trivial 4

        let t_SAS = t_of_delta delta_SAS
        let t_ASA = t_of_delta delta_ASA

        let symmetrization ?rep_ts rep_tu rep_d =
          let rep_ts =
            match rep_ts with
            | None -> rep_tu
            | Some rep_t -> rep_t in
          eq
            (rep_ts 1 2 3)
            (gluon 1 (-1) *** rep_d 2 (-2) *** rep_tu (-1) (-2) (-3) *** rep_d (-3) 3)

	let suite_symmetrization =
          "symmetrization" >:::

            [ "t6" >:: (fun () -> symmetrization t6 delta6);
              "t10" >:: (fun () -> symmetrization t10 delta10);
              "t15" >:: (fun () -> symmetrization t15 delta15);
              "t3bar" >:: (fun () -> symmetrization t3bar delta3bar);
              "t_SAS" >:: (fun () -> symmetrization t_SAS delta_SAS);
              "t_ASA" >:: (fun () -> symmetrization t_ASA delta_ASA);
              "t6'" >:: (fun () -> symmetrization ~rep_ts:t6 (t_unsymmetrized 2) delta6);
              "t10'" >:: (fun () -> symmetrization ~rep_ts:t10 (t_unsymmetrized 3) delta10);
              "t15'" >:: (fun () -> symmetrization ~rep_ts:t15 (t_unsymmetrized 4) delta15);

              "t6''" >::
                (fun () ->
                  eq
                    (t6 1 2 3)
                    (int 2 *** delta6 2 (-1) *** t6_trivial 1 (-1) (-2) *** delta6 (-2) 3));

              "t10''" >::
                (fun () ->
                  eq
                    (t10 1 2 3)
                    (int 3 *** delta10 2 (-1) *** t10_trivial 1 (-1) (-2) *** delta10 (-2) 3));

              "t15''" >::
                (fun () ->
                  eq
                    (t15 1 2 3)
                    (int 4 *** delta15 2 (-1) *** t15_trivial 1 (-1) (-2) *** delta15 (-2) 3)) ]

(* \thocwmodulesubsection{Traces} *)

(* Compute (anti-)commutators of generators in the representation~$r$,
   i.\,e.~$[r(t_a)r(t_b)]_{ij}\mp[r(t_b)r(t_a)]_{ij}$, using
   [isum<0] as summation index in the matrix products. *)
        let commutator rep_t i_sum a b i j =
          multiply [rep_t a i i_sum; rep_t b i_sum j]
          --- multiply [rep_t b i i_sum; rep_t a i_sum j]

        let anti_commutator rep_t i_sum a b i j =
          multiply [rep_t a i i_sum; rep_t b i_sum j]
          +++ multiply [rep_t b i i_sum; rep_t a i_sum j]

(* Trace of the product of three generators in the representation~$r$,
   i.\,e.~$\tr_r(r(t_a)r(t_b)r(t_c))$, using $-1,-2,-3$ as summation indices
   in the matrix products. *)
        let trace3 rep_t a b c =
          rep_t a (-1) (-2) *** rep_t b (-2) (-3) *** rep_t c (-3) (-1)

        let loop3 a b c =
          [ (LP.int 1, A.cycle (List.rev [a; b; c]));
            (LP.int (-1), (a <=> b) @ [?? c]);
            (LP.int (-1), (b <=> c) @ [?? a]);
            (LP.int (-1), (c <=> a) @ [?? b]);
            (LP.int 1, [a => a; ?? b; ?? c]);
            (LP.int 1, [?? a; b => b; ?? c]);
            (LP.int 1, [?? a; ?? b; c => c]);
            (LP.nc (-1), [?? a; ?? b; ?? c]) ]

        let suite_trace =
          "trace" >:::

            [ "tr(ttt)" >::
                (fun () -> eq (trace3 t 1 2 3) (loop3 1 2 3));

              "tr(ttt) cyclic 1" >:: (* $\tr(T_aT_bT_c)=\tr(T_bT_cT_a)$ *)
                (fun () -> eq (trace3 t 1 2 3) (trace3 t 2 3 1));

              "tr(ttt) cyclic 2" >:: (* $\tr(T_aT_bT_c)=\tr(T_cT_aT_b)$ *)
                (fun () -> eq (trace3 t 1 2 3) (trace3 t 3 1 2));

(* \begin{dubious}
     Do we expect this?
   \end{dubious} *)
              "tr(tttt)" >:: (* $\tr(T_aT_bT_cT_d)=\ldots$ *)
                (fun () ->
                  eqx
                    [(LP.int 1, A.cycle [4; 3; 2; 1])]
                    (t 1 (-1) (-2) *** t 2 (-2) (-3) *** t 3 (-3) (-4) *** t 4 (-4) (-1))) ]

        let suite_ghosts =
          "ghosts" >:::

            [ "H->gg" >::
	        (fun () ->
	          eq
                    (delta8_loop 1 2)
                    (t 1 (-1) (-2) *** t 2 (-2) (-1)));

              "H->ggg f" >::
	        (fun () ->
	          eq
                    (imag *** f 1 2 3)
                    (trace3 t 1 2 3 --- trace3 t 1 3 2));

              "H->ggg d" >::
	        (fun () ->
	          eq
                    (d 1 2 3)
                    (trace3 t 1 2 3 +++ trace3 t 1 3 2));

              "H->ggg f'" >::
	        (fun () ->
	          eq
                    (imag *** f 1 2 3)
                    (t 1 (-3) (-2) *** commutator t (-1) 2 3 (-2) (-3)));

              "H->ggg d'" >::
	        (fun () ->
	          eq
                    (d 1 2 3)
                    (t 1 (-3) (-2) *** anti_commutator t (-1) 2 3 (-2) (-3)));

              "H->ggg cyclic'" >::
	        (fun () ->
                  let trace a b c =
                    t a (-3) (-2) *** commutator t (-1) b c (-2) (-3) in
	          eq (trace 1 2 3) (trace 2 3 1)) ]

        let ff a1 a2 a3 a4 =
          [ (LP.int (-1), A.cycle [a1; a2; a3; a4]);
            (LP.int ( 1), A.cycle [a2; a1; a3; a4]);
            (LP.int ( 1), A.cycle [a1; a2; a4; a3]);
            (LP.int (-1), A.cycle [a2; a1; a4; a3]) ]

        let tf j i a b =
          [ (LP.imag ( 1), A.chain [i; a; b; j]);
            (LP.imag (-1), A.chain [i; b; a; j]) ]

        let suite_ff =
          "f*f" >:::
            [ "1" >:: (fun () -> eq (ff 1 2 3 4) (f (-1) 1 2 *** f (-1) 3 4));
              "2" >:: (fun () -> eq (ff 1 2 3 4) (f (-1) 1 2 *** f 3 4 (-1)));
              "3" >:: (fun () -> eq (ff 1 2 3 4) (f (-1) 1 2 *** f 4 (-1) 3)) ]

        let suite_tf =
          "t*f" >:::
            [ "1" >:: (fun () -> eq (tf 1 2 3 4) (t (-1) 1 2 *** f (-1) 3 4)) ]

(* \thocwmodulesubsection{Completeness Relation} *)

(* Check the completeness relation corresponding
   to $q\bar q$-scattering:
   \begin{equation}
     \parbox{38\unitlength}{%
       \fmfframe(4,2)(4,4){%
       \begin{fmfgraph*}(30,20)
         \setupFourAmp
         \fmflabel{$i$}{i2}
         \fmflabel{$j$}{i1}
         \fmflabel{$k$}{o1}
         \fmflabel{$l$}{o2}
         \fmf{fermion}{i1,v1,i2}
         \fmf{fermion}{o2,v2,o1}
         \fmf{gluon}{v1,v2}
       \end{fmfgraph*}}} =
     \parbox{38\unitlength}{%
       \fmfframe(4,2)(4,4){%
       \begin{fmfgraph*}(30,20)
         \setupFourAmp
         \fmflabel{$i$}{i2}
         \fmflabel{$j$}{i1}
         \fmflabel{$k$}{o1}
         \fmflabel{$l$}{o2}
         \fmfi{phantom_arrow}{vpath (__i1, __v1)}
         \fmfi{phantom_arrow}{vpath (__v1, __v2) sideways -thick}
         \fmfi{phantom_arrow}{vpath (__v2, __o1)}
         \fmfi{phantom_arrow}{vpath (__o2, __v2)}
         \fmfi{phantom_arrow}{reverse vpath (__v1, __v2) sideways -thick}
         \fmfi{phantom_arrow}{vpath (__v1, __i2)}
         \fmfi{plain}{vpath (__i1, __v1) join 
                      (vpath (__v1, __v2) sideways -thick) join
                      vpath (__v2, __o1)}
         \fmfi{plain}{vpath (__o2, __v2) join
                      (reverse vpath (__v1, __v2) sideways -thick) join
                      vpath (__v1, __i2)}
       \end{fmfgraph*}}} +
     \parbox{38\unitlength}{%
       \fmfframe(4,2)(4,4){%
       \begin{fmfgraph*}(30,20)
         \setupFourAmp
         \fmflabel{$i$}{i2}
         \fmflabel{$j$}{i1}
         \fmflabel{$k$}{o1}
         \fmflabel{$l$}{o2}
         \fmfi{phantom_arrow}{vpath (__i1, __v1)}
         \fmfi{phantom_arrow}{vpath (__v2, __o1)}
         \fmfi{phantom_arrow}{vpath (__o2, __v2)}
         \fmfi{phantom_arrow}{vpath (__v1, __i2)}
         \fmfi{plain}{vpath (__i1, __v1) join 
                      vpath (__v1, __i2)}
         \fmfi{plain}{vpath (__o2, __v2) join
                      vpath (__v2, __o1)}
         \fmfi{dots,label=$-1/N_C$}{vpath (__v1, __v2)}
       \end{fmfgraph*}}}
     \end{equation} *)

        (* $T_{a}^{ij} T_{a}^{kl}$ *)
        let tt i j k l =
          t (-1) i j *** t (-1) k l

        (* $ \delta^{il}\delta^{kj} - \delta^{ij}\delta^{kl}/N_C$ *)
        let tt_expected i j k l =
          [ (LP.int 1, [l => i; j => k]);
            (LP.over_nc (-1), [j => i; l => k]) ]

        let suite_tt =
          "t*t" >:::
            [ "1" >:: (* $T_{a}^{ij} T_{a}^{kl} = \delta^{il}\delta^{kj} - \delta^{ij}\delta^{kl}/N_C$ *)
	        (fun () -> eq (tt_expected 1 2 3 4) (tt 1 2 3 4)) ]

(* \thocwmodulesubsection{Lie Algebra} *)

(* Check the commutation relations $[T_a,T_b]=\ii f_{abc} T_c$
   in various representations. *)
        let lie_algebra_id rep_t =
          let lhs = imag *** f 1 2 (-1) *** t (-1) 3 4
          and rhs = commutator t (-1) 1 2 3 4 in
          eq lhs rhs

(* Check the normalization of the structure consistants
   $\mathcal{N} f_{abc} = - \ii \tr(T_a[T_b,T_c])$ *)
	let f_of_rep_id norm rep_t =
          let lhs = norm *** f 1 2 3
          and rhs = f_of_rep rep_t 1 2 3 in
          eq lhs rhs

(* \begin{dubious}
     Are the normalization factors for the traces of the higher dimensional
     representations correct?
   \end{dubious} *)
(* \begin{dubious}
     The traces don't work for the symmetrized generators
     that we need elsewhere!
   \end{dubious} *)
        let suite_lie =
          "Lie algebra relations" >:::
            [ "[t,t]=ift" >:: (fun () -> lie_algebra_id t);
              "[t8,t8]=ift8" >:: (fun () -> lie_algebra_id t8);
              "[t6,t6]=ift6" >:: (fun () -> lie_algebra_id t6);
              "[t10,t10]=ift10" >:: (fun () -> lie_algebra_id t10);
              "[t15,t15]=ift15" >:: (fun () -> lie_algebra_id t15);
              "[t3bar,t3bar]=ift3bar" >:: (fun () -> lie_algebra_id t3bar);
              "[tSAS,tSAS]=iftSAS" >:: (fun () -> lie_algebra_id t_SAS);
              "[tASA,tASA]=iftASA" >:: (fun () -> lie_algebra_id t_ASA);
              "[t6,t6]=ift6'" >:: (fun () -> lie_algebra_id (t_unsymmetrized 2));
              "[t10,t10]=ift10'" >:: (fun () -> lie_algebra_id (t_unsymmetrized 3));
              "[t15,t15]=ift15'" >:: (fun () -> lie_algebra_id (t_unsymmetrized 4));
              "[t6,t6]=ift6''" >:: (fun () -> lie_algebra_id t6_trivial);
              "[t10,t10]=ift10''" >:: (fun () -> lie_algebra_id t10_trivial);
              "[t15,t15]=ift15''" >:: (fun () -> lie_algebra_id t15_trivial);
              "if = tr(t[t,t])" >:: (fun () -> f_of_rep_id one t);
              "2n*if = tr(t8[t8,t8])" >:: (fun () -> f_of_rep_id (two *** nc) t8);
              "n*if = tr(t6[t6,t6])" >:: (fun () -> f_of_rep_id nc t6_trivial);
              "n^2*if = tr(t10[t10,t10])" >:: (fun () -> f_of_rep_id (nc *** nc) t10_trivial);
              "n^3*if = tr(t15[t15,t15])" >:: (fun () -> f_of_rep_id (nc *** nc *** nc) t15_trivial) ]

(* \thocwmodulesubsection{Ward Identities} *)

(* Testing the color part of basic Ward identities is essentially
   the same as testing the Lie algebra equations above, but with
   generators sandwiched between propagators, as in Feynman diagrams,
   where the relative signs come from the kinematic part of the
   diagrams after applying the equations of motion..   *)

        (* First the diagram with the three gluon vertex
           $\ii f_{abc} D_{cd}^{\text{gluon}} D^{ik} T_d^{kl} D^{lj}$ *)
        let ward_ft rep_t rep_d a b i j =
          imag *** f a b (-11) *** gluon (-11) (-12)
          *** rep_d i (-1) *** rep_t (-12) (-1) (-2) *** rep_d (-2) j

        (* then one diagram with two gauge couplings
           $D^{ik} T_c^{kl} D^{lm} T_c^{mn} D^{nj}$ *)
        let ward_tt1 rep_t rep_d a b i j =
          rep_d i (-1) *** rep_t a (-1) (-2) *** rep_d (-2) (-3)
          *** rep_t b (-3) (-4) *** rep_d (-4) j

        (* finally the difference of exchanged orders:
           $D^{ik} T_a^{kl} D^{lm} T_b^{mn} D^{nj}
           -D^{ik} T_b^{kl} D^{lm} T_a^{mn} D^{nj}$ *)
        let ward_tt rep_t rep_d a b i j =
          ward_tt1 rep_t rep_d a b i j --- ward_tt1 rep_t rep_d b a i j

        (* \begin{dubious}
             The optional [~fudge] factor was used for
             debugging normalizations.
           \end{dubious} *)
        let ward_id ?(fudge=one) rep_t rep_d =
          let lhs = ward_ft rep_t rep_d 1 2 3 4
          and rhs = ward_tt rep_t rep_d 1 2 3 4 in
          eq lhs (fudge *** rhs)

        let suite_ward =
          "Ward identities" >:::
            [ "fund." >:: (fun () -> ward_id t delta3);
              "adj." >:: (fun () -> ward_id t8 delta8);
              "S2" >:: (fun () -> ward_id t6 delta6);
              "S3" >:: (fun () -> ward_id t10 delta10);
              "A2" >:: (fun () -> ward_id t3bar delta3bar);
              "A3" >:: (fun () -> ward_id (t_A 3) (delta_A 3));
              "SAS" >:: (fun () -> ward_id t_SAS delta_SAS);
              "ASA" >:: (fun () -> ward_id t_ASA delta_ASA);
              "S2'" >:: (fun () -> ward_id ~fudge:two t6_trivial delta6);
              "S3'" >:: (fun () -> ward_id ~fudge:(int 3) t10_trivial delta10) ]

        let suite_ward_long =
          "Ward identities" >:::
            [ "S4" >:: (fun () -> ward_id t15 delta15);
              "S4'" >:: (fun () -> ward_id ~fudge:(int 4) t15_trivial delta15) ]

(* \thocwmodulesubsection{Jacobi Identities} *)

        (* $T_aT_bT_c$ *)
        let prod3 rep_t a b c i j =
          rep_t a i (-1) *** rep_t b (-1) (-2) *** rep_t c (-2) j

        (* $[T_a,[T_b,T_c]]$ *)
        let jacobi1 rep_t a b c i j =
          (prod3 rep_t a b c i j --- prod3 rep_t a c b i j)
          --- (prod3 rep_t b c a i j --- prod3 rep_t c b a i j)

        (* sum of cyclic permutations of $[T_a,[T_b,T_c]]$ *)
        let jacobi rep_t =
          sum [jacobi1 rep_t 1 2 3 4 5;
               jacobi1 rep_t 2 3 1 4 5;
               jacobi1 rep_t 3 1 2 4 5]

        let jacobi_id rep_t =
          assert_zero_vertex (jacobi rep_t)

        let suite_jacobi =
          "Jacobi identities" >:::
            [ "fund." >:: (fun () -> jacobi_id t);
              "adj." >:: (fun () -> jacobi_id f);
              "S2" >:: (fun () -> jacobi_id t6);
              "S3" >:: (fun () -> jacobi_id t10);
              "A2" >:: (fun () -> jacobi_id (t_A 2));
              "A3" >:: (fun () -> jacobi_id (t_A 3));
              "SAS" >:: (fun () -> jacobi_id t_SAS);
              "ASA" >:: (fun () -> jacobi_id t_ASA);
              "S2'" >:: (fun () -> jacobi_id t6_trivial);
              "S3'" >:: (fun () -> jacobi_id t10_trivial) ]

        let suite_jacobi_long =
          "Jacobi identities" >:::
            [ "S4" >:: (fun () -> jacobi_id t15);
              "S4'" >:: (fun () -> jacobi_id t15_trivial) ]

(* \thocwmodulesubsection{Casimir Operators}
   \label{pg:casimir-tests} *)

        (* We can read of the eigenvalues of the Casimir operators for
           the adjoint, totally symmetric and totally antisymmetric
           representations of~$\mathrm{SU}(N)$ from table~II of
           \texttt{hep-ph/0611341}
           \begin{subequations}
             \begin{align}
               C_2(\text{adj}) &= 2N \\
               C_2(S_n) &= \frac{n(N-1)(N+n)}{N} \\
               C_2(A_n) &= \frac{n(N-n)(N+1)}{N}
          \end{align}
           \end{subequations}
           adjusted for our normalization.
           Also from \texttt{arxiv:1912.13302}
           \begin{equation}
               C_3(S_1) =(N^2-1)(N^2-4)/N^2=\frac{N_C^4-5N_C^2+4}{N_C^2}
           \end{equation} *)

        (* Building blocks $n/N_C$ and $N_C+n$ *)
        let n_over_nc n = const (LP.ints [ (n, -1) ])
        let nc_plus n = const (LP.ints [ (1, 1); (n,0) ])

        (* $C_2(S_n) = n/N_C(N_C-1)(N_C+n)$ *)
        let c2_S n = n_over_nc n *** nc_plus (-1) *** nc_plus n

        (* $C_2(A_n) = n/N_C(N_C-n)(N_C+1)$ *)
        let c2_A n = n_over_nc n *** nc_plus (-n) *** nc_plus 1
          
        let casimir_tt i j = c2_S 1 *** delta3 i j
        let casimir_t6t6 i j = c2_S 2 *** delta6 i j
        let casimir_t10t10 i j = c2_S 3 *** delta10 i j
        let casimir_t15t15 i j = c2_S 4 *** delta15 i j
        let casimir_t3bart3bar i j = c2_A 2 *** delta3bar i j
        let casimir_tA3tA3 i j = c2_A 3 *** delta_A 3 i j

        (* $C_2(\text{adj})=2N_C$ *)
        let ca = LP.ints [(2, 1)]
        let casimir_ff a b = [(ca, 1 <=> 2); (LP.int (-2), [1=>1; 2=>2])]

        (* $C_3(S_1)=N_C^2-5+4/N_C^2$ *)
        let c3f = LP.ints [(1, 2); (-5, 0); (4, -2)]
        let casimir_ttt i j = const c3f *** delta3 i j

        let suite_casimir =
          "Casimir operators" >:::

            [ "t*t" >::
	        (fun () ->
	          eq
                    (casimir_tt 1 2)
                    (t (-1) 1 (-2) *** t (-1) (-2) 2));

              "t*t*t" >::
	        (fun () ->
	          eq
                    (casimir_ttt 1 2)
                    (d (-1) (-2) (-3) ***
                       t (-1) 1 (-4) *** t (-2) (-4) (-5) *** t (-3) (-5) 2));

              "f*f" >::
	        (fun () ->
	          eq
                    (casimir_ff 1 2)
                    (minus *** f (-1) 1 (-2) *** f (-1) (-2) 2));

              "t6*t6" >::
	        (fun () ->
	          eq
                    (casimir_t6t6 1 2)
                    (t6 (-1) 1 (-2) *** t6 (-1) (-2) 2));

              "t3bar*t3bar" >::
	        (fun () ->
	          eq
                    (casimir_t3bart3bar 1 2)
                    (t3bar (-1) 1 (-2) *** t3bar (-1) (-2) 2));

              "tA3*tA3" >::
	        (fun () ->
	          eq
                    (casimir_tA3tA3 1 2)
                    (t_A 3 (-1) 1 (-2) *** t_A 3 (-1) (-2) 2));

              "t_SAS*t_SAS" >::
	        (fun () ->
	          eq
                    (const (LP.ints [(3,1); (-9,-1)]) *** delta_SAS 1 2)
                    (t_SAS (-1) 1 (-2) *** t_SAS (-1) (-2) 2));

              "t_ASA*t_ASA" >::
	        (fun () ->
	          eq
                    (const (LP.ints [(3,1); (-9,-1)]) *** delta_ASA 1 2)
                    (t_ASA (-1) 1 (-2) *** t_ASA (-1) (-2) 2));

              "t10*t10" >::
	        (fun () ->
	          eq
                    (casimir_t10t10 1 2)
                    (t10 (-1) 1 (-2) *** t10 (-1) (-2) 2)) ]

        let suite_casimir_long =
          "Casimir operators" >:::

            [ "t15*t15" >::
	        (fun () ->
	          eq
                    (casimir_t15t15 1 2)
                    (t15 (-1) 1 (-2) *** t15 (-1) (-2) 2)) ]

(* \thocwmodulesubsection{Color Sums} *)

        let suite_colorsums =
          "(squared) color sums" >:::

            [ "gluon normalization" >::
	        (fun () ->
	          eq
                    (delta8 1 2)
                    (delta8 1 (-1) *** gluon (-1) (-2) *** delta8 (-2) 2));

              "f*f" >::
	        (fun () ->
                  let sum_ff =
                    multiply [ f (-11) (-12) (-13);
                               f (-21) (-22) (-23);
                               gluon (-11) (-21);
                               gluon (-12) (-22);
                               gluon (-13) (-23) ]
                  and expected = ints [(2, 3); (-2, 1)] in
	          eq expected sum_ff);

              "d*d" >::
	        (fun () ->
                  let sum_dd =
                    multiply [ d (-11) (-12) (-13);
                               d (-21) (-22) (-23);
                               gluon (-11) (-21);
                               gluon (-12) (-22);
                               gluon (-13) (-23) ]
                  and expected = ints [(2, 3); (-10, 1); (8, -1)] in
	          eq expected sum_dd);

              "f*d" >::
	        (fun () ->
                  let sum_fd =
                    multiply [ f (-11) (-12) (-13);
                               d (-21) (-22) (-23);
                               gluon (-11) (-21);
                               gluon (-12) (-22);
                               gluon (-13) (-23) ] in
	          assert_zero_vertex sum_fd);

              "Hgg" >::
	        (fun () ->
                  let sum_hgg =
                    multiply [ delta8_loop (-11) (-12);
                               delta8_loop (-21) (-22);
                               gluon (-11) (-21);
                               gluon (-12) (-22) ]
                  and expected = ints [(1, 2); (-1, 0)] in
	          eq expected sum_hgg) ]

        let suite =
          "Color.SU3" >:::
	    [suite_sum;
	     suite_diff;
	     suite_times;
	     suite_normalization;
	     suite_symmetrization;
	     suite_ghosts;
	     suite_propagators;
	     suite_trace;
	     suite_ff;
	     suite_tf;
	     suite_tt;
             suite_lie;
             suite_ward;
             suite_jacobi;
	     suite_casimir;
             suite_colorsums]

        let suite_long =
          "Color.SU3 long" >:::
	    [suite_ward_long;
             suite_jacobi_long;
             suite_casimir_long]

      end

  end

(* \thocwmodulesection{$\mathrm{U}(N_C)$} *)

(* \begin{dubious}
     This must not be used, because it has not yet been updated
     to the correctly symmetrized version!
   \end{dubious} *)

module U3 : SU3 =
  struct

    module A = Arrow
    open Arrow.Infix

    module B = Birdtracks
    type t = B.t
    let canonicalize = B.canonicalize
    let to_string = B.to_string
    let pp = B.pp
    let trivial = B.trivial
    let is_null = B.is_null
    let null = B.null
    let const = B.const
    let one = B.one
    let two = B.two
    let int = B.int
    let half = B.half
    let third = B.third
    let fraction = B.fraction
    let nc = B.nc
    let over_nc = B.over_nc
    let minus = B.minus
    let imag = B.imag
    let ints = B.ints
    let sum = B.sum
    let diff = B.diff
    let scale = B.scale
    let times = B.times
    let multiply = B.multiply
    let relocate = B.relocate
    let fuse = B.fuse
    let f_of_rep = B.f_of_rep
    let d_of_rep = B.d_of_rep
    module Infix = B.Infix

    let delta3 i j =
      [(LP.int 1, j ==> i)]

    let delta8 a b =
      [(LP.int 1, a <=> b)]

    let delta8_loop = delta8

    let gluon a b =
      delta8 a b

    let delta6 n m =
      [ (LP.fraction 2, [(m, 0) >=>> (n, 0); (m, 1) >=>> (n, 1)]);
        (LP.fraction 2, [(m, 0) >=>> (n, 1); (m, 1) >=>> (n, 0)]) ]

    let triples =
      [(0, 1, 2); (1, 2, 0); (2, 0, 1);
       (2, 1, 0); (0, 2, 1); (1, 0, 2)]

    let delta10 n m =
      List.map
        (fun (i, j, k) ->
          (LP.fraction 6, [(m, 0) >=>> (n, i);
                           (m, 1) >=>> (n, j);
                           (m, 2) >=>> (n, k)]))
        triples

    let t a i j =
      [ (LP.int 1, [j => a; a => i]) ]

    let f a b c =
      [ (LP.imag ( 1), A.cycle [a; b; c]);
        (LP.imag (-1), A.cycle [a; c; b]) ]

    let t8 a b c =
      Birdtracks.Infix.( minus *** imag *** f a b c )

    let d a b c =
      [ (LP.int 1, A.cycle [a; b; c]);
        (LP.int 1, A.cycle [a; c; b]) ]

    let incomplete tensor =
      failwith ("Color.Vertex: " ^ tensor ^ " not supported yet!")

    let experimental tensor =
      Printf.eprintf
        "Color.Vertex: %s support still experimental and untested!\n"
        tensor

    let epsilon tips = incomplete "epsilon-tensor"
    let epsilon_bar tails = incomplete "epsilon-tensor"

    let t6 a m n =
      [ (LP.int ( 1), [(n, 0) >=> a; a =>> (m, 0); (n, 1) >=>> (m, 1)]);
        (LP.int ( 1), [(n, 1) >=> a; a =>> (m, 0); (n, 0) >=>> (m, 1)]) ]

    let t10 a m n =
      [ (LP.int ( 1), [(n, 0) >=> a; a =>> (m, 0);
                       (n, 1) >=>> (m, 1);
                       (n, 2) >=>> (m, 2)]);
        (LP.int (-1), [(n, 0) >=>> (m, 0);
                       (n, 1) >=>> (m, 1);
                       (n, 2) >=>> (m, 2)]) ]

    let k6 m i j =
      experimental "k6-tensor";
      [ (LP.int 1, [(m, 0) >=> i; (m, 1) >=> j]);
        (LP.int 1, [(m, 1) >=> i; (m, 0) >=> j]) ]

    let k6bar m i j =
      experimental "k6-tensor";
      [ (LP.int 1, [i =>> (m, 0); j =>> (m, 1)]);
        (LP.int 1, [i =>> (m, 1); j =>> (m, 0)]) ]

    let delta_of_tableau t i j =
      incomplete "delta_of_tableau"

    let t_of_tableau tableau a k l =
      incomplete "t_of_tableau"

    (* \thocwmodulesubsection{Unit Tests} *)

    module Test : Test =
      struct

        open OUnit
        open Birdtracks
        open Infix

        let suite_lie =
          "Lie algebra relations" >:::

            [ "if = tr(t[t,t])" >::
	        (fun () -> eq (f 1 2 3) (f_of_rep t 1 2 3)) ]

        (* $N_C=N_C^2/N_C$ *)
        let cf = LP.ints [(1, 1)]

        let casimir_tt i j =
          [(cf, i ==> j)]

        let suite_casimir =
          "Casimir operators" >:::

            [ "t*t" >::
	        (fun () ->
	          eq (casimir_tt 2 1) (t (-1) (-2) 2 *** t (-1) 1 (-2))) ]

        let suite =
          "Color.U3" >:::
	    [suite_lie;
             suite_casimir]

        let suite_long =
          "Color.U3 long" >:::
	    []

      end

  end

module Vertex = SU3
