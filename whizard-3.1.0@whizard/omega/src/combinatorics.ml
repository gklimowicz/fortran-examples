(* combinatorics.ml --

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

type 'a seq = 'a list

(* \thocwmodulesection{Simple Combinatorial Functions} *)

let rec factorial' fn n =
  if n < 1 then
    fn
  else
    factorial' (n * fn) (pred n)

let factorial n =
  let result = factorial' 1 n in
  if result < 0 then
    invalid_arg "Combinatorics.factorial overflow"
  else
    result

(* \begin{multline}
     \binom{n}{k} = \frac{n!}{k!(n-k)!}
         = \frac{n(n-1)\cdots(n-k+1)}{k(k-1)\cdots1} \\
         = \frac{n(n-1)\cdots(k+1)}{(n-k)(n-k-1)\cdots1} =
       \begin{cases}
         B_{n-k+1}(n,k) & \text{for $k \le \lfloor n/2 \rfloor$} \\
         B_{k+1}(n,n-k) & \text{for $k >   \lfloor n/2 \rfloor$}
       \end{cases}
   \end{multline}
   where
   \begin{equation}
     B_{n_{\min}}(n,k) =
       \begin{cases}
         n B_{n_{\min}}(n-1,k)           & \text{for $n \ge n_{\min}$} \\
         \frac{1}{k} B_{n_{\min}}(n,k-1) & \text{for $k > 1$} \\
         1                               & \text{otherwise}
       \end{cases}
   \end{equation} *)

let rec binomial' n_min n k acc =
  if n >= n_min then
    binomial' n_min (pred n) k (n * acc)
  else if k > 1 then
    binomial' n_min n (pred k) (acc / k)
  else
    acc

let binomial n k =
  if k > n / 2 then
    binomial' (k + 1) n (n - k) 1
  else
    binomial' (n - k + 1) n k 1

(* Overflows later, but takes much more time:
   \begin{equation}
     \binom{n}{k} = \binom{n-1}{k} + \binom{n-1}{k-1}
   \end{equation} *)

let rec slow_binomial n k =
  if n < 0 || k < 0 then
    invalid_arg "Combinatorics.binomial"
  else if k = 0 || k = n then
    1
  else
    slow_binomial (pred n) k + slow_binomial (pred n) (pred k)

let multinomial n_list =
  List.fold_left (fun acc n -> acc / (factorial n))
    (factorial (List.fold_left (+) 0 n_list)) n_list

let symmetry l =
  List.fold_left (fun s (n, _) -> s * factorial n) 1 (ThoList.classify l)

(* \thocwmodulesection{Partitions} *)

(* The inner steps of the recursion (i.\,e.~$n=1$) are expanded as follows
   \begin{multline}
     \ocwlowerid{split'}(1,\lbrack p_k;p_{k-1};\ldots;p_1\rbrack,
                           \lbrack x_l;x_{l-1};\ldots;x_1\rbrack,
                           \lbrack x_{l+1};x_{l+2};\ldots;x_m\rbrack ) = \\
       \lbrack (\lbrack p_1;\ldots;p_k;x_{l+1}\rbrack,
                \lbrack x_1;\ldots;x_l;x_{l+2};\ldots;x_m\rbrack); \qquad\qquad\qquad\\
               (\lbrack p_1;\ldots;p_k;x_{l+2}\rbrack,
                \lbrack x_1;\ldots;x_l;x_{l+1};x_{l+3}\ldots;x_m\rbrack);
            \ldots; \\
               (\lbrack p_1;\ldots;p_k;x_m\rbrack,
                \lbrack x_1;\ldots;x_l;x_{l+1};\ldots;x_{m-1}\rbrack) \rbrack
   \end{multline}
   while the outer steps (i.\,e.~$n>1$) perform the same with one element
   moved from the last argument to the first argument.  At the $n$th level we have
   \begin{multline}
     \ocwlowerid{split'}(n,\lbrack p_k;p_{k-1};\ldots;p_1\rbrack,
                           \lbrack x_l;x_{l-1};\ldots;x_1\rbrack,
                           \lbrack x_{l+1};x_{l+2};\ldots;x_m\rbrack ) = \\
       \lbrack (\lbrack p_1;\ldots;p_k;x_{l+1};x_{l+2};\ldots;x_{l+n}\rbrack,
                \lbrack x_1;\ldots;x_l;x_{l+n+1};\ldots;x_m\rbrack); \ldots; \qquad\\
               (\lbrack p_1;\ldots;p_k;x_{m-n+1};x_{m-n+2};\ldots;x_{m}\rbrack,
                \lbrack x_1;\ldots;x_l;x_{l+1};\ldots;x_{m-n}\rbrack) \rbrack
   \end{multline}
   where the order of the~$\lbrack x_1;x_2;\ldots;x_m\rbrack$ is maintained in
   the partitions.  Variations on this multiple recursion idiom are used many
   times below.  *)

let rec split' n rev_part rev_head = function
  | [] -> []
  | x :: tail ->
      let rev_part' = x :: rev_part
      and parts = split' n rev_part (x :: rev_head) tail in
      if n < 1 then
        failwith "Combinatorics.split': can't happen"
      else if n = 1 then
	(List.rev rev_part', List.rev_append rev_head tail) :: parts
      else
	split' (pred n) rev_part' rev_head tail @ parts

(* Kick off the recursion for $0<n<|l|$ and handle the cases $n\in\{0,|l|\}$
   explicitely.  Use reflection symmetry for a small optimization. *)

let ordered_split_unsafe n abs_l l =
  let abs_l = List.length l in
  if n = 0 then
    [[], l]
  else if n = abs_l then
    [l, []]
  else if n <= abs_l / 2 then
    split' n [] [] l
  else
    List.rev_map (fun (a, b) -> (b, a)) (split' (abs_l - n) [] [] l)

(* Check the arguments and call the workhorse: *)

let ordered_split n l =
  let abs_l = List.length l in
  if n < 0 || n > abs_l then
    invalid_arg "Combinatorics.ordered_split"
  else
    ordered_split_unsafe n abs_l l

(* Handle equipartitions specially: *)

let split n l =
  let abs_l = List.length l in
  if n < 0 || n > abs_l then
    invalid_arg "Combinatorics.split"
  else begin
    if 2 * n = abs_l then
      match l with
      | [] -> failwith "Combinatorics.split: can't happen"
      | x :: tail ->
          List.map (fun (p1, p2) -> (x :: p1, p2)) (split' (pred n) [] [] tail)
    else
      ordered_split_unsafe n abs_l l
  end

(* If we chop off parts repeatedly, we can either keep permutations or
   suppress them.  Generically, [attach_to_fst] has type
   \begin{quote}
     [('a * 'b) list -> 'a list -> ('a list * 'b) list -> ('a list * 'b) list]
   \end{quote}
   and semantics
   \begin{multline}
     \ocwlowerid{attach\_to\_fst}
       (\lbrack (a_1,b_1),(a_2,b_2),\ldots,(a_m,b_m)\rbrack,
        \lbrack a'_1,a'_2,\ldots\rbrack) = \\
        \lbrack (\lbrack a_1,a'_1,\ldots\rbrack, b_1),
                (\lbrack a_2,a'_1,\ldots\rbrack, b_2),\ldots,
                (\lbrack a_m,a'_1,\ldots\rbrack, b_m)\rbrack
   \end{multline}
   (where some of the result can be filtered out), assumed to be
   prepended to the final argument. *)

let rec multi_split' attach_to_fst n size splits =
  if n <= 0 then
    splits
  else
    multi_split' attach_to_fst (pred n) size
      (List.fold_left (fun acc (parts, tail) ->
        attach_to_fst (ordered_split size tail) parts acc) [] splits)

let attach_to_fst_unsorted splits parts acc =
  List.fold_left (fun acc' (p, rest) -> (p :: parts, rest) :: acc') acc splits

(* Similarly, if the secod argument is a list of lists: *)

let prepend_to_fst_unsorted splits parts acc =
  List.fold_left (fun acc' (p, rest) -> (p @ parts, rest) :: acc') acc splits

let attach_to_fst_sorted splits parts acc =
  match parts with
  | [] -> List.fold_left (fun acc' (p, rest) -> ([p], rest) :: acc') acc splits
  | p :: _ as parts ->
      List.fold_left (fun acc' (p', rest) ->
        if p' > p then
          (p' :: parts, rest) :: acc'
        else
          acc') acc splits

let multi_split n size l =
  multi_split' attach_to_fst_sorted n size [([], l)]

let ordered_multi_split n size l =
  multi_split' attach_to_fst_unsorted n size [([], l)]

let rec partitions' splits = function
  | [] -> List.map (fun (h, r) -> (List.rev h, r)) splits
  | (1, size) :: more ->
      partitions'
        (List.fold_left (fun acc (parts, rest) ->
          attach_to_fst_unsorted (split size rest) parts acc)
           [] splits) more
  | (n, size) :: more ->
      partitions'
        (List.fold_left (fun acc (parts, rest) ->
          prepend_to_fst_unsorted (multi_split n size rest) parts acc)
           [] splits) more

let partitions multiplicities l =
  if List.fold_left (+) 0 multiplicities <> List.length l then
    invalid_arg "Combinatorics.partitions"
  else
    List.map fst (partitions' [([], l)]
                    (ThoList.classify (List.sort compare multiplicities)))

let rec ordered_partitions' splits = function
  | [] -> List.map (fun (h, r) -> (List.rev h, r)) splits
  | size :: more ->
      ordered_partitions'
        (List.fold_left (fun acc (parts, rest) ->
          attach_to_fst_unsorted (ordered_split size rest) parts acc)
           [] splits) more

let ordered_partitions multiplicities l =
  if List.fold_left (+) 0 multiplicities <> List.length l then
    invalid_arg "Combinatorics.ordered_partitions"
  else
    List.map fst (ordered_partitions' [([], l)] multiplicities)


let hdtl = function
  | [] -> invalid_arg "Combinatorics.hdtl"
  | h :: t -> (h, t)

let factorized_partitions multiplicities l =
  ThoList.factorize (List.map hdtl (partitions multiplicities l))

(* In order to construct keystones (cf.~chapter~\ref{sec:topology}), we
   must eliminate reflectionsc consistently.  For this to work, the lengths
   of the parts \emph{must not} be reordered arbitrarily.  Ordering with
   monotonously fallings lengths would be incorrect however, because
   then some remainders could fake a reflection symmetry and partitions
   would be dropped erroneously.  Therefore we put the longest first and
   order the remaining with rising lengths: *)

let longest_first l =
  match ThoList.classify (List.sort (fun n1 n2 -> compare n2 n1) l) with
  | [] -> []
  | longest :: rest -> longest :: List.rev rest

let keystones multiplicities l =
  if List.fold_left (+) 0 multiplicities <> List.length l then
    invalid_arg "Combinatorics.keystones"
  else
    List.map fst (partitions' [([], l)] (longest_first multiplicities))

let factorized_keystones multiplicities l =
  ThoList.factorize (List.map hdtl (keystones multiplicities l))

(* \thocwmodulesection{Choices} *)

(* The implementation is very similar to [split'], but here we don't
   have to keep track of the complements of the chosen sets. *)

let rec choose' n rev_choice = function
  | [] -> []
  | x :: tail ->
      let rev_choice' = x :: rev_choice
      and choices = choose' n rev_choice tail in
      if n < 1 then
        failwith "Combinatorics.choose': can't happen"
      else if n = 1 then
	List.rev rev_choice' :: choices
      else
	choose' (pred n) rev_choice' tail @ choices

(* [choose n] is equivalent to $(\ocwlowerid{List.map}\,\ocwlowerid{fst})\circ
   (\ocwlowerid{split\_ordered}\,\ocwlowerid{n})$, but more efficient.  *)

let choose n l =
  let abs_l = List.length l in
  if n < 0 then
    invalid_arg "Combinatorics.choose"
  else if n > abs_l then
    []
  else if n = 0 then
    [[]]
  else if n = abs_l then
    [l]
  else
    choose' n [] l

let multi_choose n size l =
  List.map fst (multi_split n size l)

let ordered_multi_choose n size l =
  List.map fst (ordered_multi_split n size l)

(* \thocwmodulesection{Permutations} *)

let rec insert x = function
  | [] -> [[x]]
  | h :: t as l ->
      (x :: l) :: List.rev_map (fun l' -> h :: l') (insert x t)

let permute l =
  List.fold_left (fun acc x -> ThoList.rev_flatmap (insert x) acc) [[]] l

(* \thocwmodulesubsection{Graded Permutations} *)

let rec insert_signed x = function
  | (eps, []) -> [(eps, [x])]
  | (eps, h :: t) -> (eps, x :: h :: t) ::
      (List.map (fun (eps', l') -> (-eps', h :: l')) (insert_signed x (eps, t)))

let rec permute_signed' = function
  | (eps, []) -> [(eps, [])]
  | (eps, h :: t) -> ThoList.flatmap (insert_signed h) (permute_signed' (eps, t))

let permute_signed l =
  permute_signed' (1, l)

(* The following are wasting at most a factor of two and there's probably
   no point in improving on this \ldots *)

let filter_sign s l =
  List.map snd (List.filter (fun (eps, _) -> eps = s) l)

let permute_even l =
  filter_sign 1 (permute_signed l)

let permute_odd l =
  filter_sign (-1) (permute_signed l)

(* \begin{dubious}
     We have a slight inconsistency here:
     [permute [] = [[]]], while
     [permute_cyclic [] = []].
     I don't know if it is worth fixing.
   \end{dubious} *)

let permute_cyclic l =
  let rec permute_cyclic' acc before = function
    | [] -> List.rev acc
    | x :: rest as after ->
       permute_cyclic' ((after @ List.rev before) :: acc) (x :: before) rest
  in
  permute_cyclic' [] [] l

(* Algorithm: toggle the signs and at the end map all signs to $+1$,
   iff the last sign is positive, i.\,e.~there's an odd number of elements. *)
let permute_cyclic_signed l =
  let rec permute_cyclic_signed' eps acc before = function
    | [] ->
       if eps > 0 then
         List.rev_map (fun (_, p) -> (1, p)) acc
       else
         List.rev acc
    | x :: rest as after ->
       let eps' = - eps in
       permute_cyclic_signed' eps' ((eps', after @ List.rev before) :: acc) (x :: before) rest
  in
  permute_cyclic_signed' (-1) [] [] l

(* \thocwmodulesubsection{Tensor Products of Permutations} *)

let permute_tensor ll =
  Product.list (fun l -> l) (List.map permute ll)

let join_signs l =
  let el, pl = List.split l in
  (List.fold_left (fun acc x -> x * acc) 1 el, pl)

let permute_tensor_signed ll =
  Product.list join_signs (List.map permute_signed ll)

let permute_tensor_even l =
  filter_sign 1 (permute_tensor_signed l)

let permute_tensor_odd l =
  filter_sign (-1) (permute_tensor_signed l)

(* \thocwmodulesubsection{Sorting} *)

let insert_inorder_signed order x (eps, l) =
  let rec insert eps' accu = function
    | [] -> (eps * eps', List.rev_append accu [x])
    | h :: t ->
        if order x h = 0 then
          invalid_arg
            "Combinatorics.insert_inorder_signed: identical elements"
        else if order x h < 0 then
          (eps * eps', List.rev_append accu (x :: h :: t))
        else
          insert (-eps') (h::accu) t
  in
  insert 1 [] l

let sort_signed ?(cmp=pcompare) l =
  List.fold_right (insert_inorder_signed cmp) l (1, [])

let sign ?(cmp=pcompare) l =
  let eps, _ = sort_signed ~cmp l in
  eps

let sign2 ?(cmp=pcompare) l =
  let a = Array.of_list l in
  let eps = ref 1 in
  for j = 0 to Array.length a - 1 do
    for i = 0 to j - 1 do
      if cmp a.(i) a.(j) > 0 then
        eps := - !eps
    done
  done;
  !eps

module Test =
  struct

    open OUnit

    let to_string =
      ThoList.to_string (ThoList.to_string string_of_int)

    let assert_equal_perms =
      assert_equal ~printer:to_string

    let count_permutations n =
      let factorial_n = factorial n
      and range = ThoList.range 1 n in
      let sorted = List.sort compare (permute range) in
      (* Verify the count \ldots *)
      assert_equal factorial_n (List.length sorted);
      (* \ldots{} check that they're all different \ldots *)
      assert_equal factorial_n (List.length (ThoList.uniq sorted));
      (* \ldots{} make sure that they a all permutations. *)
      assert_equal_perms
        [range] (ThoList.uniq (List.map (List.sort compare) sorted))

    let suite_permute =
      "permute" >:::
	[ "permute []" >::
	    (fun () ->
              assert_equal_perms [[]] (permute []));
          "permute [1]" >::
	    (fun () ->
              assert_equal_perms [[1]] (permute [1]));
          "permute [1;2;3]" >::
	    (fun () ->
              assert_equal_perms
                [ [2; 3; 1]; [2; 1; 3]; [3; 2; 1];
                  [1; 3; 2]; [1; 2; 3]; [3; 1; 2] ]
                (permute [1; 2; 3]));
          "permute [1;2;3;4]" >::
	    (fun () ->
              assert_equal_perms
                [ [3; 4; 1; 2]; [3; 1; 2; 4]; [3; 1; 4; 2];
                  [4; 3; 1; 2]; [1; 4; 2; 3]; [1; 2; 3; 4];
                  [1; 2; 4; 3]; [4; 1; 2; 3]; [1; 4; 3; 2];
                  [1; 3; 2; 4]; [1; 3; 4; 2]; [4; 1; 3; 2];
                  [3; 4; 2; 1]; [3; 2; 1; 4]; [3; 2; 4; 1];
                  [4; 3; 2; 1]; [2; 4; 1; 3]; [2; 1; 3; 4];
                  [2; 1; 4; 3]; [4; 2; 1; 3]; [2; 4; 3; 1];
                  [2; 3; 1; 4]; [2; 3; 4; 1]; [4; 2; 3; 1] ]
                (permute [1; 2; 3; 4]));
          "count permute 5" >::
            (fun () -> count_permutations 5);
          "count permute 6" >::
            (fun () -> count_permutations 6);
          "count permute 7" >::
            (fun () -> count_permutations 7);
          "count permute 8" >::
            (fun () -> count_permutations 8);
          "cyclic []" >::
	    (fun () ->
              assert_equal_perms [] (permute_cyclic []));
          "cyclic [1]" >::
	    (fun () ->
              assert_equal_perms [[1]] (permute_cyclic [1]));
          "cyclic [1;2;3]" >::
	    (fun () ->
	      assert_equal_perms
                [[1;2;3]; [2;3;1]; [3;1;2]]
                (permute_cyclic [1;2;3]));
          "cyclic [1;2;3;4]" >::
	    (fun () ->
	      assert_equal_perms
                [[1;2;3;4]; [2;3;4;1]; [3;4;1;2]; [4;1;2;3]]
                (permute_cyclic [1;2;3;4]));
          "cyclic [1;2;3] signed" >::
	    (fun () ->
	      assert_equal
                [(1,[1;2;3]); (1,[2;3;1]); (1,[3;1;2])]
                (permute_cyclic_signed [1;2;3]));
          "cyclic [1;2;3;4] signed" >::
	    (fun () ->
	      assert_equal
                [(1,[1;2;3;4]); (-1,[2;3;4;1]); (1,[3;4;1;2]); (-1,[4;1;2;3])]
                (permute_cyclic_signed [1;2;3;4]))]

    let sort_signed_not_unique =
      "not unique" >::
	(fun () ->
	  assert_raises
            (Invalid_argument
               "Combinatorics.insert_inorder_signed: identical elements")
            (fun () -> sort_signed [1;2;3;4;2]))
        
    let sort_signed_even =
      "even" >::
	(fun () ->
	  assert_equal (1, [1;2;3;4;5;6])
            (sort_signed [1;2;4;3;6;5]))

    let sort_signed_odd =
      "odd" >::
	(fun () ->
	  assert_equal (-1, [1;2;3;4;5;6])
            (sort_signed [2;3;1;5;4;6]))

    let sort_signed_all =
      "all" >::
      (fun () ->
        let l = ThoList.range 1 8 in
        assert_bool "all signed permutations"
          (List.for_all
             (fun (eps, p) ->
               let eps', p' = sort_signed p in
               eps' = eps && p' = l)
             (permute_signed l)))

    let sign_sign2 =
      "sign/sign2" >::
      (fun () ->
        let l = ThoList.range 1 8 in
          assert_bool "all permutations"
          (List.for_all
             (fun p -> sign p = sign2 p)
             (permute l)))

    let suite_sort_signed =
      "sort_signed" >:::
	[sort_signed_not_unique;
         sort_signed_even;
         sort_signed_odd;
         sort_signed_all;
         sign_sign2]

    let suite =
      "Combinatorics" >:::
	[suite_permute;
         suite_sort_signed]

  end

(*i
 *  Local Variables:
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
