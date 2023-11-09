(* linalg.ml --

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

(* This is not a functional implementations, but uses imperative
   array in Fotran style for maximimum speed. *)

exception Singular
exception Not_Square

let copy_matrix a =
  Array.init (Array.length a)
    (fun i -> Array.copy a.(i))

let matmul a b =
  let ni = Array.length a
  and nj = Array.length b.(0)
  and n = Array.length b in
  let ab = Array.make_matrix ni nj 0.0 in
  for i = 0 to pred ni do
    for j = 0 to pred nj do
      for k = 0 to pred n do
	ab.(i).(j) <- ab.(i).(j) +. a.(i).(k) *. b.(k).(j)
      done
    done
  done;
  ab

let matmulv a v =
  let na = Array.length a in
  let nv = Array.length v in
  let v' = Array.make na 0.0 in
  for i = 0 to pred na do
    for j = 0 to pred nv do
      v'.(i) <- v'.(i) +. a.(i).(j) *. v.(j)
    done
  done;
  v'

(*i
let maxval = Array.fold_left max 0.0]

let maxval a : float =
  let x = ref a.(0) in
  for i = 1 to Array.length a - 1 do
    x := max !x a.(i)
  done;
  !x
i*)

let maxabsval a : float =
  let x = ref (abs_float a.(0)) in
  for i = 1 to Array.length a - 1 do
    x := max !x (abs_float a.(i))
  done;
  !x

(*i
let minval = Array.fold_left min 0.0

let minval a : float =
  let x = ref a.(0) in
  for i = 1 to Array.length a - 1 do
    x := min !x a.(i)
  done;
  !x

let maxloc (a : float array) n =
  let n' = ref n
  and max_a : float ref = ref a.(n) in
  for i = succ n to Array.length a - 1 do
    let a_i = a.(i) in
    if a_i > !max_a then begin
      n' := i;
      max_a := a_i
    end
  done;
  !n'

let minloc (a : float array) n =
  let n' = ref n
  and min_a : float ref = ref a.(n) in
  for i = succ n to Array.length a - 1 do
    let a_i = a.(i) in
    if a_i < !min_a then begin
      n' := i;
      min_a := a_i
    end
  done;
  !n'

let rec any' f (a : float array) i =
  if i < 0 then
    false
  else if f a.(i) then
    true
  else
    any' f a (pred i)

let any f a = any' f a (Array.length a - 1)
i*)

(* \thocwmodulesection{$LU$ Decomposition}
   \begin{subequations}
   \label{eq:LU}
   \begin{equation}
     A = LU
   \end{equation}
   In more detail
   \begin{multline}
     \begin{pmatrix}
       a_{00} & a_{01} & \ldots & a_{0(n-1)} \\
       a_{10} & a_{11} & \ldots & a_{1(n-1)} \\
       \vdots & \vdots & \vdots & \vdots \\
       a_{(n-1)0} & a_{(n-1)1} & \ldots & a_{(n-1)(n-1)}
     \end{pmatrix}
     = \\
     \begin{pmatrix}
       1      & 0      & \ldots & 0      \\
       l_{10} & 1      & \ldots & 0      \\
       \vdots & \vdots & \vdots & \vdots \\
       l_{(n-1)0} & l_{(n-1)1} & \ldots & 1
     \end{pmatrix}
     \begin{pmatrix}
       u_{00} & u_{01} & \ldots & u_{0(n-1)} \\
       0      & u_{11} & \ldots & u_{1(n-1)} \\
       \vdots & \vdots & \vdots & \vdots \\
       0      & 0      & \ldots & u_{(n-1)(n-1)}
     \end{pmatrix}
   \end{multline}
   \end{subequations}
   Rewriting~(\ref{eq:LU}) in block matrix notation
   \begin{equation}
     \begin{pmatrix}
       a_{00}     & a_{0\cdot} \\
       a_{\cdot0} & A
     \end{pmatrix}
     =
     \begin{pmatrix}
       1          & 0 \\
       l_{\cdot0} & L
     \end{pmatrix}
     \begin{pmatrix}
       u_{00} & u_{0\cdot} \\
       0      & U
     \end{pmatrix}
      =
     \begin{pmatrix}
       u_{00}            & u_{0\cdot} \\
       l_{\cdot0} u_{00} & l_{\cdot0} \otimes u_{0\cdot} + LU
     \end{pmatrix}
   \end{equation}
   we can solve it easily
   \begin{subequations}
   \begin{align}
     u_{00}     &= a_{00} \\
     u_{0\cdot} &= a_{0\cdot} \\
   \label{eq:LU1}
     l_{\cdot0} &= \frac{a_{\cdot0}}{a_{00}} \\
   \label{eq:LU2}
     LU         &= A - \frac{a_{\cdot0} \otimes a_{0\cdot}}{a_{00}}
   \end{align}
   \end{subequations}
   and~(\ref{eq:LU1}) and~(\ref{eq:LU2}) define a simple iterative
   algorithm if we work from the outside in.  It just remains to add
   pivoting. *)

let swap a i j =
  let a_i = a.(i) in
  a.(i) <- a.(j);
  a.(j) <- a_i

let pivot_column v a n =
  let n' = ref n
  and max_va = ref (v.(n) *. (abs_float a.(n).(n))) in
  for i = succ n to Array.length v - 1 do
    let va_i = v.(i) *. (abs_float a.(i).(n)) in
    if va_i > !max_va then begin
      n' := i;
      max_va := va_i
    end
  done;
  !n'

let lu_decompose_in_place a =
  let n = Array.length a in
  let eps = ref 1
  and pivots = Array.make n 0
  and v =
    try
      Array.init n (fun i ->
	let a_i = a.(i) in
	if Array.length a_i <> n then
	  raise Not_Square;
	1.0 /. (maxabsval a_i))
    with
    | Division_by_zero -> raise Singular in
  for i = 0 to pred n do
    let pivot = pivot_column v a i in
    if pivot <> i then begin
      swap a pivot i;
      eps := - !eps;
      v.(pivot) <- v.(i)
    end;
    pivots.(i) <- pivot;
    let inv_a_ii =
      try 1.0 /. a.(i).(i) with Division_by_zero -> raise Singular in
    for j = succ i to pred n do
      a.(j).(i) <- inv_a_ii *. a.(j).(i)
    done;
    for j = succ i to pred n do
      for k = succ i to pred n do
	a.(j).(k) <- a.(j).(k) -. a.(j).(i) *. a.(i).(k)
      done
    done
  done;
  (pivots, !eps)

let lu_decompose_split a pivots =
  let n = Array.length pivots in
  let l = Array.make_matrix n n 0.0 in
  let u = Array.make_matrix n n 0.0 in
  for i = 0 to pred n do
    l.(i).(i) <- 1.0;
    for j = succ i to pred n do
      l.(j).(i) <- a.(j).(i)
    done
  done;
  for i = pred n downto 0 do
    swap l i pivots.(i)
  done;
  for i = 0 to pred n do
    for j = 0 to i do
      u.(j).(i) <- a.(j).(i)
    done
  done;
  (l, u)

let lu_decompose a =
  let a = copy_matrix a in
  let pivots, _ = lu_decompose_in_place a in
  lu_decompose_split a pivots

let lu_backsubstitute a pivots b =
  let n = Array.length a in
  let nonzero = ref (-1) in
  let b = Array.copy b in
  for i = 0 to pred n do
    let ll = pivots.(i) in
    let b_i = ref (b.(ll)) in
    b.(ll) <- b.(i);
    if !nonzero >= 0 then
      for j = !nonzero to pred i do
	b_i := !b_i -. a.(i).(j) *. b.(j)
      done
    else if !b_i <> 0.0 then
      nonzero := i;
    b.(i) <- !b_i
  done;
  for i = pred n downto 0 do
    let b_i = ref (b.(i)) in
    for j = succ i to pred n do
      b_i := !b_i -. a.(i).(j) *. b.(j)
    done;
    b.(i) <- !b_i /. a.(i).(i)
  done;
  b

let solve_destructive a b =
  let pivot, _ = lu_decompose_in_place a in
  lu_backsubstitute a pivot b

let solve_many_destructive a bs =
  let pivot, _ = lu_decompose_in_place a in
  List.map (lu_backsubstitute a pivot) bs

let solve a b =
  solve_destructive (copy_matrix a) b

let solve_many a bs =
  solve_many_destructive (copy_matrix a) bs

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
