(* circe2/histogram.ml --  *)
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

open Printf

type t =
    { n_bins : int;
      n_bins_float : float;
      x_min : float;
      x_max : float;
      x_min_eps : float;
      x_max_eps : float;
      mutable n_underflow : int;
      mutable underflow : float;
      mutable underflow2 : float;
      mutable n_overflow : int;
      mutable overflow : float;
      mutable overflow2 : float;
      n : int array;
      w : float array;
      w2 : float array }

let create n_bins x_min x_max =
  let eps = 100. *. Float.Double.epsilon *. abs_float (x_max -. x_min) in
  { n_bins = n_bins;
    n_bins_float = float n_bins;
    x_min = x_min;
    x_max = x_max;
    x_min_eps = x_min -. eps;
    x_max_eps = x_max +. eps;
    n_underflow = 0;
    underflow = 0.0;
    underflow2 = 0.0;
    n_overflow = 0;
    overflow = 0.0;
    overflow2 = 0.0;
    n = Array.make n_bins 0;
    w = Array.make n_bins 0.0;
    w2 = Array.make n_bins 0.0 }

let record h x f =
  let i =
    truncate
      (floor (h.n_bins_float *. (x -. h.x_min) /. (h.x_max -. h.x_min))) in
  let i =
    if i < 0 && x > h.x_min_eps then
      0
    else if i >= h.n_bins - 1 && x < h.x_max_eps then
      h.n_bins - 1
    else
      i in
  if i < 0 then begin
    h.n_underflow <- h.n_underflow + 1;
    h.underflow <- h.underflow +. f;
    h.underflow2 <- h.underflow2 +. f *. f
  end else if i >= h.n_bins then begin
    h.n_overflow <- h.n_overflow + 1;
    h.overflow <- h.overflow +. f;
    h.overflow2 <- h.overflow2 +. f *. f
  end else begin
    h.n.(i) <- h.n.(i) + 1;
    h.w.(i) <- h.w.(i) +. f;
    h.w2.(i) <- h.w2.(i) +. f *. f
  end

let normalize h =
  let sum_w = Array.fold_left (+.) (h.underflow +. h.overflow) h.w in
  let sum_w2 = sum_w *. sum_w in
  { n_bins = h.n_bins;
    n_bins_float = h.n_bins_float;
    x_min = h.x_min;
    x_max = h.x_max;
    x_min_eps = h.x_min_eps;
    x_max_eps = h.x_max_eps;
    n_underflow = h.n_underflow;
    underflow = h.underflow /. sum_w;
    underflow2 = h.underflow2 /. sum_w2;
    n_overflow = h.n_overflow;
    overflow = h.overflow /. sum_w;
    overflow2 = h.overflow2 /. sum_w2;
    n = Array.copy h.n;
    w = Array.map (fun w' -> w' /. sum_w) h.w;
    w2 = Array.map (fun w2' -> w2' /. sum_w2) h.w2 }

let to_channel oc h =
  for i = 0 to h.n_bins - 1 do
    let x_mid = h.x_min
        +. (h.x_max -. h.x_min) *. (float i +. 0.5) /. h.n_bins_float in
    if h.n.(i) > 1 then
      let n = float h.n.(i) in
      (* [let var1 = (h.w2.(i) /. n -. (h.w.(i) /. n) ** 2.0) /. (n -. 1.0)] *)
      let var2 = h.w.(i) ** 2.0 /. (n *. (n -. 1.0)) in
      let var = var2 in
      fprintf oc " %.17E %.17E %.17E\n" x_mid h.w.(i) (sqrt var)
    else if h.n.(i) = 1 then
      fprintf oc " %.17E %.17E %.17E\n" x_mid h.w.(i) h.w.(i)
    else
      fprintf oc " %.17E %.17E\n" x_mid h.w.(i)
  done

let as_bins_to_channel oc h =
  for i = 0 to h.n_bins - 1 do
    let x_min = h.x_min
        +. (h.x_max -. h.x_min) *. (float i) /. h.n_bins_float
    and x_max = h.x_min
        +. (h.x_max -. h.x_min) *. (float i +. 1.0) /. h.n_bins_float in
    fprintf oc " %.17e %.17e\n" x_min h.w.(i);
    fprintf oc " %.17e %.17e\n" x_max h.w.(i)
  done

(*i
let to_channel oc h =
  for i = 0 to h.n_bins - 1 do
    let x_min = h.x_min
        +. (h.x_max -. h.x_min) *. (float i) /. h.n_bins_float
    and x_max = h.x_min
        +. (h.x_max -. h.x_min) *. (float i +. 1.0) /. h.n_bins_float in
    fprintf oc " %.17E 0\n" x_min;
    fprintf oc " %.17E %.17E\n" x_min h.w.(i);
    fprintf oc " %.17E %.17E\n" x_max h.w.(i);
    fprintf oc " %.17E 0\n" x_max
  done
i*)

let to_file name h =
  let oc = open_out name in
  to_channel oc h;
  close_out oc

let as_bins_to_file name h =
  let oc = open_out name in
  as_bins_to_channel oc h;
  close_out oc

(* \subsection{Naive Linear Regression} *)

type regression_moments =
    { mutable n : int;
      mutable x : float;
      mutable y : float;
      mutable xx : float;
      mutable xy : float }

let init_regression_moments =
    { n = 0;
      x = 0.0;
      y = 0.0;
      xx = 0.0;
      xy = 0.0 }

let record_regression m x y =
  m.n <- m.n + 1;
  m.x <- m.x +. x;
  m.y <- m.y +. y;
  m.xx <- m.xx +. x *. x;
  m.xy <- m.xy +. x *. y
  
(* Minimize
   \begin{equation}
     f(a,b) = \sum_{i} w_{i} (ax_{i}+b-y_{i})^2 = \langle(ax+b-y)^2\rangle
   \end{equation}
   i.\,e.
   \begin{subequations}
   \begin{align}
     \frac{1}{2}\frac{\partial f}{\partial a}(a,b)
        &= \langle x(ax+b-y) \rangle
         = a\langle x^2 \rangle + b\langle x \rangle - \langle xy \rangle = 0 \\
     \frac{1}{2}\frac{\partial f}{\partial b}(a,b)
        &= \langle ax+b-y \rangle
         = a\langle x \rangle + b - \langle y \rangle = 0
   \end{align}
   \end{subequations}
   and
   \begin{subequations}
   \begin{align}
     a &= \frac{\langle xy \rangle - \langle x \rangle \langle y \rangle}%
               {\langle x^2 \rangle - \langle x \rangle^2} \\
     b &= \langle y \rangle - a\langle x \rangle
   \end{align}
   \end{subequations} *)

let linear_regression m =
  let n = float m.n in
  let x = m.x /. n
  and y = m.y /. n
  and xx = m.xx /. n
  and xy = m.xy /. n in
  let a = (xy -. x *. y) /. (xx -. x *. x) in
  let b = y -. a *. x in
  (a, b)
    
let regression h chi fx fy =
  let m = init_regression_moments in
  for i = 0 to h.n_bins - 1 do
    let x_mid = h.x_min
        +. (h.x_max -. h.x_min) *. (float i +. 0.5) /. h.n_bins_float in
    if chi x_mid then
      record_regression m (fx x_mid) (fy h.w.(i))
  done;
  linear_regression m

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)




