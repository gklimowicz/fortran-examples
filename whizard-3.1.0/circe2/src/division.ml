(* circe2/division.ml --  *)
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

let epsilon_100 = 100.0 *. Float.Double.epsilon

let equidistant n x_min x_max =
  if n <= 0 then
    invalid_arg "Division.equidistant: n <= 0"
  else
    let delta = (x_max -. x_min) /. (float n) in
    Array.init (n + 1) (fun i -> x_min +. delta *. float i)

exception Above_max of float * (float * float) * int
exception Below_min of float * (float * float) * int
exception Out_of_range of float * (float * float)
exception Rebinning_failure of string

let find_raw d x =
  let n_max = Array.length d - 1 in
  let eps = epsilon_100 *. (d.(n_max) -. d.(0)) in
  let rec find' a b =
    if b <= a + 1 then
      a
    else
      let m = (a + b) / 2 in
      if x < d.(m) then
	find' a m
      else
	find' m b in
  if x < d.(0) -. eps then
    raise (Below_min (x, (d.(0), d.(n_max)), 0))
  else if x > d.(n_max) +. eps then
    raise (Above_max (x, (d.(0), d.(n_max)), n_max - 1))
  else if x <= d.(0) then
    0
  else if x >= d.(n_max) then
    n_max - 1
  else
    find' 0 n_max

module type T =
  sig

    type t

    val copy : t -> t

    val find : t -> float -> int
    val record : t -> float -> float -> unit

    val rebin : ?power:float -> ?fixed_min:bool -> ?fixed_max:bool -> t -> t

    val caj : t -> float -> float

    val n_bins : t -> int
    val bins : t -> float array
    val to_channel : out_channel -> t -> unit

  end

(* \subsubsection{Primary Divisions} *)

module type Mono =
  sig
    include T
    val create : ?bias:(float -> float) -> int -> float -> float -> t
  end

module Mono (* [: T] *) =
  struct

    type t =
        { x : float array;
          mutable x_min : float;
          mutable x_max : float;
          n : int array;
          w : float array;
          w2 : float array;
          bias : float -> float } 

    let copy d =
      { x = Array.copy d.x;
        x_min = d.x_min;
        x_max = d.x_max;
        n = Array.copy d.n;
        w = Array.copy d.w;
        w2 = Array.copy d.w2;
        bias = d.bias }

    let create ?(bias = fun x -> 1.0) n x_min x_max =
      { x = equidistant n x_min x_max;
        x_min = x_max;
        x_max = x_min;
        n = Array.make n 0;
        w = Array.make n 0.0;
        w2 = Array.make n 0.0;
        bias = bias }

    let bins d = d.x
    let n_bins d = Array.length d.x - 1

    let find d = find_raw d.x

    let normal_float x =
      match classify_float x with
      | FP_normal | FP_subnormal | FP_zero -> true
      | FP_infinite | FP_nan -> false

    let report_denormal x f b what =
      eprintf
        "circe2: Division.record: ignoring %s (x=%g, f=%g, b=%g)\n"
        what x f b;
      flush stderr
      
(*i let report_denormal x f b what = () i*)

    let caj d x = 1.0

    let record d x f =
      if x < d.x_min then
        d.x_min <- x;
      if x > d.x_max then
        d.x_max <- x;
      let i = find d x in
      d.n.(i) <- succ d.n.(i);
      let b = d.bias x in
      let w = f *. b in
      match classify_float w with
      | FP_normal | FP_subnormal | FP_zero ->
          d.w.(i) <- d.w.(i) +. w;
          let w2 = f *. w in
          begin match classify_float w2 with
          | FP_normal | FP_subnormal | FP_zero ->
              d.w2.(i) <- d.w2.(i) +. w2
          | FP_infinite -> report_denormal x f b "w2 = [inf]"
          | FP_nan -> report_denormal x f b "w2 = [nan]"
          end
      | FP_infinite -> report_denormal x f b "w2 = [inf]"
      | FP_nan -> report_denormal x f b "w2 = [nan]"

    (* \begin{equation}
         \begin{aligned}
           d_1     &\to \frac{1}{2}(d_1+d_2) \\
           d_2     &\to \frac{1}{3}(d_1+d_2+d_3) \\
                   &\ldots\\
           d_{n-1} &\to \frac{1}{3}(d_{n-2}+d_{n-1}+d_n) \\
           d_n     &\to \frac{1}{2}(d_{n-1}+d_n)
         \end{aligned}
       \end{equation} *)

    let smooth3 f =
      match Array.length f with
      | 0 -> f
      | 1 -> Array.copy f
      | 2 -> Array.make 2 ((f.(0) +. f.(1)) /. 2.0)
      | n ->
          let f' = Array.make n 0.0 in
          f'.(0) <- (f.(0) +. f.(1)) /. 2.0;
          for i = 1 to n - 2 do
            f'.(i) <- (f.(i-1) +. f.(i) +. f.(i+1)) /. 3.0
          done;
          f'.(n-1) <- (f.(n-2) +. f.(n-1)) /. 2.0;
          f'

    (* \begin{equation}
         m_i = \left(
           \frac{\frac{\bar f_i \Delta x_i}{\sum_j\bar f_j \Delta x_j}-1}
                {\ln\left(\frac{\bar f_i \Delta x_i}{\sum_j\bar f_j \Delta x_j}\right)}
             \right)^\alpha
       \end{equation} *)

    let rebinning_weights' power fs =
      let sum_f = Array.fold_left (+.) 0.0 fs in
      if sum_f <= 0.0 then
        Array.make (Array.length fs) 1.0
      else
        Array.map (fun f ->
          let f' = f /. sum_f in
          if f' < 1.0e-12 then
            0.
          else
            ((f' -. 1.0) /. (log f')) ** power) fs

(*i\begin{figure}
     \begin{center}
       \empuse{weights}
       \begin{empgraph}(70,30)
         randomseed := 720.251;
         pickup pencircle scaled 0.7pt;
         path m[], g[];
         numeric pi; pi = 180;
         numeric dx; dx = 0.05;
         numeric dg; dg = -0.04;
         vardef adap_fct (expr x) = (x + sind(4*x*pi)/16) enddef;
         autogrid (,);
         frame.bot;
         setrange (0, 0, 1, 1.2);
         for x = 0 step dx until 1+dx/2:
           numeric r;
           r = 1 + normaldeviate/10;
           augment.m[x] (adap_fct (x), r);
           augment.m[x] (adap_fct (x+dx), r);
           augment.m[x] (adap_fct (x+dx), 0);
           augment.m[x] (adap_fct (x), 0);
           augment.g[x] (adap_fct (x), 0);
           augment.g[x] (adap_fct (x), dg);
         endfor
         for x = 0 step dx until 1-dx/2:
           gfill m[x] -- cycle withcolor 0.7white;
           gdraw m[x] -- cycle;
         endfor
         for x = 0 step dx until 1+dx/2:
           gdraw g[x];
         endfor
         glabel.bot (btex $x_0$     etex, (adap_fct (0*dx), dg));
         glabel.bot (btex $x_1$     etex, (adap_fct (1*dx), dg));
         glabel.bot (btex $x_2$     etex, (adap_fct (2*dx), dg));
         glabel.bot (btex $x_{n-1}$ etex, (adap_fct (1-dx), dg));
         glabel.bot (btex $x_n$     etex, (adap_fct (1), dg));
         glabel.lft (btex $\displaystyle
                      \bar f_i\approx\frac{m_i}{\Delta x_i}$ etex, OUT);
       \end{empgraph}
     \end{center}
     \caption{\label{fig:rebin}%
       Typical weights used in the rebinning algorithm.}
   \end{figure}
i*)

    (* The nested loops can be turned into recursions, of course.  But arrays aren't
       purely functional anyway \ldots *)

    let rebin' m x =
      let n = Array.length x - 1 in
      let x' = Array.make (n + 1) 0.0 in
      let sum_m = Array.fold_left (+.) 0.0 m in
      if sum_m <= 0.0 then
        Array.copy x
      else begin
        let step = sum_m /. (float n) in
        let k = ref 0
        and delta = ref 0.0 in
        x'.(0) <- x.(0);
        for i = 1 to n - 1 do

          (* We increment~$k$ until another $\Delta$ (a.\,k.\,a.~[step]) of the
             integral has been accumulated.  % (cf.~figure~\ref{fig:rebin}). 
	  *)
          while !delta < step do
            incr k;
            delta := !delta +. m.(!k-1)
          done;
    
          (* Correct the mismatch. *)
          delta := !delta -. step;

          (* Linearly interpolate the next bin boundary.  *)
          x'.(i) <- x.(!k) -. (x.(!k) -. x.(!k-1)) *. !delta /. m.(!k-1);

          if x'.(i) < x'.(i-1) then
            raise (Rebinning_failure
                     (sprintf "x(%d)=%g < x(%d)=%g" i x'.(i) (i-1) x'.(i-1)))

        done;
        x'.(n) <- x.(n);
        x'
      end

    (* \begin{dubious}
         Check that [x_min] and [x_max] are implemented correctly!!!!
       \end{dubious} *)

    (* \begin{dubious}
         One known problem is that the second outermost bins hinder the
         outermost bins from moving.
       \end{dubious} *)

    let rebin ?(power = 1.5) ?(fixed_min = false) ?(fixed_max = false) d =
      let n = Array.length d.w in
      let x = rebin' (rebinning_weights' power (smooth3 d.w2)) d.x in
      if not fixed_min then
        x.(0) <- (x.(0) +. min d.x_min x.(1)) /. 2.;
      if not fixed_max then
        x.(n) <- (x.(n) +. max d.x_max x.(n-1)) /. 2.;
      { x = x;
        x_min = d.x_min;
        x_max = d.x_max;
        n = Array.make n 0;
        w = Array.make n 0.0;
        w2 = Array.make n 0.0;
        bias = d.bias }

    let to_channel oc d =
      Array.iter (fun x ->
        fprintf oc " %s 0 1 0 0 1 1\n" (Float.Double.to_string x)) d.x

  end

(* \subsubsection{Polydivisions} *)

module type Poly =
  sig

    module M : Diffmaps.Real

    include T

    val create : ?bias:(float -> float) ->
      (int * M.t) list -> int -> float -> float -> t

  end

module Make_Poly (M : Diffmaps.Real) (* [: Poly] *) =
  struct

    module M = M

    type t =
        { x : float array;
          d : Mono.t array;
          n_bins : int;
          ofs : int array;
          maps : M.t array;
          n : int array;
          w : float array;
          w2 : float array }

    let copy pd =
      { x = Array.copy pd.x;
        d = Array.map Mono.copy pd.d;
        n_bins = pd.n_bins;
        ofs = Array.copy pd.ofs;
        maps = Array.copy pd.maps;
        n = Array.copy pd.n;
        w = Array.copy pd.w;
        w2 = Array.copy pd.w2 }

    let n_bins pd = pd.n_bins

    let find pd y =
      let i = find_raw pd.x y in
      let x = M.ihp pd.maps.(i) y in
      pd.ofs.(i) + Mono.find pd.d.(i) x

    let bins pd =
      let a = Array.make (pd.n_bins + 1) 0.0 in
      let bins0 = Mono.bins pd.d.(0) in
      let len = Array.length bins0 in
      Array.blit bins0 0 a 0 len;
      let ofs = ref len in
      for i = 1 to Array.length pd.d - 1 do
        let len = Mono.n_bins pd.d.(i) in
        Array.blit (Mono.bins pd.d.(i)) 1 a !ofs len;
        ofs := !ofs + len
      done;
      a
        
    type interval =
        { nbin : int;
          x_min : float;
          x_max : float;
          map : M.t }

    let interval nbin map =
      { nbin = nbin;
        x_min = M.x_min map;
        x_max = M.x_max map;
        map = map }

    let id_map n y_min y_max =
      interval n (M.id ~x_min:y_min ~x_max:y_max y_min y_max)

    let sort_intervals intervals =
      List.sort (fun i1 i2 -> compare i1.x_min i2.x_min) intervals

    (* Fill the gaps between adjacent intervals, using
       [val default : int -> float -> float -> interval] to
       construct intermediate intervals. *)

    let fill_gaps default n x_min x_max intervals =
      let rec fill_gaps' prev_x_max acc = function
        | i :: rest ->
            if i.x_min = prev_x_max then
              fill_gaps' i.x_max (i :: acc) rest
            else if i.x_min > prev_x_max then
              fill_gaps' i.x_max
                (i :: (default n prev_x_max i.x_min) :: acc) rest
            else
              invalid_arg "Polydivision.fill_gaps: overlapping"
        | [] ->
            if x_max = prev_x_max then
              List.rev acc
            else if x_max > prev_x_max then
              List.rev (default n prev_x_max x_max :: acc)
            else
              invalid_arg "Polydivision.fill_gaps: sticking out" in
      match intervals with
      | i :: rest ->
          if i.x_min = x_min then
            fill_gaps' i.x_max [i] rest
          else if i.x_min > x_min then
            fill_gaps' i.x_max (i :: [default n x_min i.x_min]) rest
          else
            invalid_arg "Polydivision.fill_gaps: sticking out"
      | [] -> [default n x_min x_max]

    let create ?bias intervals n x_min x_max =
      let intervals = List.map (fun (n, m) -> interval n m) intervals in
      match fill_gaps id_map n x_min x_max (sort_intervals intervals) with
      | [] -> failwith "Division.Poly.create: impossible"
      | interval :: _ as intervals ->
          let ndiv = List.length intervals in
          let x = Array.of_list (interval.x_min ::
                                 List.map (fun i -> i.x_max) intervals) in
          let d = Array.of_list
              (List.map (fun i ->
                Mono.create ?bias i.nbin i.x_min i.x_max) intervals) in
          let ofs = Array.make ndiv 0 in
          for i = 1 to ndiv - 1 do
            ofs.(i) <- ofs.(i-1) + Mono.n_bins d.(i-1)
          done;
          let n_bins = ofs.(ndiv-1) + Mono.n_bins d.(ndiv-1) in
          { x = x;
            d = d;
            n_bins = n_bins;
            ofs = ofs;
            maps = Array.of_list (List.map (fun i -> i.map) intervals);
            n = Array.make ndiv 0;
            w = Array.make ndiv 0.0;
            w2 = Array.make ndiv 0.0 }

    (* We can safely assume that [find_raw pd.x y = find_raw pd.x x].
       \begin{equation}
         w = \frac{f}{{\displaystyle\frac{\mathrm{d}x}{\mathrm{d}y}}}
           = f\cdot\frac{\mathrm{d}y}{\mathrm{d}x}
       \end{equation}
       Here, the jacobian makes no difference for the final result, but
       it steers VEGAS/VAMP into the right direction.  *)

    let caj pd y =
      let i = find_raw pd.x y in
      let m = pd.maps.(i)
      and d = pd.d.(i) in
      let x = M.ihp m y in
      M.caj m y *. Mono.caj d x

    let record pd y f =
      let i = find_raw pd.x y in
      let m = pd.maps.(i) in
      let x = M.ihp m y in
      let w = M.jac m x *. f in
      Mono.record pd.d.(i) x w;
      pd.n.(i) <- succ pd.n.(i);
      pd.w.(i) <- pd.w.(i) +. w;
      pd.w2.(i) <- pd.w2.(i) +. w *. w

    (* Rebin the divisions, enforcing fixed boundaries for the inner
       intervals. *)
    let rebin ?(power = 1.5) ?(fixed_min = false) ?(fixed_max = false) pd =
      let ndiv = Array.length pd.d in
      let rebin_mono i d =
        if ndiv <= 1 then
          Mono.rebin ~power ~fixed_min ~fixed_max d
        else if i = 0 then
          Mono.rebin ~power ~fixed_min ~fixed_max:true d
        else if i = ndiv - 1 then
          Mono.rebin ~power ~fixed_min:true ~fixed_max d
        else
          Mono.rebin ~power ~fixed_min:true ~fixed_max:true d in
      { x = Array.copy pd.x;
        d = Array.init ndiv (fun i -> rebin_mono i pd.d.(i));
        n_bins = pd.n_bins;
        ofs = pd.ofs;
        maps = Array.copy pd.maps;
        n = Array.make ndiv 0;
        w = Array.make ndiv 0.0;
        w2 = Array.make ndiv 0.0 }

    let to_channel oc pd =
      for i = 0 to Array.length pd.d - 1 do
        let map = M.encode pd.maps.(i)
        and bins = Mono.bins pd.d.(i)
        and j0 = if i = 0 then 0 else 1 in
        for j = j0 to Array.length bins - 1 do
          fprintf oc " %s %s\n" (Float.Double.to_string bins.(j)) map;
        done
      done

  end

module Poly = Make_Poly (Diffmaps.Default)

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
