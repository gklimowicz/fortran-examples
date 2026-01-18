(* circe2/diffmap.ml --  *)
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

module type T =
  sig

    type t

    type domain
    val x_min : t -> domain
    val x_max : t -> domain

    type codomain
    val y_min : t -> codomain
    val y_max : t -> codomain

    val phi : t -> domain -> codomain
    val ihp : t -> codomain -> domain
    val jac : t -> domain -> float
    val caj : t -> codomain -> float

    val with_domain : t -> x_min:domain -> x_max:domain -> t

    val encode : t -> string

  end

module type Real = T with type domain = float and type codomain = float

(* \subsection{Testing Real Maps} *)

module type Test =
  sig
    module M : Real
    val domain : M.t -> unit
    val inverse : M.t -> unit
    val jacobian : M.t -> unit
    val all : M.t -> unit
  end

module Make_Test (M : Real) =
  struct

    module M = M

    let steps = 1000
    let epsilon = 1.0e-6

    let diff ?(tolerance = 1.0e-13) x1 x2 =
      let d = (x1 -. x2)  in
      if abs_float d < (abs_float x1 +. abs_float x2) *. tolerance then
        0.0
      else
        d

    let derive x_min x_max f x =
      let xp = min x_max (x +. epsilon)
      and xm = max x_min (x -. epsilon) in
      (f xp -. f xm) /. (xp -.xm)

    let domain m =
      let x_min = M.x_min m
      and x_max = M.x_max m
      and y_min = M.y_min m
      and y_max = M.y_max m in
      let x_min' = M.ihp m y_min
      and x_max' = M.ihp m y_max
      and y_min' = M.phi m x_min
      and y_max' = M.phi m x_max in
      printf "   f: [%g,%g] -> [%g,%g] ([%g,%g])\n"
        x_min x_max y_min' y_max' (diff y_min' y_min) (diff y_max' y_max);
      printf "f^-1: [%g,%g] -> [%g,%g] ([%g,%g])\n"
        y_min y_max x_min' x_max' (diff x_min' x_min) (diff x_max' x_max)

    let inverse m =
      let x_min = M.x_min m
      and x_max = M.x_max m
      and y_min = M.y_min m
      and y_max = M.y_max m in
      for i = 1 to steps do
        let x = x_min +. Random.float (x_max -. x_min)
        and y = y_min +. Random.float (y_max -. y_min) in
        let x' = M.ihp m y
        and y' = M.phi m x in
        let x'' = M.ihp m y'
        and y'' = M.phi m x' in
        let dx = diff x'' x
        and dy = diff y'' y in
        if dx <> 0.0 then
          printf "f^-1 o f   : %g -> %g -> %g (%g)\n" x y' x'' dx;
        if dy <> 0.0 then
          printf "   f o f^-1: %g -> %g -> %g (%g)\n" y x' y'' dy
      done

    let jacobian m =
      let x_min = M.x_min m
      and x_max = M.x_max m
      and y_min = M.y_min m
      and y_max = M.y_max m in
      for i = 1 to steps do
        let x = x_min +. Random.float (x_max -. x_min)
        and y = y_min +. Random.float (y_max -. y_min) in
        let jac_x' = derive x_min x_max (M.phi m) x
        and jac_x = M.jac m x
        and inv_jac_y' = derive y_min y_max (M.ihp m) y
        and inv_jac_y = M.caj m y in
        let dj = diff ~tolerance:1.0e-9 jac_x' jac_x
        and dij = diff ~tolerance:1.0e-9 inv_jac_y' inv_jac_y in
        if dj <> 0.0 then
          printf "dy/dx: %g -> %g (%g)\n" x jac_x' dj;
        if dij <> 0.0 then
          printf "dx/dy: %g -> %g (%g)\n" y inv_jac_y' dij
      done

    let all m =
      printf "phi(domain) = codomain and phi(codomain) = domain";
      domain m;
      printf "ihp o phi = id (domain) and phi o ihp = id(codomain)";
      inverse m;
      printf "jacobian";
      jacobian m

  end

(* \subsection{Specific Real Maps} *)

module Id =
  struct

    type domain = float
    type codomain = float

    type t =
        { x_min : domain;
          x_max : domain;
          y_min : codomain;
          y_max : codomain;
          phi : float -> float;
          ihp : float -> float;
          jac : float -> float;
          caj : float -> float }

    let encode m = "0 1 0 0 1 1"

    let closure ~x_min ~x_max ~y_min ~y_max =
      let phi x = x
      and ihp y = y
      and jac x = 1.0
      and caj y = 1.0 in
      { x_min = x_min;
        x_max = x_max;
        y_min = y_min;
        y_max = y_max;
        phi = phi;
        ihp = ihp;
        jac = jac;
        caj = caj }

    let idmap ~x_min ~x_max ~y_min ~y_max =
      if x_min <> y_min && x_max <> y_max then
        invalid_arg "Diffmap.Id.idmap"
      else
        closure ~x_min ~x_max ~y_min ~y_max

    let with_domain m ~x_min ~x_max =
      idmap ~x_min ~x_max ~y_min:m.y_min ~y_max:m.y_max

    let create ?x_min ?x_max y_min y_max =
      idmap
        ~x_min:(match x_min with Some x -> x | None -> y_min)
        ~x_max:(match x_max with Some x -> x | None -> y_max)
        ~y_min ~y_max

    let x_min m = m.x_min
    let x_max m = m.x_max
    let y_min m = m.y_min
    let y_max m = m.y_max

    let phi m = m.phi
    let ihp m = m.ihp
    let jac m = m.jac
    let caj m  = m.caj

  end

module Linear =
  struct

    type domain = float
    type codomain = float

    type t =
        { x_min : domain;
          x_max : domain;
          y_min : codomain;
          y_max : codomain;
          a : float;
          b : float;
          phi : domain -> codomain;
          ihp : codomain -> domain;
          jac : domain -> float;
          caj : codomain -> float }

    let encode m = failwith "Diffmap.Linear: not used in Circe2"

    let closure ~x_min ~x_max ~y_min ~y_max ~a ~b =

(* \begin{equation}
      x \mapsto \lambda_{a,b}(x) = ax + b
   \end{equation} *)
      let phi x = a *. x +. b

(* \begin{equation}
      y \mapsto (\lambda_{a,b})^{-1}(y) = \frac{y-b}{a}
   \end{equation} *)
      and ihp y = (y -. b) /. a

      and jac x = a
      and caj y = 1.0 /. a in

      { x_min = x_min;
        x_max = x_max;
        y_min = y_min;
        y_max = y_max;
        a = a;
        b = b;
        phi = phi;
        ihp = ihp;
        jac = jac;
        caj = caj }

    let linearmap ~x_min ~x_max ~y_min ~y_max =
      let delta_x = x_max -. x_min
      and delta_y = y_max -. y_min in
      let a = delta_y /. delta_x
      and b = (y_min *. x_max -. y_max *. x_min) /. delta_x in
      closure ~x_min ~x_max ~y_min ~y_max ~a ~b

    let with_domain m ~x_min ~x_max =
      linearmap ~x_min ~x_max ~y_min:m.y_min ~y_max:m.y_max

    let create ?x_min ?x_max y_min y_max =
      linearmap
        ~x_min:(match x_min with Some x -> x | None -> y_min)
        ~x_max:(match x_max with Some x -> x | None -> y_max)
        ~y_min ~y_max

    let x_min m = m.x_min
    let x_max m = m.x_max
    let y_min m = m.y_min
    let y_max m = m.y_max

    let phi m = m.phi
    let ihp m = m.ihp
    let jac m = m.jac
    let caj m  = m.caj

  end

module Power =
  struct

    type domain = float
    type codomain = float

    type t =
        { x_min : domain;
          x_max : domain;
          y_min : codomain;
          y_max : codomain;
          alpha : float;
          xi : float;
          eta : float;
          a : float;
          b : float;
          phi : domain -> codomain;
          ihp : codomain -> domain;
          jac : domain -> float;
          caj : codomain -> float }

    let encode m =
      sprintf "1 %s %s %s %s %s"
        (Float.Double.to_string m.alpha)
        (Float.Double.to_string m.xi)
        (Float.Double.to_string m.eta)
        (Float.Double.to_string m.a)
        (Float.Double.to_string m.b)

    let closure ~x_min ~x_max ~y_min ~y_max ~alpha ~xi ~eta ~a ~b =

(* \begin{equation}
      x \mapsto \psi_{a,b}^{\alpha,\xi,\eta}(x)
          = \frac{1}{b}(a(x-\xi))^{\alpha} + \eta
   \end{equation} *)
      let phi x =
        (a *. (x -. xi)) ** alpha /. b +. eta

(* \begin{equation}
      y \mapsto (\psi_{a,b}^{\alpha,\xi,\eta})^{-1}(y)
          = \frac{1}{a} (b(y-\eta))^{1/\alpha} + \xi
   \end{equation} *)
      and ihp y =
        (b *. (y -. eta)) ** (1.0 /. alpha) /. a +. xi

(* \begin{equation}
     \frac{\mathrm{d}y}{\mathrm{d}x} (x)
        = \frac{a\alpha}{b} (a(x-\xi))^{\alpha-1}
   \end{equation} *)
      and jac x =
        a *. alpha *. (a *. (x -. xi)) ** (alpha -. 1.0) /. b

(* \begin{equation}
     \frac{\mathrm{d}x}{\mathrm{d}y} (y)
        = \frac{b}{a\alpha} (b(y-\eta))^{1/\alpha-1}
   \end{equation} *)
      and caj y =
        b *. (b *. (y -. eta)) ** (1.0 /. alpha -. 1.0) /. (a *. alpha) in

      { x_min = x_min;
        x_max = x_max;
        y_min = y_min;
        y_max = y_max;
        alpha = alpha;
        xi = xi;
        eta = eta;
        a = a;
        b = b;
        phi = phi;
        ihp = ihp;
        jac = jac;
        caj = caj }

(* \begin{subequations}
   \begin{align}
     a_{i}   &= \frac{   (b_{i}(y_{i}-\eta_{i}))^{1/\alpha_{i}}
                       - (b_{i}(y_{i-1}-\eta_{i}))^{1/\alpha_{i}}}%
                     {x_{i} - x_{i-1}} \\
     \xi_{i} &= \frac{   x_{i-1}|y_{i}  -\eta_{i}|^{1/\alpha_{i}}
                       - x_{i}  |y_{i-1}-\eta_{i}|^{1/\alpha_{i}}}%
                     {          |y_{i}  -\eta_{i}|^{1/\alpha_{i}}
                       -        |y_{i-1}-\eta_{i}|^{1/\alpha_{i}}}
   \end{align}
   \end{subequations}
   The degeneracy~(\ref{eq:ab-semigroup}) can finally be resolved by
   demanding~$|b|=1$ in~(\ref{eq:ai}). *)

    let powermap ~x_min ~x_max ~y_min ~y_max ~alpha ~eta =
      let b =
        if eta <= y_min then
          1.
        else if eta >= y_max then
          -1.
        else
          invalid_arg "singular" in
      let pow y = (b *. (y -. eta)) ** (1. /. alpha) in
      let delta_pow = pow y_max -. pow y_min
      and delta_x = x_max -. x_min in
      let a = delta_pow /. delta_x
      and xi = (x_min *. pow y_max -. x_max *. pow y_min) /. delta_pow in
      closure ~x_min ~x_max ~y_min ~y_max ~alpha ~xi ~eta ~a ~b

    let with_domain m ~x_min ~x_max =
      powermap ~x_min ~x_max ~y_min:m.y_min ~y_max:m.y_max
        ~alpha:m.alpha ~eta:m.eta

    let create ~alpha ~eta ?x_min ?x_max y_min y_max =
      powermap
        ~x_min:(match x_min with Some x -> x | None -> y_min)
        ~x_max:(match x_max with Some x -> x | None -> y_max)
        ~y_min ~y_max ~alpha ~eta

    let x_min m = m.x_min
    let x_max m = m.x_max
    let y_min m = m.y_min
    let y_max m = m.y_max

    let phi m = m.phi
    let ihp m = m.ihp
    let jac m = m.jac
    let caj m  = m.caj

  end

module Resonance =
  struct

    type domain = float
    type codomain = float

    type t =
        { x_min : domain;
          x_max : domain;
          y_min : codomain;
          y_max : codomain;
          xi : float;
          eta : float;
          a : float;
          b : float;
          phi : domain -> codomain;
          ihp : codomain -> domain;
          jac : domain -> float;
          caj : codomain -> float }

    let encode m =
      sprintf "2 0 %s %s %s %s"
        (Float.Double.to_string m.xi)
        (Float.Double.to_string m.eta)
        (Float.Double.to_string m.a)
        (Float.Double.to_string m.b)

    let closure ~x_min ~x_max ~y_min ~y_max ~xi ~eta ~a ~b =

(* \begin{equation}
     x \mapsto \rho_{a,b}^{\xi,\eta}(x)
         = a \tan\left(\frac{a}{b^2}(x-\xi)\right) + \eta
   \end{equation} *)
      let phi x = a *. tan (a *. (x -. xi) /. (b *. b)) +. eta

(* \begin{equation}
     y \mapsto (\rho_{a,b}^{\xi,\eta})^{-1}(y)
         = \frac{b^2}{a} \textrm{atan}\left(\frac{y-\eta}{a}\right) + \xi
   \end{equation} *)
      and ihp y = b *. b *. (atan2 (y -. eta) a) /. a +. xi

(* \begin{equation}
     \frac{\mathrm{d}y}{\mathrm{d}x}(x(y))
       = \frac{1}{{\displaystyle\frac{\mathrm{d}x}{\mathrm{d}y}(y)}}
       = \left(\frac{b^2}{(y-\eta)^2+a^2}\right)^{-1}
   \end{equation} *)
      and caj y = b *. b /. ((y -. eta) ** 2.0 +. a *. a) in

      let jac x = 1.0 /. caj (phi x) in

      { x_min = x_min;
        x_max = x_max;
        y_min = y_min;
        y_max = y_max;
        xi = xi;
        eta = eta;
        a = a;
        b = b;
        phi = phi;
        ihp = ihp;
        jac = jac;
        caj = caj }

(* \begin{subequations}
   \begin{align}
     b_{i}
       &= \sqrt{a_{i}
            \frac{x_{i}-x_{i-1}}%
                 {   \textrm{atan}\left(\frac{y_{i}  -\eta_{i}}{a_{i}}\right)
                   - \textrm{atan}\left(\frac{y_{i-1}-\eta_{i}}{a_{i}}\right)}} \\
     \xi_{i}
       &= \frac{   x_{i-1}\textrm{atan}\left(\frac{y_{i}  -\eta_{i}}{a_{i}}\right)
                 - x_{i}  \textrm{atan}\left(\frac{y_{i-1}-\eta_{i}}{a_{i}}\right)}%
               {x_{i}-x_{i-1}}
   \end{align}
   \end{subequations} *)
   
   let resonancemap ~x_min ~x_max ~y_min ~y_max ~eta ~a =
     let arc y = atan2 (y -. eta) a in
     let delta_arc = arc y_max -. arc y_min
     and delta_x = x_max -. x_min in
     let b = sqrt (a *. delta_x /. delta_arc)
     and xi = (x_min *. arc y_max -. x_max *. arc y_min) /. delta_arc in
     closure ~x_min ~x_max ~y_min ~y_max ~xi ~eta ~a ~b
   
    let with_domain m ~x_min ~x_max =
      resonancemap ~x_min ~x_max ~y_min:m.y_min ~y_max:m.y_max
        ~eta:m.eta ~a:m.a

    let create ~eta ~a ?x_min ?x_max y_min y_max =
      resonancemap
        ~x_min:(match x_min with Some x -> x | None -> y_min)
        ~x_max:(match x_max with Some x -> x | None -> y_max)
        ~y_min ~y_max ~eta ~a

    let x_min m = m.x_min
    let x_max m = m.x_max
    let y_min m = m.y_min
    let y_max m = m.y_max

    let phi m = m.phi
    let ihp m = m.ihp
    let jac m = m.jac
    let caj m  = m.caj

  end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
