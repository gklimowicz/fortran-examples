(* circe2/diffmap.mli --  *)
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

module type T =
  sig

    type t

    (* An invertible differentiable map is characterized by its
       domain~$\lbrack x_{\min},x_{\max}\rbrack$ *)
    type domain
    val x_min : t -> domain
    val x_max : t -> domain

    (* and codomain~$\lbrack y_{\min},y_{\max}\rbrack$ *)
    type codomain
    val y_min : t -> codomain
    val y_max : t -> codomain

    (* the map proper
       \begin{equation}
         \begin{aligned}
           \phi : \lbrack x_{\min},x_{\max}\rbrack 
               &\to \lbrack y_{\min},y_{\max}\rbrack  \\
             x &\mapsto y = \phi(x)
         \end{aligned}
       \end{equation} *)
    val phi : t -> domain -> codomain

    (* the inverse map
       \begin{equation}
         \begin{aligned}
           \phi^{-1} : \lbrack y_{\min},y_{\max}\rbrack 
               &\to \lbrack x_{\min},x_{\max}\rbrack  \\
             y &\mapsto x = \phi^{-1}(y)
         \end{aligned}
       \end{equation} *)
    val ihp : t -> codomain -> domain

    (* the jacobian of the map
       \begin{equation}
         \begin{aligned}
           J : \lbrack x_{\min},x_{\max}\rbrack 
               &\to \mathbf{R} \\
             x &\mapsto J(x) = \frac{\mathrm{d}\phi}{\mathrm{d}x}(x)
         \end{aligned}
       \end{equation} *)
    val jac : t -> domain -> float

    (* and finally the jacobian of the inverse map
       \begin{equation}
         \begin{aligned}
           J^{*} : \lbrack y_{\min},y_{\max}\rbrack 
               &\to \mathbf{R} \\
             y &\mapsto J^{*}(y)
                  = \frac{\mathrm{d}\phi^{-1}}{\mathrm{d}y}(y)
                  = \left(\frac{\mathrm{d}\phi}{\mathrm{d}x}(\phi^{-1}(y))\right)^{-1}
         \end{aligned}
       \end{equation} *)
    val caj : t -> codomain -> float

    (* [with_domain map x_min x_max] takes the map [map] and
       returns the `same' map with the new
       domain~$\lbrack x_{\min},x_{\max}\rbrack$ *)
    val with_domain : t -> x_min:domain -> x_max:domain -> t

    (* There is also a convention for encoding the map so that
       it can be read by \KirkeTwo/: *)
    val encode : t -> string

  end

(* For the application in \KirkeTwo/, it suffices to consider real maps.
   Introducing [domain] and [codomain] does not make any difference for
   the typechecker as long as we only use [Diffmap.Real], but it provides
   documentation and keeps the door for extensions open. *)
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

module Make_Test (M : Real) : Test with module M = M

(* \subsection{Specific Real Maps} *)

module Id :
    sig
      include Real
 
      (* [create x_min x_max y_min y_max] creates an identity map
         $\lbrack x_{\min},x_{\max}\rbrack\to
          \lbrack y_{\min},y_{\max}\rbrack$.
         \begin{equation}
            \begin{aligned}
              \iota : \lbrack x_{\min},x_{\max}\rbrack
                 &\to \lbrack x_{\min},x_{\max}\rbrack  \\
               x &\mapsto \iota(x) = x
            \end{aligned}
         \end{equation}
         Default values for [x_min] and~[x_max] are [y_min] and~[y_max],
         respectively.   Indeed, they are the only
         possible values and other values raise an exception. *)
      val create :
          ?x_min:domain -> ?x_max:domain -> codomain -> codomain -> t
    end

module Linear :
    sig
      include Real

      (* [create x_min x_max y_min y_max] creates a linear map
         $\lbrack x_{\min},x_{\max}\rbrack\to\lbrack y_{\min},y_{\max}\rbrack$.
         The parameters~$a$ and~$b$ are determined from domain and codomain.
         \begin{equation}
            \begin{aligned}
              \lambda_{a,b} :
                \lbrack x_{\min},x_{\max}\rbrack  &\to \lbrack y_{\min},y_{\max}\rbrack  \\
               x &\mapsto \lambda_{a,b}(x) = ax + b
            \end{aligned}
         \end{equation}
         Default values for [x_min] and~[x_max] are [y_min] and~[y_max],
         respectively. *)
      val create :
          ?x_min:domain -> ?x_max:domain -> codomain -> codomain -> t
    end

module Power :
    sig
      include Real

      (* [create alpha eta x_min x_max y_min y_max] creates a power map
         $\lbrack x_{\min},x_{\max}\rbrack\to\lbrack y_{\min},y_{\max}\rbrack$.
         The parameters~$\xi$, $a$ and~$b$ are determined from~$\alpha$,
         $\eta$, domain and codomain. 
         \begin{equation}
            \begin{aligned}
              \psi_{a,b}^{\alpha,\xi,\eta} : \lbrack x_{\min},x_{\max}\rbrack
                 &\to \lbrack y_{\min},y_{\max}\rbrack  \\
               x &\mapsto \psi_{a,b}^{\alpha,\xi,\eta}(x)
                      = \frac{1}{b}(a(x-\xi))^{\alpha} + \eta
            \end{aligned}
         \end{equation}
         Default values for [x_min] and~[x_max] are [y_min] and~[y_max],
         respectively. *) 
      val create : alpha:float -> eta:float ->
        ?x_min:domain -> ?x_max:domain -> codomain -> codomain -> t
    end

module Resonance :
    sig
      include Real

      (* [create eta a x_min x_max y_min y_max] creates a resonance map
         $\lbrack x_{\min},x_{\max}\rbrack\to\lbrack y_{\min},y_{\max}\rbrack$.
         \begin{equation}
            \begin{aligned}
              \rho_{a,b}^{\xi,\eta}: \lbrack x_{\min},x_{\max}\rbrack 
                 &\to \lbrack y_{\min},y_{\max}\rbrack  \\
               x &\mapsto \rho_{a,b}^{\xi,\eta}(x)
                  = a \tan\left(\frac{a}{b^2}(x-\xi)\right) + \eta
            \end{aligned}
         \end{equation}
         The parameters~$\xi$ and~$b$ are determined from~$\eta$, $a$,
         domain and codomain. Default values for [x_min] and~[x_max] are
         [y_min] and~[y_max], respectively. *) 
      val create : eta:float -> a:float ->
        ?x_min:domain -> ?x_max:domain -> codomain -> codomain -> t
    end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)


