(* bundle.mli --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>

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

(* \begin{figure}
     \begin{center}
       \begin{emp}(80,80)
         ahlength := 3mm;
         ahangle := 20;
         pickup pencircle scaled 1.5pt;
         pair nw, ne, sw, se;
         nw = (.4w,.9h);
         ne = (.9w,.9h);
         sw = (.1w,.1h);
         se = (.6w,.1h);
         for i = 0 step 0.2 until 1:
           draw (i*sw+(1-i)*se){up}..{up}(i*nw+(1-i)*ne);
         endfor
         path base, fiber;
         base = (0,.5h){right}..{right}(w,.4h);
         fiber = (.6*sw+(1-.6)*se){up}..{up}(.6*nw+(1-.6)*ne);
         pickup pencircle scaled 3pt;
         draw base;
         pickup pencircle scaled 2pt;
         draw fiber;
         pickup pencircle scaled 1.5pt;
         drawarrow (.9w,.3h){up} .. {up}point .8 of base;
         label.bot (btex $B=\pi(E)$ etex, (.9w,.3h));
         drawarrow (.7w,.2h){up} .. {-1,1}(base intersectionpoint fiber);
         label.bot (btex $x\in B$ etex, (.7w,.2h));
         drawarrow (.2w,.8h){right} .. point .8 of fiber;
         label.lft (btex $\pi^{-1}(x)$ etex, (.2w,.8h));
         label.lft (btex $E = \pi^{-1}(b)$ etex, (.2w,.6h));
         setbounds currentpicture to (0,0)--(w,0)--(w,h)--(0,h)--cycle;
       \end{emp}
     \end{center}
     \caption{\label{fig:bundle}
       The bundle structure implemented by [Bundle.T]}
   \end{figure}
   \label{Bundle}

   See figure~\ref{fig:bundle} for the geometric intuition behind the bundle structure.

   \begin{dubious}
     Does the current implementation support faithful projections with a forgetful
     comparison in the base?
   \end{dubious}
*)

module type Elt_Base =
  sig
    type elt
    type base
    val compare_elt : elt -> elt -> int
    val compare_base : base -> base -> int
  end

module type Projection =
  sig
    include Elt_Base

    (* $\pi: E \to B$ *)
    val pi : elt -> base

  end

module type T =
  sig

    type t

    type elt
    type fiber = elt list
    type base

    val add : elt -> t -> t
    val of_list : elt list -> t

    (* $\pi: E \to B$ *)
    val pi : elt -> base

    (* $\pi^{-1}: B \to E$ *)
    val inv_pi : base -> t -> fiber

    val base : t -> base list

    (* $\pi^{-1}\circ\pi$ *)
    val fiber : elt -> t -> fiber

    val fibers : t -> (base * fiber) list
  end

module Make (P : Projection) : T with type elt = P.elt and type base = P.base

(* The same thing again, but with a projection that is not hardcoded, but passed
   as an argument at runtime. *)

module type Dyn =
  sig
    type t
    type elt
    type fiber = elt list
    type base
    val add : (elt -> base) -> elt -> t -> t
    val of_list : (elt -> base) -> elt list -> t
    val inv_pi : base -> t -> fiber
    val base : t -> base list
    val fiber : (elt -> base) -> elt -> t -> fiber
    val fibers : t -> (base * fiber) list
  end

module Dyn (P : Elt_Base) : Dyn with type elt = P.elt and type base = P.base

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
