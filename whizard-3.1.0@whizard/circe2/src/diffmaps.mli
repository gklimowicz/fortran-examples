(* circe2/diffmaps.mli --  *)
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

(* \subsection{Combined Differentiable Maps} *)

module type T =
  sig
    include Diffmap.T
    val id : ?x_min:domain -> ?x_max:domain -> codomain -> codomain -> t
  end

module type Real = T with type domain = float and type codomain = float

module type Default =
  sig

    include Real

    val power : alpha:float -> eta:float ->
      ?x_min:domain -> ?x_max:domain -> codomain -> codomain -> t
    val resonance : eta:float -> a:float ->
      ?x_min:domain -> ?x_max:domain -> codomain -> codomain -> t

  end

module Default : Default

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)


