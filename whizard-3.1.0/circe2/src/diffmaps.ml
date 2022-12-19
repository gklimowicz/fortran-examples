(* circe2/diffmaps.ml --  *)
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

module Default =
  struct

    type domain = float
    type codomain = float

    type t = 
        { encode : string;
          with_domain : x_min:domain -> x_max:domain -> t;
          x_min : domain;
          x_max : domain;
          y_min : codomain;
          y_max : codomain;
          phi : domain -> codomain;
          ihp : codomain -> domain;
          jac : domain -> float;
          caj : codomain -> float }

    let encode m = m.encode
    let with_domain m = m.with_domain

    let x_min m = m.x_min
    let x_max m = m.x_max
    let y_min m = m.y_min
    let y_max m = m.y_max

    let phi m = m.phi
    let ihp m = m.ihp
    let jac m = m.jac
    let caj m  = m.caj

    let rec id ?x_min ?x_max y_min y_max =
      let m = Diffmap.Id.create ?x_min ?x_max y_min y_max in
      let with_domain ~x_min ~x_max =
        id ~x_min ~x_max y_min y_max in
      { encode = Diffmap.Id.encode m;
        with_domain = with_domain;
        x_min = Diffmap.Id.x_min m;
        x_max = Diffmap.Id.x_max m;
        y_min = Diffmap.Id.y_min m;
        y_max = Diffmap.Id.y_max m;
        phi = Diffmap.Id.phi m;
        ihp = Diffmap.Id.ihp m;
        jac = Diffmap.Id.jac m;
        caj = Diffmap.Id.caj m }

    let rec power ~alpha ~eta ?x_min ?x_max y_min y_max =
      let m = Diffmap.Power.create ~alpha ~eta ?x_min ?x_max y_min y_max in
      let with_domain ~x_min ~x_max =
        power ~alpha ~eta ~x_min ~x_max y_min y_max in
      { encode = Diffmap.Power.encode m;
        with_domain = with_domain;
        x_min = Diffmap.Power.x_min m;
        x_max = Diffmap.Power.x_max m;
        y_min = Diffmap.Power.y_min m;
        y_max = Diffmap.Power.y_max m;
        phi = Diffmap.Power.phi m;
        ihp = Diffmap.Power.ihp m;
        jac = Diffmap.Power.jac m;
        caj = Diffmap.Power.caj m }

    let rec resonance ~eta ~a ?x_min ?x_max y_min y_max =
      let m = Diffmap.Resonance.create ~eta ~a ?x_min ?x_max y_min y_max  in
      let with_domain ~x_min ~x_max =
        resonance ~eta ~a ~x_min ~x_max y_min y_max in
      { encode = Diffmap.Resonance.encode m;
        with_domain = with_domain;
        x_min = Diffmap.Resonance.x_min m;
        x_max = Diffmap.Resonance.x_max m;
        y_min = Diffmap.Resonance.y_min m;
        y_max = Diffmap.Resonance.y_max m;
        phi = Diffmap.Resonance.phi m;
        ihp = Diffmap.Resonance.ihp m;
        jac = Diffmap.Resonance.jac m;
        caj = Diffmap.Resonance.caj m }

  end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
