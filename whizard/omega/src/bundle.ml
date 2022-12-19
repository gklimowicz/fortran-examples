(* bundle.ml --

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

module type Elt_Base =
  sig
    type elt
    type base
    val compare_elt : elt -> elt -> int
    val compare_base : base -> base -> int
  end

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

module Dyn (P : Elt_Base) =
  struct

    type elt = P.elt
    type base = P.base

    type fiber = elt list

    module InvPi = Map.Make (struct type t = P.base let compare = P.compare_base end)
    module Fiber = Set.Make (struct type t = P.elt let compare = P.compare_elt end)

    type t = Fiber.t InvPi.t

    let add pi element fibers =
      let base = pi element in
      let fiber =
        try InvPi.find base fibers with Not_found -> Fiber.empty in
      InvPi.add base (Fiber.add element fiber) fibers

    let of_list pi list =
      List.fold_right (add pi) list InvPi.empty

    let fibers bundle =
      InvPi.fold
        (fun base fiber acc -> (base, Fiber.elements fiber) :: acc) bundle []

    let base bundle =
      InvPi.fold
        (fun base fiber acc -> base :: acc) bundle []
      
    let inv_pi base bundle =
      try
        Fiber.elements (InvPi.find base bundle)
      with
      | Not_found -> []

    let fiber pi elt bundle =
      inv_pi (pi elt) bundle

  end

module type Projection =
  sig
    include Elt_Base
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
    val pi : elt -> base
    val inv_pi : base -> t -> fiber
    val base : t -> base list
    val fiber : elt -> t -> fiber
    val fibers : t -> (base * fiber) list
  end

module Make (P : Projection) =
  struct

    module D = Dyn (P)

    type elt = D.elt
    type base = D.base
    type fiber = D.fiber
    type t = D.t

    let pi = P.pi

    let add = D.add pi
    let of_list = D.of_list pi
    let base = D.base
    let inv_pi = D.inv_pi
    let fibers = D.fibers

    let fiber elt bundle =
      inv_pi (pi elt) bundle

  end

(*i
module Test = Make (struct
  type fiber = int
  type base = int
  let compare_fiber = compare
  let compare_base = compare
  let pi = abs
end)

let sample = [-1; -4; 7; -8; 9; 42; -137; -42; 42; 4; 1; -9]

Test.fibers (Test.classify sample);;
i*)

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
