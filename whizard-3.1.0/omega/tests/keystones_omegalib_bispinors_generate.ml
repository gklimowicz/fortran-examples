(* keystones_omegalib_bispinors_generate.ml --

   Copyright (C) 2019-2022 by

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

open Coupling
open Keystones

let vector_spinor_current tag =
  { tag = Printf.sprintf "vector_spinor_current__%s_ff" tag;
    keystones = [ { bra = (Majorana, 0);
                    name = Printf.sprintf "f_%sf" tag;
                    args = [G (0); F (Vector, 1); F (Majorana, 2)] };
                  { bra = (Vector, 1);
                    name = Printf.sprintf "%s_ff" tag;
                    args = [G (0); F (Majorana, 0); F (Majorana, 2)] } ] }

let scalar_spinor_current tag =
  { tag = Printf.sprintf "scalar_spinor_current__%s_ff" tag;
    keystones = [ { bra = (Majorana, 0);
                    name = Printf.sprintf "f_%sf" tag;
                    args = [G (0); F (Scalar, 1); F (Majorana, 2)] };
                  { bra = (Scalar, 1);
                    name = Printf.sprintf "%s_ff" tag;
                    args = [G (0); F (Majorana, 0); F (Majorana, 2)] } ] }

let vertices =
  List.concat
    [ List.map vector_spinor_current ["v"; "a"; "vl"; "vr"];
      List.map scalar_spinor_current ["s"; "p"; "sl"; "sr"] ]

let _ =
  Keystones.generate
    ~reps:1000 ~threshold:0.70
    ~program:"keystones_omegalib_bispinors" ~omega_module:"omega95_bispinors"
    vertices;
  exit 0
