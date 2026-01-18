(* keystones_omegalib_generate.ml --

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
    keystones = [ { bra = (ConjSpinor, 0);
                    name = Printf.sprintf "f_%sf" tag;
                    args = [G (0); F (Vector, 1); F (Spinor, 2)] };
                  { bra = (Vector, 1);
                    name = Printf.sprintf "%s_ff" tag;
                    args = [G (0); F (ConjSpinor, 0); F (Spinor, 2)] };
                  { bra = (Spinor, 2);
                    name = Printf.sprintf "f_f%s" tag;
                    args = [G (0); F (ConjSpinor, 0); F (Vector, 1)] } ] }

let scalar_spinor_current tag =
  { tag = Printf.sprintf "scalar_spinor_current__%s_ff" tag;
    keystones = [ { bra = (ConjSpinor, 0);
                    name = Printf.sprintf "f_%sf" tag;
                    args = [G (0); F (Scalar, 1); F (Spinor, 2)] };
                  { bra = (Scalar, 1);
                    name = Printf.sprintf "%s_ff" tag;
                    args = [G (0); F (ConjSpinor, 0); F (Spinor, 2)] };
                  { bra = (Spinor, 2);
                    name = Printf.sprintf "f_f%s" tag;
                    args = [G (0); F (ConjSpinor, 0); F (Scalar, 1)] } ] }

(* NB: the vertex is anti-symmetric in the scalars and we need to
   use a cyclic permutation. *)
let vector_scalar_current =
  { tag = "vector_scalar_current__v_ss";
    keystones = [ { bra = (Vector, 0);
                    name = "v_ss";
                    args = [G (0); F (Scalar, 1); P (1); F (Scalar, 2); P (2)] };
                  { bra = (Scalar, 2);
                    name = "s_vs";
                    args = [G (0); F (Vector, 0); P (0); F (Scalar, 1); P (1)] } ] }

let scalar_vector_current tag =
  { tag = Printf.sprintf "transversal_vector_current__s_vv_%s" tag;
    keystones = [ { bra = (Scalar, 0);
                    name = Printf.sprintf "s_vv_%s" tag;
                    args = [G (0); F (Vector, 1); P (1); F (Vector, 2); P (2)] };
                  { bra = (Vector, 1);
                    name = Printf.sprintf "v_sv_%s" tag;
                    args = [G (0); F (Scalar, 0); P (0); F (Vector, 2); P (2)] } ] }

let vertices =
  List.concat
    [ List.map vector_spinor_current ["v"; "a"; "vl"; "vr"];
      List.map scalar_spinor_current ["s"; "p"; "sl"; "sr"];
      [ vector_scalar_current ];
      List.map scalar_vector_current ["t"; "6D"; "6DP"] ]

let _ =
  Keystones.generate
    ~program:"keystones_omegalib" ~reps:1000 ~threshold:0.70 vertices;
  exit 0
