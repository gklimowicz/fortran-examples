(* omega_SYM.ml --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       Christian Speckner <cnspeckn@googlemail.com>

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


module SYM = 
  struct

    open Coupling

    let options = Options.empty
    let caveats () = []

    let nc = 3

    type flavor = 
      | Q of int | SQ of int
      | G of int | SG of int
      | Phi

    let generations = ThoList.range 1 1

    let generations_pairs =
      List.map
        (function [a;b] -> (a, b)
          | _ -> failwith "omega_SYM.generations_pairs")
        (Product.power 2 generations)

    let generations_triples =
      List.map
        (function [a;b;c] -> (a, b, c)
          | _ -> failwith "omega_SYM.generations_triples")
        (Product.power 3 generations)

    let generations_quadruples =
      List.map
        (function [a;b;c;d] -> (a, b, c, d)
          | _ -> failwith "omega_SYM.generations_quadruples")
        (Product.power 4 generations)

    let external_flavors () =
      [ "Quarks", List.map (fun i -> Q i) generations;
        "Anti-Quarks", List.map (fun i -> Q (-i)) generations;
        "SQuarks", List.map (fun i -> SQ i) generations;
        "Anti-SQuarks", List.map (fun i -> SQ (-i)) generations;
        "Gluons", List.map (fun i -> G i) generations;
        "SGluons", List.map (fun i -> SG i) generations;
        "Other", [Phi]]

    let flavors () =
      ThoList.flatmap snd (external_flavors ())

    type gauge = unit
    type constant =
      | G_saa of int * int
      | G_saaa of int * int * int
      | G3 of int * int * int
      | I_G3 of int * int * int
      | G4 of int * int * int * int

    type orders = unit
    let orders = function
      | _ -> ()

    let lorentz = function
      | Q i ->
          if i > 0 then
            Spinor
          else if i < 0 then
            ConjSpinor
          else
            invalid_arg "SYM.lorentz (Q 0)"
      | SQ _ | Phi -> Scalar
      | G _ -> Vector
      | SG _ -> Majorana

    let color = function 
      | Q i | SQ i ->
          Color.SUN (if i > 0 then nc else if i < 0 then -nc else invalid_arg "SYM.color (Q 0)")
      | G _ | SG _ -> Color.AdjSUN nc
      | Phi -> Color.Singlet

    let nc () = nc

    let propagator = function
      | Q i ->
          if i > 0 then
            Prop_Spinor
          else if i < 0 then
            Prop_ConjSpinor
          else
            invalid_arg "SYM.lorentz (Q 0)"
      | SQ _ | Phi -> Prop_Scalar
      | G _ -> Prop_Feynman
      | SG _ -> Prop_Majorana

(*i let propagator _ =
      Only_Insertion
i*)

    let width _ = Timelike
    let goldstone _ = None

    let conjugate = function
      | Q i -> Q (-i)
      | SQ i -> SQ (-i)
      | (G _ | SG _ | Phi) as p -> p

    let fermion = function
      | Q i ->
          if i > 0 then
            1
          else if i < 0 then
            -1
          else
            invalid_arg "SYM.fermion (Q 0)"
      | SQ _ | G _ | Phi -> 0
      | SG _ -> 2

    module Ch = Charges.Null
    let charges _ = ()

    module F = Modeltools.Fusions (struct
      type f = flavor
      type c = constant
      let compare = compare
      let conjugate = conjugate
    end)

    let quark_current =
      List.map
        (fun (i, j, k) ->
          ((Q (-i), G j, Q k), FBF (-1, Psibar, V, Psi), G3 (i, j, k)))
        generations_triples

    let squark_current =
      List.map
        (fun (i, j, k) ->
          ((G j, SQ i, SQ (-k)), Vector_Scalar_Scalar 1, G3 (i, j, k)))
        generations_triples

    let three_gluon =
      List.map
        (fun (i, j, k) ->
          ((G i, G j, G k), Gauge_Gauge_Gauge 1, I_G3 (i, j, k)))
        generations_triples

    let gluon2_phi =
      List.map
        (fun (i, j) ->
          ((Phi, G i, G j), Dim5_Scalar_Gauge2 1, G_saa (i, j)))
        generations_pairs

    let vertices3 =
      quark_current @ squark_current @ three_gluon @ gluon2_phi
                                                       
    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]

    let squark_seagull =
      List.map
        (fun (i, j, k, l) ->
          ((SQ i, SQ (-j), G k, G l), Scalar2_Vector2 1, G4 (i, j, k, l)))
       generations_quadruples

    let four_gluon =
      List.map
        (fun (i, j, k, l) ->
          ((G i, G j, G k, G l), gauge4, G4 (i, j, k, l)))
       generations_quadruples

(*i
    let gluon3_phi =
      List.map
        (fun (i, j, k) ->
          ((Phi, G i, G j, G k), Dim6_Scalar_Gauge3 1, G_saaa (i, j, k)))
        generations_triples
i*)
(* \begin{dubious}
     We need at least a [Dim6_Scalar_Gauge3] vertex to support this.
   \end{dubious} *)
    let gluon3_phi =
      []

    let vertices4 =
      squark_seagull @ four_gluon @ gluon3_phi

    let vertices () = 
      (vertices3, vertices4, [])

    let table = F.of_vertices (vertices ())
    let fuse2 = F.fuse2 table
    let fuse3 = F.fuse3 table
    let fuse = F.fuse table
    let max_degree () = 4
        
    let parameters () = { input = []; derived = []; derived_arrays = [] }

    let invalid_flavor s =
      invalid_arg ("omega_SYM.flavor_of_string: " ^ s)

    let flavor_of_string s =
      let l = String.length s in
      if l < 2 then
        invalid_flavor s
      else if l = 2 then
        if String.sub s 0 1 = "q" then
          Q (int_of_string (String.sub s 1 1))
        else if String.sub s 0 1 = "Q" then
          Q (- (int_of_string (String.sub s 1 1)))
        else if String.sub s 0 1 = "g" then
          G (int_of_string (String.sub s 1 1))
        else
          invalid_flavor s
      else if l = 3 then
        if s = "phi" then
          Phi
        else if String.sub s 0 2 = "sq" then
          SQ (int_of_string (String.sub s 2 1))
        else if String.sub s 0 2 = "sQ" then
          SQ (- (int_of_string (String.sub s 2 1)))
        else if String.sub s 0 2 = "sg" then
          SG (int_of_string (String.sub s 2 1))
        else
          invalid_flavor s
      else
        invalid_flavor s

    let flavor_to_string = function
      | Q i ->
          if i > 0 then
            "q" ^ string_of_int i
          else if i < 0 then
            "Q" ^ string_of_int (-i)
          else
            invalid_arg "SYM.flavor_to_string (Q 0)"
      | SQ i -> 
          if i > 0 then
            "sq" ^ string_of_int i
          else if i < 0 then
            "sQ" ^ string_of_int (-i)
          else
            invalid_arg "SYM.flavor_to_string (SQ 0)"
      | G i -> "g" ^ string_of_int i
      | SG i -> "sg" ^ string_of_int i
      | Phi -> "phi"

    let flavor_to_TeX = function
      | Q i ->
          if i > 0 then
            "q_{" ^ string_of_int i ^ "}"
          else if i < 0 then
            "{\bar q}_{" ^ string_of_int (-i) ^ "}"
          else
            invalid_arg "SYM.flavor_to_string (Q 0)"
      | SQ i -> 
          if i > 0 then
            "{\tilde q}_{" ^ string_of_int i ^ "}"
          else if i < 0 then
            "{\bar{\tilde q}}_{" ^ string_of_int (-i) ^ "}"
          else
            invalid_arg "SYM.flavor_to_string (SQ 0)"
      | G i -> "g_{" ^ string_of_int i ^ "}"
      | SG i -> "{\tilde g}_{" ^ string_of_int i ^ "}"
      | Phi -> "phi"

    let flavor_symbol = function
      | Q i ->
          if i > 0 then
            "q" ^ string_of_int i
          else if i < 0 then
            "qbar" ^ string_of_int (-i)
          else
            invalid_arg "SYM.flavor_to_string (Q 0)"
      | SQ i -> 
          if i > 0 then
            "sq" ^ string_of_int i
          else if i < 0 then
            "sqbar" ^ string_of_int (-i)
          else
            invalid_arg "SYM.flavor_to_string (SQ 0)"
      | G i -> "g" ^ string_of_int i
      | SG i -> "sg" ^ string_of_int i
      | Phi -> "phi"

    let gauge_symbol () =
      failwith "omega_SYM.gauge_symbol: internal error"

    let pdg _ = 0
    let mass_symbol _ = "0.0_default"
    let width_symbol _ = "0.0_default"

    let string_of_int_list int_list = 
      "(" ^ String.concat "," (List.map string_of_int int_list) ^ ")"

    let constant_symbol = function
      | G_saa (i, j) -> "g_saa" ^ string_of_int_list [i;j]
      | G_saaa (i, j, k) -> "g_saaa" ^ string_of_int_list [i;j;k]
      | G3 (i, j, k) -> "g3" ^ string_of_int_list [i;j;k]
      | I_G3 (i, j, k) -> "ig3" ^ string_of_int_list [i;j;k]
      | G4 (i, j, k, l) -> "g4" ^ string_of_int_list [i;j;k;l]

  end

module O = Omega.Mixed23(Targets.Fortran_Majorana)(SYM)
let _ = O.main ()


(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
