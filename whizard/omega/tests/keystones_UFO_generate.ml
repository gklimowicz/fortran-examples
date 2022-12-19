(* keystones_UFO_generate.ml --

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

(* For testing Dirac equations \&c.~\ldots *)

let pslash =
  equivalent_tensors
    [| ConjSpinor; Spinor |]
    [ ("pslash", "P(-1,2)*Gamma(-1,1,2)") ]

let qed =
  equivalent_tensors
    [| ConjSpinor; Vector; Spinor |]
    [ ("qed", "Gamma(2,1,3)") ]

let axial =
  equivalent_tensors
    [| ConjSpinor; Vector; Spinor |]
    [ ("axial1", "Gamma5(1,-1)*Gamma(2,-1,3)");
      ("axial2", "-Gamma(2,1,-3)*Gamma5(-3,3)") ]

let left =
  equivalent_tensors
    [| ConjSpinor; Vector; Spinor |]
    [ ("left1", "(Identity(1,-1)+Gamma5(1,-1))*Gamma(2,-1,3)");
      ("left2", "2*ProjP(1,-1)*Gamma(2,-1,3)");
      ("left3", "Gamma(2,1,-3)*(Identity(-3,3)-Gamma5(-3,3))");
      ("left4", "2*Gamma(2,1,-3)*ProjM(-3,3)") ]

let right =
  equivalent_tensors
    [| ConjSpinor; Vector; Spinor |]
    [ ("right1", "(Identity(1,-1)-Gamma5(1,-1))*Gamma(2,-1,3)");
      ("right2", "2*ProjM(1,-1)*Gamma(2,-1,3)");
      ("right3", "Gamma(2,1,-3)*(Identity(-3,3)+Gamma5(-3,3))");
      ("right4", "2*Gamma(2,1,-3)*ProjP(-3,3)") ]

let vector_spinor_current tag =
  { tag = Printf.sprintf "vector_spinor_current__%s_ff" tag;
    keystones =
      [ { bra = (ConjSpinor, 0);
          name = Printf.sprintf "f_%sf" tag;
          args = [G (0); F (Vector, 1); F (Spinor, 2)] };
        { bra = (Vector, 1);
          name = Printf.sprintf "%s_ff" tag;
          args = [G (0); F (ConjSpinor, 0); F (Spinor, 2)] };
        { bra = (Spinor, 2);
          name = Printf.sprintf "f_f%s" tag;
          args = [G (0); F (ConjSpinor, 0); F (Vector, 1)] } ] }

let scalar =
  equivalent_tensors
    [| ConjSpinor; Scalar; Spinor |]
    [ ("scalar_current", "Identity(1,3)") ]

let pseudo =
  equivalent_tensors
    [| ConjSpinor; Scalar; Spinor |]
    [ ("pseudo_current", "Gamma5(1,3)") ]

let left_scalar =
  equivalent_tensors
    [| ConjSpinor; Scalar; Spinor |]
    [ ("left_scalar1", "Identity(1,3)-Gamma5(1,3)");
      ("left_scalar2", "2*ProjM(1,3)") ]

let right_scalar =
  equivalent_tensors
    [| ConjSpinor; Scalar; Spinor |]
    [ ("right_scalar1", "Identity(1,3)+Gamma5(1,3)");
      ("right_scalar2", "2*ProjP(1,3)") ]

let scalar_spinor_current tag =
  { tag = Printf.sprintf "scalar_spinor_current__%s_ff" tag;
    keystones =
      [ { bra = (ConjSpinor, 0);
          name = Printf.sprintf "f_%sf" tag;
          args = [G (0); F (Scalar, 1); F (Spinor, 2)] };
        { bra = (Scalar, 1);
          name = Printf.sprintf "%s_ff" tag;
          args = [G (0); F (ConjSpinor, 0); F (Spinor, 2)] };
        { bra = (Spinor, 2);
          name = Printf.sprintf "f_f%s" tag;
          args = [G (0); F (ConjSpinor, 0); F (Scalar, 1)] } ] }

let fermi_ss =
  equivalent_tensors
    [| ConjSpinor; Spinor; ConjSpinor; Spinor |]
    [ ("fermi_ss", "Identity(1,2)*Identity(3,4)");
      ("fermi_ss_f",
       "   (1/4) * Identity(1,4)*Identity(3,2)" ^
       " + (1/4) * Gamma(-1,1,4)*Gamma(-1,3,2)" ^
       " + (1/8) * Sigma(-1,-2,1,4)*Sigma(-1,-2,3,2)" ^
       " - (1/4) * Gamma(-1,1,-4)*Gamma5(-4,4)*Gamma(-1,3,-2)*Gamma5(-2,2)" ^
       " + (1/4) * Gamma5(1,4)*Gamma5(3,2)") ]

let fermi_vv =
  equivalent_tensors
    [| ConjSpinor; Spinor; ConjSpinor; Spinor |]
    [ ("fermi_vv", "Gamma(-1,1,2)*Gamma(-1,3,4)");
      ("fermi_vv_f",
       "           Identity(1,4)*Identity(3,2)" ^
       " - (1/2) * Gamma(-1,1,4)*Gamma(-1,3,2)" ^
       " - (1/2) * Gamma(-1,1,-4)*Gamma5(-4,4)*Gamma(-1,3,-2)*Gamma5(-2,2)" ^
       " -         Gamma5(1,4)*Gamma5(3,2)") ]

let fermi_tt =
  equivalent_tensors
    [| ConjSpinor; Spinor; ConjSpinor; Spinor |]
    [ ("fermi_tt1", "   Sigma(-1,-2,1,2)*Sigma(-1,-2,3,4)");
      ("fermi_tt2", " - Sigma(-1,-2,1,2)*Sigma(-2,-1,3,4)");
      ("fermi_tt3", " - Sigma(-2,-1,1,2)*Sigma(-1,-2,3,4)");
      ("fermi_tt_f",
       "   3     * Identity(1,4)*Identity(3,2)" ^
       " - (1/2) * Sigma(-1,-2,1,4)*Sigma(-1,-2,3,2)" ^
       " + 3     * Gamma5(1,4)*Gamma5(3,2)") ]

let fermi_aa =
  equivalent_tensors
    [| ConjSpinor; Spinor; ConjSpinor; Spinor |]
    [ ("fermi_aa", "Gamma5(1,-2)*Gamma(-1,-2,2)*Gamma5(3,-3)*Gamma(-1,-3,4)");
      ("fermi_aa_f",
       " -         Identity(1,4)*Identity(3,2)" ^
       " - (1/2) * Gamma(-1,1,4)*Gamma(-1,3,2)" ^
       " - (1/2) * Gamma(-1,1,-4)*Gamma5(-4,4)*Gamma(-1,3,-2)*Gamma5(-2,2)" ^
       " +         Gamma5(1,4)*Gamma5(3,2)") ]

let fermi_pp =
  equivalent_tensors
    [| ConjSpinor; Spinor; ConjSpinor; Spinor |]
    [ ("fermi_pp", "Gamma5(1,2)*Gamma5(3,4)");
      ("fermi_pp_f",
       "   (1/4) * Identity(1,4)*Identity(3,2)" ^
       " - (1/4) * Gamma(-1,1,4)*Gamma(-1,3,2)" ^
       " + (1/8) * Sigma(-1,-2,1,4)*Sigma(-1,-2,3,2)" ^
       " + (1/4) * Gamma(-1,1,-4)*Gamma5(-4,4)*Gamma(-1,3,-2)*Gamma5(-2,2)" ^
       " + (1/4) * Gamma5(1,4)*Gamma5(3,2)") ]

let fermi_ll =
  equivalent_tensors
    [| ConjSpinor; Spinor; ConjSpinor; Spinor |]
    [ ("fermi_ll",   "   Gamma(-1,1,-2)*ProjM(-2,2)*Gamma(-1,3,-4)*ProjM(-4,4)");
      ("fermi_ll_f", " - Gamma(-1,1,-2)*ProjM(-2,4)*Gamma(-1,3,-4)*ProjM(-4,2)") ]

let fermi_va =
  equivalent_tensors
    [| ConjSpinor; Spinor; ConjSpinor; Spinor |]
    [ ("fermi_va", "Gamma(-1,1,2)*Gamma5(3,-3)*Gamma(-1,-3,4)") ]

let fermi_av =
  equivalent_tensors
    [| ConjSpinor; Spinor; ConjSpinor; Spinor |]
    [ ("fermi_av", "Gamma5(1,-2)*Gamma(-1,-2,2)*Gamma(-1,3,4)") ]

let sqed =
  equivalent_tensors
    [| Scalar; Vector; Scalar |]
    [ ("sqed1", "P(2,3)-P(2,1)");
      ("sqed2", "2*P(2,3)+P(2,2)");
      ("sqed3", "-P(2,2)-2*P(2,1)") ]

let vector_scalar_current =
  { tag = "vector_scalar_current__v_ss";
    keystones =
      [ { bra = (Vector, 1);
          name = "v_ss";
          args = [G (0); F (Scalar, 2); P (2); F (Scalar, 0); P (0)] };
        { bra = (Scalar, 0);
          name = "s_vs";
          args = [G (0); F (Vector, 1); P (1); F (Scalar, 2); P (2)] } ] }

let svv_t =
  equivalent_tensors
    [| Scalar; Vector; Vector |]
    [ ("svv_t", "P(-1,2)*P(-1,3)*Metric(2,3)-P(2,3)*P(3,2)") ]

let scalar_vector_current tag =
  { tag = Printf.sprintf "transversal_vector_current__s_vv_%s" tag;
    keystones = [ { bra = (Scalar, 0);
                    name = Printf.sprintf "s_vv_%s" tag;
                    args = [G (0); F (Vector, 1); P (1); F (Vector, 2); P (2)] };
                  { bra = (Vector, 1);
                    name = Printf.sprintf "v_sv_%s" tag;
                    args = [G (0); F (Scalar, 0); P (0); F (Vector, 2); P (2)] } ] }

let gauge =
  equivalent_tensors
    [| Vector; Vector; Vector |]
    [ ("gauge", "   Metric(1,2)*P(3,1) - Metric(1,2)*P(3,2) \
                  + Metric(3,1)*P(2,3) - Metric(3,1)*P(2,1) \
                  + Metric(2,3)*P(1,2) - Metric(2,3)*P(1,3)") ]

let gauge_omega =
  { tag = "g_gg";
    keystones =
      [ { bra = (Vector, 0);
          name = "(0,1)*g_gg";
          args = [G (0); F (Vector, 1); P (1); F (Vector, 2); P (2)] } ] }

(* Note that $C^{-1}=-C$ for the charge conjugation matrix.*)
let charge_conjugate_s =
  equivalent_tensors
    [| Scalar; ConjSpinor; Spinor |]
    [ ("gamma1",    "Identity(2,3)");
      ("gamma1_cc", "C(3,-3)*Identity(-3,-2)*(-C(-2,2))");
      ("gamma1_cx", "C(3,-1)*(-C(-1,2))") ]

(* $C \gamma_5 C^{-1} = \gamma_5^T$ *)
let charge_conjugate_p =
  equivalent_tensors
    [| Scalar; ConjSpinor; Spinor |]
    [ ("gamma5",    "Gamma5(2,3)");
      ("gamma5_cc", "C(3,-3)*Gamma5(-3,-2)*(-C(-2,2))") ]

(* $C \gamma_\mu C^{-1} = - \gamma_\mu^T$ *)
let charge_conjugate_v =
  equivalent_tensors
    [| Vector; ConjSpinor; Spinor |]
    [ ("gamma_mu",    "Gamma(1,2,3)");
      ("gamma_mu_cc", "-C(3,-3)*Gamma(1,-3,-2)*(-C(-2,2))") ]

(* $C \gamma_5\gamma_\mu C^{-1} = (\gamma_5\gamma_\mu)^T$ *)
let charge_conjugate_a =
  equivalent_tensors
    [| Vector; ConjSpinor; Spinor |]
    [ ("gamma_5mu",    "Gamma5(2,-2)*Gamma(1,-2,3)");
      ("gamma_5mu_cc", "C(3,-3)*Gamma5(-3,-1)*Gamma(1,-1,-2)*(-C(-2,2))") ]

(* $C \sigma_{\mu\nu} C^{-1} = - \sigma_{\mu\nu}^T$ *)
let charge_conjugate_t =
  equivalent_tensors
    [| Vector; Vector; ConjSpinor; Spinor |]
    [ ("sigma_munu",    "Sigma(1,2,3,4)");
      ("sigma_munu_cc", "-C(4,-4)*Sigma(1,2,-4,-3)*(-C(-3,3))") ]

(* $C \gamma_\mu \gamma_\nu C^{-1} = \gamma_\nu^T \gamma_\mu^T$ *)
let charge_conjugate_vv =
  equivalent_tensors
    [| Vector; Vector; ConjSpinor; Spinor |]
    [ ("gamma_mu_nu",    "Gamma(1,3,-1)*Gamma(2,-1,4)");
      ("gamma_mu_nu_cc", "C(4,-4)*Gamma(2,-4,-1)*Gamma(1,-1,-3)*(-C(-3,3))") ]

let empty = { tag = "empty"; keystones = [ ] }

let vertices =
  [ (pslash, empty);
    (qed, vector_spinor_current "v");
    (axial, vector_spinor_current "a");
    (left, vector_spinor_current "vl");
    (right, vector_spinor_current "vr");
    (scalar, scalar_spinor_current "s");
    (pseudo, scalar_spinor_current "p");
    (left_scalar, scalar_spinor_current "sl");
    (right_scalar, scalar_spinor_current "sr");
    (sqed, vector_scalar_current);
    (fermi_ss, empty);
    (fermi_vv, empty);
    (fermi_tt, empty);
    (fermi_aa, empty);
    (fermi_pp, empty);
    (fermi_ll, empty);
    (fermi_va, empty);
    (fermi_av, empty);
    (svv_t, scalar_vector_current "t");
    (gauge, gauge_omega);
    (charge_conjugate_s, empty);
    (charge_conjugate_p, empty);
    (charge_conjugate_v, empty);
    (charge_conjugate_a, empty);
    (charge_conjugate_t, empty);
    (charge_conjugate_vv, empty) ]

let parse_propagator (p_tag, p_omega, p_spins, numerator, denominator) =
  let p =
    UFO.Propagator.of_propagator_UFO
      { UFO.Propagator_UFO.name = p_tag;
        UFO.Propagator_UFO.numerator = UFOx.Lorentz.of_string numerator;
        UFO.Propagator_UFO.denominator = UFOx.Lorentz.of_string denominator } in
  { p_tag; p_omega; p_spins;
    p_propagator = p }

let default_denominator =
  "P('mu', id) * P('mu', id) - Mass(id) * Mass(id) \
   + complex(0,1) * Mass(id) * Width(id)"

let scalar_propagator =
  ( "scalar", "pr_phi", (Scalar, Scalar),
    "1",
    default_denominator )

let spinor_propagator =
  ( "spinor", "pr_psi", (ConjSpinor, Spinor),
    "Gamma('mu', 1, 2) * P('mu', id) + Mass(id) * Identity(1, 2)",
    default_denominator )

let conjspinor_propagator =
  ( "conjspinor", "pr_psibar", (ConjSpinor, Spinor),
    "Gamma('mu', 1, 2) * P('mu', id) + Mass(id) * Identity(1, 2)",
    default_denominator )

let feynman_propagator =
  ( "feynman", "pr_feynman", (Vector, Vector),
    " - Metric(1, 2)",
    "P('mu', id) * P('mu', id)" )

let gauge_propagator =
  ( "gauge_propagator", "pr_gauge", (Vector, Vector),
    " - Metric(1, 2) + (1 - 42) * P(1,id) * P(2,id) / " ^
      "( P('mu', id) * P('mu', id) )",
    "P('mu', id) * P('mu', id)" )

let rxi_propagator =
  ( "rxi_propagator", "pr_rxi", (Vector, Vector),
    " - Metric(1, 2) + (1 - 42) * P(1,id) * P(2,id) / " ^
      "( P('mu', id) * P('mu', id) - 42 * Mass(id)**2 )",
    default_denominator )

let unitarity_propagator =
  ( "unitarity", "pr_unitarity", (Massive_Vector, Massive_Vector),
    "- Metric(1, 2) + Metric(1,'mu')*P('mu', id)*P(2, id)/Mass(id)**2",
    default_denominator )

let tensor_propagator =    
  ( "tensor", "pr_tensor", (Tensor_2, Tensor_2),
    "  1/2 * (Metric(1001,1002) - P(1001,id)*P(1002,id)/Mass(id)**2) \
           * (Metric(2001,2002) - P(2001,id)*P(2002,id)/Mass(id)**2) \
     + 1/2 * (Metric(1001,2002) - P(1001,id)*P(2002,id)/Mass(id)**2) \
           * (Metric(2001,1002) - P(2001,id)*P(1002,id)/Mass(id)**2) \
     - 1/3 * (Metric(1001,2001) - P(1001,id)*P(2001,id)/Mass(id)**2) \
           * (Metric(1002,2002) - P(1002,id)*P(2002,id)/Mass(id)**2) ",
    default_denominator )

let tensor_propagator_51_52 =
  ( "tensor_51_52", "pr_tensor", (Tensor_2, Tensor_2),
    "  1/2 * (Metric( 1, 2) - P( 1,id)*P( 2,id)/Mass(id)**2) \
           * (Metric(51,52) - P(51,id)*P(52,id)/Mass(id)**2) \
     + 1/2 * (Metric( 1,52) - P( 1,id)*P(52,id)/Mass(id)**2) \
           * (Metric(51, 2) - P(51,id)*P( 2,id)/Mass(id)**2) \
     - 1/3 * (Metric( 1,51) - P( 1,id)*P(51,id)/Mass(id)**2) \
           * (Metric( 2,52) - P( 2,id)*P(52,id)/Mass(id)**2) ",
    default_denominator )

let propagators =
  List.map
    parse_propagator
    [ scalar_propagator;
      spinor_propagator;
      feynman_propagator;
      gauge_propagator;
      rxi_propagator;
      unitarity_propagator;
      tensor_propagator;
      tensor_propagator_51_52 ]

let conjugate_propagators =
  List.map
    (fun p -> transpose (parse_propagator p))
    [ conjspinor_propagator ]

let all_propagators = propagators @ conjugate_propagators

let _ =
  generate_ufo
    ~reps:1000 ~threshold:0.70 ~program:"keystones_UFO"
    "fusions_UFO" vertices all_propagators;
  exit 0
