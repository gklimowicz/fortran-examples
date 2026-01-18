(* omega_SM_top.ml --

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

(* \thocwmodulesection{SM with charge $4/3$ top} *)

module type SM_flags =
  sig
    val include_anomalous : bool
    val k_matrix : bool
  end

module SM_no_anomalous : SM_flags =
  struct
    let include_anomalous = false
    let k_matrix = false
  end

module SM_gluons : SM_flags =
  struct
    let include_anomalous = false
    let k_matrix = false
  end

module Anomtop (Flags : SM_flags) =
  struct

    open Coupling

    let default_width = ref Timelike
    let use_fudged_width = ref false

    let options = Options.create
      [ "constant_width", Arg.Unit (fun () -> default_width := Constant),
        "use constant width (also in t-channel)";
        "fudged_width", Arg.Set use_fudged_width,
        "use fudge factor for charge particle width";
        "custom_width", Arg.String (fun f -> default_width := Custom f),
        "use custom width";
        "cancel_widths", Arg.Unit (fun () -> default_width := Vanishing),
        "use vanishing width";
        "running_width", Arg.Unit (fun () -> default_width := Running),
        "use running width" ]
    let caveats () = []

    type matter_field = L of int | N of int | U of int | D of int 
    type gauge_boson = Ga | Wp | Wm | Z | Gl   
    type other = Phip | Phim | Phi0 | H 
    type flavor = M of matter_field | G of gauge_boson | O of other

    let matter_field f = M f
    let gauge_boson f = G f
    let other f = O f

    type field =
      | Matter of matter_field
      | Gauge of gauge_boson
      | Other of other

    let field = function
      | M f -> Matter f
      | G f -> Gauge f
      | O f -> Other f

    type gauge = unit

    let gauge_symbol () =
      failwith "Models.Anomtop.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; N n; U n; D n ]

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Gauge Bosons", List.map gauge_boson [Ga; Z; Wp; Wm; Gl];
        "Higgs", [O H];
        "Goldstone Bosons", List.map other [Phip; Phim; Phi0] ]

    let flavors () = ThoList.flatmap snd (external_flavors ())

    let spinor n =
      if n >= 0 then
        Spinor
      else
        ConjSpinor

    let lorentz = function
      | M f ->
          begin match f with
          | L n -> spinor n | N n -> spinor n
          | U n -> spinor n | D n -> spinor n
          end
      | G f ->
          begin match f with
          | Ga | Gl -> Vector
          | Wp | Wm | Z -> Massive_Vector
          end
      | O f -> 
          Scalar

    let color = function 
      | M (U n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D n) -> Color.SUN  (if n > 0 then 3 else -3)
      | G Gl -> Color.AdjSUN 3
      | _ -> Color.Singlet

    let nc () = 3

    let prop_spinor n =
      if n >= 0 then
        Prop_Spinor
      else
        Prop_ConjSpinor

    let propagator = function
      | M f ->
          begin match f with
          | L n -> prop_spinor n | N n -> prop_spinor n
          | U n -> prop_spinor n | D n -> prop_spinor n
          end
      | G f ->
          begin match f with
          | Ga | Gl -> Prop_Feynman
          | Wp | Wm | Z -> Prop_Unitarity
          end
      | O f ->
          begin match f with
          | Phip | Phim | Phi0 -> Only_Insertion
          | H  -> Prop_Scalar
          end

(* Optionally, ask for the fudge factor treatment for the widths of
   charged particles.  Currently, this only applies to $W^\pm$ and top. *)

    let width f =
      if !use_fudged_width then
        match f with
        | G Wp | G Wm | M (U 3) | M (U (-3)) -> Fudged
        | _ -> !default_width
      else
        !default_width

    let goldstone = function
      | G f ->
          begin match f with
          | Wp -> Some (O Phip, Coupling.Integer 1)
          | Wm -> Some (O Phim, Coupling.Integer 1)
          | Z -> Some (O Phi0, Coupling.Integer 1)
          | _ -> None
          end
      | _ -> None

    let conjugate = function
      | M f ->
          M (begin match f with
          | L n -> L (-n) | N n -> N (-n)
          | U n -> U (-n) | D n -> D (-n)
          end)
      | G f ->
          G (begin match f with
          | Gl -> Gl | Ga -> Ga | Z -> Z 
          | Wp -> Wm | Wm -> Wp
          end)
      | O f ->
          O (begin match f with
          | Phip -> Phim | Phim -> Phip | Phi0 -> Phi0
          | H -> H 
          end)

    let fermion = function
      | M f ->
          begin match f with
          | L n -> if n > 0 then 1 else -1
          | N n -> if n > 0 then 1 else -1
          | U n -> if n > 0 then 1 else -1
          | D n -> if n > 0 then 1 else -1
          end
      | G f ->
          begin match f with
          | Gl | Ga | Z | Wp | Wm -> 0
          end
      | O _ -> 0

    (* Electrical charge, lepton number, baryon number.  We could avoid the
       rationals altogether by multiplying the first and last by 3 \ldots *)

    module Ch = Charges.QQ
    let ( // ) = Algebra.Small_Rational.make

    let generation' = function
      |  1 -> [ 1//1;  0//1;  0//1]
      |  2 -> [ 0//1;  1//1;  0//1]
      |  3 -> [ 0//1;  0//1;  1//1]
      | -1 -> [-1//1;  0//1;  0//1]
      | -2 -> [ 0//1; -1//1;  0//1]
      | -3 -> [ 0//1;  0//1; -1//1]
      |  n -> invalid_arg ("SM_top.generation': " ^ string_of_int n)

    let generation f =
      match f with
      | M (L n | N n | U n | D n) -> generation' n
      | G _ | O _ -> [0//1; 0//1; 0//1]

    let charge = function
      | M f ->
          begin match f with
          | L n -> if n > 0 then -1//1 else  1//1
          | N n -> 0//1
          | U (1|2) -> 2//3 
          | U ((-1)|(-2)) -> -2//3
          | U 3 -> -4//3
          | U (-3) -> 4//3
          | U n -> invalid_arg ("SM_top.charge: up quark " ^ string_of_int n)
          | D n -> if n > 0 then -1//3 else  1//3
          end
      | G f ->
          begin match f with
          | Gl | Ga | Z -> 0//1
          | Wp ->  1//1
          | Wm -> -1//1
          end
      | O f ->
          begin match f with
          | H | Phi0 ->  0//1
          | Phip ->  1//1
          | Phim -> -1//1
          end

    let lepton = function
      | M f ->
          begin match f with
          | L n | N n -> if n > 0 then 1//1 else -1//1
          | U _ | D _ -> 0//1
          end
      | G _ | O _ -> 0//1

    let baryon = function
      | M f ->
          begin match f with
          | L _ | N _ -> 0//1
          | U n | D n -> if n > 0 then 1//1 else -1//1
          end
      | G _ | O _ -> 0//1

    let charges f = 
      [ charge f; lepton f; baryon f] @ generation f

    type constant =
      | Unit | Pi | Alpha_QED | Sin2thw
      | Sinthw | Costhw | E | G_weak | Vev 
      | Q_lepton | Q_up | Q_down | Q_top | G_CC 
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down | G_NC_top
      | I_Q_W | I_G_ZWW
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ
      | G_Htt | G_Hbb | G_Hcc | G_Htautau | G_H3 | G_H4
      | Gs | I_Gs | G2 
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let input_parameters =
      []

    let derived_parameters =
      []

    let derived_parameter_arrays =
      []

    let parameters () =
      { input = input_parameters;
        derived = derived_parameters;
        derived_arrays = derived_parameter_arrays }

    module F = Modeltools.Fusions (struct
      type f = flavor
      type c = constant
      let compare = compare
      let conjugate = conjugate
    end)

(* \begin{equation}
     \mathcal{L}_{\textrm{EM}} =
        - e \sum_i q_i \bar\psi_i\fmslash{A}\psi_i
   \end{equation} *)

    let mgm ((m1, g, m2), fbf, c) = ((M m1, G g, M m2), fbf, c)

    let electromagnetic_currents' n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]
     let em_up_type_currents = 
       List.map mgm
        [ ((U (-1), Ga, U 1), FBF (1, Psibar, V, Psi), Q_up);  
          ((U (-2), Ga, U 2), FBF (1, Psibar, V, Psi), Q_up);  
          ((U (-3), Ga, U 3), FBF (1, Psibar, V, Psi), Q_top)]
     let electromagnetic_currents = 
       ThoList.flatmap electromagnetic_currents' [1;2;3] @ 
       em_up_type_currents

    let color_currents n =
        List.map mgm
          [ ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs);
            ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs) ]

(* \begin{equation}
     \mathcal{L}_{\textrm{NC}} =
        - \frac{g}{2\cos\theta_W}
            \sum_i \bar\psi_i\fmslash{Z}(g_V^i-g_A^i\gamma_5)\psi_i
   \end{equation} *)

    let neutral_currents' n =
      List.map mgm
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 
     let neutral_up_type_currents = 
       List.map mgm
        [ ((U (-1), Z, U 1), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((U (-2), Z, U 2), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((U (-3), Z, U 3), FBF (1, Psibar, VA, Psi), G_NC_top) ]        
     let neutral_currents = 
       ThoList.flatmap neutral_currents' [1;2;3] @
       neutral_up_type_currents
       
(* \begin{equation}
     \mathcal{L}_{\textrm{CC}} =
        - \frac{g}{2\sqrt2} \sum_i \bar\psi_i
               (T^+\fmslash{W}^+ + T^-\fmslash{W}^-)(1-\gamma_5)\psi_i 
   \end{equation} *)

    let charged_currents' n =
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC) ]
    let charged_up_currents = 
      List.map mgm
        [ ((U (-1), Wp, D 1), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-2), Wp, D 2), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-3), Wm, D 3), FBF (1, Psibar, VL, Psi), G_CC);
          ((D (-1), Wm, U 1), FBF (1, Psibar, VL, Psi), G_CC);
          ((D (-2), Wm, U 2), FBF (1, Psibar, VL, Psi), G_CC);
          ((D (-3), Wp, U 3), FBF (1, Psibar, VL, Psi), G_CC) ]
    let charged_currents = 
      ThoList.flatmap charged_currents' [1;2;3] @      
      charged_up_currents

    let yukawa =
      [ ((M (U (-3)), O H, M (U 3)), FBF (1, Psibar, S, Psi), G_Htt);
        ((M (D (-3)), O H, M (D 3)), FBF (1, Psibar, S, Psi), G_Hbb);
        ((M (U (-2)), O H, M (U 2)), FBF (1, Psibar, S, Psi), G_Hcc);
        ((M (L (-3)), O H, M (L 3)), FBF (1, Psibar, S, Psi), G_Htautau) ]
      
(* \begin{equation}
     \mathcal{L}_{\textrm{TGC}} =
        - e \partial_\mu A_\nu W_+^\mu W_-^\nu + \ldots
        - e \cot\theta_w  \partial_\mu Z_\nu W_+^\mu W_-^\nu + \ldots
   \end{equation} *)

    let tgc ((g1, g2, g3), t, c) = ((G g1, G g2, G g3), t, c)

    let triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs)]

    let qgc ((g1, g2, g3, g4), t, c) = ((G g1, G g2, G g3, G g4), t, c)

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          (Gl, Gl, Gl, Gl), gauge4, G2 ]

    let gauge_higgs =
      [ ((O H, G Wp, G Wm), Scalar_Vector_Vector 1, G_HWW);
        ((O H, G Z, G Z), Scalar_Vector_Vector 1, G_HZZ) ]

    let gauge_higgs4 =
      [ (O H, O H, G Wp, G Wm), Scalar2_Vector2 1, G_HHWW;
        (O H, O H, G Z, G Z), Scalar2_Vector2 1, G_HHZZ ]

    let higgs =
      [ (O H, O H, O H), Scalar_Scalar_Scalar 1, G_H3 ]

    let higgs4 =
      [ (O H, O H, O H, O H), Scalar4 1, G_H4 ]

    let goldstone_vertices =
      [ ((O Phi0, G Wm, G Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phip, G Ga, G Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phip, G Z, G Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phim, G Wp, G Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phim, G Wp, G Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (electromagnetic_currents @
       ThoList.flatmap color_currents [1;2;3] @
       neutral_currents @
       charged_currents @
       yukawa @ triple_gauge @
       gauge_higgs @ higgs @ goldstone_vertices)

    let vertices4 =
      quartic_gauge @ gauge_higgs4 @ higgs4

    let vertices () = (vertices3, vertices4, [])

(* For efficiency, make sure that [F.of_vertices vertices] is
   evaluated only once. *)

    let table = F.of_vertices (vertices ())
    let fuse2 = F.fuse2 table
    let fuse3 = F.fuse3 table
    let fuse = F.fuse table
    let max_degree () = 4

    let flavor_of_string = function
      | "e-" -> M (L 1) | "e+" -> M (L (-1))
      | "mu-" -> M (L 2) | "mu+" -> M (L (-2))
      | "tau-" -> M (L 3) | "tau+" -> M (L (-3))
      | "nue" -> M (N 1) | "nuebar" -> M (N (-1))
      | "numu" -> M (N 2) | "numubar" -> M (N (-2))
      | "nutau" -> M (N 3) | "nutaubar" -> M (N (-3))
      | "u" -> M (U 1) | "ubar" -> M (U (-1))
      | "c" -> M (U 2) | "cbar" -> M (U (-2))
      | "t" -> M (U 3) | "tbar" -> M (U (-3))
      | "d" -> M (D 1) | "dbar" -> M (D (-1))
      | "s" -> M (D 2) | "sbar" -> M (D (-2))
      | "b" -> M (D 3) | "bbar" -> M (D (-3))
      | "g" | "gl" -> G Gl
      | "A" -> G Ga | "Z" | "Z0" -> G Z
      | "W+" -> G Wp | "W-" -> G Wm
      | "H" -> O H
      | _ -> invalid_arg "Models.Anomtop.flavor_of_string" 

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Models.Anomtop.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Models.Anomtop.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Models.Anomtop.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Models.Anomtop.flavor_to_string: invalid down type quark"
          end
      | G f ->
          begin match f with
          | Gl -> "g"
          | Ga -> "A" | Z -> "Z"
          | Wp -> "W+" | Wm -> "W-"
          end
      | O f ->
          begin match f with
          | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
          | H -> "H" 
          end

    let flavor_to_TeX = function
      | M f ->
          begin match f with
          | L 1 -> "e^-" | L (-1) -> "e^+"
          | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
          | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
          | L _ -> invalid_arg
                "Models.Anomtop.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Models.Anomtop.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Models.Anomtop.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Models.Anomtop.flavor_to_TeX: invalid down type quark"
          end
      | G f ->
          begin match f with
          | Gl -> "g"
          | Ga -> "\\gamma" | Z -> "Z"
          | Wp -> "W^+" | Wm -> "W^-"
          end
      | O f ->
          begin match f with
          | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
          | H -> "H" 
          end

    let flavor_symbol = function
      | M f ->
          begin match f with
          | L n when n > 0 -> "l" ^ string_of_int n
          | L n -> "l" ^ string_of_int (abs n) ^ "b"
          | N n when n > 0 -> "n" ^ string_of_int n
          | N n -> "n" ^ string_of_int (abs n) ^ "b"
          | U n when n > 0 -> "u" ^ string_of_int n
          | U n -> "u" ^ string_of_int (abs n) ^ "b"
          | D n when n > 0 ->  "d" ^ string_of_int n
          | D n -> "d" ^ string_of_int (abs n) ^ "b"
          end
      | G f ->
          begin match f with
          | Gl -> "gl"
          | Ga -> "a" | Z -> "z"
          | Wp -> "wp" | Wm -> "wm"
          end
      | O f ->
          begin match f with
          | Phip -> "pp" | Phim -> "pm" | Phi0 -> "p0" 
          | H -> "h" 
          end

    let pdg = function
      | M f ->
          begin match f with
          | L n when n > 0 -> 9 + 2*n
          | L n -> - 9 + 2*n
          | N n when n > 0 -> 10 + 2*n
          | N n -> - 10 + 2*n
          | U n when n > 0 -> 2*n
          | U n -> 2*n
          | D n when n > 0 -> - 1 + 2*n
          | D n -> 1 + 2*n
          end
      | G f ->
          begin match f with
          | Gl -> 21
          | Ga -> 22 | Z -> 23
          | Wp -> 24 | Wm -> (-24)
          end
      | O f ->
          begin match f with
          | Phip | Phim -> 27 | Phi0 -> 26
          | H -> 25 
          end

    let mass_symbol f = 
      "mass(" ^ string_of_int (abs (pdg f)) ^ ")"

    let width_symbol f =
      "width(" ^ string_of_int (abs (pdg f)) ^ ")"

    let constant_symbol = function
      | Unit -> "unit" | Pi -> "PI" 
      | Alpha_QED -> "alpha" | E -> "e" | G_weak -> "g" | Vev -> "vev"
      | Sin2thw -> "sin2thw" | Sinthw -> "sinthw" | Costhw -> "costhw"
      | Q_lepton -> "qlep" | Q_up -> "qup" | Q_down -> "qdwn" 
      | Q_top -> "qtop"
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_NC_top -> "gnctop"
      | G_CC -> "gcc"
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww"
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc"
      | G_H3 -> "gh3" | G_H4 -> "gh4"
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f
  end

module O = Omega.Mixed23(Targets.Fortran)
    (Anomtop(SM_no_anomalous))
let _ = O.main ()

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
