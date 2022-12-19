(* modellib_NoH.ml --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       Christian Speckner <cnspeckn@googlemail.com>
       Marco Sekulla <marco.sekulla@kit.edu>
       Fabian Bach <fabian.bach@t-online.de> (only parts of this file)

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

(* \thocwmodulesection{Minimal Higgsless Model (Unitarity Gauge)} *)

module type NoH_flags =
  sig
    val triple_anom : bool
    val quartic_anom : bool
    val k_matrix : bool
    val ckm_present : bool
    val top_anom : bool
    val top_anom_4f : bool
  end

module NoH_k_matrix : NoH_flags =
  struct
    let triple_anom = false
    let quartic_anom = false
    let k_matrix = true
    let ckm_present = false
    let top_anom = false
    let top_anom_4f = false
  end

(* \thocwmodulesection{Minimal Higgsless Model including unitarization} *)

module NoH (Flags : NoH_flags) =
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
        "cms_width", Arg.Unit (fun () -> default_width := Complex_Mass),
        "use complex mass scheme";
        "running_width", Arg.Unit (fun () -> default_width := Running),
        "use running width"]
    let caveats () = []

    type f_aux_top = TTGG | TBWA | TBWZ | TTWW | BBWW |   (*i top auxiliary field "flavors" *)
                     QGUG | QBUB | QW | DL | DR |
                     QUQD1L | QUQD1R | QUQD8L | QUQD8R

    type matter_field = L of int | N of int | U of int | D of int
    type gauge_boson = Ga | Wp | Wm | Z | Gl
    type other = Phip | Phim | Phi0
                 | Aux_top of int*int*int*bool*f_aux_top    (*i lorentz*color*charge*top-side*flavor *)
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
      failwith "Modellib.NoH.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; N n; U n; D n ]

    let rec aux_top_flavors (f,l,co,ch) = List.append
      ( List.map other [ Aux_top(l,co,ch/2,true,f); Aux_top(l,co,ch/2,false,f) ] )
      ( if ch > 1 then List.append
          ( List.map other [ Aux_top(l,co,-ch/2,true,f); Aux_top(l,co,-ch/2,false,f) ] )
          ( aux_top_flavors (f,l,co,(ch-2)) )
        else [] )

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Gauge Bosons", List.map gauge_boson [Ga; Z; Wp; Wm; Gl];
        "Goldstone Bosons", List.map other [Phip; Phim; Phi0] ]

    let flavors () = List.append
      ( ThoList.flatmap snd (external_flavors ()) )
      ( ThoList.flatmap aux_top_flavors
         [ (TTGG,2,1,1); (TBWA,2,0,2); (TBWZ,2,0,2); (TTWW,2,0,1); (BBWW,2,0,1);
           (QGUG,1,1,1); (QBUB,1,0,1); (QW,1,0,3); (DL,0,0,3); (DR,0,0,3);
           (QUQD1L,0,0,3); (QUQD1R,0,0,3); (QUQD8L,0,1,3); (QUQD8R,0,1,3) ] )

    let spinor n =
      if n >= 0 then
        Spinor
      else
        ConjSpinor

    let lorentz_aux = function
      | 2 -> Tensor_1
      | 1 -> Vector
      | 0 -> Scalar
      | _ -> invalid_arg ("NoH.lorentz_aux: wrong value")

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
          begin match f with
          | Aux_top (l,_,_,_,_) -> lorentz_aux l
          | _ -> Scalar
          end

    let color = function 
      | M (U n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D n) -> Color.SUN (if n > 0 then 3 else -3)
      | G Gl -> Color.AdjSUN 3
      | O (Aux_top (_,co,_,_,_)) -> if co == 0 then Color.Singlet else Color.AdjSUN 3
      | _ -> Color.Singlet

    let nc () = 3

    let prop_spinor n =
      if n >= 0 then
        Prop_Spinor
      else
        Prop_ConjSpinor

    let prop_aux = function
      | 2 -> Aux_Tensor_1
      | 1 -> Aux_Vector
      | 0 -> Aux_Scalar
      | _ -> invalid_arg ("NoH.prop_aux: wrong value")

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
          | Aux_top (l,_,_,_,_) -> prop_aux l
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
          | Aux_top (l,co,ch,n,f) -> Aux_top (l,co,(-ch),(not n),f)
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
      |  n -> invalid_arg ("NoH.generation': " ^ string_of_int n)

    let generation f =
      if Flags.ckm_present then
        []
      else
        match f with
        | M (L n | N n | U n | D n) -> generation' n
        | G _ | O _ -> [0//1; 0//1; 0//1]

    let charge = function
      | M f ->
          begin match f with
          | L n -> if n > 0 then -1//1 else  1//1
          | N n -> 0//1
          | U n -> if n > 0 then  2//3 else -2//3
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
          | Phi0 ->  0//1
          | Phip ->  1//1
          | Phim -> -1//1
          | Aux_top (_,_,ch,_,_) -> ch//1
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
      | Unit | Half | Pi | Alpha_QED | Sin2thw
      | Sinthw | Costhw | E | G_weak | I_G_weak | Vev
      | Q_lepton | Q_up | Q_down | G_CC | G_CCQ of int*int
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down 
      | G_TVA_ttA | G_TVA_bbA 
      | G_VLR_ttZ | G_TVA_ttZ | G_TVA_bbZ 
      | G_VLR_btW | G_VLR_tbW
      | G_TLR_btW | G_TRL_tbW
      | G_TLR_btWZ | G_TRL_tbWZ
      | G_TLR_btWA | G_TRL_tbWA
      | G_TVA_ttWW | G_TVA_bbWW
      | G_TVA_ttG | G_TVA_ttGG
      | G_VLR_qGuG | G_VLR_qBuB
      | G_VLR_qBuB_u | G_VLR_qBuB_d | G_VLR_qBuB_e | G_VL_qBuB_n
      | G_VL_qW | G_VL_qW_u | G_VL_qW_d
      | G_SL_DttR | G_SR_DttR | G_SL_DttL | G_SLR_DbtR | G_SL_DbtL
      | C_quqd1R_bt | C_quqd1R_tb | C_quqd1L_bt | C_quqd1L_tb
      | C_quqd8R_bt | C_quqd8R_tb | C_quqd8L_bt | C_quqd8L_tb
      | I_Q_W | I_G_ZWW
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | I_G1_AWW | I_G1_ZWW
      | I_G1_plus_kappa_plus_G4_AWW
      | I_G1_plus_kappa_plus_G4_ZWW
      | I_G1_plus_kappa_minus_G4_AWW
      | I_G1_plus_kappa_minus_G4_ZWW
      | I_G1_minus_kappa_plus_G4_AWW
      | I_G1_minus_kappa_plus_G4_ZWW
      | I_G1_minus_kappa_minus_G4_AWW
      | I_G1_minus_kappa_minus_G4_ZWW
      | I_lambda_AWW | I_lambda_ZWW
      | G5_AWW | G5_ZWW
      | I_kappa5_AWW | I_kappa5_ZWW 
      | I_lambda5_AWW | I_lambda5_ZWW
      | Alpha_WWWW0 | Alpha_ZZWW1 | Alpha_WWWW2
      | Alpha_ZZWW0 | Alpha_ZZZZ
      | D_Alpha_ZZWW0_S | D_Alpha_ZZWW0_T | D_Alpha_ZZWW1_S
      | D_Alpha_ZZWW1_T | D_Alpha_ZZWW1_U | D_Alpha_WWWW0_S
      | D_Alpha_WWWW0_T | D_Alpha_WWWW0_U | D_Alpha_WWWW2_S
      | D_Alpha_WWWW2_T | D_Alpha_ZZZZ_S | D_Alpha_ZZZZ_T
      | Gs | I_Gs | G2
      | Mass of flavor | Width of flavor
      | K_Matrix_Coeff of int | K_Matrix_Pole of int
	  
(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | Q_lepton | Q_up | Q_down | G_NC_lepton | G_NC_neutrino 
      | G_NC_up | G_NC_down | G_CC | G_CCQ _ 
      | I_Q_W 
      | I_G_ZWW | I_G1_AWW | I_G1_ZWW | I_G_weak
      | Half | Unit 
      | I_G1_plus_kappa_plus_G4_AWW 
      | I_G1_plus_kappa_plus_G4_ZWW 
      | I_G1_minus_kappa_plus_G4_AWW 
      | I_G1_minus_kappa_plus_G4_ZWW 
      | I_G1_plus_kappa_minus_G4_AWW 
      | I_G1_plus_kappa_minus_G4_ZWW
      | I_G1_minus_kappa_minus_G4_AWW 
      | I_G1_minus_kappa_minus_G4_ZWW | I_kappa5_AWW 
      | I_kappa5_ZWW | G5_AWW | G5_ZWW 
      | I_lambda_AWW | I_lambda_ZWW | I_lambda5_AWW 
      | I_lambda5_ZWW | G_TVA_ttA | G_TVA_bbA 
      | G_VLR_ttZ | G_TVA_ttZ | G_TVA_bbZ 
      | G_VLR_btW | G_VLR_tbW | G_TLR_btW | G_TRL_tbW
      | G_TLR_btWA | G_TRL_tbWA | G_TLR_btWZ | G_TRL_tbWZ	
      | G_VLR_qBuB | G_VLR_qBuB_u | G_VLR_qBuB_d
      | G_VLR_qBuB_e | G_VL_qBuB_n | G_VL_qW | G_VL_qW_u | G_VL_qW_d
      | G_SL_DttR | G_SR_DttR  | G_SL_DttL | G_SLR_DbtR | G_SL_DbtL
      | G_TVA_ttWW | G_TVA_bbWW -> (0,1)
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW  
      |	Alpha_WWWW0 | Alpha_WWWW2 | Alpha_ZZWW0 
      | Alpha_ZZWW1 | Alpha_ZZZZ 
      | D_Alpha_WWWW0_S | D_Alpha_WWWW0_T | D_Alpha_WWWW0_U
      | D_Alpha_WWWW2_S | D_Alpha_WWWW2_T | D_Alpha_ZZWW0_S 
      | D_Alpha_ZZWW0_T | D_Alpha_ZZWW1_S | D_Alpha_ZZWW1_T
      | D_Alpha_ZZWW1_U | D_Alpha_ZZZZ_S | D_Alpha_ZZZZ_T -> (0,2)
      | Gs | I_Gs | G_TVA_ttG | G_TVA_ttGG | G_VLR_qGuG 
      | C_quqd1R_bt | C_quqd1R_tb | C_quqd1L_bt | C_quqd1L_tb
      | C_quqd8R_bt | C_quqd8R_tb | C_quqd8L_bt | C_quqd8L_tb -> (1,0)
      | G2  -> (2,0)
	(* These constants are not used, hence initialized to zero. *)
      | Sinthw | Sin2thw | Costhw | Pi 
      | Alpha_QED | G_weak | K_Matrix_Coeff _ 
      | K_Matrix_Pole _ | Mass _ | Width _ | Vev | E -> (0,0)

(* \begin{dubious}
     The current abstract syntax for parameter dependencies is admittedly
     tedious. Later, there will be a parser for a convenient concrete syntax
     as a part of a concrete syntax for models.  But as these examples show,
     it should include simple functions.
   \end{dubious} *)

(* \begin{subequations}
     \begin{align}
        \alpha_{\text{QED}} &= \frac{1}{137.0359895} \\
             \sin^2\theta_w &= 0.23124
     \end{align}
   \end{subequations} *)
    let input_parameters =
      [ Alpha_QED, 1. /. 137.0359895;
        Sin2thw, 0.23124;
        Mass (G Z), 91.187;
        Mass (M (N 1)), 0.0; Mass (M (L 1)), 0.51099907e-3;
        Mass (M (N 2)), 0.0; Mass (M (L 2)), 0.105658389;
        Mass (M (N 3)), 0.0; Mass (M (L 3)), 1.77705;
        Mass (M (U 1)), 5.0e-3; Mass (M (D 1)), 3.0e-3;
        Mass (M (U 2)), 1.2; Mass (M (D 2)), 0.1;
        Mass (M (U 3)), 174.0; Mass (M (D 3)), 4.2 ]

(* \begin{subequations}
     \begin{align}
                        e &= \sqrt{4\pi\alpha} \\
             \sin\theta_w &= \sqrt{\sin^2\theta_w} \\
             \cos\theta_w &= \sqrt{1-\sin^2\theta_w} \\
                        g &= \frac{e}{\sin\theta_w} \\
                      m_W &= \cos\theta_w m_Z \\
                        v &= \frac{2m_W}{g} \\
                  g_{CC}   =
       -\frac{g}{2\sqrt2} &= -\frac{e}{2\sqrt2\sin\theta_w} \\
       Q_{\text{lepton}}   =
      -q_{\text{lepton}}e &= e \\
           Q_{\text{up}}   =
          -q_{\text{up}}e &= -\frac{2}{3}e \\
         Q_{\text{down}}   =
        -q_{\text{down}}e &= \frac{1}{3}e \\
        \ii q_We           =
        \ii g_{\gamma WW} &= \ii e \\
              \ii g_{ZWW} &= \ii g \cos\theta_w \\
              \ii g_{WWW} &= \ii g
     \end{align}
   \end{subequations} *)



    let derived_parameters =
      [ Real E, Sqrt (Prod [Integer 4; Atom Pi; Atom Alpha_QED]);
        Real Sinthw, Sqrt (Atom Sin2thw);
        Real Costhw, Sqrt (Diff (Integer 1, Atom Sin2thw));
        Real G_weak, Quot (Atom E, Atom Sinthw);
        Real (Mass (G Wp)), Prod [Atom Costhw; Atom (Mass (G Z))];
        Real Vev, Quot (Prod [Integer 2; Atom (Mass (G Wp))], Atom G_weak);
        Real Q_lepton, Atom E;
        Real Q_up, Prod [Quot (Integer (-2), Integer 3); Atom E];
        Real Q_down, Prod [Quot (Integer 1, Integer 3); Atom E];
        Real G_CC, Neg (Quot (Atom G_weak, Prod [Integer 2; Sqrt (Integer 2)]));
        Complex I_Q_W, Prod [I; Atom E];
        Complex I_G_weak, Prod [I; Atom G_weak];
        Complex I_G_ZWW, Prod [I; Atom G_weak; Atom Costhw] ]
             
(* \begin{equation}
      - \frac{g}{2\cos\theta_w}
   \end{equation} *)
    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

(* \begin{subequations}
     \begin{align}
           - \frac{g}{2\cos\theta_w} g_V
        &= - \frac{g}{2\cos\theta_w} (T_3 - 2 q \sin^2\theta_w) \\
           - \frac{g}{2\cos\theta_w} g_A
        &= - \frac{g}{2\cos\theta_w} T_3
     \end{align}
   \end{subequations} *)
    let nc_coupling c t3 q =
      (Real_Array c,
       [Prod [g_over_2_costh; Diff (t3, Prod [Integer 2; q; Atom Sin2thw])];
        Prod [g_over_2_costh; t3]])

    let half = Quot (Integer 1, Integer 2)

    let derived_parameter_arrays =
      [ nc_coupling G_NC_neutrino half (Integer 0);
        nc_coupling G_NC_lepton (Neg half) (Integer (-1));
        nc_coupling G_NC_up half (Quot (Integer 2, Integer 3));
        nc_coupling G_NC_down (Neg half) (Quot (Integer (-1), Integer 3)) ]

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
    let mom ((m1, o, m2), fbf, c) = ((M m1, O o, M m2), fbf, c)

    let electromagnetic_currents n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]
        
    let color_currents n =
      List.map mgm
        [ ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs);
          ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs) ]

(* \begin{equation}
     \mathcal{L}_{\textrm{NC}} =
        - \frac{g}{2\cos\theta_W}
            \sum_i \bar\psi_i\fmslash{Z}(g_V^i-g_A^i\gamma_5)\psi_i
   \end{equation} *)

    let neutral_currents n =
      List.map mgm
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

(* \begin{equation}
     \mathcal{L}_{\textrm{CC}} =
        - \frac{g}{2\sqrt2} \sum_i \bar\psi_i
               (T^+\fmslash{W}^+ + T^-\fmslash{W}^-)(1-\gamma_5)\psi_i 
   \end{equation} *)

    let charged_currents' n = 
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let charged_currents'' n =
      List.map mgm 
        [ ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let charged_currents_triv = 
      ThoList.flatmap charged_currents' [1;2;3] @
      ThoList.flatmap charged_currents'' [1;2;3]

    let charged_currents_ckm = 
      let charged_currents_2 n1 n2 = 
        List.map mgm 
          [ ((D (-n1), Wm, U n2), FBF (1, Psibar, VL, Psi), G_CCQ (n2,n1));
            ((U (-n1), Wp, D n2), FBF (1, Psibar, VL, Psi), G_CCQ (n1,n2)) ] in
      ThoList.flatmap charged_currents' [1;2;3] @ 
      List.flatten (Product.list2 charged_currents_2 [1;2;3] [1;2;3])

      
(* \begin{equation}
     \mathcal{L}_{\textrm{TGC}} =
        - e \partial_\mu A_\nu W_+^\mu W_-^\nu + \ldots
        - e \cot\theta_w  \partial_\mu Z_\nu W_+^\mu W_-^\nu + \ldots
   \end{equation} *)

    let tgc ((g1, g2, g3), t, c) = ((G g1, G g2, G g3), t, c)

    let standard_triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs)]

(* \begin{multline}
     \mathcal{L}_{\textrm{TGC}}(g_1,\kappa)
        =   g_1 \mathcal{L}_T(V,W^+,W^-) \\
          + \frac{\kappa+g_1}{2} \Bigl(\mathcal{L}_T(W^-,V,W^+)
                                         - \mathcal{L}_T(W^+,V,W^-)\Bigr)\\
          + \frac{\kappa-g_1}{2} \Bigl(\mathcal{L}_L(W^-,V,W^+)
                                         - \mathcal{L}_T(W^+,V,W^-)\Bigr)
   \end{multline} *)

(* \begin{dubious}
   The whole thing in the LEP2 workshop notation:
   \begin{multline}
     \ii\mathcal{L}_{\textrm{TGC},V} / g_{WWV} = \\
            g_1^V V^\mu (W^-_{\mu\nu}W^{+,\nu}-W^+_{\mu\nu}W^{-,\nu})
          + \kappa_V  W^+_\mu W^-_\nu V^{\mu\nu}
          + \frac{\lambda_V}{m_W^2} V_{\mu\nu}
               W^-_{\rho\mu} W^{+,\hphantom{\nu}\rho}_{\hphantom{+,}\nu} \\
          + \ii g_5^V \epsilon_{\mu\nu\rho\sigma}
              \left(   (\partial^\rho W^{-,\mu}) W^{+,\nu}
                     -  W^{-,\mu}(\partial^\rho W^{+,\nu}) \right) V^\sigma \\
          + \ii g_4^V W^-_\mu W^+_\nu (\partial^\mu V^\nu + \partial^\nu V^\mu)
          - \frac{\tilde\kappa_V}{2}  W^-_\mu W^+_\nu \epsilon^{\mu\nu\rho\sigma}
              V_{\rho\sigma}
          - \frac{\tilde\lambda_V}{2m_W^2}
               W^-_{\rho\mu} W^{+,\mu}_{\hphantom{+,\mu}\nu} \epsilon^{\nu\rho\alpha\beta}
                V_{\alpha\beta}
   \end{multline}
   using the conventions of Itzykson and Zuber with $\epsilon^{0123} = +1$.
   \end{dubious} *)

(* \begin{dubious}
   This is equivalent to the notation of Hagiwara et al.~\cite{HPZH87}, if we
   remember that they have opposite signs for~$g_{WWV}$:
   \begin{multline}
     \mathcal{L}_{WWV} / (-g_{WWV})  = \\
       \ii g_1^V \left( W^\dagger_{\mu\nu} W^\mu 
                         - W^\dagger_\mu W^\mu_{\hphantom{\mu}\nu} \right) V^\nu
     + \ii \kappa_V  W^\dagger_\mu W_\nu V^{\mu\nu}
     + \ii \frac{\lambda_V}{m_W^2}
          W^\dagger_{\lambda\mu} W^\mu_{\hphantom{\mu}\nu} V^{\nu\lambda} \\
     - g_4^V  W^\dagger_\mu W_\nu
          \left(\partial^\mu V^\nu + \partial^\nu V^\mu \right)
     + g_5^V \epsilon^{\mu\nu\lambda\sigma}
           \left( W^\dagger_\mu \stackrel{\leftrightarrow}{\partial_\lambda}
                  W_\nu \right) V_\sigma\\
     + \ii \tilde\kappa_V  W^\dagger_\mu W_\nu \tilde{V}^{\mu\nu}
     + \ii\frac{\tilde\lambda_V}{m_W^2}
           W^\dagger_{\lambda\mu} W^\mu_{\hphantom{\mu}\nu} \tilde{V}^{\nu\lambda}
   \end{multline}
   Here $V^\mu$ stands for either the photon or the~$Z$ field, $W^\mu$ is the
   $W^-$ field, $W_{\mu\nu} = \partial_\mu W_\nu - \partial_\nu W_\mu$,
   $V_{\mu\nu} = \partial_\mu V_\nu - \partial_\nu V_\mu$, and
   $\tilde{V}_{\mu\nu} = \frac{1}{2} \epsilon_{\mu\nu\lambda\sigma}
   V^{\lambda\sigma}$.
   \end{dubious} *)

    let anomalous_triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Dim4_Vector_Vector_Vector_T (-1),
           I_G1_AWW);
          ((Z, Wm, Wp), Dim4_Vector_Vector_Vector_T (-1),
           I_G1_ZWW);
          ((Wm, Ga, Wp), Dim4_Vector_Vector_Vector_T 1,
           I_G1_plus_kappa_minus_G4_AWW);
          ((Wm, Z, Wp), Dim4_Vector_Vector_Vector_T 1,
           I_G1_plus_kappa_minus_G4_ZWW);
          ((Wp, Ga, Wm), Dim4_Vector_Vector_Vector_T (-1),
           I_G1_plus_kappa_plus_G4_AWW);
          ((Wp, Z, Wm), Dim4_Vector_Vector_Vector_T (-1),
           I_G1_plus_kappa_plus_G4_ZWW);
          ((Wm, Ga, Wp), Dim4_Vector_Vector_Vector_L (-1),
           I_G1_minus_kappa_plus_G4_AWW);
          ((Wm, Z, Wp), Dim4_Vector_Vector_Vector_L (-1),
           I_G1_minus_kappa_plus_G4_ZWW);
          ((Wp, Ga, Wm), Dim4_Vector_Vector_Vector_L 1,
           I_G1_minus_kappa_minus_G4_AWW);
          ((Wp, Z, Wm), Dim4_Vector_Vector_Vector_L 1,
           I_G1_minus_kappa_minus_G4_ZWW);
          ((Ga, Wm, Wp), Dim4_Vector_Vector_Vector_L5 (-1),
           I_kappa5_AWW);
          ((Z, Wm, Wp), Dim4_Vector_Vector_Vector_L5 (-1),
           I_kappa5_ZWW);
          ((Ga, Wm, Wp), Dim4_Vector_Vector_Vector_T5 (-1),
           G5_AWW);
          ((Z, Wm, Wp), Dim4_Vector_Vector_Vector_T5 (-1),
           G5_ZWW);
          ((Ga, Wp, Wm), Dim6_Gauge_Gauge_Gauge (-1),
           I_lambda_AWW);
          ((Z, Wp, Wm), Dim6_Gauge_Gauge_Gauge (-1),
           I_lambda_ZWW);
          ((Ga, Wp, Wm), Dim6_Gauge_Gauge_Gauge_5 (-1),
           I_lambda5_AWW);
          ((Z, Wp, Wm), Dim6_Gauge_Gauge_Gauge_5 (-1),
           I_lambda5_ZWW) ]

    let triple_gauge =
      if Flags.triple_anom then
        anomalous_triple_gauge
      else
        standard_triple_gauge

(* \begin{equation}
     \mathcal{L}_{\textrm{QGC}} =
        - g^2 W_{+,\mu} W_{-,\nu} W_+^\mu W_-^\nu + \ldots
   \end{equation} *)

(* Actually, quartic gauge couplings are a little bit more straightforward
   using auxiliary fields.  Here we have to impose the antisymmetry manually:
   \begin{subequations}
   \begin{multline}
     (W^{+,\mu}_1 W^{-,\nu}_2 - W^{+,\nu}_1 W^{-,\mu}_2)
     (W^+_{3,\mu} W^-_{4,\nu} - W^+_{3,\nu} W^-_{4,\mu}) \\
        = 2(W^+_1W^+_3)(W^-_2W^-_4) - 2(W^+_1W^-_4)(W^-_2W^+_3)
   \end{multline}
   also ($V$ can be $A$ or $Z$)
   \begin{multline}
     (W^{+,\mu}_1 V^\nu_2 - W^{+,\nu}_1 V^\mu_2)
     (W^-_{3,\mu} V_{4,\nu} - W^-_{3,\nu} V_{4,\mu}) \\
        = 2(W^+_1W^-_3)(V_2V_4) - 2(W^+_1V_4)(V_2W^-_3)
   \end{multline}
   \end{subequations} *)

(* \begin{subequations}
   \begin{multline}
      W^{+,\mu} W^{-,\nu} W^+_\mu W^-_\nu
   \end{multline}
   \end{subequations} *)

    let qgc ((g1, g2, g3, g4), t, c) = ((G g1, G g2, G g3, G g4), t, c)

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let standard_quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          (Gl, Gl, Gl, Gl), gauge4, G2 ]

(* \begin{subequations}
   \begin{align}
     \mathcal{L}_4
       &= \alpha_4 \left(   \frac{g^4}{2}\left(   (W^+_\mu W^{-,\mu})^2
                                                + W^+_\mu W^{+,\mu} W^-_\mu W^{-,\mu}
                                               \right)\right.\notag \\
       &\qquad\qquad\qquad \left.
                          + \frac{g^4}{\cos^2\theta_w} W^+_\mu Z^\mu W^-_\nu Z^\nu
                          + \frac{g^4}{4\cos^4\theta_w} (Z_\mu Z^\mu)^2 \right) \\
     \mathcal{L}_5
       &= \alpha_5 \left(   g^4 (W^+_\mu W^{-,\mu})^2
                          + \frac{g^4}{\cos^2\theta_w}  W^+_\mu W^{-,\mu} Z_\nu Z^\nu
                          + \frac{g^4}{4\cos^4\theta_w} (Z_\mu Z^\mu)^2 \right)
   \end{align}
   \end{subequations}
   or
   \begin{multline}
     \mathcal{L}_4 + \mathcal{L}_5
       =   (\alpha_4+2\alpha_5) g^4 \frac{1}{2} (W^+_\mu W^{-,\mu})^2 \\
         + 2\alpha_4 g^4 \frac{1}{4} W^+_\mu W^{+,\mu} W^-_\mu W^{-,\mu}
         + \alpha_4 \frac{g^4}{\cos^2\theta_w} W^+_\mu Z^\mu W^-_\nu Z^\nu \\
         + 2\alpha_5 \frac{g^4}{\cos^2\theta_w} \frac{1}{2} W^+_\mu W^{-,\mu} Z_\nu Z^\nu
         + (2\alpha_4 + 2\alpha_5) \frac{g^4}{\cos^4\theta_w} \frac{1}{8} (Z_\mu Z^\mu)^2
   \end{multline}
   and therefore
   \begin{subequations}
   \begin{align}
     \alpha_{(WW)_0} &= (\alpha_4+2\alpha_5) g^4 \\
     \alpha_{(WW)_2} &= 2\alpha_4 g^4 \\
     \alpha_{(WZ)_0} &= 2\alpha_5 \frac{g^4}{\cos^2\theta_w} \\
     \alpha_{(WZ)_1} &= \alpha_4 \frac{g^4}{\cos^2\theta_w} \\
     \alpha_{ZZ} &= (2\alpha_4 + 2\alpha_5) \frac{g^4}{\cos^4\theta_w}
   \end{align}
   \end{subequations} *)

    let anomalous_quartic_gauge =
      if Flags.quartic_anom then
        List.map qgc
          [ ((Wm, Wm, Wp, Wp),
             Vector4 [(1, C_13_42); (1, C_14_23)], Alpha_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Vector4 [1, C_12_34], Alpha_WWWW2);
            ((Wm, Wp, Z, Z),
             Vector4 [1, C_12_34], Alpha_ZZWW0);
            ((Wm, Wp, Z, Z),
             Vector4 [(1, C_13_42); (1, C_14_23)], Alpha_ZZWW1);
            ((Z, Z, Z, Z),
             Vector4 [(1, C_12_34); (1, C_13_42); (1, C_14_23)], Alpha_ZZZZ) ]
      else
        []

(* In any diagonal channel~$\chi$, the scattering amplitude~$a_\chi(s)$ is
   unitary iff\footnote{%
     Trivial proof:
     \begin{equation}
       -1 = \textrm{Im}\left(\frac{1}{a_\chi(s)}\right)
          = \frac{\textrm{Im}(a_\chi^*(s))}{ |a_\chi(s)|^2 }
          = - \frac{\textrm{Im}(a_\chi(s))}{ |a_\chi(s)|^2 }
     \end{equation}
     i.\,e.~$\textrm{Im}(a_\chi(s)) = |a_\chi(s)|^2$.}
   \begin{equation}
     \textrm{Im}\left(\frac{1}{a_\chi(s)}\right) = -1
   \end{equation}
   For a real perturbative scattering amplitude~$r_\chi(s)$ this can be
   enforced easily--and arbitrarily--by
   \begin{equation}
     \frac{1}{a_\chi(s)} = \frac{1}{r_\chi(s)} - \mathrm{i}
   \end{equation} 

*)


    let k_matrix_quartic_gauge =
      if Flags.k_matrix then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_14_23)]), D_Alpha_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42)]), D_Alpha_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_jr (0,
                   [(1, C_14_23)]), D_Alpha_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42)]), D_Alpha_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42); (1, C_14_23)]), D_Alpha_WWWW2_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42); (1, C_14_23)]), D_Alpha_ZZWW0_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42)]), D_Alpha_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_14_23)]), D_Alpha_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_jr (1,
                   [(1, C_12_34)]), D_Alpha_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_jr (1,
                   [(1, C_13_42)]), D_Alpha_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_jr (1,
                   [(1, C_14_23)]), D_Alpha_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_jr (2,
                   [(1, C_12_34)]), D_Alpha_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_jr (2,
                   [(1, C_13_42)]), D_Alpha_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_jr (2,
                   [(1, C_14_23)]), D_Alpha_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42); (1, C_14_23)]), D_Alpha_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_jr (3,
                   [(1, C_14_23)]), D_Alpha_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_jr (3,
                   [(1, C_13_42); (1, C_12_34)]), D_Alpha_ZZZZ_T)]
      else
        []



(*i Thorsten's original implementation of the K matrix, which we keep since
   it still might be usefull for the future. 


    let k_matrix_quartic_gauge =
      if Flags.k_matrix then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_tho (0, [K_Matrix_Coeff 0, 
                         K_Matrix_Pole 0]), Alpha_WWWW0);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_tho (0, [K_Matrix_Coeff 2, 
                         K_Matrix_Pole 2]), Alpha_WWWW2);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_tho (0, [(K_Matrix_Coeff 0, 
                         K_Matrix_Pole 0); (K_Matrix_Coeff 2, 
                         K_Matrix_Pole 2)]), Alpha_ZZWW0);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_tho (0, [K_Matrix_Coeff 1, 
                         K_Matrix_Pole 1]), Alpha_ZZWW1);
            ((Z, Z, Z, Z), Vector4_K_Matrix_tho (0, [K_Matrix_Coeff 0, 
                         K_Matrix_Pole 0]), Alpha_ZZZZ) ]
      else
        []

i*)

    let quartic_gauge =
      standard_quartic_gauge @ anomalous_quartic_gauge @ k_matrix_quartic_gauge

(* WK's couplings (apparently, he still intends to divide by
   $\Lambda^2_{\text{EWSB}}=16\pi^2v_{\mathrm{F}}^2$):
   \begin{subequations}
   \begin{align}
     \mathcal{L}^{\tau}_4 &=
      \left\lbrack (\partial_{\mu}H)(\partial^{\mu}H)
                     + \frac{g^2v_{\mathrm{F}}^2}{4} V_{\mu} V^{\mu} \right\rbrack^2 \\
     \mathcal{L}^{\tau}_5 &=
      \left\lbrack (\partial_{\mu}H)(\partial_{\nu}H)
                     + \frac{g^2v_{\mathrm{F}}^2}{4} V_{\mu} V_{\nu} \right\rbrack^2
   \end{align}
   \end{subequations}
   with
   \begin{equation}
      V_{\mu} V_{\nu} =
        \frac{1}{2} \left( W^+_{\mu} W^-_{\nu} + W^+_{\nu} W^-_{\mu} \right)
         + \frac{1}{2\cos^2\theta_{w}} Z_{\mu} Z_{\nu}
   \end{equation}
   (note the symmetrization!), i.\,e.
   \begin{subequations}
   \begin{align}
     \mathcal{L}_4 &= \alpha_4 \frac{g^4v_{\mathrm{F}}^4}{16} (V_{\mu} V_{\nu})^2 \\
     \mathcal{L}_5 &= \alpha_5 \frac{g^4v_{\mathrm{F}}^4}{16} (V_{\mu} V^{\mu})^2
   \end{align}
   \end{subequations} *)

    let goldstone_vertices =
      [ ((O Phi0, G Wm, G Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phip, G Ga, G Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phip, G Z, G Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phim, G Wp, G Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phim, G Wp, G Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

(* Anomalous trilinear interactions $f_i f_j V$ :
   \begin{equation}
     \Delta\mathcal{L}_{tt\gamma} =
        - e \frac{\upsilon}{\Lambda^2}
            \bar{t} i\sigma^{\mu\nu} k_\nu (d_V(k^2) + i d_A(k^2) \gamma_5) t A_\mu
   \end{equation} *)

    let anomalous_ttA =
      if Flags.top_anom then
        [ ((M (U (-3)), G Ga, M (U 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_ttA) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{bb\gamma} =
        - e \frac{\upsilon}{\Lambda^2}
            \bar{b} i\sigma^{\mu\nu} k_\nu (d_V(k^2) + i d_A(k^2) \gamma_5) b A_\mu
   \end{equation} *)

    let anomalous_bbA =
      if Flags.top_anom then
        [ ((M (D (-3)), G Ga, M (D 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_bbA) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{ttg} =
        - g_s \frac{\upsilon}{\Lambda^2}
            \bar{t}\lambda^a i\sigma^{\mu\nu}k_\nu
                (d_V(k^2)+id_A(k^2)\gamma_5)tG^a_\mu
   \end{equation} *)

    let anomalous_ttG =
      if Flags.top_anom then
        [ ((M (U (-3)), G Gl, M (U 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_ttG) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{ttZ} =
        - \frac{g}{2 c_W} \frac{\upsilon^2}{\Lambda^2}\left\lbrack
              \bar{t} \fmslash{Z} (X_L(k^2) P_L + X_R(k^2) P_R) t
            + \bar{t}\frac{i\sigma^{\mu\nu}k_\nu}{m_Z}
                  (d_V(k^2)+id_A(k^2)\gamma_5)tZ_\mu\right\rbrack
   \end{equation} *)

    let anomalous_ttZ =
      if Flags.top_anom then
        [ ((M (U (-3)), G Z, M (U 3)), FBF (1, Psibar, VLRM, Psi), G_VLR_ttZ);
          ((M (U (-3)), G Z, M (U 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_ttZ) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{bbZ} =
        - \frac{g}{2 c_W} \frac{\upsilon^2}{\Lambda^2}
              \bar{b}\frac{i\sigma^{\mu\nu}k_\nu}{m_Z}
                  (d_V(k^2)+id_A(k^2)\gamma_5)bZ_\mu
   \end{equation} *)

    let anomalous_bbZ =
      if Flags.top_anom then
        [ ((M (D (-3)), G Z, M (D 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_bbZ) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{tbW} =
        - \frac{g}{\sqrt{2}} \frac{\upsilon^2}{\Lambda^2}\left\lbrack
            \bar{b}\fmslash{W}^-(V_L(k^2) P_L+V_R(k^2) P_R) t
          + \bar{b}\frac{i\sigma^{\mu\nu}k_\nu}{m_W}
                (g_L(k^2)P_L+g_R(k^2)P_R)tW^-_\mu\right\rbrack
        + \textnormal{H.c.}
   \end{equation} *)

    let anomalous_tbW =
      if Flags.top_anom then
        [ ((M (D (-3)), G Wm, M (U 3)), FBF (1, Psibar, VLRM, Psi), G_VLR_btW);
          ((M (U (-3)), G Wp, M (D 3)), FBF (1, Psibar, VLRM, Psi), G_VLR_tbW);
          ((M (D (-3)), G Wm, M (U 3)), FBF (1, Psibar, TLRM, Psi), G_TLR_btW);
          ((M (U (-3)), G Wp, M (D 3)), FBF (1, Psibar, TRLM, Psi), G_TRL_tbW) ]
      else
        []

(* quartic fermion-gauge interactions $f_i f_j V_1 V_2$ emerging from gauge-invariant
effective operators:
   \begin{equation}
     \Delta\mathcal{L}_{ttgg} =
        - \frac{g_s^2}{2} f_{abc} \frac{\upsilon}{\Lambda^2}
            \bar{t} \lambda^a \sigma^{\mu\nu}
                (d_V(k^2)+id_A(k^2)\gamma_5)t G^b_\mu G^c_\nu
   \end{equation} *)

    let anomalous_ttGG =
      if Flags.top_anom then
        [ ((M (U (-3)), O (Aux_top (2,1,0,true,TTGG)), M (U 3)), FBF (1, Psibar, TVA, Psi), G_TVA_ttGG);
          ((O (Aux_top (2,1,0,false,TTGG)), G Gl, G Gl), Aux_Gauge_Gauge 1, I_Gs) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{tbWA} =
        - i\sin\theta_w \frac{g^2}{2\sqrt{2}} \frac{\upsilon^2}{\Lambda^2}\left\lbrack
            \bar{b}\frac{\sigma^{\mu\nu}}{m_W}
                (g_L(k^2)P_L+g_R(k^2)P_R)t A_\mu W^-_\nu \right\rbrack
        + \textnormal{H.c.}
   \end{equation} *)

    let anomalous_tbWA =
      if Flags.top_anom then
        [ ((M (D (-3)), O (Aux_top (2,0,-1,true,TBWA)), M (U 3)), FBF (1, Psibar, TLR, Psi), G_TLR_btWA);
          ((O (Aux_top (2,0,1,false,TBWA)), G Ga, G Wm), Aux_Gauge_Gauge 1, I_G_weak);
          ((M (U (-3)), O (Aux_top (2,0,1,true,TBWA)), M (D 3)), FBF (1, Psibar, TRL, Psi), G_TRL_tbWA);
          ((O (Aux_top (2,0,-1,false,TBWA)), G Wp, G Ga), Aux_Gauge_Gauge 1, I_G_weak) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{tbWZ} =
        - i\cos\theta_w \frac{g^2}{2\sqrt{2}} \frac{\upsilon^2}{\Lambda^2}\left\lbrack
            \bar{b}\frac{\sigma^{\mu\nu}}{m_W}
                (g_L(k^2)P_L+g_R(k^2)P_R)t Z_\mu W^-_\nu \right\rbrack
        + \textnormal{H.c.}
   \end{equation} *)

    let anomalous_tbWZ =
      if Flags.top_anom then
        [ ((M (D (-3)), O (Aux_top (2,0,-1,true,TBWZ)), M (U 3)), FBF (1, Psibar, TLR, Psi), G_TLR_btWZ);
          ((O (Aux_top (2,0,1,false,TBWZ)), G Z, G Wm), Aux_Gauge_Gauge 1, I_G_weak);
          ((M (U (-3)), O (Aux_top (2,0,1,true,TBWZ)), M (D 3)), FBF (1, Psibar, TRL, Psi), G_TRL_tbWZ);
          ((O (Aux_top (2,0,-1,false,TBWZ)), G Wp, G Z), Aux_Gauge_Gauge 1, I_G_weak) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{ttWW} =
        - i \frac{g^2}{2} \frac{\upsilon^2}{\Lambda^2}
            \bar{t} \frac{\sigma^{\mu\nu}}{m_W}
                (d_V(k^2)+id_A(k^2)\gamma_5)t W^-_\mu W^+_\nu
   \end{equation} *)

    let anomalous_ttWW =
      if Flags.top_anom then
        [ ((M (U (-3)), O (Aux_top (2,0,0,true,TTWW)), M (U 3)), FBF (1, Psibar, TVA, Psi), G_TVA_ttWW);
          ((O (Aux_top (2,0,0,false,TTWW)), G Wm, G Wp), Aux_Gauge_Gauge 1, I_G_weak) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{bbWW} =
        - i \frac{g^2}{2} \frac{\upsilon^2}{\Lambda^2}
            \bar{b} \frac{\sigma^{\mu\nu}}{m_W}
                (d_V(k^2)+id_A(k^2)\gamma_5)b W^-_\mu W^+_\nu
   \end{equation} *)

    let anomalous_bbWW =
      if Flags.top_anom then
        [ ((M (D (-3)), O (Aux_top (2,0,0,true,BBWW)), M (D 3)), FBF (1, Psibar, TVA, Psi), G_TVA_bbWW);
          ((O (Aux_top (2,0,0,false,BBWW)), G Wm, G Wp), Aux_Gauge_Gauge 1, I_G_weak) ]
      else
        []

(* 4-fermion contact terms emerging from operator rewriting: *)

    let anomalous_top_qGuG_tt =
      [ ((M (U (-3)), O (Aux_top (1,1,0,true,QGUG)), M (U 3)), FBF (1, Psibar, VLR, Psi), G_VLR_qGuG) ]

    let anomalous_top_qGuG_ff n =
      List.map mom
        [ ((U (-n), Aux_top (1,1,0,false,QGUG), U n), FBF (1, Psibar, V, Psi), Unit);
          ((D (-n), Aux_top (1,1,0,false,QGUG), D n), FBF (1, Psibar, V, Psi), Unit) ]

    let anomalous_top_qGuG =
      if Flags.top_anom_4f then
        anomalous_top_qGuG_tt @ ThoList.flatmap anomalous_top_qGuG_ff [1;2;3]
      else
        []

    let anomalous_top_qBuB_tt =
      [ ((M (U (-3)), O (Aux_top (1,0,0,true,QBUB)), M (U 3)), FBF (1, Psibar, VLR, Psi), G_VLR_qBuB) ]

    let anomalous_top_qBuB_ff n =
      List.map mom
        [ ((U (-n), Aux_top (1,0,0,false,QBUB), U n), FBF (1, Psibar, VLR, Psi), G_VLR_qBuB_u);
          ((D (-n), Aux_top (1,0,0,false,QBUB), D n), FBF (1, Psibar, VLR, Psi), G_VLR_qBuB_d);
          ((L (-n), Aux_top (1,0,0,false,QBUB), L n), FBF (1, Psibar, VLR, Psi), G_VLR_qBuB_e);
          ((N (-n), Aux_top (1,0,0,false,QBUB), N n), FBF (1, Psibar, VL, Psi), G_VL_qBuB_n) ]

    let anomalous_top_qBuB =
      if Flags.top_anom_4f then
        anomalous_top_qBuB_tt @ ThoList.flatmap anomalous_top_qBuB_ff [1;2;3]
      else
        []

    let anomalous_top_qW_tq =
      [ ((M (U (-3)), O (Aux_top (1,0,0,true,QW)), M (U 3)), FBF (1, Psibar, VL, Psi), G_VL_qW);
        ((M (D (-3)), O (Aux_top (1,0,-1,true,QW)), M (U 3)), FBF (1, Psibar, VL, Psi), G_VL_qW);
        ((M (U (-3)), O (Aux_top (1,0,1,true,QW)), M (D 3)), FBF (1, Psibar, VL, Psi), G_VL_qW) ]

    let anomalous_top_qW_ff n =
      List.map mom
        [ ((U (-n), Aux_top (1,0,0,false,QW), U n), FBF (1, Psibar, VL, Psi), G_VL_qW_u);
          ((D (-n), Aux_top (1,0,0,false,QW), D n), FBF (1, Psibar, VL, Psi), G_VL_qW_d);
          ((N (-n), Aux_top (1,0,0,false,QW), N n), FBF (1, Psibar, VL, Psi), G_VL_qW_u);
          ((L (-n), Aux_top (1,0,0,false,QW), L n), FBF (1, Psibar, VL, Psi), G_VL_qW_d);
          ((D (-n), Aux_top (1,0,-1,false,QW), U n), FBF (1, Psibar, VL, Psi), Half);
          ((U (-n), Aux_top (1,0,1,false,QW), D n), FBF (1, Psibar, VL, Psi), Half);
          ((L (-n), Aux_top (1,0,-1,false,QW), N n), FBF (1, Psibar, VL, Psi), Half);
          ((N (-n), Aux_top (1,0,1,false,QW), L n), FBF (1, Psibar, VL, Psi), Half) ]

    let anomalous_top_qW =
      if Flags.top_anom_4f then
        anomalous_top_qW_tq @ ThoList.flatmap anomalous_top_qW_ff [1;2;3]
      else
        []

    let anomalous_top_DuDd =
      if Flags.top_anom_4f then
        [ ((M (U (-3)), O (Aux_top (0,0,0,true,DR)), M (U 3)), FBF (1, Psibar, SR, Psi), Half);
          ((M (U (-3)), O (Aux_top (0,0,0,false,DR)), M (U 3)), FBF (1, Psibar, SL, Psi), G_SL_DttR);
          ((M (D (-3)), O (Aux_top (0,0,0,false,DR)), M (D 3)), FBF (1, Psibar, SR, Psi), G_SR_DttR);
          ((M (U (-3)), O (Aux_top (0,0,0,true,DL)), M (U 3)), FBF (1, Psibar, SL, Psi), Half);
          ((M (D (-3)), O (Aux_top (0,0,0,false,DL)), M (D 3)), FBF (1, Psibar, SL, Psi), G_SL_DttL);
          ((M (D (-3)), O (Aux_top (0,0,-1,true,DR)), M (U 3)), FBF (1, Psibar, SR, Psi), Half);
          ((M (U (-3)), O (Aux_top (0,0,1,false,DR)), M (D 3)), FBF (1, Psibar, SLR, Psi), G_SLR_DbtR);
          ((M (D (-3)), O (Aux_top (0,0,-1,true,DL)), M (U 3)), FBF (1, Psibar, SL, Psi), Half);
          ((M (U (-3)), O (Aux_top (0,0,1,false,DL)), M (D 3)), FBF (1, Psibar, SL, Psi), G_SL_DbtL) ]
      else
        []

    let anomalous_top_quqd1_tq =
      [ ((M (D (-3)), O (Aux_top (0,0,-1,true,QUQD1R)), M (U 3)), FBF (1, Psibar, SR, Psi), C_quqd1R_bt);
        ((M (U (-3)), O (Aux_top (0,0, 1,true,QUQD1R)), M (D 3)), FBF (1, Psibar, SL, Psi), C_quqd1R_tb);
        ((M (D (-3)), O (Aux_top (0,0,-1,true,QUQD1L)), M (U 3)), FBF (1, Psibar, SL, Psi), C_quqd1L_bt);
        ((M (U (-3)), O (Aux_top (0,0, 1,true,QUQD1L)), M (D 3)), FBF (1, Psibar, SR, Psi), C_quqd1L_tb) ]

    let anomalous_top_quqd1_ff n =
      List.map mom
        [ ((U (-n), Aux_top (0,0, 1,false,QUQD1R), D n), FBF (1, Psibar, SR, Psi), Half);
          ((D (-n), Aux_top (0,0,-1,false,QUQD1R), U n), FBF (1, Psibar, SL, Psi), Half);
          ((U (-n), Aux_top (0,0, 1,false,QUQD1L), D n), FBF (1, Psibar, SL, Psi), Half);
          ((D (-n), Aux_top (0,0,-1,false,QUQD1L), U n), FBF (1, Psibar, SR, Psi), Half) ]

    let anomalous_top_quqd1 =
      if Flags.top_anom_4f then
        anomalous_top_quqd1_tq @ ThoList.flatmap anomalous_top_quqd1_ff [1;2;3]
      else
        []

    let anomalous_top_quqd8_tq =
      [ ((M (D (-3)), O (Aux_top (0,1,-1,true,QUQD8R)), M (U 3)), FBF (1, Psibar, SR, Psi), C_quqd8R_bt);
        ((M (U (-3)), O (Aux_top (0,1, 1,true,QUQD8R)), M (D 3)), FBF (1, Psibar, SL, Psi), C_quqd8R_tb);
        ((M (D (-3)), O (Aux_top (0,1,-1,true,QUQD8L)), M (U 3)), FBF (1, Psibar, SL, Psi), C_quqd8L_bt);
        ((M (U (-3)), O (Aux_top (0,1, 1,true,QUQD8L)), M (D 3)), FBF (1, Psibar, SR, Psi), C_quqd8L_tb) ]

    let anomalous_top_quqd8_ff n =
      List.map mom
        [ ((U (-n), Aux_top (0,1, 1,false,QUQD8R), D n), FBF (1, Psibar, SR, Psi), Half);
          ((D (-n), Aux_top (0,1,-1,false,QUQD8R), U n), FBF (1, Psibar, SL, Psi), Half);
          ((U (-n), Aux_top (0,1, 1,false,QUQD8L), D n), FBF (1, Psibar, SL, Psi), Half);
          ((D (-n), Aux_top (0,1,-1,false,QUQD8L), U n), FBF (1, Psibar, SR, Psi), Half) ]

    let anomalous_top_quqd8 =
      if Flags.top_anom_4f then
        anomalous_top_quqd8_tq @ ThoList.flatmap anomalous_top_quqd8_ff [1;2;3]
      else
        []

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @
       ThoList.flatmap color_currents [1;2;3] @
       ThoList.flatmap neutral_currents [1;2;3] @
       (if Flags.ckm_present then
         charged_currents_ckm
       else
         charged_currents_triv) @
       triple_gauge @
       goldstone_vertices @
       anomalous_ttA @ anomalous_bbA @
       anomalous_ttZ @ anomalous_bbZ @
       anomalous_tbW @ anomalous_tbWA @ anomalous_tbWZ @
       anomalous_ttWW @ anomalous_bbWW @
       anomalous_ttG @ anomalous_ttGG @
       anomalous_top_qGuG @ anomalous_top_qBuB @
       anomalous_top_qW @ anomalous_top_DuDd @
       anomalous_top_quqd1 @ anomalous_top_quqd8)

    let vertices4 =
      quartic_gauge 

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
      | "Aux_t_ttGG0" -> O (Aux_top (2,1, 0,true,TTGG)) | "Aux_ttGG0" -> O (Aux_top (2,1, 0,false,TTGG))
      | "Aux_t_tbWA+" -> O (Aux_top (2,0, 1,true,TBWA)) | "Aux_tbWA+" -> O (Aux_top (2,0, 1,false,TBWA))
      | "Aux_t_tbWA-" -> O (Aux_top (2,0,-1,true,TBWA)) | "Aux_tbWA-" -> O (Aux_top (2,0,-1,false,TBWA))
      | "Aux_t_tbWZ+" -> O (Aux_top (2,0, 1,true,TBWZ)) | "Aux_tbWZ+" -> O (Aux_top (2,0, 1,false,TBWZ))
      | "Aux_t_tbWZ-" -> O (Aux_top (2,0,-1,true,TBWZ)) | "Aux_tbWZ-" -> O (Aux_top (2,0,-1,false,TBWZ))
      | "Aux_t_ttWW0" -> O (Aux_top (2,0, 0,true,TTWW)) | "Aux_ttWW0" -> O (Aux_top (2,0, 0,false,TTWW))
      | "Aux_t_bbWW0" -> O (Aux_top (2,0, 0,true,BBWW)) | "Aux_bbWW0" -> O (Aux_top (2,0, 0,false,BBWW))
      | "Aux_t_qGuG0" -> O (Aux_top (1,1, 0,true,QGUG)) | "Aux_qGuG0" -> O (Aux_top (1,1, 0,false,QGUG))
      | "Aux_t_qBuB0" -> O (Aux_top (1,0, 0,true,QBUB)) | "Aux_qBuB0" -> O (Aux_top (1,0, 0,false,QBUB))
      | "Aux_t_qW0"   -> O (Aux_top (1,0, 0,true,QW))   | "Aux_qW0"   -> O (Aux_top (1,0, 0,false,QW))
      | "Aux_t_qW+"   -> O (Aux_top (1,0, 1,true,QW))   | "Aux_qW+"   -> O (Aux_top (1,0, 1,false,QW))
      | "Aux_t_qW-"   -> O (Aux_top (1,0,-1,true,QW))   | "Aux_qW-"   -> O (Aux_top (1,0,-1,false,QW))
      | "Aux_t_dL0"   -> O (Aux_top (0,0, 0,true,DL))   | "Aux_dL0"   -> O (Aux_top (0,0, 0,false,DL))
      | "Aux_t_dL+"   -> O (Aux_top (0,0, 1,true,DL))   | "Aux_dL+"   -> O (Aux_top (0,0, 1,false,DL))
      | "Aux_t_dL-"   -> O (Aux_top (0,0,-1,true,DL))   | "Aux_dL-"   -> O (Aux_top (0,0,-1,false,DL))
      | "Aux_t_dR0"   -> O (Aux_top (0,0, 0,true,DR))   | "Aux_dR0"   -> O (Aux_top (0,0, 0,false,DR))
      | "Aux_t_dR+"   -> O (Aux_top (0,0, 1,true,DR))   | "Aux_dR+"   -> O (Aux_top (0,0, 1,false,DR))
      | "Aux_t_dR-"   -> O (Aux_top (0,0,-1,true,DR))   | "Aux_dR-"   -> O (Aux_top (0,0,-1,false,DR))
      | "Aux_t_quqd1L+" -> O (Aux_top (0,0, 1,true,QUQD1L)) | "Aux_quqd1L+" -> O (Aux_top (0,0, 1,false,QUQD1L))
      | "Aux_t_quqd1L-" -> O (Aux_top (0,0,-1,true,QUQD1L)) | "Aux_quqd1L-" -> O (Aux_top (0,0,-1,false,QUQD1L))
      | "Aux_t_quqd1R+" -> O (Aux_top (0,0, 1,true,QUQD1R)) | "Aux_quqd1R+" -> O (Aux_top (0,0, 1,false,QUQD1R))
      | "Aux_t_quqd1R-" -> O (Aux_top (0,0,-1,true,QUQD1R)) | "Aux_quqd1R-" -> O (Aux_top (0,0,-1,false,QUQD1R))
      | "Aux_t_quqd8L+" -> O (Aux_top (0,1, 1,true,QUQD8L)) | "Aux_quqd8L+" -> O (Aux_top (0,1, 1,false,QUQD8L))
      | "Aux_t_quqd8L-" -> O (Aux_top (0,1,-1,true,QUQD8L)) | "Aux_quqd8L-" -> O (Aux_top (0,1,-1,false,QUQD8L))
      | "Aux_t_quqd8R+" -> O (Aux_top (0,1, 1,true,QUQD8R)) | "Aux_quqd8R+" -> O (Aux_top (0,1, 1,false,QUQD8R))
      | "Aux_t_quqd8R-" -> O (Aux_top (0,1,-1,true,QUQD8R)) | "Aux_quqd8R-" -> O (Aux_top (0,1,-1,false,QUQD8R))
      | _ -> invalid_arg "Modellib.NoH.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib.NoH.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib.NoH.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib.NoH.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib.NoH.flavor_to_string: invalid down type quark"
          end
      | G f ->
          begin match f with
          | Gl -> "gl"
          | Ga -> "A" | Z -> "Z"
          | Wp -> "W+" | Wm -> "W-"
          end
      | O f ->
          begin match f with
          | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
          | Aux_top (_,_,ch,n,v) -> "Aux_" ^ (if n then "t_" else "") ^ (
              begin match v with
              | TTGG -> "ttGG" | TBWA -> "tbWA" | TBWZ -> "tbWZ"
              | TTWW -> "ttWW" | BBWW -> "bbWW"
              | QGUG -> "qGuG" | QBUB -> "qBuB"
              | QW   -> "qW"   | DL   -> "dL"   | DR   -> "dR"
              | QUQD1L -> "quqd1L" | QUQD1R -> "quqd1R"
              | QUQD8L -> "quqd8L" | QUQD8R -> "quqd8R"
              end ) ^ ( if ch > 0 then "+" else if ch < 0 then "-" else "0" )
          end

    let flavor_to_TeX = function
      | M f ->
          begin match f with
          | L 1 -> "e^-" | L (-1) -> "e^+"
          | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
          | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
          | L _ -> invalid_arg
                "Modellib.NoH.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib.NoH.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib.NoH.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Modellib.NoH.flavor_to_TeX: invalid down type quark"
          end
      | G f ->
          begin match f with
          | Gl -> "g"
          | Ga -> "\\gamma" | Z -> "Z"
          | Wp -> "W^+" | Wm -> "W^-"
          end
      | O f ->
          begin match f with
          | Phip -> "\\phi^+" | Phim -> "\\phi^-" | Phi0 -> "\\phi^0" 
          | Aux_top (_,_,ch,n,v) -> "\\textnormal{Aux_" ^ (if n then "t_" else "") ^ (
              begin match v with
              | TTGG -> "ttGG" | TBWA -> "tbWA" | TBWZ -> "tbWZ"
              | TTWW -> "ttWW" | BBWW -> "bbWW"
              | QGUG -> "qGuG" | QBUB -> "qBuB"
              | QW   -> "qW"   | DL   -> "dL"   | DR   -> "dR"
              | QUQD1L -> "quqd1L" | QUQD1R -> "quqd1R"
              | QUQD8L -> "quqd8L" | QUQD8R -> "quqd8R"
              end ) ^ ( if ch > 0 then "^+" else if ch < 0 then "^-" else "^0" ) ^ "}"
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
          | Aux_top (_,_,ch,n,v) -> "aux_" ^ (if n then "t_" else "") ^ (
              begin match v with
              | TTGG -> "ttgg" | TBWA -> "tbwa" | TBWZ -> "tbwz"
              | TTWW -> "ttww" | BBWW -> "bbww"
              | QGUG -> "qgug" | QBUB -> "qbub"
              | QW   -> "qw"   | DL   -> "dl"   | DR   -> "dr"
              | QUQD1L -> "quqd1l" | QUQD1R -> "quqd1r"
              | QUQD8L -> "quqd8l" | QUQD8R -> "quqd8r"
              end ) ^ "_" ^ ( if ch > 0 then "p" else if ch < 0 then "m" else "0" )
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
          | Aux_top (_,_,ch,t,f) -> let n =
            begin match f with
            | QW -> 0
            | QUQD1R -> 1 | QUQD1L -> 2
            | QUQD8R -> 3 | QUQD8L -> 4
            | _ -> 5
            end
            in (602 + 3*n - ch) * ( if t then (1) else (-1) )
          end

    let mass_symbol f = 
      "mass(" ^ string_of_int (abs (pdg f)) ^ ")"

    let width_symbol f =
      "width(" ^ string_of_int (abs (pdg f)) ^ ")"

    let constant_symbol = function
      | Unit -> "unit" | Half -> "half" | Pi -> "PI"
      | Alpha_QED -> "alpha" | E -> "e" | G_weak -> "g" | Vev -> "vev"
      | I_G_weak -> "ig" 
      | Sin2thw -> "sin2thw" | Sinthw -> "sinthw" | Costhw -> "costhw"
      | Q_lepton -> "qlep" | Q_up -> "qup" | Q_down -> "qdwn"
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_TVA_ttA -> "gtva_tta" | G_TVA_bbA -> "gtva_bba" 
      | G_VLR_ttZ -> "gvlr_ttz" | G_TVA_ttZ -> "gtva_ttz" | G_TVA_bbZ -> "gtva_bbz"
      | G_VLR_btW -> "gvlr_btw" | G_VLR_tbW -> "gvlr_tbw"
      | G_TLR_btW -> "gtlr_btw" | G_TRL_tbW -> "gtrl_tbw"
      | G_TLR_btWA -> "gtlr_btwa" | G_TRL_tbWA -> "gtrl_tbwa"
      | G_TLR_btWZ -> "gtlr_btwz" | G_TRL_tbWZ -> "gtrl_tbwz"
      | G_TVA_ttWW -> "gtva_ttww" | G_TVA_bbWW -> "gtva_bbww"
      | G_TVA_ttG -> "gtva_ttg" | G_TVA_ttGG -> "gtva_ttgg"
      | G_VLR_qGuG -> "gvlr_qgug"
      | G_VLR_qBuB -> "gvlr_qbub"
      | G_VLR_qBuB_u -> "gvlr_qbub_u" | G_VLR_qBuB_d -> "gvlr_qbub_d"
      | G_VLR_qBuB_e -> "gvlr_qbub_e" | G_VL_qBuB_n -> "gvl_qbub_n"
      | G_VL_qW -> "gvl_qw"
      | G_VL_qW_u -> "gvl_qw_u" | G_VL_qW_d -> "gvl_qw_d"
      | G_SL_DttR -> "gsl_dttr" | G_SR_DttR -> "gsr_dttr" | G_SL_DttL -> "gsl_dttl"
      | G_SLR_DbtR -> "gslr_dbtr" | G_SL_DbtL -> "gsl_dbtl"
      | C_quqd1R_bt -> "c_quqd1_1" | C_quqd1R_tb -> "conjg(c_quqd1_1)"
      | C_quqd1L_bt -> "conjg(c_quqd1_2)" | C_quqd1L_tb -> "c_quqd1_2"
      | C_quqd8R_bt -> "c_quqd8_1" | C_quqd8R_tb -> "conjg(c_quqd8_1)"
      | C_quqd8L_bt -> "conjg(c_quqd8_2)" | C_quqd8L_tb -> "c_quqd8_2"
      | G_CC -> "gcc"
      | G_CCQ (n1,n2) -> "gccq" ^ string_of_int n1 ^ string_of_int n2
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | I_G1_AWW -> "ig1a" | I_G1_ZWW -> "ig1z"
      | I_G1_plus_kappa_plus_G4_AWW -> "ig1pkpg4a"
      | I_G1_plus_kappa_plus_G4_ZWW -> "ig1pkpg4z"
      | I_G1_plus_kappa_minus_G4_AWW -> "ig1pkmg4a"
      | I_G1_plus_kappa_minus_G4_ZWW -> "ig1pkmg4z"
      | I_G1_minus_kappa_plus_G4_AWW -> "ig1mkpg4a"
      | I_G1_minus_kappa_plus_G4_ZWW -> "ig1mkpg4z"
      | I_G1_minus_kappa_minus_G4_AWW -> "ig1mkmg4a"
      | I_G1_minus_kappa_minus_G4_ZWW -> "ig1mkmg4z"
      | I_lambda_AWW -> "ila"
      | I_lambda_ZWW -> "ilz"
      | G5_AWW -> "rg5a"
      | G5_ZWW -> "rg5z"
      | I_kappa5_AWW -> "ik5a"
      | I_kappa5_ZWW -> "ik5z"
      | I_lambda5_AWW -> "il5a" | I_lambda5_ZWW -> "il5z"
      | Alpha_WWWW0 -> "alww0" | Alpha_WWWW2 -> "alww2"
      | Alpha_ZZWW0 -> "alzw0" | Alpha_ZZWW1 -> "alzw1"
      | Alpha_ZZZZ  -> "alzz"
      | D_Alpha_ZZWW0_S -> "dalzz0_s(gkm,mkm,"
      | D_Alpha_ZZWW0_T -> "dalzz0_t(gkm,mkm,"
      | D_Alpha_ZZWW1_S -> "dalzz1_s(gkm,mkm,"
      | D_Alpha_ZZWW1_T -> "dalzz1_t(gkm,mkm,"
      | D_Alpha_ZZWW1_U -> "dalzz1_u(gkm,mkm,"
      | D_Alpha_WWWW0_S -> "dalww0_s(gkm,mkm,"
      | D_Alpha_WWWW0_T -> "dalww0_t(gkm,mkm,"
      | D_Alpha_WWWW0_U -> "dalww0_u(gkm,mkm,"
      | D_Alpha_WWWW2_S -> "dalww2_s(gkm,mkm,"
      | D_Alpha_WWWW2_T -> "dalww2_t(gkm,mkm,"
      | D_Alpha_ZZZZ_S -> "dalz4_s(gkm,mkm,"
      | D_Alpha_ZZZZ_T -> "dalz4_t(gkm,mkm,"
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f
      | K_Matrix_Coeff i -> "kc" ^ string_of_int i
      | K_Matrix_Pole i -> "kp" ^ string_of_int i

  end

(* \thocwmodulesection{Minimal Higgsless Model including additional Resonances} *)

module AltH (Flags : NoH_flags) =
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
        "cms_width", Arg.Unit (fun () -> default_width := Complex_Mass),
        "use complex mass scheme"]
    let caveats () = []

    type f_aux_top = TTGG | TBWA | TBWZ | TTWW | BBWW |   (*i top auxiliary field "flavors" *)
                     QGUG | QBUB | QW | DL | DR

    type matter_field = L of int | N of int | U of int | D of int
    type gauge_boson = Ga | Wp | Wm | Z | Gl
    type other = Phip | Phim | Phi0 
                 | Rsigma | Rphin | Rphip | Rphim | Rphipp | Rphimm 
                 | Rf | Rtn | Rtp | Rtm | Rtpp | Rtmm 
                 | Aux_top of int*int*int*bool*f_aux_top    (*i lorentz*color*charge*top-side*flavor *)
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
      failwith "Modellib_NoH.AltH.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; N n; U n; D n ]

    let rec aux_top_flavors (f,l,co,ch) = List.append
      ( List.map other [ Aux_top(l,co,ch/2,true,f); Aux_top(l,co,ch/2,false,f) ] )
      ( if ch > 1 then List.append
          ( List.map other [ Aux_top(l,co,-ch/2,true,f); Aux_top(l,co,-ch/2,false,f) ] )
          ( aux_top_flavors (f,l,co,(ch-2)) )
        else [] )

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Gauge Bosons", List.map gauge_boson [Ga; Z; Wp; Wm; Gl];
	"Scalar Resonances", List.map other [Rsigma; Rphin; Rphip; Rphim; Rphipp; Rphimm];
	"Tensor Resonances", List.map other [Rf; Rtn; Rtp; Rtm; Rtpp; Rtmm];
        "Goldstone Bosons", List.map other [Phip; Phim; Phi0] ]

    let flavors () = List.append
      ( ThoList.flatmap snd (external_flavors ()) )
      ( ThoList.flatmap aux_top_flavors
         [ (TTGG,2,1,1); (TBWA,2,0,2); (TBWZ,2,0,2); (TTWW,2,0,1); (BBWW,2,0,1);
           (QGUG,1,1,1); (QBUB,1,0,1); (QW,1,0,3); (DL,0,0,3); (DR,0,0,3) ] )

    let spinor n =
      if n >= 0 then
        Spinor
      else
        ConjSpinor

    let lorentz_aux = function
      | 2 -> Tensor_1
      | 1 -> Vector
      | 0 -> Scalar
      | _ -> invalid_arg ("SM.lorentz_aux: wrong value")

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
          begin match f with
          | Aux_top (l,_,_,_,_) -> lorentz_aux l
          | Rf | Rtn | Rtp | Rtm | Rtpp | Rtmm -> Tensor_2
          | _ -> Scalar
          end

    let color = function 
      | M (U n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D n) -> Color.SUN (if n > 0 then 3 else -3)
      | G Gl -> Color.AdjSUN 3
      | O (Aux_top (_,co,_,_,_)) -> if co == 0 then Color.Singlet else Color.AdjSUN 3
      | _ -> Color.Singlet

    let nc () = 3

    let prop_spinor n =
      if n >= 0 then
        Prop_Spinor
      else
        Prop_ConjSpinor

    let prop_aux = function
      | 2 -> Aux_Tensor_1
      | 1 -> Aux_Vector
      | 0 -> Aux_Scalar
      | _ -> invalid_arg ("SM.prop_aux: wrong value")

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
          | Rsigma -> Prop_Scalar
          | Rphin | Rphip | Rphim | Rphipp | Rphimm -> Prop_Scalar
          | Rf -> Prop_Tensor_2
          | Rtn | Rtp | Rtm | Rtpp | Rtmm -> Prop_Tensor_2
          | Aux_top (l,_,_,_,_) -> prop_aux l
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
          |  Rsigma -> Rsigma
          | Rphin -> Rphin | Rphip -> Rphim | Rphim -> Rphip
          | Rphipp -> Rphimm | Rphimm -> Rphipp
          | Rf -> Rf
          | Rtn -> Rtn | Rtp -> Rtm | Rtm -> Rtp
          | Rtpp -> Rtmm | Rtmm -> Rtpp
          | Aux_top (l,co,ch,n,f) -> Aux_top (l,co,(-ch),(not n),f)
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
      |  n -> invalid_arg ("SM.generation': " ^ string_of_int n)

    let generation f =
      if Flags.ckm_present then
        []
      else
        match f with
        | M (L n | N n | U n | D n) -> generation' n
        | G _ | O _ -> [0//1; 0//1; 0//1]

    let charge = function
      | M f ->
          begin match f with
          | L n -> if n > 0 then -1//1 else  1//1
          | N n -> 0//1
          | U n -> if n > 0 then  2//3 else -2//3
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
          | Rsigma | Phi0 | Rphin | Rf | Rtn ->  0//1
          | Phip | Rphip | Rtp ->  1//1
          | Phim | Rphim | Rtm -> -1//1
          | Rphipp | Rtpp ->  2//1
          | Rphimm | Rtmm -> -2//1
          | Aux_top (_,_,ch,_,_) -> ch//1
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
      | Unit | Half | Pi | Alpha_QED | Sin2thw
      | Sinthw | Costhw | E | G_weak | I_G_weak | Vev
      | Q_lepton | Q_up | Q_down | G_CC | G_CCQ of int*int
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down 
      | G_TVA_ttA | G_TVA_bbA 
      | G_VLR_ttZ | G_TVA_ttZ | G_TVA_bbZ 
      | G_VLR_btW | G_VLR_tbW
      | G_TLR_btW | G_TRL_tbW
      | G_TLR_btWZ | G_TRL_tbWZ
      | G_TLR_btWA | G_TRL_tbWA
      | G_TVA_ttWW | G_TVA_bbWW
      | G_TVA_ttG | G_TVA_ttGG
      | G_VLR_qGuG | G_VLR_qBuB
      | G_VLR_qBuB_u | G_VLR_qBuB_d | G_VLR_qBuB_e | G_VL_qBuB_n
      | G_VL_qW | G_VL_qW_u | G_VL_qW_d
      | G_SL_DttR | G_SR_DttR | G_SL_DttL | G_SLR_DbtR | G_SL_DbtL
      | I_Q_W | I_G_ZWW
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | I_G1_AWW | I_G1_ZWW
      | I_G1_plus_kappa_plus_G4_AWW
      | I_G1_plus_kappa_plus_G4_ZWW
      | I_G1_plus_kappa_minus_G4_AWW
      | I_G1_plus_kappa_minus_G4_ZWW
      | I_G1_minus_kappa_plus_G4_AWW
      | I_G1_minus_kappa_plus_G4_ZWW
      | I_G1_minus_kappa_minus_G4_AWW
      | I_G1_minus_kappa_minus_G4_ZWW
      | I_lambda_AWW | I_lambda_ZWW
      | G5_AWW | G5_ZWW
      | I_kappa5_AWW | I_kappa5_ZWW 
      | I_lambda5_AWW | I_lambda5_ZWW
      | Alpha_WWWW0 | Alpha_ZZWW1 | Alpha_WWWW2
      | Alpha_ZZWW0 | Alpha_ZZZZ
      | D_Alpha_ZZWW0_S | D_Alpha_ZZWW0_T | D_Alpha_ZZWW1_S
      | D_Alpha_ZZWW1_T | D_Alpha_ZZWW1_U | D_Alpha_WWWW0_S
      | D_Alpha_WWWW0_T | D_Alpha_WWWW0_U | D_Alpha_WWWW2_S
      | D_Alpha_WWWW2_T | D_Alpha_ZZZZ_S | D_Alpha_ZZZZ_T
      | G_SWW | G_SWW_T | G_SSWW | G_SZZ | G_SZZ_T | G_SSZZ      
      | G_PNWW | G_PNZZ | G_PWZ | G_PWW
      | G_FWW | G_FZZ | G_FWW_T | G_FZZ_T
      | G_TNWW | G_TNZZ | G_TWZ | G_TWW
      | Gs | I_Gs | G2
      | Mass of flavor | Width of flavor
      | K_Matrix_Coeff of int | K_Matrix_Pole of int

(* \begin{dubious}
     The current abstract syntax for parameter dependencies is admittedly
     tedious. Later, there will be a parser for a convenient concrete syntax
     as a part of a concrete syntax for models.  But as these examples show,
     it should include simple functions.
   \end{dubious} *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

(* \begin{subequations}
     \begin{align}
        \alpha_{\text{QED}} &= \frac{1}{137.0359895} \\
             \sin^2\theta_w &= 0.23124
     \end{align}
   \end{subequations} *)
    let input_parameters =
      [ Alpha_QED, 1. /. 137.0359895;
        Sin2thw, 0.23124;
        Mass (G Z), 91.187;
        Mass (M (N 1)), 0.0; Mass (M (L 1)), 0.51099907e-3;
        Mass (M (N 2)), 0.0; Mass (M (L 2)), 0.105658389;
        Mass (M (N 3)), 0.0; Mass (M (L 3)), 1.77705;
        Mass (M (U 1)), 5.0e-3; Mass (M (D 1)), 3.0e-3;
        Mass (M (U 2)), 1.2; Mass (M (D 2)), 0.1;
        Mass (M (U 3)), 174.0; Mass (M (D 3)), 4.2 ]

(* \begin{subequations}
     \begin{align}
                        e &= \sqrt{4\pi\alpha} \\
             \sin\theta_w &= \sqrt{\sin^2\theta_w} \\
             \cos\theta_w &= \sqrt{1-\sin^2\theta_w} \\
                        g &= \frac{e}{\sin\theta_w} \\
                      m_W &= \cos\theta_w m_Z \\
                        v &= \frac{2m_W}{g} \\
                  g_{CC}   =
       -\frac{g}{2\sqrt2} &= -\frac{e}{2\sqrt2\sin\theta_w} \\
       Q_{\text{lepton}}   =
      -q_{\text{lepton}}e &= e \\
           Q_{\text{up}}   =
          -q_{\text{up}}e &= -\frac{2}{3}e \\
         Q_{\text{down}}   =
        -q_{\text{down}}e &= \frac{1}{3}e \\
        \ii q_We           =
        \ii g_{\gamma WW} &= \ii e \\
              \ii g_{ZWW} &= \ii g \cos\theta_w \\
              \ii g_{WWW} &= \ii g
     \end{align}
   \end{subequations} *)

    let derived_parameters =
      [ Real E, Sqrt (Prod [Integer 4; Atom Pi; Atom Alpha_QED]);
        Real Sinthw, Sqrt (Atom Sin2thw);
        Real Costhw, Sqrt (Diff (Integer 1, Atom Sin2thw));
        Real G_weak, Quot (Atom E, Atom Sinthw);
        Real (Mass (G Wp)), Prod [Atom Costhw; Atom (Mass (G Z))];
        Real Vev, Quot (Prod [Integer 2; Atom (Mass (G Wp))], Atom G_weak);
        Real Q_lepton, Atom E;
        Real Q_up, Prod [Quot (Integer (-2), Integer 3); Atom E];
        Real Q_down, Prod [Quot (Integer 1, Integer 3); Atom E];
        Real G_CC, Neg (Quot (Atom G_weak, Prod [Integer 2; Sqrt (Integer 2)]));
        Complex I_Q_W, Prod [I; Atom E];
        Complex I_G_weak, Prod [I; Atom G_weak];
        Complex I_G_ZWW, Prod [I; Atom G_weak; Atom Costhw] ]
             
(* \begin{equation}
      - \frac{g}{2\cos\theta_w}
   \end{equation} *)
    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

(* \begin{subequations}
     \begin{align}
           - \frac{g}{2\cos\theta_w} g_V
        &= - \frac{g}{2\cos\theta_w} (T_3 - 2 q \sin^2\theta_w) \\
           - \frac{g}{2\cos\theta_w} g_A
        &= - \frac{g}{2\cos\theta_w} T_3
     \end{align}
   \end{subequations} *)
    let nc_coupling c t3 q =
      (Real_Array c,
       [Prod [g_over_2_costh; Diff (t3, Prod [Integer 2; q; Atom Sin2thw])];
        Prod [g_over_2_costh; t3]])

    let half = Quot (Integer 1, Integer 2)

    let derived_parameter_arrays =
      [ nc_coupling G_NC_neutrino half (Integer 0);
        nc_coupling G_NC_lepton (Neg half) (Integer (-1));
        nc_coupling G_NC_up half (Quot (Integer 2, Integer 3));
        nc_coupling G_NC_down (Neg half) (Quot (Integer (-1), Integer 3)) ]

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
    let mom ((m1, o, m2), fbf, c) = ((M m1, O o, M m2), fbf, c)

    let electromagnetic_currents n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]
        
    let color_currents n =
      List.map mgm
        [ ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs);
          ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs) ]

(* \begin{equation}
     \mathcal{L}_{\textrm{NC}} =
        - \frac{g}{2\cos\theta_W}
            \sum_i \bar\psi_i\fmslash{Z}(g_V^i-g_A^i\gamma_5)\psi_i
   \end{equation} *)

    let neutral_currents n =
      List.map mgm
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

(* \begin{equation}
     \mathcal{L}_{\textrm{CC}} =
        - \frac{g}{2\sqrt2} \sum_i \bar\psi_i
               (T^+\fmslash{W}^+ + T^-\fmslash{W}^-)(1-\gamma_5)\psi_i 
   \end{equation} *)

    let charged_currents' n = 
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let charged_currents'' n =
      List.map mgm 
        [ ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let charged_currents_triv = 
      ThoList.flatmap charged_currents' [1;2;3] @
      ThoList.flatmap charged_currents'' [1;2;3]

    let charged_currents_ckm = 
      let charged_currents_2 n1 n2 = 
        List.map mgm 
          [ ((D (-n1), Wm, U n2), FBF (1, Psibar, VL, Psi), G_CCQ (n2,n1));
            ((U (-n1), Wp, D n2), FBF (1, Psibar, VL, Psi), G_CCQ (n1,n2)) ] in
      ThoList.flatmap charged_currents' [1;2;3] @ 
      List.flatten (Product.list2 charged_currents_2 [1;2;3] [1;2;3])
      
(* \begin{equation}
     \mathcal{L}_{\textrm{TGC}} =
        - e \partial_\mu A_\nu W_+^\mu W_-^\nu + \ldots
        - e \cot\theta_w  \partial_\mu Z_\nu W_+^\mu W_-^\nu + \ldots
   \end{equation} *)

    let tgc ((g1, g2, g3), t, c) = ((G g1, G g2, G g3), t, c)

    let standard_triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs)]

(* \begin{multline}
     \mathcal{L}_{\textrm{TGC}}(g_1,\kappa)
        =   g_1 \mathcal{L}_T(V,W^+,W^-) \\
          + \frac{\kappa+g_1}{2} \Bigl(\mathcal{L}_T(W^-,V,W^+)
                                         - \mathcal{L}_T(W^+,V,W^-)\Bigr)\\
          + \frac{\kappa-g_1}{2} \Bigl(\mathcal{L}_L(W^-,V,W^+)
                                         - \mathcal{L}_T(W^+,V,W^-)\Bigr)
   \end{multline} *)

(* \begin{dubious}
   The whole thing in the LEP2 workshop notation:
   \begin{multline}
     \ii\mathcal{L}_{\textrm{TGC},V} / g_{WWV} = \\
            g_1^V V^\mu (W^-_{\mu\nu}W^{+,\nu}-W^+_{\mu\nu}W^{-,\nu})
          + \kappa_V  W^+_\mu W^-_\nu V^{\mu\nu}
          + \frac{\lambda_V}{m_W^2} V_{\mu\nu}
               W^-_{\rho\mu} W^{+,\hphantom{\nu}\rho}_{\hphantom{+,}\nu} \\
          + \ii g_5^V \epsilon_{\mu\nu\rho\sigma}
              \left(   (\partial^\rho W^{-,\mu}) W^{+,\nu}
                     -  W^{-,\mu}(\partial^\rho W^{+,\nu}) \right) V^\sigma \\
          + \ii g_4^V W^-_\mu W^+_\nu (\partial^\mu V^\nu + \partial^\nu V^\mu)
          - \frac{\tilde\kappa_V}{2}  W^-_\mu W^+_\nu \epsilon^{\mu\nu\rho\sigma}
              V_{\rho\sigma}
          - \frac{\tilde\lambda_V}{2m_W^2}
               W^-_{\rho\mu} W^{+,\mu}_{\hphantom{+,\mu}\nu} \epsilon^{\nu\rho\alpha\beta}
                V_{\alpha\beta}
   \end{multline}
   using the conventions of Itzykson and Zuber with $\epsilon^{0123} = +1$.
   \end{dubious} *)

(* \begin{dubious}
   This is equivalent to the notation of Hagiwara et al.~\cite{HPZH87}, if we
   remember that they have opposite signs for~$g_{WWV}$:
   \begin{multline}
     \mathcal{L}_{WWV} / (-g_{WWV})  = \\
       \ii g_1^V \left( W^\dagger_{\mu\nu} W^\mu 
                         - W^\dagger_\mu W^\mu_{\hphantom{\mu}\nu} \right) V^\nu
     + \ii \kappa_V  W^\dagger_\mu W_\nu V^{\mu\nu}
     + \ii \frac{\lambda_V}{m_W^2}
          W^\dagger_{\lambda\mu} W^\mu_{\hphantom{\mu}\nu} V^{\nu\lambda} \\
     - g_4^V  W^\dagger_\mu W_\nu
          \left(\partial^\mu V^\nu + \partial^\nu V^\mu \right)
     + g_5^V \epsilon^{\mu\nu\lambda\sigma}
           \left( W^\dagger_\mu \stackrel{\leftrightarrow}{\partial_\lambda}
                  W_\nu \right) V_\sigma\\
     + \ii \tilde\kappa_V  W^\dagger_\mu W_\nu \tilde{V}^{\mu\nu}
     + \ii\frac{\tilde\lambda_V}{m_W^2}
           W^\dagger_{\lambda\mu} W^\mu_{\hphantom{\mu}\nu} \tilde{V}^{\nu\lambda}
   \end{multline}
   Here $V^\mu$ stands for either the photon or the~$Z$ field, $W^\mu$ is the
   $W^-$ field, $W_{\mu\nu} = \partial_\mu W_\nu - \partial_\nu W_\mu$,
   $V_{\mu\nu} = \partial_\mu V_\nu - \partial_\nu V_\mu$, and
   $\tilde{V}_{\mu\nu} = \frac{1}{2} \epsilon_{\mu\nu\lambda\sigma}
   V^{\lambda\sigma}$.
   \end{dubious} *)

    let anomalous_triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Dim4_Vector_Vector_Vector_T (-1),
           I_G1_AWW);
          ((Z, Wm, Wp), Dim4_Vector_Vector_Vector_T (-1),
           I_G1_ZWW);
          ((Wm, Ga, Wp), Dim4_Vector_Vector_Vector_T 1,
           I_G1_plus_kappa_minus_G4_AWW);
          ((Wm, Z, Wp), Dim4_Vector_Vector_Vector_T 1,
           I_G1_plus_kappa_minus_G4_ZWW);
          ((Wp, Ga, Wm), Dim4_Vector_Vector_Vector_T (-1),
           I_G1_plus_kappa_plus_G4_AWW);
          ((Wp, Z, Wm), Dim4_Vector_Vector_Vector_T (-1),
           I_G1_plus_kappa_plus_G4_ZWW);
          ((Wm, Ga, Wp), Dim4_Vector_Vector_Vector_L (-1),
           I_G1_minus_kappa_plus_G4_AWW);
          ((Wm, Z, Wp), Dim4_Vector_Vector_Vector_L (-1),
           I_G1_minus_kappa_plus_G4_ZWW);
          ((Wp, Ga, Wm), Dim4_Vector_Vector_Vector_L 1,
           I_G1_minus_kappa_minus_G4_AWW);
          ((Wp, Z, Wm), Dim4_Vector_Vector_Vector_L 1,
           I_G1_minus_kappa_minus_G4_ZWW);
          ((Ga, Wm, Wp), Dim4_Vector_Vector_Vector_L5 (-1),
           I_kappa5_AWW);
          ((Z, Wm, Wp), Dim4_Vector_Vector_Vector_L5 (-1),
           I_kappa5_ZWW);
          ((Ga, Wm, Wp), Dim4_Vector_Vector_Vector_T5 (-1),
           G5_AWW);
          ((Z, Wm, Wp), Dim4_Vector_Vector_Vector_T5 (-1),
           G5_ZWW);
          ((Ga, Wp, Wm), Dim6_Gauge_Gauge_Gauge (-1),
           I_lambda_AWW);
          ((Z, Wp, Wm), Dim6_Gauge_Gauge_Gauge (-1),
           I_lambda_ZWW);
          ((Ga, Wp, Wm), Dim6_Gauge_Gauge_Gauge_5 (-1),
           I_lambda5_AWW);
          ((Z, Wp, Wm), Dim6_Gauge_Gauge_Gauge_5 (-1),
           I_lambda5_ZWW) ]

    let triple_gauge =
      if Flags.triple_anom then
        anomalous_triple_gauge
      else
        standard_triple_gauge

(* \begin{equation}
     \mathcal{L}_{\textrm{QGC}} =
        - g^2 W_{+,\mu} W_{-,\nu} W_+^\mu W_-^\nu + \ldots
   \end{equation} *)

(* Actually, quartic gauge couplings are a little bit more straightforward
   using auxiliary fields.  Here we have to impose the antisymmetry manually:
   \begin{subequations}
   \begin{multline}
     (W^{+,\mu}_1 W^{-,\nu}_2 - W^{+,\nu}_1 W^{-,\mu}_2)
     (W^+_{3,\mu} W^-_{4,\nu} - W^+_{3,\nu} W^-_{4,\mu}) \\
        = 2(W^+_1W^+_3)(W^-_2W^-_4) - 2(W^+_1W^-_4)(W^-_2W^+_3)
   \end{multline}
   also ($V$ can be $A$ or $Z$)
   \begin{multline}
     (W^{+,\mu}_1 V^\nu_2 - W^{+,\nu}_1 V^\mu_2)
     (W^-_{3,\mu} V_{4,\nu} - W^-_{3,\nu} V_{4,\mu}) \\
        = 2(W^+_1W^-_3)(V_2V_4) - 2(W^+_1V_4)(V_2W^-_3)
   \end{multline}
   \end{subequations} *)

(* \begin{subequations}
   \begin{multline}
      W^{+,\mu} W^{-,\nu} W^+_\mu W^-_\nu
   \end{multline}
   \end{subequations} *)

    let qgc ((g1, g2, g3, g4), t, c) = ((G g1, G g2, G g3, G g4), t, c)

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let standard_quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          (Gl, Gl, Gl, Gl), gauge4, G2 ]

(* \begin{subequations}
   \begin{align}
     \mathcal{L}_4
       &= \alpha_4 \left(   \frac{g^4}{2}\left(   (W^+_\mu W^{-,\mu})^2
                                                + W^+_\mu W^{+,\mu} W^-_\mu W^{-,\mu}
                                               \right)\right.\notag \\
       &\qquad\qquad\qquad \left.
                          + \frac{g^4}{\cos^2\theta_w} W^+_\mu Z^\mu W^-_\nu Z^\nu
                          + \frac{g^4}{4\cos^4\theta_w} (Z_\mu Z^\mu)^2 \right) \\
     \mathcal{L}_5
       &= \alpha_5 \left(   g^4 (W^+_\mu W^{-,\mu})^2
                          + \frac{g^4}{\cos^2\theta_w}  W^+_\mu W^{-,\mu} Z_\nu Z^\nu
                          + \frac{g^4}{4\cos^4\theta_w} (Z_\mu Z^\mu)^2 \right)
   \end{align}
   \end{subequations}
   or
   \begin{multline}
     \mathcal{L}_4 + \mathcal{L}_5
       =   (\alpha_4+2\alpha_5) g^4 \frac{1}{2} (W^+_\mu W^{-,\mu})^2 \\
         + 2\alpha_4 g^4 \frac{1}{4} W^+_\mu W^{+,\mu} W^-_\mu W^{-,\mu}
         + \alpha_4 \frac{g^4}{\cos^2\theta_w} W^+_\mu Z^\mu W^-_\nu Z^\nu \\
         + 2\alpha_5 \frac{g^4}{\cos^2\theta_w} \frac{1}{2} W^+_\mu W^{-,\mu} Z_\nu Z^\nu
         + (2\alpha_4 + 2\alpha_5) \frac{g^4}{\cos^4\theta_w} \frac{1}{8} (Z_\mu Z^\mu)^2
   \end{multline}
   and therefore
   \begin{subequations}
   \begin{align}
     \alpha_{(WW)_0} &= (\alpha_4+2\alpha_5) g^4 \\
     \alpha_{(WW)_2} &= 2\alpha_4 g^4 \\
     \alpha_{(WZ)_0} &= 2\alpha_5 \frac{g^4}{\cos^2\theta_w} \\
     \alpha_{(WZ)_1} &= \alpha_4 \frac{g^4}{\cos^2\theta_w} \\
     \alpha_{ZZ} &= (2\alpha_4 + 2\alpha_5) \frac{g^4}{\cos^4\theta_w}
   \end{align}
   \end{subequations} *)

    let anomalous_quartic_gauge =
      if Flags.quartic_anom then
        List.map qgc
          [ ((Wm, Wm, Wp, Wp),
             Vector4 [(1, C_13_42); (1, C_14_23)], Alpha_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Vector4 [1, C_12_34], Alpha_WWWW2);
            ((Wm, Wp, Z, Z),
             Vector4 [1, C_12_34], Alpha_ZZWW0);
            ((Wm, Wp, Z, Z),
             Vector4 [(1, C_13_42); (1, C_14_23)], Alpha_ZZWW1);
            ((Z, Z, Z, Z),
             Vector4 [(1, C_12_34); (1, C_13_42); (1, C_14_23)], Alpha_ZZZZ) ]
      else
        []

(* In any diagonal channel~$\chi$, the scattering amplitude~$a_\chi(s)$ is
   unitary iff\footnote{%
     Trivial proof:
     \begin{equation}
       -1 = \textrm{Im}\left(\frac{1}{a_\chi(s)}\right)
          = \frac{\textrm{Im}(a_\chi^*(s))}{ |a_\chi(s)|^2 }
          = - \frac{\textrm{Im}(a_\chi(s))}{ |a_\chi(s)|^2 }
     \end{equation}
     i.\,e.~$\textrm{Im}(a_\chi(s)) = |a_\chi(s)|^2$.}
   \begin{equation}
     \textrm{Im}\left(\frac{1}{a_\chi(s)}\right) = -1
   \end{equation}
   For a real perturbative scattering amplitude~$r_\chi(s)$ this can be
   enforced easily--and arbitrarily--by
   \begin{equation}
     \frac{1}{a_\chi(s)} = \frac{1}{r_\chi(s)} - \mathrm{i}
   \end{equation} 

*)


    let k_matrix_quartic_gauge =
      if Flags.k_matrix then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_14_23)]), D_Alpha_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42)]), D_Alpha_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_jr (0,
                   [(1, C_14_23)]), D_Alpha_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42)]), D_Alpha_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42); (1, C_14_23)]), D_Alpha_WWWW2_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42); (1, C_14_23)]), D_Alpha_ZZWW0_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42)]), D_Alpha_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_14_23)]), D_Alpha_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_jr (1,
                   [(1, C_12_34)]), D_Alpha_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_jr (1,
                   [(1, C_13_42)]), D_Alpha_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_jr (1,
                   [(1, C_14_23)]), D_Alpha_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_jr (2,
                   [(1, C_12_34)]), D_Alpha_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_jr (2,
                   [(1, C_13_42)]), D_Alpha_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_jr (2,
                   [(1, C_14_23)]), D_Alpha_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_12_34)]), D_Alpha_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_jr (0,
                   [(1, C_13_42); (1, C_14_23)]), D_Alpha_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_jr (3,
                   [(1, C_14_23)]), D_Alpha_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_jr (3,
                   [(1, C_13_42); (1, C_12_34)]), D_Alpha_ZZZZ_T)]
      else
        []



(*i Thorsten's original implementation of the K matrix, which we keep since
   it still might be usefull for the future. 


    let k_matrix_quartic_gauge =
      if Flags.k_matrix then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_tho (0, [K_Matrix_Coeff 0, 
                         K_Matrix_Pole 0]), Alpha_WWWW0);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_tho (0, [K_Matrix_Coeff 2, 
                         K_Matrix_Pole 2]), Alpha_WWWW2);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_tho (0, [(K_Matrix_Coeff 0, 
                         K_Matrix_Pole 0); (K_Matrix_Coeff 2, 
                         K_Matrix_Pole 2)]), Alpha_ZZWW0);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_tho (0, [K_Matrix_Coeff 1, 
                         K_Matrix_Pole 1]), Alpha_ZZWW1);
            ((Z, Z, Z, Z), Vector4_K_Matrix_tho (0, [K_Matrix_Coeff 0, 
                         K_Matrix_Pole 0]), Alpha_ZZZZ) ]
      else
        []

i*)

    let quartic_gauge =
      standard_quartic_gauge @ anomalous_quartic_gauge @ k_matrix_quartic_gauge

(* WK's couplings (apparently, he still intends to divide by
   $\Lambda^2_{\text{EWSB}}=16\pi^2v_{\mathrm{F}}^2$):
   with
   \begin{equation}
      V_{\mu} V_{\nu} =
        \frac{1}{2} \left( W^+_{\mu} W^-_{\nu} + W^+_{\nu} W^-_{\mu} \right)
         + \frac{1}{2\cos^2\theta_{w}} Z_{\mu} Z_{\nu}
   \end{equation}
   (note the symmetrization!), i.\,e.
   \begin{subequations}
   \begin{align}
     \mathcal{L}_4 &= \alpha_4 \frac{g^4v_{\mathrm{F}}^4}{16} (V_{\mu} V_{\nu})^2 \\
     \mathcal{L}_5 &= \alpha_5 \frac{g^4v_{\mathrm{F}}^4}{16} (V_{\mu} V^{\mu})^2
   \end{align}
   \end{subequations} *)

    let goldstone_vertices =
      [ ((O Phi0, G Wm, G Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phip, G Ga, G Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phip, G Z, G Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phim, G Wp, G Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phim, G Wp, G Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

(* New Resonances *)

(*
  \begin{dubious}
    There is an extra minus in the Lagrangian to have the same sign as
    HWW or HZZ vertex. 
    Effectivly this doesn't matter for VBS, because $(-1)^2=1$.
    This is only for completeness.
  \end{dubious}
  \begin{subequations}
    \begin{align}
      \mathbf{V}_\mu &= -\mathrm{i} g\mathbf{W}_\mu+\mathrm{i} g^\prime\mathbf{B}_\mu \\
      \mathbf{W}_\mu &= W_\mu^a\frac{\tau^a}{2} \\
      \mathbf{B}_\mu &= W_\mu^a\frac{\tau^3}{2} \\
      \tau^{++}&= \tau^+ \otimes \tau^+ \\
      \tau^+ &= \frac{1}{2} \left (\tau^+ \otimes \tau^3 + \tau^3+\tau^+ \right ) \\
      \tau^0 &= \frac{1}{\sqrt{6}} \left (\tau^3\otimes\tau^3 -\tau^+ \otimes \tau^- - \tau^-+\tau^+ \right ) \\
      \tau^- &= \frac{1}{2} \left (\tau^- \otimes \tau^3 + \tau^3+\tau^- \right ) \\
      \tau^{--}&= \tau^- \otimes \tau^- 
    \end{align}
  \end{subequations}  
*)

(* Scalar Isoscalar
   \begin{equation}
    \mathcal{L}_{\sigma}=
              -\frac{g_\sigma v}{2} \text{tr}
	      \left\lbrack \mathbf{V}_\mu \mathbf{V}^\mu \right\rbrack \sigma
   \end{equation}
*)
    let rsigma3 =
      [ ((O Rsigma, G Wp, G Wm), Scalar_Vector_Vector 1, G_SWW);
        ((O Rsigma, G Z, G Z), Scalar_Vector_Vector 1, G_SZZ) ]

    let rsigma3t =
      [ ((O Rsigma, G Wp, G Wm), Scalar_Vector_Vector_t 1, G_SWW_T);
        ((O Rsigma, G Z, G Z), Scalar_Vector_Vector_t 1, G_SZZ_T) ]

    let rsigma4 =
      [ (O Rsigma, O Rsigma, G Wp, G Wm), Scalar2_Vector2 1, G_SSWW;
        (O Rsigma, O Rsigma, G Z, G Z), Scalar2_Vector2 1, G_SSZZ ]

(* Scalar Isotensor
  \begin{subequations}
    \begin{align}
     \mathcal{L}_{\phi}&=
              \frac{g_\phi v}{4} \text{Tr}
	      \left \lbrack \left ( \mathbf{V}_\mu \otimes \mathbf{V}^\mu - \frac{\tau^{aa}}{6} \text{Tr} \left \lbrack \mathbf{V}_\mu \mathbf{V}^\mu \right \rbrack\right ) {\mathbf{\phi}} \right \rbrack\\
     \phi&=\sqrt{2} \left (\phi^{++}\tau^{++}+\phi^+\tau^++\phi^0\tau^0+\phi^-\tau^- + \phi^{--}\tau^{--} \right )
    \end{align}
  \end{subequations}
*)
    let rphi3 =
      [ ((O Rphin, G Wp, G Wm), Scalar_Vector_Vector 1, G_PNWW);
        ((O Rphin, G Z, G Z), Scalar_Vector_Vector 1, G_PNZZ) ;
        ((O Rphip, G Z, G Wm), Scalar_Vector_Vector 1, G_PWZ) ;
        ((O Rphipp, G Wm, G Wm), Scalar_Vector_Vector 1, G_PWW) ;
        ((O Rphim, G Wp, G Z), Scalar_Vector_Vector 1, G_PWZ) ;
        ((O Rphimm, G Wp, G Wp), Scalar_Vector_Vector 1, G_PWW) ]

(* Tensor IsoScalar
*)
    let rf3 =
      [ ((O Rf, G Wp, G Wm), Tensor_2_Vector_Vector_1 1, G_FWW);
        ((O Rf, G Z, G Z), Tensor_2_Vector_Vector_1 1, G_FZZ) ]

    let rf3t =
      [ ((O Rf, G Wp, G Wm), Tensor_2_Vector_Vector_t 1, G_FWW_T);
        ((O Rf, G Z, G Z), Tensor_2_Vector_Vector_t 1, G_FZZ_T) ]

(* Tensor Isotensor
  \begin{subequations}
    \begin{align}
     \mathcal{L}_{t}
    \end{align}
  \end{subequations}
*)
    let rt3 =
      [ ((O Rtn, G Wp, G Wm), Tensor_2_Vector_Vector_1 1, G_TNWW);
        ((O Rtn, G Z, G Z), Tensor_2_Vector_Vector_1 1, G_TNZZ) ;
        ((O Rtp, G Z, G Wm), Tensor_2_Vector_Vector_1 1, G_TWZ) ;
        ((O Rtpp, G Wm, G Wm), Tensor_2_Vector_Vector_1 1, G_TWW) ;
        ((O Rtm, G Wp, G Z), Tensor_2_Vector_Vector_1 1, G_TWZ) ;
        ((O Rtmm, G Wp, G Wp), Tensor_2_Vector_Vector_1 1, G_TWW) ]


(* Anomalous trilinear interactions $f_i f_j V$ :
   \begin{equation}
     \Delta\mathcal{L}_{tt\gamma} =
        - e \frac{\upsilon}{\Lambda^2}
            \bar{t} i\sigma^{\mu\nu} k_\nu (d_V(k^2) + i d_A(k^2) \gamma_5) t A_\mu
   \end{equation} *)

    let anomalous_ttA =
      if Flags.top_anom then
        [ ((M (U (-3)), G Ga, M (U 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_ttA) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{bb\gamma} =
        - e \frac{\upsilon}{\Lambda^2}
            \bar{b} i\sigma^{\mu\nu} k_\nu (d_V(k^2) + i d_A(k^2) \gamma_5) b A_\mu
   \end{equation} *)

    let anomalous_bbA =
      if Flags.top_anom then
        [ ((M (D (-3)), G Ga, M (D 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_bbA) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{ttg} =
        - g_s \frac{\upsilon}{\Lambda^2}
            \bar{t}\lambda^a i\sigma^{\mu\nu}k_\nu
                (d_V(k^2)+id_A(k^2)\gamma_5)tG^a_\mu
   \end{equation} *)

    let anomalous_ttG =
      if Flags.top_anom then
        [ ((M (U (-3)), G Gl, M (U 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_ttG) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{ttZ} =
        - \frac{g}{2 c_W} \frac{\upsilon^2}{\Lambda^2}\left\lbrack
              \bar{t} \fmslash{Z} (X_L(k^2) P_L + X_R(k^2) P_R) t
            + \bar{t}\frac{i\sigma^{\mu\nu}k_\nu}{m_Z}
                  (d_V(k^2)+id_A(k^2)\gamma_5)tZ_\mu\right\rbrack
   \end{equation} *)

    let anomalous_ttZ =
      if Flags.top_anom then
        [ ((M (U (-3)), G Z, M (U 3)), FBF (1, Psibar, VLRM, Psi), G_VLR_ttZ);
          ((M (U (-3)), G Z, M (U 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_ttZ) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{bbZ} =
        - \frac{g}{2 c_W} \frac{\upsilon^2}{\Lambda^2}
              \bar{b}\frac{i\sigma^{\mu\nu}k_\nu}{m_Z}
                  (d_V(k^2)+id_A(k^2)\gamma_5)bZ_\mu
   \end{equation} *)

    let anomalous_bbZ =
      if Flags.top_anom then
        [ ((M (D (-3)), G Z, M (D 3)), FBF (1, Psibar, TVAM, Psi), G_TVA_bbZ) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{tbW} =
        - \frac{g}{\sqrt{2}} \frac{\upsilon^2}{\Lambda^2}\left\lbrack
            \bar{b}\fmslash{W}^-(V_L(k^2) P_L+V_R(k^2) P_R) t
          + \bar{b}\frac{i\sigma^{\mu\nu}k_\nu}{m_W}
                (g_L(k^2)P_L+g_R(k^2)P_R)tW^-_\mu\right\rbrack
        + \textnormal{H.c.}
   \end{equation} *)

    let anomalous_tbW =
      if Flags.top_anom then
        [ ((M (D (-3)), G Wm, M (U 3)), FBF (1, Psibar, VLRM, Psi), G_VLR_btW);
          ((M (U (-3)), G Wp, M (D 3)), FBF (1, Psibar, VLRM, Psi), G_VLR_tbW);
          ((M (D (-3)), G Wm, M (U 3)), FBF (1, Psibar, TLRM, Psi), G_TLR_btW);
          ((M (U (-3)), G Wp, M (D 3)), FBF (1, Psibar, TRLM, Psi), G_TRL_tbW) ]
      else
        []

(* quartic fermion-gauge interactions $f_i f_j V_1 V_2$ emerging from gauge-invariant
effective operators:
   \begin{equation}
     \Delta\mathcal{L}_{ttgg} =
        - \frac{g_s^2}{2} f_{abc} \frac{\upsilon}{\Lambda^2}
            \bar{t} \lambda^a \sigma^{\mu\nu}
                (d_V(k^2)+id_A(k^2)\gamma_5)t G^b_\mu G^c_\nu
   \end{equation} *)

    let anomalous_ttGG =
      if Flags.top_anom then
        [ ((M (U (-3)), O (Aux_top (2,1,0,true,TTGG)), M (U 3)), FBF (1, Psibar, TVA, Psi), G_TVA_ttGG);
          ((O (Aux_top (2,1,0,false,TTGG)), G Gl, G Gl), Aux_Gauge_Gauge 1, I_Gs) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{tbWA} =
        - i\sin\theta_w \frac{g^2}{2\sqrt{2}} \frac{\upsilon^2}{\Lambda^2}\left\lbrack
            \bar{b}\frac{\sigma^{\mu\nu}}{m_W}
                (g_L(k^2)P_L+g_R(k^2)P_R)t A_\mu W^-_\nu \right\rbrack
        + \textnormal{H.c.}
   \end{equation} *)

    let anomalous_tbWA =
      if Flags.top_anom then
        [ ((M (D (-3)), O (Aux_top (2,0,-1,true,TBWA)), M (U 3)), FBF (1, Psibar, TLR, Psi), G_TLR_btWA);
          ((O (Aux_top (2,0,1,false,TBWA)), G Ga, G Wm), Aux_Gauge_Gauge 1, I_G_weak);
          ((M (U (-3)), O (Aux_top (2,0,1,true,TBWA)), M (D 3)), FBF (1, Psibar, TRL, Psi), G_TRL_tbWA);
          ((O (Aux_top (2,0,-1,false,TBWA)), G Wp, G Ga), Aux_Gauge_Gauge 1, I_G_weak) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{tbWZ} =
        - i\cos\theta_w \frac{g^2}{2\sqrt{2}} \frac{\upsilon^2}{\Lambda^2}\left\lbrack
            \bar{b}\frac{\sigma^{\mu\nu}}{m_W}
                (g_L(k^2)P_L+g_R(k^2)P_R)t Z_\mu W^-_\nu \right\rbrack
        + \textnormal{H.c.}
   \end{equation} *)

    let anomalous_tbWZ =
      if Flags.top_anom then
        [ ((M (D (-3)), O (Aux_top (2,0,-1,true,TBWZ)), M (U 3)), FBF (1, Psibar, TLR, Psi), G_TLR_btWZ);
          ((O (Aux_top (2,0,1,false,TBWZ)), G Z, G Wm), Aux_Gauge_Gauge 1, I_G_weak);
          ((M (U (-3)), O (Aux_top (2,0,1,true,TBWZ)), M (D 3)), FBF (1, Psibar, TRL, Psi), G_TRL_tbWZ);
          ((O (Aux_top (2,0,-1,false,TBWZ)), G Wp, G Z), Aux_Gauge_Gauge 1, I_G_weak) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{ttWW} =
        - i \frac{g^2}{2} \frac{\upsilon^2}{\Lambda^2}
            \bar{t} \frac{\sigma^{\mu\nu}}{m_W}
                (d_V(k^2)+id_A(k^2)\gamma_5)t W^-_\mu W^+_\nu
   \end{equation} *)

    let anomalous_ttWW =
      if Flags.top_anom then
        [ ((M (U (-3)), O (Aux_top (2,0,0,true,TTWW)), M (U 3)), FBF (1, Psibar, TVA, Psi), G_TVA_ttWW);
          ((O (Aux_top (2,0,0,false,TTWW)), G Wm, G Wp), Aux_Gauge_Gauge 1, I_G_weak) ]
      else
        []

(* \begin{equation}
     \Delta\mathcal{L}_{bbWW} =
        - i \frac{g^2}{2} \frac{\upsilon^2}{\Lambda^2}
            \bar{b} \frac{\sigma^{\mu\nu}}{m_W}
                (d_V(k^2)+id_A(k^2)\gamma_5)b W^-_\mu W^+_\nu
   \end{equation} *)

    let anomalous_bbWW =
      if Flags.top_anom then
        [ ((M (D (-3)), O (Aux_top (2,0,0,true,BBWW)), M (D 3)), FBF (1, Psibar, TVA, Psi), G_TVA_bbWW);
          ((O (Aux_top (2,0,0,false,BBWW)), G Wm, G Wp), Aux_Gauge_Gauge 1, I_G_weak) ]
      else
        []

(* 4-fermion contact terms emerging from operator rewriting: *)

    let anomalous_top_qGuG_tt =
      [ ((M (U (-3)), O (Aux_top (1,1,0,true,QGUG)), M (U 3)), FBF (1, Psibar, VLR, Psi), G_VLR_qGuG) ]

    let anomalous_top_qGuG_ff n =
      List.map mom
        [ ((U (-n), Aux_top (1,1,0,false,QGUG), U n), FBF (1, Psibar, V, Psi), Unit);
          ((D (-n), Aux_top (1,1,0,false,QGUG), D n), FBF (1, Psibar, V, Psi), Unit) ]

    let anomalous_top_qGuG =
      if Flags.top_anom_4f then
        anomalous_top_qGuG_tt @ ThoList.flatmap anomalous_top_qGuG_ff [1;2;3]
      else
        []

    let anomalous_top_qBuB_tt =
      [ ((M (U (-3)), O (Aux_top (1,0,0,true,QBUB)), M (U 3)), FBF (1, Psibar, VLR, Psi), G_VLR_qBuB) ]

    let anomalous_top_qBuB_ff n =
      List.map mom
        [ ((U (-n), Aux_top (1,0,0,false,QBUB), U n), FBF (1, Psibar, VLR, Psi), G_VLR_qBuB_u);
          ((D (-n), Aux_top (1,0,0,false,QBUB), D n), FBF (1, Psibar, VLR, Psi), G_VLR_qBuB_d);
          ((L (-n), Aux_top (1,0,0,false,QBUB), L n), FBF (1, Psibar, VLR, Psi), G_VLR_qBuB_e);
          ((N (-n), Aux_top (1,0,0,false,QBUB), N n), FBF (1, Psibar, VL, Psi), G_VL_qBuB_n) ]

    let anomalous_top_qBuB =
      if Flags.top_anom_4f then
        anomalous_top_qBuB_tt @ ThoList.flatmap anomalous_top_qBuB_ff [1;2;3]
      else
        []

    let anomalous_top_qW_tq =
      [ ((M (U (-3)), O (Aux_top (1,0,0,true,QW)), M (U 3)), FBF (1, Psibar, VL, Psi), G_VL_qW);
        ((M (D (-3)), O (Aux_top (1,0,-1,true,QW)), M (U 3)), FBF (1, Psibar, VL, Psi), G_VL_qW);
        ((M (U (-3)), O (Aux_top (1,0,1,true,QW)), M (D 3)), FBF (1, Psibar, VL, Psi), G_VL_qW) ]

    let anomalous_top_qW_ff n =
      List.map mom
        [ ((U (-n), Aux_top (1,0,0,false,QW), U n), FBF (1, Psibar, VL, Psi), G_VL_qW_u);
          ((D (-n), Aux_top (1,0,0,false,QW), D n), FBF (1, Psibar, VL, Psi), G_VL_qW_d);
          ((N (-n), Aux_top (1,0,0,false,QW), N n), FBF (1, Psibar, VL, Psi), G_VL_qW_u);
          ((L (-n), Aux_top (1,0,0,false,QW), L n), FBF (1, Psibar, VL, Psi), G_VL_qW_d);
          ((D (-n), Aux_top (1,0,-1,false,QW), U n), FBF (1, Psibar, VL, Psi), Half);
          ((U (-n), Aux_top (1,0,1,false,QW), D n), FBF (1, Psibar, VL, Psi), Half);
          ((L (-n), Aux_top (1,0,-1,false,QW), N n), FBF (1, Psibar, VL, Psi), Half);
          ((N (-n), Aux_top (1,0,1,false,QW), L n), FBF (1, Psibar, VL, Psi), Half) ]

    let anomalous_top_qW =
      if Flags.top_anom_4f then
        anomalous_top_qW_tq @ ThoList.flatmap anomalous_top_qW_ff [1;2;3]
      else
        []

    let anomalous_top_DuDd =
      if Flags.top_anom_4f then
        [ ((M (U (-3)), O (Aux_top (0,0,0,true,DR)), M (U 3)), FBF (1, Psibar, SR, Psi), Half);
          ((M (U (-3)), O (Aux_top (0,0,0,false,DR)), M (U 3)), FBF (1, Psibar, SL, Psi), G_SL_DttR);
          ((M (D (-3)), O (Aux_top (0,0,0,false,DR)), M (D 3)), FBF (1, Psibar, SR, Psi), G_SR_DttR);
          ((M (U (-3)), O (Aux_top (0,0,0,true,DL)), M (U 3)), FBF (1, Psibar, SL, Psi), Half);
          ((M (D (-3)), O (Aux_top (0,0,0,false,DL)), M (D 3)), FBF (1, Psibar, SL, Psi), G_SL_DttL);
          ((M (D (-3)), O (Aux_top (0,0,-1,true,DR)), M (U 3)), FBF (1, Psibar, SR, Psi), Half);
          ((M (U (-3)), O (Aux_top (0,0,1,false,DR)), M (D 3)), FBF (1, Psibar, SLR, Psi), G_SLR_DbtR);
          ((M (D (-3)), O (Aux_top (0,0,-1,true,DL)), M (U 3)), FBF (1, Psibar, SL, Psi), Half);
          ((M (U (-3)), O (Aux_top (0,0,1,false,DL)), M (D 3)), FBF (1, Psibar, SL, Psi), G_SL_DbtL) ]
      else
        []

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @
       ThoList.flatmap color_currents [1;2;3] @
       ThoList.flatmap neutral_currents [1;2;3] @
       (if Flags.ckm_present then
         charged_currents_ckm
       else
         charged_currents_triv) @
       triple_gauge @
       goldstone_vertices @
       rsigma3 @ rsigma3t @ rphi3 @ rf3 @ rf3t @ rt3 @
       anomalous_ttA @ anomalous_bbA @
       anomalous_ttZ @ anomalous_bbZ @
       anomalous_tbW @ anomalous_tbWA @ anomalous_tbWZ @
       anomalous_ttWW @ anomalous_bbWW @
       anomalous_ttG @ anomalous_ttGG @
       anomalous_top_qGuG @ anomalous_top_qBuB @
       anomalous_top_qW @ anomalous_top_DuDd)

    let vertices4 =
      quartic_gauge @ rsigma4

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
      | "Rsigma" -> O Rsigma
      | "Rphi0" -> O Rphin
      | "Rphi+" -> O Rphip |  "Rphi-" -> O Rphim
      | "Rphi++" -> O Rphip |  "Rphi--" -> O Rphimm
      | "Rf" -> O Rf
      | "Rt0" -> O Rtn
      | "Rt+" -> O Rtp |  "Rt-" -> O Rtm
      | "Rt++" -> O Rtp |  "Rt--" -> O Rtmm
      | "Aux_t_ttGG0" -> O (Aux_top (2,1, 0,true,TTGG)) | "Aux_ttGG0" -> O (Aux_top (2,1, 0,false,TTGG))
      | "Aux_t_tbWA+" -> O (Aux_top (2,0, 1,true,TBWA)) | "Aux_tbWA+" -> O (Aux_top (2,0, 1,false,TBWA))
      | "Aux_t_tbWA-" -> O (Aux_top (2,0,-1,true,TBWA)) | "Aux_tbWA-" -> O (Aux_top (2,0,-1,false,TBWA))
      | "Aux_t_tbWZ+" -> O (Aux_top (2,0, 1,true,TBWZ)) | "Aux_tbWZ+" -> O (Aux_top (2,0, 1,false,TBWZ))
      | "Aux_t_tbWZ-" -> O (Aux_top (2,0,-1,true,TBWZ)) | "Aux_tbWZ-" -> O (Aux_top (2,0,-1,false,TBWZ))
      | "Aux_t_ttWW0" -> O (Aux_top (2,0, 0,true,TTWW)) | "Aux_ttWW0" -> O (Aux_top (2,0, 0,false,TTWW))
      | "Aux_t_bbWW0" -> O (Aux_top (2,0, 0,true,BBWW)) | "Aux_bbWW0" -> O (Aux_top (2,0, 0,false,BBWW))
      | "Aux_t_qGuG0" -> O (Aux_top (1,1, 0,true,QGUG)) | "Aux_qGuG0" -> O (Aux_top (1,1, 0,false,QGUG))
      | "Aux_t_qBuB0" -> O (Aux_top (1,0, 0,true,QBUB)) | "Aux_qBuB0" -> O (Aux_top (1,0, 0,false,QBUB))
      | "Aux_t_qW0"   -> O (Aux_top (1,0, 0,true,QW))   | "Aux_qW0"   -> O (Aux_top (1,0, 0,false,QW))
      | "Aux_t_qW+"   -> O (Aux_top (1,0, 1,true,QW))   | "Aux_qW+"   -> O (Aux_top (1,0, 1,false,QW))
      | "Aux_t_qW-"   -> O (Aux_top (1,0,-1,true,QW))   | "Aux_qW-"   -> O (Aux_top (1,0,-1,false,QW))
      | "Aux_t_dL0"   -> O (Aux_top (0,0, 0,true,DL))   | "Aux_dL0"   -> O (Aux_top (0,0, 0,false,DL))
      | "Aux_t_dL+"   -> O (Aux_top (0,0, 1,true,DL))   | "Aux_dL+"   -> O (Aux_top (0,0, 1,false,DL))
      | "Aux_t_dL-"   -> O (Aux_top (0,0,-1,true,DL))   | "Aux_dL-"   -> O (Aux_top (0,0,-1,false,DL))
      | "Aux_t_dR0"   -> O (Aux_top (0,0, 0,true,DR))   | "Aux_dR0"   -> O (Aux_top (0,0, 0,false,DR))
      | "Aux_t_dR+"   -> O (Aux_top (0,0, 1,true,DR))   | "Aux_dR+"   -> O (Aux_top (0,0, 1,false,DR))
      | "Aux_t_dR-"   -> O (Aux_top (0,0,-1,true,DR))   | "Aux_dR-"   -> O (Aux_top (0,0,-1,false,DR))
      | _ -> invalid_arg "Modellib_NoH.AltH.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib_NoH.AltH.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib_NoH.AltH.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib_NoH.AltH.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_NoH.AltH.flavor_to_string: invalid down type quark"
          end
      | G f ->
          begin match f with
          | Gl -> "gl"
          | Ga -> "A" | Z -> "Z"
          | Wp -> "W+" | Wm -> "W-"
          end
      | O f ->
          begin match f with
          | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
          | Rsigma -> "Rsigma"
          | Rphin -> "Rphin" | Rphip -> "Rphi+" | Rphim -> "Rphi-"
          | Rphipp -> "Rphi++" | Rphimm -> "Rphi--"
          | Rf -> "Rf"
          | Rtn -> "Rtn" | Rtp -> "Rt+" | Rtm -> "Rt-"
          | Rtpp -> "Rt++" | Rtmm -> "Rt--"
          | Aux_top (_,_,ch,n,v) -> "Aux_" ^ (if n then "t_" else "") ^ (
              begin match v with
              | TTGG -> "ttGG" | TBWA -> "tbWA" | TBWZ -> "tbWZ"
              | TTWW -> "ttWW" | BBWW -> "bbWW"
              | QGUG -> "qGuG" | QBUB -> "qBuB"
              | QW   -> "qW"   | DL   -> "dL"   | DR   -> "dR"
              end ) ^ ( if ch > 0 then "+" else if ch < 0 then "-" else "0" )
          end

    let flavor_to_TeX = function
      | M f ->
          begin match f with
          | L 1 -> "e^-" | L (-1) -> "e^+"
          | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
          | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
          | L _ -> invalid_arg
                "Modellib_NoH.AltH.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib_NoH.AltH.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib_NoH.AltH.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Modellib_NoH.AltH.flavor_to_TeX: invalid down type quark"
          end
      | G f ->
          begin match f with
          | Gl -> "g"
          | Ga -> "\\gamma" | Z -> "Z"
          | Wp -> "W^+" | Wm -> "W^-"
          end
      | O f ->
          begin match f with
          | Phip -> "\\phi^+" | Phim -> "\\phi^-" | Phi0 -> "\\phi^0" 
          | Rsigma -> "\\sigma"
          | Rphip -> "\\phi^+" | Rphim -> "\\phi^-" | Rphin -> "\\phi^0" 
          | Rphipp -> "\\phi^{++}" | Rphimm -> "\\phi^{--}"
          | Rf -> "f"
          | Rtp -> "t^+" | Rtm -> "t^-" | Rtn -> "t^0" 
          | Rtpp -> "t^{++}" | Rtmm -> "t^{--}"
          | Aux_top (_,_,ch,n,v) -> "\\textnormal{Aux_" ^ (if n then "t_" else "") ^ (
              begin match v with
              | TTGG -> "ttGG" | TBWA -> "tbWA" | TBWZ -> "tbWZ"
              | TTWW -> "ttWW" | BBWW -> "bbWW"
              | QGUG -> "qGuG" | QBUB -> "qBuB"
              | QW   -> "qW"   | DL   -> "dL"   | DR   -> "dR"
              end ) ^ ( if ch > 0 then "^+" else if ch < 0 then "^-" else "^0" ) ^ "}"
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
          | Rsigma -> "rsi"
          | Rphip -> "rpp" | Rphim -> "rpm" | Rphin -> "rpn"
          | Rphipp -> "rppp" | Rphimm -> "rpmm"
          | Rf -> "rf"
          | Rtp -> "rtp" | Rtm -> "rtm" | Rtn -> "rtn"
          | Rtpp -> "rtpp" | Rtmm -> "rtmm"
          | Aux_top (_,_,ch,n,v) -> "aux_" ^ (if n then "t_" else "") ^ (
              begin match v with
              | TTGG -> "ttgg" | TBWA -> "tbwa" | TBWZ -> "tbwz"
              | TTWW -> "ttww" | BBWW -> "bbww"
              | QGUG -> "qgug" | QBUB -> "qbub"
              | QW   -> "qw"   | DL   -> "dl"   | DR   -> "dr"
              end ) ^ "_" ^ ( if ch > 0 then "p" else if ch < 0 then "m" else "0" )
          end

(* Introducing new Resonances from 45, there are no PDG values *)

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
          | Rsigma -> 45
          | Rphin -> 46 | Rphip | Rphim -> 47 
          | Rphipp | Rphimm -> 48
          | Rf -> 52
          | Rtn -> 53 | Rtp | Rtm -> 54 
          | Rtpp | Rtmm -> 55
          | Aux_top (_,_,_,_,_) -> 81
          end

    let mass_symbol f = 
      "mass(" ^ string_of_int (abs (pdg f)) ^ ")"

    let width_symbol f =
      "width(" ^ string_of_int (abs (pdg f)) ^ ")"

    let constant_symbol = function
      | Unit -> "unit" | Half -> "half" | Pi -> "PI"
      | Alpha_QED -> "alpha" | E -> "e" | G_weak -> "g" | Vev -> "vev"
      | I_G_weak -> "ig" 
      | Sin2thw -> "sin2thw" | Sinthw -> "sinthw" | Costhw -> "costhw"
      | Q_lepton -> "qlep" | Q_up -> "qup" | Q_down -> "qdwn"
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_TVA_ttA -> "gtva_tta" | G_TVA_bbA -> "gtva_bba" 
      | G_VLR_ttZ -> "gvlr_ttz" | G_TVA_ttZ -> "gtva_ttz" | G_TVA_bbZ -> "gtva_bbz"
      | G_VLR_btW -> "gvlr_btw" | G_VLR_tbW -> "gvlr_tbw"
      | G_TLR_btW -> "gtlr_btw" | G_TRL_tbW -> "gtrl_tbw"
      | G_TLR_btWA -> "gtlr_btwa" | G_TRL_tbWA -> "gtrl_tbwa"
      | G_TLR_btWZ -> "gtlr_btwz" | G_TRL_tbWZ -> "gtrl_tbwz"
      | G_TVA_ttWW -> "gtva_ttww" | G_TVA_bbWW -> "gtva_bbww"
      | G_TVA_ttG -> "gtva_ttg" | G_TVA_ttGG -> "gtva_ttgg"
      | G_VLR_qGuG -> "gvlr_qgug"
      | G_VLR_qBuB -> "gvlr_qbub"
      | G_VLR_qBuB_u -> "gvlr_qbub_u" | G_VLR_qBuB_d -> "gvlr_qbub_d"
      | G_VLR_qBuB_e -> "gvlr_qbub_e" | G_VL_qBuB_n -> "gvl_qbub_n"
      | G_VL_qW -> "gvl_qw"
      | G_VL_qW_u -> "gvl_qw_u" | G_VL_qW_d -> "gvl_qw_d"
      | G_SL_DttR -> "gsl_dttr" | G_SR_DttR -> "gsr_dttr" | G_SL_DttL -> "gsl_dttl"
      | G_SLR_DbtR -> "gslr_dbtr" | G_SL_DbtL -> "gsl_dbtl"
      | G_CC -> "gcc"
      | G_CCQ (n1,n2) -> "gccq" ^ string_of_int n1 ^ string_of_int n2
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | I_G1_AWW -> "ig1a" | I_G1_ZWW -> "ig1z"
      | I_G1_plus_kappa_plus_G4_AWW -> "ig1pkpg4a"
      | I_G1_plus_kappa_plus_G4_ZWW -> "ig1pkpg4z"
      | I_G1_plus_kappa_minus_G4_AWW -> "ig1pkmg4a"
      | I_G1_plus_kappa_minus_G4_ZWW -> "ig1pkmg4z"
      | I_G1_minus_kappa_plus_G4_AWW -> "ig1mkpg4a"
      | I_G1_minus_kappa_plus_G4_ZWW -> "ig1mkpg4z"
      | I_G1_minus_kappa_minus_G4_AWW -> "ig1mkmg4a"
      | I_G1_minus_kappa_minus_G4_ZWW -> "ig1mkmg4z"
      | I_lambda_AWW -> "ila"
      | I_lambda_ZWW -> "ilz"
      | G5_AWW -> "rg5a"
      | G5_ZWW -> "rg5z"
      | I_kappa5_AWW -> "ik5a"
      | I_kappa5_ZWW -> "ik5z"
      | I_lambda5_AWW -> "il5a" | I_lambda5_ZWW -> "il5z"
      | Alpha_WWWW0 -> "alww0" | Alpha_WWWW2 -> "alww2"
      | Alpha_ZZWW0 -> "alzw0" | Alpha_ZZWW1 -> "alzw1"
      | Alpha_ZZZZ  -> "alzz"
      | D_Alpha_ZZWW0_S -> "dalzz0_s(gkm,mkm,"
      | D_Alpha_ZZWW0_T -> "dalzz0_t(gkm,mkm,"
      | D_Alpha_ZZWW1_S -> "dalzz1_s(gkm,mkm,"
      | D_Alpha_ZZWW1_T -> "dalzz1_t(gkm,mkm,"
      | D_Alpha_ZZWW1_U -> "dalzz1_u(gkm,mkm,"
      | D_Alpha_WWWW0_S -> "dalww0_s(gkm,mkm,"
      | D_Alpha_WWWW0_T -> "dalww0_t(gkm,mkm,"
      | D_Alpha_WWWW0_U -> "dalww0_u(gkm,mkm,"
      | D_Alpha_WWWW2_S -> "dalww2_s(gkm,mkm,"
      | D_Alpha_WWWW2_T -> "dalww2_t(gkm,mkm,"
      | D_Alpha_ZZZZ_S -> "dalz4_s(gkm,mkm,"
      | D_Alpha_ZZZZ_T -> "dalz4_t(gkm,mkm,"
      | G_SWW -> "gsww" | G_SZZ -> "gszz"
      | G_SWW_T -> "gswwt" | G_SZZ_T -> "gszzt"
      | G_PNWW -> "gpnww" | G_PNZZ -> "gpnzz"
      | G_PWZ -> "gpwz" | G_PWW -> "gpww"
      | G_FWW -> "gfww" | G_FZZ -> "gfzz"
      | G_FWW_T -> "gfwwt" | G_FZZ_T -> "gfzzt"
      | G_TNWW -> "gtnww" | G_TNZZ -> "gtnzz"
      | G_TWZ -> "gtwz" | G_TWW -> "gtww"
      | G_SSWW -> "gssww" | G_SSZZ -> "gsszz"
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f
      | K_Matrix_Coeff i -> "kc" ^ string_of_int i
      | K_Matrix_Pole i -> "kp" ^ string_of_int i

  end
