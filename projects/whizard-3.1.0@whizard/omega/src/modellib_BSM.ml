(* modellib_BSM.ml --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       cf. main AUTHORS file

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

(* \thocwmodulesection{Littlest Higgs Model} *)

module type BSM_flags =
  sig
    val u1_gauged      : bool
    val anom_ferm_ass  : bool
  end

module BSM_bsm : BSM_flags =
  struct
    let u1_gauged         = true
    let anom_ferm_ass     = false
  end

module BSM_ungauged : BSM_flags =
  struct
    let u1_gauged         = false
    let anom_ferm_ass     = false
  end

module BSM_anom : BSM_flags = 
  struct
    let u1_gauged         = false
    let anom_ferm_ass     = true
  end

module Littlest (Flags : BSM_flags) =
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
        "use complex mass scheme" ;
        "running_width", Arg.Unit (fun () -> default_width := Running),
        "use running width" ]
    let caveats () = []

    let gauge_symbol () =
      failwith "Modellib_BSM.Littlest.gauge_symbol: internal error"

    type matter_field = L of int | N of int | U of int | D of int 
        | TopH | TopHb
    type gauge_boson = Ga | Wp | Wm | Z | Gl | WHp | WHm 
        | ZH | AH
    type other = Phip | Phim | Phi0 | H | Eta | Psi0 
        | Psi1 | Psip | Psim | Psipp | Psimm

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
      failwith "Modellib_BSM.Littlest.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; N n; U n; D n ]

(* Since [Phi] already belongs to the EW Goldstone bosons we use [Psi]
   for the TeV scale complex triplet. *)

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Heavy Quarks", List.map matter_field [TopH; TopHb];
        "Heavy Scalars", List.map other 
          [Psi0; Psi1; Psip; Psim; Psipp; Psimm];
        "Gauge Bosons", List.map gauge_boson 
          (if Flags.u1_gauged then
           [Ga; Z; Wp; Wm; Gl; WHp; WHm; ZH; AH]
              else
           [Ga; Z; Wp; Wm; Gl; WHp; WHm; ZH]);
        "Higgs", List.map other
          (if Flags.u1_gauged then [H] 
              else [H; Eta]);
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
          | TopH -> Spinor | TopHb -> ConjSpinor 
          end
      | G f ->
          begin match f with
          | Ga | Gl -> Vector
          | Wp | Wm | Z | WHp | WHm | ZH | AH -> Massive_Vector
          end
      | O f ->
          begin match f with
          | Phip | Phim | Phi0 | H | Eta | Psi0 
          | Psi1 | Psip | Psim | Psipp | Psimm -> Scalar
          end

    let color = function 
      | M (U n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D n) -> Color.SUN  (if n > 0 then 3 else -3)
      | M TopH -> Color.SUN 3 | M TopHb -> Color.SUN (-3)
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
          | TopH -> Prop_Spinor | TopHb -> Prop_ConjSpinor 
          end
      | G f ->
          begin match f with
          | Ga | Gl -> Prop_Feynman
          | Wp | Wm | Z | WHp | WHm | ZH | AH -> Prop_Unitarity
          end
      | O f ->
          begin match f with
          | Phip | Phim | Phi0 -> Only_Insertion
          | H | Eta | Psi0 | Psi1 | Psip | Psim 
          | Psipp | Psimm -> Prop_Scalar
          end

(* Optionally, ask for the fudge factor treatment for the widths of
   charged particles.  Currently, this only applies to $W^\pm$ and top. *)

    let width f =
      if !use_fudged_width then
        match f with
        | G Wp | G Wm | M (U 3) | M (U (-3))
        | G WHp | G WHm | G ZH | G AH
        | M TopH | M TopHb -> Fudged
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
          | TopH -> TopHb | TopHb -> TopH
          end)
      | G f ->
          G (begin match f with
          | Gl -> Gl | Ga -> Ga | Z -> Z
          | Wp -> Wm | Wm -> Wp | WHm -> WHp
          | WHp -> WHm | ZH -> ZH | AH -> AH
          end)
      | O f ->
          O (begin match f with
          | Psi0 -> Psi0 | Psi1 -> Psi1 | Psip -> Psim
          | Psim -> Psip | Psipp -> Psimm | Psimm -> Psipp
          | Phip -> Phim | Phim -> Phip | Phi0 -> Phi0
          | H -> H | Eta -> Eta
          end)

    let fermion = function
      | M f ->
          begin match f with
          | L n -> if n > 0 then 1 else -1
          | N n -> if n > 0 then 1 else -1
          | U n -> if n > 0 then 1 else -1
          | D n -> if n > 0 then 1 else -1
          | TopH -> 1 | TopHb -> -1
          end
      | G f ->
          begin match f with
          | Gl | Ga | Z | Wp | Wm | WHp 
          | WHm | AH | ZH -> 0
          end
      | O f ->
          begin match f with
          | Psi0 | Psi1 | Psip | Psim | Psipp | Psimm 
          | Phip | Phim | Phi0 | H | Eta -> 0
          end

(* This model does NOT have a conserved generation charge
   even in absence of CKM mixing because of the heavy top
   admixture. *)

    module Ch = Charges.QQ
    let ( // ) = Algebra.Small_Rational.make

    let charge = function
      | M f ->
          begin match f with
          | L n -> if n > 0 then -1//1 else  1//1
          | N n -> 0//1
          | U n -> if n > 0 then  2//3 else -2//3
          | D n -> if n > 0 then -1//3 else  1//3
          | TopH -> 2//3
          | TopHb -> -2//3
          end
      | G f ->
          begin match f with
          | Gl | Ga | Z | AH | ZH -> 0//1
          | Wp | WHp ->  1//1
          | Wm | WHm -> -1//1
          end
      | O f ->
          begin match f with
          | H | Phi0 | Eta | Psi1 | Psi0 ->  0//1
          | Phip | Psip ->  1//1
          | Phim | Psim -> -1//1
          | Psipp -> 2//1
          | Psimm -> -2//1
          end

    let lepton = function
      | M f ->
          begin match f with
          | L n | N n -> if n > 0 then 1//1 else -1//1
          | U _ | D _ | _ -> 0//1
          end
      | G _ | O _ -> 0//1

    let baryon = function
      | M f ->
          begin match f with
          | L _ | N _ -> 0//1
          | U n | D n -> if n > 0 then 1//1 else -1//1
          | TopH -> 1//1
          | TopHb -> -1//1
          end
      | G _ | O _ -> 0//1

    let charges f = 
      [ charge f; lepton f; baryon f]

    type constant =
      | Unit | Pi | Alpha_QED | Sin2thw
      | Sinthw | Costhw | E | G_weak | Vev | VHeavy
      | Supp | Supp2
      | Sinpsi | Cospsi | Atpsi | Sccs  (* Mixing angles of SU(2) *)
      | Q_lepton | Q_up | Q_down | Q_Z_up | G_CC | G_CCtop
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down | G_NC_heavy
      | G_NC_h_neutrino | G_NC_h_lepton | G_NC_h_up | G_NC_h_down 
      | G_CC_heavy | G_ZHTHT | G_ZTHT | G_AHTHTH | G_AHTHT | G_AHTT
      | G_CC_WH | G_CC_W
      | I_Q_W | I_G_ZWW | I_G_WWW
      | I_G_AHWW | I_G_ZHWW | I_G_ZWHW | I_G_AHWHWH | I_G_ZHWHWH
      | I_G_AHWHW | I_Q_H  
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | G_WH4 | G_WHWHWW | G_WHWWW | G_WH3W 
      | G_WWAAH | G_WWAZH | G_WWZZH | G_WWZAH | G_WHWHAAH 
      | G_WHWHAZH | G_WHWHZZH | G_WHWHZAH | G_WWZHAH 
      | G_WHWHZHAH | G_WHWZZ | G_WHWAZ | G_WHWAAH | G_WHWZAH
      | G_WHWZHZH | G_WHWZHAH | G_WHWAZH | G_WHWZZH
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ
      | G_PsiWW | G_PsiWHW | G_PsiZZ | G_PsiZHZH 
      | G_PsiZHZ | G_PsiZAH | G_PsiZHAH | G_PsiAHAH
      | G_PsiZW | G_PsiZWH | G_PsiAHW | G_PsiAHWH 
      | G_PsiZHW | G_PsiZHWH
      | G_PsippWW | G_PsippWHW | G_PsippWHWH
      | G_PsiHW | G_PsiHWH | G_Psi0W | G_Psi0WH 
      | G_Psi1W | G_Psi1WH | G_PsiPPW | G_PsiPPWH
      | G_Psi1HAH | G_Psi01AH | G_AHPsip | G_Psi1HZ
      | G_Psi1HZH | G_Psi01Z | G_Psi01ZH | G_ZPsip | G_ZPsipp | G_ZHPsipp
      | G_HHAA | G_HHWHW | G_HHZHZ | G_HHAHZ | G_HHZHAH
      | G_HPsi0WW | G_HPsi0WHW | G_HPsi0ZZ
      | G_HPsi0ZHZH | G_HPsi0ZHZ | G_HPsi0AHAH | G_HPsi0ZAH | G_HPsi0ZHAH 
      | G_HPsipWA | G_HPsipWHA | G_HPsipWZ | G_HPsipWHZ | G_HPsipWAH
      | G_HPsipWHAH | G_HPsipWZH | G_HPsipWHZH | G_HPsippWW | G_HPsippWHWH
      | G_HPsippWHW | G_Psi00ZH | G_Psi00AH | G_Psi00ZHAH
      | G_Psi0pWA | G_Psi0pWHA | G_Psi0pWZ | G_Psi0pWHZ | G_Psi0pWAH
      | G_Psi0pWHAH | G_Psi0pWZH | G_Psi0pWHZH | G_Psi0ppWW | G_Psi0ppWHWH
      | G_Psi0ppWHW | I_G_Psi0pWA | I_G_Psi0pWHA | I_G_Psi0pWZ | I_G_Psi0pWHZ 
      | I_G_Psi0pWAH | I_G_Psi0pWHAH | I_G_Psi0pWZH | I_G_Psi0pWHZH 
      | I_G_Psi0ppWW | I_G_Psi0ppWHWH | I_G_Psi0ppWHW 
      | G_PsippZZ | G_PsippZHZH | G_PsippAZ | G_PsippAAH | G_PsippZAH
      | G_PsippWA | G_PsippWHA | G_PsippWZ | G_PsippWHZ | G_PsippWAH
      | G_PsippWHAH | G_PsippWZH | G_PsippWHZH
      | G_PsiccZZ | G_PsiccAZ | G_PsiccAAH | G_PsiccZZH | G_PsiccAZH
      | G_PsiccZAH
      | G_Htt | G_Hbb | G_Hcc | G_Htautau | G_H3 | G_H4
      | G_Hthth | G_Htht | G_Ethth | G_Etht | G_Ett
      | G_HHtt | G_HHthth | G_HHtht
      | G_Psi0tt | G_Psi0bb | G_Psi0cc | G_Psi0tautau
      | G_Psi1tt | G_Psi1bb | G_Psi1cc | G_Psi1tautau
      | G_Psipq3 | G_Psipq2 | G_Psipl3 | G_Psi0tth | G_Psi1tth
      | G_Psipbth | G_Ebb 
      | G_HGaGa | G_HGaZ | G_EGaGa | G_EGaZ | G_EGlGl
      | Gs | I_Gs | G2
      | G_HWHW | G_HWHWH | G_HAHAH | G_HZHZ | G_HZHAH | G_HAHZ
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let input_parameters =
      []

    let derived_parameters =
      []

    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

    let nc_coupling c t3 q =
      (Real_Array c,
       [Prod [g_over_2_costh; Diff (t3, Prod [Integer 2; q; Atom Sin2thw])];
        Prod [g_over_2_costh; t3]])

    let half = Quot (Integer 1, Integer 2)

    let derived_parameter_arrays =
      [ nc_coupling G_NC_neutrino half (Integer 0);
        nc_coupling G_NC_lepton (Neg half) (Integer (-1));
        nc_coupling G_NC_up half (Quot (Integer 2, Integer 3));
        nc_coupling G_NC_down (Neg half) (Quot (Integer (-1), Integer 3));
        nc_coupling G_NC_h_neutrino half (Integer 0);
        nc_coupling G_NC_h_lepton (Neg half) (Integer (-1));
        nc_coupling G_NC_h_up half (Quot (Integer 2, Integer 3));
        nc_coupling G_NC_h_down (Neg half) (Quot (Integer (-1), Integer 3)) ]

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

    let mgm ((m1, g, m2), fbf, c) = ((M m1, G g, M m2), fbf, c)
    let mhm ((m1, h, m2), fbf, c) = ((M m1, O h, M m2), fbf, c)
    let tgc ((g1, g2, g3), t, c) = ((G g1, G g2, G g3), t, c)
    let qgc ((g1, g2, g3, g4), t, c) = ((G g1, G g2, G g3, G g4), t, c)
    let hgg ((h, g1, g2), coup, c) = ((O h, G g1, G g2), coup, c)
    let ghh ((g, h1, h2), coup, c) = ((G g, O h1, O h2), coup, c)
    let hhgg ((h1, h2, g1, g2), coup, c) = ((O h1, O h2, G g1, G g2), coup, c)

    let electromagnetic_currents n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);  
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]

    let neutral_currents n =
      List.map mgm
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

(* The sign of this coupling is just the one of the T3, being -(1/2) for
   leptons and down quarks, and +(1/2) for neutrinos and up quarks. *)

    let neutral_heavy_currents n =
      List.map mgm
       ([ ((L (-n), ZH, L n), FBF ((-1), Psibar, VL, Psi), G_NC_heavy);
          ((N (-n), ZH, N n), FBF (1, Psibar, VL, Psi), G_NC_heavy);
          ((U (-n), ZH, U n), FBF (1, Psibar, VL, Psi), G_NC_heavy);
          ((D (-n), ZH, D n), FBF ((-1), Psibar, VL, Psi), G_NC_heavy)]
        @
          (if Flags.u1_gauged then
        [ ((L (-n), AH, L n), FBF (1, Psibar, VA, Psi), G_NC_h_lepton);
          ((N (-n), AH, N n), FBF (1, Psibar, VA, Psi), G_NC_h_neutrino);
          ((D (-n), AH, D n), FBF (1, Psibar, VA, Psi), G_NC_h_down)]
          else
            []))

    let color_currents n =
      List.map mgm
      [ ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs);
        ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs)]
	
   let heavy_top_currents = 
     List.map mgm
       ([ ((TopHb, Ga, TopH), FBF (1, Psibar, V, Psi), Q_up);
          ((TopHb, Z, TopH), FBF (1, Psibar, V, Psi), Q_Z_up);
	  ((TopHb, Gl, TopH), FBF (1, Psibar, V, Psi), Gs);
          ((TopHb, Z, U 3), FBF (1, Psibar, VL, Psi), G_ZTHT);
          ((U (-3), Z, TopH), FBF (1, Psibar, VL, Psi), G_ZTHT);
          ((TopHb, ZH, U 3), FBF (1, Psibar, VL, Psi), G_ZHTHT);
          ((U (-3), ZH, TopH), FBF (1, Psibar, VL, Psi), G_ZHTHT);
          ((U (-3), Wp, D 3), FBF (1, Psibar, VL, Psi), G_CCtop);
          ((D (-3), Wm, U 3), FBF (1, Psibar, VL, Psi), G_CCtop);
          ((TopHb, WHp, D 3), FBF (1, Psibar, VL, Psi), G_CC_WH);
          ((D (-3), WHm, TopH), FBF (1, Psibar, VL, Psi), G_CC_WH);
          ((TopHb, Wp, D 3), FBF (1, Psibar, VL, Psi), G_CC_W);
          ((D (-3), Wm, TopH), FBF (1, Psibar, VL, Psi), G_CC_W)] 
          @ 
            (if Flags.u1_gauged then
        [ ((U (-3), AH, U 3), FBF (1, Psibar, VA, Psi), G_AHTT);
          ((TopHb, AH, TopH), FBF (1, Psibar, VA, Psi), G_AHTHTH);
          ((TopHb, AH, U 3), FBF (1, Psibar, VR, Psi), G_AHTHT);
          ((U (-3), AH, TopH), FBF (1, Psibar, VR, Psi), G_AHTHT)]
            else
              []))


(* \begin{equation}
     \mathcal{L}_{\textrm{CC}} =
        - \frac{g}{2\sqrt2} \sum_i \bar\psi_i
               (T^+\fmslash{W}^+ + T^-\fmslash{W}^-)(1-\gamma_5)\psi_i 
   \end{equation} *)

    let charged_currents n =
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC);
          ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ]

    let charged_heavy_currents n =
      List.map mgm 
       ([ ((L (-n), WHm, N n), FBF (1, Psibar, VL, Psi), G_CC_heavy);
          ((N (-n), WHp, L n), FBF (1, Psibar, VL, Psi), G_CC_heavy);
          ((D (-n), WHm, U n), FBF (1, Psibar, VL, Psi), G_CC_heavy);
          ((U (-n), WHp, D n), FBF (1, Psibar, VL, Psi), G_CC_heavy)]
          @
            (if Flags.u1_gauged then
        [ ((U (-n), AH, U n), FBF (1, Psibar, VA, Psi), G_NC_h_up)]
            else
              []))
  

(* We specialize the third generation since there is an additional shift 
   coming from the admixture of the heavy top quark. The universal shift, 
   coming from the mixing in the non-Abelian gauge boson sector is 
   unobservable. (Redefinition of coupling constants by measured ones. *)

    let yukawa =
      List.map mhm
        [ ((U (-3), H, U 3), FBF (1, Psibar, S, Psi), G_Htt);
          ((D (-3), H, D 3), FBF (1, Psibar, S, Psi), G_Hbb);
          ((U (-2), H, U 2), FBF (1, Psibar, S, Psi), G_Hcc);
          ((L (-3), H, L 3), FBF (1, Psibar, S, Psi), G_Htautau)] 
        
    let yukawa_add' = 
      List.map mhm
        [ ((TopHb, H, TopH), FBF (1, Psibar, S, Psi), G_Hthth);
          ((TopHb, H, U 3), FBF (1, Psibar, SLR, Psi), G_Htht);
          ((U (-3), H, TopH), FBF (1, Psibar, SLR, Psi), G_Htht);
          ((U (-3), Psi0, U 3), FBF (1, Psibar, S, Psi), G_Psi0tt);
          ((D (-3), Psi0, D 3), FBF (1, Psibar, S, Psi), G_Psi0bb);
          ((U (-2), Psi0, U 2), FBF (1, Psibar, S, Psi), G_Psi0cc);
          ((L (-3), Psi0, L 3), FBF (1, Psibar, S, Psi), G_Psi0tautau);
          ((U (-3), Psi1, U 3), FBF (1, Psibar, P, Psi), G_Psi1tt);
          ((D (-3), Psi1, D 3), FBF (1, Psibar, P, Psi), G_Psi1bb);
          ((U (-2), Psi1, U 2), FBF (1, Psibar, P, Psi), G_Psi1cc);
          ((L (-3), Psi1, L 3), FBF (1, Psibar, P, Psi), G_Psi1tautau);
          ((U (-3), Psip, D 3), FBF (1, Psibar, SLR, Psi), G_Psipq3);
          ((U (-2), Psip, D 2), FBF (1, Psibar, SLR, Psi), G_Psipq2);
          ((N (-3), Psip, L 3), FBF (1, Psibar, SR, Psi), G_Psipl3);
          ((D (-3), Psim, U 3), FBF (1, Psibar, SLR, Psi), G_Psipq3);
          ((D (-2), Psim, U 2), FBF (1, Psibar, SLR, Psi), G_Psipq2);
          ((L (-3), Psim, N 3), FBF (1, Psibar, SL, Psi), G_Psipl3);
          ((TopHb, Psi0, U 3), FBF (1, Psibar, SL, Psi), G_Psi0tth);
          ((U (-3), Psi0, TopH), FBF (1, Psibar, SR, Psi), G_Psi0tth);
          ((TopHb, Psi1, U 3), FBF (1, Psibar, SL, Psi), G_Psi1tth);
          ((U (-3), Psi1, TopH), FBF (1, Psibar, SR, Psi), G_Psi1tth);
          ((TopHb, Psip, D 3), FBF (1, Psibar, SL, Psi), G_Psipbth);
          ((D (-3), Psim, TopH), FBF (1, Psibar, SR, Psi), G_Psipbth)]
 
    let yukawa_add = 
        if Flags.u1_gauged then
          yukawa_add'
        else
          yukawa_add' @
          List.map mhm
          [ ((U (-3), Eta, U 3), FBF (1, Psibar, P, Psi), G_Ett);
            ((TopHb, Eta, U 3), FBF (1, Psibar, SLR, Psi), G_Etht);
            ((D (-3), Eta, D 3), FBF (1, Psibar, P, Psi), G_Ebb);
            ((U (-3), Eta, TopH), FBF (1, Psibar, SLR, Psi), G_Etht)]
      
(* \begin{equation}
     \mathcal{L}_{\textrm{TGC}} =
        - e \partial_\mu A_\nu W_+^\mu W_-^\nu + \ldots
        - e \cot\theta_w  \partial_\mu Z_\nu W_+^\mu W_-^\nu + \ldots
   \end{equation} *)

    let standard_triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs) ]

    let heavy_triple_gauge =
      List.map tgc
       ([ ((Ga, WHm, WHp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, WHm, WHp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((ZH, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZHWW);
          ((Z, WHm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWHW);
          ((Z, Wm, WHp), Gauge_Gauge_Gauge (-1), I_G_ZWHW);
          ((ZH, WHm, Wp), Gauge_Gauge_Gauge 1, I_G_WWW);
          ((ZH, Wm, WHp), Gauge_Gauge_Gauge (-1), I_G_WWW);
          ((ZH, WHm, WHp), Gauge_Gauge_Gauge (-1), I_G_ZHWHWH)]          
          @
            (if Flags.u1_gauged then
        [ ((AH, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_AHWW);
          ((AH, WHm, Wp), Gauge_Gauge_Gauge 1, I_G_AHWHW);
          ((AH, Wm, WHp), Gauge_Gauge_Gauge (-1), I_G_AHWHW);
          ((AH, WHm, WHp), Gauge_Gauge_Gauge 1, I_G_AHWHWH)]
            else
              []))

    let triple_gauge =
        standard_triple_gauge @ heavy_triple_gauge

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let standard_quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          (Gl, Gl, Gl, Gl), gauge4, G2 ]

    let heavy_quartic_gauge = 
      List.map qgc
       ([ (WHm, Wp, WHm, Wp), gauge4, G_WWWW;
          (Wm, WHp, Wm, WHp), gauge4, G_WWWW;
          (WHm, WHp, WHm, WHp), gauge4, G_WH4;
          (Wm, Wp, WHm, WHp), gauge4, G_WHWHWW;        
          (Wm, Wp, Wm, WHp), gauge4, G_WHWWW;
          (Wm, Wp, WHm, Wp), gauge4, G_WHWWW;
          (WHm, WHp, Wm, WHp), gauge4, G_WH3W;
          (WHm, WHp, WHm, Wp), gauge4, G_WH3W;
          (WHm, Z, WHp, Z), minus_gauge4, G_ZZWW;
          (WHm, Z, WHp, Ga), minus_gauge4, G_AZWW;
          (WHm, Ga, WHp, ZH), minus_gauge4, G_AAWW;
          (WHm, Z, WHp, ZH), minus_gauge4, G_ZZWW;
          (Wm, ZH, Wp, ZH), minus_gauge4, G_WWWW;
          (Wm, Ga, Wp, ZH), minus_gauge4, G_WWAZH;
          (Wm, Z, Wp, ZH), minus_gauge4, G_WWZZH;
          (WHm, Ga, WHp, ZH), minus_gauge4, G_WHWHAZH;
          (WHm, Z, WHp, ZH), minus_gauge4, G_WHWHZZH;
          (WHm, ZH, WHp, ZH), minus_gauge4, G_WH4;
          (WHm, Z, Wp, Z), minus_gauge4, G_WHWZZ;
          (Wm, Z, WHp, Z), minus_gauge4, G_WHWZZ;
          (WHm, Ga, Wp, Z), minus_gauge4, G_WHWAZ;
          (Wm, Ga, WHp, Z), minus_gauge4, G_WHWAZ;
          (WHm, ZH, Wp, ZH), minus_gauge4, G_WHWZHZH;
          (Wm, ZH, WHp, ZH), minus_gauge4, G_WHWZHZH;
          (WHm, Ga, Wp, ZH), minus_gauge4, G_WHWAZH;
          (Wm, Ga, WHp, ZH), minus_gauge4, G_WHWAZH;
          (WHm, Z, Wp, ZH), minus_gauge4, G_WHWZZH;
          (Wm, Z, WHp, ZH), minus_gauge4, G_WHWZZH]
      @ 
        (if Flags.u1_gauged then
          [ (Wm, Ga, Wp, AH), minus_gauge4, G_WWAAH;            
            (Wm, Z, Wp, AH), minus_gauge4, G_WWZAH;
            (WHm, Ga, WHp, AH), minus_gauge4, G_WHWHAAH;
            (WHm, Z, WHp, AH), minus_gauge4, G_WHWHZAH;
            (Wm, ZH, Wp, AH), minus_gauge4, G_WWZHAH;
            (WHm, ZH, WHp, AH), minus_gauge4, G_WHWHZHAH;
            (WHm, Ga, Wp, AH), minus_gauge4, G_WHWAAH;
            (Wm, Ga, WHp, AH), minus_gauge4, G_WHWAAH;        
            (WHm, Z, Wp, AH), minus_gauge4, G_WHWZAH;
            (Wm, Z, WHp, AH), minus_gauge4, G_WHWZAH;
            (WHm, ZH, Wp, AH), minus_gauge4, G_WHWZHAH;
            (Wm, ZH, WHp, AH), minus_gauge4, G_WHWZHAH]
        else
          []))

    let quartic_gauge =
      standard_quartic_gauge @ heavy_quartic_gauge

    let standard_gauge_higgs' =
      List.map hgg
        [ ((H, Wp, Wm), Scalar_Vector_Vector 1, G_HWW);
          ((H, Z, Z), Scalar_Vector_Vector 1, G_HZZ) ]

    let heavy_gauge_higgs = 
      List.map hgg
       ([ ((H, Wp, WHm), Scalar_Vector_Vector 1, G_HWHW);
          ((H, WHp, Wm), Scalar_Vector_Vector 1, G_HWHW);
          ((H, WHp, WHm), Scalar_Vector_Vector 1, G_HWHWH);
          ((H, ZH, ZH), Scalar_Vector_Vector 1, G_HWHWH);
          ((H, ZH, Z), Scalar_Vector_Vector 1, G_HZHZ);
          ((H, Wp, Wm), Scalar_Vector_Vector 1, G_HZHAH)]
      @ 
        (if Flags.u1_gauged then
          [((H, AH, AH), Scalar_Vector_Vector 1, G_HAHAH);            
           ((H, Z, AH), Scalar_Vector_Vector 1, G_HAHZ)]
        else
          []))
          
    let triplet_gauge_higgs = 
      List.map hgg
       ([ ((Psi0, Wp, Wm), Scalar_Vector_Vector 1, G_PsiWW);
          ((Psi0, WHp, WHm), Scalar_Vector_Vector (-1), G_PsiWW);
          ((Psi0, WHp, Wm), Scalar_Vector_Vector 1, G_PsiWHW);
          ((Psi0, WHm, Wp), Scalar_Vector_Vector 1, G_PsiWHW);
          ((Psi0, Z, Z), Scalar_Vector_Vector 1, G_PsiZZ);        
          ((Psi0, ZH, ZH), Scalar_Vector_Vector 1, G_PsiZHZH);        
          ((Psi0, ZH, Z), Scalar_Vector_Vector 1, G_PsiZHZ);        
          ((Psim, Wp, Z), Scalar_Vector_Vector 1, G_PsiZW);
          ((Psip, Wm, Z), Scalar_Vector_Vector 1, G_PsiZW);
          ((Psim, WHp, Z), Scalar_Vector_Vector 1, G_PsiZWH);
          ((Psip, WHm, Z), Scalar_Vector_Vector 1, G_PsiZWH);
          ((Psim, Wp, ZH), Scalar_Vector_Vector 1, G_PsiZHW);
          ((Psip, Wm, ZH), Scalar_Vector_Vector 1, G_PsiZHW);
          ((Psim, WHp, ZH), Scalar_Vector_Vector 1, G_PsiZHWH);
          ((Psip, WHm, ZH), Scalar_Vector_Vector 1, G_PsiZHWH);
          ((Psimm, Wp, Wp), Scalar_Vector_Vector 1, G_PsippWW);
          ((Psipp, Wm, Wm), Scalar_Vector_Vector 1, G_PsippWW);
          ((Psimm, WHp, Wp), Scalar_Vector_Vector 1, G_PsippWHW);
          ((Psipp, WHm, Wm), Scalar_Vector_Vector 1, G_PsippWHW);
          ((Psimm, WHp, WHp), Scalar_Vector_Vector 1, G_PsippWHWH);
          ((Psipp, WHm, WHm), Scalar_Vector_Vector 1, G_PsippWHWH)] 
      @
        (if Flags.u1_gauged then
          [((Psi0, AH, Z), Scalar_Vector_Vector 1, G_PsiZAH);        
           ((Psi0, AH, ZH), Scalar_Vector_Vector 1, G_PsiZHAH);        
           ((Psi0, AH, AH), Scalar_Vector_Vector 1, G_PsiAHAH);
           ((Psim, Wp, AH), Scalar_Vector_Vector 1, G_PsiAHW);
           ((Psip, Wm, AH), Scalar_Vector_Vector 1, G_PsiAHW);
           ((Psim, WHp, AH), Scalar_Vector_Vector 1, G_PsiAHWH);
           ((Psip, WHm, AH), Scalar_Vector_Vector 1, G_PsiAHWH)]            
        else
          []))

    let triplet_gauge2_higgs = 
      List.map ghh
       ([ ((Wp, H, Psim), Vector_Scalar_Scalar 1, G_PsiHW);
          ((Wm, H, Psip), Vector_Scalar_Scalar 1, G_PsiHW);
          ((WHp, H, Psim), Vector_Scalar_Scalar 1, G_PsiHWH);
          ((WHm, H, Psip), Vector_Scalar_Scalar 1, G_PsiHWH);
          ((Wp, Psi0, Psim), Vector_Scalar_Scalar 1, G_Psi0W);
          ((Wm, Psi0, Psip), Vector_Scalar_Scalar 1, G_Psi0W);
          ((WHp, Psi0, Psim), Vector_Scalar_Scalar 1, G_Psi0WH);
          ((WHm, Psi0, Psip), Vector_Scalar_Scalar 1, G_Psi0WH);
          ((Wp, Psi1, Psim), Vector_Scalar_Scalar 1, G_Psi1W);
          ((Wm, Psi1, Psip), Vector_Scalar_Scalar (-1), G_Psi1W);
          ((WHp, Psi1, Psim), Vector_Scalar_Scalar 1, G_Psi1WH);
          ((WHm, Psi1, Psip), Vector_Scalar_Scalar (-1), G_Psi1WH);
          ((Wp, Psip, Psimm), Vector_Scalar_Scalar 1, G_PsiPPW);
          ((Wm, Psim, Psipp), Vector_Scalar_Scalar 1, G_PsiPPW);
          ((WHp, Psip, Psimm), Vector_Scalar_Scalar 1, G_PsiPPWH);
          ((WHm, Psim, Psipp), Vector_Scalar_Scalar 1, G_PsiPPWH);
          ((Ga, Psip, Psim), Vector_Scalar_Scalar 1, Q_lepton);
          ((Ga, Psipp, Psimm), Vector_Scalar_Scalar 2, Q_lepton);
          ((Z, H, Psi1), Vector_Scalar_Scalar 1, G_Psi1HZ);
          ((ZH, H, Psi1), Vector_Scalar_Scalar 1, G_Psi1HZH);
          ((Z, Psi0, Psi1), Vector_Scalar_Scalar 1, G_Psi01Z);
          ((ZH, Psi0, Psi1), Vector_Scalar_Scalar 1, G_Psi01ZH);
          ((Z, Psip, Psim), Vector_Scalar_Scalar 1, G_ZPsip);
          ((Z, Psipp, Psimm), Vector_Scalar_Scalar 2, G_ZPsipp);
          ((ZH, Psipp, Psimm), Vector_Scalar_Scalar 2, G_ZHPsipp)]        
        @ 
          (if Flags.u1_gauged then
            [((AH, H, Psi1), Vector_Scalar_Scalar 1, G_Psi1HAH);
             ((AH, Psi0, Psi1), Vector_Scalar_Scalar 1, G_Psi01AH); 
             ((AH, Psip, Psim), Vector_Scalar_Scalar 1, G_AHPsip);
             ((AH, Psipp, Psimm), Vector_Scalar_Scalar 2, G_AHPsip)]
          else []))
           
    let standard_gauge_higgs = 
      standard_gauge_higgs' @ heavy_gauge_higgs @ triplet_gauge_higgs @
      triplet_gauge2_higgs

    let standard_gauge_higgs4 =      
      List.map hhgg
      [ (H, H, Wp, Wm), Scalar2_Vector2 1, G_HHWW;
        (H, H, Z, Z), Scalar2_Vector2 1, G_HHZZ ]

    let littlest_gauge_higgs4 = 
      List.map hhgg
        ([ (H, H, WHp, WHm), Scalar2_Vector2 (-1), G_HHWW;
           (H, H, ZH, ZH), Scalar2_Vector2 (-1), G_HHWW;
           (H, H, Wp, WHm), Scalar2_Vector2 1, G_HHWHW;
           (H, H, WHp, Wm), Scalar2_Vector2 1, G_HHWHW;
           (H, H, ZH, Z), Scalar2_Vector2 (-1), G_HHZHZ;
           (H, Psi0, Wp, Wm), Scalar2_Vector2 1, G_HPsi0WW;
           (H, Psi0, WHp, WHm), Scalar2_Vector2 (-1), G_HPsi0WW;
           (H, Psi0, WHp, Wm), Scalar2_Vector2 1, G_HPsi0WHW;
           (H, Psi0, Wp, WHm), Scalar2_Vector2 1, G_HPsi0WHW;
           (H, Psi0, Z, Z), Scalar2_Vector2 1, G_HPsi0ZZ;
           (H, Psi0, ZH, ZH), Scalar2_Vector2 1, G_HPsi0ZHZH;
           (H, Psi0, ZH, Z), Scalar2_Vector2 1, G_HPsi0ZHZ;
           (H, Psim, Wp, Ga), Scalar2_Vector2 1, G_HPsipWA;
           (H, Psip, Wm, Ga), Scalar2_Vector2 1, G_HPsipWA;
           (H, Psim, WHp, Ga), Scalar2_Vector2 1, G_HPsipWHA;
           (H, Psip, WHm, Ga), Scalar2_Vector2 1, G_HPsipWHA;
           (H, Psim, Wp, Z), Scalar2_Vector2 1, G_HPsipWZ;
           (H, Psip, Wm, Z), Scalar2_Vector2 1, G_HPsipWZ;
           (H, Psim, WHp, Z), Scalar2_Vector2 1, G_HPsipWHZ;
           (H, Psip, WHm, Z), Scalar2_Vector2 1, G_HPsipWHZ;
           (H, Psim, Wp, ZH), Scalar2_Vector2 1, G_HPsipWZH;
           (H, Psip, Wm, ZH), Scalar2_Vector2 1, G_HPsipWZH;
           (H, Psim, WHp, ZH), Scalar2_Vector2 1, G_HPsipWHZH;
           (H, Psip, WHm, ZH), Scalar2_Vector2 1, G_HPsipWHZH;
           (H, Psimm, Wp, Wp), Scalar2_Vector2 1, G_HPsippWW;
           (H, Psipp, Wm, Wm), Scalar2_Vector2 1, G_HPsippWW;
           (H, Psimm, WHp, WHp), Scalar2_Vector2 1, G_HPsippWHWH;
           (H, Psipp, WHm, WHm), Scalar2_Vector2 1, G_HPsippWHWH;
           (H, Psimm, WHp, Wp), Scalar2_Vector2 1, G_HPsippWHW;
           (H, Psipp, WHm, Wm), Scalar2_Vector2 1, G_HPsippWHW;
           (Psi0, Psi0, Wp, Wm), Scalar2_Vector2 2, G_HHWW;
           (Psi0, Psi0, WHp, WHm), Scalar2_Vector2 (-2), G_HHWW;
           (Psi0, Psi0, Z, Z), Scalar2_Vector2 4, G_HHZZ;
           (Psi0, Psi0, ZH, ZH), Scalar2_Vector2 1, G_Psi00ZH;
           (Psi0, Psi0, WHp, Wm), Scalar2_Vector2 2, G_HHWHW;
           (Psi0, Psi0, Wp, WHm), Scalar2_Vector2 2, G_HHWHW;        
           (Psi0, Psi0, Z, ZH), Scalar2_Vector2 4, G_HHZHZ;
           (Psi0, Psim, Wp, Ga), Scalar2_Vector2 1, G_Psi0pWA;
           (Psi0, Psip, Wm, Ga), Scalar2_Vector2 1, G_Psi0pWA;
           (Psi0, Psim, WHp, Ga), Scalar2_Vector2 1, G_Psi0pWHA;
           (Psi0, Psip, WHm, Ga), Scalar2_Vector2 1, G_Psi0pWHA;
           (Psi0, Psim, Wp, Z), Scalar2_Vector2 1, G_Psi0pWZ;
           (Psi0, Psip, Wm, Z), Scalar2_Vector2 1, G_Psi0pWZ;
           (Psi0, Psim, WHp, Z), Scalar2_Vector2 1, G_Psi0pWHZ;
           (Psi0, Psip, WHm, Z), Scalar2_Vector2 1, G_Psi0pWHZ;
           (Psi0, Psim, Wp, ZH), Scalar2_Vector2 1, G_Psi0pWZH;
           (Psi0, Psip, Wm, ZH), Scalar2_Vector2 1, G_Psi0pWZH;
           (Psi0, Psim, WHp, ZH), Scalar2_Vector2 1, G_Psi0pWHZH;
           (Psi0, Psip, WHm, ZH), Scalar2_Vector2 1, G_Psi0pWHZH;
           (Psi0, Psimm, Wp, Wp), Scalar2_Vector2 1, G_Psi0ppWW;
           (Psi0, Psipp, Wm, Wm), Scalar2_Vector2 1, G_Psi0ppWW;
           (Psi0, Psimm, WHp, WHp), Scalar2_Vector2 1, G_Psi0ppWHWH;
           (Psi0, Psipp, WHm, WHm), Scalar2_Vector2 1, G_Psi0ppWHWH;
           (Psi0, Psimm, WHp, Wp), Scalar2_Vector2 1, G_Psi0ppWHW;
           (Psi0, Psipp, WHm, Wm), Scalar2_Vector2 1, G_Psi0ppWHW;
           (Psi1, Psi1, Wp, Wm), Scalar2_Vector2 2, G_HHWW;
           (Psi1, Psi1, WHp, WHm), Scalar2_Vector2 (-2), G_HHWW;
           (Psi1, Psi1, Z, Z), Scalar2_Vector2 4, G_HHZZ;
           (Psi1, Psi1, ZH, ZH), Scalar2_Vector2 1, G_Psi00ZH;
           (Psi1, Psi1, WHp, Wm), Scalar2_Vector2 2, G_HHWHW;
           (Psi1, Psi1, Wp, WHm), Scalar2_Vector2 2, G_HHWHW;       
           (Psi1, Psi1, Z, ZH), Scalar2_Vector2 4, G_HHZHZ;
           (Psi1, Psim, Wp, Ga), Scalar2_Vector2 1, I_G_Psi0pWA;
           (Psi1, Psip, Wm, Ga), Scalar2_Vector2 (-1), I_G_Psi0pWA;
           (Psi1, Psim, WHp, Ga), Scalar2_Vector2 1, I_G_Psi0pWHA;
           (Psi1, Psip, WHm, Ga), Scalar2_Vector2 (-1), I_G_Psi0pWHA;
           (Psi1, Psim, Wp, Z), Scalar2_Vector2 1, I_G_Psi0pWZ;
           (Psi1, Psip, Wm, Z), Scalar2_Vector2 (-1), I_G_Psi0pWZ;
           (Psi1, Psim, WHp, Z), Scalar2_Vector2 1, I_G_Psi0pWHZ;
           (Psi1, Psip, WHm, Z), Scalar2_Vector2 (-1), I_G_Psi0pWHZ;
           (Psi1, Psim, Wp, ZH), Scalar2_Vector2 1, I_G_Psi0pWZH;
           (Psi1, Psip, Wm, ZH), Scalar2_Vector2 (-1), I_G_Psi0pWZH;
           (Psi1, Psim, WHp, ZH), Scalar2_Vector2 1, I_G_Psi0pWHZH;
           (Psi1, Psip, WHm, ZH), Scalar2_Vector2 (-1), I_G_Psi0pWHZH;
           (Psi1, Psimm, Wp, Wp), Scalar2_Vector2 1, I_G_Psi0ppWW;
           (Psi1, Psipp, Wm, Wm), Scalar2_Vector2 (-1), I_G_Psi0ppWW;
           (Psi1, Psimm, WHp, WHp), Scalar2_Vector2 1, I_G_Psi0ppWHWH;
           (Psi1, Psipp, WHm, WHm), Scalar2_Vector2 (-1), I_G_Psi0ppWHWH;
           (Psi1, Psimm, WHp, Wp), Scalar2_Vector2 1, I_G_Psi0ppWHW;
           (Psi1, Psipp, WHm, Wm), Scalar2_Vector2 (-1), I_G_Psi0ppWHW;
           (Psip, Psim, Wp, Wm), Scalar2_Vector2 4, G_HHWW;
           (Psip, Psim, WHp, WHm), Scalar2_Vector2 1, G_Psi00ZH;
           (Psip, Psim, WHp, Wm), Scalar2_Vector2 4, G_HHWHW;
           (Psip, Psim, Wp, WHm), Scalar2_Vector2 4, G_HHWHW;        
           (Psip, Psim, Z, Z), Scalar2_Vector2 1, G_PsippZZ;
           (Psip, Psim, Ga, Ga), Scalar2_Vector2 2, G_AAWW;
           (Psip, Psim, ZH, ZH), Scalar2_Vector2 1, G_PsippZHZH;
           (Psip, Psim, Ga, Z), Scalar2_Vector2 4, G_PsippAZ;
           (Psip, Psimm, Wp, Ga), Scalar2_Vector2 1, G_PsippWA;
           (Psim, Psipp, Wm, Ga), Scalar2_Vector2 1, G_PsippWA;
           (Psip, Psimm, WHp, Ga), Scalar2_Vector2 1, G_PsippWHA;
           (Psim, Psipp, WHm, Ga), Scalar2_Vector2 1, G_PsippWHA;
           (Psip, Psimm, Wp, Z), Scalar2_Vector2 1, G_PsippWZ;
           (Psim, Psipp, Wm, Z), Scalar2_Vector2 1, G_PsippWZ;
           (Psip, Psimm, WHp, Z), Scalar2_Vector2 1, G_PsippWHZ;
           (Psim, Psipp, WHm, Z), Scalar2_Vector2 1, G_PsippWHZ;
           (Psip, Psimm, Wp, ZH), Scalar2_Vector2 1, G_PsippWZH;
           (Psim, Psipp, Wm, ZH), Scalar2_Vector2 1, G_PsippWZH;
           (Psip, Psimm, WHp, ZH), Scalar2_Vector2 1, G_PsippWHZH;
           (Psim, Psipp, WHm, ZH), Scalar2_Vector2 1, G_PsippWHZH;
           (Psipp, Psimm, Wp, Wm), Scalar2_Vector2 2, G_HHWW;
           (Psipp, Psimm, WHp, WHm), Scalar2_Vector2 (-2), G_HHWW;
           (Psipp, Psimm, WHp, Wm), Scalar2_Vector2 2, G_HHWHW;
           (Psipp, Psimm, Wp, WHm), Scalar2_Vector2 2, G_HHWHW;        
           (Psipp, Psimm, Z, Z), Scalar2_Vector2 1, G_PsiccZZ;
           (Psipp, Psimm, Ga, Ga), Scalar2_Vector2 8, G_AAWW;
           (Psipp, Psimm, ZH, ZH), Scalar2_Vector2 1, G_Psi00ZH;
           (Psipp, Psimm, Ga, Z), Scalar2_Vector2 1, G_PsiccAZ;
           (Psipp, Psimm, Z, ZH), Scalar2_Vector2 4, G_PsiccZZH;
           (Psipp, Psimm, Ga, ZH), Scalar2_Vector2 4, G_PsiccAZH]
         @
           (if Flags.u1_gauged then
             [(H, H, AH, AH), Scalar2_Vector2 1, G_HHAA;            
              (H, H, AH, Z), Scalar2_Vector2 (-1), G_HHAHZ;
              (H, H, ZH, AH), Scalar2_Vector2 (-1), G_HHZHAH;
              (H, Psi0, AH, AH), Scalar2_Vector2 1, G_HPsi0AHAH;
              (H, Psi0, Z, AH), Scalar2_Vector2 1, G_HPsi0ZAH;
              (H, Psi0, ZH, AH), Scalar2_Vector2 1, G_HPsi0ZHAH;
              (H, Psim, Wp, AH), Scalar2_Vector2 1, G_HPsipWAH;
              (H, Psip, Wm, AH), Scalar2_Vector2 1, G_HPsipWAH;
              (H, Psim, WHp, AH), Scalar2_Vector2 1, G_HPsipWHAH;
              (H, Psip, WHm, AH), Scalar2_Vector2 1, G_HPsipWHAH;
              (Psi0, Psi0, AH, AH), Scalar2_Vector2 1, G_Psi00AH;
              (Psi0, Psi0, Z, AH), Scalar2_Vector2 4, G_HHAHZ;
              (Psi0, Psi0, AH, ZH), Scalar2_Vector2 1, G_Psi00ZHAH;
              (Psi0, Psim, Wp, AH), Scalar2_Vector2 1, G_Psi0pWAH;
              (Psi0, Psip, Wm, AH), Scalar2_Vector2 1, G_Psi0pWAH;
              (Psi0, Psim, WHp, AH), Scalar2_Vector2 1, G_Psi0pWHAH;
              (Psi0, Psip, WHm, AH), Scalar2_Vector2 1, G_Psi0pWHAH;
              (Psi1, Psi1, AH, AH), Scalar2_Vector2 1, G_Psi00AH;
              (Psi1, Psi1, Z, AH), Scalar2_Vector2 4, G_HHAHZ;
              (Psi1, Psi1, AH, ZH), Scalar2_Vector2 1, G_Psi00ZHAH;
              (Psi1, Psim, Wp, AH), Scalar2_Vector2 1, I_G_Psi0pWAH;
              (Psi1, Psip, Wm, AH), Scalar2_Vector2 (-1), I_G_Psi0pWAH;
              (Psi1, Psim, WHp, AH), Scalar2_Vector2 1, I_G_Psi0pWHAH;
              (Psi1, Psip, WHm, AH), Scalar2_Vector2 (-1), I_G_Psi0pWHAH;
              (Psip, Psim, AH, AH), Scalar2_Vector2 1, G_Psi00AH;
              (Psip, Psim, Ga, AH), Scalar2_Vector2 4, G_PsippAAH;
              (Psip, Psim, Z, AH), Scalar2_Vector2 4, G_PsippZAH;
              (Psip, Psimm, Wp, AH), Scalar2_Vector2 1, G_PsippWAH;
              (Psim, Psipp, Wm, AH), Scalar2_Vector2 1, G_PsippWAH;
              (Psip, Psimm, WHp, AH), Scalar2_Vector2 1, G_PsippWHAH;
              (Psim, Psipp, WHm, AH), Scalar2_Vector2 1, G_PsippWHAH;
              (Psipp, Psimm, AH, AH), Scalar2_Vector2 1, G_Psi00AH;
              (Psipp, Psimm, AH, ZH), Scalar2_Vector2 (-1), G_Psi00ZHAH;
              (Psipp, Psimm, Ga, AH), Scalar2_Vector2 4, G_PsiccAAH;
              (Psipp, Psimm, Z, AH), Scalar2_Vector2 4, G_PsiccZAH]
           else []))

    let standard_higgs =
      [ (O H, O H, O H), Scalar_Scalar_Scalar 1, G_H3 ]
        
   let anomaly_higgs = 
     List.map hgg
      [ (Eta, Gl, Gl), Dim5_Scalar_Gauge2_Skew 1, G_EGlGl;
        (Eta, Ga, Ga), Dim5_Scalar_Gauge2_Skew 1, G_EGaGa; 
        (Eta, Ga, Z), Dim5_Scalar_Gauge2_Skew 1, G_EGaZ] 
(*    @ [ (H, Ga, Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
        (H, Ga, Z), Dim5_Scalar_Gauge2 1, G_HGaZ ]           *)

    let standard_higgs4 =
      [ (O H, O H, O H, O H), Scalar4 1, G_H4 ]

    let gauge_higgs =
        standard_gauge_higgs

    let gauge_higgs4 =
        standard_gauge_higgs4

    let higgs =
        standard_higgs

    let higgs4 =
        standard_higgs4

    let top_quartic = 
      [ ((M (U (-3)), O H, O H, M (U 3)), GBBG (1, Psibar, S2, Psi), G_HHtt);
   ((M (TopHb), O H, O H, M TopH), GBBG (1, Psibar, S2, Psi), G_HHthth);
   ((M (U (-3)), O H, O H, M TopH), GBBG (1, Psibar, S2LR, Psi), G_HHtht);
   ((M (TopHb), O H, O H, M (U 3)), GBBG (1, Psibar, S2LR, Psi), G_HHtht)]

    let goldstone_vertices =
      List.map hgg
      [ ((Phi0, Wm, Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((Phip, Ga, Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((Phip, Z, Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((Phim, Wp, Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((Phim, Wp, Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @ 
       ThoList.flatmap color_currents [1;2;3] @
       ThoList.flatmap neutral_currents [1;2;3] @
       ThoList.flatmap neutral_heavy_currents [1;2;3] @       
       ThoList.flatmap charged_currents [1;2;3] @
       ThoList.flatmap charged_heavy_currents [1;2;3] @
       heavy_top_currents @ 
       (if Flags.u1_gauged then []
           else anomaly_higgs) @
       yukawa @ yukawa_add @ triple_gauge @ 
       gauge_higgs @ higgs @ goldstone_vertices)

    let vertices4 =
      quartic_gauge @ gauge_higgs4 @ higgs4 @ top_quartic

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
      | "th" -> M TopH  | "thbar" -> M TopHb
      | "g" | "gl" -> G Gl
      | "A" -> G Ga | "Z" | "Z0" -> G Z
      | "AH" | "AH0" | "Ah" | "Ah0" -> G AH
      | "ZH" | "ZH0" | "Zh" | "Zh0" -> G ZH
      | "W+" -> G Wp | "W-" -> G Wm
      | "WH+" -> G WHp | "WH-" -> G WHm
      | "H" | "h" -> O H | "eta" | "Eta" -> O Eta
      | "Psi" | "Psi0" | "psi" | "psi0" -> O Psi0
      | "Psi1" | "psi1" -> O Psi1
      | "Psi+" | "psi+" | "Psip" | "psip" -> O Psip
      | "Psi-" | "psi-" | "Psim" | "psim" -> O Psim
      | "Psi++" | "psi++" | "Psipp" | "psipp" -> O Psipp
      | "Psi--" | "psi--" | "Psimm" | "psimm" -> O Psimm
      | _ -> invalid_arg "Modellib_BSM.Littlest.flavor_of_string" 


    let flavor_to_string = function
      | M f -> 
          begin match f with 
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg "Modellib_BSM.Littlest.flavor_to_string" 
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg "Modellib_BSM.Littlest.flavor_to_string" 
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg "Modellib_BSM.Littlest.flavor_to_string" 
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg "Modellib_BSM.Littlest.flavor_to_string" 
          | TopH -> "th" | TopHb -> "thbar"
          end
      | G f -> 
          begin match f with
          | Gl -> "g"
          | Ga -> "A" | Z -> "Z"
          | Wp -> "W+" | Wm -> "W-"
          | ZH -> "ZH" | AH -> "AH" | WHp -> "WHp" | WHm -> "WHm"
          end
      | O f ->
          begin match f with
          | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
          | H -> "H" | Eta -> "Eta"
          | Psi0 -> "Psi0" | Psi1 -> "Psi1" | Psip -> "Psi+" 
          | Psim -> "Psi-" | Psipp -> "Psi++" | Psimm -> "Psi--"
          end                

    let flavor_to_TeX = function
      | M f -> 
          begin match f with 
          | L 1 -> "e^-" | L (-1) -> "e^+"
          | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
          | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
          | L _ -> invalid_arg "Modellib_BSM.Littlest.flavor_to_TeX" 
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg "Modellib_BSM.Littlest.flavor_to_TeX" 
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg "Modellib_BSM.Littlest.flavor_to_TeX" 
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg "Modellib_BSM.Littlest.flavor_to_TeX" 
          | TopH -> "T" | TopHb -> "\\bar{T}"
          end
      | G f -> 
          begin match f with
          | Gl -> "g"
          | Ga -> "\\gamma" | Z -> "Z"
          | Wp -> "W^+" | Wm -> "W^-"
          | ZH -> "Z_H" | AH -> "\\gamma_H" | WHp -> "W_H^+" | WHm -> "W_H^-"
          end
      | O f ->
          begin match f with
          | Phip -> "\\Phi^+" | Phim -> "\\Phi^-" | Phi0 -> "\\Phi^0" 
          | H -> "H" | Eta -> "\\eta"
          | Psi0 -> "\\Psi_S" | Psi1 -> "\\Psi_P" | Psip -> "\\Psi^+" 
          | Psim -> "\\Psi^-" | Psipp -> "\\Psi^{++}" | Psimm -> "\\Psi^{--}"
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
          | TopH -> "th" | TopHb -> "thb" 
          end
      | G f -> 
          begin match f with
          | Gl -> "gl"
          | Ga -> "a" | Z -> "z"
          | Wp -> "wp" | Wm -> "wm"
          | ZH -> "zh" | AH -> "ah" | WHp -> "whp" | WHm -> "whm"
          end
      | O f ->
          begin match f with
          | Phip -> "pp" | Phim -> "pm" | Phi0 -> "p0" 
          | H -> "h" | Eta -> "eta"
          | Psi0 -> "psi0" | Psi1 -> "psi1" | Psip -> "psip" 
          | Psim -> "psim" | Psipp -> "psipp" | Psimm -> "psimm"
          end

(* There are PDG numbers for Z', Z'', W', 32-34, respectively.
   We just introduce a number 38 for Y0 as a Z'''.
   As well, there is the number 8 for a t'. But we cheat a little bit and 
   take the number 35 which is reserved for a heavy scalar Higgs for the 
   Eta scalar.
   For the heavy Higgs states we take 35 and 36 for the neutral ones, 37 for 
   the charged and 38 for the doubly-charged. 
   The pseudoscalar gets the 39.
*)  

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
          | TopH -> 8 | TopHb -> (-8)
          end
      | G f ->
          begin match f with
          | Gl -> 21
          | Ga -> 22 | Z -> 23
          | Wp -> 24 | Wm -> (-24)
          | AH -> 32 | ZH -> 33 | WHp -> 34 | WHm -> (-34) 
          end
      | O f ->
          begin match f with
          | Phip | Phim -> 27 | Phi0 -> 26
          | Psi0 -> 35 | Psi1 -> 36 | Psip -> 37 | Psim -> (-37)
          | Psipp -> 38 | Psimm -> (-38)
          | H -> 25 | Eta -> 39
          end

    let mass_symbol f = 
      "mass(" ^ string_of_int (abs (pdg f)) ^ ")"

    let width_symbol f =
      "width(" ^ string_of_int (abs (pdg f)) ^ ")"

    let constant_symbol = function
      | Unit -> "unit" | Pi -> "PI" | VHeavy -> "vheavy"
      | Alpha_QED -> "alpha" | E -> "e" | G_weak -> "g" | Vev -> "vev"
      | Sin2thw -> "sin2thw" | Sinthw -> "sinthw" | Costhw -> "costhw"
      | Sinpsi -> "sinpsi" | Cospsi -> "cospsi" 
      | Atpsi -> "atpsi" | Sccs -> "sccs"
      | Supp -> "vF" | Supp2 -> "v2F2" 
      | Q_lepton -> "qlep" | Q_up -> "qup" | Q_down -> "qdwn"
      | Q_Z_up -> "qzup" 
      | G_ZHTHT -> "gzhtht" | G_ZTHT -> "gztht"
      | G_AHTHTH -> "gahthth" | G_AHTHT -> "gahtht" | G_AHTT -> "gahtt"
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_CC -> "gcc" | G_CCtop -> "gcctop" | G_CC_heavy -> "gcch" 
      | G_CC_WH -> "gccwh" | G_CC_W -> "gccw" 
      | G_NC_h_lepton -> "gnchlep" | G_NC_h_neutrino -> "gnchneu"
      | G_NC_h_up -> "gnchup" | G_NC_h_down -> "gnchdwn"   
      | G_NC_heavy -> "gnch"  
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" | I_G_WWW -> "igwww"
      | I_G_AHWW -> "igahww" | I_G_ZHWW -> "igzhww" | I_G_ZWHW -> "igzwhw"
      | I_G_AHWHWH -> "igahwhwh" | I_G_ZHWHWH -> "igzhwhwh"
      | I_G_AHWHW -> "igahwhw"
      | I_Q_H -> "iqh" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_WH4 -> "gwh4" | G_WHWHWW -> "gwhwhww" | G_WHWWW -> "gwhwww" 
      | G_WH3W -> "gwh3w"
      | G_WWAAH -> "gwwaah" | G_WWAZH -> "gwwazh" | G_WWZZH -> "gwwzzh" 
      | G_WWZAH -> "gwwzah" | G_WHWHAAH -> "gwhwhaah"  
      | G_WHWHAZH -> "gwhwhazh" | G_WHWHZZH -> "gwhwhzzh" 
      | G_WHWHZAH -> "gwhwhzah" 
      | G_WWZHAH -> "gwwzhah" | G_WHWHZHAH -> "gwhwhzhah"
      | G_WHWZZ -> "gwhwzz" | G_WHWAZ -> "gwhwaz" 
      | G_WHWAAH -> "gwhwaah" | G_WHWZAH -> "gwhwzah"
      | G_WHWZHZH -> "gwhwzhzh" | G_WHWZHAH -> "gwhwzhah"
      | G_WHWAZH -> "gwhwazh" | G_WHWZZH -> "gwhwzzh"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_HWHW -> "ghwhw" | G_HWHWH -> "ghwhwh" | G_HAHAH -> "ghahah"
      | G_HZHZ -> "ghzhz" | G_HZHAH -> "ghzhah"
      | G_HAHZ -> "ghahz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc"
      | G_Hthth -> "ghthth" | G_Htht -> "ghtht"
      | G_HHtt -> "ghhtt" | G_HHthth -> "ghhthth" | G_HHtht -> "ghhtht"
      | G_Psi0tt -> "gpsi0tt" | G_Psi0bb -> "gpsi0bb" 
      | G_Psi0cc -> "gpsi0cc" | G_Psi0tautau -> "gpsi0tautau"
      | G_Psi1tt -> "gpsi1tt" | G_Psi1bb -> "gpsi1bb" 
      | G_Psi1cc -> "gpsi1cc" | G_Psi1tautau -> "gpsi1tautau"
      | G_Psipq3 -> "gpsipq3" | G_Psipq2 -> "gpsipq2" | G_Psipl3 -> "gpsipl3"
      | G_Psi0tth -> "gpsi0tth" | G_Psi1tth -> "gpsi1tth"
      | G_Psipbth -> "gpsipbth"
      | G_Ethth -> "gethth" | G_Etht -> "getht"
      | G_Ett -> "gett" | G_Ebb -> "gebb"
      | G_HGaGa -> "ghgaga" | G_HGaZ -> "ghgaz"
      | G_EGaGa -> "geaa" | G_EGaZ -> "geaz" | G_EGlGl -> "gegg"
      | G_H3 -> "gh3" | G_H4 -> "gh4"
      | G_PsiWW -> "gpsiww" | G_PsiWHW -> "gpsiwhw" 
      | G_PsiZZ -> "gpsizz" | G_PsiZHZH -> "gpsizhzh" 
      | G_PsiZHZ -> "gpsizhz" | G_PsiZAH -> "gpsizah" 
      | G_PsiZHAH -> "gpsizhah" | G_PsiAHAH -> "gpsiahah"
      | G_PsiZW -> "gpsizw" | G_PsiZWH -> "gpsizwh" | G_PsiAHW -> "gpsiahw"
      | G_PsiAHWH -> "gpsiahwh" | G_PsiZHW -> "gpsizhw" 
      | G_PsiZHWH -> "gpsizhwh"
      | G_PsippWW -> "gpsippww" | G_PsippWHW -> "gpsippwhw" 
      | G_PsippWHWH -> "gpsippwhwh"
      | Gs -> "gs" | G2 -> "gs**2" | I_Gs -> "igs"
      | G_PsiHW -> "gpsihw" | G_PsiHWH -> "gpsihwh" 
      | G_Psi0W -> "gpsi0w" | G_Psi0WH -> "gpsi0wh"
      | G_Psi1W -> "gpsi1w" | G_Psi1WH -> "gpsi1wh" 
      | G_PsiPPW -> "gpsippw" | G_PsiPPWH -> "gpsippwh"
      | G_Psi1HAH -> "gpsihah" | G_Psi01AH -> "gpsi0ah" 
      | G_AHPsip -> "gahpsip" | G_Psi1HZ -> "gpsi1hz"
      | G_Psi1HZH -> "gpsi1hzh" | G_Psi01Z -> "gpsi01z" 
      | G_Psi01ZH -> "gpsi01zh" | G_ZPsip -> "gzpsip" 
      | G_ZPsipp -> "gzpsipp" | G_ZHPsipp -> "gzhpsipp"
      | G_HHAA -> "ghhaa" | G_HHWHW -> "ghhwhw" | G_HHZHZ -> "ghhzhz"
      | G_HHAHZ -> "ghhahz" | G_HHZHAH -> "ghhzhah"
      | G_HPsi0WW -> "ghpsi0ww" | G_HPsi0WHW -> "ghpsi0whw" 
      | G_HPsi0ZZ -> "ghpsi0zz" | G_HPsi0ZHZH -> "ghpsi0zhzh" 
      | G_HPsi0ZHZ -> "ghpsi0zhz" | G_HPsi0AHAH -> "ghpsi0ahah"
      | G_HPsi0ZAH -> "ghpsi0zah" | G_HPsi0ZHAH -> "ghpsi0zhah"
      | G_HPsipWA -> "ghpsipwa" | G_HPsipWHA -> "ghpsipwha" 
      | G_HPsipWZ -> "ghpsipwz" | G_HPsipWHZ -> "ghpsiwhz" 
      | G_HPsipWAH -> "ghpsipwah" | G_HPsipWHAH -> "ghpsipwhah" 
      | G_HPsipWZH -> "ghpsipwzh" | G_HPsipWHZH -> "ghpsipwhzh" 
      | G_HPsippWW -> "ghpsippww" | G_HPsippWHWH -> "ghpsippwhwh"
      | G_HPsippWHW -> "ghpsippwhw" | G_Psi00ZH -> "gpsi00zh" 
      | G_Psi00AH -> "gpsi00ah" | G_Psi00ZHAH -> "gpsi00zhah"
      | G_Psi0pWA -> "gpsi0pwa" | G_Psi0pWHA -> "gpsi0pwha" 
      | G_Psi0pWZ -> "gpsi0pwz" | G_Psi0pWHZ -> "gpsi0pwhz" 
      | G_Psi0pWAH -> "gpsi0pwah" | G_Psi0pWHAH -> "gpsi0pwhah" 
      | G_Psi0pWZH -> "gpsi0pwzh" | G_Psi0pWHZH -> "gpsi0pwhzh" 
      | G_Psi0ppWW -> "gpsi0ppww" | G_Psi0ppWHWH -> "gpsi0ppwhwh"
      | G_Psi0ppWHW -> "gpsi0ppwhw"
      | I_G_Psi0pWA -> "i_gpsi0pwa" | I_G_Psi0pWHA -> "i_gpsi0pwha" 
      | I_G_Psi0pWZ -> "i_gpsi0pwz" | I_G_Psi0pWHZ -> "i_gpsi0pwhz" 
      | I_G_Psi0pWAH -> "i_gpsi0pwah" | I_G_Psi0pWHAH -> "i_gpsi0pwhah" 
      | I_G_Psi0pWZH -> "i_gpsi0pwzh" | I_G_Psi0pWHZH -> "i_gpsi0pwhzh" 
      | I_G_Psi0ppWW -> "i_gpsi0ppww" | I_G_Psi0ppWHWH -> "i_gpsi0ppwhwh"
      | I_G_Psi0ppWHW -> "i_gpsi0ppwhw" 
      | G_PsippZZ -> "gpsippzz" | G_PsippZHZH -> "gpsippzhzh"
      | G_PsippAZ -> "gpsippaz" | G_PsippAAH -> "gpsippaah" 
      | G_PsippZAH -> "gpsippzah"
      | G_PsippWA -> "gpsippwa" | G_PsippWHA -> "gpsippwha" 
      | G_PsippWZ -> "gpsippwz" | G_PsippWHZ -> "gpsippwhz" 
      | G_PsippWAH -> "gpsippwah" | G_PsippWHAH -> "gpsippwhah"
      | G_PsippWZH -> "gpsippwzh" | G_PsippWHZH -> "gpsippwhzh"
      | G_PsiccZZ -> "gpsicczz" | G_PsiccAZ -> "gpsiccaz"
      | G_PsiccAAH -> "gpsiccaah" | G_PsiccZZH -> "gpsicczzh" 
      | G_PsiccAZH -> "gpsiccazh" | G_PsiccZAH -> "gpsicczah"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f

  end

module Littlest_Tpar (Flags : BSM_flags) =
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

    type flavor = L of int | N of int | U of int | D of int 
        | Topp | Toppb 
        | Ga | Wp | Wm | Z | Gl | Lodd of int | Nodd of int 
        | Uodd of int | Dodd of int
        | WHp | WHm | ZH | AH | Phip | Phim | Phi0 | H | Eta | Psi0 
        | Psi1 | Psip | Psim | Psipp | Psimm

    type gauge = unit

    let gauge_symbol () =
      failwith "Modellib_BSM.Littlest_Tpar.gauge_symbol: internal error"

    let family n = [ L n; N n; U n; D n; Dodd n; Nodd n; Lodd n; Uodd n ]

(* Since [Phi] already belongs to the EW Goldstone bosons we use [Psi]
   for the TeV scale complex triplet. 

   We use the notation Todd1 = Uodd 3, Todd2 = Uodd 4.
*)

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Heavy Quarks", [Topp; Toppb; Uodd 4; Uodd (-4)];
        "Heavy Scalars", [Psi0; Psi1; Psip; Psim; Psipp; Psimm];
        "Gauge Bosons", if Flags.u1_gauged then
          [Ga; Z; Wp; Wm; Gl; WHp; WHm; ZH; AH]
            else
          [Ga; Z; Wp; Wm; Gl; WHp; WHm; ZH];
        "Higgs", if Flags.u1_gauged then [H] 
        else [H; Eta];
        "Goldstone Bosons", [Phip; Phim; Phi0] ]

    let flavors () = ThoList.flatmap snd (external_flavors ()) 

    let spinor n =
      if n >= 0 then
        Spinor
      else
        ConjSpinor

    let lorentz = function
      | L n -> spinor n | N n -> spinor n
      | U n -> spinor n | D n -> spinor n
      | Topp -> Spinor | Toppb -> ConjSpinor 
      | Ga | Gl -> Vector
      | Wp | Wm | Z | WHp | WHm | ZH | AH -> Massive_Vector
      | _ -> Scalar

    let color = function 
      | U n -> Color.SUN (if n > 0 then 3 else -3)
      | Uodd n -> Color.SUN (if n > 0 then 3 else -3)
      | D n -> Color.SUN  (if n > 0 then 3 else -3)
      | Dodd n -> Color.SUN (if n > 0 then 3 else -3)
      | Topp -> Color.SUN 3 | Toppb -> Color.SUN (-3)
      | Gl -> Color.AdjSUN 3
      | _ -> Color.Singlet

    let nc () = 3

    let prop_spinor n =
      if n >= 0 then
        Prop_Spinor
      else
        Prop_ConjSpinor

    let propagator = function
      | L n -> prop_spinor n | N n -> prop_spinor n
      | Lodd n -> prop_spinor n | Nodd n -> prop_spinor n
      | U n -> prop_spinor n | D n -> prop_spinor n
      | Uodd n -> prop_spinor n | Dodd n -> prop_spinor n
      | Topp -> Prop_Spinor | Toppb -> Prop_ConjSpinor 
      | Ga | Gl -> Prop_Feynman
      | Wp | Wm | Z | WHp | WHm | ZH | AH -> Prop_Unitarity
      | Phip | Phim | Phi0 -> Only_Insertion
      | H | Eta | Psi0 | Psi1 | Psip | Psim | Psipp | Psimm -> Prop_Scalar

(* Optionally, ask for the fudge factor treatment for the widths of
   charged particles.  Currently, this only applies to $W^\pm$ and top. *)

    let width f =
      if !use_fudged_width then
        match f with
        | Wp | Wm | U 3 | U (-3) 
        | WHp | WHm | ZH | AH 
        | Uodd _ | Dodd _ | Nodd _ | Lodd _
        | Topp | Toppb -> Fudged
        | _ -> !default_width
      else
        !default_width

    let goldstone = function
      | Wp -> Some (Phip, Coupling.Integer 1)
      | Wm -> Some (Phim, Coupling.Integer 1)
      | Z -> Some (Phi0, Coupling.Integer 1)
      | _ -> None

    let conjugate = function
      | L n -> L (-n) | N n -> N (-n)
      | Lodd n -> L (-n) | Nodd n -> N (-n)
      | U n -> U (-n) | D n -> D (-n)
      | Uodd n -> U (-n) | Dodd n -> D (-n)
      | Topp -> Toppb | Toppb -> Topp
      | Gl -> Gl | Ga -> Ga | Z -> Z 
      | Wp -> Wm | Wm -> Wp | WHm -> WHp
      | WHp -> WHm | ZH -> ZH | AH -> AH
      | Psi0 -> Psi0 | Psi1 -> Psi1 | Psip -> Psim
      | Psim -> Psip | Psipp -> Psimm | Psimm -> Psipp
      | Phip -> Phim | Phim -> Phip | Phi0 -> Phi0
      | H -> H | Eta -> Eta

    let fermion = function
      | L n -> if n > 0 then 1 else -1
      | N n -> if n > 0 then 1 else -1
      | U n -> if n > 0 then 1 else -1
      | D n -> if n > 0 then 1 else -1
      | Lodd n -> if n > 0 then 1 else -1
      | Nodd n -> if n > 0 then 1 else -1
      | Uodd n -> if n > 0 then 1 else -1
      | Dodd n -> if n > 0 then 1 else -1
      | Topp -> 1 | Toppb -> -1
      | Gl | Ga | Z | Wp | Wm | WHp | WHm | AH | ZH -> 0
      | _ -> 0

    module Ch = Charges.QQ
    let ( // ) = Algebra.Small_Rational.make

    let charge = function
       | L n | Lodd n -> if n > 0 then -1//1 else  1//1
       | N n | Nodd n -> 0//1
       | U n | Uodd n -> if n > 0 then  2//3 else -2//3
       | D n | Dodd n -> if n > 0 then -1//3 else  1//3
       | Topp -> 2//3
       | Toppb -> -2//3
       | Gl | Ga | Z | AH | ZH -> 0//1
       | Wp | WHp ->  1//1
       | Wm | WHm -> -1//1
       | H | Phi0 | Eta | Psi1 | Psi0 ->  0//1
       | Phip | Psip ->  1//1
       | Phim | Psim -> -1//1
       | Psipp -> 2//1
       | Psimm -> -2//1

    let lepton = function
       | L n | N n | Lodd n | Nodd n 
          -> if n > 0 then 1//1 else -1//1
       | U _ | D _ | _ -> 0//1

    let baryon = function
       | L _ | N _ -> 0//1
       | U n | D n | Uodd n | Dodd n 
          -> if n > 0 then 1//1 else -1//1
       | Topp -> 1//1
       | Toppb -> -1//1
       | _ -> 0//1

    let charges f = 
      [ charge f; lepton f; baryon f]

    type constant =
      | Unit | Pi | Alpha_QED | Sin2thw
      | Sinthw | Costhw | E | G_weak | Vev | VHeavy
      | Supp | Supp2
      | Sinpsi | Cospsi | Atpsi | Sccs  (* Mixing angles of SU(2) *)
      | Q_lepton | Q_up | Q_down | Q_Z_up | G_CC | G_CCtop
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down | G_NC_heavy
      | G_NC_h_neutrino | G_NC_h_lepton | G_NC_h_up | G_NC_h_down 
      | G_CC_heavy | G_ZHTHT | G_ZTHT | G_AHTHTH | G_AHTHT | G_AHTT
      | G_CC_WH | G_CC_W 
      | Gs | I_Gs | G2
      | I_Q_W | I_G_ZWW | I_G_WWW
      | I_G_AHWW | I_G_ZHWW | I_G_ZWHW | I_G_AHWHWH | I_G_ZHWHWH
      | I_G_AHWHW | I_Q_H  
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | G_WH4 | G_WHWHWW | G_WHWWW | G_WH3W 
      | G_WWAAH | G_WWAZH | G_WWZZH | G_WWZAH | G_WHWHAAH 
      | G_WHWHAZH | G_WHWHZZH | G_WHWHZAH | G_WWZHAH 
      | G_WHWHZHAH | G_WHWZZ | G_WHWAZ | G_WHWAAH | G_WHWZAH
      | G_WHWZHZH | G_WHWZHAH | G_WHWAZH | G_WHWZZH
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ
      | G_PsiWW | G_PsiWHW | G_PsiZZ | G_PsiZHZH 
      | G_PsiZHZ | G_PsiZAH | G_PsiZHAH | G_PsiAHAH
      | G_PsiZW | G_PsiZWH | G_PsiAHW | G_PsiAHWH 
      | G_PsiZHW | G_PsiZHWH
      | G_PsippWW | G_PsippWHW | G_PsippWHWH
      | G_PsiHW | G_PsiHWH | G_Psi0W | G_Psi0WH 
      | G_Psi1W | G_Psi1WH | G_PsiPPW | G_PsiPPWH
      | G_Psi1HAH | G_Psi01AH | G_AHPsip | G_Psi1HZ
      | G_Psi1HZH | G_Psi01Z | G_Psi01ZH | G_ZPsip | G_ZPsipp | G_ZHPsipp
      | G_HHAA | G_HHWHW | G_HHZHZ | G_HHAHZ | G_HHZHAH
      | G_HPsi0WW | G_HPsi0WHW | G_HPsi0ZZ
      | G_HPsi0ZHZH | G_HPsi0ZHZ | G_HPsi0AHAH | G_HPsi0ZAH | G_HPsi0ZHAH 
      | G_HPsipWA | G_HPsipWHA | G_HPsipWZ | G_HPsipWHZ | G_HPsipWAH
      | G_HPsipWHAH | G_HPsipWZH | G_HPsipWHZH | G_HPsippWW | G_HPsippWHWH
      | G_HPsippWHW | G_Psi00ZH | G_Psi00AH | G_Psi00ZHAH
      | G_Psi0pWA | G_Psi0pWHA | G_Psi0pWZ | G_Psi0pWHZ | G_Psi0pWAH
      | G_Psi0pWHAH | G_Psi0pWZH | G_Psi0pWHZH | G_Psi0ppWW | G_Psi0ppWHWH
      | G_Psi0ppWHW | I_G_Psi0pWA | I_G_Psi0pWHA | I_G_Psi0pWZ | I_G_Psi0pWHZ 
      | I_G_Psi0pWAH | I_G_Psi0pWHAH | I_G_Psi0pWZH | I_G_Psi0pWHZH 
      | I_G_Psi0ppWW | I_G_Psi0ppWHWH | I_G_Psi0ppWHW 
      | G_PsippZZ | G_PsippZHZH | G_PsippAZ | G_PsippAAH | G_PsippZAH
      | G_PsippWA | G_PsippWHA | G_PsippWZ | G_PsippWHZ | G_PsippWAH
      | G_PsippWHAH | G_PsippWZH | G_PsippWHZH
      | G_PsiccZZ | G_PsiccAZ | G_PsiccAAH | G_PsiccZZH | G_PsiccAZH
      | G_PsiccZAH
      | G_Htt | G_Hbb | G_Hcc | G_Htautau | G_H3 | G_H4
      | G_Hthth | G_Htht | G_Ethth | G_Etht | G_Ett
      | G_HHtt | G_HHthth | G_HHtht
      | G_Psi0tt | G_Psi0bb | G_Psi0cc | G_Psi0tautau
      | G_Psi1tt | G_Psi1bb | G_Psi1cc | G_Psi1tautau
      | G_Psipq3 | G_Psipq2 | G_Psipl3 | G_Psi0tth | G_Psi1tth
      | G_Psipbth | G_Ebb 
      | G_HGaGa | G_HGaZ | G_EGaGa | G_EGaZ | G_EGlGl
      | G_HWHW | G_HWHWH | G_HAHAH | G_HZHZ | G_HZHAH | G_HAHZ
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let input_parameters =
      []

    let derived_parameters =
      []

    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

    let nc_coupling c t3 q =
      (Real_Array c,
       [Prod [g_over_2_costh; Diff (t3, Prod [Integer 2; q; Atom Sin2thw])];
        Prod [g_over_2_costh; t3]])

    let half = Quot (Integer 1, Integer 2)

    let derived_parameter_arrays =
      [ nc_coupling G_NC_neutrino half (Integer 0);
        nc_coupling G_NC_lepton (Neg half) (Integer (-1));
        nc_coupling G_NC_up half (Quot (Integer 2, Integer 3));
        nc_coupling G_NC_down (Neg half) (Quot (Integer (-1), Integer 3));
        nc_coupling G_NC_h_neutrino half (Integer 0);
        nc_coupling G_NC_h_lepton (Neg half) (Integer (-1));
        nc_coupling G_NC_h_up half (Quot (Integer 2, Integer 3));
        nc_coupling G_NC_h_down (Neg half) (Quot (Integer (-1), Integer 3)) ]

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

    let electromagnetic_currents n =
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);  
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down);
	  ((Lodd (-n), Ga, Lodd n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((Uodd (-n), Ga, Uodd n), FBF (1, Psibar, V, Psi), Q_up);  
          ((Dodd (-n), Ga, Dodd n), FBF (1, Psibar, V, Psi), Q_down) ]

    let color_currents n =
        [ ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs);  
          ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs);
	  ((Uodd (-n), Gl, Uodd n), FBF ((-1), Psibar, V, Psi), Gs);  
          ((Dodd (-n), Gl, Dodd n), FBF ((-1), Psibar, V, Psi), Gs) ]  

    let neutral_currents n =
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

(* The sign of this coupling is just the one of the T3, being -(1/2) for
   leptons and down quarks, and +(1/2) for neutrinos and up quarks. *)

    let neutral_heavy_currents n =
       ([ ((L (-n), ZH, L n), FBF ((-1), Psibar, VL, Psi), G_NC_heavy);
          ((N (-n), ZH, N n), FBF (1, Psibar, VL, Psi), G_NC_heavy);
          ((U (-n), ZH, U n), FBF (1, Psibar, VL, Psi), G_NC_heavy);
          ((D (-n), ZH, D n), FBF ((-1), Psibar, VL, Psi), G_NC_heavy)]
        @
          (if Flags.u1_gauged then
        [ ((L (-n), AH, L n), FBF (1, Psibar, VA, Psi), G_NC_h_lepton);
          ((N (-n), AH, N n), FBF (1, Psibar, VA, Psi), G_NC_h_neutrino);
          ((D (-n), AH, D n), FBF (1, Psibar, VA, Psi), G_NC_h_down)]
          else
            []))

   let heavy_top_currents = 
       ([ ((Toppb, Ga, Topp), FBF (1, Psibar, V, Psi), Q_up);
          ((Toppb, Z, Topp), FBF (1, Psibar, V, Psi), Q_Z_up);
	  ((Toppb, Gl, Topp), FBF (1, Psibar, V, Psi), Gs);
          ((Toppb, Z, U 3), FBF (1, Psibar, VL, Psi), G_ZTHT);
          ((U (-3), Z, Topp), FBF (1, Psibar, VL, Psi), G_ZTHT);
          ((Toppb, ZH, U 3), FBF (1, Psibar, VL, Psi), G_ZHTHT);
          ((U (-3), ZH, Topp), FBF (1, Psibar, VL, Psi), G_ZHTHT);
          ((U (-3), Wp, D 3), FBF (1, Psibar, VL, Psi), G_CCtop);
          ((D (-3), Wm, U 3), FBF (1, Psibar, VL, Psi), G_CCtop);
          ((Toppb, WHp, D 3), FBF (1, Psibar, VL, Psi), G_CC_WH);
          ((D (-3), WHm, Topp), FBF (1, Psibar, VL, Psi), G_CC_WH);
          ((Toppb, Wp, D 3), FBF (1, Psibar, VL, Psi), G_CC_W);
          ((D (-3), Wm, Topp), FBF (1, Psibar, VL, Psi), G_CC_W)] 
          @ 
            (if Flags.u1_gauged then
        [ ((U (-3), AH, U 3), FBF (1, Psibar, VA, Psi), G_AHTT);
          ((Toppb, AH, Topp), FBF (1, Psibar, VA, Psi), G_AHTHTH);
          ((Toppb, AH, U 3), FBF (1, Psibar, VR, Psi), G_AHTHT);
          ((U (-3), AH, Topp), FBF (1, Psibar, VR, Psi), G_AHTHT)]
            else
              []))


(* \begin{equation}
     \mathcal{L}_{\textrm{CC}} =
        - \frac{g}{2\sqrt2} \sum_i \bar\psi_i
               (T^+\fmslash{W}^+ + T^-\fmslash{W}^-)(1-\gamma_5)\psi_i 
   \end{equation} *)

    let charged_currents n =
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC);
          ((L (-n), WHm, N n), FBF (1, Psibar, VL, Psi), G_CC_heavy);
          ((N (-n), WHp, L n), FBF (1, Psibar, VL, Psi), G_CC_heavy);
          ((D (-n), WHm, U n), FBF (1, Psibar, VL, Psi), G_CC_heavy);
          ((U (-n), WHp, D n), FBF (1, Psibar, VL, Psi), G_CC_heavy)]

    let quark_currents n =
       ([ ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC)]
          @
            (if Flags.u1_gauged then
        [ ((U (-n), AH, U n), FBF (1, Psibar, VA, Psi), G_NC_h_up)]
            else
              []))
  

(* We specialize the third generation since there is an additional shift 
   coming from the admixture of the heavy top quark. The universal shift, 
   coming from the mixing in the non-Abelian gauge boson sector is 
   unobservable. (Redefinition of coupling constants by measured ones. *)

    let yukawa =
      [ ((U (-3), H, U 3), FBF (1, Psibar, S, Psi), G_Htt);
        ((D (-3), H, D 3), FBF (1, Psibar, S, Psi), G_Hbb);
        ((U (-2), H, U 2), FBF (1, Psibar, S, Psi), G_Hcc);
        ((L (-3), H, L 3), FBF (1, Psibar, S, Psi), G_Htautau)] 

    let yukawa_add' = 
      [ ((Toppb, H, Topp), FBF (1, Psibar, S, Psi), G_Hthth);
        ((Toppb, H, U 3), FBF (1, Psibar, SLR, Psi), G_Htht);
        ((U (-3), H, Topp), FBF (1, Psibar, SLR, Psi), G_Htht);
        ((U (-3), Psi0, U 3), FBF (1, Psibar, S, Psi), G_Psi0tt);
        ((D (-3), Psi0, D 3), FBF (1, Psibar, S, Psi), G_Psi0bb);
        ((U (-2), Psi0, U 2), FBF (1, Psibar, S, Psi), G_Psi0cc);
        ((L (-3), Psi0, L 3), FBF (1, Psibar, S, Psi), G_Psi0tautau);
        ((U (-3), Psi1, U 3), FBF (1, Psibar, P, Psi), G_Psi1tt);
        ((D (-3), Psi1, D 3), FBF (1, Psibar, P, Psi), G_Psi1bb);
        ((U (-2), Psi1, U 2), FBF (1, Psibar, P, Psi), G_Psi1cc);
        ((L (-3), Psi1, L 3), FBF (1, Psibar, P, Psi), G_Psi1tautau);
        ((U (-3), Psip, D 3), FBF (1, Psibar, SLR, Psi), G_Psipq3);
        ((U (-2), Psip, D 2), FBF (1, Psibar, SLR, Psi), G_Psipq2);
        ((N (-3), Psip, L 3), FBF (1, Psibar, SR, Psi), G_Psipl3);
        ((D (-3), Psim, U 3), FBF (1, Psibar, SLR, Psi), G_Psipq3);
        ((D (-2), Psim, U 2), FBF (1, Psibar, SLR, Psi), G_Psipq2);
        ((L (-3), Psim, N 3), FBF (1, Psibar, SL, Psi), G_Psipl3);
        ((Toppb, Psi0, U 3), FBF (1, Psibar, SL, Psi), G_Psi0tth);
        ((U (-3), Psi0, Topp), FBF (1, Psibar, SR, Psi), G_Psi0tth);
        ((Toppb, Psi1, U 3), FBF (1, Psibar, SL, Psi), G_Psi1tth);
        ((U (-3), Psi1, Topp), FBF (1, Psibar, SR, Psi), G_Psi1tth);
        ((Toppb, Psip, D 3), FBF (1, Psibar, SL, Psi), G_Psipbth);
        ((D (-3), Psim, Topp), FBF (1, Psibar, SR, Psi), G_Psipbth)]
 
    let yukawa_add = 
      if Flags.u1_gauged then
        yukawa_add'
      else
        yukawa_add' @
      [ ((U (-3), Eta, U 3), FBF (1, Psibar, P, Psi), G_Ett);
        ((Toppb, Eta, U 3), FBF (1, Psibar, SLR, Psi), G_Etht);
        ((D (-3), Eta, D 3), FBF (1, Psibar, P, Psi), G_Ebb);
        ((U (-3), Eta, Topp), FBF (1, Psibar, SLR, Psi), G_Etht)]
      
(* \begin{equation}
     \mathcal{L}_{\textrm{TGC}} =
        - e \partial_\mu A_\nu W_+^\mu W_-^\nu + \ldots
        - e \cot\theta_w  \partial_\mu Z_\nu W_+^\mu W_-^\nu + \ldots
   \end{equation} *)

(* Check. *)

    let standard_triple_gauge =
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);           
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs) ]

    let heavy_triple_gauge =
       ([ ((Ga, WHm, WHp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, WHm, WHp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((ZH, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZHWW);
          ((Z, WHm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWHW);
          ((Z, Wm, WHp), Gauge_Gauge_Gauge (-1), I_G_ZWHW);
          ((ZH, WHm, Wp), Gauge_Gauge_Gauge 1, I_G_WWW);
          ((ZH, Wm, WHp), Gauge_Gauge_Gauge (-1), I_G_WWW);
          ((ZH, WHm, WHp), Gauge_Gauge_Gauge (-1), I_G_ZHWHWH)]          
          @
            (if Flags.u1_gauged then
        [ ((AH, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_AHWW);
          ((AH, WHm, Wp), Gauge_Gauge_Gauge 1, I_G_AHWHW);
          ((AH, Wm, WHp), Gauge_Gauge_Gauge (-1), I_G_AHWHW);
          ((AH, WHm, WHp), Gauge_Gauge_Gauge 1, I_G_AHWHWH)]
            else
              []))

    let triple_gauge =
        standard_triple_gauge @ heavy_triple_gauge

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let standard_quartic_gauge =
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          (Gl, Gl, Gl, Gl), gauge4, G2]

    let heavy_quartic_gauge = 
     [ (WHm, Wp, WHm, Wp), gauge4, G_WWWW;
        (Wm, WHp, Wm, WHp), gauge4, G_WWWW;
        (WHm, WHp, WHm, WHp), gauge4, G_WH4;
        (Wm, Wp, WHm, WHp), gauge4, G_WHWHWW;        
        (Wm, Wp, Wm, WHp), gauge4, G_WHWWW;
        (Wm, Wp, WHm, Wp), gauge4, G_WHWWW;
        (WHm, WHp, Wm, WHp), gauge4, G_WH3W;
        (WHm, WHp, WHm, Wp), gauge4, G_WH3W;
        (WHm, Z, WHp, Z), minus_gauge4, G_ZZWW;
        (WHm, Z, WHp, Ga), minus_gauge4, G_AZWW;
        (WHm, Ga, WHp, ZH), minus_gauge4, G_AAWW;
        (WHm, Z, WHp, ZH), minus_gauge4, G_ZZWW;
        (Wm, ZH, Wp, ZH), minus_gauge4, G_WWWW;
        (Wm, Ga, Wp, ZH), minus_gauge4, G_WWAZH;
        (Wm, Z, Wp, ZH), minus_gauge4, G_WWZZH;
        (WHm, Ga, WHp, ZH), minus_gauge4, G_WHWHAZH;
        (WHm, Z, WHp, ZH), minus_gauge4, G_WHWHZZH;
        (WHm, ZH, WHp, ZH), minus_gauge4, G_WH4;
        (WHm, Z, Wp, Z), minus_gauge4, G_WHWZZ;
        (Wm, Z, WHp, Z), minus_gauge4, G_WHWZZ;
        (WHm, Ga, Wp, Z), minus_gauge4, G_WHWAZ;
        (Wm, Ga, WHp, Z), minus_gauge4, G_WHWAZ;
        (WHm, ZH, Wp, ZH), minus_gauge4, G_WHWZHZH;
        (Wm, ZH, WHp, ZH), minus_gauge4, G_WHWZHZH;
        (WHm, Ga, Wp, ZH), minus_gauge4, G_WHWAZH;
        (Wm, Ga, WHp, ZH), minus_gauge4, G_WHWAZH;
        (WHm, Z, Wp, ZH), minus_gauge4, G_WHWZZH;
        (Wm, Z, WHp, ZH), minus_gauge4, G_WHWZZH]
        @ 
          (if Flags.u1_gauged then
      [ (Wm, Ga, Wp, AH), minus_gauge4, G_WWAAH;            
        (Wm, Z, Wp, AH), minus_gauge4, G_WWZAH;
        (WHm, Ga, WHp, AH), minus_gauge4, G_WHWHAAH;
        (WHm, Z, WHp, AH), minus_gauge4, G_WHWHZAH;
        (Wm, ZH, Wp, AH), minus_gauge4, G_WWZHAH;
        (WHm, ZH, WHp, AH), minus_gauge4, G_WHWHZHAH;
        (WHm, Ga, Wp, AH), minus_gauge4, G_WHWAAH;
        (Wm, Ga, WHp, AH), minus_gauge4, G_WHWAAH;        
        (WHm, Z, Wp, AH), minus_gauge4, G_WHWZAH;
        (Wm, Z, WHp, AH), minus_gauge4, G_WHWZAH;
        (WHm, ZH, Wp, AH), minus_gauge4, G_WHWZHAH;
        (Wm, ZH, WHp, AH), minus_gauge4, G_WHWZHAH]
          else
            [])

    let quartic_gauge =
      standard_quartic_gauge @ heavy_quartic_gauge

    let standard_gauge_higgs' =
      [ ((H, Wp, Wm), Scalar_Vector_Vector 1, G_HWW);
        ((H, Z, Z), Scalar_Vector_Vector 1, G_HZZ) ]

    let heavy_gauge_higgs = 
      [ ((H, Wp, WHm), Scalar_Vector_Vector 1, G_HWHW);
        ((H, WHp, Wm), Scalar_Vector_Vector 1, G_HWHW);
        ((H, WHp, WHm), Scalar_Vector_Vector 1, G_HWHWH);
        ((H, ZH, ZH), Scalar_Vector_Vector 1, G_HWHWH);
        ((H, ZH, Z), Scalar_Vector_Vector 1, G_HZHZ);
        ((H, Wp, Wm), Scalar_Vector_Vector 1, G_HZHAH)]
        @ 
          (if Flags.u1_gauged then
       [((H, AH, AH), Scalar_Vector_Vector 1, G_HAHAH);            
        ((H, Z, AH), Scalar_Vector_Vector 1, G_HAHZ)]
          else
            [])

    let triplet_gauge_higgs = 
      [ ((Psi0, Wp, Wm), Scalar_Vector_Vector 1, G_PsiWW);
        ((Psi0, WHp, WHm), Scalar_Vector_Vector (-1), G_PsiWW);
        ((Psi0, WHp, Wm), Scalar_Vector_Vector 1, G_PsiWHW);
        ((Psi0, WHm, Wp), Scalar_Vector_Vector 1, G_PsiWHW);
        ((Psi0, Z, Z), Scalar_Vector_Vector 1, G_PsiZZ);        
        ((Psi0, ZH, ZH), Scalar_Vector_Vector 1, G_PsiZHZH);        
        ((Psi0, ZH, Z), Scalar_Vector_Vector 1, G_PsiZHZ);        
        ((Psim, Wp, Z), Scalar_Vector_Vector 1, G_PsiZW);
        ((Psip, Wm, Z), Scalar_Vector_Vector 1, G_PsiZW);
        ((Psim, WHp, Z), Scalar_Vector_Vector 1, G_PsiZWH);
        ((Psip, WHm, Z), Scalar_Vector_Vector 1, G_PsiZWH);
        ((Psim, Wp, ZH), Scalar_Vector_Vector 1, G_PsiZHW);
        ((Psip, Wm, ZH), Scalar_Vector_Vector 1, G_PsiZHW);
        ((Psim, WHp, ZH), Scalar_Vector_Vector 1, G_PsiZHWH);
        ((Psip, WHm, ZH), Scalar_Vector_Vector 1, G_PsiZHWH);
        ((Psimm, Wp, Wp), Scalar_Vector_Vector 1, G_PsippWW);
        ((Psipp, Wm, Wm), Scalar_Vector_Vector 1, G_PsippWW);
        ((Psimm, WHp, Wp), Scalar_Vector_Vector 1, G_PsippWHW);
        ((Psipp, WHm, Wm), Scalar_Vector_Vector 1, G_PsippWHW);
        ((Psimm, WHp, WHp), Scalar_Vector_Vector 1, G_PsippWHWH);
        ((Psipp, WHm, WHm), Scalar_Vector_Vector 1, G_PsippWHWH)] 
        @
          (if Flags.u1_gauged then
       [((Psi0, AH, Z), Scalar_Vector_Vector 1, G_PsiZAH);        
        ((Psi0, AH, ZH), Scalar_Vector_Vector 1, G_PsiZHAH);        
        ((Psi0, AH, AH), Scalar_Vector_Vector 1, G_PsiAHAH);
        ((Psim, Wp, AH), Scalar_Vector_Vector 1, G_PsiAHW);
        ((Psip, Wm, AH), Scalar_Vector_Vector 1, G_PsiAHW);
        ((Psim, WHp, AH), Scalar_Vector_Vector 1, G_PsiAHWH);
        ((Psip, WHm, AH), Scalar_Vector_Vector 1, G_PsiAHWH)]            
          else
            [])

    let triplet_gauge2_higgs = 
      [ ((Wp, H, Psim), Vector_Scalar_Scalar 1, G_PsiHW);
        ((Wm, H, Psip), Vector_Scalar_Scalar 1, G_PsiHW);
        ((WHp, H, Psim), Vector_Scalar_Scalar 1, G_PsiHWH);
        ((WHm, H, Psip), Vector_Scalar_Scalar 1, G_PsiHWH);
        ((Wp, Psi0, Psim), Vector_Scalar_Scalar 1, G_Psi0W);
        ((Wm, Psi0, Psip), Vector_Scalar_Scalar 1, G_Psi0W);
        ((WHp, Psi0, Psim), Vector_Scalar_Scalar 1, G_Psi0WH);
        ((WHm, Psi0, Psip), Vector_Scalar_Scalar 1, G_Psi0WH);
        ((Wp, Psi1, Psim), Vector_Scalar_Scalar 1, G_Psi1W);
        ((Wm, Psi1, Psip), Vector_Scalar_Scalar (-1), G_Psi1W);
        ((WHp, Psi1, Psim), Vector_Scalar_Scalar 1, G_Psi1WH);
        ((WHm, Psi1, Psip), Vector_Scalar_Scalar (-1), G_Psi1WH);
        ((Wp, Psip, Psimm), Vector_Scalar_Scalar 1, G_PsiPPW);
        ((Wm, Psim, Psipp), Vector_Scalar_Scalar 1, G_PsiPPW);
        ((WHp, Psip, Psimm), Vector_Scalar_Scalar 1, G_PsiPPWH);
        ((WHm, Psim, Psipp), Vector_Scalar_Scalar 1, G_PsiPPWH);
        ((Ga, Psip, Psim), Vector_Scalar_Scalar 1, Q_lepton);
        ((Ga, Psipp, Psimm), Vector_Scalar_Scalar 2, Q_lepton);
        ((Z, H, Psi1), Vector_Scalar_Scalar 1, G_Psi1HZ);
        ((ZH, H, Psi1), Vector_Scalar_Scalar 1, G_Psi1HZH);
        ((Z, Psi0, Psi1), Vector_Scalar_Scalar 1, G_Psi01Z);
        ((ZH, Psi0, Psi1), Vector_Scalar_Scalar 1, G_Psi01ZH);
        ((Z, Psip, Psim), Vector_Scalar_Scalar 1, G_ZPsip);
        ((Z, Psipp, Psimm), Vector_Scalar_Scalar 2, G_ZPsipp);
        ((ZH, Psipp, Psimm), Vector_Scalar_Scalar 2, G_ZHPsipp)]        
        @ 
          (if Flags.u1_gauged then
       [((AH, H, Psi1), Vector_Scalar_Scalar 1, G_Psi1HAH);
        ((AH, Psi0, Psi1), Vector_Scalar_Scalar 1, G_Psi01AH); 
        ((AH, Psip, Psim), Vector_Scalar_Scalar 1, G_AHPsip);
        ((AH, Psipp, Psimm), Vector_Scalar_Scalar 2, G_AHPsip)]
          else [])
       
    let standard_gauge_higgs = 
      standard_gauge_higgs' @ heavy_gauge_higgs @ triplet_gauge_higgs @
      triplet_gauge2_higgs

    let standard_gauge_higgs4 =
      [ (H, H, Wp, Wm), Scalar2_Vector2 1, G_HHWW;
        (H, H, Z, Z), Scalar2_Vector2 1, G_HHZZ ]

    let littlest_gauge_higgs4 = 
      [ (H, H, WHp, WHm), Scalar2_Vector2 (-1), G_HHWW;
        (H, H, ZH, ZH), Scalar2_Vector2 (-1), G_HHWW;
        (H, H, Wp, WHm), Scalar2_Vector2 1, G_HHWHW;
        (H, H, WHp, Wm), Scalar2_Vector2 1, G_HHWHW;
        (H, H, ZH, Z), Scalar2_Vector2 (-1), G_HHZHZ;
        (H, Psi0, Wp, Wm), Scalar2_Vector2 1, G_HPsi0WW;
        (H, Psi0, WHp, WHm), Scalar2_Vector2 (-1), G_HPsi0WW;
        (H, Psi0, WHp, Wm), Scalar2_Vector2 1, G_HPsi0WHW;
        (H, Psi0, Wp, WHm), Scalar2_Vector2 1, G_HPsi0WHW;
        (H, Psi0, Z, Z), Scalar2_Vector2 1, G_HPsi0ZZ;
        (H, Psi0, ZH, ZH), Scalar2_Vector2 1, G_HPsi0ZHZH;
        (H, Psi0, ZH, Z), Scalar2_Vector2 1, G_HPsi0ZHZ;
        (H, Psim, Wp, Ga), Scalar2_Vector2 1, G_HPsipWA;
        (H, Psip, Wm, Ga), Scalar2_Vector2 1, G_HPsipWA;
        (H, Psim, WHp, Ga), Scalar2_Vector2 1, G_HPsipWHA;
        (H, Psip, WHm, Ga), Scalar2_Vector2 1, G_HPsipWHA;
        (H, Psim, Wp, Z), Scalar2_Vector2 1, G_HPsipWZ;
        (H, Psip, Wm, Z), Scalar2_Vector2 1, G_HPsipWZ;
        (H, Psim, WHp, Z), Scalar2_Vector2 1, G_HPsipWHZ;
        (H, Psip, WHm, Z), Scalar2_Vector2 1, G_HPsipWHZ;
        (H, Psim, Wp, ZH), Scalar2_Vector2 1, G_HPsipWZH;
        (H, Psip, Wm, ZH), Scalar2_Vector2 1, G_HPsipWZH;
        (H, Psim, WHp, ZH), Scalar2_Vector2 1, G_HPsipWHZH;
        (H, Psip, WHm, ZH), Scalar2_Vector2 1, G_HPsipWHZH;
        (H, Psimm, Wp, Wp), Scalar2_Vector2 1, G_HPsippWW;
        (H, Psipp, Wm, Wm), Scalar2_Vector2 1, G_HPsippWW;
        (H, Psimm, WHp, WHp), Scalar2_Vector2 1, G_HPsippWHWH;
        (H, Psipp, WHm, WHm), Scalar2_Vector2 1, G_HPsippWHWH;
        (H, Psimm, WHp, Wp), Scalar2_Vector2 1, G_HPsippWHW;
        (H, Psipp, WHm, Wm), Scalar2_Vector2 1, G_HPsippWHW;
        (Psi0, Psi0, Wp, Wm), Scalar2_Vector2 2, G_HHWW;
        (Psi0, Psi0, WHp, WHm), Scalar2_Vector2 (-2), G_HHWW;
        (Psi0, Psi0, Z, Z), Scalar2_Vector2 4, G_HHZZ;
        (Psi0, Psi0, ZH, ZH), Scalar2_Vector2 1, G_Psi00ZH;
        (Psi0, Psi0, WHp, Wm), Scalar2_Vector2 2, G_HHWHW;
        (Psi0, Psi0, Wp, WHm), Scalar2_Vector2 2, G_HHWHW;        
        (Psi0, Psi0, Z, ZH), Scalar2_Vector2 4, G_HHZHZ;
        (Psi0, Psim, Wp, Ga), Scalar2_Vector2 1, G_Psi0pWA;
        (Psi0, Psip, Wm, Ga), Scalar2_Vector2 1, G_Psi0pWA;
        (Psi0, Psim, WHp, Ga), Scalar2_Vector2 1, G_Psi0pWHA;
        (Psi0, Psip, WHm, Ga), Scalar2_Vector2 1, G_Psi0pWHA;
        (Psi0, Psim, Wp, Z), Scalar2_Vector2 1, G_Psi0pWZ;
        (Psi0, Psip, Wm, Z), Scalar2_Vector2 1, G_Psi0pWZ;
        (Psi0, Psim, WHp, Z), Scalar2_Vector2 1, G_Psi0pWHZ;
        (Psi0, Psip, WHm, Z), Scalar2_Vector2 1, G_Psi0pWHZ;
        (Psi0, Psim, Wp, ZH), Scalar2_Vector2 1, G_Psi0pWZH;
        (Psi0, Psip, Wm, ZH), Scalar2_Vector2 1, G_Psi0pWZH;
        (Psi0, Psim, WHp, ZH), Scalar2_Vector2 1, G_Psi0pWHZH;
        (Psi0, Psip, WHm, ZH), Scalar2_Vector2 1, G_Psi0pWHZH;
        (Psi0, Psimm, Wp, Wp), Scalar2_Vector2 1, G_Psi0ppWW;
        (Psi0, Psipp, Wm, Wm), Scalar2_Vector2 1, G_Psi0ppWW;
        (Psi0, Psimm, WHp, WHp), Scalar2_Vector2 1, G_Psi0ppWHWH;
        (Psi0, Psipp, WHm, WHm), Scalar2_Vector2 1, G_Psi0ppWHWH;
        (Psi0, Psimm, WHp, Wp), Scalar2_Vector2 1, G_Psi0ppWHW;
        (Psi0, Psipp, WHm, Wm), Scalar2_Vector2 1, G_Psi0ppWHW;
        (Psi1, Psi1, Wp, Wm), Scalar2_Vector2 2, G_HHWW;
        (Psi1, Psi1, WHp, WHm), Scalar2_Vector2 (-2), G_HHWW;
        (Psi1, Psi1, Z, Z), Scalar2_Vector2 4, G_HHZZ;
        (Psi1, Psi1, ZH, ZH), Scalar2_Vector2 1, G_Psi00ZH;
        (Psi1, Psi1, WHp, Wm), Scalar2_Vector2 2, G_HHWHW;
        (Psi1, Psi1, Wp, WHm), Scalar2_Vector2 2, G_HHWHW;       
        (Psi1, Psi1, Z, ZH), Scalar2_Vector2 4, G_HHZHZ;
        (Psi1, Psim, Wp, Ga), Scalar2_Vector2 1, I_G_Psi0pWA;
        (Psi1, Psip, Wm, Ga), Scalar2_Vector2 (-1), I_G_Psi0pWA;
        (Psi1, Psim, WHp, Ga), Scalar2_Vector2 1, I_G_Psi0pWHA;
        (Psi1, Psip, WHm, Ga), Scalar2_Vector2 (-1), I_G_Psi0pWHA;
        (Psi1, Psim, Wp, Z), Scalar2_Vector2 1, I_G_Psi0pWZ;
        (Psi1, Psip, Wm, Z), Scalar2_Vector2 (-1), I_G_Psi0pWZ;
        (Psi1, Psim, WHp, Z), Scalar2_Vector2 1, I_G_Psi0pWHZ;
        (Psi1, Psip, WHm, Z), Scalar2_Vector2 (-1), I_G_Psi0pWHZ;
        (Psi1, Psim, Wp, ZH), Scalar2_Vector2 1, I_G_Psi0pWZH;
        (Psi1, Psip, Wm, ZH), Scalar2_Vector2 (-1), I_G_Psi0pWZH;
        (Psi1, Psim, WHp, ZH), Scalar2_Vector2 1, I_G_Psi0pWHZH;
        (Psi1, Psip, WHm, ZH), Scalar2_Vector2 (-1), I_G_Psi0pWHZH;
        (Psi1, Psimm, Wp, Wp), Scalar2_Vector2 1, I_G_Psi0ppWW;
        (Psi1, Psipp, Wm, Wm), Scalar2_Vector2 (-1), I_G_Psi0ppWW;
        (Psi1, Psimm, WHp, WHp), Scalar2_Vector2 1, I_G_Psi0ppWHWH;
        (Psi1, Psipp, WHm, WHm), Scalar2_Vector2 (-1), I_G_Psi0ppWHWH;
        (Psi1, Psimm, WHp, Wp), Scalar2_Vector2 1, I_G_Psi0ppWHW;
        (Psi1, Psipp, WHm, Wm), Scalar2_Vector2 (-1), I_G_Psi0ppWHW;
        (Psip, Psim, Wp, Wm), Scalar2_Vector2 4, G_HHWW;
        (Psip, Psim, WHp, WHm), Scalar2_Vector2 1, G_Psi00ZH;
        (Psip, Psim, WHp, Wm), Scalar2_Vector2 4, G_HHWHW;
        (Psip, Psim, Wp, WHm), Scalar2_Vector2 4, G_HHWHW;        
        (Psip, Psim, Z, Z), Scalar2_Vector2 1, G_PsippZZ;
        (Psip, Psim, Ga, Ga), Scalar2_Vector2 2, G_AAWW;
        (Psip, Psim, ZH, ZH), Scalar2_Vector2 1, G_PsippZHZH;
        (Psip, Psim, Ga, Z), Scalar2_Vector2 4, G_PsippAZ;
        (Psip, Psimm, Wp, Ga), Scalar2_Vector2 1, G_PsippWA;
        (Psim, Psipp, Wm, Ga), Scalar2_Vector2 1, G_PsippWA;
        (Psip, Psimm, WHp, Ga), Scalar2_Vector2 1, G_PsippWHA;
        (Psim, Psipp, WHm, Ga), Scalar2_Vector2 1, G_PsippWHA;
        (Psip, Psimm, Wp, Z), Scalar2_Vector2 1, G_PsippWZ;
        (Psim, Psipp, Wm, Z), Scalar2_Vector2 1, G_PsippWZ;
        (Psip, Psimm, WHp, Z), Scalar2_Vector2 1, G_PsippWHZ;
        (Psim, Psipp, WHm, Z), Scalar2_Vector2 1, G_PsippWHZ;
        (Psip, Psimm, Wp, ZH), Scalar2_Vector2 1, G_PsippWZH;
        (Psim, Psipp, Wm, ZH), Scalar2_Vector2 1, G_PsippWZH;
        (Psip, Psimm, WHp, ZH), Scalar2_Vector2 1, G_PsippWHZH;
        (Psim, Psipp, WHm, ZH), Scalar2_Vector2 1, G_PsippWHZH;
        (Psipp, Psimm, Wp, Wm), Scalar2_Vector2 2, G_HHWW;
        (Psipp, Psimm, WHp, WHm), Scalar2_Vector2 (-2), G_HHWW;
        (Psipp, Psimm, WHp, Wm), Scalar2_Vector2 2, G_HHWHW;
        (Psipp, Psimm, Wp, WHm), Scalar2_Vector2 2, G_HHWHW;        
        (Psipp, Psimm, Z, Z), Scalar2_Vector2 1, G_PsiccZZ;
        (Psipp, Psimm, Ga, Ga), Scalar2_Vector2 8, G_AAWW;
        (Psipp, Psimm, ZH, ZH), Scalar2_Vector2 1, G_Psi00ZH;
        (Psipp, Psimm, Ga, Z), Scalar2_Vector2 1, G_PsiccAZ;
        (Psipp, Psimm, Z, ZH), Scalar2_Vector2 4, G_PsiccZZH;
        (Psipp, Psimm, Ga, ZH), Scalar2_Vector2 4, G_PsiccAZH]
        @
          (if Flags.u1_gauged then
       [(H, H, AH, AH), Scalar2_Vector2 1, G_HHAA;            
        (H, H, AH, Z), Scalar2_Vector2 (-1), G_HHAHZ;
        (H, H, ZH, AH), Scalar2_Vector2 (-1), G_HHZHAH;
        (H, Psi0, AH, AH), Scalar2_Vector2 1, G_HPsi0AHAH;
        (H, Psi0, Z, AH), Scalar2_Vector2 1, G_HPsi0ZAH;
        (H, Psi0, ZH, AH), Scalar2_Vector2 1, G_HPsi0ZHAH;
        (H, Psim, Wp, AH), Scalar2_Vector2 1, G_HPsipWAH;
        (H, Psip, Wm, AH), Scalar2_Vector2 1, G_HPsipWAH;
        (H, Psim, WHp, AH), Scalar2_Vector2 1, G_HPsipWHAH;
        (H, Psip, WHm, AH), Scalar2_Vector2 1, G_HPsipWHAH;
        (Psi0, Psi0, AH, AH), Scalar2_Vector2 1, G_Psi00AH;
        (Psi0, Psi0, Z, AH), Scalar2_Vector2 4, G_HHAHZ;
        (Psi0, Psi0, AH, ZH), Scalar2_Vector2 1, G_Psi00ZHAH;
        (Psi0, Psim, Wp, AH), Scalar2_Vector2 1, G_Psi0pWAH;
        (Psi0, Psip, Wm, AH), Scalar2_Vector2 1, G_Psi0pWAH;
        (Psi0, Psim, WHp, AH), Scalar2_Vector2 1, G_Psi0pWHAH;
        (Psi0, Psip, WHm, AH), Scalar2_Vector2 1, G_Psi0pWHAH;
        (Psi1, Psi1, AH, AH), Scalar2_Vector2 1, G_Psi00AH;
        (Psi1, Psi1, Z, AH), Scalar2_Vector2 4, G_HHAHZ;
        (Psi1, Psi1, AH, ZH), Scalar2_Vector2 1, G_Psi00ZHAH;
        (Psi1, Psim, Wp, AH), Scalar2_Vector2 1, I_G_Psi0pWAH;
        (Psi1, Psip, Wm, AH), Scalar2_Vector2 (-1), I_G_Psi0pWAH;
        (Psi1, Psim, WHp, AH), Scalar2_Vector2 1, I_G_Psi0pWHAH;
        (Psi1, Psip, WHm, AH), Scalar2_Vector2 (-1), I_G_Psi0pWHAH;
        (Psip, Psim, AH, AH), Scalar2_Vector2 1, G_Psi00AH;
        (Psip, Psim, Ga, AH), Scalar2_Vector2 4, G_PsippAAH;
        (Psip, Psim, Z, AH), Scalar2_Vector2 4, G_PsippZAH;
        (Psip, Psimm, Wp, AH), Scalar2_Vector2 1, G_PsippWAH;
        (Psim, Psipp, Wm, AH), Scalar2_Vector2 1, G_PsippWAH;
        (Psip, Psimm, WHp, AH), Scalar2_Vector2 1, G_PsippWHAH;
        (Psim, Psipp, WHm, AH), Scalar2_Vector2 1, G_PsippWHAH;
        (Psipp, Psimm, AH, AH), Scalar2_Vector2 1, G_Psi00AH;
        (Psipp, Psimm, AH, ZH), Scalar2_Vector2 (-1), G_Psi00ZHAH;
        (Psipp, Psimm, Ga, AH), Scalar2_Vector2 4, G_PsiccAAH;
        (Psipp, Psimm, Z, AH), Scalar2_Vector2 4, G_PsiccZAH]
          else [])

    let standard_higgs =
      [ (H, H, H), Scalar_Scalar_Scalar 1, G_H3 ]
        
   let anomaly_higgs = 
      [ (Eta, Gl, Gl), Dim5_Scalar_Gauge2_Skew 1, G_EGlGl;
        (Eta, Ga, Ga), Dim5_Scalar_Gauge2_Skew 1, G_EGaGa; 
        (Eta, Ga, Z), Dim5_Scalar_Gauge2_Skew 1, G_EGaZ] 
(*    @ [ (H, Ga, Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
        (H, Ga, Z), Dim5_Scalar_Gauge2 1, G_HGaZ ]           *)

    let standard_higgs4 =
      [ (H, H, H, H), Scalar4 1, G_H4 ]

    let gauge_higgs =
        standard_gauge_higgs

    let gauge_higgs4 =
        standard_gauge_higgs4

    let higgs =
        standard_higgs

    let higgs4 =
        standard_higgs4

    let top_quartic = 
      [ ((U (-3), H, H, U 3), GBBG (1, Psibar, S2, Psi), G_HHtt);
   ((Toppb, H, H, Topp), GBBG (1, Psibar, S2, Psi), G_HHthth);
   ((U (-3), H, H, Topp), GBBG (1, Psibar, S2LR, Psi), G_HHtht);
   ((Toppb, H, H, U 3), GBBG (1, Psibar, S2LR, Psi), G_HHtht)]

    let goldstone_vertices =
      [ ((Phi0, Wm, Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((Phip, Ga, Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((Phip, Z, Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((Phim, Wp, Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((Phim, Wp, Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @ 
       ThoList.flatmap neutral_currents [1;2;3] @
       ThoList.flatmap color_currents [1;2;3] @
       ThoList.flatmap neutral_heavy_currents [1;2;3] @       
       ThoList.flatmap charged_currents [1;2;3] @
       ThoList.flatmap quark_currents [1;2] @       
       heavy_top_currents @ 
       (if Flags.u1_gauged then []
           else anomaly_higgs) @
       yukawa @ yukawa_add @ triple_gauge @ 
       gauge_higgs @ higgs @ goldstone_vertices)

    let vertices4 =
      quartic_gauge @ gauge_higgs4 @ higgs4 @ top_quartic

    let vertices () = (vertices3, vertices4, [])

(* For efficiency, make sure that [F.of_vertices vertices] is
   evaluated only once. *)

    let table = F.of_vertices (vertices ())
    let fuse2 = F.fuse2 table
    let fuse3 = F.fuse3 table
    let fuse = F.fuse table
    let max_degree () = 4

    let flavor_of_string = function
      | "e-" -> L 1 | "e+" -> L (-1)
      | "mu-" -> L 2 | "mu+" -> L (-2)
      | "tau-" -> L 3 | "tau+" -> L (-3)
      | "nue" -> N 1 | "nuebar" -> N (-1)
      | "numu" -> N 2 | "numubar" -> N (-2)
      | "nutau" -> N 3 | "nutaubar" -> N (-3)
      | "u" -> U 1 | "ubar" -> U (-1)
      | "c" -> U 2 | "cbar" -> U (-2)
      | "t" -> U 3 | "tbar" -> U (-3)
      | "d" -> D 1 | "dbar" -> D (-1)
      | "s" -> D 2 | "sbar" -> D (-2)
      | "b" -> D 3 | "bbar" -> D (-3)
      | "tp" -> Topp  | "tpbar" -> Toppb
      | "g" -> Gl
      | "A" -> Ga | "Z" | "Z0" -> Z
      | "AH" | "AH0" | "Ah" | "Ah0" -> AH
      | "ZH" | "ZH0" | "Zh" | "Zh0" -> ZH
      | "W+" -> Wp | "W-" -> Wm
      | "WH+" -> WHp | "WH-" -> WHm
      | "H" | "h" -> H | "eta" | "Eta" -> Eta
      | "Psi" | "Psi0" | "psi" | "psi0" -> Psi0
      | "Psi1" | "psi1" -> Psi1
      | "Psi+" | "psi+" | "Psip" | "psip" -> Psip
      | "Psi-" | "psi-" | "Psim" | "psim" -> Psim
      | "Psi++" | "psi++" | "Psipp" | "psipp" -> Psipp
      | "Psi--" | "psi--" | "Psimm" | "psimm" -> Psimm
      | _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_of_string" 

    let flavor_to_string = function
      | L 1 -> "e-" | L (-1) -> "e+"
      | L 2 -> "mu-" | L (-2) -> "mu+"
      | L 3 -> "tau-" | L (-3) -> "tau+"
      | L _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_string" 
      | N 1 -> "nue" | N (-1) -> "nuebar"
      | N 2 -> "numu" | N (-2) -> "numubar"
      | N 3 -> "nutau" | N (-3) -> "nutaubar"
      | N _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_string" 
      | Lodd 1 -> "l1odd-" | Lodd (-1) -> "l1odd+"
      | Lodd 2 -> "l2odd-" | Lodd (-2) -> "l2odd+"
      | Lodd 3 -> "l3odd-" | Lodd (-3) -> "l3odd+"
      | Lodd _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_string" 
      | Nodd 1 -> "n1odd" | Nodd (-1) -> "n1oddbar"
      | Nodd 2 -> "n2odd" | Nodd (-2) -> "n2oddbar"
      | Nodd 3 -> "n3odd" | Nodd (-3) -> "n3oddbar"
      | Nodd _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_string" 
      | U 1 -> "u" | U (-1) -> "ubar"
      | U 2 -> "c" | U (-2) -> "cbar"
      | U 3 -> "t" | U (-3) -> "tbar"
      | U _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_string" 
      | D 1 -> "d" | D (-1) -> "dbar"
      | D 2 -> "s" | D (-2) -> "sbar"
      | D 3 -> "b" | D (-3) -> "bbar"
      | D _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_string" 
      | Uodd 1 -> "uodd" | Uodd (-1) -> "uoddbar"
      | Uodd 2 -> "codd" | Uodd (-2) -> "coddbar"
      | Uodd 3 -> "t1odd" | Uodd (-3) -> "t1oddbar"
      | Uodd 4 -> "t2odd" | Uodd (-4) -> "t2oddbar"
      | Uodd _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_string" 
      | Dodd 1 -> "dodd" | Dodd (-1) -> "doddbar"
      | Dodd 2 -> "sodd" | Dodd (-2) -> "soddbar"
      | Dodd 3 -> "bodd" | Dodd (-3) -> "boddbar"
      | Dodd _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_string" 
      | Topp -> "tp" | Toppb -> "tpbar"
      | Gl -> "g"
      | Ga -> "A" | Z -> "Z"
      | Wp -> "W+" | Wm -> "W-"
      | ZH -> "ZH" | AH -> "AH" | WHp -> "WHp" | WHm -> "WHm"
      | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
      | H -> "H" | Eta -> "Eta"
      | Psi0 -> "Psi0" | Psi1 -> "Psi1" | Psip -> "Psi+" 
      | Psim -> "Psi-" | Psipp -> "Psi++" | Psimm -> "Psi--"

    let flavor_to_TeX = function
      | L 1 -> "e^-" | L (-1) -> "e^+"
      | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
      | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
      | L _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_TeX" 
      | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
      | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
      | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
      | N _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_TeX" 
      | Lodd 1 -> "L_1^-" | Lodd (-1) -> "L_1^+"
      | Lodd 2 -> "L_2^-" | Lodd (-2) -> "L_2^+"
      | Lodd 3 -> "L_3^-" | Lodd (-3) -> "L_3^+"
      | Lodd _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_TeX" 
      | Nodd 1 -> "N_1" | Nodd (-1) -> "\\bar{N}_1"
      | Nodd 2 -> "N_2" | Nodd (-2) -> "\\bar{N}_2"
      | Nodd 3 -> "N_3" | Nodd (-3) -> "\\bar{N}_3"
      | Nodd _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_TeX" 
      | U 1 -> "u" | U (-1) -> "\\bar{u}"
      | U 2 -> "c" | U (-2) -> "\\bar{c}"
      | U 3 -> "t" | U (-3) -> "\\bar{t}"
      | U _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_TeX" 
      | D 1 -> "d" | D (-1) -> "\\bar{d}"
      | D 2 -> "s" | D (-2) -> "\\bar{s}"
      | D 3 -> "b" | D (-3) -> "\\bar{b}"
      | D _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_TeX" 
      | Uodd 1 -> "U" | Uodd (-1) -> "\\bar{U}"
      | Uodd 2 -> "C" | Uodd (-2) -> "\\bar{C}"
      | Uodd 3 -> "T_1" | Uodd (-3) -> "\\bar{T}_1"
      | Uodd 4 -> "T_2" | Uodd (-4) -> "\\bar{T}_2"
      | Uodd _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_TeX" 
      | Dodd 1 -> "D" | Dodd (-1) -> "\\bar{D}"
      | Dodd 2 -> "S" | Dodd (-2) -> "\\bar{S}"
      | Dodd 3 -> "B" | Dodd (-3) -> "\\bar{B}"
      | Dodd _ -> invalid_arg "Modellib_BSM.Littlest_Tpar.flavor_to_TeX" 
      | Topp -> "T^\\prime" | Toppb -> "\\bar{T}^\\prime"
      | Gl -> "g"
      | Ga -> "\\gamma" | Z -> "Z"
      | Wp -> "W^+" | Wm -> "W^-"
      | ZH -> "Z_H" | AH -> "\\gamma_H" | WHp -> "W_H^+" | WHm -> "W_H^-"
      | Phip -> "\\Phi^+" | Phim -> "\\Phi^-" | Phi0 -> "\\Phi^0" 
      | H -> "H" | Eta -> "\\eta"
      | Psi0 -> "\\Psi_S" | Psi1 -> "\\Psi_P" | Psip -> "\\Psi^+" 
      | Psim -> "\\Psi^-" | Psipp -> "\\Psi^{++}" | Psimm -> "\\Psi^{--}"

    let flavor_symbol = function
      | L n when n > 0 -> "l" ^ string_of_int n
      | L n -> "l" ^ string_of_int (abs n) ^ "b"
      | Lodd n when n > 0 -> "lodd" ^ string_of_int n
      | Lodd n -> "lodd" ^ string_of_int (abs n) ^ "b"
      | N n when n > 0 -> "n" ^ string_of_int n
      | N n -> "n" ^ string_of_int (abs n) ^ "b"
      | Nodd n when n > 0 -> "nodd" ^ string_of_int n
      | Nodd n -> "nodd" ^ string_of_int (abs n) ^ "b"
      | U n when n > 0 -> "u" ^ string_of_int n
      | U n -> "u" ^ string_of_int (abs n) ^ "b"
      | D n when n > 0 ->  "d" ^ string_of_int n
      | D n -> "d" ^ string_of_int (abs n) ^ "b"
      | Uodd n when n > 0 -> "uodd" ^ string_of_int n
      | Uodd n -> "uodd" ^ string_of_int (abs n) ^ "b"
      | Dodd n when n > 0 ->  "dodd" ^ string_of_int n
      | Dodd n -> "dodd" ^ string_of_int (abs n) ^ "b"
      | Topp -> "tp" | Toppb -> "tpb" 
      | Gl -> "gl"
      | Ga -> "a" | Z -> "z"
      | Wp -> "wp" | Wm -> "wm"
      | ZH -> "zh" | AH -> "ah" | WHp -> "whp" | WHm -> "whm"
      | Phip -> "pp" | Phim -> "pm" | Phi0 -> "p0" 
      | H -> "h" | Eta -> "eta"
      | Psi0 -> "psi0" | Psi1 -> "psi1" | Psip -> "psip" 
      | Psim -> "psim" | Psipp -> "psipp" | Psimm -> "psimm"

(* There are PDG numbers for Z', Z'', W', 32-34, respectively.
   We just introduce a number 38 for Y0 as a Z'''.
   As well, there is the number 8 for a t'. But we cheat a little bit and 
   take the number 35 which is reserved for a heavy scalar Higgs for the 
   Eta scalar.
   For the heavy Higgs states we take 35 and 36 for the neutral ones, 37 for 
   the charged and 38 for the doubly-charged. 
   The pseudoscalar gets the 39.
   For the odd fermions we add 40 to the values for the SM particles.
*)  

    let pdg = function
      | L n when n > 0 -> 9 + 2*n
      | L n -> - 9 + 2*n
      | N n when n > 0 -> 10 + 2*n
      | N n -> - 10 + 2*n
      | U n when n > 0 -> 2*n
      | U n -> 2*n
      | D n when n > 0 -> - 1 + 2*n
      | D n -> 1 + 2*n
      | Lodd n when n > 0 -> 49 + 2*n
      | Lodd n -> - 49 + 2*n
      | Nodd n when n > 0 -> 50 + 2*n
      | Nodd n -> - 50 + 2*n
      | Uodd n when n > 0 -> 40 + 2*n
      | Uodd n -> -40 + 2*n
      | Dodd n when n > 0 -> 39 + 2*n
      | Dodd n -> -39 + 2*n
      | Topp -> 8 | Toppb -> (-8)
      | Gl -> 21
      | Ga -> 22 | Z -> 23
      | Wp -> 24 | Wm -> (-24)
      | AH -> 32 | ZH -> 33 | WHp -> 34 | WHm -> (-34) 
      | Phip | Phim -> 27 | Phi0 -> 26
      | Psi0 -> 35 | Psi1 -> 36 | Psip -> 37 | Psim -> (-37)
      | Psipp -> 38 | Psimm -> (-38)
      | H -> 25 | Eta -> 39

    let mass_symbol f = 
      "mass(" ^ string_of_int (abs (pdg f)) ^ ")"

    let width_symbol f =
      "width(" ^ string_of_int (abs (pdg f)) ^ ")"

    let constant_symbol = function
      | Unit -> "unit" | Pi -> "PI" | VHeavy -> "vheavy"
      | Alpha_QED -> "alpha" | E -> "e" | G_weak -> "g" | Vev -> "vev"
      | Sin2thw -> "sin2thw" | Sinthw -> "sinthw" | Costhw -> "costhw"
      | Sinpsi -> "sinpsi" | Cospsi -> "cospsi" 
      | Atpsi -> "atpsi" | Sccs -> "sccs"
      | Supp -> "vF" | Supp2 -> "v2F2" 
      | Q_lepton -> "qlep" | Q_up -> "qup" | Q_down -> "qdwn"
      | Q_Z_up -> "qzup" 
      | G_ZHTHT -> "gzhtht" | G_ZTHT -> "gztht"
      | G_AHTHTH -> "gahthth" | G_AHTHT -> "gahtht" | G_AHTT -> "gahtt"
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_CC -> "gcc" | G_CCtop -> "gcctop" | G_CC_heavy -> "gcch" 
      | G_CC_WH -> "gccwh" | G_CC_W -> "gccw" 
      | G_NC_h_lepton -> "gnchlep" | G_NC_h_neutrino -> "gnchneu"
      | G_NC_h_up -> "gnchup" | G_NC_h_down -> "gnchdwn"   
      | G_NC_heavy -> "gnch"  
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" | I_G_WWW -> "igwww"
      | I_G_AHWW -> "igahww" | I_G_ZHWW -> "igzhww" | I_G_ZWHW -> "igzwhw"
      | I_G_AHWHWH -> "igahwhwh" | I_G_ZHWHWH -> "igzhwhwh"
      | I_G_AHWHW -> "igahwhw"
      | I_Q_H -> "iqh" 
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_WH4 -> "gwh4" | G_WHWHWW -> "gwhwhww" | G_WHWWW -> "gwhwww" 
      | G_WH3W -> "gwh3w"
      | G_WWAAH -> "gwwaah" | G_WWAZH -> "gwwazh" | G_WWZZH -> "gwwzzh" 
      | G_WWZAH -> "gwwzah" | G_WHWHAAH -> "gwhwhaah"  
      | G_WHWHAZH -> "gwhwhazh" | G_WHWHZZH -> "gwhwhzzh" 
      | G_WHWHZAH -> "gwhwhzah" 
      | G_WWZHAH -> "gwwzhah" | G_WHWHZHAH -> "gwhwhzhah"
      | G_WHWZZ -> "gwhwzz" | G_WHWAZ -> "gwhwaz" 
      | G_WHWAAH -> "gwhwaah" | G_WHWZAH -> "gwhwzah"
      | G_WHWZHZH -> "gwhwzhzh" | G_WHWZHAH -> "gwhwzhah"
      | G_WHWAZH -> "gwhwazh" | G_WHWZZH -> "gwhwzzh"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_HWHW -> "ghwhw" | G_HWHWH -> "ghwhwh" | G_HAHAH -> "ghahah"
      | G_HZHZ -> "ghzhz" | G_HZHAH -> "ghzhah"
      | G_HAHZ -> "ghahz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc"
      | G_Hthth -> "ghthth" | G_Htht -> "ghtht"
      | G_HHtt -> "ghhtt" | G_HHthth -> "ghhthth" | G_HHtht -> "ghhtht"
      | G_Psi0tt -> "gpsi0tt" | G_Psi0bb -> "gpsi0bb" 
      | G_Psi0cc -> "gpsi0cc" | G_Psi0tautau -> "gpsi0tautau"
      | G_Psi1tt -> "gpsi1tt" | G_Psi1bb -> "gpsi1bb" 
      | G_Psi1cc -> "gpsi1cc" | G_Psi1tautau -> "gpsi1tautau"
      | G_Psipq3 -> "gpsipq3" | G_Psipq2 -> "gpsipq2" | G_Psipl3 -> "gpsipl3"
      | G_Psi0tth -> "gpsi0tth" | G_Psi1tth -> "gpsi1tth"
      | G_Psipbth -> "gpsipbth"
      | G_Ethth -> "gethth" | G_Etht -> "getht"
      | G_Ett -> "gett" | G_Ebb -> "gebb"
      | G_HGaGa -> "ghgaga" | G_HGaZ -> "ghgaz"
      | G_EGaGa -> "geaa" | G_EGaZ -> "geaz" | G_EGlGl -> "gegg"
      | G_H3 -> "gh3" | G_H4 -> "gh4"
      | G_PsiWW -> "gpsiww" | G_PsiWHW -> "gpsiwhw" 
      | G_PsiZZ -> "gpsizz" | G_PsiZHZH -> "gpsizhzh" 
      | G_PsiZHZ -> "gpsizhz" | G_PsiZAH -> "gpsizah" 
      | G_PsiZHAH -> "gpsizhah" | G_PsiAHAH -> "gpsiahah"
      | G_PsiZW -> "gpsizw" | G_PsiZWH -> "gpsizwh" | G_PsiAHW -> "gpsiahw"
      | G_PsiAHWH -> "gpsiahwh" | G_PsiZHW -> "gpsizhw" 
      | G_PsiZHWH -> "gpsizhwh"
      | G_PsippWW -> "gpsippww" | G_PsippWHW -> "gpsippwhw" 
      | G_PsippWHWH -> "gpsippwhwh"
      | G_PsiHW -> "gpsihw" | G_PsiHWH -> "gpsihwh" 
      | G_Psi0W -> "gpsi0w" | G_Psi0WH -> "gpsi0wh"
      | G_Psi1W -> "gpsi1w" | G_Psi1WH -> "gpsi1wh" 
      | G_PsiPPW -> "gpsippw" | G_PsiPPWH -> "gpsippwh"
      | G_Psi1HAH -> "gpsihah" | G_Psi01AH -> "gpsi0ah" 
      | G_AHPsip -> "gahpsip" | G_Psi1HZ -> "gpsi1hz"
      | G_Psi1HZH -> "gpsi1hzh" | G_Psi01Z -> "gpsi01z" 
      | G_Psi01ZH -> "gpsi01zh" | G_ZPsip -> "gzpsip" 
      | G_ZPsipp -> "gzpsipp" | G_ZHPsipp -> "gzhpsipp"
      | G_HHAA -> "ghhaa" | G_HHWHW -> "ghhwhw" | G_HHZHZ -> "ghhzhz"
      | G_HHAHZ -> "ghhahz" | G_HHZHAH -> "ghhzhah"
      | G_HPsi0WW -> "ghpsi0ww" | G_HPsi0WHW -> "ghpsi0whw" 
      | G_HPsi0ZZ -> "ghpsi0zz" | G_HPsi0ZHZH -> "ghpsi0zhzh" 
      | G_HPsi0ZHZ -> "ghpsi0zhz" | G_HPsi0AHAH -> "ghpsi0ahah"
      | G_HPsi0ZAH -> "ghpsi0zah" | G_HPsi0ZHAH -> "ghpsi0zhah"
      | G_HPsipWA -> "ghpsipwa" | G_HPsipWHA -> "ghpsipwha" 
      | G_HPsipWZ -> "ghpsipwz" | G_HPsipWHZ -> "ghpsiwhz" 
      | G_HPsipWAH -> "ghpsipwah" | G_HPsipWHAH -> "ghpsipwhah" 
      | G_HPsipWZH -> "ghpsipwzh" | G_HPsipWHZH -> "ghpsipwhzh" 
      | G_HPsippWW -> "ghpsippww" | G_HPsippWHWH -> "ghpsippwhwh"
      | G_HPsippWHW -> "ghpsippwhw" | G_Psi00ZH -> "gpsi00zh" 
      | G_Psi00AH -> "gpsi00ah" | G_Psi00ZHAH -> "gpsi00zhah"
      | G_Psi0pWA -> "gpsi0pwa" | G_Psi0pWHA -> "gpsi0pwha" 
      | G_Psi0pWZ -> "gpsi0pwz" | G_Psi0pWHZ -> "gpsi0pwhz" 
      | G_Psi0pWAH -> "gpsi0pwah" | G_Psi0pWHAH -> "gpsi0pwhah" 
      | G_Psi0pWZH -> "gpsi0pwzh" | G_Psi0pWHZH -> "gpsi0pwhzh" 
      | G_Psi0ppWW -> "gpsi0ppww" | G_Psi0ppWHWH -> "gpsi0ppwhwh"
      | G_Psi0ppWHW -> "gpsi0ppwhw"
      | I_G_Psi0pWA -> "i_gpsi0pwa" | I_G_Psi0pWHA -> "i_gpsi0pwha" 
      | I_G_Psi0pWZ -> "i_gpsi0pwz" | I_G_Psi0pWHZ -> "i_gpsi0pwhz" 
      | I_G_Psi0pWAH -> "i_gpsi0pwah" | I_G_Psi0pWHAH -> "i_gpsi0pwhah" 
      | I_G_Psi0pWZH -> "i_gpsi0pwzh" | I_G_Psi0pWHZH -> "i_gpsi0pwhzh" 
      | I_G_Psi0ppWW -> "i_gpsi0ppww" | I_G_Psi0ppWHWH -> "i_gpsi0ppwhwh"
      | I_G_Psi0ppWHW -> "i_gpsi0ppwhw" 
      | G_PsippZZ -> "gpsippzz" | G_PsippZHZH -> "gpsippzhzh"
      | G_PsippAZ -> "gpsippaz" | G_PsippAAH -> "gpsippaah" 
      | G_PsippZAH -> "gpsippzah"
      | G_PsippWA -> "gpsippwa" | G_PsippWHA -> "gpsippwha" 
      | G_PsippWZ -> "gpsippwz" | G_PsippWHZ -> "gpsippwhz" 
      | G_PsippWAH -> "gpsippwah" | G_PsippWHAH -> "gpsippwhah"
      | G_PsippWZH -> "gpsippwzh" | G_PsippWHZH -> "gpsippwhzh"
      | G_PsiccZZ -> "gpsicczz" | G_PsiccAZ -> "gpsiccaz"
      | G_PsiccAAH -> "gpsiccaah" | G_PsiccZZH -> "gpsicczzh" 
      | G_PsiccAZH -> "gpsiccazh" | G_PsiccZAH -> "gpsicczah"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f
  end

module Simplest (Flags : BSM_flags) =
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

(* We do not introduce the Goldstones for the heavy vectors here. The heavy
   quarks are simply numerated by their generation, the assignments whether
   they are up- or down-type will be defined by the model. *)

    type flavor = L of int | N of int | U of int | D of int | QH of int      
        | NH of int | Wp | Wm | Ga | Z | Xp | Xm | X0 | Y0 | ZH 
        | Phip | Phim | Phi0 | H | Eta | Gl

    type gauge = unit

    let gauge_symbol () =
      failwith "Modellib_BSM.Simplest.gauge_symbol: internal error"

    let family n = [ L n; N n; U n; D n; QH n; NH n ]

(* Note that we add all heavy quarks, [U], [D], [C], [S], in order to have 
   both embeddings included. *)

    let external_flavors () =
      [ "1st Generation (incl. heavy)", ThoList.flatmap family [1; -1];
        "2nd Generation (incl. heavy)", ThoList.flatmap family [2; -2];
        "3rd Generation (incl. heavy)", ThoList.flatmap family [3; -3];
        "Gauge Bosons", [Ga; Z; Wp; Wm; Gl; Xp; Xm; X0; Y0; ZH];
        "Higgs", [H; Eta];
        "Goldstone Bosons", [Phip; Phim; Phi0] ]

    let flavors () = ThoList.flatmap snd (external_flavors ())

    let spinor n =
      if n >= 0 then
        Spinor
      else
        ConjSpinor

    let lorentz = function
      | L n -> spinor n | N n -> spinor n
      | U n -> spinor n | D n -> spinor n
      | QH n -> spinor n | NH n -> spinor n
      | Ga | Gl -> Vector
      | Wp | Wm | Z | Xp | Xm | X0 | Y0 | ZH -> Massive_Vector
      | _ -> Scalar

    let color = function 
      | U n -> Color.SUN (if n > 0 then 3 else -3)
      | D n -> Color.SUN  (if n > 0 then 3 else -3)
      | QH n -> Color.SUN  (if n > 0 then 3 else -3)
      | Gl -> Color.AdjSUN 3 
      | _ -> Color.Singlet

    let nc () = 3

    let prop_spinor n =
      if n >= 0 then
        Prop_Spinor
      else
        Prop_ConjSpinor

    let propagator = function
      | L n -> prop_spinor n | N n -> prop_spinor n
      | U n -> prop_spinor n | D n -> prop_spinor n
      | QH n -> prop_spinor n | NH n -> prop_spinor n
      | Ga | Gl -> Prop_Feynman
      | Wp | Wm | Z | Xp | Xm | X0 | Y0 | ZH -> Prop_Unitarity
      | Phip | Phim | Phi0 -> Only_Insertion
      | H | Eta -> Prop_Scalar

(* Optionally, ask for the fudge factor treatment for the widths of
   charged particles.  Currently, this only applies to $W^\pm$ and top. *)

    let width f =
      if !use_fudged_width then
        match f with
        | Wp | Wm | U 3 | U (-3) | QH _ | NH _ -> Fudged
        | _ -> !default_width
      else
        !default_width

    let goldstone = function
      | Wp -> Some (Phip, Coupling.Integer 1)
      | Wm -> Some (Phim, Coupling.Integer 1)
      | Z -> Some (Phi0, Coupling.Integer 1)
      | _ -> None

    let conjugate = function
      | L n -> L (-n) | N n -> N (-n)
      | U n -> U (-n) | D n -> D (-n)
      | QH n -> QH (-n) | NH n -> NH (-n)
      | Ga -> Ga | Gl -> Gl | Z -> Z 
      | Wp -> Wm | Wm -> Wp 
      | Xp -> Xm | Xm -> Xp | X0 -> X0 | Y0 -> Y0 | ZH -> ZH
      | Phip -> Phim | Phim -> Phip | Phi0 -> Phi0
      | H -> H | Eta -> Eta

    let fermion = function
      | L n -> if n > 0 then 1 else -1
      | N n -> if n > 0 then 1 else -1
      | U n -> if n > 0 then 1 else -1
      | D n -> if n > 0 then 1 else -1
      | QH n -> if n > 0 then 1 else -1
      | NH n -> if n > 0 then 1 else -1
      | Ga | Gl | Z | Wp | Wm | Xp | Xm | X0 | Y0 | ZH -> 0
      | _ -> 0 


    module Ch = Charges.QQ
    let ( // ) = Algebra.Small_Rational.make

    let charge = function
       | L n -> if n > 0 then -1//1 else  1//1
       | N n | NH n -> 0//1
       | U n -> if n > 0 then  2//3 else -2//3
       | QH 3 -> 2//3 | QH (-3) -> -2//3
       | QH (1|2) ->
          if Flags.anom_ferm_ass then
             2//3
          else
             -1//3
       | QH ((-1)|(-2)) ->
          if Flags.anom_ferm_ass then
             -2//3
          else
             1//3
       | QH n -> invalid_arg ("Simplest.charge: QH " ^ string_of_int n)  
       | D n -> if n > 0 then -1//3 else  1//3
       | Gl | Ga | Z | ZH | X0 | Y0 -> 0//1
       | Wp | Xp ->  1//1
       | Wm | Xm -> -1//1
       | H | Phi0 | Eta ->  0//1
       | Phip ->  1//1
       | Phim -> -1//1

    let lepton = function
       | L n | N n | NH n 
          -> if n > 0 then 1//1 else -1//1
       | U _ | D _ | _ -> 0//1

    let baryon = function
       | L _ | N _ -> 0//1
       | U n | D n | QH n 
          -> if n > 0 then 1//1 else -1//1
       | _ -> 0//1

    let charges f = 
      [ charge f; lepton f; baryon f]

    type constant =
      | Unit | Pi | Alpha_QED | Sin2thw
      | Sinthw | Costhw | E | G_weak | Vev | VHeavy
      | Supp | Supp2
      | Sinpsi | Cospsi | Atpsi | Sccs  (* Mixing angles of SU(2) *)
      | Q_lepton | Q_up | Q_down | Q_Z_up | G_CC | I_G_CC
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down
      | G_NC_X | G_NC_X_t | G_NC_Y | G_NC_Y_t | G_NC_H
      | G_NC_h_neutrino | G_NC_h_lepton | G_NC_h_up | G_NC_h_down
      | G_NC_h_top | G_NC_h_bot | G_NCH_N | G_NCH_U | G_NCH_D | G_NCHt
      | G_zhthth       
      | I_Q_W | I_G_ZWW | I_G_WWW
      | I_G_Z1 | I_G_Z2 | I_G_Z3 | I_G_Z4 | I_G_Z5 | I_G_Z6
      | I_Q_H | Gs | I_Gs | G2
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | I_Q_ZH 
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ | G_HHZZH 
      | G_heavy_HVV | G_heavy_HWW | G_heavy_HZZ | G_HHthth
      | G_Htt | G_Hbb | G_Hcc | G_Htautau | G_H3 | G_H4
      | G_Hthth | G_Htht | G_Ethth | G_Etht | G_Ett | G_Hqhq
      | G_Ebb | G_ZEH | G_ZHEH | G_Hgg
      | G_HGaGa | G_HGaZ | G_EGaGa | G_EGaZ | G_EGlGl 
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

(* \begin{dubious}
     The current abstract syntax for parameter dependencies is admittedly
     tedious. Later, there will be a parser for a convenient concrete syntax
     as a part of a concrete syntax for models.  But as these examples show,
     it should include simple functions.
   \end{dubious} *)


    let input_parameters =
      []

    let derived_parameters =
      [] 

    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

    let nc_coupling c t3 q =
      (Real_Array c,
       [Prod [g_over_2_costh; Diff (t3, Prod [Integer 2; q; Atom Sin2thw])];
        Prod [g_over_2_costh; t3]])

    let half = Quot (Integer 1, Integer 2)

    let derived_parameter_arrays =
      [ nc_coupling G_NC_neutrino half (Integer 0);
        nc_coupling G_NC_lepton (Neg half) (Integer (-1));
        nc_coupling G_NC_up half (Quot (Integer 2, Integer 3));
        nc_coupling G_NC_down (Neg half) (Quot (Integer (-1), Integer 3));
        nc_coupling G_NC_h_neutrino half (Integer 0);
        nc_coupling G_NC_h_lepton (Neg half) (Integer (-1));
        nc_coupling G_NC_h_up half (Quot (Integer 2, Integer 3));
        nc_coupling G_NC_h_down (Neg half) (Quot (Integer (-1), Integer 3)) ]


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

    let electromagnetic_currents n =
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);  
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]

    let color_currents n =
      [ ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs);
        ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs);
        ((QH (-n), Gl, QH n), FBF ((-1), Psibar, V, Psi), Gs)]

    let neutral_currents n =
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

    let xy_currents = 
      ThoList.flatmap 
        (fun n -> [ ((N (-n), X0, N n), FBF ((-1), Psibar, VL, Psi), G_NC_X);
                    ((L (-n), Xm, N n), FBF ((-1), Psibar, VL, Psi), G_NC_X);
                    ((N (-n), Xp, L n), FBF ((-1), Psibar, VL, Psi), G_NC_X);
                    ((N (-n), Y0, N n), FBF ((-1), Psibar, VL, Psi), G_NC_Y);
                    ((NH (-n), X0, N n), FBF ((-1), Psibar, VL, Psi), G_CC);
                    ((N (-n), X0, NH n), FBF ((-1), Psibar, VL, Psi), G_CC);
                    ((NH (-n), Y0, N n), FBF ((-1), Psibar, VL, Psi), I_G_CC);
                    ((N (-n), Y0, NH n), FBF ((-1), Psibar, VL, Psi), I_G_CC);
                    ((L (-n), Xm, NH n), FBF ((-1), Psibar, VL, Psi), G_CC);
                    ((NH (-n), Xp, L n), FBF ((-1), Psibar, VL, Psi), G_CC)])
        [1;2;3]
      @ 
        [ ((U (-3), X0, U 3), FBF (1, Psibar, VL, Psi), G_NC_X_t);
          ((U (-3), Y0, U 3), FBF (1, Psibar, VL, Psi), G_NC_Y_t);
          ((U (-3), X0, QH 3), FBF (1, Psibar, VL, Psi), G_CC);
          ((QH (-3), X0, U 3), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-3), Y0, QH 3), FBF (1, Psibar, VL, Psi), I_G_CC);
          ((QH (-3), Y0, U 3), FBF (1, Psibar, VL, Psi), I_G_CC);
          ((D (-3), Xm, U 3), FBF (1, Psibar, VL, Psi), G_NC_X_t);
          ((U (-3), Xp, D 3), FBF (1, Psibar, VL, Psi), G_NC_X_t);
          ((D (-3), Xm, QH 3), FBF (1, Psibar, VL, Psi), G_CC);
          ((QH (-3), Xp, D 3), FBF (1, Psibar, VL, Psi), G_CC);
          ((QH (-3), Wp, D 3), FBF (1, Psibar, VL, Psi), G_NC_X_t);
          ((D (-3), Wm, QH 3), FBF (1, Psibar, VL, Psi), G_NC_X_t);
          ((QH (-3), Z, U 3), FBF (1, Psibar, VL, Psi), G_NCHt);
          ((U (-3), Z, QH 3), FBF (1, Psibar, VL, Psi), G_NCHt)]
      @
        ThoList.flatmap
          (fun n -> 
            if Flags.anom_ferm_ass then
              [ ((U (-n), X0, U n), FBF ((-1), Psibar, VL, Psi), G_NC_X);
                ((U (-n), Y0, U n), FBF ((-1), Psibar, VL, Psi), G_NC_Y);
                ((D (-n), Xm, U n), FBF ((-1), Psibar, VL, Psi), G_NC_X);
                ((U (-n), Xp, D n), FBF ((-1), Psibar, VL, Psi), G_NC_X);
                ((QH (-n), X0, U n), FBF ((-1), Psibar, VL, Psi), G_CC);
                ((U (-n), X0, QH n), FBF ((-1), Psibar, VL, Psi), G_CC);
                ((QH (-n), Y0, U n), FBF ((-1), Psibar, VL, Psi), I_G_CC);
                ((U (-n), Y0, QH n), FBF ((-1), Psibar, VL, Psi), I_G_CC);
                ((D (-n), Xm, QH n), FBF ((-1), Psibar, VL, Psi), G_CC);
                ((QH (-n), Xp, D n), FBF ((-1), Psibar, VL, Psi), G_CC);
                ((QH (-n), Wp, D n), FBF ((-1), Psibar, VL, Psi), G_NC_X);
                ((D (-n), Wm, QH n), FBF ((-1), Psibar, VL, Psi), G_NC_X);
                ((QH (-n), Z, U n), FBF (1, Psibar, VL, Psi), G_NC_H);
                ((U (-n), Z, QH n), FBF (1, Psibar, VL, Psi), G_NC_H)]
            else
              [ ((D (-n), X0, D n), FBF (1, Psibar, VL, Psi), G_NC_X);
                ((D (-n), Y0, D n), FBF (1, Psibar, VL, Psi), G_NC_Y);
                ((D (-n), Xm, U n), FBF (1, Psibar, VL, Psi), G_NC_X);
                ((U (-n), Xp, D n), FBF (1, Psibar, VL, Psi), G_NC_X);
                ((QH (-n), X0, D n), FBF ((-1), Psibar, VL, Psi), G_CC);
                ((D (-n), X0, QH n), FBF ((-1), Psibar, VL, Psi), G_CC);
                ((QH (-n), Y0, D n), FBF ((-1), Psibar, VL, Psi), I_G_CC);
                ((D (-n), Y0, QH n), FBF ((-1), Psibar, VL, Psi), I_G_CC);
                ((QH (-n), Xm, U n), FBF (1, Psibar, VL, Psi), G_CC);
                ((U (-n), Xp, QH n), FBF (1, Psibar, VL, Psi), G_CC);
                ((QH (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_NC_X);
                ((U (-n), Wp, QH n), FBF (1, Psibar, VL, Psi), G_NC_X);
                ((QH (-n), Z, D n), FBF (1, Psibar, VL, Psi), G_NC_H);
                ((D (-n), Z, QH n), FBF (1, Psibar, VL, Psi), G_NC_H)])
          [1; 2]
         

(* The sign of this coupling is just the one of the T3, being -(1/2) for
   leptons and down quarks, and +(1/2) for neutrinos and up quarks. *)

    let neutral_heavy_currents n =
        [ ((L (-n), ZH, L n), FBF (1, Psibar, VLR, Psi), G_NC_h_lepton);
          ((N (-n), ZH, N n), FBF ((-1), Psibar, VLR, Psi), G_NC_h_neutrino);
          ((U (-n), ZH, U n), FBF ((-1), Psibar, VLR, Psi), (if n = 3 then
                                   G_NC_h_top else G_NC_h_up));
          ((D (-n), ZH, D n), FBF (1, Psibar, VLR, Psi), (if n = 3 then 
                                   G_NC_h_bot else G_NC_h_down));
          ((NH (-n), ZH, NH n), FBF (1, Psibar, VLR, Psi), G_NCH_N);
          ((QH (-n), ZH, QH n), FBF (1, Psibar, VLR, Psi), (if n = 3 then
             G_NCH_U else if Flags.anom_ferm_ass then G_NCH_U else G_NCH_D))]
                                    

    let heavy_currents n = 
      [ ((QH (-n), Ga, QH n), FBF (1, Psibar, V, Psi), (if n=3 then Q_up else
            if Flags.anom_ferm_ass then Q_up else Q_down))]

    let charged_currents n =
      [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
        ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC);
        ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
        ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ] 
        
    let yukawa =
      [ ((U (-3), H, U 3), FBF (1, Psibar, S, Psi), G_Htt);
        ((D (-3), H, D 3), FBF (1, Psibar, S, Psi), G_Hbb);
        ((U (-2), H, U 2), FBF (1, Psibar, S, Psi), G_Hcc);
        ((L (-3), H, L 3), FBF (1, Psibar, S, Psi), G_Htautau) ]

    let yukawa_add = 
      [ ((QH (-3), H, U 3), FBF (1, Psibar, SL, Psi), G_Htht);
        ((U (-3), H, QH 3), FBF (1, Psibar, SR, Psi), G_Htht);
        ((QH (-3), Eta, U 3), FBF (1, Psibar, SR, Psi), G_Etht);
        ((U (-3), Eta, QH 3), FBF (1, Psibar, SL, Psi), G_Etht);
        ((D (-3), Eta, D 3), FBF (1, Psibar, P, Psi), G_Ebb);
        ((U (-3), Eta, U 3), FBF (1, Psibar, P, Psi), G_Ett)]
        @ 
          ThoList.flatmap
            (fun n -> 
          if Flags.anom_ferm_ass then
      [ ((QH (-n), H, U n), FBF (1, Psibar, SL, Psi), G_Hqhq);
        ((U (-n), H, QH n), FBF (1, Psibar, SR, Psi), G_Hqhq)]
          else
      [ ((QH (-n), H, D n), FBF (1, Psibar, SL, Psi), G_Hqhq);
        ((D (-n), H, QH n), FBF (1, Psibar, SR, Psi), G_Hqhq)])
            [1;2]


    let standard_triple_gauge =
      [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
        ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
        ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs)]
        
    let heavy_triple_gauge =
      [ ((Ga, Xm, Xp), Gauge_Gauge_Gauge 1, I_Q_W);
        ((Z, Xm, Xp), Gauge_Gauge_Gauge 1,  I_Q_ZH);
        ((Z, X0, Y0), Gauge_Gauge_Gauge 1,  I_G_Z1);
        ((ZH, X0, Y0), Gauge_Gauge_Gauge 1, I_G_Z2);
        ((Y0, Wm, Xp), Gauge_Gauge_Gauge 1, I_G_Z3);
        ((Y0, Wp, Xm), Gauge_Gauge_Gauge (-1), I_G_Z3);
        ((X0, Wm, Xp), Gauge_Gauge_Gauge 1, I_G_Z4);
        ((X0, Wp, Xm), Gauge_Gauge_Gauge 1, I_G_Z4);
        ((ZH, Xm, Xp), Gauge_Gauge_Gauge 1, I_G_Z5);
        ((ZH, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_Z6)]

    let triple_gauge =
      standard_triple_gauge @ heavy_triple_gauge
                                
    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let standard_quartic_gauge =
      [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
        (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
        (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
        (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
        (Gl, Gl, Gl, Gl), gauge4, G2]
        
    let heavy_quartic_gauge = 
      []
        

    let quartic_gauge =
      standard_quartic_gauge @ heavy_quartic_gauge

    let standard_gauge_higgs' =
      [ ((H, Wp, Wm), Scalar_Vector_Vector 1, G_HWW);
        ((H, Z, Z), Scalar_Vector_Vector 1, G_HZZ) ]

    let heavy_gauge_higgs = 
      [ ((H, Wp, Xm), Scalar_Vector_Vector 1, G_heavy_HWW);
        ((H, Wm, Xp), Scalar_Vector_Vector 1, G_heavy_HWW);
        ((H, Z, X0),  Scalar_Vector_Vector 1,  G_heavy_HVV);
        ((H, ZH, X0), Scalar_Vector_Vector 1, G_heavy_HVV)]

    let standard_gauge_higgs = 
      standard_gauge_higgs' @ heavy_gauge_higgs 

    let standard_gauge_higgs4 =
      [ (H, H, Wp, Wm), Scalar2_Vector2 1, G_HHWW;
        (H, H, Z, Z), Scalar2_Vector2 1, G_HHZZ ]

    let heavy_gauge_higgs4 =
      [ (H, H, Z, ZH), Scalar2_Vector2 1, G_HHZZH;
        (H, H, Xp, Xm), Scalar2_Vector2 (-1), G_HHWW;
        (H, H, ZH, ZH), Scalar2_Vector2 (-1), G_HHZZ ]

    let standard_higgs =
      [ (H, H, H), Scalar_Scalar_Scalar 1, G_H3 ]
        
   let anomaly_higgs = 
      [ (Eta, Gl, Gl), Dim5_Scalar_Gauge2_Skew 1, G_EGlGl;
        (Eta, Ga, Ga), Dim5_Scalar_Gauge2_Skew 1, G_EGaGa; 
        (Eta, Ga, Z), Dim5_Scalar_Gauge2_Skew 1, G_EGaZ  ] 
(*    @ [ (H, Ga, Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
        (H, Ga, Z), Dim5_Scalar_Gauge2 1, G_HGaZ ]   *)

    let standard_higgs4 =
      [ (H, H, H, H), Scalar4 1, G_H4 ]

    let gauge_higgs =
        standard_gauge_higgs

    let gauge_higgs4 =
        standard_gauge_higgs4 @ heavy_gauge_higgs4

    let higgs =
        standard_higgs

    let eta_higgs_gauge =
      [ (Z, Eta, H), Vector_Scalar_Scalar 1, G_ZEH;
        (ZH, Eta, H), Vector_Scalar_Scalar 1, G_ZHEH;
        (X0, Eta, H), Vector_Scalar_Scalar 1, G_CC ]    

    let top_quartic = 
      [ ((QH (-3), H, H, QH 3), GBBG (1, Psibar, S2, Psi), G_HHthth)]

    let higgs4 =
        standard_higgs4

    let goldstone_vertices =
      [ ((Phi0, Wm, Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((Phip, Ga, Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((Phip, Z, Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((Phim, Wp, Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((Phim, Wp, Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @ 
       ThoList.flatmap color_currents [1;2;3] @ 
       ThoList.flatmap neutral_currents [1;2;3] @
       ThoList.flatmap neutral_heavy_currents [1;2;3] @       
       ThoList.flatmap heavy_currents [1;2;3] @       
       ThoList.flatmap charged_currents [1;2;3] @
       xy_currents @ anomaly_higgs @ 
       eta_higgs_gauge @
       yukawa @ yukawa_add @ 
       triple_gauge @ 
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
      | "e-" -> L 1 | "e+" -> L (-1)
      | "mu-" -> L 2 | "mu+" -> L (-2)
      | "tau-" -> L 3 | "tau+" -> L (-3)
      | "nue" -> N 1 | "nuebar" -> N (-1)
      | "numu" -> N 2 | "numubar" -> N (-2)
      | "nutau" -> N 3 | "nutaubar" -> N (-3)
      | "nh1" -> NH 1 | "nh1bar" -> NH (-1)
      | "nh2" -> NH 2 | "nh2bar" -> NH (-2)
      | "nh3" -> NH 3 | "nh3bar" -> NH (-3)
      | "u" -> U 1 | "ubar" -> U (-1)
      | "c" -> U 2 | "cbar" -> U (-2)
      | "t" -> U 3 | "tbar" -> U (-3)
      | "d" -> D 1 | "dbar" -> D (-1)
      | "s" -> D 2 | "sbar" -> D (-2)
      | "b" -> D 3 | "bbar" -> D (-3)            
      | "uh" -> if Flags.anom_ferm_ass then QH 1 else invalid_arg
          "Modellib_BSM.Simplest.flavor_of_string"
      | "dh" -> if Flags.anom_ferm_ass then invalid_arg
            "Modellib_BSM.Simplest.flavor_of_string" else QH 1
      | "uhbar" -> if Flags.anom_ferm_ass then QH (-1) else invalid_arg
          "Modellib_BSM.Simplest.flavor_of_string"
      | "dhbar" -> if Flags.anom_ferm_ass then invalid_arg
            "Modellib_BSM.Simplest.flavor_of_string" else QH (-1)
      | "ch" -> if Flags.anom_ferm_ass then QH 2 else invalid_arg
          "Modellib_BSM.Simplest.flavor_of_string"
      | "sh" -> if Flags.anom_ferm_ass then invalid_arg
            "Modellib_BSM.Simplest.flavor_of_string" else QH 2
      | "chbar" -> if Flags.anom_ferm_ass then QH (-2) else invalid_arg
          "Modellib_BSM.Simplest.flavor_of_string"
      | "shbar" -> if Flags.anom_ferm_ass then invalid_arg
            "Modellib_BSM.Simplest.flavor_of_string" else QH (-2)
      | "th" -> QH 3 | "thbar" -> QH (-3)
      | "eta" | "Eta" -> Eta
      | "A" -> Ga | "Z" | "Z0" -> Z | "g" | "gl" -> Gl
      | "ZH" | "ZH0" | "Zh" | "Zh0" -> ZH
      | "W+" -> Wp | "W-" -> Wm
      | "X+" -> Xp | "X-" -> Xm
      | "X0" -> X0 | "Y0" -> Y0
      | "H" -> H
      | _ -> invalid_arg "Modellib_BSM.Simplest.flavor_of_string" 

    let flavor_to_string = function
      | L 1 -> "e-" | L (-1) -> "e+"
      | L 2 -> "mu-" | L (-2) -> "mu+"
      | L 3 -> "tau-" | L (-3) -> "tau+"
      | L _ -> invalid_arg
            "Modellib_BSM.Simplest.flavor_to_string: invalid lepton"
      | N 1 -> "nue" | N (-1) -> "nuebar"
      | N 2 -> "numu" | N (-2) -> "numubar"
      | N 3 -> "nutau" | N (-3) -> "nutaubar"
      | N _ -> invalid_arg
            "Modellib_BSM.Simplest.flavor_to_string: invalid neutrino"
      | U 1 -> "u" | U (-1) -> "ubar"
      | U 2 -> "c" | U (-2) -> "cbar"
      | U 3 -> "t" | U (-3) -> "tbar"
      | U _ -> invalid_arg
            "Modellib_BSM.Simplest.flavor_to_string: invalid up type quark"
      | D 1 -> "d" | D (-1) -> "dbar"
      | D 2 -> "s" | D (-2) -> "sbar"
      | D 3 -> "b" | D (-3) -> "bbar"
      | D _ -> invalid_arg
            "Modellib_BSM.Simplest.flavor_to_string: invalid down type quark"
      | QH 1 -> if Flags.anom_ferm_ass then "uh" else "dh"
      | QH 2 -> if Flags.anom_ferm_ass then "ch" else "sh"
      | QH 3 -> "th" 
      | QH (-1) -> if Flags.anom_ferm_ass then "uhbar" else "dhbar"
      | QH (-2) -> if Flags.anom_ferm_ass then "chbar" else "shbar"
      | QH (-3) -> "thbar"
      | QH _ -> invalid_arg 
            "Modellib_BSM.Simplest.flavor_to_string: invalid heavy quark"
      | NH n when n > 0 -> "nh" ^ string_of_int n
      | NH n -> "nh" ^ string_of_int (abs n) ^ "bar" 
      | Ga -> "A" | Z -> "Z" | Gl -> "gl"
      | Wp -> "W+" | Wm -> "W-"
      | Xp -> "X+" | Xm -> "X-" | X0 -> "X0" | Y0 -> "Y0" | ZH -> "ZH" 
      | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
      | H -> "H" | Eta -> "Eta"

    let flavor_to_TeX = function
      | L 1 -> "e^-" | L (-1) -> "\\e^+"
      | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
      | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
      | L _ -> invalid_arg
            "Modellib_BSM.Simplest.flavor_to_TeX: invalid lepton"
      | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
      | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
      | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
      | N _ -> invalid_arg
            "Modellib_BSM.Simplest.flavor_to_TeX: invalid neutrino"
      | U 1 -> "u" | U (-1) -> "\\bar{u}"
      | U 2 -> "c" | U (-2) -> "\\bar{c}"
      | U 3 -> "t" | U (-3) -> "\\bar{t}"
      | U _ -> invalid_arg
            "Modellib_BSM.Simplest.flavor_to_TeX: invalid up type quark"
      | D 1 -> "d" | D (-1) -> "\\bar{d}"
      | D 2 -> "s" | D (-2) -> "\\bar{s}"
      | D 3 -> "b" | D (-3) -> "\\bar{b}"
      | D _ -> invalid_arg
            "Modellib_BSM.Simplest.flavor_to_TeX: invalid down type quark"
      | QH 1 -> if Flags.anom_ferm_ass then "U" else "D"
      | QH 2 -> if Flags.anom_ferm_ass then "C" else "S"
      | QH 3 -> "T" 
      | QH (-1) -> if Flags.anom_ferm_ass then "\\bar{U}" else "\\bar{D}"
      | QH (-2) -> if Flags.anom_ferm_ass then "\\bar{C}" else "\\bar{S}"
      | QH (-3) -> "thbar"
      | QH _ -> invalid_arg 
            "Modellib_BSM.Simplest.flavor_to_TeX: invalid heavy quark"
      | NH n when n > 0 -> "N_" ^ string_of_int n
      | NH n -> "\\bar{N}_" ^ string_of_int (abs n)
      | Ga -> "\\gamma" | Z -> "Z" | Gl -> "g"
      | Wp -> "W^+" | Wm -> "W^-"
      | Xp -> "X^+" | Xm -> "X^-" | X0 -> "X^0" | Y0 -> "Y^0" | ZH -> "Z_H" 
      | Phip -> "\\phi^+" | Phim -> "\\phi^-" | Phi0 -> "\\phi^0" 
      | H -> "H" | Eta -> "\\eta"

    let flavor_symbol = function
      | L n when n > 0 -> "l" ^ string_of_int n
      | L n -> "l" ^ string_of_int (abs n) ^ "b"
      | N n when n > 0 -> "n" ^ string_of_int n
      | N n -> "n" ^ string_of_int (abs n) ^ "b"
      | U n when n > 0 -> "u" ^ string_of_int n
      | U n -> "u" ^ string_of_int (abs n) ^ "b"
      | D n when n > 0 ->  "d" ^ string_of_int n
      | D n -> "d" ^ string_of_int (abs n) ^ "b"
      | NH n when n > 0 -> "nh" ^ string_of_int n
      | NH n -> "nh" ^ string_of_int (abs n) ^ "b"
      | QH n when n > 0 -> "qh" ^ string_of_int n
      | QH n -> "qh" ^ string_of_int (abs n) ^ "b"
      | Ga -> "a" | Z -> "z" | Gl -> "gl"
      | Wp -> "wp" | Wm -> "wm"
      | Xp -> "xp" | Xm -> "xm" | X0 -> "x0" | Y0 -> "y0" | ZH -> "zh"
      | Phip -> "pp" | Phim -> "pm" | Phi0 -> "p0" 
      | H -> "h" | Eta -> "eta"

(* There are PDG numbers for Z', Z'', W', 32-34, respectively.
   We just introduce a number 38 for Y0 as a Z'''.
   As well, there is the number 8 for a t'. But we cheat a little bit and 
   take the number 35 which is reserved for a heavy scalar Higgs for the 
   Eta scalar. 

   We abuse notation for the heavy quarks and take the PDG code for their 
   SUSY partners!!! (What about an update of the PDG numbering scheme?)
   Thereby we take only those for up-type (s)quarks. The heavy neutrinos get 
   the numbers of the sneutrinos. 
*)  

    let pdg = function
      | L n when n > 0 -> 9 + 2*n
      | L n -> - 9 + 2*n
      | N n when n > 0 -> 10 + 2*n
      | N n -> - 10 + 2*n
      | U n when n > 0 -> 2*n
      | U n -> 2*n
      | D n when n > 0 -> - 1 + 2*n
      | D n -> 1 + 2*n
      | NH n when n > 0 -> 1000010 + 2*n 
      | NH n -> - 1000010 + 2*n
      | QH 3 -> 1000006 
      | QH (-3) -> - 1000006
      | QH n when n > 0 -> if Flags.anom_ferm_ass then 
          1000000 + 2*n   else   999999 + 2*n
      | QH n -> if Flags.anom_ferm_ass then
          - 1000000 + 2*n   else   - 999999 + 2*n
      | Gl -> 21 
      | Ga -> 22 | Z -> 23
      | Wp -> 24 | Wm -> (-24)
      | Xp -> 34 | Xm -> (-34) | ZH -> 32 | X0 -> 33 | Y0 -> 38
      | Phip | Phim -> 27 | Phi0 -> 26
      | H -> 25 | Eta -> 36

(* As in the case of SUSY we introduce an internal dummy pdf code in order
   to have manageable arrays. Heavy neutrinos get numbers 41,43,45, while the
   heavy quarks have the numbers 40,42,44. I take them all as up type
   here. 
 *)

    let pdg_mw = function
      | L n when n > 0 -> 9 + 2*n
      | L n -> - 9 + 2*n
      | N n when n > 0 -> 10 + 2*n
      | N n -> - 10 + 2*n
      | U n when n > 0 -> 2*n
      | U n -> 2*n
      | D n when n > 0 -> - 1 + 2*n
      | D n -> 1 + 2*n
      | NH n when n > 0 -> 39 + 2*n 
      | NH n -> - 39 + 2*n
      | QH n when n > 0 -> 38 + 2*n
      | QH n -> - 38 + 2*n
      | Gl -> 21
      | Ga -> 22 | Z -> 23
      | Wp -> 24 | Wm -> (-24)
      | Xp -> 34 | Xm -> (-34) | ZH -> 32 | X0 -> 33 | Y0 -> 38
      | Phip | Phim -> 27 | Phi0 -> 26
      | H -> 25 | Eta -> 36

    let mass_symbol f = 
      "mass(" ^ string_of_int (abs (pdg_mw f)) ^ ")"

    let width_symbol f =
      "width(" ^ string_of_int (abs (pdg_mw f)) ^ ")"

    let constant_symbol = function
      | Unit -> "unit" | Pi -> "PI" | VHeavy -> "vheavy"
      | Alpha_QED -> "alpha" | E -> "e" | G_weak -> "g" | Vev -> "vev"
      | Sin2thw -> "sin2thw" | Sinthw -> "sinthw" | Costhw -> "costhw"
      | Sinpsi -> "sinpsi" | Cospsi -> "cospsi" 
      | Atpsi -> "atpsi" | Sccs -> "sccs"
      | Supp -> "vF" | Supp2 -> "v2F2" 
      | Q_lepton -> "qlep" | Q_up -> "qup" | Q_down -> "qdwn"
      | Q_Z_up -> "qzup" 
      | G_zhthth -> "gzhthth" 
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_NC_X -> "gncx" | G_NC_X_t -> "gncxt"
      | G_NC_Y -> "gncy" | G_NC_Y_t -> "gncyt" | G_NC_H -> "gnch"
      | G_CC -> "gcc" | I_G_CC -> "i_gcc"
      | G_NC_h_lepton -> "gnchlep" | G_NC_h_neutrino -> "gnchneu"
      | G_NC_h_up -> "gnchup" | G_NC_h_down -> "gnchdwn"   
      | G_NC_h_top -> "gnchtop" | G_NC_h_bot -> "gnchbot"
      | G_NCH_N -> "gnchn" | G_NCH_U -> "gnchu" | G_NCH_D -> "gnchd"
      | G_NCHt -> "gncht"
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" | I_G_WWW -> "igwww"
      | I_Q_H -> "iqh" | I_Q_ZH -> "iqzh"
      | I_G_Z1 -> "igz1" | I_G_Z2 -> "igz2" | I_G_Z3 -> "igz3" 
      | I_G_Z4 -> "igz4" | I_G_Z5 -> "igz5" | I_G_Z6 -> "igz6"
      | G_HHthth -> "ghhthth" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_heavy_HVV -> "ghyhvv"
      | G_heavy_HWW -> "ghyhww"
      | G_heavy_HZZ -> "ghyhzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_HHZZH -> "ghhzzh" 
      | G_Hgg -> "ghgg"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc"
      | G_Hthth -> "ghthth" | G_Htht -> "ghtht"
      | G_Hqhq -> "ghqhq"
      | G_Ethth -> "gethth" | G_Etht -> "getht"
      | G_Ett -> "gett" | G_Ebb -> "gebb"
      | G_HGaGa -> "ghgaga" | G_HGaZ -> "ghgaz"
      | G_EGaGa -> "geaa" | G_EGaZ -> "geaz" | G_EGlGl -> "gegg"
      | G_ZEH -> "gzeh" | G_ZHEH -> "gzheh" 
      | G_H3 -> "gh3" | G_H4 -> "gh4"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
  end

module Xdim (Flags : BSM_flags) =
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

    type matter_field = L of int | N of int | U of int | D of int
    type gauge_boson = Ga | Wp | Wm | Z | Gl
    type other = Phip | Phim | Phi0 | H | Grav
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
      failwith "Modellib_BSM.Xdim.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; N n; U n; D n ]

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Gauge Bosons", List.map gauge_boson [Ga; Z; Wp; Wm; Gl];
        "Higgs", List.map other [H];
        "Graviton", List.map other [Grav];
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
          begin match f with 
          | Grav -> Tensor_2
          | _ -> Scalar
          end

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
          | H -> Prop_Scalar
          | Grav -> Prop_Tensor_2
          end

(* Optionally, ask for the fudge factor treatment for the widths of
   charged particles.  Currently, this only applies to $W^\pm$ and top. *)

    let width f =
      if !use_fudged_width then
        match f with
        | G Wp | G Wm | M (U 3) | M (U (-3)) | O Grav -> Fudged
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
          | H -> H | Grav -> Grav
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

    module Ch = Charges.QQ
    let ( // ) = Algebra.Small_Rational.make

    let generation' = function
      |  1 -> [ 1//1;  0//1;  0//1]
      |  2 -> [ 0//1;  1//1;  0//1]
      |  3 -> [ 0//1;  0//1;  1//1]
      | -1 -> [-1//1;  0//1;  0//1]
      | -2 -> [ 0//1; -1//1;  0//1]
      | -3 -> [ 0//1;  0//1; -1//1]
      |  n -> invalid_arg ("Xdim.generation': " ^ string_of_int n)

    let generation f =
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
          | H | Phi0 | Grav ->  0//1
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
      | Q_lepton | Q_up | Q_down | G_CC | G_CCQ of int*int
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down
      | Gs | I_Gs | G2
      | I_Q_W | I_G_ZWW
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ
      | G_Htt | G_Hbb | G_Hcc | G_Htautau | G_H3 | G_H4
      | G_HGaZ | G_HGaGa | G_Hgg | G_Grav
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let input_parameters =
      []

    let derived_parameters =
      []

    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

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

    let mgm ((m1, g, m2), fbf, c) = ((M m1, G g, M m2), fbf, c)
    let mom ((m1, o, m2), fbf, c) = ((M m1, O o, M m2), fbf, c)

    let electromagnetic_currents n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]
       
    let neutral_currents n =
      List.map mgm
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

    let color_currents n =
      List.map mgm
        [ ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs);
          ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs) ]

    let charged_currents n = 
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC);
          ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let gravity_currents n = 
      List.map mom
        [ ((L (-n), Grav, L n), Graviton_Spinor_Spinor 1, G_Grav);
          ((N (-n), Grav, N n), Graviton_Spinor_Spinor 1, G_Grav);
          ((U (-n), Grav, U n), Graviton_Spinor_Spinor 1, G_Grav);
          ((D (-n), Grav, D n), Graviton_Spinor_Spinor 1, G_Grav) ] 

    let yukawa =
      List.map mom 
        [ ((U (-3), H, U 3), FBF (1, Psibar, S, Psi), G_Htt);
          ((D (-3), H, D 3), FBF (1, Psibar, S, Psi), G_Hbb);
          ((U (-2), H, U 2), FBF (1, Psibar, S, Psi), G_Hcc);
          ((L (-3), H, L 3), FBF (1, Psibar, S, Psi), G_Htautau) ]

    let tgc ((g1, g2, g3), t, c) = ((G g1, G g2, G g3), t, c)

    let standard_triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs) ]

    let triple_gauge =
        standard_triple_gauge

    let qgc ((g1, g2, g3, g4), t, c) = ((G g1, G g2, G g3, G g4), t, c)

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let standard_quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          (Gl, Gl, Gl, Gl), gauge4, G2] 

    let quartic_gauge =
      standard_quartic_gauge

    let gravity_gauge = 
      [ (O Grav, G Z, G Z), Graviton_Vector_Vector 1, G_Grav;
        (O Grav, G Wp, G Wm), Graviton_Vector_Vector 1, G_Grav;
        (O Grav, G Ga, G Ga), Graviton_Vector_Vector 1, G_Grav;
        (O Grav, G Gl, G Gl), Graviton_Vector_Vector 1, G_Grav ]

    let standard_gauge_higgs =
      [ ((O H, G Wp, G Wm), Scalar_Vector_Vector 1, G_HWW);
        ((O H, G Z, G Z), Scalar_Vector_Vector 1, G_HZZ) ]

    let standard_gauge_higgs4 =
      [ (O H, O H, G Wp, G Wm), Scalar2_Vector2 1, G_HHWW;
        (O H, O H, G Z, G Z), Scalar2_Vector2 1, G_HHZZ ]
       
    let standard_higgs =
      [ (O H, O H, O H), Scalar_Scalar_Scalar 1, G_H3 ]

    let standard_higgs4 =
      [ (O H, O H, O H, O H), Scalar4 1, G_H4 ]

    let gravity_higgs = 
      [ (O Grav, O H, O H), Graviton_Scalar_Scalar 1, G_Grav]

    let anomalous_gauge_higgs =
      []

    let anomalous_gauge_higgs4 =
      []

    let anomalous_higgs =
      []

    let anomaly_higgs = 
      [ (O H, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
        (O H, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_HGaZ;
        (O H, G Gl, G Gl), Dim5_Scalar_Gauge2 1, G_Hgg ]

    let anomalous_higgs4 =
      []

    let gauge_higgs =
        standard_gauge_higgs

    let gauge_higgs4 =
        standard_gauge_higgs4

    let higgs =
        standard_higgs @ gravity_higgs

    let higgs4 =
        standard_higgs4

    let goldstone_vertices =
      [ ((O Phi0, G Wm, G Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phip, G Ga, G Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phip, G Z, G Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phim, G Wp, G Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phim, G Wp, G Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @
       ThoList.flatmap neutral_currents [1;2;3] @
       ThoList.flatmap color_currents [1;2;3] @
       ThoList.flatmap charged_currents [1;2;3] @
       ThoList.flatmap gravity_currents [1;2;3] @
       yukawa @ triple_gauge @ gravity_gauge @
       gauge_higgs @ higgs @ anomaly_higgs 
       @ goldstone_vertices)

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
      | "GG" -> O Grav
      | _ -> invalid_arg "Modellib_BSM.Xdim.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib_BSM.Xdim.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib_BSM.Xdim.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib_BSM.Xdim.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_BSM.Xdim.flavor_to_string: invalid down type quark"
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
          | H -> "H" | Grav -> "GG"
          end

    let flavor_to_TeX = function
      | M f ->
          begin match f with
          | L 1 -> "e^-" | L (-1) -> "e^+"
          | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
          | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
          | L _ -> invalid_arg
                "Modellib_BSM.Xdim.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib_BSM.Xdim.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib_BSM.Xdim.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Modellib_BSM.Xdim.flavor_to_TeX: invalid down type quark"
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
          | H -> "H" | Grav -> "G"
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
          | H -> "h" | Grav -> "gv"
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
          | H -> 25 | Grav -> 39
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
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | G_CC -> "gcc"
      | G_CCQ (n1,n2) -> "gccq" ^ string_of_int n1 ^ string_of_int n2
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc"
      | G_HGaZ -> "ghgaz" | G_HGaGa -> "ghgaga" | G_Hgg -> "ghgg"
      | G_H3 -> "gh3" | G_H4 -> "gh4" | G_Grav -> "ggrav"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f

  end
 
module UED (Flags : BSM_flags) =
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

    type matter_field = L of int | N of int | U of int | D of int 
          | L_K1_L of int | L_K1_R of int | N_K1 of int
          | L_K2_L of int | L_K2_R of int | N_K2 of int
          | U_K1_L of int | U_K2_L of int | D_K1_L of int | D_K2_L of int
          | U_K1_R of int | U_K2_R of int | D_K1_R of int | D_K2_R of int
    type gauge_boson = Ga | Wp | Wm | Z | Gl | Gl_K1 | Gl_K2
          | B1 | B2 | Z1 | Z2 | Wp1 | Wm1 | Wp2 | Wm2
    type other = Phip | Phim | Phi0 | H | H1up | H1um 
          | H1dp | H1dm | H2up |H2um | H2dp |H2dm  
          | Grav
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
      failwith "Modellib_BSM.UED.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; N n; U n; D n; L_K1_L n;
      L_K1_R n; L_K2_L n; L_K2_R n; N_K1 n; N_K2 n; U_K1_L n; U_K2_L n;
      D_K1_L n; D_K2_L n; U_K1_R n; U_K2_R n; D_K1_R n; D_K2_R n]

(* We don't introduce a special index for the higher excitations but make
   them parts of the particles' names. *)

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Gauge Bosons", List.map gauge_boson [Ga; Z; Wp; Wm; Gl; 
              Gl_K1; Gl_K2; B1; B2; Z1; Z2; Wp1 ; Wm1; Wp2; Wm2];
        "Higgs", List.map other [H; H1up; H1um; H1dp; H1dm;
              H2up; H2um; H2dp; H2dm];
        "Graviton", List.map other [Grav];
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
          | L_K1_L n -> spinor n | L_K1_R n -> spinor n
          | L_K2_L n -> spinor n | L_K2_R n -> spinor n
          | N_K1 n -> spinor n | N_K2 n -> spinor n
          | U_K1_L n -> spinor n | U_K1_R n -> spinor n
          | U_K2_L n -> spinor n | U_K2_R n -> spinor n
          | D_K1_L n -> spinor n | D_K1_R n -> spinor n
          | D_K2_L n -> spinor n | D_K2_R n -> spinor n
          end
      | G f ->
          begin match f with
          | Ga | Gl -> Vector                
          | Wp | Wm | Z | Gl_K1 | Gl_K2 | B1 | B2
          | Z1 | Z2 | Wp1 | Wm1 | Wp2 | Wm2 -> Massive_Vector
          end
      | O f -> 
          begin match f with 
          | Grav -> Tensor_2
          | _ -> Scalar
          end

    let color = function 
      | M (U n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D n) -> Color.SUN  (if n > 0 then 3 else -3)
      | M (U_K1_L n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D_K1_L n) -> Color.SUN  (if n > 0 then 3 else -3)
      | M (U_K1_R n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D_K1_R n) -> Color.SUN  (if n > 0 then 3 else -3)
      | M (U_K2_L n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D_K2_L n) -> Color.SUN  (if n > 0 then 3 else -3)
      | M (U_K2_R n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D_K2_R n) -> Color.SUN  (if n > 0 then 3 else -3)
      | G Gl | G Gl_K1 | G Gl_K2 -> Color.AdjSUN 3
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
          | L_K1_L n -> prop_spinor n | L_K1_R n -> prop_spinor n
          | L_K2_L n -> prop_spinor n | L_K2_R n -> prop_spinor n
          | N_K1 n -> prop_spinor n | N_K2 n -> prop_spinor n
          | U_K1_L n -> prop_spinor n | U_K1_R n -> prop_spinor n
          | U_K2_L n -> prop_spinor n | U_K2_R n -> prop_spinor n
          | D_K1_L n -> prop_spinor n | D_K1_R n -> prop_spinor n
          | D_K2_L n -> prop_spinor n | D_K2_R n -> prop_spinor n
          end
      | G f ->
          begin match f with
          | Ga | Gl -> Prop_Feynman
          | Wp | Wm | Z | Gl_K1 | Gl_K2 | B1 | B2 
          | Z1 | Z2 | Wp1 | Wm1 | Wp2 | Wm2 -> Prop_Unitarity
          end
      | O f ->
          begin match f with
          | Phip | Phim | Phi0 -> Only_Insertion
          | H | H1up | H1um | H1dp | H1dm | H2up 
          | H2um | H2dp | H2dm -> Prop_Scalar
          | Grav -> Prop_Tensor_2
          end

(* Optionally, ask for the fudge factor treatment for the widths of
   charged particles.  Currently, this only applies to $W^\pm$ and top. *)

    let width f =
      if !use_fudged_width then
        match f with
        | G Wp | G Wm | M (U 3) | M (U (-3)) | O Grav -> Fudged
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
          | L_K1_L n -> L_K1_L (-n) | L_K1_R n -> L_K1_R (-n)
          | L_K2_L n -> L_K2_L (-n) | L_K2_R n -> L_K2_R (-n)
          | N_K1 n -> N_K1 (-n) | N_K2 n -> N_K2 (-n)
          | U_K1_L n -> U_K1_L (-n) | U_K1_R n -> U_K1_R (-n)
          | U_K2_L n -> U_K2_L (-n) | U_K2_R n -> U_K2_R (-n)
          | D_K1_L n -> D_K1_L (-n) | D_K1_R n -> D_K1_R (-n)
          | D_K2_L n -> D_K2_L (-n) | D_K2_R n -> D_K2_R (-n)
          end)
      | G f ->
          G (begin match f with
          | Gl -> Gl | Ga -> Ga | Z -> Z
          | Wp -> Wm | Wm -> Wp 
          | Gl_K1 -> Gl_K1 | Gl_K2 -> Gl_K2 | B1 -> B1 | B2 -> B2
          | Z1 -> Z1 | Z2 -> Z2 | Wp1 -> Wm1 | Wm1 -> Wp1 
          | Wp2 -> Wm2 | Wm2 -> Wp2 
          end)
      | O f ->
          O (begin match f with
          | Phip -> Phim | Phim -> Phip | Phi0 -> Phi0
          | H -> H | H1up -> H1um | H1um -> H1up 
          | H1dp -> H1dm | H1dm -> H1dp 
          | H2up -> H2um | H2um -> H2up 
          | H2dp -> H2dm | H2dm -> H2dp 
          | Grav -> Grav
          end)

    let fermion = function
      | M f ->
          begin match f with
          | L n -> if n > 0 then 1 else -1
          | N n -> if n > 0 then 1 else -1
          | U n -> if n > 0 then 1 else -1
          | D n -> if n > 0 then 1 else -1
          | L_K1_L n -> if n > 0 then 1 else -1
          | L_K2_L n -> if n > 0 then 1 else -1
          | L_K1_R n -> if n > 0 then 1 else -1
          | L_K2_R n -> if n > 0 then 1 else -1
          | U_K1_L n -> if n > 0 then 1 else -1
          | U_K2_L n -> if n > 0 then 1 else -1
          | U_K1_R n -> if n > 0 then 1 else -1
          | U_K2_R n -> if n > 0 then 1 else -1
          | D_K1_L n -> if n > 0 then 1 else -1
          | D_K2_L n -> if n > 0 then 1 else -1
          | D_K1_R n -> if n > 0 then 1 else -1
          | D_K2_R n -> if n > 0 then 1 else -1
          | N_K1 n -> if n > 0 then 1 else -1
          | N_K2 n -> if n > 0 then 1 else -1
          end
      | G f ->
          begin match f with
          | Gl | Ga | Z | Wp | Wm | Gl_K1 | Gl_K2 
          | B1 | B2 | Z1 | Z2 | Wp1 | Wm1 | Wp2 
          | Wm2 -> 0
          end
      | O _ -> 0

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
      match f with
      | M (L n | N n | U n | D n | L_K1_L n | L_K2_L n
         | L_K1_R n | L_K2_R n | N_K1 n | N_K2 n | U_K1_L n
         | U_K2_L n | U_K1_R n | U_K2_R n | D_K1_L n | D_K2_L n
         | D_K1_R n | D_K2_R n ) -> generation' n
      | G _ | O _ -> [0//1; 0//1; 0//1]

    let charge = function
      | M f ->
          begin match f with
          | L n | L_K1_L n | L_K2_L n | L_K1_R n 
          | L_K2_R n -> if n > 0 then -1//1 else  1//1
          | N n | N_K1 n | N_K2 n -> 0//1
          | U n | U_K1_L n | U_K2_L n | U_K1_R n 
          | U_K2_R n -> if n > 0 then  2//3 else -2//3
          | D n | D_K1_L n | D_K2_L n | D_K1_R n 
          | D_K2_R n -> if n > 0 then -1//3 else  1//3
          end
      | G f ->
          begin match f with
          | Gl | Gl_K1 | Gl_K2 | Ga | Z 
          | B1 | B2 | Z1 | Z2 -> 0//1
          | Wp | Wp1 | Wp2 ->  1//1
          | Wm | Wm1 | Wm2 -> -1//1
          end
      | O f ->
          begin match f with
          | H | Phi0 | Grav ->  0//1
          | H1up | H1dp | H2up | H2dp | Phip ->  1//1
          | H1um | H1dm | H2um | H2dm | Phim -> -1//1
          end

    let lepton = function
      | M f ->
          begin match f with
          | L n | N n | L_K1_L n | L_K1_R n | L_K2_L n 
          | L_K2_R n | N_K1 n | N_K2 n -> if n > 0 then 1//1 else -1//1
          | U _ | D _ | _ -> 0//1
          end
      | G _ | O _ -> 0//1

    let baryon = function
      | M f ->
          begin match f with
          | U n | D n | U_K1_L n | U_K1_R n | U_K2_L n 
          | U_K2_R n | D_K1_L n | D_K1_R n | D_K2_L n 
          | D_K2_R n -> if n > 0 then 1//1 else -1//1
          | L _ | N _ | _ -> 0//1
          end
      | G _ | O _ -> 0//1

    let charges f = 
      [ charge f; lepton f; baryon f] @ generation f

    type constant =
      | Unit | Pi | Alpha_QED | Sin2thw
      | Sinthw | Costhw | E | G_weak | Vev
      | Q_lepton | Q_up | Q_down | G_CC | G_CCQ of int*int
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down
      | I_Q_W | I_G_ZWW | I_Q_W_K | I_G_ZWW_K1 | I_G_ZWW_K2 
      | I_G_ZWW_K3
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ
      | G_Htt | G_Hbb | G_Hcc | G_Htautau | G_H3 | G_H4
      | G_HGaZ | G_HGaGa | G_Hgg
      | Gs | I_Gs | I_GsRt2 | G2 | G22 | G_Grav
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let input_parameters =
      []

    let derived_parameters =
      []

    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

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

    let mgm ((m1, g, m2), fbf, c) = ((M m1, G g, M m2), fbf, c)
    let mom ((m1, o, m2), fbf, c) = ((M m1, O o, M m2), fbf, c)

    let electromagnetic_currents n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]
       
    let neutral_currents n =
      List.map mgm
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

    let charged_currents n = 
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC);
          ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let color_currents n =
        List.map mgm
          [ ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs);
            ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs) ]

    let gravity_currents n = 
      List.map mom
        [ ((L (-n), Grav, L n), Graviton_Spinor_Spinor 1, G_Grav);
          ((N (-n), Grav, N n), Graviton_Spinor_Spinor 1, G_Grav);
          ((U (-n), Grav, U n), Graviton_Spinor_Spinor 1, G_Grav);
          ((D (-n), Grav, D n), Graviton_Spinor_Spinor 1, G_Grav) ] 

    let yukawa =
      List.map mom 
        [ ((U (-3), H, U 3), FBF (1, Psibar, S, Psi), G_Htt);
          ((D (-3), H, D 3), FBF (1, Psibar, S, Psi), G_Hbb);
          ((U (-2), H, U 2), FBF (1, Psibar, S, Psi), G_Hcc);
          ((L (-3), H, L 3), FBF (1, Psibar, S, Psi), G_Htautau) ]

    let tgc ((g1, g2, g3), t, c) = ((G g1, G g2, G g3), t, c)

(* Gluons should be included in just that way. *)

    let standard_triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Ga, Wm1, Wp1), Gauge_Gauge_Gauge 1, I_Q_W_K);
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((Z, Wm1, Wp1), Gauge_Gauge_Gauge 1, I_G_ZWW_K1);
          ((Z1, Wm, Wp1), Gauge_Gauge_Gauge 1, I_G_ZWW_K2);
          ((Z1, Wm1, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW_K2);
          ((Z2, Wm1, Wp2), Gauge_Gauge_Gauge 1, I_G_ZWW_K3);
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs);
          ((Gl, Gl_K2, Gl_K2), Gauge_Gauge_Gauge (-1), I_Gs);
          ((Gl, Gl_K1, Gl_K1), Gauge_Gauge_Gauge 1, I_Gs);
          ((Gl_K2, Gl_K1, Gl_K1), Gauge_Gauge_Gauge 1, I_GsRt2)]

    let triple_gauge =
        standard_triple_gauge

    let qgc ((g1, g2, g3, g4), t, c) = ((G g1, G g2, G g3, G g4), t, c)

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let standard_quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          ((Gl, Gl, Gl, Gl), gauge4, G2);
          ((Gl, Gl, Gl_K1, Gl_K1), gauge4, G2);
          ((Gl, Gl, Gl_K2, Gl_K2), gauge4, G2);
          ((Gl_K1, Gl_K1, Gl_K2, Gl_K2), gauge4, G2);
          ((Gl_K2, Gl_K2, Gl_K2, Gl_K2), gauge4, G22)]

    let quartic_gauge =
      standard_quartic_gauge

    let gravity_gauge = 
      [ (O Grav, G Z, G Z), Graviton_Vector_Vector 1, G_Grav;
        (O Grav, G Wp, G Wm), Graviton_Vector_Vector 1, G_Grav;
        (O Grav, G Ga, G Ga), Graviton_Vector_Vector 1, G_Grav;
        (O Grav, G Gl, G Gl), Graviton_Vector_Vector 1, G_Grav ]

    let standard_gauge_higgs =
      [ ((O H, G Wp, G Wm), Scalar_Vector_Vector 1, G_HWW);
        ((O H, G Z, G Z), Scalar_Vector_Vector 1, G_HZZ) ]

    let standard_gauge_higgs4 =
      [ (O H, O H, G Wp, G Wm), Scalar2_Vector2 1, G_HHWW;
        (O H, O H, G Z, G Z), Scalar2_Vector2 1, G_HHZZ ]
       
    let standard_higgs =
      [ (O H, O H, O H), Scalar_Scalar_Scalar 1, G_H3 ]

    let standard_higgs4 =
      [ (O H, O H, O H, O H), Scalar4 1, G_H4 ]

    let gravity_higgs = 
      [ (O Grav, O H, O H), Graviton_Scalar_Scalar 1, G_Grav]

    let anomalous_gauge_higgs =
      []

    let anomalous_gauge_higgs4 =
      []

    let anomalous_higgs =
      []

    let anomaly_higgs = 
      [ (O H, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
        (O H, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_HGaZ;
        (O H, G Gl, G Gl), Dim5_Scalar_Gauge2 1, G_Hgg ]

    let anomalous_higgs4 =
      []

    let gauge_higgs =
        standard_gauge_higgs

    let gauge_higgs4 =
        standard_gauge_higgs4

    let higgs =
        standard_higgs @ gravity_higgs

    let higgs4 =
        standard_higgs4

    let goldstone_vertices =
      [ ((O Phi0, G Wm, G Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phip, G Ga, G Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phip, G Z, G Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phim, G Wp, G Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phim, G Wp, G Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @
       ThoList.flatmap neutral_currents [1;2;3] @
       ThoList.flatmap charged_currents [1;2;3] @
       ThoList.flatmap color_currents [1;2;3] @
       ThoList.flatmap gravity_currents [1;2;3] @
       yukawa @ triple_gauge @ gravity_gauge @
       gauge_higgs @ higgs @ anomaly_higgs 
       @ goldstone_vertices)

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
      | "uk1l" -> M (U_K1_L 1) | "uk1lbar" -> M (U_K1_L (-1))
      | "ck1l" -> M (U_K1_L 2) | "ck1lbar" -> M (U_K1_L (-2))
      | "tk1l" -> M (U_K1_L 3) | "tk1lbar" -> M (U_K1_L (-3))
      | "dk1l" -> M (D_K1_L 1) | "dk1lbar" -> M (D_K1_L (-1))
      | "sk1l" -> M (D_K1_L 2) | "sk1lbar" -> M (D_K1_L (-2))
      | "bk1l" -> M (D_K1_L 3) | "bk1lbar" -> M (D_K1_L (-3))
      | "uk1r" -> M (U_K1_R 1) | "uk1rbar" -> M (U_K1_R (-1))
      | "ck1r" -> M (U_K1_R 2) | "ck1rbar" -> M (U_K1_R (-2))
      | "tk1r" -> M (U_K1_R 3) | "tk1rbar" -> M (U_K1_R (-3))
      | "dk1r" -> M (D_K1_R 1) | "dk1rbar" -> M (D_K1_R (-1))
      | "sk1r" -> M (D_K1_R 2) | "sk1rbar" -> M (D_K1_R (-2))
      | "bk1r" -> M (D_K1_R 3) | "bk1rbar" -> M (D_K1_R (-3))
      | "uk2l" -> M (U_K2_L 1) | "uk2lbar" -> M (U_K2_L (-1))
      | "ck2l" -> M (U_K2_L 2) | "ck2lbar" -> M (U_K2_L (-2))
      | "tk2l" -> M (U_K2_L 3) | "tk2lbar" -> M (U_K2_L (-3))
      | "dk2l" -> M (D_K2_L 1) | "dk2lbar" -> M (D_K2_L (-1))
      | "sk2l" -> M (D_K2_L 2) | "sk2lbar" -> M (D_K2_L (-2))
      | "bk2l" -> M (D_K2_L 3) | "bk2lbar" -> M (D_K2_L (-3))
      | "uk2r" -> M (U_K2_R 1) | "uk2rbar" -> M (U_K2_R (-1))
      | "ck2r" -> M (U_K2_R 2) | "ck2rbar" -> M (U_K2_R (-2))
      | "tk2r" -> M (U_K2_R 3) | "tk2rbar" -> M (U_K2_R (-3))
      | "dk2r" -> M (D_K2_R 1) | "dk2rbar" -> M (D_K2_R (-1))
      | "sk2r" -> M (D_K2_R 2) | "sk2rbar" -> M (D_K2_R (-2))
      | "bk2r" -> M (D_K2_R 3) | "bk2rbar" -> M (D_K2_R (-3))
      | "g" | "gl" -> G Gl
      | "g_k1" | "gl_k1" -> G Gl_K1
      | "g_k2" | "gl_k2" -> G Gl_K2
      | "b1" -> G B1 | "b2" -> G B2 | "z1" -> G Z1 | "z2" -> G Z2
      | "W1+" -> G Wp1 | "W1-" -> G Wm1 
      | "W2+" -> G Wp2 | "W2-" -> G Wm2 
      | "A" -> G Ga | "Z" | "Z0" -> G Z
      | "W+" -> G Wp | "W-" -> G Wm
      | "H" -> O H | "H1u+" -> O H1up | "H1u-" -> O H1um
      | "H1d+" -> O H1dp | "H1d-" -> O H1dm
      | "H2u+" -> O H2up | "H2u-" -> O H2um
      | "H2d+" -> O H2dp | "H2d-" -> O H2dm
      | "GG" -> O Grav
      | "ek1l-" -> M (L_K1_L 1) | "ek1l+" -> M (L_K1_L (-1))
      | "muk1l-" -> M (L_K1_L 2) | "mu1l+" -> M (L_K1_L (-2))
      | "tauk1l-" -> M (L_K1_L 3) | "tauk1l+" -> M (L_K1_L (-3))
      | "ek1r-" -> M (L_K1_R 1) | "ek1r+" -> M (L_K1_R (-1))
      | "muk1r-" -> M (L_K1_R 2) | "mu1r+" -> M (L_K1_R (-2))
      | "tau1r-" -> M (L_K1_R 3) | "tauk1r+" -> M (L_K1_R (-3))
      | "ek2l-" -> M (L_K2_L 1) | "ek2l+" -> M (L_K2_L (-1))
      | "muk2l-" -> M (L_K2_L 2) | "mu2l+" -> M (L_K2_L (-2))
      | "tauk2l-" -> M (L_K2_L 3) | "tauk2l+" -> M (L_K2_L (-3))
      | "ek2r-" -> M (L_K2_R 1) | "ek2r+" -> M (L_K2_R (-1))
      | "muk2r-" -> M (L_K2_R 2) | "mu2r+" -> M (L_K2_R (-2))
      | "tau2r-" -> M (L_K2_R 3) | "tauk2r+" -> M (L_K2_R (-3))
      | "nuek1" -> M (N_K1 1) | "nuek1bar" -> M (N_K1 (-1))
      | "numuk1" -> M (N_K1 2) | "numuk1bar" -> M (N_K1 (-2))
      | "nutauk1" -> M (N_K1 3) | "nutauk1bar" -> M (N_K1 (-3))
      | "nuek2" -> M (N_K2 1) | "nuek2bar" -> M (N_K2 (-1))
      | "numuk2" -> M (N_K2 2) | "numuk2bar" -> M (N_K2 (-2))
      | "nutauk2" -> M (N_K2 3) | "nutauk2bar" -> M (N_K2 (-3))
      | _ -> invalid_arg "Modellib_BSM.UED.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid down type quark"
          | U_K1_L 1 -> "uk1l" | U_K1_L (-1) -> "uk1lbar"
          | U_K1_L 2 -> "ck1l" | U_K1_L (-2) -> "ck1lbar"
          | U_K1_L 3 -> "tk1l" | U_K1_L (-3) -> "tk1lbar"
          | U_K1_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid up type quark"
          | D_K1_L 1 -> "dk1l" | D_K1_L (-1) -> "dk1lbar"
          | D_K1_L 2 -> "sk1l" | D_K1_L (-2) -> "sk1lbar"
          | D_K1_L 3 -> "bk1l" | D_K1_L (-3) -> "bk1lbar"
          | D_K1_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid down type quark"
          | U_K1_R 1 -> "uk1r" | U_K1_R (-1) -> "uk1rbar"
          | U_K1_R 2 -> "ck1r" | U_K1_R (-2) -> "ck1rbar"
          | U_K1_R 3 -> "tk1r" | U_K1_R (-3) -> "tk1rbar"
          | U_K1_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid up type quark"
          | D_K1_R 1 -> "dk1r" | D_K1_R (-1) -> "dk1rbar"
          | D_K1_R 2 -> "sk1r" | D_K1_R (-2) -> "sk1rbar"
          | D_K1_R 3 -> "bk1r" | D_K1_R (-3) -> "bk1rbar"
          | D_K1_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid down type quark"
          | U_K2_L 1 -> "uk2l" | U_K2_L (-1) -> "uk2lbar"
          | U_K2_L 2 -> "ck2l" | U_K2_L (-2) -> "ck2lbar"
          | U_K2_L 3 -> "tk2l" | U_K2_L (-3) -> "tk2lbar"
          | U_K2_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid up type quark"
          | D_K2_L 1 -> "dk2l" | D_K2_L (-1) -> "dk2lbar"
          | D_K2_L 2 -> "sk2l" | D_K2_L (-2) -> "sk2lbar"
          | D_K2_L 3 -> "bk2l" | D_K2_L (-3) -> "bk2lbar"
          | D_K2_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid down type quark"
          | U_K2_R 1 -> "uk2r" | U_K2_R (-1) -> "uk2rbar"
          | U_K2_R 2 -> "ck2r" | U_K2_R (-2) -> "ck2rbar"
          | U_K2_R 3 -> "tk2r" | U_K2_R (-3) -> "tk2rbar"
          | U_K2_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid up type quark"
          | D_K2_R 1 -> "dk2r" | D_K2_R (-1) -> "dk2rbar"
          | D_K2_R 2 -> "sk2r" | D_K2_R (-2) -> "sk2rbar"
          | D_K2_R 3 -> "bk2r" | D_K2_R (-3) -> "bk2rbar"
          | D_K2_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid down type quark"
          | L_K1_L 1 -> "ek1l-" | L_K1_L (-1) -> "ek1l+"
          | L_K1_L 2 -> "muk1l-" | L_K1_L (-2) -> "muk1l+"
          | L_K1_L 3 -> "tauk1l-" | L_K1_L (-3) -> "tauk1l+"
          | L_K1_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid lepton"
          | L_K1_R 1 -> "ek1r-" | L_K1_R (-1) -> "ek1r+"
          | L_K1_R 2 -> "muk1r-" | L_K1_R (-2) -> "muk1r+"
          | L_K1_R 3 -> "tauk1r-" | L_K1_R (-3) -> "tauk1r+"
          | L_K1_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid lepton"
          | L_K2_L 1 -> "ek2l-" | L_K2_L (-1) -> "ek2l+"
          | L_K2_L 2 -> "muk2l-" | L_K2_L (-2) -> "muk2l+"
          | L_K2_L 3 -> "tauk2l-" | L_K2_L (-3) -> "tauk2l+"
          | L_K2_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid lepton"
          | L_K2_R 1 -> "ek2r-" | L_K2_R (-1) -> "ek2r+"
          | L_K2_R 2 -> "muk2r-" | L_K2_R (-2) -> "muk2r+"
          | L_K2_R 3 -> "tauk2r-" | L_K2_R (-3) -> "tauk2r+"
          | L_K2_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid lepton"
          | N_K1 1 -> "nuek1" | N_K1 (-1) -> "nuek1bar"
          | N_K1 2 -> "numuk1" | N_K1 (-2) -> "numuk1bar"
          | N_K1 3 -> "nutauk1" | N_K1 (-3) -> "nutauk1bar"
          | N_K1 _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid neutrino"
          | N_K2 1 -> "nuek2" | N_K2 (-1) -> "nuek2bar"
          | N_K2 2 -> "numuk2" | N_K2 (-2) -> "numuk2bar"
          | N_K2 3 -> "nutauk2" | N_K2 (-3) -> "nutauk2bar"
          | N_K2 _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_string: invalid neutrino"
          end
      | G f ->
          begin match f with
          | Gl -> "g" 
          | Ga -> "A" | Z -> "Z"
          | Wp -> "W+" | Wm -> "W-"
          | Gl_K1 -> "gk1" | Gl_K2 -> "gk2"
          | B1 -> "b1" | B2 -> "b2" | Z1 -> "z1" | Z2 -> "z2"
          | Wp1 -> "W1+" | Wm1 -> "W1-" 
          | Wp2 -> "W2+" | Wm2 -> "W2-" 
          end
      | O f ->
          begin match f with
          | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
          | H -> "H" | H1up -> "H1u+" | H1um -> "H1u-"
          | H1dp -> "H1d+" | H1dm -> "H1d-"
          | H2up -> "H2u+" | H2um -> "H2u-"
          | H2dp -> "H2d+" | H2dm -> "H2d-"
          | Grav -> "GG"
          end

    let flavor_to_TeX = function
      | M f ->
          begin match f with
          | L 1 -> "e^-" | L (-1) -> "e^+"
          | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
          | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
          | L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid down type quark"
          | U_K1_L 1 -> "u^\\prime_L" | U_K1_L (-1) -> "\\bar{u}^\\prime_L"
          | U_K1_L 2 -> "c^\\prime_L" | U_K1_L (-2) -> "\\bar{c}^\\prime_L"
          | U_K1_L 3 -> "t^\\prime_L" | U_K1_L (-3) -> "\\bar{t}^\\prime_L"
          | U_K1_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid up type quark"
          | D_K1_L 1 -> "d^\\prime_L" | D_K1_L (-1) -> "\\bar{d}^\\prime_L"
          | D_K1_L 2 -> "s^\\prime_L" | D_K1_L (-2) -> "\\bar{s}^\\prime_L"
          | D_K1_L 3 -> "b^\\prime_L" | D_K1_L (-3) -> "\\bar{b}^\\prime_L"
          | D_K1_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid down type quark"
          | U_K1_R 1 -> "u^\\prime_R" | U_K1_R (-1) -> "\\bar{u}^\\prime_R"
          | U_K1_R 2 -> "c^\\prime_R" | U_K1_R (-2) -> "\\bar{c}^\\prime_R"
          | U_K1_R 3 -> "t^\\prime_R" | U_K1_R (-3) -> "\\bar{t}^\\prime_R"
          | U_K1_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid up type quark"
          | D_K1_R 1 -> "d^\\prime_R" | D_K1_R (-1) -> "\\bar{d}^\\prime_R"
          | D_K1_R 2 -> "s^\\prime_R" | D_K1_R (-2) -> "\\bar{s}^\\prime_R"
          | D_K1_R 3 -> "b^\\prime_R" | D_K1_R (-3) -> "\\bar{b}^\\prime_R"
          | D_K1_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid down type quark"
          | U_K2_L 1 -> "u^{\\prime\\prime}_L" | U_K2_L (-1) -> "\\bar{u}^{\\prime\\prime}_L"
          | U_K2_L 2 -> "c^{\\prime\\prime}_L" | U_K2_L (-2) -> "\\bar{c}^{\\prime\\prime}_L"
          | U_K2_L 3 -> "t^{\\prime\\prime}_L" | U_K2_L (-3) -> "\\bar{t}^{\\prime\\prime}_L"
          | U_K2_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid up type quark"
          | D_K2_L 1 -> "d^{\\prime\\prime}_L" | D_K2_L (-1) -> "\\bar{d}^{\\prime\\prime}_L"
          | D_K2_L 2 -> "s^{\\prime\\prime}_L" | D_K2_L (-2) -> "\\bar{s}^{\\prime\\prime}_L"
          | D_K2_L 3 -> "b^{\\prime\\prime}_L" | D_K2_L (-3) -> "\\bar{b}^{\\prime\\prime}_L"
          | D_K2_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid down type quark"
          | U_K2_R 1 -> "u^{\\prime\\prime}_R" | U_K2_R (-1) -> "\\bar{u}^{\\prime\\prime}_R"
          | U_K2_R 2 -> "c^{\\prime\\prime}_R" | U_K2_R (-2) -> "\\bar{c}^{\\prime\\prime}_R"
          | U_K2_R 3 -> "t^{\\prime\\prime}_R" | U_K2_R (-3) -> "\\bar{t}^{\\prime\\prime}_R"
          | U_K2_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid up type quark"
          | D_K2_R 1 -> "d^\\prime_R" | D_K2_R (-1) -> "\\bar{d}^{\\prime\\prime}_R"
          | D_K2_R 2 -> "s^\\prime_R" | D_K2_R (-2) -> "\\bar{s}^{\\prime\\prime}_R"
          | D_K2_R 3 -> "b^\\prime_R" | D_K2_R (-3) -> "\\bar{b}^{\\prime\\prime}_R"
          | D_K2_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid down type quark"
          | L_K1_L 1 -> "e_L^{\\prime,,-}" | L_K1_L (-1) -> "\\bar{e}_L^{\\prime,,+}"
          | L_K1_L 2 -> "\\mu_L^{\\prime,,-}" | L_K1_L (-2) -> "\\bar{\\mu}_L^{{\\prime,,+}"
          | L_K1_L 3 -> "\\tau_L^{\\prime,,-}" | L_K1_L (-3) -> "\\bar{\\tau}_L^{\\prime,,+}"
          | L_K1_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid lepton"
          | L_K1_R 1 -> "e_R^{\\prime,,-}" | L_K1_R (-1) -> "\\bar{e}_R^{\\prime,,+}"
          | L_K1_R 2 -> "\\mu_R{\\prime,,-}" | L_K1_R (-2) -> "\\bar{\\mu}_R^{\\prime,,+}"
          | L_K1_R 3 -> "\\tau_R¬{\\prime,,-}" | L_K1_R (-3) -> "\\bar{\\tau}_R¬{\\prime,,+}"
          | L_K1_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid lepton"
          | L_K2_L 1 -> "e^{\\prime\\prime,,-}_L" | L_K2_L (-1) -> "\\bar{e}_L^{\\prime\\prime,,+}"
          | L_K2_L 2 -> "\\mu_L^{\\prime\\prime,,-}" | L_K2_L (-2) -> "\\bar{\\mu}_L^{\\prime\\prime,,+}"
          | L_K2_L 3 -> "\\tau_L^{\\prime\\prime,,-}" | L_K2_L (-3) -> "\\bar{\\tau}_L^{\\prime\\prime,,+}"
          | L_K2_L _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid lepton"
          | L_K2_R 1 -> "e_R^{\\prime\\prime,,-}" | L_K2_R (-1) -> "\\bar{e}_R^{\\prime\\prime,,+}"
          | L_K2_R 2 -> "\\mu_R^{\\prime\\prime,,-}" | L_K2_R (-2) -> "\\bar{\\mu}_R^{\\prime\\prime,,+}"
          | L_K2_R 3 -> "\\tau_R{\\prime\\prime,,-}" | L_K2_R (-3) -> "\\bar{\\tau}_R^{\\prime\\prime,,+}"
          | L_K2_R _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid lepton"
          | N_K1 1 -> "\\nu_e^\\prime" | N_K1 (-1) -> "\\bar{\\nu}_e^\\prime"
          | N_K1 2 -> "\\nu_\\mu^\\prime" | N_K1 (-2) -> "\\bar{\\nu}_\\mu^\\prime"
          | N_K1 3 -> "\\nu_\\tau^\\prime" | N_K1 (-3) -> "\\bar{\\nu}_\\tau^\\prime"
          | N_K1 _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid neutrino"
          | N_K2 1 -> "\\nu_e^{\\prime\\prime}" | N_K2 (-1) -> "\\bar{\\nu}_e^{\\prime\\prime}"
          | N_K2 2 -> "\\nu_\\mu^{\\prime\\prime}" | N_K2 (-2) -> "\\bar{\\nu}_\\mu^{\\prime\\prime}"
          | N_K2 3 -> "\\nu_\\tau^{\\prime\\prime}" | N_K2 (-3) -> "\\bar{\\nu}_\\tau^{\\prime\\prime}"
          | N_K2 _ -> invalid_arg
                "Modellib_BSM.UED.flavor_to_TeX: invalid neutrino"
          end
      | G f ->
          begin match f with
          | Gl -> "g" 
          | Ga -> "\\gamma" | Z -> "Z"
          | Wp -> "W^+" | Wm -> "W^-"
          | Gl_K1 -> "g^\\prime" | Gl_K2 -> "g^{\\prime\\prime}"
          | B1 -> "B^\\prime" | B2 -> "B^{\\prime\\prime}" 
          | Z1 -> "Z^\\prime" | Z2 -> "Z^{\\prime\\prime}"
          | Wp1 -> "W^{\\prime,,+}" | Wm1 -> "W^{\\prime,,-}" 
          | Wp2 -> "W^{\\prime\\prime,,+}" | Wm2 -> "W^{\\prime\\prime,,-}" 
          end
      | O f ->
          begin match f with
          | Phip -> "\\phi^+" | Phim -> "\\phi^-" | Phi0 -> "\\phi^0" 
          | H -> "H" | H1up -> "H1u+" | H1um -> "H1u-"
          | H1dp -> "H1d+" | H1dm -> "H1d-"
          | H2up -> "H2u+" | H2um -> "H2u-"
          | H2dp -> "H2d+" | H2dm -> "H2d-"
          | Grav -> "G^\\prime"
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
          | L_K1_L n when n > 0 -> "lk1l" ^ string_of_int n
          | L_K1_L n -> "lk1l" ^ string_of_int (abs n) ^ "b"
          | L_K1_R n when n > 0 -> "lk1r" ^ string_of_int n
          | L_K1_R n -> "lk1r" ^ string_of_int (abs n) ^ "b"
          | L_K2_L n when n > 0 -> "lk2l" ^ string_of_int n
          | L_K2_L n -> "lk2l" ^ string_of_int (abs n) ^ "b"
          | L_K2_R n when n > 0 -> "lk2r" ^ string_of_int n
          | L_K2_R n -> "lk2r" ^ string_of_int (abs n) ^ "b"
          | U_K1_L n when n > 0 -> "uk1l" ^ string_of_int n
          | U_K1_L n -> "uk1l" ^ string_of_int (abs n) ^ "b"
          | U_K1_R n when n > 0 -> "uk1r" ^ string_of_int n
          | U_K1_R n -> "uk1r" ^ string_of_int (abs n) ^ "b"
          | U_K2_L n when n > 0 -> "uk2l" ^ string_of_int n
          | U_K2_L n -> "uk2l" ^ string_of_int (abs n) ^ "b"
          | U_K2_R n when n > 0 -> "uk2r" ^ string_of_int n
          | U_K2_R n -> "uk2r" ^ string_of_int (abs n) ^ "b"
          | D_K1_L n when n > 0 -> "dk1l" ^ string_of_int n
          | D_K1_L n -> "dk1l" ^ string_of_int (abs n) ^ "b"
          | D_K1_R n when n > 0 -> "dk1r" ^ string_of_int n
          | D_K1_R n -> "dk1r" ^ string_of_int (abs n) ^ "b"
          | D_K2_L n when n > 0 -> "dk2l" ^ string_of_int n
          | D_K2_L n -> "dk2l" ^ string_of_int (abs n) ^ "b"
          | D_K2_R n when n > 0 -> "dk2r" ^ string_of_int n
          | D_K2_R n -> "dk2r" ^ string_of_int (abs n) ^ "b"
          | N_K1 n when n > 0 -> "nk1" ^ string_of_int n
          | N_K1 n -> "nk1" ^ string_of_int (abs n) ^ "b"
          | N_K2 n when n > 0 -> "nk2" ^ string_of_int n
          | N_K2 n -> "nk2" ^ string_of_int (abs n) ^ "b"
          end
      | G f ->
          begin match f with
          | Gl -> "gl"
          | Ga -> "a" | Z -> "z"
          | Wp -> "wp" | Wm -> "wm"
          | Gl_K1 -> "gk1" | Gl_K2 -> "gk2"
          | B1 -> "b1" | B2 -> "b2" | Z1 -> "z1" | Z2 -> "z2"
          | Wp1 -> "wp1" | Wm1 -> "wm1" 
          | Wp2 -> "wp2" | Wm2 -> "wm2" 
          end
      | O f ->
          begin match f with
          | Phip -> "pp" | Phim -> "pm" | Phi0 -> "p0" 
          | H -> "h" | H1up -> "h1up" | H1um -> "h1um"
          | H1dp -> "h1dp" | H1dm -> "h1dm"
          | H2up -> "h2up" | H2um -> "h2um"
          | H2dp -> "h2dp" | H2dm -> "h2dm"
          | Grav -> "gv"
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
          | U_K1_L n when n > 0 -> 4000000 + 2*n
          | U_K1_L n -> - 4000000 + 2*n
          | D_K1_L n when n > 0 -> 3999999 + 2*n
          | D_K1_L n -> - 3999999 + 2*n
          | U_K1_R n when n > 0 -> 5000000 + 2*n
          | U_K1_R n -> - 5000000 + 2*n
          | D_K1_R n when n > 0 -> 4999999 + 2*n
          | D_K1_R n -> - 4999999 + 2*n
          | U_K2_L n when n > 0 -> 6000000 + 2*n
          | U_K2_L n -> - 6000000 + 2*n
          | D_K2_L n when n > 0 -> 5999999 + 2*n
          | D_K2_L n -> - 5999999 + 2*n
          | U_K2_R n when n > 7000000 -> 2*n
          | U_K2_R n -> - 7000000 + 2*n
          | D_K2_R n when n > 0 -> 6999999 + 2*n
          | D_K2_R n -> - 6999999 + 2*n
          | L_K1_L n when n > 0 -> 4000009 + 2*n
          | L_K1_L n -> - 4000009 + 2*n
          | L_K1_R n when n > 0 -> 5000009 + 2*n
          | L_K1_R n -> - 5000009 + 2*n
          | L_K2_L n when n > 0 -> 6000009 + 2*n
          | L_K2_L n -> - 6000009 + 2*n
          | L_K2_R n when n > 0 -> 7000009 + 2*n
          | L_K2_R n -> - 7000009 + 2*n
          | N_K1 n when n > 0 -> 4000010 + 2*n
          | N_K1 n -> - 4000010 + 2*n
          | N_K2 n when n > 0 -> 6000010 + 2*n
          | N_K2 n -> - 6000010 + 2*n
          end
      | G f ->
          begin match f with
          | Gl -> 21
          | Ga -> 22 | Z -> 23
          | Wp -> 24 | Wm -> (-24)
          | Gl_K1 -> 4000021 | Gl_K2 -> 6000021
          | B1 -> 4000022 | B2 -> 6000022
          | Z1 -> 4000023 | Z2 -> 6000024
          | Wp1 -> 4000024 | Wm1 -> (-4000024)
          | Wp2 -> 6000024 | Wm2 -> (-6000024)
          end
      | O f ->
          begin match f with
          | Phip | Phim -> 27 | Phi0 -> 26
          | H -> 25 | H1up -> 4000036 | H1um -> (-4000036)
          | H1dp -> 4000037 | H1dm -> (-4000037)
          | H2up -> 6000036 | H2um -> (-6000036)
          | H2dp -> 6000037 | H2dm -> (-6000037)
          | Grav -> 39
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
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_CC -> "gcc"
      | G_CCQ (n1,n2) -> "gccq" ^ string_of_int n1 ^ string_of_int n2
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" 
      | I_Q_W_K -> "iqwk" | I_G_ZWW_K1 -> "igzwwk1" 
      | I_G_ZWW_K2 -> "igzwwk2" | I_G_ZWW_K3 -> "igzwwk3"  
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc"
      | G_HGaZ -> "ghgaz" | G_HGaGa -> "ghgaga" | G_Hgg -> "ghgg"
      | G_H3 -> "gh3" | G_H4 -> "gh4"
      | G2 -> "gs**2" | Gs -> "gs" | I_Gs -> "igs" | I_GsRt2 -> "igs/sqrt(2.0_default)"
      | G22 -> "gs**2/2.0_default"
      | G_Grav -> "ggrav"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f

  end

module GravTest (Flags : BSM_flags) =
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

    type matter_field = L of int | N of int | U of int | D of int | SL of int
    type gauge_boson = Ga | Wp | Wm | Z | Gl | Phino
    type other = Phip | Phim | Phi0 | H | Grino
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
      failwith "Modellib_BSM.SM.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; SL n; N n; U n; D n ]

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Gauge Bosons", List.map gauge_boson [Ga; Z; Wp; Wm; Gl; Phino];
        "Higgs", List.map other [H];
        "Gravitino", List.map other [Grino];
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
          | SL _ -> Scalar 
          end
      | G f ->
          begin match f with
          | Ga | Gl -> Vector
          | Wp | Wm | Z -> Massive_Vector
          | Phino -> Majorana
          end
      | O f -> 
          begin match f with 
          | Grino -> Vectorspinor
          | _ -> Scalar
          end

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
          | SL n -> Prop_Scalar
          end
      | G f ->
          begin match f with
          | Ga | Gl -> Prop_Feynman
          | Wp | Wm | Z -> Prop_Unitarity 
          | Phino -> Prop_Majorana
          end
      | O f ->
          begin match f with
          | Phip | Phim | Phi0 -> Only_Insertion
          | H -> Prop_Scalar
          | Grino -> Prop_Vectorspinor
          end

(* Optionally, ask for the fudge factor treatment for the widths of
   charged particles.  Currently, this only applies to $W^\pm$ and top. *)

    let width f =
      if !use_fudged_width then
        match f with
        | G Wp | G Wm | M (U 3) | M (U (-3)) | O Grino -> Fudged
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
          | SL n -> SL (-n)
          end)
      | G f ->
          G (begin match f with
          | Gl -> Gl | Ga -> Ga | Z -> Z
          | Wp -> Wm | Wm -> Wp | Phino -> Phino 
          end)
      | O f ->
          O (begin match f with
          | Phip -> Phim | Phim -> Phip | Phi0 -> Phi0
          | H -> H | Grino -> Grino
          end)

    let fermion = function
      | M f ->
          begin match f with
          | L n -> if n > 0 then 1 else -1
          | N n -> if n > 0 then 1 else -1
          | U n -> if n > 0 then 1 else -1
          | D n -> if n > 0 then 1 else -1
          | SL _ -> 0
          end
      | G f ->
          begin match f with
          | Gl | Ga | Z | Wp | Wm -> 0
          | Phino -> 2
          end
      | O f -> 
          begin match f with
          | Grino -> 2 
          | _ -> 0
          end

    module Ch = Charges.QQ

    let ( // ) = Algebra.Small_Rational.make

    let generation' = function
      |  1 -> [ 1//1;  0//1;  0//1]
      |  2 -> [ 0//1;  1//1;  0//1]
      |  3 -> [ 0//1;  0//1;  1//1]
      | -1 -> [-1//1;  0//1;  0//1]
      | -2 -> [ 0//1; -1//1;  0//1]
      | -3 -> [ 0//1;  0//1; -1//1]
      |  n -> invalid_arg ("SM3.generation': " ^ string_of_int n)

    let generation f =
      match f with
      | M (L n | N n | U n | D n | SL n) -> generation' n
      | G _ | O _ -> [0//1; 0//1; 0//1]

    let charge = function
      | M f ->
          begin match f with
          | L n -> if n > 0 then -1//1 else  1//1
          | SL n -> if n > 0 then -1//1 else  1//1
          | N n -> 0//1
          | U n -> if n > 0 then  2//3 else -2//3
          | D n -> if n > 0 then -1//3 else  1//3
          end
      | G f ->
          begin match f with
          | Gl | Ga | Z | Phino -> 0//1
          | Wp ->  1//1
          | Wm -> -1//1
          end
      | O f ->
          begin match f with
          | H | Phi0 | Grino ->  0//1
          | Phip ->  1//1
          | Phim -> -1//1
          end

    let lepton = function
      | M f ->
          begin match f with
          | L n | N n | SL n -> if n > 0 then 1//1 else -1//1
          | U _ | D _ -> 0//1
          end
      | G _ | O _ -> 0//1

    let baryon = function
      | M f ->
          begin match f with
          | L _ | N _ | SL _ -> 0//1
          | U n | D n -> if n > 0 then 1//1 else -1//1
          end
      | G _ | O _ -> 0//1

    let charges f =
      [ charge f; lepton f; baryon f] @ generation f

    type constant =
      | Unit | Pi | Alpha_QED | Sin2thw
      | Sinthw | Costhw | E | G_weak | Vev
      | Q_lepton | Q_up | Q_down | G_CC | G_CCQ of int*int
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down
      | I_Q_W | I_G_ZWW
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ
      | G_Htt | G_Hbb | G_Hcc | G_Htautau | G_H3 | G_H4
      | G_HGaZ | G_HGaGa | G_Hgg
      | G_strong | G_Grav
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let input_parameters =
      []

    let derived_parameters =
      []

    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

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

    let mgm ((m1, g, m2), fbf, c) = ((M m1, G g, M m2), fbf, c)
    let mom ((m1, o, m2), fbf, c) = ((M m1, O o, M m2), fbf, c)

    let electromagnetic_currents n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]
       
    let neutral_currents n =
      List.map mgm
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

    let charged_currents n = 
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC);
          ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let yukawa =
      List.map mom 
        [ ((U (-3), H, U 3), FBF (1, Psibar, S, Psi), G_Htt);
          ((D (-3), H, D 3), FBF (1, Psibar, S, Psi), G_Hbb);
          ((U (-2), H, U 2), FBF (1, Psibar, S, Psi), G_Hcc);
          ((L (-3), H, L 3), FBF (1, Psibar, S, Psi), G_Htautau) ]

    let tgc ((g1, g2, g3), t, c) = ((G g1, G g2, G g3), t, c)

    let standard_triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW) ]

    let triple_gauge =
        standard_triple_gauge

    let qgc ((g1, g2, g3, g4), t, c) = ((G g1, G g2, G g3, G g4), t, c)

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let standard_quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW ]

    let quartic_gauge =
      standard_quartic_gauge

    let standard_gauge_higgs =
      [ ((O H, G Wp, G Wm), Scalar_Vector_Vector 1, G_HWW);
        ((O H, G Z, G Z), Scalar_Vector_Vector 1, G_HZZ) ]

    let standard_gauge_higgs4 =
      [ (O H, O H, G Wp, G Wm), Scalar2_Vector2 1, G_HHWW;
        (O H, O H, G Z, G Z), Scalar2_Vector2 1, G_HHZZ ]
       
    let standard_higgs =
      [ (O H, O H, O H), Scalar_Scalar_Scalar 1, G_H3 ]

    let standard_higgs4 =
      [ (O H, O H, O H, O H), Scalar4 1, G_H4 ]

    let anomalous_gauge_higgs =
      []

    let anomalous_gauge_higgs4 =
      []

    let anomalous_higgs =
      []

    let anomaly_higgs = 
      [ (O H, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
        (O H, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_HGaZ;
        (O H, G Gl, G Gl), Dim5_Scalar_Gauge2 1, G_Hgg ]

    let gravitino_coup n = 
      [ (O Grino, M (SL (-n)), M (L n)), GBG (1, Gravbar, POT, Psi), G_Grav;
        (M (L (-n)), M (SL n), O Grino), GBG (1, Psibar, POT, Grav), G_Grav]

    let gravitino_gauge = 
      [ (O Grino, G Ga, G Phino), GBG (1, Gravbar, V, Chi), G_Grav ]


    let anomalous_higgs4 =
      []

    let gauge_higgs =
        standard_gauge_higgs

    let gauge_higgs4 =
        standard_gauge_higgs4

    let higgs =
        standard_higgs 

    let higgs4 =
        standard_higgs4

    let goldstone_vertices =
      [ ((O Phi0, G Wm, G Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phip, G Ga, G Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phip, G Z, G Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phim, G Wp, G Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phim, G Wp, G Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @
       ThoList.flatmap neutral_currents [1;2;3] @
       ThoList.flatmap charged_currents [1;2;3] @
       ThoList.flatmap gravitino_coup [1;2;3] @
       gravitino_gauge @
       yukawa @ triple_gauge @
       gauge_higgs @ higgs @ anomaly_higgs 
       @ goldstone_vertices)

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
      | "se-" -> M (SL 1) | "se+" -> M (SL (-1))
      | "smu-" -> M (SL 2) | "smu+" -> M (SL (-2))
      | "stau-" -> M (SL 3) | "stau+" -> M (SL (-3))
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
      | "GG" -> O Grino
      | "phino" | "Phino" -> G Phino
      | _ -> invalid_arg "Modellib_BSM.GravTest.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib_BSM.GravTest.flavor_to_string: invalid lepton"
          | SL 1 -> "se-" | SL (-1) -> "se+"
          | SL 2 -> "smu-" | SL (-2) -> "smu+"
          | SL 3 -> "stau-" | SL (-3) -> "stau+"
          | SL _ -> invalid_arg
                "Modellib_BSM.GravTest.flavor_to_string: invalid slepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib_BSM.GravTest.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib_BSM.SM.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_BSM.GravTest.flavor_to_string: invalid down type quark"
          end
      | G f ->
          begin match f with
          | Gl -> "g"
          | Ga -> "A" | Z -> "Z"
          | Wp -> "W+" | Wm -> "W-"
          | Phino -> "phino"
          end
      | O f ->
          begin match f with
          | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
          | H -> "H" | Grino -> "GG"
          end

    let flavor_to_TeX = function
      | M f ->
          begin match f with
          | L 1 -> "e^-" | L (-1) -> "e^+"
          | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
          | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
          | L _ -> invalid_arg
                "Modellib_BSM.GravTest.flavor_to_TeX: invalid lepton"
          | SL 1 -> "\\tilde{e}^-" | SL (-1) -> "\\tilde{e}^+"
          | SL 2 -> "\\tilde{\\mu}^-" | SL (-2) -> "\\tilde{\\mu}^+"
          | SL 3 -> "\\tilde{\\tau}^-" | SL (-3) -> "\\tilde{\\tau}^+"
          | SL _ -> invalid_arg
                "Modellib_BSM.GravTest.flavor_to_TeX: invalid slepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib_BSM.GravTest.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib_BSM.SM.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Modellib_BSM.GravTest.flavor_to_TeX: invalid down type quark"
          end
      | G f ->
          begin match f with
          | Gl -> "g"
          | Ga -> "\\gamma" | Z -> "Z"
          | Wp -> "W^+" | Wm -> "W^-"
          | Phino -> "\\tilde{\\phi}"
          end
      | O f ->
          begin match f with
          | Phip -> "\\phi^+" | Phim -> "\\phi^-" | Phi0 -> "\\phi^0" 
          | H -> "H" | Grino -> "\\tilde{G}"
          end

    let flavor_symbol = function
      | M f ->
          begin match f with
          | L n when n > 0 -> "l" ^ string_of_int n
          | L n -> "l" ^ string_of_int (abs n) ^ "b"
          | SL n when n > 0 -> "sl" ^ string_of_int n
          | SL n -> "sl" ^ string_of_int (abs n) ^ "b"
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
          | Phino -> "phino"
          end
      | O f ->
          begin match f with
          | Phip -> "pp" | Phim -> "pm" | Phi0 -> "p0" 
          | H -> "h" | Grino -> "gv" 
          end

    let pdg = function
      | M f ->
          begin match f with
          | L n when n > 0 -> 9 + 2*n
          | L n -> - 9 + 2*n
          | SL n when n > 0 -> 39 + 2*n
          | SL n -> - 39 + 2*n
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
          | Phino -> 46
          end
      | O f ->
          begin match f with
          | Phip | Phim -> 27 | Phi0 -> 26
          | H -> 25 | Grino -> 39
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
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_CC -> "gcc"
      | G_CCQ (n1,n2) -> "gccq" ^ string_of_int n1 ^ string_of_int n2
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc"
      | G_HGaZ -> "ghgaz" | G_HGaGa -> "ghgaga" | G_Hgg -> "ghgg"
      | G_H3 -> "gh3" | G_H4 -> "gh4"
      | G_strong -> "gs" | G_Grav -> "ggrav"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f

  end

module Template (Flags : BSM_flags) =
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
      failwith "Modellib_BSM.Template.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; N n; U n; D n ]

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Gauge Bosons", List.map gauge_boson [Ga; Z; Wp; Wm; Gl];
        "Higgs", List.map other [H];
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
      | O f -> Scalar

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
          | H -> Prop_Scalar
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

    module Ch = Charges.QQ

    let ( // ) = Algebra.Small_Rational.make

    let generation' = function
      |  1 -> [ 1//1;  0//1;  0//1]
      |  2 -> [ 0//1;  1//1;  0//1]
      |  3 -> [ 0//1;  0//1;  1//1]
      | -1 -> [-1//1;  0//1;  0//1]
      | -2 -> [ 0//1; -1//1;  0//1]
      | -3 -> [ 0//1;  0//1; -1//1]
      |  n -> invalid_arg ("Template.generation': " ^ string_of_int n)

    let generation f =
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
      | Q_lepton | Q_up | Q_down | G_CC
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down
      | I_Q_W | I_G_ZWW
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ 
      | G_Htt | G_Hbb | G_Hcc | G_Hmm | G_Htautau | G_H3 | G_H4
      | G_HGaZ | G_HGaGa | G_Hgg
      | Gs | I_Gs | G2
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let input_parameters = []

    let derived_parameters = [] 

    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

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

    let mgm ((m1, g, m2), fbf, c) = ((M m1, G g, M m2), fbf, c)

    let electromagnetic_currents n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]
        
    let color_currents n =
        List.map mgm
          [ ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs);
            ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs) ]

    let neutral_currents n =
      List.map mgm
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

    let charged_currents n = 
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC);  
          ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let yukawa =
      [ ((M (U (-3)), O H, M (U 3)), FBF (1, Psibar, S, Psi), G_Htt);
        ((M (D (-3)), O H, M (D 3)), FBF (1, Psibar, S, Psi), G_Hbb);
        ((M (U (-2)), O H, M (U 2)), FBF (1, Psibar, S, Psi), G_Hcc);
        ((M (L (-2)), O H, M (L 2)), FBF (1, Psibar, S, Psi), G_Hmm);
        ((M (L (-3)), O H, M (L 3)), FBF (1, Psibar, S, Psi), G_Htautau) ]

    let tgc ((g1, g2, g3), t, c) = ((G g1, G g2, G g3), t, c)

    let triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs) ]

    let qgc ((g1, g2, g3, g4), t, c) = ((G g1, G g2, G g3, G g4), t, c)

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          (Gl, Gl, Gl, Gl), gauge4, G2]

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

    let anomaly_higgs = 
      []
(*      [ (O H, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
        (O H, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_HGaZ; 
        (O H, G Gl, G Gl), Dim5_Scalar_Gauge2 1, G_Hgg]  *)

    let goldstone_vertices =
      [ ((O Phi0, G Wm, G Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phip, G Ga, G Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phip, G Z, G Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phim, G Wp, G Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phim, G Wp, G Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @
       ThoList.flatmap color_currents [1;2;3] @
       ThoList.flatmap neutral_currents [1;2;3] @
       ThoList.flatmap charged_currents [1;2;3] @
       yukawa @ triple_gauge @ gauge_higgs @ higgs @ 
       anomaly_higgs @ goldstone_vertices)

    let vertices4 =
      quartic_gauge @ gauge_higgs4 @ higgs4

    let vertices () = (vertices3, vertices4, [])

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
      | _ -> invalid_arg "Modellib_BSM.Template.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib_BSM.Template.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib_BSM.Template.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib_BSM.Template.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_BSM.Template.flavor_to_string: invalid down type quark"
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
                "Modellib_BSM.Template.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib_BSM.Template.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib_BSM.Template.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Modellib_BSM.Template.flavor_to_TeX: invalid down type quark"
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
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_CC -> "gcc"
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc" | G_Hmm -> "ghmm"
      | G_HGaZ -> "ghgaz" | G_HGaGa -> "ghgaga" | G_Hgg -> "ghgg"
      | G_H3 -> "gh3" | G_H4 -> "gh4" 
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f

  end

module HSExt (Flags : BSM_flags) =
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

    type matter_field = L of int | N of int | U of int | D of int
    type gauge_boson = Ga | Wp | Wm | Z | Gl 
    type other = Phip | Phim | Phi0 | H | S
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
      failwith "Modellib_BSM.HSExt.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; N n; U n; D n ]

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Gauge Bosons", List.map gauge_boson [Ga; Z; Wp; Wm; Gl];
        "Higgs", List.map other [H; S];
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
      | O f -> Scalar

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
          | H | S -> Prop_Scalar
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
          | H -> H | S -> S
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

    module Ch = Charges.QQ

    let ( // ) = Algebra.Small_Rational.make

    let generation' = function
      |  1 -> [ 1//1;  0//1;  0//1]
      |  2 -> [ 0//1;  1//1;  0//1]
      |  3 -> [ 0//1;  0//1;  1//1]
      | -1 -> [-1//1;  0//1;  0//1]
      | -2 -> [ 0//1; -1//1;  0//1]
      | -3 -> [ 0//1;  0//1; -1//1]
      |  n -> invalid_arg ("HSExt.generation': " ^ string_of_int n)

    let generation f =
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
          | H | Phi0 | S ->  0//1
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
      | Q_lepton | Q_up | Q_down | G_CC
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down
      | I_Q_W | I_G_ZWW
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ 
      | G_SWW | G_SSWW | G_SZZ | G_SSZZ | G_HSWW | G_HSZZ
      | G_Htt | G_Hbb | G_Hcc | G_Hmm | G_Htautau | G_H3 | G_H4_1
      | G_H4_2 | G_H4_3 | G_H4_4 | G_H4_5
      | G_Stt | G_Sbb | G_Scc | G_Smm | G_Stautau | G_HSS | G_HHS
      | G_HGaZ | G_HGaGa | G_Hgg | G_SGaZ | G_SGaGa | G_Sgg
      | Gs | I_Gs | G2
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let input_parameters = []

    let derived_parameters = [] 

    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

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

    let mgm ((m1, g, m2), fbf, c) = ((M m1, G g, M m2), fbf, c)

    let electromagnetic_currents n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]
        
    let color_currents n =
        List.map mgm
          [ ((U (-n), Gl, U n), FBF ((-1), Psibar, V, Psi), Gs);
            ((D (-n), Gl, D n), FBF ((-1), Psibar, V, Psi), Gs) ]

    let neutral_currents n =
      List.map mgm
        [ ((L (-n), Z, L n), FBF (1, Psibar, VA, Psi), G_NC_lepton);
          ((N (-n), Z, N n), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
          ((U (-n), Z, U n), FBF (1, Psibar, VA, Psi), G_NC_up);
          ((D (-n), Z, D n), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

    let charged_currents n = 
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC);  
          ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let yukawa =
      [ ((M (U (-3)), O H, M (U 3)), FBF (1, Psibar, Coupling.S, Psi), G_Htt);
        ((M (D (-3)), O H, M (D 3)), FBF (1, Psibar, Coupling.S, Psi), G_Hbb);
        ((M (U (-2)), O H, M (U 2)), FBF (1, Psibar, Coupling.S, Psi), G_Hcc);
        ((M (L (-2)), O H, M (L 2)), FBF (1, Psibar, Coupling.S, Psi), G_Hmm);
        ((M (L (-3)), O H, M (L 3)), 
	      FBF (1, Psibar, Coupling.S, Psi), G_Htautau);
	((M (U (-3)), O S, M (U 3)), FBF (1, Psibar, Coupling.S, Psi), G_Stt);
        ((M (D (-3)), O S, M (D 3)), FBF (1, Psibar, Coupling.S, Psi), G_Sbb);
        ((M (U (-2)), O S, M (U 2)), FBF (1, Psibar, Coupling.S, Psi), G_Scc);
        ((M (L (-2)), O S, M (L 2)), FBF (1, Psibar, Coupling.S, Psi), G_Smm);
        ((M (L (-3)), O S, M (L 3)), 
	      FBF (1, Psibar, Coupling.S, Psi), G_Stautau) ]

    let tgc ((g1, g2, g3), t, c) = ((G g1, G g2, G g3), t, c)

    let triple_gauge =
      List.map tgc
        [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);
          ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs) ]

    let qgc ((g1, g2, g3, g4), t, c) = ((G g1, G g2, G g3, G g4), t, c)

    let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
    let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]
    let quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          (Gl, Gl, Gl, Gl), gauge4, G2]

    let gauge_higgs =
      [ ((O H, G Wp, G Wm), Scalar_Vector_Vector 1, G_HWW);
        ((O H, G Z, G Z), Scalar_Vector_Vector 1, G_HZZ);
	((O S, G Wp, G Wm), Scalar_Vector_Vector 1, G_SWW);
        ((O S, G Z, G Z), Scalar_Vector_Vector 1, G_SZZ) ]

    let gauge_higgs4 =
      [ (O H, O H, G Wp, G Wm), Scalar2_Vector2 1, G_HHWW;
        (O H, O H, G Z, G Z), Scalar2_Vector2 1, G_HHZZ;
	(O H, O S, G Wp, G Wm), Scalar2_Vector2 1, G_HSWW;
        (O H, O S, G Z, G Z), Scalar2_Vector2 1, G_HSZZ;
	(O S, O S, G Wp, G Wm), Scalar2_Vector2 1, G_SSWW;
        (O S, O S, G Z, G Z), Scalar2_Vector2 1, G_SSZZ ]
       
    let higgs =
      [ (O H, O H, O H), Scalar_Scalar_Scalar 1, G_H3;
	(O S, O H, O H), Scalar_Scalar_Scalar 1, G_HHS;
	(O S, O S, O H), Scalar_Scalar_Scalar 1, G_HSS ]

    let higgs4 =
      [ (O H, O H, O H, O H), Scalar4 1, G_H4_1;
	(O H, O H, O H, O S), Scalar4 1, G_H4_2;
	(O H, O H, O S, O S), Scalar4 1, G_H4_3;
	(O H, O S, O S, O S), Scalar4 1, G_H4_4;
	(O S, O S, O S, O S), Scalar4 1, G_H4_5 ]

    let anomaly_higgs = 
      [ (O H, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
        (O H, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_HGaZ; 
        (O H, G Gl, G Gl), Dim5_Scalar_Gauge2 1, G_Hgg;
	(O S, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_SGaGa;
        (O S, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_SGaZ; 
        (O S, G Gl, G Gl), Dim5_Scalar_Gauge2 1, G_Sgg ]

    let goldstone_vertices =
      [ ((O Phi0, G Wm, G Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phip, G Ga, G Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phip, G Z, G Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phim, G Wp, G Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phim, G Wp, G Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @
       ThoList.flatmap color_currents [1;2;3] @
       ThoList.flatmap neutral_currents [1;2;3] @
       ThoList.flatmap charged_currents [1;2;3] @
       yukawa @ triple_gauge @ gauge_higgs @ higgs @ 
       anomaly_higgs @ goldstone_vertices)

    let vertices4 =
      quartic_gauge @ gauge_higgs4 @ higgs4

    let vertices () = (vertices3, vertices4, [])

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
      | "S" -> O S
      | _ -> invalid_arg "Modellib_BSM.HSExt.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib_BSM.HSExt.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib_BSM.HSExt.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib_BSM.HSExt.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_BSM.HSExt.flavor_to_string: invalid down type quark"
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
          | H -> "H" | S -> "S"
          end

    let flavor_to_TeX = function
      | M f ->
          begin match f with
          | L 1 -> "e^-" | L (-1) -> "e^+"
          | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
          | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
          | L _ -> invalid_arg
                "Modellib_BSM.HSExt.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib_BSM.HSExt.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib_BSM.HSExt.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Modellib_BSM.HSExt.flavor_to_TeX: invalid down type quark"
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
          | H -> "H" | S -> "S"
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
          | H -> "h" | S -> "s"
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
          | H -> 25 | S -> 35
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
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_CC -> "gcc"
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_SWW -> "gsww" | G_SZZ -> "gszz"
      | G_SSWW -> "gssww" | G_SSZZ -> "gsszz"
      | G_HSWW -> "ghsww" | G_HSZZ -> "ghszz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc" | G_Hmm -> "ghmm"
      | G_Stt -> "gstt" | G_Sbb -> "gsbb"
      | G_Stautau -> "gstautau" | G_Scc -> "gscc" | G_Smm -> "gsmm"
      | G_HGaZ -> "ghgaz" | G_HGaGa -> "ghgaga" | G_Hgg -> "ghgg"
      | G_SGaZ -> "gsgaz" | G_SGaGa -> "gsgaga" | G_Sgg -> "gsgg"
      | G_H3 -> "gh3" | G_H4_1 -> "gh4_1" | G_H4_2 -> "gh4_2" 
      | G_H4_3 -> "gh4_3" | G_H4_4 -> "gh4_4" | G_H4_5 -> "gh4_5" 
      | G_HHS -> "ghhs" | G_HSS -> "ghss"
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f

  end

(* \thocwmodulesection{Three-Site Higgsless Model} *)

module type Threeshl_options =
   sig
      val include_ckm: bool
      val include_hf: bool
      val diet: bool
   end

module Threeshl_no_ckm: Threeshl_options =
   struct
      let include_ckm = false
      let include_hf = true
      let diet = false
   end

module Threeshl_ckm: Threeshl_options =
   struct
      let include_ckm = true
      let include_hf = true
      let diet = false
   end

module Threeshl_no_ckm_no_hf: Threeshl_options =
   struct
      let include_ckm = false
      let include_hf = false
      let diet = false
   end

module Threeshl_ckm_no_hf: Threeshl_options =
   struct
      let include_ckm = true
      let include_hf = false
      let diet = false
   end

module Threeshl_diet_no_hf: Threeshl_options = 
   struct
      let include_ckm = false
      let include_hf = false
      let diet = true
   end

module Threeshl_diet: Threeshl_options =
   struct
      let include_ckm = false
      let include_hf = true
      let diet = true
   end

(* We use one generic implementation of the model and implement different features via option
modules given to a functor *)
module Threeshl (Module_options: Threeshl_options) =
   struct

      open Coupling

      let modname = "Modellib_BSM.Threeshl"

      (* Shamelessly stolen from Modellib.SM3, but with no support for fudged width yet *)
      let default_width = ref Timelike

      (* If this flag is set true, all gauge bosons are assumed to be massless and are assigned
      feynman gauge propagators. This in conjunction with the unbroken three site model is intended for
      checking gauge invariance via the ward identites. *)
      let all_feynman = ref false

      let options = Options.create [
         "constant_width", Arg.Unit (fun _ -> default_width := Constant),
            "use constant width (also in t-channel)";
         "custom_width", Arg.String (fun x -> default_width := Custom x),
            "use custom width";
         "cancel_widths", Arg.Unit (fun _ -> default_width := Vanishing),
            "use vanishing width";
         "all_feynman", Arg.Unit (fun _ -> all_feynman := true),
            "assign feynman gauge propagators to all gauge bosons\n"
            ^ "\t(for checking the ward identities); use only if you *really* know\n"
            ^ "\twhat you are doing"]
       let caveats () = []

      (* The quantum numbers that are carried by the particles. \verb$csign$ is \emph{not} the charge
      carried by the particle, but differentiates between particles (\verb$Pos$) and antiparticles
      (\verb$Neg$) *)
      type kkmode = Light | Heavy
      type generation = Gen0 | Gen1 | Gen2
      type csign = Pos | Neg
      type isospin = Iso_up | Iso_down

      (* Necessary to represent the indices of the couplings defined in FORTRAN *)
      type kk2 = Light2 | Heavy2 | Light_Heavy

      (* Map the different types to the constants used in the FORTRAN module *)
      let fspec_of_kkmode = function Light -> "l_mode" | Heavy -> "h_mode"
      let fspec_of_kk2 = function
         Light2 -> "l_mode" | Heavy2 -> "h_mode" | Light_Heavy -> "lh_mode"
      let fspec_of_gen = function Gen0 -> "gen_0" | Gen1 -> "gen_1" | Gen2 -> "gen_2"
      let fspec_of_iso = function Iso_up -> "iso_up" | Iso_down -> "iso_down"

      (* Covert the ``charge sign'' into a numeric sign (used e.g. in the determination of the MCID
      codes) *)
      let int_of_csign = function Pos -> 1 | Neg -> -1

      (* Convert the generation into an integer (dito) *)
      let int_of_gen = function Gen0 -> 1 | Gen1 -> 2 | Gen2 -> 3

      (* The type \verb$flavor$ is implemented as a variant. Fermions are implemented as a variant
      differentating between leptons and quarks (seemed the most natural way as this is also the way
      in which the FORTRAN code is structured). Bosons are implemented as a variant the
      differentiates between $W$, $Z$ and $A$. All other quantum numbers that are required for
      identifying the particles are carried by the variant constructors. *)
      type fermion = 
         | Lepton of (kkmode * csign * generation * isospin)
         | Quark of (kkmode * csign * generation * isospin)

      type boson =
         | W of (kkmode * csign)
         | Z of kkmode
         | A
         | G

      type flavor = Fermion of fermion | Boson of boson
   
      (* Helpers to construct particles from quantum numbers *)
      let lepton kk cs gen iso = Lepton (kk, cs, gen, iso)
      let quark kk cs gen iso = Quark (kk, cs, gen, iso)
      let w kk cs = W (kk, cs)
      let z kk = Z kk
      let flavor_of_f x = Fermion x
      let flavor_of_b x = Boson x

      (* Map a list of functions to the list (partially) applied to a value *)
      let revmap funs v = List.map (fun x -> x v) funs

      (* The same for a list of values; the result is flattened *)
      let revmap2 funs vals = ThoList.flatmap (revmap funs) vals

      (* Functions to loop the constructors over quantum numbers for list creation purposes *)
      let loop_kk flist = revmap2 flist [Light; Heavy]
      let loop_cs flist = revmap2 flist [Pos; Neg]
      let loop_gen flist = revmap2 flist [Gen0; Gen1; Gen2]
      let loop_iso flist = revmap2 flist [Iso_up; Iso_down]
      let loop_kk2 flist = revmap2 flist [Light2; Heavy2; Light_Heavy]

      (* Conditional looping over kk modes depending on whether to include heavy fermions *)
      let cloop_kk flist = match Module_options.include_hf with
         | true -> loop_kk flist
         | false -> revmap flist Light
      let cloop_kk2 flist = match Module_options.include_hf with
         | true -> loop_kk2 flist
         | false -> revmap flist Light2

      (* Having defined the necessary helpers, the magic of currying makes building lists of
      particles as easy as nesting the loop functions in the correct order... *)
      let all_leptons = loop_iso (loop_gen (loop_cs (cloop_kk [lepton] )))
      let all_quarks = loop_iso( loop_gen (loop_cs (cloop_kk [quark] )))
      let all_bosons = (loop_cs (loop_kk [w] )) @ [Z Light; Z Heavy; A; G]
      
      (* Converts a flavor spec to the BCD identifier defined in the FORTRAN module. Splitting the
      function into two parts \verb$prefix$ and \verb$rump$ removes a lot of redundancy. *)
      let bcdi_of_flavor = 
      let prefix = function
         | Fermion (Lepton (Heavy, _, _, _)) | Fermion (Quark (Heavy, _, _, _))
         | Boson (W (Heavy, _)) | Boson (Z Heavy) -> "h"
         | _ -> ""
      in let rump = function
         | Fermion (Lepton spec) -> (match spec with
            | (_, _, Gen0, Iso_up) -> "nue"
            | (_, _, Gen0, Iso_down) -> "e"
            | (_, _, Gen1, Iso_up) -> "numu"
            | (_, _, Gen1, Iso_down) -> "mu"
            | (_, _, Gen2, Iso_up) -> "nutau"
            | (_, _, Gen2, Iso_down) -> "tau")
         | Fermion (Quark spec) -> (match spec with
            | (_, _, Gen0, Iso_up) -> "u"
            | (_, _, Gen0, Iso_down) -> "d"
            | (_, _, Gen1, Iso_up) -> "c"
            | (_, _, Gen1, Iso_down) -> "s"
            | (_, _, Gen2, Iso_up) -> "t"
            | (_, _, Gen2, Iso_down) -> "b")
         | Boson (W _) -> "w" | Boson (Z _) -> "z"
         | Boson A -> invalid_arg (modname ^ ".bcd_of_flavor: no bcd for photon!")
         | Boson G -> invalid_arg (modname ^ ".bcd_of_flavor: no bcd for gluon!")
      in function x -> (prefix x) ^ (rump x) ^ "_bcd"

      (* The function defined in the model signature which returns the colour representation of a
      particle *)
      let color =
      let quarkrep = function
         | (_, Pos, _, _) -> Color.SUN 3
         | (_, Neg, _, _) -> Color.SUN (-3)
      in function
         | Fermion (Quark x) -> quarkrep x
         | Boson G -> Color.AdjSUN 3
         | _ -> Color.Singlet
      
    let nc () = 3

      (* Function for calculating the MCID code of a particle. Convenctions have been choosen such
      that the heavy modes are identified by the same numbers as the light ones, prefixed with
      \verb$99$. This is supposedly in accord with the conventions for adding new particles to the list
      of MCID codes. This function is required by the signature. *)
      let pdg =
      let iso_delta = function Iso_down -> 0 | Iso_up -> 1
      in let gen_delta = function Gen0 -> 0 | Gen1 -> 2 | Gen2 -> 4
      in let kk_delta = function Light -> 0 | Heavy -> 9900
      in function
         | Fermion ( Lepton (kk, cs, gen, iso)) ->
            (int_of_csign cs) * (11 + (gen_delta gen) + (iso_delta iso) + (kk_delta kk))
         | Fermion ( Quark (kk, cs, gen, iso)) -> 
            (int_of_csign cs) * (1 + (gen_delta gen) + (iso_delta iso)+ (kk_delta kk))
         | Boson (W (kk, cs)) -> (int_of_csign cs) * (24 + (kk_delta kk))
         | Boson (Z kk) -> 23 + (kk_delta kk)
         | Boson A -> 22
         | Boson G -> 21

      (* Returns the lorentz representation of a particle; required by the signature. *)
      let lorentz = 
      let spinor = function
         | (_, Pos, _, _) -> Spinor
         | (_, Neg, _, _) -> ConjSpinor
      in function
         | Fermion (Lepton x) | Fermion (Quark x) -> spinor x
         | Boson (W _) | Boson (Z _) -> Massive_Vector
         | Boson A -> Vector
         | Boson G -> Vector

      (* O'Mega supports models that allow different gauges; however, we only implement unitary
      gauge and therefore stub this (SM3 does the same thing). The \verb$gauge$ type as well as
      \verb$gauge_symbol$ are required by the signature. *)
      type gauge = unit

      let gauge_symbol () =
         failwith (modname ^ ".gauge_symbol: internal error")

      (* Returns the propagator for a given particle type. Required by signature. *)
      let propagator =
      let spinorprop = function
         | (_, Pos, _, _) -> Prop_Spinor
         | (_, Neg, _, _) -> Prop_ConjSpinor
      in function 
         | Fermion (Lepton x) | Fermion (Quark x) -> spinorprop x
         | Boson (W _) | Boson (Z _) ->
            (match !all_feynman with false -> Prop_Unitarity | true -> Prop_Feynman)
         | Boson A -> Prop_Feynman
         | Boson G -> Prop_Feynman

      (* Return the width of a particle, required by signature. \\
      \emph{TODO:} Refine such that stable particles always are treade via vanishing width, as this
      might speed up the generated code a bit. *)
      let width _ = !default_width

      (* Returns the conjugate particle; required by signature. *)
      let conjugate =
      let conj_csign = function
         | Pos -> Neg
         | Neg -> Pos
      in function
         | Fermion (Lepton (kk, cs, gen, iso)) -> Fermion (Lepton (kk, conj_csign cs, gen, iso))
         | Fermion (Quark (kk, cs, gen, iso)) -> Fermion (Quark (kk, conj_csign cs, gen, iso))
         | Boson (W (kk, cs)) -> Boson (W (kk, conj_csign cs))
         |  x -> x

      (* Tells the diagram generator whether a particle is a fermion, a conjugate fermion or a
      boson. Required by signature *)
      let fermion = function
         | Fermion (Lepton (_, cs, _, _)) | Fermion (Quark (_, cs, _, _)) -> int_of_csign cs
         | Boson _ -> 0

      (* Charges are: charge, lepton number, baryon number, generation. Required by signature *)
      module Ch = Charges.QQ
      let ( // ) = Algebra.Small_Rational.make

      let qn_charge = function
         | Boson b -> (match b with
            | W (_, c) -> (int_of_csign (c)) // 1
            | _ -> 0//1)
         | Fermion f -> (match f with
            | Lepton (_, c, _, Iso_up) -> 0//1
            | Lepton (_, c, _, Iso_down) -> (-1 * int_of_csign (c)) // 1
            | Quark (_, c, _, Iso_up) -> (2 * int_of_csign (c)) // 3
            | Quark (_, c, _, Iso_down) -> (-1 * int_of_csign (c)) // 3)

      let qn_lepton = function
         | Fermion (Lepton (_, c, _, _)) -> int_of_csign (c) // 1
         | _ -> 0//1

      let qn_baryon = function
         | Fermion (Quark (_, c, _, _)) -> int_of_csign (c) // 1
         | _ -> 0//1

      (* Generation is conditional: if we enable the nontrivial CKM matrix, all particles carry
      generation [0; 0; 0] *)
      let qn_generation x =
         let qn cs gen =
            let c = int_of_csign (cs) in
            match gen with
               | Gen0 -> [c//1; 0//1; 0//1]
               | Gen1 -> [0//1; c//1; 0//1]
               | Gen2 -> [0//1; 0//1; c//1]
         in
         if Module_options.include_ckm then
            [0//1; 0//1; 0//1]
         else
            match x with
               | Fermion (Lepton (_, c, g, _)) -> qn c g
               | Fermion (Quark (_, c, g, _)) -> qn c g
               | _ -> [0//1; 0//1; 0//1]

      let charges x =
         [qn_charge x; qn_lepton x; qn_baryon x] @ (qn_generation x)

      (* A variant to represent the different coupling constants, choosen to mimic the FORTRAN part.
      Required by signature. *)
      type constant =
         | G_a_lep | G_a_quark of isospin
         | G_aww | G_aaww
         | G_w_lep of (kkmode * kkmode * generation * kkmode * generation)
         | G_w_quark of (kkmode * kkmode * generation * kkmode * generation)
         | G_z_lep of (kkmode * kk2 * generation * isospin)
         | G_z_quark of (kkmode * kk2 * generation * isospin)
         | G_wwz of (kk2 * kkmode)
         | G_wwzz of (kk2 * kk2)
         | G_wwza of (kk2 * kkmode)
         | G_wwww of int
         | G_s
         | IG_s
         | G_s2

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

      (* Functions for the construction of constants from indices *)
      let g_a_quark x = G_a_quark x
      let g_w_lep kk1 kk2 gen1 kk3 gen2 = G_w_lep (kk1, kk2, gen1, kk3, gen2)
      let g_w_quark kk1 kk2 gen1 kk3 gen2 = G_w_quark (kk1, kk2, gen1, kk3, gen2)
      let g_z_lep kk1 kk2 gen iso = G_z_lep (kk1, kk2, gen, iso)
      let g_z_quark kk1 kk2 gen iso = G_z_quark (kk1, kk2, gen, iso)
      let g_wwz kk1 kk2 = G_wwz (kk1, kk2)
      let g_wwzz kk1 kk2 = G_wwzz (kk1, kk2)
      let g_wwza kk1 kk2 = G_wwza (kk1, kk2)
      let g_wwww nhw = if (nhw >= 0) && (nhw <= 4) then G_wwww nhw
            else failwith (modname ^ ".g_wwww: invalid integer, very bad")

      (* Build a list of the different constants *)
      let clist = [G_a_lep; G_aww; G_aaww] @ (loop_iso [g_a_quark]) @
         (loop_gen (cloop_kk (loop_gen (cloop_kk (loop_kk [g_w_lep] ))))) @
         (loop_gen (cloop_kk (loop_gen (cloop_kk (loop_kk [g_w_quark] ))))) @
         (loop_iso (loop_gen (cloop_kk2 (loop_kk [g_z_lep] )))) @
         (loop_iso (loop_gen (cloop_kk2 (loop_kk [g_z_quark] )))) @
         (loop_kk (loop_kk2 [g_wwz] )) @ (loop_kk2 (loop_kk2 [g_wwzz] )) @
         (loop_kk (loop_kk2 [g_wwza] )) @ (List.map g_wwww [0; 1; 2; 3; 4])

      (* Maximum number of lines meeting at a vertex, required by signature. *)
      let max_degree () = 4

      (* Transform a pair of kk identifiers into a kk2 identifier *)
      let get_kk2 = function (Light, Light) -> Light2 | (Heavy, Heavy) -> Heavy2
         | (Light, Heavy) | (Heavy, Light) -> Light_Heavy

      (* Flip isospin *)
      let conj_iso = function Iso_up -> Iso_down | Iso_down -> Iso_up

      (* Below, lists of couplings are generated which ultimately are joined into a list of all
      couplings in the model. The generated lists can be viewed using the \verb$dump.ml$ script in the
      O'Mega toplevel directory. \\
      The individual couplings are defined as 5-tupels resp. 6-tupels consisting in this
      order of the particles meeting at the vertex, the coupling type (see \verb$couplings.ml$) and the
      coupling constant. *)

      (* List of $llA$ type vertices *)
      let vertices_all =
      let vgen kk gen =
         ((Fermion (Lepton (kk, Neg, gen, Iso_down)), Boson A, Fermion (Lepton (kk, Pos, gen,
            Iso_down))), FBF(1, Psibar, V, Psi), G_a_lep)
      in loop_gen (cloop_kk [vgen])

      (* List of $qqA$ type vertices *)
      let vertices_aqq =
      let vgen kk gen iso =
         ((Fermion (Quark (kk, Neg, gen, iso)), Boson A, Fermion (Quark (kk, Pos, gen,
            iso))), FBF(1, Psibar, V, Psi), G_a_quark iso)
      in loop_iso (loop_gen (cloop_kk [vgen]))


      (* List of $\nu lW$ type vertices *)
      let vertices_wll =
      let vgen kkw kk_f kk_fbar iso_f gen =
         ((Fermion (Lepton (kk_fbar, Neg, gen, conj_iso iso_f)),
            Boson (W (kkw, (match iso_f with Iso_up -> Neg | _ -> Pos))),
            Fermion (Lepton (kk_f, Pos, gen, iso_f))),
            FBF (1, Psibar, VA2, Psi),
            G_w_lep (kkw, (match iso_f with Iso_up -> kk_f | _ -> kk_fbar), gen,
               (match iso_f with Iso_up -> kk_fbar | _ -> kk_f), gen) )
      in loop_gen (loop_iso (cloop_kk (cloop_kk (loop_kk [vgen] ))))

      (* The same list, but without couplings between the $W^\prime$ and light fermions *)
      let vertices_wll_diet =
      let filter = function
         | ((Fermion (Lepton (Light, _, _, _)), Boson (W (Heavy, _)),
            Fermion (Lepton (Light, _, _, _))), _, _) -> false
         | _ -> true
      in List.filter filter vertices_wll

      (* List of $udW$ type vertices, flavor-diagonal *)
      let vertices_wqq_no_ckm =
      let vgen kkw kk_f kk_fbar iso_f gen =
         ((Fermion (Quark (kk_fbar, Neg, gen, conj_iso iso_f)),
            Boson (W (kkw, (match iso_f with Iso_up -> Neg | _ -> Pos))),
            Fermion (Quark (kk_f, Pos, gen, iso_f))),
            FBF (1, Psibar, VA2, Psi),
            G_w_quark (kkw, (match iso_f with Iso_up -> kk_f | _ -> kk_fbar), gen,
               (match iso_f with Iso_up -> kk_fbar | _ -> kk_f), gen) )
      in loop_gen (loop_iso (cloop_kk (cloop_kk (loop_kk [vgen] ))))

      (* The same list, but without couplings between the $W^\prime$ and the first two generations
      of quarks *)
      let vertices_wqq_no_ckm_diet =
      let filter = function
         | ((Fermion (Quark (Light, _, gen, _)), Boson (W (Heavy, _)),
            Fermion (Quark (Light, _, _, _))), _, _) -> 
               (match gen with Gen2 -> true | _ -> false)
         | _ -> true
      in List.filter filter vertices_wqq_no_ckm

      (* List of $udW$ type vertices, including non flavor-diagonal couplings *)
      let vertices_wqq =
      let vgen kkw kk_f gen_f kk_fbar gen_fbar iso_f =
         ((Fermion (Quark (kk_fbar, Neg, gen_fbar, conj_iso iso_f)),
            Boson (W (kkw, (match iso_f with Iso_up -> Neg | _ -> Pos))),
            Fermion (Quark (kk_f, Pos, gen_f, iso_f))),
            FBF (1, Psibar, VA2, Psi),
            G_w_quark (match iso_f with
               | Iso_up -> (kkw, kk_f, gen_f, kk_fbar, gen_fbar)
               | Iso_down -> (kkw, kk_fbar, gen_fbar, kk_f, gen_f)))
      in loop_iso (loop_gen (cloop_kk (loop_gen (cloop_kk (loop_kk [vgen] )))))


      (* List of $llZ$ / $\nu\nu Z$ type vertices *)
      let vertices_zll =
      let vgen kkz kk_f kk_fbar gen iso =
         ((Fermion (Lepton (kk_fbar, Neg, gen, iso)), Boson (Z kkz),
            Fermion (Lepton (kk_f, Pos, gen, iso))),
            FBF (1, Psibar, VA2, Psi),
            G_z_lep (kkz, get_kk2 (kk_f, kk_fbar), gen, iso))
      in loop_iso (loop_gen (cloop_kk (cloop_kk (loop_kk [vgen] ))))

      (* List of $qqZ$ type vertices *)
      let vertices_zqq =
      let vgen kkz kk_f kk_fbar gen iso =
         ((Fermion (Quark (kk_fbar, Neg, gen, iso)), Boson (Z kkz),
            Fermion (Quark (kk_f, Pos, gen, iso))),
            FBF (1, Psibar, VA2, Psi),
            G_z_quark (kkz, get_kk2 (kk_f, kk_fbar), gen, iso))
      in loop_iso (loop_gen (cloop_kk (cloop_kk (loop_kk [vgen] ))))

      (* $gq\bar{q}$ *)
      let vertices_gqq =
      let vgen kk gen iso =
         ((Fermion (Quark (kk, Neg, gen, iso)), Boson G, Fermion (Quark (kk, Pos, gen, iso))),
            FBF (1, Psibar, V, Psi), G_s)
      in loop_iso (loop_gen (cloop_kk [vgen]))

      (* AWW *)
      let vertices_aww =
      let vgen kk =
         ( (Boson A, Boson (W (kk, Pos)), Boson (W (kk, Neg))), Gauge_Gauge_Gauge 1, G_aww)
      in loop_kk [vgen]

      (* ZWW *)
      let vertices_zww =
      let vgen kkz kkwp kkwm =
         ((Boson (Z kkz), Boson (W (kkwp, Pos)), Boson (W (kkwm, Neg))), Gauge_Gauge_Gauge 1, 
            G_wwz (get_kk2 (kkwp, kkwm), kkz))
      in loop_kk (loop_kk (loop_kk [vgen]))

      (* $ggg$ *)
      let vertices_ggg = [(Boson G, Boson G, Boson G), Gauge_Gauge_Gauge (-1), IG_s]

      (* Stolen from Modellib.SM; the signs seem to be OK. See \verb$couplings.ml$ for more docs. *)
      let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
      let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)]

      (* AAWW *)
      let vertices_aaww =
      let vgen kk =
         ((Boson A, Boson (W (kk, Pos)), Boson A, Boson (W (kk, Neg))), minus_gauge4, G_aaww)
      in loop_kk [vgen]

      (* WWZZ *)
      let vertices_wwzz =
      let vgen kkwp kkwm kk2z =
         ((Boson (Z (match kk2z with Heavy2 -> Heavy | Light2 | Light_Heavy -> Light)),
            Boson (W (kkwp, Pos)),
            Boson (Z (match kk2z with Heavy2 | Light_Heavy -> Heavy | Light2 -> Light)),
            Boson (W (kkwm, Neg))), minus_gauge4, G_wwzz (get_kk2 (kkwp, kkwm), kk2z))
      in loop_kk2 (loop_kk (loop_kk [vgen]))

      (* WWZA *)
      let vertices_wwza =
      let vgen kkwp kkwm kkz =
         ((Boson A, Boson (W (kkwp, Pos)), Boson (Z kkz), Boson (W (kkwm, Neg))),
            minus_gauge4, G_wwza (get_kk2 (kkwp, kkwm), kkz))
      in loop_kk (loop_kk (loop_kk [vgen]))

      (* WWWW *)
      let vertices_wwww =
      let count = function Light2 -> 0 | Light_Heavy -> 1 | Heavy2 -> 2
      in let vgen kk2wp kk2wm =
         ((Boson (W ((match kk2wp with Heavy2 -> Heavy | Light2 | Light_Heavy -> Light), Pos)),
            Boson (W ((match kk2wm with Heavy2 -> Heavy | Light2 | Light_Heavy -> Light), Neg)),
            Boson (W ((match kk2wp with Heavy2 | Light_Heavy -> Heavy | Light2 -> Light), Pos)),
            Boson (W ((match kk2wm with Heavy2 | Light_Heavy -> Heavy | Light2 -> Light), Neg))),
            gauge4, G_wwww ((count kk2wp) + (count kk2wm)))
      in loop_kk2 (loop_kk2 [vgen])

      (* gggg *)
      let vertices_gggg = [(Boson G, Boson G, Boson G, Boson G), gauge4, G_s2]

      (* The list of couplings is transformed into the fusion lists required by the generator by
      the Model.Fusions functor. *)

      (* This is copy\& paste from the other models; check again with Thorsten if it is correct *)
      module F = Modeltools.Fusions (struct
         type f = flavor
         type c = constant
         let compare = compare
         let conjugate = conjugate
      end )

      (* Not sure yet whether F.fusex also creates the conjugate vertices; by looking at the
      implementation of the other models, I assume it doesn't. Still, better ask Thorsten to be
      sure!!!\\
      \emph{Update:} Still didn't get to ask, but since the results are consistent, I suspect my assertion
      is correct. \\
      The stuff below is required by the signature. *)

      let vertices () = (vertices_all @ vertices_aqq @ 
         (match Module_options.diet with
            | false -> vertices_wll
            | true -> vertices_wll_diet) @
         (match (Module_options.include_ckm, Module_options.diet) with
            | (true, false) -> vertices_wqq
            | (false, false) -> vertices_wqq_no_ckm
            | (false, true) -> vertices_wqq_no_ckm_diet
            | (true, true) -> raise (Failure
               ("Modules4.Threeshl.vertices: CKM matrix together with option diet is not" ^
               " implemented yet!"))) @
         vertices_zll @ vertices_zqq @ vertices_aww @ vertices_zww @ vertices_gqq @ vertices_ggg,
         vertices_aaww @ vertices_wwzz @ vertices_wwza @ vertices_wwww @ vertices_gggg
         , [])
      let table = F.of_vertices (vertices ())
      let fuse2 = F.fuse2 table
      let fuse3 = F.fuse3 table
      let fuse = F.fuse table

      (* A function that returns a list of a flavours known to the model, required by the signature.
      *)
      let flavors () = (List.map flavor_of_f (all_leptons @ all_quarks)) @
         (List.map flavor_of_b all_bosons)

      (* dito, external flavours, also required. *)
      let external_flavors () = [
         "light leptons", List.map flavor_of_f (loop_iso (loop_gen( loop_cs [lepton Light])));
         "light quarks", List.map flavor_of_f (loop_iso (loop_gen( loop_cs [quark Light])));
         "light gauge bosons", List.map flavor_of_b [W (Light, Pos); W (Light, Neg); Z Light; A];
         "heavy gauge bosons", List.map flavor_of_b [W (Heavy, Pos); W (Heavy, Neg); Z Heavy]] @
         (match Module_options.include_hf with
            | true -> [
               "heavy leptons", List.map flavor_of_f (loop_iso (loop_gen( loop_cs [lepton Heavy])));
               "heavy quarks", List.map flavor_of_f (loop_iso (loop_gen( loop_cs [quark Heavy])))]
            | false -> [] ) @ ["gluons", [Boson G]]

      (* Which of the particles are goldstones? $\rightarrow$ none. Required by the signature. *)
      let goldstone x = None

      (* This is wrong but handy for debugging the constant identifier generation via -params.
      Usually, this function would return a record consisting of the parameters as well as
      expression for the dependent quantities that can be used to generate FORTRAN code for
      calculating them. However, we have a seperate module for the threeshl, so we can abuse this
      for debugging. Required by signature. *)
      let parameters () = {input = List.map (fun x -> (x, 0.)) clist;
         derived = []; derived_arrays = []}

      (* Convert a flavour into a ID string with which it will be referred by the user interface of
      the compiled generator. Required by signature *)
      let flavor_to_string =
      let prefix = function
         | Fermion (Lepton (Heavy, _, _, _)) | Fermion (Quark (Heavy, _, _, _))
         | Boson (W (Heavy, _)) | Boson (Z Heavy) -> "H"
         | _ -> ""
      in let postfix = function
         | Fermion (Lepton (_, cs, _, Iso_down)) -> (match cs with Pos -> "-" | Neg -> "+")
         | Fermion (Quark (_, Neg, _, _)) | Fermion (Lepton (_, Neg, _, Iso_up)) -> "bar"
         | Boson (W (_, cs)) -> (match cs with Pos -> "+" | Neg -> "-")
         | _ -> ""
      in let rump = function
         | Fermion (Lepton desc) -> (match desc with
            | (_, _, Gen0, Iso_up) -> "nue"
            | (_, _, Gen0, Iso_down) -> "e"
            | (_, _, Gen1, Iso_up) -> "numu"
            | (_, _, Gen1, Iso_down) -> "mu"
            | (_, _, Gen2, Iso_up) -> "nutau"
            | (_, _, Gen2, Iso_down) -> "tau")
         | Fermion (Quark desc) -> (match desc with
            | (_, _, Gen0, Iso_up) -> "u"
            | (_, _, Gen0, Iso_down) -> "d"
            | (_, _, Gen1, Iso_up) -> "c"
            | (_, _, Gen1, Iso_down) -> "s"
            | (_, _, Gen2, Iso_up) -> "t"
            | (_, _, Gen2, Iso_down) -> "b")
         | Boson (W _) -> "W" | Boson (Z _) -> "Z" | Boson A -> "A" | Boson G -> "gl"
      in function x -> (prefix x) ^ (rump x) ^ (postfix x)

      (* Conversion of the ID string into a particle flavor. Instead of going through all cases
      again, we generate a ``dictionary'' of flavor / ID pairs which we use to identify the correct
      flavor. Required by signature. *)
      let flavor_of_string x =
      let dict = List.map (fun x -> (x, flavor_to_string x)) (flavors ())
      in let get_ident = function (x, _) -> x
      in try
            get_ident (List.find (fun (_, y) -> (x = y)) dict)
         with
            Not_found -> (match x with
               | "g" -> Boson G
               | _ -> invalid_arg (modname ^ ".flavor_of_string")
            )

      (* Converts a flavor into a symbol used as identification in the generated FORTRAN code (has
      to comply to the conventions of valid FORTRAN identifiers therefore). We stick to the same
      convenctions as SM3, prefixing heavy modes with a \verb$H$. Required by signature. *)
      let flavor_symbol =
      let prefix = function
         | Fermion (Lepton (Heavy, _, _, _)) | Fermion (Quark (Heavy, _, _, _))
         | Boson (W (Heavy, _)) | Boson (Z Heavy) -> "H"
         | _ -> ""
      in let postfix = function
         | Fermion (Lepton (_, Neg, _, _)) | Fermion (Quark (_, Neg, _, _)) -> "b"
         | _ -> ""
      in let rump = function
         | Fermion spec -> (match spec with
            | Lepton (_, _, gen, Iso_up) -> "n" ^ (string_of_int (int_of_gen gen))
            | Lepton (_, _, gen, Iso_down) -> "l" ^ (string_of_int (int_of_gen gen))
            | Quark (_, _, gen, Iso_up) -> "u" ^ (string_of_int (int_of_gen gen))
            | Quark (_, _, gen, Iso_down) -> "d"^ (string_of_int (int_of_gen gen)))
         | Boson spec -> (match spec with
            | W (_, Pos) -> "wp" | W (_, Neg) -> "wm"
            | Z _ -> "z" | A -> "a" | G -> "gl" )
      in function
         x -> (prefix x) ^ (rump x) ^ (postfix x)

      (* Generate TeX for a flavor *)
      let flavor_to_TeX =
      let bar x y = match  x with Neg -> "\\overline{" ^ y ^ "}" | Pos -> y
      in let pm x y = match x with Neg -> "{" ^ y ^ "}^+" | Pos -> "{" ^ y ^ "}^-"
      in let prime x y = match x with Light -> y | Heavy -> "{" ^ y ^ "}^\\prime"
      in function
         | Fermion (Lepton desc) -> (match desc with
            | (kk, cs, gen, Iso_up) -> prime kk (bar cs (match gen with
               | Gen0 -> "\\nu_e"
               | Gen1 -> "\\nu_\\mu"
               | Gen2 -> "\\nu_\\tau"))
            | (kk, cs, gen, Iso_down) -> prime kk (pm cs (match gen with
               | Gen0 -> "e" | Gen1 -> "\\mu" | Gen2 -> "\\tau")))
         | Fermion (Quark (kk, cs, gen, iso)) -> prime kk (bar cs (match (gen, iso) with
            | (Gen0, Iso_up) -> "u"
            | (Gen0, Iso_down) -> "d"
            | (Gen1, Iso_up) -> "c"
            | (Gen1, Iso_down) -> "s"
            | (Gen2, Iso_up) -> "t"
            | (Gen2, Iso_down) -> "b"))
         | Boson spec -> (match spec with
            | W (kk, cs) -> prime kk (pm (match cs with Pos -> Neg | Neg -> Pos) "W")
            | Z kk -> prime kk "Z"
            | A -> "A" | G -> "g")
         
      (* Returns the string referring to the particle mass in the generated FORTRAN code. Required
      by signature. *)
      let mass_symbol = function
         | Boson A | Boson G-> "0._default"
         | x -> "mass_array(" ^ (bcdi_of_flavor x) ^ ")"

      (* Dito, for width. Required by signature. *)
      let width_symbol = function
         | Boson A | Boson G -> "0._default"
         | x -> "width_array(" ^ (bcdi_of_flavor x) ^ ")"
      
      (* Determines the string referring to a coupling constant in the generated FORTRAN code.
      Required by signature. *)
      let constant_symbol =
      let c = ", "
      in let g_w_ferm = function
         (kk1, kk2, gen1, kk3, gen2) ->
            ":, " ^ (fspec_of_kkmode kk1) ^ c ^ (fspec_of_kkmode kk2) ^ c ^ (fspec_of_gen gen1) ^ c ^
            (fspec_of_kkmode kk3) ^ c ^ (fspec_of_gen gen2)
      in let g_z_ferm = function
         (kk1, kk2, gen, iso) ->
            ":, " ^ (fspec_of_kkmode kk1) ^ c ^ (fspec_of_kk2 kk2) ^ c ^ (fspec_of_gen gen) ^ c ^
            (fspec_of_iso iso)
      in function
         | G_a_lep -> "g_a_lep"
         | G_s -> "g_s_norm"
         | IG_s -> "ig_s_norm"
         | G_s2 -> "g_s_norm2"
         | G_a_quark iso -> "g_a_quark(" ^ (fspec_of_iso iso) ^ ")"
         | G_aww -> "ig_aww"
         | G_aaww -> "g_aaww"
         | G_w_lep spec -> "g_w_lep_va(" ^ (g_w_ferm spec) ^ ")"
         | G_w_quark spec -> "g_w_quark_va(" ^ (g_w_ferm spec) ^ ")"
         | G_z_lep spec -> "g_z_lep_va(" ^ (g_z_ferm spec) ^ ")"
         | G_z_quark spec -> "g_z_quark_va(" ^ (g_z_ferm spec) ^ ")"
         | G_wwz (kk1, kk2) -> "ig_wwz(" ^ (fspec_of_kk2 kk1) ^ c ^
            (fspec_of_kkmode kk2) ^ ")"
         | G_wwzz (kk1, kk2) -> "g_wwzz(" ^ (fspec_of_kk2 kk1) ^ c ^
            (fspec_of_kk2 kk2) ^ ")"
         | G_wwza (kk1, kk2) -> "g_wwza(" ^(fspec_of_kk2 kk1) ^ c ^
            (fspec_of_kkmode kk2) ^ ")"
         | G_wwww nhw -> if (0 <= nhw) && (nhw <= 4) then
            "g_wwww(" ^ (string_of_int nhw) ^ ")"
            else failwith "Modules4.Threeshl.constant_symbol: invalid int for G_wwww; very bad"

   end

(* \thocwmodulesection{THDM with and without non-trivial flavor structure} *)

module type THDM_flags =
  sig
    val ckm_present   : bool
  end

module THDM : THDM_flags = 
  struct
    let ckm_present = false
  end

module THDM_CKM : THDM_flags =
  struct
    let ckm_present = true
  end

module TwoHiggsDoublet (Flags : THDM_flags) =
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

    type matter_field = L of int | N of int | U of int | D of int
    type gauge_boson = Ga | Wp | Wm | Z | Gl
    type other = Phip | Phim | Phi0 | Hh | HA | HH | Hp | Hm
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
      failwith "Modellib_BSM.TwoHiggsDoublet.gauge_symbol: internal error"

    let family n = List.map matter_field [ L n; N n; U n; D n ]

    let external_flavors () =
      [ "1st Generation", ThoList.flatmap family [1; -1];
        "2nd Generation", ThoList.flatmap family [2; -2];
        "3rd Generation", ThoList.flatmap family [3; -3];
        "Gauge Bosons", List.map gauge_boson [Ga; Z; Wp; Wm; Gl];
        "Higgs", List.map other [Hh; HH; HA; Hp; Hm];
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
      | O f -> Scalar

    let color = function 
      | M (U n) -> Color.SUN (if n > 0 then 3 else -3)
      | M (D n) -> Color.SUN  (if n > 0 then 3 else -3)
      | G Gl  -> Color.AdjSUN 3
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
          | Hh | HH | HA | Hp | Hm -> Prop_Scalar
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
          | Hh -> Hh | HH -> HH | HA -> HA
	  | Hp -> Hm | Hm -> Hp
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
      |  n -> invalid_arg ("Modellib_BSM.TwoHiggsDoublet.generation': " ^ string_of_int n)

    let generation f =
      if Flags.ckm_present then
	[]
      else
	match f with
	  | M (L n | N n | U n | D n) -> generation' n
	  | G _ | O _ -> [ 0//1; 0//1; 0//1] 

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
          | Hh | HH | HA | Phi0 ->  0//1
          | Hp | Phip ->  1//1
          | Hm | Phim -> -1//1
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
      | Q_lepton | Q_up | Q_down | G_CC
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down
      | I_Q_W | I_G_ZWW | I_G_WWW
      | G_WWWW | G_ZZWW | G_AZWW | G_AAWW
      | G_htt | G_hbb | G_hcc | G_htautau | G_hmumu
      | G_Htt | G_Hbb | G_Hcc | G_Htautau | G_Hmumu
      | I_G_Att | I_G_Abb | I_G_Acc | I_G_Atautau | I_G_Amumu
      | G_Htb | G_Hcs | G_Htaunu | G_Hmunu
      | Gs | I_Gs | G2
      | G_AHpHm | G_ZHpHm | G_Zh1h2 | G_Zh1h3 | G_Zh2h3 
      | G_WpHmh1 | G_WpHmh2 | G_WpHmh3 | G_WmHph1 | G_WmHph2 
      | G_WmHph3 | G_h1ZZ | G_h2ZZ | G_h3ZZ | G_h1WpWm | G_h2WpWm 
      | G_h3WpWm | G_hhWpWm | G_hhZZ | G_HpHmAA | G_HpHmZZ | G_HpHmAZ 
      | G_HpHmWpWm | G_h1HpAWm | G_h2HpAWm | G_h3HpAWm | G_h1HpZWm 
      | G_h2HpZWm | G_h3HpZWm | G_h1HpAWmC | G_h2HpAWmC | G_h3HpAWmC 
      | G_h1HpZWmC | G_h2HpZWmC | G_h3HpZWmC | G_h1HpHm | G_h2HpHm | G_h3HpHm 
      | G_h111 | G_h112 | G_h113 | G_h221 | G_h222 | G_h223 | G_h331 | G_h332 
      | G_h333 | G_h123 | G_HpHmHpHm | G_HpHm11 | G_HpHm12 | G_HpHm13 
      | G_HpHm22 | G_HpHm23 | G_HpHm33 | G_h1111 | G_h1112 | G_h1113 
      | G_h1122 | G_h1123 | G_h1133 | G_h1222 | G_h1223 | G_h1233 | G_h1333 
      | G_h2222 | G_h2223 | G_h2233 | G_h2333 | G_h3333 
      | G_h1uu | G_h2uu | G_h3uu | G_h1uc | G_h2uc | G_h3uc | G_h1ut
      | G_h2ut | G_h3ut | G_h1cu | G_h2cu | G_h3cu | G_h1cc | G_h2cc
      | G_h3cc | G_h1ct | G_h2ct | G_h3ct | G_h1tu | G_h2tu | G_h3tu
      | G_h1tc | G_h2tc | G_h3tc | G_h1tt | G_h2tt | G_h3tt
      | G_h1dd | G_h2dd | G_h3dd | G_h1ds | G_h2ds | G_h3ds | G_h1db
      | G_h2db | G_h3db | G_h1sd | G_h2sd | G_h3sd | G_h1ss | G_h2ss
      | G_h3ss | G_h1sb | G_h2sb | G_h3sb | G_h1bd | G_h2bd | G_h3bd
      | G_h1bs | G_h2bs | G_h3bs | G_h1bb | G_h2bb | G_h3bb
      | G_hud | G_hus | G_hub | G_hcd | G_hcs | G_hcb | G_htd | G_hts | G_htb
      | G_hdu | G_hdc | G_hdt | G_hsu | G_hsc | G_hst | G_hbu | G_hbc | G_hbt
      | G_he1n1 | G_he1n2 | G_he1n3 | G_he2n1 | G_he2n2 | G_he2n3 | G_he3n1
      | G_he3n2 | G_he3n3 | G_hn1e1 | G_hn1e2 | G_hn1e3 | G_hn2e1 | G_hn2e2
      | G_hn2e3 | G_hn3e1 | G_hn3e2 | G_hn3e3
      | G_h1e1e1 | G_h2e1e1 | G_h3e1e1 | G_h1e1e2 | G_h2e1e2 | G_h3e1e2
      | G_h1e1e3 | G_h2e1e3 | G_h3e1e3 | G_h1e2e1 | G_h2e2e1 | G_h3e2e1
      | G_h1e2e2 | G_h2e2e2 | G_h3e2e2 | G_h1e2e3 | G_h2e2e3 | G_h3e2e3
      | G_h1e3e1 | G_h2e3e1 | G_h3e3e1 | G_h1e3e2 | G_h2e3e2 | G_h3e3e2
      | G_h1e3e3 | G_h2e3e3 | G_h3e3e3
      | Mass of flavor | Width of flavor

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let g_over_2_costh =
      Quot (Neg (Atom G_weak), Prod [Integer 2; Atom Costhw])

    let nc_coupling c t3 q =
      (Real_Array c,
       [Prod [g_over_2_costh; Diff (t3, Prod [Integer 2; q; Atom Sin2thw])];
        Prod [g_over_2_costh; t3]])

    let half = Quot (Integer 1, Integer 2)

    let array_list = [G_h1uu; G_h2uu; G_h3uu; G_h1uc; G_h2uc; G_h3uc; G_h1ut;
        G_h2ut; G_h3ut; G_h1cu; G_h2cu; G_h3cu; G_h1cc; G_h2cc;
        G_h3cc; G_h1ct; G_h2ct; G_h3ct; G_h1tu; G_h2tu; G_h3tu;
        G_h1tc; G_h2tc; G_h3tc; G_h1tt; G_h2tt; G_h3tt;
        G_h1dd; G_h2dd; G_h3dd; G_h1ds; G_h2ds; G_h3ds; G_h1db;
        G_h2db; G_h3db; G_h1sd; G_h2sd; G_h3sd; G_h1ss; G_h2ss;
        G_h3ss; G_h1sb; G_h2sb; G_h3sb; G_h1bd; G_h2bd; G_h3bd;
        G_h1bs; G_h2bs; G_h3bs; G_h1bb; G_h2bb; G_h3bb;
        G_hud; G_hus; G_hub; G_hcd; G_hcs; G_hcb; G_htd; G_hts; G_htb;
        G_hdu; G_hdc; G_hdt; G_hsu; G_hsc; G_hst; G_hbu; G_hbc; G_hbt;
        G_he1n1; G_he1n2; G_he1n3; G_he2n1; G_he2n2; G_he2n3; G_he3n1;
        G_he3n2; G_he3n3; G_hn1e1; G_hn1e2; G_hn1e3; G_hn2e1; G_hn2e2;
        G_hn2e3; G_hn3e1; G_hn3e2; G_hn3e3;
        G_h1e1e1; G_h2e1e1; G_h3e1e1; G_h1e1e2; G_h2e1e2; G_h3e1e2;
        G_h1e1e3; G_h2e1e3; G_h3e1e3; G_h1e2e1; G_h2e2e1; G_h3e2e1;
        G_h1e2e2; G_h2e2e2; G_h3e2e2; G_h1e2e3; G_h2e2e3; G_h3e2e3;
        G_h1e3e1; G_h2e3e1; G_h3e3e1; G_h1e3e2; G_h2e3e2; G_h3e3e2;
        G_h1e3e3; G_h2e3e3; G_h3e3e3]

    let add_complex_array_tag c = (Complex_Array c, [Integer 0; Integer 0])

    let derived_parameter_arrays =
      [ nc_coupling G_NC_neutrino half (Integer 0);
        nc_coupling G_NC_lepton (Neg half) (Integer (-1));
        nc_coupling G_NC_up half (Quot (Integer 2, Integer 3));
        nc_coupling G_NC_down (Neg half) (Quot (Integer (-1), Integer 3));
       ] @ (List.map add_complex_array_tag array_list)

    let parameters () =
      { input = []; derived = []; derived_arrays = derived_parameter_arrays}

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

    let electromagnetic_currents n =
      List.map mgm
        [ ((L (-n), Ga, L n), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-n), Ga, U n), FBF (1, Psibar, V, Psi), Q_up);
          ((D (-n), Ga, D n), FBF (1, Psibar, V, Psi), Q_down) ]

    let color_currents n =
      List.map mgm
        [ ((U (-n), Gl, U n), FBF (1, Psibar, V, Psi), Gs);
          ((D (-n), Gl, D n), FBF (1, Psibar, V, Psi), Gs) ]

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

    let charged_currents n =
      List.map mgm
        [ ((L (-n), Wm, N n), FBF (1, Psibar, VL, Psi), G_CC);
          ((N (-n), Wp, L n), FBF (1, Psibar, VL, Psi), G_CC);
          ((D (-n), Wm, U n), FBF (1, Psibar, VL, Psi), G_CC);
          ((U (-n), Wp, D n), FBF (1, Psibar, VL, Psi), G_CC) ] 

    let yukawa =
      [ ((M (U (-1)),O Hh,M (U 1)), FBF (1,Psibar,SP,Psi), G_h1uu);
	((M (U (-1)),O HH,M (U 1)), FBF (1,Psibar,SP,Psi), G_h2uu);
	((M (U (-1)),O HA,M (U 1)), FBF (1,Psibar,SP,Psi), G_h3uu);
	((M (U (-2)),O Hh,M (U 2)), FBF (1,Psibar,SP,Psi), G_h1cc);
	((M (U (-2)),O HH,M (U 2)), FBF (1,Psibar,SP,Psi), G_h2cc);
	((M (U (-2)),O HA,M (U 2)), FBF (1,Psibar,SP,Psi), G_h3cc);
	((M (U (-3)),O Hh,M (U 3)), FBF (1,Psibar,SP,Psi), G_h1tt);
	((M (U (-3)),O HH,M (U 3)), FBF (1,Psibar,SP,Psi), G_h2tt);
	((M (U (-3)),O HA,M (U 3)), FBF (1,Psibar,SP,Psi), G_h3tt); 

	((M (D (-1)),O Hh,M (D 1)), FBF (1,Psibar,SP,Psi), G_h1dd); 
	((M (D (-1)),O HH,M (D 1)), FBF (1,Psibar,SP,Psi), G_h2dd); 
	((M (D (-1)),O HA,M (D 1)), FBF (1,Psibar,SP,Psi), G_h3dd); 
	((M (D (-2)),O Hh,M (D 2)), FBF (1,Psibar,SP,Psi), G_h1ss); 
	((M (D (-2)),O HH,M (D 2)), FBF (1,Psibar,SP,Psi), G_h2ss); 
	((M (D (-2)),O HA,M (D 2)), FBF (1,Psibar,SP,Psi), G_h3ss); 
	((M (D (-3)),O Hh,M (D 3)), FBF (1,Psibar,SP,Psi), G_h1bb); 
	((M (D (-3)),O HH,M (D 3)), FBF (1,Psibar,SP,Psi), G_h2bb); 
	((M (D (-3)),O HA,M (D 3)), FBF (1,Psibar,SP,Psi), G_h3bb); 
				    			   
	((M (U (-1)),O Hp,M (D 1)), FBF (1,Psibar,SP,Psi), G_hud); 
	((M (U (-2)),O Hp,M (D 2)), FBF (1,Psibar,SP,Psi), G_hcs); 
	((M (U (-3)),O Hp,M (D 3)), FBF (1,Psibar,SP,Psi), G_htb);
				    			   
	((M (D (-1)),O Hm,M (U 1)), FBF (1,Psibar,SP,Psi), G_hdu); 
	((M (D (-2)),O Hm,M (U 2)), FBF (1,Psibar,SP,Psi), G_hsc); 
	((M (D (-3)),O Hm,M (U 3)), FBF (1,Psibar,SP,Psi), G_hbt);

        ((M (L (-1)),O Hh,M (L 1)), FBF (1,Psibar,SP,Psi), G_h1e1e1); 
        ((M (L (-1)),O HH,M (L 1)), FBF (1,Psibar,SP,Psi), G_h2e1e1); 
        ((M (L (-1)),O HA,M (L 1)), FBF (1,Psibar,SP,Psi), G_h3e1e1); 
        ((M (L (-2)),O Hh,M (L 2)), FBF (1,Psibar,SP,Psi), G_h1e2e2); 
        ((M (L (-2)),O HH,M (L 2)), FBF (1,Psibar,SP,Psi), G_h2e2e2); 
        ((M (L (-2)),O HA,M (L 2)), FBF (1,Psibar,SP,Psi), G_h3e2e2); 
        ((M (L (-3)),O Hh,M (L 3)), FBF (1,Psibar,SP,Psi), G_h1e3e3); 
        ((M (L (-3)),O HH,M (L 3)), FBF (1,Psibar,SP,Psi), G_h2e3e3); 
        ((M (L (-3)),O HA,M (L 3)), FBF (1,Psibar,SP,Psi), G_h3e3e3) 
				    			   
(*i	((M (N (-1)),O Hp,M (L 1)), FBF (1,Psibar,SR,Psi), G_hn1e1); 
	((M (N (-2)),O Hp,M (L 2)), FBF (1,Psibar,SR,Psi), G_hn2e2); 
	((M (N (-3)),O Hp,M (L 3)), FBF (1,Psibar,SR,Psi), G_hn3e3); 
	((M (L (-1)),O Hm,M (N 1)), FBF (1,Psibar,SL,Psi), G_he1n1); 
	((M (L (-2)),O Hm,M (N 2)), FBF (1,Psibar,SL,Psi), G_he2n2); 
	((M (L (-3)),O Hm,M (N 3)), FBF (1,Psibar,SL,Psi), G_he3n3); i*)
				    			   
 ]	@
       if Flags.ckm_present then
       [((M (U (-1)),O Hh, M (U 2)), FBF (1,Psibar,SP,Psi), G_h1uc);
	((M (U (-1)),O Hh, M (U 3)), FBF (1,Psibar,SP,Psi), G_h1ut);
	((M (U (-2)),O Hh,M (U 1)), FBF (1,Psibar,SP,Psi), G_h1cu);
	((M (U (-2)),O Hh,M (U 3)), FBF (1,Psibar,SP,Psi), G_h1ct);
	((M (U (-1)),O HH,M (U 2)), FBF (1,Psibar,SP,Psi), G_h2uc);
	((M (U (-1)),O HH,M (U 3)), FBF (1,Psibar,SP,Psi), G_h2ut);
	((M (U (-1)),O HA,M (U 2)), FBF (1,Psibar,SP,Psi), G_h3uc);
	((M (U (-1)),O HA,M (U 3)), FBF (1,Psibar,SP,Psi), G_h3ut);
	((M (U (-2)),O HH,M (U 1)), FBF (1,Psibar,SP,Psi), G_h2cu);
	((M (U (-2)),O HH,M (U 3)), FBF (1,Psibar,SP,Psi), G_h2ct);
	((M (U (-2)),O HA,M (U 1)), FBF (1,Psibar,SP,Psi), G_h3cu); 
	((M (U (-2)),O HA,M (U 3)), FBF (1,Psibar,SP,Psi), G_h3ct);
	((M (U (-3)),O Hh,M (U 1)), FBF (1,Psibar,SP,Psi), G_h1tu);
	((M (U (-3)),O Hh,M (U 2)), FBF (1,Psibar,SP,Psi), G_h1tc); 
	((M (U (-3)),O HH,M (U 1)), FBF (1,Psibar,SP,Psi), G_h2tu);
	((M (U (-3)),O HH,M (U 2)), FBF (1,Psibar,SP,Psi), G_h2tc);
	((M (U (-3)),O HA,M (U 1)), FBF (1,Psibar,SP,Psi), G_h3tu);
	((M (U (-3)),O HA,M (U 2)), FBF (1,Psibar,SP,Psi), G_h3tc); 

	((M (D (-1)),O Hh,M (D 2)), FBF (1,Psibar,SP,Psi), G_h1ds); 
	((M (D (-1)),O Hh,M (D 3)), FBF (1,Psibar,SP,Psi), G_h1db); 
	((M (D (-1)),O HH,M (D 2)), FBF (1,Psibar,SP,Psi), G_h2ds); 
	((M (D (-1)),O HH,M (D 3)), FBF (1,Psibar,SP,Psi), G_h2db); 
	((M (D (-1)),O HA,M (D 2)), FBF (1,Psibar,SP,Psi), G_h3ds); 
	((M (D (-1)),O HA,M (D 3)), FBF (1,Psibar,SP,Psi), G_h3db); 
	((M (D (-2)),O Hh,M (D 1)), FBF (1,Psibar,SP,Psi), G_h1sd); 
	((M (D (-2)),O Hh,M (D 3)), FBF (1,Psibar,SP,Psi), G_h1sb); 
	((M (D (-2)),O HH,M (D 1)), FBF (1,Psibar,SP,Psi), G_h2sd); 
	((M (D (-2)),O HH,M (D 3)), FBF (1,Psibar,SP,Psi), G_h2sb); 
	((M (D (-2)),O HA,M (D 1)), FBF (1,Psibar,SP,Psi), G_h3sd); 

(*i	((M (N (-1)),O Hp,M (L 2)), FBF (1,Psibar,SR,Psi), G_hn1e2); 
        ((M (N (-1)),O Hp,M (L 3)), FBF (1,Psibar,SR,Psi), G_hn1e3); 
        ((M (N (-2)),O Hp,M (L 1)), FBF (1,Psibar,SR,Psi), G_hn2e1); 
        ((M (N (-2)),O Hp,M (L 3)), FBF (1,Psibar,SR,Psi), G_hn2e3); 
        ((M (N (-3)),O Hp,M (L 1)), FBF (1,Psibar,SR,Psi), G_hn3e1); 
        ((M (N (-3)),O Hp,M (L 2)), FBF (1,Psibar,SR,Psi), G_hn3e2); 
        ((M (L (-1)),O Hm,M (N 2)), FBF (1,Psibar,SL,Psi), G_he1n2); 
        ((M (L (-1)),O Hm,M (N 3)), FBF (1,Psibar,SL,Psi), G_he1n3); 
        ((M (L (-2)),O Hm,M (N 1)), FBF (1,Psibar,SL,Psi), G_he2n1); 
        ((M (L (-2)),O Hm,M (N 3)), FBF (1,Psibar,SL,Psi), G_he2n3); 
        ((M (L (-3)),O Hm,M (N 1)), FBF (1,Psibar,SL,Psi), G_he3n1); 
        ((M (L (-3)),O Hm,M (N 2)), FBF (1,Psibar,SL,Psi), G_he3n2); i*)

        ((M (L (-1)),O Hh,M (L 2)), FBF (1,Psibar,SP,Psi), G_h1e1e2); 
        ((M (L (-1)),O Hh,M (L 3)), FBF (1,Psibar,SP,Psi), G_h1e1e3); 
        ((M (L (-1)),O HH,M (L 2)), FBF (1,Psibar,SP,Psi), G_h2e1e2); 
        ((M (L (-1)),O HH,M (L 3)), FBF (1,Psibar,SP,Psi), G_h2e1e3); 
        ((M (L (-1)),O HA,M (L 2)), FBF (1,Psibar,SP,Psi), G_h3e1e2); 
        ((M (L (-1)),O HA,M (L 3)), FBF (1,Psibar,SP,Psi), G_h3e1e3); 
        ((M (L (-2)),O Hh,M (L 1)), FBF (1,Psibar,SP,Psi), G_h1e2e1); 
        ((M (L (-2)),O Hh,M (L 3)), FBF (1,Psibar,SP,Psi), G_h1e2e3); 
        ((M (L (-2)),O HH,M (L 1)), FBF (1,Psibar,SP,Psi), G_h2e2e1); 
        ((M (L (-2)),O HH,M (L 3)), FBF (1,Psibar,SP,Psi), G_h2e2e3); 
        ((M (L (-2)),O HA,M (L 1)), FBF (1,Psibar,SP,Psi), G_h3e2e1); 
        ((M (L (-2)),O HA,M (L 3)), FBF (1,Psibar,SP,Psi), G_h3e2e3); 
        ((M (L (-3)),O Hh,M (L 1)), FBF (1,Psibar,SP,Psi), G_h1e3e1); 
        ((M (L (-3)),O Hh,M (L 2)), FBF (1,Psibar,SP,Psi), G_h1e3e2); 
        ((M (L (-3)),O HH,M (L 1)), FBF (1,Psibar,SP,Psi), G_h2e3e1); 
        ((M (L (-3)),O HH,M (L 2)), FBF (1,Psibar,SP,Psi), G_h2e3e2); 
        ((M (L (-3)),O HA,M (L 1)), FBF (1,Psibar,SP,Psi), G_h3e3e1); 
        ((M (L (-3)),O HA,M (L 2)), FBF (1,Psibar,SP,Psi), G_h3e3e2) ]
	else
         [] 

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
          ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_Gs) ]

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
    let quartic_gauge =
      List.map qgc
        [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
          (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
          (Wm, Z, Wp, Ga), minus_gauge4, G_AZWW;
          (Wm, Ga, Wp, Ga), minus_gauge4, G_AAWW;
          (Gl, Gl, Gl, Gl), gauge4, G2]

    let gauge_higgs =
      [ (G Ga, O Hp, O Hm), Vector_Scalar_Scalar 1, G_AHpHm;
        (G Z, O Hp, O Hm), Vector_Scalar_Scalar 1, G_ZHpHm;
        (G Z, O Hh, O HH), Vector_Scalar_Scalar 1, G_Zh1h2;
        (G Z, O Hh, O HA), Vector_Scalar_Scalar 1, G_Zh1h3;
        (G Z, O HH, O HA), Vector_Scalar_Scalar 1, G_Zh2h3;
        (G Wp, O Hm, O Hh), Vector_Scalar_Scalar 1, G_WpHmh1;
        (G Wp, O Hm, O HH), Vector_Scalar_Scalar 1, G_WpHmh2;
        (G Wp, O Hm, O HA), Vector_Scalar_Scalar 1, G_WpHmh3;
        (G Wm, O Hp, O Hh), Vector_Scalar_Scalar 1, G_WmHph1;
        (G Wm, O Hp, O HH), Vector_Scalar_Scalar 1, G_WmHph2;
        (G Wm, O Hp, O HA), Vector_Scalar_Scalar 1, G_WmHph3;
        (O Hh, G Z, G Z), Scalar_Vector_Vector 1, G_h1ZZ;
        (O HH, G Z, G Z), Scalar_Vector_Vector 1, G_h2ZZ;
        (O HA, G Z, G Z), Scalar_Vector_Vector 1, G_h3ZZ;
        (O Hh, G Wp, G Wm), Scalar_Vector_Vector 1, G_h1WpWm;
        (O HH, G Wp, G Wm), Scalar_Vector_Vector 1, G_h2WpWm;
        (O HA, G Wp, G Wm), Scalar_Vector_Vector 1, G_h3WpWm ]

    let gauge_higgs4 =
      [ (O Hh, O Hh, G Wp, G Wm), Scalar2_Vector2 1, G_hhWpWm;
        (O HH, O HH, G Wp, G Wm), Scalar2_Vector2 1, G_hhWpWm;
        (O HA, O HA, G Wp, G Wm), Scalar2_Vector2 1, G_hhWpWm;
        (O Hh, O Hh, G Z, G Z), Scalar2_Vector2 1, G_hhZZ;
        (O HH, O HH, G Z, G Z), Scalar2_Vector2 1, G_hhZZ;
        (O HA, O HA, G Z, G Z), Scalar2_Vector2 1, G_hhZZ;        
        (O Hp, O Hm, G Ga, G Ga), Scalar2_Vector2 1, G_HpHmAA;
        (O Hp, O Hm, G Z, G Z), Scalar2_Vector2 1, G_HpHmZZ;
        (O Hp, O Hm, G Ga, G Z), Scalar2_Vector2 1, G_HpHmAZ;
        (O Hp, O Hm, G Wp, G Wm), Scalar2_Vector2 1, G_HpHmWpWm;        
        (O Hh, O Hp, G Ga, G Wm), Scalar2_Vector2 1, G_h1HpAWm;
        (O HH, O Hp, G Ga, G Wm), Scalar2_Vector2 1, G_h2HpAWm;
        (O HA, O Hp, G Ga, G Wm), Scalar2_Vector2 1, G_h3HpAWm;
        (O Hh, O Hp, G Z, G Wm), Scalar2_Vector2 1, G_h1HpZWm;
        (O HH, O Hp, G Z, G Wm), Scalar2_Vector2 1, G_h2HpZWm;
        (O HA, O Hp, G Z, G Wm), Scalar2_Vector2 1, G_h3HpZWm;
        (O Hh, O Hm, G Ga, G Wp), Scalar2_Vector2 1, G_h1HpAWmC; 
        (O HH, O Hm, G Ga, G Wp), Scalar2_Vector2 1, G_h2HpAWmC;
        (O HA, O Hm, G Ga, G Wp), Scalar2_Vector2 1, G_h3HpAWmC;
        (O Hh, O Hm, G Z, G Wp), Scalar2_Vector2 1, G_h1HpZWmC;
        (O HH, O Hm, G Z, G Wp), Scalar2_Vector2 1, G_h2HpZWmC;
        (O HA, O Hm, G Z, G Wp), Scalar2_Vector2 1, G_h3HpZWmC ]
       
    let higgs =
      [ (O Hh, O Hp, O Hm), Scalar_Scalar_Scalar 1, G_h1HpHm;
	(O HH, O Hp, O Hm), Scalar_Scalar_Scalar 1, G_h2HpHm;
	(O HA, O Hp, O Hm), Scalar_Scalar_Scalar 1, G_h3HpHm;
	(O Hh, O Hh, O Hh), Scalar_Scalar_Scalar 1, G_h111;
	(O Hh, O Hh, O HH), Scalar_Scalar_Scalar 1, G_h112;
	(O Hh, O Hh, O HA), Scalar_Scalar_Scalar 1, G_h113;
	(O HH, O HH, O Hh), Scalar_Scalar_Scalar 1, G_h221;
	(O HH, O HH, O HH), Scalar_Scalar_Scalar 1, G_h222;
	(O HH, O HH, O HA), Scalar_Scalar_Scalar 1, G_h223;
	(O HA, O HA, O Hh), Scalar_Scalar_Scalar 1, G_h331;
	(O HA, O HA, O HH), Scalar_Scalar_Scalar 1, G_h332;
	(O HA, O HA, O HA), Scalar_Scalar_Scalar 1, G_h333;
	(O Hh, O HH, O HA), Scalar_Scalar_Scalar 1, G_h123 ]

    let higgs4 =
      [ (O Hp, O Hm, O Hp, O Hm), Scalar4 1, G_HpHmHpHm;
	(O Hp, O Hm, O Hh, O Hh), Scalar4 1, G_HpHm11;
	(O Hp, O Hm, O Hh, O HH), Scalar4 1, G_HpHm12;
	(O Hp, O Hm, O Hh, O HA), Scalar4 1, G_HpHm13;
	(O Hp, O Hm, O HH, O HH), Scalar4 1, G_HpHm22;
	(O Hp, O Hm, O HH, O HA), Scalar4 1, G_HpHm23;
	(O Hp, O Hm, O HA, O HA), Scalar4 1, G_HpHm33; 
	(O Hh, O Hh, O Hh, O Hh), Scalar4 1, G_h1111;
	(O Hh, O Hh, O Hh, O HH), Scalar4 1, G_h1112;
	(O Hh, O Hh, O Hh, O HA), Scalar4 1, G_h1113;
	(O Hh, O Hh, O HH, O HH), Scalar4 1, G_h1122;
	(O Hh, O Hh, O HH, O HA), Scalar4 1, G_h1123;
	(O Hh, O Hh, O HA, O HA), Scalar4 1, G_h1133;
	(O Hh, O HH, O HH, O HH), Scalar4 1, G_h1222;
	(O Hh, O HH, O HH, O HA), Scalar4 1, G_h1223;
	(O Hh, O HH, O HA, O HA), Scalar4 1, G_h1233;
	(O Hh, O HA, O HA, O HA), Scalar4 1, G_h1333;
	(O HH, O HH, O HH, O HH), Scalar4 1, G_h2222;
	(O HH, O HH, O HH, O HA), Scalar4 1, G_h2223;
	(O HH, O HH, O HA, O HA), Scalar4 1, G_h2233;
	(O HH, O HA, O HA, O HA), Scalar4 1, G_h2333;
	(O HA, O HA, O HA, O HA), Scalar4 1, G_h3333 ]

    let goldstone_vertices =
      [ ((O Phi0, G Wm, G Wp), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phip, G Ga, G Wm), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phip, G Z, G Wm), Scalar_Vector_Vector 1, I_G_ZWW);
        ((O Phim, G Wp, G Ga), Scalar_Vector_Vector 1, I_Q_W);
        ((O Phim, G Wp, G Z), Scalar_Vector_Vector 1, I_G_ZWW) ]

    let vertices3 =
      (ThoList.flatmap electromagnetic_currents [1;2;3] @
       ThoList.flatmap color_currents [1;2;3] @
       ThoList.flatmap neutral_currents [1;2;3] @
       ThoList.flatmap charged_currents [1;2;3] @
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
      | "g" -> G Gl
      | "A" -> G Ga | "Z" | "Z0" -> G Z
      | "W+" -> G Wp | "W-" -> G Wm
      | "h0" -> O Hh
      | "H0" -> O HH
      | "A0" -> O HA
      | "Hp" | "H+" -> O Hp
      | "Hm" | "H-" -> O Hm
      | _ -> invalid_arg "Modellib_BSM.TwoHiggsDoublet.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib_BSM.TwoHiggsDoublet.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib_BSM.TwoHiggsDoublet.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib_BSM.TwoHiggsDoublet.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_BSM.TwoHiggsDoublet.flavor_to_string: invalid down type quark"
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
          | Hh -> "h0" | HH -> "H0" | HA -> "A0"
	  | Hp -> "H+" | Hm -> "H-"
          end

    let flavor_to_TeX = function
      | M f ->
          begin match f with
          | L 1 -> "e^-" | L (-1) -> "e^+"
          | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
          | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
          | L _ -> invalid_arg
                "Modellib_BSM.TwoHiggsDoublet.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib_BSM.TwoHiggsDoublet.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib_BSM.TwoHiggsDoublet.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Modellib_BSM.TwoHiggsDoublet.flavor_to_TeX: invalid down type quark"
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
          | Hh -> "h^0" | HH -> "H^0" | HA -> "A^0"
	  | Hp -> "H^+" | Hm -> "H^-"
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
          | Hh -> "h" | HH -> "h0" | HA -> "a0"
          | Hp -> "hp" | Hm -> "hm"
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
          | Phip -> 27 | Phim -> -27 | Phi0 -> 26
          | Hh -> 25
          | HH -> 35
          | HA -> 36
          | Hp -> 37
          | Hm -> -37
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
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu"
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_CC -> "gcc"
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" | I_G_WWW -> "igwww"
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_AZWW -> "gazww" | G_AAWW -> "gaaww"
      | G_htt -> "ghtt" | G_hbb -> "ghbb" | G_hcc -> "ghcc"
      | G_Htt -> "gh0tt" | G_Hbb -> "gh0bb" | G_Hcc -> "gh0cc" 
      | I_G_Att -> "iga0tt" | I_G_Abb -> "iga0bb" | I_G_Acc -> "iga0cc"
      | G_htautau -> "ghtautau" | G_hmumu -> "ghmumu"
      | G_Htautau -> "gh0tautau" | G_Hmumu -> "gh0mumu"
      | I_G_Atautau -> "iga0tautau" | I_G_Amumu -> "iga0mumu"
      | G_Htb -> "ghptb" | G_Hcs -> "ghpcs"
      | G_Htaunu -> "ghptaunu" | G_Hmunu -> "ghpmunu"
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | G_AHpHm -> "gAHpHm" | G_ZHpHm -> "gZHpHm" 
      | G_Zh1h2 -> "gZh1h2" | G_Zh1h3 -> "gZh1h3" 
      | G_Zh2h3 -> "gZh2h3" | G_WpHmh1 -> "gWpHmh1" 
      | G_WpHmh2 -> "gWpHmh2" | G_WpHmh3 -> "gWpHmh3" 
      | G_WmHph1 -> "gWmHph1" | G_WmHph2 -> "gWmHph2" 
      | G_WmHph3 -> "gWmHph3" | G_h1ZZ -> "gh1ZZ" 
      | G_h2ZZ -> "gh2ZZ" | G_h3ZZ -> "gh3ZZ" 
      | G_h1WpWm -> "gh1WpWm" | G_h2WpWm -> "gh2WpWm" 
      | G_h3WpWm -> "gh3WpWm" | G_hhWpWm -> "ghhWpWm" 
      | G_hhZZ -> "ghhZZ" | G_HpHmAA -> "gHpHmAA" 
      | G_HpHmZZ -> "gHpHmZZ" | G_HpHmAZ -> "gHpHmAZ" 
      | G_HpHmWpWm -> "gHpHmWpWm" | G_h1HpAWm -> "gh1HpAWm" 
      | G_h2HpAWm -> "gh2HpAWm" | G_h3HpAWm -> "gh3HpAWm" 
      | G_h1HpZWm -> "gh1HpZWm" | G_h2HpZWm -> "gh2HpZWm" 
      | G_h3HpZWm -> "gh3HpZWm" | G_h1HpAWmC -> "gh1HpAWmC" 
      | G_h2HpAWmC -> "gh2HpAWmC" | G_h3HpAWmC -> "gh3HpAWmC" 
      | G_h1HpZWmC -> "gh1HpZWmC" | G_h2HpZWmC -> "gh2HpZWmC" 
      | G_h3HpZWmC -> "gh3HpZWmC" | G_h1HpHm -> "gh1HpHm" 
      | G_h2HpHm -> "gh2HpHm" | G_h3HpHm -> "gh3HpHm" 
      | G_h111 -> "gh111" | G_h112 -> "gh112" | G_h113 -> "gh113" 
      | G_h221 -> "gh221" | G_h222 -> "gh222" | G_h223 -> "gh223" 
      | G_h331 -> "gh331" | G_h332 -> "gh332" 
      | G_h333 -> "gh333" | G_h123 -> "gh123" 
      | G_HpHmHpHm -> "gHpHmHpHm" | G_HpHm11 -> "gHpHm11" 
      | G_HpHm12 -> "gHpHm12" | G_HpHm13 -> "gHpHm13" 
      | G_HpHm22 -> "gHpHm22" | G_HpHm23 -> "gHpHm23" 
      | G_HpHm33 -> "gHpHm33" | G_h1111 -> "gh1111" 
      | G_h1112 -> "gh1112" | G_h1113 -> "gh1113" 
      | G_h1122 -> "gh1122" | G_h1123 -> "gh1123" 
      | G_h1133 -> "gh1133" | G_h1222 -> "gh1222" 
      | G_h1223 -> "gh1223" | G_h1233 -> "gh1233" 
      | G_h1333 -> "gh1333" | G_h2222 -> "gh2222" 
      | G_h2223 -> "gh2223" | G_h2233 -> "gh2233" 
      | G_h2333 -> "gh2333" | G_h3333 -> "gh3333" 
      | G_h1uu -> "gh1uu" | G_h2uu -> "gh2uu" 
      | G_h3uu -> "gh3uu" | G_h1uc -> "gh1uc" 
      | G_h2uc -> "gh2uc" | G_h3uc -> "gh3uc" 
      | G_h1ut -> "gh1ut" | G_h2ut -> "gh2ut" 
      | G_h3ut -> "gh3ut" | G_h1cu -> "gh1cu" 
      | G_h2cu -> "gh2cu" | G_h3cu -> "gh3cu" 
      | G_h1cc -> "gh1cc" | G_h2cc -> "gh2cc" 
      | G_h3cc -> "gh3cc" | G_h1ct -> "gh1ct" 
      | G_h2ct -> "gh2ct" | G_h3ct -> "gh3ct" 
      | G_h1tu -> "gh1tu" | G_h2tu -> "gh2tu" 
      | G_h3tu -> "gh3tu" | G_h1tc -> "gh1tc" 
      | G_h2tc -> "gh2tc" | G_h3tc -> "gh3tc" 
      | G_h1tt -> "gh1tt" | G_h2tt -> "gh2tt" 
      | G_h3tt -> "gh3tt"
      | G_h1dd -> "gh1dd" | G_h2dd -> "gh2dd" 
      | G_h3dd -> "gh3dd" | G_h1ds -> "gh1ds" 
      | G_h2ds -> "gh2ds" | G_h3ds -> "gh3ds" 
      | G_h1db -> "gh1db" | G_h2db -> "gh2db" 
      | G_h3db -> "gh3db" | G_h1sd -> "gh1sd" 
      | G_h2sd -> "gh2sd" | G_h3sd -> "gh3sd" 
      | G_h1ss -> "gh1ss" | G_h2ss -> "gh2ss" 
      | G_h3ss -> "gh3ss" | G_h1sb -> "gh1sb" 
      | G_h2sb -> "gh2sb" | G_h3sb -> "gh3sb" 
      | G_h1bd -> "gh1bd" | G_h2bd -> "gh2bd" 
      | G_h3bd -> "gh3bd" | G_h1bs -> "gh1bs" 
      | G_h2bs -> "gh2bs" | G_h3bs -> "gh3bs" 
      | G_h1bb -> "gh1bb" | G_h2bb -> "gh2bb" 
      | G_h3bb -> "gh3bb"
      | G_hud -> "ghud" | G_hus -> "ghus" 
      | G_hub -> "ghub" | G_hcd -> "ghcd" 
      | G_hcs -> "ghcs" | G_hcb -> "ghcb" 
      | G_htd -> "ghtd" | G_hts -> "ghts" 
      | G_htb -> "ghtb" 
      | G_hdu -> "ghdu" | G_hdc -> "ghdc" 
      | G_hdt -> "ghdt" | G_hsu -> "ghsu" 
      | G_hsc -> "ghsc" | G_hst -> "ghst" 
      | G_hbu -> "ghbu" | G_hbc -> "ghbc" 
      | G_hbt -> "ghbt" 
      | G_he1n1 -> "ghe1n1" | G_he1n2 -> "ghe1n2" 
      | G_he1n3 -> "ghe1n3" | G_he2n1 -> "ghe2n1" 
      | G_he2n2 -> "ghe2n2" | G_he2n3 -> "ghe2n3" 
      | G_he3n1 -> "ghe3n1" | G_he3n2 -> "ghe3n2" 
      | G_he3n3 -> "ghe3n3" | G_hn1e1 -> "ghn1e1" 
      | G_hn1e2 -> "ghn1e2" | G_hn1e3 -> "ghn1e3" 
      | G_hn2e1 -> "ghn2e1" | G_hn2e2 -> "ghn2e2" 
      | G_hn2e3 -> "ghn2e3" | G_hn3e1 -> "ghn3e1" 
      | G_hn3e2 -> "ghn3e2" | G_hn3e3 -> "ghn3e3" 
      | G_h1e1e1 -> "gh1e1e1" | G_h2e1e1 -> "gh2e1e1" 
      | G_h3e1e1 -> "gh3e1e1" | G_h1e1e2 -> "gh1e1e2" 
      | G_h2e1e2 -> "gh2e1e2" | G_h3e1e2 -> "gh3e1e2" 
      | G_h1e1e3 -> "gh1e1e3" | G_h2e1e3 -> "gh2e1e3" 
      | G_h3e1e3 -> "gh3e1e3" | G_h1e2e1 -> "gh1e2e1" 
      | G_h2e2e1 -> "gh2e2e1" | G_h3e2e1 -> "gh3e2e1" 
      | G_h1e2e2 -> "gh1e2e2" | G_h2e2e2 -> "gh2e2e2" 
      | G_h3e2e2 -> "gh3e2e2" | G_h1e2e3 -> "gh1e2e3" 
      | G_h2e2e3 -> "gh2e2e3" | G_h3e2e3 -> "gh3e2e3" 
      | G_h1e3e1 -> "gh1e3e1" | G_h2e3e1 -> "gh2e3e1" 
      | G_h3e3e1 -> "gh3e3e1" | G_h1e3e2 -> "gh1e3e2" 
      | G_h2e3e2 -> "gh2e3e2" | G_h3e3e2 -> "gh3e3e2" 
      | G_h1e3e3 -> "gh1e3e3" | G_h2e3e3 -> "gh2e3e3" 
      | G_h3e3e3 -> "gh3e3e3" 
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f

  end

module type SSC_flags =
  sig
    val higgs_triangle : bool (* $H\gamma\gamma$, $Hg\gamma$ and $Hgg$ couplings *)
    val higgs_hmm : bool    
    val triple_anom : bool
    val quartic_anom : bool
    val higgs_anom : bool
    val k_matrix : bool
    val k_matrix_tm : bool      
    val ckm_present : bool
    val top_anom : bool
    val top_anom_4f : bool
    val cf_arbitrary : bool
    val higgs_matrix : bool
  end


module SSC_kmatrix: SSC_flags =
  struct
    let higgs_triangle = false
    let higgs_hmm = false
    let triple_anom = false
    let quartic_anom = true
    let higgs_anom = false
    let k_matrix = true
    let k_matrix_tm = false
    let ckm_present = false
    let top_anom = false
    let top_anom_4f = false
    let cf_arbitrary = false
    let higgs_matrix = false
  end

module SSC_kmatrix_2: SSC_flags =
  struct
    let higgs_triangle = false
    let higgs_hmm = false
    let triple_anom = false
    let quartic_anom = true
    let higgs_anom = false
    let k_matrix = true
    let k_matrix_tm = true      
    let ckm_present = false
    let top_anom = false
    let top_anom_4f = false
    let cf_arbitrary = true
    let higgs_matrix = true
  end


(* \thocwmodulesection{Complete Minimal Standard Model including additional Resonances} *)

module SSC (Flags : SSC_flags) =
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
                     QGUG | QBUB | QW | DL | DR

    type matter_field = L of int | N of int | U of int | D of int
    type gauge_boson = Ga | Wp | Wm | Z | Gl
    type other = Phip | Phim | Phi0 | H
                 | Rsigma | Rphin | Rphisn | Rphip | Rphim | Rphipp | Rphimm 
                 | Rf | Rtn | Rtsn | Rtp | Rtm | Rtpp | Rtmm 
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
      failwith "Modellib_BSM.SSC.gauge_symbol: internal error"

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
        "Higgs", List.map other [H];
	"Scalar Resonances", List.map other [Rsigma; Rphin; Rphisn; Rphip; Rphim; Rphipp; Rphimm];
	"Tensor Resonances", List.map other [Rf; Rtn; Rtsn; Rtp; Rtm; Rtpp; Rtmm];
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
          | Rf | Rtn | Rtsn | Rtp | Rtm | Rtpp | Rtmm -> Tensor_2
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
          | H | Rsigma -> Prop_Scalar
          | Rphin | Rphisn | Rphip | Rphim | Rphipp | Rphimm -> Prop_Scalar
          | Rf -> Prop_Tensor_2
          | Rtn | Rtsn | Rtp | Rtm | Rtpp | Rtmm -> Prop_Tensor_2
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
          | H -> H |  Rsigma -> Rsigma
          | Rphin -> Rphin | Rphisn-> Rphisn | Rphip -> Rphim | Rphim -> Rphip
          | Rphipp -> Rphimm | Rphimm -> Rphipp
          | Rf -> Rf
          | Rtn -> Rtn | Rtsn -> Rtsn | Rtp -> Rtm | Rtm -> Rtp
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
          | H | Rsigma | Phi0 | Rphin | Rphisn 
          | Rf | Rtn | Rtsn ->  0//1
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
      | G_SP_ttH
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
      | FS0_HHWW | FS0_HHZZ
      | FS1_HHWW | FS1_HHZZ
      | FM0_HHWW | FM0_HHZZ 
      | FM1_HHWW | FM1_HHZZ
      | FM7_HHWW | FM7_HHZZ
      | Alpha_WWWW0 | Alpha_ZZWW1 | Alpha_WWWW2
      | Alpha_ZZWW0 | Alpha_ZZZZ
      | FT0_WWWW0 | FT0_WWWW2
      | FT0_ZZWW0 | FT0_ZZWW1
      | FT0_ZZZZ  | FT0_AAAA
      | FT0_AAWW0 | FT0_AAWW1
      | FT0_AAZZ  
      | FT0_AZWW0 | FT0_AZWW1
      | FT0_AAAZ  | FT0_AZZZ
      | FT1_WWWW0 | FT1_WWWW2
      | FT1_ZZWW0 | FT1_ZZWW1
      | FT1_ZZZZ  | FT1_AAAA 
      | FT1_AAWW0 | FT1_AAWW1
      | FT1_AAZZ  
      | FT1_AZWW0 | FT1_AZWW1
      | FT1_AAAZ  | FT1_AZZZ
      | FT2_WWWW0 | FT2_WWWW2
      | FT2_ZZWW0 | FT2_ZZWW1
      | FT2_ZZZZ  | FT2_AAAA
      | FT2_AAWW0 | FT2_AAWW1
      | FT2_AAZZ  
      | FT2_AZWW0 | FT2_AZWW1
      | FT2_AAAZ  | FT2_AZZZ
      | FM0_WWWW0 | FM0_WWWW2
      | FM0_ZZWW0 | FM0_ZZWW1
      | FM0_ZZZZ  
      | FM1_WWWW0 | FM1_WWWW2
      | FM1_ZZWW0 | FM1_ZZWW1
      | FM1_ZZZZ 
      | FM7_WWWW0 | FM7_WWWW2
      | FM7_ZZWW0 | FM7_ZZWW1
      | FM7_ZZZZ
      | D_Alpha_ZZWW0_S | D_Alpha_ZZWW0_T | D_Alpha_ZZWW1_S
      | D_Alpha_ZZWW1_T | D_Alpha_ZZWW1_U | D_Alpha_WWWW0_S
      | D_Alpha_WWWW0_T | D_Alpha_WWWW0_U | D_Alpha_WWWW2_S
      | D_Alpha_WWWW2_T | D_Alpha_ZZZZ_S | D_Alpha_ZZZZ_T
      | D_FT0_ZZWW0_S | D_FT0_ZZWW0_T | D_FT0_ZZWW0_U 
      | D_FT0_ZZWW1_S | D_FT0_ZZWW1_T | D_FT0_ZZWW1_U 
      | D_FT0_WWWW0_S | D_FT0_WWWW0_T | D_FT0_WWWW0_U 
      | D_FT0_WWWW2_S | D_FT0_WWWW2_T | D_FT0_WWWW2_U    
      | D_FT0_ZZZZ_S  | D_FT0_ZZZZ_T  | D_FT0_ZZZZ_U
      | D_FT0_AAAA_S  | D_FT0_AAAA_T  | D_FT0_AAAA_U
      | D_FT0_AAWW0_S | D_FT0_AAWW0_T | D_FT0_AAWW0_U 
      | D_FT0_AAWW1_S | D_FT0_AAWW1_T | D_FT0_AAWW1_U
      | D_FT0_AAZZ_S  | D_FT0_AAZZ_T  | D_FT0_AAZZ_U  
      | D_FT0_AZWW0_S | D_FT0_AZWW0_T | D_FT0_AZWW0_U
      | D_FT0_AZWW1_S | D_FT0_AZWW1_T | D_FT0_AZWW1_U
      | D_FT0_AAAZ_S  | D_FT0_AAAZ_T  | D_FT0_AAAZ_U
      | D_FT0_AZZZ_S  | D_FT0_AZZZ_T  | D_FT0_AZZZ_U     
      | D_FT1_ZZWW0_S | D_FT1_ZZWW0_T | D_FT1_ZZWW0_U 
      | D_FT1_ZZWW1_S | D_FT1_ZZWW1_T | D_FT1_ZZWW1_U 
      | D_FT1_WWWW0_S | D_FT1_WWWW0_T | D_FT1_WWWW0_U 
      | D_FT1_WWWW2_S | D_FT1_WWWW2_T | D_FT1_WWWW2_U    
      | D_FT1_ZZZZ_S  | D_FT1_ZZZZ_T  | D_FT1_ZZZZ_U 
      | D_FT1_AAAA_S  | D_FT1_AAAA_T  | D_FT1_AAAA_U
      | D_FT1_AAWW0_S | D_FT1_AAWW0_T | D_FT1_AAWW0_U 
      | D_FT1_AAWW1_S | D_FT1_AAWW1_T | D_FT1_AAWW1_U
      | D_FT1_AAZZ_S  | D_FT1_AAZZ_T  | D_FT1_AAZZ_U  
      | D_FT1_AZWW0_S | D_FT1_AZWW0_T | D_FT1_AZWW0_U
      | D_FT1_AZWW1_S | D_FT1_AZWW1_T | D_FT1_AZWW1_U   
      | D_FT1_AAAZ_S  | D_FT1_AAAZ_T  | D_FT1_AAAZ_U
      | D_FT1_AZZZ_S  | D_FT1_AZZZ_T  | D_FT1_AZZZ_U      
      | D_FT2_ZZWW0_S | D_FT2_ZZWW0_T | D_FT2_ZZWW0_U 
      | D_FT2_ZZWW1_S | D_FT2_ZZWW1_T | D_FT2_ZZWW1_U 
      | D_FT2_WWWW0_S | D_FT2_WWWW0_T | D_FT2_WWWW0_U 
      | D_FT2_WWWW2_S | D_FT2_WWWW2_T | D_FT2_WWWW2_U    
      | D_FT2_ZZZZ_S  | D_FT2_ZZZZ_T  | D_FT2_ZZZZ_U 
      | D_FT2_AAAA_S  | D_FT2_AAAA_T  | D_FT2_AAAA_U
      | D_FT2_AAWW0_S | D_FT2_AAWW0_T | D_FT2_AAWW0_U 
      | D_FT2_AAWW1_S | D_FT2_AAWW1_T | D_FT2_AAWW1_U
      | D_FT2_AAZZ_S  | D_FT2_AAZZ_T  | D_FT2_AAZZ_U 
      | D_FT2_AZWW0_S | D_FT2_AZWW0_T | D_FT2_AZWW0_U
      | D_FT2_AZWW1_S | D_FT2_AZWW1_T | D_FT2_AZWW1_U 
      | D_FT2_AAAZ_S  | D_FT2_AAAZ_T  | D_FT2_AAAZ_U
      | D_FT2_AZZZ_S  | D_FT2_AZZZ_T  | D_FT2_AZZZ_U 
      | D_FTrsi_ZZWW0_S | D_FTrsi_ZZWW0_T | D_FTrsi_ZZWW0_U 
      | D_FTrsi_ZZWW1_S | D_FTrsi_ZZWW1_T | D_FTrsi_ZZWW1_U 
      | D_FTrsi_WWWW0_S | D_FTrsi_WWWW0_T | D_FTrsi_WWWW0_U 
      | D_FTrsi_WWWW2_S | D_FTrsi_WWWW2_T | D_FTrsi_WWWW2_U    
      | D_FTrsi_ZZZZ_S  | D_FTrsi_ZZZZ_T  | D_FTrsi_ZZZZ_U 
      | D_FTrsi_AAAA_S  | D_FTrsi_AAAA_T  | D_FTrsi_AAAA_U
      | D_FTrsi_AAWW0_S | D_FTrsi_AAWW0_T | D_FTrsi_AAWW0_U 
      | D_FTrsi_AAWW1_S | D_FTrsi_AAWW1_T | D_FTrsi_AAWW1_U
      | D_FTrsi_AAZZ_S  | D_FTrsi_AAZZ_T  | D_FTrsi_AAZZ_U 
      | D_FTrsi_AZWW0_S | D_FTrsi_AZWW0_T | D_FTrsi_AZWW0_U
      | D_FTrsi_AZWW1_S | D_FTrsi_AZWW1_T | D_FTrsi_AZWW1_U 
      | D_FTrsi_AAAZ_S  | D_FTrsi_AAAZ_T  | D_FTrsi_AAAZ_U
      | D_FTrsi_AZZZ_S  | D_FTrsi_AZZZ_T  | D_FTrsi_AZZZ_U       
      | D_FM0_ZZWW0_S | D_FM0_ZZWW0_T | D_FM0_ZZWW0_U 
      | D_FM0_ZZWW1_S | D_FM0_ZZWW1_T | D_FM0_ZZWW1_U 
      | D_FM0_WWWW0_S | D_FM0_WWWW0_T | D_FM0_WWWW0_U 
      | D_FM0_WWWW2_S | D_FM0_WWWW2_T | D_FM0_WWWW2_U    
      | D_FM0_ZZZZ_S  | D_FM0_ZZZZ_T  | D_FM0_ZZZZ_U       
      | D_FM1_ZZWW0_S | D_FM1_ZZWW0_T | D_FM1_ZZWW0_U 
      | D_FM1_ZZWW1_S | D_FM1_ZZWW1_T | D_FM1_ZZWW1_U 
      | D_FM1_WWWW0_S | D_FM1_WWWW0_T | D_FM1_WWWW0_U 
      | D_FM1_WWWW2_S | D_FM1_WWWW2_T | D_FM1_WWWW2_U    
      | D_FM1_ZZZZ_S  | D_FM1_ZZZZ_T  | D_FM1_ZZZZ_U 
      | D_FM7_ZZWW0_S | D_FM7_ZZWW0_T | D_FM7_ZZWW0_U 
      | D_FM7_ZZWW1_S | D_FM7_ZZWW1_T | D_FM7_ZZWW1_U 
      | D_FM7_WWWW0_S | D_FM7_WWWW0_T | D_FM7_WWWW0_U 
      | D_FM7_WWWW2_S | D_FM7_WWWW2_T | D_FM7_WWWW2_U    
      | D_FM7_ZZZZ_S  | D_FM7_ZZZZ_T  | D_FM7_ZZZZ_U
      | D_Alpha_HHHH_S  | D_Alpha_HHHH_T 
      | D_Alpha_HHZZ0_S | D_Alpha_HHWW0_S 
      | D_Alpha_HHZZ0_T | D_Alpha_HHWW0_T
      | D_Alpha_HHZZ1_S | D_Alpha_HHWW1_S 
      | D_Alpha_HHZZ1_T | D_Alpha_HHWW1_T
      | D_Alpha_HHZZ1_U | D_Alpha_HHWW1_U
      | D_FM0_HHZZ0_S | D_FM0_HHWW0_S 
      | D_FM0_HHZZ0_T | D_FM0_HHWW0_T
      | D_FM0_HHZZ0_U | D_FM0_HHWW0_U
      | D_FM0_HHZZ1_S | D_FM0_HHWW1_S 
      | D_FM0_HHZZ1_T | D_FM0_HHWW1_T
      | D_FM0_HHZZ1_U | D_FM0_HHWW1_U 
      | D_FM1_HHZZ0_S | D_FM1_HHWW0_S 
      | D_FM1_HHZZ0_T | D_FM1_HHWW0_T
      | D_FM1_HHZZ0_U | D_FM1_HHWW0_U
      | D_FM1_HHZZ1_S | D_FM1_HHWW1_S 
      | D_FM1_HHZZ1_T | D_FM1_HHWW1_T
      | D_FM1_HHZZ1_U | D_FM1_HHWW1_U 
      | D_FM7_HHZZ0_S | D_FM7_HHWW0_S 
      | D_FM7_HHZZ0_T | D_FM7_HHWW0_T
      | D_FM7_HHZZ0_U | D_FM7_HHWW0_U
      | D_FM7_HHZZ1_S | D_FM7_HHWW1_S 
      | D_FM7_HHZZ1_T | D_FM7_HHWW1_T
      | D_FM7_HHZZ1_U | D_FM7_HHWW1_U
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ
      | G_SWW | G_SWW_T | G_SSWW | G_SZZ 
      | G_SZZ_T | G_SSZZ | G_SHH
      | G_SAA_T | G_SAZ_T 
      | G_PNWW | G_PNZZ | G_PWZ | G_PWW
      | G_PSNWW | G_PSNZZ | G_PSNHH
      | G_FWW | G_FZZ | G_FWW_CF | G_FZZ_CF 
      | G_FWW_T | G_FZZ_T | G_FHH | G_FHH_CF
      | G_TNWW | G_TNZZ | G_TSNWW | G_TSNZZ | G_TWZ | G_TWW
      | G_TNWW_CF | G_TNZZ_CF | G_TSNWW_CF | G_TSNZZ_CF
      | G_TWZ_CF | G_TWW_CF
      | G_Htt | G_Hbb | G_Hcc | G_Hmm | G_Htautau | G_H3 | G_H4
      | FS_H4
      | G_HGaZ | G_HGaGa | G_Hgg
      | G_HGaZ_anom | G_HGaGa_anom | G_HZZ_anom | G_HWW_anom  
      | G_HGaZ_u | G_HZZ_u | G_HWW_u
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

(* \begin{dubious}
   \ldots{} to be continued \ldots{}
   The quartic couplings can't be correct, because the dimensions are wrong!
   \begin{subequations}
     \begin{align}
                  g_{HWW} &= g m_W = 2 \frac{m_W^2}{v}\\
                 g_{HHWW} &= 2 \frac{m_W^2}{v^2} = \frac{g^2}{2} \\
                  g_{HZZ} &= \frac{g}{\cos\theta_w}m_Z \\
                 g_{HHZZ} &= 2 \frac{m_Z^2}{v^2} = \frac{g^2}{2\cos\theta_w} \\
                  g_{Htt} &= \lambda_t \\
                  g_{Hbb} &= \lambda_b=\frac{m_b}{m_t}\lambda_t \\
                  g_{H^3} &= - \frac{3g}{2}\frac{m_H^2}{m_W} = - 3 \frac{m_H^2}{v} 
                  g_{H^4} &= - \frac{3g^2}{4} \frac{m_W^2}{v^2} = -3 \frac{m_H^2}{v^2}  
     \end{align}
   \end{subequations}
   \end{dubious} *)

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

    let yukawa =
      [ ((M (U (-3)), O H, M (U 3)), FBF (1, Psibar, S, Psi), G_Htt);
        ((M (D (-3)), O H, M (D 3)), FBF (1, Psibar, S, Psi), G_Hbb);
        ((M (U (-2)), O H, M (U 2)), FBF (1, Psibar, S, Psi), G_Hcc);
        ((M (L (-3)), O H, M (L 3)), FBF (1, Psibar, S, Psi), G_Htautau) ] @
      if Flags.higgs_hmm then
      [ ((M (L (-2)), O H, M (L 2)), FBF (1, Psibar, S, Psi), G_Hmm)]
          else
      []

      
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
            ((Z, Z, Z, Z),
             Vector4 [(1, C_12_34); (1, C_13_42); (1, C_14_23)], Alpha_ZZZZ);
            ((Wm, Wp, Z, Z),
             Vector4 [1, C_12_34], Alpha_ZZWW0);
            ((Wm, Wp, Z, Z),
             Vector4 [(1, C_13_42); (1, C_14_23)], Alpha_ZZWW1)]
	    @
	  (if Flags.k_matrix_tm then
	      List.map qgc
	   [((Wm, Wm, Wp, Wp),		  
             Dim8_Vector4_t_0 [1, C_13_42], FT0_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_WWWW2); 
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_WWWW2); 
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_WWWW2); 
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_0 [1, C_13_42], FM0_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_0 [1, C_14_23], FM0_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_0 [1, C_12_34], FM0_WWWW2); 
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_1 [1, C_13_42], FM1_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_1 [1, C_14_23], FM1_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_1 [1, C_12_34], FM1_WWWW2);  
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_7 [1, C_13_42], FM7_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_7 [1, C_14_23], FM7_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_7 [1, C_12_34], FM7_WWWW2);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_ZZWW1);      
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_ZZWW1);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_0 [1, C_12_34], FM0_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_0 [1, C_13_42], FM0_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_0 [1, C_14_23], FM0_ZZWW1);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_1 [1, C_12_34], FM1_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_1 [1, C_13_42], FM1_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_1 [1, C_14_23], FM1_ZZWW1); 
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_7 [1, C_12_34], FM7_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_7 [1, C_13_42], FM7_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_7 [1, C_14_23], FM7_ZZWW1);
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_ZZZZ);
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_ZZZZ);             
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_0 [1, C_12_34], FM0_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_0 [1, C_13_42], FM0_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_0 [1, C_14_23], FM0_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_1 [1, C_12_34], FM1_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_1 [1, C_13_42], FM1_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_1 [1, C_14_23], FM1_ZZZZ);
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_7 [1, C_12_34], FM7_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_7 [1, C_13_42], FM7_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_7 [1, C_14_23], FM7_ZZZZ);
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AAAA);  
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AAAA); 
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AAAA);
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AAAA);  
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AAAA); 
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AAAA);
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AAAA);  
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AAAA); 
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AAAA);
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AAWW0);   
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AAWW1);
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AAWW1); 
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AAWW0);   
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AAWW1);
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AAWW1);              
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AAWW0);   
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AAWW1);
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AAWW1);
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AAZZ);  
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AAZZ); 
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AAZZ);
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AAZZ);  
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AAZZ); 
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AAZZ);             
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AAZZ);  
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AAZZ); 
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AAZZ);
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AZWW0);  
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AZWW1); 
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AZWW1);
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AZWW0);  
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AZWW1); 
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AZWW1);             
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AZWW0);  
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AZWW1); 
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AZWW1);
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AAAZ);  
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AAAZ); 
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AAAZ);
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AAAZ);  
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AAAZ); 
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AAAZ);             
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AAAZ);  
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AAAZ); 
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AAAZ);
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AZZZ);  
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AZZZ); 
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AZZZ);
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AZZZ);  
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AZZZ); 
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AZZZ);             
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AZZZ);  
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AZZZ); 
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AZZZ)]
      else
	[] )
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
        
    let k_matrix_quartic_gauge_t_0 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_12_34)]), D_FT0_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_13_42)]), D_FT0_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_14_23)]), D_FT0_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_ZZZZ_U);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAAA_T); 
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAAA_U);        
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AAAA_T);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AAAA_U);                   
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAWW0_S);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAWW0_T);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAWW0_U);                   
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAWW1_S);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAWW1_T);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AAWW1_S);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AAWW1_T); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_12_34)]), D_FT0_AAWW1_S); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_13_42)]), D_FT0_AAWW1_U);
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_14_23)]), D_FT0_AAWW1_T);
            ((Ga, Z, Ga, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAZZ_S);
            ((Ga, Z, Ga, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAZZ_T);
            ((Ga, Z, Ga, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAZZ_U);
            ((Z, Ga, Z, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAZZ_S); 
            ((Z, Ga, Z, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAZZ_T);
            ((Z, Ga, Z, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAZZ_U); 
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAZZ_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAZZ_U); 
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AZWW0_S);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AZWW0_T);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AZWW0_U); 
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AZWW1_S);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AZWW1_T);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AZWW1_U); 
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AZWW1_S);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AZWW1_T);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AZWW1_U); 
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AZWW1_S);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AZWW1_T);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AZWW1_U); 
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AZWW1_S);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AZWW1_T);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AZWW1_U); 
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AAAZ_S);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AAAZ_T);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AAAZ_U); 
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AAAZ_S);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AAAZ_T);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AAAZ_U); 
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AAAZ_S);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AAAZ_T);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AAAZ_U); 
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AZZZ_S);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AZZZ_T);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AZZZ_U); 
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AZZZ_S);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AZZZ_T);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AZZZ_U); 
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AZZZ_S);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AZZZ_T);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AZZZ_U)]        
      else
        []
        
    let k_matrix_quartic_gauge_t_1 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_12_34)]), D_FT1_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_13_42)]), D_FT1_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_14_23)]), D_FT1_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_ZZZZ_U);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AAAA_T); 
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AAAA_U);        
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAAA_T);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAAA_U);                   
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_AAWW0_S);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AAWW0_T);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AAWW0_U);                   
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_AAWW1_S);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AAWW1_T);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AAWW1_S);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AAWW1_T); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_12_34)]), D_FT1_AAWW1_S); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_13_42)]), D_FT1_AAWW1_U);
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_14_23)]), D_FT1_AAWW1_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AAZZ_T); 
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AAZZ_U);        
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAZZ_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAZZ_U); 
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_AZWW0_S);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AZWW0_T);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AZWW0_U); 
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AZWW1_S);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AZWW1_T);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AZWW1_U); 
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AZWW1_S);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AZWW1_T);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AZWW1_U); 
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AZWW1_S);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AZWW1_T);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AZWW1_U); 
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AZWW1_S);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AZWW1_T);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AZWW1_U); 
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAAZ_S);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAAZ_T);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAAZ_U); 
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAAZ_S);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAAZ_T);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAAZ_U); 
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAAZ_S);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAAZ_T);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAAZ_U); 
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AZZZ_S);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AZZZ_T);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AZZZ_U); 
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AZZZ_S);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AZZZ_T);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AZZZ_U); 
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AZZZ_S);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AZZZ_T);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AZZZ_U)]        
      else
        []        

    let k_matrix_quartic_gauge_t_2 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_12_34)]), D_FT2_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_13_42)]), D_FT2_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_14_23)]), D_FT2_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_ZZZZ_U);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AAAA_T); 
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AAAA_U);        
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAAA_T);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAAA_U);                   
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AAWW0_S);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AAWW0_T);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AAWW0_U);                   
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AAWW1_S);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AAWW1_T);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AAWW1_S);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AAWW1_T); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_12_34)]), D_FT2_AAWW1_S); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_13_42)]), D_FT2_AAWW1_U);
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_14_23)]), D_FT2_AAWW1_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AAZZ_T); 
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AAZZ_U);        
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAZZ_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAZZ_U); 
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AZWW0_S);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AZWW0_T);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AZWW0_U); 
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AZWW1_S);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AZWW1_T);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AZWW1_U); 
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AZWW1_S);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AZWW1_T);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AZWW1_U); 
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AZWW1_S);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AZWW1_T);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AZWW1_U); 
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AZWW1_S);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AZWW1_T);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AZWW1_U); 
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAAZ_S);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAAZ_T);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAAZ_U); 
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAAZ_S);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAAZ_T);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAAZ_U); 
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAAZ_S);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAAZ_T);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAAZ_U); 
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AZZZ_S);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AZZZ_T);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AZZZ_U); 
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AZZZ_S);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AZZZ_T);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AZZZ_U); 
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AZZZ_S);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AZZZ_T);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AZZZ_U)]        
      else
        []
        
    let k_matrix_quartic_gauge_t_rsi =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_12_34)]), D_FTrsi_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_13_42)]), D_FTrsi_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_14_23)]), D_FTrsi_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_ZZZZ_U);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAAA_T); 
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAAA_U);        
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AAAA_T);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AAAA_U);                   
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAWW0_S);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAWW0_T);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAWW0_U);                   
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAWW1_S);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAWW1_T);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AAWW1_S);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AAWW1_T); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_12_34)]), D_FTrsi_AAWW1_S); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_13_42)]), D_FTrsi_AAWW1_U);
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_14_23)]), D_FTrsi_AAWW1_T);
            ((Ga, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAZZ_S);
            ((Ga, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAZZ_T);
            ((Ga, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAZZ_U);
            ((Z, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAZZ_S); 
            ((Z, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAZZ_T);
            ((Z, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAZZ_U); 
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAZZ_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAZZ_U); 
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AZWW0_S);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AZWW0_T);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AZWW0_U); 
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AZWW1_S);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AZWW1_T);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AZWW1_U); 
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AZWW1_S);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AZWW1_T);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AZWW1_U); 
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AZWW1_S);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AZWW1_T);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AZWW1_U); 
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AZWW1_S);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AZWW1_T);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AZWW1_U); 
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AAAZ_S);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AAAZ_T);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AAAZ_U); 
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AAAZ_S);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AAAZ_T);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AAAZ_U); 
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AAAZ_S);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AAAZ_T);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AAAZ_U); 
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AZZZ_S);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AZZZ_T);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AZZZ_U); 
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AZZZ_S);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AZZZ_T);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AZZZ_U); 
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AZZZ_S);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AZZZ_T);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AZZZ_U)]        
      else
        []        
        
    let k_matrix_quartic_gauge_m_0 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m0 (1,
                   [(1, C_12_34)]), D_FM0_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m0 (1,
                   [(1, C_13_42)]), D_FM0_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m0 (1,
                   [(1, C_14_23)]), D_FM0_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m0 (2,
                   [(1, C_12_34)]), D_FM0_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m0 (2,
                   [(1, C_13_42)]), D_FM0_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m0 (2,
                   [(1, C_14_23)]), D_FM0_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (3,
                   [(1, C_14_23)]), D_FM0_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (3,
                   [(1, C_13_42)]), D_FM0_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (3,
                   [(1, C_12_34)]), D_FM0_ZZZZ_U)]        
      else
        []

    let k_matrix_quartic_gauge_m_1 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m1 (1,
                   [(1, C_12_34)]), D_FM1_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m1 (1,
                   [(1, C_13_42)]), D_FM1_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m1 (1,
                   [(1, C_14_23)]), D_FM1_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m1 (2,
                   [(1, C_12_34)]), D_FM1_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m1 (2,
                   [(1, C_13_42)]), D_FM1_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m1 (2,
                   [(1, C_14_23)]), D_FM1_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (3,
                   [(1, C_14_23)]), D_FM1_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (3,
                   [(1, C_13_42)]), D_FM1_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (3,
                   [(1, C_12_34)]), D_FM1_ZZZZ_U)]        
      else
        []
        
    let k_matrix_quartic_gauge_m_7 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m7 (1,
                   [(1, C_12_34)]), D_FM7_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m7 (1,
                   [(1, C_13_42)]), D_FM7_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m7 (1,
                   [(1, C_14_23)]), D_FM7_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m7 (2,
                   [(1, C_12_34)]), D_FM7_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m7 (2,
                   [(1, C_13_42)]), D_FM7_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m7 (2,
                   [(1, C_14_23)]), D_FM7_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (3,
                   [(1, C_14_23)]), D_FM7_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (3,
                   [(1, C_13_42)]), D_FM7_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (3,
                   [(1, C_12_34)]), D_FM7_ZZZZ_U)]        
      else
        []    

    let k_matrix_2scalar_2gauge =
      if Flags.k_matrix_tm then
        if Flags.higgs_matrix then
            [ ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms 
                   (0,  [(1, C_12_34)]), D_Alpha_HHZZ0_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (0,  [(1, C_13_42)]), D_Alpha_HHZZ0_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (0,  [(1, C_14_23)]), D_Alpha_HHZZ0_T);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms 
                   (3,  [(1, C_14_23)]), D_Alpha_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (3,  [(1, C_13_42)]), D_Alpha_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (3,  [(1, C_12_34)]), D_Alpha_HHZZ1_U);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms 
                   (6,  [(1, C_13_42)]), D_Alpha_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (6,  [(1, C_12_34)]), D_Alpha_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (6,  [(1, C_14_23)]), D_Alpha_HHZZ1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms 
                   (0,  [(1, C_12_34)]), D_Alpha_HHWW0_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (2,  [(1, C_13_42)]), D_Alpha_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (1,  [(1, C_14_23)]), D_Alpha_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (1,  [(1, C_13_42)]), D_Alpha_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (2,  [(1, C_14_23)]), D_Alpha_HHWW0_T);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms 
                   (3,  [(1, C_14_23)]), D_Alpha_HHWW1_S);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms 
                   (6,  [(1, C_13_42)]), D_Alpha_HHWW1_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (4,  [(1, C_13_42)]), D_Alpha_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (5,  [(1, C_12_34)]), D_Alpha_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (8,  [(1, C_14_23)]), D_Alpha_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (7,  [(1, C_12_34)]), D_Alpha_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (5,  [(1, C_13_42)]), D_Alpha_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (4,  [(1, C_12_34)]), D_Alpha_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (7,  [(1, C_14_23)]), D_Alpha_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (8,  [(1, C_12_34)]), D_Alpha_HHWW1_U) ]
        else
            []
      else
          []
          
    let k_matrix_2scalar_2gauge_m =
      if Flags.k_matrix_tm then
        if Flags.higgs_matrix then
            [ ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM0_HHZZ0_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (0,  [(1, C_13_42)]), D_FM0_HHZZ0_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (0,  [(1, C_14_23)]), D_FM0_HHZZ0_U);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf 
                   (3,  [(1, C_14_23)]), D_FM0_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (3,  [(1, C_13_42)]), D_FM0_HHZZ1_U);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (3,  [(1, C_12_34)]), D_FM0_HHZZ1_T);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM0_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (6,  [(1, C_12_34)]), D_FM0_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (6,  [(1, C_14_23)]), D_FM0_HHZZ1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM0_HHWW0_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (2,  [(1, C_13_42)]), D_FM0_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (1,  [(1, C_14_23)]), D_FM0_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (1,  [(1, C_13_42)]), D_FM0_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (2,  [(1, C_14_23)]), D_FM0_HHWW0_T);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (3,  [(1, C_14_23)]), D_FM0_HHWW1_S);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM0_HHWW1_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (4,  [(1, C_13_42)]), D_FM0_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (5,  [(1, C_12_34)]), D_FM0_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (8,  [(1, C_14_23)]), D_FM0_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (7,  [(1, C_12_34)]), D_FM0_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (5,  [(1, C_13_42)]), D_FM0_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (4,  [(1, C_12_34)]), D_FM0_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (7,  [(1, C_14_23)]), D_FM0_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (8,  [(1, C_12_34)]), D_FM0_HHWW1_T) ]
        else
            []
      else
          [] 
          
    let k_matrix_2scalar_2gauge_m_1 =
      if Flags.k_matrix_tm then
        if Flags.higgs_matrix then
            [ ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM1_HHZZ0_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (0,  [(1, C_13_42)]), D_FM1_HHZZ0_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (0,  [(1, C_14_23)]), D_FM1_HHZZ0_U);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf 
                   (3,  [(1, C_14_23)]), D_FM1_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (3,  [(1, C_13_42)]), D_FM1_HHZZ1_U);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (3,  [(1, C_12_34)]), D_FM1_HHZZ1_T);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM1_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (6,  [(1, C_12_34)]), D_FM1_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (6,  [(1, C_14_23)]), D_FM1_HHZZ1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM1_HHWW0_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (2,  [(1, C_13_42)]), D_FM1_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (1,  [(1, C_14_23)]), D_FM1_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (1,  [(1, C_13_42)]), D_FM1_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (2,  [(1, C_14_23)]), D_FM1_HHWW0_T);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (3,  [(1, C_14_23)]), D_FM1_HHWW1_S);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM1_HHWW1_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (4,  [(1, C_13_42)]), D_FM1_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (5,  [(1, C_12_34)]), D_FM1_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (8,  [(1, C_14_23)]), D_FM1_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (7,  [(1, C_12_34)]), D_FM1_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (5,  [(1, C_13_42)]), D_FM1_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (4,  [(1, C_12_34)]), D_FM1_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (7,  [(1, C_14_23)]), D_FM1_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (8,  [(1, C_12_34)]), D_FM1_HHWW1_T) ]
        else
            []
      else
          []
          
    let k_matrix_2scalar_2gauge_m_7 =
      if Flags.k_matrix_tm then
        if Flags.higgs_matrix then
            [ ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM7_HHZZ0_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (0,  [(1, C_13_42)]), D_FM7_HHZZ0_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (0,  [(1, C_14_23)]), D_FM7_HHZZ0_U);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf 
                   (3,  [(1, C_14_23)]), D_FM7_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (3,  [(1, C_13_42)]), D_FM7_HHZZ1_U);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (3,  [(1, C_12_34)]), D_FM7_HHZZ1_T);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM7_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (6,  [(1, C_12_34)]), D_FM7_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (6,  [(1, C_14_23)]), D_FM7_HHZZ1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM7_HHWW0_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (2,  [(1, C_13_42)]), D_FM7_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (1,  [(1, C_14_23)]), D_FM7_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (1,  [(1, C_13_42)]), D_FM7_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (2,  [(1, C_14_23)]), D_FM7_HHWW0_T);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (3,  [(1, C_14_23)]), D_FM7_HHWW1_S);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM7_HHWW1_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (4,  [(1, C_13_42)]), D_FM7_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (5,  [(1, C_12_34)]), D_FM7_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (8,  [(1, C_14_23)]), D_FM7_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (7,  [(1, C_12_34)]), D_FM7_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (5,  [(1, C_13_42)]), D_FM7_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (4,  [(1, C_12_34)]), D_FM7_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (7,  [(1, C_14_23)]), D_FM7_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (8,  [(1, C_12_34)]), D_FM7_HHWW1_T) ]
        else
            []
      else
          []      

    let k_matrix_4scalar =
      if Flags.k_matrix then
        if Flags.higgs_matrix then
            [ ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
	           (0,  [(1, C_12_34)]), D_Alpha_HHHH_S);
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
	           (0, [(1, C_13_42)]), D_Alpha_HHHH_T); 
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
                   (0, [(1, C_14_23)]), D_Alpha_HHHH_T); 
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
                   (3, [(1, C_14_23)]), D_Alpha_HHHH_S);
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
                   (3, [(1, C_13_42)]), D_Alpha_HHHH_T);
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
                   (3, [(1, C_12_34)]), D_Alpha_HHHH_T) ]
        else
            []
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
      @ k_matrix_quartic_gauge_t_0 @ k_matrix_quartic_gauge_t_1 @ k_matrix_quartic_gauge_t_2
      @ k_matrix_quartic_gauge_t_rsi
      @ k_matrix_quartic_gauge_m_0 @ k_matrix_quartic_gauge_m_1 @ k_matrix_quartic_gauge_m_7

    let standard_gauge_higgs =
      [ ((O H, G Wp, G Wm), Scalar_Vector_Vector 1, G_HWW);
        ((O H, G Z, G Z), Scalar_Vector_Vector 1, G_HZZ) ]

    let standard_gauge_higgs4 =
      [ (O H, O H, G Wp, G Wm), Scalar2_Vector2 1, G_HHWW;
        (O H, O H, G Z, G Z), Scalar2_Vector2 1, G_HHZZ ]

    let dim8_gauge_higgs4 =
      [ (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_1 1, FS0_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_1 1, FS0_HHZZ;
        (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_2 1, FS1_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_2 1, FS1_HHZZ ]
    
    let dim8_gauge_higgs4_m =
      [ (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_m_0 1, FM0_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_m_0 1, FM0_HHZZ;
        (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_m_1 1, FM1_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_m_1 1, FM1_HHZZ;
        (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_m_7 1, FM7_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_m_7 1, FM7_HHZZ]
       
    let standard_higgs =
      [ (O H, O H, O H), Scalar_Scalar_Scalar 1, G_H3 ]

    let standard_higgs4 =
      [ (O H, O H, O H, O H), Scalar4 1, G_H4 ]

    let fs_higgs4 =
      [ (O H, O H, O H, O H), Dim8_Scalar4 1, FS_H4 ]



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

(* Breaking thinks up
   \begin{subequations}
   \begin{align}
     \mathcal{L}^{\tau,H^4}_4 &=
       \left\lbrack (\partial_{\mu}H)(\partial^{\mu}H) \right\rbrack^2 \\
     \mathcal{L}^{\tau,H^4}_5 &=
       \left\lbrack (\partial_{\mu}H)(\partial^{\mu}H) \right\rbrack^2
   \end{align}
   \end{subequations}
   and
   \begin{subequations}
   \begin{align}
     \mathcal{L}^{\tau,H^2V^2}_4 &= \frac{g^2v_{\mathrm{F}}^2}{2}
              (\partial_{\mu}H)(\partial^{\mu}H) V_{\mu}V^{\mu}   \\
     \mathcal{L}^{\tau,H^2V^2}_5 &= \frac{g^2v_{\mathrm{F}}^2}{2}
              (\partial_{\mu}H)(\partial_{\nu}H) V_{\mu}V_{\nu}
   \end{align}
   \end{subequations}
   i.\,e.
   \begin{subequations}
   \begin{align}
     \mathcal{L}^{\tau,H^2V^2}_4 &=
        \frac{g^2v_{\mathrm{F}}^2}{2}
          \left\lbrack
              (\partial_{\mu}H)(\partial^{\mu}H) W^+_{\nu}W^{-,\nu}
            + \frac{1}{2\cos^2\theta_{w}} (\partial_{\mu}H)(\partial^{\mu}H) Z_{\nu} Z^{\nu}
          \right\rbrack \\
     \mathcal{L}^{\tau,H^2V^2}_5 &=
          \frac{g^2v_{\mathrm{F}}^2}{2}
          \left\lbrack
              (W^{+,\mu}\partial_{\mu}H) (W^{-,\nu}\partial_{\nu}H)
            + \frac{1}{2\cos^2\theta_{w}} (Z^{\mu}\partial_{\mu}H)(Z^{\nu}\partial_{\nu}H)
          \right\rbrack
   \end{align}
   \end{subequations} *)

(* \begin{multline}
     \tau^4_8 \mathcal{L}^{\tau,H^2V^2}_4 + \tau^5_8 \mathcal{L}^{\tau,H^2V^2}_5 = \\
       - \frac{g^2v_{\mathrm{F}}^2}{2} \Biggl\lbrack
            2\tau^4_8
              \frac{1}{2}(\ii\partial_{\mu}H)(\ii\partial^{\mu}H) W^+_{\nu}W^{-,\nu}
          + \tau^5_8
              (W^{+,\mu}\ii\partial_{\mu}H) (W^{-,\nu}\ii\partial_{\nu}H) \\
          + \frac{2\tau^4_8}{\cos^2\theta_{w}}
              \frac{1}{4} (\ii\partial_{\mu}H)(\ii\partial^{\mu}H) Z_{\nu} Z^{\nu}
          + \frac{\tau^5_8}{\cos^2\theta_{w}}
              \frac{1}{2} (Z^{\mu}\ii\partial_{\mu}H)(Z^{\nu}\ii\partial_{\nu}H)
          \Biggr\rbrack
   \end{multline}
   where the two powers of $\ii$ make the sign conveniently negative,
   i.\,e.
   \begin{subequations}
   \begin{align}
     \alpha_{(\partial H)^2W^2}^2 &= \tau^4_8 g^2v_{\mathrm{F}}^2\\
     \alpha_{(\partial HW)^2}^2 &= \frac{\tau^5_8 g^2v_{\mathrm{F}}^2}{2}  \\
     \alpha_{(\partial H)^2Z^2}^2 &= \frac{\tau^4_8 g^2v_{\mathrm{F}}^2}{\cos^2\theta_{w}} \\ 
     \alpha_{(\partial HZ)^2}^2 &=\frac{\tau^5_8 g^2v_{\mathrm{F}}^2}{2\cos^2\theta_{w}}
   \end{align}
   \end{subequations} *)

    let anomalous_gauge_higgs =
      [ (O H, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_HGaGa_anom;
        (O H, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_HGaZ_anom;
        (O H, G Z, G Z), Dim5_Scalar_Gauge2 1, G_HZZ_anom;
        (O H, G Wp, G Wm), Dim5_Scalar_Gauge2 1, G_HWW_anom;
        (O H, G Ga, G Z), Dim5_Scalar_Vector_Vector_U 1, G_HGaZ_u;
        (O H, G Z, G Z), Dim5_Scalar_Vector_Vector_U 1, G_HZZ_u;
        (O H, G Wp, G Wm), Dim5_Scalar_Vector_Vector_U 1, G_HWW_u;
        (O H, G Wm, G Wp), Dim5_Scalar_Vector_Vector_U 1, G_HWW_u
      ]

    let anomalous_gauge_higgs4 =
      []

    let anomalous_higgs =
      []

    let higgs_triangle_vertices = 
      if Flags.higgs_triangle then
        [ (O H, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
          (O H, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_HGaZ;
          (O H, G Gl, G Gl), Dim5_Scalar_Gauge2 1, G_Hgg ]
      else
        []

    let anomalous_higgs4 =
      []

    let gauge_higgs =
      if Flags.higgs_anom then
        standard_gauge_higgs @ anomalous_gauge_higgs
      else
        standard_gauge_higgs


    let gauge_higgs4 =
      ( if Flags.higgs_anom then
          standard_gauge_higgs4 @ anomalous_gauge_higgs4
        else
          standard_gauge_higgs4 ) @
      ( if Flags.higgs_matrix then
          (dim8_gauge_higgs4 @ dim8_gauge_higgs4_m @ k_matrix_2scalar_2gauge 
           @ k_matrix_2scalar_2gauge_m @ k_matrix_2scalar_2gauge_m_1 @ k_matrix_2scalar_2gauge_m_7)
	 else
	   [] )

    let higgs =
      if Flags.higgs_anom then
        standard_higgs @ anomalous_higgs
      else
        standard_higgs

    let higgs4 =
      ( if Flags.higgs_anom then
          standard_higgs4 @ anomalous_higgs4
        else
          standard_higgs4 ) @ 
      ( if Flags.higgs_matrix then
          (fs_higgs4 @ k_matrix_4scalar )
	 else
	   [] )


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
    Effectivly this doesn't matter for SSC, because $(-1)^2=1$.
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
   Old representation
   \begin{equation}
    \mathcal{L}_{\sigma}=
              -\frac{g_\sigma v}{2} \text{tr}
	      \left\lbrack \mathbf{V}_\mu \mathbf{V}^\mu \right\rbrack \sigma
   \end{equation}
*)

(* \begin{dubious}
   Transversal couplings like rsigma3t and rf3t are to be calculated in the new
   higgs matrix representation.
   \end{dubious} *)


    let rsigma3 =
      [ ((O Rsigma, G Wp, G Wm), Scalar_Vector_Vector 1, G_SWW);
        ((O Rsigma, G Z, G Z), Scalar_Vector_Vector 1, G_SZZ) ]

    let rsigma3h =
      [ ((O Rsigma, O H, O H), Dim5_Scalar_Scalar2 1, G_SHH) ]

    let rsigma3t =
      [ ((O Rsigma, G Wp, G Wm), Scalar_Vector_Vector_t 1, G_SWW_T);
        ((O Rsigma, G Z, G Z), Scalar_Vector_Vector_t 1, G_SZZ_T);
        ((O Rsigma, G Ga, G Ga), Scalar_Vector_Vector_t 1, G_SAA_T);
        ((O Rsigma, G Ga, G Z), Scalar_Vector_Vector_t 1, G_SAZ_T) ]

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
        ((O Rphisn, G Wp, G Wm), Scalar_Vector_Vector 1, G_PSNWW);
        ((O Rphisn, G Z, G Z), Scalar_Vector_Vector 1, G_PSNZZ) ;
        ((O Rphip, G Z, G Wm), Scalar_Vector_Vector 1, G_PWZ) ;
        ((O Rphipp, G Wm, G Wm), Scalar_Vector_Vector 1, G_PWW) ;
        ((O Rphim, G Wp, G Z), Scalar_Vector_Vector 1, G_PWZ) ;
        ((O Rphimm, G Wp, G Wp), Scalar_Vector_Vector 1, G_PWW) ]

    let rphi3h =
      [ ((O Rphisn, O H, O H), Dim5_Scalar_Scalar2 1, G_PSNHH) ]

(* Tensor IsoScalar
*)
    let rf3 =
      [ ((O Rf, G Wp, G Wm), Tensor_2_Vector_Vector_1 1, G_FWW);
        ((O Rf, G Z, G Z), Tensor_2_Vector_Vector_1 1, G_FZZ) ]
    
    let rf3cf =
      [ ((O Rf, G Wp, G Wm), Tensor_2_Vector_Vector 1, G_FWW);
        ((O Rf, G Z, G Z), Tensor_2_Vector_Vector 1, G_FZZ);
        ((O Rf, G Wp, G Wm), Tensor_2_Vector_Vector_cf 1, G_FWW_CF);
        ((O Rf, G Z, G Z), Tensor_2_Vector_Vector_cf 1, G_FZZ_CF) ]

    let rf3h =
      [ ((O Rf, O H, O H), Tensor_2_Scalar_Scalar 1, G_FHH);
        ((O Rf, O H, O H), Tensor_2_Scalar_Scalar_cf 1, G_FHH_CF) ]
     
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
        ((O Rtsn, G Wp, G Wm), Tensor_2_Vector_Vector_1 1, G_TSNWW);
        ((O Rtsn, G Z, G Z), Tensor_2_Vector_Vector_1 1, G_TSNZZ) ;
        ((O Rtp, G Z, G Wm), Tensor_2_Vector_Vector_1 1, G_TWZ) ;
        ((O Rtpp, G Wm, G Wm), Tensor_2_Vector_Vector_1 1, G_TWW) ;
        ((O Rtm, G Wp, G Z), Tensor_2_Vector_Vector_1 1, G_TWZ) ;
        ((O Rtmm, G Wp, G Wp), Tensor_2_Vector_Vector_1 1, G_TWW) ]

    let rt3cf =
      [ ((O Rtn, G Wp, G Wm), Tensor_2_Vector_Vector 1, G_TNWW);
        ((O Rtn, G Z, G Z), Tensor_2_Vector_Vector 1, G_TNZZ) ;
        ((O Rtsn, G Wp, G Wm), Tensor_2_Vector_Vector 1, G_TSNWW);
        ((O Rtsn, G Z, G Z), Tensor_2_Vector_Vector 1, G_TSNZZ) ;
        ((O Rtp, G Z, G Wm), Tensor_2_Vector_Vector 1, G_TWZ) ;
        ((O Rtpp, G Wm, G Wm), Tensor_2_Vector_Vector 1, G_TWW) ;
        ((O Rtm, G Wp, G Z), Tensor_2_Vector_Vector 1, G_TWZ) ;
        ((O Rtmm, G Wp, G Wp), Tensor_2_Vector_Vector 1, G_TWW);
        ((O Rtn, G Wp, G Wm), Tensor_2_Vector_Vector_cf 1, G_TNWW_CF);
        ((O Rtn, G Z, G Z), Tensor_2_Vector_Vector_cf 1, G_TNZZ_CF) ;
        ((O Rtsn, G Wp, G Wm), Tensor_2_Vector_Vector_cf 1, G_TSNWW_CF);
        ((O Rtsn, G Z, G Z), Tensor_2_Vector_Vector_cf 1, G_TSNZZ_CF) ;
        ((O Rtp, G Z, G Wm), Tensor_2_Vector_Vector_cf 1, G_TWZ_CF) ;
        ((O Rtpp, G Wm, G Wm), Tensor_2_Vector_Vector_cf 1, G_TWW_CF) ;
        ((O Rtm, G Wp, G Z), Tensor_2_Vector_Vector_cf 1, G_TWZ_CF) ;
        ((O Rtmm, G Wp, G Wp), Tensor_2_Vector_Vector_cf 1, G_TWW_CF) ]


(* Anomalous trilinear interactions $f_i f_j V$ and $ttH$:
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

(* \begin{equation}
     \Delta\mathcal{L}_{ttH} =
        - \frac{1}{\sqrt{2}} \bar{t} (Y_V(k^2)+iY_A(k^2)\gamma_5)t H
   \end{equation} *)

    let anomalous_ttH =
      if Flags.top_anom then
        [ ((M (U (-3)), O H, M (U 3)), FBF (1, Psibar, SPM, Psi), G_SP_ttH) ]
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
       yukawa @ triple_gauge @
       gauge_higgs @ higgs @ higgs_triangle_vertices 
       @ goldstone_vertices @
       rsigma3 @ rsigma3t @ rphi3 @
       ( if Flags.cf_arbitrary then
	    (rf3cf @ rt3cf)
	 else
	    (rf3 @ rt3) ) @
       rf3t @ 
       ( if Flags.higgs_matrix then
	    (rsigma3h @ rphi3h @ rf3h )
	 else
	    [] ) @
       anomalous_ttA @ anomalous_bbA @
       anomalous_ttZ @ anomalous_bbZ @
       anomalous_tbW @ anomalous_tbWA @ anomalous_tbWZ @
       anomalous_ttWW @ anomalous_bbWW @
       anomalous_ttG @ anomalous_ttGG @
       anomalous_ttH @
       anomalous_top_qGuG @ anomalous_top_qBuB @
       anomalous_top_qW @ anomalous_top_DuDd)

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
      | "H" -> O H | "Rsigma" -> O Rsigma
      | "Rphi0" -> O Rphin
      | "Rphis0" -> O Rphisn
      | "Rphi+" -> O Rphip |  "Rphi-" -> O Rphim
      | "Rphi++" -> O Rphip |  "Rphi--" -> O Rphimm
      | "Rf" -> O Rf
      | "Rt0" -> O Rtn
      | "Rts0" -> O Rtsn
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
      | _ -> invalid_arg "Modellib_BSM.SSC.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib_BSM.SSC.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib_BSM.SSC.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib_BSM.SSC.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_BSM.SSC.flavor_to_string: invalid down type quark"
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
          | H -> "H" | Rsigma -> "Rsigma"
          | Rphin -> "Rphin"  | Rphisn -> "Rphisn"
          | Rphip -> "Rphi+" | Rphim -> "Rphi-"
          | Rphipp -> "Rphi++" | Rphimm -> "Rphi--"
          | Rf -> "Rf"
          | Rtn -> "Rtn" | Rtsn -> "Rtsn" | Rtp -> "Rt+" | Rtm -> "Rt-"
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
                "Modellib_BSM.SSC.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib_BSM.SSC.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib_BSM.SSC.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Modellib_BSM.SSC.flavor_to_TeX: invalid down type quark"
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
          | H -> "H" | Rsigma -> "\\sigma"
          | Rphip -> "\\phi^+" | Rphim -> "\\phi^-" | Rphin -> "\\phi^0" 
          | Rphisn -> "\\phi_s^0" 
          | Rphipp -> "\\phi^{++}" | Rphimm -> "\\phi^{--}"
          | Rf -> "f"
          | Rtp -> "t^+" | Rtm -> "t^-" | Rtn -> "t^0" | Rtsn -> "t_s^0" 
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
          | H -> "h" | Rsigma -> "rsi"
          | Rphip -> "rpp" | Rphim -> "rpm" | Rphin -> "rpn"
          | Rphisn -> "rpsn"
          | Rphipp -> "rppp" | Rphimm -> "rpmm"
          | Rf -> "rf"
          | Rtp -> "rtp" | Rtm -> "rtm" | Rtn -> "rtn"
          | Rtsn -> "rtsn"
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
          | H -> 25 | Rsigma -> 45
          | Rphin -> 46 | Rphip | Rphim -> 47 
          | Rphipp | Rphimm -> 48
          | Rphisn -> 49   
          | Rf -> 52
          | Rtn -> 53 | Rtp | Rtm -> 54 
          | Rtpp | Rtmm -> 55
          | Rtsn -> 59
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
      | G_SP_ttH -> "gsp_tth"
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
      | FT0_WWWW0 -> "at0ww0" | FT0_WWWW2 -> "at0ww2"
      | FT0_ZZWW0 -> "at0zw0" | FT0_ZZWW1 -> "at0zw1"
      | FT0_ZZZZ  -> "at0zz"  | FT0_AAAA  -> "at0aa"
      | FT0_AAWW0 -> "at0aw0" | FT0_AAWW1 -> "at0aw1"
      | FT0_AAZZ -> "at0az"   
      | FT0_AZWW0 -> "at0azw0" | FT0_AZWW1 -> "at0azw1"
      | FT0_AAAZ  -> "at03az"  | FT0_AZZZ  -> "at0a3z"
      | FT1_WWWW0 -> "at1ww0" | FT1_WWWW2 -> "at1ww2"
      | FT1_ZZWW0 -> "at1zw0" | FT1_ZZWW1 -> "at1zw1"
      | FT1_ZZZZ  -> "at1zz"  | FT1_AAAA  -> "at1aa" 
      | FT1_AAWW0 -> "at1aw0" | FT1_AAWW1 -> "at1aw1"
      | FT1_AAZZ -> "at1az"   
      | FT1_AZWW0 -> "at1azw0" | FT1_AZWW1 -> "at1azw1"
      | FT1_AAAZ  -> "at13az"  | FT1_AZZZ  -> "at1a3z"
      | FT2_WWWW0 -> "at2ww0" | FT2_WWWW2 -> "at2ww2"
      | FT2_ZZWW0 -> "at2zw0" | FT2_ZZWW1 -> "at2zw1"
      | FT2_ZZZZ  -> "at2zz"  | FT2_AAAA  -> "at2aa"
      | FT2_AAWW0 -> "at2aw0" | FT2_AAWW1 -> "at2aw1"
      | FT2_AAZZ -> "at2az"   
      | FT2_AZWW0 -> "at2azw0" | FT2_AZWW1 -> "at2azw1"
      | FT2_AAAZ  -> "at23az"  | FT2_AZZZ  -> "at2a3z"
      | FM0_WWWW0 -> "am0ww0,am0ww0" | FM0_WWWW2 -> "am0ww2,am0ww2"
      | FM0_ZZWW0 -> "am0zw0/costhw**2,am0zw0*costhw**2" | FM0_ZZWW1 -> "am0zw1/costhw**2,am0zw1*costhw**2"
      | FM0_ZZZZ  -> "am0zz,am0zz" 
      | FM1_WWWW0 -> "am1ww0,am1ww0" | FM1_WWWW2 -> "am1ww2,am1ww2"
      | FM1_ZZWW0 -> "am1zw0/costhw**2,am1zw0*costhw**2" | FM1_ZZWW1 -> "am1zw1/costhw**2,am1zw1*costhw**2"
      | FM1_ZZZZ  -> "am1zz,am1zz"  
      | FM7_WWWW0 -> "am7ww0,am7ww0,am7ww0" | FM7_WWWW2 -> "am7ww2,am7ww2,am7ww2"
      | FM7_ZZWW0 -> "am7zw0/costhw**2,am7zw0,am7zw0*costhw**2" | FM7_ZZWW1 -> "am7zw1/costhw**2,am7zw1,am7zw1*costhw**2"
      | FM7_ZZZZ  -> "am7zz,am7zz,am7zz"
      | FS0_HHWW -> "fs0hhww" | FS0_HHZZ -> "fs0hhzz"
      | FS1_HHWW -> "fs1hhww" | FS1_HHZZ -> "fs1hhzz"
      | FS_H4 -> "fsh4"
      | FM0_HHWW -> "fm0hhww" | FM0_HHZZ -> "fm0hhzz" 
      | FM1_HHWW -> "fm1hhww" | FM1_HHZZ -> "fm1hhzz"  
      | FM7_HHWW -> "fm7hhww" | FM7_HHZZ -> "fm7hhzz"
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
      | D_Alpha_ZZZZ_S  -> "dalz4_s(gkm,mkm,"
      | D_Alpha_ZZZZ_T  -> "dalz4_t(gkm,mkm,"
      | D_FT0_ZZWW0_S -> "datzz0_s_0(gkm,mkm,"
      | D_FT0_ZZWW0_T -> "datzz0_t_0(gkm,mkm,"
      | D_FT0_ZZWW0_U -> "datzz0_u_0(gkm,mkm,"
      | D_FT0_ZZWW1_S -> "datzz1_s_0(gkm,mkm,"
      | D_FT0_ZZWW1_T -> "datzz1_t_0(gkm,mkm,"
      | D_FT0_ZZWW1_U -> "datzz1_u_0(gkm,mkm,"
      | D_FT0_WWWW0_S -> "datww0_s_0(gkm,mkm,"
      | D_FT0_WWWW0_T -> "datww0_t_0(gkm,mkm,"
      | D_FT0_WWWW0_U -> "datww0_u_0(gkm,mkm,"
      | D_FT0_WWWW2_S -> "datww2_s_0(gkm,mkm,"
      | D_FT0_WWWW2_T -> "datww2_t_0(gkm,mkm,"
      | D_FT0_WWWW2_U -> "datww2_u_0(gkm,mkm,"
      | D_FT0_ZZZZ_S  -> "datz4_s_0(gkm,mkm,"
      | D_FT0_ZZZZ_T  -> "datz4_t_0(gkm,mkm,"
      | D_FT0_ZZZZ_U  -> "datz4_u_0(gkm,mkm,"
      | D_FT0_AAAA_S  -> "data4_s_0(gkm,mkm,"
      | D_FT0_AAAA_T  -> "data4_t_0(gkm,mkm,"
      | D_FT0_AAAA_U  -> "data4_u_0(gkm,mkm," 
      | D_FT0_AAWW0_S -> "dataw0_s_0(gkm,mkm,"
      | D_FT0_AAWW0_T -> "dataw0_t_0(gkm,mkm,"
      | D_FT0_AAWW0_U -> "dataw0_u_0(gkm,mkm,"
      | D_FT0_AAWW1_S -> "dataw1_s_0(gkm,mkm,"
      | D_FT0_AAWW1_T -> "dataw1_t_0(gkm,mkm,"
      | D_FT0_AAWW1_U -> "dataw1_u_0(gkm,mkm,"
      | D_FT0_AAZZ_S  -> "dataz_s_0(gkm,mkm,"
      | D_FT0_AAZZ_T  -> "dataz_t_0(gkm,mkm,"
      | D_FT0_AAZZ_U  -> "dataz_u_0(gkm,mkm," 
      | D_FT0_AZWW0_S -> "datazw0_s_0(gkm,mkm,"
      | D_FT0_AZWW0_T -> "datazw0_t_0(gkm,mkm,"
      | D_FT0_AZWW0_U -> "datazw0_u_0(gkm,mkm,"
      | D_FT0_AZWW1_S -> "datazw1_s_0(gkm,mkm,"
      | D_FT0_AZWW1_T -> "datazw1_t_0(gkm,mkm,"
      | D_FT0_AZWW1_U -> "datazw1_u_0(gkm,mkm," 
      | D_FT0_AAAZ_S -> "dat3az_s_0(gkm,mkm,"
      | D_FT0_AAAZ_T -> "dat3az_t_0(gkm,mkm,"
      | D_FT0_AAAZ_U -> "dat3az_u_0(gkm,mkm," 
      | D_FT0_AZZZ_S -> "data3z_s_0(gkm,mkm,"
      | D_FT0_AZZZ_T -> "data3z_t_0(gkm,mkm,"
      | D_FT0_AZZZ_U -> "data3z_u_0(gkm,mkm,"            
      | D_FT1_ZZWW0_S -> "datzz0_s_1(gkm,mkm,"
      | D_FT1_ZZWW0_T -> "datzz0_t_1(gkm,mkm,"
      | D_FT1_ZZWW0_U -> "datzz0_u_1(gkm,mkm,"
      | D_FT1_ZZWW1_S -> "datzz1_s_1(gkm,mkm,"
      | D_FT1_ZZWW1_T -> "datzz1_t_1(gkm,mkm,"
      | D_FT1_ZZWW1_U -> "datzz1_u_1(gkm,mkm,"
      | D_FT1_WWWW0_S -> "datww0_s_1(gkm,mkm,"
      | D_FT1_WWWW0_T -> "datww0_t_1(gkm,mkm,"
      | D_FT1_WWWW0_U -> "datww0_u_1(gkm,mkm,"
      | D_FT1_WWWW2_S -> "datww2_s_1(gkm,mkm,"
      | D_FT1_WWWW2_T -> "datww2_t_1(gkm,mkm,"
      | D_FT1_WWWW2_U -> "datww2_u_1(gkm,mkm,"
      | D_FT1_ZZZZ_S  -> "datz4_s_1(gkm,mkm,"
      | D_FT1_ZZZZ_T  -> "datz4_t_1(gkm,mkm,"
      | D_FT1_ZZZZ_U  -> "datz4_u_1(gkm,mkm,"
      | D_FT1_AAAA_S  -> "data4_s_1(gkm,mkm,"
      | D_FT1_AAAA_T  -> "data4_t_1(gkm,mkm,"
      | D_FT1_AAAA_U  -> "data4_u_1(gkm,mkm," 
      | D_FT1_AAWW0_S -> "dataw0_s_1(gkm,mkm,"
      | D_FT1_AAWW0_T -> "dataw0_t_1(gkm,mkm,"
      | D_FT1_AAWW0_U -> "dataw0_u_1(gkm,mkm,"
      | D_FT1_AAWW1_S -> "dataw1_s_1(gkm,mkm,"
      | D_FT1_AAWW1_T -> "dataw1_t_1(gkm,mkm,"
      | D_FT1_AAWW1_U -> "dataw1_u_1(gkm,mkm,"
      | D_FT1_AAZZ_S  -> "dataz_s_1(gkm,mkm,"
      | D_FT1_AAZZ_T  -> "dataz_t_1(gkm,mkm,"
      | D_FT1_AAZZ_U  -> "dataz_u_1(gkm,mkm,"
      | D_FT1_AZWW0_S -> "datazw0_s_1(gkm,mkm,"
      | D_FT1_AZWW0_T -> "datazw0_t_1(gkm,mkm,"
      | D_FT1_AZWW0_U -> "datazw0_u_1(gkm,mkm,"
      | D_FT1_AZWW1_S -> "datazw1_s_1(gkm,mkm,"
      | D_FT1_AZWW1_T -> "datazw1_t_1(gkm,mkm,"
      | D_FT1_AZWW1_U -> "datazw1_u_1(gkm,mkm,"
      | D_FT1_AAAZ_S -> "dat3az_s_1(gkm,mkm,"
      | D_FT1_AAAZ_T -> "dat3az_t_1(gkm,mkm,"
      | D_FT1_AAAZ_U -> "dat3az_u_1(gkm,mkm," 
      | D_FT1_AZZZ_S -> "data3z_s_1(gkm,mkm,"
      | D_FT1_AZZZ_T -> "data3z_t_1(gkm,mkm,"
      | D_FT1_AZZZ_U -> "data3z_u_1(gkm,mkm,"      
      | D_FT2_ZZWW0_S -> "datzz0_s_2(gkm,mkm,"
      | D_FT2_ZZWW0_T -> "datzz0_t_2(gkm,mkm,"
      | D_FT2_ZZWW0_U -> "datzz0_u_2(gkm,mkm,"
      | D_FT2_ZZWW1_S -> "datzz1_s_2(gkm,mkm,"
      | D_FT2_ZZWW1_T -> "datzz1_t_2(gkm,mkm,"
      | D_FT2_ZZWW1_U -> "datzz1_u_2(gkm,mkm,"
      | D_FT2_WWWW0_S -> "datww0_s_2(gkm,mkm,"
      | D_FT2_WWWW0_T -> "datww0_t_2(gkm,mkm,"
      | D_FT2_WWWW0_U -> "datww0_u_2(gkm,mkm,"
      | D_FT2_WWWW2_S -> "datww2_s_2(gkm,mkm,"
      | D_FT2_WWWW2_T -> "datww2_t_2(gkm,mkm,"
      | D_FT2_WWWW2_U -> "datww2_u_2(gkm,mkm,"
      | D_FT2_ZZZZ_S  -> "datz4_s_2(gkm,mkm,"
      | D_FT2_ZZZZ_T  -> "datz4_t_2(gkm,mkm,"
      | D_FT2_ZZZZ_U  -> "datz4_u_2(gkm,mkm,"
      | D_FT2_AAAA_S  -> "data4_s_2(gkm,mkm,"
      | D_FT2_AAAA_T  -> "data4_t_2(gkm,mkm,"
      | D_FT2_AAAA_U  -> "data4_u_2(gkm,mkm," 
      | D_FT2_AAWW0_S -> "dataw0_s_2(gkm,mkm,"
      | D_FT2_AAWW0_T -> "dataw0_t_2(gkm,mkm,"
      | D_FT2_AAWW0_U -> "dataw0_u_2(gkm,mkm,"
      | D_FT2_AAWW1_S -> "dataw1_s_2(gkm,mkm,"
      | D_FT2_AAWW1_T -> "dataw1_t_2(gkm,mkm,"
      | D_FT2_AAWW1_U -> "dataw1_u_2(gkm,mkm,"
      | D_FT2_AAZZ_S  -> "dataz_s_2(gkm,mkm,"
      | D_FT2_AAZZ_T  -> "dataz_t_2(gkm,mkm,"
      | D_FT2_AAZZ_U  -> "dataz_u_2(gkm,mkm,"
      | D_FT2_AZWW0_S -> "datazw0_s_2(gkm,mkm,"
      | D_FT2_AZWW0_T -> "datazw0_t_2(gkm,mkm,"
      | D_FT2_AZWW0_U -> "datazw0_u_2(gkm,mkm,"
      | D_FT2_AZWW1_S -> "datazw1_s_2(gkm,mkm,"
      | D_FT2_AZWW1_T -> "datazw1_t_2(gkm,mkm,"
      | D_FT2_AZWW1_U -> "datazw1_u_2(gkm,mkm,"
      | D_FT2_AAAZ_S -> "dat3az_s_2(gkm,mkm,"
      | D_FT2_AAAZ_T -> "dat3az_t_2(gkm,mkm,"
      | D_FT2_AAAZ_U -> "dat3az_u_2(gkm,mkm," 
      | D_FT2_AZZZ_S -> "data3z_s_2(gkm,mkm,"
      | D_FT2_AZZZ_T -> "data3z_t_2(gkm,mkm,"
      | D_FT2_AZZZ_U -> "data3z_u_2(gkm,mkm,"
      | D_FTrsi_ZZWW0_S -> "datzz0_s_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW0_T -> "datzz0_t_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW0_U -> "datzz0_u_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW1_S -> "datzz1_s_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW1_T -> "datzz1_t_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW1_U -> "datzz1_u_rsi(gkm,mkm,"
      | D_FTrsi_WWWW0_S -> "datww0_s_rsi(gkm,mkm,"
      | D_FTrsi_WWWW0_T -> "datww0_t_rsi(gkm,mkm,"
      | D_FTrsi_WWWW0_U -> "datww0_u_rsi(gkm,mkm,"
      | D_FTrsi_WWWW2_S -> "datww2_s_rsi(gkm,mkm,"
      | D_FTrsi_WWWW2_T -> "datww2_t_rsi(gkm,mkm,"
      | D_FTrsi_WWWW2_U -> "datww2_u_rsi(gkm,mkm,"
      | D_FTrsi_ZZZZ_S  -> "datz4_s_rsi(gkm,mkm,"
      | D_FTrsi_ZZZZ_T  -> "datz4_t_rsi(gkm,mkm,"
      | D_FTrsi_ZZZZ_U  -> "datz4_u_rsi(gkm,mkm," 
      | D_FTrsi_AAAA_S  -> "data4_s_rsi(gkm,mkm,"
      | D_FTrsi_AAAA_T  -> "data4_t_rsi(gkm,mkm,"
      | D_FTrsi_AAAA_U  -> "data4_u_rsi(gkm,mkm,"  
      | D_FTrsi_AAWW0_S -> "dataw0_s_rsi(gkm,mkm,"
      | D_FTrsi_AAWW0_T -> "dataw0_t_rsi(gkm,mkm,"
      | D_FTrsi_AAWW0_U -> "dataw0_u_rsi(gkm,mkm,"
      | D_FTrsi_AAWW1_S -> "dataw1_s_rsi(gkm,mkm,"
      | D_FTrsi_AAWW1_T -> "dataw1_t_rsi(gkm,mkm,"
      | D_FTrsi_AAWW1_U -> "dataw1_u_rsi(gkm,mkm,"
      | D_FTrsi_AAZZ_S  -> "dataz_s_rsi(gkm,mkm,"
      | D_FTrsi_AAZZ_T  -> "dataz_t_rsi(gkm,mkm,"
      | D_FTrsi_AAZZ_U  -> "dataz_u_rsi(gkm,mkm,"    
      | D_FTrsi_AZWW0_S -> "datazw0_s_rsi(gkm,mkm,"
      | D_FTrsi_AZWW0_T -> "datazw0_t_rsi(gkm,mkm,"
      | D_FTrsi_AZWW0_U -> "datazw0_u_rsi(gkm,mkm,"
      | D_FTrsi_AZWW1_S -> "datazw1_s_rsi(gkm,mkm,"
      | D_FTrsi_AZWW1_T -> "datazw1_t_rsi(gkm,mkm,"
      | D_FTrsi_AZWW1_U -> "datazw1_u_rsi(gkm,mkm,"
      | D_FTrsi_AAAZ_S -> "dat3az_s_rsi(gkm,mkm,"
      | D_FTrsi_AAAZ_T -> "dat3az_t_rsi(gkm,mkm,"
      | D_FTrsi_AAAZ_U -> "dat3az_u_rsi(gkm,mkm," 
      | D_FTrsi_AZZZ_S -> "data3z_s_rsi(gkm,mkm,"
      | D_FTrsi_AZZZ_T -> "data3z_t_rsi(gkm,mkm,"
      | D_FTrsi_AZZZ_U -> "data3z_u_rsi(gkm,mkm,"       
      | D_FM0_ZZWW0_S -> "damzz0_s_0(gkm,mkm,"
      | D_FM0_ZZWW0_T -> "damzz0_t_0(gkm,mkm,"
      | D_FM0_ZZWW0_U -> "damzz0_u_0(gkm,mkm,"
      | D_FM0_ZZWW1_S -> "damzz1_s_0(gkm,mkm,"
      | D_FM0_ZZWW1_T -> "damzz1_t_0(gkm,mkm,"
      | D_FM0_ZZWW1_U -> "damzz1_u_0(gkm,mkm,"
      | D_FM0_WWWW0_S -> "damww0_s_0(gkm,mkm,"
      | D_FM0_WWWW0_T -> "damww0_t_0(gkm,mkm,"
      | D_FM0_WWWW0_U -> "damww0_u_0(gkm,mkm,"
      | D_FM0_WWWW2_S -> "damww2_s_0(gkm,mkm,"
      | D_FM0_WWWW2_T -> "damww2_t_0(gkm,mkm,"
      | D_FM0_WWWW2_U -> "damww2_u_0(gkm,mkm,"
      | D_FM0_ZZZZ_S  -> "damz4_s_0(gkm,mkm,"
      | D_FM0_ZZZZ_T  -> "damz4_t_0(gkm,mkm,"
      | D_FM0_ZZZZ_U  -> "damz4_u_0(gkm,mkm,"
      | D_FM1_ZZWW0_S -> "damzz0_s_1(gkm,mkm,"
      | D_FM1_ZZWW0_T -> "damzz0_t_1(gkm,mkm,"
      | D_FM1_ZZWW0_U -> "damzz0_u_1(gkm,mkm,"
      | D_FM1_ZZWW1_S -> "damzz1_s_1(gkm,mkm,"
      | D_FM1_ZZWW1_T -> "damzz1_t_1(gkm,mkm,"
      | D_FM1_ZZWW1_U -> "damzz1_u_1(gkm,mkm,"
      | D_FM1_WWWW0_S -> "damww0_s_1(gkm,mkm,"
      | D_FM1_WWWW0_T -> "damww0_t_1(gkm,mkm,"
      | D_FM1_WWWW0_U -> "damww0_u_1(gkm,mkm,"
      | D_FM1_WWWW2_S -> "damww2_s_1(gkm,mkm,"
      | D_FM1_WWWW2_T -> "damww2_t_1(gkm,mkm,"
      | D_FM1_WWWW2_U -> "damww2_u_1(gkm,mkm,"
      | D_FM1_ZZZZ_S  -> "damz4_s_1(gkm,mkm,"
      | D_FM1_ZZZZ_T  -> "damz4_t_1(gkm,mkm,"
      | D_FM1_ZZZZ_U  -> "damz4_u_1(gkm,mkm,"
      | D_FM7_ZZWW0_S -> "damzz0_s_7(gkm,mkm,"
      | D_FM7_ZZWW0_T -> "damzz0_t_7(gkm,mkm,"
      | D_FM7_ZZWW0_U -> "damzz0_u_7(gkm,mkm,"
      | D_FM7_ZZWW1_S -> "damzz1_s_7(gkm,mkm,"
      | D_FM7_ZZWW1_T -> "damzz1_t_7(gkm,mkm,"
      | D_FM7_ZZWW1_U -> "damzz1_u_7(gkm,mkm,"
      | D_FM7_WWWW0_S -> "damww0_s_7(gkm,mkm,"
      | D_FM7_WWWW0_T -> "damww0_t_7(gkm,mkm,"
      | D_FM7_WWWW0_U -> "damww0_u_7(gkm,mkm,"
      | D_FM7_WWWW2_S -> "damww2_s_7(gkm,mkm,"
      | D_FM7_WWWW2_T -> "damww2_t_7(gkm,mkm,"
      | D_FM7_WWWW2_U -> "damww2_u_7(gkm,mkm,"
      | D_FM7_ZZZZ_S  -> "damz4_s_7(gkm,mkm,"
      | D_FM7_ZZZZ_T  -> "damz4_t_7(gkm,mkm,"
      | D_FM7_ZZZZ_U  -> "damz4_u_7(gkm,mkm,"
      | D_Alpha_HHHH_S  -> "dalh4_s(gkm,mkm,"
      | D_Alpha_HHHH_T  -> "dalh4_t(gkm,mkm,"
      | D_Alpha_HHWW0_S -> "dalhw0_s(gkm,mkm,"
      | D_Alpha_HHWW0_T -> "dalhw0_t(gkm,mkm,"
      | D_Alpha_HHZZ0_S -> "dalhz0_s(gkm,mkm,"
      | D_Alpha_HHZZ0_T -> "dalhz0_t(gkm,mkm,"
      | D_Alpha_HHWW1_S -> "dalhw1_s(gkm,mkm,"
      | D_Alpha_HHWW1_T -> "dalhw1_t(gkm,mkm,"
      | D_Alpha_HHWW1_U -> "dalhw1_u(gkm,mkm,"
      | D_Alpha_HHZZ1_S -> "dalhz1_s(gkm,mkm,"
      | D_Alpha_HHZZ1_T -> "dalhz1_t(gkm,mkm,"
      | D_Alpha_HHZZ1_U -> "dalhz1_u(gkm,mkm,"
      | D_FM0_HHWW0_S -> "damhw0_s_0(gkm,mkm,"
      | D_FM0_HHWW0_T -> "damhw0_t_0(gkm,mkm,"
      | D_FM0_HHWW0_U -> "damhw0_u_0(gkm,mkm,"
      | D_FM0_HHZZ0_S -> "damhz0_s_0(gkm,mkm,"
      | D_FM0_HHZZ0_T -> "damhz0_t_0(gkm,mkm,"
      | D_FM0_HHZZ0_U -> "damhz0_u_0(gkm,mkm,"
      | D_FM0_HHWW1_S -> "damhw1_s_0(gkm,mkm,"
      | D_FM0_HHWW1_T -> "damhw1_t_0(gkm,mkm,"
      | D_FM0_HHWW1_U -> "damhw1_u_0(gkm,mkm,"
      | D_FM0_HHZZ1_S -> "damhz1_s_0(gkm,mkm,"
      | D_FM0_HHZZ1_T -> "damhz1_t_0(gkm,mkm,"
      | D_FM0_HHZZ1_U -> "damhz1_u_0(gkm,mkm,"  
      | D_FM1_HHWW0_S -> "damhw0_s_1(gkm,mkm,"
      | D_FM1_HHWW0_T -> "damhw0_t_1(gkm,mkm,"
      | D_FM1_HHWW0_U -> "damhw0_u_1(gkm,mkm,"
      | D_FM1_HHZZ0_S -> "damhz0_s_1(gkm,mkm,"
      | D_FM1_HHZZ0_T -> "damhz0_t_1(gkm,mkm,"
      | D_FM1_HHZZ0_U -> "damhz0_u_1(gkm,mkm,"
      | D_FM1_HHWW1_S -> "damhw1_s_1(gkm,mkm,"
      | D_FM1_HHWW1_T -> "damhw1_t_1(gkm,mkm,"
      | D_FM1_HHWW1_U -> "damhw1_u_1(gkm,mkm,"
      | D_FM1_HHZZ1_S -> "damhz1_s_1(gkm,mkm,"
      | D_FM1_HHZZ1_T -> "damhz1_t_1(gkm,mkm,"
      | D_FM1_HHZZ1_U -> "damhz1_u_1(gkm,mkm," 
      | D_FM7_HHWW0_S -> "damhw0_s_7(gkm,mkm,"
      | D_FM7_HHWW0_T -> "damhw0_t_7(gkm,mkm,"
      | D_FM7_HHWW0_U -> "damhw0_u_7(gkm,mkm,"
      | D_FM7_HHZZ0_S -> "damhz0_s_7(gkm,mkm,"
      | D_FM7_HHZZ0_T -> "damhz0_t_7(gkm,mkm,"
      | D_FM7_HHZZ0_U -> "damhz0_u_7(gkm,mkm,"
      | D_FM7_HHWW1_S -> "damhw1_s_7(gkm,mkm,"
      | D_FM7_HHWW1_T -> "damhw1_t_7(gkm,mkm,"
      | D_FM7_HHWW1_U -> "damhw1_u_7(gkm,mkm,"
      | D_FM7_HHZZ1_S -> "damhz1_s_7(gkm,mkm,"
      | D_FM7_HHZZ1_T -> "damhz1_t_7(gkm,mkm,"
      | D_FM7_HHZZ1_U -> "damhz1_u_7(gkm,mkm,"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_SWW -> "gsww" | G_SZZ -> "gszz"
      | G_SHH -> "gshh"
      | G_SWW_T -> "gswwt" | G_SZZ_T -> "gszzt"
      | G_SAA_T -> "gsaat" | G_SAZ_T -> "gsazt"
      | G_PNWW -> "gpnww" | G_PNZZ -> "gpnzz"
      | G_PSNWW -> "gpsnww" | G_PSNZZ -> "gpsnzz"
      | G_PSNHH -> "gpsnhh"
      | G_PWZ -> "gpwz" | G_PWW -> "gpww"
      | G_FWW -> "gfww" | G_FZZ -> "gfzz"
      | G_FWW_CF -> "gfwwcf" | G_FZZ_CF -> "gfzzcf"
      | G_FHH -> "gfhh" | G_FHH_CF -> "gfhhcf"
      | G_FWW_T -> "gfwwt" | G_FZZ_T -> "gfzzt"
      | G_TNWW -> "gtnww" | G_TNZZ -> "gtnzz"
      | G_TNWW_CF -> "gtnwwcf" | G_TNZZ_CF -> "gtnzzcf"
      | G_TSNWW -> "gtsnww" | G_TSNZZ -> "gtsnzz"
      | G_TSNWW_CF -> "gtsnwwcf" | G_TSNZZ_CF -> "gtsnzzcf"
      | G_TWZ -> "gtwz" | G_TWW -> "gtww"
      | G_TWZ_CF -> "gtwzcf" | G_TWW_CF -> "gtwwcf"
      | G_SSWW -> "gssww" | G_SSZZ -> "gsszz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc" | G_Hmm -> "ghmm"
      | G_HGaZ -> "ghgaz" | G_HGaGa -> "ghgaga" | G_Hgg -> "ghgg"
      | G_HGaGa_anom -> "ghgaga_ac" | G_HGaZ_anom -> "ghgaz_ac"
      | G_HZZ_anom -> "ghzz_ac" | G_HWW_anom -> "ghww_ac"
      | G_HGaZ_u -> "ghgaz_u" | G_HZZ_u -> "ghzz_u" 
      | G_HWW_u -> "ghww_u" 
      | G_H3 -> "gh3" | G_H4 -> "gh4"
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f
      | K_Matrix_Coeff i -> "kc" ^ string_of_int i
      | K_Matrix_Pole i -> "kp" ^ string_of_int i

  end

(* \thocwmodulesection{Complete Minimal Standard Model including additional Resonances (alternate Tensor)} *)

module SSC_AltT (Flags : SSC_flags) =
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
                     QGUG | QBUB | QW | DL | DR

    type matter_field = L of int | N of int | U of int | D of int
    type gauge_boson = Ga | Wp | Wm | Z | Gl
    type other = Phip | Phim | Phi0 | H
                 | Rsigma | Rphin | Rphisn | Rphip | Rphim | Rphipp | Rphimm 
                 | Rf | Rtn | Rtsn | Rtp | Rtm | Rtpp | Rtmm
                 | Rff | Rfv | Rfphi
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
      failwith "Modellib_BSM.SSC_AltT.gauge_symbol: internal error"

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
        "Higgs", List.map other [H];
	"Scalar Resonances", List.map other [Rsigma; Rphin; Rphisn; Rphip; Rphim; Rphipp; Rphimm];
	"Tensor Resonances", List.map other [Rf; Rtn; Rtsn; Rtp; Rtm; Rtpp; Rtmm];
	"Alternate Tensor", List.map other [Rff; Rfv; Rfphi];
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
          | Rf | Rtn | Rtsn | Rtp | Rtm | Rtpp | Rtmm -> Tensor_2
          | Rff -> Tensor_2
          | Rfv -> Vector
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
          | H | Rsigma -> Prop_Scalar
          | Rphin | Rphisn  | Rphip | Rphim | Rphipp | Rphimm -> Prop_Scalar
          | Rf -> Prop_Tensor_2
	  | Rff -> Prop_Tensor_pure
	  | Rfv -> Prop_Vector_pure
          | Rfphi -> Prop_Scalar
          | Rtn | Rtsn | Rtp | Rtm | Rtpp | Rtmm -> Prop_Tensor_2
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
          | H -> H |  Rsigma -> Rsigma
          | Rphin -> Rphin | Rphip -> Rphim | Rphim -> Rphip
          | Rphisn -> Rphisn
          | Rphipp -> Rphimm | Rphimm -> Rphipp
          | Rf -> Rf
          | Rff -> Rff | Rfv -> Rfv | Rfphi -> Rfphi
          | Rtn -> Rtn | Rtsn -> Rtsn | Rtp -> Rtm | Rtm -> Rtp
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
          | H | Rsigma | Phi0 | Rphin | Rphisn 
          | Rf | Rff | Rfv | Rfphi | Rtn | Rtsn ->  0//1
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
      | G_SP_ttH
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
      | FS0_HHWW | FS0_HHZZ
      | FS1_HHWW | FS1_HHZZ
      | FM0_HHWW | FM0_HHZZ 
      | FM1_HHWW | FM1_HHZZ   
      | FM7_HHWW | FM7_HHZZ
      | Alpha_WWWW0 | Alpha_ZZWW1 | Alpha_WWWW2
      | Alpha_ZZWW0 | Alpha_ZZZZ
      | FT0_WWWW0 | FT0_WWWW2
      | FT0_ZZWW0 | FT0_ZZWW1
      | FT0_ZZZZ  | FT0_AAAA
      | FT0_AAWW0 | FT0_AAWW1
      | FT0_AAZZ  
      | FT0_AZWW0 | FT0_AZWW1
      | FT0_AAAZ  | FT0_AZZZ
      | FT1_WWWW0 | FT1_WWWW2
      | FT1_ZZWW0 | FT1_ZZWW1
      | FT1_ZZZZ  | FT1_AAAA  
      | FT1_AAWW0 | FT1_AAWW1
      | FT1_AAZZ  
      | FT1_AZWW0 | FT1_AZWW1 
      | FT1_AAAZ  | FT1_AZZZ
      | FT2_WWWW0 | FT2_WWWW2
      | FT2_ZZWW0 | FT2_ZZWW1
      | FT2_ZZZZ  | FT2_AAAA
      | FT2_AAWW0 | FT2_AAWW1
      | FT2_AAZZ  
      | FT2_AZWW0 | FT2_AZWW1 
      | FT2_AAAZ  | FT2_AZZZ
      | FM0_WWWW0 | FM0_WWWW2
      | FM0_ZZWW0 | FM0_ZZWW1
      | FM0_ZZZZ  
      | FM1_WWWW0 | FM1_WWWW2
      | FM1_ZZWW0 | FM1_ZZWW1
      | FM1_ZZZZ     
      | FM7_WWWW0 | FM7_WWWW2
      | FM7_ZZWW0 | FM7_ZZWW1
      | FM7_ZZZZ
      | D_Alpha_ZZWW0_S | D_Alpha_ZZWW0_T | D_Alpha_ZZWW1_S
      | D_Alpha_ZZWW1_T | D_Alpha_ZZWW1_U | D_Alpha_WWWW0_S
      | D_Alpha_WWWW0_T | D_Alpha_WWWW0_U | D_Alpha_WWWW2_S
      | D_Alpha_WWWW2_T | D_Alpha_ZZZZ_S | D_Alpha_ZZZZ_T
      | D_FT0_ZZWW0_S | D_FT0_ZZWW0_T | D_FT0_ZZWW0_U 
      | D_FT0_ZZWW1_S | D_FT0_ZZWW1_T | D_FT0_ZZWW1_U 
      | D_FT0_WWWW0_S | D_FT0_WWWW0_T | D_FT0_WWWW0_U 
      | D_FT0_WWWW2_S | D_FT0_WWWW2_T | D_FT0_WWWW2_U    
      | D_FT0_ZZZZ_S  | D_FT0_ZZZZ_T  | D_FT0_ZZZZ_U
      | D_FT0_AAAA_S  | D_FT0_AAAA_T  | D_FT0_AAAA_U  
      | D_FT0_AAWW0_S | D_FT0_AAWW0_T | D_FT0_AAWW0_U 
      | D_FT0_AAWW1_S | D_FT0_AAWW1_T | D_FT0_AAWW1_U
      | D_FT0_AAZZ_S  | D_FT0_AAZZ_T  | D_FT0_AAZZ_U 
      | D_FT0_AZWW0_S | D_FT0_AZWW0_T | D_FT0_AZWW0_U
      | D_FT0_AZWW1_S | D_FT0_AZWW1_T | D_FT0_AZWW1_U 
      | D_FT0_AAAZ_S  | D_FT0_AAAZ_T  | D_FT0_AAAZ_U
      | D_FT0_AZZZ_S  | D_FT0_AZZZ_T  | D_FT0_AZZZ_U           
      | D_FT1_ZZWW0_S | D_FT1_ZZWW0_T | D_FT1_ZZWW0_U 
      | D_FT1_ZZWW1_S | D_FT1_ZZWW1_T | D_FT1_ZZWW1_U 
      | D_FT1_WWWW0_S | D_FT1_WWWW0_T | D_FT1_WWWW0_U 
      | D_FT1_WWWW2_S | D_FT1_WWWW2_T | D_FT1_WWWW2_U    
      | D_FT1_ZZZZ_S  | D_FT1_ZZZZ_T  | D_FT1_ZZZZ_U 
      | D_FT1_AAAA_S  | D_FT1_AAAA_T  | D_FT1_AAAA_U  
      | D_FT1_AAWW0_S | D_FT1_AAWW0_T | D_FT1_AAWW0_U 
      | D_FT1_AAWW1_S | D_FT1_AAWW1_T | D_FT1_AAWW1_U
      | D_FT1_AAZZ_S  | D_FT1_AAZZ_T  | D_FT1_AAZZ_U  
      | D_FT1_AZWW0_S | D_FT1_AZWW0_T | D_FT1_AZWW0_U
      | D_FT1_AZWW1_S | D_FT1_AZWW1_T | D_FT1_AZWW1_U
      | D_FT1_AAAZ_S  | D_FT1_AAAZ_T  | D_FT1_AAAZ_U
      | D_FT1_AZZZ_S  | D_FT1_AZZZ_T  | D_FT1_AZZZ_U       
      | D_FT2_ZZWW0_S | D_FT2_ZZWW0_T | D_FT2_ZZWW0_U 
      | D_FT2_ZZWW1_S | D_FT2_ZZWW1_T | D_FT2_ZZWW1_U 
      | D_FT2_WWWW0_S | D_FT2_WWWW0_T | D_FT2_WWWW0_U 
      | D_FT2_WWWW2_S | D_FT2_WWWW2_T | D_FT2_WWWW2_U    
      | D_FT2_ZZZZ_S  | D_FT2_ZZZZ_T  | D_FT2_ZZZZ_U 
      | D_FT2_AAAA_S  | D_FT2_AAAA_T  | D_FT2_AAAA_U  
      | D_FT2_AAWW0_S | D_FT2_AAWW0_T | D_FT2_AAWW0_U 
      | D_FT2_AAWW1_S | D_FT2_AAWW1_T | D_FT2_AAWW1_U
      | D_FT2_AAZZ_S  | D_FT2_AAZZ_T  | D_FT2_AAZZ_U  
      | D_FT2_AZWW0_S | D_FT2_AZWW0_T | D_FT2_AZWW0_U
      | D_FT2_AZWW1_S | D_FT2_AZWW1_T | D_FT2_AZWW1_U  
      | D_FT2_AAAZ_S  | D_FT2_AAAZ_T  | D_FT2_AAAZ_U
      | D_FT2_AZZZ_S  | D_FT2_AZZZ_T  | D_FT2_AZZZ_U  
      | D_FTrsi_ZZWW0_S | D_FTrsi_ZZWW0_T | D_FTrsi_ZZWW0_U 
      | D_FTrsi_ZZWW1_S | D_FTrsi_ZZWW1_T | D_FTrsi_ZZWW1_U 
      | D_FTrsi_WWWW0_S | D_FTrsi_WWWW0_T | D_FTrsi_WWWW0_U 
      | D_FTrsi_WWWW2_S | D_FTrsi_WWWW2_T | D_FTrsi_WWWW2_U    
      | D_FTrsi_ZZZZ_S  | D_FTrsi_ZZZZ_T  | D_FTrsi_ZZZZ_U 
      | D_FTrsi_AAAA_S  | D_FTrsi_AAAA_T  | D_FTrsi_AAAA_U
      | D_FTrsi_AAWW0_S | D_FTrsi_AAWW0_T | D_FTrsi_AAWW0_U 
      | D_FTrsi_AAWW1_S | D_FTrsi_AAWW1_T | D_FTrsi_AAWW1_U
      | D_FTrsi_AAZZ_S  | D_FTrsi_AAZZ_T  | D_FTrsi_AAZZ_U 
      | D_FTrsi_AZWW0_S | D_FTrsi_AZWW0_T | D_FTrsi_AZWW0_U
      | D_FTrsi_AZWW1_S | D_FTrsi_AZWW1_T | D_FTrsi_AZWW1_U 
      | D_FTrsi_AAAZ_S  | D_FTrsi_AAAZ_T  | D_FTrsi_AAAZ_U
      | D_FTrsi_AZZZ_S  | D_FTrsi_AZZZ_T  | D_FTrsi_AZZZ_U        
      | D_FM0_ZZWW0_S | D_FM0_ZZWW0_T | D_FM0_ZZWW0_U 
      | D_FM0_ZZWW1_S | D_FM0_ZZWW1_T | D_FM0_ZZWW1_U 
      | D_FM0_WWWW0_S | D_FM0_WWWW0_T | D_FM0_WWWW0_U 
      | D_FM0_WWWW2_S | D_FM0_WWWW2_T | D_FM0_WWWW2_U    
      | D_FM0_ZZZZ_S | D_FM0_ZZZZ_T | D_FM0_ZZZZ_U       
      | D_FM1_ZZWW0_S | D_FM1_ZZWW0_T | D_FM1_ZZWW0_U 
      | D_FM1_ZZWW1_S | D_FM1_ZZWW1_T | D_FM1_ZZWW1_U 
      | D_FM1_WWWW0_S | D_FM1_WWWW0_T | D_FM1_WWWW0_U 
      | D_FM1_WWWW2_S | D_FM1_WWWW2_T | D_FM1_WWWW2_U    
      | D_FM1_ZZZZ_S | D_FM1_ZZZZ_T | D_FM1_ZZZZ_U 
      | D_FM7_ZZWW0_S | D_FM7_ZZWW0_T | D_FM7_ZZWW0_U 
      | D_FM7_ZZWW1_S | D_FM7_ZZWW1_T | D_FM7_ZZWW1_U 
      | D_FM7_WWWW0_S | D_FM7_WWWW0_T | D_FM7_WWWW0_U 
      | D_FM7_WWWW2_S | D_FM7_WWWW2_T | D_FM7_WWWW2_U    
      | D_FM7_ZZZZ_S | D_FM7_ZZZZ_T | D_FM7_ZZZZ_U
      | D_Alpha_HHHH_S  | D_Alpha_HHHH_T 
      | D_Alpha_HHZZ0_S | D_Alpha_HHWW0_S 
      | D_Alpha_HHZZ0_T | D_Alpha_HHWW0_T
      | D_Alpha_HHZZ1_S | D_Alpha_HHWW1_S 
      | D_Alpha_HHZZ1_T | D_Alpha_HHWW1_T
      | D_Alpha_HHZZ1_U | D_Alpha_HHWW1_U
      | D_FM0_HHZZ0_S | D_FM0_HHWW0_S 
      | D_FM0_HHZZ0_T | D_FM0_HHWW0_T
      | D_FM0_HHZZ0_U | D_FM0_HHWW0_U      
      | D_FM0_HHZZ1_S | D_FM0_HHWW1_S 
      | D_FM0_HHZZ1_T | D_FM0_HHWW1_T
      | D_FM0_HHZZ1_U | D_FM0_HHWW1_U   
      | D_FM1_HHZZ0_S | D_FM1_HHWW0_S 
      | D_FM1_HHZZ0_T | D_FM1_HHWW0_T
      | D_FM1_HHZZ0_U | D_FM1_HHWW0_U      
      | D_FM1_HHZZ1_S | D_FM1_HHWW1_S 
      | D_FM1_HHZZ1_T | D_FM1_HHWW1_T
      | D_FM1_HHZZ1_U | D_FM1_HHWW1_U
      | D_FM7_HHZZ0_S | D_FM7_HHWW0_S 
      | D_FM7_HHZZ0_T | D_FM7_HHWW0_T
      | D_FM7_HHZZ0_U | D_FM7_HHWW0_U      
      | D_FM7_HHZZ1_S | D_FM7_HHWW1_S 
      | D_FM7_HHZZ1_T | D_FM7_HHWW1_T
      | D_FM7_HHZZ1_U | D_FM7_HHWW1_U
      | G_HWW | G_HHWW | G_HZZ | G_HHZZ
      | G_SWW | G_SWW_T | G_SSWW | G_SZZ 
      | G_SZZ_T | G_SSZZ | G_SHH
      | G_SAA_T | G_SAZ_T 
      | G_PNWW | G_PNZZ | G_PWZ | G_PWW
      | G_PSNWW | G_PSNZZ | G_PSNHH
      | G_FWW | G_FZZ | G_FWW_CF | G_FZZ_CF 
      | G_FWW_T | G_FZZ_T | G_FHH | G_FHH_CF
      | G_FFWW | G_FFZZ | G_FFWW_CF | G_FFZZ_CF 
      | G_FFHH | G_FFHH_CF
      | G_FVWW | G_FVZZ | G_FVHH | G_FVWW_CF | G_FVZZ_CF | G_FVHH_CF
      | G_FDDSWW | G_FDDSZZ | G_FDDSHH | G_FDDSWW_CF 
      | G_FDDSZZ_CF | G_FDDSHH_CF | G_FSWW | G_FSZZ | G_FSHH 
      | G_TNWW | G_TNZZ | G_TSNWW | G_TSNZZ | G_TWZ | G_TWW
      | G_TNWW_CF | G_TNZZ_CF | G_TSNWW_CF | G_TSNZZ_CF
      | G_TWZ_CF | G_TWW_CF
      | G_Htt | G_Hbb | G_Hcc | G_Hmm | G_Htautau | G_H3 | G_H4
      | FS_H4
      | G_HGaZ | G_HGaGa | G_Hgg
      | G_HGaZ_anom | G_HGaGa_anom | G_HZZ_anom | G_HWW_anom  
      | G_HGaZ_u | G_HZZ_u | G_HWW_u
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

(* \begin{dubious}
   \ldots{} to be continued \ldots{}
   The quartic couplings can't be correct, because the dimensions are wrong!
   \begin{subequations}
     \begin{align}
                  g_{HWW} &= g m_W = 2 \frac{m_W^2}{v}\\
                 g_{HHWW} &= 2 \frac{m_W^2}{v^2} = \frac{g^2}{2} \\
                  g_{HZZ} &= \frac{g}{\cos\theta_w}m_Z \\
                 g_{HHZZ} &= 2 \frac{m_Z^2}{v^2} = \frac{g^2}{2\cos\theta_w} \\
                  g_{Htt} &= \lambda_t \\
                  g_{Hbb} &= \lambda_b=\frac{m_b}{m_t}\lambda_t \\
                  g_{H^3} &= - \frac{3g}{2}\frac{m_H^2}{m_W} = - 3 \frac{m_H^2}{v} 
                  g_{H^4} &= - \frac{3g^2}{4} \frac{m_W^2}{v^2} = -3 \frac{m_H^2}{v^2}  
     \end{align}
   \end{subequations}
   \end{dubious} *)

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

    let yukawa =
      [ ((M (U (-3)), O H, M (U 3)), FBF (1, Psibar, S, Psi), G_Htt);
        ((M (D (-3)), O H, M (D 3)), FBF (1, Psibar, S, Psi), G_Hbb);
        ((M (U (-2)), O H, M (U 2)), FBF (1, Psibar, S, Psi), G_Hcc);
        ((M (L (-3)), O H, M (L 3)), FBF (1, Psibar, S, Psi), G_Htautau) ] @
      if Flags.higgs_hmm then
      [ ((M (L (-2)), O H, M (L 2)), FBF (1, Psibar, S, Psi), G_Hmm)]
          else
      []

      
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
            ((Z, Z, Z, Z),
             Vector4 [(1, C_12_34); (1, C_13_42); (1, C_14_23)], Alpha_ZZZZ);
            ((Wm, Wp, Z, Z),
             Vector4 [1, C_12_34], Alpha_ZZWW0);
            ((Wm, Wp, Z, Z),
             Vector4 [(1, C_13_42); (1, C_14_23)], Alpha_ZZWW1)]	    
	@
	  (if Flags.k_matrix_tm then
	      List.map qgc
           [((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_WWWW2);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_WWWW2);             
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_WWWW2); 
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_0 [1, C_13_42], FM0_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_0 [1, C_14_23], FM0_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_0 [1, C_12_34], FM0_WWWW2); 
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_1 [1, C_13_42], FM1_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_1 [1, C_14_23], FM1_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_1 [1, C_12_34], FM1_WWWW2);  
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_7 [1, C_13_42], FM7_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_7 [1, C_14_23], FM7_WWWW0);
            ((Wm, Wm, Wp, Wp),
             Dim8_Vector4_m_7 [1, C_12_34], FM7_WWWW2);             
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_ZZWW1); 
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_ZZWW1);              
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_ZZWW1);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_0 [1, C_12_34], FM0_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_0 [1, C_13_42], FM0_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_0 [1, C_14_23], FM0_ZZWW1);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_1 [1, C_12_34], FM1_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_1 [1, C_13_42], FM1_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_1 [1, C_14_23], FM1_ZZWW1); 
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_7 [1, C_12_34], FM7_ZZWW0);   
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_7 [1, C_13_42], FM7_ZZWW1);
            ((Wm, Wp, Z, Z),
             Dim8_Vector4_m_7 [1, C_14_23], FM7_ZZWW1);              
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_ZZZZ);
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_ZZZZ);             
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_0 [1, C_12_34], FM0_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_0 [1, C_13_42], FM0_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_0 [1, C_14_23], FM0_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_1 [1, C_12_34], FM1_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_1 [1, C_13_42], FM1_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_1 [1, C_14_23], FM1_ZZZZ);
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_7 [1, C_12_34], FM7_ZZZZ);  
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_7 [1, C_13_42], FM7_ZZZZ); 
            ((Z, Z, Z, Z),
             Dim8_Vector4_m_7 [1, C_14_23], FM7_ZZZZ);
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AAAA);  
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AAAA); 
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AAAA);
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AAAA);  
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AAAA); 
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AAAA);             
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AAAA);  
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AAAA); 
            ((Ga, Ga, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AAAA);
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AAWW0);   
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AAWW1);
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AAWW1); 
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AAWW0);   
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AAWW1);
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AAWW1);              
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AAWW0);   
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AAWW1);
            ((Wm, Wp, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AAWW1);
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AAZZ);  
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AAZZ); 
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AAZZ);
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AAZZ);  
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AAZZ); 
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AAZZ);             
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AAZZ);  
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AAZZ); 
            ((Z, Z, Ga, Ga),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AAZZ);
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AZWW0);  
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AZWW1); 
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AZWW1);
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AZWW0);  
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AZWW1); 
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AZWW1);             
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AZWW0);  
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AZWW1); 
            ((Ga, Z, Wp, Wm),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AZWW1);
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AAAZ);  
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AAAZ); 
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AAAZ);
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AAAZ);  
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AAAZ); 
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AAAZ);             
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AAAZ);  
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AAAZ); 
            ((Ga, Ga, Ga, Z),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AAAZ);
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_12_34], FT0_AZZZ);  
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_13_42], FT0_AZZZ); 
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_0 [1, C_14_23], FT0_AZZZ);
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_12_34], FT1_AZZZ);  
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_13_42], FT1_AZZZ); 
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_1 [1, C_14_23], FT1_AZZZ);             
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_12_34], FT2_AZZZ);  
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_13_42], FT2_AZZZ); 
            ((Ga, Z, Z, Z),
             Dim8_Vector4_t_2 [1, C_14_23], FT2_AZZZ)]
      else
        [] )
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
        
    let k_matrix_quartic_gauge_t_0 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_12_34)]), D_FT0_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_13_42)]), D_FT0_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_14_23)]), D_FT0_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_ZZZZ_U);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAAA_T); 
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAAA_U);        
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AAAA_T);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AAAA_U);                   
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAWW0_S);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAWW0_T);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAWW0_U);                   
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAWW1_S);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAWW1_T);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AAWW1_S);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AAWW1_T); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_12_34)]), D_FT0_AAWW1_S); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_13_42)]), D_FT0_AAWW1_U);
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t0 (2,
                   [(1, C_14_23)]), D_FT0_AAWW1_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AAZZ_T); 
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AAZZ_U);        
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AAZZ_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AAZZ_U); 
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_12_34)]), D_FT0_AZWW0_S);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_13_42)]), D_FT0_AZWW0_T);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t0 (0,
                   [(1, C_14_23)]), D_FT0_AZWW0_U); 
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AZWW1_S);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AZWW1_T);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AZWW1_U); 
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AZWW1_S);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AZWW1_T);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AZWW1_U); 
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AZWW1_S);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AZWW1_T);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AZWW1_U); 
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_12_34)]), D_FT0_AZWW1_S);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_13_42)]), D_FT0_AZWW1_T);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t0 (1,
                   [(1, C_14_23)]), D_FT0_AZWW1_U); 
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AAAZ_S);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AAAZ_T);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AAAZ_U); 
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AAAZ_S);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AAAZ_T);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AAAZ_U); 
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AAAZ_S);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AAAZ_T);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AAAZ_U); 
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AZZZ_S);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AZZZ_T);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AZZZ_U); 
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AZZZ_S);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AZZZ_T);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AZZZ_U); 
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_12_34)]), D_FT0_AZZZ_S);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_13_42)]), D_FT0_AZZZ_T);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t0 (3,
                   [(1, C_14_23)]), D_FT0_AZZZ_U)]        
      else
        []
        
    let k_matrix_quartic_gauge_t_1 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_12_34)]), D_FT1_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_13_42)]), D_FT1_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_14_23)]), D_FT1_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_ZZZZ_U);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1(0,
                   [(1, C_12_34)]), D_FT1_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AAAA_T); 
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AAAA_U);        
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAAA_T);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAAA_U);                   
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_AAWW0_S);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AAWW0_T);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AAWW0_U);                   
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_AAWW1_S);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AAWW1_T);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AAWW1_S);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AAWW1_T); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_12_34)]), D_FT1_AAWW1_S); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_13_42)]), D_FT1_AAWW1_U);
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t1 (2,
                   [(1, C_14_23)]), D_FT1_AAWW1_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AAZZ_T); 
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AAZZ_U);        
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAZZ_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAZZ_U); 
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_12_34)]), D_FT1_AZWW0_S);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_13_42)]), D_FT1_AZWW0_T);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t1 (0,
                   [(1, C_14_23)]), D_FT1_AZWW0_U); 
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AZWW1_S);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AZWW1_T);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AZWW1_U); 
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AZWW1_S);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AZWW1_T);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AZWW1_U); 
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AZWW1_S);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AZWW1_T);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AZWW1_U); 
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_12_34)]), D_FT1_AZWW1_S);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_13_42)]), D_FT1_AZWW1_T);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t1 (1,
                   [(1, C_14_23)]), D_FT1_AZWW1_U); 
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAAZ_S);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAAZ_T);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAAZ_U); 
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAAZ_S);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAAZ_T);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAAZ_U); 
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AAAZ_S);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AAAZ_T);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AAAZ_U); 
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AZZZ_S);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AZZZ_T);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AZZZ_U); 
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AZZZ_S);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AZZZ_T);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AZZZ_U); 
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_12_34)]), D_FT1_AZZZ_S);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_13_42)]), D_FT1_AZZZ_T);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t1 (3,
                   [(1, C_14_23)]), D_FT1_AZZZ_U)]        
      else
        []        
        
    let k_matrix_quartic_gauge_t_2 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_12_34)]), D_FT2_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_13_42)]), D_FT2_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_14_23)]), D_FT2_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_ZZZZ_U);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AAAA_T); 
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AAAA_U);        
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAAA_T);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAAA_U);                   
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AAWW0_S);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AAWW0_T);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AAWW0_U);                   
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AAWW1_S);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AAWW1_T);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AAWW1_S);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AAWW1_T); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_12_34)]), D_FT2_AAWW1_S); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_13_42)]), D_FT2_AAWW1_U);
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t2 (2,
                   [(1, C_14_23)]), D_FT2_AAWW1_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AAZZ_T); 
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AAZZ_U);        
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAZZ_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAZZ_U); 
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_12_34)]), D_FT2_AZWW0_S);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_13_42)]), D_FT2_AZWW0_T);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t2 (0,
                   [(1, C_14_23)]), D_FT2_AZWW0_U); 
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AZWW1_S);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AZWW1_T);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AZWW1_U); 
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AZWW1_S);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AZWW1_T);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AZWW1_U); 
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AZWW1_S);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AZWW1_T);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AZWW1_U); 
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_12_34)]), D_FT2_AZWW1_S);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_13_42)]), D_FT2_AZWW1_T);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t2 (1,
                   [(1, C_14_23)]), D_FT2_AZWW1_U); 
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAAZ_S);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAAZ_T);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAAZ_U); 
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAAZ_S);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAAZ_T);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAAZ_U); 
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AAAZ_S);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AAAZ_T);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AAAZ_U); 
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AZZZ_S);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AZZZ_T);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AZZZ_U); 
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AZZZ_S);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AZZZ_T);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AZZZ_U); 
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_12_34)]), D_FT2_AZZZ_S);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_13_42)]), D_FT2_AZZZ_T);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t2 (3,
                   [(1, C_14_23)]), D_FT2_AZZZ_U)]        
      else
        []
        
    let k_matrix_quartic_gauge_t_rsi =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_12_34)]), D_FTrsi_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_13_42)]), D_FTrsi_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_14_23)]), D_FTrsi_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_ZZZZ_U);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAAA_T); 
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAAA_U);        
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AAAA_S);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AAAA_T);
            ((Ga, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AAAA_U);                   
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAWW0_S);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAWW0_T);
            ((Wm, Wp, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAWW0_U);                   
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAWW1_S);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAWW1_T);
            ((Wm, Ga, Wp, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AAWW1_S);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AAWW1_U);
            ((Wp, Ga, Ga, Wm), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AAWW1_T); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_12_34)]), D_FTrsi_AAWW1_S); 
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_13_42)]), D_FTrsi_AAWW1_U);
            ((Ga, Wp, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (2,
                   [(1, C_14_23)]), D_FTrsi_AAWW1_T);
            ((Ga, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAZZ_S);
            ((Ga, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAZZ_T);
            ((Ga, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAZZ_U);
            ((Z, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAZZ_S); 
            ((Z, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAZZ_T);
            ((Z, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAZZ_U); 
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AAZZ_S);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AAZZ_T);
            ((Ga, Ga, Z, Z), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AAZZ_U); 
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_12_34)]), D_FTrsi_AZWW0_S);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_13_42)]), D_FTrsi_AZWW0_T);
            ((Ga, Z, Wp, Wm), Vector4_K_Matrix_cf_t_rsi (0,
                   [(1, C_14_23)]), D_FTrsi_AZWW0_U); 
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AZWW1_S);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AZWW1_T);
            ((Wp, Ga, Wm, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AZWW1_U); 
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AZWW1_S);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AZWW1_T);
            ((Wm, Ga, Wp, Z), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AZWW1_U); 
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AZWW1_S);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AZWW1_T);
            ((Z, Wm, Ga, Wp), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AZWW1_U); 
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_12_34)]), D_FTrsi_AZWW1_S);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_13_42)]), D_FTrsi_AZWW1_T);
            ((Wp, Z, Wm, Ga), Vector4_K_Matrix_cf_t_rsi (1,
                   [(1, C_14_23)]), D_FTrsi_AZWW1_U); 
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AAAZ_S);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AAAZ_T);
            ((Ga, Ga, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AAAZ_U); 
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AAAZ_S);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AAAZ_T);
            ((Z, Ga, Ga, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AAAZ_U); 
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AAAZ_S);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AAAZ_T);
            ((Ga, Ga, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AAAZ_U); 
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AZZZ_S);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AZZZ_T);
            ((Z, Z, Z, Ga), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AZZZ_U); 
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AZZZ_S);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AZZZ_T);
            ((Ga, Z, Z, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AZZZ_U); 
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_12_34)]), D_FTrsi_AZZZ_S);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_13_42)]), D_FTrsi_AZZZ_T);
            ((Z, Z, Ga, Z), Vector4_K_Matrix_cf_t_rsi (3,
                   [(1, C_14_23)]), D_FTrsi_AZZZ_U)]        
      else
        []        
        
    let k_matrix_quartic_gauge_m_0 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m0 (1,
                   [(1, C_12_34)]), D_FM0_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m0 (1,
                   [(1, C_13_42)]), D_FM0_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m0 (1,
                   [(1, C_14_23)]), D_FM0_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m0 (2,
                   [(1, C_12_34)]), D_FM0_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m0 (2,
                   [(1, C_13_42)]), D_FM0_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m0 (2,
                   [(1, C_14_23)]), D_FM0_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_12_34)]), D_FM0_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_13_42)]), D_FM0_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (0,
                   [(1, C_14_23)]), D_FM0_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (3,
                   [(1, C_14_23)]), D_FM0_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (3,
                   [(1, C_13_42)]), D_FM0_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m0 (3,
                   [(1, C_12_34)]), D_FM0_ZZZZ_U)]        
      else
        []

    let k_matrix_quartic_gauge_m_1 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m1 (1,
                   [(1, C_12_34)]), D_FM1_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m1 (1,
                   [(1, C_13_42)]), D_FM1_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m1 (1,
                   [(1, C_14_23)]), D_FM1_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m1 (2,
                   [(1, C_12_34)]), D_FM1_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m1 (2,
                   [(1, C_13_42)]), D_FM1_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m1 (2,
                   [(1, C_14_23)]), D_FM1_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_12_34)]), D_FM1_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_13_42)]), D_FM1_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (0,
                   [(1, C_14_23)]), D_FM1_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (3,
                   [(1, C_14_23)]), D_FM1_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (3,
                   [(1, C_13_42)]), D_FM1_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m1 (3,
                   [(1, C_12_34)]), D_FM1_ZZZZ_U)]        
      else
        []
        
    let k_matrix_quartic_gauge_m_7 =
      if Flags.k_matrix_tm then
        List.map qgc
          [ ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_WWWW0_S);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_WWWW0_T);
            ((Wm, Wp, Wm, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_WWWW0_U);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_WWWW0_S); 
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_WWWW0_T);
            ((Wp, Wm, Wp, Wm), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_WWWW0_U); 
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_WWWW2_S);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_WWWW2_T);
            ((Wm, Wm, Wp, Wp), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_WWWW2_U);                   
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_ZZWW0_S);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_ZZWW0_T);
            ((Wm, Wp, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_ZZWW0_U);                   
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_ZZWW1_S);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_ZZWW1_T);
            ((Wm, Z, Wp, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m7 (1,
                   [(1, C_12_34)]), D_FM7_ZZWW1_S);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m7 (1,
                   [(1, C_13_42)]), D_FM7_ZZWW1_U);
            ((Wp, Z, Z, Wm), Vector4_K_Matrix_cf_m7 (1,
                   [(1, C_14_23)]), D_FM7_ZZWW1_T); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m7 (2,
                   [(1, C_12_34)]), D_FM7_ZZWW1_S); 
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m7 (2,
                   [(1, C_13_42)]), D_FM7_ZZWW1_U);
            ((Z, Wp, Wm, Z), Vector4_K_Matrix_cf_m7 (2,
                   [(1, C_14_23)]), D_FM7_ZZWW1_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_12_34)]), D_FM7_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_13_42)]), D_FM7_ZZZZ_T); 
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (0,
                   [(1, C_14_23)]), D_FM7_ZZZZ_U);        
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (3,
                   [(1, C_14_23)]), D_FM7_ZZZZ_S);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (3,
                   [(1, C_13_42)]), D_FM7_ZZZZ_T);
            ((Z, Z, Z, Z), Vector4_K_Matrix_cf_m7 (3,
                   [(1, C_12_34)]), D_FM7_ZZZZ_U)]        
      else
        []    

    let k_matrix_2scalar_2gauge =
      if Flags.k_matrix_tm then
        if Flags.higgs_matrix then
            [ ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms 
                   (0,  [(1, C_12_34)]), D_Alpha_HHZZ0_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (0,  [(1, C_13_42)]), D_Alpha_HHZZ0_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (0,  [(1, C_14_23)]), D_Alpha_HHZZ0_T);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms 
                   (3,  [(1, C_14_23)]), D_Alpha_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (3,  [(1, C_13_42)]), D_Alpha_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (3,  [(1, C_12_34)]), D_Alpha_HHZZ1_U);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms 
                   (6,  [(1, C_13_42)]), D_Alpha_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (6,  [(1, C_12_34)]), D_Alpha_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_K_Matrix_ms
                   (6,  [(1, C_14_23)]), D_Alpha_HHZZ1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms 
                   (0,  [(1, C_12_34)]), D_Alpha_HHWW0_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (2,  [(1, C_13_42)]), D_Alpha_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (1,  [(1, C_14_23)]), D_Alpha_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (1,  [(1, C_13_42)]), D_Alpha_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (2,  [(1, C_14_23)]), D_Alpha_HHWW0_T);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms 
                   (3,  [(1, C_14_23)]), D_Alpha_HHWW1_S);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms 
                   (6,  [(1, C_13_42)]), D_Alpha_HHWW1_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (4,  [(1, C_13_42)]), D_Alpha_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (5,  [(1, C_12_34)]), D_Alpha_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (8,  [(1, C_14_23)]), D_Alpha_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (7,  [(1, C_12_34)]), D_Alpha_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (5,  [(1, C_13_42)]), D_Alpha_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (4,  [(1, C_12_34)]), D_Alpha_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (7,  [(1, C_14_23)]), D_Alpha_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_K_Matrix_ms
                   (8,  [(1, C_12_34)]), D_Alpha_HHWW1_U) ]
        else
            []
      else
          []
          
    let k_matrix_2scalar_2gauge_m =
      if Flags.k_matrix_tm then
        if Flags.higgs_matrix then
            [ ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM0_HHZZ0_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (0,  [(1, C_13_42)]), D_FM0_HHZZ0_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (0,  [(1, C_14_23)]), D_FM0_HHZZ0_U);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf 
                   (3,  [(1, C_14_23)]), D_FM0_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (3,  [(1, C_13_42)]), D_FM0_HHZZ1_U);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (3,  [(1, C_12_34)]), D_FM0_HHZZ1_T);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM0_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (6,  [(1, C_12_34)]), D_FM0_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_0_K_Matrix_cf
                   (6,  [(1, C_14_23)]), D_FM0_HHZZ1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM0_HHWW0_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (2,  [(1, C_13_42)]), D_FM0_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (1,  [(1, C_14_23)]), D_FM0_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (1,  [(1, C_13_42)]), D_FM0_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (2,  [(1, C_14_23)]), D_FM0_HHWW0_T);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (3,  [(1, C_14_23)]), D_FM0_HHWW1_S);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM0_HHWW1_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (4,  [(1, C_13_42)]), D_FM0_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (5,  [(1, C_12_34)]), D_FM0_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (8,  [(1, C_14_23)]), D_FM0_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (7,  [(1, C_12_34)]), D_FM0_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (5,  [(1, C_13_42)]), D_FM0_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (4,  [(1, C_12_34)]), D_FM0_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (7,  [(1, C_14_23)]), D_FM0_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_0_K_Matrix_cf
                   (8,  [(1, C_12_34)]), D_FM0_HHWW1_T) ]
        else
            []
      else
          [] 
          
    let k_matrix_2scalar_2gauge_m_1 =
      if Flags.k_matrix_tm then
        if Flags.higgs_matrix then
            [ ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM1_HHZZ0_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (0,  [(1, C_13_42)]), D_FM1_HHZZ0_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (0,  [(1, C_14_23)]), D_FM1_HHZZ0_U);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf 
                   (3,  [(1, C_14_23)]), D_FM1_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (3,  [(1, C_13_42)]), D_FM1_HHZZ1_U);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (3,  [(1, C_12_34)]), D_FM1_HHZZ1_T);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM1_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (6,  [(1, C_12_34)]), D_FM1_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_1_K_Matrix_cf
                   (6,  [(1, C_14_23)]), D_FM1_HHZZ1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM1_HHWW0_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (2,  [(1, C_13_42)]), D_FM1_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (1,  [(1, C_14_23)]), D_FM1_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (1,  [(1, C_13_42)]), D_FM1_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (2,  [(1, C_14_23)]), D_FM1_HHWW0_T);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (3,  [(1, C_14_23)]), D_FM1_HHWW1_S);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM1_HHWW1_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (4,  [(1, C_13_42)]), D_FM1_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (5,  [(1, C_12_34)]), D_FM1_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (8,  [(1, C_14_23)]), D_FM1_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (7,  [(1, C_12_34)]), D_FM1_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (5,  [(1, C_13_42)]), D_FM1_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (4,  [(1, C_12_34)]), D_FM1_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (7,  [(1, C_14_23)]), D_FM1_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_1_K_Matrix_cf
                   (8,  [(1, C_12_34)]), D_FM1_HHWW1_T) ]
        else
            []
      else
          [] 
          
    let k_matrix_2scalar_2gauge_m_7 =
      if Flags.k_matrix_tm then
        if Flags.higgs_matrix then
            [ ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM7_HHZZ0_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (0,  [(1, C_13_42)]), D_FM7_HHZZ0_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (0,  [(1, C_14_23)]), D_FM7_HHZZ0_U);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf 
                   (3,  [(1, C_14_23)]), D_FM7_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (3,  [(1, C_13_42)]), D_FM7_HHZZ1_U);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (3,  [(1, C_12_34)]), D_FM7_HHZZ1_T);
	      ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM7_HHZZ1_S);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (6,  [(1, C_12_34)]), D_FM7_HHZZ1_T);
              ((O H,O H,G Z,G Z), DScalar2_Vector2_m_7_K_Matrix_cf
                   (6,  [(1, C_14_23)]), D_FM7_HHZZ1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf 
                   (0,  [(1, C_12_34)]), D_FM7_HHWW0_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (2,  [(1, C_13_42)]), D_FM7_HHWW0_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (1,  [(1, C_14_23)]), D_FM7_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (1,  [(1, C_13_42)]), D_FM7_HHWW0_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (2,  [(1, C_14_23)]), D_FM7_HHWW0_T);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (3,  [(1, C_14_23)]), D_FM7_HHWW1_S);
	      ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (6,  [(1, C_13_42)]), D_FM7_HHWW1_S);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (4,  [(1, C_13_42)]), D_FM7_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (5,  [(1, C_12_34)]), D_FM7_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (8,  [(1, C_14_23)]), D_FM7_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (7,  [(1, C_12_34)]), D_FM7_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (5,  [(1, C_13_42)]), D_FM7_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (4,  [(1, C_12_34)]), D_FM7_HHWW1_T);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (7,  [(1, C_14_23)]), D_FM7_HHWW1_U);
              ((O H,O H,G Wp,G Wm), DScalar2_Vector2_m_7_K_Matrix_cf
                   (8,  [(1, C_12_34)]), D_FM7_HHWW1_T) ]
        else
            []
      else
          []      

    let k_matrix_4scalar =
      if Flags.k_matrix then
        if Flags.higgs_matrix then
            [ ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
	           (0,  [(1, C_12_34)]), D_Alpha_HHHH_S);
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
	           (0, [(1, C_13_42)]), D_Alpha_HHHH_T); 
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
                   (0, [(1, C_14_23)]), D_Alpha_HHHH_T); 
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
                   (3, [(1, C_14_23)]), D_Alpha_HHHH_S);
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
                   (3, [(1, C_13_42)]), D_Alpha_HHHH_T);
              ((O H,O H,O H,O H), DScalar4_K_Matrix_ms
                   (3, [(1, C_12_34)]), D_Alpha_HHHH_T) ]
        else
            []
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
      @ k_matrix_quartic_gauge_t_0 @ k_matrix_quartic_gauge_t_1 @ k_matrix_quartic_gauge_t_2
      @ k_matrix_quartic_gauge_t_rsi
      @ k_matrix_quartic_gauge_m_0 @ k_matrix_quartic_gauge_m_1 @ k_matrix_quartic_gauge_m_7

    let standard_gauge_higgs =
      [ ((O H, G Wp, G Wm), Scalar_Vector_Vector 1, G_HWW);
        ((O H, G Z, G Z), Scalar_Vector_Vector 1, G_HZZ) ]

    let standard_gauge_higgs4 =
      [ (O H, O H, G Wp, G Wm), Scalar2_Vector2 1, G_HHWW;
        (O H, O H, G Z, G Z), Scalar2_Vector2 1, G_HHZZ ]

    let dim8_gauge_higgs4 =
      [ (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_1 1, FS0_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_1 1, FS0_HHZZ;
        (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_2 1, FS1_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_2 1, FS1_HHZZ ]
        
    let dim8_gauge_higgs4_m =
      [ (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_m_0 1, FM0_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_m_0 1, FM0_HHZZ;
        (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_m_1 1, FM1_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_m_1 1, FM1_HHZZ;
        (O H, O H, G Wp, G Wm), Dim8_Scalar2_Vector2_m_7 1, FM7_HHWW;
        (O H, O H, G Z, G Z), Dim8_Scalar2_Vector2_m_7 1, FM7_HHZZ]    
       
    let standard_higgs =
      [ (O H, O H, O H), Scalar_Scalar_Scalar 1, G_H3 ]

    let standard_higgs4 =
      [ (O H, O H, O H, O H), Scalar4 1, G_H4 ]

    let fs_higgs4 =
      [ (O H, O H, O H, O H), Dim8_Scalar4 1, FS_H4 ]



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

(* Breaking thinks up
   \begin{subequations}
   \begin{align}
     \mathcal{L}^{\tau,H^4}_4 &=
       \left\lbrack (\partial_{\mu}H)(\partial^{\mu}H) \right\rbrack^2 \\
     \mathcal{L}^{\tau,H^4}_5 &=
       \left\lbrack (\partial_{\mu}H)(\partial^{\mu}H) \right\rbrack^2
   \end{align}
   \end{subequations}
   and
   \begin{subequations}
   \begin{align}
     \mathcal{L}^{\tau,H^2V^2}_4 &= \frac{g^2v_{\mathrm{F}}^2}{2}
              (\partial_{\mu}H)(\partial^{\mu}H) V_{\mu}V^{\mu}   \\
     \mathcal{L}^{\tau,H^2V^2}_5 &= \frac{g^2v_{\mathrm{F}}^2}{2}
              (\partial_{\mu}H)(\partial_{\nu}H) V_{\mu}V_{\nu}
   \end{align}
   \end{subequations}
   i.\,e.
   \begin{subequations}
   \begin{align}
     \mathcal{L}^{\tau,H^2V^2}_4 &=
        \frac{g^2v_{\mathrm{F}}^2}{2}
          \left\lbrack
              (\partial_{\mu}H)(\partial^{\mu}H) W^+_{\nu}W^{-,\nu}
            + \frac{1}{2\cos^2\theta_{w}} (\partial_{\mu}H)(\partial^{\mu}H) Z_{\nu} Z^{\nu}
          \right\rbrack \\
     \mathcal{L}^{\tau,H^2V^2}_5 &=
          \frac{g^2v_{\mathrm{F}}^2}{2}
          \left\lbrack
              (W^{+,\mu}\partial_{\mu}H) (W^{-,\nu}\partial_{\nu}H)
            + \frac{1}{2\cos^2\theta_{w}} (Z^{\mu}\partial_{\mu}H)(Z^{\nu}\partial_{\nu}H)
          \right\rbrack
   \end{align}
   \end{subequations} *)

(* \begin{multline}
     \tau^4_8 \mathcal{L}^{\tau,H^2V^2}_4 + \tau^5_8 \mathcal{L}^{\tau,H^2V^2}_5 = \\
       - \frac{g^2v_{\mathrm{F}}^2}{2} \Biggl\lbrack
            2\tau^4_8
              \frac{1}{2}(\ii\partial_{\mu}H)(\ii\partial^{\mu}H) W^+_{\nu}W^{-,\nu}
          + \tau^5_8
              (W^{+,\mu}\ii\partial_{\mu}H) (W^{-,\nu}\ii\partial_{\nu}H) \\
          + \frac{2\tau^4_8}{\cos^2\theta_{w}}
              \frac{1}{4} (\ii\partial_{\mu}H)(\ii\partial^{\mu}H) Z_{\nu} Z^{\nu}
          + \frac{\tau^5_8}{\cos^2\theta_{w}}
              \frac{1}{2} (Z^{\mu}\ii\partial_{\mu}H)(Z^{\nu}\ii\partial_{\nu}H)
          \Biggr\rbrack
   \end{multline}
   where the two powers of $\ii$ make the sign conveniently negative,
   i.\,e.
   \begin{subequations}
   \begin{align}
     \alpha_{(\partial H)^2W^2}^2 &= \tau^4_8 g^2v_{\mathrm{F}}^2\\
     \alpha_{(\partial HW)^2}^2 &= \frac{\tau^5_8 g^2v_{\mathrm{F}}^2}{2}  \\
     \alpha_{(\partial H)^2Z^2}^2 &= \frac{\tau^4_8 g^2v_{\mathrm{F}}^2}{\cos^2\theta_{w}} \\ 
     \alpha_{(\partial HZ)^2}^2 &=\frac{\tau^5_8 g^2v_{\mathrm{F}}^2}{2\cos^2\theta_{w}}
   \end{align}
   \end{subequations} *)

    let anomalous_gauge_higgs =
      [ (O H, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_HGaGa_anom;
        (O H, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_HGaZ_anom;
        (O H, G Z, G Z), Dim5_Scalar_Gauge2 1, G_HZZ_anom;
        (O H, G Wp, G Wm), Dim5_Scalar_Gauge2 1, G_HWW_anom;
        (O H, G Ga, G Z), Dim5_Scalar_Vector_Vector_U 1, G_HGaZ_u;
        (O H, G Z, G Z), Dim5_Scalar_Vector_Vector_U 1, G_HZZ_u;
        (O H, G Wp, G Wm), Dim5_Scalar_Vector_Vector_U 1, G_HWW_u;
        (O H, G Wm, G Wp), Dim5_Scalar_Vector_Vector_U 1, G_HWW_u
      ]

    let anomalous_gauge_higgs4 =
      []

    let anomalous_higgs =
      []

    let higgs_triangle_vertices = 
      if Flags.higgs_triangle then
        [ (O H, G Ga, G Ga), Dim5_Scalar_Gauge2 1, G_HGaGa;
          (O H, G Ga, G Z), Dim5_Scalar_Gauge2 1, G_HGaZ;
          (O H, G Gl, G Gl), Dim5_Scalar_Gauge2 1, G_Hgg ]
      else
        []

    let anomalous_higgs4 =
      []

    let gauge_higgs =
      if Flags.higgs_anom then
        standard_gauge_higgs @ anomalous_gauge_higgs
      else
        standard_gauge_higgs


    let gauge_higgs4 =
      ( if Flags.higgs_anom then
          standard_gauge_higgs4 @ anomalous_gauge_higgs4
        else
          standard_gauge_higgs4 ) @
      ( if Flags.higgs_matrix then
          (dim8_gauge_higgs4 @ dim8_gauge_higgs4_m @ k_matrix_2scalar_2gauge 
           @ k_matrix_2scalar_2gauge_m @ k_matrix_2scalar_2gauge_m_1 @ k_matrix_2scalar_2gauge_m_7)
	 else
	   [] )

    let higgs =
      if Flags.higgs_anom then
        standard_higgs @ anomalous_higgs
      else
        standard_higgs

    let higgs4 =
      ( if Flags.higgs_anom then
          standard_higgs4 @ anomalous_higgs4
        else
          standard_higgs4 ) @ 
      ( if Flags.higgs_matrix then
          (fs_higgs4 @ k_matrix_4scalar )
	 else
	   [] )


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
    Effectivly this doesn't matter for SSC, because $(-1)^2=1$.
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
   Old representation
   \begin{equation}
    \mathcal{L}_{\sigma}=
              -\frac{g_\sigma v}{2} \text{tr}
	      \left\lbrack \mathbf{V}_\mu \mathbf{V}^\mu \right\rbrack \sigma
   \end{equation}
*)

(* \begin{dubious}
   Transversal couplings like rsigma3t and rf3t are to be calculated in the new
   higgs matrix representation.
   \end{dubious} *)


    let rsigma3 =
      [ ((O Rsigma, G Wp, G Wm), Scalar_Vector_Vector 1, G_SWW);
        ((O Rsigma, G Z, G Z), Scalar_Vector_Vector 1, G_SZZ) ]

    let rsigma3h =
      [ ((O Rsigma, O H, O H), Dim5_Scalar_Scalar2 1, G_SHH) ]

    let rsigma3t =
      [ ((O Rsigma, G Wp, G Wm), Scalar_Vector_Vector_t 1, G_SWW_T);
        ((O Rsigma, G Z, G Z), Scalar_Vector_Vector_t 1, G_SZZ_T);
        ((O Rsigma, G Ga, G Ga), Scalar_Vector_Vector_t 1, G_SAA_T);
        ((O Rsigma, G Ga, G Z), Scalar_Vector_Vector_t 1, G_SAZ_T) ]

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
        ((O Rphisn, G Wp, G Wm), Scalar_Vector_Vector 1, G_PSNWW);
        ((O Rphisn, G Z, G Z), Scalar_Vector_Vector 1, G_PSNZZ) ;
        ((O Rphip, G Z, G Wm), Scalar_Vector_Vector 1, G_PWZ) ;
        ((O Rphipp, G Wm, G Wm), Scalar_Vector_Vector 1, G_PWW) ;
        ((O Rphim, G Wp, G Z), Scalar_Vector_Vector 1, G_PWZ) ;
        ((O Rphimm, G Wp, G Wp), Scalar_Vector_Vector 1, G_PWW) ]

    let rphi3h =
      [ ((O Rphisn, O H, O H), Dim5_Scalar_Scalar2 1, G_PSNHH) ]

(* Tensor IsoScalar
*)
    let rf3 =
      [ ((O Rf, G Wp, G Wm), Tensor_2_Vector_Vector_1 1, G_FWW);
        ((O Rf, G Z, G Z), Tensor_2_Vector_Vector_1 1, G_FZZ) ]

    let rf3cf =
      [ ((O Rf, G Wp, G Wm), Tensor_2_Vector_Vector 1, G_FWW);
        ((O Rf, G Z, G Z), Tensor_2_Vector_Vector 1, G_FZZ);
        ((O Rf, G Wp, G Wm), Tensor_2_Vector_Vector_cf 1, G_FWW_CF);
        ((O Rf, G Z, G Z), Tensor_2_Vector_Vector_cf 1, G_FZZ_CF) ]

    let rff3cf =
      [ ((O Rff, G Wp, G Wm), Tensor_2_Vector_Vector 1, G_FFWW);
        ((O Rff, G Z, G Z), Tensor_2_Vector_Vector 1, G_FFZZ);
        ((O Rff, G Wp, G Wm), Tensor_2_Vector_Vector_cf 1, G_FFWW_CF);
        ((O Rff, G Z, G Z), Tensor_2_Vector_Vector_cf 1, G_FFZZ_CF) ]

    let rfv3cf =
      [ ((O Rfv, G Wp, G Wm), TensorVector_Vector_Vector 1, G_FVWW);
        ((O Rfv, G Z, G Z), TensorVector_Vector_Vector 1, G_FVZZ);
        ((O Rfv, G Wp, G Wm), TensorVector_Vector_Vector_cf 1, G_FVWW_CF);
        ((O Rfv, G Z, G Z), TensorVector_Vector_Vector_cf 1, G_FVZZ_CF) ]

    let rfddphi3cf =
      [ ((O Rfphi, G Wp, G Wm), TensorScalar_Vector_Vector 1, G_FDDSWW);
        ((O Rfphi, G Z, G Z), TensorScalar_Vector_Vector 1, G_FDDSZZ);
        ((O Rfphi, G Wp, G Wm), TensorScalar_Vector_Vector_cf 1, G_FDDSWW_CF);
        ((O Rfphi, G Z, G Z), TensorScalar_Vector_Vector_cf 1, G_FDDSZZ_CF) ]

    let rfphi3cf =
      [ ((O Rfphi, G Wp, G Wm), Scalar_Vector_Vector 1, G_FSWW);
        ((O Rfphi, G Z, G Z), Scalar_Vector_Vector 1, G_FSZZ) ]

    let rf3h =
      [ ((O Rf, O H, O H), Tensor_2_Scalar_Scalar 1, G_FHH);
        ((O Rf, O H, O H), Tensor_2_Scalar_Scalar_cf 1, G_FHH_CF) ]

    let rff3h =
      [ ((O Rff, O H, O H), Tensor_2_Scalar_Scalar 1, G_FFHH);
        ((O Rff, O H, O H), Tensor_2_Scalar_Scalar_cf 1, G_FFHH_CF);
        ((O Rfv, O H, O H), TensorVector_Scalar_Scalar 1, G_FVHH);
        ((O Rfv, O H, O H), TensorVector_Scalar_Scalar_cf 1, G_FVHH_CF);
        ((O Rfphi, O H, O H), TensorScalar_Scalar_Scalar 1, G_FDDSHH);
        ((O Rfphi, O H, O H), TensorScalar_Scalar_Scalar_cf 1, G_FDDSHH_CF);
        ((O Rfphi, O H, O H), Dim5_Scalar_Scalar2 1, G_FSHH) ]
     
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
        ((O Rtsn, G Wp, G Wm), Tensor_2_Vector_Vector_1 1, G_TSNWW);
        ((O Rtsn, G Z, G Z), Tensor_2_Vector_Vector_1 1, G_TSNZZ) ;
        ((O Rtp, G Z, G Wm), Tensor_2_Vector_Vector_1 1, G_TWZ) ;
        ((O Rtpp, G Wm, G Wm), Tensor_2_Vector_Vector_1 1, G_TWW) ;
        ((O Rtm, G Wp, G Z), Tensor_2_Vector_Vector_1 1, G_TWZ) ;
        ((O Rtmm, G Wp, G Wp), Tensor_2_Vector_Vector_1 1, G_TWW) ]

    let rt3cf =
      [ ((O Rtn, G Wp, G Wm), Tensor_2_Vector_Vector 1, G_TNWW);
        ((O Rtn, G Z, G Z), Tensor_2_Vector_Vector 1, G_TNZZ) ;
        ((O Rtsn, G Wp, G Wm), Tensor_2_Vector_Vector 1, G_TSNWW);
        ((O Rtsn, G Z, G Z), Tensor_2_Vector_Vector 1, G_TSNZZ) ;
        ((O Rtp, G Z, G Wm), Tensor_2_Vector_Vector 1, G_TWZ) ;
        ((O Rtpp, G Wm, G Wm), Tensor_2_Vector_Vector 1, G_TWW) ;
        ((O Rtm, G Wp, G Z), Tensor_2_Vector_Vector 1, G_TWZ) ;
        ((O Rtmm, G Wp, G Wp), Tensor_2_Vector_Vector 1, G_TWW);
        ((O Rtn, G Wp, G Wm), Tensor_2_Vector_Vector_cf 1, G_TNWW_CF);
        ((O Rtn, G Z, G Z), Tensor_2_Vector_Vector_cf 1, G_TNZZ_CF) ;
        ((O Rtsn, G Wp, G Wm), Tensor_2_Vector_Vector_cf 1, G_TSNWW_CF);
        ((O Rtsn, G Z, G Z), Tensor_2_Vector_Vector_cf 1, G_TSNZZ_CF) ;
        ((O Rtp, G Z, G Wm), Tensor_2_Vector_Vector_cf 1, G_TWZ_CF) ;
        ((O Rtpp, G Wm, G Wm), Tensor_2_Vector_Vector_cf 1, G_TWW_CF) ;
        ((O Rtm, G Wp, G Z), Tensor_2_Vector_Vector_cf 1, G_TWZ_CF) ;
        ((O Rtmm, G Wp, G Wp), Tensor_2_Vector_Vector_cf 1, G_TWW_CF) ]


(* Anomalous trilinear interactions $f_i f_j V$ and $ttH$:
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

(* \begin{equation}
     \Delta\mathcal{L}_{ttH} =
        - \frac{1}{\sqrt{2}} \bar{t} (Y_V(k^2)+iY_A(k^2)\gamma_5)t H
   \end{equation} *)

    let anomalous_ttH =
      if Flags.top_anom then
        [ ((M (U (-3)), O H, M (U 3)), FBF (1, Psibar, SPM, Psi), G_SP_ttH) ]
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
       yukawa @ triple_gauge @
       gauge_higgs @ higgs @ higgs_triangle_vertices 
       @ goldstone_vertices @
       rsigma3 @ rsigma3t @ rphi3 @
       ( if Flags.cf_arbitrary then
	    ( rt3cf @ rff3cf @
             rfv3cf @ rfphi3cf @ rfddphi3cf )
	 else
	    (rf3 @ rt3 ) ) @
       rf3t @ 
       ( if Flags.higgs_matrix then
	    (rsigma3h @ rff3h )
	 else
	    [] ) @
       anomalous_ttA @ anomalous_bbA @
       anomalous_ttZ @ anomalous_bbZ @
       anomalous_tbW @ anomalous_tbWA @ anomalous_tbWZ @
       anomalous_ttWW @ anomalous_bbWW @
       anomalous_ttG @ anomalous_ttGG @
       anomalous_ttH @
       anomalous_top_qGuG @ anomalous_top_qBuB @
       anomalous_top_qW @ anomalous_top_DuDd)

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
      | "H" -> O H | "Rsigma" -> O Rsigma
      | "Rphi0" -> O Rphin
      | "Rphis0" -> O Rphisn
      | "Rphi+" -> O Rphip |  "Rphi-" -> O Rphim
      | "Rphi++" -> O Rphip |  "Rphi--" -> O Rphimm
      | "Rf" -> O Rf
      | "Rff" -> O Rff
      | "Rfv" -> O Rfv
      | "Rfphi" -> O Rfphi
      | "Rt0" -> O Rtn
      | "Rts0" -> O Rtsn
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
      | _ -> invalid_arg "Modellib_BSM.SSC_AltT.flavor_of_string"

    let flavor_to_string = function
      | M f ->
          begin match f with
          | L 1 -> "e-" | L (-1) -> "e+"
          | L 2 -> "mu-" | L (-2) -> "mu+"
          | L 3 -> "tau-" | L (-3) -> "tau+"
          | L _ -> invalid_arg
                "Modellib_BSM.SSC_AltT.flavor_to_string: invalid lepton"
          | N 1 -> "nue" | N (-1) -> "nuebar"
          | N 2 -> "numu" | N (-2) -> "numubar"
          | N 3 -> "nutau" | N (-3) -> "nutaubar"
          | N _ -> invalid_arg
                "Modellib_BSM.SSC_AltT.flavor_to_string: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "ubar"
          | U 2 -> "c" | U (-2) -> "cbar"
          | U 3 -> "t" | U (-3) -> "tbar"
          | U _ -> invalid_arg
                "Modellib_BSM.SSC_AltT.flavor_to_string: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "dbar"
          | D 2 -> "s" | D (-2) -> "sbar"
          | D 3 -> "b" | D (-3) -> "bbar"
          | D _ -> invalid_arg
                "Modellib_BSM.SSC_AltT.flavor_to_string: invalid down type quark"
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
          | H -> "H" | Rsigma -> "Rsigma"
          | Rphin -> "Rphin" | Rphip -> "Rphi+" | Rphim -> "Rphi-"
          | Rphipp -> "Rphi++" | Rphimm -> "Rphi--"
          | Rphisn -> "Rphisn"
          | Rf -> "Rf" 
          | Rff -> "Rff" | Rfv -> "Rfv" | Rfphi -> "Rfphi"
          | Rtn -> "Rtn" | Rtsn -> "Rtsn" | Rtp -> "Rt+" | Rtm -> "Rt-"
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
                "Modellib_BSM.SSC_AltT.flavor_to_TeX: invalid lepton"
          | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
          | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
          | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
          | N _ -> invalid_arg
                "Modellib_BSM.SSC_AltT.flavor_to_TeX: invalid neutrino"
          | U 1 -> "u" | U (-1) -> "\\bar{u}"
          | U 2 -> "c" | U (-2) -> "\\bar{c}"
          | U 3 -> "t" | U (-3) -> "\\bar{t}"
          | U _ -> invalid_arg
                "Modellib_BSM.SSC_AltT.flavor_to_TeX: invalid up type quark"
          | D 1 -> "d" | D (-1) -> "\\bar{d}"
          | D 2 -> "s" | D (-2) -> "\\bar{s}"
          | D 3 -> "b" | D (-3) -> "\\bar{b}"
          | D _ -> invalid_arg
                "Modellib_BSM.SSC_AltT.flavor_to_TeX: invalid down type quark"
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
          | H -> "H" | Rsigma -> "\\sigma"
          | Rphip -> "\\phi^+" | Rphim -> "\\phi^-" | Rphin -> "\\phi^0" 
          | Rphisn -> "\\phi_s^0" 
          | Rphipp -> "\\phi^{++}" | Rphimm -> "\\phi^{--}"
          | Rf -> "f"
          | Rff -> "f^f" | Rfv -> "f^v" | Rfphi -> "f^s"
          | Rtp -> "t^+" | Rtm -> "t^-" | Rtn -> "t^0" | Rtsn -> "t_s^0" 
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
          | H -> "h" | Rsigma -> "rsi"
          | Rphip -> "rpp" | Rphim -> "rpm" | Rphin -> "rpn"
          | Rphisn -> "rpsn"
          | Rphipp -> "rppp" | Rphimm -> "rpmm"
          | Rf -> "rf"
          | Rff -> "rff" | Rfv -> "rfv" | Rfphi -> "rfphi" 
          | Rtp -> "rtp" | Rtm -> "rtm" | Rtn -> "rtn" 
          | Rtsn -> "rtsn"
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
          | H -> 25 | Rsigma -> 45
          | Rphin -> 46 | Rphip | Rphim -> 47 
          | Rphipp | Rphimm -> 48
          | Rphisn -> 49
          | Rf -> 52
          | Rtn -> 53 | Rtp | Rtm -> 54 
          | Rtpp | Rtmm -> 55
          | Rff -> 56 | Rfv -> 57 | Rfphi -> 58
          | Rtsn -> 59
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
      | G_SP_ttH -> "gsp_tth"
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
      | FT0_WWWW0 -> "at0ww0" | FT0_WWWW2 -> "at0ww2"
      | FT0_ZZWW0 -> "at0zw0" | FT0_ZZWW1 -> "at0zw1"
      | FT0_ZZZZ  -> "at0zz"  | FT0_AAAA  -> "at0aa"
      | FT0_AAWW0 -> "at0aw0" | FT0_AAWW1 -> "at0aw1"
      | FT0_AAZZ -> "at0az" 
      | FT0_AZWW0 -> "at0azw0" | FT0_AZWW1 -> "at0azw1"
      | FT0_AAAZ  -> "at03az"  | FT0_AZZZ  -> "at0a3z"
      | FT1_WWWW0 -> "at1ww0" | FT1_WWWW2 -> "at1ww2"
      | FT1_ZZWW0 -> "at1zw0" | FT1_ZZWW1 -> "at1zw1"
      | FT1_ZZZZ  -> "at1zz"  | FT1_AAAA  -> "at1aa"
      | FT1_AAWW0 -> "at1aw0" | FT1_AAWW1 -> "at1aw1"
      | FT1_AAZZ -> "at1az"   
      | FT1_AZWW0 -> "at1azw0" | FT1_AZWW1 -> "at1azw1"
      | FT1_AAAZ  -> "at13az"  | FT1_AZZZ  -> "at1a3z"
      | FT2_WWWW0 -> "at2ww0" | FT2_WWWW2 -> "at2ww2"
      | FT2_ZZWW0 -> "at2zw0" | FT2_ZZWW1 -> "at2zw1"
      | FT2_ZZZZ  -> "at2zz"  | FT2_AAAA  -> "at2aa"
      | FT2_AAWW0 -> "at2aw0" | FT2_AAWW1 -> "at2aw1"
      | FT2_AAZZ -> "at2az"   
      | FT2_AZWW0 -> "at2azw0" | FT2_AZWW1 -> "at2azw1"
      | FT2_AAAZ  -> "at23az"  | FT2_AZZZ  -> "at2a3z"
      | FM0_WWWW0 -> "am0ww0,am0ww0" | FM0_WWWW2 -> "am0ww2,am0ww2"
      | FM0_ZZWW0 -> "am0zw0/costhw**2,am0zw0*costhw**2" | FM0_ZZWW1 -> "am0zw1/costhw**2,am0zw1*costhw**2"
      | FM0_ZZZZ  -> "am0zz,am0zz" 
      | FM1_WWWW0 -> "am1ww0,am1ww0" | FM1_WWWW2 -> "am1ww2,am1ww2"
      | FM1_ZZWW0 -> "am1zw0/costhw**2,am1zw0*costhw**2" | FM1_ZZWW1 -> "am1zw1/costhw**2,am1zw1*costhw**2"
      | FM1_ZZZZ  -> "am1zz,am1zz"  
      | FM7_WWWW0 -> "am7ww0,am7ww0,am7ww0" | FM7_WWWW2 -> "am7ww2,am7ww2,am7ww2"
      | FM7_ZZWW0 -> "am7zw0/costhw**2,am7zw0,am7zw0*costhw**2" | FM7_ZZWW1 -> "am7zw1/costhw**2,am7zw1,am7zw1*costhw**2"
      | FM7_ZZZZ  -> "am7zz,am7zz,am7zz"
      | FS0_HHWW -> "fs0hhww" | FS0_HHZZ -> "fs0hhzz"
      | FS1_HHWW -> "fs1hhww" | FS1_HHZZ -> "fs1hhzz"
      | FS_H4 -> "fsh4"
      | FM0_HHWW -> "fm0hhww" | FM0_HHZZ -> "fm0hhzz"
      | FM1_HHWW -> "fm1hhww" | FM1_HHZZ -> "fm1hhzz"   
      | FM7_HHWW -> "fm7hhww" | FM7_HHZZ -> "fm7hhzz"
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
      | D_Alpha_ZZZZ_S  -> "dalz4_s(gkm,mkm,"
      | D_Alpha_ZZZZ_T  -> "dalz4_t(gkm,mkm,"
      | D_FT0_ZZWW0_S -> "datzz0_s_0(gkm,mkm,"
      | D_FT0_ZZWW0_T -> "datzz0_t_0(gkm,mkm,"
      | D_FT0_ZZWW0_U -> "datzz0_u_0(gkm,mkm,"
      | D_FT0_ZZWW1_S -> "datzz1_s_0(gkm,mkm,"
      | D_FT0_ZZWW1_T -> "datzz1_t_0(gkm,mkm,"
      | D_FT0_ZZWW1_U -> "datzz1_u_0(gkm,mkm,"
      | D_FT0_WWWW0_S -> "datww0_s_0(gkm,mkm,"
      | D_FT0_WWWW0_T -> "datww0_t_0(gkm,mkm,"
      | D_FT0_WWWW0_U -> "datww0_u_0(gkm,mkm,"
      | D_FT0_WWWW2_S -> "datww2_s_0(gkm,mkm,"
      | D_FT0_WWWW2_T -> "datww2_t_0(gkm,mkm,"
      | D_FT0_WWWW2_U -> "datww2_u_0(gkm,mkm,"
      | D_FT0_ZZZZ_S  -> "datz4_s_0(gkm,mkm,"
      | D_FT0_ZZZZ_T  -> "datz4_t_0(gkm,mkm,"
      | D_FT0_ZZZZ_U  -> "datz4_u_0(gkm,mkm,"
      | D_FT0_AAAA_S  -> "data4_s_0(gkm,mkm,"
      | D_FT0_AAAA_T  -> "data4_t_0(gkm,mkm,"
      | D_FT0_AAAA_U  -> "data4_u_0(gkm,mkm,"  
      | D_FT0_AAWW0_S -> "dataw0_s_0(gkm,mkm,"
      | D_FT0_AAWW0_T -> "dataw0_t_0(gkm,mkm,"
      | D_FT0_AAWW0_U -> "dataw0_u_0(gkm,mkm,"
      | D_FT0_AAWW1_S -> "dataw1_s_0(gkm,mkm,"
      | D_FT0_AAWW1_T -> "dataw1_t_0(gkm,mkm,"
      | D_FT0_AAWW1_U -> "dataw1_u_0(gkm,mkm,"
      | D_FT0_AAZZ_S  -> "dataz_s_0(gkm,mkm,"
      | D_FT0_AAZZ_T  -> "dataz_t_0(gkm,mkm,"
      | D_FT0_AAZZ_U  -> "dataz_u_0(gkm,mkm,"  
      | D_FT0_AZWW0_S -> "datazw0_s_0(gkm,mkm,"
      | D_FT0_AZWW0_T -> "datazw0_t_0(gkm,mkm,"
      | D_FT0_AZWW0_U -> "datazw0_u_0(gkm,mkm,"
      | D_FT0_AZWW1_S -> "datazw0_s_1(gkm,mkm,"
      | D_FT0_AZWW1_T -> "datazw0_t_1(gkm,mkm,"
      | D_FT0_AZWW1_U -> "datazw0_u_1(gkm,mkm," 
      | D_FT0_AAAZ_S -> "dat3az_s_0(gkm,mkm,"
      | D_FT0_AAAZ_T -> "dat3az_t_0(gkm,mkm,"
      | D_FT0_AAAZ_U -> "dat3az_u_0(gkm,mkm," 
      | D_FT0_AZZZ_S -> "data3z_s_0(gkm,mkm,"
      | D_FT0_AZZZ_T -> "data3z_t_0(gkm,mkm,"
      | D_FT0_AZZZ_U -> "data3z_u_0(gkm,mkm,"             
      | D_FT1_ZZWW0_S -> "datzz0_s_1(gkm,mkm,"
      | D_FT1_ZZWW0_T -> "datzz0_t_1(gkm,mkm,"
      | D_FT1_ZZWW0_U -> "datzz0_u_1(gkm,mkm,"
      | D_FT1_ZZWW1_S -> "datzz1_s_1(gkm,mkm,"
      | D_FT1_ZZWW1_T -> "datzz1_t_1(gkm,mkm,"
      | D_FT1_ZZWW1_U -> "datzz1_u_1(gkm,mkm,"
      | D_FT1_WWWW0_S -> "datww0_s_1(gkm,mkm,"
      | D_FT1_WWWW0_T -> "datww0_t_1(gkm,mkm,"
      | D_FT1_WWWW0_U -> "datww0_u_1(gkm,mkm,"
      | D_FT1_WWWW2_S -> "datww2_s_1(gkm,mkm,"
      | D_FT1_WWWW2_T -> "datww2_t_1(gkm,mkm,"
      | D_FT1_WWWW2_U -> "datww2_u_1(gkm,mkm,"
      | D_FT1_ZZZZ_S  -> "datz4_s_1(gkm,mkm,"
      | D_FT1_ZZZZ_T  -> "datz4_t_1(gkm,mkm,"
      | D_FT1_ZZZZ_U  -> "datz4_u_1(gkm,mkm,"
      | D_FT1_AAAA_S  -> "data4_s_1(gkm,mkm,"
      | D_FT1_AAAA_T  -> "data4_t_1(gkm,mkm,"
      | D_FT1_AAAA_U  -> "data4_u_1(gkm,mkm,"  
      | D_FT1_AAWW0_S -> "dataw0_s_1(gkm,mkm,"
      | D_FT1_AAWW0_T -> "dataw0_t_1(gkm,mkm,"
      | D_FT1_AAWW0_U -> "dataw0_u_1(gkm,mkm,"
      | D_FT1_AAWW1_S -> "dataw1_s_1(gkm,mkm,"
      | D_FT1_AAWW1_T -> "dataw1_t_1(gkm,mkm,"
      | D_FT1_AAWW1_U -> "dataw1_u_1(gkm,mkm,"
      | D_FT1_AAZZ_S  -> "dataz_s_1(gkm,mkm,"
      | D_FT1_AAZZ_T  -> "dataz_t_1(gkm,mkm,"
      | D_FT1_AAZZ_U  -> "dataz_u_1(gkm,mkm,"
      | D_FT1_AZWW0_S -> "datazw0_s_1(gkm,mkm,"
      | D_FT1_AZWW0_T -> "datazw0_t_1(gkm,mkm,"
      | D_FT1_AZWW0_U -> "datazw0_u_1(gkm,mkm,"
      | D_FT1_AZWW1_S -> "datazw1_s_1(gkm,mkm,"
      | D_FT1_AZWW1_T -> "datazw1_t_1(gkm,mkm,"
      | D_FT1_AZWW1_U -> "datazw1_u_1(gkm,mkm," 
      | D_FT1_AAAZ_S -> "dat3az_s_1(gkm,mkm,"
      | D_FT1_AAAZ_T -> "dat3az_t_1(gkm,mkm,"
      | D_FT1_AAAZ_U -> "dat3az_u_1(gkm,mkm," 
      | D_FT1_AZZZ_S -> "data3z_s_1(gkm,mkm,"
      | D_FT1_AZZZ_T -> "data3z_t_1(gkm,mkm,"
      | D_FT1_AZZZ_U -> "data3z_u_1(gkm,mkm,"      
      | D_FT2_ZZWW0_S -> "datzz0_s_2(gkm,mkm,"
      | D_FT2_ZZWW0_T -> "datzz0_t_2(gkm,mkm,"
      | D_FT2_ZZWW0_U -> "datzz0_u_2(gkm,mkm,"
      | D_FT2_ZZWW1_S -> "datzz1_s_2(gkm,mkm,"
      | D_FT2_ZZWW1_T -> "datzz1_t_2(gkm,mkm,"
      | D_FT2_ZZWW1_U -> "datzz1_u_2(gkm,mkm,"
      | D_FT2_WWWW0_S -> "datww0_s_2(gkm,mkm,"
      | D_FT2_WWWW0_T -> "datww0_t_2(gkm,mkm,"
      | D_FT2_WWWW0_U -> "datww0_u_2(gkm,mkm,"
      | D_FT2_WWWW2_S -> "datww2_s_2(gkm,mkm,"
      | D_FT2_WWWW2_T -> "datww2_t_2(gkm,mkm,"
      | D_FT2_WWWW2_U -> "datww2_u_2(gkm,mkm,"
      | D_FT2_ZZZZ_S  -> "datz4_s_2(gkm,mkm,"
      | D_FT2_ZZZZ_T  -> "datz4_t_2(gkm,mkm,"
      | D_FT2_ZZZZ_U  -> "datz4_u_2(gkm,mkm," 
      | D_FT2_AAAA_S  -> "data4_s_2(gkm,mkm,"
      | D_FT2_AAAA_T  -> "data4_t_2(gkm,mkm,"
      | D_FT2_AAAA_U  -> "data4_u_2(gkm,mkm,"  
      | D_FT2_AAWW0_S -> "dataw0_s_2(gkm,mkm,"
      | D_FT2_AAWW0_T -> "dataw0_t_2(gkm,mkm,"
      | D_FT2_AAWW0_U -> "dataw0_u_2(gkm,mkm,"
      | D_FT2_AAWW1_S -> "dataw1_s_2(gkm,mkm,"
      | D_FT2_AAWW1_T -> "dataw1_t_2(gkm,mkm,"
      | D_FT2_AAWW1_U -> "dataw1_u_2(gkm,mkm,"
      | D_FT2_AAZZ_S  -> "dataz_s_2(gkm,mkm,"
      | D_FT2_AAZZ_T  -> "dataz_t_2(gkm,mkm,"
      | D_FT2_AAZZ_U  -> "dataz_u_2(gkm,mkm,"    
      | D_FT2_AZWW0_S -> "datazw0_s_2(gkm,mkm,"
      | D_FT2_AZWW0_T -> "datazw0_t_2(gkm,mkm,"
      | D_FT2_AZWW0_U -> "datazw0_u_2(gkm,mkm,"
      | D_FT2_AZWW1_S -> "datazw1_s_2(gkm,mkm,"
      | D_FT2_AZWW1_T -> "datazw1_t_2(gkm,mkm,"
      | D_FT2_AZWW1_U -> "datazw1_u_2(gkm,mkm,"
      | D_FT2_AAAZ_S -> "dat3az_s_2(gkm,mkm,"
      | D_FT2_AAAZ_T -> "dat3az_t_2(gkm,mkm,"
      | D_FT2_AAAZ_U -> "dat3az_u_2(gkm,mkm," 
      | D_FT2_AZZZ_S -> "data3z_s_2(gkm,mkm,"
      | D_FT2_AZZZ_T -> "data3z_t_2(gkm,mkm,"
      | D_FT2_AZZZ_U -> "data3z_u_2(gkm,mkm,"
      | D_FTrsi_ZZWW0_S -> "datzz0_s_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW0_T -> "datzz0_t_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW0_U -> "datzz0_u_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW1_S -> "datzz1_s_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW1_T -> "datzz1_t_rsi(gkm,mkm,"
      | D_FTrsi_ZZWW1_U -> "datzz1_u_rsi(gkm,mkm,"
      | D_FTrsi_WWWW0_S -> "datww0_s_rsi(gkm,mkm,"
      | D_FTrsi_WWWW0_T -> "datww0_t_rsi(gkm,mkm,"
      | D_FTrsi_WWWW0_U -> "datww0_u_rsi(gkm,mkm,"
      | D_FTrsi_WWWW2_S -> "datww2_s_rsi(gkm,mkm,"
      | D_FTrsi_WWWW2_T -> "datww2_t_rsi(gkm,mkm,"
      | D_FTrsi_WWWW2_U -> "datww2_u_rsi(gkm,mkm,"
      | D_FTrsi_ZZZZ_S  -> "datz4_s_rsi(gkm,mkm,"
      | D_FTrsi_ZZZZ_T  -> "datz4_t_rsi(gkm,mkm,"
      | D_FTrsi_ZZZZ_U  -> "datz4_u_rsi(gkm,mkm," 
      | D_FTrsi_AAAA_S  -> "data4_s_rsi(gkm,mkm,"
      | D_FTrsi_AAAA_T  -> "data4_t_rsi(gkm,mkm,"
      | D_FTrsi_AAAA_U  -> "data4_u_rsi(gkm,mkm,"  
      | D_FTrsi_AAWW0_S -> "dataw0_s_rsi(gkm,mkm,"
      | D_FTrsi_AAWW0_T -> "dataw0_t_rsi(gkm,mkm,"
      | D_FTrsi_AAWW0_U -> "dataw0_u_rsi(gkm,mkm,"
      | D_FTrsi_AAWW1_S -> "dataw1_s_rsi(gkm,mkm,"
      | D_FTrsi_AAWW1_T -> "dataw1_t_rsi(gkm,mkm,"
      | D_FTrsi_AAWW1_U -> "dataw1_u_rsi(gkm,mkm,"
      | D_FTrsi_AAZZ_S  -> "dataz_s_rsi(gkm,mkm,"
      | D_FTrsi_AAZZ_T  -> "dataz_t_rsi(gkm,mkm,"
      | D_FTrsi_AAZZ_U  -> "dataz_u_rsi(gkm,mkm,"    
      | D_FTrsi_AZWW0_S -> "datazw0_s_rsi(gkm,mkm,"
      | D_FTrsi_AZWW0_T -> "datazw0_t_rsi(gkm,mkm,"
      | D_FTrsi_AZWW0_U -> "datazw0_u_rsi(gkm,mkm,"
      | D_FTrsi_AZWW1_S -> "datazw1_s_rsi(gkm,mkm,"
      | D_FTrsi_AZWW1_T -> "datazw1_t_rsi(gkm,mkm,"
      | D_FTrsi_AZWW1_U -> "datazw1_u_rsi(gkm,mkm,"
      | D_FTrsi_AAAZ_S -> "dat3az_s_rsi(gkm,mkm,"
      | D_FTrsi_AAAZ_T -> "dat3az_t_rsi(gkm,mkm,"
      | D_FTrsi_AAAZ_U -> "dat3az_u_rsi(gkm,mkm," 
      | D_FTrsi_AZZZ_S -> "data3z_s_rsi(gkm,mkm,"
      | D_FTrsi_AZZZ_T -> "data3z_t_rsi(gkm,mkm,"
      | D_FTrsi_AZZZ_U -> "data3z_u_rsi(gkm,mkm,"      
      | D_FM0_ZZWW0_S -> "damzz0_s_0(gkm,mkm,"
      | D_FM0_ZZWW0_T -> "damzz0_t_0(gkm,mkm,"
      | D_FM0_ZZWW0_U -> "damzz0_u_0(gkm,mkm,"
      | D_FM0_ZZWW1_S -> "damzz1_s_0(gkm,mkm,"
      | D_FM0_ZZWW1_T -> "damzz1_t_0(gkm,mkm,"
      | D_FM0_ZZWW1_U -> "damzz1_u_0(gkm,mkm,"
      | D_FM0_WWWW0_S -> "damww0_s_0(gkm,mkm,"
      | D_FM0_WWWW0_T -> "damww0_t_0(gkm,mkm,"
      | D_FM0_WWWW0_U -> "damww0_u_0(gkm,mkm,"
      | D_FM0_WWWW2_S -> "damww2_s_0(gkm,mkm,"
      | D_FM0_WWWW2_T -> "damww2_t_0(gkm,mkm,"
      | D_FM0_WWWW2_U -> "damww2_u_0(gkm,mkm,"
      | D_FM0_ZZZZ_S  -> "damz4_s_0(gkm,mkm,"
      | D_FM0_ZZZZ_T  -> "damz4_t_0(gkm,mkm,"
      | D_FM0_ZZZZ_U  -> "damz4_u_0(gkm,mkm,"
      | D_FM1_ZZWW0_S -> "damzz0_s_1(gkm,mkm,"
      | D_FM1_ZZWW0_T -> "damzz0_t_1(gkm,mkm,"
      | D_FM1_ZZWW0_U -> "damzz0_u_1(gkm,mkm,"
      | D_FM1_ZZWW1_S -> "damzz1_s_1(gkm,mkm,"
      | D_FM1_ZZWW1_T -> "damzz1_t_1(gkm,mkm,"
      | D_FM1_ZZWW1_U -> "damzz1_u_1(gkm,mkm,"
      | D_FM1_WWWW0_S -> "damww0_s_1(gkm,mkm,"
      | D_FM1_WWWW0_T -> "damww0_t_1(gkm,mkm,"
      | D_FM1_WWWW0_U -> "damww0_u_1(gkm,mkm,"
      | D_FM1_WWWW2_S -> "damww2_s_1(gkm,mkm,"
      | D_FM1_WWWW2_T -> "damww2_t_1(gkm,mkm,"
      | D_FM1_WWWW2_U -> "damww2_u_1(gkm,mkm,"
      | D_FM1_ZZZZ_S  -> "damz4_s_1(gkm,mkm,"
      | D_FM1_ZZZZ_T  -> "damz4_t_1(gkm,mkm,"
      | D_FM1_ZZZZ_U  -> "damz4_u_1(gkm,mkm,"
      | D_FM7_ZZWW0_S -> "damzz0_s_7(gkm,mkm,"
      | D_FM7_ZZWW0_T -> "damzz0_t_7(gkm,mkm,"
      | D_FM7_ZZWW0_U -> "damzz0_u_7(gkm,mkm,"
      | D_FM7_ZZWW1_S -> "damzz1_s_7(gkm,mkm,"
      | D_FM7_ZZWW1_T -> "damzz1_t_7(gkm,mkm,"
      | D_FM7_ZZWW1_U -> "damzz1_u_7(gkm,mkm,"
      | D_FM7_WWWW0_S -> "damww0_s_7(gkm,mkm,"
      | D_FM7_WWWW0_T -> "damww0_t_7(gkm,mkm,"
      | D_FM7_WWWW0_U -> "damww0_u_7(gkm,mkm,"
      | D_FM7_WWWW2_S -> "damww2_s_7(gkm,mkm,"
      | D_FM7_WWWW2_T -> "damww2_t_7(gkm,mkm,"
      | D_FM7_WWWW2_U -> "damww2_u_7(gkm,mkm,"
      | D_FM7_ZZZZ_S  -> "damz4_s_7(gkm,mkm,"
      | D_FM7_ZZZZ_T  -> "damz4_t_7(gkm,mkm,"
      | D_FM7_ZZZZ_U  -> "damz4_u_7(gkm,mkm,"
      | D_Alpha_HHHH_S  -> "dalh4_s(gkm,mkm,"
      | D_Alpha_HHHH_T  -> "dalh4_t(gkm,mkm,"
      | D_Alpha_HHWW0_S -> "dalhw0_s(gkm,mkm,"
      | D_Alpha_HHWW0_T -> "dalhw0_t(gkm,mkm,"
      | D_Alpha_HHZZ0_S -> "dalhz0_s(gkm,mkm,"
      | D_Alpha_HHZZ0_T -> "dalhz0_t(gkm,mkm,"
      | D_Alpha_HHWW1_S -> "dalhw1_s(gkm,mkm,"
      | D_Alpha_HHWW1_T -> "dalhw1_t(gkm,mkm,"
      | D_Alpha_HHWW1_U -> "dalhw1_u(gkm,mkm,"
      | D_Alpha_HHZZ1_S -> "dalhz1_s(gkm,mkm,"
      | D_Alpha_HHZZ1_T -> "dalhz1_t(gkm,mkm,"
      | D_Alpha_HHZZ1_U -> "dalhz1_u(gkm,mkm,"
      | D_FM0_HHWW0_S -> "damhw0_s_0(gkm,mkm,"
      | D_FM0_HHWW0_T -> "damhw0_t_0(gkm,mkm,"
      | D_FM0_HHWW0_U -> "damhw0_u_0(gkm,mkm,"
      | D_FM0_HHZZ0_S -> "damhz0_s_0(gkm,mkm,"
      | D_FM0_HHZZ0_T -> "damhz0_t_0(gkm,mkm,"
      | D_FM0_HHZZ0_U -> "damhz0_u_0(gkm,mkm,"
      | D_FM0_HHWW1_S -> "damhw1_s_0(gkm,mkm,"
      | D_FM0_HHWW1_T -> "damhw1_t_0(gkm,mkm,"
      | D_FM0_HHWW1_U -> "damhw1_u_0(gkm,mkm,"
      | D_FM0_HHZZ1_S -> "damhz1_s_0(gkm,mkm,"
      | D_FM0_HHZZ1_T -> "damhz1_t_0(gkm,mkm,"
      | D_FM0_HHZZ1_U -> "damhz1_u_0(gkm,mkm,"   
      | D_FM1_HHWW0_S -> "damhw0_s_1(gkm,mkm,"
      | D_FM1_HHWW0_T -> "damhw0_t_1(gkm,mkm,"
      | D_FM1_HHWW0_U -> "damhw0_u_1(gkm,mkm,"
      | D_FM1_HHZZ0_S -> "damhz0_s_1(gkm,mkm,"
      | D_FM1_HHZZ0_T -> "damhz0_t_1(gkm,mkm,"
      | D_FM1_HHZZ0_U -> "damhz0_u_1(gkm,mkm,"
      | D_FM1_HHWW1_S -> "damhw1_s_1(gkm,mkm,"
      | D_FM1_HHWW1_T -> "damhw1_t_1(gkm,mkm,"
      | D_FM1_HHWW1_U -> "damhw1_u_1(gkm,mkm,"
      | D_FM1_HHZZ1_S -> "damhz1_s_1(gkm,mkm,"
      | D_FM1_HHZZ1_T -> "damhz1_t_1(gkm,mkm,"
      | D_FM1_HHZZ1_U -> "damhz1_u_1(gkm,mkm," 
      | D_FM7_HHWW0_S -> "damhw0_s_1(gkm,mkm,"
      | D_FM7_HHWW0_T -> "damhw0_t_1(gkm,mkm,"
      | D_FM7_HHWW0_U -> "damhw0_u_1(gkm,mkm,"
      | D_FM7_HHZZ0_S -> "damhz0_s_1(gkm,mkm,"
      | D_FM7_HHZZ0_T -> "damhz0_t_1(gkm,mkm,"
      | D_FM7_HHZZ0_U -> "damhz0_u_1(gkm,mkm,"
      | D_FM7_HHWW1_S -> "damhw1_s_1(gkm,mkm,"
      | D_FM7_HHWW1_T -> "damhw1_t_1(gkm,mkm,"
      | D_FM7_HHWW1_U -> "damhw1_u_1(gkm,mkm,"
      | D_FM7_HHZZ1_S -> "damhz1_s_1(gkm,mkm,"
      | D_FM7_HHZZ1_T -> "damhz1_t_1(gkm,mkm,"
      | D_FM7_HHZZ1_U -> "damhz1_u_1(gkm,mkm,"
      | G_HWW -> "ghww" | G_HZZ -> "ghzz"
      | G_HHWW -> "ghhww" | G_HHZZ -> "ghhzz"
      | G_SWW -> "gsww" | G_SZZ -> "gszz"
      | G_SHH -> "gshh"
      | G_SWW_T -> "gswwt" | G_SZZ_T -> "gszzt"
      | G_SAA_T -> "gsaat" | G_SAZ_T -> "gsazt"
      | G_PNWW -> "gpnww" | G_PNZZ -> "gpnzz"
      | G_PSNWW -> "gpsnww" | G_PSNZZ -> "gpsnzz"
      | G_PSNHH -> "gpsnhh"
      | G_PWZ -> "gpwz" | G_PWW -> "gpww"
      | G_FWW -> "gfww" | G_FZZ -> "gfzz"
      | G_FWW_CF -> "gfwwcf" | G_FZZ_CF -> "gfzzcf"
      | G_FHH -> "gfhh" | G_FHH_CF -> "gfhhcf"
      | G_FWW_T -> "gfwwt" | G_FZZ_T -> "gfzzt"
      | G_FFWW -> "gffww" | G_FFZZ -> "gffzz"
      | G_FFWW_CF -> "gffwwcf" | G_FFZZ_CF -> "gffzzcf"
      | G_FFHH -> "gffhh" | G_FFHH_CF -> "gffhhcf"
      | G_FVWW -> "gfvww" | G_FVZZ -> "gfvzz"
      | G_FVWW_CF -> "gfvwwcf" | G_FVZZ_CF -> "gfvzzcf"
      | G_FVHH -> "gfvhh" | G_FVHH_CF -> "gfvhhcf"
      | G_FDDSWW -> "gfddsww" | G_FDDSZZ -> "gfddszz"
      | G_FDDSWW_CF -> "gfddswwcf" | G_FDDSZZ_CF -> "gfddszzcf"
      | G_FDDSHH -> "gfddshh" | G_FDDSHH_CF -> "gfddshhcf"
      | G_FSWW -> "gfsww" | G_FSZZ -> "gfszz"
      | G_FSHH -> "gfshh"
      | G_TNWW -> "gtnww" | G_TNZZ -> "gtnzz"
      | G_TNWW_CF -> "gtnwwcf" | G_TNZZ_CF -> "gtnzzcf"
      | G_TSNWW -> "gtsnww" | G_TSNZZ -> "gtsnzz"
      | G_TSNWW_CF -> "gtsnwwcf" | G_TSNZZ_CF -> "gtsnzzcf"
      | G_TWZ -> "gtwz" | G_TWW -> "gtww"
      | G_TWZ_CF -> "gtwzcf" | G_TWW_CF -> "gtwwcf"
      | G_SSWW -> "gssww" | G_SSZZ -> "gsszz"
      | G_Htt -> "ghtt" | G_Hbb -> "ghbb"
      | G_Htautau -> "ghtautau" | G_Hcc -> "ghcc" | G_Hmm -> "ghmm"
      | G_HGaZ -> "ghgaz" | G_HGaGa -> "ghgaga" | G_Hgg -> "ghgg"
      | G_HGaGa_anom -> "ghgaga_ac" | G_HGaZ_anom -> "ghgaz_ac"
      | G_HZZ_anom -> "ghzz_ac" | G_HWW_anom -> "ghww_ac"
      | G_HGaZ_u -> "ghgaz_u" | G_HZZ_u -> "ghzz_u" 
      | G_HWW_u -> "ghww_u" 
      | G_H3 -> "gh3" | G_H4 -> "gh4"
      | Gs -> "gs" | I_Gs -> "igs" | G2 -> "gs**2"
      | Mass f -> "mass" ^ flavor_symbol f
      | Width f -> "width" ^ flavor_symbol f
      | K_Matrix_Coeff i -> "kc" ^ string_of_int i
      | K_Matrix_Pole i -> "kp" ^ string_of_int i

  end
