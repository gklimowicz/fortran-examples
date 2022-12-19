(* modellib_PSSSM.ml --

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


(* \thocwmodulesection{Extended Supersymmetric Standard Model(s)} *)

(* This is based on the NMSSM implementation by Felix Braam, and extended to
   the exotica -- leptoquarks, leptoquarkinos, additional Higgses etc. -- by
   Daniel Wiesler. Note that for the Higgs sector vertices the conventions of
   the Franke/Fraas paper have been used. *)

module type extMSSM_flags =
  sig
    val ckm_present       : bool
  end

module PSSSM : extMSSM_flags =
  struct 
    let ckm_present       = false
  end

module PSSSM_QCD : extMSSM_flags =
  struct 
    let ckm_present       = false
  end

module ExtMSSM (Flags : extMSSM_flags) = 
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

(*additional combinatorics *)
(* yields a list of tuples consistig of the off-diag combinations of the elements in "set" *)
    let choose2 set =
      List.map (function [x;y] -> (x,y) | _ -> failwith "choose2")
        (Combinatorics.choose 2 set)


(* [pairs] *)
(* [pairs] appends the diagonal combinations to [choose2] *)    	
    let rec diag = function
      | [] -> []
      | x1 :: rest -> (x1, x1) :: diag rest

    let pairs l = choose2 l @ diag l

(* [triples] *)
(* rank 3 generalization of [pairs] *)
   let rec cloop set i j k =
     if i > ((List.length set)-1) then []
     else if j > i then cloop set (succ i) (j-i-1) (j-i-1)    
     else if k > j then cloop set i (succ j) (k-j-1)  
     else (List.nth set i, List.nth set j, List.nth set k) :: cloop set i j (succ k)

    let triples set = cloop set 0 0 0
(* [two_and_one] *)

    let rec two_and_one' l1 z n =
       if n < 0 then []
       else
       ((fst (List.nth (pairs l1) n)),(snd (List.nth (pairs l1) n)), z):: two_and_one' l1 z (pred n) 

    let two_and_one l1 l2 = 
       let f z = two_and_one' l1 z ((List.length (pairs l1))-1)
       in
       List.flatten ( List.map f l2 ) 


    type gen = 
      | G of int | GG of gen*gen

    let rec string_of_gen = function
      | G n when n > 0  -> string_of_int n
      | G n -> string_of_int (abs n) ^ "c" 
      | GG (g1,g2) -> string_of_gen g1 ^ "_" ^ string_of_gen g2

(* With this we distinguish the flavour. *)

    type sff = 
      | SL | SN | SU | SD

    let string_of_sff = function
      | SL -> "sl" | SN -> "sn" | SU -> "su" | SD -> "sd"         

(* With this we distinguish the mass eigenstates. At the moment we have to cheat 
   a little bit for the sneutrinos. Because we are dealing with massless 
   neutrinos there is only one sort of sneutrino. *)

    type sfm =
      | M1 | M2

    let string_of_sfm = function 
      | M1 -> "1" | M2 -> "2"

(* We also introduce special types for the charginos and neutralinos. *)

    type char = 
      | C1 | C2 | C1c | C2c | C3 | C3c | C4 | C4c

    type neu =
      | N1 | N2 | N3 | N4 | N5 | N6 | N7 | N8 | N9 | N10 | N11 

    let int_of_char = function
      | C1 -> 1 | C2 -> 2 | C1c -> -1 | C2c -> -2
      | C3 -> 3 | C4 -> 4 | C3c -> -3 | C4c -> -4

    let string_of_char c = string_of_int (int_of_char c)

    let conj_char = function
      | C1 -> C1c | C2 -> C2c | C1c -> C1 | C2c -> C2
      | C3 -> C3c | C4 -> C4c | C3c -> C3 | C4c -> C4

    let string_of_neu = function
      | N1 -> "1" | N2 -> "2" | N3 -> "3" | N4 -> "4" | N5 -> "5" | N6 -> "6" 
      | N7 -> "7" | N8 -> "8" | N9 -> "9" | N10 -> "10"| N11 -> "11"

(* For NMSSM-like the Higgs bosons, we follow the conventions of 
   Franke/Fraas. Daniel Wiesler: extended to E6 models. *)

    type shiggs =
      | S1 | S2 | S3 | S4 | S5 | S6 | S7 | S8 | S9

    type phiggs =
      | P1 | P2 | P3 | P4 | P5 | P6 | P7 

(* [HCx] is always the $H^+$, [HCxc] the $H^-$. *)

    type chiggs = 
      | HC1 | HC2 | HC3 | HC4 | HC5 | HC1c | HC2c | HC3c | HC4c | HC5c 

    let conj_chiggs = function
      | HC1 -> HC1c | HC2 -> HC2c | HC1c -> HC1 | HC2c -> HC2
      | HC3 -> HC3c | HC4 -> HC4c | HC3c -> HC3 | HC4c -> HC4
      | HC5 -> HC5c | HC5c -> HC5

    let string_of_shiggs = function
      | S1 -> "1" | S2 -> "2" | S3 -> "3" | S4 -> "4" | S5 -> "5" 
      | S6 -> "6" | S7 -> "7" | S8 -> "8" | S9 -> "9"

    let string_of_phiggs = function
      | P1 -> "1" | P2 -> "2" | P3 -> "3" | P4 -> "4" | P5 -> "5" 
      | P6 -> "6" | P7 -> "7" 

    let nlist = [ N1; N2; N3; N4; N5; N6; N7; N8; N9; N10; N11 ]
    let slist = [ S1; S2; S3; S4; S5; S6; S7; S8; S9 ]
    let plist = [ P1; P2; P3; P4; P5; P6; P7 ]
    let clist = [ HC1; HC2; HC3; HC4; HC5; HC1c; HC2c; HC3c; HC4c; HC5c ]
    let charlist = [ C1; C2; C3; C4; C1c; C2c; C3c; C4c ]

    type flavor =
      | L of int | N of int
      | U of int | D of int
      | Sup of sfm*int | Sdown of sfm*int 
      | Ga | Wp | Wm | Z | Gl 
      | Slepton of sfm*int | Sneutrino of int 
      | Neutralino of neu | Chargino of char 
      | Gluino | SHiggs of shiggs | PHiggs of phiggs 
      | CHiggs of chiggs
      |	LQ of sfm*int 
      | LQino of int 
      
    let string_of_fermion_type = function
      | L _ -> "l" | U _ -> "u" | D _ -> "d" | N _ -> "n"
      | _ -> failwith 
            "Modellib_PSSSM.ExtMSSM.string_of_fermion_type: invalid fermion type"

    let string_of_fermion_gen = function
      | L g | U g | D g | N g -> string_of_int (abs (g))
      | _ -> failwith 
            "Modellib_PSSSM.ExtMSSM.string_of_fermion_gen: invalid fermion type"
            
    type gauge = unit

    let gauge_symbol () =
      failwith "Modellib_PSSSM.ExtMSSM.gauge_symbol: internal error"       

(* At this point we will forget graviton and -ino. *) 

    let family g = [ L g; N g; Slepton (M1,g); 
                     Slepton (M2,g); Sneutrino g;
                     U g; D g; Sup (M1,g); Sup (M2,g);
                     Sdown (M1,g); Sdown (M2,g); 
                     LQ (M1,g); LQ (M2,g); LQino g ]

    let external_flavors () = 
        [ "1st Generation matter", ThoList.flatmap family [1; -1];
          "2nd Generation matter", ThoList.flatmap family [2; -2];
          "3rd Generation matter", ThoList.flatmap family [3; -3];
          "Gauge Bosons", [Ga; Z; Wp; Wm; Gl];
          "Charginos", List.map (fun a -> Chargino a) charlist;
          "Neutralinos", List.map (fun a -> Neutralino a) nlist;     
          "Higgs Bosons", List.map (fun a -> SHiggs a) slist @ 
                          List.map (fun a -> PHiggs a) plist @
                          List.map (fun a -> CHiggs a) clist;
          "Gluino", [Gluino]]
          
    let flavors () = ThoList.flatmap snd (external_flavors ())

    let spinor n m =
      if n >= 0 && m >= 0 then
        Spinor
      else if
        n <= 0 && m <=0 then
        ConjSpinor
      else
        invalid_arg "Modellib_PSSSM.ExtMSSM.spinor: internal error"

    let lorentz = function
      | L g -> spinor g 0 | N g -> spinor g 0
      | U g -> spinor g 0 | D g -> spinor g 0 
      | LQino g -> spinor g 0
      | Chargino c -> spinor (int_of_char c) 0 
      | Ga | Gl -> Vector
      | Wp | Wm | Z -> Massive_Vector
      | SHiggs _ | PHiggs _ | CHiggs _  
      | Sup _ | Sdown _ | Slepton _ | Sneutrino _ | LQ _ -> Scalar 
      | Neutralino _ | Gluino -> Majorana 

    let color = function
      | U g -> Color.SUN (if g > 0 then 3 else -3)
      | Sup (m,g) -> Color.SUN (if g > 0 then 3 else -3)
      | D g -> Color.SUN (if g > 0 then 3 else -3)
      | Sdown (m,g) -> Color.SUN  (if g > 0 then 3 else -3)
      | LQ (m,g) -> Color.SUN  (if g > 0 then 3 else -3)
      | LQino g -> Color.SUN  (if g > 0 then 3 else -3)
      | Gl | Gluino -> Color.AdjSUN 3
      | _ -> Color.Singlet   

    let nc () = 3

    let prop_spinor n m =
      if n >= 0 && m >=0 then
        Prop_Spinor
      else if 
        n <=0 && m <=0 then
        Prop_ConjSpinor
      else 
        invalid_arg "Modellib_PSSSM.ExtMSSM.prop_spinor: internal error"

    let propagator = function
      | L g -> prop_spinor g 0 | N g -> prop_spinor g 0
      | U g -> prop_spinor g 0 | D g -> prop_spinor g 0
      | LQino g -> prop_spinor g 0
      | Chargino c -> prop_spinor (int_of_char c) 0 
      | Ga | Gl -> Prop_Feynman
      | Wp | Wm | Z -> Prop_Unitarity
      | SHiggs _ | PHiggs _ | CHiggs _ -> Prop_Scalar
      | Sup _ | Sdown _ | Slepton _ | Sneutrino _ -> Prop_Scalar
      | LQ _ -> Prop_Scalar
      | Gluino -> Prop_Majorana 
      | Neutralino _ -> Prop_Majorana

(* Optionally, ask for the fudge factor treatment for the widths of
   charged particles.  Currently, this only applies to $W^\pm$ and top. *)
                                                                     
    let width f =
      if !use_fudged_width then
        match f with
        | Wp | Wm | U 3 | U (-3) -> Fudged
        | _ -> !default_width
      else
        !default_width 

    let goldstone _ = None

    let conjugate = function
      | L g -> L (-g) | N g -> N (-g)
      | U g -> U (-g) | D g -> D (-g)
      | Sup (m,g) -> Sup (m,-g) 
      | Sdown (m,g) -> Sdown (m,-g) 
      | Slepton (m,g) -> Slepton (m,-g) 
      | Sneutrino g -> Sneutrino (-g)
      | Gl -> Gl | Ga -> Ga | Z -> Z
      | Wp -> Wm | Wm -> Wp
      | SHiggs s -> SHiggs s 
      | PHiggs p -> PHiggs p 
      | CHiggs c -> CHiggs (conj_chiggs c)
      | Gluino -> Gluino 
      | Neutralino n -> Neutralino n | Chargino c -> Chargino (conj_char c)
      | LQino g -> LQino (-g)
      | LQ (m,g) -> LQ (m,-g)

   let fermion = function
     | L g -> if g > 0 then 1 else -1
     | N g -> if g > 0 then 1 else -1
     | U g -> if g > 0 then 1 else -1
     | D g -> if g > 0 then 1 else -1
     | Gl | Ga | Z | Wp | Wm -> 0 
     | SHiggs _ | PHiggs _ | CHiggs _ -> 0       
     | Neutralino _ -> 2
     | Chargino c -> if (int_of_char c) > 0 then 1 else -1
     | Sup _ -> 0 | Sdown _ -> 0 
     | Slepton _ -> 0 | Sneutrino _ -> 0          
     | Gluino -> 2 
     | LQ _ -> 0
     | LQino g -> if g > 0 then 1 else -1

(* This model does NOT have a conserved generation quantum number. *)

    module Ch = Charges.QQ

    let ( // ) = Algebra.Small_Rational.make

    let charge = function
      | L n -> if n > 0 then -1//1 else  1//1
      | Slepton (_,n) -> if n > 0 then -1//1 else  1//1
      | N n -> 0//1
      | Sneutrino n -> 0//1
      | U n -> if n > 0 then  2//3 else -2//3
      | Sup (_,n) -> if n > 0 then  2//3 else -2//3
      | D n | LQ (_,n) | LQino n -> if n > 0 then -1//3 else  1//3          
      | Sdown (_,n) -> if n > 0 then -1//3 else  1//3          
      | Gl | Ga | Z | Neutralino _ | Gluino -> 0//1
      | Wp ->  1//1
      | Wm -> -1//1
      | SHiggs _ | PHiggs _  ->  0//1
      | CHiggs (HC1|HC2|HC3|HC4|HC5) ->  1//1
      | CHiggs (HC1c|HC2c|HC3c|HC4c|HC5c) -> -1//1
      | Chargino (C1|C2|C3|C4) -> 1//1 
      | Chargino (C1c|C2c|C3c|C4c) -> -1//1 

    let lepton = function
      | L n | N n -> if n > 0 then 1//1 else -1//1
      | Slepton (_,n) 
      | Sneutrino n -> if n > 0 then 1//1 else -1//1
      | LQ (_,n) | LQino n -> if n > 0 then 1//1 else -1//1
      | _ -> 0//1

    let baryon = function
      | U n | D n -> if n > 0 then 1//1 else -1//1
      | Sup (_,n) | Sdown (_,n) -> if n > 0 then 1//1 else -1//1
      | LQ (_,n) | LQino n -> if n > 0 then 1//1 else -1//1
      | _ -> 0//1

    let charges f =
      [ charge f; lepton f; baryon f] 

(* We introduce a Boolean type vc as a pseudonym for Vertex Conjugator to 
   distinguish between vertices containing complex mixing matrices like the 
   CKM--matrix or the sfermion or neutralino/chargino--mixing matrices, which 
   have to become complex conjugated. The true--option stands for the conjugated 
   vertex, the false--option for the unconjugated vertex. *)

    type vc = bool

    type constant =
      | E | G 
      | Q_lepton | Q_up | Q_down | Q_charg           
      | G_Z | G_CC | G_CCQ of vc*int*int
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down 
      | I_Q_W | I_G_ZWW | G_WWWW | G_ZZWW | G_PZWW | G_PPWW 
      | G_strong | G_SS | I_G_S 
      | Gs
      | G_NZN of neu*neu | G_CZC of char*char 
      | G_YUK_FFS of flavor*flavor*shiggs
      | G_YUK_FFP of flavor*flavor*phiggs
      | G_YUK_LCN of int
      | G_YUK_UCD of int*int | G_YUK_DCU of int*int 
      | G_NHC of vc*neu*char 
      | G_YUK_C of vc*flavor*char*sff*sfm
      | G_YUK_Q of vc*int*flavor*char*sff*sfm
      | G_YUK_N of vc*flavor*neu*sff*sfm
      | G_YUK_G of vc*flavor*sff*sfm
      | G_NWC of neu*char | G_CWN of char*neu
      | G_CSC of char*char*shiggs	
      | G_CPC of char*char*phiggs	
      | G_WSQ of vc*int*int*sfm*sfm
      | G_SLSNW of vc*int*sfm 
      | G_ZSF of sff*int*sfm*sfm
      | G_CICIS of neu*neu*shiggs
      | G_CICIP of neu*neu*phiggs
      | G_GH_WPC of phiggs   | G_GH_WSC of shiggs
      | G_GH_ZSP of shiggs*phiggs | G_GH_WWS of shiggs
      | G_GH_ZZS of shiggs | G_GH_ZCC 
      | G_GH_GaCC  
      | G_GH4_ZZPP of phiggs*phiggs
      | G_GH4_ZZSS of shiggs*shiggs
      | G_GH4_ZZCC  | G_GH4_GaGaCC
      | G_GH4_ZGaCC | G_GH4_WWCC
      | G_GH4_WWPP of phiggs*phiggs
      | G_GH4_WWSS of shiggs*shiggs
      | G_GH4_ZWSC of shiggs
      | G_GH4_GaWSC of shiggs
      | G_GH4_ZWPC of phiggs
      | G_GH4_GaWPC of phiggs
      | G_WWSFSF of sff*int*sfm*sfm 
      | G_WPSLSN of vc*int*sfm
      | G_H3_SCC of shiggs
      | G_H3_SSS of shiggs*shiggs*shiggs
      | G_H3_SPP of shiggs*phiggs*phiggs
      | G_SFSFS of shiggs*sff*int*sfm*sfm
      | G_SFSFP of phiggs*sff*int*sfm*sfm
      | G_HSNSL of vc*int*sfm  
      | G_HSUSD of vc*sfm*sfm*int*int 
      | G_WPSUSD of vc*sfm*sfm*int*int  
      | G_WZSUSD of vc*sfm*sfm*int*int  
      | G_WZSLSN of vc*int*sfm | G_GlGlSQSQ
      | G_PPSFSF of sff 
      | G_ZZSFSF of sff*int*sfm*sfm | G_ZPSFSF of sff*int*sfm*sfm 
      | G_GlZSFSF of sff*int*sfm*sfm | G_GlPSQSQ 
      | G_GlWSUSD of vc*sfm*sfm*int*int
      | G_YUK_LQ_S of int*shiggs*int
      | G_YUK_LQ_P of int*phiggs*int
      | G_LQ_NEU of sfm*int*int*neu
      | G_LQ_EC_UC of vc*sfm*int*int*int
      | G_LQ_GG of sfm*int*int
      | G_LQ_SSU of sfm*sfm*sfm*int*int*int
      | G_LQ_SSD of sfm*sfm*int*int*int
      | G_LQ_S of sfm*sfm*int*shiggs*int
      | G_LQ_P of sfm*sfm*int*phiggs*int
      | G_ZLQ of int*sfm*sfm
      | G_ZZLQLQ | G_ZPLQLQ | G_PPLQLQ | G_ZGlLQLQ | G_PGlLQLQ | G_NLQC | G_GlGlLQLQ

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

(* \begin{subequations}
     \begin{align}
        \alpha_{\text{QED}} &= \frac{1}{137.0359895} \\
             \sin^2\theta_w &= 0.23124
     \end{align}
   \end{subequations}

Here we must perhaps allow for complex input parameters. So split them
into their modulus and their phase. At first, we leave them real; the 
generalization to complex parameters is obvious. *)


    let parameters () =
      { input = [];
        derived = [];
        derived_arrays = [] }   
      
    module F = Modeltools.Fusions (struct
      type f = flavor
      type c = constant
      let compare = compare
      let conjugate = conjugate
    end)


(* For the couplings there are generally two possibilities concerning the
   sign of the covariant derivative. 
   \begin{equation} 
   {\rm CD}^\pm = \partial_\mu \pm \ii g T^a A^a_\mu 
   \end{equation} 
   The particle data group defines the signs consistently to be positive. 
   Since the convention for that signs also influence the phase definitions 
   of the gaugino/higgsino fields via the off-diagonal entries in their
   mass matrices it would be the best to adopt that convention. *)

(*** REVISED: Compatible with CD+.  FB ***)
    let electromagnetic_currents_3 g =
        [ ((L (-g), Ga, L g), FBF (1, Psibar, V, Psi), Q_lepton);
          ((U (-g), Ga, U g), FBF (1, Psibar, V, Psi), Q_up);
          ((D (-g), Ga, D g), FBF (1, Psibar, V, Psi), Q_down)]
        
(*** REVISED: Compatible with CD+. FB***)
    let electromagnetic_sfermion_currents g m =
        [ ((Ga, Slepton (m,-g), Slepton (m,g)), Vector_Scalar_Scalar 1, Q_lepton);
          ((Ga, Sup (m,-g), Sup (m,g)), Vector_Scalar_Scalar 1, Q_up);
          ((Ga, Sdown (m,-g), Sdown (m,g)), Vector_Scalar_Scalar 1, Q_down)]       

(*** REVISED: Compatible with CD+. FB***)
    let electromagnetic_currents_2 c =
      let cc = conj_char c in
      [ ((Chargino cc, Ga, Chargino c), FBF (1, Psibar, V, Psi), Q_charg) ]

(*** REVISED: Compatible with CD+. FB***)
    let neutral_currents g =
      [ ((L (-g), Z, L g), FBF (1, Psibar, VA, Psi), G_NC_lepton);
        ((N (-g), Z, N g), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
        ((U (-g), Z, U g), FBF (1, Psibar, VA, Psi), G_NC_up);
        ((D (-g), Z, D g), FBF (1, Psibar, VA, Psi), G_NC_down)]

(* \begin{equation}
     \mathcal{L}_{\textrm{CC}} =
        \mp \frac{g}{2\sqrt2} \sum_i \bar\psi_i \gamma^\mu
               (1-\gamma_5)(T^+W^+_\mu+T^-W^-_\mu)\psi_i ,
   \end{equation}
   where the sign corresponds to $\text{CD}_\pm$, respectively.  *)

(*** REVISED: Compatible with CD+. ***)
        (* Remark: The definition with the other sign compared to the SM files
           comes from the fact that $g_{cc} = 1/(2\sqrt{2})$ is used 
           overwhelmingly often in the SUSY Feynman rules, so that JR 
           decided to use a different definiton for [g_cc] in SM and MSSM. *)
(**    FB         **)
    let charged_currents g =
      [ ((L (-g), Wm, N g), FBF ((-1), Psibar, VL, Psi), G_CC);
        ((N (-g), Wp, L g), FBF ((-1), Psibar, VL, Psi), G_CC) ]

(* The quark with the inverted generation (the antiparticle) is the outgoing 
   one, the other the incoming. The vertex attached to the outgoing up-quark 
   contains the CKM matrix element {\em not} complex conjugated, while the 
   vertex with the outgoing down-quark has the conjugated CKM matrix 
   element. *)

(*** REVISED: Compatible with CD+. FB ***)
    let charged_quark_currents g h = 
        [ ((D (-g), Wm, U h), FBF ((-1), Psibar, VL, Psi), G_CCQ (true,g,h));
          ((U (-g), Wp, D h), FBF ((-1), Psibar, VL, Psi), G_CCQ (false,h,g))] 

(*** REVISED: Compatible with CD+.FB ***)
    let charged_chargino_currents n c =
      let cc = conj_char c in 
      [ ((Chargino cc, Wp, Neutralino n), 
                    FBF (1, Psibar, VLR, Chi), G_CWN (c,n));
        ((Neutralino n, Wm, Chargino c), 
                    FBF (1, Chibar, VLR, Psi), G_NWC (n,c)) ]

(*** REVISED: Compatible with CD+. FB***)
    let charged_slepton_currents g m =
      [ ((Wm, Slepton (m,-g), Sneutrino g), Vector_Scalar_Scalar (-1), G_SLSNW 
           (true,g,m));
        ((Wp, Slepton (m,g), Sneutrino (-g)), Vector_Scalar_Scalar 1, G_SLSNW 
           (false,g,m)) ]
 
(*** REVISED: Compatible with CD+. FB***)
    let charged_squark_currents' g h m1 m2 =
      [ ((Wm, Sup (m1,g), Sdown (m2,-h)), Vector_Scalar_Scalar (-1), G_WSQ 
             (true,g,h,m1,m2));
          ((Wp, Sup (m1,-g), Sdown (m2,h)), Vector_Scalar_Scalar 1, G_WSQ 
             (false,g,h,m1,m2)) ]
    let charged_squark_currents g h = 
    List.flatten (Product.list2 (charged_squark_currents' g h) [M1;M2] [M1;M2] ) 

(*** REVISED: Compatible with CD+. FB ***)
    let neutral_sfermion_currents' g m1 m2 =
      [ ((Z, Slepton (m1,-g), Slepton (m2,g)), Vector_Scalar_Scalar (-1), 
            G_ZSF (SL,g,m1,m2));
        ((Z, Sup (m1,-g), Sup (m2,g)), Vector_Scalar_Scalar (-1), 
            G_ZSF (SU,g,m1,m2));
        ((Z, Sdown (m1,-g), Sdown (m2,g)), Vector_Scalar_Scalar (-1), 
            G_ZSF (SD,g,m1,m2))]
    let neutral_sfermion_currents g = 
      List.flatten (Product.list2 (neutral_sfermion_currents'
                  g) [M1;M2] [M1;M2]) @
      [ ((Z, Sneutrino (-g), Sneutrino g), Vector_Scalar_Scalar (-1), 
            G_ZSF (SN,g,M1,M1)) ]

(* The reality of the coupling of the Z-boson to two identical neutralinos 
   makes the vector part of the coupling vanish. So we distinguish them not 
   by the name but by the structure of the couplings. *)  

(*** REVISED: Compatible with CD+. FB***)
    let neutral_Z (n,m) =  
      [ ((Neutralino n, Z, Neutralino m), FBF (1, Chibar, VA, Chi), 
              (G_NZN (n,m))) ]

(*** REVISED: Compatible with CD+. FB***)
    let charged_Z c1 c2 =
      let cc1 = conj_char c1 in
      ((Chargino cc1, Z, Chargino c2), FBF ((-1), Psibar, VA , Psi), 
               G_CZC (c1,c2)) 

(*** REVISED: Compatible with CD+. 
   Remark: This is pure octet. FB***)        
    
    let yukawa_v =
      [ (Gluino, Gl, Gluino), FBF (1, Chibar, V, Chi), Gs]

(*** REVISED: Independent of the sign of CD. ***)
(*** REVISED: Felix Braam: Compact version using new COMBOS + FF-Couplings *)
    let yukawa_higgs_FFS f s   = 
        [((conjugate f, SHiggs s, f ), FBF (1, Psibar, S, Psi), 
           G_YUK_FFS (conjugate f, f, s))]          
    let yukawa_higgs_FFP f p   =  
        [((conjugate f, PHiggs p, f), FBF (1, Psibar, P, Psi), 
           G_YUK_FFP (conjugate f ,f , p))] 

(* JR: Only the first charged Higgs. *)
    let yukawa_higgs_NLC g = 
      [ ((N (-g), CHiggs HC1, L g), FBF (1, Psibar, Coupling.SR, Psi), 
            G_YUK_LCN g);
        ((L (-g), CHiggs HC1c, N g), FBF (1, Psibar, Coupling.SL, Psi), 
            G_YUK_LCN g)]
    
    let yukawa_higgs g = 
       yukawa_higgs_NLC g @
       List.flatten ( Product.list2 yukawa_higgs_FFS  [L g; U g; D g] [S1; S2; S3]) @ 
       List.flatten ( Product.list2 yukawa_higgs_FFP  [L g; U g; D g] [P1; P2]) 


(* JR: Only the first charged Higgs. *)   
(*** REVISED: Independent of the sign of CD. FB***)
    let yukawa_higgs_quark (g,h) =
      [ ((U (-g), CHiggs HC1, D h), FBF (1, Psibar, SLR, Psi), 
            G_YUK_UCD (g, h)); 
        ((D (-h), CHiggs HC1c, U g), FBF (1, Psibar, SLR, Psi), 
            G_YUK_DCU (g, h))  ]

(*** REVISED: Compatible with CD+.FB*)
(*** REVISED: Compact version using new COMBOS*)
    let yukawa_shiggs_2 c1 c2 s =
      let cc1 = conj_char c1 in
       ((Chargino cc1, SHiggs s, Chargino c2), FBF (1, Psibar, SLR, Psi), 
           G_CSC (c1,c2,s))  

    let yukawa_phiggs_2 c1 c2 p =
      let cc1 = conj_char c1 in
      ((Chargino cc1, PHiggs p, Chargino c2), FBF (1, Psibar, SLR, Psi), 
           G_CPC (c1,c2,p))  

    let yukawa_higgs_2 = 
      Product.list3 yukawa_shiggs_2 [C1;C2] [C1;C2] [S1;S2;S3] @ 
      Product.list3 yukawa_phiggs_2 [C1;C2] [C1;C2] [P1;P2] 

(* JR: Only the first charged Higgs. *)
(*** REVISED: Compatible with CD+.FB ***)
    let higgs_charg_neutr n c =
      let cc = conj_char c in
      [ ((Neutralino n, CHiggs HC1c, Chargino c), FBF (-1, Chibar, SLR, Psi), 
                   G_NHC (false,n,c));
        ((Chargino cc, CHiggs HC1, Neutralino n), FBF (-1, Psibar, SLR, Chi), 
                   G_NHC (true,n,c)) ]

(*** REVISED: Compatible with CD+. FB***)    
(*** REVISED: Compact version using new COMBOS*)    
    let shiggs_neutr (n,m,s)  =
       ((Neutralino n, SHiggs s, Neutralino m), FBF (1, Chibar, SP, Chi), 
           G_CICIS (n,m,s)) 
    let phiggs_neutr (n,m,p) =
       ((Neutralino n, PHiggs p, Neutralino m), FBF (1, Chibar, SP, Chi), 
           G_CICIP (n,m,p)) 
    
    let higgs_neutr = 						
      List.map shiggs_neutr (two_and_one [N1;N2;N3;N4;N5] [S1;S2;S3]) @ 
      List.map phiggs_neutr (two_and_one [N1;N2;N3;N4;N5] [P1;P2]) 

(*** REVISED: Compatible with CD+. FB***)
       let yukawa_n_2 n m g = 
         [ ((Neutralino n, Slepton (m,-g), L g), FBF (1, Chibar, SLR, Psi),  
                G_YUK_N (true,L g,n,SL,m));
           ((L (-g), Slepton (m,g), Neutralino n), FBF (1, Psibar, SLR, Chi), 
                G_YUK_N (false,L g,n,SL,m));
           ((Neutralino n, Sup (m,-g), U g), FBF (1, Chibar, SLR, Psi), 
                G_YUK_N (true,U g,n,SU,m));
           ((U (-g), Sup (m,g), Neutralino n), FBF (1, Psibar, SLR, Chi), 
                G_YUK_N (false,U g,n,SU,m));
           ((Neutralino n, Sdown (m,-g), D g), FBF (1, Chibar, SLR, Psi), 
                G_YUK_N (true,D g,n,SD,m));
           ((D (-g), Sdown (m,g), Neutralino n), FBF (1, Psibar, SLR, Chi), 
                G_YUK_N (false,D g,n,SD,m)) ]
     let yukawa_n_3 n g =
         [ ((Neutralino n, Sneutrino (-g), N g), FBF (1, Chibar, SLR, Psi), 
                G_YUK_N (true,N g,n,SN,M1));
           ((N (-g), Sneutrino g, Neutralino n), FBF (1, Psibar, SLR, Chi), 
                G_YUK_N (false,N g, n,SN,M1)) ]

    let yukawa_n_5 g m =
          [ ((U (-g), Sup (m,g), Gluino), FBF (1, Psibar, SLR, Chi), 
                G_YUK_G (false,U g,SU,m));
           ((D (-g), Sdown (m,g), Gluino), FBF (1, Psibar, SLR, Chi), 
                G_YUK_G (false,D g,SD,m));
           ((Gluino, Sup (m,-g), U g), FBF (1, Chibar, SLR, Psi), 
                G_YUK_G (true,U g,SU,m));
           ((Gluino, Sdown (m,-g), D g), FBF (1, Chibar, SLR, Psi), 
                G_YUK_G (true,D g,SD,m))]
    let yukawa_n =
      List.flatten (Product.list3 yukawa_n_2 [N1;N2;N3;N4;N5] [M1;M2] [1;2;3]) @
      List.flatten (Product.list2 yukawa_n_3 [N1;N2;N3;N4;N5] [1;2;3]) @
      List.flatten (Product.list2 yukawa_n_5 [1;2;3] [M1;M2]) 
      

(*** REVISED: Compatible with CD+.FB ***)
    let yukawa_c_2 c g  = 
         let cc = conj_char c in
         [ ((L (-g), Sneutrino g, Chargino cc), BBB (1, Psibar, SLR, 
              Psibar), G_YUK_C (true,L g,c,SN,M1));
           ((Chargino c, Sneutrino (-g), L g), PBP (1, Psi, SLR, Psi), 
              G_YUK_C (false,L g,c,SN,M1)) ]
    let yukawa_c_3 c m g =
         let cc = conj_char c in
         [ ((N (-g), Slepton (m,g), Chargino c), FBF (1, Psibar, SLR, 
              Psi), G_YUK_C (true,N g,c,SL,m));
           ((Chargino cc, Slepton (m,-g), N g), FBF (1, Psibar, SLR, 
              Psi), G_YUK_C (false,N g,c,SL,m)) ]
    let yukawa_c c = 
      ThoList.flatmap (yukawa_c_2 c) [1;2;3] @ 
      List.flatten (Product.list2 (yukawa_c_3 c) [M1;M2] [1;2;3]) 


(*** REVISED: Compatible with CD+. FB***)
   let yukawa_cq' c (g,h) m = 
       let cc = conj_char c in
         [ ((Chargino c, Sup (m,-g), D h), PBP (1, Psi, SLR, Psi), 
            G_YUK_Q (false,g,D h,c,SU,m));
           ((D (-h), Sup (m,g), Chargino cc), BBB (1, Psibar, SLR, Psibar), 
            G_YUK_Q (true,g,D h,c,SU,m));
           ((Chargino cc, Sdown (m,-g), U h), FBF (1, Psibar, SLR, Psi), 
            G_YUK_Q (true,g,U h,c,SD,m));
           ((U (-h), Sdown (m,g), Chargino c), FBF (1, Psibar, SLR, Psi), 
            G_YUK_Q (false,g,U h,c,SD,m)) ]               
    let yukawa_cq c =      
     if Flags.ckm_present then
       List.flatten (Product.list2 (yukawa_cq' c) [(1,1);(1,2);(2,1);(2,2);(1,3);(2,3);(3,3);(3,2);(3,1)] [M1;M2]) 
     else
       List.flatten (Product.list2 (yukawa_cq' c) [(1,1);(2,2);(3,3)] [M1;M2]) 


(*** REVISED: Compatible with CD+. 
   Remark: Singlet and octet gluon exchange. The coupling is divided by
   sqrt(2) to account for the correct normalization of the Lie algebra
   generators.
**FB*)         
    let col_currents g =
      [ ((D (-g), Gl, D g), FBF ((-1), Psibar, V, Psi), Gs);
        ((U (-g), Gl, U g), FBF ((-1), Psibar, V, Psi), Gs)]

(*** REVISED: Compatible with CD+. 
   Remark: Singlet and octet gluon exchange. The coupling is divided by
   sqrt(2) to account for the correct normalization of the Lie algebra
   generators.
**FB*)

(** LQ-coupl.  **DW**)
   

   let chg = function
     | M1 -> M2 | M2 -> M1


(** LQ - Yuk's **)


   let yuk_lqino_se_uc1' g1 g2 g3 m =
     let cm = chg m in
       [ ((U (-g3), Slepton (m,-g2), LQino g1), FBF (1, Psibar, SLR, Psi), 
             G_LQ_EC_UC (true,cm,g1,g2,g3)) ]

   let yuk_lqino_se_uc1 g1 g2 g3 =
      ThoList.flatmap (yuk_lqino_se_uc1' g1 g2 g3) [M1;M2] 

   let yuk_lqino_se_uc2' g1 g2 g3 m =
     let cm = chg m in 
       [ ((LQino (-g1), Slepton (m,g2), U g3), FBF (1, Psibar, SLR, Psi), 
             G_LQ_EC_UC (false,cm,g1,g2,g3)) ]

   let yuk_lqino_se_uc2 g1 g2 g3 =
      ThoList.flatmap (yuk_lqino_se_uc2' g1 g2 g3) [M1;M2] 


   let yuk_lqino_sn_dc1 g1 g2 g3 =
     [ ((D (-g3), Sneutrino (-g2), LQino g1), FBF (-1, Psibar, SLR, Psi), 
           G_LQ_EC_UC (true,M2,g1,g2,g3)) ]

   let yuk_lqino_sn_dc2 g1 g2 g3 =
     [ ((LQino (-g1), Sneutrino g2, D g3), FBF (-1, Psibar, SLR, Psi), 
           G_LQ_EC_UC (false,M2,g1,g2,g3)) ]

   let yuk_lqino_ec_su1' g1 g2 g3 m =
     let cm = chg m in
       [ ((LQino (-g1), Sup (m,g3), L g2), FBF (1, Psibar, SLR, Psi), 
             G_LQ_EC_UC (true,cm,g1,g2,g3)) ]

   let yuk_lqino_ec_su1 g1 g2 g3 =
      ThoList.flatmap (yuk_lqino_ec_su1' g1 g2 g3) [M1;M2]

   let yuk_lqino_ec_su2' g1 g2 g3 m =
     let cm = chg m in
       [ ((L (-g2), Sup (m,-g3), LQino (g1)), FBF (1, Psibar, SLR, Psi), 
             G_LQ_EC_UC (false,cm,g1,g2,g3)) ]

   let yuk_lqino_ec_su2 g1 g2 g3 =
      ThoList.flatmap (yuk_lqino_ec_su2' g1 g2 g3) [M1;M2]

   let yuk_lqino_nc_sd1 g1 g2 g3  =
     [ ((LQino (-g1), Sdown (M1,g3), N g2), FBF (-1, Psibar, SLR, Psi), 
           G_LQ_EC_UC (true,M2,g1,g2,g3)) ]

   let yuk_lqino_nc_sd2 g1 g2 g3  =
     [ ((N (-g2), Sdown (M1,-g3), LQino (g1)), FBF (-1, Psibar, SLR, Psi), 
           G_LQ_EC_UC (false,M2,g1,g2,g3)) ]

   let yuk_lq_ec_uc' g1 g2 g3 m =
     [ ((L (-g2), LQ (m,g1), U (-g3)), BBB (1, Psibar, SLR, Psibar), 
           G_LQ_EC_UC (false,m,g1,g2,g3)) ]
   
   let yuk_lq_ec_uc g1 g2 g3 =
      ThoList.flatmap (yuk_lq_ec_uc' g1 g2 g3) [M1;M2] 

   let yuk_lq_ec_uc2' g1 g2 g3 m =
     [ ((L (g2), LQ (m,-g1), U (g3)), PBP (1, Psi, SLR, Psi), 
           G_LQ_EC_UC (true,m,g1,g2,g3)) ]
   
   let yuk_lq_ec_uc2 g1 g2 g3 =
      ThoList.flatmap (yuk_lq_ec_uc2' g1 g2 g3) [M1;M2] 

   let yuk_lq_nc_dc g1 g2 g3 =
     [ ((N (-g2), LQ (M2,g1), D (-g3)), BBB (-1, Psibar, SLR, Psibar), 
           G_LQ_EC_UC (false,M2,g1,g2,g3)) ]
   
   let yuk_lq_nc_dc2 g1 g2 g3 =
     [ ((N (g2), LQ (M2,-g1), D (g3)), PBP (-1, Psi, SLR, Psi), 
           G_LQ_EC_UC (true,M2,g1,g2,g3)) ]
   
(*** Daniel Wiesler: LQ - F-Term w/ vev ***)

   let lq_se_su' g1 g2 g3 m1 m2 m3 =
     [ ((LQ (m1,g1), Slepton (m2,-g2), Sup (m3,-g3)), Scalar_Scalar_Scalar 1, 
           G_LQ_SSU (m1,m2,m3,g1,g2,g3)) ]
   
   let lq_se_su g1 g2 g3 =
      List.flatten (Product.list3 (lq_se_su' g1 g2 g3) [M1;M2] [M1;M2] [M1;M2] )

   let lq_snu_sd' g1 g2 g3 m1 m2  =
     [ ((LQ (m1,g1), Sdown (m2,-g2), Sneutrino (-g3)), Scalar_Scalar_Scalar 1, 
           G_LQ_SSD (m1,m2,g1,g2,g3)) ]
   
   let lq_snu_sd g1 g2 g3 =
      List.flatten (Product.list2 (lq_snu_sd' g1 g2 g3) [M1;M2] [M1;M2] )

(*** Daniel Wiesler: LQ - Higgs ***)

   let lq_shiggs' g1 s g2 m1 m2 =
     [ ((LQ (m1,g1), SHiggs s, LQ (m2,-g2)), Scalar_Scalar_Scalar 1, G_LQ_S (m1,m2,g1,s,g2))]

   let lq_shiggs g1 s g2 =
      List.flatten ( Product.list2 (lq_shiggs' g1 s g2) [M1;M2] [M1;M2]) 

   let lq_phiggs' g1 p g2 m1 m2 =
     [ ((LQ (m1,g1), PHiggs p, LQ (m2,-g2)), Scalar_Scalar_Scalar 1, G_LQ_P (m1,m2,g1,p,g2))]

   let lq_phiggs g1 p g2 =
      List.flatten ( Product.list2 (lq_phiggs' g1 p g2) [M1;M2] [M1;M2]) 

   let yuk_lqino_shiggs g1 s g2 =
      [ ((LQino (-g1), SHiggs s, LQino g2), FBF (1, Psibar, SLR, Psi), 
            G_YUK_LQ_S (g1,s,g2)) ]

   let yuk_lqino_phiggs g1 p g2 =
      [ ((LQino (-g1), PHiggs p, LQino g2), FBF (1, Psibar, SLR, Psi), 
            G_YUK_LQ_P (g1,p,g2)) ]

(*** Daniel Wiesler: LQ - Neutralinos. ***)

   let lqino_lq_neu' n g1 g2 m =
      [ ((Neutralino n, LQ (m,-g1), LQino g2), FBF (1, Chibar, SLR, Psi), 
            G_LQ_NEU (m,g1,g2,n)) ]

   let lqino_lq_neu n g1 g2 =
      ThoList.flatmap (lqino_lq_neu' n g1 g2) [M1;M2]

   let lqino_lq_neu2' n g1 g2 m =
      [ ((LQino (-g2), LQ (m,g1), Neutralino n), FBF (1, Psibar, SLR, Chi), 
            G_LQ_NEU (m,g1,g2,n)) ]

   let lqino_lq_neu2 n g1 g2 =
      ThoList.flatmap (lqino_lq_neu2' n g1 g2) [M1;M2]

(*** Daniel Wiesler: LQ-LQino-Gluino ***)
   let lqino_lq_gg' g1 g2 m =
     [ ((Gluino, LQ (m,-g1), LQino g2), FBF (1, Chibar, SLR, Psi), 
           G_LQ_GG (m,g1,g2)) ]

   let lqino_lq_gg  g1 g2 =
      ThoList.flatmap (lqino_lq_gg' g1 g2) [M1;M2]

(*** Daniel Wiesler: LQ - Gauge ***)

   let col_lqino_currents g =
      [ ((LQino (-g), Gl, LQino g), FBF ((-1), Psibar, V, Psi), Gs)]

   let neutr_lqino_current g = 
      [ ((LQino (-g), Z, LQino g), FBF (1, Psibar, V, Psi), G_NLQC)]

   let col_lq_currents m g =
      [ ((Gl, LQ (m,-g), LQ (m,g)), Vector_Scalar_Scalar (-1), Gs)]

   let lq_neutr_Z g m1 m2 =
      [ ((Z, LQ (m1,-g), LQ (m2,g)), Vector_Scalar_Scalar (-1), G_ZLQ (g,m1,m2))]

   let em_lq_currents g m = 
      [ ((Ga, LQ (m,-g), LQ (m,g)), Vector_Scalar_Scalar 1, Q_down)]       
   
   let em_lqino_currents g = 
      [ ((LQino (-g), Ga, LQino g), FBF (1, Psibar, V, Psi), Q_down)]       

   let gluon2_lq2' g m = 
      [ ((LQ (m,g), LQ (m,-g), Gl, Gl), Scalar2_Vector2 2, G_GlGlLQLQ)]
   let gluon2_lq2 g = 
      ThoList.flatmap (gluon2_lq2' g) [M1;M2] 

   let lq_gauge4' g m = 
      [ ((Z, Z, LQ (m,g), LQ (m,-g)), Scalar2_Vector2 1, G_ZZLQLQ);
        ((Z, Ga, LQ (m,g), LQ (m,-g)), Scalar2_Vector2 1, G_ZPLQLQ);
        ((Ga, Ga, LQ (m,g), LQ (m,-g)), Scalar2_Vector2 1, G_PPLQLQ)]

   let lq_gauge4 g = 
      ThoList.flatmap (lq_gauge4' g) [M1;M2] 

   let lq_gg_gauge2' g m = 
      [ ((Z, Gl, LQ (m,g), LQ (m,-g)), Scalar2_Vector2 1, G_ZGlLQLQ);
        ((Ga, Gl, LQ (m,g), LQ (m,-g)), Scalar2_Vector2 1, G_PGlLQLQ)]

   let lq_gg_gauge2 g = 
      ThoList.flatmap (lq_gg_gauge2' g) [M1;M2] 


   
   let col_sfermion_currents g m = 
      [ ((Gl, Sup (m,-g), Sup (m,g)), Vector_Scalar_Scalar (-1), Gs);
        ((Gl, Sdown (m,-g), Sdown (m,g)), Vector_Scalar_Scalar (-1), Gs)]

(*** REVISED: Compatible with CD+. **FB*)
   let triple_gauge =
      [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);  
        ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
        ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_G_S)]

(*** REVISED: Independent of the sign of CD. **FB*) 
   let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
   let minus_gauge4 = Vector4 [(-2, C_13_42); (1, C_12_34); (1, C_14_23)] 
   let quartic_gauge =
      [ (Wm, Wp, Wm, Wp), gauge4, G_WWWW;
        (Wm, Z, Wp, Z), minus_gauge4, G_ZZWW;
        (Wm, Z, Wp, Ga), minus_gauge4, G_PZWW;
        (Wm, Ga, Wp, Ga), minus_gauge4, G_PPWW;
        (Gl, Gl, Gl, Gl), gauge4, G_SS]

(* The [Scalar_Vector_Vector] couplings do not depend on the choice of the
   sign of the covariant derivative since they are quadratic in the
   gauge couplings. *)

(* JR: Only the first charged Higgs. *)
(*** REVISED: Compatible with CD+. ***)
(*** Revision: 2005-03-10: first two vertices corrected. ***)
(*** REVISED: Felix Braam: Compact version using new COMBOS*)
(*** REVISED: Felix Braam: Couplings adjusted to FF-convention*)
     let gauge_higgs_WPC p=
      [ ((Wm, CHiggs HC1, PHiggs p), Vector_Scalar_Scalar 1, G_GH_WPC p);
        ((Wp, CHiggs HC1c, PHiggs p), Vector_Scalar_Scalar 1, G_GH_WPC p)]
     let gauge_higgs_WSC s=
       [((Wm, CHiggs HC1, SHiggs s),Vector_Scalar_Scalar 1, G_GH_WSC s);
        ((Wp, CHiggs HC1c, SHiggs s),Vector_Scalar_Scalar (-1), G_GH_WSC s)]
     let gauge_higgs_ZSP s p =
        [((Z, SHiggs s, PHiggs p),Vector_Scalar_Scalar 1, G_GH_ZSP (s,p))]
     let gauge_higgs_WWS s=
        ((SHiggs s, Wp, Wm),Scalar_Vector_Vector 1, G_GH_WWS s)
     let gauge_higgs_ZZS s=
        ((SHiggs s, Z, Z), Scalar_Vector_Vector 1, G_GH_ZZS s)
     let gauge_higgs_ZCC =
        ((Z, CHiggs HC1, CHiggs HC1c),Vector_Scalar_Scalar 1, G_GH_ZCC )
     let gauge_higgs_GaCC =
        ((Ga, CHiggs HC1, CHiggs HC1c),Vector_Scalar_Scalar 1, G_GH_GaCC )

     let gauge_higgs =
       ThoList.flatmap gauge_higgs_WPC [P1;P2] @
       ThoList.flatmap gauge_higgs_WSC [S1;S2;S3] @
       List.flatten (Product.list2 gauge_higgs_ZSP [S1;S2;S3] [P1;P2]) @
       List.map gauge_higgs_WWS [S1;S2;S3] @
       List.map gauge_higgs_ZZS [S1;S2;S3] @
       [gauge_higgs_ZCC] @ [gauge_higgs_GaCC] 

(*** REVISED: Compact version using new COMBOS*)
(*** REVISED: Couplings adjusted to FF-convention*)

     let gauge_higgs4_ZZPP (p1,p2) = 
       ((PHiggs p1, PHiggs p2, Z, Z), Scalar2_Vector2 1, G_GH4_ZZPP (p1,p2))

     let gauge_higgs4_ZZSS (s1,s2) = 
        ((SHiggs s1, SHiggs s2 , Z, Z), Scalar2_Vector2 1, G_GH4_ZZSS (s1,s2))

(* JR: Only the first charged Higgs. *)
     let gauge_higgs4_ZZCC =
        ((CHiggs HC1, CHiggs HC1c, Z, Z), Scalar2_Vector2 1, G_GH4_ZZCC)

     let gauge_higgs4_GaGaCC =
        ((CHiggs HC1, CHiggs HC1c, Ga, Ga), Scalar2_Vector2 1, G_GH4_GaGaCC)

     let gauge_higgs4_ZGaCC =
        ((CHiggs HC1, CHiggs HC1c, Ga, Z), Scalar2_Vector2 1, G_GH4_ZGaCC )

     let gauge_higgs4_WWCC =
        ((CHiggs HC1, CHiggs HC1c, Wp, Wm), Scalar2_Vector2 1, G_GH4_WWCC )

     let gauge_higgs4_WWPP (p1,p2) =
        ((PHiggs p1, PHiggs p2, Wp, Wm), Scalar2_Vector2 1, G_GH4_WWPP (p1,p2))

     let gauge_higgs4_WWSS (s1,s2) =
        ((SHiggs s1, SHiggs s2, Wp, Wm), Scalar2_Vector2 1, G_GH4_WWSS (s1,s2))  

(* JR: Only the first charged Higgs. *)
     let gauge_higgs4_ZWSC s =
       [ ((CHiggs HC1, SHiggs s, Wm, Z), Scalar2_Vector2 1, G_GH4_ZWSC s); 
         ((CHiggs HC1c, SHiggs s, Wp, Z), Scalar2_Vector2 1, G_GH4_ZWSC s)]

     let gauge_higgs4_GaWSC s =
       [ ((CHiggs HC1, SHiggs s, Wm, Ga), Scalar2_Vector2 1, G_GH4_GaWSC s); 
         ((CHiggs HC1c, SHiggs s, Wp, Ga), Scalar2_Vector2 1, G_GH4_GaWSC s) ]

     let gauge_higgs4_ZWPC p =
       [ ((CHiggs HC1, PHiggs p, Wm, Z), Scalar2_Vector2 1, G_GH4_ZWPC p); 
         ((CHiggs HC1c, PHiggs p, Wp, Z), Scalar2_Vector2 (-1), G_GH4_ZWPC p)]

     let gauge_higgs4_GaWPC p =
       [ ((CHiggs HC1, PHiggs p, Wm, Ga), Scalar2_Vector2 1, G_GH4_GaWPC p); 
         ((CHiggs HC1c, PHiggs p, Wp, Ga), Scalar2_Vector2 (-1), 
             G_GH4_GaWPC p) ]
         
     let gauge_higgs4 = 
       List.map gauge_higgs4_ZZPP (pairs [P1;P2]) @
       List.map gauge_higgs4_ZZSS (pairs [S1;S2;S3]) @
       [gauge_higgs4_ZZCC] @ [gauge_higgs4_GaGaCC] @
       [gauge_higgs4_ZGaCC] @ [gauge_higgs4_WWCC] @
       List.map gauge_higgs4_WWPP (pairs [P1;P2]) @
       List.map gauge_higgs4_WWSS (pairs [S1;S2;S3]) @
       ThoList.flatmap gauge_higgs4_ZWSC [S1;S2;S3] @
       ThoList.flatmap gauge_higgs4_GaWSC [S1;S2;S3] @
       ThoList.flatmap gauge_higgs4_ZWPC [P1;P2] @
       ThoList.flatmap gauge_higgs4_GaWPC [P1;P2] 

(*** Added by Felix Braam. ***)
    let gauge_sfermion4' g m1 m2 =
       [ ((Wp, Wm, Slepton (m1,g), Slepton (m2,-g)), Scalar2_Vector2 1, 
            G_WWSFSF (SL,g,m1,m2));
        ((Z, Ga, Slepton (m1,g), Slepton (m2,-g)), Scalar2_Vector2 1, 
           G_ZPSFSF (SL,g,m1,m2));
        ((Z, Z, Slepton (m1,g), Slepton (m2,-g)), Scalar2_Vector2 1, 
           G_ZZSFSF (SL,g,m1,m2)); 
        ((Wp, Wm, Sup (m1,g), Sup (m2,-g)), Scalar2_Vector2 1, G_WWSFSF 
           (SU,g,m1,m2));
        ((Wp, Wm, Sdown (m1,g), Sdown (m2,-g)), Scalar2_Vector2 1, G_WWSFSF 
           (SD,g,m1,m2));
        ((Z, Z, Sup (m1,g), Sup (m2,-g)), Scalar2_Vector2 1, G_ZZSFSF 
           (SU,g,m1,m2));
        ((Z, Z, Sdown (m1,g), Sdown (m2,-g)), Scalar2_Vector2 1, G_ZZSFSF 
           (SD,g,m1,m2));
        ((Z, Ga, Sup (m1,g), Sup (m2,-g)), Scalar2_Vector2 1, G_ZPSFSF 
           (SU,g,m1,m2));
        ((Z, Ga, Sdown (m1,g), Sdown (m2,-g)), Scalar2_Vector2 1, G_ZPSFSF 
           (SD,g,m1,m2)) ]


    let gauge_sfermion4'' g m =
      [ ((Wp, Ga, Slepton (m,g), Sneutrino (-g)), Scalar2_Vector2 1, 
           G_WPSLSN (false,g,m));
        ((Wm, Ga, Slepton (m,-g), Sneutrino g), Scalar2_Vector2 1, 
           G_WPSLSN (true,g,m));
        ((Wp, Z, Slepton (m,g), Sneutrino (-g)), Scalar2_Vector2 1, 
           G_WZSLSN (false,g,m));
        ((Wm, Z, Slepton (m,-g), Sneutrino g), Scalar2_Vector2 1,
           G_WZSLSN (true,g,m));
        ((Ga, Ga, Slepton (m,g), Slepton (m,-g)), Scalar2_Vector2 1, G_PPSFSF SL); 
        ((Ga, Ga, Sup (m,g), Sup (m,-g)), Scalar2_Vector2 1, G_PPSFSF SU);
        ((Ga, Ga, Sdown (m,g), Sdown (m,-g)), Scalar2_Vector2 1, G_PPSFSF SD)]


    let gauge_sfermion4 g =
      List.flatten (Product.list2 (gauge_sfermion4' g) [M1;M2] [M1;M2]) @
      ThoList.flatmap (gauge_sfermion4'' g) [M1;M2] @
      [ ((Wp, Wm, Sneutrino g, Sneutrino (-g)), Scalar2_Vector2 1, G_WWSFSF 
           (SN,g,M1,M1));
        ((Z, Z, Sneutrino g, Sneutrino (-g)), Scalar2_Vector2 1, G_ZZSFSF 
           (SN,g,M1,M1)) ]


(*** Modified by Felix Braam. ***)
    let gauge_squark4'' g h m1 m2 = 
      [ ((Wp, Ga, Sup (m1,-g), Sdown (m2,h)), Scalar2_Vector2 1, G_WPSUSD 
           (false,m1,m2,g,h));
        ((Wm, Ga, Sup (m1,g), Sdown (m2,-h)), Scalar2_Vector2 1, G_WPSUSD 
           (true,m1,m2,g,h));
        ((Wp, Z, Sup (m1,-g), Sdown (m2,h)), Scalar2_Vector2 1, G_WZSUSD 
           (false,m1,m2,g,h));
        ((Wm, Z, Sup (m1,g), Sdown (m2,-h)), Scalar2_Vector2 1, G_WZSUSD 
           (true,m1,m2,g,h)) ]
    let gauge_squark4' g h = List.flatten (Product.list2 (gauge_squark4'' g h) 
                                              [M1;M2] [M1;M2])
    let gauge_squark4 =
      if Flags.ckm_present then
        List.flatten (Product.list2 gauge_squark4' [1;2;3] [1;2;3]) 
      else
        ThoList.flatmap (fun g -> gauge_squark4' g g) [1;2;3]

    let gluon_w_squark'' g h m1 m2 =
      [ ((Gl, Wp, Sup (m1,-g), Sdown (m2,h)), 
            Scalar2_Vector2 1, G_GlWSUSD (false,m1,m2,g,h));
        ((Gl, Wm, Sup (m1,g), Sdown (m2,-h)), 
            Scalar2_Vector2 1, G_GlWSUSD (true,m1,m2,g,h)) ]
    let gluon_w_squark' g h = 
      List.flatten (Product.list2 (gluon_w_squark'' g h) [M1;M2] [M1;M2])
    let gluon_w_squark = 
      if Flags.ckm_present then
        List.flatten (Product.list2 gluon_w_squark' [1;2;3] [1;2;3]) 
      else
        ThoList.flatmap (fun g -> gluon_w_squark' g g) [1;2;3]

(*** Modified by Felix Braam. ***)
    let gluon_gauge_squark' g m1 m2 =
      [ ((Gl, Z, Sup (m1,g), Sup (m2,-g)), 
            Scalar2_Vector2 2, G_GlZSFSF (SU,g,m1,m2));
        ((Gl, Z, Sdown (m1,g), Sdown (m2,-g)), 
            Scalar2_Vector2 2, G_GlZSFSF (SD,g,m1,m2)) ]
    let gluon_gauge_squark'' g m =
      [ ((Gl, Ga, Sup (m,g), Sup (m,-g)), Scalar2_Vector2 2, G_GlPSQSQ);
        ((Gl, Ga, Sdown (m,g), Sdown (m,-g)), Scalar2_Vector2 (-1), G_GlPSQSQ) ]

(*** Modified by Felix Braam. ***)
    let gluon_gauge_squark g =
      List.flatten (Product.list2 (gluon_gauge_squark' g) [M1;M2] [M1;M2]) @
      ThoList.flatmap (gluon_gauge_squark'' g) [M1;M2]

    let gluon2_squark2' g m = 
      [ ((Gl, Gl, Sup (m,g), Sup (m,-g)), Scalar2_Vector2 2, G_GlGlSQSQ);
        ((Gl, Gl, Sdown (m,g), Sdown (m,-g)), Scalar2_Vector2 2, G_GlGlSQSQ) ] 
    let gluon2_squark2 g = 
      ThoList.flatmap (gluon2_squark2' g) [M1;M2] 


(* JR: Only the first charged Higgs. *)
(*** REVISED: Independent of the sign of CD. ***)
(*** REVISED: Felix Braam: Compact version using new COMBOS *)
(*** REVISED: Felix Braam: Couplings adjusted to FF-convention *)
    let higgs_SCC s =
       ((CHiggs HC1, CHiggs HC1c, SHiggs s), Scalar_Scalar_Scalar 1, 
           G_H3_SCC s )
    let higgs_SSS (s1,s2,s3)=
        ((SHiggs s1, SHiggs s2, SHiggs s3), Scalar_Scalar_Scalar 1, 
            G_H3_SSS (s1,s2,s3))
    let higgs_SPP (p1,p2,s) =
        ((SHiggs s, PHiggs p1, PHiggs p2), Scalar_Scalar_Scalar 1, 
            G_H3_SPP (s,p1,p2))

    let higgs =
       List.map higgs_SCC [S1;S2;S3]@
       List.map higgs_SSS (triples [S1;S2;S3])@
       List.map higgs_SPP (two_and_one [P1;P2] [S1;S2;S3])


    let higgs4 = []
(* The vertices of the type Higgs - Sfermion - Sfermion are independent of 
   the choice of the CD sign since they are quadratic in the gauge 
   coupling. *) 

(* JR: Only the first charged Higgs. *)
(*** REVISED: Independent of the sign of CD. ***)
    let higgs_sneutrino' s g =
       ((SHiggs s, Sneutrino g, Sneutrino (-g)), Scalar_Scalar_Scalar 1, 
           G_SFSFS (s,SN,g,M1,M1))
      let higgs_sneutrino'' g m = 
        [((CHiggs HC1, Sneutrino (-g), Slepton (m,g)), 
              Scalar_Scalar_Scalar 1, G_HSNSL (false,g,m)); 
         ((CHiggs HC1c, Sneutrino g, Slepton (m,-g)), Scalar_Scalar_Scalar 1, 
              G_HSNSL (true,g,m))] 
      let higgs_sneutrino = 
        Product.list2 higgs_sneutrino' [S1;S2;S3] [1;2;3] @
        List.flatten ( Product.list2  higgs_sneutrino'' [1;2;3] [M1;M2] )   
        

(* Under the assumption that there is no mixing between the left- and
   right-handed sfermions for the first two generations there is only a 
   coupling of the form Higgs - sfermion1 - sfermion2 for the third 
   generation. All the others are suppressed by $m_f/M_W$. *)

(*** REVISED: Independent of the sign of CD. ***)
      let higgs_sfermion_S s g m1 m2 =
        [ ((SHiggs s, Slepton (m1,g), Slepton (m2,-g)), Scalar_Scalar_Scalar 1,
              G_SFSFS (s,SL,g,m1,m2));
          ((SHiggs s, Sup (m1,g), Sup (m2,-g)), Scalar_Scalar_Scalar 1, 
              G_SFSFS (s,SU,g,m1,m2));
          ((SHiggs s, Sdown (m1,g), Sdown (m2,-g)), Scalar_Scalar_Scalar 1, 
              G_SFSFS (s,SD,g,m1,m2))]

    let higgs_sfermion' g m1 m2 =
         (higgs_sfermion_S S1 g m1 m2) @ (higgs_sfermion_S S2 g m1 m2) @ (higgs_sfermion_S S3 g m1 m2)  
 
    let higgs_sfermion_P p g m1 m2 = 
        [ ((PHiggs p, Slepton (m1,g), Slepton (m2,-g)), Scalar_Scalar_Scalar 1,
              G_SFSFP (p,SL,g,m1,m2));
          ((PHiggs p, Sup (m1,g), Sup (m2,-g)), Scalar_Scalar_Scalar 1, 
              G_SFSFP (p,SU,g,m1,m2));
          ((PHiggs p, Sdown (m1,g), Sdown (m2,-g)), Scalar_Scalar_Scalar 1, 
              G_SFSFP (p,SD,g,m1,m2)) ]

    let higgs_sfermion'' g m1 m2 =
         (higgs_sfermion_P P1 g m1 m2) @ (higgs_sfermion_P P2 g m1 m2)   
    let higgs_sfermion = List.flatten (Product.list3 higgs_sfermion' [1;2;3] [M1;M2] [M1;M2])  @ 
        List.flatten (Product.list3 higgs_sfermion'' [1;2;3] [M1;M2] [M1;M2]) 

(* JR: Only the first charged Higgs. *)
(*** REVISED: Independent of the sign of CD. ***)
    let higgs_squark' g h m1 m2 =
      [ ((CHiggs HC1, Sup (m1,-g), Sdown (m2,h)), Scalar_Scalar_Scalar 1, 
              G_HSUSD (false,m1,m2,g,h)); 
        ((CHiggs HC1c, Sup (m1,g), Sdown (m2,-h)), Scalar_Scalar_Scalar 1, 
              G_HSUSD (true,m1,m2,g,h)) ]
    let higgs_squark_a g h = higgs_squark' g h M1 M1 
    let higgs_squark_b (g,h) = List.flatten (Product.list2 (higgs_squark' g h)
                                             [M1;M2] [M1;M2]) 
    let higgs_squark =          
      if Flags.ckm_present then
        List.flatten (Product.list2 higgs_squark_a [1;2] [1;2]) @ 
        ThoList.flatmap higgs_squark_b [(1,3);(2,3);(3,3);(3,1);(3,2)] 
      else
        higgs_squark_a 1 1 @ higgs_squark_a 2 2 @ higgs_squark_b (3,3)

    let vertices3 = 
        (ThoList.flatmap electromagnetic_currents_3 [1;2;3] @
         ThoList.flatmap electromagnetic_currents_2 [C1;C2] @
         List.flatten (Product.list2 electromagnetic_sfermion_currents [1;2;3] 
                         [M1;M2]) @ 
         ThoList.flatmap neutral_currents [1;2;3] @
         ThoList.flatmap neutral_sfermion_currents [1;2;3] @  
         ThoList.flatmap charged_currents [1;2;3] @
         List.flatten (Product.list2 charged_slepton_currents [1;2;3] 
                         [M1;M2]) @ 
         (if Flags.ckm_present then 
           List.flatten (Product.list2 charged_quark_currents [1;2;3] 
                           [1;2;3]) @
           List.flatten (Product.list2 charged_squark_currents [1;2;3] 
                           [1;2;3]) @ 
           ThoList.flatmap yukawa_higgs_quark [(1,3);(2,3);(3,3);(3,1);(3,2)]
         else
           charged_quark_currents 1 1 @
           charged_quark_currents 2 2 @
           charged_quark_currents 3 3 @
           charged_squark_currents 1 1 @
           charged_squark_currents 2 2 @
           charged_squark_currents 3 3 @ 
           ThoList.flatmap yukawa_higgs_quark [(3,3)]) @ 
(*i         ThoList.flatmap yukawa_higgs [1;2;3] @  i*)
         yukawa_higgs 3 @ yukawa_n @ 
         ThoList.flatmap yukawa_c [C1;C2] @ 
         ThoList.flatmap yukawa_cq [C1;C2] @ 
         List.flatten (Product.list2 charged_chargino_currents [N1;N2;N3;N4;N5] 
                         [C1;C2]) @ triple_gauge @ 
         ThoList.flatmap neutral_Z (pairs [N1;N2;N3;N4;N5]) @         
         Product.list2 charged_Z [C1;C2] [C1;C2] @ 
         gauge_higgs @ higgs @ yukawa_higgs_2 @ 
(*i         List.flatten (Product.list2 yukawa_higgs_quark [1;2;3] [1;2;3]) @  i*)
         List.flatten (Product.list2 higgs_charg_neutr [N1;N2;N3;N4;N5] [C1;C2]) @ 
         higgs_neutr @ higgs_sneutrino @ higgs_sfermion @ 
         higgs_squark @ yukawa_v @
         ThoList.flatmap col_currents [1;2;3] @
         List.flatten (Product.list2 col_sfermion_currents [1;2;3] [M1;M2])) @
         List.flatten (Product.list2 col_lq_currents [M1;M2] [1;2;3]) @
         ThoList.flatmap col_lqino_currents [1;2;3] @
         ThoList.flatmap em_lqino_currents [1;2;3] @
         ThoList.flatmap neutr_lqino_current [1;2;3] @
         List.flatten (Product.list3 yuk_lqino_se_uc1 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lqino_se_uc2 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lqino_ec_su1 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lqino_ec_su2 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lqino_sn_dc1 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lqino_sn_dc2 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lqino_nc_sd1 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lqino_nc_sd2 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lq_ec_uc [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lq_ec_uc2 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lq_nc_dc [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 yuk_lq_nc_dc2 [1;2;3] [1;2;3] [1;2;3]) @ 
         List.flatten (Product.list3 lq_neutr_Z [1;2;3] [M1;M2] [M1;M2]) @
         List.flatten (Product.list2 em_lq_currents [1;2;3] [M1;M2]) @ 
         List.flatten (Product.list3 lq_shiggs [1;2;3] [S1;S2;S3;S4;S5;S6;S7;S8;S9] [1;2;3]) @
         List.flatten (Product.list3 lq_phiggs [1;2;3] [P1;P2;P3;P4;P5;P6;P7] [1;2;3]) @
         List.flatten (Product.list3 yuk_lqino_shiggs [1;2;3] [S1;S2;S3;S4;S5;S6;S7;S8;S9] [1;2;3]) @
         List.flatten (Product.list3 yuk_lqino_phiggs [1;2;3] [P1;P2;P3;P4;P5;P6;P7] [1;2;3]) @
         List.flatten (Product.list3 lqino_lq_neu nlist [1;2;3] [1;2;3]) @
         List.flatten (Product.list3 lqino_lq_neu2 nlist [1;2;3] [1;2;3]) @
         List.flatten (Product.list3 lq_se_su [1;2;3] [1;2;3] [1;2;3]) @
         List.flatten (Product.list3 lq_snu_sd [1;2;3] [1;2;3] [1;2;3]) @
         List.flatten (Product.list2 lqino_lq_gg [1;2;3] [1;2;3])

    let vertices4 =
       (quartic_gauge @ higgs4 @ gauge_higgs4 @ 
        ThoList.flatmap gauge_sfermion4 [1;2;3] @
        gauge_squark4 @ gluon_w_squark @
        ThoList.flatmap gluon2_squark2  [1;2;3] @
        ThoList.flatmap gluon_gauge_squark [1;2;3] @
        ThoList.flatmap gluon2_lq2  [1;2;3] @            
        ThoList.flatmap lq_gauge4 [1;2;3] @
        ThoList.flatmap lq_gg_gauge2 [1;2;3])
        
    let vertices () = (vertices3, vertices4, [])

    let table = F.of_vertices (vertices ())
    let fuse2 = F.fuse2 table
    let fuse3 = F.fuse3 table
    let fuse = F.fuse table                                    
    let max_degree () = 4


(* SLHA2-Nomenclature for neutral Higgses *)
    let flavor_of_string s = 
      match s with
          | "e-" -> L 1 | "e+" -> L (-1)
          | "mu-" -> L 2 | "mu+" -> L (-2)
          | "tau-" -> L 3 | "tau+" -> L (-3)
          | "nue" -> N 1 | "nuebar" -> N (-1)
          | "numu" -> N 2 | "numubar" -> N (-2)
          | "nutau" -> N 3 | "nutaubar" -> N (-3)
          | "se1-" -> Slepton (M1,1) | "se1+" -> Slepton (M1,-1)
          | "smu1-" -> Slepton (M1,2) | "smu1+" -> Slepton (M1,-2)
          | "stau1-" -> Slepton (M1,3) | "stau1+" -> Slepton (M1,-3)
          | "se2-" -> Slepton (M2,1) | "se2+" -> Slepton (M2,-1)
          | "smu2-" -> Slepton (M2,2) | "smu2+" -> Slepton (M2,-2)
          | "stau2-" -> Slepton (M2,3) | "stau2+" -> Slepton (M2,-3)
          | "snue" -> Sneutrino 1 | "snue*" -> Sneutrino (-1)
          | "snumu" -> Sneutrino 2 | "snumu*" -> Sneutrino (-2)
          | "snutau" -> Sneutrino 3 | "snutau*" -> Sneutrino (-3)
          | "u" -> U 1 | "ubar" -> U (-1)
          | "c" -> U 2 | "cbar" -> U (-2)
          | "t" -> U 3 | "tbar" -> U (-3)
          | "d" -> D 1 | "dbar" -> D (-1)
          | "s" -> D 2 | "sbar" -> D (-2)
          | "b" -> D 3 | "bbar" -> D (-3)
          | "A" -> Ga | "Z" | "Z0" -> Z
          | "W+" -> Wp | "W-" -> Wm
          | "gl" | "g" -> Gl 
          | "h01" -> SHiggs S1 | "h02" -> SHiggs S2 | "h03" -> SHiggs S3 
          | "A01" -> PHiggs P1 | "A02" -> PHiggs P2 
          | "h04" -> SHiggs S4 | "h05" -> SHiggs S5 | "h06" -> SHiggs S6 
          | "A03" -> PHiggs P3 | "A04" -> PHiggs P4 
          | "h07" -> SHiggs S7 | "h08" -> SHiggs S8 | "h09" -> SHiggs S9 
          | "A05" -> PHiggs P5 | "A06" -> PHiggs P6 | "A07" -> PHiggs P7 
(* JR: Only the first charged Higgs. *)
          | "H+" -> CHiggs HC1 | "H-" -> CHiggs HC1c
          | "su1" -> Sup (M1,1) | "su1c" -> Sup (M1,-1)
          | "sc1" -> Sup (M1,2) | "sc1c" -> Sup (M1,-2)
          | "st1" -> Sup (M1,3) | "st1c" -> Sup (M1,-3)
          | "su2" -> Sup (M2,1) | "su2c" -> Sup (M2,-1)
          | "sc2" -> Sup (M2,2) | "sc2c" -> Sup (M2,-2)
          | "st2" -> Sup (M2,3) | "st2c" -> Sup (M2,-3)
          | "sgl" | "sg" -> Gluino
          | "sd1" -> Sdown (M1,1) | "sd1c" -> Sdown (M1,-1)
          | "ss1" -> Sdown (M1,2) | "ss1c" -> Sdown (M1,-2)
          | "sb1" -> Sdown (M1,3) | "sb1c" -> Sdown (M1,-3)
          | "sd2" -> Sdown (M2,1) | "sd2c" -> Sdown (M2,-1)
          | "ss2" -> Sdown (M2,2) | "ss2c" -> Sdown (M2,-2)
          | "sb2" -> Sdown (M2,3) | "sb2c" -> Sdown (M2,-3)
          | "neu1" -> Neutralino N1 | "neu2" -> Neutralino N2
          | "neu3" -> Neutralino N3 | "neu4" -> Neutralino N4      
          | "neu5" -> Neutralino N5 | "neu6" -> Neutralino N6 
          | "neu7" -> Neutralino N7 | "neu8" -> Neutralino N8 
          | "neu9" -> Neutralino N9 | "neu10" -> Neutralino N10 
          | "neu11" -> Neutralino N11
          | "ch1+" -> Chargino C1 | "ch2+" -> Chargino C2
          | "ch1-" -> Chargino C1c | "ch2-" -> Chargino C2c
          | "ch3+" -> Chargino C3 | "ch4+" -> Chargino C4
          | "ch3-" -> Chargino C3c | "ch4-" -> Chargino C4c
          | "lq11" -> LQ (M1,1) | "lq11c" -> LQ (M1,-1)
          | "lq12" -> LQ (M2,1) | "lq12c" -> LQ (M2,-1)
          | "lq21" -> LQ (M1,2) | "lq21c" -> LQ (M1,-2)
          | "lq22" -> LQ (M2,2) | "lq22c" -> LQ (M2,-2)
          | "lq31" -> LQ (M1,3) | "lq31c" -> LQ (M1,-3)
          | "lq32" -> LQ (M2,3) | "lq32c" -> LQ (M2,-3)
          | "lqino1" -> LQino 1 | "lqino1b" -> LQino (-1)
          | "lqino2" -> LQino 2 | "lqino2b" -> LQino (-2)
          | "lqino3" -> LQino 3 | "lqino3b" -> LQino (-3)
          | s -> invalid_arg ("HUBABUBA: %s Modellib_PSSSM.ExtMSSM.flavor_of_string:" ^ s)
                
    let flavor_to_string = function
      | L 1 -> "e-" | L (-1) -> "e+"
      | L 2 -> "mu-" | L (-2) -> "mu+"
      | L 3 -> "tau-" | L (-3) -> "tau+"
      | N 1 -> "nue" | N (-1) -> "nuebar"
      | N 2 -> "numu" | N (-2) -> "numubar"
      | N 3 -> "nutau" | N (-3) -> "nutaubar"
      | U 1 -> "u" | U (-1) -> "ubar"
      | U 2 -> "c" | U (-2) -> "cbar"
      | U 3 -> "t" | U (-3) -> "tbar"
      | U _ -> invalid_arg
            "Modellib_PSSSM.ExtMSSM.flavor_to_string: invalid up type quark"
      | D 1 -> "d" | D (-1) -> "dbar"
      | D 2 -> "s" | D (-2) -> "sbar"
      | D 3 -> "b" | D (-3) -> "bbar"
      | D _ -> invalid_arg
            "Modellib_PSSSM.ExtMSSM.flavor_to_string: invalid down type quark"
      | Gl -> "gl" | Gluino -> "sgl"
      | Ga -> "A" | Z -> "Z" 
      | Wp -> "W+" | Wm -> "W-"
      | SHiggs S1 -> "h01" | SHiggs S2 -> "h02" | SHiggs S3 -> "h03" 
      | PHiggs P1 -> "A01" | PHiggs P2 -> "A02"
      | SHiggs S4 -> "h04" | SHiggs S5 -> "h05" | SHiggs S6 -> "h06" 
      | PHiggs P3 -> "A03" | PHiggs P4 -> "A04"
      | SHiggs S7 -> "h07" | SHiggs S8 -> "h08" | SHiggs S9 -> "h09" 
      | PHiggs P5 -> "A05" | PHiggs P6 -> "A06" | PHiggs P7 -> "A07"
(* JR: Only the first charged Higgs. *)
      | CHiggs HC1 -> "H+" | CHiggs HC1c -> "H-"
      | CHiggs HC2 -> "HX_1+" | CHiggs HC2c -> "HX_1-"
      | CHiggs HC3 -> "HX_2+" | CHiggs HC3c -> "HX_2-"
      | CHiggs HC4 -> "HX_3+" | CHiggs HC4c -> "HX_3-"
      | CHiggs HC5 -> "HX_4+" | CHiggs HC5c -> "HX_4-"
      | Slepton (M1,1) -> "se1-" | Slepton (M1,-1) -> "se1+"
      | Slepton (M1,2) -> "smu1-" | Slepton (M1,-2) -> "smu1+"
      | Slepton (M1,3) -> "stau1-" | Slepton (M1,-3) -> "stau1+"
      | Slepton (M2,1) -> "se2-" | Slepton (M2,-1) -> "se2+"
      | Slepton (M2,2) -> "smu2-" | Slepton (M2,-2) -> "smu2+"
      | Slepton (M2,3) -> "stau2-" | Slepton (M2,-3) -> "stau2+"
      | Sneutrino 1 -> "snue" | Sneutrino (-1) -> "snue*"
      | Sneutrino 2 -> "snumu" | Sneutrino (-2) -> "snumu*"
      | Sneutrino 3 -> "snutau" | Sneutrino (-3) -> "snutau*"
      | Sup (M1,1) -> "su1" | Sup (M1,-1) -> "su1c"
      | Sup (M1,2) -> "sc1" | Sup (M1,-2) -> "sc1c"
      | Sup (M1,3) -> "st1" | Sup (M1,-3) -> "st1c"
      | Sup (M2,1) -> "su2" | Sup (M2,-1) -> "su2c"
      | Sup (M2,2) -> "sc2" | Sup (M2,-2) -> "sc2c"
      | Sup (M2,3) -> "st2" | Sup (M2,-3) -> "st2c"
      | Sdown (M1,1) -> "sd1" | Sdown (M1,-1) -> "sd1c"
      | Sdown (M1,2) -> "ss1" | Sdown (M1,-2) -> "ss1c"
      | Sdown (M1,3) -> "sb1" | Sdown (M1,-3) -> "sb1c"
      | Sdown (M2,1) -> "sd2" | Sdown (M2,-1) -> "sd2c"
      | Sdown (M2,2) -> "ss2" | Sdown (M2,-2) -> "ss2c"
      | Sdown (M2,3) -> "sb2" | Sdown (M2,-3) -> "sb2c"
      | Neutralino n -> "neu" ^ string_of_neu n
      | Chargino C1 -> "ch1+" | Chargino C1c -> "ch1-"
      | Chargino C2 -> "ch2+" | Chargino C2c -> "ch2-"
      | Chargino C3 -> "ch3+" | Chargino C3c -> "ch3-"
      | Chargino C4 -> "ch4+" | Chargino C4c -> "ch4-"
      | LQ (M1,1) -> "lq11" | LQ (M1,-1) -> "lq11c"
      | LQ (M2,1) -> "lq12" | LQ (M2,-1) -> "lq12c"
      | LQ (M1,2) -> "lq21" | LQ (M1,-2) -> "lq21c"
      | LQ (M2,2) -> "lq22" | LQ (M2,-2) -> "lq22c"
      | LQ (M1,3) -> "lq31" | LQ (M1,-3) -> "lq31c"
      | LQ (M2,3) -> "lq32" | LQ (M2,-3) -> "lq32c"
      | LQino 1 -> "lqino1" | LQino (-1) -> "lqino1b"
      | LQino 2 -> "lqino2" | LQino (-2) -> "lqino2b"
      | LQino 3 -> "lqino3" | LQino (-3) -> "lqino3b"
      | _ -> invalid_arg "Modellib_PSSSM.ExtMSSM.flavor_to_string"
                
    let flavor_to_TeX = function
      | L 1 -> "e^-" | L (-1) -> "e^+"
      | L 2 -> "\\mu^-" | L (-2) -> "\\mu^+"
      | L 3 -> "\\tau^-" | L (-3) -> "\\tau^+"
      | N 1 -> "\\nu_e" | N (-1) -> "\\bar{\\nu}_e"
      | N 2 -> "\\nu_\\mu" | N (-2) -> "\\bar{\\nu}_\\mu"
      | N 3 -> "\\nu_\\tau" | N (-3) -> "\\bar{\\nu}_\\tau"
      | U 1 -> "u" | U (-1) -> "\\bar{u}"
      | U 2 -> "c" | U (-2) -> "\\bar{c}"
      | U 3 -> "t" | U (-3) -> "\\bar{t}"
      | D 1 -> "d" | D (-1) -> "\\bar{d}"
      | D 2 -> "s" | D (-2) -> "\\bar{s}"
      | D 3 -> "b" | D (-3) -> "\\bar{b}"
      | L _ -> invalid_arg
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid lepton"
      | N _ -> invalid_arg
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid neutrino"
      | U _ -> invalid_arg
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid up type quark"
      | D _ -> invalid_arg
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid down type quark"
      | Gl -> "g" | Gluino -> "\\widetilde{g}"
      | Ga -> "\\gamma" | Z -> "Z" | Wp -> "W^+" | Wm -> "W^-"
      | SHiggs S1 -> "S_1" | SHiggs S2 -> "S_2" | SHiggs S3 -> "S_3" 
      | SHiggs S4 -> "S_4" | SHiggs S5 -> "S_5" | SHiggs S6 -> "S_6" 
      | SHiggs S7 -> "S_7" | SHiggs S8 -> "S_8" | SHiggs S9 -> "S_9" 
      | PHiggs P1 -> "P_1" | PHiggs P2 -> "P_2" | PHiggs P3 -> "P_3" 
      | PHiggs P4 -> "P_4" | PHiggs P5 -> "P_5" | PHiggs P6 -> "P_6" 
      | PHiggs P7 -> "P_7"
      | CHiggs HC1 -> "H^+" | CHiggs HC1c -> "H^-"
      | CHiggs HC2 -> "X_{H,1}^+" | CHiggs HC2c -> "X_{H,1}^-"
      | CHiggs HC3 -> "X_{H,2}^+" | CHiggs HC3c -> "X_{H,2}^-"
      | CHiggs HC4 -> "X_{H,3}^+" | CHiggs HC4c -> "X_{H,3}^-"
      | CHiggs HC5 -> "X_{H,4}^+" | CHiggs HC5c -> "X_{H,4}^-"
      | Slepton (M1,1) -> "\\widetilde{e}_1^-" 
      | Slepton (M1,-1) -> "\\widetilde{e}_1^+"
      | Slepton (M1,2) -> "\\widetilde{\\mu}_1^-" 
      | Slepton (M1,-2) -> "\\widetilde{\\mu}_1^+"
      | Slepton (M1,3) -> "\\widetilde{\\tau}_1^-" 
      | Slepton (M1,-3) -> "\\widetilde{\\tau}_1^+"
      | Slepton (M2,1) -> "\\widetilde{e}_2^-" 
      | Slepton (M2,-1) -> "\\widetilde{e}_2^+"
      | Slepton (M2,2) -> "\\widetilde{\\mu}_2^-" 
      | Slepton (M2,-2) -> "\\widetilde{\\mu}_2^+"
      | Slepton (M2,3) -> "\\widetilde{\\tau}_2^-" 
      | Slepton (M2,-3) -> "\\widetilde{\\tau}_2^+"
      | Sneutrino 1 -> "\\widetilde{\\nu}_e" 
      | Sneutrino (-1) -> "\\widetilde{\\nu}_e^*"
      | Sneutrino 2 -> "\\widetilde{\\nu}_\\mu" 
      | Sneutrino (-2) -> "\\widetilde{\\nu}_\\mu^*"
      | Sneutrino 3 -> "\\widetilde{\\nu}_\\tau" 
      | Sneutrino (-3) -> "\\widetilde{\\nu}_\\tau^*"
      | Sup (M1,1)  -> "\\widetilde{u}_1" 
      | Sup (M1,-1) -> "\\widetilde{u}_1^*"
      | Sup (M1,2)  -> "\\widetilde{c}_1" 
      | Sup (M1,-2) -> "\\widetilde{c}_1^*"
      | Sup (M1,3)  -> "\\widetilde{t}_1" 
      | Sup (M1,-3) -> "\\widetilde{t}_1^*"
      | Sup (M2,1)  -> "\\widetilde{u}_2" 
      | Sup (M2,-1) -> "\\widetilde{u}_2^*"
      | Sup (M2,2)  -> "\\widetilde{c}_2" 
      | Sup (M2,-2) -> "\\widetilde{c}_2^*"
      | Sup (M2,3)  -> "\\widetilde{t}_2" 
      | Sup (M2,-3) -> "\\widetilde{t}_2^*"
      | Sdown (M1,1)  -> "\\widetilde{d}_1" 
      | Sdown (M1,-1) -> "\\widetilde{d}_1^*"
      | Sdown (M1,2)  -> "\\widetilde{s}_1" 
      | Sdown (M1,-2) -> "\\widetilde{s}_1^*"
      | Sdown (M1,3)  -> "\\widetilde{b}_1" 
      | Sdown (M1,-3) -> "\\widetilde{b}_1^*"
      | Sdown (M2,1)  -> "\\widetilde{d}_2" 
      | Sdown (M2,-1) -> "\\widetilde{d}_2^*"
      | Sdown (M2,2)  -> "\\widetilde{s}_2" 
      | Sdown (M2,-2) -> "\\widetilde{s}_2^*"
      | Sdown (M2,3)  -> "\\widetilde{b}_2" 
      | Sdown (M2,-3) -> "\\widetilde{b}_2^*"
      | Neutralino N1 -> "\\widetilde{\\chi}^0_1"
      | Neutralino N2 -> "\\widetilde{\\chi}^0_2"
      | Neutralino N3 -> "\\widetilde{\\chi}^0_3"
      | Neutralino N4 -> "\\widetilde{\\chi}^0_4"
      | Neutralino N5 -> "\\widetilde{\\chi}^0_5"
      | Neutralino N6 -> "\\widetilde{\\chi}^0_6"
      | Neutralino N7 -> "\\widetilde{\\chi}^0_7"
      | Neutralino N8 -> "\\widetilde{\\chi}^0_8"
      | Neutralino N9 -> "\\widetilde{\\chi}^0_9"
      | Neutralino N10 -> "\\widetilde{\\chi}^0_{10}"
      | Neutralino N11 -> "\\widetilde{\\chi}^0_{11}"
      | Slepton _ -> invalid_arg
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid slepton"
      | Sneutrino _ -> invalid_arg
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid sneutrino"
      | Sup _ -> invalid_arg
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid up type squark"
      | Sdown _ -> invalid_arg 
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid down type squark"
      | Chargino C1  -> "\\widetilde{\\chi}_1^+" 
      | Chargino C1c -> "\\widetilde{\\chi}_1^-"
      | Chargino C2  -> "\\widetilde{\\chi}_2^+" 
      | Chargino C2c -> "\\widetilde{\\chi}_2^-"
      | Chargino C3  -> "\\widetilde{\\chi}_3^+" 
      | Chargino C3c -> "\\widetilde{\\chi}_3^-"
      | Chargino C4  -> "\\widetilde{\\chi}_4^+" 
      | Chargino C4c -> "\\widetilde{\\chi}_4^-"
      | LQ (M1,1) -> "D_{1,,1}" | LQ (M1,-1) -> "D_{1,,1}^*"
      | LQ (M2,1) -> "D_{1,,2}" | LQ (M2,-1) -> "D_{1,,2}^*"
      | LQ (M1,2) -> "D_{2,,1}" | LQ (M1,-2) -> "D_{2,,1}^*"
      | LQ (M2,2) -> "D_{2,,2}" | LQ (M2,-2) -> "D_{2,,2}^*"
      | LQ (M1,3) -> "D_{3,,1}" | LQ (M1,-3) -> "D_{3,,1}^*"
      | LQ (M2,3) -> "D_{3,,2}" | LQ (M2,-3) -> "D_{3,,2}^*"
      | LQino 1 -> "\\widetilde{D}_1" | LQino (-1) -> "\\bar\\widetilde{D}_1"
      | LQino 2 -> "\\widetilde{D}_2" | LQino (-2) -> "\\bar\\widetilde{D}_2"
      | LQino 3 -> "\\widetilde{D}_3" | LQino (-3) -> "\\bar\\widetilde{D}_3"
      | LQ _ -> invalid_arg 
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid leptoquark type"
      | LQino _ -> invalid_arg 
            "Modellib_PSSSM.ExtMSSM.flavor_to_TeX: invalid leptoquarkino type"

    let flavor_symbol = function
      | L g when g > 0 -> "l" ^ string_of_int g
      | L g -> "l" ^ string_of_int (abs g) ^ "b"  
      | N g when g > 0 -> "n" ^ string_of_int g
      | N g -> "n" ^ string_of_int (abs g) ^ "b"      
      | U g when g > 0 -> "u" ^ string_of_int g 
      | U g -> "u" ^ string_of_int (abs g) ^ "b"  
      | D g when g > 0 ->  "d" ^ string_of_int g 
      | D g -> "d" ^ string_of_int (abs g) ^ "b"    
      | Gl -> "gl" 
      | Ga -> "a" | Z -> "z"
      | Wp -> "wp" | Wm -> "wm"
      | Slepton (M1,g) when g > 0 -> "sl1" ^ string_of_int g 
      | Slepton (M1,g) -> "sl1c" ^ string_of_int (abs g)
      | Slepton (M2,g) when g > 0 -> "sl2" ^ string_of_int g
      | Slepton (M2,g) -> "sl2c" ^ string_of_int (abs g)
      | Sneutrino g when g > 0 -> "sn" ^ string_of_int g
      | Sneutrino g -> "snc" ^ string_of_int (abs g)
      | Sup (M1,g) when g > 0 -> "su1" ^ string_of_int g
      | Sup (M1,g) -> "su1c" ^ string_of_int (abs g)
      | Sup (M2,g) when g > 0 -> "su2" ^ string_of_int g
      | Sup (M2,g) -> "su2c" ^ string_of_int (abs g)
      | Sdown (M1,g) when g > 0 ->  "sd1" ^ string_of_int g
      | Sdown (M1,g) -> "sd1c" ^ string_of_int (abs g)
      | Sdown (M2,g) when g > 0 ->  "sd2" ^ string_of_int g
      | Sdown (M2,g) -> "sd2c" ^ string_of_int (abs g)
      | Neutralino n -> "neu" ^ (string_of_neu n)
      | Chargino c when (int_of_char c) > 0 -> "cp" ^ string_of_char c
      | Chargino c -> "cm" ^ string_of_int (abs (int_of_char c))
      | Gluino -> "sgl" 
      | SHiggs s -> "h0" ^ (string_of_shiggs s)
      | PHiggs p -> "A0" ^ (string_of_phiggs p)
      | CHiggs HC1 -> "hp" | CHiggs HC1c -> "hm" 
      | CHiggs _ -> invalid_arg "charged Higgs not yet implemented"
      | LQ (M1,g) when g > 0 -> "lq" ^ string_of_int g ^ "1"
      | LQ (M1,g) -> "lq" ^ string_of_int (abs g) ^ "1c"
      | LQ (M2,g) when g > 0 -> "lq" ^ string_of_int g ^ "2"
      | LQ (M2,g) -> "lq" ^ string_of_int (abs g) ^ "2c"
      | LQino g when g > 0 -> "lqino" ^ string_of_int g
      | LQino g -> "lqino" ^ string_of_int (abs g) ^ "b"

     let pdg = function
      | L g when g > 0 -> 9 + 2*g
      | L g -> - 9 + 2*g
      | N g when g > 0 -> 10 + 2*g
      | N g -> - 10 + 2*g
      | U g  when g > 0 -> 2*g
      | U g  -> 2*g
      | D g  when g > 0 -> - 1 + 2*g
      | D g  -> 1 + 2*g
      | Gl -> 21 
      | Ga -> 22 | Z -> 23
      | Wp -> 24 | Wm -> (-24)
      | SHiggs S1 -> 25 | SHiggs S2 -> 35 | PHiggs P1 -> 36
(* JR: Only the first charged Higgs. *)
      | CHiggs HC1 -> 37 | CHiggs HC1c -> (-37)
      | CHiggs _ -> invalid_arg "charged Higgs not yet implemented"
      | Slepton (M1,g) when g > 0 -> 1000009 + 2*g
      | Slepton (M1,g) -> - 1000009 + 2*g
      | Slepton (M2,g) when g > 0 -> 2000009 + 2*g
      | Slepton (M2,g) -> - 2000009 + 2*g            
      | Sneutrino g when g > 0 -> 1000010 + 2*g
      | Sneutrino g -> - 1000010 + 2*g            
      | Sup (M1,g) when g > 0 -> 1000000 + 2*g
      | Sup (M1,g) -> - 1000000 + 2*g
      | Sup (M2,g) when g > 0 -> 2000000 + 2*g
      | Sup (M2,g) -> - 2000000 + 2*g
      | Sdown (M1,g) when g > 0 -> 999999 + 2*g
      | Sdown (M1,g) -> - 999999 + 2*g
      | Sdown (M2,g) when g > 0 -> 1999999 + 2*g
      | Sdown (M2,g) -> - 1999999 + 2*g
      | Gluino -> 1000021
(* JR: only the first two charginos. *)
      | Chargino C1 -> 1000024 | Chargino C1c -> (-1000024)
      | Chargino C2 -> 1000037 | Chargino C2c -> (-1000037)
      | Chargino C3 -> 1000039 | Chargino C3c -> (-1000039)
      | Chargino C4 -> 1000041 | Chargino C4c -> (-1000041)
      | Neutralino N1 -> 1000022 | Neutralino N2 -> 1000023
      | Neutralino N3 -> 1000025 | Neutralino N4 -> 1000035
(* According to SLHA2 (not anymore ?!?)*)
      | Neutralino N5 -> 1000045 | Neutralino N6 -> 1000046 
      | Neutralino N7 -> 1000047 | Neutralino N8 -> 1000048 
      | Neutralino N9 -> 1000049 | Neutralino N10 -> 1000050 
      | Neutralino N11 -> 1000051
      | PHiggs P2 -> 46 | PHiggs P3 -> 47 | PHiggs P4 -> 48 
      | PHiggs P5 -> 49 | PHiggs P6 -> 50 | PHiggs P7 -> 51 
      | SHiggs S3 -> 45 | SHiggs S4 -> 52 | SHiggs S5 -> 53 
      | SHiggs S6 -> 54 | SHiggs S7 -> 55 | SHiggs S8 -> 56
      | SHiggs S9 -> 57
      | LQ (M1,g) when g > 0 -> 1000059 + g
      | LQ (M1,g) -> - 1000059 + g
      | LQ (M2,g) when g > 0 -> 2000059 + g
      | LQ (M2,g) -> - 2000059 + g
      | LQino g when g > 0 -> 59 + g
      | LQino g -> -59 + g

(* We must take care of the pdg numbers for the two different kinds of 
   sfermions in the MSSM. The particle data group in its Monte Carlo particle 
   numbering scheme takes only into account mixtures of the third generation 
   squarks and the stau. For the other sfermions we will use the number of the 
   lefthanded field for the lighter mixed state and the one for the righthanded
   for the heavier. Below are the official pdg numbers from the Particle
   Data Group. In order not to produce arrays with some million entries in 
   the Fortran code for the masses and the widths we introduce our private 
   pdg numbering scheme which only extends not too far beyond 42. 
   Our private scheme then has the following pdf numbers (for the sparticles
   the subscripts $L/R$ and $1/2$ are taken synonymously): 

   \begin{center}
      \renewcommand{\arraystretch}{1.2}
       \begin{tabular}{|r|l|l|}\hline
         $d$                    & down-quark         &      1 \\\hline
         $u$                    & up-quark           &      2 \\\hline
         $s$                    & strange-quark      &      3 \\\hline
         $c$                    & charm-quark        &      4 \\\hline
         $b$                    & bottom-quark       &      5 \\\hline
         $t$                    & top-quark          &      6 \\\hline\hline
         $e^-$                  & electron           &     11 \\\hline
         $\nu_e$                & electron-neutrino  &     12 \\\hline
         $\mu^-$                & muon               &     13 \\\hline
         $\nu_\mu$              & muon-neutrino      &     14 \\\hline
         $\tau^-$               & tau                &     15 \\\hline
         $\nu_\tau$             & tau-neutrino       &     16 \\\hline\hline
         $g$                    & gluon              & (9) 21 \\\hline
         $\gamma$               & photon             &     22 \\\hline
         $Z^0$                  & Z-boson            &     23 \\\hline
         $W^+$                  & W-boson            &     24 \\\hline\hline
         $h^0$                  & light Higgs boson  &     25 \\\hline
         $H^0$                  & heavy Higgs boson  &     35 \\\hline
         $A^0$                  & pseudoscalar Higgs &     36 \\\hline
         $H^+$                  & charged Higgs      &     37 \\\hline\hline
         $\tilde{d}_L$          & down-squark 1      &     41 \\\hline 
         $\tilde{u}_L$          & up-squark 1        &     42 \\\hline
         $\tilde{s}_L$          & strange-squark 1   &     43 \\\hline
         $\tilde{c}_L$          & charm-squark 1     &     44 \\\hline
         $\tilde{b}_L$          & bottom-squark 1    &     45 \\\hline
         $\tilde{t}_L$          & top-squark 1       &     46 \\\hline
         $\tilde{d}_R$          & down-squark 2      &     47 \\\hline 
         $\tilde{u}_R$          & up-squark 2        &     48 \\\hline
         $\tilde{s}_R$          & strange-squark 2   &     49 \\\hline
         $\tilde{c}_R$          & charm-squark 2     &     50 \\\hline
         $\tilde{b}_R$          & bottom-squark 2    &     51 \\\hline
         $\tilde{t}_R$          & top-squark 2       &     52 \\\hline\hline
         $\tilde{e}_L$          & selectron 1        &     53 \\\hline
         $\tilde{\nu}_{e,L}$    & electron-sneutrino &     54 \\\hline
         $\tilde{\mu}_L$        & smuon 1            &     55 \\\hline
         $\tilde{\nu}_{\mu,L}$  & muon-sneutrino     &     56 \\\hline
         $\tilde{\tau}_L$       & stau 1             &     57 \\\hline
         $\tilde{\nu}_{\tau,L}$ & tau-sneutrino      &     58 \\\hline
         $\tilde{e}_R$          & selectron 2        &     59 \\\hline
         $\tilde{\mu}_R$        & smuon 2            &     61 \\\hline
         $\tilde{\tau}_R$       & stau 2             &     63 \\\hline\hline
         $\tilde{g}$            & gluino             &     64 \\\hline
         $\tilde{\chi}^0_1$     & neutralino 1       &     65 \\\hline
         $\tilde{\chi}^0_2$     & neutralino 2       &     66 \\\hline
         $\tilde{\chi}^0_3$     & neutralino 3       &     67 \\\hline
         $\tilde{\chi}^0_4$     & neutralino 4       &     68 \\\hline
         $\tilde{\chi}^0_4$     & neutralino 5       &     69 \\\hline
         $\tilde{\chi4}^+_1$    & chargino 1         &     70 \\\hline
         $\tilde{\chi}^+_2$     & chargino 2         &     71 \\\hline\hline
         $a$                    & pseudoscalar       &     72 \\\hline
         $s$                    & scalar singlet     &     73 \\\hline
         $\tilde{G}$            & gravitino          &     -- \\\hline\hline 
     \end{tabular}
   \end{center}   *)

    let pdg_mw = function
      | L g when g > 0 -> 9 + 2*g
      | L g -> - 9 + 2*g
      | N g when g > 0 -> 10 + 2*g
      | N g -> - 10 + 2*g
      | U g when g > 0 -> 2*g
      | U g -> 2*g
      | D g when g > 0 -> - 1 + 2*g
      | D g -> 1 + 2*g
      | Gl -> 21
      | Ga -> 22 | Z -> 23
      | Wp -> 24 | Wm -> (-24)
      | SHiggs S1 -> 25 | SHiggs S2 -> 35 | PHiggs P1 -> 36
(* JR: Only the first charged Higgs. *)
      | CHiggs HC1 -> 37 | CHiggs HC1c -> (-37)
      | CHiggs _ -> invalid_arg "charged Higgs not yet implemented"
      | Sup (M1,g) when g > 0 -> 40 + 2*g
      | Sup (M1,g) -> - 40 + 2*g
      | Sup (M2,g) when g > 0 -> 46 + 2*g
      | Sup (M2,g) -> - 46 + 2*g
      | Sdown (M1,g) when g > 0 -> 39 + 2*g
      | Sdown (M1,g) -> - 39 + 2*g
      | Sdown (M2,g) when g > 0 -> 45 + 2*g
      | Sdown (M2,g) -> - 45 + 2*g           
      | Slepton (M1,g) when g > 0 -> 51 + 2*g
      | Slepton (M1,g) -> - 51 + 2*g
      | Slepton (M2,g) when g > 0 -> 57 + 2*g
      | Slepton (M2,g) -> - 57 + 2*g            
      | Sneutrino g when g > 0 ->  52 + 2*g
      | Sneutrino g -> - 52 + 2*g            
      | Gluino -> 64
(* JR: Only the first two charginos. *)
      | Chargino C1 -> 70 | Chargino C1c -> (-70)
      | Chargino C2 -> 71 | Chargino C2c -> (-71)
      | Chargino C3 -> 106 | Chargino C3c -> (-106)
      | Chargino C4 -> 107 | Chargino C4c -> (-107)
      | Neutralino N1 -> 65 | Neutralino N2 -> 66
      | Neutralino N3 -> 67 | Neutralino N4 -> 68 
      | Neutralino N5 -> 69 | Neutralino N6 -> 100 
      | Neutralino N7 -> 101 | Neutralino N8 -> 102
      | Neutralino N9 -> 103 | Neutralino N10 -> 104 
      | Neutralino N11 -> 105
      | PHiggs P2 -> 72 | PHiggs P3 -> 89 | PHiggs P4 -> 90 
      | PHiggs P5 -> 91 | PHiggs P6 -> 92 | PHiggs P7 -> 93   
      | SHiggs S3 -> 73 | SHiggs S4 -> 94 | SHiggs S5 -> 95 
      | SHiggs S6 -> 96 | SHiggs S7 -> 97 | SHiggs S8 -> 98 
      | SHiggs S9 -> 99             
      | LQ (M1,g) when g > 0 -> 78 + 2*g
      | LQ (M1,g) -> - 78 + 2*g
      | LQ (M2,g) when g > 0 -> 79 + 2*g
      | LQ (M2,g) -> - 79 + 2*g
      | LQino g when g > 0 -> 85 + g
      | LQino g -> - 85 + g

    let mass_symbol f =
      "mass(" ^ string_of_int (abs (pdg_mw f)) ^ ")"  

    let width_symbol f =
      "width(" ^ string_of_int (abs (pdg_mw f)) ^ ")"  

    let conj_symbol = function
      | false, str -> str
      | true, str -> str ^ "_c"

    let constant_symbol = function
      | E -> "e" | G -> "g" | G_Z -> "gz"
      | Q_lepton -> "qlep" | Q_up -> "qup" | Q_down -> "qdwn"
      | Q_charg -> "qchar"
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu" 
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_CC -> "gcc"
      | G_CCQ (vc,g1,g2) -> conj_symbol (vc, "g_ccq" ) ^ "(" 
          ^ string_of_int g1 ^ "," ^ string_of_int g2 ^ ")"
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_PZWW -> "gpzww" | G_PPWW -> "gppww"   
      | G_GH4_ZZPP (p1,p2) -> "g_ZZA0A0(" ^ string_of_phiggs p1 ^ "," ^ 
          string_of_phiggs p2 ^ ")" 
      | G_GH4_ZZSS (s1,s2) -> "g_ZZh0h0(" ^ string_of_shiggs s1 ^ "," ^ 
          string_of_shiggs s2 ^ ")"
      | G_GH4_ZZCC  -> "g_zzhphm"
      | G_GH4_GaGaCC -> "g_AAhphm"
      | G_GH4_ZGaCC -> "g_zAhphm"
      | G_GH4_WWCC -> "g_wwhphm"
      | G_GH4_WWPP (p1,p2) -> "g_WWA0A0(" ^ string_of_phiggs p1 ^ "," ^ 
          string_of_phiggs p2 ^ ")"
      | G_GH4_WWSS (s1,s2) -> "g_WWh0h0(" ^ string_of_shiggs s1 ^ "," ^ 
          string_of_shiggs s2 ^ ")"
      | G_GH4_ZWSC s -> "g_ZWhph0(" ^ string_of_shiggs s ^")"
      | G_GH4_GaWSC s -> "g_AWhph0(" ^ string_of_shiggs s ^")"
      | G_GH4_ZWPC p -> "g_ZWhpA0(" ^ string_of_phiggs p ^")"
      | G_GH4_GaWPC p -> "g_AWhpA0(" ^ string_of_phiggs p ^")"             
      | G_CICIS (n1,n2,s) -> "g_neuneuh0(" ^ string_of_neu n1 ^ "," ^ 
          string_of_neu n2 ^ "," ^ string_of_shiggs s ^ ")"
      | G_CICIP (n1,n2,p) ->  "g_neuneuA0(" ^ string_of_neu n1 ^ "," ^ 
          string_of_neu n2 ^ "," ^ string_of_phiggs p ^ ")" 
      | G_H3_SCC s -> "g_h0hphm(" ^ string_of_shiggs s ^ ")"
      | G_H3_SPP (s,p1,p2) -> "g_h0A0A0(" ^ string_of_shiggs s ^ "," ^ 
          string_of_phiggs p1 ^ "," ^ string_of_phiggs p2 ^ ")"
      | G_H3_SSS (s1,s2,s3) -> "g_h0h0h0(" ^ string_of_shiggs s1 ^ "," ^ 
          string_of_shiggs s2 ^ "," ^ string_of_shiggs s3 ^ ")"
      | G_CSC (c1,c2,s) -> "g_chchh0(" ^ string_of_char c1 ^ "," ^ 
          string_of_char c2 ^ "," ^ string_of_shiggs s ^ ")"  
      | G_CPC (c1,c2,p) ->  "g_chchA0(" ^ string_of_char c1 ^ "," ^ 
          string_of_char c2 ^ "," ^ string_of_phiggs p ^")"  
      | G_YUK_FFS (f1,f2,s) -> "g_yuk_h0_" ^ string_of_fermion_type f1 ^ 
          string_of_fermion_type f2 ^ "(" ^ string_of_shiggs s ^ "," ^ 
          string_of_fermion_gen f1 ^ ")"
      | G_YUK_FFP (f1,f2,p) -> "g_yuk_A0_" ^ string_of_fermion_type f1 ^ 
          string_of_fermion_type f2 ^ "(" ^ string_of_phiggs p ^ "," ^ 
          string_of_fermion_gen f1 ^ ")"
      | G_YUK_LCN g -> "g_yuk_hp_ln(" ^ string_of_int g ^ ")"
      | G_NWC (n,c) -> "g_nwc(" ^ string_of_char c ^ "," ^ string_of_neu n ^ ")" 
      | G_CWN (c,n) -> "g_cwn(" ^ string_of_char c ^ "," ^ string_of_neu n ^ ")" 
      | G_SLSNW (vc,g,m) -> conj_symbol (vc, "g_wslsn") ^ "(" ^ string_of_int g 
          ^ "," ^ string_of_sfm m ^ ")"
      | G_NZN (n1,n2) -> "g_zneuneu(" ^ string_of_neu n1 ^ "," 
          ^ string_of_neu n2 ^ ")"
      | G_CZC (c1,c2) -> "g_zchch(" ^ string_of_char c1 ^ "," 
          ^ string_of_char c2 ^ ")" 
      | Gs -> "gs"
      | G_YUK_UCD (n,m) -> "g_yuk_hp_ud(" ^ string_of_int n ^ "," ^ 
          string_of_int m ^ ")" 
      | G_YUK_DCU (n,m) -> "g_yuk_hm_du(" ^ string_of_int n ^ "," ^ 
          string_of_int m ^ ")" 
      | G_YUK_N (vc,f,n,sf,m) -> conj_symbol (vc, "g_yuk_neu_" ^ 
          string_of_fermion_type f ^ string_of_sff sf) ^ "(" ^ 
          string_of_fermion_gen f ^ "," ^ string_of_neu n ^ "," ^ 
          string_of_sfm m ^ ")" 
      | G_YUK_G (vc,f,sf,m) -> conj_symbol (vc, "g_yuk_gluino_" ^ 
          string_of_fermion_type f ^ string_of_sff sf) ^ "(" ^ 
          string_of_fermion_gen f  ^ "," ^ string_of_sfm m ^ ")"
      | G_YUK_C (vc,f,c,sf,m) -> conj_symbol (vc, "g_yuk_char_" ^ 
          string_of_fermion_type f ^ string_of_sff sf) ^ "(" ^ 
          string_of_fermion_gen f ^ "," ^ string_of_char c ^ "," ^ 
          string_of_sfm m ^ ")" 
      | G_YUK_Q (vc,g1,f,c,sf,m) -> conj_symbol (vc, "g_yuk_char_" ^ 
          string_of_fermion_type f ^ string_of_sff sf) ^"("^string_of_int g1 ^ 
          "," ^ string_of_fermion_gen f ^ "," ^ string_of_char c ^ "," ^ 
          string_of_sfm m ^ ")"
      | G_WPSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "g_wA_susd") ^ "(" ^ 
          string_of_int g1 ^ "," ^ string_of_int g2 ^ "," ^ string_of_sfm m1 ^ 
          "," ^ string_of_sfm m2 ^ ")" 
      | G_WZSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "g_wz_susd") ^ "(" ^ 
          string_of_int g1 ^ "," ^ string_of_int g2 ^ "," ^ string_of_sfm m1 
          ^ "," ^ string_of_sfm m2 ^ ")" 
      (* 3vertex: Higgs-Gauge a la Franke-Fraas *)

(* Nomenclature consistent with [flavor_of_string] *)     
      | G_GH_ZSP (s,p) -> "g_zh0a0(" ^ string_of_shiggs s ^ "," ^
          string_of_phiggs p ^ ")"
      | G_GH_WSC s -> "g_Whph0(" ^ string_of_shiggs s ^ ")"
      | G_GH_WPC p -> "g_WhpA0(" ^ string_of_phiggs p^ ")"        
      | G_GH_ZZS s -> "g_ZZh0(" ^ string_of_shiggs s ^ ")"  
      | G_GH_WWS s -> "g_WWh0(" ^ string_of_shiggs s ^ ")"
      | G_GH_ZCC -> "g_Zhmhp"
      | G_GH_GaCC -> "g_Ahmhp"
      | G_ZSF (f,g,m1,m2) -> "g_z" ^ string_of_sff f ^ string_of_sff f ^ "(" ^ 
          string_of_int g ^ "," ^ string_of_sfm m1 ^ "," ^ string_of_sfm m2 ^ ")" 
      | G_HSNSL (vc,g,m) -> conj_symbol (vc, "g_hp_sl" ^ string_of_sfm m ^ "sn1" )
          ^ "(" ^ string_of_int g ^ ")"
      | G_GlGlSQSQ -> "g_gg_sqsq" 
      | G_PPSFSF f -> "g_AA_" ^ string_of_sff f ^ string_of_sff f 
      | G_ZZSFSF (f,g,m1,m2) -> "g_zz_" ^ string_of_sff f ^string_of_sff f ^ "(" ^ 
          string_of_int g ^ "," ^ string_of_sfm m1 ^ "," ^ string_of_sfm m2 ^ ")" 
      | G_ZPSFSF (f,g,m1,m2) -> "g_zA_" ^ string_of_sff f ^string_of_sff f ^ "(" 
          ^ string_of_int g ^","^ string_of_sfm m1 ^ "," ^ string_of_sfm m2 ^ ")" 
      | G_GlPSQSQ -> "g_gA_sqsq" 
      | G_GlZSFSF (f,g,m1,m2) -> "g_gz_" ^ string_of_sff f ^ string_of_sff f ^ 
          "("  ^ string_of_int g ^ "," ^ string_of_sfm m1 ^ "," ^ string_of_sfm 
          m2 ^ ")"
      | G_GlWSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "g_gw_susd") ^ "(" ^ 
          string_of_int g1 ^ "," ^string_of_int g2 ^ "," ^ string_of_sfm m1 ^ "," 
          ^ string_of_sfm m2 ^ ")" 
      | G_strong -> "gs" | G_SS -> "gs**2" 
      | I_G_S -> "igs"           
      | G_NHC (vc,n,c) -> conj_symbol(vc,"g_neuhmchar") ^ "(" ^ 
          string_of_neu n ^ "," ^ string_of_char c ^ ")"
      | G_WWSFSF (f,g,m1,m2) -> "g_ww_" ^ string_of_sff f ^ string_of_sff f ^ 
          "(" ^ string_of_int g ^ "," ^ string_of_sfm m1 ^ "," ^ string_of_sfm m2 
          ^ ")"
      | G_WPSLSN (vc,g,m) -> conj_symbol (vc, "g_wA_slsn") ^"("^ string_of_int g 
          ^ "," ^ string_of_sfm m ^ ")" 
      | G_WZSLSN (vc,g,m) -> conj_symbol (vc, "g_wz_slsn") ^ "(" ^ string_of_int 
          g ^ "," ^ string_of_sfm m ^ ")" 
      | G_SFSFS (s,f,g,m1,m2) -> "g_h0_"^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "(" ^ string_of_shiggs s ^ "," ^
          string_of_int g ^ ")"   
      | G_SFSFP (p,f,g,m1,m2) -> "g_A0_"^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "(" ^ string_of_phiggs p ^ "," ^
          string_of_int g ^ ")"
      | G_HSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "g_hp_su" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ) ^ "(" ^ string_of_int g1 ^ "," 
          ^ string_of_int g2 ^ ")"
      | G_WSQ (vc,g1,g2,m1,m2) -> conj_symbol (vc, "g_wsusd") ^ "(" 
          ^ string_of_int g1 ^ "," ^ string_of_int g2 ^ "," ^ string_of_sfm m1 ^
          "," ^ string_of_sfm m2 ^ ")"
      | G_YUK_LQ_S (g1,s,g3) -> "g_yuk_lq_s(" ^ string_of_int g1 ^ "," ^ 
          string_of_shiggs s ^"," ^ string_of_int g3 ^")"       
      | G_YUK_LQ_P (g1,p,g3) -> "g_yuk_lq_p(" ^ string_of_int g1 ^ "," ^ 
          string_of_phiggs p ^ "," ^ string_of_int g3 ^ ")"       
      | G_LQ_NEU (m,g1,g2,n) -> "g_lq_neu(" ^ string_of_sfm m ^ ","  ^ 
          string_of_int g1 ^ "," ^ string_of_int g2 ^ "," ^ string_of_neu n ^ ")"
      | G_LQ_GG (m,g1,g2) -> "g_lq_gg(" ^ string_of_sfm m ^ ","  ^ 
          string_of_int g1 ^ "," ^ string_of_int g2 ^ ")"       
      | G_LQ_EC_UC (vc,m,g1,g2,g3) -> conj_symbol(vc,"g_lq_ec_uc") ^ "("  ^ 
          string_of_sfm m ^ "," ^ string_of_int g1 ^ "," ^ string_of_int g2 ^ "," 
          ^ string_of_int g3 ^ ")"       
      | G_LQ_SSU (m1,m2,m3,g1,g2,g3) -> "g_lq_sst(" ^ string_of_sfm m1 ^ ","  ^ 
          string_of_sfm m2 ^ ","  ^ string_of_sfm m3 ^ "," ^ string_of_int g1 ^ "," ^ 
          string_of_int g2 ^ "," ^ string_of_int g3 ^ ")"       
      | G_LQ_SSD (m1,m2,g1,g2,g3) -> "g_lq_ssta(" ^ string_of_sfm m1 ^ ","  ^ 
          string_of_sfm m2 ^ ","  ^ string_of_int g1 ^ "," ^ string_of_int g2 ^ 
          "," ^ string_of_int g3 ^ ")"       
      | G_LQ_S (m1,m2,g1,s,g2) -> "g_lq_s(" ^ string_of_sfm m1 ^ "," ^ string_of_sfm m2 ^ ","
          ^ string_of_int g1 ^ "," ^ string_of_shiggs s ^ "," ^ string_of_int g2 ^ ")"
      | G_LQ_P (m1,m2,g1,p,g2) -> "g_lq_s(" ^ string_of_sfm m1 ^ ","  ^ string_of_sfm m2 ^ ","
          ^ string_of_int g1 ^ "," ^ string_of_phiggs p ^ "," ^ string_of_int g2 ^ ")"
      | G_ZLQ (g,m1,m2) -> "g_zlqlq(" ^ string_of_int g ^ "," ^ 
          string_of_sfm m1 ^ "," ^ string_of_sfm m2 ^ ")"
      | G_ZZLQLQ -> "g_zz_lqlq"
      | G_ZPLQLQ -> "g_zA_lqlq"
      | G_PPLQLQ -> "g_AA_lqlq"
      | G_ZGlLQLQ -> "g_zg_lqlq"
      | G_PGlLQLQ -> "g_Ag_lqlq"
      | G_GlGlLQLQ -> "g_gg_lqlq"
      | G_NLQC -> "g_nlqc"
      
  end
