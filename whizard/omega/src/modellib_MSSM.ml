(* modellib_MSSM.ml --

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

(* modellib_MSSM.ml -- *)

(* \thocwmodulesection{Minimal Supersymmetric Standard Model} *)

module type MSSM_flags =
  sig
    val include_goldstone : bool
    val include_four      : bool
    val ckm_present       : bool
    val gravitino         : bool
    val higgs_triangle    : bool
  end

module MSSM_no_goldstone : MSSM_flags =
  struct
    let include_goldstone = false
    let include_four      = true
    let ckm_present       = false
    let gravitino         = false
    let higgs_triangle    = false
  end

module MSSM_goldstone : MSSM_flags =
  struct
    let include_goldstone = true
    let include_four      = true
    let ckm_present       = false
    let gravitino         = false
    let higgs_triangle    = false    
  end

module MSSM_no_4 : MSSM_flags = 
  struct 
    let include_goldstone = false
    let include_four      = false
    let ckm_present       = false
    let gravitino         = false
    let higgs_triangle    = false
  end

module MSSM_no_4_ckm : MSSM_flags = 
  struct 
    let include_goldstone = false
    let include_four      = false
    let ckm_present       = true
    let gravitino         = false
    let higgs_triangle    = false
  end

module MSSM_Grav : MSSM_flags = 
  struct 
    let include_goldstone = false
    let include_four      = false
    let ckm_present       = false
    let gravitino         = true
    let higgs_triangle    = false
  end

module MSSM_Hgg : MSSM_flags = 
  struct 
    let include_goldstone = false
    let include_four      = false
    let ckm_present       = false
    let gravitino         = false
    let higgs_triangle    = true
  end


module MSSM (Flags : MSSM_flags) = 
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
      | C1 | C2 | C1c | C2c

    type neu =
      | N1 | N2 | N3 | N4 

    let int_of_char = function
      | C1 -> 1 | C2 -> 2 | C1c -> -1 | C2c -> -2

    let string_of_char = function
      | C1 -> "1" | C2 -> "2" | C1c -> "-1" | C2c -> "-2"

    let conj_char = function
      | C1 -> C1c | C2 -> C2c | C1c -> C1 | C2c -> C2

    let string_of_neu = function
      | N1 -> "1" | N2 -> "2" | N3 -> "3" | N4 -> "4" 

(* Also we need types to distinguish the Higgs bosons. We follow the 
   conventions of Kuroda, which means
   \begin{align}
   \label{eq:higgs3}
   H_1 &= 
      \begin{pmatrix} 
          \frac{1}{\sqrt{2}} 
          \bigl( 
          v_1 + H^0 \cos\alpha - h^0
          \sin\alpha + \ii A^0 \sin\beta - \ii \phi^0 \cos\beta 
          \bigr) \\ 
          H^- \sin\beta - \phi^- \cos\beta 
       \end{pmatrix},
       \\ & \notag \\
   H_2 & = 
      \begin{pmatrix} 
          H^+ \cos\beta + \phi^+ \sin\beta \\
          \frac{1}{\sqrt{2}}
          \bigl( 
           v_2 + H^0 \sin\alpha + h^0 \cos\alpha + \ii A^0 \cos\beta + 
           \ii \phi^0 \sin\beta 
           \bigr)
      \end{pmatrix} 
      \label{eq:higgs4}
   \end{align}
   This is a different sign convention compared to, e.g., 
   Weinberg's volume iii. We will refer to it as [GS+]. 
*)

    type higgs =
      | H1   (* the light scalar Higgs *)
      | H2   (* the heavy scalar Higgs *)
      | H3   (* the pseudoscalar Higgs *)
      | H4   (* the charged Higgs *)
      | H5   (* the neutral Goldstone boson *)
      | H6   (* the charged Goldstone boson *)
      | DH of higgs*higgs

    let rec string_of_higgs = function
      | H1 -> "h1" | H2 -> "h2" | H3 -> "h3" | H4 -> "h4" 
      | H5 -> "p1" | H6 -> "p2" 
      | DH (h1,h2) -> string_of_higgs h1 ^ string_of_higgs h2

    type flavor =
      | L of int | N of int
      | U of int | D of int
      | Sup of sfm*int | Sdown of sfm*int 
      | Ga | Wp | Wm | Z | Gl
      | Slepton of sfm*int | Sneutrino of int 
      | Neutralino of neu | Chargino of char 
      | Gluino | Grino 
      | Phip | Phim | Phi0 | H_Heavy | H_Light | Hp | Hm | A 

    type gauge = unit

    let gauge_symbol () =
      failwith "Modellib_MSSM.MSSM.gauge_symbol: internal error"       

(* At this point we will forget graviton and -tino. *) 
        
    let lep_family g = [ L g; N g; Slepton (M1,g); 
                         Slepton (M2,g); Sneutrino g ] 
    let family g = 
      [ L g; N g; Slepton (M1,g); Slepton (M2,g); Sneutrino g;
        U g; D g; Sup (M1,g); Sup (M2,g); Sdown (M1,g); 
        Sdown (M2,g)]

    let external_flavors'' = 
        [ "1st Generation", ThoList.flatmap family [1; -1];
          "2nd Generation", ThoList.flatmap family [2; -2];
          "3rd Generation", ThoList.flatmap family [3; -3];
          "Gauge Bosons", [Ga; Z; Wp; Wm; Gl];
          "Charginos", [Chargino C1; Chargino C2; Chargino C1c; Chargino C2c];
          "Neutralinos", [Neutralino N1; Neutralino N2; Neutralino N3; 
                          Neutralino N4];
          "Higgs Bosons", [H_Heavy; H_Light; Hp; Hm; A];
          "Gluinos", [Gluino]]
    let external_flavors' =
      if Flags.gravitino then external_flavors'' @ ["Gravitino", [Grino]]
      else
        external_flavors''
    let external_flavors () =
      if Flags.include_goldstone then external_flavors' @ ["Goldstone Bosons", 
                                                           [Phip; Phim; Phi0]]
      else
        external_flavors' 

          
    let flavors () = ThoList.flatmap snd (external_flavors ())

    let spinor n =
      if n >= 0 then
        Spinor
      else if
        n <= 0 then
        ConjSpinor
      else
        invalid_arg "Modellib_MSSM.MSSM.spinor: internal error"

    let lorentz = function
      | L g -> spinor g | N g -> spinor g
      | U g -> spinor g | D g -> spinor g
      | Chargino c -> spinor (int_of_char c) 
      | Ga -> Vector  
(*i      | Ga -> Ward_Vector i*)
      | Gl -> Vector
      | Wp | Wm | Z -> Massive_Vector
      | H_Heavy | H_Light | Hp | Hm | A -> Scalar
      | Phip | Phim | Phi0 -> Scalar
      | Sup _ | Sdown _ | Slepton _ | Sneutrino _ -> Scalar 
      | Neutralino _ -> Majorana 
      | Gluino -> Majorana
      | Grino -> Vectorspinor

    let color = function
      | U g -> Color.SUN (if g > 0 then 3 else -3)
      | Sup (m,g) -> Color.SUN (if g > 0 then 3 else -3)
      | D g -> Color.SUN (if g > 0 then 3 else -3)
      | Sdown (m,g) -> Color.SUN  (if g > 0 then 3 else -3)
      | Gl | Gluino -> Color.AdjSUN 3
      | _ -> Color.Singlet   

    let nc () = 3

    let prop_spinor n =
      if n >= 0 then
        Prop_Spinor
      else if 
        n <=0 then
        Prop_ConjSpinor
      else 
        invalid_arg "Modellib_MSSM.MSSM.prop_spinor: internal error"

    let propagator = function
      | L g -> prop_spinor g | N g -> prop_spinor g
      | U g -> prop_spinor g | D g -> prop_spinor g
      | Chargino c -> prop_spinor (int_of_char c)
      | Ga | Gl -> Prop_Feynman
      | Wp | Wm | Z -> Prop_Unitarity
      | H_Heavy | H_Light | Hp | Hm | A -> Prop_Scalar
      | Phip | Phim | Phi0 -> if Flags.include_goldstone then Prop_Scalar 
                                     else Only_Insertion  
      | Slepton _ | Sneutrino _ | Sup _ | Sdown _ -> Prop_Scalar
      | Gluino -> Prop_Majorana | Neutralino _ -> Prop_Majorana
      | Grino -> Only_Insertion

(* Note, that we define the gravitino only as an insertion since when using propagators
   we are effectively going to a higher order in the gravitational coupling. This would
   enforce us to also include higher-dimensional vertices with two gravitinos for 
   a consistent power counting in $1/M_{\text{Planck}}$. *)

(*i      | Grino -> Prop_Vectorspinor   i*)

(* Optionally, ask for the fudge factor treatment for the widths of
   charged particles.  Currently, this only applies to $W^\pm$ and top. *)
                                                                     
    let width f =
      if !use_fudged_width then
        match f with
        | Wp | Wm | U 3 | U (-3) -> Fudged
        | _ -> !default_width
      else
        !default_width 

(* For the Goldstone bosons we adopt the conventions of the Kuroda paper. 
    \begin{subequations}
    \begin{equation}
        H_1 \equiv \begin{pmatrix} \left( v_1 + H^0 \cos\alpha - h^0 \sin
        \alpha + \ii A^0 \sin\beta - \ii \cos\beta \phi^0 \right) / \sqrt{2} \\
        H^- \sin\beta - \phi^- \cos\beta \end{pmatrix}
    \end{equation}
    \begin{equation}
        H_2 \equiv \begin{pmatrix} H^+ \cos\beta + \phi^+ \sin\beta \\ \left( 
        v_2 + H^0 \sin\alpha + h^0 \cos\alpha + \ii A^0 \cos\beta + \ii 
        \phi^0 \sin\beta \right) / \sqrt{2} \end{pmatrix}        
   \end{equation}
   \end{subequations}
*)

    let goldstone = function
      | Wp -> Some (Phip, Coupling.Integer 1)
      | Wm -> Some (Phim, Coupling.Integer 1)
      | Z -> Some (Phi0, Coupling.Integer 1)
      | _ -> None

    let conjugate = function
      | L g -> L (-g) | N g -> N (-g)
      | U g -> U (-g) | D g -> D (-g)
      | Sup (m,g) -> Sup (m,-g) 
      | Sdown (m,g) -> Sdown (m,-g) 
      | Slepton (m,g) -> Slepton (m,-g) 
      | Sneutrino g -> Sneutrino (-g)
      | Gl -> Gl (* | Gl0 -> Gl0 *)
      | Ga -> Ga | Z -> Z | Wp -> Wm | Wm -> Wp
      | H_Heavy -> H_Heavy | H_Light -> H_Light | A -> A
      | Hp -> Hm | Hm -> Hp 
      | Phip -> Phim | Phim -> Phip | Phi0 -> Phi0
      | Gluino -> Gluino 
      | Grino -> Grino
      | Neutralino n -> Neutralino n | Chargino c -> Chargino (conj_char c)

   let fermion = function
     | L g -> if g > 0 then 1 else -1
     | N g -> if g > 0 then 1 else -1
     | U g -> if g > 0 then 1 else -1
     | D g -> if g > 0 then 1 else -1
     | Gl | Ga | Z | Wp | Wm -> 0    (* | Gl0 -> 0 *)
     | H_Heavy | H_Light | Hp | Hm | A -> 0       
     | Phip | Phim | Phi0 -> 0            
     | Neutralino _ -> 2
     | Chargino c -> if (int_of_char c) > 0 then 1 else -1
     | Sup _ -> 0 | Sdown _ -> 0 
     | Slepton _ -> 0 | Sneutrino _ -> 0          
     | Gluino | Grino -> 2 

(* Because the O'Caml compiler only allows 248 constructors we must divide the 
   constants into subgroups of constants, e.g. for the Higgs couplings. In the 
   MSSM there are a lot of angles among the parameters, the Weinberg-angle, the 
   angle describing the Higgs vacuum structure, the mixing angle of the real 
   parts of the Higgs dubletts, the mixing angles of the sfermions. Therefore we 
   are going to define the trigonometric functions of those angles not as 
   constants but as functors of the angels. Sums and differences of angles are 
   only used as arguments for the $\alpha$ and $\beta$ angles, so it makes no 
   sense to define special functions for differences and sums of angles. *)

    type angle = 
      | Thw | Al | Be | Th_SF of sff*int | Delta | CKM_12 | CKM_13 | CKM_23

    let string_of_angle = function
      | Thw -> "thw" | Al -> "al" | Be -> "be" | Delta -> "d" 
      | CKM_12 -> "ckm12" | CKM_13 -> "ckm13" | CKM_23 -> "ckm23"
      | Th_SF (f,g) -> "th" ^ string_of_sff f ^ string_of_int g          

(* We introduce a Boolean type vc as a pseudonym for Vertex Conjugator to 
   distinguish between vertices containing complex mixing matrices like the 
   CKM--matrix or the sfermion or neutralino/chargino--mixing matrices, which 
   have to become complex conjugated. The true--option stands for the conjugated 
   vertex, the false--option for the unconjugated vertex. *)

    type vc = bool

    type constant =
      | Unit | Pi | Alpha_QED | Sin2thw
      | Sin of angle | Cos of angle | E | G | Vev | Tanb | Tana 
      | Cos2be | Cos2al | Sin2be | Sin2al | Sin4al | Sin4be | Cos4be
      | Cosapb | Cosamb | Sinapb | Sinamb | Cos2am2b | Sin2am2b       
      | Eidelta
      | Mu | AU of int | AD of int | AL of int 
      | V_CKM of int*int | M_SF of sff*int*sfm*sfm 
      | M_V of char*char  (* left chargino mixing matrix *)
      | M_U of char*char  (* right chargino mixing matrix *)
      | M_N of neu*neu   (* neutralino mixing matrix *)
      | V_0 of neu*neu | A_0 of neu*neu | V_P of char*char | A_P of char*char
      | L_CN of char*neu | R_CN of char*neu | L_NC of neu*char | R_NC of neu*char
(*i      | L_NF of neu*sff*sfm | R_NF of neu*sff*sfm i*)
      | S_NNH1 of neu*neu | P_NNH1 of neu*neu
      | S_NNH2 of neu*neu | P_NNH2 of neu*neu 
      | S_NNA of neu*neu | P_NNA of neu*neu
      | S_NNG of neu*neu | P_NNG of neu*neu
      | L_CNG of char*neu | R_CNG of char*neu
      | L_NCH of neu*char | R_NCH of neu*char
      | Q_lepton | Q_up | Q_down | Q_charg           
      | G_Z | G_CC | G_CCQ of vc*int*int
      | G_NC_neutrino | G_NC_lepton | G_NC_up | G_NC_down 
      | I_Q_W | I_G_ZWW | G_WWWW | G_ZZWW | G_PZWW | G_PPWW 
      | G_strong | G_SS | I_G_S | G_S_Sqrt 
      | Gs
      | M of flavor | W of flavor    
      | G_NZN of neu*neu | G_CZC of char*char | G_NNA
      | G_YUK of int*int
      | G_YUK_1 of int*int | G_YUK_2 of int*int | G_YUK_3 of int*int 
      | G_YUK_4 of int*int | G_NHC of neu*char | G_CHN of char*neu
      | G_YUK_C of vc*int*char*sff*sfm
      | G_YUK_Q of vc*int*int*char*sff*sfm
      | G_YUK_N of vc*int*neu*sff*sfm
      | G_YUK_G of vc*int*sff*sfm
      | G_NGC of neu*char | G_CGN of char*neu 
      | SUM_1 
      | G_NWC of neu*char | G_CWN of char*neu
      | G_CH1C of char*char | G_CH2C of char*char | G_CAC of char*char
      | G_CGC of char*char
      | G_SWS of vc*int*int*sfm*sfm
      | G_SLSNW of vc*int*sfm 
      | G_ZSF of sff*int*sfm*sfm
      | G_CICIH1 of neu*neu | G_CICIH2 of neu*neu | G_CICIA of neu*neu
      | G_CICIG of neu*neu 
      | G_GH of int | G_GHGo of int
      | G_GLGLH | G_GLGLHH | G_GLGLA | G_PPH | G_PPHH | G_PPA
      | G_WWSFSF of sff*int*sfm*sfm 
      | G_WPSLSN of vc*int*sfm
      | G_H3 of int | G_H4 of int
      | G_HGo3 of int | G_HGo4 of int | G_GG4 of int
      | G_H1SFSF of sff*int*sfm*sfm | G_H2SFSF of sff*int*sfm*sfm 
      | G_ASFSF of sff*int*sfm*sfm 
      | G_HSNSL of vc*int*sfm  
      | G_GoSFSF of sff*int*sfm*sfm 
      | G_GoSNSL of vc*int*sfm 
      | G_HSUSD of vc*sfm*sfm*int*int | G_GSUSD of vc*sfm*sfm*int*int 
      | G_WPSUSD of vc*sfm*sfm*int*int  
      | G_WZSUSD of vc*sfm*sfm*int*int  
      | G_WZSLSN of vc*int*sfm | G_GlGlSQSQ
      | G_PPSFSF of sff 
      | G_ZZSFSF of sff*int*sfm*sfm | G_ZPSFSF of sff*int*sfm*sfm 
      | G_GlZSFSF of sff*int*sfm*sfm | G_GlPSQSQ 
      | G_GlWSUSD of vc*sfm*sfm*int*int
      | G_GH4 of int | G_GHGo4 of int 
      | G_H1H2SFSF of sff*sfm*sfm*int 
      | G_H1H1SFSF of sff*sfm*sfm*int 
      | G_H2H2SFSF of sff*sfm*sfm*int 
      | G_HHSFSF of sff*sfm*sfm*int 
      | G_AASFSF of sff*sfm*sfm*int 
      | G_HH1SLSN of vc*sfm*int | G_HH2SLSN of vc*sfm*int 
      | G_HASLSN of vc*sfm*int   
      | G_HH1SUSD of vc*sfm*sfm*int*int 
      | G_HH2SUSD of vc*sfm*sfm*int*int 
      | G_HASUSD of vc*sfm*sfm*int*int 
      | G_AG0SFSF of sff*sfm*sfm*int 
      | G_HGSFSF of sff*sfm*sfm*int 
      | G_GGSFSF of sff*sfm*sfm*int 
      | G_G0G0SFSF of sff*sfm*sfm*int
      | G_HGSNSL of vc*sfm*int | G_H1GSNSL of vc*sfm*int 
      | G_H2GSNSL of vc*sfm*int | G_AGSNSL of vc*sfm*int 
      | G_GGSNSL of vc*sfm*int 
      | G_HGSUSD of vc*sfm*sfm*int*int 
      | G_H1GSUSD of vc*sfm*sfm*int*int 
      | G_H2GSUSD of vc*sfm*sfm*int*int 
      | G_AGSUSD of vc*sfm*sfm*int*int 
      | G_GGSUSD of vc*sfm*sfm*int*int 
      | G_SN4 of int*int
      | G_SN2SL2_1 of sfm*sfm*int*int | G_SN2SL2_2 of sfm*sfm*int*int
      | G_SF4 of sff*sff*sfm*sfm*sfm*sfm*int*int
      | G_SF4_3 of sff*sff*sfm*sfm*sfm*sfm*int*int*int
      | G_SF4_4 of sff*sff*sfm*sfm*sfm*sfm*int*int*int*int
      | G_SL4 of sfm*sfm*sfm*sfm*int
      | G_SL4_2 of sfm*sfm*sfm*sfm*int*int
      | G_SN2SQ2 of sff*sfm*sfm*int*int
      | G_SL2SQ2 of sff*sfm*sfm*sfm*sfm*int*int
      | G_SUSDSNSL of vc*sfm*sfm*sfm*int*int*int
      | G_SU4 of sfm*sfm*sfm*sfm*int
      | G_SU4_2 of sfm*sfm*sfm*sfm*int*int
      | G_SD4 of sfm*sfm*sfm*sfm*int
      | G_SD4_2 of sfm*sfm*sfm*sfm*int*int
      | G_SU2SD2 of sfm*sfm*sfm*sfm*int*int*int*int
      | G_HSF31 of higgs*int*sfm*sfm*sff*sff
      | G_HSF32 of higgs*int*int*sfm*sfm*sff*sff
      | G_HSF41 of higgs*int*sfm*sfm*sff*sff
      | G_HSF42 of higgs*int*int*sfm*sfm*sff*sff
      | G_Grav | G_Gr_Ch of char | G_Gr_Z_Neu of neu
      | G_Gr_A_Neu of neu | G_Gr4_Neu of neu 
      | G_Gr4_A_Ch of char | G_Gr4_Z_Ch of char
      | G_Grav_N | G_Grav_U of int*sfm | G_Grav_D of int*sfm 
      | G_Grav_L of int*sfm | G_Grav_Uc of int*sfm | G_Grav_Dc of int*sfm 
      | G_Grav_Lc of int*sfm | G_GravGl
      | G_Gr_H_Ch of char | G_Gr_H1_Neu of neu
      | G_Gr_H2_Neu of neu | G_Gr_H3_Neu of neu
      | G_Gr4A_Sl of int*sfm | G_Gr4A_Slc of int*sfm 
      | G_Gr4A_Su of int*sfm | G_Gr4A_Suc of int*sfm 
      | G_Gr4A_Sd of int*sfm | G_Gr4A_Sdc of int*sfm 
      | G_Gr4Z_Sn | G_Gr4Z_Snc
      | G_Gr4Z_Sl of int*sfm | G_Gr4Z_Slc of int*sfm 
      | G_Gr4Z_Su of int*sfm | G_Gr4Z_Suc of int*sfm 
      | G_Gr4Z_Sd of int*sfm | G_Gr4Z_Sdc of int*sfm 
      | G_Gr4W_Sl of int*sfm | G_Gr4W_Slc of int*sfm 
      | G_Gr4W_Su of int*sfm | G_Gr4W_Suc of int*sfm 
      | G_Gr4W_Sd of int*sfm | G_Gr4W_Sdc of int*sfm 
      | G_Gr4W_Sn | G_Gr4W_Snc
      | G_Gr4Gl_Su of int*sfm | G_Gr4Gl_Suc of int*sfm 
      | G_Gr4Gl_Sd of int*sfm | G_Gr4Gl_Sdc of int*sfm
      | G_Gr4_Z_H1 of neu | G_Gr4_Z_H2 of neu | G_Gr4_Z_H3 of neu
      | G_Gr4_W_H of neu | G_Gr4_W_Hc of neu | G_Gr4_H_A of char
      | G_Gr4_H_Z of char

(* Two integer counters for the QCD and EW order of the couplings. *)

    type orders = int * int

    let orders = function 
      | _ -> (0,0)

    let ferm_of_sff = function
      | SL, g -> (L g) | SN, g -> (N g) 
      | SU, g -> (U g) | SD, g -> (D g)

(* \begin{subequations}
     \begin{align}
        \alpha_{\text{QED}} &= \frac{1}{137.0359895} \\
             \sin^2\theta_w &= 0.23124
     \end{align}
   \end{subequations}

Here we must perhaps allow for complex input parameters. So split them
into their modulus and their phase. At first, we leave them real; the 
generalization to complex parameters is obvious. *)

    module Ch = Charges.QQ

    let ( // ) = Algebra.Small_Rational.make

    let generation' = function
      |  1 -> [ 1//1;  0//1;  0//1]
      |  2 -> [ 0//1;  1//1;  0//1]
      |  3 -> [ 0//1;  0//1;  1//1]
      | -1 -> [-1//1;  0//1;  0//1]
      | -2 -> [ 0//1; -1//1;  0//1]
      | -3 -> [ 0//1;  0//1; -1//1]
      |  n -> invalid_arg ("MSSM.generation': " ^ string_of_int n)

    let generation f =
      if Flags.ckm_present then
        []
      else
        match f with
        | L n | N n | U n | D n | Sup (_,n) 
        | Sdown (_,n) | Slepton (_,n) 
        | Sneutrino n -> generation' n
        | _ -> [0//1; 0//1; 0//1]

    let charge = function
      | L n -> if n > 0 then -1//1 else  1//1
      | Slepton (_,n) -> if n > 0 then -1//1 else  1//1
      | N n -> 0//1
      | Sneutrino n -> 0//1
      | U n -> if n > 0 then  2//3 else -2//3
      | Sup (_,n) -> if n > 0 then  2//3 else -2//3
      | D n -> if n > 0 then -1//3 else  1//3          
      | Sdown (_,n) -> if n > 0 then -1//3 else  1//3          
      | Gl | Ga | Z | Neutralino _ | Gluino -> 0//1
      | Wp ->  1//1
      | Wm -> -1//1
      | H_Heavy | H_Light | Phi0 ->  0//1
      | Hp | Phip ->  1//1
      | Hm | Phim -> -1//1
      | Chargino (C1 | C2) -> 1//1 
      | Chargino (C1c | C2c) -> -1//1 
      | _ -> 0//1

    let lepton = function
      | L n | N n -> if n > 0 then 1//1 else -1//1
      | Slepton (_,n) 
      | Sneutrino n -> if n > 0 then 1//1 else -1//1
      | _ -> 0//1

    let baryon = function
      | U n | D n -> if n > 0 then 1//1 else -1//1
      | Sup (_,n) | Sdown (_,n) -> if n > 0 then 1//1 else -1//1
      | _ -> 0//1

    let charges f =
      [ charge f; lepton f; baryon f] @ generation f

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

(*** REVISED: Compatible with CD+. ***)
    let electromagnetic_currents_3 g = 
     [((U (-g), Ga, U g), FBF (1, Psibar, V, Psi), Q_up);
      ((D (-g), Ga, D g), FBF (1, Psibar, V, Psi), Q_down);
      ((L (-g), Ga, L g), FBF (1, Psibar, V, Psi), Q_lepton) ]
        
(*** REVISED: Compatible with CD+. ***)
    let electromagnetic_sfermion_currents g m =
        [ ((Ga, Slepton (m,-g), Slepton (m,g)), Vector_Scalar_Scalar 1, Q_lepton);
          ((Ga, Sup (m,-g), Sup (m,g)), Vector_Scalar_Scalar 1, Q_up);
          ((Ga, Sdown (m,-g), Sdown (m,g)), Vector_Scalar_Scalar 1, Q_down) ]

(*** REVISED: Compatible with CD+. ***)
    let electromagnetic_currents_2 c =
      let cc = conj_char c in
      [ ((Chargino cc, Ga, Chargino c), FBF (1, Psibar, V, Psi), Q_charg) ]

(*** REVISED: Compatible with CD+. ***)
    let neutral_currents g =
      [ ((L (-g), Z, L g), FBF (1, Psibar, VA, Psi), G_NC_lepton);
        ((N (-g), Z, N g), FBF (1, Psibar, VA, Psi), G_NC_neutrino);
        ((U (-g), Z, U g), FBF (1, Psibar, VA, Psi), G_NC_up);
        ((D (-g), Z, D g), FBF (1, Psibar, VA, Psi), G_NC_down) ] 

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
    let charged_currents g =
      [ ((L (-g), Wm, N g), FBF ((-1), Psibar, VL, Psi), G_CC);
        ((N (-g), Wp, L g), FBF ((-1), Psibar, VL, Psi), G_CC) ]

(* The quark with the inverted generation (the antiparticle) is the outgoing 
   one, the other the incoming. The vertex attached to the outgoing up-quark 
   contains the CKM matrix element {\em not} complex conjugated, while the 
   vertex with the outgoing down-quark has the conjugated CKM matrix 
   element. *)

(*** REVISED: Compatible with CD+. ***)
    let charged_quark_currents g h = 
      [ ((D (-g), Wm, U h), FBF ((-1), Psibar, VL, Psi), G_CCQ (true,g,h));
        ((U (-g), Wp, D h), FBF ((-1), Psibar, VL, Psi), G_CCQ (false,h,g))] 

(*** REVISED: Compatible with CD+. ***)
    let charged_chargino_currents n c =
      let cc = conj_char c in 
      [ ((Chargino cc, Wp, Neutralino n), 
                    FBF (1, Psibar, VLR, Chi), G_CWN (c,n));
        ((Neutralino n, Wm, Chargino c), 
                    FBF (1, Chibar, VLR, Psi), G_NWC (n,c)) ]

(*** REVISED: Compatible with CD+. ***)
    let charged_slepton_currents g m =
      [ ((Wm, Slepton (m,-g), Sneutrino g), Vector_Scalar_Scalar (-1), G_SLSNW 
           (true,g,m));
        ((Wp, Slepton (m,g), Sneutrino (-g)), Vector_Scalar_Scalar 1, G_SLSNW 
           (false,g,m)) ]
 
(*** REVISED: Compatible with CD+. ***)
    let charged_squark_currents' g h m1 m2 =
      [ ((Wm, Sup (m1,g), Sdown (m2,-h)), Vector_Scalar_Scalar (-1), G_SWS 
           (true,g,h,m1,m2));
        ((Wp, Sup (m1,-g), Sdown (m2,h)), Vector_Scalar_Scalar 1, G_SWS 
           (false,g,h,m1,m2)) ]
    let charged_squark_currents g h = List.flatten (Product.list2 
            (charged_squark_currents' g h) [M1;M2] [M1;M2]) 

(*** REVISED: Compatible with CD+. ***)
    let neutral_sfermion_currents' g m1 m2 =
      [ ((Z, Slepton (m1,-g), Slepton (m2,g)), Vector_Scalar_Scalar (-1), G_ZSF 
           (SL,g,m1,m2));
        ((Z, Sup (m1,-g), Sup (m2,g)), Vector_Scalar_Scalar (-1), G_ZSF 
           (SU,g,m1,m2));
        ((Z, Sdown (m1,-g), Sdown (m2,g)), Vector_Scalar_Scalar (-1), G_ZSF 
           (SD,g,m1,m2)) ]
    let neutral_sfermion_currents g = 
      List.flatten (Product.list2 (neutral_sfermion_currents'
                  g) [M1;M2] [M1;M2]) @
      [ ((Z, Sneutrino (-g), Sneutrino g), Vector_Scalar_Scalar (-1), G_ZSF 
           (SN,g,M1,M1)) ]

(* The reality of the coupling of the Z-boson to two identical neutralinos 
   makes the vector part of the coupling vanish. So we distinguish them not 
   by the name but by the structure of the couplings. *)  

(*** REVISED: Compatible with CD+. ***)
    let neutral_Z_1 (n,m) =  
      [ ((Neutralino n, Z, Neutralino m), FBF (1, Chibar, VA, Chi), 
              (G_NZN (n,m))) ]
(*** REVISED: Compatible with CD+. ***)
    let neutral_Z_2 n =
      [ ((Neutralino n, Z, Neutralino n), FBF (1, Chibar, Coupling.A, Chi), 
         (G_NZN (n,n)) )]

(* For very compressed spectra, radiative decays of the next-to-lightest neutralino 
   become important. The formula can be found Haber/Wyler, 1989. In abuse, we 
   include this loop-induced coupling together in the same model variant with the
   triangle Higgs couplings. *)
    let neutral_A =
      if Flags.higgs_triangle then
        [ ((Neutralino N2, Ga, Neutralino N1), FBF (1, Chibar, TVAM, Chi), G_NNA) ]
      else
	[]

(*** REVISED: Compatible with CD+. ***)
    let charged_Z c1 c2 =
      let cc1 = conj_char c1 in
      ((Chargino cc1, Z, Chargino c2), FBF ((-1), Psibar, VA, Psi), 
               G_CZC (c1,c2)) 

(*** REVISED: Compatible with CD+. ***)        
    let yukawa_v =
      [ ((Gluino, Gl, Gluino), FBF (1, Chibar, V, Chi), Gs) ]

(*** REVISED: Independent of the sign of CD. ***)
    let yukawa_higgs g = 
      [ ((N (-g), Hp, L g), FBF (1, Psibar, Coupling.SR, Psi), G_YUK (6,g));
        ((L (-g), Hm, N g), FBF (1, Psibar, Coupling.SL, Psi), G_YUK (6,g));
        ((L (-g), H_Heavy, L g), FBF (1, Psibar, S, Psi), G_YUK (7,g));
        ((L (-g), H_Light, L g), FBF (1, Psibar, S, Psi), G_YUK (8,g));
        ((L (-g), A, L g), FBF (1, Psibar, P, Psi), G_YUK (9,g));
        ((U (-g), H_Heavy, U g), FBF (1, Psibar, S, Psi), G_YUK (10,g));
        ((U (-g), H_Light, U g), FBF (1, Psibar, S, Psi), G_YUK (11,g));
        ((U (-g), A, U g), FBF (1, Psibar, P, Psi), G_YUK (12,g));
        ((D (-g), H_Heavy, D g), FBF (1, Psibar, S, Psi), G_YUK (13,g));
        ((D (-g), H_Light, D g), FBF (1, Psibar, S, Psi), G_YUK (14,g));
        ((D (-g), A, D g), FBF (1, Psibar, P, Psi), G_YUK (15,g)) ]

(*** REVISED: Compatible with CD+ and GS+. ***)
    let yukawa_goldstone g = 
      [ ((N (-g), Phip, L g), FBF (1, Psibar, Coupling.SR, Psi), G_YUK (19,g));
        ((L (-g), Phim, N g), FBF (1, Psibar, Coupling.SL, Psi), G_YUK (19,g));
        ((L (-g), Phi0, L g), FBF (1, Psibar, P, Psi), G_YUK (16,g));
        ((U (-g), Phi0, U g), FBF (1, Psibar, P, Psi), G_YUK (17,g));
        ((D (-g), Phi0, D g), FBF (1, Psibar, P, Psi), G_YUK (18,g)) ]
        
(*** REVISED: Independent of the sign of CD. ***)
    let yukawa_higgs_quark (g,h) =
      [ ((U (-g), Hp, D h), FBF (1, Psibar, SLR, Psi), G_YUK_1 (g, h)); 
        ((D (-h), Hm, U g), FBF (1, Psibar, SLR, Psi), G_YUK_2 (g, h))  ]

(*** REVISED: Compatible with CD+ and GS+. ***)
    let yukawa_goldstone_quark g h =
        [ ((U (-g), Phip, D h), FBF (1, Psibar, SLR, Psi), G_YUK_3 (g, h)); 
          ((D (-h), Phim, U g), FBF (1, Psibar, SLR, Psi), G_YUK_4 (g, h)) ]

(*** REVISED: Compatible with CD+. *)
    let yukawa_higgs_2' (c1,c2) =
      let cc1 = conj_char c1 in
      [ ((Chargino cc1, H_Heavy, Chargino c2), FBF (1, Psibar, SLR, Psi), 
                G_CH2C (c1,c2));
        ((Chargino cc1, H_Light, Chargino c2), FBF (1, Psibar, SLR, Psi),
                G_CH1C (c1,c2));
        ((Chargino cc1, A, Chargino c2), FBF (1, Psibar, SLR, Psi), 
                G_CAC (c1,c2)) ]
    let yukawa_higgs_2'' c =
      let cc = conj_char c in
      [ ((Chargino cc, H_Heavy, Chargino c), FBF (1, Psibar, S, Psi), 
                G_CH2C (c,c));
        ((Chargino cc, H_Light, Chargino c), FBF (1, Psibar, S, Psi),
                G_CH1C (c,c));
        ((Chargino cc, A, Chargino c), FBF (1, Psibar, P, Psi), 
                G_CAC (c,c)) ]
    let yukawa_higgs_2 = 
      ThoList.flatmap yukawa_higgs_2'  [(C1,C2);(C2,C1)] @ 
      ThoList.flatmap yukawa_higgs_2'' [C1;C2] 

(*** REVISED: Compatible with CD+ and GS+. ***)
    let yukawa_goldstone_2' (c1,c2) = 
      let cc1 = conj_char c1 in
      [ ((Chargino cc1, Phi0, Chargino c2), FBF (1, Psibar, SLR, Psi), 
                G_CGC (c1,c2)) ]
    let yukawa_goldstone_2'' c = 
      let cc = conj_char c in 
      [ ((Chargino cc, Phi0, Chargino c), FBF (1, Psibar, P, Psi), 
                G_CGC (c,c)) ]
    let yukawa_goldstone_2 = 
      ThoList.flatmap yukawa_goldstone_2' [(C1,C2);(C2,C1)] @
      ThoList.flatmap yukawa_goldstone_2'' [C1;C2] 

(*** REVISED: Compatible with CD+. ***)
    let higgs_charg_neutr n c =
      let cc = conj_char c in
      [ ((Neutralino n, Hm, Chargino c), FBF (-1, Chibar, SLR, Psi), 
                   G_NHC (n,c));
        ((Chargino cc, Hp, Neutralino n), FBF (-1, Psibar, SLR, Chi), 
                   G_CHN (c,n)) ]

(*** REVISED: Compatible with CD+ and GS+. ***) 
    let goldstone_charg_neutr n c =
      let cc = conj_char c in 
      [ ((Neutralino n, Phim, Chargino c), FBF (1, Chibar, SLR, Psi), 
              G_NGC (n,c));
        ((Chargino cc, Phip, Neutralino n), FBF (1, Psibar, SLR, Chi), 
                   G_CGN (c,n)) ]

(*** REVISED: Compatible with CD+. ***)        
    let higgs_neutr' (n,m) =
      [ ((Neutralino n, H_Heavy, Neutralino m), FBF (1, Chibar, SP, Chi), 
               G_CICIH2 (n,m));
        ((Neutralino n, H_Light, Neutralino m), FBF (1, Chibar, SP, Chi), 
               G_CICIH1 (n,m));
        ((Neutralino n, A, Neutralino m), FBF (1, Chibar, SP, Chi), 
               G_CICIA (n,m)) ]
    let higgs_neutr'' n =
      [ ((Neutralino n, H_Heavy, Neutralino n), FBF (1, Chibar, S, Chi), 
               G_CICIH2 (n,n));
        ((Neutralino n, H_Light, Neutralino n), FBF (1, Chibar, S, Chi), 
               G_CICIH1 (n,n));
        ((Neutralino n, A, Neutralino n), FBF (1, Chibar, P, Chi), 
               G_CICIA (n,n)) ]
    let higgs_neutr = 
      ThoList.flatmap higgs_neutr'  [(N1,N2);(N1,N3);(N1,N4);
                                     (N2,N3);(N2,N4);(N3,N4)] @ 
      ThoList.flatmap higgs_neutr'' [N1;N2;N3;N4] 

(*** REVISED: Compatible with CD+ and GS+. ***) 
    let goldstone_neutr' (n,m) =
      [ ((Neutralino n, Phi0, Neutralino m), FBF (1, Chibar, SP, Chi), 
               G_CICIG (n,m)) ]
    let goldstone_neutr'' n = 
      [ ((Neutralino n, Phi0, Neutralino n), FBF (1, Chibar, P, Chi), 
               G_CICIG (n,n)) ]
    let goldstone_neutr = 
      ThoList.flatmap goldstone_neutr'  [(N1,N2);(N1,N3);(N1,N4);
                                     (N2,N3);(N2,N4);(N3,N4)] @ 
      ThoList.flatmap goldstone_neutr'' [N1;N2;N3;N4] 
              

(*** REVISED: Compatible with CD+. ***)
     let yukawa_n_1 n g = 
         [ ((Neutralino n, Slepton (M1,-g), L g), FBF (1, Chibar, Coupling.SL,
              Psi), G_YUK_N (true,g,n,SL,M1));
           ((Neutralino n, Slepton (M2,-g), L g), FBF (1, Chibar, SR, Psi), 
              G_YUK_N (true,g,n,SL,M2));
           ((L (-g), Slepton (M1,g), Neutralino n), FBF (1, Psibar, SR, Chi), 
              G_YUK_N (false,g,n,SL,M1));
           ((L (-g), Slepton (M2,g), Neutralino n), FBF (1, Psibar, Coupling.SL,
              Chi), G_YUK_N (false,g,n,SL,M2));
           ((Neutralino n, Sup (M1,-g), U g), FBF (1, Chibar, Coupling.SL, 
              Psi), G_YUK_N (true,g,n,SU,M1));
           ((Neutralino n, Sup (M2,-g), U g), FBF (1, Chibar, SR, Psi), 
              G_YUK_N (true,g,n,SU,M2));
           ((U (-g), Sup (M1,g), Neutralino n), FBF (1, Psibar, SR, Chi), 
              G_YUK_N (false,g,n,SU,M1));
           ((U (-g), Sup (M2,g), Neutralino n), FBF (1, Psibar, Coupling.SL, 
              Chi), G_YUK_N (false,g,n,SU,M2));
           ((Neutralino n, Sdown (M1,-g), D g), FBF (1, Chibar, Coupling.SL, 
              Psi), G_YUK_N (true,g,n,SD,M1));
           ((Neutralino n, Sdown (M2,-g), D g), FBF (1, Chibar, SR, Psi), 
              G_YUK_N (true,g,n,SD,M2));
           ((D (-g), Sdown (M1,g), Neutralino n), FBF (1, Psibar, SR, Chi), 
              G_YUK_N (false,g,n,SD,M1));
           ((D (-g), Sdown (M2,g), Neutralino n), FBF (1, Psibar, Coupling.SL, 
              Chi), G_YUK_N (false,g,n,SD,M2)) ]
     let yukawa_n_2 n m = 
         [ ((Neutralino n, Slepton (m,-3), L 3), FBF (1, Chibar, SLR, Psi), 
              G_YUK_N (true,3,n,SL,m));
           ((L (-3), Slepton (m,3), Neutralino n), FBF (1, Psibar, SLR, Chi), 
              G_YUK_N (false,3,n,SL,m));
           ((Neutralino n, Sup (m,-3), U 3), FBF (1, Chibar, SLR, Psi), 
              G_YUK_N (true,3,n,SU,m));
           ((U (-3), Sup (m,3), Neutralino n), FBF (1, Psibar, SLR, Chi), 
              G_YUK_N (false,3,n,SU,m));
           ((Neutralino n, Sdown (m,-3), D 3), FBF (1, Chibar, SLR, Psi), 
              G_YUK_N (true,3,n,SD,m));
           ((D (-3), Sdown (m,3), Neutralino n), FBF (1, Psibar, SLR, Chi), 
              G_YUK_N (false,3,n,SD,m)) ]
     let yukawa_n_3 n g =
         [ ((Neutralino n, Sneutrino (-g), N g), FBF (1, Chibar, Coupling.SL, 
              Psi), G_YUK_N (true,g,n,SN,M1));
           ((N (-g), Sneutrino g, Neutralino n), FBF (1, Psibar, SR, Chi), 
              G_YUK_N (false,g,n,SN,M1)) ]
     let yukawa_n_4 g =
         [ ((U (-g), Sup (M1,g), Gluino), FBF ((-1), Psibar, SR, Chi), G_S_Sqrt);
           ((D (-g), Sdown (M1,g), Gluino), FBF ((-1), Psibar, SR, Chi), G_S_Sqrt);
           ((Gluino, Sup (M1,-g), U g), FBF ((-1), Chibar, Coupling.SL, Psi), G_S_Sqrt);
           ((Gluino, Sdown (M1,-g), D g), FBF ((-1), Chibar, Coupling.SL, Psi), G_S_Sqrt);
           ((U (-g), Sup (M2,g), Gluino), FBF (1, Psibar, Coupling.SL, Chi), G_S_Sqrt);
           ((D (-g), Sdown (M2,g), Gluino), FBF (1, Psibar, Coupling.SL, Chi), G_S_Sqrt);
           ((Gluino, Sup (M2,-g), U g), FBF (1, Chibar, SR, Psi), G_S_Sqrt);
           ((Gluino, Sdown (M2,-g), D g), FBF (1, Chibar, SR, Psi), G_S_Sqrt)]
    let yukawa_n_5 m =
          [ ((U (-3), Sup (m,3), Gluino), FBF (1, Psibar, SLR, Chi), 
                  G_YUK_G (false,3,SU,m));
            ((D (-3), Sdown (m,3), Gluino), FBF (1, Psibar, SLR, Chi), 
                  G_YUK_G (false,3,SD,m));
            ((Gluino, Sup (m,-3), U 3), FBF (1, Chibar, SLR, Psi), 
                  G_YUK_G (true,3,SU,m));
            ((Gluino, Sdown (m,-3), D 3), FBF (1, Chibar, SLR, Psi), 
                  G_YUK_G (true,3,SD,m))]
    let yukawa_n =
      List.flatten (Product.list2 yukawa_n_1 [N1;N2;N3;N4] [1;2]) @
      List.flatten (Product.list2 yukawa_n_2 [N1;N2;N3;N4] [M1;M2]) @
      List.flatten (Product.list2 yukawa_n_3 [N1;N2;N3;N4] [1;2;3]) @
      ThoList.flatmap yukawa_n_4 [1;2] @ 
      ThoList.flatmap yukawa_n_5 [M1;M2]      

(*** REVISED: Compatible with CD+. ***)
    let yukawa_c_1 c g =
         let cc = conj_char c in
         [ ((L (-g), Sneutrino g, Chargino cc), BBB (1, Psibar, Coupling.SR, 
              Psibar), G_YUK_C (true,g,c,SN,M1));
           ((Chargino c, Sneutrino (-g), L g), PBP (1, Psi, Coupling.SL, Psi), 
              G_YUK_C (false,g,c,SN,M1)) ]
    let yukawa_c_2 c = 
         let cc = conj_char c in
         [ ((L (-3), Sneutrino 3, Chargino cc), BBB (1, Psibar, SLR, 
              Psibar), G_YUK_C (true,3,c,SN,M1));
           ((Chargino c, Sneutrino (-3), L 3), PBP (1, Psi, SLR, Psi), 
              G_YUK_C (false,3,c,SN,M1)) ]
    let yukawa_c_3 c m g =
         let cc = conj_char c in
         [ ((N (-g), Slepton (m,g), Chargino c), FBF (1, Psibar, Coupling.SR, 
              Psi), G_YUK_C (true,g,c,SL,m));
           ((Chargino cc, Slepton (m,-g), N g), FBF (1, Psibar, Coupling.SL, 
              Psi), G_YUK_C (false,g,c,SL,m)) ]
    let yukawa_c c = 
      ThoList.flatmap (yukawa_c_1 c) [1;2] @ 
      yukawa_c_2 c @
      List.flatten (Product.list2 (yukawa_c_3 c) [M1] [1;2]) @
      List.flatten (Product.list2 (yukawa_c_3 c) [M1;M2] [3]) 

(*** REVISED: Compatible with CD+. ***)
   let yukawa_cq' c (g,h) m = 
       let cc = conj_char c in
         [ ((Chargino c, Sup (m,-g), D h), PBP (1, Psi, SLR, Psi), 
            G_YUK_Q (false,g,h,c,SU,m));
           ((D (-h), Sup (m,g), Chargino cc), BBB (1, Psibar, SLR, Psibar), 
            G_YUK_Q (true,g,h,c,SU,m));
           ((Chargino cc, Sdown (m,-h), U g), FBF (1, Psibar, SLR, Psi), 
            G_YUK_Q (true,g,h,c,SD,m));
           ((U (-g), Sdown (m,h), Chargino c), FBF (1, Psibar, SLR, Psi), 
            G_YUK_Q (false,g,h,c,SD,m)) ]               
      let yukawa_cq'' c (g,h) =
        let cc = conj_char c in
          [ ((Chargino c, Sup (M1,-g), D h), PBP (1, Psi, Coupling.SL, Psi), 
                G_YUK_Q (false,g,h,c,SU,M1));
            ((D (-h), Sup (M1,g), Chargino cc), 
                BBB (1, Psibar, Coupling.SR, Psibar), G_YUK_Q (true,g,h,c,SU,M1));
            ((Chargino cc, Sdown (M1,-h), U g), 
                FBF (1, Psibar, Coupling.SL, Psi), G_YUK_Q (true,g,h,c,SD,M1));
            ((U (-g), Sdown (M1,h), Chargino c), 
                FBF (1, Psibar, Coupling.SR, Psi), G_YUK_Q (false,g,h,c,SD,M1)) ]
   let yukawa_cq c =      
     if Flags.ckm_present then
       List.flatten (Product.list2 (yukawa_cq' c) [(1,3);(2,3);(3,3);
                                                   (3,2);(3,1)] [M1;M2]) @
       ThoList.flatmap (yukawa_cq'' c) [(1,1);(1,2);(2,1);(2,2)] 
     else
       ThoList.flatmap (yukawa_cq' c (3,3)) [M1;M2] @
       ThoList.flatmap (yukawa_cq'' c) [(1,1);(2,2)]


(*** REVISED: Compatible with CD+. 
   Remark: Singlet and octet gluon exchange. The coupling is divided by
   sqrt(2) to account for the correct normalization of the Lie algebra
   generators.
***)         
    let col_currents g =
      [ ((D (-g), Gl, D g), FBF ((-1), Psibar, V, Psi), Gs);
        ((U (-g), Gl, U g), FBF ((-1), Psibar, V, Psi), Gs)]

(*** REVISED: Compatible with CD+. 
   Remark: Singlet and octet gluon exchange. The coupling is divided by
   sqrt(2) to account for the correct normalization of the Lie algebra
   generators.
***)

   let col_sfermion_currents g m = 
      [ ((Gl, Sup (m,-g), Sup (m,g)), Vector_Scalar_Scalar (-1), Gs);
        ((Gl, Sdown (m,-g), Sdown (m,g)), Vector_Scalar_Scalar (-1), Gs)]


(* The gravitino coupling is generically $1/(4 M_{Pl.})$ *)

(*** Triple vertices containing graivitinos. ***)
    let triple_gravitino' g = 
      [ ((Grino, Sneutrino (-g), N g), GBG (1, Gravbar, Coupling.SL, Psi), G_Grav_N);
        ((N (-g), Sneutrino g, Grino), GBG (1, Psibar, Coupling.SL, Grav), G_Grav_N)]

    let triple_gravitino'' g m = 
      [ ((Grino, Slepton (m, -g), L g), GBG (1, Gravbar, SLR, Psi), G_Grav_L (g,m));
        ((L (-g), Slepton (m, g), Grino), GBG (1, Psibar, SLR, Grav), G_Grav_Lc (g,m));
        ((Grino, Sup (m, -g), U g), GBG (1, Gravbar, SLR, Psi), G_Grav_U (g,m));
        ((U (-g), Sup (m, g), Grino), GBG (1, Psibar, SLR, Grav), G_Grav_Uc (g,m));
        ((Grino, Sdown (m, -g), D g), GBG (1, Gravbar, SLR, Psi), G_Grav_D (g,m));
        ((D (-g), Sdown (m, g), Grino), GBG (1, Psibar, SLR, Grav), G_Grav_Dc (g,m)) ]

    let higgs_ch_gravitino c =
      let cc = conj_char c in      
      [ ((Grino, Hm, Chargino c), GBG (1, Gravbar, SLR, Psi), G_Gr_H_Ch c); 
        ((Chargino cc, Hp, Grino), GBG (1, Psibar, SLR, Grav), G_Gr_H_Ch cc) ]

    let higgs_neu_gravitino n = 
      [ ((Grino, H_Light, Neutralino n), GBG (1, Gravbar, SLR, Chi), G_Gr_H1_Neu n);
        ((Grino, H_Heavy, Neutralino n), GBG (1, Gravbar, SLR, Chi), G_Gr_H2_Neu n);
        ((Grino, A, Neutralino n), GBG (1, Gravbar, SLR, Chi), G_Gr_H3_Neu n) ]

    let gravitino_gaugino_3 = 
      [ ((Grino, Gl, Gluino), GBG (1, Gravbar, V, Chi), G_Grav);
        ((Gluino, Gl, Grino), GBG (1, Chibar, V, Grav), G_Grav);
        ((Chargino C1c, Wp, Grino), GBG (1, Psibar, VLR, Grav), G_Gr_Ch C1);
        ((Chargino C2c, Wp, Grino), GBG (1, Psibar, VLR, Grav), G_Gr_Ch C2);
        ((Grino, Wm, Chargino C1), GBG (1, Gravbar, VLR, Psi), G_Gr_Ch C1c);
        ((Grino, Wm, Chargino C2), GBG (1, Gravbar, VLR, Psi), G_Gr_Ch C2c); 
        ((Grino, Z, Neutralino N1), GBG (1, Gravbar, VLR, Chi), G_Gr_Z_Neu N1);
        ((Grino, Z, Neutralino N2), GBG (1, Gravbar, VLR, Chi), G_Gr_Z_Neu N2);
        ((Grino, Z, Neutralino N3), GBG (1, Gravbar, VLR, Chi), G_Gr_Z_Neu N3);
        ((Grino, Z, Neutralino N4), GBG (1, Gravbar, VLR, Chi), G_Gr_Z_Neu N4);
        ((Grino, Ga, Neutralino N1), GBG (1, Gravbar, VLR, Chi), G_Gr_A_Neu N1);
        ((Grino, Ga, Neutralino N2), GBG (1, Gravbar, VLR, Chi), G_Gr_A_Neu N2);
        ((Grino, Ga, Neutralino N3), GBG (1, Gravbar, VLR, Chi), G_Gr_A_Neu N3);
        ((Grino, Ga, Neutralino N4), GBG (1, Gravbar, VLR, Chi), G_Gr_A_Neu N4) ]
        
    let triple_gravitino = 
      ThoList.flatmap triple_gravitino' [1;2;3] @  
      List.flatten (Product.list2 triple_gravitino'' [1;2;3] [M1; M2]) @  
      ThoList.flatmap higgs_ch_gravitino [C1; C2] @  
      ThoList.flatmap higgs_neu_gravitino [N1; N2; N3; N4] @ 
      gravitino_gaugino_3 


(*** REVISED: Compatible with CD+. ***)
   let triple_gauge =
      [ ((Ga, Wm, Wp), Gauge_Gauge_Gauge 1, I_Q_W);  
        ((Z, Wm, Wp), Gauge_Gauge_Gauge 1, I_G_ZWW);
        ((Gl, Gl, Gl), Gauge_Gauge_Gauge 1, I_G_S)]

(*** REVISED: Independent of the sign of CD. ***) 
   let gauge4 = Vector4 [(2, C_13_42); (-1, C_12_34); (-1, C_14_23)]
   let gluon4 = Vector4 [(-1, C_13_42); (-1, C_12_34); (-1, C_14_23)]       
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

(*** REVISED: Compatible with CD+. ***)
(*** Revision: 2005-03-10: first two vertices corrected. ***)
    let gauge_higgs =
      [ ((Wm, Hp, A), Vector_Scalar_Scalar 1, G_GH 1);
        ((Wp, Hm, A), Vector_Scalar_Scalar 1, G_GH 1);
        ((Z, H_Heavy, A), Vector_Scalar_Scalar 1, G_GH 3);
        ((Z, H_Light, A), Vector_Scalar_Scalar 1, G_GH 2);
        ((H_Heavy, Wp, Wm), Scalar_Vector_Vector 1, G_GH 5);
        ((H_Light, Wp, Wm), Scalar_Vector_Vector 1, G_GH 4);
        ((Wm, Hp, H_Heavy), Vector_Scalar_Scalar 1, G_GH 7);
        ((Wp, Hm, H_Heavy), Vector_Scalar_Scalar (-1), G_GH 7);
        ((Wm, Hp, H_Light), Vector_Scalar_Scalar 1, G_GH 6);
        ((Wp, Hm, H_Light), Vector_Scalar_Scalar (-1), G_GH 6);        
        ((H_Heavy, Z, Z), Scalar_Vector_Vector 1, G_GH 9);
        ((H_Light, Z, Z), Scalar_Vector_Vector 1, G_GH 8);
        ((Z, Hp, Hm), Vector_Scalar_Scalar 1, G_GH 10);
        ((Ga, Hp, Hm), Vector_Scalar_Scalar 1, G_GH 11) ] @
      (if Flags.higgs_triangle then
       [((H_Light, Gl, Gl), Dim5_Scalar_Gauge2 1, G_GLGLH);
        ((H_Heavy, Gl, Gl), Dim5_Scalar_Gauge2 1, G_GLGLHH);
        ((A, Gl, Gl), Dim5_Scalar_Gauge2_Skew 1, G_GLGLA);
        ((H_Light, Ga, Ga), Dim5_Scalar_Gauge2 1, G_PPH);
        ((H_Heavy, Ga, Ga), Dim5_Scalar_Gauge2 1, G_PPHH);
        ((A, Ga, Ga), Dim5_Scalar_Gauge2 1, G_PPA)]
       else
         [])        

(*** REVISED: Compatible with CD+ and GS+. ***)
    let gauge_higgs_gold =
      [ ((Wp, Phi0, Phim), Vector_Scalar_Scalar 1, G_GH 1);
        ((Wm, Phi0, Phip), Vector_Scalar_Scalar 1, G_GH 1);
        ((Z, H_Heavy, Phi0), Vector_Scalar_Scalar 1, G_GH 2);
        ((Z, H_Light, Phi0), Vector_Scalar_Scalar (-1), G_GH 3);
        ((Wp, H_Heavy, Phim), Vector_Scalar_Scalar 1, G_GH 6);
        ((Wm, H_Heavy, Phip), Vector_Scalar_Scalar (-1), G_GH 6);
        ((Wp, H_Light, Phim), Vector_Scalar_Scalar (-1), G_GH 7);
        ((Wm, H_Light, Phip), Vector_Scalar_Scalar 1, G_GH 7);        
        ((Phim, Wp, Ga), Scalar_Vector_Vector 1, G_GHGo 1);
        ((Phip, Wm, Ga), Scalar_Vector_Vector 1, G_GHGo 1);
        ((Phim, Wp, Z), Scalar_Vector_Vector 1, G_GHGo 2);
        ((Phip, Wm, Z), Scalar_Vector_Vector 1, G_GHGo 2);
        ((Z, Phip, Phim), Vector_Scalar_Scalar 1, G_GH 10);
        ((Ga, Phip, Phim), Vector_Scalar_Scalar 1, G_GH 11) ]

    let gauge_higgs4 = 
      [ ((A, A, Z, Z), Scalar2_Vector2 1, G_GH4 1);
        ((H_Heavy, H_Heavy, Z, Z), Scalar2_Vector2 1, G_GH4 3);
        ((H_Light, H_Light, Z, Z), Scalar2_Vector2 1, G_GH4 2);
        ((Hp, Hm, Z, Z), Scalar2_Vector2 1, G_GH4 4);
        ((Hp, Hm, Ga, Ga), Scalar2_Vector2 1, G_GH4 5);
        ((Hp, Hm, Ga, Z), Scalar2_Vector2 1, G_GH4 6);
        ((Hp, H_Heavy, Wm, Z), Scalar2_Vector2 1, G_GH4 8);
        ((Hm, H_Heavy, Wp, Z), Scalar2_Vector2 1, G_GH4 8);
        ((Hp, H_Light, Wm, Z), Scalar2_Vector2 1, G_GH4 7);
        ((Hm, H_Light, Wp, Z), Scalar2_Vector2 1, G_GH4 7);
        ((Hp, H_Heavy, Wm, Ga), Scalar2_Vector2 1, G_GH4 10);
        ((Hm, H_Heavy, Wp, Ga), Scalar2_Vector2 1, G_GH4 10);
        ((Hp, H_Light, Wm, Ga), Scalar2_Vector2 1, G_GH4 9);
        ((Hm, H_Light, Wp, Ga), Scalar2_Vector2 1, G_GH4 9);
        ((A, A, Wp, Wm), Scalar2_Vector2 1, G_GH4 11); 
        ((H_Heavy, H_Heavy, Wp, Wm), Scalar2_Vector2 1, G_GH4 13);
        ((H_Light, H_Light, Wp, Wm), Scalar2_Vector2 1, G_GH4 12);
        ((Hp, Hm, Wp, Wm), Scalar2_Vector2 1, G_GH4 14);
        ((Hp, A, Wm, Z), Scalar2_Vector2 1, G_GH4 15);
        ((Hm, A, Wp, Z), Scalar2_Vector2 (-1), G_GH4 15);
        ((Hp, A, Wm, Ga), Scalar2_Vector2 1, G_GH4 16);
        ((Hm, A, Wp, Ga), Scalar2_Vector2 (-1), G_GH4 16) ]

    let gauge_higgs_gold4 =
      [ ((Z, Z, Phi0, Phi0), Scalar2_Vector2 1, G_GHGo4 1);
        ((Z, Z, Phip, Phim), Scalar2_Vector2 1, G_GHGo4 2);
        ((Ga, Ga, Phip, Phim), Scalar2_Vector2 1, G_GHGo4 3);
        ((Z, Ga, Phip, Phim), Scalar2_Vector2 1, G_GHGo4 4);
        ((Wp, Wm, Phip, Phim), Scalar2_Vector2 1, G_GHGo4 5);
        ((Wp, Wm, Phi0, Phi0), Scalar2_Vector2 1, G_GHGo4 5);
        ((Wp, Z, Phim, Phi0), Scalar2_Vector2 1, G_GHGo4 6);
        ((Wm, Z, Phip, Phi0), Scalar2_Vector2 (-1), G_GHGo4 6);
        ((Wp, Ga, Phim, Phi0), Scalar2_Vector2 1, G_GHGo4 7);
        ((Wm, Ga, Phip, Phi0), Scalar2_Vector2 (-1), G_GHGo4 7);
        ((Wp, Z, Phim, H_Heavy), Scalar2_Vector2 1, G_GHGo4 9);
        ((Wm, Z, Phip, H_Heavy), Scalar2_Vector2 1, G_GHGo4 9);
        ((Wp, Ga, Phim, H_Heavy), Scalar2_Vector2 1, G_GHGo4 11);
        ((Wm, Ga, Phip, H_Heavy), Scalar2_Vector2 1, G_GHGo4 11);
        ((Wp, Z, Phim, H_Light), Scalar2_Vector2 1, G_GHGo4 8);
        ((Wm, Z, Phip, H_Light), Scalar2_Vector2 1, G_GHGo4 8);
        ((Wp, Ga, Phim, H_Light), Scalar2_Vector2 1, G_GHGo4 10);
        ((Wm, Ga, Phip, H_Light), Scalar2_Vector2 1, G_GHGo4 10) ]

    let gauge_sfermion4' g m1 m2 =
       [ ((Wp, Wm, Slepton (m1,g), Slepton (m2,-g)), Scalar2_Vector2 1, 
            G_WWSFSF (SL,g,m1,m2));
        ((Z, Ga, Slepton (m1,g), Slepton (m2,-g)), Scalar2_Vector2 1, 
           G_ZPSFSF (SL,g,m1,m2));
        ((Z, Z, Slepton (m1,g), Slepton (m2,-g)), Scalar2_Vector2 1, G_ZZSFSF 
           (SL,g,m1,m2));
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
      [ ((Wp, Ga, Slepton (m,g), Sneutrino (-g)), Scalar2_Vector2 1, G_WPSLSN 
           (false,g,m));
        ((Wm, Ga, Slepton (m,-g), Sneutrino g), Scalar2_Vector2 1, 
           G_WPSLSN (true,g,m));
        ((Wp, Z, Slepton (m,g), Sneutrino (-g)), Scalar2_Vector2 1, G_WZSLSN 
           (false,g,m));
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

    let gluon_gauge_squark' g m1 m2 =
      [ ((Gl, Z, Sup (m1,g), Sup (m2,-g)), Scalar2_Vector2 2, G_GlZSFSF (SU,g,m1,m2));
        ((Gl, Z, Sdown (m1,g), Sdown (m2,-g)), Scalar2_Vector2 2, G_GlZSFSF (SD,g,m1,m2)) ]
    let gluon_gauge_squark'' g m =
      [ ((Gl, Ga, Sup (m,g), Sup (m,-g)), Scalar2_Vector2 2, G_GlPSQSQ);
        ((Gl, Ga, Sdown (m,g), Sdown (m,-g)), Scalar2_Vector2 (-1), G_GlPSQSQ) ]
    let gluon_gauge_squark g =
      List.flatten (Product.list2 (gluon_gauge_squark' g) [M1;M2] [M1;M2]) @
      ThoList.flatmap (gluon_gauge_squark'' g) [M1;M2]

    let gluon2_squark2 g m =
      [ ((Gl, Gl, Sup (m,g), Sup (m,-g)), Scalar2_Vector2 1, G_GlGlSQSQ);
        ((Gl, Gl, Sdown (m,g), Sdown (m,-g)), Scalar2_Vector2 1, G_GlGlSQSQ)]

(*** REVISED: Independent of the sign of CD. ***)
    let higgs =
      [ ((Hp, Hm, H_Heavy), Scalar_Scalar_Scalar 1, G_H3 1);
        ((Hp, Hm, H_Light), Scalar_Scalar_Scalar 1, G_H3 2);
        ((H_Heavy, H_Heavy, H_Light), Scalar_Scalar_Scalar 1, G_H3 3);
        ((H_Heavy, H_Heavy, H_Heavy), Scalar_Scalar_Scalar 1, G_H3 4);
        ((H_Light, H_Light, H_Light), Scalar_Scalar_Scalar 1, G_H3 5);
        ((H_Heavy, H_Light, H_Light), Scalar_Scalar_Scalar 1, G_H3 6);
        ((H_Heavy, A, A), Scalar_Scalar_Scalar 1, G_H3 7);
        ((H_Light, A, A), Scalar_Scalar_Scalar 1, G_H3 8) ]

(*** REVISED: Compatible with GS+, independent of the sign of CD. ***)
    let higgs_gold = 
      [ ((H_Heavy, A, Phi0), Scalar_Scalar_Scalar 1, G_HGo3 1);
        ((H_Light, A, Phi0), Scalar_Scalar_Scalar 1, G_HGo3 2);
        ((H_Heavy, Hp, Phim), Scalar_Scalar_Scalar 1, G_HGo3 3);
        ((H_Heavy, Hm, Phip), Scalar_Scalar_Scalar 1, G_HGo3 3);
        ((H_Light, Hp, Phim), Scalar_Scalar_Scalar 1, G_HGo3 4);
        ((H_Light, Hm, Phip), Scalar_Scalar_Scalar 1, G_HGo3 4);
        ((A, Hp, Phim), Scalar_Scalar_Scalar (-1), G_HGo3 5);
        ((A, Hm, Phip), Scalar_Scalar_Scalar 1, G_HGo3 5);   
        ((H_Heavy, Phi0, Phi0), Scalar_Scalar_Scalar (-1), G_H3 7);
        ((H_Heavy, Phip, Phim), Scalar_Scalar_Scalar (-1), G_H3 7);
        ((H_Light, Phi0, Phi0), Scalar_Scalar_Scalar (-1), G_H3 8);
        ((H_Light, Phip, Phim), Scalar_Scalar_Scalar (-1), G_H3 8) ]

(* Here follow purely scalar quartic vertices which are only available for the
   no-Whizard colored version. *)

(*** REVISED: Independent of the sign of CD. ***)
    let higgs4 =
      [ ((Hp, Hm, Hp, Hm), Scalar4 1, G_H4 1);
        ((Hp, Hm, H_Heavy, H_Heavy), Scalar4 1, G_H4 2);
        ((Hp, Hm, H_Light, H_Light), Scalar4 1, G_H4 3);
        ((Hp, Hm, H_Heavy, H_Light), Scalar4 1, G_H4 4);
        ((Hp, Hm, A, A), Scalar4 1, G_H4 5);
        ((H_Heavy, H_Heavy, H_Heavy, H_Heavy), Scalar4 1, G_H4 6);
        ((H_Light, H_Light, H_Light, H_Light), Scalar4 1, G_H4 6);
        ((H_Heavy, H_Heavy, H_Light, H_Light), Scalar4 1, G_H4 7);
        ((H_Heavy, H_Light, H_Light, H_Light), Scalar4 1, G_H4 8);
        ((H_Heavy, H_Heavy, H_Heavy, H_Light), Scalar4 (-1), G_H4 8);
        ((H_Heavy, H_Heavy, A, A), Scalar4 1, G_H4 9);
        ((H_Light, H_Light, A, A), Scalar4 (-1), G_H4 9);
        ((H_Heavy, H_Light, A, A), Scalar4 1, G_H4 10);
        ((A, A, A, A), Scalar4 1, G_H4 11) ]
        
(*** REVISED: Compatible with GS+, independent of the sign of CD. ***)
    let higgs_gold4 =
      [ ((H_Heavy, H_Heavy, A, Phi0), Scalar4 1, G_HGo4 1);
        ((H_Heavy, H_Light, A, Phi0), Scalar4 1, G_HGo4 2);
        ((H_Light, H_Light, A, Phi0), Scalar4 (-1), G_HGo4 1);
        ((A, A, A, Phi0), Scalar4 3, G_HGo4 3);
        ((Hp, Hm, A, Phi0), Scalar4 1, G_HGo4 3);
        ((H_Heavy, H_Heavy, Hp, Phim), Scalar4 1, G_HGo4 4);
        ((H_Heavy, H_Heavy, Hm, Phip), Scalar4 1, G_HGo4 4);
        ((H_Heavy, H_Light, Hp, Phim), Scalar4 1, G_HGo4 5);
        ((H_Heavy, H_Light, Hm, Phip), Scalar4 1, G_HGo4 5);
        ((H_Light, H_Light, Hp, Phim), Scalar4 (-1), G_HGo4 4);
        ((H_Light, H_Light, Hm, Phip), Scalar4 (-1), G_HGo4 4);
        ((A, A, Hp, Phim), Scalar4 1, G_HGo4 6);
        ((A, A, Hm, Phip), Scalar4 1, G_HGo4 6);
        ((H_Heavy, A, Hp, Phim), Scalar4 1, G_HGo4 7);
        ((H_Heavy, A, Hm, Phip), Scalar4 (-1), G_HGo4 7);
        ((H_Light, A, Hp, Phim), Scalar4 1, G_HGo4 8);
        ((H_Light, A, Hm, Phip), Scalar4 (-1), G_HGo4 8);
        ((Hp, Hm, Hp, Phim), Scalar4 2, G_HGo4 6);
        ((Hp, Hm, Hm, Phip), Scalar4 2, G_HGo4 6);
        ((H_Heavy, H_Heavy, Phi0, Phi0), Scalar4 (-1), G_H4 9);
        ((H_Heavy, H_Light, Phi0, Phi0), Scalar4 (-1), G_H4 10);
        ((H_Light, H_Light, Phi0, Phi0), Scalar4 1, G_H4 9);
        ((A, A, Phi0, Phi0), Scalar4 1, G_HGo4 9);
        ((Hp, Hm, Phi0, Phi0), Scalar4 1, G_HGo4 10);
        ((H_Heavy, Hp, Phim, Phi0), Scalar4 1, G_HGo4 8);
        ((H_Heavy, Hm, Phip, Phi0), Scalar4 (-1), G_HGo4 8);
        ((H_Light, Hp, Phim, Phi0), Scalar4 (-1), G_HGo4 7);
        ((H_Light, Hm, Phip, Phi0), Scalar4 1, G_HGo4 7);
        ((A, Hp, Phim, Phi0), Scalar4 1, G_HGo4 11);
        ((A, Hm, Phip, Phi0), Scalar4 1, G_HGo4 11);
        ((H_Heavy, H_Heavy, Phip, Phim), Scalar4 1, G_HGo4 12);
        ((H_Heavy, H_Light, Phip, Phim), Scalar4 1, G_HGo4 13);
        ((H_Light, H_Light, Phip, Phim), Scalar4 1, G_HGo4 14);
        ((A, A, Phip, Phim), Scalar4 1, G_HGo4 15);
        ((Hp, Hm, Phip, Phim), Scalar4 1, G_HGo4 16);
        ((Hp, Hp, Phim, Phim), Scalar4 1, G_HGo4 17);
        ((Hm, Hm, Phip, Phip), Scalar4 1, G_HGo4 17);
        ((Hp, Phim, Phi0, Phi0), Scalar4 (-1), G_HGo4 6);
        ((Hm, Phip, Phi0, Phi0), Scalar4 (-1), G_HGo4 6);
        ((A, Phi0, Phi0, Phi0), Scalar4 (-3), G_HGo4 6);
        ((A, Phi0, Phip, Phim), Scalar4 (-1), G_HGo4 6);
        ((Hp, Phim, Phip, Phim), Scalar4 (-2), G_HGo4 6);
        ((Hm, Phip, Phip, Phim), Scalar4 (-2), G_HGo4 6) ]

(*** REVISED: Independent of the sign of CD and GS. ***)
    let goldstone4 =
      [ ((Phi0, Phi0, Phi0, Phi0), Scalar4 1, G_GG4 1);
        ((Phip, Phim, Phi0, Phi0), Scalar4 1, G_GG4 2);
        ((Phip, Phim, Phip, Phim), Scalar4 1, G_GG4 3) ]
    
(* The vertices of the type Higgs - Sfermion - Sfermion are independent of 
   the choice of the CD sign since they are quadratic in the gauge 
   coupling. *) 

(*** REVISED: Independent of the sign of CD. ***)
    let higgs_sneutrino' g =
      [ ((H_Heavy, Sneutrino g, Sneutrino (-g)), Scalar_Scalar_Scalar 1, 
                       G_H2SFSF (SN,g,M1,M1));
        ((H_Light, Sneutrino g, Sneutrino (-g)), Scalar_Scalar_Scalar 1, 
                       G_H1SFSF (SN,g,M1,M1));
        ((Hp, Sneutrino (-g), Slepton (M1,g)), Scalar_Scalar_Scalar 1, 
              G_HSNSL (false,g,M1)); 
        ((Hm, Sneutrino g, Slepton (M1,-g)), Scalar_Scalar_Scalar 1, 
              G_HSNSL (true,g,M1)) ]
      let higgs_sneutrino'' = 
      [ ((Hp, Sneutrino (-3), Slepton (M2,3)), Scalar_Scalar_Scalar 1, 
              G_HSNSL (false,3,M2)); 
        ((Hm, Sneutrino 3, Slepton (M2,-3)), Scalar_Scalar_Scalar 1, 
              G_HSNSL (false,3,M2)) ]
      let higgs_sneutrino = 
        ThoList.flatmap higgs_sneutrino' [1;2;3] @ higgs_sneutrino''  
        

(* Under the assumption that there is no mixing between the left- and
   right-handed sfermions for the first two generations there is only a 
   coupling of the form Higgs - sfermion1 - sfermion2 for the third 
   generation. All the others are suppressed by $m_f/M_W$. *)

(*** REVISED: Independent of the sign of CD. ***)
    let higgs_sfermion' g m1 m2 =
      [ ((H_Heavy, Slepton (m1,g), Slepton (m2,-g)), Scalar_Scalar_Scalar 1,
            G_H2SFSF (SL,g,m1,m2));
        ((H_Light, Slepton (m1,g), Slepton (m2,-g)), Scalar_Scalar_Scalar 1,
            G_H1SFSF (SL,g,m1,m2));
        ((H_Heavy, Sup (m1,g), Sup (m2,-g)), Scalar_Scalar_Scalar 1, 
            G_H2SFSF (SU,g,m1,m2));
        ((H_Heavy, Sdown (m1,g), Sdown (m2,-g)), Scalar_Scalar_Scalar 1, 
            G_H2SFSF (SD,g,m1,m2));
        ((H_Light, Sup (m1,g), Sup (m2,-g)), Scalar_Scalar_Scalar 1, 
            G_H1SFSF (SU,g,m1,m2));
        ((H_Light, Sdown (m1,g), Sdown (m2,-g)), Scalar_Scalar_Scalar 1, 
            G_H1SFSF (SD,g,m1,m2)) ]
      let higgs_sfermion'' m1 m2 = 
        [ ((A, Slepton (m1,3), Slepton (m2,-3)), Scalar_Scalar_Scalar 1,
           G_ASFSF (SL,3,m1,m2));
          ((A, Sup (m1,3), Sup (m2,-3)), Scalar_Scalar_Scalar 1, 
            G_ASFSF (SU,3,m1,m2));
          ((A, Sdown (m1,3), Sdown (m2,-3)), Scalar_Scalar_Scalar 1, 
            G_ASFSF (SD,3,m1,m2)) ]
    let higgs_sfermion = List.flatten (Product.list2 (higgs_sfermion' 3) 
                                           [M1;M2] [M1;M2]) @ 
        (higgs_sfermion' 1 M1 M1) @ (higgs_sfermion' 1 M2 M2) @
        (higgs_sfermion' 2 M1 M1) @ (higgs_sfermion' 2 M2 M2) @ 
        List.flatten (Product.list2 higgs_sfermion'' [M1;M2] [M1;M2]) 
   
(*i    let higgs_sfermion g = List.flatten (Product.list2 (higgs_sfermion' g) 
                                           [M1;M2] [M1;M2])  i*)
     
(*** REVISED: Independent of the sign of CD, compatible with GS+. ***)
    let goldstone_sfermion' g m1 m2 =
      [ ((Phi0, Slepton (m1,g), Slepton (m2,-g)), Scalar_Scalar_Scalar 1,
                       G_GoSFSF (SL,g,m1,m2));
        ((Phi0, Sup (m1,g), Sup (m2,-g)), Scalar_Scalar_Scalar 1, 
                       G_GoSFSF (SU,g,m1,m2));
        ((Phi0, Sdown (m1,g), Sdown (m2,-g)), Scalar_Scalar_Scalar 1, 
                       G_GoSFSF (SD,g,m1,m2))]
    let goldstone_sfermion'' g =
      [ ((Phip, Sneutrino (-g), Slepton (M1,g)), Scalar_Scalar_Scalar 1, 
                       G_GoSNSL (false,g,M1)); 
        ((Phim, Sneutrino g, Slepton (M1,-g)), Scalar_Scalar_Scalar 1, 
                       G_GoSNSL (true,g,M1)) ]
    let goldstone_sfermion''' g = 
      [ ((Phip, Sneutrino (-g), Slepton (M2,g)), Scalar_Scalar_Scalar 1, 
                       G_GoSNSL (false,g,M2)); 
        ((Phim, Sneutrino g, Slepton (M2,-g)), Scalar_Scalar_Scalar 1, 
                       G_GoSNSL (true,g,M2))]
    let goldstone_sfermion = 
      List.flatten (Product.list2 (goldstone_sfermion' 3) [M1;M2] [M1;M2]) @
      ThoList.flatmap goldstone_sfermion'' [1;2;3] @ 
      goldstone_sfermion''' 3 

(*** REVISED: Independent of the sign of CD. ***)
    let higgs_squark' g h m1 m2 =
      [ ((Hp, Sup (m1,-g), Sdown (m2,h)), Scalar_Scalar_Scalar 1, 
              G_HSUSD (false,m1,m2,g,h)); 
        ((Hm, Sup (m1,g), Sdown (m2,-h)), Scalar_Scalar_Scalar 1, 
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

(*** REVISED: Independent of the sign of CD, compatible with GS+. ***)
    let goldstone_squark' g h m1 m2 =
      [ ((Phip, Sup (m1,-g), Sdown (m2,h)), Scalar_Scalar_Scalar 1, 
              G_GSUSD (false,m1,m2,g,h)); 
        ((Phim, Sup (m1,g), Sdown (m2,-h)), Scalar_Scalar_Scalar 1, 
              G_GSUSD (true,m1,m2,g,h)) ]
    let goldstone_squark_a g h = goldstone_squark' g h M1 M1 
    let goldstone_squark_b (g,h) = List.flatten (Product.list2 
            (goldstone_squark' g h) [M1;M2] [M1;M2]) 
    let goldstone_squark =          
         List.flatten (Product.list2 goldstone_squark_a [1;2] [1;2]) @ 
         ThoList.flatmap goldstone_squark_b [(1,3);(2,3);(3,3);(3,1);(3,2)] 

(* BAUSTELLE: For the quartic scalar coupligs we does not allow [whiz_col]. *)

    let higgs_sneutrino4' g m =
      [ ((Hp, H_Heavy, Slepton (m,g), Sneutrino (-g)), Scalar4 1, 
              G_HH2SLSN (false,m,g));
        ((Hm, H_Heavy, Slepton (m,-g), Sneutrino g), Scalar4 1, 
              G_HH2SLSN (true,m,g));
        ((Hp, H_Light, Slepton (m,g), Sneutrino (-g)), Scalar4 1, 
              G_HH1SLSN (false,m,g));
        ((Hm, H_Light, Slepton (m,-g), Sneutrino g), Scalar4 1, 
              G_HH1SLSN (true,m,g));
        ((Hp, A, Slepton (m,g), Sneutrino (-g)), Scalar4 1, 
              G_HASLSN (false,m,g));
        ((Hm, A, Slepton (m,-g), Sneutrino g), Scalar4 1, 
              G_HASLSN (true,m,g)) ]
    let higgs_sneutrino4 g = 
      ThoList.flatmap (higgs_sneutrino4' g) [M1;M2] @
      [ ((H_Heavy, H_Heavy, Sneutrino g, Sneutrino (-g)), Scalar4 1, 
            G_H2H2SFSF (SN,M1,M1,g));
        ((H_Heavy, H_Light, Sneutrino g, Sneutrino (-g)), Scalar4 1, 
            G_H1H2SFSF (SN,M1,M1,g));
        ((H_Light, H_Light, Sneutrino g, Sneutrino (-g)), Scalar4 1, 
            G_H1H1SFSF (SN,M1,M1,g));
        ((Hp, Hm, Sneutrino g, Sneutrino (-g)), Scalar4 1, G_HHSFSF (SN,M1,M1,g)) ]
        
    let higgs_sfermion4' g m1 m2 =
      [ ((H_Heavy, H_Heavy, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
            G_H2H2SFSF (SL,m1,m2,g));
        ((H_Heavy, H_Light, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
            G_H1H2SFSF (SL,m1,m2,g));
        ((H_Light, H_Light, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
            G_H1H1SFSF (SL,m1,m2,g));
        ((A, A, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
            G_AASFSF (SL,m1,m2,g));
        ((Hp, Hm, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
            G_HHSFSF (SL,m1,m2,g));
        ((H_Heavy, H_Heavy, Sup (m1,g), Sup (m2,-g)), Scalar4 1, 
            G_H2H2SFSF (SU,m1,m2,g));
        ((H_Heavy, H_Heavy, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, 
            G_H2H2SFSF (SD,m1,m2,g));
        ((H_Light, H_Light, Sup (m1,g), Sup (m2,-g)), Scalar4 1, 
            G_H1H1SFSF (SU,m1,m2,g));
        ((H_Light, H_Light, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, 
            G_H1H1SFSF (SD,m1,m2,g));
        ((H_Light, H_Heavy, Sup (m1,g), Sup (m2,-g)), Scalar4 1, 
            G_H1H2SFSF (SU,m1,m2,g));
        ((H_Light, H_Heavy, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, 
            G_H1H2SFSF (SD,m1,m2,g));
        ((Hp, Hm, Sup (m1,g), Sup (m2,-g)), Scalar4 1, G_HHSFSF (SU,m1,m2,g));
        ((Hp, Hm, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, G_HHSFSF (SD,m1,m2,g));
        ((A, A, Sup (m1,g), Sup (m2,-g)), Scalar4 1, G_AASFSF (SU,m1,m2,g));
        ((A, A, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, G_AASFSF (SD,m1,m2,g)) ]
    let higgs_sfermion4 g = List.flatten (Product.list2 (higgs_sfermion4' g) 
                                                [M1;M2] [M1;M2])

    let higgs_squark4' g h m1 m2 =
       [ ((Hp, H_Light, Sup (m1,-g), Sdown (m2,h)), Scalar4 1, 
            G_HH1SUSD (false,m1,m2,g,h));
         ((Hm, H_Light, Sup (m1,g), Sdown (m2,-h)), Scalar4 1, 
            G_HH1SUSD (true,m1,m2,g,h));
         ((Hp, H_Heavy, Sup (m1,-g), Sdown (m2,h)), Scalar4 1, 
            G_HH2SUSD (false,m1,m2,g,h));
         ((Hm, H_Heavy, Sup (m1,g), Sdown (m2,-h)), Scalar4 1, 
            G_HH2SUSD (true,m1,m2,g,h));
         ((Hp, A, Sup (m1,-g), Sdown (m2,h)), Scalar4 1, 
            G_HASUSD (false,m1,m2,g,h));
         ((Hm, A, Sup (m1,g), Sdown (m2,-h)), Scalar4 1, 
            G_HASUSD (true,m1,m2,g,h)) ]
    let higgs_squark4 g h = List.flatten (Product.list2 (higgs_squark4' g h) 
                                                [M1;M2] [M1;M2])
         
    let higgs_gold_sneutrino' g m =
      [ ((Hp, Phi0, Sneutrino (-g), Slepton (m,g)), Scalar4 1, G_HGSNSL (false,m,g));
        ((Hm, Phi0, Sneutrino g, Slepton (m,-g)), Scalar4 1, G_HGSNSL (true,m,g));
        ((H_Heavy, Phip, Sneutrino (-g), Slepton (m,g)), Scalar4 1, 
            G_H2GSNSL (false,m,g));
        ((H_Heavy, Phim, Sneutrino g, Slepton (m,-g)), Scalar4 1, 
            G_H2GSNSL (true,m,g));
        ((H_Light, Phip, Sneutrino (-g), Slepton (m,g)), Scalar4 1, 
            G_H1GSNSL (false,m,g));
        ((H_Light, Phim, Sneutrino g, Slepton (m,-g)), Scalar4 1, 
            G_H1GSNSL (true,m,g));
        ((A, Phip, Sneutrino (-g), Slepton (m,g)), Scalar4 1, G_AGSNSL (false,m,g));
        ((A, Phim, Sneutrino g, Slepton (m,-g)), Scalar4 1, G_AGSNSL (true,m,g));
        ((Phi0, Phip, Sneutrino (-g), Slepton (m,g)), Scalar4 1, G_GGSNSL (false,m,g));
        ((Phi0, Phim, Sneutrino g, Slepton (m,-g)), Scalar4 1, G_GGSNSL (true,m,g))]
    let higgs_gold_sneutrino g =
      ThoList.flatmap (higgs_gold_sneutrino' g) [M1;M2] @
       [ ((A, Phi0, Sneutrino g, Sneutrino (-g)), Scalar4 1, 
             G_AG0SFSF (SN,M1,M1,g));
         ((Hp, Phim, Sneutrino g, Sneutrino (-g)), Scalar4 1, 
             G_HGSFSF (SN,M1,M1,g));
         ((Hm, Phip, Sneutrino g, Sneutrino (-g)), Scalar4 1, 
             G_HGSFSF (SN,M1,M1,g));
         ((Phip, Phim, Sneutrino g, Sneutrino (-g)), Scalar4 1, 
             G_GGSFSF (SN,M1,M1,g));
         ((Phi0, Phi0, Sneutrino g, Sneutrino (-g)), Scalar4 1, 
             G_G0G0SFSF (SN,M1,M1,g)) ]
 
    let higgs_gold_sfermion' g m1 m2 =
       [ ((A, Phi0, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
             G_AG0SFSF (SL,m1,m2,g));
         ((Hp, Phim, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
             G_HGSFSF (SL,m1,m2,g));
         ((Hm, Phip, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
             G_HGSFSF (SL,m1,m2,g));
         ((Phip, Phim, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
             G_GGSFSF (SL,m1,m2,g));
         ((Phi0, Phi0, Slepton (m1,g), Slepton (m2,-g)), Scalar4 1, 
             G_G0G0SFSF (SL,m1,m2,g));
         ((A, Phi0, Sup (m1,g), Sup (m2,-g)), Scalar4 1, G_AG0SFSF (SU,m1,m2,g));
         ((A, Phi0, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, 
             G_AG0SFSF (SD,m1,m2,g));
         ((Hp, Phim, Sup (m1,g), Sup (m2,-g)), Scalar4 1, G_HGSFSF (SU,m1,m2,g));
         ((Hm, Phip, Sup (m1,g), Sup (m2,-g)), Scalar4 1, G_HGSFSF (SU,m1,m2,g));
         ((Hp, Phim, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, 
             G_HGSFSF (SD,m1,m2,g));
         ((Hm, Phip, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, 
             G_HGSFSF (SD,m1,m2,g));
         ((Phip, Phim, Sup (m1,g), Sup (m2,-g)), Scalar4 1, 
             G_GGSFSF (SU,m1,m2,g));
         ((Phip, Phim, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, 
             G_GGSFSF (SD,m1,m2,g));
         ((Phi0, Phi0, Sup (m1,g), Sup (m2,-g)), Scalar4 1, 
             G_G0G0SFSF (SU,m1,m2,g));
         ((Phi0, Phi0, Sdown (m1,g), Sdown (m2,-g)), Scalar4 1, 
             G_G0G0SFSF (SD,m1,m2,g)) ]
    let higgs_gold_sfermion g = List.flatten (Product.list2 
             (higgs_gold_sfermion' g) [M1;M2] [M1;M2])

    let higgs_gold_squark' g h m1 m2 = 
      [ ((Hp, Phi0, Sup (m1,-g), Sdown (m2,h)), Scalar4 1, 
            G_HGSUSD (false,m1,m2,g,h));
        ((Hm, Phi0, Sup (m1,g), Sdown (m2,-h)), Scalar4 1, 
            G_HGSUSD (true,m1,m2,g,h));
        ((H_Heavy, Phip, Sup (m1,-g), Sdown (m2,h)), Scalar4 1, 
            G_H2GSUSD (false,m1,m2,g,h));
        ((H_Heavy, Phim, Sup (m1,g), Sdown (m2,-h)), Scalar4 1, 
            G_H2GSUSD (true,m1,m2,g,h));
        ((H_Light, Phip, Sup (m1,-g), Sdown (m2,h)), Scalar4 1, 
            G_H1GSUSD (false,m1,m2,g,h));
        ((H_Light, Phim, Sup (m1,g), Sdown (m2,-h)), Scalar4 1, 
            G_H1GSUSD (true,m1,m2,g,h));
        ((A, Phip, Sup (m1,-g), Sdown (m2,h)), Scalar4 1, 
            G_AGSUSD (false,m1,m2,g,h));
        ((A, Phim, Sup (m1,g), Sdown (m2,-h)), Scalar4 1, 
            G_AGSUSD (true,m1,m2,g,h));
        ((Phi0, Phip, Sup (m1,-g), Sdown (m2,h)), Scalar4 1, 
            G_GGSUSD (false,m1,m2,g,h));
        ((Phi0, Phim, Sup (m1,g), Sdown (m2,-h)), Scalar4 1, 
            G_GGSUSD (true,m1,m2,g,h)) ]
    let higgs_gold_squark g h = List.flatten (Product.list2 (higgs_gold_squark' 
                g h) [M1;M2] [M1;M2])

    let sneutrino4' (g,h) =
      [ ((Sneutrino g, Sneutrino h, Sneutrino (-g), Sneutrino (-h)), Scalar4 1, 
            G_SN4 (g,h))]
    let sneutrino4 = ThoList.flatmap sneutrino4' 
                   [(1,1);(1,2);(1,3);(2,2);(2,3);(3,3)]

    let sneu2_slep2_1' g h m1 m2 = 
       ((Sneutrino (-g), Sneutrino g, Slepton (m1,-h), Slepton (m2,h)), Scalar4 1, 
              G_SN2SL2_1 (m1,m2,g,h))
    let sneu2_slep2_2' (g,h) m1 m2 =
        ((Sneutrino g, Sneutrino (-h), Slepton (m1,-g), Slepton (m2,h)), Scalar4 1,
              G_SN2SL2_2 (m1,m2,g,h)) 
    let sneu2_slep2_1 g h = Product.list2 (sneu2_slep2_1' g h) [M1;M2] [M1;M2]
    let sneu2_slep2_2 (g,h) = Product.list2 (sneu2_slep2_2' (g,h)) [M1;M2] [M1;M2]

(* The 4-slepton-vertices have the following structure: The sleptons come up in 
   pairs of a positive and a negative slepton of the same generation; there is 
   no vertex with e.g. two negative selectrons and two positive smuons, that of 
   course would be a contradiction to the conservation of the separate slepton 
   numbers of each generation which is not implemented in the MSSM. Because there 
   is no CKM-mixing for the sleptons (in case of massless neutrinos) we maximally 
   have two different generations of sleptons in a 4-slepton-vertex. *)

    let slepton4_1gen' g (m1,m2,m3,m4) =
      [ ((Slepton (m1,-g), Slepton (m2,g), Slepton (m3,-g), Slepton (m4,g)), 
            Scalar4 1, G_SL4 (m1,m2,m3,m4,g)) ]
    let slepton4_1gen g = ThoList.flatmap (slepton4_1gen' g) [(M1,M1,M1,M1);
        (M1,M1,M1,M2); (M1,M1,M2,M1); (M1,M1,M2,M2); (M1,M2,M1,M2); (M1,M2,M2,M1);
        (M1,M2,M2,M2); (M2,M1,M2,M2); (M2,M2,M2,M2) ]      
    let slepton4_2gen' (g,h) (m1,m2) (m3,m4) =
       ((Slepton (m1,-g), Slepton (m2,g), Slepton (m3,-h), Slepton (m4,h)), 
           Scalar4 1, G_SL4_2 (m1,m2,m3,m4,g,h)) 
    let slepton4_2gen (g,h) = 
      Product.list2 (slepton4_2gen' (g,h)) [(M1,M1);(M1,M2);(M2,M1);(M2,M2)] 
        [(M1,M1);(M1,M2);(M2,M1);(M2,M2)]
                                                        
    let sneu2_squark2' g h m1 m2 =
       [ ((Sneutrino (-g), Sneutrino g, Sup (m1,-h), Sup (m2,h)), Scalar4 1, 
               G_SN2SQ2 (SU,m1,m2,g,h)); 
         ((Sneutrino (-g), Sneutrino g, Sdown (m1,-h), Sdown (m2,h)), Scalar4 1, 
               G_SN2SQ2 (SD,m1,m2,g,h)) ]
    let sneu2_squark2 g h = List.flatten (Product.list2 (sneu2_squark2' g h) 
                                            [M1;M2] [M1;M2]) 

    let slepton2_squark2'' g h m1 m2 m3 m4 = 
       [ ((Slepton (m1,-g), Slepton (m2,g), Sup (m3,-h), Sup (m4,h)), Scalar4 1, 
               G_SL2SQ2 (SU,m1,m2,m3,m4,g,h));
         ((Slepton (m1,-g), Slepton (m2,g), Sdown (m3,-h), Sdown (m4,h)), 
               Scalar4 1, G_SL2SQ2 (SD,m1,m2,m3,m4,g,h)) ]
    let slepton2_squark2' g h m1 m2 =
      List.flatten (Product.list2 (slepton2_squark2'' g h m1 m2) [M1;M2] [M1;M2])
    let slepton2_squark2 g h = 
      List.flatten (Product.list2 (slepton2_squark2' g h) [M1;M2] [M1;M2])

    let slep_sneu_squark2'' g1 g2 g3 m1 m2 m3 = 
      [ ((Sup (m1,-g1), Sdown (m2,g2), Slepton (m3,-g3), Sneutrino g3), 
              Scalar4 1, G_SUSDSNSL (false,m1,m2,m3,g1,g2,g3));
        ((Sup (m1,g1), Sdown (m2,-g2), Slepton (m3,g3), Sneutrino (-g3)), 
              Scalar4 1, G_SUSDSNSL (true,m1,m2,m3,g1,g2,g3)) ]
    let slep_sneu_squark2' g1 g2 g3 m1 =
      List.flatten (Product.list2 (slep_sneu_squark2'' g1 g2 g3 m1) 
                      [M1;M2] [M1;M2])
    let slep_sneu_squark2 g1 g2 =
      List.flatten (Product.list2 (slep_sneu_squark2' g1 g2) [1;2;3] [M1;M2])
        
(* There are three kinds of 4-squark-vertices: Four up-Squarks, four down-squarks          
   or two up- and two down-squarks. *)

    let sup4_1gen' g (m1,m2,m3,m4) =
      [ ((Sup (m1,-g), Sup (m2,g), Sup (m3,-g), Sup (m4,g)), Scalar4 1, 
              G_SU4 (m1,m2,m3,m4,g)) ]
    let sup4_1gen g = ThoList.flatmap (sup4_1gen' g) [(M1,M1,M1,M1);
        (M1,M1,M1,M2); (M1,M1,M2,M1); (M1,M1,M2,M2); (M1,M2,M1,M2); (M1,M2,M2,M1);
        (M1,M2,M2,M2); (M2,M1,M2,M2); (M2,M2,M2,M2) ]      
    let sup4_2gen' (g,h) (m1,m2) (m3,m4) =
       ((Sup (m1,-g), Sup (m2,g), Sup (m3,-h), Sup (m4,h)), Scalar4 1,
              G_SU4_2 (m1,m2,m3,m4,g,h)) 
    let sup4_2gen (g,h) = 
      Product.list2 (sup4_2gen' (g,h)) [(M1,M1);(M1,M2);(M2,M1);(M2,M2)] 
        [(M1,M1);(M1,M2);(M2,M1);(M2,M2)]    

    let sdown4_1gen' g (m1,m2,m3,m4) =
      [ ((Sdown (m1,-g), Sdown (m2,g), Sdown (m3,-g), Sdown (m4,g)), Scalar4 1, 
              G_SD4 (m1,m2,m3,m4,g)) ]
    let sdown4_1gen g = ThoList.flatmap (sdown4_1gen' g) [(M1,M1,M1,M1);
        (M1,M1,M1,M2); (M1,M1,M2,M1); (M1,M1,M2,M2); (M1,M2,M1,M2); (M1,M2,M2,M1);
        (M1,M2,M2,M2); (M2,M1,M2,M2); (M2,M2,M2,M2) ]      
    let sdown4_2gen' (g,h) (m1,m2) (m3,m4) =
       ((Sdown (m1,-g), Sdown (m2,g), Sdown (m3,-h), Sdown (m4,h)), Scalar4 1,
              G_SD4_2 (m1,m2,m3,m4,g,h)) 
    let sdown4_2gen (g,h) = 
      Product.list2 (sdown4_2gen' (g,h)) [(M1,M1);(M1,M2);(M2,M1);(M2,M2)] 
        [(M1,M1);(M1,M2);(M2,M1);(M2,M2)]  

    let sup2_sdown2_3 g1 g2 g3 g4 m1 m2 m3 m4 =
        ((Sup (m1,-g1), Sup (m2,g2), Sdown (m3,-g3), Sdown 
            (m4,g4)), Scalar4 1, G_SU2SD2 (m1,m2,m3,m4,g1,g2,g3,g4))
    let sup2_sdown2_2 g1 g2 g3 g4 m1 m2 =
      Product.list2 (sup2_sdown2_3 g1 g2 g3 g4 m1 m2) [M1;M2] [M1;M2]
    let sup2_sdown2_1 g1 g2 g3 g4 =
      List.flatten (Product.list2 (sup2_sdown2_2 g1 g2 g3 g4) [M1;M2] [M1;M2])
    let sup2_sdown2 g1 g2 =
      List.flatten (Product.list2 (sup2_sdown2_1 g1 g2) [1;2;3] [1;2;3])

    let quartic_grav_gauge g m = 
      [ ((Grino, Slepton (m, -g), Ga, L g), GBBG (1, Gravbar, SLRV, Psi), G_Gr4A_Sl (g,m));
        ((L (-g), Slepton (m, g), Ga, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4A_Slc (g,m));
        ((Grino, Sup (m, -g), Ga, U g), GBBG (1, Gravbar, SLRV, Psi), G_Gr4A_Su (g,m));
        ((U (-g), Sup (m, g), Ga, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4A_Suc (g,m));
        ((Grino, Sdown (m, -g), Ga, D g), GBBG (1, Gravbar, SLRV, Psi), G_Gr4A_Sd (g,m));
        ((D (-g), Sdown (m, g), Ga, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4A_Sdc (g,m));
        ((Grino, Slepton (m, -g), Z, L g), GBBG (1, Gravbar, SLRV, Psi), G_Gr4Z_Sl (g,m));
        ((L (-g), Slepton (m, g), Z, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4Z_Slc (g,m));
        ((Grino, Sup (m, -g), Z, U g), GBBG (1, Gravbar, SLRV, Psi), G_Gr4Z_Su (g,m));
        ((U (-g), Sup (m, g), Z, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4Z_Suc (g,m));
        ((Grino, Sdown (m, -g), Z, D g), GBBG (1, Gravbar, SLRV, Psi), G_Gr4Z_Sd (g,m));
        ((D (-g), Sdown (m, g), Z, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4Z_Sdc (g,m));
        ((Grino, Sup (m, -g), Gl, U g), GBBG (1, Gravbar, SLRV, Psi), G_Gr4Gl_Su (g,m));
        ((U (-g), Sup (m, g), Gl, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4Gl_Suc (g,m));
        ((Grino, Sdown (m, -g), Gl, D g), GBBG (1, Gravbar, SLRV, Psi), G_Gr4Gl_Sd (g,m));
        ((D (-g), Sdown (m, g), Gl, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4Gl_Sdc (g,m));
        ((Grino, Slepton (m, -g), Wm, N g), GBBG (1, Gravbar, SLV, Psi), G_Gr4W_Sl (g,m));
        ((N (-g), Slepton (m, g), Wp, Grino), GBBG (1, Psibar, SLV, Grav), G_Gr4Z_Slc (g,m));
        ((Grino, Sup (m, -g), Wp, D g), GBBG (1, Gravbar, SLV, Psi), G_Gr4W_Su (g,m));
        ((D (-g), Sup (m, g), Wm, Grino), GBBG (1, Psibar, SLV, Grav), G_Gr4W_Suc (g,m));
        ((Grino, Sdown (m, -g), Wm, U g), GBBG (1, Gravbar, SLV, Psi), G_Gr4W_Sd (g,m));
        ((U (-g), Sdown (m, g), Wp, Grino), GBBG (1, Psibar, SLV, Grav), G_Gr4W_Sdc (g,m)) ]

    let quartic_grav_sneutrino g = 
      [ ((Grino, Sneutrino (-g), Z, N g), GBBG (1, Gravbar, SLV, Psi), G_Gr4Z_Sn);
        ((N (-g), Sneutrino g, Z, Grino), GBBG (1, Psibar, SLV, Grav), G_Gr4Z_Snc);
        ((Grino, Sneutrino (-g), Wp, L g), GBBG (1, Gravbar, SLV, Psi), G_Gr4W_Sn);
        ((L (-g), Sneutrino g, Wm, Grino), GBBG (1, Psibar, SLV, Grav), G_Gr4W_Snc) ]

    let quartic_grav_neu n = 
      [ ((Grino, Wp, Wm, Neutralino n), GBBG (1, Gravbar, V2LR, Chi), G_Gr4_Neu n);
        ((Grino, H_Light, Z, Neutralino n), GBBG (1, Gravbar, SLRV, Chi), G_Gr4_Z_H1 n);
        ((Grino, H_Heavy, Z, Neutralino n), GBBG (1, Gravbar, SLRV, Chi), G_Gr4_Z_H2 n);
        ((Grino, A, Z, Neutralino n), GBBG (1, Gravbar, SLRV, Chi), G_Gr4_Z_H3 n);
        ((Grino, Hm, Wp, Neutralino n), GBBG (1, Gravbar, SLRV, Chi), G_Gr4_W_H n);
        ((Grino, Hp, Wm, Neutralino n), GBBG (1, Gravbar, SLRV, Chi), G_Gr4_W_Hc n) ]

    let quartic_grav_char c = 
      let cc = conj_char c in
      [ ((Grino, Wm, Ga, Chargino c), GBBG (1, Gravbar, V2LR, Psi), G_Gr4_A_Ch c);
        ((Grino, Wm, Z, Chargino c), GBBG (1, Gravbar, V2LR, Psi), G_Gr4_Z_Ch c);
        ((Chargino cc, Wp, Ga, Grino), GBBG ((-1), Psibar, V2LR, Grav), G_Gr4_A_Ch cc);
        ((Chargino cc, Wp, Z, Grino), GBBG ((-1), Psibar, V2LR, Grav), G_Gr4_Z_Ch cc);
        ((Grino, Hm, Ga, Chargino c), GBBG (1, Gravbar, SLRV, Psi), G_Gr4_H_A c); 
        ((Chargino cc, Hp, Ga, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4_H_A cc);
        ((Grino, Hm, Z, Chargino c), GBBG (1, Gravbar, SLRV, Psi), G_Gr4_H_Z c); 
        ((Chargino cc, Hp, Z, Grino), GBBG (1, Psibar, SLRV, Grav), G_Gr4_H_Z cc)]
       
    let quartic_gravitino =
      [ ((Grino, Gl, Gl, Gluino), GBBG (1, Gravbar, V2, Chi), G_GravGl)] @
      ThoList.flatmap quartic_grav_neu [N1; N2; N3; N4] @
      ThoList.flatmap quartic_grav_char [C1; C2] @
      List.flatten (Product.list2 quartic_grav_gauge [1; 2; 3] [M1; M2]) @
      ThoList.flatmap quartic_grav_sneutrino [1; 2; 3]
                
    let vertices3'' = 
      if Flags.ckm_present then
        (ThoList.flatmap electromagnetic_currents_3 [1;2;3] @
         ThoList.flatmap electromagnetic_currents_2 [C1;C2] @
         List.flatten (Product.list2 
                electromagnetic_sfermion_currents [1;2;3] [M1;M2]) @ 
         ThoList.flatmap neutral_currents [1;2;3] @
         ThoList.flatmap neutral_sfermion_currents [1;2;3] @  
         ThoList.flatmap charged_currents [1;2;3] @
         List.flatten (Product.list2 charged_slepton_currents [1;2;3] 
                         [M1;M2]) @ 
         List.flatten (Product.list2 charged_quark_currents [1;2;3] 
                         [1;2;3]) @
         List.flatten (Product.list2 charged_squark_currents [1;2;3] 
                         [1;2;3]) @ 
         ThoList.flatmap yukawa_higgs_quark [(1,3);(2,3);(3,3);(3,1);(3,2)] @
         yukawa_higgs 3 @ yukawa_n @ 
         ThoList.flatmap yukawa_c [C1;C2] @ 
         ThoList.flatmap yukawa_cq [C1;C2] @ 
         List.flatten (Product.list2 charged_chargino_currents [N1;N2;N3;N4] 
                         [C1;C2]) @ triple_gauge @ 
         ThoList.flatmap neutral_Z_1 [(N1,N2);(N1,N3);(N1,N4);(N2,N3);(N2,N4);
                                    (N3,N4)] @
         ThoList.flatmap neutral_Z_2 [N1;N2;N3;N4] @ neutral_A @
         Product.list2 charged_Z [C1;C2] [C1;C2] @ 
         gauge_higgs @ higgs @ yukawa_higgs_2 @ 
         List.flatten (Product.list2 higgs_charg_neutr [N1;N2;N3;N4] [C1;C2]) @ 
         higgs_neutr @ higgs_sneutrino @ higgs_sfermion @ 
         higgs_squark @ yukawa_v  @
         ThoList.flatmap col_currents [1;2;3] @ 
         List.flatten (Product.list2 col_sfermion_currents [1;2;3] [M1;M2]))
      else
        (ThoList.flatmap electromagnetic_currents_3 [1;2;3] @
         ThoList.flatmap electromagnetic_currents_2 [C1;C2] @
         List.flatten (Product.list2 
                electromagnetic_sfermion_currents [1;2;3] [M1;M2]) @ 
         ThoList.flatmap neutral_currents [1;2;3] @
         ThoList.flatmap neutral_sfermion_currents [1;2;3] @  
         ThoList.flatmap charged_currents [1;2;3] @
         List.flatten (Product.list2 charged_slepton_currents [1;2;3] 
                         [M1;M2]) @ 
         charged_quark_currents 1 1 @
         charged_quark_currents 2 2 @
         charged_quark_currents 3 3 @
         charged_squark_currents 1 1 @
         charged_squark_currents 2 2 @
         charged_squark_currents 3 3 @ 
         ThoList.flatmap yukawa_higgs_quark [(3,3)] @
         yukawa_higgs 3 @ yukawa_n @ 
         ThoList.flatmap yukawa_c [C1;C2] @ 
         ThoList.flatmap yukawa_cq [C1;C2] @ 
         List.flatten (Product.list2 charged_chargino_currents [N1;N2;N3;N4] 
                         [C1;C2]) @ triple_gauge @ 
         ThoList.flatmap neutral_Z_1 [(N1,N2);(N1,N3);(N1,N4);(N2,N3);(N2,N4);
                                    (N3,N4)] @
         ThoList.flatmap neutral_Z_2 [N1;N2;N3;N4] @ neutral_A @
         Product.list2 charged_Z [C1;C2] [C1;C2] @ 
         gauge_higgs @ higgs @ yukawa_higgs_2 @ 
         List.flatten (Product.list2 higgs_charg_neutr [N1;N2;N3;N4] [C1;C2]) @ 
         higgs_neutr @ higgs_sneutrino @ higgs_sfermion @ 
         higgs_squark @ yukawa_v  @
         ThoList.flatmap col_currents [1;2;3] @ 
         List.flatten (Product.list2 col_sfermion_currents [1;2;3] [M1;M2]))

    let vertices3' = 
      if Flags.gravitino then (vertices3'' @ triple_gravitino) 
      else vertices3''
    let vertices3 = 
      if Flags.include_goldstone then
        (vertices3' @ yukawa_goldstone 3 @
         gauge_higgs_gold @ higgs_gold @ yukawa_goldstone_2 @ 
         (if Flags.ckm_present then
           List.flatten (Product.list2 yukawa_goldstone_quark [1;2;3] 
                           [1;2;3]) @
           List.flatten (Product.list2 goldstone_charg_neutr [N1;N2;N3;N4] 
                           [C1;C2])
         else
           yukawa_goldstone_quark 1 1 @
           yukawa_goldstone_quark 2 2 @
           yukawa_goldstone_quark 3 3) @
         goldstone_neutr @ goldstone_sfermion @ goldstone_squark)
      else vertices3'
        
    let vertices4''' = 
      (quartic_gauge @ higgs4 @ gauge_higgs4 @ 
       ThoList.flatmap gauge_sfermion4 [1;2;3] @
       gauge_squark4 @ gluon_w_squark @
       List.flatten (Product.list2 gluon2_squark2  [1;2;3] [M1;M2]) @
       ThoList.flatmap gluon_gauge_squark [1;2;3])
    let vertices4'' =
      if Flags.gravitino then (vertices4''' @ quartic_gravitino)          
        else vertices4'''
   let vertices4' =
    if Flags.include_four then
       (vertices4'' @  
        ThoList.flatmap higgs_sfermion4 [1;2;3] @
        ThoList.flatmap higgs_sneutrino4 [1;2;3] @
        List.flatten (Product.list2 higgs_squark4 [1;2;3] [1;2;3]) @ 
        sneutrino4 @ 
        List.flatten (Product.list2 sneu2_slep2_1 [1;2;3] [1;2;3]) @
        ThoList.flatmap sneu2_slep2_2 [(1,2);(1,3);(2,3);(2,1);(3,1);(3,2)] @ 
        ThoList.flatmap slepton4_1gen [1;2;3] @
        ThoList.flatmap slepton4_2gen [(1,2);(1,3);(2,3)] @
        List.flatten (Product.list2 sneu2_squark2 [1;2;3] [1;2;3]) @
        List.flatten (Product.list2 slepton2_squark2 [1;2;3] [1;2;3]) @
        List.flatten (Product.list2 slep_sneu_squark2 [1;2;3] [1;2;3]) @
        ThoList.flatmap sup4_1gen [1;2;3] @
        ThoList.flatmap sup4_2gen [(1,2);(1,3);(2,3)] @
        ThoList.flatmap sdown4_1gen [1;2;3] @
        ThoList.flatmap sdown4_2gen [(1,2);(1,3);(2,3)] @
        List.flatten (Product.list2 sup2_sdown2 [1;2;3] [1;2;3]))
        else
      vertices4''
   let vertices4 =  
     if Flags.include_goldstone then 
       (vertices4' @ higgs_gold4 @ gauge_higgs_gold4 @ goldstone4 @  
        ThoList.flatmap higgs_gold_sneutrino [1;2;3] @ 
        ThoList.flatmap higgs_gold_sfermion [1;2;3] @  
        List.flatten (Product.list2 higgs_gold_squark [1;2;3] [1;2;3]))
     else 
       vertices4' 

    let vertices () = (vertices3, vertices4, [])

    let table = F.of_vertices (vertices ())
    let fuse2 = F.fuse2 table
    let fuse3 = F.fuse3 table
    let fuse = F.fuse table                                    
    let max_degree () = 4

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
      | "H" -> H_Heavy | "h" -> H_Light | "A0" -> A 
      | "H+" -> Hp | "H-" -> Hm
      | "phi0" -> Phi0 | "phi+" -> Phip | "phim" -> Phim
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
      | "ch1+" -> Chargino C1 | "ch2+" -> Chargino C2
      | "ch1-" -> Chargino C1c | "ch2-" -> Chargino C2c
      | "GR" -> Grino
      | _ -> invalid_arg "Modellib_MSSM.MSSM.flavor_of_string"

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
      | D 1 -> "d" | D (-1) -> "dbar"
      | D 2 -> "s" | D (-2) -> "sbar"
      | D 3 -> "b" | D (-3) -> "bbar"
      | L _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_string: invalid lepton"
      | N _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_string: invalid neutrino"
      | U _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_string: invalid up type quark"
      | D _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_string: invalid down type quark"
      | Gl -> "gl" | Gluino -> "sgl"
      | Ga -> "A" | Z -> "Z" | Wp -> "W+" | Wm -> "W-"
      | Phip -> "phi+" | Phim -> "phi-" | Phi0 -> "phi0" 
      | H_Heavy -> "H" | H_Light -> "h" | A -> "A0" 
      | Hp -> "H+" | Hm -> "H-"
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
      | Neutralino N1 -> "neu1"
      | Neutralino N2 -> "neu2"
      | Neutralino N3 -> "neu3"
      | Neutralino N4 -> "neu4"
      | Slepton _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_string: invalid slepton"
      | Sneutrino _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_string: invalid sneutrino"
      | Sup _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_string: invalid up type squark"
      | Sdown _ -> invalid_arg 
            "Modellib_MSSM.MSSM.flavor_to_string: invalid down type squark"
      | Chargino C1 -> "ch1+" | Chargino C1c -> "ch1-"
      | Chargino C2 -> "ch2+" | Chargino C2c -> "ch2-"
      | Grino -> "GR"

    let flavor_symbol = function
      | L g when g > 0 -> "l" ^ string_of_int g
      | L g -> "l" ^ string_of_int (abs g) ^ "b"  
      | N g when g > 0 -> "n" ^ string_of_int g
      | N g -> "n" ^ string_of_int (abs g) ^ "b"      
      | U g when g > 0 -> "u" ^ string_of_int g 
      | U g -> "u" ^ string_of_int (abs g) ^ "b"  
      | D g when g > 0 ->  "d" ^ string_of_int g 
      | D g -> "d" ^ string_of_int (abs g) ^ "b"    
      | Gl -> "gl" | Ga -> "a" | Z -> "z"
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
      | Gluino -> "sgl" | Phip -> "pp" | Phim -> "pm" | Phi0 -> "p0"
      | H_Heavy -> "h0h" | H_Light -> "h0l" | A -> "a0"
      | Hp -> "hp" | Hm -> "hm" | Grino -> "gv"
                
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
            "Modellib_MSSM.MSSM.flavor_to_TeX: invalid lepton"
      | N _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_TeX: invalid neutrino"
      | U _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_TeX: invalid up type quark"
      | D _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_TeX: invalid down type quark"
      | Gl -> "g" | Gluino -> "\\widetilde{g}"
      | Ga -> "\\gamma" | Z -> "Z" | Wp -> "W^+" | Wm -> "W^-"
      | Phip -> "\\phi^+" | Phim -> "\\phi^-" | Phi0 -> "\\phi^0" 
      | H_Heavy -> "H^0" | H_Light -> "h^0" | A -> "A^0" 
      | Hp -> "H^+" | Hm -> "H^-"
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
      | Slepton _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_TeX: invalid slepton"
      | Sneutrino _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_TeX: invalid sneutrino"
      | Sup _ -> invalid_arg
            "Modellib_MSSM.MSSM.flavor_to_TeX: invalid up type squark"
      | Sdown _ -> invalid_arg 
            "Modellib_MSSM.MSSM.flavor_to_TeX: invalid down type squark"
      | Chargino C1  -> "\\widetilde{\\chi}_1^+" 
      | Chargino C1c -> "\\widetilde{\\chi}_1^-"
      | Chargino C2  -> "\\widetilde{\\chi}_2^+" 
      | Chargino C2c -> "\\widetilde{\\chi}_2^-"
      | Grino -> "\\widetilde{G}"

    let pdg = function
      | L g when g > 0 -> 9 + 2*g
      | L g -> - 9 + 2*g
      | N g when g > 0 -> 10 + 2*g
      | N g -> - 10 + 2*g
      | U g when g > 0 -> 2*g
      | U g -> 2*g
      | D g when g > 0 -> - 1 + 2*g
      | D g -> 1 + 2*g
      | Gl -> 21 | Ga -> 22 | Z -> 23
      | Wp -> 24 | Wm -> (-24)
      | H_Light -> 25 | H_Heavy -> 35 | A -> 36
      | Hp -> 37 | Hm -> (-37)
      | Phip | Phim -> 27 | Phi0 -> 26              
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
      | Grino -> 1000039
      | Chargino C1 -> 1000024 | Chargino C1c -> (-1000024)
      | Chargino C2 -> 1000037 | Chargino C2c -> (-1000037)
      | Neutralino N1 -> 1000022 | Neutralino N2 -> 1000023
      | Neutralino N3 -> 1000025 | Neutralino N4 -> 1000035


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
         $\widetilde{\psi}_\mu$     & gravitino          &     39 \\\hline\hline
         $\widetilde{d}_L$          & down-squark 1      &     41 \\\hline 
         $\widetilde{u}_L$          & up-squark 1        &     42 \\\hline
         $\widetilde{s}_L$          & strange-squark 1   &     43 \\\hline
         $\widetilde{c}_L$          & charm-squark 1     &     44 \\\hline
         $\widetilde{b}_L$          & bottom-squark 1    &     45 \\\hline
         $\widetilde{t}_L$          & top-squark 1       &     46 \\\hline
         $\widetilde{d}_R$          & down-squark 2      &     47 \\\hline 
         $\widetilde{u}_R$          & up-squark 2        &     48 \\\hline
         $\widetilde{s}_R$          & strange-squark 2   &     49 \\\hline
         $\widetilde{c}_R$          & charm-squark 2     &     50 \\\hline
         $\widetilde{b}_R$          & bottom-squark 2    &     51 \\\hline
         $\widetilde{t}_R$          & top-squark 2       &     52 \\\hline\hline
         $\widetilde{e}_L$          & selectron 1        &     53 \\\hline
         $\widetilde{\nu}_{e,L}$    & electron-sneutrino &     54 \\\hline
         $\widetilde{\mu}_L$        & smuon 1            &     55 \\\hline
         $\widetilde{\nu}_{\mu,L}$  & muon-sneutrino     &     56 \\\hline
         $\widetilde{\tau}_L$       & stau 1             &     57 \\\hline
         $\widetilde{\nu}_{\tau,L}$ & tau-sneutrino      &     58 \\\hline
         $\widetilde{e}_R$          & selectron 2        &     59 \\\hline
         $\widetilde{\mu}_R$        & smuon 2            &     61 \\\hline
         $\widetilde{\tau}_R$       & stau 2             &     63 \\\hline\hline
         $\widetilde{g}$            & gluino             &     64 \\\hline
         $\widetilde{\chi}^0_1$     & neutralino 1       &     65 \\\hline
         $\widetilde{\chi}^0_2$     & neutralino 2       &     66 \\\hline
         $\widetilde{\chi}^0_3$     & neutralino 3       &     67 \\\hline
         $\widetilde{\chi}^0_4$     & neutralino 4       &     68 \\\hline
         $\widetilde{\chi}^+_1$     & chargino 1         &     69 \\\hline
         $\widetilde{\chi}^+_2$     & chargino 2         &     70 \\\hline\hline
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
      | Gl -> 21 | Ga -> 22 | Z -> 23
      | Wp -> 24 | Wm -> (-24)
      | H_Light -> 25 | H_Heavy -> 35 | A -> 36
      | Hp -> 37 | Hm -> (-37)
      | Phip | Phim -> 27 | Phi0 -> 26              
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
      | Grino -> 39
      | Gluino -> 64
      | Chargino C1 -> 69 | Chargino C1c -> (-69)
      | Chargino C2 -> 70 | Chargino C2c -> (-70)
      | Neutralino N1 -> 65 | Neutralino N2 -> 66
      | Neutralino N3 -> 67 | Neutralino N4 -> 68 

    let mass_symbol f =
      "mass(" ^ string_of_int (abs (pdg_mw f)) ^ ")"  

    let width_symbol f =
      "width(" ^ string_of_int (abs (pdg_mw f)) ^ ")"  

    let conj_symbol = function
      | false, str -> str
      | true, str -> str ^ "_c"

    let constant_symbol = function
      | Unit -> "unit" | Pi -> "PI"
      | Alpha_QED -> "alpha" | E -> "e" | G -> "g" | Vev -> "vev"
      | Sin2thw -> "sin2thw" | Eidelta -> "eidelta" | Mu -> "mu" | G_Z -> "gz"
      | Sin a -> "sin" ^ string_of_angle a | Cos a -> "cos" ^ string_of_angle a 
      | Sin2am2b -> "sin2am2b" | Cos2am2b -> "cos2am2b" | Sinamb -> "sinamb"
      | Sinapb -> "sinapb" | Cosamb -> "cosamb" | Cosapb -> "cosapb" 
      | Cos4be -> "cos4be" | Sin4be -> "sin4be" | Sin4al -> "sin4al" 
      | Sin2al -> "sin2al" | Cos2al -> "cos2al" | Sin2be -> "sin2be" 
      | Cos2be -> "cos2be" | Tana -> "tana" | Tanb -> "tanb"
      | Q_lepton -> "qlep" | Q_up -> "qup" | Q_down -> "qdwn"
      | Q_charg -> "qchar"
      | V_CKM (g1,g2) -> "vckm_" ^ string_of_int g1 ^ string_of_int g2
      | M_SF (f,g,m1,m2) -> "mix_" ^ string_of_sff f ^ string_of_int g
          ^ string_of_sfm m1 ^ string_of_sfm m2
      | AL g -> "al_" ^ string_of_int g 
      | AD g -> "ad_" ^ string_of_int g 
      | AU g -> "au_" ^ string_of_int g 
      | A_0 (n1,n2) -> "a0_" ^ string_of_neu n1 ^ string_of_neu n2
      | A_P (c1,c2) -> "ap_" ^ string_of_char c1 ^ string_of_char c2
      | V_0 (n1,n2) -> "v0_" ^ string_of_neu n1 ^ string_of_neu n2
      | V_P (c1,c2) -> "vp_" ^ string_of_char c1 ^ string_of_char c2
      | M_N (n1,n2) -> "mn_" ^ string_of_neu n1 ^ string_of_neu n2
      | M_U (c1,c2) -> "mu_" ^ string_of_char c1 ^ string_of_char c2
      | M_V (c1,c2) -> "mv_" ^ string_of_char c1 ^ string_of_char c2
      | L_NC (n,c) -> "lnc_" ^ string_of_neu n ^ string_of_char c
      | R_NC (n,c) -> "rnc_" ^ string_of_neu n ^ string_of_char c
      | L_CN (c,n) -> "lcn_" ^ string_of_char c ^ string_of_neu n
      | R_CN (c,n) -> "rcn_" ^ string_of_char c ^ string_of_neu n
      | L_NCH (n,c) -> "lnch_" ^ string_of_neu n ^ string_of_char c
      | R_NCH (n,c) -> "rnch_" ^ string_of_neu n ^ string_of_char c
      | L_CNG (c,n) -> "lcng_" ^ string_of_char c ^ string_of_neu n 
      | R_CNG (c,n) -> "rcng_" ^ string_of_char c ^ string_of_neu n 
      | S_NNA (n1,n2) -> "snna_" ^ string_of_neu n1 ^ string_of_neu n2
      | P_NNA (n1,n2) -> "pnna_" ^ string_of_neu n1 ^ string_of_neu n2
      | S_NNG (n1,n2) -> "snng_" ^ string_of_neu n1 ^ string_of_neu n2
      | P_NNG (n1,n2) -> "pnng_" ^ string_of_neu n1 ^ string_of_neu n2
      | S_NNH1 (n1,n2) -> "snnh1_" ^ string_of_neu n1 ^ string_of_neu n2
      | P_NNH1 (n1,n2) -> "pnnh1_" ^ string_of_neu n1 ^ string_of_neu n2
      | S_NNH2 (n1,n2) -> "snnh2_" ^ string_of_neu n1 ^ string_of_neu n2
      | P_NNH2 (n1,n2) -> "pnnh2_" ^ string_of_neu n1 ^ string_of_neu n2
      | G_NC_lepton -> "gnclep" | G_NC_neutrino -> "gncneu" 
      | G_NC_up -> "gncup" | G_NC_down -> "gncdwn"
      | G_CC -> "gcc"
      | G_CCQ (vc,g1,g2) -> conj_symbol (vc, "gccq_" ^ string_of_int g1 ^ "_" 
          ^ string_of_int g2)
      | I_Q_W -> "iqw" | I_G_ZWW -> "igzww" 
      | G_WWWW -> "gw4" | G_ZZWW -> "gzzww"
      | G_PZWW -> "gpzww" | G_PPWW -> "gppww"   
      | G_GH 1 -> "ghaw" 
      | G_GH 2 -> "gh1az" | G_GH 3 -> "gh2az"
      | G_GH 4 -> "gh1ww" | G_GH 5 -> "gh2ww" 
      | G_GH 6 -> "ghh1w" | G_GH 7 -> "ghh2w" 
      | G_GH 8 -> "gh1zz" | G_GH 9 -> "gh2zz" 
      | G_GH 10 -> "ghhz" | G_GH 11 -> "ghhp"            
      | G_GH _ ->  failwith "this G_GH coupling is not available"
      | G_GLGLH -> "gglglh" | G_GLGLHH -> "gglglhh" 
      | G_GLGLA -> "gglgla" | G_PPH -> "gpph"
      | G_PPHH -> "gpphh" | G_PPA -> "gppa"
      | G_GHGo n -> "g_hgh(" ^ string_of_int n ^ ")"  
      | G_GH4 1 -> "gaazz" | G_GH4 2 -> "gh1h1zz" | G_GH4 3 -> "gh2h2zz"
      | G_GH4 4 -> "ghphmzz" | G_GH4 5 -> "ghphmpp" | G_GH4 6 -> "ghphmpz"
      | G_GH4 7 -> "ghh1wz" | G_GH4 8 -> "ghh2wz"
      | G_GH4 9 -> "ghh1wp" | G_GH4 10 -> "ghh2wp" 
      | G_GH4 11 -> "gaaww" | G_GH4 12 -> "gh1h1ww" | G_GH4 13 -> "gh2h2ww" 
      | G_GH4 14 -> "ghhww" | G_GH4 15 -> "ghawz" | G_GH4 16 -> "ghawp" 
      | G_GH4 _ ->  failwith "this G_GH4 coupling is not available"
      | G_CICIH1 (n1,n2) -> "gcicih1_" ^ string_of_neu n1 ^ "_" 
          ^ string_of_neu n2
      | G_CICIH2 (n1,n2) -> "gcicih2_" ^ string_of_neu n1 ^ "_" 
          ^ string_of_neu n2 
      | G_CICIA (n1,n2) -> "gcicia_" ^ string_of_neu n1 ^ "_" 
          ^ string_of_neu n2
      | G_CICIG (n1,n2) -> "gcicig_" ^ string_of_neu n1 ^ "_" 
          ^ string_of_neu n2 
      | G_H3 n -> "gh3_" ^ string_of_int n
      | G_H4 n -> "gh4_" ^ string_of_int n 
      | G_HGo3 n -> "ghg3_" ^ string_of_int n 
      | G_HGo4 n -> "ghg4_" ^ string_of_int n 
      | G_GG4 n -> "ggg4_" ^ string_of_int n
      | G_strong -> "gs" | G_SS -> "gs**2" 
      | Gs -> "gs"
      | I_G_S -> "igs"
      | G_S_Sqrt -> "gssq"
      | G_NWC (n,c) -> "gnwc_" ^ string_of_neu n ^ "_" ^ string_of_char c 
      | G_CWN (c,n) -> "gcwn_" ^ string_of_char c ^ "_" ^ string_of_neu n 
      | G_CH1C (c1,c2) -> "gch1c_" ^ string_of_char c1 ^ "_" ^ string_of_char c2
      | G_CH2C (c1,c2) -> "gch2c_" ^ string_of_char c1 ^ "_" ^ string_of_char c2
      | G_CAC (c1,c2) -> "gcac_" ^ string_of_char c1 ^ "_" ^ string_of_char c2
      | G_CGC (c1,c2) -> "gcgc_" ^ string_of_char c1 ^ "_" ^ string_of_char c2
      | G_YUK (i,g) -> "g_yuk" ^ string_of_int i ^ "_" ^ string_of_int g
      | G_NZN (n1,n2) -> "gnzn_" ^ string_of_neu n1 ^ "_" ^ string_of_neu n2
      | G_NNA -> "gnna"
      | G_CZC (c1,c2) -> "gczc_" ^ string_of_char c1 ^ "_" ^ string_of_char 
          c2 
      | G_YUK_1 (n,m) -> "g_yuk1_" ^ string_of_int n ^ "_" ^ string_of_int m 
      | G_YUK_2 (n,m) -> "g_yuk2_" ^ string_of_int n ^ "_" ^ string_of_int m 
      | G_YUK_3 (n,m) -> "g_yuk3_" ^ string_of_int n ^ "_" ^ string_of_int m 
      | G_YUK_4 (n,m) -> "g_yuk4_" ^ string_of_int n ^ "_" ^ string_of_int m 
      | G_YUK_C (vc,g,c,sf,m) -> conj_symbol (vc, "g_yuk_ch" ^ string_of_char c  
          ^ "_" ^ string_of_sff sf ^ string_of_sfm m ^ "_" ^ string_of_int g ) 
      | G_YUK_N (vc,g,n,sf,m) -> conj_symbol (vc, "g_yuk_n" ^ string_of_neu n
          ^ "_" ^ string_of_sff sf ^ string_of_sfm m ^ "_" ^ string_of_int g )
      | G_YUK_G (vc,g,sf,m) -> conj_symbol (vc, "g_yuk_g" ^ string_of_sff sf 
          ^ string_of_sfm m ^ "_" ^ string_of_int g)
      | G_YUK_Q (vc,g1,g2,c,sf,m) -> conj_symbol (vc, "g_yuk_ch" ^ string_of_char c  
          ^ "_" ^ string_of_sff sf ^ string_of_sfm m ^ "_" ^ string_of_int g1 
          ^ "_" ^ string_of_int g2) 
      | G_NHC (n,c) -> "g_nhc_" ^ string_of_neu n ^ "_" ^ string_of_char c 
      | G_CHN (c,n) -> "g_chn_" ^ string_of_neu n ^ "_" ^ string_of_char c
      | G_NGC (n,c) -> "g_ngc_" ^ string_of_neu n ^ string_of_char c
      | G_CGN (c,n) -> "g_cgn_" ^ string_of_char c ^ string_of_neu n
      | SUM_1 -> "sum1"
      | G_SLSNW (vc,g,m) -> conj_symbol (vc, "gsl" ^ string_of_sfm m ^ "_" 
          ^ string_of_int g ^ "snw") 
      | G_ZSF (f,g,m1,m2) -> "g" ^ string_of_sff f ^ string_of_sfm m1 ^ "z" 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int g 
      | G_WWSFSF (f,g,m1,m2) -> "gww" ^ string_of_sff f ^ string_of_sfm m1   
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int g 
      | G_WPSLSN (vc,g,m) -> conj_symbol (vc, "gpwsl" ^ string_of_sfm m 
          ^ "sn_" ^ string_of_int g)
      | G_WZSLSN (vc,g,m) -> conj_symbol (vc, "gwzsl" ^ string_of_sfm m 
          ^ "sn_" ^ string_of_int g)
      | G_H1SFSF (f,g,m1,m2) -> "gh1" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int g    
      | G_H2SFSF (f,g,m1,m2) -> "gh2" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int g  
      | G_ASFSF (f,g,m1,m2) -> "ga" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int g
      | G_HSNSL (vc,g,m) -> conj_symbol (vc, "ghsnsl" ^ string_of_sfm m ^ "_" 
          ^ string_of_int g)
      | G_GoSFSF (f,g,m1,m2) -> "ggo" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int g 
      | G_GoSNSL (vc,g,m) -> conj_symbol (vc, "ggosnsl" ^ string_of_sfm m ^ "_" 
          ^ string_of_int g) 
      | G_HSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "ghsu" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 ^ "_" 
          ^ string_of_int g2)
      | G_GSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "ggsu" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 ^ "_" 
          ^ string_of_int g2)
      | G_WPSUSD (vc,m1,m2,n,m) -> conj_symbol (vc, "gpwpsu" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int n ^ "_" 
          ^ string_of_int m)
      | G_WZSUSD (vc,m1,m2,n,m) -> conj_symbol (vc, "gzwpsu" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int n ^ "_" 
          ^ string_of_int m)
      | G_SWS (vc,g1,g2,m1,m2) -> conj_symbol (vc, "gs" ^ string_of_sfm m1 ^ "ws" 
          ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 ^ "_" ^ string_of_int g2)
      | G_GlGlSQSQ -> "gglglsqsq" 
      | G_PPSFSF f -> "gpp" ^ string_of_sff f ^ string_of_sff f 
      | G_ZZSFSF (f,g,m1,m2) -> "gzz" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int g 
      | G_ZPSFSF (f,g,m1,m2) -> "gzp" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int g 
      | G_GlPSQSQ -> "gglpsqsq" 
      | G_GlZSFSF (f,g,m1,m2) -> "ggl" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int g 
      | G_GlWSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "gglwsu" 
          ^ string_of_sfm m1 ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 
          ^ "_" ^ string_of_int g2)
      | G_GHGo4 1 -> "gzzg0g0" | G_GHGo4 2 -> "gzzgpgm" 
      | G_GHGo4 3 -> "gppgpgm" | G_GHGo4 4 -> "gzpgpgm" 
      | G_GHGo4 5 -> "gwwgpgm" | G_GHGo4 6 -> "gwwg0g0"
      | G_GHGo4 7 -> "gwzg0g" | G_GHGo4 8 -> "gwzg0g" 
      | G_GHGo4 9 -> "gwzh1g" | G_GHGo4 10 -> "gwzh2g"
      | G_GHGo4 11 -> "gwph1g" | G_GHGo4 12 -> "gwph2g"   
      | G_GHGo4 _ -> failwith "Coupling G_GHGo4 is not available"
      | G_HSF31 (h,g,m1,m2,f1,f2) -> "g_" ^ string_of_higgs h ^ 
          string_of_int g ^ string_of_sfm m1 ^ string_of_sfm m2 ^
          string_of_sff f1 ^ string_of_sff f2
      | G_HSF32 (h,g1,g2,m1,m2,f1,f2) -> "g_" ^ string_of_higgs h ^ 
          string_of_int g1 ^ "_" ^ string_of_int g2 ^ string_of_sfm m1 ^ 
          string_of_sfm m2 ^ string_of_sff f1 ^ string_of_sff f2
      | G_HSF41 (h,g,m1,m2,f1,f2) -> "g_" ^ string_of_higgs h ^ 
          string_of_int g ^ string_of_sfm m1 ^ string_of_sfm m2 ^
          string_of_sff f1 ^ string_of_sff f2
      | G_HSF42 (h,g1,g2,m1,m2,f1,f2) -> "g_" ^ string_of_higgs h ^ 
          string_of_int g1 ^ "_" ^ string_of_int g2 ^ string_of_sfm m1 ^ 
          string_of_sfm m2 ^ string_of_sff f1 ^ string_of_sff f2
      | G_H1H1SFSF (f,m1,m2,n) -> "gh1h1" ^ string_of_sff f ^ string_of_sfm 
          m1 ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int n      
      | G_H1H2SFSF (f,m1,m2,n) -> "gh1h2" ^ string_of_sff f ^ string_of_sfm 
          m1 ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int n   
      | G_H2H2SFSF (f,m1,m2,n) -> "gh2h2" ^ string_of_sff f ^ string_of_sfm 
          m1 ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int n  
      | G_HHSFSF (f,m1,m2,n) -> "ghh" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int n  
      | G_AASFSF (f,m1,m2,n) -> "gaa" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int n 
      | G_HH1SUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "ghh1su" 
          ^ string_of_sfm m1 ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 
          ^ "_" ^ string_of_int g2)
      | G_HH2SUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "ghh2su" 
          ^ string_of_sfm m1 ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 
          ^ "_" ^ string_of_int g2)
      | G_HASUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "ghasu" 
          ^ string_of_sfm m1 ^ "sd" ^ string_of_sfm m2 ^ "_" 
          ^ string_of_int g1 ^ "_" ^ string_of_int g2 ^ "_c")
      | G_HH1SLSN (vc,m,g) -> conj_symbol (vc, "ghh1sl" ^ string_of_sfm m 
                                           ^ "sn_" ^ string_of_int g)
      | G_HH2SLSN (vc,m,g) -> conj_symbol (vc, "ghh2sl" ^ string_of_sfm m 
                                           ^ "sn_" ^ string_of_int g)
      | G_HASLSN (vc,m,g) -> conj_symbol (vc, "ghasl" ^ string_of_sfm m  
                                           ^ "sn_" ^ string_of_int g)
      | G_AG0SFSF (f,m1,m2,n) -> "gag0" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int n  
      | G_HGSFSF (f,m1,m2,n) -> "ghg" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m1 ^ "_" ^ string_of_int n  
      | G_GGSFSF (f,m1,m2,n) -> "ggg" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int n  
      | G_G0G0SFSF (f,m1,m2,n) -> "gg0g0" ^ string_of_sff f ^ string_of_sfm m1 
          ^ string_of_sff f ^ string_of_sfm m2 ^ "_" ^ string_of_int n  
      | G_HGSNSL (vc,m,n) -> conj_symbol (vc, "ghgsnsl" ^ string_of_sfm m ^ "_" 
          ^ string_of_int n)
      | G_H1GSNSL (vc,m,n) -> conj_symbol (vc, "gh1gsnsl" ^ string_of_sfm m ^ "_" 
          ^ string_of_int n)
      | G_H2GSNSL (vc,m,n) -> conj_symbol (vc, "gh2gsnsl" ^ string_of_sfm m ^ "_" 
          ^ string_of_int n) 
      | G_AGSNSL (vc,m,n) -> conj_symbol (vc, "gagsnsl" ^ string_of_sfm m ^ "_" 
          ^ string_of_int n) 
      | G_GGSNSL (vc,m,n) -> conj_symbol (vc, "gggsnsl" ^ string_of_sfm m ^ "_" 
          ^ string_of_int n) 
      | G_HGSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "gghpsu" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 ^ "_" 
          ^ string_of_int g2)
      | G_H1GSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "gh1gpsu" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 ^ "_" 
          ^ string_of_int g2)
      | G_H2GSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "gh2gpsu" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 ^ "_" 
          ^ string_of_int g2)
      | G_AGSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "gagpsu" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 ^ "_" 
          ^ string_of_int g2)
      | G_GGSUSD (vc,m1,m2,g1,g2) -> conj_symbol (vc, "gggpsu" ^ string_of_sfm m1 
          ^ "sd" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 ^ "_" 
          ^ string_of_int g2)
      | G_SN4 (g1,g2) -> "gsn4_" ^ string_of_int g1 ^ "_" ^ string_of_int g2 
      | G_SN2SL2_1 (m1,m2,g1,g2) -> "gsl_" ^ string_of_int g1 ^ "_sl_" 
          ^ string_of_int g1 ^ "_sl" ^ string_of_sfm m1 ^ "_" ^ string_of_int g2 
          ^ "_sl" ^ string_of_sfm m2 ^ "_" ^ string_of_int g2 
      | G_SN2SL2_2 (m1,m2,g1,g2) -> "gsl_" ^ string_of_int g1 ^ "_sl_" 
          ^ string_of_int g2 ^ "_sl" ^ string_of_sfm m1 ^ "_" ^ string_of_int g1 
          ^ "_sl" ^ string_of_sfm m2 ^ "_" ^ string_of_int g2 ^ "_mix"
      | G_SF4 (f1,f2,m1,m2,m3,m4,g1,g2) -> "gsf" ^ string_of_sff f1 ^ 
          string_of_sff f2 ^ string_of_sfm m1 ^ string_of_sfm m2 ^ 
          string_of_sfm m3 ^ string_of_sfm m4 ^ string_of_int g1 ^
          string_of_int g2 
      | G_SF4_3 (f1,f2,m1,m2,m3,m4,g1,g2,g3) -> "gsf" ^ string_of_sff f1 ^ 
          string_of_sff f2 ^ string_of_sfm m1 ^ string_of_sfm m2 ^ 
          string_of_sfm m3 ^ string_of_sfm m4 ^ string_of_int g1 ^
          string_of_int g2 ^ "_" ^ string_of_int g3 
      | G_SF4_4 (f1,f2,m1,m2,m3,m4,g1,g2,g3,g4) -> "gsf" ^ string_of_sff f1 ^ 
          string_of_sff f2 ^ string_of_sfm m1 ^ string_of_sfm m2 ^ 
          string_of_sfm m3 ^ string_of_sfm m4 ^ string_of_int g1 ^ "_" ^ 
          string_of_int g2 ^ string_of_int g3 ^ "_" ^ string_of_int g4
      | G_SL4 (m1,m2,m3,m4,g) -> "gsl" ^ string_of_sfm m1 ^ "_" 
          ^ "sl" ^ string_of_sfm m2 ^ "_" ^ "sl" ^ string_of_sfm m3 ^ "_" 
          ^ "sl" ^ string_of_sfm m4 ^ "_" ^ string_of_int g 
      | G_SL4_2 (m1,m2,m3,m4,g1,g2) -> "gsl" ^ string_of_sfm m1 ^ "_" 
          ^ "sl" ^ string_of_sfm m2 ^ "_" ^ "sl" ^ string_of_sfm m3 ^ "_" 
          ^ "sl" ^ string_of_sfm m4 ^ "_" ^ string_of_int g1 ^ "_" ^
          string_of_int g2
      | G_SN2SQ2 (f,m1,m2,g1,g2) -> "gsn_" ^ string_of_int g1 ^ "_sn_" 
          ^ string_of_int g1 ^ "_" ^ string_of_sff f ^ string_of_sfm m1 ^ "_" 
          ^ string_of_int g2 ^ "_" ^ string_of_sff f ^ string_of_sfm m2 ^ "_" 
          ^ string_of_int g2
      | G_SL2SQ2 (f,m1,m2,m3,m4,g1,g2) -> "gsl" ^ string_of_sfm m1 ^ "_" 
          ^ string_of_int g1 ^ "_sl" ^ string_of_sfm m2 ^ "_" ^ string_of_int g1 
          ^ "_" ^ string_of_sff f ^ string_of_sfm m3 ^ "_" ^ string_of_int g2 
          ^ "_" ^ string_of_sff f ^ string_of_sfm m4 ^ "_" ^ string_of_int g2 
      | G_SUSDSNSL (vc,m1,m2,m3,g1,g2,g3) -> conj_symbol (vc, "gsl" 
          ^ string_of_sfm m3 ^ "_" ^ string_of_int g3 ^ "_sn_" ^ string_of_int g3 
          ^ "_su" ^ string_of_sfm m1 ^ "_" ^ string_of_int g1 ^ "_sd" 
          ^ string_of_sfm m2 ^ "_" ^ string_of_int g2) 
      | G_SU4 (m1,m2,m3,m4,g) -> "gsu" ^ string_of_sfm m1 ^ "_" 
          ^ "_su" ^ string_of_sfm m2 ^ "_" ^ "_su" ^ string_of_sfm m3 ^ "_" ^ 
          "_su" ^ string_of_sfm m4 ^ "_" ^ string_of_int g 
      | G_SU4_2 (m1,m2,m3,m4,g1,g2) -> "gsu" ^ string_of_sfm m1 ^ "_" 
          ^ "_su" ^ string_of_sfm m2 ^ "_" ^ "_su" ^ string_of_sfm m3 ^ "_" ^ 
          "_su" ^ string_of_sfm m4 ^ "_" ^ string_of_int g1 ^ "_" ^ 
          string_of_int g2
      | G_SD4 (m1,m2,m3,m4,g) -> "gsd" ^ string_of_sfm m1 ^ "_" 
          ^ "_sd" ^ string_of_sfm m2 ^ "_" ^ "_sd" ^ string_of_sfm m3 ^ "_" 
          ^ "_sd" ^ string_of_sfm m4 ^ "_" ^ string_of_int g 
      | G_SD4_2 (m1,m2,m3,m4,g1,g2) -> "gsd" ^ string_of_sfm m1 ^ "_" 
          ^ "_sd" ^ string_of_sfm m2 ^ "_" ^ "_sd" ^ string_of_sfm m3 ^ "_" 
          ^ "_sd" ^ string_of_sfm m4 ^ "_" ^ string_of_int g1 ^ "_" ^ 
          string_of_int g2
      | G_SU2SD2 (m1,m2,m3,m4,g1,g2,g3,g4) -> "gsu" ^ string_of_sfm m1 
          ^ "_" ^ string_of_int g1 ^ "_su" ^ string_of_sfm m2 ^ "_" 
          ^ string_of_int g2 ^ "_sd" ^ string_of_sfm m3 ^ "_" ^ string_of_int g3 
          ^ "_sd" ^ string_of_sfm m4 ^ "_" ^ string_of_int g4   
      | M f -> "mass" ^ flavor_symbol f
      | W f -> "width" ^ flavor_symbol f 
      | G_Grav -> "ggrav" | G_Gr_Ch C1 -> "ggrch1" | G_Gr_Ch C2 -> "ggrch2"
      | G_Gr_Ch C1c -> "ggrch1c" | G_Gr_Ch C2c  -> "ggrch2c"
      | G_Gr_Z_Neu n -> "ggrzneu" ^ string_of_neu n 
      | G_Gr_A_Neu n -> "ggraneu" ^ string_of_neu n 
      | G_Gr4_Neu n -> "ggr4neu" ^ string_of_neu n
      | G_Gr4_A_Ch C1 -> "ggr4ach1" | G_Gr4_A_Ch C2 -> "ggr4ach2" 
      | G_Gr4_A_Ch C1c -> "ggr4ach1c" | G_Gr4_A_Ch C2c -> "ggr4ach2c" 
      | G_Gr4_Z_Ch C1 -> "ggr4zch1" | G_Gr4_Z_Ch C2 -> "ggr4zch2"
      | G_Gr4_Z_Ch C1c -> "ggr4zch1c" | G_Gr4_Z_Ch C2c -> "ggr4zch2c"
      | G_Grav_N -> "ggravn" 
      | G_GravGl -> "gs * ggrav"
      | G_Grav_L  (g,m) -> "ggravl" ^ string_of_int g ^ string_of_sfm m 
      | G_Grav_Lc (g,m) -> "ggravl" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Grav_U  (g,m) -> "ggravu" ^ string_of_int g ^ string_of_sfm m
      | G_Grav_Uc (g,m) -> "ggravu" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Grav_D  (g,m) -> "ggravd" ^ string_of_int g ^ string_of_sfm m
      | G_Grav_Dc (g,m) -> "ggravd" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr_H_Ch C1 -> "ggrhch1" | G_Gr_H_Ch C2 -> "ggrhch2" 
      | G_Gr_H_Ch C1c -> "ggrhch1c" | G_Gr_H_Ch C2c -> "ggrhch2c" 
      | G_Gr_H1_Neu n -> "ggrh1neu" ^ string_of_neu n
      | G_Gr_H2_Neu n -> "ggrh2neu" ^ string_of_neu n
      | G_Gr_H3_Neu n -> "ggrh3neu" ^ string_of_neu n
      | G_Gr4A_Sl (g,m) -> "ggr4asl" ^ string_of_int g ^ string_of_sfm m
      | G_Gr4A_Slc (g,m) -> "ggr4asl" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4A_Su (g,m) -> "ggr4asu" ^ string_of_int g ^ string_of_sfm m
      | G_Gr4A_Suc (g,m) -> "ggr4asu" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4A_Sd (g,m) -> "ggr4asd" ^ string_of_int g ^ string_of_sfm m
      | G_Gr4A_Sdc (g,m) -> "ggr4asd" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4Z_Sn -> "ggr4zsn" | G_Gr4Z_Snc -> "ggr4zsnc"
      | G_Gr4Z_Sl (g,m) -> "ggr4zsl" ^ string_of_int g ^ string_of_sfm m
      | G_Gr4Z_Slc (g,m) -> "ggr4zsl" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4Z_Su (g,m) -> "ggr4zsu" ^ string_of_int g ^ string_of_sfm m
      | G_Gr4Z_Suc (g,m) -> "ggr4zsu" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4Z_Sd (g,m) -> "ggr4zsd" ^ string_of_int g ^ string_of_sfm m
      | G_Gr4Z_Sdc (g,m) -> "ggr4zsd" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4W_Sl (g,m) -> "ggr4wsl" ^ string_of_int g ^ string_of_sfm m
      | G_Gr4W_Slc (g,m) -> "ggr4wsl" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4W_Su (g,m) -> "ggr4wsu" ^ string_of_int g ^ string_of_sfm m
      | G_Gr4W_Suc (g,m) -> "ggr4wsu" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4W_Sd (g,m) -> "ggr4wsd" ^ string_of_int g ^ string_of_sfm m
      | G_Gr4W_Sdc (g,m) -> "ggr4wsd" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4Gl_Su (g,m) -> "ggr4glsu" ^ string_of_int g ^ string_of_sfm m 
      | G_Gr4Gl_Suc (g,m) -> "ggr4glsu" ^ string_of_int g ^ string_of_sfm m ^ "c" 
      | G_Gr4Gl_Sd (g,m) -> "ggr4glsd" ^ string_of_int g ^ string_of_sfm m 
      | G_Gr4Gl_Sdc (g,m) -> "ggr4glsd" ^ string_of_int g ^ string_of_sfm m ^ "c"
      | G_Gr4_Z_H1 n -> "ggr4zh1_" ^ string_of_neu n
      | G_Gr4_Z_H2 n -> "ggr4zh2_" ^ string_of_neu n
      | G_Gr4_Z_H3 n -> "ggr4zh3_" ^ string_of_neu n
      | G_Gr4_W_H n -> "ggr4wh_" ^ string_of_neu n
      | G_Gr4_W_Hc n -> "ggr4whc_" ^ string_of_neu n
      | G_Gr4_H_A C1 -> "ggr4ha1" | G_Gr4_H_A C2 -> "ggr4ha2" 
      | G_Gr4_H_A C1c -> "ggr4ha1c" | G_Gr4_H_A C2c -> "ggr4ha2c" 
      | G_Gr4_H_Z C1 -> "ggr4hz1" | G_Gr4_H_Z C2 -> "ggr4hz2" 
      | G_Gr4_H_Z C1c -> "ggr4hz1c" | G_Gr4_H_Z C2c -> "ggr4hz2c" 
      | G_Gr4W_Sn -> "ggr4wsn"
      | G_Gr4W_Snc -> "ggr4wsnc"
      
  end

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
