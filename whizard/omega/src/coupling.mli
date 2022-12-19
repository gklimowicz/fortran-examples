(* coupling.mli --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       Christian Speckner <cnspeckn@googlemail.com>
       Marco Sekulla <marco.sekulla@kit.edu>
       So Young Shim <soyoung.shim@desy.de> (only parts of this file)

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

(* The enumeration types used for communication from [Models]
   to [Targets].  On the physics side, the modules in [Models]
   must implement the Feynman rules according to the conventions
   set up here.  On the numerics side, the modules in [Targets]
   must handle all cases according to the same conventions.  *)

(* \thocwmodulesection{Propagators}
   The Lorentz representation of the particle. NB:  O'Mega
   treats all lines as \emph{outgoing} and particles are therefore
   transforming as [ConjSpinor] and antiparticles as [Spinor]. *)
type lorentz =
  | Scalar
  | Spinor (* $\psi$ *)
  | ConjSpinor (* $\bar\psi$ *)
  | Majorana (* $\chi$ *) 
  | Maj_Ghost (* SUSY ghosts *)
  | Vector
(*i  | Ward_Vector  i*)
  | Massive_Vector
  | Vectorspinor (* supersymmetric currents and gravitinos *)
  | Tensor_1
  | Tensor_2 (* massive gravitons (large extra dimensions) *)
  | BRS of lorentz

type lorentz3 = lorentz * lorentz * lorentz
type lorentz4 = lorentz * lorentz * lorentz * lorentz
type lorentzn = lorentz list

type fermion_lines = (int * int) list

(* \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{2.2}
       \begin{tabular}{|r|l|l|}\hline
                    & only Dirac fermions & incl.~Majorana fermions \\\hline
           [Prop_Scalar]
         & \multicolumn{2}{ l |}{%
             $\displaystyle\phi(p)\leftarrow
              \frac{\ii}{p^2-m^2+\ii m\Gamma}\phi(p)$} \\\hline
           [Prop_Spinor]
         & $\displaystyle\psi(p)\leftarrow
            \frac{\ii(-\fmslash{p}+m)}{p^2-m^2+\ii m\Gamma}\psi(p)$
         & $\displaystyle\psi(p)\leftarrow
            \frac{\ii(-\fmslash{p}+m)}{p^2-m^2+\ii m\Gamma}\psi(p)$ \\\hline
           [Prop_ConjSpinor]
         & $\displaystyle\bar\psi(p)\leftarrow
            \bar\psi(p)\frac{\ii(\fmslash{p}+m)}{p^2-m^2+\ii m\Gamma}$
         & $\displaystyle\psi(p)\leftarrow
            \frac{\ii(-\fmslash{p}+m)}{p^2-m^2+\ii m\Gamma}\psi(p)$ \\\hline
           [Prop_Majorana]
         & \multicolumn{1}{ c |}{N/A}
         & $\displaystyle\chi(p)\leftarrow
            \frac{\ii(-\fmslash{p}+m)}{p^2-m^2+\ii m\Gamma}\chi(p)$ \\\hline
           [Prop_Unitarity]
         & \multicolumn{2}{ l |}{%
             $\displaystyle\epsilon_\mu(p)\leftarrow
              \frac{\ii}{p^2-m^2+\ii m\Gamma}
              \left(-g_{\mu\nu}+\frac{p_\mu p_\nu}{m^2}\right)\epsilon^\nu(p)$} \\\hline
           [Prop_Feynman]
         & \multicolumn{2}{ l |}{%
             $\displaystyle\epsilon^\nu(p)\leftarrow
              \frac{-\ii}{p^2-m^2+\ii m\Gamma}\epsilon^\nu(p)$} \\\hline
           [Prop_Gauge]
         & \multicolumn{2}{ l |}{%
             $\displaystyle\epsilon_\mu(p)\leftarrow
              \frac{\ii}{p^2}
              \left(-g_{\mu\nu}+(1-\xi)\frac{p_\mu p_\nu}{p^2}\right)\epsilon^\nu(p)$} \\\hline
           [Prop_Rxi]
         & \multicolumn{2}{ l |}{%
             $\displaystyle\epsilon_\mu(p)\leftarrow
              \frac{\ii}{p^2-m^2+\ii m\Gamma}
              \left(-g_{\mu\nu}+(1-\xi)\frac{p_\mu p_\nu}{p^2-\xi m^2}\right)
              \epsilon^\nu(p)$} \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:propagators} Propagators.  NB: The sign of the
        momenta in the spinor propagators comes about because O'Mega
        treats all momenta as \emph{outgoing} and the charge flow for
        [Spinor] is therefore opposite to the momentum, while the charge
        flow for [ConjSpinor] is parallel to the momentum.}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.5}
       \begin{tabular}{|r|l|}\hline
           [Aux_Scalar]
         & $\displaystyle\phi(p)\leftarrow\ii\phi(p)$ \\\hline
           [Aux_Spinor]
         & $\displaystyle\psi(p)\leftarrow\ii\psi(p)$ \\\hline
           [Aux_ConjSpinor]
         & $\displaystyle\bar\psi(p)\leftarrow\ii\bar\psi(p)$ \\\hline
           [Aux_Vector]
         & $\displaystyle\epsilon^\mu(p)\leftarrow\ii\epsilon^\mu(p)$ \\\hline
           [Aux_Tensor_1]
         & $\displaystyle T^{\mu\nu}(p)\leftarrow\ii T^{\mu\nu}(p)$ \\\hline
           [Only_Insertion]
         & \multicolumn{1}{ c |}{N/A} \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:aux-propagators} Auxiliary and non propagating fields}
   \end{table}
   If there were no vectors or auxiliary fields, we could deduce the propagator from
   the Lorentz representation.  While we're at it, we can introduce
   ``propagators'' for the contact interactions of auxiliary fields
   as well.  [Prop_Gauge] and [Prop_Feynman] are redundant as special
   cases of [Prop_Rxi].

   The special case [Only_Insertion] corresponds to operator insertions
   that do not correspond to a propagating field all.  These are used
   for checking Slavnov-Taylor identities
   \begin{equation}
      \partial_\mu\Braket{\text{out}|W^\mu(x)|\text{in}}
        = m_W\Braket{\text{out}|\phi(x)|\text{in}}
   \end{equation}
   of gauge theories in unitarity gauge where the Goldstone bosons are
   not propagating. Numerically, it would suffice to use a vanishing
   propagator, but then superflous fusions would be calculated in
   production code in which the Slavnov-Taylor identities are not tested. *)

type 'a propagator =
  | Prop_Scalar | Prop_Ghost 
  | Prop_Spinor | Prop_ConjSpinor | Prop_Majorana 
  | Prop_Unitarity | Prop_Feynman | Prop_Gauge of 'a | Prop_Rxi of 'a
  | Prop_Tensor_2 | Prop_Tensor_pure | Prop_Vector_pure
  | Prop_Vectorspinor
  | Prop_Col_Scalar | Prop_Col_Feynman | Prop_Col_Majorana 
  | Prop_Col_Unitarity 
  | Aux_Scalar | Aux_Vector | Aux_Tensor_1
  | Aux_Col_Scalar | Aux_Col_Vector | Aux_Col_Tensor_1
  | Aux_Spinor | Aux_ConjSpinor | Aux_Majorana
  | Only_Insertion
  | Prop_UFO of string

(* \begin{JR}
   We don't need different fermionic propagators as supposed by the variable
   names [Prop_Spinor], [Prop_ConjSpinor] or [Prop_Majorana]. The 
   propagator in all cases has to be multiplied on the left hand side of the
   spinor out of which a new one should be built. All momenta are treated as 
   \emph{outgoing}, so for the propagation of the different fermions the 
   following table arises, in which the momentum direction is always downwards
   and the arrows show whether the momentum and the fermion line,
   respectively are parallel or antiparallel to the direction of calculation:
   \begin{center}
   \begin{tabular}{|l|c|c|c|c|}\hline
     Fermion type & fermion arrow & mom. & calc. & sign \\\hline\hline 
     Dirac fermion &     $\uparrow$  &  $\uparrow~\downarrow$ & 
       $\uparrow~\uparrow$ & negative \\\hline
     Dirac antifermion & $\downarrow$ & $\downarrow~\downarrow$ & 
       $\uparrow~\downarrow$ & negative   \\\hline
     Majorana fermion  & - & $\uparrow~\downarrow$ & - & negative \\\hline
   \end{tabular}
   \end{center}
   So the sign of the momentum is always negative and no further distinction
   is needed. 
   \end{JR} *)

type width =
  | Vanishing
  | Constant
  | Timelike
  | Running
  | Fudged
  | Complex_Mass
  | Custom of string

(* \thocwmodulesection{Vertices}
   The combined $S-P$ and $V-A$ couplings (see
   tables~\ref{tab:dim4-fermions-SP}, \ref{tab:dim4-fermions-VA},
   \ref{tab:dim4-fermions-SPVA-maj} and~\ref{tab:dim4-fermions-SPVA-maj2})
   are redundant, of course, but they allow some targets to create
   more efficient numerical code.\footnote{An additional benefit
   is that the counting of Feynman diagrams is not upset by a splitting
   of the vectorial and axial pieces of gauge bosons.} Choosing VA2 over
	VA will cause the FORTRAN backend to pass the coupling as a whole array *)
type fermion = Psi | Chi | Grav 
type fermionbar = Psibar | Chibar | Gravbar 
type boson =
  | SP | SPM | S | P | SL | SR | SLR | VA | V | A | VL | VR | VLR | VLRM | VAM
  | TVA | TLR | TRL | TVAM | TLRM | TRLM
  | POT | MOM | MOM5 | MOML | MOMR | LMOM | RMOM | VMOM | VA2 | VA3 | VA3M
type boson2 = S2 | P2 | S2P | S2L | S2R | S2LR 
  | SV | PV | SLV | SRV | SLRV | V2 | V2LR


(* The integer is an additional coefficient that multiplies the respective
   coupling constant.  This allows to reduce the number of required coupling
   constants in manifestly symmetrc cases.  Most of times it will be equal
   unity, though. *)

(* The two vertex types [PBP] and [BBB] for the couplings of two fermions or 
   two antifermions ("clashing arrows") is unavoidable in supersymmetric 
   theories.
   \begin{dubious}
     \ldots{} tho doesn't like the names and has promised to find a better
     mnemonics! 
   \end{dubious} *) 

type 'a vertex3 =
  | FBF of int * fermionbar * boson * fermion
  | PBP of int * fermion * boson * fermion
  | BBB of int * fermionbar * boson * fermionbar
  | GBG of int * fermionbar * boson * fermion (* gravitino-boson-fermion *)
  | Gauge_Gauge_Gauge of int | Aux_Gauge_Gauge of int
  | I_Gauge_Gauge_Gauge of int
  | Scalar_Vector_Vector of int
  | Aux_Vector_Vector of int | Aux_Scalar_Vector of int
  | Scalar_Scalar_Scalar of int | Aux_Scalar_Scalar of int
  | Vector_Scalar_Scalar of int
  | Graviton_Scalar_Scalar of int
  | Graviton_Vector_Vector of int
  | Graviton_Spinor_Spinor of int
  | Dim4_Vector_Vector_Vector_T of int
  | Dim4_Vector_Vector_Vector_L of int
  | Dim4_Vector_Vector_Vector_T5 of int
  | Dim4_Vector_Vector_Vector_L5 of int
  | Dim6_Gauge_Gauge_Gauge of int
  | Dim6_Gauge_Gauge_Gauge_5 of int
  | Aux_DScalar_DScalar of int | Aux_Vector_DScalar of int
  | Dim5_Scalar_Gauge2 of int (* %
      $\frac12 \phi F_{1,\mu\nu} F_2^{\mu\nu} = - \frac12
      \phi (\ii \partial_{[\mu,} V_{1,\nu]})(\ii \partial^{[\mu,} V_2^{\nu]})$ *)
  | Dim5_Scalar_Gauge2_Skew of int 
      (* %
      $\frac14 \phi F_{1,\mu\nu} \tilde{F}_2^{\mu\nu} = -
      \phi (\ii \partial_\mu V_{1,\nu})(\ii \partial_\rho V_{2,\sigma})\epsilon^{\mu\nu\rho\sigma}$ *) 
  | Dim5_Scalar_Scalar2 of int (* %
    $\phi_1 \partial_\mu \phi_2 \partial^\mu \phi_3$ *)
  | Dim5_Scalar_Vector_Vector_T of int (* %
      $\phi(\ii\partial_\mu V_1^\nu)(\ii\partial_\nu V_2^\mu)$ *)
  | Dim5_Scalar_Vector_Vector_TU of int (* %
      $(\ii\partial_\nu\phi) (\ii\partial_\mu V_1^\nu) V_2^\mu$ *)
  | Dim5_Scalar_Vector_Vector_U of int (* %
      $(\ii\partial_\nu\phi) (\ii\partial_\mu V^\nu) V^\mu$ *)
  | Scalar_Vector_Vector_t of int (* %
      $ ( \partial_\mu V_\nu-\partial_\nu V_\mu )^2 $ *)
  | Dim6_Vector_Vector_Vector_T of int (* %
      $V_1^\mu ((\ii\partial_\nu V_2^\rho) %
       \ii\overleftrightarrow{\partial_\mu}(\ii\partial_\rho V_3^\nu))$ *)
  | Tensor_2_Vector_Vector of int (* %
      $T^{\mu\nu} (V_{1,\mu}V_{2,\nu} + V_{1,\nu}V_{2,\mu})$ *)
  | Tensor_2_Vector_Vector_1 of int (* % 
      $T^{\mu\nu} (V_{1,\mu}V_{2,\nu} + V_{1,\nu}V_{2,\mu} - g_{\mu,\nu}V_1^\rho V_{2,\rho} )$ *) 
  | Tensor_2_Vector_Vector_cf of int (* % 
      $T^{\mu\nu} ( %
        - \frac{c_f}{2} g_{\mu,\nu}V_1^\rho V_{2,\rho} )$ *) 
  | Tensor_2_Scalar_Scalar of int (* % 
      $T^{\mu\nu} (\partial_{\mu}\phi_1\partial_{\nu}\phi_2 + %
      		  \partial_{\nu}\phi_1\partial_{\mu}\phi_2 )$ *) 
  | Tensor_2_Scalar_Scalar_cf of int (* % 
      $T^{\mu\nu} ( - \frac{c_f}{2} g_{\mu,\nu} %
		  \partial_{\rho}\phi_1\partial_{\rho}\phi_2 )$ *) 
  | Tensor_2_Vector_Vector_t of int (* %
    $T^{\mu\nu} (V_{1,\mu}V_{2,\nu} + V_{1,\nu}V_{2,\mu} - g_{\mu,\nu}V_1^\rho V_{2,\rho} )$ *) 
  | Dim5_Tensor_2_Vector_Vector_1 of int (* %
      $T^{\alpha\beta} (V_1^\mu
         \ii\overleftrightarrow\partial_\alpha
         \ii\overleftrightarrow\partial_\beta V_{2,\mu}$ *)
  | Dim5_Tensor_2_Vector_Vector_2 of int
      (* %
      $T^{\alpha\beta}
       (   V_1^\mu \ii\overleftrightarrow\partial_\beta (\ii\partial_\mu V_{2,\alpha})
         + V_1^\mu \ii\overleftrightarrow\partial_\alpha (\ii\partial_\mu V_{2,\beta}))$ *)
  | Dim7_Tensor_2_Vector_Vector_T of int (* %
      $T^{\alpha\beta} ((\ii\partial^\mu V_1^\nu)
                          \ii\overleftrightarrow\partial_\alpha
                          \ii\overleftrightarrow\partial_\beta
					    (\ii\partial_\nu V_{2,\mu})) $ *)
  | Dim6_Scalar_Vector_Vector_D of int
    (* %
       $\ii \phi ( - (\partial^\mu \partial^\nu W^{-}_{\mu})W^{+}_{\nu} 
                     - (\partial^\mu \partial^\nu W^{+}_{\nu})W^{-}_{\mu}
					    \\ \mbox{} \qquad
                     + ( (\partial^\rho \partial_\rho W^{-}_{\mu})W^{+}_{\nu}
                        + (\partial^\rho \partial_\rho W^{+}_{\nu})W^{-}_{\mu})
					    g^{\mu\nu}) $ *)
  | Dim6_Scalar_Vector_Vector_DP of int
    (* %
       $\ii ( (\partial^\mu H)(\partial^\nu W^{-}_{\mu})W^{+}_{\nu}
               + (\partial^\nu H)(\partial^\mu W^{+}_{\nu})W^{-}_{\mu}
					     \\ \mbox{} \qquad
               - ((\partial^\rho H)(\partial_\rho W^{-}_{\mu})W^{+}_{\nu}
                   (\partial^\rho H)(\partial^\rho W^{+}_{\nu})W^{-}_{\mu})
					     g^{\mu\nu})  $*)
  | Dim6_HAZ_D of int   (* %
      $\ii ((\partial^\mu \partial^\nu A_{\mu})Z_{\nu} 
      + (\partial^\rho \partial_\rho A_{\mu})Z_{\nu}g^{\mu\nu} )$ *)
  | Dim6_HAZ_DP of int   (* %
      $\ii ((\partial^{\nu} A_{\mu})(\partial^{\mu} H)Z_{\nu} 
      - (\partial^{\rho} A_{\mu})(\partial_{\rho} H)Z_{\nu} g^{\mu\nu})$ *)
  | Dim6_AWW_DP of int   (* % 
      $\ii ((\partial^{\rho} A_{\mu}) W^{-}_{\nu} W^{+}_{\rho} g^{\mu\nu} 
      - (\partial^{\nu} A_{\mu}) W^{-}_{\nu} W^{+}_{\rho} g^{\mu\rho}) $ *)
  | Dim6_AWW_DW of int
    (*%
      $\ii [ (3(\partial^\rho A_{\mu})W^{-}_{\nu}W^{+}_{\rho} 
      - (\partial^\rho W^{-}_{\nu})A_{\mu}W^{+}_{\rho}
      + (\partial^\rho W^{+}_{\rho})A_{\mu} W^{-}_{\nu})g^{\mu\nu}
      \\ \mbox{} \qquad
      +(-3(\partial^\nu A_{\mu})W^{-}_{\nu}W^{+}_{\rho} 
      -  (\partial^\nu W^{-}_{\nu})A_{\mu}W^{+}_{\rho}
      + (\partial^\nu W^{+}_{\rho})A_{\mu}W^{-}_{\nu})g^{\mu\rho}
      \\ \mbox{} \qquad
      +(2(\partial^\mu W^{-}_{\nu})A_{\mu}W^{+}_{\rho} 
      - 2(\partial^\mu W^{+}_{\rho})A_{\mu}W^{-}_{\nu})g^{\nu\rho} ]$ 
    *)
  | Dim6_HHH of int   (*% 
      $\ii(-(\partial^{\mu}H_1)(\partial_{\mu}H_2)H_3 
      - (\partial^{\mu}H_1)H_2(\partial_{\mu}H_3) 
      - H_1(\partial^{\mu}H_2)(\partial_{\mu}H_3) )$ *)
  | Dim6_Gauge_Gauge_Gauge_i of int
    (*% 
      $\ii 
      (-(\partial^{\nu}V_{\mu})(\partial^{\rho}V_{\nu})(\partial^{\mu}V_{\rho})
      + (\partial^{\rho}V_{\mu})(\partial^{\mu}V_{\nu})(\partial^{\nu}V_{\rho})
      \\ \mbox{} \qquad
      + (-\partial^{\nu}V_{\rho} g^{\mu\rho} 
      + \partial^{\mu}V_{\rho} g^{\nu\rho})
      (\partial^{\sigma}V_{\mu})(\partial_{\sigma}V_{\nu})
      + (\partial^{\rho}V_{\nu} g^{\mu\nu} - \partial^{\mu}V_{\nu} g^{\nu\rho})
      (\partial^{\sigma}V_{\mu})(\partial_{\sigma}V_{\rho})
      \\ \mbox{} \qquad
      + (-\partial^{\rho}V_{\mu} g^{\mu\nu} + \partial^{\mu}V_{\mu} g^{\mu\rho})
      (\partial^{\sigma}V_{\nu})(\partial_{\sigma}V_{\rho}) )$ *)
  | Gauge_Gauge_Gauge_i of int
  | Dim6_GGG of int
  | Dim6_WWZ_DPWDW of int
    (* %
       $\ii( ((\partial^\rho V_{\mu})V_{\nu}V_{\rho} 
       - (\partial^{\rho}V_{\nu})V_{\mu}V_{\rho})g^{\mu\nu}
       - (\partial^{\nu}V_{\mu})V_{\nu}V_{\rho}g^{\mu\rho} 
       + (\partial^{\mu}V_{\nu})V_{\mu}V_{\rho})g^{\rho\nu} )$ *)
  | Dim6_WWZ_DW of int
    (* % 
       $\ii( ((\partial^\mu V_{\mu})V_{\nu}V_{\rho} 
       + V_{\mu}(\partial^\mu V_{\nu})V_{\rho})g^{\nu\rho}
       - ((\partial^\nu V_{\mu})V_{\nu}V_{\rho} 
       + V_{\mu}(\partial^\nu V_{\nu})V_{\rho})g^{\mu\rho})$ *)
  | Dim6_WWZ_D of int   (* %
      $\ii (  V_{\mu})V_{\nu}(\partial^{\nu}V_{\rho})g^{\mu\rho} 
      + V_{\mu}V_{\nu}(\partial^{\mu}V_{\rho})g^{\nu\rho})$ 
    *)
  | TensorVector_Vector_Vector of int
  | TensorVector_Vector_Vector_cf of int
  | TensorVector_Scalar_Scalar of int
  | TensorVector_Scalar_Scalar_cf of int
  | TensorScalar_Vector_Vector of int
  | TensorScalar_Vector_Vector_cf of int
  | TensorScalar_Scalar_Scalar of int
  | TensorScalar_Scalar_Scalar_cf of int

(* As long as we stick to renormalizable couplings, there are only
   three types of quartic couplings: [Scalar4], [Scalar2_Vector2]
   and [Vector4].  However, there are three inequivalent contractions
   for the latter and the general vertex will be a linear combination
   with integer coefficients:
   \begin{subequations}
   \begin{align}
     \ocwupperid{Scalar4}\,1 :&\;\;\;\;\;
       \phi_1 \phi_2 \phi_3 \phi_4 \\
     \ocwupperid{Scalar2\_Vector2}\,1 :&\;\;\;\;\;
       \phi_1^{\vphantom{\mu}} \phi_2^{\vphantom{\mu}}
        V_3^\mu V_{4,\mu}^{\vphantom{\mu}} \\
     \ocwupperid{Vector4}\,\lbrack 1, \ocwupperid{C\_12\_34} \rbrack :&\;\;\;\;\;
        V_1^\mu V_{2,\mu}^{\vphantom{\mu}}
        V_3^\nu V_{4,\nu}^{\vphantom{\mu}} \\
     \ocwupperid{Vector4}\,\lbrack 1, \ocwupperid{C\_13\_42} \rbrack :&\;\;\;\;\;
        V_1^\mu V_2^\nu
        V_{3,\mu}^{\vphantom{\mu}} V_{4,\nu}^{\vphantom{\mu}} \\
     \ocwupperid{Vector4}\,\lbrack 1, \ocwupperid{C\_14\_23} \rbrack :&\;\;\;\;\;
        V_1^\mu V_2^\nu
        V_{3,\nu}^{\vphantom{\mu}} V_{4,\mu}^{\vphantom{\mu}}
   \end{align}
   \end{subequations} *)

type contract4 = C_12_34 | C_13_42 | C_14_23

(*i\begin{dubious}
     CS objected to the polymorphic [type 'a vertex4], since it broke the
     implementation of some of his extensions.  Is there another way of
     getting coupling constants into [Vector4_K_Matrix], besides the brute
     force solution of declaring the possible coupling constants here?
     \textit{I'd like to put the blame on CS for two reasons: it's not clear
     that the brute force solution will actually work and everytime a new
     vertex that depends non-linearly on coupling contanst pops up, the
     problem will make another appearance.}
   \end{dubious}i*)

type 'a vertex4 =
  | Scalar4 of int 
  | Scalar2_Vector2 of int
  | Vector4 of (int * contract4) list
  | DScalar4 of (int * contract4) list
  | DScalar2_Vector2 of (int * contract4) list
  | Dim8_Scalar2_Vector2_1 of int
  | Dim8_Scalar2_Vector2_2 of int
  | Dim8_Scalar2_Vector2_m_0 of int
  | Dim8_Scalar2_Vector2_m_1 of int  
  | Dim8_Scalar2_Vector2_m_7 of int
  | Dim8_Scalar4 of int
  | Dim8_Vector4_t_0 of (int * contract4) list
  | Dim8_Vector4_t_1 of (int * contract4) list  
  | Dim8_Vector4_t_2 of (int * contract4) list
  | Dim8_Vector4_m_0 of (int * contract4) list  
  | Dim8_Vector4_m_1 of (int * contract4) list  
  | Dim8_Vector4_m_7 of (int * contract4) list
  | GBBG of int * fermionbar * boson2 * fermion

(* In some applications, we have to allow for contributions outside of
   perturbation theory.  The most prominent example is heavy gauge boson
   scattering at very high energies, where the perturbative expression
   violates unitarity. *)
 
(* One solution is the `$K$-matrix' ansatz.  Such unitarizations typically
   introduce effective propagators and/or vertices that violate crossing
   symmetry and vanish in the $t$-channel.  This can be taken care of in
   [Fusion] by filtering out vertices that have the wrong momenta.  *)

(* In this case the ordering of the fields in a vertex of the Feynman
   rules becomes significant. In particular, we assume that $(V_1,V_2,V_3,V_4)$
   implies
   \begin{equation}
     \parbox{25mm}{\fmfframe(2,3)(2,3){\begin{fmfgraph*}(20,20)
       \fmfleft{v1,v2}
       \fmfright{v4,v3}
       \fmflabel{$V_1$}{v1}
       \fmflabel{$V_2$}{v2}
       \fmflabel{$V_3$}{v3}
       \fmflabel{$V_4$}{v4}
       \fmf{plain}{v,v1}
       \fmf{plain}{v,v2}
       \fmf{plain}{v,v3}
       \fmf{plain}{v,v4}
       \fmfblob{.2w}{v}
     \end{fmfgraph*}}}
       \qquad\Longrightarrow\qquad
     \parbox{45mm}{\fmfframe(2,3)(2,3){\begin{fmfgraph*}(40,20)
       \fmfleft{v1,v2}
       \fmfright{v4,v3}
       \fmflabel{$V_1$}{v1}
       \fmflabel{$V_2$}{v2}
       \fmflabel{$V_3$}{v3}
       \fmflabel{$V_4$}{v4}
       \fmf{plain}{v1,v12,v2}
       \fmf{plain}{v3,v34,v4}
       \fmf{dots,label=$\Theta((p_1+p_2)^2)$,tension=0.7}{v12,v34}
       \fmfdot{v12,v34}
     \end{fmfgraph*}}}
   \end{equation}
   The list of pairs of parameters denotes the location and strengths
   of the poles in the $K$-matrix ansatz:
   \begin{equation}
     (c_1,a_1,c_2,a_2,\ldots,c_n,a_n) \Longrightarrow
       f(s) = \sum_{i=1}^{n} \frac{c_i}{s-a_i}
   \end{equation} *)
  | Vector4_K_Matrix_tho of int * ('a * 'a) list    
  | Vector4_K_Matrix_jr of int * (int * contract4) list
  | Vector4_K_Matrix_cf_t0 of int * (int * contract4) list 
  | Vector4_K_Matrix_cf_t1 of int * (int * contract4) list
  | Vector4_K_Matrix_cf_t2 of int * (int * contract4) list
  | Vector4_K_Matrix_cf_t_rsi of int * (int * contract4) list
  | Vector4_K_Matrix_cf_m0 of int * (int * contract4) list  
  | Vector4_K_Matrix_cf_m1 of int * (int * contract4) list   
  | Vector4_K_Matrix_cf_m7 of int * (int * contract4) list
  | DScalar2_Vector2_K_Matrix_ms of int * (int * contract4) list
  | DScalar2_Vector2_m_0_K_Matrix_cf of int * (int * contract4) list
  | DScalar2_Vector2_m_1_K_Matrix_cf of int * (int * contract4) list
  | DScalar2_Vector2_m_7_K_Matrix_cf of int * (int * contract4) list
  | DScalar4_K_Matrix_ms of int * (int * contract4) list
  | Dim6_H4_P2 of int
    (* %
       $\ii( -(\partial^{\mu}H_1)(\partial_{\mu}H_2) H_3 H_4 
       - (\partial^{\mu}H_1)H_2(\partial_{\mu}H_3) H_4 
       -(\partial^{\mu}H_1)H_2 H_3 (\partial_{mu}H_4)
       \\ \mbox{} \qquad
       - H_1(\partial^{\mu}H_2)(\partial_{\mu}H_3) H_4
       - H_1(\partial^{\mu}H_2) H_3(\partial_{\mu} H_4) 
       - H_1 H_2 (\partial^{\mu}H_3)(\partial_{\mu} H_4) )$ *)
  | Dim6_AHWW_DPB of int   (* %
      $\ii H ( (\partial^{\rho} A_{\mu}) W_{\nu}W_{\rho} g^{\mu\nu} 
       - (\partial^{\nu}A_{\mu})W_{\nu}W_{\rho}g^{\mu\rho})$ *)
  | Dim6_AHWW_DPW of int
    (* %
      $\ii ( ((\partial^{\rho}A_{\mu})W_{\nu}W_{\rho} 
      - (\partial^{\rho} H)A_{\mu}W_{\nu}W_{\rho})g^{\mu\nu}
       \\ \mbox{} \qquad
      (-(\partial^{\nu}A_{\mu})W_{\nu}W_{\rho} 
      + (\partial^{\nu} H)A_{\mu}W_{\nu}W_{\rho})g^{\mu\rho})$ 
    *)
  | Dim6_AHWW_DW of int
    (* %
       $\ii H( (3(\partial^{\rho}A_{\mu})W_{\nu}W_{\rho} 
       - A_{\mu}(\partial^{\rho}W_{\nu})W_{\rho} 
       + A_{\mu}W_{\nu}(\partial^{\rho}W_{\rho})) g^{\mu\nu} 
       \\ \mbox{} \qquad
       + (-3(\partial^{\nu}A_{\mu})W_{\nu}W_{\rho} 
       - A_{\mu}(\partial^{\nu}W_{\nu})W_{\rho} 
       + A_{\mu}W_{\nu}(\partial^{\nu}W_{\rho})) g^{\mu\rho} 
       \\ \mbox{} \qquad
       + 2(A_{\mu}(\partial^{\mu}W_{\nu})W_{\rho} 
       + A_{\mu}W_{\nu}(\partial^{\mu}W_{\rho}))) g^{\nu\rho}) $ 
    *)
  | Dim6_Vector4_DW of int   (*%
      $\ii ( -V_{1,\mu}V_{2,\nu}V^{3,\nu}V^{4,\mu} 
      - V_{1,\mu}V_{2,\nu}V^{3,\mu}V^{4,\nu} \\
      \mbox{} \qquad
      + 2V_{1,\mu}V^{2,\mu}V_{3,\nu}V^{4,\nu} $
    *)
  | Dim6_Vector4_W of int
    (* %
      $\ii (((\partial^{\rho}V_{1,\mu})V_{2}^{\mu}
       (\partial^{\sigma}V_{3,\rho})V_{4,\sigma} 
       + V_{1,\mu}(\partial^{\rho}V_{2}^{\mu})
       (\partial^{\sigma}V_{3,\rho})V_{4,\sigma} 
       \\ \mbox{} \qquad
       + (\partial^{\sigma}V_{1,\mu})V_{2}^{\mu}V_{3,\rho}
       (\partial^{\rho}V_{4,\sigma}) 
       + V_{1,\mu}(\partial^{\sigma}V_{2}^{\mu})V_{3,\rho}
       (\partial^{\rho}V_{4,\sigma}))
       \\ \mbox{} \qquad
       + ((\partial^{\sigma}V_{1,\mu})V_{2,\nu}
       (\partial^{\nu}V_{3}^{\mu})V_{4,\sigma} 
       - V_{1,\mu}(\partial^{\sigma}V_{2,\nu})
       (\partial^{\nu}V_{3}^{\mu})V_{4,\sigma} 
       \\ \mbox{} \qquad
       - (\partial^{\nu}V_{1}^{\mu})V_{2,\nu}
       (\partial^{\sigma}V_{3,\mu})V_{4,\sigma} 
       - (\partial^{\sigma}V_{1,\mu})V_{2,\nu}V_{3}^{\mu}
       (\partial^{\nu}V_{4,\sigma}))
       \\ \mbox{} \qquad
       + ( -(\partial^{\rho}V_{1,\mu})V_{2,\nu}
       (\partial^{\nu}V_{3,\rho})V_{4}^{\mu} 
       + (\partial^{\rho}V_{1,\mu})V_{2,\nu}V_{3,\rho}
       (\partial^{\nu}V_{4}^{\mu}) 
       \\ \mbox{} \qquad
       - V_{1,\mu}(\partial^{\rho}V_{2,\nu})V_{3,\rho}
       (\partial^{\nu}V_{4}^{\mu}) 
       - (\partial^{\nu}V_{1,\mu})V_{2,\nu}V_{3,\rho}
       (\partial^{\rho}V_{4}^{\mu}) ) 
       \\ \mbox{} \qquad
       +( -(\partial^{\sigma}V_{1,\mu})V_{2,\nu}
       (\partial^{\mu}V_{3}^{\nu})V_{4,\sigma} 
       + V_{1,\mu}(\partial^{\sigma}V_{2,\nu})
       (\partial^{\mu}V_{3}^{\nu})V_{4,\sigma} 
       \\ \mbox{} \qquad
       - V_{1,\mu}(\partial^{\mu}V_{2,\nu})
       (\partial^{\sigma}V_{3}^{\nu})V_{4,\sigma} 
       - V_{1,\mu}(\partial^{\sigma}V_{2,\nu})V_{3}^{\nu}
       (\partial^{\mu}V_{4,\sigma})
       \\ \mbox{} \qquad
       + ( -V_{1,\mu}(\partial^{\rho}V_{2,\nu})
       (\partial^{\mu}V_{3,\rho})V_{4}^{\nu} 
       - (\partial^{\rho}V_{1,\mu})V_{2,\nu}V_{3,\rho}
       (\partial^{\mu}V_{4}^{\nu})
       \\ \mbox{} \qquad
       + V_{1,\mu}(\partial^{\rho}V_{2,\nu})V_{3,\rho}
       (\partial^{\mu}V_{4}^{\nu}) 
       - V_{1,\mu}(\partial^{\mu}V_{2,\nu})V_{3,\rho}
       (\partial^{\rho}V_{4}^{\nu}) )
       \\ \mbox{} \qquad
       + ((\partial^{\nu}V_{1,\mu})V_{2,\nu}
       (\partial^{\mu}V_{3,\rho})V_{4}^{\rho} 
       + V_{1,\mu}(\partial^{\mu}V_{2,\nu})
       (\partial^{\nu}V_{3,\rho})V_{4}^{\rho}
       \\ \mbox{} \qquad 
       + (\partial^{\nu}V_{1,\mu})V_{2,\nu}V_{3,\rho}
       (\partial^{\mu}V_{4}^{\rho})
       + V_{1,\mu}(\partial^{\mu}V_{2,\nu})V_{3,\rho}
       (\partial^{\nu}V_{4}^{\rho})) 
       \\ \mbox{} \qquad
       + (\partial^{\rho}V_{1,\mu})V_{2,\nu}V_{3}^{\mu}
       (\partial_{\rho}V_{4}^{\nu})
       - (\partial^{\rho}V_{1,\mu})V_{2}^{\mu}V_{3,\nu}
       (\partial_{\rho}V_{4}^{\nu})
       \\ \mbox{} \qquad
       + V_{1,\mu}(\partial^{\rho}V_{2,\nu})
       (\partial_{\rho}V_{3}^{\mu})V_{4}^{\nu}
       - V_{1,\mu}(\partial^{\rho}V_{2}^{\mu})
       (\partial_{\rho}V_{3,\nu})V_{4}^{\nu}
       \\ \mbox{} \qquad
       + (\partial^{\rho}V_{1,\mu})V_{2,\nu}
       (\partial_{\rho}V_{3}^{\nu})V_{4}^{\mu}
       -  (\partial^{\rho}V_{1,\mu})V_{2}^{\mu}
       (\partial_{\rho}V_{3, \nu})V_{4}^{\nu}
       \\ \mbox{} \qquad
       + V_{1,\mu}(\partial^{\rho}V_{2,\nu})V_{3}^{\nu}
       (\partial_{\rho}V_{4}^{\mu})
       - V_{1,\mu}(\partial^{\rho}V_{2}^{\mu})V_{3,\nu}
       (\partial_{\rho}V_{4}^{\nu}) )$  
    *)
  | Dim6_Scalar2_Vector2_D of int
    (*%
      $\ii H_1 H_2 (-(\partial^{\mu}\partial^{\nu}V_{3,\mu})V_{4,\nu} 
      + (\partial^{\mu}\partial_{\mu}V_{3,\nu})V_{4}^{\nu} \\
      \mbox{}\qquad
      - V_{3,\mu}(\partial^{\mu}\partial^{\nu}V_{4,\nu}) 
      + V_{3,\mu}(\partial^{\nu}\partial_{\nu}V_{4}^{\mu}))$ 
    *)
  | Dim6_Scalar2_Vector2_DP of int
    (*%
      $\ii ((\partial^{\mu}H_1)H_2(\partial^{\nu}V_{3,\mu})V_{4,\nu} 
      - (\partial^{\nu}H_1)H_2(\partial_{\nu}V_{3,\mu})V^{4,\mu} 
      + H_1(\partial^{\mu}H_2)(\partial^{\nu}V_{3,\mu})V_{4,\nu} \\
      \mbox{} \qquad 
      -  H_1(\partial^{\nu}H_2)(\partial_{\nu}V_{3,\mu})V^{4,\mu}
      + (\partial^{\nu}H_1)H_2V_{3,\mu}(\partial^{\mu}V_{4,\nu}) 
      - (\partial^{\nu}H_1)H_2V_{3,\mu}(\partial_{\nu}V^{4,\mu}) \\
      \mbox{} \qquad
      + H_1(\partial^{\nu}H_2)V_{3,\mu}(\partial^{\mu}V_{4,\nu}) 
      - H_1(\partial^{\nu}H_2)V_{3,\mu}(\partial_{\nu}V^{4,\mu})) $
    *)
  | Dim6_Scalar2_Vector2_PB of int   
    (*%
      $\ii (H_1H_2(\partial^{\nu}V_{3,\mu})(\partial^{\mu}V_{4,\nu}) 
      - H_1H_2(\partial^{\nu}V_{3,\mu})(\partial_{\nu}V^{4,\mu})) $
    *)
  | Dim6_HHZZ_T of int   (*% 
    $\ii H_1H_2V_{3,\mu}V^{4,\mu}$ *)
  | Dim6_HWWZ_DW of int
    (* %
       $\ii( H_1(\partial^{\rho}W_{2,\mu})W^{3,\mu}Z_{4,\rho} 
       - H_1W_{2,\mu}(\partial^{\rho}W^{3,\mu})Z_{4,\rho} 
       - 2H_1(\partial^{\nu}W_{2,\mu})W_{3,\nu}Z^{4,\mu} \\
       \mbox{} \qquad
       - H_1W_{2,\mu}(\partial^{\nu}W_{3,\nu})Z^{4,\mu}
       + H_1(\partial^{\mu}W_{2,\mu})W_{3,\nu}Z^{4,\nu} 
       + 2H_1W_{2,\mu}(\partial^{\mu}W_{3,\nu})Z^{4,\nu})$
    *)
  | Dim6_HWWZ_DPB of int
    (* %
      $\ii ( - H_1W_{2,\mu}W_{3,\nu}(\partial^{\nu}Z^{4,\mu}) + 
      H_1W_{2,\mu}W_{3,\nu}(\partial^{\mu}Z^{4,\nu}))$ *)
  | Dim6_HWWZ_DDPW of int
    (* % 
       $ \ii(H_1(\partial^{\nu}W_{2,\mu})W^{3,\mu}Z_{4,\nu} 
       - H_1W_{2,\mu}(\partial^{\nu}W^{3,\mu})Z_{4,\nu} 
       - H_1(\partial^{\nu}W_{2,\mu})W_{3,\nu}Z^{4,\mu} \\
       \mbox{} \qquad
       + H_1W_{2,\mu}W_{3,\nu}(\partial^{\nu}Z^{4,\mu})
       + H_1W_{2,\mu}(\partial^{\mu}W_{3,\nu})Z^{4,\nu} 
       - H_1W_{2,\mu}W_{3,\nu}(\partial^{\mu}Z^{4,\nu}))$ *)   
  | Dim6_HWWZ_DPW of int
    (* %
       $\ii ( H_1(\partial^{\nu}W_{2,\mu})W^{3,\mu}Z_{4,\nu} 
       - H_1W_{2,\mu}(\partial^{\nu}W^{3,\mu})Z_{4,\nu}
       + (\partial^{\nu}H_1)W_{2,\mu}W_{3,\nu}Z^{4,\mu} \\
       \mbox{} \qquad
       - H_1(\partial^{\nu}W_{2,\mu})W_{3,\nu}Z^{4,\mu}
       - (\partial^{\mu}H_1)W_{2,\mu}W_{3,\nu}Z^{4,\nu}
       + H_1W_{2,\mu}(\partial^{\mu}W_{3,\nu})Z^{4,\nu} )$ *)
  | Dim6_AHHZ_D of int
    (* % 
      $\ii (H_1H_2(\partial^{\mu}\partial^{\nu}A_{\mu})Z_{\nu} - 
      H_1H_2(\partial^{\nu}\partial_{\nu}A_{\mu})Z^{\mu})$ *)
  | Dim6_AHHZ_DP of int
    (* %
       $\ii ((\partial^{\mu}H_1)H_2(\partial^{\nu}A_{\mu})Z_{\nu} 
       + H_1(\partial^{\mu}H_2)(\partial^{\nu}A_{\mu})Z_{\nu} \\
       \mbox{} \qquad
       - (\partial^{\nu}H_1)H_2(\partial_{\nu}A_{\mu})Z^{\mu} - 
       H_1(\partial^{\nu}H_2)(\partial_{\nu}A_{\mu})Z^{\mu} ) $ *)
  | Dim6_AHHZ_PB of int
    (* %
       $\ii (H_1H_2(\partial^{\nu}A_{\mu})(\partial_{\nu}Z^{\mu}) - 
       H_1H_2(\partial^{\nu}A_{\mu})(\partial^{\mu}Z_{\nu}))$ *)

type 'a vertexn =
  | UFO of Algebra.QC.t * string * lorentzn * fermion_lines * Color.Vertex.t

(*  An obvious candidate for addition to [boson] is [T], of course. *)

(* \begin{dubious}
     This list is sufficient for the minimal standard model, but not comprehensive
     enough for most of its extensions, supersymmetric or otherwise.
     In particular, we need a \emph{general} parameterization for all trilinear
     vertices.  One straightforward possibility are polynomials in the momenta for
     each combination of fields.
   \end{dubious}
   \begin{JR}
   Here we use the rules which can be found in~\cite{Denner:Majorana}
   and are more properly described in [Targets] where the performing of the fusion
   rules in analytical expressions is encoded.
   \end{JR}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.2}
       \begin{tabular}{|r|l|l|}\hline
                    & only Dirac fermions & incl.~Majorana fermions \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, S, Psi)]:
                              $\mathcal{L}_I=g_S\bar\psi_1 S\psi_2$}\\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot g_S\bar\psi_1 S$
                    & $\psi_2\leftarrow\ii\cdot g_S\psi_1 S$ \\\hline
              [F21] & $\bar\psi_2\leftarrow\ii\cdot g_S S \bar\psi_1$
                    & $\psi_2\leftarrow\ii\cdot g_SS\psi_1$ \\\hline
              [F13] & $S\leftarrow\ii\cdot g_S\bar\psi_1\psi_2$
                    & $S\leftarrow\ii\cdot g_S\psi_1^T{\mathrm{C}}\psi_2$ \\\hline
              [F31] & $S\leftarrow\ii\cdot g_S\psi_{2,\alpha}\bar\psi_{1,\alpha}$
                    & $S\leftarrow\ii\cdot g_S\psi_2^T{\mathrm{C}} \psi_1$\\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot g_SS\psi_2$
                    & $\psi_1\leftarrow\ii\cdot g_SS\psi_2$ \\\hline
              [F32] & $\psi_1\leftarrow\ii\cdot g_S\psi_2 S$
                    & $\psi_1\leftarrow\ii\cdot g_S\psi_2 S$ \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, P, Psi)]:
                              $\mathcal{L}_I=g_P\bar\psi_1 P\gamma_5\psi_2$} \\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot g_P\bar\psi_1\gamma_5 P$
                    & $\psi_2\leftarrow\ii\cdot g_P \gamma_5\psi_1 P$ \\\hline
              [F21] & $\bar\psi_2\leftarrow\ii\cdot g_P P\bar\psi_1\gamma_5$
                    & $\psi_2\leftarrow\ii\cdot g_P P\gamma_5\psi_1$ \\\hline
              [F13] & $P\leftarrow\ii\cdot g_P\bar\psi_1\gamma_5\psi_2$
                    & $P\leftarrow\ii\cdot g_P\psi_1^T {\mathrm{C}}\gamma_5\psi_2$ \\\hline
              [F31] & $P\leftarrow\ii\cdot g_P[\gamma_5\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $P\leftarrow\ii\cdot g_P\psi_2^T {\mathrm{C}}\gamma_5\psi_1$ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot g_P P\gamma_5\psi_2$
                    & $\psi_1\leftarrow\ii\cdot g_P P\gamma_5\psi_2$ \\\hline
              [F32] & $\psi_1\leftarrow\ii\cdot g_P \gamma_5\psi_2 P$
                    & $\psi_1\leftarrow\ii\cdot g_P \gamma_5\psi_2 P$ \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, V, Psi)]:
                              $\mathcal{L}_I=g_V\bar\psi_1\fmslash{V}\psi_2$} \\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot g_V\bar\psi_1\fmslash{V}$
                    & $\psi_{2,\alpha}\leftarrow\ii\cdot
                       (-g_V)\psi_{1,\beta}\fmslash{V}_{\alpha\beta}$ \\\hline
              [F21] & $\bar\psi_{2,\beta}\leftarrow\ii\cdot
                       g_V\fmslash{V}_{\alpha\beta} \bar\psi_{1,\alpha}$
                    & $\psi_2\leftarrow\ii\cdot (-g_V)\fmslash{V}\psi_1$  \\\hline
              [F13] & $V_\mu\leftarrow\ii\cdot g_V\bar\psi_1\gamma_\mu\psi_2$
                    & $V_\mu\leftarrow\ii\cdot
                       g_V (\psi_1)^T {\mathrm{C}}\gamma_{\mu}\psi_2$ \\\hline
              [F31] & $V_\mu\leftarrow\ii\cdot g_V[\gamma_\mu\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $V_\mu\leftarrow\ii\cdot
                       (-g_V)(\psi_2)^T {\mathrm{C}}\gamma_{\mu}\psi_1$ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot g_V\fmslash{V}\psi_2$
                    & $\psi_1\leftarrow\ii\cdot g_V\fmslash{V}\psi_2$ \\\hline
              [F32] & $\psi_{1,\alpha}\leftarrow\ii\cdot
                       g_V\psi_{2,\beta}\fmslash{V}_{\alpha\beta}$
                    & $\psi_{1,\alpha}\leftarrow\ii\cdot
                       g_V\psi_{2,\beta}\fmslash{V}_{\alpha\beta}$ \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, A, Psi)]:
                              $\mathcal{L}_I=g_A\bar\psi_1\gamma_5\fmslash{A}\psi_2$}  \\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot g_A\bar\psi_1\gamma_5\fmslash{A}$
                    & $\psi_{2,\alpha}\leftarrow\ii\cdot
                       g_A\psi_{\beta}[\gamma_5\fmslash{A}]_{\alpha\beta}$ \\\hline
              [F21] & $\bar\psi_{2,\beta}\leftarrow\ii\cdot g_A
                       [\gamma_5\fmslash{A}]_{\alpha\beta} \bar\psi_{1,\alpha}$
                    & $\psi_2\leftarrow\ii\cdot g_A \gamma_5\fmslash{A}\psi$ \\\hline
              [F13] & $A_\mu\leftarrow\ii\cdot g_A\bar\psi_1\gamma_5\gamma_\mu\psi_2$
                    & $A_\mu\leftarrow\ii\cdot
                       g_A \psi_1^T {\textrm{C}}\gamma_5\gamma_{\mu}\psi_2$ \\\hline
              [F31] & $A_\mu\leftarrow\ii\cdot
                       g_A[\gamma_5\gamma_\mu\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $A_\mu\leftarrow\ii\cdot
                       g_A \psi_2^T {\textrm{C}}\gamma_5\gamma_{\mu}\psi_1$ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot g_A\gamma_5\fmslash{A}\psi_2$
                    & $\psi_1\leftarrow\ii\cdot g_A\gamma_5\fmslash{A}\psi_2$ \\\hline
              [F32] & $\psi_{1,\alpha}\leftarrow\ii\cdot g_A
                       \psi_{2,\beta}[\gamma_5\fmslash{A}]_{\alpha\beta}$
                    & $\psi_{1,\alpha}\leftarrow\ii\cdot
                       g_A\psi_{2,\beta}[\gamma_5\fmslash{A}]_{\alpha\beta}$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions} Dimension-4 trilinear fermionic couplings.
       The momenta are unambiguous, because there are no derivative couplings
       and all participating fields are different.}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|r|l|l|}\hline
                    & only Dirac fermions & incl.~Majorana fermions \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, T, Psi)]:
                              $\mathcal{L}_I=g_TT_{\mu\nu}\bar\psi_1
                               [\gamma^\mu,\gamma^\nu]_-\psi_2$}\\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot g_T
                       \bar\psi_1[\gamma^\mu,\gamma^\nu]_-T_{\mu\nu}$
                    & $\bar\psi_2\leftarrow\ii\cdot g_T \cdots$ \\\hline
              [F21] & $\bar\psi_2\leftarrow\ii\cdot g_T T_{\mu\nu}
                       \bar\psi_1[\gamma^\mu,\gamma^\nu]_-$
                    & $\bar\psi_2\leftarrow\ii\cdot g_T \cdots$ \\\hline
              [F13] & $T_{\mu\nu}\leftarrow\ii\cdot g_T\bar\psi_1[\gamma_\mu,\gamma_\nu]_-\psi_2$
                    & $T_{\mu\nu}\leftarrow\ii\cdot g_T \cdots $ \\\hline
              [F31] & $T_{\mu\nu}\leftarrow\ii\cdot g_T
                       [[\gamma_\mu,\gamma_\nu]_-\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $T_{\mu\nu}\leftarrow\ii\cdot g_T \cdots $ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot g_T T_{\mu\nu}[\gamma^\mu,\gamma^\nu]_-\psi_2$
                    & $\psi_1\leftarrow\ii\cdot g_T \cdots$ \\\hline
              [F32] & $\psi_1\leftarrow\ii\cdot g_T [\gamma^\mu,\gamma^\nu]_-\psi_2 T_{\mu\nu}$
                    & $\psi_1\leftarrow\ii\cdot g_T \cdots$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-fermions} Dimension-5 trilinear fermionic couplings
       (NB: the coefficients and signs are not fixed yet).
       The momenta are unambiguous, because there are no derivative couplings
       and all participating fields are different.}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|r|l|l|}\hline
                    & only Dirac fermions & incl.~Majorana fermions \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, SP, Psi)]:
                              $\mathcal{L}_I=\bar\psi_1\phi(g_S+g_P\gamma_5)\psi_2$}\\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot\bar\psi_1(g_S+g_P\gamma_5)\phi$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F21] & $\bar\psi_2\leftarrow\ii\cdot\phi\bar\psi_1(g_S+g_P\gamma_5)$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F13] & $\phi\leftarrow\ii\cdot\bar\psi_1(g_S+g_P\gamma_5)\psi_2$
                    & $\phi\leftarrow\ii\cdot\cdots$ \\\hline
              [F31] & $\phi\leftarrow\ii\cdot[(g_S+g_P\gamma_5)\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $\phi\leftarrow\ii\cdot\cdots$ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot \phi(g_S+g_P\gamma_5)\psi_2$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
              [F32] & $\psi_1\leftarrow\ii\cdot(g_S+g_P\gamma_5)\psi_2\phi$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, SL, Psi)]:
                              $\mathcal{L}_I=g_L\bar\psi_1\phi(1-\gamma_5)\psi_2$}\\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot g_L\bar\psi_1(1-\gamma_5)\phi$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F21] & $\bar\psi_2\leftarrow\ii\cdot g_L\phi\bar\psi_1(1-\gamma_5)$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F13] & $\phi\leftarrow\ii\cdot g_L\bar\psi_1(1-\gamma_5)\psi_2$
                    & $\phi\leftarrow\ii\cdot\cdots$ \\\hline
              [F31] & $\phi\leftarrow\ii\cdot g_L[(1-\gamma_5)\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $\phi\leftarrow\ii\cdot\cdots$ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot g_L\phi(1-\gamma_5)\psi_2$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
              [F32] & $\psi_1\leftarrow\ii\cdot g_L(1-\gamma_5)\psi_2\phi$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, SR, Psi)]:
                              $\mathcal{L}_I=g_R\bar\psi_1\phi(1+\gamma_5)\psi_2$}\\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot g_R\bar\psi_1(1+\gamma_5)\phi$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F21] & $\bar\psi_2\leftarrow\ii\cdot g_R\phi\bar\psi_1(1+\gamma_5)$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F13] & $\phi\leftarrow\ii\cdot g_R\bar\psi_1(1+\gamma_5)\psi_2$
                    & $\phi\leftarrow\ii\cdot\cdots$ \\\hline
              [F31] & $\phi\leftarrow\ii\cdot g_R[(1+\gamma_5)\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $\phi\leftarrow\ii\cdot\cdots$ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot g_R\phi(1+\gamma_5)\psi_2$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
              [F32] & $\psi_1\leftarrow\ii\cdot g_R(1+\gamma_5)\psi_2\phi$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, SLR, Psi)]:
                              $\mathcal{L}_I=g_L\bar\psi_1\phi(1-\gamma_5)\psi_2
                                            +g_R\bar\psi_1\phi(1+\gamma_5)\psi_2$}\\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions-SP} Combined dimension-4 trilinear fermionic couplings.}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|r|l|l|}\hline
                    & only Dirac fermions & incl.~Majorana fermions \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, VA, Psi)]:
                              $\mathcal{L}_I=\bar\psi_1\fmslash{Z}(g_V-g_A\gamma_5)\psi_2$}\\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot\bar\psi_1\fmslash{Z}(g_V-g_A\gamma_5)$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F21] & $\bar\psi_{2,\beta}\leftarrow\ii\cdot
                       [\fmslash{Z}(g_V-g_A\gamma_5)]_{\alpha\beta}\bar\psi_{1,\alpha}$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F13] & $Z_\mu\leftarrow\ii\cdot\bar\psi_1\gamma_\mu(g_V-g_A\gamma_5)\psi_2$
                    & $Z_\mu\leftarrow\ii\cdot \cdots$ \\\hline
              [F31] & $Z_\mu\leftarrow\ii\cdot
                       [\gamma_\mu(g_V-g_A\gamma_5)\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $Z_\mu\leftarrow\ii\cdot \cdots$ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot\fmslash{Z}(g_V-g_A\gamma_5)\psi_2$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
              [F32] & $\psi_{1,\alpha}\leftarrow\ii\cdot
                       \psi_{2,\beta}[\fmslash{Z}(g_V-g_A\gamma_5)]_{\alpha\beta}$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, VL, Psi)]:
                              $\mathcal{L}_I=g_L\bar\psi_1\fmslash{Z}(1-\gamma_5)\psi_2$}\\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot g_L\bar\psi_1\fmslash{Z}(1-\gamma_5)$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F21] & $\bar\psi_{2,\beta}\leftarrow\ii\cdot
                       g_L[\fmslash{Z}(1-\gamma_5)]_{\alpha\beta}\bar\psi_{1,\alpha}$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F13] & $Z_\mu\leftarrow\ii\cdot g_L\bar\psi_1\gamma_\mu(1-\gamma_5)\psi_2$
                    & $Z_\mu\leftarrow\ii\cdot \cdots$ \\\hline
              [F31] & $Z_\mu\leftarrow\ii\cdot
                       g_L[\gamma_\mu(1-\gamma_5)\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $Z_\mu\leftarrow\ii\cdot \cdots$ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot g_L\fmslash{Z}(1-\gamma_5)\psi_2$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
              [F32] & $\psi_{1,\alpha}\leftarrow\ii\cdot
                       g_L\psi_{2,\beta}[\fmslash{Z}(1-\gamma_5)]_{\alpha\beta}$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, VR, Psi)]:
                              $\mathcal{L}_I=g_R\bar\psi_1\fmslash{Z}(1+\gamma_5)\psi_2$}\\\hline
              [F12] & $\bar\psi_2\leftarrow\ii\cdot g_R\bar\psi_1\fmslash{Z}(1+\gamma_5)$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F21] & $\bar\psi_{2,\beta}\leftarrow\ii\cdot
                       g_R[\fmslash{Z}(1+\gamma_5)]_{\alpha\beta}\bar\psi_{1,\alpha}$
                    & $\psi_2\leftarrow\ii\cdot \cdots$ \\\hline
              [F13] & $Z_\mu\leftarrow\ii\cdot g_R\bar\psi_1\gamma_\mu(1+\gamma_5)\psi_2$
                    & $Z_\mu\leftarrow\ii\cdot \cdots$ \\\hline
              [F31] & $Z_\mu\leftarrow\ii\cdot
                       g_R[\gamma_\mu(1+\gamma_5)\psi_2]_\alpha\bar\psi_{1,\alpha}$
                    & $Z_\mu\leftarrow\ii\cdot \cdots$ \\\hline
              [F23] & $\psi_1\leftarrow\ii\cdot g_R\fmslash{Z}(1+\gamma_5)\psi_2$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
              [F32] & $\psi_{1,\alpha}\leftarrow\ii\cdot
                       g_R\psi_{2,\beta}[\fmslash{Z}(1+\gamma_5)]_{\alpha\beta}$
                    & $\psi_1\leftarrow\ii\cdot\cdots$ \\\hline
         \multicolumn{3}{|l|}{[FBF (Psibar, VLR, Psi)]:
                              $\mathcal{L}_I=g_L\bar\psi_1\fmslash{Z}(1-\gamma_5)\psi_2
                                            +g_R\bar\psi_1\fmslash{Z}(1+\gamma_5)\psi_2$}\\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions-VA} Combined dimension-4 trilinear
       fermionic couplings continued.}
   \end{table}
  \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[FBF (Psibar, S, Chi)]: $\bar\psi S\chi$}\\\hline
              [F12] & $\chi\leftarrow\psi S$
            & [F21] & $\chi\leftarrow S \psi$ \\\hline
              [F13] & $S\leftarrow \psi^T{\rm C}\chi$
            & [F31] & $S\leftarrow \chi^T {\rm C}\psi$ \\\hline
              [F23] & $\psi\leftarrow S\chi$
            & [F32] & $\psi\leftarrow\chi S$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Psibar, P, Chi)]: $\bar\psi P\gamma_5\chi$}\\\hline
              [F12] & $\chi\leftarrow \gamma_5 \psi P$
            & [F21] & $\chi\leftarrow P \gamma_5 \psi$ \\\hline
              [F13] & $P\leftarrow \psi^T {\rm C}\gamma_5\chi$
            & [F31] & $P\leftarrow \chi^T {\rm C}\gamma_5\psi$ \\\hline
              [F23] & $\psi\leftarrow P\gamma_5\chi$
            & [F32] & $\psi\leftarrow\gamma_5\chi P$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Psibar, V, Chi)]: $\bar\psi\fmslash{V}\chi$}\\\hline
              [F12] & $\chi_{\alpha}\leftarrow-\psi_{\beta}\fmslash{V}_{\alpha\beta}$
            & [F21] & $\chi\leftarrow-\fmslash{V}\psi$ \\\hline
              [F13] & $V_{\mu}\leftarrow \psi^T {\rm C}\gamma_{\mu}\chi$ 
            & [F31] & $V_{\mu}\leftarrow \chi^T {\rm C}(-\gamma_{\mu}\psi)$ \\\hline
              [F23] & $\psi\leftarrow\fmslash{V}\chi$
            & [F32] & $\psi_\alpha\leftarrow\chi_\beta\fmslash{V}_{\alpha\beta}$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Psibar, A, Chi)]: $\bar\psi\gamma^5\fmslash{A}\chi$}\\\hline
              [F12] & $\chi_{\alpha}\leftarrow\psi_{\beta}\lbrack \gamma^5 \fmslash{A} \rbrack_{\alpha\beta}$
            & [F21] & $\chi\leftarrow\gamma^5\fmslash{A}\psi$ \\\hline
              [F13] & $A_{\mu}\leftarrow \psi^T {\rm C}\gamma^5\gamma_{\mu}\chi$ 
            & [F31] & $A_{\mu}\leftarrow \chi^T {\rm C}(\gamma^5 \gamma_{\mu}\psi)$ \\\hline
              [F23] & $\psi\leftarrow\gamma^5\fmslash{A}\chi$
            & [F32] & $\psi_\alpha\leftarrow\chi_\beta\lbrack \gamma^5 \fmslash{A} \rbrack_{\alpha\beta}$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions-maj} Dimension-4 trilinear couplings
       including one Dirac and one Majorana fermion}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[FBF (Psibar, SP, Chi)]:
                              $\bar\psi\phi(g_S+g_P\gamma_5)\chi$}\\\hline
              [F12] & $\chi \leftarrow (g_S+g_P\gamma_5)\psi \phi$
            & [F21] & $\chi\leftarrow\phi(g_S+g_P\gamma_5)\psi$ \\\hline
              [F13] & $\phi\leftarrow \psi^T {\rm C}(g_S+g_P\gamma_5)\chi$
            & [F31] & $\phi\leftarrow \chi^T {\rm C}(g_S+g_P\gamma_5) \chi$ \\\hline
              [F23] & $\psi\leftarrow \phi(g_S+g_P\gamma_5)\chi$
            & [F32] & $\psi\leftarrow(g_S+g_P\gamma_5)\chi\phi$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Psibar, VA, Chi)]:
                              $\bar\psi\fmslash{Z}(g_V - g_A\gamma_5)\chi$}\\\hline
              [F12] & $\chi_\alpha\leftarrow
                       \psi_\beta[\fmslash{Z}(-g_V-g_A\gamma_5)]_{\alpha\beta}$
            & [F21] & $\chi\leftarrow\fmslash{Z}(-g_V-g_A\gamma_5)]
                       \psi$ \\\hline
              [F13] & $Z_\mu\leftarrow \psi^T {\rm C}\gamma_\mu(g_V-g_A\gamma_5)\chi$
            & [F31] & $Z_\mu\leftarrow \chi^T {\rm C}\gamma_\mu(-g_V-g_A\gamma_5)\psi$ \\\hline
              [F23] & $\psi\leftarrow\fmslash{Z}(g_V-g_A\gamma_5)\chi$
            & [F32] & $\psi_\alpha\leftarrow
                       \chi_\beta[\fmslash{Z}(g_V-g_A\gamma_5)]_{\alpha\beta}$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions-SPVA-maj} Combined dimension-4 trilinear
       fermionic couplings including one Dirac and one Majorana fermion.}
   \end{table}
  \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, S, Psi)]: $\bar\chi S\psi$}\\\hline
              [F12] & $\psi\leftarrow\chi S$ 
            & [F21] & $\psi\leftarrow S\chi$  \\\hline
              [F13] & $S\leftarrow \chi^T {\rm C}\psi$ 
            & [F31] & $S\leftarrow \psi^T {\rm C}\chi$ \\\hline
              [F23] & $\chi\leftarrow S \psi$
            & [F32] & $\chi\leftarrow\psi S$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, P, Psi)]: $\bar\chi P\gamma_5\psi$}\\\hline
              [F12] & $\psi\leftarrow\gamma_5\chi P$ 
            & [F21] & $\psi\leftarrow P\gamma_5\chi$  \\\hline
              [F13] & $P\leftarrow \chi^T {\rm C}\gamma_5\psi$ 
            & [F31] & $P\leftarrow \psi^T {\rm C}\gamma_5\chi$ \\\hline
              [F23] & $\chi\leftarrow P \gamma_5 \psi$
            & [F32] & $\chi\leftarrow \gamma_5 \psi P$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, V, Psi)]: $\bar\chi\fmslash{V}\psi$}\\\hline
              [F12] & $\psi_\alpha\leftarrow-\chi_\beta\fmslash{V}_{\alpha\beta}$ 
            & [F21] & $\psi\leftarrow-\fmslash{V}\chi$  \\\hline
              [F13] & $V_{\mu}\leftarrow \chi^T {\rm C}\gamma_{\mu}\psi$  
            & [F31] & $V_{\mu}\leftarrow \psi^T {\rm C}(-\gamma_{\mu}\chi)$ \\\hline
              [F23] & $\chi\leftarrow\fmslash{V}\psi$
            & [F32] & $\chi_{\alpha}\leftarrow\psi_{\beta}\fmslash{V}_{\alpha\beta}$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, A, Psi)]: $\bar\chi\gamma^5\fmslash{A}\psi$}\\\hline
              [F12] & $\psi_\alpha\leftarrow\chi_\beta\lbrack\gamma^5\fmslash{A} \rbrack_{\alpha\beta}$ 
            & [F21] & $\psi\leftarrow\gamma^5\fmslash{A}\chi$  \\\hline
              [F13] & $A_{\mu}\leftarrow \chi^T {\rm C}(\gamma^5\gamma_{\mu}\psi)$ 
            & [F31] & $A_{\mu}\leftarrow \psi^T {\rm C}\gamma^5\gamma_{\mu}\chi$  \\\hline
              [F23] & $\chi\leftarrow\gamma^5\fmslash{A}\psi$
            & [F32] & $\chi_{\alpha}\leftarrow\psi_{\beta}\lbrack\gamma^5\fmslash{A} \rbrack_{\alpha\beta}$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions-maj'} Dimension-4 trilinear couplings
       including one Dirac and one Majorana fermion}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, SP, Psi)]: $\bar\chi\phi(g_S+g_P\gamma_5)\psi$}\\\hline
              [F12] & $\psi\leftarrow(g_S+g_P\gamma_5)\chi\phi$ 
            & [F21] & $\psi\leftarrow \phi(g_S+g_P\gamma_5)\chi$  \\\hline
              [F13] & $\phi\leftarrow \chi^T {\rm C}(g_S+g_P\gamma_5) \psi$ 
            & [F31] & $\phi\leftarrow \psi^T {\rm C}(g_S+g_P\gamma_5)\chi$ \\\hline
              [F23] & $\chi\leftarrow\phi(g_S+g_P\gamma_5)\psi$
            & [F32] & $\chi \leftarrow (g_S+g_P\gamma_5)\psi \phi$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, VA, Psi)]:
                              $\bar\chi\fmslash{Z}(g_V - g_A\gamma_5)\psi$}\\\hline
              [F12] & $\psi_\alpha\leftarrow
                       \chi_\beta[\fmslash{Z}(-g_V-g_A\gamma_5)]_{\alpha\beta}$
            & [F21] & $\psi\leftarrow\fmslash{Z}(-g_V-g_A\gamma_5)\chi$  \\\hline
              [F13] & $Z_\mu\leftarrow \chi^T {\rm C}\gamma_\mu(g_V-g_A\gamma_5)\psi$ 
            & [F31] & $Z_\mu\leftarrow \psi^T {\rm C}\gamma_\mu(-g_V-g_A\gamma_5)\chi$ \\\hline
              [F23] & $\chi\leftarrow\fmslash{Z}(g_V-g_A\gamma_5)]
                       \psi$
            & [F32] & $\chi_\alpha\leftarrow\psi_\beta[\fmslash{Z}(g_V-g_A\gamma_5)]_{\alpha\beta}$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions-SPVA-maj'} Combined dimension-4 trilinear
       fermionic couplings including one Dirac and one Majorana fermion.}
   \end{table}
  \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, S, Chi)]: $\bar\chi_a S\chi_b$}\\\hline
              [F12] & $\chi_b\leftarrow\chi_a S$
            & [F21] & $\chi_b\leftarrow S \chi_a$ \\\hline
              [F13] & $S\leftarrow \chi^T_a {\rm C}\chi_b$
            & [F31] & $S\leftarrow \chi^T_b {\rm C}\chi_a$ \\\hline
              [F23] & $\chi_a\leftarrow S\chi_b$
            & [F32] & $\chi_a\leftarrow\chi S_b$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, P, Chi)]: $\bar\chi_a P\gamma_5\psi_b$}\\\hline
              [F12] & $\chi_b\leftarrow \gamma_5 \chi_a P$
            & [F21] & $\chi_b\leftarrow P \gamma_5 \chi_a$ \\\hline
              [F13] & $P\leftarrow \chi^T_a {\rm C}\gamma_5\chi_b$
            & [F31] & $P\leftarrow \chi^T_b {\rm C}\gamma_5\chi_a$ \\\hline
              [F23] & $\chi_a\leftarrow P\gamma_5\chi_b$
            & [F32] & $\chi_a\leftarrow\gamma_5\chi_b P$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, V, Chi)]: $\bar\chi_a\fmslash{V}\chi_b$}\\\hline
              [F12] & $\chi_{b,\alpha}\leftarrow-\chi_{a,\beta}\fmslash{V}_{\alpha\beta}$
            & [F21] & $\chi_b\leftarrow-\fmslash{V}\chi_a$ \\\hline
              [F13] & $V_{\mu}\leftarrow \chi^T_a {\rm C}\gamma_{\mu}\chi_b$ 
            & [F31] & $V_{\mu}\leftarrow - \chi^T_b {\rm C}\gamma_{\mu}\chi_a$ \\\hline
              [F23] & $\chi_a\leftarrow\fmslash{V}\chi_b$
            & [F32] & $\chi_{a,\alpha}\leftarrow\chi_{b,\beta}\fmslash{V}_{\alpha\beta}$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, A, Chi)]: $\bar\chi_a\gamma^5\fmslash{A}\chi_b$}\\\hline
              [F12] & $\chi_{b,\alpha}\leftarrow\chi_{a,\beta}\lbrack\gamma^5\fmslash{A} \rbrack_{\alpha\beta}$
            & [F21] & $\chi_b\leftarrow\gamma^5\fmslash{A}\chi_a$ \\\hline
              [F13] & $A_{\mu}\leftarrow \chi^T_a {\rm C}\gamma^5\gamma_{\mu}\chi_b$ 
            & [F31] & $A_{\mu}\leftarrow \chi^T_b {\rm C}(\gamma^5\gamma_{\mu}\chi_a)$ \\\hline
              [F23] & $\chi_a\leftarrow\gamma^5\fmslash{A}\chi_b$
            & [F32] & $\chi_{a,\alpha}\leftarrow\chi_{b,\beta}\lbrack\gamma^5\fmslash{A} \rbrack_{\alpha\beta}$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions-maj2} Dimension-4 trilinear couplings
       of two Majorana fermions}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, SP, Chi)]:
                              $\bar\chi\phi_a(g_S+g_P\gamma_5)\chi_b$}\\\hline
              [F12] & $\chi_b \leftarrow (g_S+g_P\gamma_5)\chi_a \phi$
            & [F21] & $\chi_b\leftarrow\phi(g_S+g_P\gamma_5)\chi_a$ \\\hline
              [F13] & $\phi\leftarrow \chi^T_a {\rm C}(g_S+g_P\gamma_5)\chi_b$
            & [F31] & $\phi\leftarrow \chi^T_b {\rm C}(g_S+g_P\gamma_5) \chi_a$ \\\hline
              [F23] & $\chi_a\leftarrow \phi(g_S+g_P\gamma_5)\chi_b$
            & [F32] & $\chi_a\leftarrow(g_S+g_P\gamma_5)\chi_b\phi$ \\\hline
         \multicolumn{4}{|l|}{[FBF (Chibar, VA, Chi)]:
                              $\bar\chi_a\fmslash{Z}(g_V-g_A\gamma_5)\chi_b$}\\\hline
              [F12] & $\chi_{b,\alpha}\leftarrow\chi_{a,\beta}[\fmslash{Z}(-g_V-g_A\gamma_5)]_{\alpha\beta}$
            & [F21] & $\chi_b\leftarrow\fmslash{Z}(-g_V-g_A\gamma_5)]\chi_a$ \\\hline
              [F13] & $Z_\mu\leftarrow \chi^T_a {\rm C}\gamma_\mu(g_V-g_A\gamma_5)\chi_b$
            & [F31] & $Z_\mu\leftarrow \chi^T_b {\rm C}\gamma_\mu(-g_V-g_A\gamma_5)\chi_a$ \\\hline
              [F23] & $\chi_a\leftarrow\fmslash{Z}(g_V-g_A\gamma_5)\chi_b$
            & [F32] & $\chi_{a,\alpha}\leftarrow
                       \chi_{b,\beta}[\fmslash{Z}(g_V-g_A\gamma_5)]_{\alpha\beta}$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions-SPVA-maj2} Combined dimension-4 trilinear
       fermionic couplings of two Majorana fermions.}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Gauge_Gauge_Gauge]:
                              $\mathcal{L}_I=gf_{abc}
                               A_a^\mu A_b^\nu\partial_\mu A_{c,\nu}$}\\\hline
              [_]  & $A_a^\mu\leftarrow\ii\cdot
                      (-\ii g/2)\cdot C_{abc}^{\mu\rho\sigma}(-k_2-k_3,k_2,k_3)
                      A^b_\rho A^c_\sigma$\\\hline
         \multicolumn{2}{|l|}{[Aux_Gauge_Gauge]:
                              $\mathcal{L}_I=gf_{abc}X_{a,\mu\nu}(k_1)
                               ( A_b^{\mu}(k_2)A_c^{\nu}(k_3)
                                -A_b^{\nu}(k_2)A_c^{\mu}(k_3))$}\\\hline
              [F23]$\lor$[F32] & $X_a^{\mu\nu}(k_2+k_3)\leftarrow\ii\cdot
                                  gf_{abc}( A_b^\mu(k_2)A_c^\nu(k_3)
                                   -A_b^\nu(k_2)A_c^\mu(k_3))$ \\\hline
              [F12]$\lor$[F13] & $A_{a,\mu}(k_1+k_{2/3})\leftarrow\ii\cdot
                                  gf_{abc}X_{b,\nu\mu}(k_1)A_c^\nu(k_{2/3})$ \\\hline
              [F21]$\lor$[F31] & $A_{a,\mu}(k_{2/3}+k_1)\leftarrow\ii\cdot
                                  gf_{abc}A_b^\nu(k_{2/3}) X_{c,\mu\nu}(k_1)$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-bosons} Dimension-4 Vector Boson couplings with
       \emph{outgoing} momenta.
       See~(\ref{eq:C123}) and~(\ref{eq:C123'}) for the definition of the
       antisymmetric tensor $C^{\mu_1\mu_2\mu_3}(k_1,k_2,k_3)$.}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[Scalar_Vector_Vector]:
                              $\mathcal{L}_I=g\phi V_1^\mu V_{2,\mu}$}\\\hline
              [F13] & $\leftarrow\ii\cdot g\cdots$
            & [F31] & $\leftarrow\ii\cdot g\cdots$ \\\hline
              [F12] & $\leftarrow\ii\cdot g\cdots$
            & [F21] & $\leftarrow\ii\cdot g\cdots$ \\\hline
              [F23] & $\phi\leftarrow\ii\cdot g V_1^\mu V_{2,\mu}$
            & [F32] & $\phi\leftarrow\ii\cdot g V_{2,\mu} V_1^\mu$ \\\hline
         \multicolumn{4}{|l|}{[Aux_Vector_Vector]:
                              $\mathcal{L}_I=gX V_1^\mu V_{2,\mu}$}\\\hline
              [F13] & $\leftarrow\ii\cdot g\cdots$
            & [F31] & $\leftarrow\ii\cdot g\cdots$ \\\hline
              [F12] & $\leftarrow\ii\cdot g\cdots$
            & [F21] & $\leftarrow\ii\cdot g\cdots$ \\\hline
              [F23] & $X\leftarrow\ii\cdot g V_1^\mu V_{2,\mu}$
            & [F32] & $X\leftarrow\ii\cdot g V_{2,\mu} V_1^\mu$ \\\hline
         \multicolumn{4}{|l|}{[Aux_Scalar_Vector]:
                              $\mathcal{L}_I=gX^\mu \phi V_\mu$}\\\hline
              [F13] & $\leftarrow\ii\cdot g\cdots$
            & [F31] & $\leftarrow\ii\cdot g\cdots$ \\\hline
              [F12] & $\leftarrow\ii\cdot g\cdots$
            & [F21] & $\leftarrow\ii\cdot g\cdots$ \\\hline
              [F23] & $\leftarrow\ii\cdot g\cdots$
            & [F32] & $\leftarrow\ii\cdot g\cdots$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:scalar-vector}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[Scalar_Scalar_Scalar]:
                              $\mathcal{L}_I=g\phi_1\phi_2\phi_3$}\\\hline
              [F13] & $\phi_2\leftarrow\ii\cdot g \phi_1\phi_3$
            & [F31] & $\phi_2\leftarrow\ii\cdot g \phi_3\phi_1$ \\\hline
              [F12] & $\phi_3\leftarrow\ii\cdot g \phi_1\phi_2$
            & [F21] & $\phi_3\leftarrow\ii\cdot g \phi_2\phi_1$ \\\hline
              [F23] & $\phi_1\leftarrow\ii\cdot g \phi_2\phi_3$
            & [F32] & $\phi_1\leftarrow\ii\cdot g \phi_3\phi_2$ \\\hline
         \multicolumn{4}{|l|}{[Aux_Scalar_Scalar]:
                              $\mathcal{L}_I=gX\phi_1\phi_2$}\\\hline
              [F13] & $\leftarrow\ii\cdot g\cdots$
            & [F31] & $\leftarrow\ii\cdot g\cdots$ \\\hline
              [F12] & $\leftarrow\ii\cdot g\cdots$
            & [F21] & $\leftarrow\ii\cdot g\cdots$ \\\hline
              [F23] & $X\leftarrow\ii\cdot g \phi_1\phi_2$
            & [F32] & $X\leftarrow\ii\cdot g \phi_2\phi_1$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:scalars}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Vector_Scalar_Scalar]:
                              $\mathcal{L}_I=gV^\mu\phi_1
                               \ii\overleftrightarrow{\partial_\mu}\phi_2$}\\\hline
              [F23] & $V^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g(k_2^\mu-k_3^\mu)\phi_1(k_2)\phi_2(k_3)$ \\\hline
              [F32] & $V^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g(k_2^\mu-k_3^\mu)\phi_2(k_3)\phi_1(k_2)$ \\\hline
              [F12] & $\phi_2(k_1+k_2)\leftarrow\ii\cdot
                       g(k_1^\mu+2k_2^\mu)V_\mu(k_1)\phi_1(k_2)$ \\\hline
              [F21] & $\phi_2(k_1+k_2)\leftarrow\ii\cdot
                       g(k_1^\mu+2k_2^\mu)\phi_1(k_2)V_\mu(k_1)$ \\\hline
              [F13] & $\phi_1(k_1+k_3)\leftarrow\ii\cdot
                       g(-k_1^\mu-2k_3^\mu)V_\mu(k_1)\phi_2(k_3)$ \\\hline
              [F31] & $\phi_1(k_1+k_3)\leftarrow\ii\cdot
                       g(-k_1^\mu-2k_3^\mu)\phi_2(k_3)V_\mu(k_1)$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:scalar-current}
       \ldots}
   \end{table} *)
(* \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Aux_DScalar_DScalar]:
                              $\mathcal{L}_I=g\chi
                               (\ii\partial_\mu\phi_1)(\ii\partial^\mu\phi_2)$}\\\hline
              [F23] & $\chi(k_2+k_3)\leftarrow\ii\cdot
                       g (k_2\cdot k_3) \phi_1(k_2) \phi_2(k_3) $ \\\hline
              [F32] & $\chi(k_2+k_3)\leftarrow\ii\cdot
                       g (k_3\cdot k_2) \phi_2(k_3) \phi_1(k_2) $ \\\hline
              [F12] & $\phi_2(k_1+k_2)\leftarrow\ii\cdot
                       g ((-k_1-k_2) \cdot k_2) \chi(k_1) \phi_1(k_2) $ \\\hline
              [F21] & $\phi_2(k_1+k_2)\leftarrow\ii\cdot
                       g (k_2 \cdot (-k_1-k_2)) \phi_1(k_2) \chi(k_1) $ \\\hline
              [F13] & $\phi_1(k_1+k_3)\leftarrow\ii\cdot
                       g ((-k_1-k_3) \cdot k_3) \chi(k_1) \phi_2(k_3) $ \\\hline
              [F31] & $\phi_1(k_1+k_3)\leftarrow\ii\cdot
                       g (k_3 \cdot (-k_1-k_3)) \phi_2(k_3) \chi(k_1) $ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dscalar-dscalar}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Aux_Vector_DScalar]:
                              $\mathcal{L}_I=g\chi V_\mu (\ii\partial^\mu\phi)$}\\\hline
              [F23] & $\chi(k_2+k_3)\leftarrow\ii\cdot
                       g k_3^\mu V_\mu(k_2) \phi(k_3) $ \\\hline
              [F32] & $\chi(k_2+k_3)\leftarrow\ii\cdot
                       g \phi(k_3) k_3^\mu V_\mu(k_2) $ \\\hline
              [F12] & $\phi(k_1+k_2)\leftarrow\ii\cdot
                       g \chi(k_1) (-k_1-k_2)^\mu V_\mu(k_2) $ \\\hline
              [F21] & $\phi(k_1+k_2)\leftarrow\ii\cdot
                       g (-k_1-k_2)^\mu V_\mu(k_2) \chi(k_1) $ \\\hline
              [F13] & $V_\mu(k_1+k_3)\leftarrow\ii\cdot
                       g (-k_1-k_3)_\mu \chi(k_1) \phi(k_3) $ \\\hline
              [F31] & $V_\mu(k_1+k_3)\leftarrow\ii\cdot
                       g (-k_1-k_3)_\mu \phi(k_3) \chi(k_1) $ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:vector-dscalar}
       \ldots}
   \end{table} 
*)




(* Signify which two of three fields are fused: *)
type fuse2 = F23 | F32 | F31 | F13 | F12 | F21

(* Signify which three of four fields are fused: *)
type fuse3 =
  | F123 | F231 | F312 | F132 | F321 | F213
  | F124 | F241 | F412 | F142 | F421 | F214
  | F134 | F341 | F413 | F143 | F431 | F314
  | F234 | F342 | F423 | F243 | F432 | F324

(* Explicit enumeration types make no sense for higher degrees.  *)
type fusen = int list

(* The third member of the triplet will contain the coupling constant: *)
type 'a t =
  | V3 of 'a vertex3 * fuse2 * 'a
  | V4 of 'a vertex4 * fuse3 * 'a
  | Vn of 'a vertexn * fusen * 'a

(* \thocwmodulesection{Gauge Couplings}
   Dimension-4 trilinear vector boson couplings
   \begin{subequations}
   \begin{multline}
     f_{abc}\partial^{\mu}A^{a,\nu}A^b_{\mu}A^c_{\nu} \rightarrow
        \ii f_{abc}k_1^\mu A^{a,\nu}(k_1)A^b_{\mu}(k_2)A^c_{\nu}(k_3) \\
       = -\frac{\ii}{3!} f_{a_1a_2a_3} C^{\mu_1\mu_2\mu_3}(k_1,k_2,k_3)
            A^{a_1}_{\mu_1}(k_1)A^{a_2}_{\mu_2}(k_2)A^{a_3}_{\mu_3}(k_3)
   \end{multline}
   with the totally antisymmetric tensor (under simultaneous permutations
   of all quantum numbers $\mu_i$ and $k_i$) and all momenta \emph{outgoing}
   \begin{equation}
   \label{eq:C123}
     C^{\mu_1\mu_2\mu_3}(k_1,k_2,k_3) =
             (   g^{\mu_1\mu_2} (k_1^{\mu_3}-k_2^{\mu_3})
               + g^{\mu_2\mu_3} (k_2^{\mu_1}-k_3^{\mu_1})
               + g^{\mu_3\mu_1} (k_3^{\mu_2}-k_1^{\mu_2}) )
   \end{equation}
   \end{subequations}
   Since~$f_{a_1a_2a_3}C^{\mu_1\mu_2\mu_3}(k_1,k_2,k_3)$ is totally symmetric
   (under simultaneous permutations of all quantum numbers $a_i$, $\mu_i$ and $k_i$),
   it is easy to take the partial derivative 
   \begin{subequations}
   \label{eq:AofAA}
   \begin{equation}
     A^{a,\mu}(k_2+k_3) =
       - \frac{\ii}{2!} f_{abc}C^{\mu\rho\sigma}(-k_2-k_3,k_2,k_3) A^b_\rho(k_2)A^c_\sigma(k_3)
   \end{equation}
   with
   \begin{equation}
   \label{eq:C123'}
      C^{\mu\rho\sigma}(-k_2-k_3,k_2,k_3) =
             (   g^{\rho\sigma} ( k_2^{\mu}   -k_3^{\mu}   )
               + g^{\mu\sigma}  (2k_3^{\rho}  +k_2^{\rho}  )
               - g^{\mu\rho}    (2k_2^{\sigma}+k_3^{\sigma}) )
   \end{equation}
   i.\,e.
   \begin{multline}
   \label{eq:fuse-gauge}
     A^{a,\mu}(k_2+k_3) = - \frac{\ii}{2!} f_{abc}
         \bigl(   (k_2^{\mu}-k_3^{\mu})A^b(k_2) \cdot A^c(k_3) \\
                + (2k_3+k_2)\cdot A^b(k_2)A^{c,\mu}(k_3)
                - A^{b,\mu}(k_2)A^c(k_3)\cdot(2k_2+k_3) \bigr)
   \end{multline}
   \end{subequations}
   \begin{dubious}
      Investigate the rearrangements proposed in~\cite{HELAS} for improved
      numerical stability.
   \end{dubious} *)

(* \thocwmodulesubsection{Non-Gauge Vector Couplings}
   As a basis for the dimension-4 couplings of three vector bosons, we
   choose ``transversal'' and ``longitudinal'' (with respect to the first
   vector field) tensors that are odd and even under permutation of the
   second and third argument
   \begin{subequations}
   \begin{align}
      \mathcal{L}_T(V_1,V_2,V_3)
        &= V_1^\mu (V_{2,\nu}\ii\overleftrightarrow{\partial_\mu}V_3^\nu)
         = - \mathcal{L}_T(V_1,V_3,V_2) \\
      \mathcal{L}_L(V_1,V_2,V_3)
        &= (\ii\partial_\mu V_1^\mu) V_{2,\nu}V_3^\nu
         = \mathcal{L}_L(V_1,V_3,V_2)
   \end{align}
   \end{subequations}
   Using partial integration in~$\mathcal{L}_L$, we find the
   convenient combinations
   \begin{subequations}
   \begin{align}
     \mathcal{L}_T(V_1,V_2,V_3) + \mathcal{L}_L(V_1,V_2,V_3)
        &= - 2 V_1^\mu \ii\partial_\mu V_{2,\nu} V_3^\nu \\
     \mathcal{L}_T(V_1,V_2,V_3) - \mathcal{L}_L(V_1,V_2,V_3)
        &=   2 V_1^\mu V_{2,\nu} \ii\partial_\mu V_3^\nu
   \end{align}
   \end{subequations}
   As an important example, we can rewrite the dimension-4 ``anomalous'' triple
   gauge couplings
   \begin{multline}
     \ii\mathcal{L}_{\textrm{TGC}}(g_1,\kappa,g_4)/g_{VWW}
        = g_1 V^\mu (W^-_{\mu\nu} W^{+,\nu} - W^+_{\mu\nu} W^{-,\nu}) \\
          + \kappa W^+_\mu W^-_\nu V^{\mu\nu}
          + g_4 W^+_\mu W^-_\nu (\partial^\mu V^\nu + \partial^\nu V^\mu)
   \end{multline}
   as
   \begin{multline}
     \mathcal{L}_{\textrm{TGC}}(g_1,\kappa,g_4)
        =   g_1 \mathcal{L}_T(V,W^-,W^+) \\
          - \frac{\kappa+g_1-g_4}{2} \mathcal{L}_T(W^-,V,W^+)
          + \frac{\kappa+g_1+g_4}{2} \mathcal{L}_T(W^+,V,W^-) \\
          - \frac{\kappa-g_1-g_4}{2} \mathcal{L}_L(W^-,V,W^+)
          + \frac{\kappa-g_1+g_4}{2} \mathcal{L}_L(W^+,V,W^-)
   \end{multline}
   \thocwmodulesubsection{$CP$ Violation}
   \begin{subequations}
   \begin{align}
      \mathcal{L}_{\tilde T}(V_1,V_2,V_3)
        &= V_{1,\mu}(V_{2,\rho}\ii\overleftrightarrow{\partial_\nu}
                     V_{3,\sigma})\epsilon^{\mu\nu\rho\sigma}
         = + \mathcal{L}_T(V_1,V_3,V_2) \\
      \mathcal{L}_{\tilde L}(V_1,V_2,V_3)
        &= (\ii\partial_\mu V_{1,\nu})
             V_{2,\rho}V_{3,\sigma}\epsilon^{\mu\nu\rho\sigma}
         = - \mathcal{L}_L(V_1,V_3,V_2)
   \end{align}
   \end{subequations}
   Here the notations~$\tilde T$ and~$\tilde L$ are clearly
   \textit{abuse de langage}, because
   $\mathcal{L}_{\tilde L}(V_1,V_2,V_3)$ is actually the
   transversal combination, due to the antisymmetry of~$\epsilon$.
   Using partial integration in~$\mathcal{L}_{\tilde L}$, we could again find
   combinations
   \begin{subequations}
   \begin{align}
     \mathcal{L}_{\tilde T}(V_1,V_2,V_3) + \mathcal{L}_{\tilde L}(V_1,V_2,V_3)
        &= - 2 V_{1,\mu} V_{2,\nu} \ii\partial_\rho V_{3,\sigma}
                 \epsilon^{\mu\nu\rho\sigma} \\
     \mathcal{L}_{\tilde T}(V_1,V_2,V_3) - \mathcal{L}_{\tilde L}(V_1,V_2,V_3)
        &= - 2 V_{1,\mu} \ii\partial_\nu V_{2,\rho} V_{3,\sigma}
                 \epsilon^{\mu\nu\rho\sigma}
   \end{align}
   \end{subequations}
   but we don't need them, since
   \begin{multline}
     \ii\mathcal{L}_{\textrm{TGC}}(g_5,\tilde\kappa)/g_{VWW}
       = g_5 \epsilon_{\mu\nu\rho\sigma}
              (W^{+,\mu} \ii\overleftrightarrow{\partial^\rho} W^{-,\nu}) V^\sigma \\
          - \frac{\tilde\kappa_V}{2}  W^-_\mu W^+_\nu \epsilon^{\mu\nu\rho\sigma}
              V_{\rho\sigma}
   \end{multline}
   is immediately recognizable as
   \begin{equation}
     \mathcal{L}_{\textrm{TGC}}(g_5,\tilde\kappa) / g_{VWW}
       = - \ii g_5 \mathcal{L}_{\tilde L}(V,W^-,W^+)
          + \tilde\kappa \mathcal{L}_{\tilde T}(V,W^-,W^+)
   \end{equation}
%%% #procedure decl
%%%   symbol g1, kappa;
%%%   vector V, Wp, Wm, k0, kp, km;
%%%   vector v, V1, V2, V3, k1, k2, k3;
%%%   index mu, nu;
%%% #endprocedure
%%% 
%%% #call decl
%%% 
%%% global L_T(k1,V1,k2,V2,k3,V3)
%%%   = (V1.k2 - V1.k3) * V2.V3;
%%% 
%%% global L_L(k1,V1,k2,V2,k3,V3)
%%%   = - V1.k1 * V2.V3;
%%% 
%%% global L_g1(k1,V1,k2,V2,k3,V3)
%%%   = - V1(mu) * (   (k2(mu)*V2(nu) - k2(nu)*V2(mu)) * V3(nu)
%%%                  - (k3(mu)*V3(nu) - k3(nu)*V3(mu)) * V2(nu) );
%%% 
%%% global L_kappa(k1,V1,k2,V2,k3,V3)
%%%   = (k1(mu)*V1(nu) - k1(nu)*V1(mu)) * V2(mu) * V3(nu);
%%% 
%%% print;
%%% .sort
%%% .store
%%% 
%%% #call decl
%%% 
%%% local lp = L_T(k1,V1,k2,V2,k3,V3) + L_L(k1,V1,k2,V2,k3,V3);
%%% local lm = L_T(k1,V1,k2,V2,k3,V3) - L_L(k1,V1,k2,V2,k3,V3);
%%% print;
%%% .sort
%%% id k1.v? = - k2.v - k3.v;
%%% print;
%%% .sort
%%% .store
%%% 
%%% #call decl
%%% 
%%% local [sum(TL)-g1] = - L_g1(k0,V,km,Wm,kp,Wp)
%%%    + L_T(k0,V,kp,Wp,km,Wm)
%%%    + (L_T(km,Wm,k0,V,kp,Wp) - L_T(kp,Wp,k0,V,km,Wm)) / 2
%%%    - (L_L(km,Wm,k0,V,kp,Wp) - L_L(kp,Wp,k0,V,km,Wm)) / 2;
%%% 
%%% local [sum(TL)-kappa] = - L_kappa(k0,V,km,Wm,kp,Wp)
%%%    + (L_T(km,Wm,k0,V,kp,Wp) - L_T(kp,Wp,k0,V,km,Wm)) / 2
%%%    + (L_L(km,Wm,k0,V,kp,Wp) - L_L(kp,Wp,k0,V,km,Wm)) / 2;
%%% 
%%% local delta =
%%%    - (g1 * L_g1(k0,V,km,Wm,kp,Wp) + kappa * L_kappa(k0,V,km,Wm,kp,Wp))
%%%    + g1 * L_T(k0,V,kp,Wp,km,Wm)
%%%    + (  g1 + kappa) / 2 * (L_T(km,Wm,k0,V,kp,Wp) - L_T(kp,Wp,k0,V,km,Wm))
%%%    + (- g1 + kappa) / 2 * (L_L(km,Wm,k0,V,kp,Wp) - L_L(kp,Wp,k0,V,km,Wm));
%%% 
%%% print;
%%% .sort
%%% 
%%% id k0.v? = - kp.v - km.v;
%%% print;
%%% .sort
%%% .store
%%% 
%%% .end *)

(* \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dim4_Vector_Vector_Vector_T]:
                              $\mathcal{L}_I=gV_1^\mu
                               V_{2,\nu}\ii\overleftrightarrow{\partial_\mu}V_3^\nu$}\\\hline
              [F23] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g(k_2^\mu-k_3^\mu)V_{2,\nu}(k_2)V_3^\nu(k_3)$ \\\hline
              [F32] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g(k_2^\mu-k_3^\mu)V_3^\nu(k_3)V_{2,\nu}(k_2)$ \\\hline
              [F12] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot
                       g(2k_2^\nu+k_1^\nu)V_{1,\nu}(k_1)V_2^\mu(k_2)$ \\\hline
              [F21] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot
                       g(2k_2^\nu+k_1^\nu)V_2^\mu(k_2)V_{1,\nu}(k_1)$ \\\hline
              [F13] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot
                       g(-k_1^\nu-2k_3^\nu)V_1^\nu(k_1)V_3^\mu(k_3)$ \\\hline
              [F31] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot
                       g(-k_1^\nu-2k_3^\nu)V_3^\mu(k_3)V_1^\nu(k_1)$ \\\hline
         \multicolumn{2}{|l|}{[Dim4_Vector_Vector_Vector_L]:
                              $\mathcal{L}_I=g\ii\partial_\mu V_1^\mu
                               V_{2,\nu}V_3^\nu$}\\\hline
              [F23] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g(k_2^\mu+k_3^\mu)V_{2,\nu}(k_2)V_3^\nu(k_3)$ \\\hline
              [F32] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g(k_2^\mu+k_3^\mu)V_3^\nu(k_3)V_{2,\nu}(k_2)$ \\\hline
              [F12] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot
                       g(-k_1^\nu)V_{1,\nu}(k_1)V_2^\mu(k_2)$ \\\hline
              [F21] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot
                       g(-k_1^\nu)V_2^\mu(k_2)V_{1,\nu}(k_1)$ \\\hline
              [F13] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot
                       g(-k_1^\nu)V_1^\nu(k_1)V_3^\mu(k_3)$ \\\hline
              [F31] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot
                       g(-k_1^\nu)V_3^\mu(k_3)V_1^\nu(k_1)$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-TGC}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dim4_Vector_Vector_Vector_T5]:
                              $\mathcal{L}_I=gV_{1,\mu}
                               V_{2,\rho}\ii\overleftrightarrow{\partial_\nu}
                               V_{3,\sigma}\epsilon^{\mu\nu\rho\sigma}$}\\\hline
              [F23] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(k_{2,\nu}-k_{3,\nu})
                       V_{2,\rho}(k_2)V_{3,\sigma}(k_3)$ \\\hline
              [F32] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(k_{2,\nu}-k_{3,\nu})
                       V_{3,\sigma}(k_3)V_{2,\rho}(k_2)$ \\\hline
              [F12] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(2k_{2,\nu}+k_{1,\nu})
                       V_{1,\rho}(k_1)V_{2,\sigma}(k_2)$ \\\hline
              [F21] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(2k_{2,\nu}+k_{1,\nu})
                       V_{2,\sigma}(k_2)V_{1,\rho}(k_1)$ \\\hline
              [F13] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(-k_{1,\nu}-2k_{3,\nu})
                       V_{1,\rho}(k_1)V_{3,\sigma}(k_3)$ \\\hline
              [F31] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(-k_{1,\nu}-2k_{3,\nu})
                       V_{3,\sigma}(k_3)V_{1,\rho}(k_1)$ \\\hline
         \multicolumn{2}{|l|}{[Dim4_Vector_Vector_Vector_L5]:
                              $\mathcal{L}_I=g\ii\partial_\mu V_{1,\nu}
                               V_{2,\nu}V_{3,\sigma}\epsilon^{\mu\nu\rho\sigma}$}\\\hline
              [F23] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(k_{2,\nu}+k_{3,\nu})
                       V_{2,\rho}(k_2)V_{3,\sigma}(k_3)$ \\\hline
              [F32] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(k_{2,\nu}+k_{3,\nu})
                       V_{2,\rho}(k_2)V_{3,\sigma}(k_3)$ \\\hline
              [F12] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(-k_{1,\nu})
                       V_{1,\rho}(k_1)V_{2,\sigma}(k_2)$ \\\hline
              [F21] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(-k_{1,\nu})
                       V_{2,\sigma}(k_2)V_{1,\rho}(k_1)$ \\\hline
              [F13] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(-k_{1,\nu})
                       V_{1,\rho}(k_1)V_{3,\sigma}(k_3)$ \\\hline
              [F31] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot
                       g\epsilon^{\mu\nu\rho\sigma}(-k_{1,\nu})
                       V_{3,\sigma}(k_3)V_{1,\rho}(k_1)$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-TGC5}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dim6_Gauge_Gauge_Gauge]:
                              $\mathcal{L}_I=gF_1^{\mu\nu}F_{2,\nu\rho}
                               F_{3,\hphantom{\rho}\mu}^{\hphantom{3,}\rho}$}\\\hline
              [_]  & $A_1^\mu(k_2+k_3)\leftarrow-\ii\cdot
                      \Lambda^{\mu\rho\sigma}(-k_2-k_3,k_2,k_3)
                      A_{2,\rho} A_{c,\sigma}$\\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim6-TGC}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dim6_Gauge_Gauge_Gauge_5]:
                              $\mathcal{L}_I=g/2\cdot\epsilon^{\mu\nu\lambda\tau}
                               F_{1,\mu\nu}F_{2,\tau\rho}
                               F_{3,\hphantom{\rho}\lambda}^{\hphantom{3,}\rho}$}\\\hline
            [F23]  & $A_1^\mu(k_2+k_3)\leftarrow-\ii\cdot
                      \Lambda_5^{\mu\rho\sigma}(-k_2-k_3,k_2,k_3)
                      A_{2,\rho} A_{3,\sigma}$\\\hline
            [F32]  & $A_1^\mu(k_2+k_3)\leftarrow-\ii\cdot
                      \Lambda_5^{\mu\rho\sigma}(-k_2-k_3,k_2,k_3)
                      A_{3,\sigma} A_{2,\rho}$\\\hline
            [F12]  & $A_3^\mu(k_1+k_2)\leftarrow-\ii\cdot$\\\hline
            [F21]  & $A_3^\mu(k_1+k_2)\leftarrow-\ii\cdot$\\\hline
            [F13]  & $A_2^\mu(k_1+k_3)\leftarrow-\ii\cdot$\\\hline
            [F31]  & $A_2^\mu(k_1+k_3)\leftarrow-\ii\cdot$\\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim6-TGC5}
       \ldots}
   \end{table} *)

(* \thocwmodulesection{$\textrm{SU}(2)$ Gauge Bosons}
   An important special case for table~\ref{tab:dim4-bosons} are the two
   usual coordinates of~$\textrm{SU}(2)$
   \begin{equation}
     W_\pm = \frac{1}{\sqrt2} \left(W_1 \mp \ii W_2\right)
   \end{equation}
   i.\,e.
   \begin{subequations}
   \begin{align}
     W_1 &= \frac{1}{\sqrt2} \left(W_+ + W_-\right) \\
     W_2 &= \frac{\ii}{\sqrt2} \left(W_+ - W_-\right)
   \end{align}
   \end{subequations}
   and
   \begin{equation}
     W_1^\mu W_2^\nu - W_2^\mu W_1^\nu
       = \ii\left(W_-^\mu W_+^\nu - W_+^\mu W_-^\nu\right)
   \end{equation}
   Thus the symmtry remains after the change of basis:
   \begin{multline}
      \epsilon^{abc} W_a^{\mu_1}W_b^{\mu_2}W_c^{\mu_3}
         = \ii W_-^{\mu_1} (W_+^{\mu_2}W_3^{\mu_3} - W_3^{\mu_2}W_+^{\mu_3}) \\
        + \ii W_+^{\mu_1} (W_3^{\mu_2}W_-^{\mu_3} - W_-^{\mu_2}W_3^{\mu_3})
        + \ii W_3^{\mu_1} (W_-^{\mu_2}W_+^{\mu_3} - W_+^{\mu_2}W_-^{\mu_3})
   \end{multline} *)

(* \thocwmodulesection{Quartic Couplings and Auxiliary Fields}
   Quartic couplings can be replaced by cubic couplings to a non-propagating
   auxiliary field. The quartic term should get a negative sign so that it the
   energy is bounded from below for identical fields. In the language of
   functional integrals
   \begin{subequations}
   \label{eq:quartic-aux}
   \begin{multline}
     \mathcal{L}_{\phi^4} = - g^2\phi_1\phi_2\phi_3\phi_4
       \Longrightarrow \\
     \mathcal{L}_{X\phi^2}
       = X^*X \pm gX\phi_1\phi_2 \pm gX^*\phi_3\phi_4
       = (X^* \pm g\phi_1\phi_2)(X \pm g\phi_3\phi_4)
         - g^2\phi_1\phi_2\phi_3\phi_4
   \end{multline}
   and in the language of Feynman diagrams
   \begin{equation}
     \parbox{21mm}{\begin{fmfgraph*}(20,20)
       \fmfleft{e1,e2}
       \fmfright{e3,e4}
       \fmf{plain}{v,e1}
       \fmf{plain}{v,e2}
       \fmf{plain}{v,e3}
       \fmf{plain}{v,e4}
       \fmfv{d.sh=circle,d.si=dot_size,label=$-\ii g^2$}{v}  
     \end{fmfgraph*}}
       \qquad\Longrightarrow\qquad
     \parbox{21mm}{\begin{fmfgraph*}(20,20)
       \fmfleft{e1,e2}
       \fmfright{e3,e4}
       \fmf{plain}{v12,e1}
       \fmf{plain}{v12,e2}
       \fmf{plain}{v34,e3}
       \fmf{plain}{v34,e4}
       \fmf{dashes,label=$+\ii$}{v12,v34}
       \fmfv{d.sh=circle,d.si=dot_size,label=$\pm\ii g$}{v12}  
       \fmfv{d.sh=circle,d.si=dot_size,label=$\pm\ii g$}{v34}  
     \end{fmfgraph*}}
   \end{equation}
   \end{subequations}
   The other choice of signs
   \begin{equation}
     \mathcal{L}_{X\phi^2}'
       = - X^*X \pm gX\phi_1\phi_2 \mp gX^*\phi_3\phi_4
       = - (X^* \pm g\phi_1\phi_2)(X \mp g\phi_3\phi_4)
         - g^2\phi_1\phi_2\phi_3\phi_4
   \end{equation}
   can not be extended easily to identical particles and is therefore
   not used.  For identical particles we have
   \begin{multline}
     \mathcal{L}_{\phi^4} = - \frac{g^2}{4!}\phi^4
       \Longrightarrow \\
     \mathcal{L}_{X\phi^2}
       = \frac{1}{2}X^2 \pm \frac{g}{2}X\phi^2 \pm \frac{g}{2}X\phi^2
       = \frac{1}{2}\left(X \pm \frac{g}{2}\phi^2\right)
                      \left(X \pm \frac{g}{2}\phi^2\right)
         - \frac{g^2}{4!}\phi^4
   \end{multline}
   \begin{dubious}
     Explain the factor~$1/3$ in the functional setting and its
     relation to the three diagrams in the graphical setting?
   \end{dubious}

   \thocwmodulesubsection{Quartic Gauge Couplings}
   \begin{figure}
     \begin{subequations}
     \label{eq:Feynman-QCD}
     \begin{align}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,24)
         \threeexternal{k,,\mu,,a}{p}{p'}
         \fmf{gluon}{v,e1}
         \fmf{fermion}{e2,v,e3}
         \fmfdot{v}  \end{fmfgraph*}}} \,&=
       \begin{split}
           \mbox{} + & \ii g\gamma_\mu T_a
       \end{split} \\
     \label{eq:TGV}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,24)
         \threeexternal{1}{2}{3}
         \fmf{gluon}{v,e1}
         \fmf{gluon}{v,e2}
         \fmf{gluon}{v,e3}
         \threeoutgoing
       \end{fmfgraph*}}} \,&= 
       \begin{split}
           & g f_{a_1a_2a_3} C^{\mu_1\mu_2\mu_3} (k_1,k_2,k_3)
       \end{split} \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,24)
         \fmfsurround{d1,e1,d2,e2,d3,e3,d4,e4}
         \fmf{gluon}{v,e1}
         \fmf{gluon}{v,e2}
         \fmf{gluon}{v,e3}
         \fmf{gluon}{v,e4}
         \fmflabel{1}{e1}
         \fmflabel{2}{e2}
         \fmflabel{3}{e3}
         \fmflabel{4}{e4}
         \fmfdot{v}
         \fmffreeze
         \fmf{warrow_right}{v,e1}
         \fmf{warrow_right}{v,e2}
         \fmf{warrow_right}{v,e3}
         \fmf{warrow_right}{v,e4}
       \end{fmfgraph*}}} \,&= 
       \begin{split}
           \mbox{} - & \ii g^2 f_{a_1a_2b}f_{a_3a_4b}
                       (g_{\mu_1\mu_3} g_{\mu_4\mu_2} - g_{\mu_1\mu_4} g_{\mu_2\mu_3}) \\
           \mbox{} - & \ii g^2 f_{a_1a_3b}f_{a_4a_2b}
                       (g_{\mu_1\mu_4} g_{\mu_2\mu_3} - g_{\mu_1\mu_2} g_{\mu_3\mu_4}) \\
           \mbox{} - & \ii g^2 f_{a_1a_4b}f_{a_2a_3b}
                       (g_{\mu_1\mu_2} g_{\mu_3\mu_4} - g_{\mu_1\mu_3} g_{\mu_4\mu_2})
       \end{split}
     \end{align}
     \end{subequations}
     \caption{\label{fig:gauge-feynman-rules} Gauge couplings.
       See~(\ref{eq:C123}) for the definition of the antisymmetric
       tensor $C^{\mu_1\mu_2\mu_3}(k_1,k_2,k_3)$.}
   \end{figure}
   \begin{figure}
     \begin{equation}
       \label{eq:Feynman-QCD'}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,24)
         \fmfsurround{d1,e1,d2,e2,d3,e3,d4,e4}
         \fmf{gluon}{v12,e1}
         \fmf{gluon}{v12,e2}
         \fmf{gluon}{v34,e3}
         \fmf{gluon}{v34,e4}
         \fmf{dashes}{v12,v34}
         \fmflabel{1}{e1}
         \fmflabel{2}{e2}
         \fmflabel{3}{e3}
         \fmflabel{4}{e4}
         \fmfdot{v12,v34}
         \fmffreeze
         \fmf{warrow_right}{v12,e1}
         \fmf{warrow_right}{v12,e2}
         \fmf{warrow_right}{v34,e3}
         \fmf{warrow_right}{v34,e4}
       \end{fmfgraph*}}} \,= 
           \mbox{} - \ii g^2 f_{a_1a_2b}f_{a_3a_4b}
             (g_{\mu_1\mu_3} g_{\mu_4\mu_2} - g_{\mu_1\mu_4} g_{\mu_2\mu_3})
     \end{equation}
     \caption{\label{fig:gauge-feynman-rules'} Gauge couplings.}
   \end{figure}
   The three crossed versions of
   figure~\ref{fig:gauge-feynman-rules'} reproduces the quartic coupling in
   figure~\ref{fig:gauge-feynman-rules}, because
   \begin{multline}
     - \ii g^2 f_{a_1a_2b}f_{a_3a_4b}
         (g_{\mu_1\mu_3} g_{\mu_4\mu_2} - g_{\mu_1\mu_4} g_{\mu_2\mu_3}) \\
       = (\ii g f_{a_1a_2b} T_{\mu_1\mu_2,\nu_1\nu_2})
         \left(\frac{\ii g^{\nu_1\nu_3} g^{\nu_2\nu_4}}{2}\right)
         (\ii g f_{a_3a_4b} T_{\mu_3\mu_4,\nu_3\nu_4})
   \end{multline}
   with $T_{\mu_1\mu_2,\mu_3\mu_4} = 
     g_{\mu_1\mu_3}g_{\mu_4\mu_2}-g_{\mu_1\mu_4}g_{\mu_2\mu_3}$. *)

(* \thocwmodulesection{Gravitinos and supersymmetric currents}
   In supergravity theories there is a fermionic partner of the graviton, the
   gravitino. Therefore we have introduced the Lorentz type [Vectorspinor]. 
*)

(*    \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[GBG (Fermbar, MOM, Ferm)]:
                              $\bar\psi_1(\ii\fmslash{\partial}\pm m)\phi\psi_2$}\\\hline
              [F12] & $\psi_2\leftarrow-(\fmslash{k}\mp m)\psi_1S$
            & [F21] & $\psi_2\leftarrow-S(\fmslash{k}\mp m)\psi_1$ \\\hline
              [F13] & $S\leftarrow \psi^T_1 {\rm C}(\fmslash{k}\pm m)\psi_2$ 
            & [F31] & $S\leftarrow \psi^T_2 {\rm C}(-(\fmslash{k}\mp m)\psi_1)$ \\\hline
              [F23] & $\psi_1\leftarrow S(\fmslash{k}\pm m)\psi_2$
            & [F32] & $\psi_1\leftarrow(\fmslash{k}\pm m)\psi_2 S$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Fermbar, MOM5, Ferm)]:
                              $\bar\psi_1(\ii\fmslash{\partial}\pm m)\phi\gamma^5\psi_2$}\\\hline
              [F12] & $\psi_2\leftarrow(\fmslash{k}\pm m)\gamma^5\psi_1P$
            & [F21] & $\psi_2\leftarrow P(\fmslash{k}\pm m)\gamma^5\psi_1$ \\\hline
              [F13] & $P\leftarrow \psi^T_1 {\rm C}(\fmslash{k}\pm m)\gamma^5\psi_2$ 
            & [F31] & $P\leftarrow \psi^T_2 {\rm C}(\fmslash{k}\pm m)\gamma^5\psi_1$ \\\hline
              [F23] & $\psi_1\leftarrow P(\fmslash{k}\pm m)\gamma^5\psi_2$
            & [F32] & $\psi_1\leftarrow(\fmslash{k}\pm m)\gamma^5\psi_2 P$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Fermbar, MOML, Ferm)]:
                              $\bar\psi_1 (\ii\fmslash{\partial}\pm m)\phi(1-\gamma^5)\psi_2$}\\\hline
              [F12] & $\psi_2\leftarrow-(1-\gamma^5)(\fmslash{k}\mp m)\psi_1\phi$
            & [F21] & $\psi_2\leftarrow-\phi(1-\gamma^5)(\fmslash{k}\mp m)\psi_1$ \\\hline
              [F13] & $\phi\leftarrow \psi^T_1 {\rm C}(\fmslash{k}\pm m)(1-\gamma^5)\psi_2$ 
            & [F31] & $\phi\leftarrow \psi^T_2 {\rm C}(1-\gamma^5)(-(\fmslash{k}\mp m)\psi_1)$ \\\hline
              [F23] & $\psi_1\leftarrow\phi(\fmslash{k}\pm m)(1-\gamma^5)\psi_2$
            & [F32] & $\psi_1\leftarrow(\fmslash{k}\pm m)(1-\gamma^5)\psi_2 \phi$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Fermbar, LMOM, Ferm)]:
                              $\bar\psi_1 \phi(1-\gamma^5)(\ii\fmslash{\partial}\pm m)\psi_2$}\\\hline
              [F12] & $\psi_2\leftarrow-(\fmslash{k}\mp m)\psi_1(1-\gamma^5)\phi$
            & [F21] & $\psi_2\leftarrow-\phi(\fmslash{k}\mp m)(1-\gamma^5)\psi_1$ \\\hline
              [F13] & $\phi\leftarrow \psi^T_1 {\rm C}(1-\gamma^5)(\fmslash{k}\pm m)\psi_2$ 
            & [F31] & $\phi\leftarrow \psi^T_2 {\rm C}(-(\fmslash{k}\mp m)(1-\gamma^5)\psi_1)$ \\\hline
              [F23] & $\psi_1\leftarrow\phi(1-\gamma^5)(\fmslash{k}\pm m)\psi_2$
            & [F32] & $\psi_1\leftarrow(1-\gamma^5)(\fmslash{k}\pm m)\psi_2 \phi$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Fermbar, VMOM, Ferm)]:
                              $\bar\psi_1 \ii\fmslash{\partial}_\alpha V_\beta \lbrack \gamma^\alpha, \gamma^\beta \rbrack \psi_2$}\\\hline
              [F12] & $\psi_2\leftarrow-\lbrack\fmslash{k},\gamma^\alpha\rbrack\psi_1 V_\alpha$
            & [F21] & $\psi_2\leftarrow-\lbrack\fmslash{k},\fmslash{V}\rbrack\psi_1$ \\\hline
              [F13] & $V_\alpha\leftarrow \psi^T_1 {\rm C}\lbrack\fmslash{k},\gamma_\alpha\rbrack\psi_2$ 
            & [F31] & $V_\alpha\leftarrow \psi^T_2 {\rm C}(-\lbrack\fmslash{k}, \gamma_\alpha\rbrack\psi_1)$ \\\hline
              [F23] & $\psi_1\leftarrow\rbrack\fmslash{k},\fmslash{V}\rbrack\psi_2$
            & [F32] & $\psi_1\leftarrow\lbrack\fmslash{k},\gamma^\alpha\rbrack\psi_2 V_\alpha$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim4-fermions-MOM} Combined dimension-4 trilinear
       fermionic couplings including a momentum. $Ferm$ stands for $Psi$ and
       $Chi$. The case of $MOMR$ is identical to $MOML$ if one substitutes
       $1+\gamma^5$ for $1-\gamma^5$, as well as for $LMOM$ and $RMOM$. The 
       mass term forces us to keep the chiral projector always on the left 
       after "inverting the line" for $MOML$ while on the right for $LMOM$.} 
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{2}{|l|}{[GBBG (Fermbar, S2LR, Ferm)]: $\bar\psi_1 S_1 S_2
(g_L P_L + g_R P_R) \psi_2$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_2\leftarrow S_1 S_2 (g_R P_L + g_L P_R) \psi_1$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_1 \leftarrow S_1 S_2 (g_L P_L + g_R P_R) \psi_2$ \\ \hline
              [F134] [F143] [F314] & $S_1 \leftarrow \psi^T_1 C S_2 (g_L P_L + g_R P_R) \psi_2$ \\ \hline
              [F124] [F142] [F214] & $S_2 \leftarrow \psi^T_1 C S_1 (g_L P_L + g_R P_R) \psi_2$ \\ \hline
              [F413] [F431] [F341] & $S_1 \leftarrow \psi^T_2 C S_2 (g_R P_L + g_L P_R) \psi_1$ \\ \hline
              [F412] [F421] [F241] & $S_2 \leftarrow \psi^T_2 C S_1 (g_R P_L + g_L P_R) \psi_1$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Fermbar, S2, Ferm)]: $\bar\psi_1 S_1 S_2
\gamma^5 \psi_2$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_2\leftarrow S_1 S_2 \gamma^5 \psi_1$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_1 \leftarrow S_1 S_2 \gamma^5 \psi_2$ \\ \hline
              [F134] [F143] [F314] & $S_1 \leftarrow \psi^T_1 C S_2 \gamma^5 \psi_2$ \\ \hline
              [F124] [F142] [F214] & $S_2 \leftarrow \psi^T_1 C S_1 \gamma^5 \psi_2$ \\ \hline
              [F413] [F431] [F341] & $S_1 \leftarrow \psi^T_2 C S_2 \gamma^5 \psi_1$ \\ \hline
              [F412] [F421] [F241] & $S_2 \leftarrow \psi^T_2 C S_1 \gamma^5 \psi_1$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Fermbar, V2, Ferm)]: $\bar\psi_1 \lbrack \fmslash{V}_1 , \fmslash{V}_2 \rbrack \psi_2$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_2\leftarrow - \lbrack \fmslash{V}_1 , \fmslash{V}_2 \rbrack \psi_1$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_1 \leftarrow \lbrack \fmslash{V}_1 , \fmslash{V}_2 \rbrack \psi_2$ \\ \hline
              [F134] [F143] [F314] & $V_{1\:\alpha} \leftarrow \psi^T_1 C \lbrack \gamma_\alpha , \fmslash{V}_2 \rbrack \psi_2$ \\ \hline
              [F124] [F142] [F214] & $V_{2\:\alpha} \leftarrow \psi^T_1 C (-\lbrack \gamma_\alpha , \fmslash{V}_1 \rbrack) \psi_2$ \\ \hline
              [F413] [F431] [F341] & $V_{1\:\alpha} \leftarrow \psi^T_2 C (-\lbrack \gamma_\alpha , \fmslash{V}_2 \rbrack) \psi_1$ \\ \hline
              [F412] [F421] [F241] & $V_{2\:\alpha} \leftarrow \psi^T_2 C \lbrack \gamma_\alpha , \fmslash{V}_1 \rbrack \psi_1$ \\ \hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-mom2} Vertices with two fermions ($Ferm$ stands 
   for $Psi$ and $Chi$, but not for $Grav$) and two bosons (two scalars, 
   scalar/vector, two vectors) for the BRST transformations. Part I}
   \end{table}      
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{2}{|l|}{[GBBG (Fermbar, SV, Ferm)]: $\bar\psi_1 \fmslash{V} S \psi_2$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_2\leftarrow - \fmslash{V} S \psi_1$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_1 \leftarrow \fmslash{V} S \psi_2$ \\ \hline
              [F134] [F143] [F314] & $V_\alpha \leftarrow \psi^T_1 C \gamma_\alpha S \psi_2$ \\ \hline
              [F124] [F142] [F214] & $S \leftarrow \psi^T_1 C \fmslash{V} \psi_2$ \\ \hline
              [F413] [F431] [F341] & $V_\alpha \leftarrow \psi^T_2 C (- \gamma_\alpha S \psi_1)$ \\ \hline
              [F412] [F421] [F241] & $S \leftarrow \psi^T_2 C (- \fmslash{V} \psi_1)$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Fermbar, PV, Ferm)]: $\bar\psi_1 \fmslash{V} \gamma^5 P \psi_2$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_2\leftarrow \fmslash{V} \gamma^5 P \psi_1$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_1 \leftarrow \fmslash{V} \gamma^5 P \psi_2$ \\ \hline
              [F134] [F143] [F314] & $V_\alpha \leftarrow \psi^T_1 C \gamma_\alpha \gamma^5 P \psi_2$ \\ \hline
              [F124] [F142] [F214] & $P \leftarrow \psi^T_1 C \fmslash{V} \gamma^5 \psi_2$ \\ \hline
              [F413] [F431] [F341] & $V_\alpha \leftarrow \psi^T_2 C \gamma_\alpha \gamma^5 P \psi_1$ \\ \hline
              [F412] [F421] [F241] & $P \leftarrow \psi^T_2 C \fmslash{V} \gamma^5 \psi_1$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Fermbar, S(L/R)V, Ferm)]: $\bar\psi_1 \fmslash{V} (1 \mp\gamma^5) \phi \psi_2$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_2\leftarrow - \fmslash{V} (1\pm\gamma^5) \phi \psi_1$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_1 \leftarrow \fmslash{V} (1\mp\gamma^5) \phi \psi_2$ \\ \hline
              [F134] [F143] [F314] & $V_\alpha \leftarrow \psi^T_1 C \gamma_\alpha (1\mp\gamma^5) \phi \psi_2$ \\ \hline
              [F124] [F142] [F214] & $\phi \leftarrow \psi^T_1 C \fmslash{V} (1\mp\gamma^5) \psi_2$ \\ \hline
              [F413] [F431] [F341] & $V_\alpha \leftarrow \psi^T_2 C \gamma_\alpha (-(1\pm\gamma^5) \phi \psi_1)$ \\ \hline
              [F412] [F421] [F241] & $\phi \leftarrow \psi^T_2 C \fmslash{V} (-(1\pm\gamma^5) \psi_1)$ \\ \hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-mom2} Vertices with two fermions ($Ferm$ stands 
   for $Psi$ and $Chi$, but not for $Grav$) and two bosons (two scalars, 
   scalar/vector, two vectors) for the BRST transformations. Part II}
   \end{table}      
  \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[GBG (Gravbar, POT, Psi)]: $\bar\psi_\mu S \gamma^\mu \psi$}\\\hline
              [F12] & $\psi\leftarrow - \gamma^\mu \psi_\mu S$
            & [F21] & $\psi\leftarrow - S\gamma^\mu \psi_\mu$ \\\hline
              [F13] & $S\leftarrow \psi^T_\mu {\rm C} \gamma^\mu \psi$
            & [F31] & $S\leftarrow \psi^T{\rm C} (-\gamma^\mu)\psi_\mu$ \\\hline
              [F23] & $\psi_\mu\leftarrow S\gamma_\mu\psi$
            & [F32] & $\psi_\mu\leftarrow \gamma_\mu \psi S$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Gravbar, S, Psi)]: $\bar\psi_\mu \fmslash{k}_S S \gamma^\mu \psi$}\\\hline
              [F12] & $\psi\leftarrow \gamma^\mu \fmslash{k}_S \psi_\mu S$
            & [F21] & $\psi\leftarrow S\gamma^\mu \fmslash{k}_S \psi_\mu$ \\\hline
              [F13] & $S\leftarrow \psi^T_\mu {\rm C} \fmslash{k}_S \gamma^\mu \psi$
            & [F31] & $S\leftarrow \psi^T{\rm C}\gamma^\mu\fmslash{k}_S \psi_\mu$ \\\hline
              [F23] & $\psi_\mu\leftarrow S\fmslash{k}_S\gamma_\mu\psi$
            & [F32] & $\psi_\mu\leftarrow \fmslash{k}_S \gamma_\mu \psi S$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Gravbar, P, Psi)]: $\bar\psi_\mu \fmslash{k}_P P \gamma^\mu \gamma_5 \psi$}\\\hline
              [F12] & $\psi\leftarrow \gamma^\mu\fmslash{k}_P\gamma_5\psi_\mu P$
            & [F21] & $\psi\leftarrow P\gamma^\mu\fmslash{k}_P\gamma_5\psi_\mu$ \\\hline
              [F13] & $P\leftarrow \psi^T_\mu {\rm C}\fmslash{k}_P\gamma^\mu\gamma_5\psi$
            & [F31] & $P\leftarrow \psi^T {\rm C}\gamma^\mu\fmslash{k}_P\gamma_5\psi_\mu$ \\\hline
              [F23] & $\psi_\mu\leftarrow P\fmslash{k}_P \gamma_\mu \gamma_5 \psi$
            & [F32] & $\psi_\mu\leftarrow \fmslash{k}_P \gamma_\mu \gamma_5 \psi P$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Gravbar, V, Psi)]: $\bar\psi_\mu\lbrack\fmslash{k}_V,\fmslash{V}\rbrack\gamma^\mu\gamma^5\psi$}\\\hline
              [F12] & $\psi\leftarrow \gamma^5\gamma^\mu \lbrack \fmslash{k}_V , \gamma^\alpha \rbrack \psi_\mu V_\alpha$ 
            & [F21] & $\psi\leftarrow \gamma^5\gamma^\mu \lbrack \fmslash{k}_V , \fmslash{V} \rbrack\psi_\mu$ \\\hline
              [F13] & $V_{\mu}\leftarrow \psi^T_\rho {\rm C} \lbrack \fmslash{k}_V , \gamma_\mu \rbrack \gamma^\rho \gamma^5 \psi$
            & [F31] & $V_{\mu}\leftarrow \psi^T {\rm C} \gamma^5 \gamma^{\rho} \lbrack \fmslash{k}_V , \gamma_\mu \rbrack \psi_\rho$ \\\hline
              [F23] & $\psi_\mu\leftarrow\lbrack \fmslash{k}_V , \fmslash{V} \rbrack \gamma_\mu \gamma^5 \psi $
            & [F32] & $\psi_\mu\leftarrow\lbrack \fmslash{k}_V , \gamma^\alpha \rbrack \gamma_\mu \gamma^5 \psi V_\alpha$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-fermions-gravdirac} Dimension-5 trilinear 
       couplings including one Dirac, one Gravitino fermion and one additional particle.The option [POT] is for the coupling of the supersymmetric current to the derivative of the quadratic terms in the superpotential.}
   \end{table}
  \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[GBG (Psibar, POT, Grav)]: $\bar\psi \gamma^\mu S \psi_\mu$}\\\hline
              [F12] & $\psi_\mu\leftarrow - \gamma_\mu \psi S$
            & [F21] & $\psi_\mu\leftarrow - S \gamma_\mu\psi$ \\\hline
              [F13] & $S\leftarrow \psi^T{\rm C}\gamma^\mu\psi_\mu$
            & [F31] & $S\leftarrow \psi^T_\mu {\rm C} (-\gamma^\mu) \psi$ \\\hline
              [F23] & $\psi\leftarrow S\gamma^\mu\psi_\mu$
            & [F32] & $\psi\leftarrow \gamma^\mu\psi_\mu S$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Psibar, S, Grav)]: $\bar\psi \gamma^\mu \fmslash{k}_S S \psi_\mu$}\\\hline
              [F12] & $\psi_\mu\leftarrow \fmslash{k}_S \gamma_\mu \psi S$
            & [F21] & $\psi_\mu\leftarrow S \fmslash{k}_S \gamma_\mu\psi$ \\\hline
              [F13] & $S\leftarrow \psi^T{\rm C}\gamma^\mu\fmslash{k}_S \psi_\mu$
            & [F31] & $S\leftarrow \psi^T_\mu {\rm C} \fmslash{k}_S \gamma^\mu \psi$ \\\hline
              [F23] & $\psi\leftarrow S\gamma^\mu\fmslash{k}_S\psi_\mu$
            & [F32] & $\psi\leftarrow \gamma^\mu\fmslash{k}_S\psi_\mu S$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Psibar, P, Grav)]: $\bar\psi \gamma^\mu\gamma^5 P\fmslash{k}_P \psi_\mu$}\\\hline
              [F12] & $\psi_\mu\leftarrow -\fmslash{k}_P \gamma_\mu \gamma^5 \psi P$
            & [F21] & $\psi_\mu\leftarrow -P\fmslash{k}_P \gamma_\mu \gamma^5 \psi$ \\\hline
              [F13] & $P\leftarrow \psi^T {\rm C}\gamma^\mu\gamma^5\fmslash{k}_P\psi_\mu$
            & [F31] & $P\leftarrow -\psi^T_\mu {\rm C}\fmslash{k}_P\gamma^\mu\gamma_5\psi$ \\\hline
              [F23] & $\psi\leftarrow P \gamma^\mu\gamma^5\fmslash{k}_P\psi_\mu$
            & [F32] & $\psi\leftarrow \gamma^\mu\gamma^5\fmslash{k}_P\psi_\mu P$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Psibar, V, Grav)]: $\bar\psi\gamma^5\gamma^\mu\lbrack\fmslash{k}_V,\fmslash{V}\rbrack\psi_\mu$}\\\hline
              [F12] & $\psi_\mu\leftarrow \lbrack \fmslash{k}_V , \gamma^\alpha \rbrack \gamma_\mu \gamma^5 \psi V_\alpha$ 
            & [F21] & $\psi_\mu\leftarrow \lbrack \fmslash{k}_V , \fmslash{V} \rbrack \gamma_\mu \gamma^5 \psi$ \\\hline
              [F13] & $V_{\mu}\leftarrow \psi^T {\rm C} \gamma^5 \gamma^\rho \lbrack \fmslash{k}_V , \gamma_\mu \rbrack \psi_\rho$ 
            & [F31] & $V_{\mu}\leftarrow \psi^T_\rho {\rm C} \lbrack \fmslash{k}_V , \gamma_\mu \rbrack \gamma^\rho \gamma^5 \psi$ \\\hline
              [F23] & $\psi\leftarrow\gamma^5\gamma^\mu\lbrack \fmslash{k}_V , \fmslash{V} \rbrack\psi_\mu$
            & [F32] & $\psi\leftarrow\gamma^5\gamma^\mu\lbrack \fmslash{k}_V , \gamma^\alpha \rbrack\psi_\mu V_\alpha$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-fermions-diracgrav} Dimension-5 trilinear 
       couplings including one conjugated Dirac, one Gravitino fermion and one additional particle.}
   \end{table}
  \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[GBG (Gravbar, POT, Chi)]: $\bar\psi_\mu S \gamma^\mu \chi$}\\\hline
              [F12] & $\chi\leftarrow - \gamma^\mu \psi_\mu S$
            & [F21] & $\chi\leftarrow - S\gamma^\mu \psi_\mu$ \\\hline
              [F13] & $S\leftarrow \psi^T_\mu {\rm C} \gamma^\mu \chi$
            & [F31] & $S\leftarrow \chi^T{\rm C} (-\gamma^\mu)\psi_\mu$ \\\hline
              [F23] & $\psi_\mu\leftarrow S\gamma_\mu\chi$
            & [F32] & $\psi_\mu\leftarrow \gamma_\mu \chi S$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Gravbar, S, Chi)]: $\bar\psi_\mu \fmslash{k}_S S \gamma^\mu \chi$}\\\hline
              [F12] & $\chi\leftarrow \gamma^\mu \fmslash{k}_S \psi_\mu S$
            & [F21] & $\chi\leftarrow S\gamma^\mu \fmslash{k}_S \psi_\mu$ \\\hline
              [F13] & $S\leftarrow \psi^T_\mu {\rm C} \fmslash{k}_S \gamma^\mu \chi$
            & [F31] & $S\leftarrow \chi^T{\rm C}\gamma^\mu\fmslash{k}_S \psi_\mu$ \\\hline
              [F23] & $\psi_\mu\leftarrow S\fmslash{k}_S\gamma_\mu\chi$
            & [F32] & $\psi_\mu\leftarrow \fmslash{k}_S \gamma_\mu \chi S$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Gravbar, P, Chi)]: $\bar\psi_\mu \fmslash{k}_P P \gamma^\mu \gamma_5 \chi$}\\\hline
              [F12] & $\chi\leftarrow \gamma^\mu\fmslash{k}_P\gamma_5\psi_\mu P$
            & [F21] & $\chi\leftarrow P\gamma^\mu\fmslash{k}_P\gamma_5\psi_\mu$ \\\hline
              [F13] & $P\leftarrow \psi^T_\mu {\rm C}\fmslash{k}_P\gamma^\mu\gamma_5\chi$
            & [F31] & $P\leftarrow \chi^T {\rm C}\gamma^\mu\fmslash{k}_P\gamma_5\psi_\mu$ \\\hline
              [F23] & $\psi_\mu\leftarrow P\fmslash{k}_P \gamma_\mu \gamma_5 \chi$
            & [F32] & $\psi_\mu\leftarrow \fmslash{k}_P \gamma_\mu \gamma_5 \chi P$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Gravbar, V, Chi)]: $\bar\psi_\mu\lbrack\fmslash{k}_V,\fmslash{V}\rbrack\gamma^\mu\gamma^5\chi$}\\\hline
              [F12] & $\chi\leftarrow \gamma^5\gamma^\mu \lbrack \fmslash{k}_V , \gamma^\alpha \rbrack \psi_\mu V_\alpha$ 
            & [F21] & $\chi\leftarrow \gamma^5\gamma^\mu \lbrack \fmslash{k}_V , \fmslash{V} \rbrack\psi_\mu$ \\\hline
              [F13] & $V_{\mu}\leftarrow \psi^T_\rho {\rm C} \lbrack \fmslash{k}_V , \gamma_\mu \rbrack \gamma^\rho \gamma^5 \chi$
            & [F31] & $V_{\mu}\leftarrow \chi^T {\rm C} \gamma^5 \gamma^{\rho} \lbrack \fmslash{k}_V , \gamma_\mu \rbrack \psi_\rho$ \\\hline
              [F23] & $\psi_\mu\leftarrow\lbrack \fmslash{k}_V , \fmslash{V} \rbrack \gamma_\mu \gamma^5 \chi $
            & [F32] & $\psi_\mu\leftarrow\lbrack \fmslash{k}_V , \gamma^\alpha \rbrack \gamma_\mu \gamma^5 \chi V_\alpha$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-fermions-gravmajo} Dimension-5 trilinear 
       couplings including one Majorana, one Gravitino fermion and one 
       additional particle. The table is essentially the same as the one 
       with the Dirac fermion and only written for the sake of completeness.}
   \end{table}
  \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{4}{|l|}{[GBG (Chibar, POT, Grav)]: $\bar\chi \gamma^\mu S \psi_\mu$}\\\hline
              [F12] & $\psi_\mu\leftarrow - \gamma_\mu \chi S$
            & [F21] & $\psi_\mu\leftarrow - S \gamma_\mu\chi$ \\\hline
              [F13] & $S\leftarrow \chi^T{\rm C}\gamma^\mu\psi_\mu$
            & [F31] & $S\leftarrow \psi^T_\mu {\rm C} (-\gamma^\mu) \chi$ \\\hline
              [F23] & $\chi\leftarrow S\gamma^\mu\psi_\mu$
            & [F32] & $\chi\leftarrow \gamma^\mu\psi_\mu S$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Chibar, S, Grav)]: $\bar\chi \gamma^\mu \fmslash{k}_S S \psi_\mu$}\\\hline
              [F12] & $\psi_\mu\leftarrow \fmslash{k}_S \gamma_\mu \chi S$
            & [F21] & $\psi_\mu\leftarrow S \fmslash{k}_S \gamma_\mu\chi$ \\\hline
              [F13] & $S\leftarrow \chi^T{\rm C}\gamma^\mu\fmslash{k}_S \psi_\mu$
            & [F31] & $S\leftarrow \psi^T_\mu {\rm C} \fmslash{k}_S \gamma^\mu \chi$ \\\hline
              [F23] & $\chi\leftarrow S\gamma^\mu\fmslash{k}_S\psi_\mu$
            & [F32] & $\chi\leftarrow \gamma^\mu\fmslash{k}_S\psi_\mu S$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Chibar, P, Grav)]: $\bar\chi \gamma^\mu\gamma^5 P\fmslash{k}_P \psi_\mu$}\\\hline
              [F12] & $\psi_\mu\leftarrow -\fmslash{k}_P \gamma_\mu \gamma^5 \chi P$
            & [F21] & $\psi_\mu\leftarrow -P\fmslash{k}_P \gamma_\mu \gamma^5 \chi$ \\\hline
              [F13] & $P\leftarrow \chi^T {\rm C}\gamma^\mu\gamma^5\fmslash{k}_P\psi_\mu$
            & [F31] & $P\leftarrow -\psi^T_\mu {\rm C}\fmslash{k}_P\gamma^\mu\gamma_5\chi$ \\\hline
              [F23] & $\chi\leftarrow P \gamma^\mu\gamma^5\fmslash{k}_P\psi_\mu$
            & [F32] & $\chi\leftarrow \gamma^\mu\gamma^5\fmslash{k}_P\psi_\mu P$ \\\hline
         \multicolumn{4}{|l|}{[GBG (Chibar, V, Grav)]: $\bar\chi\gamma^5\gamma^\mu\lbrack\fmslash{k}_V,\fmslash{V}\rbrack\psi_\mu$}\\\hline
              [F12] & $\psi_\mu\leftarrow \lbrack \fmslash{k}_V , \gamma^\alpha \rbrack \gamma_\mu \gamma^5 \chi V_\alpha$ 
            & [F21] & $\psi_\mu\leftarrow \lbrack \fmslash{k}_V , \fmslash{V} \rbrack \gamma_\mu \gamma^5 \chi$ \\\hline
              [F13] & $V_{\mu}\leftarrow \chi^T {\rm C} \gamma^5 \gamma^\rho \lbrack \fmslash{k}_V , \gamma_\mu \rbrack \psi_\rho$ 
            & [F31] & $V_{\mu}\leftarrow \psi^T_\rho {\rm C} \lbrack \fmslash{k}_V , \gamma_\mu \rbrack \gamma^\rho \gamma^5 \chi$ \\\hline
              [F23] & $\chi\leftarrow\gamma^5\gamma^\mu\lbrack \fmslash{k}_V , \fmslash{V} \rbrack\psi_\mu$
            & [F32] & $\chi\leftarrow\gamma^5\gamma^\mu\lbrack \fmslash{k}_V , \gamma^\alpha \rbrack\psi_\mu V_\alpha$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-fermions-majograv} Dimension-5 trilinear 
       couplings including one conjugated Majorana, one Gravitino fermion and 
       one additional particle. This table is not only the same as the one 
       with the conjugated Dirac fermion but also the same part of the 
       Lagrangian density as the one with the Majorana particle on the right 
       of the gravitino.}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{2}{|l|}{[GBBG (Gravbar, S2, Psi)]: $\bar\psi_\mu S_1 S_2
\gamma^\mu \psi$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi\leftarrow - \gamma^\mu S_1 S_2 \psi_\mu$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_\mu \leftarrow \gamma_\mu S_1 S_2 \psi$ \\ \hline
              [F134] [F143] [F314] & $S_1 \leftarrow \psi^T_\mu C S_2 \gamma^\mu \psi$ \\ \hline
              [F124] [F142] [F214] & $S_2 \leftarrow \psi^T_\mu C S_1 \gamma^\mu \psi$ \\ \hline
              [F413] [F431] [F341] & $S_1 \leftarrow - \psi^T C S_2 \gamma^\mu \psi_\mu$ \\ \hline
              [F412] [F421] [F241] & $S_2 \leftarrow - \psi^T C S_1 \gamma^\mu \psi_\mu$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Gravbar, SV, Psi)]: $\bar\psi_\mu S \fmslash{V} \gamma^\mu \gamma^5 \psi$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi\leftarrow \gamma^5 \gamma^\mu S \fmslash{V} \psi_\mu$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_\mu \leftarrow \fmslash{V} S \gamma_\mu \gamma^5 \psi$ \\ \hline
              [F134] [F143] [F314] & $S \leftarrow \psi^T_\mu C \fmslash{V} \gamma^\mu \gamma^5 \psi$ \\ \hline
              [F124] [F142] [F214] & $V_\mu \leftarrow \psi^T_\rho C S \gamma_\mu \gamma^\rho \gamma^5 \psi$ \\ \hline
              [F413] [F431] [F341] & $S \leftarrow \psi^T C \gamma^5 \gamma^\mu \fmslash{V} \psi_\mu$ \\ \hline
              [F412] [F421] [F241] & $V_\mu \leftarrow \psi^T C S \gamma^5 \gamma^\rho \gamma_\mu \psi_\rho$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Gravbar, PV, Psi)]: $\bar\psi_\mu P \fmslash{V} \gamma^\mu \psi$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi\leftarrow \gamma^\mu P \fmslash{V} \psi_\mu$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_\mu \leftarrow \fmslash{V} P \gamma_\mu \psi$ \\ \hline
              [F134] [F143] [F314] & $P \leftarrow \psi^T_\mu C \fmslash{V} \gamma^\mu \psi$ \\ \hline
              [F124] [F142] [F214] & $V_\mu \leftarrow \psi^T_\rho C P \gamma_\mu \gamma^\rho \psi$ \\ \hline
              [F413] [F431] [F341] & $P \leftarrow \psi^T C \gamma^\mu \fmslash{V} \psi_\mu$ \\ \hline
              [F412] [F421] [F241] & $V_\mu \leftarrow \psi^T C P \gamma^\rho \gamma_\mu \psi_\rho$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Gravbar, V2, Psi)]: $\bar\psi_\mu f_{abc} \lbrack \fmslash{V}^a , \fmslash{V}^b \rbrack\gamma^\mu \gamma^5 \psi$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi\leftarrow f_{abc} \gamma^5 \gamma^\mu \lbrack \fmslash{V}^a , \fmslash{V}^b \rbrack \psi_\mu$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi_\mu \leftarrow f_{abc} \lbrack \fmslash{V}^a , \fmslash{V}^b \rbrack \gamma_\mu \gamma^5 \psi$ \\ \hline
              [F134] [F143] [F314] [F124] [F142] [F214] & $V_\mu^a \leftarrow\psi^T_\rho C f_{abc} \lbrack \gamma_\mu , \fmslash{V}^b \rbrack \gamma^\rho \gamma^5 \psi$ \\ \hline
              [F413] [F431] [F341] [F412] [F421] [F241] & $V_\mu^a \leftarrow\psi^T C f_{abc} \gamma^5 \gamma^\rho\lbrack \gamma_\mu , \fmslash{V}^b \rbrack \psi_\rho$ \\ \hline 
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-gravferm2boson} Dimension-5 trilinear
       couplings including one Dirac, one Gravitino fermion and two additional bosons. In each lines we list the fusion possibilities with the same order of the fermions, but the order of the bosons is arbitrary (of course, one has to take care of this order in the mapping of the wave functions in [fusion]).}
   \end{table}      
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|r<{:}l|}\hline
         \multicolumn{2}{|l|}{[GBBG (Psibar, S2, Grav)]: $\bar\psi S_1 S_2
\gamma^\mu \psi_\mu$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_\mu\leftarrow - \gamma_\mu S_1 S_2 \psi$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi \leftarrow \gamma^\mu S_1 S_2 \psi_\mu$ \\ \hline
              [F134] [F143] [F314] & $S_1 \leftarrow \psi^T C S_2 \gamma^\mu \psi_\mu$ \\ \hline
              [F124] [F142] [F214] & $S_2 \leftarrow \psi^T C S_1 \gamma^\mu \psi_\mu$ \\ \hline
              [F413] [F431] [F341] & $S_1 \leftarrow - \psi^T_\mu C S_2 \gamma^\mu \psi$ \\ \hline
              [F412] [F421] [F241] & $S_2 \leftarrow - \psi^T_\mu C S_1 \gamma^\mu \psi$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Psibar, SV, Grav)]: $\bar\psi S \gamma^\mu \gamma^5 \fmslash{V} \psi_\mu$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_\mu\leftarrow \fmslash{V} S \gamma^5 \gamma^\mu \psi$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi\leftarrow \gamma^\mu\gamma^5 S\fmslash{V}\psi_\mu$ \\ \hline
              [F134] [F143] [F314] & $S \leftarrow \psi^T C \gamma^\mu \gamma^5 \fmslash{V}\psi$ \\ \hline
              [F124] [F142] [F214] & $V_\mu \leftarrow \psi^T C \gamma^\rho \gamma^5 S \gamma_\mu \psi_\rho$ \\ \hline
              [F413] [F431] [F341] & $S \leftarrow \psi^T_\mu C \fmslash{V} \gamma^5 \gamma^\mu \psi$ \\ \hline
              [F412] [F421] [F241] & $V_\mu \leftarrow \psi^T_\rho C S \gamma_\mu \gamma^5 \gamma^\rho \psi$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Psibar, PV, Grav)]: $\bar\psi P \gamma^\mu \fmslash{V} \psi_\mu$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_\mu\leftarrow \fmslash{V}\gamma_\mu P \psi$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi\leftarrow \gamma^\mu\fmslash{V} P\psi_\mu$ \\ \hline
              [F134] [F143] [F314] & $P \leftarrow \psi^T C \gamma^\mu\fmslash{V}\psi_\mu$ \\ \hline
              [F124] [F142] [F214] & $V_\mu \leftarrow \psi^T C P \gamma^\rho \gamma_\mu \psi_\rho$ \\ \hline
              [F413] [F431] [F341] & $P \leftarrow \psi^T_\mu C \fmslash{V}\gamma^\mu \psi$ \\ \hline
              [F412] [F421] [F241] & $V_\mu \leftarrow \psi^T_\rho C P \gamma_\mu \gamma^\rho \psi$ \\ \hline
         \multicolumn{2}{|l|}{[GBBG (Psibar, V2, Grav)]: $\bar\psi f_{abc} \gamma^5 \gamma^\mu \lbrack \fmslash{V}^a , \fmslash{V}^b \rbrack\psi_\mu$}\\\hline
              [F123] [F213] [F132] [F231] [F312] [F321] & $\psi_\mu\leftarrow f_{abc} \lbrack \fmslash{V}^a , \fmslash{V}^b \rbrack \gamma_\mu \gamma^5 \psi$ \\ \hline
              [F423] [F243] [F432] [F234] [F342] [F324] & $\psi\leftarrow f_{abc} \gamma^5\gamma^\mu\lbrack \fmslash{V}^a , \fmslash{V}^b \rbrack\psi_\mu$ \\ \hline
              [F134] [F143] [F314] [F124] [F142] [F214] & $V_\mu^a \leftarrow\psi^T C f_{abc} \gamma^5\gamma^\rho\lbrack \gamma_\mu , \fmslash{V}^b \rbrack\psi_\rho$ \\ \hline
              [F413] [F431] [F341] [F412] [F421] [F241] & $V_\mu^a \leftarrow\psi^T_\rho C f_{abc}\lbrack \gamma_\mu , \fmslash{V}^b \rbrack\gamma^\rho\gamma^5 \psi$ \\ \hline 
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-gravferm2boson2} Dimension-5 trilinear
       couplings including one conjugated Dirac, one Gravitino fermion and two additional bosons. The couplings of Majorana fermions to the gravitino and two bosons are essentially the same as for Dirac fermions and they are omitted here.}
   \end{table}    
*)

(* \thocwmodulesection{Perturbative Quantum Gravity and Kaluza-Klein Interactions}
   The gravitational coupling constant and the relative strength of
   the dilaton coupling are abbreviated as
   \begin{subequations}
   \begin{align}
     \kappa &= \sqrt{16\pi G_N} \\
     \omega &= \sqrt{\frac{2}{3(n+2)}} = \sqrt{\frac{2}{3(d-2)}}\,,
   \end{align}
   \end{subequations}
   where~$n=d-4$ is the number of extra space dimensions. *)

(* In~(\ref{eq:graviton-feynman-rules3}-\ref{eq:dilaton-feynman-rules5}),
   we use the notation of~\cite{Han/Lykken/Zhang:1999:Kaluza-Klein}:
   \begin{subequations}
   \begin{equation}
     C_{\mu\nu,\rho\sigma} =
       g_{\mu\rho} g_{\nu\sigma} + g_{\mu\sigma} g_{\nu\rho}
         - g_{\mu\nu} g_{\rho\sigma}
   \end{equation}
   \begin{multline}
     D_{\mu\nu,\rho\sigma}(k_1,k_2) =
       g_{\mu\nu} k_{1,\sigma} k_{2,\rho} \\
         \mbox{}
             - (   g_{\mu\sigma} k_{1,\nu} k_{2,\rho}
                 + g_{\mu\rho} k_{1,\sigma} k_{2,\nu}
                 - g_{\rho\sigma} k_{1,\mu} k_{2,\nu}
                 + (\mu\leftrightarrow\nu))
   \end{multline}
   \begin{multline}
     E_{\mu\nu,\rho\sigma}(k_1,k_2) =
       g_{\mu\nu} (k_{1,\rho} k_{1,\sigma}
         + k_{2,\rho} k_{2,\sigma} + k_{1,\rho} k_{2,\sigma}) \\
           \mbox{}
             - (   g_{\nu\sigma} k_{1,\mu} k_{1,\rho}
                 + g_{\nu\rho} k_{2,\mu} k_{2,\sigma}
                 + (\mu\leftrightarrow\nu))
   \end{multline}
   \begin{multline}
     F_{\mu\nu,\rho\sigma\lambda}(k_1,k_2,k_3) = \\
       g_{\mu\rho} g_{\sigma\lambda} (k_2 - k_3)_{\nu}
         + g_{\mu\sigma} g_{\lambda\rho} (k_3 - k_1)_{\nu}
         + g_{\mu\lambda} g_{\rho\sigma} (k_1 - k_2)_{\nu}
             + (\mu\leftrightarrow\nu)
   \end{multline}
   \begin{multline}
     G_{\mu\nu,\rho\sigma\lambda\delta} =
       g_{\mu\nu} (g_{\rho\sigma}g_{\lambda\delta} - g_{\rho\delta}g_{\lambda\sigma})
         \\ \mbox{}
         + ( g_{\mu\rho}g_{\nu\delta}g_{\lambda\sigma}
               + g_{\mu\lambda}g_{\nu\sigma}g_{\rho\delta}
               - g_{\mu\rho}g_{\nu\sigma}g_{\lambda\delta}
               - g_{\mu\lambda}g_{\nu\delta}g_{\rho\sigma}
             + (\mu\leftrightarrow\nu) )
   \end{multline}
   \end{subequations} *)

(* \begin{figure}
     \begin{subequations}
     \label{eq:graviton-feynman-rules3}
     \begin{align}
       \label{eq:graviton-scalar-scalar}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Threeexternal{1}{2}{h_{\mu\nu}}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{dbl_dots}{v,e3}
         \threeoutgoing
       \end{fmfgraph*}}} \,&= 
         \begin{split}
           \mbox{} & - \ii \frac{\kappa}{2} g_{\mu\nu} m^2
             + \ii \frac{\kappa}{2} C_{\mu\nu,\mu_1\mu_2}k^{\mu_1}_1k^{\mu_2}_2
         \end{split} \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Threeexternal{1}{2}{h_{\mu\nu}}
         \fmf{photon}{v,e1}
         \fmf{photon}{v,e2}
         \fmf{dbl_dots}{v,e3}
         \threeoutgoing
       \end{fmfgraph*}}} \,&= 
       \begin{split}
           \mbox{} - \ii \frac{\kappa}{2} m^2 C_{\mu\nu,\mu_1\mu_2}
              - \ii \frac{\kappa}{2}
             (& k_1k_2 C_{\mu\nu,\mu_1\mu_2} \\
              &\mbox{} + D_{\mu\nu,\mu_1\mu_2}(k_1,k_2) \\
              &\mbox{} + \xi^{-1} E_{\mu\nu,\mu_1\mu_2}(k_1,k_2))
       \end{split} \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Threeexternal{p}{p'}{h_{\mu\nu}}
         \fmf{fermion}{e1,v,e2}
         \fmf{dbl_dots}{v,e3}
         \fmfdot{v}
       \end{fmfgraph*}}} \,&= 
       \begin{split}
           \mbox{} - \ii \frac{\kappa}{2} m g_{\mu\nu}
               - \ii \frac{\kappa}{8}
             (& \gamma_{\mu}(p+p')_{\nu} + \gamma_{\nu}(p+p')_{\mu} \\
              & \mbox{} - 2 g_{\mu\nu} (\fmslash{p}+\fmslash{p}') )
       \end{split}
     \end{align}
     \end{subequations}
     \caption{\label{fig:graviton-feynman-rules3} Three-point graviton couplings.}
   \end{figure}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Graviton_Scalar_Scalar]:
                              $h_{\mu\nu} C^{\mu\nu}_{0}(k_1,k_2)\phi_1\phi_2$}\\\hline
              [F12|F21]
                    & $\phi_2 \leftarrow \ii\cdot
                       h_{\mu\nu} C^{\mu\nu}_{0} (k_1, -k-k_1)\phi_1 $ \\\hline
              [F13|F31]
                    & $\phi_1 \leftarrow \ii\cdot
                       h_{\mu\nu} C^{\mu\nu}_{0} (-k-k_2, k_2)\phi_2 $ \\\hline
              [F23|F32]
                    & $h^{\mu\nu} \leftarrow \ii\cdot
                       C^{\mu\nu}_0 (k_1,k_2)\phi_1\phi_2 $ \\\hline
         \multicolumn{2}{|l|}{[Graviton_Vector_Vector]:
                              $h_{\mu\nu} C^{\mu\nu,\mu_1\mu_2}_1(k_1,k_2,\xi)
                               V_{\mu_1}V_{\mu_2} $}\\\hline
              [F12|F21] & $ V^\mu_2 \leftarrow \ii\cdot h_{\kappa\lambda} 
                        C^{\kappa\lambda,\mu\nu}_1(-k-k_1,k_1\xi) V_{1,\nu}$ \\\hline
              [F13|F31] & $ V^\mu_1 \leftarrow \ii\cdot h_{\kappa\lambda} 
                        C^{\kappa\lambda,\mu\nu}_1(-k-k_2,k_2,\xi) V_{2,\nu}$ \\\hline
              [F23|F32]
                    & $h^{\mu\nu} \leftarrow \ii\cdot
                       C^{\mu\nu,\mu_1\mu_2}_1(k_1,k_2,\xi)
                       V_{1,\mu_1}V_{2,\mu_2} $ \\\hline
         \multicolumn{2}{|l|}{[Graviton_Spinor_Spinor]:
                              $h_{\mu\nu} \bar\psi_1
                               C^{\mu\nu}_{\frac{1}{2}}(k_1,k_2)\psi_2 $}\\\hline
              [F12] & $ \bar\psi_2 \leftarrow \ii\cdot
                        h_{\mu\nu} \bar\psi_1 C^{\mu\nu}_{\frac{1}{2}}(k_1,-k-k_1) $ \\\hline
              [F21] & $ \bar\psi_2 \leftarrow \ii\cdot\ldots $ \\\hline
              [F13] & $ \psi_1 \leftarrow \ii\cdot
                        h_{\mu\nu}C^{\mu\nu}_{\frac{1}{2}}(-k-k_2,k_2)\psi_2$ \\\hline
              [F31] & $ \psi_1 \leftarrow \ii\cdot\ldots $ \\\hline
              [F23] & $ h^{\mu\nu} \leftarrow \ii\cdot
                        \bar\psi_1 C^{\mu\nu}_{\frac{1}{2}}(k_1,k_2)\psi_2 $ \\\hline
              [F32] & $ h^{\mu\nu} \leftarrow \ii\cdot\ldots $ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:graviton-three-point} \ldots}
   \end{table}
   Derivation of~(\ref{eq:graviton-scalar-scalar})
   \begin{subequations}
   \begin{align}
     L &= \frac{1}{2} (\partial_\mu \phi) (\partial^\mu \phi) - \frac{m^2}{2} \phi^2 \\
     (\partial_\mu\phi) \frac{\partial L}{\partial(\partial^\nu\phi)}
       &= (\partial_\mu\phi)(\partial_\nu\phi) \\
     T_{\mu\nu} &= -g_{\mu\nu} L +
       (\partial_\mu\phi) \frac{\partial L}{\partial(\partial^\nu\phi)}
       + 
   \end{align}
   \end{subequations}
   \begin{subequations}
   \begin{align}
     C^{\mu\nu}_{0}(k_1,k_2)
       &= C^{\mu\nu,\mu_1\mu_2} k_{1,\mu_1} k_{2,\mu_2} \\
     C^{\mu\nu,\mu_1\mu_2}_1(k_1,k_2,\xi)
       &= k_1k_2 C^{\mu\nu,\mu_1\mu_2}
            + D^{\mu\nu,\mu_1\mu_2}(k_1,k_2)
            + \xi^{-1} E^{\mu\nu,\mu_1\mu_2}(k_1,k_2) \\
     C^{\mu\nu}_{\frac{1}{2},\alpha\beta}(p,p')
       &=   \gamma^{\mu}_{\alpha\beta}(p+p')^{\nu}
          + \gamma^{\nu}_{\alpha\beta}(p+p')^{\mu}
          - 2 g^{\mu\nu} (\fmslash{p}+\fmslash{p}')_{\alpha\beta}
   \end{align}
   \end{subequations} *)

(* \begin{figure}
     \begin{subequations}
     \label{eq:dilaton-feynman-rules3}
     \begin{align}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Threeexternal{1}{2}{\phi(k)}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{dots}{v,e3}
         \threeoutgoing
       \end{fmfgraph*}}} \,&=
          - \ii \omega \kappa 2m^2 - \ii \omega \kappa k_1k_2 \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Threeexternal{1}{2}{\phi(k)}
         \fmf{photon}{v,e1}
         \fmf{photon}{v,e2}
         \fmf{dots}{v,e3}
         \threeoutgoing
       \end{fmfgraph*}}} \,&= 
             - \ii \omega \kappa g_{\mu_1\mu_2}m^2
             - \ii \omega \kappa
                 \xi^{-1} (k_{1,\mu_1}k_{\mu_2} + k_{2,\mu_2}k_{\mu_1}) \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Threeexternal{p}{p'}{\phi(k)}
         \fmf{fermion}{e1,v,e2}
         \fmf{dots}{v,e3}
         \fmfdot{v}
       \end{fmfgraph*}}} \,&= 
           - \ii \omega \kappa 2m
           + \ii \omega \kappa \frac{3}{4}(\fmslash{p}+\fmslash{p}')
     \end{align}
     \end{subequations}
     \caption{\label{fig:dilaton-feynman-rules3} Three-point dilaton couplings.}
   \end{figure}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.4}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dilaton_Scalar_Scalar]:
                              $\phi \ldots k_1k_2\phi_1\phi_2 $}\\\hline
              [F12|F21] & $ \phi_2 \leftarrow \ii\cdot k_1(-k-k_1)\phi\phi_1 $ \\\hline
              [F13|F31] & $ \phi_1 \leftarrow \ii\cdot (-k-k_2)k_2\phi\phi_2 $ \\\hline
              [F23|F32] & $ \phi   \leftarrow \ii\cdot k_1k_2\phi_1\phi_2 $ \\\hline
         \multicolumn{2}{|l|}{[Dilaton_Vector_Vector]:
                              $\phi \ldots $}\\\hline
              [F12] & $ V_{2,\mu} \leftarrow \ii\cdot\ldots $ \\\hline
              [F21] & $ V_{2,\mu} \leftarrow \ii\cdot\ldots $ \\\hline
              [F13] & $ V_{1,\mu} \leftarrow \ii\cdot\ldots $ \\\hline
              [F31] & $ V_{1,\mu} \leftarrow \ii\cdot\ldots $ \\\hline
              [F23] & $ \phi \leftarrow \ii\cdot\ldots $ \\\hline
              [F32] & $ \phi \leftarrow \ii\cdot\ldots $ \\\hline
         \multicolumn{2}{|l|}{[Dilaton_Spinor_Spinor]:
                              $\phi \ldots $}\\\hline
              [F12] & $ \bar\psi_2 \leftarrow \ii\cdot\ldots $ \\\hline
              [F21] & $ \bar\psi_2 \leftarrow \ii\cdot\ldots $ \\\hline
              [F13] & $ \psi_1 \leftarrow \ii\cdot\ldots $ \\\hline
              [F31] & $ \psi_1 \leftarrow \ii\cdot\ldots $ \\\hline
              [F23] & $ \phi \leftarrow \ii\cdot\ldots $ \\\hline
              [F32] & $ \phi \leftarrow \ii\cdot\ldots $ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dilaton-three-point} \ldots}
   \end{table} *)

(* \begin{figure}
     \begin{subequations}
     \label{eq:graviton-feynman-rules4}
     \begin{align}
       \label{eq:graviton-scalar-scalar-scalar}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{h_{\mu\nu}}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{plain}{v,e3}
         \fmf{dbl_dots}{v,e4}
         \fouroutgoing
       \end{fmfgraph*}}} \,&=
         \begin{split}
           \mbox{} & ??? 
         \end{split} \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{h_{\mu\nu}}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{dbl_dots}{v,e4}
         \fouroutgoing
       \end{fmfgraph*}}} \,&=
         \begin{split}
           \mbox{} & 
              - \ii g\frac{\kappa}{2} C_{\mu\nu,\mu_3\rho}(k_1-k_2)^{\rho} T^{a_3}_{n_2n_1}
         \end{split} \\
       \label{eq:graviton-scalar-vector-vector}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{h_{\mu\nu}}
         \fmf{plain}{v,e1}
         \fmf{photon}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{dbl_dots}{v,e4}
         \fouroutgoing
       \end{fmfgraph*}}} \,&=
         \begin{split}
           \mbox{} & ??? 
         \end{split} \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{h_{\mu\nu}}
         \fmf{photon}{v,e1}
         \fmf{photon}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{dbl_dots}{v,e4}
         \fouroutgoing
       \end{fmfgraph*}}} \,&=
         \begin{split}
           \mbox{} - g \frac{\kappa}{2} f^{a_1a_2a_3}
               (&           C_{\mu\nu,\mu_1\mu_2} (k_1-k_2)_{\mu_3} \\
                & \mbox{} + C_{\mu\nu,\mu_2\mu_3} (k_2-k_3)_{\mu_1} \\
                & \mbox{} + C_{\mu\nu,\mu_3\mu_1} (k_3-k_1)_{\mu_2} \\
                & \mbox{} + F_{\mu\nu,\mu_1\mu_2\mu_3}(k_1,k_2,k_3) )
         \end{split} \\
       \label{eq:graviton-yukawa}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{h_{\mu\nu}}
         \fmf{fermion}{e1,v,e2}
         \fmf{plain}{v,e3}
         \fmf{dbl_dots}{v,e4}
         \fmfdot{v}
         \fmffreeze
         \fmf{warrow_right}{v,e3}
         \fmf{warrow_right}{v,e4}
       \end{fmfgraph*}}} \,&=
         \begin{split}
           \mbox{} & ???
         \end{split} \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{h_{\mu\nu}}
         \fmf{fermion}{e1,v,e2}
         \fmf{photon}{v,e3}
         \fmf{dbl_dots}{v,e4}
         \fmfdot{v}
         \fmffreeze
         \fmf{warrow_right}{v,e3}
         \fmf{warrow_right}{v,e4}
       \end{fmfgraph*}}} \,&=
         \begin{split}
           \mbox{} & \ii g\frac{\kappa}{4}
              (C_{\mu\nu,\mu_3\rho} - g_{\mu\nu}g_{\mu_3\rho})
              \gamma^{\rho} T^{a_3}_{n_2n_1}
         \end{split}
     \end{align}
     \end{subequations}
     \caption{\label{fig:graviton-feynman-rules4} Four-point graviton couplings.
       (\ref{eq:graviton-scalar-scalar-scalar}),
       (\ref{eq:graviton-scalar-vector-vector}),
       and~(\ref{eq:graviton-yukawa)} are missing
       in~\cite{Han/Lykken/Zhang:1999:Kaluza-Klein}, but should be generated
       by standard model Higgs selfcouplings, Higgs-gaugeboson couplings, and
       Yukawa couplings.}
   \end{figure} *)

(* \begin{figure}
     \begin{subequations}
     \label{eq:dilaton-feynman-rules4}
     \begin{align}
       \label{eq:dilaton-scalar-scalar-scalar}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{\phi(k)}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{plain}{v,e3}
         \fmf{dots}{v,e4}
         \fouroutgoing
       \end{fmfgraph*}}} \,&= ??? \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{\phi(k)}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{dots}{v,e4}
         \fouroutgoing
       \end{fmfgraph*}}} \,&=
          - \ii \omega \kappa (k_1 + k_2)_{\mu_3} T^{a_3}_{n_1,n_2} \\
       \label{eq:dilaton-scalar-vector-vector}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{\phi(k)}
         \fmf{plain}{v,e1}
         \fmf{photon}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{dots}{v,e4}
         \fouroutgoing
       \end{fmfgraph*}}} \,&= ??? \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{\phi(k)}
         \fmf{photon}{v,e1}
         \fmf{photon}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{dots}{v,e4}
         \fouroutgoing
       \end{fmfgraph*}}} \,&= 0 \\
       \label{eq:dilaton-yukawa}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{h_{\mu\nu}}
         \fmf{fermion}{e1,v,e2}
         \fmf{plain}{v,e3}
         \fmf{dots}{v,e4}
         \fmfdot{v}
         \fmffreeze
         \fmf{warrow_right}{v,e3}
         \fmf{warrow_right}{v,e4}
       \end{fmfgraph*}}} \,&= ??? \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fourexternal{1}{2}{3}{\phi(k)}
         \fmf{fermion}{e1,v,e2}
         \fmf{photon}{v,e3}
         \fmf{dots}{v,e4}
         \fmfdot{v}
         \fmffreeze
         \fmf{warrow_right}{v,e3}
         \fmf{warrow_right}{v,e4}
       \end{fmfgraph*}}} \,&=
         - \ii \frac{3}{2} \omega g \kappa \gamma_{\mu_3} T^{a_3}_{n_1n_2}
     \end{align}
     \end{subequations}
     \caption{\label{fig:dilaton-feynman-rules4} Four-point dilaton couplings.
       (\ref{eq:dilaton-scalar-scalar-scalar}),
       (\ref{eq:dilaton-scalar-vector-vector})
       and~(\ref{eq:dilaton-yukawa}) are missing
       in~\cite{Han/Lykken/Zhang:1999:Kaluza-Klein}, but could be generated
       by standard model Higgs selfcouplings, Higgs-gaugeboson couplings,
       and Yukawa couplings.}
   \end{figure} *)

(* \begin{figure}
     \begin{subequations}
     \label{eq:graviton-feynman-rules5}
     \begin{align}
       \label{eq:graviton-scalar-scalar-scalar-scalar}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fiveexternal{1}{2}{3}{4}{h_{\mu\nu}}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{plain}{v,e3}
         \fmf{plain}{v,e4}
         \fmf{dots}{v,e5}
         \fiveoutgoing
       \end{fmfgraph*}}} \,&=
         \begin{split}
           \mbox{} & ???
         \end{split} \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fiveexternal{1}{2}{3}{4}{h_{\mu\nu}}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{photon}{v,e4}
         \fmf{dots}{v,e5}
         \fiveoutgoing
       \end{fmfgraph*}}} \,&=
         \begin{split}
           \mbox{} & - \ii g^2 \frac{\kappa}{2} C_{\mu\nu,\mu_3\mu_4}
                (T^{a_3}T^{a_4} + T^{a_4}T^{a_3})_{n_2n_1}
         \end{split} \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fiveexternal{1}{2}{3}{4}{h_{\mu\nu}}
         \fmf{photon}{v,e1}
         \fmf{photon}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{photon}{v,e4}
         \fmf{dots}{v,e5}
         \fiveoutgoing
       \end{fmfgraph*}}} \,&=
         \begin{split}
           \mbox{} - \ii g^2 \frac{\kappa}{2} 
              (& f^{ba_1a_3} f^{ba_2a_4} G_{\mu\nu,\mu_1\mu_2\mu_3\mu_4} \\
               & \mbox + f^{ba_1a_2} f^{ba_3a_4} G_{\mu\nu,\mu_1\mu_3\mu_2\mu_4} \\
               & \mbox + f^{ba_1a_4} f^{ba_2a_3} G_{\mu\nu,\mu_1\mu_2\mu_4\mu_3} )
         \end{split}
     \end{align}
     \end{subequations}
     \caption{\label{fig:graviton-feynman-rules5} Five-point graviton couplings.
       (\ref{eq:graviton-scalar-scalar-scalar-scalar}) is missing
       in~\cite{Han/Lykken/Zhang:1999:Kaluza-Klein}, but should be generated
       by standard model Higgs selfcouplings.}
   \end{figure} *)

(* \begin{figure}
     \begin{subequations}
     \label{eq:dilaton-feynman-rules5}
     \begin{align}
       \label{eq:dilaton-scalar-scalar-scalar-scalar}
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fiveexternal{1}{2}{3}{4}{\phi(k)}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{plain}{v,e3}
         \fmf{plain}{v,e4}
         \fmf{dots}{v,e5}
         \fiveoutgoing
       \end{fmfgraph*}}} \,&= ??? \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fiveexternal{1}{2}{3}{4}{\phi(k)}
         \fmf{plain}{v,e1}
         \fmf{plain}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{photon}{v,e4}
         \fmf{dots}{v,e5}
         \fiveoutgoing
       \end{fmfgraph*}}} \,&=
          \ii \omega g^2 \kappa g_{\mu_3\mu_4}
                (T^{a_3}T^{a_4} + T^{a_4}T^{a_3})_{n_2n_1} \\
       \parbox{28mm}{\fmfframe(2,2)(2,1){\begin{fmfgraph*}(24,22)
         \Fiveexternal{1}{2}{3}{4}{\phi(k)}
         \fmf{photon}{v,e1}
         \fmf{photon}{v,e2}
         \fmf{photon}{v,e3}
         \fmf{photon}{v,e4}
         \fmf{dots}{v,e5}
         \fiveoutgoing
       \end{fmfgraph*}}} \,&= 0
     \end{align}
     \end{subequations}
     \caption{\label{fig:dilaton-feynman-rules5} Five-point dilaton couplings.
       (\ref{eq:dilaton-scalar-scalar-scalar-scalar}) is missing
       in~\cite{Han/Lykken/Zhang:1999:Kaluza-Klein}, but could be generated
       by standard model Higgs selfcouplings.}
   \end{figure} *)

(* \thocwmodulesection{Dependent Parameters}
   This is a simple abstract syntax for parameter dependencies.
   Later, there will be a parser for a convenient concrete syntax
   as a part of a concrete syntax for models.  There is no intention
   to do \emph{any} symbolic manipulation with this.  The expressions
   will be translated directly by [Targets] to the target language.  *)

type 'a expr =
  | I
  | Integer of int
  | Float of float
  | Atom of 'a
  | Sum of 'a expr list
  | Diff of 'a expr * 'a expr 
  | Neg of 'a expr
  | Prod of 'a expr list
  | Quot of 'a expr * 'a expr 
  | Rec of 'a expr
  | Pow of 'a expr * int
  | PowX of 'a expr * 'a expr
  | Sqrt of 'a expr
  | Sin of 'a expr
  | Cos of 'a expr
  | Tan of 'a expr
  | Cot of 'a expr
  | Asin of 'a expr
  | Acos of 'a expr
  | Atan of 'a expr
  | Atan2 of 'a expr * 'a expr
  | Sinh of 'a expr
  | Cosh of 'a expr
  | Tanh of 'a expr
  | Exp of 'a expr
  | Log of 'a expr
  | Log10 of 'a expr
  | Conj of 'a expr
  | Abs of 'a expr

type 'a variable = Real of 'a | Complex of 'a
type 'a variable_array = Real_Array of 'a | Complex_Array of 'a

type 'a parameters =
    { input : ('a * float) list;
      derived : ('a variable * 'a expr) list;
      derived_arrays : ('a variable_array * 'a expr list) list }

(* \thocwmodulesection{More Exotic Couplings}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dim5_Scalar_Vector_Vector_T]:
                              $\mathcal{L}_I=g\phi
                               (\ii\partial_\mu V_1^\nu)(\ii\partial_\nu V_2^\mu)$}\\\hline
              [F23] & $\phi(k_2+k_3)\leftarrow\ii\cdot g
                         k_3^\mu V_{1,\mu}(k_2) k_2^\nu V_{2,\nu}(k_3)$ \\\hline
              [F32] & $\phi(k_2+k_3)\leftarrow\ii\cdot g
                         k_2^\mu V_{2,\mu}(k_3) k_3^\nu V_{1,\nu}(k_2)$ \\\hline
              [F12] & $V_2^\mu(k_1+k_2)\leftarrow\ii\cdot g
                         k_2^\mu \phi(k_1) (-k_1^\nu-k_2^\nu) V_{1,\nu}(k_2)$ \\\hline
              [F21] & $V_2^\mu(k_1+k_2)\leftarrow\ii\cdot g
                         k_2^\mu (-k_1^\nu-k_2^\nu)V_{1,\nu}(k_2) \phi(k_1)$ \\\hline
              [F13] & $V_1^\mu(k_1+k_3)\leftarrow\ii\cdot g
                         k_3^\mu \phi(k_1) (-k_1^\nu-k_3^\nu)V_{2,\nu}(k_3)$ \\\hline
              [F31] & $V_1^\mu(k_1+k_3)\leftarrow\ii\cdot g
                         k_3^\mu (-k_1^\nu-k_3^\nu)V_{2,\nu}(k_3) \phi(k_1)$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-scalar-vector-vector}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dim6_Vector_Vector_Vector_T]:
                              $\mathcal{L}_I=gV_1^\mu
                               ((\ii\partial_\nu V_2^\rho)%
                                \ii\overleftrightarrow{\partial_\mu}
                                (\ii\partial_\rho V_3^\nu))$}\\\hline
              [F23] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot g
                        (k_2^\mu - k_3^\mu) k_3^\nu  V_{2,\nu} (k_2)
                                            k_2^\rho V_{3,\rho}(k_3)$ \\\hline
              [F32] & $V_1^\mu(k_2+k_3)\leftarrow\ii\cdot g
                        (k_2^\mu - k_3^\mu) k_2^\nu  V_{3,\nu} (k_3)
                                            k_3^\rho V_{2,\rho}(k_2)$ \\\hline
              [F12] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot g
                        k_2^\mu (k_1^\nu+2k_2^\nu)   V_{1,\nu} (k_1)
                                (-k_1^\rho-k_2^\rho) V_{2,\rho}(k_2)$ \\\hline
              [F21] & $V_3^\mu(k_1+k_2)\leftarrow\ii\cdot g
                        k_2^\mu (-k_1^\rho-k_2^\rho) V_{2,\rho}(k_2)
                                (k_1^\nu+2k_2^\nu)   V_{1,\nu} (k_1)$ \\\hline
              [F13] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot g
                        k_3^\mu (k_1^\nu+2k_3^\nu)   V_{1,\nu} (k_1)
                                (-k_1^\rho-k_3^\rho) V_{3,\rho}(k_3)$ \\\hline
              [F31] & $V_2^\mu(k_1+k_3)\leftarrow\ii\cdot g
                        k_3^\mu (-k_1^\rho-k_3^\rho) V_{3,\rho}(k_3)
                                (k_1^\nu+2k_3^\nu)   V_{1,\nu} (k_1)$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim6-vector-vector-vector}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Tensor_2_Vector_Vector]:
                              $\mathcal{L}_I=gT^{\mu\nu}
                               (V_{1,\mu}V_{2,\nu} + V_{1,\nu}V_{2,\mu})$}\\\hline
              [F23] & $T^{\mu\nu}(k_2+k_3)\leftarrow\ii\cdot g
                       (V_{1,\mu}(k_2) V_{2,\nu}(k_3) + V_{1,\nu}(k_2) V_{2,\mu}(k_3))$ \\\hline
              [F32] & $T^{\mu\nu}(k_2+k_3)\leftarrow\ii\cdot g
                       (V_{2,\nu}(k_3) V_{1,\mu}(k_2) + V_{2,\mu}(k_3) V_{1,\nu}(k_2))$ \\\hline
              [F12] & $V_2^\mu(k_1+k_2)\leftarrow\ii\cdot g
                       (T^{\mu\nu}(k_1) + T^{\nu\mu}(k_1)) V_{1,\nu}(k_2)$ \\\hline
              [F21] & $V_2^\mu(k_1+k_2)\leftarrow\ii\cdot g
                       V_{1,\nu}(k_2)(T^{\mu\nu}(k_1) + T^{\nu\mu}(k_1))$ \\\hline
              [F13] & $V_1^\mu(k_1+k_3)\leftarrow\ii\cdot g
                       (T^{\mu\nu}(k_1) + T^{\nu\mu}(k_1)) V_{2,\nu}(k_3)$ \\\hline
              [F31] & $V_1^\mu(k_1+k_3)\leftarrow\ii\cdot g
                       V_{2,\nu}(k_3) (T^{\mu\nu}(k_1) + T^{\nu\mu}(k_1))$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:tensor2-vector-vector}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dim5_Tensor_2_Vector_Vector_1]:
                              $\mathcal{L}_I=gT^{\alpha\beta} 
                                (V_1^\mu
                                 \ii\overleftrightarrow\partial_\alpha
                                 \ii\overleftrightarrow\partial_\beta V_{2,\mu})$}\\\hline
              [F23] & $T^{\alpha\beta}(k_2+k_3)\leftarrow\ii\cdot g
                        (k_2^\alpha-k_3^\alpha)(k_2^\beta-k_3^\beta)
                        V_1^\mu(k_2)V_{2,\mu}(k_3)$ \\\hline
              [F32] & $T^{\alpha\beta}(k_2+k_3)\leftarrow\ii\cdot g
                        (k_2^\alpha-k_3^\alpha)(k_2^\beta-k_3^\beta)
                        V_{2,\mu}(k_3)V_1^\mu(k_2)$ \\\hline
              [F12] & $V_2^\mu(k_1+k_2)\leftarrow\ii\cdot g
                        (k_1^\alpha+2k_2^\alpha) (k_1^\beta+2k_2^\beta)
                         T_{\alpha\beta}(k_1) V_1^\mu(k_2)$ \\\hline
              [F21] & $V_2^\mu(k_1+k_2)\leftarrow\ii\cdot g
                        (k_1^\alpha+2k_2^\alpha) (k_1^\beta+2k_2^\beta)
                         V_1^\mu(k_2) T_{\alpha\beta}(k_1)$ \\\hline
              [F13] & $V_1^\mu(k_1+k_3)\leftarrow\ii\cdot g
                        (k_1^\alpha+2k_3^\alpha) (k_1^\beta+2k_3^\beta)
                         T_{\alpha\beta}(k_1) V_2^\mu(k_3)$ \\\hline
              [F31] & $V_1^\mu(k_1+k_3)\leftarrow\ii\cdot g
                        (k_1^\alpha+2k_3^\alpha) (k_1^\beta+2k_3^\beta)
                         V_2^\mu(k_3) T_{\alpha\beta}(k_1)$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-tensor2-vector-vector-1}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dim5_Tensor_2_Vector_Vector_2]:
                              $\mathcal{L}_I=gT^{\alpha\beta}
            (   V_1^\mu \ii\overleftrightarrow\partial_\beta (\ii\partial_\mu V_{2,\alpha})
              + V_1^\mu \ii\overleftrightarrow\partial_\alpha (\ii\partial_\mu V_{2,\beta}))
                 $}\\\hline
              [F23] & $T^{\alpha\beta}(k_2+k_3)\leftarrow\ii\cdot g
                       (k_3^\beta-k_2^\beta) k_3^\mu V_{1,\mu}(k_2)V_2^\alpha(k_3)
                         + (\alpha\leftrightarrow\beta)$ \\\hline
              [F32] & $T^{\alpha\beta}(k_2+k_3)\leftarrow\ii\cdot g
                       (k_3^\beta-k_2^\beta) V_2^\alpha(k_3) k_3^\mu V_{1,\mu}(k_2)
                         + (\alpha\leftrightarrow\beta)$ \\\hline
              [F12] & $V_2^\alpha(k_1+k_2)\leftarrow\ii\cdot g
                         (k_1^\beta+2k_2^\beta) 
                         (T^{\alpha\beta}(k_1)+T^{\beta\alpha}(k_1))
                         (k_1^\mu+k_2^\mu) V_{1,\mu}(k_2)$ \\\hline
              [F21] & $V_2^\alpha(k_1+k_2)\leftarrow\ii\cdot g
                         (k_1^\mu+k_2^\mu) V_{1,\mu}(k_2)
                         (k_1^\beta+2k_2^\beta)
                         (T^{\alpha\beta}(k_1)+T^{\beta\alpha}(k_1))$ \\\hline
              [F13] & $V_1^\alpha(k_1+k_3)\leftarrow\ii\cdot g
                         (k_1^\beta+2k_3^\beta) 
                         (T^{\alpha\beta}(k_1)+T^{\beta\alpha}(k_1))
                         (k_1^\mu+k_3^\mu) V_{2,\mu}(k_3)$ \\\hline
              [F31] & $V_1^\alpha(k_1+k_3)\leftarrow\ii\cdot g
                         (k_1^\mu+k_3^\mu) V_{2,\mu}(k_3)
                         (k_1^\beta+2k_3^\beta)
                         (T^{\alpha\beta}(k_1)+T^{\beta\alpha}(k_1))$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim5-tensor2-vector-vector-1'}
       \ldots}
   \end{table}
   \begin{table}
     \begin{center}
       \renewcommand{\arraystretch}{1.3}
       \begin{tabular}{|>{\qquad}r<{:}l|}\hline
         \multicolumn{2}{|l|}{[Dim7_Tensor_2_Vector_Vector_T]:
                              $\mathcal{L}_I=gT^{\alpha\beta}
                                ((\ii\partial^\mu V_1^\nu)
                                 \ii\overleftrightarrow\partial_\alpha
                                 \ii\overleftrightarrow\partial_\beta
                                 (\ii\partial_\nu V_{2,\mu}))$}\\\hline
              [F23] & $T^{\alpha\beta}(k_2+k_3)\leftarrow\ii\cdot g
                        (k_2^\alpha-k_3^\alpha)(k_2^\beta-k_3^\beta)
                        k_3^\mu V_{1,\mu}(k_2) k_2^\nu V_{2,\nu}(k_3)$ \\\hline
              [F32] & $T^{\alpha\beta}(k_2+k_3)\leftarrow\ii\cdot g
                        (k_2^\alpha-k_3^\alpha)(k_2^\beta-k_3^\beta)
                        k_2^\nu V_{2,\nu}(k_3) k_3^\mu V_{1,\mu}(k_2)$ \\\hline
              [F12] & $V_2^\mu(k_1+k_2)\leftarrow\ii\cdot g
                        k_2^\mu
                        (k_1^\alpha+2k_2^\alpha) (k_1^\beta+2k_2^\beta)
                         T_{\alpha\beta}(k_1) (-k_1^\nu-k_2^\nu)V_{1,\nu}(k_2)$ \\\hline
              [F21] & $V_2^\mu(k_1+k_2)\leftarrow\ii\cdot g
                        k_2^\mu (-k_1^\nu-k_2^\nu)V_{1,\nu}(k_2) 
                        (k_1^\alpha+2k_2^\alpha) (k_1^\beta+2k_2^\beta)
                         T_{\alpha\beta}(k_1)$ \\\hline
              [F13] & $V_1^\mu(k_1+k_3)\leftarrow\ii\cdot g
                        k_3^\mu
                        (k_1^\alpha+2k_3^\alpha) (k_1^\beta+2k_3^\beta)
                         T_{\alpha\beta}(k_1) (-k_1^\nu-k_3^\nu) V_{2,\nu}(k_3)$ \\\hline
              [F31] & $V_1^\mu(k_1+k_3)\leftarrow\ii\cdot g
                        k_3^\mu (-k_1^\nu-k_3^\nu) V_{2,\nu}(k_3) 
                        (k_1^\alpha+2k_3^\alpha) (k_1^\beta+2k_3^\beta)
                         T_{\alpha\beta}(k_1)$ \\\hline
       \end{tabular}
     \end{center}
     \caption{\label{tab:dim7-tensor2-vector-vector-T}
       \ldots}
   \end{table} *)
