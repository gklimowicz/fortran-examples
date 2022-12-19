(* modellib_SM.mli --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
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

(* \thocwmodulesection{Hardcoded Models} *)

module Phi3 : Model.T with module Ch = Charges.Null
module Phi4 : Model.T with module Ch = Charges.Null
module QED : Model.T with module Ch = Charges.ZZ
module QCD : Model.T with module Ch = Charges.ZZ

module type SM_flags =
  sig
    val higgs_triangle : bool (* $H\gamma\gamma$, $Hg\gamma$ and $Hgg couplings *)
    val higgs_hmm : bool
    val triple_anom : bool
    val quartic_anom : bool
    val higgs_anom : bool
    val dim6 : bool  
    val k_matrix : bool
    val ckm_present : bool   
    val top_anom : bool
    val top_anom_4f : bool
    val tt_threshold : bool
  end

module SM_no_anomalous : SM_flags
module SM_anomalous : SM_flags
module SM_k_matrix : SM_flags
module SM_no_anomalous_ckm : SM_flags
module SM_anomalous_ckm : SM_flags
module SM_Higgs : SM_flags
module SM_Higgs_CKM : SM_flags
module SM_anomalous_top : SM_flags
module SM_tt_threshold : SM_flags
module SM_dim6 : SM_flags
  
module SM : functor (F : SM_flags) -> Model.Gauge with module Ch = Charges.QQ

module SM_Rxi : Model.T with module Ch = Charges.QQ

module Groves : functor (M : Model.Gauge) -> Model.Gauge with module Ch = M.Ch
module SM_clones : Model.Gauge with module Ch = Charges.QQ
