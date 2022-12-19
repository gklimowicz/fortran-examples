(* modellib_BSM.mli --

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


(* \thocwmodulesection{More Hardcoded BSM Models} *)

module type BSM_flags = 
  sig 
    val u1_gauged         : bool
    val anom_ferm_ass     : bool
  end

module type THDM_flags =
  sig
    val ckm_present       : bool
  end

module BSM_bsm : BSM_flags
module BSM_ungauged : BSM_flags
module BSM_anom : BSM_flags
module Littlest : functor (F: BSM_flags) -> Model.Gauge with module Ch = Charges.QQ
module Littlest_Tpar : functor (F: BSM_flags) -> Model.T with module Ch = Charges.QQ
module Simplest : functor (F: BSM_flags) -> Model.T with module Ch = Charges.QQ
module Xdim : functor (F: BSM_flags) -> Model.Gauge with module Ch = Charges.QQ
module UED : functor (F: BSM_flags) -> Model.Gauge with module Ch = Charges.QQ
module GravTest : functor (F: BSM_flags) -> Model.Gauge with module Ch = Charges.QQ
module THDM : THDM_flags
module THDM_CKM : THDM_flags
module TwoHiggsDoublet : functor (F : THDM_flags) -> Model.Gauge with module Ch = Charges.QQ
module Template : functor (F : BSM_flags) -> Model.Gauge with module Ch = Charges.QQ
module HSExt : functor (F : BSM_flags) -> Model.Gauge with module Ch = Charges.QQ

module type Threeshl_options =
  sig
    val include_ckm: bool
    val include_hf: bool
    val diet: bool
  end

module Threeshl_no_ckm: Threeshl_options
module Threeshl_ckm: Threeshl_options
module Threeshl_no_ckm_no_hf: Threeshl_options
module Threeshl_ckm_no_hf: Threeshl_options
module Threeshl_diet_no_hf: Threeshl_options 
module Threeshl_diet: Threeshl_options
module Threeshl: functor (Module_options: Threeshl_options) ->
  Model.T with module Ch = Charges.QQ

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

module SSC_kmatrix: SSC_flags

module SSC_kmatrix_2: SSC_flags

module SSC: functor (F : SSC_flags) -> Model.Gauge with module Ch = Charges.QQ

module SSC_AltT: functor (F : SSC_flags) -> Model.Gauge with module Ch = Charges.QQ
