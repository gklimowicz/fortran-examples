(* modellib_PSSSM.mli --

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

(* \thocwmodulesection{Extended Supersymmetric Models} *)

(* We do not introduce the possibility here of using four point couplings 
   or not. We simply add the relevant and leave the rest out. No 
   possibility for Goldstone bosons is given. But we allow for CKM mixing.
*)

module type extMSSM_flags = 
  sig 
    val ckm_present       : bool
  end

module PSSSM : extMSSM_flags
module PSSSM_QCD : extMSSM_flags
module ExtMSSM : functor (F: extMSSM_flags) -> Model.T with module Ch = Charges.QQ

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
