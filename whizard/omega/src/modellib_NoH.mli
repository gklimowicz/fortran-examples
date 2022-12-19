(* modellib_NoH.mli --

   Copyright (C) 1999-2022 by

       Wolfgang Kilian <kilian@physik.uni-siegen.de>
       Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
       Juergen Reuter <juergen.reuter@desy.de>
       with contributions from
       Christian Speckner <cnspeckn@googlemail.com>
       Marco Sekulla <marco.sekulla@kit.edu>

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

module type NoH_flags =
  sig
    val triple_anom : bool
    val quartic_anom : bool
    val k_matrix : bool
    val ckm_present : bool
    val top_anom : bool
    val top_anom_4f : bool
  end

module NoH_k_matrix : NoH_flags

module NoH : functor (F : NoH_flags) -> Model.Gauge with module Ch = Charges.QQ
module AltH: functor (F : NoH_flags) -> Model.Gauge with module Ch = Charges.QQ
