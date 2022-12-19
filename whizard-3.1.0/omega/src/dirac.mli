(* dirac.mli --

   Copyright (C) 1999-2017 by

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

(* \thocwmodulesection{Dirac $\gamma$-matrices} *)

module type T =
  sig

    (* Matrices with complex rational entries. *)
    type qc = Algebra.QC.t
    type t = qc array array

    (* Complex rational constants. *)
    val zero : qc
    val one : qc
    val minus_one : qc
    val i : qc
    val minus_i : qc

    (* Basic $\gamma$-matrices. *)
    val unit : t
    val null : t
    val gamma0 : t
    val gamma1 : t
    val gamma2 : t
    val gamma3 : t
    val gamma5 : t

    (* $(\gamma_0,\gamma_1,\gamma_2,\gamma_3)$ *)
    val gamma : t array

    (* Charge conjugation *)
    val cc : t

    (* Algebraic operations on $\gamma$-matrices *)
    val neg : t -> t
    val add : t -> t -> t
    val sub : t -> t -> t
    val mul : t -> t -> t
    val times : qc -> t -> t
    val transpose : t -> t
    val adjoint : t -> t
    val conj : t -> t
    val product : t list -> t

    (* Toplevel *)
    val pp : Format.formatter -> t -> unit

    (* Unit tests *)
    val test_suite : OUnit.test
  end

module Chiral : T
module Dirac : T
module Majorana : T
