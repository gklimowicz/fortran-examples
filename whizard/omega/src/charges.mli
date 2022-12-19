(* charges.mli --

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

(* \thocwmodulesection{Abstract Type} *)

module type T = 
  sig

    (* The abstract type of the set of conserved charges or
       additive quantum numbers.  *)
    type t

    (* Add the quantum numbers of a pair or a list of particles. *)
    val add : t -> t -> t
    val sum : t list -> t

    (* Test the charge conservation.  *)
    val is_null : t -> bool

  end


(* \thocwmodulesection{Trivial Realisation} *)
module Null : T with type t = unit

(* \thocwmodulesection{Nontrivial Realisations} *)

(* \thocwmodulesubsection{$\mathbf{Z}$} *)
module Z : T with type t = int

(* \thocwmodulesubsection{$\mathbf{Z}\times\mathbf{Z}\times\cdots\times\mathbf{Z}$} *)
module ZZ : T with type t = int list

(* \thocwmodulesubsection{$\mathbf{Q}$} *)
module Q : T with type t = Algebra.Small_Rational.t

(* \thocwmodulesubsection{$\mathbf{Q}\times\mathbf{Q}\times\cdots\times\mathbf{Q}$} *)
module QQ : T with type t = Algebra.Small_Rational.t list

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
