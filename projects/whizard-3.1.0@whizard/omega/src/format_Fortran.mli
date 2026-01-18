(* format_Fortran.mli -- Fortran90+ continuation lines etc.

   Copyright (C) 2019-2022 by

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

(* Mimic parts of the [Format] API with support for Fortran
   style line continuation. *)

type formatter

val std_formatter : formatter

val fprintf : formatter -> ('a, Format.formatter, unit) format -> 'a
val printf : ('a, Format.formatter, unit) format -> 'a

(* Start a new line, \emph{not} a continuation! *)
val pp_newline : formatter -> unit -> unit
val newline : unit -> unit

val pp_flush : formatter -> unit -> unit
val flush : unit -> unit

val formatter_of_out_channel : ?width:int -> out_channel -> formatter
val formatter_of_buffer : ?width:int -> Buffer.t -> formatter

val pp_set_formatter_out_channel : formatter -> ?width:int -> out_channel -> unit
val set_formatter_out_channel : ?width:int -> out_channel -> unit

(* This must be exposed for the benefit of
   [Targets.Make_Fortran().print_interface],
   because somebody decided to use it for the $K$-matrix
   support. Is this really necessary? *)
val pp_switch_line_continuation : formatter -> bool -> unit
val switch_line_continuation : bool -> unit

module Test : sig val suite : OUnit.test end

