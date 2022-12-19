(* target.mli --

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

module type T =
  sig
    type amplitudes

    val options : Options.t
    type diagnostic = All | Arguments | Momenta | Gauge

(* Format the amplitudes as a sequence of strings. *)
    val amplitudes_to_channel : string -> out_channel ->
      (diagnostic * bool) list -> amplitudes -> unit

    val parameters_to_channel : out_channel -> unit

  end

module type Maker =
    functor (F : Fusion.Maker) ->
      functor (P : Momentum.T) -> functor (M : Model.T) ->
        T with type amplitudes = Fusion.Multi(F)(P)(M).amplitudes

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
