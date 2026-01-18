(* circe2/histogram.mli --  *)
(* Copyright (C) 2001-2022 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
   Circe2 is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by 
   the Free Software Foundation; either version 2, or (at your option)
   any later version.
   Circe2 is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  *)  

type t
val create : int -> float -> float -> t
val record : t -> float -> float -> unit
val normalize : t -> t
val to_channel : out_channel -> t -> unit
val to_file : string -> t -> unit
val as_bins_to_channel : out_channel -> t -> unit
val as_bins_to_file : string -> t -> unit

val regression : t -> (float -> bool) ->
  (float -> float) -> (float -> float) -> float * float

(*i
 *  Local Variables:
 *  mode:caml
 *  indent-tabs-mode:nil
 *  page-delimiter:"^(\\* .*\n"
 *  End:
i*)
